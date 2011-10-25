! $Id$
!
! PIERNIK Code Copyright (C) 2006 Michal Hanasz
!
!    This file is part of PIERNIK code.
!
!    PIERNIK is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    PIERNIK is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with PIERNIK.  If not, see <http://www.gnu.org/licenses/>.
!
!    Initial implementation of PIERNIK code was based on TVD split MHD code by
!    Ue-Li Pen
!        see: Pen, Arras & Wong (2003) for algorithm and
!             http://www.cita.utoronto.ca/~pen/MHD
!             for original source code "mhd.f90"
!
!    For full list of developers see $PIERNIK_HOME/license/pdt.txt
!

#include "piernik.h"

!>
!! \brief Module that contains HDF5 I/O routines for writing single-precision data dumps
!<
module data_hdf5

! pulled by ANY

   implicit none

   private
   public :: init_data, write_hdf5

   interface
      subroutine h5_write
         implicit none
      end subroutine h5_write
   end interface

   procedure(h5_write), pointer :: write_hdf5 => NULL() !h5_write_to_single_file

contains

   subroutine init_data

      use dataio_pub,  only: multiple_h5files

      implicit none

      if (multiple_h5files) then
         write_hdf5 => h5_write_to_multiple_files
      else
         write_hdf5 => h5_write_to_single_file
      endif


   end subroutine init_data

!>
!! \brief Routine calculating quantities for .hdf files
!<
   subroutine datafields_hdf5(var, tab, ierrh, cg)

      use common_hdf5, only: common_shortcuts
      use constants,   only: varlen, half, xdim, ydim, zdim
      use fluidindex,  only: flind
      use fluidtypes,  only: component_fluid
      use grid_cont,   only: grid_container

      implicit none

      character(len=varlen), intent(in) :: var
      real(kind=4), dimension(:,:,:)    :: tab
      integer, intent(out)              :: ierrh
      type(grid_container), pointer, intent(in) :: cg
      type(component_fluid), pointer :: fl_dni
      integer :: i_xyz
#ifdef COSM_RAYS
      integer :: i
      integer, parameter    :: auxlen = varlen - 1
      character(len=auxlen) :: aux
#endif /* COSM_RAYS */
#define RNG cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke

      call common_shortcuts(var, fl_dni, i_xyz)

      ierrh = 0
      tab = 0.0

      select case (var)
#ifdef COSM_RAYS
         case ("cr1" : "cr9")
            read(var,'(A2,I1)') aux, i !> \deprecated BEWARE 0 <= i <= 9, no other indices can be dumped to hdf file
            tab(:,:,:) = real(cg%u(flind%crs%beg+i-1, RNG), kind=4)
#endif /* COSM_RAYS */
         case ("dend", "deni", "denn")
            tab(:,:,:) = real(cg%u(fl_dni%idn, RNG), kind=4)
         case ("vlxd", "vlxn", "vlxi", "vlyd", "vlyn", "vlyi", "vlzd", "vlzn", "vlzi")
            tab(:,:,:) = real(cg%u(fl_dni%imx + i_xyz, RNG) / cg%u(fl_dni%idn, RNG), kind=4)
         case ("enen", "enei")
#ifdef ISO
            tab(:,:,:) = real(half *( cg%u(fl_dni%imx, RNG)**2 + &
                 &                   cg%u(fl_dni%imy, RNG)**2 + &
                 &                   cg%u(fl_dni%imz, RNG)**2 ) / cg%u(fl_dni%idn, RNG), kind=4)
#else /* !ISO */
            tab(:,:,:) = real(cg%u(fl_dni%ien, RNG), kind=4)
#endif /* !ISO */
#ifdef NEUTRAL
         case ("pren")
#ifndef ISO
            tab(:,:,:) = real( cg%u(flind%neu%ien, RNG) - half * ( &
                 &             cg%u(flind%neu%imx, RNG)**2 + &
                 &             cg%u(flind%neu%imy, RNG)**2 + &
                 &             cg%u(flind%neu%imz, RNG)**2 ) / cg%u(flind%neu%idn, RNG), kind=4) * real(flind%neu%gam_1, kind=4)
#endif /* !ISO */
#endif /* NEUTRAL */
         case ("prei")
#ifndef ISO
            tab(:,:,:) = real( cg%u(flind%ion%ien, RNG) - half *( &
                 &             cg%u(flind%ion%imx, RNG)**2 + &
                 &             cg%u(flind%ion%imy, RNG)**2 + &
                 &             cg%u(flind%ion%imz, RNG)**2 ) / cg%u(flind%ion%idn, RNG), kind=4) * real(flind%ion%gam_1, kind=4) - &
                 &       real( half*(flind%ion%gam_1)*(cg%b(xdim, RNG)**2 + cg%b(ydim, RNG)**2 + cg%b(zdim, RNG)**2), kind=4)
#endif /* !ISO */
         case ("magx", "magy", "magz")
            tab(:,:,:) = real(cg%b(xdim + i_xyz, RNG), kind=4)
         case ("gpot")
            if (associated(cg%gpot)) tab(:,:,:) = real(cg%gpot(RNG), kind=4)
         case ("mgso")
            if (associated(cg%sgp))  tab(:,:,:) = real(cg%sgp(RNG),  kind=4)
         case default
            ierrh = -1
      end select

#undef RNG

   end subroutine datafields_hdf5

!
! ------------------------------------------------------------------------------------
!
   subroutine h5_write_to_single_file

      use common_hdf5, only: nhdf_vars, hdf_vars, set_common_attributes
      use constants,   only: ndims, cwdlen, I_ONE
      use dataio_pub,  only: printio, msg, die, nhdf, problem_name, run_id, nhdf
      use dataio_user, only: user_vars_hdf5
      use domain,      only: dom, is_multicg !, is_uneven
      use grid,        only: all_cg
      use gc_list,     only: cg_list_element
      use grid_cont,   only: grid_container
      use hdf5,        only: HID_T, HSIZE_T, H5FD_MPIO_COLLECTIVE_F, H5P_DATASET_CREATE_F, H5P_DATASET_XFER_F, &
           &                 H5S_SELECT_SET_F, H5T_NATIVE_REAL, H5F_ACC_TRUNC_F, H5P_FILE_ACCESS_F, H5P_DEFAULT_F, &
           &                 h5dwrite_f, h5screate_simple_f, h5pcreate_f, h5dcreate_f, h5sclose_f, h5dget_space_f, h5sselect_hyperslab_f, &
           &                 h5pset_dxpl_mpio_f, h5dclose_f, h5open_f, h5close_f, h5fcreate_f, h5fclose_f, h5pclose_f, h5pset_fapl_mpio_f !, h5pset_chunk_f
      use mpisetup,    only: comm, ierr, master, FIRST
      use mpi,         only: MPI_CHARACTER, MPI_INFO_NULL

      implicit none

      integer(HID_T)          :: file_id       ! File identifier
      integer(HID_T)          :: plist_id, plist_idf ! Property list identifier
      integer                 :: ierrh, i
      integer(kind=4)         :: error
      logical                 :: ok_var
      character(len=cwdlen)   :: fname
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg
      real(kind=4), allocatable :: data (:,:,:)  ! Data to write
      integer(kind=4), parameter :: rank = ndims        !< Dataset rank = 3
      integer(HID_T) :: dset_id                 !< Dataset identifier
      integer(HID_T) :: filespace               !< Dataspace identifier in file
      integer(HID_T) :: memspace                !< Dataspace identifier in memory
      integer(HSIZE_T), dimension(rank) :: count, offset, stride, block, dimsf, chunk_dims

      ! Initialize HDF5 library and Fortran interfaces.
      !
      if (master) then
         write(fname, '(a,a1,a3,a1,i4.4,a3)') trim(problem_name),"_", trim(run_id),"_", nhdf,".h5"
         write(msg,'(3a)') 'Writing datafile ', trim(fname), " ... "
         call printio(msg, .true.)
      endif
      call MPI_Bcast(fname, cwdlen, MPI_CHARACTER, FIRST, comm, ierr)

      call h5open_f(error)
      !
      ! Setup file access property list with parallel I/O access.
      !
      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_idf, error)
      call h5pset_fapl_mpio_f(plist_idf, comm, MPI_INFO_NULL, error)
      !
      ! Create the file collectively.
      !
      call h5fcreate_f(trim(fname), H5F_ACC_TRUNC_F, file_id, error, creation_prp = H5P_DEFAULT_F, access_prp = plist_idf)
      call h5pclose_f(plist_idf, error)

      dimsf  = dom%n_d(:)    ! Dataset dimensions
      !
      ! Create the data space for the  dataset.
      !
      call h5screate_simple_f(rank, dimsf, filespace, error)

      do i = 1, nhdf_vars

         ! Create chunked dataset.
         call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)

         ! Cannot use in multiblock
         ! if (.not. is_uneven) call h5pset_chunk_f(plist_id, rank, chunk_dims, error) !> \todo check how much performance it gives (massively parallel I/O is required)

         call h5dcreate_f(file_id, hdf_vars(i), H5T_NATIVE_REAL, filespace, dset_id, error, plist_id)
         call h5sclose_f(filespace, error)

         call h5dget_space_f(dset_id, filespace, error)

         ! Create property list for collective dataset write
         call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
         if (.not. is_multicg) call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

         cgl => all_cg%first
         do while (associated(cgl))
            cg => cgl%cg

            if (.not.allocated(data)) allocate(data(cg%nxb, cg%nyb, cg%nzb))
            ierrh = 0; ok_var = .false.
            call datafields_hdf5(hdf_vars(i), data, ierrh, cg)
            if (associated(user_vars_hdf5) .and. ierrh /= 0) call user_vars_hdf5(hdf_vars(i), data, ierrh, cg)
            if (ierrh>=0) ok_var = .true.
            if (.not.ok_var) then
               write(msg,'(3a)') "[data_hdf5:h5_write_to_single_file]: ", hdf_vars(i)," is not defined in datafields_hdf5, neither in user_vars_hdf5."
               call die(msg)
            endif

            chunk_dims = cg%n_b(:) ! Chunks dimensions
            call h5screate_simple_f(rank, chunk_dims, memspace, error)

            ! Each process defines dataset in memory and writes it to the hyperslab in the file.
            stride(:) = 1
            count(:) =  1
            block(:) = chunk_dims(:)
            offset(:) = cg%off(:)

            ! Select hyperslab in the file.
            call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error, stride, block)

            ! Write the dataset collectively.
            call h5dwrite_f(dset_id, H5T_NATIVE_REAL, data, dimsf, error, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

            ! Close dataspaces.
            call h5sclose_f(memspace, error)

            if (allocated(data)) deallocate(data)

            cgl => cgl%nxt
         enddo

         !call h5pclose_f(plist_id, error) ! does it matter?

         ! Close the dataset.
         call h5dclose_f(dset_id, error)
      enddo

      call h5sclose_f(filespace, error)

      ! Close the property list.
      call h5fclose_f(file_id, error)
      call h5close_f(error)

      call set_common_attributes(fname)

      nhdf = nhdf + I_ONE

   end subroutine h5_write_to_single_file

   function h5_filename() result(f)
      use constants,  only: fnamelen
      use dataio_pub, only: problem_name, run_id, nhdf
      use mpisetup,   only: proc
      implicit none
      character(len=fnamelen) :: f
      write(f, '(a,"_",a3,i4.4,".cpu",i5.5,".h5")') trim(problem_name), trim(run_id), nhdf, proc
   end function h5_filename

   subroutine h5_write_to_multiple_files

      use constants,       only: dsetnamelen, fnamelen, xdim, ydim, zdim, I_ONE
      use common_hdf5,     only: nhdf_vars, hdf_vars
      use dataio_pub,      only: die, msg, printio
      use dataio_user,     only: user_vars_hdf5
      use gc_list,         only: cg_list_element
      use grid,            only: all_cg
      use grid_cont,       only: grid_container
      use h5lt,            only: h5ltmake_dataset_float_f, h5ltmake_dataset_double_f
      use hdf5,            only: H5F_ACC_TRUNC_F, h5fcreate_f, h5open_f, h5fclose_f, h5close_f, HID_T, h5gcreate_f, &
           &                     h5gclose_f, HSIZE_T
      use mpisetup,        only: master

      implicit none

      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg
      integer(kind=4), parameter :: rank = 3
      integer(kind=4) :: error, i
      integer :: error8
      integer(HID_T) :: file_id, grp_id
      integer(kind=8) :: ngc           !> current grid index
      integer(HSIZE_T), dimension(rank) :: dims
      character(len=dsetnamelen) :: gname
      character(len=fnamelen) :: fname
      logical :: ok_var
      real(kind=4), allocatable :: data (:,:,:)  ! Data to write

      fname = h5_filename()
      if (master) then
         write(msg,'(3a)') 'Writing datafile ', trim(fname), " ... "
         call printio(msg, .true.)
      endif

      call h5open_f(error)
      call h5fcreate_f(fname, H5F_ACC_TRUNC_F, file_id, error)
      cgl => all_cg%first
      ngc = 0
      do while (associated(cgl))
         cg => cgl%cg

         write(gname,'("grid",i8.8)') ngc
         call h5gcreate_f(file_id, gname, grp_id, error)

         ! set attributes here
         call h5ltmake_dataset_double_f(grp_id, "fbnd", int(2,kind=4), [integer(kind=HSIZE_T):: shape(cg%fbnd)], &
                                      & cg%fbnd, error)

         if (.not.allocated(data)) allocate(data(cg%n_b(xdim),cg%n_b(ydim),cg%n_b(zdim)))
         dims = cg%n_b(:)
         do i = I_ONE, int(nhdf_vars, kind=4)
            error = 0; ok_var = .false.
            call datafields_hdf5(hdf_vars(i), data, error8, cg)
            if (associated(user_vars_hdf5) .and. error8 /= 0) call user_vars_hdf5(hdf_vars(i), data, error8, cg)
            if (error8>=0) ok_var = .true.
            if (.not.ok_var) then
               write(msg,'(3a)') "[data_hdf5:h5_write_to_multiple_files]: Neither datafields_hdf5", &
                               & " nor user_vars_hdf5 defines ", hdf_vars(i)
               call die(msg)
            endif
            call h5ltmake_dataset_float_f(grp_id, hdf_vars(i), rank, dims, data(:,:,:), error)
         enddo
         if (allocated(data)) deallocate(data)
         call h5gclose_f(grp_id, error)
         ngc = ngc + 1
         cgl => cgl%nxt
      enddo
      call h5fclose_f(file_id, error)
      call h5close_f(error)

      if (master) then
         write(msg,'(a)') 'done'
         call printio(msg)
      endif

   end subroutine h5_write_to_multiple_files

end module data_hdf5
