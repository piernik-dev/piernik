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

   function datafields_descr(var) result(f)

      use gdf,    only: gdf_field_type
      use units,  only: cm, gram, sek

      implicit none

      character(len=*), intent(in) :: var
      type(gdf_field_type) :: f

      f%f2cgs = 1.0
      f%stag  = 0
      f%fn    = trim(var)
      f%fu    = 'fixme'
      select case (trim(var))
         case ("dend", "deni", "denn")
            f%fu = "\rm{g}/\rm{cm}^3"
            f%f2cgs = gram/cm**3
         case ("vlxd", "vlxn", "vlxi", "vlyd", "vlyn", "vlyi", "vlzd", "vlzn", "vlzi")
            f%fu = "\rm{cm}/\rm{s}"
            f%f2cgs = cm/sek
         case ("enen", "enei")
            f%fu = "\rm{g}*\rm{cm}^2/\rm{s}^2"
            f%f2cgs = gram*cm**2/sek**2
         case ("pren", "prei")
            f%fu = "\rm{g}/\rm{cm}/\rm{s}^2"
            f%f2cgs = gram/cm/sek**2
         case ("magx", "magy", "magz")
            f%stag = 1
         case ("cr1" : "cr9")
         case ("mgso")
         case ("gpot")
         case ("trcr")
      end select
   end function datafields_descr

   subroutine create_datafields_descrs(place)

      use common_hdf5,  only: hdf_vars
      use gdf,          only: gdf_field_type, fmax
      use hdf5,         only: HID_T, h5gcreate_f, h5gclose_f
      use helpers_hdf5, only: create_attribute

      implicit none

      integer(HID_T), intent(in) :: place

      integer         :: i
      integer(kind=4) :: error
      integer(HID_T)  :: g_id
      type(gdf_field_type), target :: f
      integer(kind=8), pointer, dimension(:) :: ibuf
      character(len=fmax), pointer :: sbuf

      allocate(ibuf(1))
      do i = lbound(hdf_vars,1), ubound(hdf_vars,1)
         f = datafields_descr(hdf_vars(i))
         call h5gcreate_f(place, hdf_vars(i), g_id, error)
         call create_attribute(g_id, 'field_to_cgs', [f%f2cgs])
         ibuf = f%stag
         call create_attribute(g_id, 'staggering',   ibuf)
         sbuf => f%fu
         call create_attribute(g_id, 'field_units',  sbuf)
         sbuf => f%fn
         call create_attribute(g_id, 'field_name',   sbuf)
         call h5gclose_f(g_id, error)
      enddo
      deallocate(ibuf)
   end subroutine create_datafields_descrs
!>
!! \brief Routine calculating quantities for .hdf files
!<
   subroutine datafields_hdf5(var, tab, ierrh, cg)

      use common_hdf5, only: common_shortcuts
      use constants,   only: varlen, xdim
      use fluidtypes,  only: component_fluid
      use func,        only: ekin, emag
      use grid_cont,   only: grid_container
#ifndef ISO
      use constants,   only: ydim, zdim
#endif /* !ISO */
#if defined(COSM_RAYS) || defined(TRACER) || !defined(ISO)
      use fluidindex,  only: flind
#endif /* COSM_RAYS || TRACER || !ISO */

      implicit none

      character(len=varlen),          intent(in)  :: var
      real(kind=4), dimension(:,:,:)              :: tab
      integer,                        intent(out) :: ierrh
      type(grid_container),  pointer, intent(in)  :: cg
      class(component_fluid), pointer              :: fl_dni
      integer(kind=4)                             :: i_xyz
#ifdef COSM_RAYS
      integer                                     :: i
      integer, parameter                          :: auxlen = varlen - 1
      character(len=auxlen)                       :: aux
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
#ifdef TRACER
         case ("trcr")
            tab(:,:,:) = real(cg%u(flind%trc%beg, RNG),4)
#endif /* TRACER */
         case ("dend", "deni", "denn")
            tab(:,:,:) = real(cg%u(fl_dni%idn, RNG), kind=4)
         case ("vlxd", "vlxn", "vlxi", "vlyd", "vlyn", "vlyi", "vlzd", "vlzn", "vlzi")
            tab(:,:,:) = real(cg%u(fl_dni%imx + i_xyz, RNG) / cg%u(fl_dni%idn, RNG), kind=4)
         case ("enen", "enei")
#ifdef ISO
            tab(:,:,:) = real(ekin(cg%u(fl_dni%imx, RNG), cg%u(fl_dni%imy, RNG), cg%u(fl_dni%imz, RNG), cg%u(fl_dni%idn, RNG)), kind=4)
#else /* !ISO */
            tab(:,:,:) = real(cg%u(fl_dni%ien, RNG), kind=4)
#endif /* !ISO */
         case ("pren")
#ifndef ISO
            tab(:,:,:) = real(flind%neu%gam_1, kind=4) * real( cg%u(flind%neu%ien, RNG) - &
                 &       ekin(cg%u(flind%neu%imx, RNG), cg%u(flind%neu%imy, RNG), cg%u(flind%neu%imz, RNG), cg%u(flind%neu%idn, RNG)), kind=4)
#endif /* !ISO */
         case ("prei")
#ifndef ISO
            tab(:,:,:) = real(flind%ion%gam_1, kind=4) * real( cg%u(flind%ion%ien, RNG) - &
                 &       ekin(cg%u(flind%ion%imx, RNG), cg%u(flind%ion%imy, RNG), cg%u(flind%ion%imz, RNG), cg%u(flind%ion%idn, RNG)), kind=4) - &
                 &       real(flind%ion%gam_1*emag(cg%b(xdim, RNG), cg%b(ydim, RNG), cg%b(zdim, RNG)), kind=4)
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

      use mpi,         only: MPI_CHARACTER
      use common_hdf5, only: set_common_attributes
      use constants,   only: cwdlen, I_ONE
      use dataio_pub,  only: printio, printinfo, nhdf, thdf, tmr_hdf, wd_wr, piernik_hdf5_version, piernik_hdf5_version2, &
         & msg, run_id, problem_name, use_v2_io, last_hdf_time
      use mpisetup,    only: comm, mpi_err, master, FIRST
      use timer,       only: set_timer

      implicit none

      character(len=cwdlen)             :: fname
      real :: phv

      thdf = set_timer(tmr_hdf,.true.)
      nhdf = nhdf + I_ONE
      ! Initialize HDF5 library and Fortran interfaces.
      !
      phv = piernik_hdf5_version ; if (use_v2_io) phv = piernik_hdf5_version2
      if (master) then
         write(fname, '(2a,a1,a3,a1,i4.4,a3)') trim(wd_wr), trim(problem_name),"_", trim(run_id),"_", nhdf,".h5" !> \todo: merge with function restart_fname()
         write(msg,'(a,es23.16,a,f5.2,1x,2a)') 'ordered t ',last_hdf_time,': Writing datafile v', phv, trim(fname), " ... "
         call printio(msg, .true.)
      endif
      call MPI_Bcast(fname, cwdlen, MPI_CHARACTER, FIRST, comm, mpi_err)

      call set_common_attributes(fname)
      if (use_v2_io) then
         call h5_write_to_single_file_v2(fname)
      else
         call h5_write_to_single_file_v1(fname)
      endif

      thdf = set_timer(tmr_hdf)
      if (master) then
         write(msg,'(a6,f10.2,a2)') ' done ', thdf, ' s'
         call printinfo(msg, .true.)
      endif

   end subroutine h5_write_to_single_file

   subroutine h5_write_to_single_file_v2(fname)
      use common_hdf5, only: write_to_hdf5_v2, O_OUT
      use gdf,         only: gdf_create_field_types
      use mpisetup,    only: comm, mpi_err, master

      implicit none

      character(len=*), intent(in)      :: fname

      call write_to_hdf5_v2(fname, O_OUT, create_empty_cg_datasets_in_output, write_cg_to_output)

      if (master) call gdf_create_field_types(fname,create_datafields_descrs)
      call MPI_Barrier(comm, mpi_err)

   end subroutine h5_write_to_single_file_v2

!> \brief Write all grid containers to the file

   subroutine write_cg_to_output(cgl_g_id, cg_n, cg_all_n_b)

      use constants,   only: xdim, ydim, zdim, ndims
      use common_hdf5, only: get_nth_cg, hdf_vars, cg_output, hdf_vars
      use dataio_pub,  only: die, nproc_io, can_i_write
      use cg_list,     only: cg_list_element
      use grid_cont,   only: grid_container
      use grid,        only: leaves
      use hdf5,        only: HID_T, HSIZE_T, H5T_NATIVE_REAL, h5sclose_f, h5dwrite_f, h5sselect_none_f, h5screate_simple_f
      use mpi,         only: MPI_REAL, MPI_STATUS_IGNORE
      use mpisetup,    only: master, FIRST, proc, comm, mpi_err

      implicit none

      integer(HID_T),                           intent(in) :: cgl_g_id    !> cg group identifier
      integer(kind=4), dimension(:),   pointer, intent(in) :: cg_n        !> offset for cg group numbering
      integer(kind=4), dimension(:,:), pointer, intent(in) :: cg_all_n_b  !> all cg sizes

      integer(HID_T)                              :: filespace_id, memspace_id
      integer(kind=4)                             :: error
      integer(kind=4), parameter        :: rank = 3
      integer(HSIZE_T), dimension(:), allocatable :: dims
      integer                                     :: i, ncg, n
      type(grid_container), pointer               :: cg
      type(cg_list_element), pointer              :: cgl
      real(kind=4), dimension(:,:,:), pointer     :: data
      type(cg_output) :: cg_desc

      call cg_desc%init(cgl_g_id, cg_n, nproc_io, hdf_vars)

      if (cg_desc%tot_cg_n < 1) call die("[data_hdf5:write_cg_to_output] no cg available!")

      ! all arrays are rank 3 here
      allocate(dims(ndims))
      ! Allocate data with the size of first cg
      allocate( data(cg_all_n_b(xdim, 1), cg_all_n_b(ydim, 1), cg_all_n_b(zdim, 1)) )

      if (nproc_io == 1) then ! perform serial write
         ! write all cg, one by one
         do ncg = 1, cg_desc%tot_cg_n
            dims = [ cg_all_n_b(xdim, ncg), cg_all_n_b(ydim, ncg), cg_all_n_b(zdim, ncg) ]
            call recycle_data(dims, cg_all_n_b, ncg, data)

            if (master) then
               if (.not. can_i_write) call die("[data_hdf5:write_cg_to_output] Master can't write")

               do i = lbound(hdf_vars,1), ubound(hdf_vars,1)
                  if (cg_desc%cg_src_p(ncg) == proc) then
                     cg => get_nth_cg(cg_desc%cg_src_n(ncg))
                     call get_data_from_cg(hdf_vars(i), cg, data)
                  else
                     call MPI_Recv(data(1,1,1), size(data), MPI_REAL, cg_desc%cg_src_p(ncg), ncg + cg_desc%tot_cg_n*i, comm, MPI_STATUS_IGNORE, mpi_err)
                  endif
                  call h5dwrite_f(cg_desc%dset_id(ncg, i), H5T_NATIVE_REAL, data, dims, error, xfer_prp = cg_desc%xfer_prp)
               enddo
            else
               if (can_i_write) call die("[data_hdf5:write_cg_to_output] Slave can write")
               if (cg_desc%cg_src_p(ncg) == proc) then
                  cg => get_nth_cg(cg_desc%cg_src_n(ncg))
                  do i = lbound(hdf_vars,1), ubound(hdf_vars,1)
                     call get_data_from_cg(hdf_vars(i), cg, data)
                     call MPI_Send(data(1,1,1), size(data), MPI_REAL, FIRST, ncg + cg_desc%tot_cg_n*i, comm, mpi_err)
                  enddo
               endif
            endif
         enddo
      else ! perform parallell write
         ! This piece will be a generalization of the serial case. It should work correctly also for nproc_io == 1 so it should replace the serial code
         if (can_i_write) then
            ! write own
            n = 0
            cgl => leaves%first

            do while (associated(cgl))
               n = n + 1
               ncg = cg_desc%offsets(proc) + n
               dims = [ cg_all_n_b(xdim, ncg), cg_all_n_b(ydim, ncg), cg_all_n_b(zdim, ncg) ]
               call recycle_data(dims, cg_all_n_b, ncg, data)
               cg => cgl%cg

               do i = lbound(hdf_vars,1), ubound(hdf_vars,1)
                  call get_data_from_cg(hdf_vars(i), cg, data)
                  call h5dwrite_f(cg_desc%dset_id(ncg, i), H5T_NATIVE_REAL, data, dims, error, xfer_prp = cg_desc%xfer_prp)
               enddo

               cgl => cgl%nxt
            enddo

            ! Behold the MAGIC in its purest form!
            ! Following block of code does exactly *nothing*, yet it's necessary for collective calls of PHDF5
            !>
            !! \deprecated BEWARE, we assume that at least 1 cg exist on a given proc
            !! \todo there should be something like H5S_NONE as a contradiction to H5S_ALL, yet I cannot find it...
            !<

            dims = [ cg_all_n_b(xdim, 1), cg_all_n_b(ydim, 1), cg_all_n_b(zdim, 1) ]

            ! completely bogus values, only to make HDF5 happy
            call h5screate_simple_f(rank, dims, filespace_id, error)
            call h5screate_simple_f(rank, dims, memspace_id, error)
            ! empty filespace
            call h5sselect_none_f(filespace_id, error)
            ! empty memoryspace
            call h5sselect_none_f(memspace_id, error)

            call recycle_data(dims, cg_all_n_b, 1, data)
            do ncg = 1, maxval(cg_n)-n
               do i = lbound(hdf_vars, 1), ubound(hdf_vars, 1)
                  call h5dwrite_f(cg_desc%dset_id(1, i), H5T_NATIVE_REAL, data, dims, error, &
                     xfer_prp = cg_desc%xfer_prp, file_space_id = filespace_id, mem_space_id = memspace_id)
               enddo
            enddo

            call h5sclose_f(memspace_id, error)
            call h5sclose_f(filespace_id, error)
            ! receive (from whom?)
         else
            call die("[data_hdf5:write_cg_to_output] nproc != nproc_io not implemented yet")
            ! send (where?)
         endif
         !call die("[data_hdf5:write_cg_to_output] Parallel v2 I/O not implemented yet")
      endif

      ! clean up
      if (allocated(dims)) deallocate(dims)
      if (associated(data)) deallocate(data)
      call cg_desc%clean()

      contains
         ! Try to avoid pointless data reallocation for every cg if shape doesn't change
         subroutine recycle_data(dims, cg_all_n_b, i, data)
            use constants, only: xdim, ydim, zdim
            use hdf5,      only: HSIZE_T
            implicit none
            integer(HSIZE_T), dimension(:) :: dims
            integer(kind=4), dimension(:,:), pointer, intent(in) :: cg_all_n_b  !> all cg sizes
            integer, intent(in) :: i
            real(kind=4), dimension(:,:,:), pointer :: data

            if (associated(data)) then
               if ( any(dims /= shape(data)) ) then
                  deallocate(data)
                  allocate(data(cg_all_n_b(xdim, i), cg_all_n_b(ydim, i), cg_all_n_b(zdim, i)))
               endif
            endif
         end subroutine recycle_data

   end subroutine write_cg_to_output

   subroutine get_data_from_cg(hdf_var, cg, tab)

      use dataio_pub,  only: die, msg
      use dataio_user, only: user_vars_hdf5
      use grid_cont,   only: grid_container

      implicit none

      character(len=*), intent(in)                           :: hdf_var
      type(grid_container), pointer, intent(in)              :: cg
      real(kind=4), dimension(:,:,:), pointer, intent(inout) :: tab

      integer :: ierrh
      logical :: ok_var

      ierrh = 0
      ok_var = .false.
      call datafields_hdf5(hdf_var, tab, ierrh, cg)
      if (associated(user_vars_hdf5) .and. ierrh /= 0) call user_vars_hdf5(hdf_var, tab, ierrh, cg)
      if (ierrh>=0) ok_var = .true.
      if (.not.ok_var) then
         write(msg,'(3a)') "[data_hdf5:get_data_from_cg]: ", hdf_var," is not defined in datafields_hdf5, neither in user_vars_hdf5."
         call die(msg)
      endif
   end subroutine get_data_from_cg

   subroutine create_empty_cg_datasets_in_output(cg_g_id, cg_n_b, Z_avail, g)

      use common_hdf5, only: create_empty_cg_dataset, hdf_vars
      use hdf5,        only: HID_T, HSIZE_T

      implicit none

      integer(HID_T), intent(in)                           :: cg_g_id
      integer(kind=4), dimension(:,:), pointer, intent(in) :: cg_n_b
      logical(kind=4), intent(in)                          :: Z_avail
      integer, intent(in)                                  :: g

      integer :: i

      do i = lbound(hdf_vars,1), ubound(hdf_vars,1)
         call create_empty_cg_dataset(cg_g_id, hdf_vars(i), int(cg_n_b(g, :), kind=HSIZE_T), Z_avail)
      enddo
   end subroutine create_empty_cg_datasets_in_output

   subroutine h5_write_to_single_file_v1(fname)

      use common_hdf5, only: nhdf_vars, hdf_vars
      use constants,   only: ndims
      use domain,      only: is_multicg !, is_uneven
      use grid,        only: finest
      use cg_list,     only: cg_list_element
      use grid_cont,   only: grid_container
      use hdf5,        only: HID_T, HSIZE_T, H5FD_MPIO_COLLECTIVE_F, H5P_DATASET_CREATE_F, H5P_DATASET_XFER_F, &
           &                 H5S_SELECT_SET_F, H5T_NATIVE_REAL, H5F_ACC_RDWR_F, H5P_FILE_ACCESS_F, &
           &                 h5dwrite_f, h5screate_simple_f, h5pcreate_f, h5dcreate_f, h5sclose_f, h5dget_space_f, h5sselect_hyperslab_f, &
           &                 h5pset_dxpl_mpio_f, h5dclose_f, h5open_f, h5close_f, h5fopen_f, h5fclose_f, h5pclose_f, h5pset_fapl_mpio_f !, h5pset_chunk_f
      use mpisetup,    only: comm
      use mpi,         only: MPI_INFO_NULL

      implicit none

      character(len=*), intent(in)      :: fname
      integer(HID_T)                    :: file_id                 !< File identifier
      integer(HID_T)                    :: plist_id, plist_idf     !< Property list identifier
      integer                           :: i
      integer(kind=4)                   :: error
      type(cg_list_element), pointer    :: cgl
      type(grid_container),  pointer    :: cg
      real(kind=4), pointer             :: data (:,:,:)            !< Data to write
      integer(kind=4), parameter        :: rank = ndims            !< Dataset rank = 3
      integer(HID_T)                    :: dset_id                 !< Dataset identifier
      integer(HID_T)                    :: filespace               !< Dataspace identifier in file
      integer(HID_T)                    :: memspace                !< Dataspace identifier in memory
      integer(HSIZE_T), dimension(rank) :: count, offset, stride, block, dimsf, chunk_dims

      call h5open_f(error)
      !
      ! Setup file access property list with parallel I/O access.
      !
      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_idf, error)
      call h5pset_fapl_mpio_f(plist_idf, comm, MPI_INFO_NULL, error)
      !
      ! Create the file collectively.
      !
      call h5fopen_f(trim(fname), H5F_ACC_RDWR_F, file_id, error, access_prp = plist_idf)
      call h5pclose_f(plist_idf, error)

      !! \todo check if finest is complete, if not then find finest complete level
      dimsf  = finest%n_d(:)    ! Dataset dimensions
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

         !! \todo if there are fine levels restrict the data first and write warning that v2 should be used instead
         cgl => finest%first
         do while (associated(cgl))
            cg => cgl%cg

            if (.not.associated(data)) allocate(data(cg%nxb, cg%nyb, cg%nzb))
            call get_data_from_cg(hdf_vars(i), cg, data)

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

            if (associated(data)) deallocate(data)

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

   end subroutine h5_write_to_single_file_v1

   function h5_filename() result(f)
      use constants,  only: fnamelen
      use dataio_pub, only: problem_name, run_id, nhdf, wd_wr
      use mpisetup,   only: proc
      implicit none
      character(len=fnamelen) :: f
      write(f, '(2a,"_",a3,i4.4,".cpu",i5.5,".h5")') trim(wd_wr), trim(problem_name), trim(run_id), nhdf, proc
   end function h5_filename

   subroutine h5_write_to_multiple_files

      use constants,       only: dsetnamelen, fnamelen, xdim, ydim, zdim, I_ONE
      use common_hdf5,     only: nhdf_vars, hdf_vars
      use dataio_pub,      only: msg, printio, printinfo, tmr_hdf, thdf, last_hdf_time, piernik_hdf5_version
      use cg_list,         only: cg_list_element
      use grid,            only: leaves
      use grid_cont,       only: grid_container
      use h5lt,            only: h5ltmake_dataset_float_f, h5ltmake_dataset_double_f
      use hdf5,            only: H5F_ACC_TRUNC_F, h5fcreate_f, h5open_f, h5fclose_f, h5close_f, HID_T, h5gcreate_f, &
           &                     h5gclose_f, HSIZE_T
      use mpisetup,        only: master
      use timer,           only: set_timer

      implicit none

      type(cg_list_element), pointer    :: cgl
      type(grid_container),  pointer    :: cg
      integer(kind=4), parameter        :: rank = 3
      integer(kind=4)                   :: error, i
      integer(HID_T)                    :: file_id, grp_id
      integer(kind=8)                   :: ngc              !> current grid index
      integer(HSIZE_T), dimension(rank) :: dims
      character(len=dsetnamelen)        :: gname
      character(len=fnamelen)           :: fname
      real(kind=4), pointer             :: data (:,:,:)     !> Data to write

      thdf = set_timer(tmr_hdf,.true.)
      fname = h5_filename()
      if (master) then
         write(msg,'(a,es23.16,a,f5.2,1x,2a)') 'ordered t ',last_hdf_time,': Writing datafile v', piernik_hdf5_version, trim(fname), " ... "
         call printio(msg, .true.)
      endif

      call h5open_f(error)
      call h5fcreate_f(fname, H5F_ACC_TRUNC_F, file_id, error)
      cgl => leaves%first
      ngc = 0
      do while (associated(cgl))
         cg => cgl%cg

         write(gname,'("grid",i8.8)') ngc-1
         call h5gcreate_f(file_id, gname, grp_id, error)

         ! set attributes here
         call h5ltmake_dataset_double_f(grp_id, "fbnd", int(2,kind=4), shape(cg%fbnd,kind=HSIZE_T), &
                                      & cg%fbnd, error)

         if (.not.associated(data)) allocate(data(cg%n_b(xdim),cg%n_b(ydim),cg%n_b(zdim)))
         dims = cg%n_b(:)
         do i = I_ONE, int(nhdf_vars, kind=4)
            call get_data_from_cg(hdf_vars(i), cg, data)
            call h5ltmake_dataset_float_f(grp_id, hdf_vars(i), rank, dims, data(:,:,:), error)
         enddo
         if (associated(data)) deallocate(data)
         call h5gclose_f(grp_id, error)
         ngc = ngc + 1
         cgl => cgl%nxt
      enddo
      call h5fclose_f(file_id, error)
      call h5close_f(error)

      thdf = set_timer(tmr_hdf)
      if (master) then
         write(msg,'(a6,f10.2,a2)') ' done ', thdf, ' s'
         call printinfo(msg, .true.)
      endif

   end subroutine h5_write_to_multiple_files

end module data_hdf5
