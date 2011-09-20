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

#define RNG cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke
#include "piernik.h"

!>
!! \brief (KK) Module that contains I/O routines using HDF5 library
!!
!! Modules contains routines for creating HDF5 output such as
!! plots, snapshots, restart files.
!!
!<
module restart_hdf5

! pulled by ANY

   use constants, only: FLUID, MAG, dsetnamelen

   implicit none

   private
   public :: read_restart_hdf5, write_restart_hdf5, write_arr_to_restart, read_arr_from_restart

   !< dataset names for restart files
   enum, bind(C)
      enumerator :: WA=MAG+1, GPOT, HGPOT, GP, SGP, SGPM, CS_ISO2, B0, U0
   end enum
   character(len=dsetnamelen), dimension(FLUID:U0), parameter :: dname = [ "fluid  ", "mag    ", "wa     ", "gpot   ", "hgpot  ", "gp     ", "sgp    ", &
        &                                                                  "sgpm   ", "cs_iso2", "b0     ", "u0     "  ]

   interface write_arr_to_restart
      module procedure write_4darr_to_restart, write_3darr_to_restart
   end interface

   interface read_arr_from_restart
      module procedure read_4darr_from_restart, read_3darr_from_restart
   end interface

contains

!>
!! \brief Routine to set parameters and dimensions of arrays in restart file
!! \param area_type case name; possibilities:
!!   AT_ALL_B - whole domain with mpi boundaries,
!!   AT_OUT_B - physical domain with outer boundaries,
!    AT_NO_B  - only physical domain without any boundaries
!! \param area grid dimensions in the file
!! \param chnk dimensions of data array dumped by this process
!! \param lleft left limits of data from array to be dumped
!! \param lright right limits of data from array to be dumped
!! \loffs offset in area for this process
!<
   subroutine set_dims_to_write(area_type, area, chnk, lleft, lright, loffs, cg)

      use constants,  only: ndims, AT_ALL_B, AT_OUT_B, AT_NO_B, AT_USER, LO, HI
      use dataio_pub, only: die, warn
      use domain,     only: dom, has_dir, cdd, is_uneven, is_mpi_noncart
      use grid_cont,  only: grid_container
      use mpi,        only: MPI_COMM_NULL
      use user_hooks, only: at_user_settings

      implicit none

      integer(kind=4),                   intent(in)  :: area_type
      integer,         dimension(ndims), intent(out) :: area, lleft, lright, chnk
      integer(kind=8), dimension(ndims), intent(out) :: loffs
      type(grid_container), pointer,     intent(in)  :: cg

      select case (area_type)
         case (AT_ALL_B)                           ! whole domain with mpi boundaries
            if (is_mpi_noncart) call die("[restart_hdf5:set_dims_to_write] allbnd dump is too hard to implement with noncartesian domain division") !psize, pcoords
            if (cdd%comm3d == MPI_COMM_NULL) call die("[restart_hdf5:set_dims_to_write] allbnd dump requires cdd%comm3d and cdd%")
            if (is_uneven) call warn("[restart_hdf5:set_dims_to_write] allbnd dump with uneven domain division")
            chnk(:)   = cg%n_
            area(:)   = dom%n_d(:) + 2 * cg%nb * cdd%psize(:) ! \todo invent something better
            lleft(:)  = 1
            lright(:) = chnk
            loffs(:)  = cg%off(:) + 2 * cg%nb * cdd%pcoords(:) !\todo invent something better
         case (AT_OUT_B)                                   ! physical domain with outer boundaries
            area(:)   = dom%n_t(:)
            lleft(:)  = cg%ijkse(:, LO)
            lright(:) = cg%ijkse(:, HI)
            chnk(:)   = cg%n_b(:)
            where (cg%off(:) == 0 .and. has_dir(:))
               lleft(:)  = lleft(:)  - cg%nb
               chnk(:)   = chnk(:)   + cg%nb
            endwhere
            where (cg%h_cor1(:) == dom%n_d(:) .and. has_dir(:))
               lright(:) = lright(:) + cg%nb
               chnk(:)   = chnk(:)   + cg%nb
            endwhere
            loffs(:)  = cg%off(:)
            where (loffs(:)>0) loffs(:) = loffs(:) + cg%nb ! the block adjacent to the left boundary are cg%nb cells wider than cg%n_b(:)
         case (AT_NO_B)                                    ! only physical domain without any boundaries
            area(:)   = dom%n_d(:)
            lleft(:)  = cg%ijkse(:, LO)
            lright(:) = cg%ijkse(:, HI)
            chnk(:)   = cg%n_b(:)
            loffs(:)  = cg%off(:)
         case (AT_USER)                                    ! user defined domain (with no reference to simulations domain)
            if (associated(at_user_settings)) then
               call at_user_settings(area, lleft, lright, chnk, loffs)
            else
               call die("[restart_hdf5:set_dims_to_write] Routine at_user_settings not associated")
            endif
         case default
            call die("[restart_hdf5:set_dims_to_write] Non-recognized area_type.")
            area(:) = 0 ! suppres compiler warnings
      end select

   end subroutine set_dims_to_write

!>
!! \brief This routine writes restart dump and updates restart counter
!<

   subroutine write_restart_hdf5(debug_res)

      use constants,   only: cwdlen, AT_IGNORE, AT_ALL_B, AT_OUT_B, AT_NO_B, I_ONE, FLUID, MAG
      use dataio_hdf5, only: set_common_attributes
      use dataio_pub,  only: chdf, nres, set_container_chdf, problem_name, run_id, msg, printio, hdf
      use global,      only: nstep
      use grid,        only: all_cg
      use gc_list,     only: cg_list_element
      use grid_cont,   only: grid_container
      use hdf5,        only: HID_T, H5P_FILE_ACCESS_F, H5F_ACC_TRUNC_F, h5open_f, h5close_f, h5fcreate_f, h5fclose_f, h5pcreate_f, h5pclose_f, h5pset_fapl_mpio_f
      !, H5P_DATASET_XFER_F, h5pset_preserve_f
      use dataio_user, only: problem_write_restart
      use mpi,         only: MPI_CHARACTER
      use mpisetup,    only: comm, info, ierr, master, FIRST

      implicit none

      logical, optional, intent(in) :: debug_res

      integer               :: i
      integer, parameter    :: extlen = 4
      character(len=extlen) :: file_extension
      character(len=cwdlen) :: filename  !> HDF File name
      integer(HID_T)        :: file_id       !> File identifier
      integer(HID_T)        :: plist_id      !> Property list identifier
      integer(kind=4)       :: area_type
      integer(kind=4)       :: error
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      ! Construct the filename
      if (present(debug_res)) then
         file_extension = '.dmp'
      else
         file_extension = '.res'
      endif

      if (master) then
         write(filename, '(a,a1,a3,a1,i4.4,a4)') trim(problem_name), '_', run_id, '_', nres, file_extension
         write(msg,'(3a)') 'Writing restart ', trim(filename), " ... "
         call printio(msg, .true.)
      endif
      call MPI_Bcast(filename, cwdlen, MPI_CHARACTER, FIRST, comm, ierr)
      call set_container_chdf(nstep)

      ! Set up a new HDF5 file for parallel write
      call h5open_f(error)
      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, error)
      call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)

      ! Write all data in parallel
      cgl => all_cg%first
      do while (associated(cgl))
         cg => cgl%cg

         !> \todo where (cg%q(:)%restart), write cg%q(:)%arr automatically, elsewhere write just names
         if (associated(problem_write_restart)) call problem_write_restart(file_id, cg)

         do i = lbound(cg%q(:), dim=1), ubound(cg%q(:), dim=1)
            if (cg%q(i)%restart_mode /= AT_IGNORE) call write_arr_to_restart(file_id, cg%q(i)%arr, cg%q(i)%restart_mode, cg%q(i)%name, cg)
         enddo

         ! Write fluids
         area_type = AT_NO_B
         if (present(debug_res)) area_type = AT_ALL_B
         if (associated(cg%u%arr)) call write_arr_to_restart(file_id, cg%u%arr, area_type, dname(FLUID), cg)

         ! Write magnetic field
         area_type = AT_OUT_B ! unlike fluids, we need magnetic field boundaries values. Then chunks might be non-uniform
         if (present(debug_res)) area_type = AT_ALL_B
         if (associated(cg%b%arr)) call write_arr_to_restart(file_id, cg%b%arr, area_type, dname(MAG), cg)
         cgl => cgl%nxt
      enddo

!     \todo writing axes using collective I/O takes order of magnitude more than
!        dumping U and B arrays alltogether, since XYZ-axis is not even read
!        back during restart I'm commenting this out. Rewrite or punt.
!      call write_axes_to_restart(file_id)

      ! End of parallel writing (close the HDF5 file stuff)
      call h5pclose_f(plist_id, error)

      ! dump cg
      !call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      !call h5pset_preserve_f(plist_id, .true., error)
      !call write_grid_containter(cg, file_id, plist_id)
      !call h5pclose_f(plist_id, error)
      !
      call h5fclose_f(file_id, error)
      call h5close_f(error)

      ! Write some global variables
      call set_common_attributes(filename, chdf)

      nres = nres + I_ONE

   end subroutine write_restart_hdf5

   !----------------------------------------------------------------------------------
   ! Write fluid, mag or other variables (4-D and 3-D arrays)
   !
   ! Having both rank-3 array pointer and rank-4 array pointer doesn;t look elegant, but works.
   ! Is there a way to pass only one, "universal" array pointer in Fortran?

   subroutine prep_arr_write(rank, ir, area_type, loffs, chunk_dims, dimsf, file_id, dname, memspace, plist_id, filespace, dset_id, dplist_id, dfilespace)

      use hdf5,       only: HID_T, HSIZE_T, H5T_NATIVE_DOUBLE, &
           &                H5P_DATASET_CREATE_F, H5S_SELECT_SET_F, H5P_DATASET_XFER_F, H5FD_MPIO_INDEPENDENT_F, H5FD_MPIO_COLLECTIVE_F, &
           &                h5screate_simple_f, h5pcreate_f, h5dcreate_f, h5dget_space_f, &
           &                h5pset_chunk_f, h5pset_dxpl_mpio_f, h5sselect_hyperslab_f
      use domain,     only: is_uneven
      use constants,  only: ndims, AT_OUT_B, LONG

      implicit none

      integer, parameter                             :: rank4 = 1 + ndims
      integer(kind=4), intent(in)                    :: rank
      integer, intent(in)                            :: ir
      integer(kind=4), intent(in)                    :: area_type !> no boundaries, only outer boundaries or all boundaries
      integer(HSIZE_T), dimension(rank4), intent(in) :: chunk_dims, dimsf
      integer(kind=8),  dimension(ndims), intent(in) :: loffs
      integer(HID_T), intent(in)                     :: file_id   !> File identifier
      character(len=*), intent(in)                   :: dname

      integer(HID_T), intent(out) :: dset_id    !> Dataset identifier
      integer(HID_T), intent(out) :: filespace  !>
      integer(HID_T), intent(out) :: dfilespace !>
      integer(HID_T), intent(out) :: memspace   !> Dataspace identifier in memory
      integer(HID_T), intent(out) :: plist_id   !>
      integer(HID_T), intent(out) :: dplist_id  !>

      integer(HSIZE_T), dimension(rank4) :: count, offset, stride, block
      integer(kind=4) :: error

      ! Create the file space for the dataset and make it chunked if possible
      call h5screate_simple_f(rank, dimsf(ir:), dfilespace, error)
      call h5pcreate_f(H5P_DATASET_CREATE_F, dplist_id, error)
      if (.not. (is_uneven .or. (area_type == AT_OUT_B))) call h5pset_chunk_f(dplist_id, rank, chunk_dims(ir:), error)
      call h5dcreate_f(file_id, dname, H5T_NATIVE_DOUBLE, dfilespace, dset_id, error, dplist_id)

      ! Each process defines dataset in memory and writes it to the hyperslab in the file.
      call h5dget_space_f(dset_id, filespace, error)

      stride(:) = 1
      count(:)  = 1
      block(:)  = chunk_dims(:)
      offset(:) = [0_LONG, loffs(:)]
      call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset(ir:), count(ir:), error, stride(ir:), block(ir:))

      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)

#ifdef INDEPENDENT_ATOUTB
      !> \warning collective write may be unproper (some fields may stay unwritten)
      !! seen for: Intel compiler version 12.0.4, HDF5 1.8.5 , OpenMPI 1.4.1
      !!    n_d = [250, 250, 200], nb = 5, psize = [ 5,  5, 8]
      !! or n_d = [500, 500, 200], nb = 5, psize = [10, 10, 4]
      !! \warning restart write may take tens of minutes for settings like above
      !! \todo check if similar problem exists for is_uneven
      !! \todo remove when problem is well resolved with newer versions of HDF5
      !<
      if (area_type == AT_OUT_B) then
         call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
      else
#endif /* INDEPENDENT_ATOUTB */
         call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
#ifdef INDEPENDENT_ATOUTB
      endif
#endif /* INDEPENDENT_ATOUTB */

      call h5screate_simple_f(rank, chunk_dims(ir:), memspace, error)
      return
   end subroutine prep_arr_write

   subroutine clean_arr_write(memspace, plist_id, filespace, dset_id, dplist_id, dfilespace)
      use hdf5,       only: HID_T, h5sclose_f, h5pclose_f, h5dclose_f
      implicit none
      integer(HID_T), intent(inout) :: memspace, plist_id, filespace, dset_id, dplist_id, dfilespace
      integer(kind=4) :: error

      call h5sclose_f(memspace, error)
      call h5pclose_f(plist_id, error)
      call h5sclose_f(filespace, error)
      call h5dclose_f(dset_id, error)
      call h5pclose_f(dplist_id, error)
      call h5sclose_f(dfilespace, error)
   end subroutine clean_arr_write

   subroutine write_4darr_to_restart(file_id, pa4d, area_type, dname, cg)

      use constants,  only: xdim, ydim, zdim, ndims
      use dataio_pub, only: die
      use grid_cont,  only: grid_container
      use hdf5,       only: HID_T, HSIZE_T, H5T_NATIVE_DOUBLE, h5dwrite_f

      implicit none

      integer(HID_T), intent(in)                    :: file_id   !> File identifier
      real, pointer, dimension(:,:,:,:), intent(in) :: pa4d      !> 4-D array pointer
      integer(kind=4), intent(in)                   :: area_type !> no boundaries, only outer boundaries or all boundaries
      character(len=*), intent(in)                  :: dname
      type(grid_container), pointer, intent(in)     :: cg

      integer, parameter :: rank4 = 1 + ndims
      integer(HSIZE_T), dimension(rank4) :: dimsf, chunk_dims
      integer(HID_T) :: dset_id               !> Dataset identifier
      integer(HID_T) :: dplist_id, plist_id   !> Property list identifiers
      integer(HID_T) :: dfilespace, filespace !> Dataspace identifiers in file
      integer(HID_T) :: memspace              !> Dataspace identifier in memory

      integer,         dimension(ndims) :: area, lleft, lright, chnk
      integer(kind=8), dimension(ndims) :: loffs

      integer(kind=4) :: rank
      integer(kind=4) :: error

      call set_dims_to_write(area_type, area, chnk, lleft, lright, loffs, cg)

      if (.not.associated(pa4d)) call die("[restart_hdf5:write_4darr_to_restart] Null pointer given.")
      dimsf      = [size(pa4d,1), area(:)]   ! Dataset dimensions
      chunk_dims = [size(pa4d,1), chnk(:)]   ! Chunks dimensions
      rank = rank4

      call prep_arr_write(rank, 1, area_type, loffs, chunk_dims, dimsf, file_id, dname, memspace, plist_id, filespace, dset_id, dplist_id, dfilespace)

      ! write data
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, pa4d(:, lleft(xdim):lright(xdim), lleft(ydim):lright(ydim), lleft(zdim):lright(zdim)), &
           &          dimsf(:), error, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

      call clean_arr_write(memspace, plist_id, filespace, dset_id, dplist_id, dfilespace)

   end subroutine write_4darr_to_restart

   subroutine write_3darr_to_restart(file_id, pa3d, area_type, dname, cg)

      use constants,  only: xdim, ydim, zdim, ndims
      use dataio_pub, only: die
      use grid_cont,  only: grid_container
      use hdf5,       only: HID_T, HSIZE_T, H5T_NATIVE_DOUBLE, h5dwrite_f

      implicit none

      integer(HID_T), intent(in)                    :: file_id   !> File identifier
      real, pointer, dimension(:,:,:), intent(inout):: pa3d      !> \deprecated 3-D array pointer; not required for cg%q(:) arrays
      integer(kind=4), intent(in)                   :: area_type !> no boundaries, only outer boundaries or all boundaries
      character(len=*), intent(in)                  :: dname
      type(grid_container), pointer, intent(in)     :: cg

      integer, parameter :: rank4 = 1 + ndims
      integer(HSIZE_T), dimension(rank4) :: dimsf, chunk_dims
      integer(HID_T) :: dset_id               !> Dataset identifier
      integer(HID_T) :: dplist_id, plist_id   !> Property list identifiers
      integer(HID_T) :: dfilespace, filespace !> Dataspace identifiers in file
      integer(HID_T) :: memspace              !> Dataspace identifier in memory

      integer,         dimension(ndims) :: area, lleft, lright, chnk
      integer(kind=8), dimension(ndims) :: loffs

      integer(kind=4) :: rank
      integer(kind=4) :: error

      if (.not. associated(pa3d)) then
         pa3d => cg%get_na_ptr(dname)
         if (.not. associated(pa3d)) call die("[restart_hdf5:write_3darr_to_restart] Null pointer given.")
      endif

      call set_dims_to_write(area_type, area, chnk, lleft, lright, loffs, cg)

      dimsf = [1, area(:)]      ! Dataset dimensions
      chunk_dims = [1, chnk(:)] ! Chunks dimensions

      rank = ndims

      call prep_arr_write(rank, 2, area_type, loffs, chunk_dims, dimsf, file_id, dname, memspace, plist_id, filespace, dset_id, dplist_id, dfilespace)

      ! write data
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, pa3d(lleft(xdim):lright(xdim), lleft(ydim):lright(ydim), lleft(zdim):lright(zdim)), &
           &          dimsf(2:), error, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

      call clean_arr_write(memspace, plist_id, filespace, dset_id, dplist_id, dfilespace)

   end subroutine write_3darr_to_restart

!!$   !----------------------------------------------------------------------------------
!!$   !
!!$   !  WRITE Axes
!!$   !
!!$   ! AMR: the axes should be associated with fluid and datasets
!!$   subroutine write_axes_to_restart(file_id)
!!$
!!$      use constants,  only: xdim, ydim, zdim, ndims
!!$      use hdf5,       only: HSIZE_T, HID_T, H5T_NATIVE_DOUBLE, H5FD_MPIO_COLLECTIVE_F, H5P_DATASET_CREATE_F, H5P_DATASET_XFER_F, H5S_SELECT_SET_F, &
!!$           &                h5screate_simple_f, h5pcreate_f, h5pset_chunk_f, h5pset_dxpl_mpio_f, h5dcreate_f, h5dget_space_f, h5dwrite_f, &
!!$           &                h5sclose_f, h5pclose_f, h5dclose_f, h5sselect_hyperslab_f
!!$      use domain,     only: dom, is_uneven
!!$      use grid,       only: cg
!!$
!!$      implicit none
!!$
!!$      integer(HID_T), intent(in) :: file_id !> File identifier
!!$
!!$      integer :: dir
!!$      integer :: error
!!$      integer, parameter :: rank=1
!!$      integer(HSIZE_T), dimension(rank) :: offset, count, stride, block, dimsf, chunk_dims
!!$      integer(HID_T) :: dset_id                !> Dataset identifier
!!$      integer(HID_T) :: dplist_id, plist_id    !> Property list identifiers
!!$      integer(HID_T) :: dfilespace, filespace  !> Dataspace identifiers in file
!!$      integer(HID_T) :: memspace               !> Dataspace identifier in memory
!!$      integer, parameter :: asis_n_len = 6
!!$      character(len=asis_n_len) :: dset_axis_n !> Dataspace name
!!$      character(len=ndims), parameter :: axis_n = "XYZ"
!!$
!!$      do dir = xdim, zdim
!!$
!!$         write(dset_axis_n,'(a,"-axis")')axis_n(dir:dir)
!!$         dimsf      = [dom%n_d(dir)] ! Dataset dimensions
!!$         chunk_dims = [cg%n_b(dir)]  ! Chunk dimensions
!!$
!!$         ! Create the file space for the dataset and make it chunked if possible
!!$         call h5screate_simple_f(rank, dimsf, dfilespace, error)
!!$         call h5pcreate_f(H5P_DATASET_CREATE_F, dplist_id, error)
!!$         if (.not. is_uneven) call h5pset_chunk_f(dplist_id, rank, chunk_dims, error)
!!$         call h5dcreate_f(file_id, dset_axis_n, H5T_NATIVE_DOUBLE, dfilespace, dset_id, error, dplist_id)
!!$
!!$         ! It is sufficient that only one CPU writes each piece of axis data.
!!$         ! The other CPUs also need to call at least everything up to h5dcreate_f, and then just close the stuff.
!!$         if (cg%off(mod(dir, ndims)+1) == 0 .and. cg%off(mod(dir+ndims-2, ndims)+1) == 0) then
!!$
!!$            call h5dget_space_f(dset_id, filespace, error)
!!$
!!$            ! Each contributing process defines dataset in memory and writes it to the hyperslab in the file.
!!$            stride(:) = 1
!!$            count(:)  = 1
!!$            block(:)  = chunk_dims(:)
!!$            offset(:) = [ cg%off(dir) ]
!!$            call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error, stride, block)
!!$
!!$            call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
!!$            call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
!!$
!!$            ! Create the memory space for the dataset.
!!$            call h5screate_simple_f(rank, chunk_dims, memspace, error)
!!$            select case (dir)
!!$               case (xdim)
!!$                  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, cg%x(cg%is:cg%ie), dimsf, error, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
!!$               case (ydim)
!!$                  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, cg%y(cg%js:cg%je), dimsf, error, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
!!$               case (zdim)
!!$                  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, cg%z(cg%ks:cg%ke), dimsf, error, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
!!$            end select
!!$            call h5sclose_f(memspace, error)
!!$            call h5pclose_f(plist_id, error)
!!$            call h5sclose_f(filespace, error)
!!$
!!$         endif
!!$
!!$         call h5dclose_f(dset_id, error)
!!$         call h5pclose_f(dplist_id, error)
!!$         call h5sclose_f(dfilespace, error)
!!$
!!$      enddo
!!$
!!$   end subroutine write_axes_to_restart

! This routine reads only interior cells of a pointed array from the restart file.
! Boundary cells are exchanged with the neughbours. Corner boundary cells are not guaranteed to be correct (area_type = AT_ALL_B not implemented yet).
! External boundary cells are not stored in the restart file and thus all of them are lost (area_type = AT_OUT_B not implemented yet).

   subroutine prep_arr_read(rank, ir, loffs, chunk_dims, file_id, dname, memspace, plist_id, filespace, dset_id)

      use hdf5,       only: HID_T, HSIZE_T, H5S_SELECT_SET_F, H5P_DATASET_XFER_F, H5FD_MPIO_COLLECTIVE_F, &
           &                h5dopen_f, h5sget_simple_extent_ndims_f, h5dget_space_f, &
           &                h5pcreate_f, h5pset_dxpl_mpio_f, h5sselect_hyperslab_f, h5screate_simple_f
      use constants,  only: ndims, LONG
      use dataio_pub, only: msg, die

      implicit none

      integer, parameter                             :: rank4 = 1 + ndims
      integer(kind=4), intent(in)                    :: rank
      integer, intent(in)                            :: ir
      integer(HSIZE_T), dimension(rank4), intent(in) :: chunk_dims
      integer(kind=8),  dimension(ndims), intent(in) :: loffs
      integer(HID_T), intent(in)                     :: file_id   !> File identifier
      character(len=*), intent(in)                   :: dname

      integer(HID_T), intent(out) :: dset_id    !> Dataset identifier
      integer(HID_T), intent(out) :: filespace  !>
      integer(HID_T), intent(out) :: memspace   !> Dataspace identifier in memory
      integer(HID_T), intent(out) :: plist_id   !>

      integer(HSIZE_T), dimension(rank4) :: count, offset, stride, block
      integer(kind=4) :: error
      integer(kind=4) :: rankf

      ! Create dataset.and filespace
      call h5dopen_f(file_id, dname, dset_id, error)
      if (error /= 0) then
         write(msg, '(3a)') "[restart_hdf5:prep_arr_read] Opening dataset '", dname,"' failed."
         call die(msg)
      endif

      call h5dget_space_f(dset_id, filespace, error)
      call h5sget_simple_extent_ndims_f (filespace, rankf, error)
      if (rank /= rankf) then
         write(msg,'(3a,2(i2,a))')"[restart_hdf5:prep_arr_read] Rank mismatch in array '", dname, "' (", rank, " /= ", rankf, ")"
         call die(msg)
      endif

      ! Select hyperslab in the file.
      stride(:) = 1
      count(:)  = 1
      block(:)  = chunk_dims(:)
      offset(:) = [ 0_LONG, loffs(:) ]
      call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset(ir:), count(ir:), error, stride(ir:), block(ir:))

      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
      call h5screate_simple_f(rank, chunk_dims(ir:), memspace, error)

   end subroutine prep_arr_read

   subroutine clean_arr_read(memspace, plist_id, filespace, dset_id)
      use hdf5,       only: HID_T, h5sclose_f, h5pclose_f, h5dclose_f
      implicit none
      integer(HID_T), intent(inout) :: memspace, plist_id, filespace, dset_id
      integer(kind=4) :: error

      call h5sclose_f(memspace, error)
      call h5pclose_f(plist_id, error)
      call h5sclose_f(filespace, error)
      call h5dclose_f(dset_id, error)

   end subroutine clean_arr_read

   subroutine read_4darr_from_restart(file_id, pa4d, area_type, dname, cg)

      use constants,  only: xdim, ydim, zdim, ndims
      use dataio_pub, only: msg, die
      use grid_cont,  only: grid_container
      use hdf5,       only: HID_T, HSIZE_T, SIZE_T, H5T_NATIVE_DOUBLE, h5dread_f

      implicit none

      integer(HID_T), intent(in)                       :: file_id   ! File identifier
      real, dimension(:,:,:,:), pointer, intent(inout) :: pa4d      ! pointer to (:, 1:cg%nx, 1:cg%ny, 1:cg%nz)-sized array
      integer(kind=4), intent(in)                      :: area_type !> no boundaries, only outer boundaries or all boundaries
      character(len=*), intent(in)                     :: dname
      type(grid_container), pointer, intent(in)        :: cg

      integer(HID_T)        :: dset_id       ! Dataset identifier
      integer(HID_T)        :: plist_id      ! Property list identifier
      integer(HID_T)        :: filespace     ! Dataspace identifier in file
      integer(HID_T)        :: memspace      ! Dataspace identifier in memory

      integer, parameter :: rank4 = 1 + ndims
      integer(HSIZE_T), dimension(rank4) :: dimsf, chunk_dims

      integer,         dimension(ndims) :: area, lleft, lright, chnk
      integer(kind=8), dimension(ndims) :: loffs
      integer :: ir
      integer(kind=4) :: rank
      integer(kind=4) :: error

      call set_dims_to_write(area_type, area, chnk, lleft, lright, loffs, cg)

      dimsf = [size(pa4d,1), area(:)]      ! Dataset dimensions
      chunk_dims = [size(pa4d,1), chnk(:)] ! Chunks dimensions
      rank = rank4
      ir = rank4 - rank + 1 ! 1 for 4-D arrays, 2 for 3-D arrays (to simplify use of count(:), offset(:), stride(:), block(:), dimsf(:) and chunk_dims(:)

      call prep_arr_read(rank, ir, loffs, chunk_dims, file_id, dname, memspace, plist_id, filespace, dset_id)

      ! Read the array
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, pa4d(:, lleft(xdim):lright(xdim), lleft(ydim):lright(ydim), lleft(zdim):lright(zdim)), &
           &         dimsf(ir:), error, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

      if (error /= 0) then
         write(msg, '(3a)') "[restart_hdf5:read_4darr_from_restart] Reading dataset '", dname,"' failed."
         call die(msg)
      endif

      call clean_arr_read(memspace, plist_id, filespace, dset_id)
      ! rank-4 arrays (cg%u%arr(:,:,:,:) and b(:,:,:,:)) have their own guardcell-exchange routines, which can also be called here

   end subroutine read_4darr_from_restart

   subroutine read_3darr_from_restart(file_id, pa3d, area_type, dname, cg)

      use constants,    only: xdim, ydim, zdim, ndims
      use dataio_pub,   only: msg, die
      use internal_bnd, only: arr3d_boundaries
      use grid_cont,    only: grid_container
      use hdf5,         only: HID_T, HSIZE_T, SIZE_T, h5dread_f, H5T_NATIVE_DOUBLE

      implicit none

      integer(HID_T), intent(in)                       :: file_id   ! File identifier
      real, dimension(:,:,:), pointer, intent(inout)   :: pa3d      ! \deprecated pointer to (1:cg%nx, 1:cg%ny, 1:cg%nz)-sized array; not required for cg%q(:) arrays
      integer(kind=4), intent(in)                      :: area_type !> no boundaries, only outer boundaries or all boundaries
      character(len=*), intent(in)                     :: dname
      type(grid_container), pointer, intent(in)        :: cg

      integer(HID_T)        :: dset_id       ! Dataset identifier
      integer(HID_T)        :: plist_id      ! Property list identifier
      integer(HID_T)        :: filespace     ! Dataspace identifier in file
      integer(HID_T)        :: memspace      ! Dataspace identifier in memory

      integer, parameter :: rank4 = 1 + ndims
      integer(HSIZE_T), dimension(rank4) :: dimsf, chunk_dims

      integer,         dimension(ndims) :: area, lleft, lright, chnk
      integer(kind=8), dimension(ndims) :: loffs
      integer :: ir
      integer(kind=4) :: rank
      integer(kind=4) :: error

      call set_dims_to_write(area_type, area, chnk, lleft, lright, loffs, cg)

      dimsf = [1, area(:)]      ! Dataset dimensions
      chunk_dims = [1, chnk(:)] ! Chunks dimensions
      if (.not. associated(pa3d)) then
         if (cg%exists(dname)) call die("[restart_hdf5:read_3darr_from_restart] Already read.")
         call cg%add_na(dname)
         pa3d => cg%get_na_ptr(dname)
         if (.not. associated(pa3d)) call die("[restart_hdf5:read_3darr_from_restart] Null pointer given.")
      endif
      rank = ndims
      ir = rank4 - rank + 1 ! 1 for 4-D arrays, 2 for 3-D arrays (to simplify use of count(:), offset(:), stride(:), block(:), dimsf(:) and chunk_dims(:)

      call prep_arr_read(rank, ir, loffs, chunk_dims, file_id, dname, memspace, plist_id, filespace, dset_id)

      ! Read the array
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, pa3d(lleft(xdim):lright(xdim), lleft(ydim):lright(ydim), lleft(zdim):lright(zdim)), &
           &         dimsf(ir:), error, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

      if (error /= 0) then
         write(msg, '(3a)') "[restart_hdf5:read_3darr_from_restart] Reading dataset '", dname,"' failed."
         call die(msg)
      endif

      call clean_arr_read(memspace, plist_id, filespace, dset_id)

      ! Originally the pa3d array was written with the guardcells. The internal guardcells will be exchanged but the external ones are lost.
      call arr3d_boundaries(pa3d, area_type=area_type, dname=dname)

   end subroutine read_3darr_from_restart

   subroutine read_restart_hdf5(chdf)

      use constants,   only: cwdlen, cbuff_len, domlen, idlen, xdim, ydim, zdim, AT_IGNORE, AT_NO_B, AT_OUT_B, LO, HI, I_ONE, FLUID, MAG
      use dataio_pub,  only: msg, printio, warn, die, require_init_prob, problem_name, run_id, piernik_hdf5_version, hdf
      use dataio_user, only: problem_read_restart
      use domain,      only: dom, has_dir
      use fluidindex,  only: flind
      use func,        only: fix_string
      use global,      only: magic_mass, t, dt
      use grid,        only: all_cg
      use gc_list,     only: cg_list_element
      use grid_cont,   only: grid_container
      use hdf5,        only: HID_T, SIZE_T, H5P_FILE_ACCESS_F, H5F_ACC_RDONLY_F, &
           &                 h5open_f, h5pcreate_f, h5pset_fapl_mpio_f, h5fopen_f, h5pclose_f, h5fclose_f, h5close_f
      use h5lt,        only: h5ltget_attribute_double_f, h5ltget_attribute_int_f, h5ltget_attribute_string_f
      use mpi,         only: MPI_CHARACTER, MPI_INTEGER, MPI_DOUBLE_PRECISION
      use mpisetup,    only: comm, ierr, info, comm, master, FIRST

      implicit none

      type(hdf)             :: chdf
      integer               :: nu, i
      character(len=cwdlen) :: filename  ! File name

      integer(HID_T)        :: file_id       ! File identifier
      integer(HID_T)        :: plist_id      ! Property list identifier
      integer(SIZE_T)       :: bufsize

      integer(kind=4)       :: error
      logical               :: file_exist

      real, dimension(1)    :: rbuf
      integer(kind=4), dimension(1) :: ibuf

      real                  :: restart_hdf5_version
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      nu = flind%all

      if (master) then
         write(filename,'(a,a1,a3,a1,i4.4,a4)') trim(problem_name),'_', run_id,'_', chdf%nres,'.res'
         write(msg, '(2a)') 'Reading restart file: ', trim(filename)
         call printio(msg)
      endif
      call MPI_Bcast(filename, cwdlen, MPI_CHARACTER, FIRST, comm, ierr)

      inquire(file = filename, exist = file_exist)
      if (.not. file_exist) then
         write(msg,'(3a)') '[restart_hdf5:read_restart_hdf5]: Restart file: ', trim(filename),' does not exist'
         call die(msg)
      endif

      call h5open_f(error)
      if (master) then
         call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, error)

         call h5ltget_attribute_double_f(file_id,"/","piernik", rbuf, error)
         if (error /= 0) call die("[restart_hdf5:read_restart_hdf5] Cannot read 'piernik' attribute from the restart file. The file may be either damaged or incompatible")
         if (rbuf(1) > piernik_hdf5_version) then
            write(msg,'(2(a,f5.2))')"[restart_hdf5:read_restart_hdf5] Cannot read future versions of the restart file: ", rbuf(1)," > ", piernik_hdf5_version
            call die(msg)
         else if (int(rbuf(1)) < int(piernik_hdf5_version)) then
            write(msg,'(2(a,f5.2))')"[restart_hdf5:read_restart_hdf5] The restart file is too ancient. It is unlikely that it could work correctly: ", rbuf(1)," << ", piernik_hdf5_version
            call die(msg)
         else if (rbuf(1) < piernik_hdf5_version) then
            write(msg,'(2(a,f5.2))')"[restart_hdf5:read_restart_hdf5] Old versions of the restart file may not always work fully correctly: ", rbuf(1)," < ", piernik_hdf5_version
            call warn(msg)
         endif

         call h5ltget_attribute_int_f(file_id,"/","nxd", ibuf, error)
         if (ibuf(1) /= dom%n_d(xdim) .or. error /= 0) call die("[restart_hdf5:read_restart_hdf5] nxd does not match")
         if (has_dir(xdim)) then
            call h5ltget_attribute_double_f(file_id,"/","xmin", rbuf, error)
            if (rbuf(1) /= dom%edge(xdim, LO) .or. error /= 0) call die("[restart_hdf5:read_restart_hdf5] xmin does not match")
            call h5ltget_attribute_double_f(file_id,"/","xmax", rbuf, error)
            if (rbuf(1) /= dom%edge(xdim, HI) .or. error /= 0) call die("[restart_hdf5:read_restart_hdf5] xmax does not match")
         endif

         call h5ltget_attribute_int_f(file_id,"/","nyd", ibuf, error)
         if (ibuf(1) /= dom%n_d(ydim) .or. error /= 0) call die("[restart_hdf5:read_restart_hdf5] nyd does not match")
         if (has_dir(ydim)) then
            call h5ltget_attribute_double_f(file_id,"/","ymin", rbuf, error)
            if (rbuf(1) /= dom%edge(ydim, LO) .or. error /= 0) call die("[restart_hdf5:read_restart_hdf5] ymin does not match")
            call h5ltget_attribute_double_f(file_id,"/","ymax", rbuf, error)
            if (rbuf(1) /= dom%edge(ydim, HI) .or. error /= 0) call die("[restart_hdf5:read_restart_hdf5] ymax does not match")
         endif

         call h5ltget_attribute_int_f(file_id,"/","nzd", ibuf, error)
         if (ibuf(1) /= dom%n_d(zdim) .or. error /= 0) call die("[restart_hdf5:read_restart_hdf5] nzd does not match")
         if (has_dir(zdim)) then
            call h5ltget_attribute_double_f(file_id,"/","zmin", rbuf, error)
            if (rbuf(1) /= dom%edge(zdim, LO) .or. error /= 0) call die("[restart_hdf5:read_restart_hdf5] zmin does not match")
            call h5ltget_attribute_double_f(file_id,"/","zmax", rbuf, error)
            if (rbuf(1) /= dom%edge(zdim, HI) .or. error /= 0) call die("[restart_hdf5:read_restart_hdf5] zmax does not match")
         endif

         call h5fclose_f(file_id, error)
      endif

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, error)

      call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, error, access_prp = plist_id)
      call h5pclose_f(plist_id, error)

      cgl => all_cg%first
      do while (associated(cgl))
         cg => cgl%cg

         !> \todo read existing cg%q(:)%arr automatically, create fresh cg%q(:)%arr where (.not. cg%q(:)%restart)
         if (associated(problem_read_restart)) call problem_read_restart(file_id, cg)

         do i = lbound(cg%q(:), dim=1), ubound(cg%q(:), dim=1)
            if (cg%q(i)%restart_mode /= AT_IGNORE) call read_arr_from_restart(file_id, cg%q(i)%arr, cg%q(i)%restart_mode, cg%q(i)%name, cg)
         enddo

         !  READ FLUID VARIABLES
         if (associated(cg%u%arr)) call read_arr_from_restart(file_id, cg%u%arr, AT_NO_B, dname(FLUID), cg)

         !  READ MAG VARIABLES
         if (associated(cg%b%arr)) call read_arr_from_restart(file_id, cg%b%arr, AT_OUT_B, dname(MAG), cg)
         cgl => cgl%nxt
      enddo

      call h5fclose_f(file_id, error)

      if (master) then
         call h5fopen_f (filename, H5F_ACC_RDONLY_F, file_id, error)
         bufsize = 1
         call h5ltget_attribute_double_f(file_id,"/","piernik", rbuf, error)
         restart_hdf5_version = rbuf(1)
         call h5ltget_attribute_double_f(file_id,"/","time", rbuf, error)
         t = rbuf(1)
         call h5ltget_attribute_double_f(file_id,"/","timestep", rbuf, error)
         dt = rbuf(1)
         call h5ltget_attribute_double_f(file_id,"/","magic_mass", rbuf, error)
         magic_mass = rbuf(1)
         call h5ltget_attribute_int_f(file_id,"/","nstep", ibuf, error)
         chdf%nstep = ibuf(1)
         call h5ltget_attribute_int_f(file_id,"/","nres", ibuf, error)
         chdf%nres = ibuf(1)
         call h5ltget_attribute_int_f(file_id,"/","nhdf", ibuf, error)
         chdf%nhdf = ibuf(1)
         call h5ltget_attribute_int_f(file_id,"/","step_res", ibuf, error)
         chdf%step_res = ibuf(1)
         call h5ltget_attribute_int_f(file_id,"/","step_hdf", ibuf, error)
         chdf%step_hdf = ibuf(1)
         call h5ltget_attribute_double_f(file_id,"/","next_t_tsl", rbuf, error)
         chdf%next_t_tsl = rbuf(1)
         call h5ltget_attribute_double_f(file_id,"/","next_t_log", rbuf, error)
         chdf%next_t_log = rbuf(1)
         call h5ltget_attribute_double_f(file_id,"/","last_hdf_time", rbuf, error)
         chdf%last_hdf_time = rbuf(1)

         call h5ltget_attribute_string_f(file_id,"/","problem_name", problem_name, error)
         call h5ltget_attribute_string_f(file_id,"/","domain", chdf%domain_dump, error)
         call h5ltget_attribute_string_f(file_id,"/","run_id", chdf%new_id, error)

         if (restart_hdf5_version > 1.11) then
            call h5ltget_attribute_int_f(file_id,"/","require_init_prob", ibuf, error)
            require_init_prob = ibuf(1)
         endif

         problem_name = fix_string(problem_name)   !> \deprecated BEWARE: >=HDF5-1.8.4 has weird issues with strings
         chdf%new_id  = fix_string(chdf%new_id)    !> \deprecated   this bit hacks it around
         chdf%domain_dump  = fix_string(chdf%domain_dump)

         call h5fclose_f(file_id, error)

         write(msg,'(2a)') 'Done reading restart file: ', trim(filename)
         call printio(msg)
      endif
      call h5close_f(error)

      call MPI_Bcast(restart_hdf5_version,    I_ONE, MPI_DOUBLE_PRECISION, FIRST, comm, ierr)

      call MPI_Bcast(chdf%nstep,    I_ONE, MPI_INTEGER, FIRST, comm, ierr)
      call MPI_Bcast(chdf%nres,     I_ONE, MPI_INTEGER, FIRST, comm, ierr)
      call MPI_Bcast(chdf%nhdf,     I_ONE, MPI_INTEGER, FIRST, comm, ierr)
      call MPI_Bcast(chdf%step_res, I_ONE, MPI_INTEGER, FIRST, comm, ierr)
      call MPI_Bcast(chdf%step_hdf, I_ONE, MPI_INTEGER, FIRST, comm, ierr)
      if (restart_hdf5_version > 1.11) call MPI_Bcast(require_init_prob, I_ONE, MPI_INTEGER, FIRST, comm, ierr)

      call MPI_Bcast(chdf%next_t_tsl,    I_ONE, MPI_DOUBLE_PRECISION, FIRST, comm, ierr)
      call MPI_Bcast(chdf%next_t_log,    I_ONE, MPI_DOUBLE_PRECISION, FIRST, comm, ierr)
      call MPI_Bcast(chdf%last_hdf_time, I_ONE, MPI_DOUBLE_PRECISION, FIRST, comm, ierr)
      call MPI_Bcast(t,                  I_ONE, MPI_DOUBLE_PRECISION, FIRST, comm, ierr)
      call MPI_Bcast(dt,                 I_ONE, MPI_DOUBLE_PRECISION, FIRST, comm, ierr)

      call MPI_Bcast(problem_name, cbuff_len, MPI_CHARACTER, FIRST, comm, ierr)
      call MPI_Bcast(chdf%domain_dump,domlen, MPI_CHARACTER, FIRST, comm, ierr)
      call MPI_Bcast(chdf%new_id,  idlen,     MPI_CHARACTER, FIRST, comm, ierr)

   end subroutine read_restart_hdf5

end module restart_hdf5
