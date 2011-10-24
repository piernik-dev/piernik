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
!! \brief Module that contains HDF5 I/O routines for reading and writing restart files.
!<
module restart_hdf5

! pulled by ANY

   implicit none

   private
   public :: read_restart_hdf5, write_restart_hdf5, read_arr_from_restart

   integer, parameter :: STAT_OK = 0

contains

!>
!! \brief Wrapper routine that writes a version 2.x or version 1.x restart file, depending on use_v2_io switch
!<

   subroutine write_restart_hdf5

      use dataio_pub, only: use_v2_io

      implicit none

      if (use_v2_io) then
         call write_restart_hdf5_v2
      else
         call write_restart_hdf5_v1
      endif

   end subroutine write_restart_hdf5

!>
!! \brief Wrapper routine that reads a version 2.x or version 1.x restart file, depending on use_v2_io switch. On v2 failure it falls back to v1.
!<

   subroutine read_restart_hdf5(chdf)

      use common_hdf5, only: hdf
      use dataio_pub, only: use_v2_io

      implicit none

      type(hdf), intent(inout) :: chdf
      integer :: status_v2

      if (use_v2_io) call read_restart_hdf5_v2(status_v2)
      if (status_v2 /= STAT_OK .or. .not. use_v2_io) call read_restart_hdf5_v1(chdf)

   end subroutine read_restart_hdf5

!>
!! \brief Routine to set dimensions of arrays related to grid containers in restart file
!! \param area_type case name; possibilities:
!!   AT_OUT_B - physical domain with outer boundaries,
!!   AT_NO_B  - only physical domain without any boundaries
!!   AT_USER  - user defined domain
!! \param chnk dimensions of data array dumped by this process
!! \param lleft left limits of data from array to be dumped
!! \param lright right limits of data from array to be dumped
!! \loffs offset in area for this process
!<
   subroutine set_dims_for_restart(area_type, chnk, lleft, lright, loffs, cg)

      use constants,  only: ndims, AT_OUT_B, AT_NO_B, AT_USER, LO, HI
      use dataio_pub, only: die
      use domain,     only: dom
      use grid_cont,  only: grid_container
      use user_hooks, only: at_user_dims_settings

      implicit none

      integer(kind=4),                   intent(in)  :: area_type
      integer,         dimension(ndims), intent(out) :: lleft, lright, chnk
      integer(kind=8), dimension(ndims), intent(out) :: loffs
      type(grid_container), pointer,     intent(in)  :: cg

      ! only physical domain without any boundaries
      lleft(:)  = cg%ijkse(:, LO)
      lright(:) = cg%ijkse(:, HI)
      chnk(:)   = cg%n_b(:)
      loffs(:)  = cg%off(:)

      select case (area_type)
         case (AT_OUT_B)                                   ! physical domain with outer boundaries
            where (cg%off(:) == 0 .and. dom%has_dir(:))
               lleft(:)  = lleft(:)  - cg%nb
               chnk(:)   = chnk(:)   + cg%nb
            endwhere
            where (cg%h_cor1(:) == dom%n_d(:) .and. dom%has_dir(:))
               lright(:) = lright(:) + cg%nb
               chnk(:)   = chnk(:)   + cg%nb
            endwhere
            where (loffs(:)>0) loffs(:) = loffs(:) + cg%nb ! the block adjacent to the left boundary are cg%nb cells wider than cg%n_b(:)
         case (AT_NO_B)                                    ! only physical domain without any boundaries
            ! Nothing special
         case (AT_USER)                                    ! user defined domain (with no reference to simulations domain)
            if (associated(at_user_dims_settings)) then
               call at_user_dims_settings(lleft, lright, chnk, loffs)
            else
               call die("[restart_hdf5:set_dims_to_write] Routine at_user_dims_settings not associated")
            endif
         case default
            call die("[restart_hdf5:set_dims_for_restart] Non-recognized area_type.")
      end select

   end subroutine set_dims_for_restart

!>
!! \brief Routine to set dimensions of arrays related to domain in restart file
!! \param area_type case name; possibilities:
!!   AT_OUT_B - physical domain with outer boundaries,
!    AT_NO_B  - only physical domain without any boundaries
!! \param area grid dimensions in the file
!<

   subroutine set_area_for_restart(area_type, area)

      use constants,  only: ndims, AT_OUT_B, AT_NO_B, AT_USER
      use dataio_pub, only: die
      use domain,     only: dom
      use user_hooks, only: at_user_area_settings

      implicit none

      integer(kind=4),           intent(in)  :: area_type
      integer, dimension(ndims), intent(out) :: area

      select case (area_type)
         case (AT_OUT_B)                                   ! physical domain with outer boundaries
            area(:) = dom%n_t(:)
         case (AT_NO_B)                                    ! only physical domain without any boundaries
            area(:) = dom%n_d(:)
         case (AT_USER)                                    ! user defined domain (with no reference to simulations domain)
            if (associated(at_user_area_settings)) then
               call at_user_area_settings(area)
            else
               call die("[restart_hdf5:set_area_for_restart] Routine at_user_area_settings not associated")
            endif
         case default
            call die("[restart_hdf5:set_area_for_restart] Non-recognized area_type.")
            area(:) = 0 ! suppress compiler warnings
      end select

   end subroutine set_area_for_restart

!>
!! \brief This routine writes restart dump and updates restart counter
!<

   subroutine write_restart_hdf5_v1

      use common_hdf5, only: set_common_attributes, chdf, set_container_chdf, hdf
      use constants,   only: cwdlen, I_ONE
      use dataio_pub,  only: nres, problem_name, run_id, msg, printio
      use global,      only: nstep
      use grid,        only: all_cg
      use grid_cont,   only: grid_container
      use hdf5,        only: HID_T, H5P_FILE_ACCESS_F, H5F_ACC_TRUNC_F, h5open_f, h5close_f, h5fcreate_f, h5fclose_f, h5pcreate_f, h5pclose_f, h5pset_fapl_mpio_f
      !, H5P_DATASET_XFER_F, h5pset_preserve_f
      use dataio_user, only: problem_write_restart
      use mpi,         only: MPI_CHARACTER, MPI_INFO_NULL
      use mpisetup,    only: comm, ierr, master, FIRST

      implicit none

      integer(kind=4)       :: i
      integer, parameter    :: extlen = 4
      character(len=extlen), parameter :: file_extension = '.res'
      character(len=cwdlen) :: filename  !> HDF File name
      integer(HID_T)        :: file_id       !> File identifier
      integer(HID_T)        :: plist_id      !> Property list identifier
      integer(kind=4)       :: error
      type(grid_container), pointer :: fcg

      ! Construct the filename

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
      call h5pset_fapl_mpio_f(plist_id, comm, MPI_INFO_NULL, error)
      call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)

      ! Write scalar fields that were declared to be required on restart
      fcg  => all_cg%first%cg
      if (allocated(fcg%q)) then
         do i = lbound(fcg%q(:), dim=1, kind=4), ubound(fcg%q(:), dim=1, kind=4)
            call write_arr_to_restart(file_id, i, .true.)
         enddo
      endif

      if (allocated(fcg%w)) then
         do i = lbound(fcg%w(:), dim=1, kind=4), ubound(fcg%w(:), dim=1, kind=4)
            call write_arr_to_restart(file_id, i, .false.)
         enddo
      endif

      !> problem-specific restart writes. Everything that was not written by the above write_arr_to_restart calls
      if (associated(problem_write_restart)) call problem_write_restart(file_id)

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

   end subroutine write_restart_hdf5_v1

   !----------------------------------------------------------------------------------
   ! Write fluid, mag or other variables (4-D and 3-D arrays)

   subroutine write_arr_to_restart(file_id, ind, tgt3d)

      use constants,  only: xdim, ydim, zdim, ndims, AT_OUT_B, AT_IGNORE, LONG, dsetnamelen
      use dataio_pub, only: die
      use domain,     only: is_multicg
      use gc_list,    only: cg_list_element
      use grid,       only: all_cg
      use grid_cont,  only: grid_container
      use hdf5,       only: HID_T, HSIZE_T, H5T_NATIVE_DOUBLE, h5dwrite_f, h5sclose_f, h5pclose_f, h5dclose_f, &
           &                H5P_DATASET_CREATE_F, H5S_SELECT_SET_F, H5P_DATASET_XFER_F, H5FD_MPIO_INDEPENDENT_F, H5FD_MPIO_COLLECTIVE_F, &
           &                h5screate_simple_f, h5pcreate_f, h5dcreate_f, h5dget_space_f, h5pset_dxpl_mpio_f, h5sselect_hyperslab_f

      implicit none

      integer(HID_T), intent(in)                    :: file_id   !> File identifier
      integer(kind=4), intent(in)                   :: ind       !> index of cg%q(:) or cg%w(:) array
      logical, intent(in)                           :: tgt3d     !> .true. for 3D array, .false. otherwise

      real, pointer, dimension(:,:,:)               :: pa3d     !> pointer to specified scope of pa3d
      real, pointer, dimension(:,:,:,:)             :: pa4d     !> pointer to specified scope of pa4d
      integer, parameter :: rank4 = 1 + ndims
      integer(HSIZE_T), dimension(rank4) :: dimsf, chunk_dims
      integer(HID_T) :: dset_id               !> Dataset identifier
      integer(HID_T) :: dplist_id, plist_id   !> Property list identifiers
      integer(HID_T) :: dfilespace, filespace !> Dataspace identifiers in file
      integer(HID_T) :: memspace              !> Dataspace identifier in memory

      integer,         dimension(ndims) :: area, lleft, lright, chnk
      integer(kind=8), dimension(ndims) :: loffs
      integer(HSIZE_T), dimension(rank4) :: cnt, offset, stride
      integer(kind=4) :: rank, error, area_type
      integer :: ir, dim1
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg
      character(len=dsetnamelen) ::dname

      if (tgt3d) then
         if (ind < lbound(all_cg%first%cg%q(:), dim=1) .or. ind > ubound(all_cg%first%cg%q(:), dim=1)) call die("[restart_hdf5:write_arr_to_restart] Invalid 3D array")
         dim1 = 1
         rank = ndims
         dname = all_cg%first%cg%q(ind)%name
         area_type = all_cg%first%cg%q(ind)%restart_mode
      else
         if (ind < lbound(all_cg%first%cg%w(:), dim=1) .or. ind > ubound(all_cg%first%cg%w(:), dim=1)) call die("[restart_hdf5:write_arr_to_restart] Invalid 4D array")
         dim1 = size(all_cg%first%cg%w(ind)%arr, dim=1)
         rank = rank4
         dname = all_cg%first%cg%w(ind)%name
         area_type = all_cg%first%cg%w(ind)%restart_mode
      endif
      ir = rank4 - rank + 1 ! 1 for 4-D arrays, 2 for 3-D arrays (to simplify use of count(:), offset(:), stride(:), block(:), dimsf(:) and chunk_dims(:)

      if (area_type == AT_IGNORE) return !> \todo write a list of unsaved arrays?
      call set_area_for_restart(area_type, area)

      dimsf      = [dim1, area(:)] ! Dataset dimensions
      ! Create the file space for the dataset and make it chunked if possible
      call h5screate_simple_f(rank, dimsf(ir:), dfilespace, error)
      call h5pcreate_f(H5P_DATASET_CREATE_F, dplist_id, error)

      !> \todo figure out a neater condition that allows for chunked I/O. Perhaps something computed in domain module?
      !if (.not. (is_uneven .or. (area_type == AT_OUT_B))) call h5pset_chunk_f(dplist_id, rank, chunk_dims(ir:), error)
      call h5dcreate_f(file_id, dname, H5T_NATIVE_DOUBLE, dfilespace, dset_id, error, dplist_id)

      ! Each process defines dataset in memory and writes it to the hyperslab in the file.
      call h5dget_space_f(dset_id, filespace, error)

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
         if (.not. is_multicg) call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
#ifdef INDEPENDENT_ATOUTB
      endif
#endif /* INDEPENDENT_ATOUTB */

      stride(:) = 1
      cnt(:)  = 1

      cgl => all_cg%first
      do while (associated(cgl))
         cg => cgl%cg

         call set_dims_for_restart(area_type, chnk, lleft, lright, loffs, cg)

         chunk_dims = [dim1, chnk(:)] ! Chunks dimensions
         offset(:) = [0_LONG, loffs(:)]
         call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset(ir:), cnt(ir:), error, stride(ir:), chunk_dims(ir:))

         call h5screate_simple_f(rank, chunk_dims(ir:), memspace, error)

         ! write data
         if (tgt3d) then
            pa3d => cg%q(ind)%arr(lleft(xdim):lright(xdim), lleft(ydim):lright(ydim), lleft(zdim):lright(zdim))
            if (associated(pa3d)) then
               call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, pa3d, dimsf(ir:), error, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
            else
               call die("[restart_hdf5:write_arr_to_restart] unassociated 3D array pointer")
            endif
         else
            pa4d => cg%w(ind)%arr(:, lleft(xdim):lright(xdim), lleft(ydim):lright(ydim), lleft(zdim):lright(zdim))
            if (associated(pa4d)) then
               call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, pa4d, dimsf(:), error, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
            else
               call die("[restart_hdf5:write_arr_to_restart] unassociated 4D array pointer")
            endif
         endif

         call h5sclose_f(memspace, error)

         cgl => cgl%nxt
      enddo

      call h5pclose_f(plist_id, error)
      call h5sclose_f(filespace, error)
      call h5dclose_f(dset_id, error)
      call h5pclose_f(dplist_id, error)
      call h5sclose_f(dfilespace, error)

   end subroutine write_arr_to_restart

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
!!$            if (.not. is_multicg) call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
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
! Boundary cells are exchanged with the neughbours. Corner boundary cells are not guaranteed to be correct.
! External boundary cells are not stored in the restart file and thus all of them are lost (area_type = AT_OUT_B not implemented yet).

   subroutine read_arr_from_restart(file_id, ind, tgt3d, alt_area_type, alt_name)

      use constants,    only: xdim, ydim, zdim, ndims, LONG, AT_IGNORE, dsetnamelen
      use dataio_pub,   only: msg, die
      use domain,       only: is_multicg
      use grid,         only: all_cg
      use gc_list,      only: cg_list_element
      use grid_cont,    only: grid_container
      use hdf5,         only: HID_T, HSIZE_T, SIZE_T, H5T_NATIVE_DOUBLE, h5dread_f, h5sclose_f, h5pclose_f, h5dclose_f, &
           &                  H5S_SELECT_SET_F, H5P_DATASET_XFER_F, H5FD_MPIO_COLLECTIVE_F, &
           &                  h5dopen_f, h5sget_simple_extent_ndims_f, h5dget_space_f, &
           &                  h5pcreate_f, h5pset_dxpl_mpio_f, h5sselect_hyperslab_f, h5screate_simple_f
      use external_bnd, only: arr3d_boundaries

      implicit none

      integer(HID_T),             intent(in) :: file_id            !> File identifier
      integer(kind=4),            intent(in) :: ind                !> index of cg%q(:) or cg%w(:) arrays
      logical,                    intent(in) :: tgt3d              !> .true. for 3D arrays, .false. otherwise
      integer(kind=4),  optional, intent(in) :: alt_area_type
      character(len=*), optional, intent(in) :: alt_name           !> used only in galdisk* setups

      real, pointer, dimension(:,:,:)        :: pa3d               !> pointer to specified scope of pa3d
      real, pointer, dimension(:,:,:,:)      :: pa4d               !> pointer to specified scope of pa4d
      integer, parameter                     :: rank4 = 1 + ndims
      integer(HSIZE_T), dimension(rank4)     :: dimsf, chunk_dims
      integer(HID_T)                         :: dset_id            !> Dataset identifier
      integer(HID_T)                         :: plist_id           !> Property list identifier
      integer(HID_T)                         :: filespace          !> Dataspace identifier in file
      integer(HID_T)                         :: memspace           !> Dataspace identifier in memory
      integer(HSIZE_T), dimension(rank4)     :: cnt, offset, stride
      integer,          dimension(ndims)     :: area, lleft, lright, chnk
      integer(kind=8),  dimension(ndims)     :: loffs
      integer(kind=4)                        :: rank, rankf, error, area_type
      integer                                :: ir, dim1
      type(cg_list_element), pointer         :: cgl
      type(grid_container),  pointer         :: cg
      character(len=dsetnamelen)             :: dname, cgname

      !> \deprecated Duplicated code
      if (tgt3d) then
         if (ind < lbound(all_cg%first%cg%q(:), dim=1) .or. ind > ubound(all_cg%first%cg%q(:), dim=1)) call die("[restart_hdf5:read_arr_from_restart] Invalid 3D array")
         dim1 = 1
         rank = ndims
         cgname = all_cg%first%cg%q(ind)%name
         area_type = all_cg%first%cg%q(ind)%restart_mode
      else
         if (ind < lbound(all_cg%first%cg%w(:), dim=1) .or. ind > ubound(all_cg%first%cg%w(:), dim=1)) call die("[restart_hdf5:read_arr_from_restart] Invalid 4D array")
         dim1 = size(all_cg%first%cg%w(ind)%arr, dim=1)
         rank = rank4
         cgname = all_cg%first%cg%w(ind)%name
         area_type = all_cg%first%cg%w(ind)%restart_mode
      endif
      dname = cgname
      ir = rank4 - rank + 1 ! 1 for 4-D arrays, 2 for 3-D arrays (to simplify use of count(:), offset(:), stride(:), block(:), dimsf(:) and chunk_dims(:)

      if (present(alt_area_type)) area_type = alt_area_type
      if (area_type == AT_IGNORE) return !> \todo write a list of unsaved arrays?
      call set_area_for_restart(area_type, area)

      dimsf = [dim1, area(:)]      ! Dataset dimensions

      ! Create dataset.and filespace
      if (present(alt_name)) dname = alt_name
      call h5dopen_f(file_id, dname, dset_id, error)
      if (error /= 0) then
         write(msg, '(3a)') "[restart_hdf5:read_arr_from_restart] Opening dataset '", dname,"' failed."
         call die(msg)
      endif

      call h5dget_space_f(dset_id, filespace, error)
      call h5sget_simple_extent_ndims_f (filespace, rankf, error)
      if (rank /= rankf) then
         write(msg,'(3a,2(i2,a))')"[restart_hdf5:read_arr_from_restart] Rank mismatch in array '", dname, "' (", rank, " /= ", rankf, ")"
         call die(msg)
      endif

      ! Select hyperslab in the file.
      stride(:) = 1
      cnt(:)  = 1

      cgl => all_cg%first
      do while (associated(cgl))
         cg => cgl%cg

         call set_dims_for_restart(area_type, chnk, lleft, lright, loffs, cg)

         chunk_dims = [dim1, chnk(:)] ! Chunks dimensions

         offset(:) = [ 0_LONG, loffs(:) ]
         call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset(ir:), cnt(ir:), error, stride(ir:), chunk_dims(ir:))

         call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
         if (.not. is_multicg) call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
         call h5screate_simple_f(rank, chunk_dims(ir:), memspace, error)

         ! Read the array
         if (tgt3d) then
            pa3d => cg%get_na_ptr(cgname)
            if (associated(pa3d)) then
               call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, pa3d(lleft(xdim):lright(xdim), lleft(ydim):lright(ydim), lleft(zdim):lright(zdim)), &
                    &         dimsf(ir:), error, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
            else
               write(msg,'(5a)')"Unassociated '",trim(cgname),"' 3D array while reading dataset '",trim(dname),"' from the restart file"
               call die(msg)
            endif
         else
            pa4d => cg%get_na_ptr_4d(cgname)
            if (associated(pa4d)) then
               call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, pa4d(:, lleft(xdim):lright(xdim), lleft(ydim):lright(ydim), lleft(zdim):lright(zdim)), &
                    &         dimsf(ir:), error, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
            else
               write(msg,'(5a)')"Unassociated '",trim(cgname),"' 4D array while reading dataset '",trim(dname),"' from the restart file"
               call die(msg)
            endif
         endif

         if (error /= 0) then
            write(msg, '(3a)') "[restart_hdf5:read_arr_from_restart] Reading dataset '", dname,"' failed."
            call die(msg)
         endif

         call h5sclose_f(memspace, error)

         cgl => cgl%nxt
      enddo

      call h5pclose_f(plist_id, error)
      call h5sclose_f(filespace, error)
      call h5dclose_f(dset_id, error)

      if (tgt3d) call arr3d_boundaries(cg%get_na_ind(dname), area_type=area_type)
         ! Originally the pa3d array was written with the guardcells. The internal guardcells will be exchanged but the external ones are lost.

      ! rank-4 arrays (cg%u(:,:,:,:) and b(:,:,:,:)) have their own guardcell-exchange routines, which can also be called here
      !> \todo consider also calling 4D boundaries

   end subroutine read_arr_from_restart

   subroutine read_restart_hdf5_v1(chdf)

      use common_hdf5, only: hdf
      use constants,   only: cwdlen, cbuff_len, domlen, idlen, xdim, ydim, zdim, LO, HI, I_ONE
      use dataio_pub,  only: msg, printio, warn, die, require_init_prob, problem_name, run_id, piernik_hdf5_version, fix_string
      use dataio_user, only: problem_read_restart
      use domain,      only: dom
      use fluidindex,  only: flind
      use global,      only: magic_mass, t, dt
      use grid,        only: all_cg
      use grid_cont,   only: grid_container
      use hdf5,        only: HID_T, SIZE_T, H5P_FILE_ACCESS_F, H5F_ACC_RDONLY_F, &
           &                 h5open_f, h5pcreate_f, h5pset_fapl_mpio_f, h5fopen_f, h5pclose_f, h5fclose_f, h5close_f
      use h5lt,        only: h5ltget_attribute_double_f, h5ltget_attribute_int_f, h5ltget_attribute_string_f
      use mpi,         only: MPI_CHARACTER, MPI_INTEGER, MPI_DOUBLE_PRECISION, MPI_INFO_NULL
      use mpisetup,    only: comm, ierr, comm, master, FIRST

      implicit none

      type(hdf), intent(inout) :: chdf

      integer               :: nu
      integer(kind=4)       :: i
      character(len=cwdlen) :: filename  ! File name

      integer(HID_T)        :: file_id       ! File identifier
      integer(HID_T)        :: plist_id      ! Property list identifier
      integer(SIZE_T)       :: bufsize

      integer(kind=4)       :: error
      logical               :: file_exist

      real, dimension(1)    :: rbuf
      integer(kind=4), dimension(1) :: ibuf

      real                  :: restart_hdf5_version
      type(grid_container), pointer :: fcg

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
         if (error /= 0) call die("[restart_hdf5:read_restart_hdf5_v1] Cannot read 'piernik' attribute from the restart file. The file may be either damaged or incompatible")
         if (rbuf(1) > piernik_hdf5_version) then
            write(msg,'(2(a,f5.2))')"[restart_hdf5:read_restart_hdf5_v1] Cannot read future versions of the restart file: ", rbuf(1)," > ", piernik_hdf5_version
            call die(msg)
         else if (int(rbuf(1)) < int(piernik_hdf5_version)) then
            write(msg,'(2(a,f5.2))')"[restart_hdf5:read_restart_hdf5_v1] The restart file is too ancient. It is unlikely that it could work correctly: ", rbuf(1)," << ", piernik_hdf5_version
            call die(msg)
         else if (rbuf(1) < piernik_hdf5_version) then
            write(msg,'(2(a,f5.2))')"[restart_hdf5:read_restart_hdf5_v1] Old versions of the restart file may not always work fully correctly: ", rbuf(1)," < ", piernik_hdf5_version
            call warn(msg)
         endif

         call h5ltget_attribute_int_f(file_id,"/","nxd", ibuf, error)
         if (ibuf(1) /= dom%n_d(xdim) .or. error /= 0) call die("[restart_hdf5:read_restart_hdf5_v1] nxd does not match")
         if (dom%has_dir(xdim)) then
            call h5ltget_attribute_double_f(file_id,"/","xmin", rbuf, error)
            if (rbuf(1) /= dom%edge(xdim, LO) .or. error /= 0) call die("[restart_hdf5:read_restart_hdf5_v1] xmin does not match")
            call h5ltget_attribute_double_f(file_id,"/","xmax", rbuf, error)
            if (rbuf(1) /= dom%edge(xdim, HI) .or. error /= 0) call die("[restart_hdf5:read_restart_hdf5_v1] xmax does not match")
         endif

         call h5ltget_attribute_int_f(file_id,"/","nyd", ibuf, error)
         if (ibuf(1) /= dom%n_d(ydim) .or. error /= 0) call die("[restart_hdf5:read_restart_hdf5_v1] nyd does not match")
         if (dom%has_dir(ydim)) then
            call h5ltget_attribute_double_f(file_id,"/","ymin", rbuf, error)
            if (rbuf(1) /= dom%edge(ydim, LO) .or. error /= 0) call die("[restart_hdf5:read_restart_hdf5_v1] ymin does not match")
            call h5ltget_attribute_double_f(file_id,"/","ymax", rbuf, error)
            if (rbuf(1) /= dom%edge(ydim, HI) .or. error /= 0) call die("[restart_hdf5:read_restart_hdf5_v1] ymax does not match")
         endif

         call h5ltget_attribute_int_f(file_id,"/","nzd", ibuf, error)
         if (ibuf(1) /= dom%n_d(zdim) .or. error /= 0) call die("[restart_hdf5:read_restart_hdf5_v1] nzd does not match")
         if (dom%has_dir(zdim)) then
            call h5ltget_attribute_double_f(file_id,"/","zmin", rbuf, error)
            if (rbuf(1) /= dom%edge(zdim, LO) .or. error /= 0) call die("[restart_hdf5:read_restart_hdf5_v1] zmin does not match")
            call h5ltget_attribute_double_f(file_id,"/","zmax", rbuf, error)
            if (rbuf(1) /= dom%edge(zdim, HI) .or. error /= 0) call die("[restart_hdf5:read_restart_hdf5_v1] zmax does not match")
         endif

         call h5fclose_f(file_id, error)
      endif

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
      call h5pset_fapl_mpio_f(plist_id, comm, MPI_INFO_NULL, error)

      call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, error, access_prp = plist_id)
      call h5pclose_f(plist_id, error)

      ! set up things such as register user rank-3 and rank-4 arrays to be read by read_arr_from_restart. Read also anything that is not read by all read_arr_from_restart calls
      if (associated(problem_read_restart)) call problem_read_restart(file_id)

      ! read auxiliary variables
      fcg => all_cg%first%cg
      if (allocated(fcg%q)) then
         do i = lbound(fcg%q(:), dim=1, kind=4), ubound(fcg%q(:), dim=1, kind=4)
            call read_arr_from_restart(file_id, i, .true.)
         enddo
      endif

      if (allocated(fcg%w)) then
         do i = lbound(fcg%w(:), dim=1, kind=4), ubound(fcg%w(:), dim=1, kind=4)
            call read_arr_from_restart(file_id, i, .false.)
         enddo
      endif

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

   end subroutine read_restart_hdf5_v1

!-------------------------------------------------------------------------------
! Routines for multi-file, multi-domain restarts

!>
!! \description The format of Piernik v2 restart and data files
!!
!! * A restart point or a data dump may consist of one or a multiple files
!!
!! The following data should be present in all files in identical copies:
!! * file version
!! * physical time
!! * number of timestep
!! * problem.par
!! * env
!! * base domain description
!!   * grid size
!!   * physical size
!!   * boundary condition types
!! * refinement topology description:
!!   * a list of all patches (boxes)
!!     * size
!!     * offset
!!     * refinement level
!! * some other attributes that are present in v1 restart files
!!
!! The following data should be unique for each file:
!! * list of grid containers
!!
!! A grid container contains only essential, non-redundant data and consists of:
!! * size
!! * offset
!! * refinement level
!! * rank-4 arrays (restart only)
!!   * "fluid"
!!   * "mag"
!!   * other cg%w arrays marked for restart
!! * rank-3 arrays
!!   * cg%q arrays marked for restart (restart only)
!!   * named quantities for data dumps
!!
!! \todo Check if it is possible to filter the data through shuffle and gzip -9  during write
!<

!>
!! \brief Write a multi-file, multi-domain restart file
!!
!! \warning Not implemented yet
!<

   subroutine write_restart_hdf5_v2

      use dataio_pub, only: die

      implicit none

      call die("[restart_hdf5:write_restart_hdf5_v2] Not implemented yet")

   end subroutine write_restart_hdf5_v2

!>
!! \brief Read a multi-file, multi-domain restart file
!!
!! \warning Not implemented yet
!!
!! \details Reading a v2 restart file should proceed as follows
!! * Identify the restart point file(s) on master
!! * Send the identification to the slaves
!! * On each processor check what files are present if any
!! * Check consistency of the data that should be present in all files in identical copies
!! * Prepare list of available grid containers
!! * If the list does not fill completely the boxes listed in refinement topology description on any process, communicate lists to the master
!! * If the compiled list does not fill completely the boxes listed in refinement topology description then call die()
!! * Compute domain decomposition for current refinement topology description (different number of processors and different decomposition method is allowed)
!! * Calculate on master which parts (if any) cannot be read directly and should be communicated, send the lists to the involved slaves
!! * Read the data (own and requested by others) from the file(s)
!! * Receive anything that had to be read on onther processes
!! * Check for local completness
!!
!! Optionally:
!! * When a single file is visible by multiple processes and the fileststem does not tolerate massively concurrent I/O load,
!!   choose one or few process to do the I/O and make the others to wait for the data.
!! * It should be possible to read runtime parameters from data stored in restart point rather than from actual problem.par file(s).
!!   A way to override some of them (in most cases  END_CONTROL::tend and nend) should be provided.
!<

   subroutine read_restart_hdf5_v2(status_v2)

      use constants,  only: INVALID
      use dataio_pub, only: warn

      implicit none

      integer, intent(out) :: status_v2

      call warn("[restart_hdf5:read_restart_hdf5_v2] Not implemented yet")

      status_v2 = INVALID

   end subroutine read_restart_hdf5_v2

end module restart_hdf5
