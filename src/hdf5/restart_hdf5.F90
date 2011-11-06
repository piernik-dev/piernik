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
   use constants,   only: dsetnamelen

   implicit none

   private
   public :: read_restart_hdf5, write_restart_hdf5, read_arr_from_restart

   integer,                      parameter :: STAT_OK = 0
   character(len=dsetnamelen/2), parameter :: cg_gname = "cg" ! leave the other half of dsetnamelen for id number

!> \brief Add an attribute to the given group and initialize its value
   interface create_attribute
      module procedure create_int_attribute
      module procedure create_real_attribute
   end interface

!> \brief Compare a 1D array across all processed. Die if any difference is found
   interface compare_array1D
      module procedure compare_int_array1D
      module procedure compare_real_array1D
   end interface

contains

!>
!! \brief Wrapper routine that writes a version 2.x or version 1.x restart file, depending on use_v2_io switch
!<

   subroutine write_restart_hdf5

      use constants,   only: I_ONE, cwdlen, WR
      use common_hdf5, only: set_common_attributes
      use dataio_pub,  only: msg, printio, printinfo, tmr_hdf, thdf, use_v2_io, nres
      use mpisetup,    only: comm, ierr, master
      use timer,       only: set_timer

      implicit none
      character(len=cwdlen)         :: filename  ! File name

      thdf = set_timer(tmr_hdf,.true.)

      filename = restart_fname(WR)
      if (master) then
         write(msg,'(3a)') 'Writing restart ', trim(filename), " ... "
         call printio(msg, .true.)
      endif
      call set_common_attributes(filename)

      if (use_v2_io) then
         call write_restart_hdf5_v2(filename)
      else
         call write_restart_hdf5_v1(filename)
      endif
      call MPI_Barrier(comm, ierr)

      thdf = set_timer(tmr_hdf)
      if (master) then
         write(msg,'(a6,f10.2,a2)') ' done ', thdf, ' s'
         call printinfo(msg, .true.)
      endif

      nres = nres + I_ONE

   end subroutine write_restart_hdf5

!>
!! \brief Wrapper routine that reads a version 2.x or version 1.x restart file, depending on use_v2_io switch. On v2 failure it falls back to v1.
!<

   subroutine read_restart_hdf5

      use dataio_pub, only: use_v2_io

      implicit none

      integer :: status_v2

      if (use_v2_io) call read_restart_hdf5_v2(status_v2)
      if (status_v2 /= STAT_OK .or. .not. use_v2_io) call read_restart_hdf5_v1

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
!! \brief Construct a filename for the restart file
!!
!! \todo add a part number for multifile restarts
!! \todo sanitize problem_name and run_id
!<

   function restart_fname(wr_rd) result(filename)

      use constants,   only: cwdlen, RD, WR
      use dataio_pub,  only: problem_name, run_id, nres, wd_wr, wd_rd
      use mpi,         only: MPI_CHARACTER
      use mpisetup,    only: comm, ierr, master, FIRST

      implicit none

      character(len=cwdlen) :: filename, temp  ! File name
      integer(kind=4), intent(in)   :: wr_rd

      if (master) then
         write(temp,'(a,a1,a3,a1,i4.4,a4)') trim(problem_name),'_', run_id,'_', nres,'.res'
         select case (wr_rd)
            case (RD)
               write(filename,'(2a)') trim(wd_rd),trim(temp)
            case (WR)
               write(filename,'(2a)') trim(wd_wr),trim(temp)
            case default
               write(filename,'(2a)') './',trim(temp)
         end select
      endif

      call MPI_Bcast(filename, cwdlen, MPI_CHARACTER, FIRST, comm, ierr)

   end function restart_fname

!> \brief Write restart dump v1.x

   subroutine write_restart_hdf5_v1(filename)

      use constants,   only: cwdlen
      use dataio_user, only: problem_write_restart
      use grid,        only: all_cg
      use grid_cont,   only: grid_container
      use hdf5,        only: HID_T, H5P_FILE_ACCESS_F, H5F_ACC_RDWR_F, h5open_f, h5close_f, h5fopen_f, h5fclose_f, h5pcreate_f, h5pclose_f, h5pset_fapl_mpio_f
      !, H5P_DATASET_XFER_F, h5pset_preserve_f
      use mpi,         only: MPI_INFO_NULL
      use mpisetup,    only: comm, FIRST

      implicit none
      character(len=cwdlen), intent(in) :: filename      !> HDF File name

      integer(kind=4)               :: i
      integer(HID_T)                :: file_id       !> File identifier
      integer(HID_T)                :: plist_id      !> Property list identifier
      integer(kind=4)               :: error
      type(grid_container), pointer :: fcg

      ! Set up a new HDF5 file for parallel write
      call h5open_f(error)
      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
      call h5pset_fapl_mpio_f(plist_id, comm, MPI_INFO_NULL, error)
      call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error, access_prp = plist_id)

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

      ! problem-specific restart writes. Everything that was not written by the above write_arr_to_restart calls
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

   end subroutine write_restart_hdf5_v1

!> \brief Write fluid, mag or other variables (4-D and 3-D arrays) to v1 restart dump

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

      integer(HID_T),  intent(in)        :: file_id                   !> File identifier
      integer(kind=4), intent(in)        :: ind                       !> index of cg%q(:) or cg%w(:) array
      logical,         intent(in)        :: tgt3d                     !> .true. for 3D array, .false. otherwise

      real, pointer, dimension(:,:,:)    :: pa3d                      !> pointer to specified scope of pa3d
      real, pointer, dimension(:,:,:,:)  :: pa4d                      !> pointer to specified scope of pa4d
      integer, parameter                 :: rank4 = 1 + ndims
      integer(HSIZE_T), dimension(rank4) :: dimsf, chunk_dims
      integer(HID_T)                     :: dset_id                   !> Dataset identifier
      integer(HID_T)                     :: dplist_id, plist_id       !> Property list identifiers
      integer(HID_T)                     :: dfilespace, filespace     !> Dataspace identifiers in file
      integer(HID_T)                     :: memspace                  !> Dataspace identifier in memory

      integer,         dimension(ndims)  :: area, lleft, lright, chnk
      integer(kind=8), dimension(ndims)  :: loffs
      integer(HSIZE_T), dimension(rank4) :: cnt, offset, stride
      integer(kind=4)                    :: rank, error, area_type
      integer                            :: ir, dim1
      type(cg_list_element), pointer     :: cgl
      type(grid_container),  pointer     :: cg
      character(len=dsetnamelen)         :: dname

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

!>
!! \brief Read interior cells of a specified array from the restart file.
!!
!! \details  Boundary cells are exchanged with the neughbours. Corner boundary cells are not guaranteed to be correct.
!<

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
            pa3d => cg%q(cg%get_na_ind(cgname))%arr(lleft(xdim):lright(xdim), lleft(ydim):lright(ydim), lleft(zdim):lright(zdim))
            if (associated(pa3d)) then
               call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, pa3d, dimsf(ir:), error, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
            else
               write(msg,'(5a)')"Unassociated '",trim(cgname),"' 3D array while reading dataset '",trim(dname),"' from the restart file"
               call die(msg)
            endif
         else
            pa4d => cg%w(cg%get_na_ind_4d(cgname))%arr(:, lleft(xdim):lright(xdim), lleft(ydim):lright(ydim), lleft(zdim):lright(zdim))
            if (associated(pa4d)) then
               call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, pa4d, dimsf(ir:), error, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
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

!> \brief Read a v1 restart point.

   subroutine read_restart_hdf5_v1

      use constants,   only: cwdlen, cbuff_len, domlen, idlen, xdim, ydim, zdim, LO, HI, I_ONE, RD
      use dataio_pub,  only: msg, warn, die, printio, require_init_prob, problem_name, run_id, piernik_hdf5_version, fix_string, &
           &                 domain_dump, last_hdf_time, next_t_log, next_t_tsl, nhdf, nres, step_hdf, new_id, step_res
      use dataio_user, only: problem_read_restart
      use domain,      only: dom
      use fluidindex,  only: flind
      use global,      only: magic_mass, t, dt, nstep
      use grid,        only: all_cg
      use grid_cont,   only: grid_container
      use hdf5,        only: HID_T, SIZE_T, H5P_FILE_ACCESS_F, H5F_ACC_RDONLY_F, &
           &                 h5open_f, h5pcreate_f, h5pset_fapl_mpio_f, h5fopen_f, h5pclose_f, h5fclose_f, h5close_f
      use h5lt,        only: h5ltget_attribute_double_f, h5ltget_attribute_int_f, h5ltget_attribute_string_f
      use mpi,         only: MPI_CHARACTER, MPI_INTEGER, MPI_DOUBLE_PRECISION, MPI_INFO_NULL
      use mpisetup,    only: comm, ierr, master, FIRST

      implicit none

      integer                       :: nu
      integer(kind=4)               :: i
      character(len=cwdlen)         :: filename      !> File name

      integer(HID_T)                :: file_id       !> File identifier
      integer(HID_T)                :: plist_id      !> Property list identifier
      integer(SIZE_T)               :: bufsize

      integer(kind=4)               :: error
      logical                       :: file_exist

      real,            dimension(1) :: rbuf
      integer(kind=4), dimension(1) :: ibuf

      real                          :: restart_hdf5_version
      type(grid_container), pointer :: fcg

      nu = flind%all

      filename = restart_fname(RD)

      if (master) then
         write(msg, '(2a)') 'Reading restart file: ', trim(filename)
         call printio(msg)
      endif

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
         nstep = ibuf(1)
         call h5ltget_attribute_int_f(file_id,"/","nres", ibuf, error)
         nres = ibuf(1)
         call h5ltget_attribute_int_f(file_id,"/","nhdf", ibuf, error)
         nhdf = ibuf(1)
         call h5ltget_attribute_int_f(file_id,"/","step_res", ibuf, error)
         step_res = ibuf(1)
         call h5ltget_attribute_int_f(file_id,"/","step_hdf", ibuf, error)
         step_hdf = ibuf(1)
         call h5ltget_attribute_double_f(file_id,"/","next_t_tsl", rbuf, error)
         next_t_tsl = rbuf(1)
         call h5ltget_attribute_double_f(file_id,"/","next_t_log", rbuf, error)
         next_t_log = rbuf(1)
         call h5ltget_attribute_double_f(file_id,"/","last_hdf_time", rbuf, error)
         last_hdf_time = rbuf(1)

         call h5ltget_attribute_string_f(file_id,"/","problem_name", problem_name, error)
         call h5ltget_attribute_string_f(file_id,"/","domain", domain_dump, error)
         call h5ltget_attribute_string_f(file_id,"/","run_id", new_id, error)

         if (restart_hdf5_version > 1.11) then
            call h5ltget_attribute_int_f(file_id,"/","require_init_prob", ibuf, error)
            require_init_prob = ibuf(1)
         endif

         problem_name = fix_string(problem_name)   !> \deprecated BEWARE: >=HDF5-1.8.4 has weird issues with strings
         new_id  = fix_string(new_id)    !> \deprecated   this bit hacks it around
         domain_dump  = fix_string(domain_dump)

         call h5fclose_f(file_id, error)

         write(msg,'(2a)') 'Done reading restart file: ', trim(filename)
         call printio(msg)
      endif
      call h5close_f(error)

      call MPI_Bcast(restart_hdf5_version,    I_ONE, MPI_DOUBLE_PRECISION, FIRST, comm, ierr)

      call MPI_Bcast(nstep,    I_ONE, MPI_INTEGER, FIRST, comm, ierr)
      call MPI_Bcast(nres,     I_ONE, MPI_INTEGER, FIRST, comm, ierr)
      call MPI_Bcast(nhdf,     I_ONE, MPI_INTEGER, FIRST, comm, ierr)
      call MPI_Bcast(step_res, I_ONE, MPI_INTEGER, FIRST, comm, ierr)
      call MPI_Bcast(step_hdf, I_ONE, MPI_INTEGER, FIRST, comm, ierr)
      if (restart_hdf5_version > 1.11) call MPI_Bcast(require_init_prob, I_ONE, MPI_INTEGER, FIRST, comm, ierr)

      call MPI_Bcast(next_t_tsl,    I_ONE, MPI_DOUBLE_PRECISION, FIRST, comm, ierr)
      call MPI_Bcast(next_t_log,    I_ONE, MPI_DOUBLE_PRECISION, FIRST, comm, ierr)
      call MPI_Bcast(last_hdf_time, I_ONE, MPI_DOUBLE_PRECISION, FIRST, comm, ierr)
      call MPI_Bcast(t,             I_ONE, MPI_DOUBLE_PRECISION, FIRST, comm, ierr)
      call MPI_Bcast(dt,            I_ONE, MPI_DOUBLE_PRECISION, FIRST, comm, ierr)

      call MPI_Bcast(problem_name, cbuff_len, MPI_CHARACTER, FIRST, comm, ierr)
      call MPI_Bcast(domain_dump,  domlen,    MPI_CHARACTER, FIRST, comm, ierr)
      call MPI_Bcast(new_id,       idlen,     MPI_CHARACTER, FIRST, comm, ierr)

   end subroutine read_restart_hdf5_v1

!-------------------------------------------------------------------------------
! Routines for multi-file, multi-domain restarts

!>
!! \note The format of Piernik v2 restart and data files
!!
!! A restart point or a data dump may consist of one or a multiple files
!!
!! The following data should be present in all files in identical copies:
!! - file version
!! - physical time
!! - number of timestep
!! - problem.par
!! - env
!! - base domain description
!!   - grid size
!!   - physical size
!!   - boundary condition types
!! - refinement topology description:
!!   - a list of all patches (boxes)
!!     - size
!!     - offset
!!     - refinement level
!! - some other attributes that are present in v1 restart files
!! - user attributes
!!
!! The following data should be unique for each file:
!! - list of grid containers
!!
!! A grid container contains only essential, non-redundant data and consists of:
!! - size
!! - offset
!! - refinement level
!! - rank-4 arrays (restart only)
!!   - "fluid"
!!   - "mag"
!!   - other cg%w arrays marked for restart
!! - rank-3 arrays
!!   - cg%q arrays marked for restart (restart only)
!!   - named quantities for data dumps
!!
!! At the moment we assume that all problem-specific data is either registeres in cg%q and cg%w arrays or can be written to global attributes.
!!
!! It will be possible to provide support for problem-specific data that should be appended to grid containers when anyone will need it
!!
!! \todo Check if it is possible to filter the data through shuffle and gzip -9  during write
!<

!> \brief Generate numbered cg group name
   function n_cg_name(g)
      implicit none
      integer, intent(in)        :: g !< group number
      character(len=dsetnamelen) :: n_cg_name
      write(n_cg_name,'(2a,i8.8)')trim(cg_gname), "_", g
   end function n_cg_name

!>
!! \brief Write a multi-file, multi-domain restart file
!!
!! \details There are three approaches to be implemented:
!! - Single-file, serial I/O. The easiest way. Master writes everything, slaves send their data to the master. Does not take advantage of parallel filesystems. Best choice for non-parallel filesystems.
!! - Multi-file, serial I/O. An extension of the above approach. Selected processes (can be all) write to their files, other processes send them their data.
!!   Can take advantage of parallel filesystems. Can use local scratch filesystems. Requires additional work on reading.
!! - Single-file, parallel I/O. The most ambitious approach. Selected processes (can be all) write to the files, other processes send them their data.
!!   Can take advantage of parallel filesystems. Currently does not allow for compression during write.
!!   Requires a lot of pseudo-collective operations. The "flexible PHDF5" would simplify the code, but it needs to be implemented.first.
!!
!! \warning Partial implementation: Single-file, serial I/O works for non-AMR setups.
!<

   subroutine write_restart_hdf5_v2(filename)

      use constants,   only: cwdlen, dsetnamelen, singlechar, xdim, zdim, ndims, I_ONE, I_TWO, AT_IGNORE, INT4
      use dataio_pub,  only: die, nproc_io, can_i_write
      use dataio_user, only: problem_write_restart
      use domain,      only: dom
      use gc_list,     only: cg_list_element
      use grid,        only: all_cg
      use grid_cont,   only: grid_container
      use hdf5,        only: HID_T, HSIZE_T, H5F_ACC_RDWR_F, H5P_FILE_ACCESS_F, H5Z_FILTER_DEFLATE_F, &
           &                 h5open_f, h5close_f, h5fopen_f, h5fclose_f, h5gcreate_f, h5gopen_f, h5gclose_f, &
           &                 h5pcreate_f, h5pclose_f, h5pset_fapl_mpio_f, h5zfilter_avail_f
      use mpi,         only: MPI_INFO_NULL, MPI_INTEGER, MPI_INTEGER8, MPI_STATUS_IGNORE
      use mpisetup,    only: comm, proc, FIRST, LAST, master

      implicit none

      character(len=cwdlen), intent(in)             :: filename
      integer(HID_T)                                :: file_id                                  !> File identifier
      integer(HID_T)                                :: plist_id                                 !> Property list identifier
      integer(HID_T)                                :: cgl_g_id,  cg_g_id                       !> cg list and cg group identifiers
      integer(HID_T)                                :: doml_g_id, dom_g_id                      !> domain list and domain group identifiers
      integer(kind=4)                               :: error, cg_cnt
      integer                                       :: g, p, i
      integer(kind=4)                               :: drank
      integer, parameter                            :: tag = I_ONE
      integer(HSIZE_T), dimension(:),   allocatable :: ddims
      integer(kind=4),  dimension(:),   pointer     :: cg_n                                     !> offset for cg group numbering
      integer(kind=4),  dimension(:,:), pointer     :: cg_all_n_b                               !> sizes of all cg
      integer(kind=4),  dimension(:),   allocatable :: cg_rl                                    !> list of refinement levels from all cgs/procs
      integer(kind=4),  dimension(:,:), allocatable :: cg_n_b                                   !> list of n_b from all cgs/procs
      integer(kind=8),  dimension(:,:), allocatable :: cg_off                                   !> list of offsets from all cgs/procs
      type(cg_list_element), pointer                :: cgl
      type(grid_container),  pointer                :: fcg
      logical(kind=4)                               :: Z_avail                                  !> .true. if HDF5 was compiled with zlib support
      character(len=singlechar), dimension(ndims), parameter :: d_pref = [ "x", "y", "z" ]
      character(len=dsetnamelen)                    :: d_label

      ! Create a new file and initialize it

      ! Prepare groups and datasets for grid containers on the master
      allocate(cg_n(FIRST:LAST))
      call MPI_Allgather(all_cg%cnt, I_ONE, MPI_INTEGER, cg_n, I_ONE, MPI_INTEGER, comm, error)
      cg_cnt = sum(cg_n(:))
      allocate(cg_all_n_b(cg_cnt, ndims))

      if (master) then

         ! Open the HDF5 file only in master process and create all groups required for cg storage.
         ! Create also all related datasets and attributes. Do not write big datasets yet.

         call h5open_f(error)
         call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)

         call h5gcreate_f(file_id, cg_gname, cgl_g_id, error)

         call create_attribute(cgl_g_id, "cg_count", [ cg_cnt ])

         Z_avail = .false.
         if (nproc_io == 1) call h5zfilter_avail_f(H5Z_FILTER_DEFLATE_F, Z_avail, error)
         !> \todo test it thoroughly before enabling for > 1

         ! Do not assume that the master knows all the lists
         do p = FIRST, LAST
            allocate(cg_rl(cg_n(p)), cg_n_b(cg_n(p), ndims), cg_off(cg_n(p), ndims))
            if (p == FIRST) then
               g = 1
               cgl => all_cg%first
               do while (associated(cgl))
                  cg_rl(g)     = 1  !> \deprecated change this once we introduce many levels
                  cg_n_b(g, :) = cgl%cg%n_b(:)
                  cg_off(g, :) = cgl%cg%off(:)
                  g = g + 1
                  cgl => cgl%nxt
               enddo
            else
               call MPI_Recv(cg_rl,  size(cg_rl),  MPI_INTEGER,  p, tag,       comm, MPI_STATUS_IGNORE, error)
               call MPI_Recv(cg_n_b, size(cg_n_b), MPI_INTEGER,  p, tag+I_ONE, comm, MPI_STATUS_IGNORE, error)
               call MPI_Recv(cg_off, size(cg_off), MPI_INTEGER8, p, tag+I_TWO, comm, MPI_STATUS_IGNORE, error)
            endif

            do g = 1, cg_n(p)
               call h5gcreate_f(cgl_g_id, n_cg_name(sum(cg_n(:p))-cg_n(p)+g), cg_g_id, error)

               call create_attribute(cg_g_id, "level", [ cg_rl(g) ] )
               call create_attribute(cg_g_id, "n_b", cg_n_b(g, :))
               call create_attribute(cg_g_id, "off", int(cg_off(g, :), kind=4))

               cg_all_n_b(sum(cg_n(:p))-cg_n(p)+g, :) = cg_n_b(g, :)

               if (any(cg_off(g, :) > 2.**31)) call die("[restart_hdf5:write_restart_hdf5_v2] large offsets require better treatment")

               !! -- this bit is restart specific
               fcg  => all_cg%first%cg

               if (allocated(fcg%q)) then
                  drank = ndims
                  allocate(ddims(drank))
                  do i = lbound(fcg%q(:), dim=1, kind=4), ubound(fcg%q(:), dim=1, kind=4)
                     if (fcg%q(i)%restart_mode /= AT_IGNORE) &
                          call create_empty_cg_dataset(cg_g_id, fcg%q(i)%name, int(cg_n_b(g, :), kind=HSIZE_T), Z_avail)
                  enddo
                  deallocate(ddims)
               endif

               if (allocated(fcg%w)) then
                  drank = ndims + I_ONE
                  allocate(ddims(drank))
                  do i = lbound(fcg%w(:), dim=1, kind=4), ubound(fcg%w(:), dim=1, kind=4)
                     if (fcg%w(i)%restart_mode /= AT_IGNORE) &
                          call create_empty_cg_dataset(cg_g_id, fcg%w(i)%name, int([ size(fcg%w(i)%arr, dim=1, kind=4), cg_n_b(g, :) ], kind=HSIZE_T), Z_avail)
                  enddo
                  deallocate(ddims)
               endif
               !! -------------------------------

               call h5gclose_f(cg_g_id, error)
            enddo

            deallocate(cg_rl, cg_n_b, cg_off)
         enddo

         call h5gclose_f(cgl_g_id, error)

         ! describe_domains
         call h5gcreate_f(file_id, "domains", doml_g_id, error)

         call h5gcreate_f(doml_g_id, "base", dom_g_id, error)
         call create_attribute(dom_g_id, "n_d", dom%n_d(:))
         do i = xdim, zdim
            write(d_label, '(2a)') d_pref(i), "-edge_position"
            call create_attribute(dom_g_id, d_label, dom%edge(i, :))
            write(d_label, '(2a)') d_pref(i), "-boundary_type"
            call create_attribute(dom_g_id, d_label, int(dom%bnd(i, :), kind=4))
         enddo

         call h5gclose_f(dom_g_id, error)

         call create_attribute(cgl_g_id, "fine_count", [ 0_INT4 ]) ! we have only base domain at the moment

         !> \todo add here all fine domains
         ! name "fine_00000001"
         ! attributes: n_d(:), off(:), refinement
         call h5gclose_f(doml_g_id, error)

         call h5fclose_f(file_id, error)
         call h5close_f(error)

      else ! send all the necessary information to the master
         allocate(cg_rl(all_cg%cnt), cg_n_b(all_cg%cnt, ndims), cg_off(all_cg%cnt, ndims))
         g = 1
         cgl => all_cg%first
         do while (associated(cgl))
            cg_rl(g)     = 1  !> \deprecated change this once we introduce many levels
            cg_n_b(g, :) = cgl%cg%n_b(:)
            cg_off(g, :) = cgl%cg%off(:)
            g = g + 1
            cgl => cgl%nxt
         enddo
         call MPI_Send(cg_rl,  size(cg_rl),  MPI_INTEGER,  FIRST, tag,       comm, error)
         call MPI_Send(cg_n_b, size(cg_n_b), MPI_INTEGER,  FIRST, tag+I_ONE, comm, error)
         call MPI_Send(cg_off, size(cg_off), MPI_INTEGER8, FIRST, tag+I_TWO, comm, error)
         deallocate(cg_rl, cg_n_b, cg_off)
      endif

      call MPI_Bcast(cg_all_n_b, size(cg_all_n_b), MPI_INTEGER, FIRST, comm, error)

      call MPI_Barrier(comm, error)
      ! Reopen the HDF5 file for parallel write
      call h5open_f(error)
      if (can_i_write) then
         if (nproc_io > 1) then
            call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
            ! when nproc_io < nproc we'll probably need another communicator for subset of processes that have can_i_write flag set
            call h5pset_fapl_mpio_f(plist_id, comm, MPI_INFO_NULL, error)
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error, access_prp = plist_id)
            call h5pclose_f(plist_id, error)
         else
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)
         endif
         call h5gopen_f(file_id, cg_gname, cgl_g_id, error)
      endif

      !! -- this bit is restart specific
      call write_cg_to_restart(cgl_g_id, cg_n(:), cg_all_n_b(:,:))
      !! -------------------------------

      if (can_i_write) then
         call h5gclose_f(cgl_g_id, error)
         if (associated(problem_write_restart)) call problem_write_restart(file_id) ! may be called by not all processes
         call h5fclose_f(file_id, error)  ! Close the file
      endif

      call h5close_f(error)            ! Close HDF5 stuff

      deallocate(cg_n, cg_all_n_b)

   end subroutine write_restart_hdf5_v2

!> \brief Create an empty double precision dataset of given dimensions. Use compression if available.

   subroutine create_empty_cg_dataset(cg_g_id, name, ddims, Z_avail)

     use constants,  only: I_NINE
     use dataio_pub, only: enable_compression
     use hdf5,       only: HID_T, HSIZE_T, H5P_DATASET_CREATE_F, H5T_NATIVE_DOUBLE, &
          &                h5dcreate_f, h5dclose_f, h5screate_simple_f, h5sclose_f, h5pcreate_f, h5pclose_f, h5pset_deflate_f, h5pset_shuffle_f, h5pset_chunk_f

     implicit none

     integer(HID_T), intent(in)                 :: cg_g_id !< group id where to create the dataset
     character(len=*), intent(in)               :: name    !< name
     integer(HSIZE_T), dimension(:), intent(in) :: ddims   !< dimensionality
     logical(kind=4), intent(in)                        :: Z_avail !< can use compression?

     integer(HID_T)  :: prp_id, filespace, dset_id
     integer(kind=4) :: error

     call h5pcreate_f(H5P_DATASET_CREATE_F, prp_id, error)
     if (enable_compression .and. Z_avail) then
        call h5pset_shuffle_f(prp_id, error)
        call h5pset_deflate_f(prp_id, I_NINE, error)
        call h5pset_chunk_f(prp_id, size(ddims, kind=4), ddims, error)
     endif

     call h5screate_simple_f(size(ddims, kind=4), ddims, filespace, error)
     call h5dcreate_f(cg_g_id, name, H5T_NATIVE_DOUBLE, filespace, dset_id, error, dcpl_id = prp_id)
     call h5dclose_f(dset_id, error)
     call h5sclose_f(filespace, error)
     call h5pclose_f(prp_id, error)

   end subroutine create_empty_cg_dataset

!> \brief Attach an 32-bit integer attribute (scalar or rank-1 small array) to the given group.

   subroutine create_int_attribute(g_id, name, int_array)

     use constants, only: I_ONE
     use hdf5,      only: H5T_NATIVE_INTEGER, HID_T, HSIZE_T, &
          &               h5acreate_f, h5aclose_f, h5awrite_f, h5screate_simple_f, h5sclose_f

     implicit none

     integer(HID_T), intent(in)                :: g_id      !< group id where to create the attribute
     character(len=*), intent(in)              :: name      !< name
     integer(kind=4), dimension(:), intent(in) :: int_array !< the data

     integer(HID_T)  :: aspace_id, attr_id
     integer(kind=4) :: error

     call h5screate_simple_f(I_ONE, [ size(int_array, kind=HSIZE_T) ], aspace_id, error)
     call h5acreate_f(g_id, name, H5T_NATIVE_INTEGER, aspace_id, attr_id, error)
     call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, int_array, [ size(int_array, kind=HSIZE_T) ], error)
     call h5aclose_f(attr_id, error)
     call h5sclose_f(aspace_id, error)

   end subroutine create_int_attribute

!> \brief Attach an 64-bit real attribute (scalar or rank-1 small array) to the given group.

   subroutine create_real_attribute(g_id, name, real_array)

     use constants, only: I_ONE
     use hdf5,      only: H5T_NATIVE_DOUBLE, HID_T, HSIZE_T, &
          &               h5acreate_f, h5aclose_f, h5awrite_f, h5screate_simple_f, h5sclose_f

     implicit none

     integer(HID_T), intent(in)     :: g_id       !< group id where to create the attribute
     character(len=*), intent(in)   :: name       !< name
     real, dimension(:), intent(in) :: real_array !< the data

     integer(HID_T)  :: aspace_id, attr_id
     integer(kind=4) :: error

     call h5screate_simple_f(I_ONE, [ size(real_array, kind=HSIZE_T) ], aspace_id, error)
     call h5acreate_f(g_id, name, H5T_NATIVE_DOUBLE, aspace_id, attr_id, error)
     call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, real_array, [ size(real_array, kind=HSIZE_T) ], error)
     call h5aclose_f(attr_id, error)
     call h5sclose_f(aspace_id, error)

   end subroutine create_real_attribute

!> \brief Write all grid containers to the file

   subroutine write_cg_to_restart(cgl_g_id, cg_n, cg_all_n_b)

      use constants,  only: xdim, ydim, zdim, ndims, AT_IGNORE
      use dataio_pub, only: die, nproc_io, can_i_write
      use grid,       only: all_cg
      use grid_cont,  only: grid_container
      use hdf5,       only: HID_T, HSIZE_T, H5P_DATASET_XFER_F, H5FD_MPIO_INDEPENDENT_F, H5T_NATIVE_DOUBLE, &
           &                h5dopen_f, h5dclose_f, h5dwrite_f, h5gopen_f, h5gclose_f, &
           &                h5pcreate_f, h5pclose_f, h5pset_dxpl_mpio_f
      use mpi,        only: MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE
      use mpisetup,   only: master, nproc, FIRST, LAST, proc, comm

      implicit none

      integer(HID_T),                           intent(in) :: cgl_g_id    !> cg group identifier
      integer(kind=4), dimension(:),   pointer, intent(in) :: cg_n        !> offset for cg group numbering
      integer(kind=4), dimension(:,:), pointer, intent(in) :: cg_all_n_b  !> all cg sizes

      integer(HID_T)                              :: cg_g_id            !> cg group identifier
      integer(HID_T)                              :: dset_id, plist_id
      integer(kind=4)                             :: error
      integer(HSIZE_T), dimension(:), allocatable :: dims
      real, pointer,    dimension(:,:,:)          :: pa3d
      real, pointer,    dimension(:,:,:,:)        :: pa4d
      integer                                     :: i, ncg
      type(grid_container), pointer               :: cg
      integer, allocatable, dimension(:)          :: cg_src_p, cg_src_n
      integer, allocatable, dimension(:)          :: q_lst, w_lst

      ! construct source addresses of the cg to be written
      allocate(cg_src_p(1:sum(cg_n(:))), cg_src_n(1:sum(cg_n(:))))
      do i = FIRST, LAST
         cg_src_p(sum(cg_n(:i))-cg_n(i)+1:sum(cg_n(:i))) = i
         do ncg = 1, cg_n(i)
            cg_src_n(sum(cg_n(:i))-cg_n(i)+ncg) = ncg
         enddo
      enddo

      !> \todo Do a consistency check

      ! find out which fields need to be written
      if (allocated(all_cg%first%cg%q)) then
         do i = lbound(all_cg%first%cg%q(:), dim=1, kind=4), ubound(all_cg%first%cg%q(:), dim=1, kind=4)
            if (all_cg%first%cg%q(i)%restart_mode /= AT_IGNORE) call append_int_to_array(q_lst, i)
         enddo
      endif
      if (.not.allocated(q_lst)) allocate(q_lst(0))  ! without it intrinsics like size, ubound, lbound return bogus values

      if (allocated(all_cg%first%cg%w)) then
         do i = lbound(all_cg%first%cg%w(:), dim=1, kind=4), ubound(all_cg%first%cg%w(:), dim=1, kind=4)
            if (all_cg%first%cg%w(i)%restart_mode /= AT_IGNORE) call append_int_to_array(w_lst, i)
         enddo
      endif
      if (.not.allocated(w_lst)) allocate(w_lst(0))  ! without it intrinsics like size, ubound, lbound return bogus values

      ! write all cg, one by one
      do ncg = 1, sum(cg_n(:))

         if (can_i_write) then
            call h5gopen_f(cgl_g_id, n_cg_name(ncg), cg_g_id, error)
            call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
         endif

         if (nproc_io > 1) call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)    ! \todo move property to problem.par

         if (nproc_io == 1) then ! perform serial write
            if (master) then
               if (.not. can_i_write) call die("[restart_hdf5:write_cg_to_restart] Master can't write")

               if (size(q_lst) > 0) then
                  allocate(dims(ndims))
                  do i = lbound(q_lst, dim=1), ubound(q_lst, dim=1)
                     call h5dopen_f(cg_g_id, all_cg%first%cg%q(q_lst(i))%name, dset_id, error)
                     if (cg_src_p(ncg) == proc) then
                        cg => get_nth_cg(cg_src_n(ncg))
                        pa3d => cg%q(q_lst(i))%arr(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) !< \todo use set_dims_for_restart
                        dims(:) = cg%n_b
                     else
                        allocate(pa3d(cg_all_n_b(ncg, xdim), cg_all_n_b(ncg, ydim), cg_all_n_b(ncg, zdim)))
                        call MPI_Recv(pa3d(:,:,:), size(pa3d(:,:,:)), MPI_DOUBLE_PRECISION, cg_src_p(ncg), ncg + sum(cg_n(:))*i, comm, MPI_STATUS_IGNORE, error)
                        dims(:) = shape(pa3d)
                     endif
                     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, pa3d, dims, error, xfer_prp = plist_id)
                     call h5dclose_f(dset_id, error)
                     if (cg_src_p(ncg) /= proc) deallocate(pa3d)
                  enddo
                  deallocate(dims)
               endif

               if (size(w_lst) > 0) then
                  allocate(dims(ndims+1))
                  do i = lbound(w_lst, dim=1), ubound(w_lst, dim=1)
                     call h5dopen_f(cg_g_id, all_cg%first%cg%w(w_lst(i))%name, dset_id, error)
                     if (cg_src_p(ncg) == proc) then
                        cg => get_nth_cg(cg_src_n(ncg))
                        pa4d => cg%w(w_lst(i))%arr(:, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) !< \todo use set_dims_for_restart
                        dims(:) = [ size(pa4d, dim=1, kind=4), cg%n_b ]
                     else
                        allocate(pa4d(size(all_cg%first%cg%w(w_lst(i))%arr(:,:,:,:), dim=1), cg_all_n_b(ncg, xdim), cg_all_n_b(ncg, ydim), cg_all_n_b(ncg, zdim)))
                        call MPI_Recv(pa4d(:,:,:,:), size(pa4d(:,:,:,:)), MPI_DOUBLE_PRECISION, cg_src_p(ncg), ncg + sum(cg_n(:))*(size(q_lst)+i), comm, MPI_STATUS_IGNORE, error)
                        dims(:) = shape(pa4d)
                     endif
                     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, pa4d, dims, error, xfer_prp = plist_id)
                     call h5dclose_f(dset_id, error)
                     if (cg_src_p(ncg) /= proc) deallocate(pa4d)
                  enddo
                  deallocate(dims)
               endif
            else
               if (can_i_write) call die("[restart_hdf5:write_cg_to_restart] Slave can write")
               if (cg_src_p(ncg) == proc) then
                  cg => get_nth_cg(cg_src_n(ncg))
                  if (size(q_lst) > 0) then
                     do i = lbound(q_lst, dim=1), ubound(q_lst, dim=1)
                        allocate(pa3d(cg%n_b(xdim), cg%n_b(ydim), cg%n_b(zdim)))
                        pa3d(:,:,:) = cg%q(q_lst(i))%arr(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke)
                        call MPI_Send(pa3d(:,:,:), size(pa3d(:,:,:)), MPI_DOUBLE_PRECISION, FIRST, ncg + sum(cg_n(:))*i, comm, error)
                        deallocate(pa3d)
                     enddo
                  endif
                  if (size(w_lst) > 0) then
                     do i = lbound(w_lst, dim=1), ubound(w_lst, dim=1)
                        allocate(pa4d(size(cg%w(w_lst(i))%arr(:,:,:,:), dim=1), cg%n_b(xdim), cg%n_b(ydim), cg%n_b(zdim)))
                        pa4d(:,:,:,:) = cg%w(w_lst(i))%arr(:, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke)
                        call MPI_Send(pa4d(:,:,:,:), size(pa4d(:,:,:,:)), MPI_DOUBLE_PRECISION, FIRST, ncg + sum(cg_n(:))*(size(q_lst)+i), comm, error)
                        deallocate(pa4d)
                     enddo
                  endif
               endif
            endif
         else ! perform parallell write
            ! This piece will be a generalization of the serial case. It should work correctly also for nproc_io == 1 so it should replace the serial code
            if (can_i_write) then
               ! write own
               ! receive (from whom?)
            else
               ! send (where?)
            endif
            call die("[restart_hdf5:write_cg_to_restart] Parallel v2 I/O not implemented yet")
         endif

         if (can_i_write) then
            call h5pclose_f(plist_id, error)
            call h5gclose_f(cg_g_id, error)
         endif

      enddo

      ! clean up
      deallocate(cg_src_p, cg_src_n)
      if (allocated(q_lst)) deallocate(q_lst)
      if (allocated(w_lst)) deallocate(w_lst)

   end subroutine write_cg_to_restart

!> \brief Expand given integer array by one and store the value i ni the last cell

   subroutine append_int_to_array(arr, i)

      implicit none

      integer, dimension(:), allocatable, intent(inout) :: arr
      integer, intent(in) :: i

      integer, allocatable, dimension(:) :: tmp

      if (allocated(arr)) then
         allocate(tmp(lbound(arr(:), dim=1):ubound(arr(:), dim=1) + 1))
         tmp(:ubound(arr(:), dim=1)) = arr(:)
         call move_alloc(from=tmp, to=arr)
      else
         allocate(arr(1))
      endif
      arr(ubound(arr(:))) = i

   end subroutine append_int_to_array

!> \brief find a n-th grid container on the cg_all list

   function get_nth_cg(n) result(cg)

      use dataio_pub, only: die
      use gc_list,    only: cg_list_element
      use grid,       only: all_cg
      use grid_cont,  only: grid_container

      implicit none

      integer, intent(in) :: n

      type(grid_container), pointer :: cg
      type(cg_list_element), pointer :: cgl

      integer :: i

      nullify(cg)
      i = 1
      cgl => all_cg%first
      do while (associated(cgl))
         if (i == n) then
            cg => cgl%cg
            exit
         endif
         i = i + 1
         cgl => cgl%nxt
      enddo

      if (.not. associated(cg)) call die("[restart_hdf5:get_nth_cg] cannot find n-th cg")

   end function get_nth_cg

!>
!! \brief Read a multi-file, multi-domain restart file
!!
!! \warning Not implemented yet
!!
!! \details Reading a v2 restart file should proceed as follows
!! - Identify the restart point file(s) on master
!! - Send the identification to the slaves
!! - On each processor check what files are present if any
!! - Check consistency of the data that should be present in all files in identical copies
!! - Prepare list of available grid containers
!! - If the list does not fill completely the boxes listed in refinement topology description on any process, communicate lists to the master
!! - If the compiled list does not fill completely the boxes listed in refinement topology description then call die()
!! - Compute domain decomposition for current refinement topology description (different number of processors and different decomposition method is allowed)
!! - Calculate on master which parts (if any) cannot be read directly and should be communicated, send the lists to the involved slaves
!! - Read the data (own and requested by others) from the file(s)
!! - Receive anything that had to be read on onther processes
!! - Check for local completness
!!
!! Optionally:
!! - When a single file is visible by multiple processes and the fileststem does not tolerate massively concurrent I/O load,
!!   choose one or few process to do the I/O and make the others to wait for the data.
!! - It should be possible to read runtime parameters from data stored in restart point rather than from actual problem.par file(s).
!!   A way to override some of them (in most cases  END_CONTROL::tend and nend) should be provided.
!<

   subroutine read_restart_hdf5_v2(status_v2)

      use constants,  only: cwdlen, cbuff_len, INVALID, RD
      use dataio_pub, only: die, warn, printio, msg, piernik_hdf5_version2
      use hdf5,       only: HID_T, H5F_ACC_RDONLY_F, h5open_f, h5close_f, h5fopen_f, h5fclose_f
      use h5lt,       only: h5ltget_attribute_double_f, h5ltget_attribute_int_f
      use mpisetup,   only: master

      implicit none

      integer, intent(out)  :: status_v2

      integer(HID_T) :: file_id       !> File identifier
      character(len=cwdlen) :: filename
      logical :: file_exist
      integer(kind=4) :: error
      real, dimension(:), allocatable :: rbuf
      integer(kind=4), dimension(:), allocatable :: ibuf
      character(len=cbuff_len), dimension(*), parameter :: real_attrs = [ "time         ", "timestep     ", "last_hdf_time",  &
           &                                                              "magic_mass   ", "next_t_tsl   ", "next_t_log   " ]
      character(len=cbuff_len), dimension(*), parameter :: int_attrs = [ "nstep            ", "nres             ", "nhdf             ", &
           &                                                             "step_res         ", "step_hdf         ", "require_init_prob" ]
      !> \deprecated same strings are used independently in set_common_attributes*
      integer :: ia

      call warn("[restart_hdf5:read_restart_hdf5_v2] Not implemented yet")

      filename = restart_fname(RD)
      if (master) then
         write(msg, '(2a)') 'Reading restart file: ', trim(filename)
         call printio(msg)
      endif

      inquire(file = filename, exist = file_exist)
      if (.not. file_exist) then
         write(msg,'(3a)') '[restart_hdf5:read_restart_hdf5_v2]: Restart file: ', trim(filename),' does not exist'
         call die(msg)
      endif

      status_v2 = INVALID

      call h5open_f(error)

      call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, error)

      ! Check file format version and compare attributes in the root of the restart point file"
      allocate(rbuf(1))

      call h5ltget_attribute_double_f(file_id, "/", "piernik", rbuf, error) !> \deprecated: magic string across multiple files
      if (error /= 0) call die("[restart_hdf5:read_restart_hdf5_v2] Cannot read 'piernik' attribute from the restart file. The file may be either damaged or incompatible")
      if (rbuf(1) > piernik_hdf5_version2) then
         write(msg,'(2(a,f5.2))')"[restart_hdf5:read_restart_hdf5_v2] Cannot read future versions of the restart file: ", rbuf(1)," > ", piernik_hdf5_version2
         call die(msg)
      else if (int(rbuf(1)) < int(piernik_hdf5_version2)) then
         write(msg,'(2(a,f5.2))')"[restart_hdf5:read_restart_hdf5_v2] The restart file is too ancient. It is unlikely that it could work correctly: ", rbuf(1)," << ", piernik_hdf5_version2
         call warn(msg)
         call h5fclose_f(file_id, error)
         call h5close_f(error)
         return
      else if (rbuf(1) < piernik_hdf5_version2) then
         write(msg,'(2(a,f5.2))')"[restart_hdf5:read_restart_hdf5_v2] Old versions of the restart file may not always work fully correctly: ", rbuf(1)," < ", piernik_hdf5_version2
         call warn(msg)
      endif
      call compare_real_array1D(rbuf(:))

      do ia = lbound(real_attrs, dim=1), ubound(real_attrs, dim=1)
         call h5ltget_attribute_double_f(file_id, "/", trim(real_attrs(ia)), rbuf, error)
         call compare_array1D(rbuf(:))
      enddo

      deallocate(rbuf)

      allocate(ibuf(1))
      do ia = lbound(int_attrs, dim=1), ubound(int_attrs, dim=1)
         call h5ltget_attribute_int_f(file_id, "/", trim(int_attrs(ia)), ibuf, error)
         call compare_array1D(ibuf(:))
      enddo
      deallocate(ibuf)

      call h5fclose_f(file_id, error)
      call h5close_f(error)

   end subroutine read_restart_hdf5_v2

!>
!! \brief Compare a 1D real array across all processed. Die if any difference is found
!!
!! \todo Add an optional parameter (string), which would identify, which comparison failed
!<

   subroutine compare_real_array1D(arr)

      use constants,  only: I_ONE
      use dataio_pub, only: die
      use mpi,        only: MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE
      use mpisetup,   only: slave, LAST, comm, proc

      implicit none

      real, dimension(:), intent(in) :: arr

      real, dimension(:), allocatable :: aux
      integer(kind=4), parameter :: tag = 10
      integer(kind=4) :: error

      allocate(aux(size(arr(:))))

      if (proc /= LAST) call MPI_Send(arr(:), size(arr(:), kind=4), MPI_DOUBLE_PRECISION, proc+I_ONE, tag, comm, error)
      if (slave) then
         call MPI_Recv(aux(:), size(aux(:), kind=4), MPI_DOUBLE_PRECISION, proc-I_ONE, tag, comm, MPI_STATUS_IGNORE, error)
         if (any(aux(:) /= arr(:))) call die("[restart_hdf5:compare_real_array1D] Inconsistency found.")
      endif

      deallocate(aux)

   end subroutine compare_real_array1D

!> \brief Compare a 1D integer array across all processed. Die if any difference is found

   subroutine compare_int_array1D(arr)

      use constants,  only: I_ONE
      use dataio_pub, only: die
      use mpi,        only: MPI_INTEGER, MPI_STATUS_IGNORE
      use mpisetup,   only: slave, LAST, comm, proc

      implicit none

      integer(kind=4), dimension(:), intent(in) :: arr

      integer(kind=4), dimension(:), allocatable :: aux
      integer(kind=4), parameter :: tag = 10
      integer(kind=4) :: error

      allocate(aux(size(arr(:))))

      if (proc /= LAST) call MPI_Send(arr(:), size(arr(:), kind=4), MPI_INTEGER, proc+I_ONE, tag, comm, error)
      if (slave) then
         call MPI_Recv(aux(:), size(aux(:), kind=4), MPI_INTEGER, proc-I_ONE, tag, comm, MPI_STATUS_IGNORE, error)
         if (any(aux(:) /= arr(:))) call die("[restart_hdf5:compare_int_array1D] Inconsistency found.")
      endif

      deallocate(aux)

   end subroutine compare_int_array1D

end module restart_hdf5
