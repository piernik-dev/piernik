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
!! \brief Module that contains HDF5 I/O routines for reading and writing old restart files. Provided for compatibility with old archives only.
!<
module restart_hdf5_v1
! pulled by HDF5

   implicit none

   private
   public :: read_restart_hdf5_v1, write_restart_hdf5_v1, read_arr_from_restart

contains

!>
!! \brief Routine to set dimensions of arrays related to grid containers in restart file
!! \param area_type case name; possibilities:
!!   AT_OUT_B - physical domain with outer boundaries,
!!   AT_NO_B  - only physical domain without any boundaries
!!   AT_USER  - user defined domain
!! \param chnk dimensions of data array dumped by this process
!! \param lleft left limits of data from array to be dumped
!! \param lright right limits of data from array to be dumped
!! \param loffs offset in area for this process
!<
   subroutine set_dims_for_restart(area_type, chnk, lleft, lright, loffs, cg)

      use constants,  only: ndims, AT_OUT_B, AT_NO_B, AT_USER, LO, HI
      use dataio_pub, only: die
      use domain,     only: dom
      use grid_cont,  only: grid_container
      use user_hooks, only: at_user_dims_settings

      implicit none

      integer(kind=4),                   intent(in)  :: area_type
      integer(kind=4), dimension(ndims), intent(out) :: lleft, lright
      integer,         dimension(ndims), intent(out) :: chnk
      integer(kind=8), dimension(ndims), intent(out) :: loffs
      type(grid_container), pointer,     intent(in)  :: cg !< current grid container

      ! only physical domain without any boundaries
      lleft(:)  = cg%ijkse(:, LO)
      lright(:) = cg%ijkse(:, HI)
      chnk(:)   = cg%n_b(:)
      loffs(:)  = cg%my_se(:, LO)

      select case (area_type)
         case (AT_OUT_B)                                   ! physical domain with outer boundaries
            where (cg%my_se(:, LO) == 0 .and. dom%has_dir(:))
               lleft(:)  = lleft(:)  - dom%nb
               chnk(:)   = chnk(:)   + dom%nb
            endwhere
            where (cg%h_cor1(:) == dom%n_d(:) .and. dom%has_dir(:)) !! \warning this should be checked against level%n_d
               lright(:) = lright(:) + dom%nb
               chnk(:)   = chnk(:)   + dom%nb
            endwhere
            where (loffs(:)>0) loffs(:) = loffs(:) + dom%nb ! the block adjacent to the left boundary are dom%nb cells wider than cg%n_b(:)
         case (AT_NO_B)                                    ! only physical domain without any boundaries
            ! Nothing special
         case (AT_USER)                                    ! user defined domain (with no reference to simulations domain)
            if (associated(at_user_dims_settings)) then
               call at_user_dims_settings(lleft, lright, chnk, loffs)
            else
               call die("[restart_hdf5_v1:set_dims_to_write] Routine at_user_dims_settings not associated")
            endif
         case default
            call die("[restart_hdf5_v1:set_dims_for_restart] Non-recognized area_type.")
      end select

   end subroutine set_dims_for_restart

!>
!! \brief Routine to set dimensions of arrays related to domain in restart file
!! \param area_type case name; possibilities:
!!   AT_OUT_B - physical domain with outer boundaries,
!!   AT_NO_B  - only physical domain without any boundaries
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
               call die("[restart_hdf5_v1:set_area_for_restart] Routine at_user_area_settings not associated")
            endif
         case default
            call die("[restart_hdf5_v1:set_area_for_restart] Non-recognized area_type.")
            area(:) = 0 ! suppress compiler warnings
      end select

   end subroutine set_area_for_restart

!> \brief Write restart dump v1.x

   subroutine write_restart_hdf5_v1(filename)

      use constants,        only: cwdlen
      use hdf5,             only: HID_T, H5P_FILE_ACCESS_F, H5F_ACC_RDWR_F, h5open_f, h5close_f, h5fopen_f, h5fclose_f, h5pcreate_f, h5pclose_f, h5pset_fapl_mpio_f
      !, H5P_DATASET_XFER_F, h5pset_preserve_f
      use mpi,              only: MPI_INFO_NULL
      use mpisetup,         only: comm
      use named_array_list, only: qna, wna

      implicit none
      character(len=cwdlen), intent(in) :: filename      !< HDF File name

      integer(kind=4)                   :: i
      integer(HID_T)                    :: file_id       !< File identifier
      integer(HID_T)                    :: plist_id      !< Property list identifier
      integer(kind=4)                   :: error

      ! Set up a new HDF5 file for parallel write
      call h5open_f(error)
      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
      call h5pset_fapl_mpio_f(plist_id, comm, MPI_INFO_NULL, error)
      call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error, access_prp = plist_id)

      ! Write scalar fields that were declared to be required on restart
      if (allocated(qna%lst)) then
         do i = lbound(qna%lst(:), dim=1, kind=4), ubound(qna%lst(:), dim=1, kind=4)
            call write_arr_to_restart(file_id, i, .true.)
         enddo
      endif

      if (allocated(wna%lst)) then
         do i = lbound(wna%lst(:), dim=1, kind=4), ubound(wna%lst(:), dim=1, kind=4)
            call write_arr_to_restart(file_id, i, .false.)
         enddo
      endif

!     \todo writing axes using collective I/O takes order of magnitude more than
!        dumping U and B arrays altogether, since XYZ-axis is not even read
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

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: ndims, AT_IGNORE, LONG, dsetnamelen
      use dataio_pub,       only: die
      use domain,           only: is_multicg
      use grid_cont,        only: grid_container
      use hdf5,             only: HID_T, HSIZE_T, H5T_NATIVE_DOUBLE, h5dwrite_f, h5sclose_f, h5pclose_f, h5dclose_f, &
           &                      H5P_DATASET_CREATE_F, H5S_SELECT_SET_F, H5P_DATASET_XFER_F, H5FD_MPIO_COLLECTIVE_F, &
           &                      h5screate_simple_f, h5pcreate_f, h5dcreate_f, h5dget_space_f, h5pset_dxpl_mpio_f, h5sselect_hyperslab_f
      use named_array_list, only: qna, wna
#ifdef INDEPENDENT_ATOUTB
      use constants,        only: AT_OUT_B
      use hdf5,             only: H5FD_MPIO_INDEPENDENT_F
#endif /* INDEPENDENT_ATOUTB */

      implicit none

      integer(HID_T),  intent(in)        :: file_id                   !< File identifier
      integer(kind=4), intent(in)        :: ind                       !< index of cg%q(:) or cg%w(:) array
      logical,         intent(in)        :: tgt3d                     !< .true. for 3D array, .false. otherwise

      real, pointer, dimension(:,:,:)    :: pa3d                      !< pointer to specified scope of pa3d
      real, pointer, dimension(:,:,:,:)  :: pa4d                      !< pointer to specified scope of pa4d
      integer, parameter                 :: rank4 = 1 + ndims
      integer(HSIZE_T), dimension(rank4) :: dimsf, chunk_dims
      integer(HID_T)                     :: dset_id                   !< Dataset identifier
      integer(HID_T)                     :: dplist_id, plist_id       !< Property list identifiers
      integer(HID_T)                     :: dfilespace, filespace     !< Dataspace identifiers in file
      integer(HID_T)                     :: memspace                  !< Dataspace identifier in memory

      integer, dimension(ndims)          :: area, chnk
      integer(kind=4), dimension(ndims)  :: lleft, lright
      integer(kind=8), dimension(ndims)  :: loffs
      integer(HSIZE_T), dimension(rank4) :: cnt, offset, stride
      integer(kind=4)                    :: rank, error, area_type
      integer                            :: ir, dim1
      type(cg_list_element), pointer     :: cgl
      type(grid_container),  pointer     :: cg
      character(len=dsetnamelen)         :: dname

      if (tgt3d) then
         if (ind < lbound(qna%lst(:), dim=1) .or. ind > ubound(qna%lst(:), dim=1)) call die("[restart_hdf5_v1:write_arr_to_restart] Invalid 3D array")
         dim1 = 1
         rank = ndims
         dname = qna%lst(ind)%name
         area_type = qna%lst(ind)%restart_mode
      else
         if (ind < lbound(wna%lst(:), dim=1) .or. ind > ubound(wna%lst(:), dim=1)) call die("[restart_hdf5_v1:write_arr_to_restart] Invalid 4D array")
         dim1 = wna%lst(ind)%dim4
         rank = rank4
         dname = wna%lst(ind)%name
         area_type = wna%lst(ind)%restart_mode
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
      !> \warning collective write may be improper (some fields may stay unwritten)
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

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         call set_dims_for_restart(area_type, chnk, lleft, lright, loffs, cg)

         chunk_dims = [dim1, chnk(:)] ! Chunks dimensions
         offset(:) = [0_LONG, loffs(:)]
         call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset(ir:), cnt(ir:), error, stride(ir:), chunk_dims(ir:))

         call h5screate_simple_f(rank, chunk_dims(ir:), memspace, error)

         ! write data
         if (tgt3d) then
            pa3d => cg%q(ind)%span(lleft, lright)
            if (associated(pa3d)) then
               call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, pa3d, dimsf(ir:), error, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
            else
               call die("[restart_hdf5_v1:write_arr_to_restart] unassociated 3D array pointer")
            endif
         else
            pa4d => cg%w(ind)%span(lleft, lright)
            if (associated(pa4d)) then
               call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, pa4d, dimsf(:), error, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
            else
               call die("[restart_hdf5_v1:write_arr_to_restart] unassociated 4D array pointer")
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
!! \details  Boundary cells are exchanged with the neighbours. Corner boundary cells are not guaranteed to be correct.
!<

   subroutine read_arr_from_restart(file_id, ind, tgt3d, alt_area_type, alt_name)

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: ndims, LONG, AT_IGNORE, dsetnamelen
      use dataio_pub,       only: msg, die
      use domain,           only: is_multicg
      use grid_cont,        only: grid_container
      use hdf5,             only: HID_T, HSIZE_T, H5T_NATIVE_DOUBLE, h5dread_f, h5sclose_f, h5pclose_f, h5dclose_f, &
           &                      H5S_SELECT_SET_F, H5P_DATASET_XFER_F, H5FD_MPIO_COLLECTIVE_F, &
           &                      h5dopen_f, h5sget_simple_extent_ndims_f, h5dget_space_f, &
           &                      h5pcreate_f, h5pset_dxpl_mpio_f, h5sselect_hyperslab_f, h5screate_simple_f
      use named_array_list, only: qna, wna

      implicit none

      integer(HID_T),             intent(in) :: file_id            !< File identifier
      integer,                    intent(in) :: ind                !< index of cg%q(:) or cg%w(:) arrays
      logical,                    intent(in) :: tgt3d              !< .true. for 3D arrays, .false. otherwise
      integer(kind=4),  optional, intent(in) :: alt_area_type
      character(len=*), optional, intent(in) :: alt_name           !< used only in galdisk* setups

      real, pointer, dimension(:,:,:)        :: pa3d               !< pointer to specified scope of pa3d
      real, pointer, dimension(:,:,:,:)      :: pa4d               !< pointer to specified scope of pa4d
      integer, parameter                     :: rank4 = 1 + ndims
      integer(HSIZE_T), dimension(rank4)     :: dimsf, chunk_dims
      integer(HID_T)                         :: dset_id            !< Dataset identifier
      integer(HID_T)                         :: plist_id           !< Property list identifier
      integer(HID_T)                         :: filespace          !< Dataspace identifier in file
      integer(HID_T)                         :: memspace           !< Dataspace identifier in memory
      integer(HSIZE_T), dimension(rank4)     :: cnt, offset, stride
      integer, dimension(ndims)              :: area, chnk
      integer(kind=4), dimension(ndims)      :: lleft, lright
      integer(kind=8),  dimension(ndims)     :: loffs
      integer(kind=4)                        :: rank, rankf, error, area_type
      integer                                :: ir, dim1
      type(cg_list_element), pointer         :: cgl
      type(grid_container),  pointer         :: cg
      character(len=dsetnamelen)             :: dname, cgname

      !> \deprecated Duplicated code
      if (tgt3d) then
         if (ind < lbound(qna%lst(:), dim=1) .or. ind > ubound(qna%lst(:), dim=1)) call die("[restart_hdf5_v1:read_arr_from_restart] Invalid 3D array")
         dim1 = 1
         rank = ndims
         cgname = qna%lst(ind)%name
         area_type = qna%lst(ind)%restart_mode
      else
         if (ind < lbound(wna%lst(:), dim=1) .or. ind > ubound(wna%lst(:), dim=1)) call die("[restart_hdf5_v1:read_arr_from_restart] Invalid 4D array")
         dim1 = wna%lst(ind)%dim4
         rank = rank4
         cgname = wna%lst(ind)%name
         area_type = wna%lst(ind)%restart_mode
      endif
      dname = cgname
      ir = rank4 - rank + 1 ! 1 for 4-D arrays, 2 for 3-D arrays (to simplify use of count(:), offset(:), stride(:), block(:), dimsf(:) and chunk_dims(:)

      if (present(alt_area_type)) area_type = alt_area_type
      if (area_type == AT_IGNORE) return !! \todo write a list of unsaved arrays?
      call set_area_for_restart(area_type, area)

      dimsf = [dim1, area(:)]      ! Dataset dimensions

      ! Create dataset.and filespace
      if (present(alt_name)) dname = alt_name
      call h5dopen_f(file_id, dname, dset_id, error)
      if (error /= 0) then
         write(msg, '(3a)') "[restart_hdf5_v1:read_arr_from_restart] Opening dataset '", dname,"' failed."
         call die(msg)
      endif

      call h5dget_space_f(dset_id, filespace, error)
      call h5sget_simple_extent_ndims_f (filespace, rankf, error)
      if (rank /= rankf) then
         write(msg,'(3a,2(i2,a))')"[restart_hdf5_v1:read_arr_from_restart] Rank mismatch in array '", dname, "' (", rank, " /= ", rankf, ")"
         call die(msg)
      endif

      ! Select hyperslab in the file.
      stride(:) = 1
      cnt(:)  = 1

      cgl => leaves%first
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
            pa3d => cg%q(qna%ind(cgname))%span(lleft, lright)
            if (associated(pa3d)) then
               call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, pa3d, dimsf(ir:), error, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
            else
               write(msg,'(5a)')"Unassociated '",trim(cgname),"' 3D array while reading dataset '",trim(dname),"' from the restart file"
               call die(msg)
            endif
         else
            pa4d => cg%w(wna%ind(cgname))%span(lleft, lright)
            if (associated(pa4d)) then
               call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, pa4d, dimsf(ir:), error, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
            else
               write(msg,'(5a)')"Unassociated '",trim(cgname),"' 4D array while reading dataset '",trim(dname),"' from the restart file"
               call die(msg)
            endif
         endif

         if (error /= 0) then
            write(msg, '(3a)') "[restart_hdf5_v1:read_arr_from_restart] Reading dataset '", dname,"' failed."
            call die(msg)
         endif

         call h5sclose_f(memspace, error)

         cgl => cgl%nxt
      enddo

      call h5pclose_f(plist_id, error)
      call h5sclose_f(filespace, error)
      call h5dclose_f(dset_id, error)

      if (tgt3d) call leaves%leaf_arr3d_boundaries(qna%ind(dname), area_type=area_type)
      ! Originally the pa3d array was written with the guardcells. The internal guardcells will be exchanged but the external ones are lost.

      ! rank-4 arrays (cg%u(:,:,:,:) and b(:,:,:,:)) have their own guardcell-exchange routines, which can also be called here
      !> \todo consider also calling 4D boundaries

   end subroutine read_arr_from_restart

!> \brief Read a v1 restart point.

   subroutine read_restart_hdf5_v1

      use common_hdf5,      only: output_fname
      use constants,        only: cwdlen, cbuff_len, domlen, idlen, xdim, ydim, zdim, LO, HI, RD
      use dataio_pub,       only: msg, warn, die, printio, require_problem_IC, problem_name, piernik_hdf5_version, fix_string, &
           &                      domain_dump, last_hdf_time, last_res_time, last_log_time, last_tsl_time, nhdf, nres, new_id
      use dataio_user,      only: user_reg_var_restart, user_attrs_rd
      use domain,           only: dom
      use fluidindex,       only: flind
      use global,           only: t, dt, nstep
      use hdf5,             only: HID_T, H5P_FILE_ACCESS_F, H5F_ACC_RDONLY_F, &
           &                      h5open_f, h5pcreate_f, h5pset_fapl_mpio_f, h5fopen_f, h5pclose_f, h5fclose_f, h5close_f
      use h5lt,             only: h5ltget_attribute_double_f, h5ltget_attribute_int_f, h5ltget_attribute_string_f
      use mass_defect,      only: magic_mass
      use mpi,              only: MPI_INFO_NULL
      use mpisetup,         only: comm, master, piernik_MPI_Bcast, ibuff, rbuff, cbuff, slave
      use named_array_list, only: qna, wna

      implicit none

      integer                                  :: nu
      integer                                  :: i
      character(len=cwdlen)                    :: filename      !< File name

      integer(HID_T)                           :: file_id       !< File identifier
      integer(HID_T)                           :: plist_id      !< Property list identifier

      integer(kind=4)                          :: error
      logical                                  :: file_exist

      real,            dimension(1)            :: rbuf
      real,            dimension(flind%fluids) :: rbufm
      integer(kind=4), dimension(1)            :: ibuf

      real                                     :: restart_hdf5_version

      nu = flind%all

      filename = output_fname(RD, '.res', nres, bcast=.true.)

      if (master) then
         write(msg, '(2a)') 'Reading restart file v1: ', trim(filename)
         call printio(msg)
      endif

      inquire(file = trim(filename), exist = file_exist)
      if (.not. file_exist) then
         write(msg,'(3a)') '[restart_hdf5_v1:read_restart_hdf5_v1]: Restart file: ', trim(filename),' does not exist'
         call die(msg)
      endif

      call h5open_f(error)
      if (master) then
         call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, error)

         call h5ltget_attribute_double_f(file_id,"/","piernik", rbuf, error)
         if (error /= 0) call die("[restart_hdf5_v1:read_restart_hdf5_v1] Cannot read 'piernik' attribute from the restart file. The file may be either damaged or incompatible")
         if (rbuf(1) > piernik_hdf5_version) then
            write(msg,'(2(a,f5.2))')"[restart_hdf5_v1:read_restart_hdf5_v1] Cannot read future versions of the restart file: ", rbuf(1)," > ", piernik_hdf5_version
            call die(msg)
         else if (int(rbuf(1)) < int(piernik_hdf5_version)) then
            write(msg,'(2(a,f5.2))')"[restart_hdf5_v1:read_restart_hdf5_v1] The restart file is too ancient. It is unlikely that it could work correctly: ", rbuf(1)," << ", piernik_hdf5_version
            call die(msg)
         else if (rbuf(1) < piernik_hdf5_version) then
            write(msg,'(2(a,f5.2))')"[restart_hdf5_v1:read_restart_hdf5_v1] Old versions of the restart file may not always work fully correctly: ", rbuf(1)," < ", piernik_hdf5_version
            call warn(msg)
         endif

         call h5ltget_attribute_int_f(file_id,"/","nxd", ibuf, error)
         if (ibuf(1) /= dom%n_d(xdim) .or. error /= 0) call die("[restart_hdf5_v1:read_restart_hdf5_v1] nxd does not match")
         if (dom%has_dir(xdim)) then
            call h5ltget_attribute_double_f(file_id,"/","xmin", rbuf, error)
            if (rbuf(1) /= dom%edge(xdim, LO) .or. error /= 0) call die("[restart_hdf5_v1:read_restart_hdf5_v1] xmin does not match")
            call h5ltget_attribute_double_f(file_id,"/","xmax", rbuf, error)
            if (rbuf(1) /= dom%edge(xdim, HI) .or. error /= 0) call die("[restart_hdf5_v1:read_restart_hdf5_v1] xmax does not match")
         endif

         call h5ltget_attribute_int_f(file_id,"/","nyd", ibuf, error)
         if (ibuf(1) /= dom%n_d(ydim) .or. error /= 0) call die("[restart_hdf5_v1:read_restart_hdf5_v1] nyd does not match")
         if (dom%has_dir(ydim)) then
            call h5ltget_attribute_double_f(file_id,"/","ymin", rbuf, error)
            if (rbuf(1) /= dom%edge(ydim, LO) .or. error /= 0) call die("[restart_hdf5_v1:read_restart_hdf5_v1] ymin does not match")
            call h5ltget_attribute_double_f(file_id,"/","ymax", rbuf, error)
            if (rbuf(1) /= dom%edge(ydim, HI) .or. error /= 0) call die("[restart_hdf5_v1:read_restart_hdf5_v1] ymax does not match")
         endif

         call h5ltget_attribute_int_f(file_id,"/","nzd", ibuf, error)
         if (ibuf(1) /= dom%n_d(zdim) .or. error /= 0) call die("[restart_hdf5_v1:read_restart_hdf5_v1] nzd does not match")
         if (dom%has_dir(zdim)) then
            call h5ltget_attribute_double_f(file_id,"/","zmin", rbuf, error)
            if (rbuf(1) /= dom%edge(zdim, LO) .or. error /= 0) call die("[restart_hdf5_v1:read_restart_hdf5_v1] zmin does not match")
            call h5ltget_attribute_double_f(file_id,"/","zmax", rbuf, error)
            if (rbuf(1) /= dom%edge(zdim, HI) .or. error /= 0) call die("[restart_hdf5_v1:read_restart_hdf5_v1] zmax does not match")
         endif

         call h5fclose_f(file_id, error)
      endif

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
      call h5pset_fapl_mpio_f(plist_id, comm, MPI_INFO_NULL, error)

      call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, error, access_prp = plist_id)
      call h5pclose_f(plist_id, error)

      ! set up things such as register user rank-3 and rank-4 arrays to be read by read_arr_from_restart. Read also anything that is not read by all read_arr_from_restart calls
      if (associated(user_attrs_rd)) call user_attrs_rd(file_id)
      if (associated(user_reg_var_restart)) call user_reg_var_restart

      ! read auxiliary variables
      if (allocated(qna%lst)) then
         do i = lbound(qna%lst(:), dim=1, kind=4), ubound(qna%lst(:), dim=1, kind=4)
            call read_arr_from_restart(file_id, i, .true.)
         enddo
      endif

      if (allocated(wna%lst)) then
         do i = lbound(wna%lst(:), dim=1, kind=4), ubound(wna%lst(:), dim=1, kind=4)
            call read_arr_from_restart(file_id, i, .false.)
         enddo
      endif

      call h5fclose_f(file_id, error)

      if (master) then
         call h5fopen_f (filename, H5F_ACC_RDONLY_F, file_id, error)
         call h5ltget_attribute_double_f(file_id,"/","piernik",       rbuf, error) ; restart_hdf5_version = rbuf(1)
         call h5ltget_attribute_double_f(file_id,"/","time",          rbuf, error) ; t  = rbuf(1)
         call h5ltget_attribute_double_f(file_id,"/","timestep",      rbuf, error) ; dt = rbuf(1)
         call h5ltget_attribute_double_f(file_id,"/","magic_mass",   rbufm, error) ; magic_mass(:) = rbufm(:)
         call h5ltget_attribute_double_f(file_id,"/","last_log_time", rbuf, error) ; last_log_time = rbuf(1)
         call h5ltget_attribute_double_f(file_id,"/","last_tsl_time", rbuf, error) ; last_tsl_time = rbuf(1)
         call h5ltget_attribute_double_f(file_id,"/","last_hdf_time", rbuf, error) ; last_hdf_time = rbuf(1)
         call h5ltget_attribute_double_f(file_id,"/","last_res_time", rbuf, error) ; last_res_time = rbuf(1)
         call h5ltget_attribute_int_f(file_id,"/","nstep", ibuf, error) ; nstep = ibuf(1)
         call h5ltget_attribute_int_f(file_id,"/","nres",  ibuf, error) ; nres  = ibuf(1)
         call h5ltget_attribute_int_f(file_id,"/","nhdf",  ibuf, error) ; nhdf  = ibuf(1)

         call h5ltget_attribute_string_f(file_id,"/","problem_name", problem_name, error)
         call h5ltget_attribute_string_f(file_id,"/","domain",       domain_dump,  error)
         call h5ltget_attribute_string_f(file_id,"/","run_id",       new_id,       error)

         if (restart_hdf5_version > 1.11) then
            call h5ltget_attribute_int_f(file_id,"/","require_problem_IC", ibuf, error)
            require_problem_IC = ibuf(1)
         endif

         problem_name = fix_string(problem_name)   !> \deprecated BEWARE: >=HDF5-1.8.4 has weird issues with strings
         new_id  = fix_string(new_id)    !> \deprecated   this bit hacks it around
         domain_dump  = fix_string(domain_dump)

         call h5fclose_f(file_id, error)

         write(msg,'(2a)') 'Done reading restart file v1: ', trim(filename)
         call printio(msg)
      endif
      call h5close_f(error)

      call piernik_MPI_Bcast(restart_hdf5_version)

      if (master) then
         ibuff(1) = nstep
         ibuff(2) = nres
         ibuff(3) = nhdf
         if (restart_hdf5_version > 1.11) ibuff(5) = require_problem_IC

         rbuff(1) = last_log_time
         rbuff(2) = last_tsl_time
         rbuff(3) = last_hdf_time
         rbuff(4) = last_res_time
         rbuff(6) = t
         rbuff(7) = dt

         cbuff(1) = problem_name
         cbuff(2)(1:domlen) = domain_dump
         cbuff(3)(1:idlen)  = new_id
      endif

      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(rbuff)
      call piernik_MPI_Bcast(cbuff, cbuff_len)

      if (slave) then
         nstep = ibuff(1)
         nres = ibuff(2)
         nhdf = ibuff(3)
         if (restart_hdf5_version > 1.11) require_problem_IC = ibuff(5)

         last_log_time = rbuff(1)
         last_tsl_time = rbuff(2)
         last_hdf_time = rbuff(3)
         last_res_time = rbuff(4)
         t = rbuff(6)
         dt = rbuff(7)

         problem_name = cbuff(1)
         domain_dump = cbuff(2)(1:domlen)
         new_id = cbuff(3)(1:idlen)

      endif

   end subroutine read_restart_hdf5_v1

end module restart_hdf5_v1
