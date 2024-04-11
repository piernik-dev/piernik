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

!> \brief Module that contains HDF5 I/O routines for reading and writing restart files v2 (AMR).

module restart_hdf5_v2
! pulled by HDF5

   use constants, only: ndims

   implicit none

   private
   public :: read_restart_hdf5_v2, write_restart_hdf5_v2, cg_essentials

   type :: cg_essentials                            !< All vital attributes of a cg in one place
      integer(kind=4), dimension(ndims) :: n_b
      integer(kind=8), dimension(ndims) :: off
      integer(kind=4)                   :: level
   end type cg_essentials
   ! nearly duplicate with level_essentials::level_t

contains

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
!! At the moment we assume that all problem-specific data is either registered in cg%q and cg%w arrays or can be written to global attributes.
!!
!! It will be possible to provide support for problem-specific data that should be appended to grid containers when anyone will need it
!!
!! \todo Check if it is possible to filter the data through shuffle and gzip -9  during write
!<

!>
!! \brief Write a multi-file, multi-domain HDF5 file
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

      use common_hdf5, only: write_to_hdf5_v2, O_RES
      use constants,   only: cwdlen

      implicit none

      character(len=cwdlen), intent(in) :: filename

      call write_to_hdf5_v2(filename, O_RES, create_empty_cg_datasets_in_restart, write_cg_to_restart)

   end subroutine write_restart_hdf5_v2

!> \brief create empty datasets for each cg to store restart data
   subroutine create_empty_cg_datasets_in_restart(cg_g_id, cg_n_b, cg_n_o, Z_avail, n_part, st_g_id)

      use common_hdf5,      only: create_empty_cg_dataset, O_RES
      use constants,        only: AT_IGNORE, AT_OUT_B, I_ONE, PPP_IO, PPP_CG
      use dataio_pub,       only: die
      use hdf5,             only: HID_T, HSIZE_T
      use named_array_list, only: qna, wna
      use ppp,              only: ppp_main
#ifdef NBODY
      use data_hdf5,        only: gdf_translate
      use particles_io,     only: pvarn
#endif /* NBODY */

      implicit none

      integer(HID_T),                intent(in) :: cg_g_id
      integer(kind=4), dimension(:), intent(in) :: cg_n_b
      integer(kind=4), dimension(:), intent(in) :: cg_n_o
      logical(kind=4),               intent(in) :: Z_avail
      integer(kind=8),               intent(in) :: n_part
      integer(HID_T),                intent(in) :: st_g_id

      integer :: i
      integer(HSIZE_T), dimension(:), allocatable :: d_size
      character(len=*), parameter :: wrce_label = "IO_write_empty_dset"

      if (size(cg_n_b) /= size(cg_n_o)) call die("[restart_hdf5_v2:create_empty_cg_datasets_in_restart] size(cg_n_b) /= size(cg_n_o)")
      call ppp_main%start(wrce_label, PPP_IO + PPP_CG)

      allocate(d_size(size(cg_n_b)))
      if (allocated(qna%lst)) then
         do i = lbound(qna%lst(:), dim=1), ubound(qna%lst(:), dim=1)
            if (qna%lst(i)%restart_mode == AT_OUT_B) then
               d_size = int(cg_n_o, kind=HSIZE_T)
            else
               d_size = int(cg_n_b, kind=HSIZE_T)
            endif
            if (qna%lst(i)%restart_mode > AT_IGNORE) &  ! create "/data/grid_%08d/qna%lst(i)%name"
                 call create_empty_cg_dataset(cg_g_id, qna%lst(i)%name, d_size, Z_avail, O_RES)
         enddo
      endif
      deallocate(d_size)

      allocate(d_size(size(cg_n_b)+I_ONE))
      if (allocated(wna%lst)) then
         do i = lbound(wna%lst(:), dim=1), ubound(wna%lst(:), dim=1)
            if (wna%lst(i)%restart_mode == AT_OUT_B) then
               d_size = int([ wna%lst(i)%dim4, cg_n_o ], kind=HSIZE_T)
            else
               d_size = int([ wna%lst(i)%dim4, cg_n_b ], kind=HSIZE_T)
            endif
            if (wna%lst(i)%restart_mode > AT_IGNORE) &  ! create "/data/grid_%08d/wna%lst(i)%name"
                 call create_empty_cg_dataset(cg_g_id, wna%lst(i)%name, d_size, Z_avail, O_RES)
         enddo
      endif
      deallocate(d_size)

#ifdef NBODY
      do i = lbound(pvarn, 1, kind=4), ubound(pvarn, 1, kind=4)
         call create_empty_cg_dataset(st_g_id, gdf_translate(pvarn(i)), [n_part], Z_avail, O_RES)
      enddo
#endif /* NBODY */

      call ppp_main%stop(wrce_label, PPP_IO + PPP_CG)

      return
      if (.false.) i = int(n_part) + int(st_g_id)

   end subroutine create_empty_cg_datasets_in_restart

!> \brief Write all grid containers to the file

   subroutine write_cg_to_restart(cgl_g_id, cg_n, cg_all_n_b, cg_all_n_o)

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use common_hdf5,      only: get_nth_cg, cg_output
      use constants,        only: xdim, ydim, zdim, ndims, dsetnamelen, I_ONE, PPP_IO, PPP_CG
      use dataio_pub,       only: die, nproc_io, can_i_write
      use grid_cont,        only: grid_container
      use hdf5,             only: HID_T, HSIZE_T, H5T_NATIVE_DOUBLE, h5sclose_f, h5dwrite_f, h5sselect_none_f, h5screate_simple_f
      use MPIF,             only: MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, MPI_COMM_WORLD
      use MPIFUN,           only: MPI_Recv, MPI_Send
      use mpisetup,         only: master, FIRST, proc, err_mpi, tag_ub
      use named_array_list, only: qna, wna
      use ppp,              only: ppp_main
#ifdef NBODY
      use data_hdf5,        only: gdf_translate
      use particles_io,     only: parallel_nbody_datafields, serial_nbody_datafields, pvarn
#endif /* NBODY */

      implicit none

      integer(HID_T),                           intent(in)  :: cgl_g_id    !< cg group identifier
      integer(kind=4), dimension(:),   pointer, intent(in)  :: cg_n        !< offset for cg group numbering
      integer(kind=4), dimension(:,:), pointer, intent(in)  :: cg_all_n_b  !< all cg sizes
      integer(kind=4), dimension(:,:), pointer, intent(in)  :: cg_all_n_o  !< all cg sizes, expanded by external boundaries

      integer(HID_T)                                        :: filespace_id, memspace_id
      integer(kind=4)                                       :: error
      integer(kind=4), parameter                            :: rank4 = I_ONE + ndims
      integer(kind=4), parameter                            :: rank3 = ndims
      integer(HSIZE_T), dimension(:), allocatable           :: dims
      real, pointer,    dimension(:,:,:)                    :: pa3d
      real, pointer,    dimension(:,:,:,:)                  :: pa4d
      integer                                               :: tot_lst_n, ic
      integer(kind=4)                                       :: ncg, n, i, ntags
      type(grid_container),  pointer                        :: cg
      type(cg_list_element), pointer                        :: cgl
      integer, allocatable, dimension(:)                    :: qr_lst, wr_lst
      type(cg_output)                                       :: cg_desc
      character(len=dsetnamelen), dimension(:), allocatable :: dsets
      real, target, dimension(0,0,0)                        :: null_r3d
      real, target, dimension(0,0,0,0)                      :: null_r4d
      integer(kind=4), dimension(ndims)                     :: n_b
      character(len=*), parameter :: wrcg_label = "IO_write_restart_v2_all_cg", wrcg1s_label = "IO_write_restart_v2_1cg_(serial)", wrcg1p_label = "IO_write_restart_v2_1cg_(parallel)"

      call ppp_main%start(wrcg_label, PPP_IO)

      call qna%get_reslst(qr_lst)
      call wna%get_reslst(wr_lst)
      tot_lst_n = size(qr_lst) + size(wr_lst)
      allocate(dsets(tot_lst_n))
      ic = 1
      if (size(qr_lst) > 0) then
         do i = lbound(qr_lst, dim=1, kind=4), ubound(qr_lst, dim=1, kind=4)
            dsets(ic) = trim(qna%lst(qr_lst(i))%name)
            ic = ic + 1
         enddo
      endif
      if (size(wr_lst) > 0) then
         do i = lbound(wr_lst, dim=1, kind=4), ubound(wr_lst, dim=1, kind=4)
            dsets(ic) = trim(wna%lst(wr_lst(i))%name)
            ic = ic + 1
         enddo
      endif
#ifdef NBODY
      call cg_desc%init(cgl_g_id, cg_n, nproc_io, ntags, dsets, gdf_translate(pvarn))
#else /* !NBODY */
      call cg_desc%init(cgl_g_id, cg_n, nproc_io, ntags, dsets)
#endif /* !NBODY */

      if (nproc_io == 1) then ! perform serial write
         ! write all cg, one by one
         if (cg_desc%tot_cg_n * ntags > tag_ub) call die("[restart_hdf5_v2:write_cg_to_restart] this MPI implementation has too low MPI_TAG_UB attribute")
         do ncg = 1, cg_desc%tot_cg_n
            call ppp_main%start(wrcg1s_label, PPP_IO + PPP_CG)

            if (master) then
               if (.not. can_i_write) call die("[restart_hdf5_v2:write_cg_to_restart] Master can't write")

               if (size(qr_lst) > 0) then
                  allocate(dims(rank3))
                  do i = lbound(qr_lst, dim=1, kind=4), ubound(qr_lst, dim=1, kind=4)
                     if (cg_desc%cg_src_p(ncg) == proc) then
                        cg => get_nth_cg(cg_desc%cg_src_n(ncg))
                        pa3d => cg%q(qr_lst(i))%span(pick_area(cg, qna%lst(qr_lst(i))%restart_mode))
                        dims(:) = cg%n_b
                     else
                        n_b = pick_size(ncg, qna%lst(qr_lst(i))%restart_mode)
                        allocate(pa3d(n_b(xdim), n_b(ydim), n_b(zdim)))
                        call MPI_Recv(pa3d(:,:,:), size(pa3d(:,:,:), kind=4), MPI_DOUBLE_PRECISION, cg_desc%cg_src_p(ncg), ncg + cg_desc%tot_cg_n*i, MPI_COMM_WORLD, MPI_STATUS_IGNORE, err_mpi)
                        dims(:) = shape(pa3d)
                     endif
                     call h5dwrite_f(cg_desc%dset_id(ncg, i), H5T_NATIVE_DOUBLE, pa3d, dims, error, xfer_prp = cg_desc%xfer_prp)
                     if (cg_desc%cg_src_p(ncg) /= proc) deallocate(pa3d)
                  enddo
                  deallocate(dims)
               endif

               if (size(wr_lst) > 0) then
                  allocate(dims(rank4))
                  do i = lbound(wr_lst, dim=1, kind=4), ubound(wr_lst, dim=1, kind=4)
                     if (cg_desc%cg_src_p(ncg) == proc) then
                        cg => get_nth_cg(cg_desc%cg_src_n(ncg))
                        pa4d => cg%w(wr_lst(i))%span(pick_area(cg, wna%lst(wr_lst(i))%restart_mode))
                        dims(:) = [ wna%lst(wr_lst(i))%dim4, cg%n_b ]
                     else
                        n_b = pick_size(ncg, wna%lst(wr_lst(i))%restart_mode)
                        allocate(pa4d(wna%lst(wr_lst(i))%dim4, n_b(xdim), n_b(ydim), n_b(zdim)))
                        call MPI_Recv(pa4d(:,:,:,:), size(pa4d(:,:,:,:), kind=4), MPI_DOUBLE_PRECISION, cg_desc%cg_src_p(ncg), &
                           ncg + cg_desc%tot_cg_n*(size(qr_lst, kind=4)+i), MPI_COMM_WORLD, MPI_STATUS_IGNORE, err_mpi)
                        dims(:) = shape(pa4d)
                     endif
                     call h5dwrite_f(cg_desc%dset_id(ncg, i+size(qr_lst)), H5T_NATIVE_DOUBLE, pa4d, dims, error, xfer_prp = cg_desc%xfer_prp)
                     if (cg_desc%cg_src_p(ncg) /= proc) deallocate(pa4d)
                  enddo
                  deallocate(dims)
               endif
            else
               if (can_i_write) call die("[restart_hdf5_v2:write_cg_to_restart] Slave can write")
               if (cg_desc%cg_src_p(ncg) == proc) then
                  cg => get_nth_cg(cg_desc%cg_src_n(ncg))
                  if (size(qr_lst) > 0) then
                     do i = lbound(qr_lst, dim=1, kind=4), ubound(qr_lst, dim=1, kind=4)
                        pa3d => cg%q(qr_lst(i))%span(pick_area(cg, qna%lst(qr_lst(i))%restart_mode))
                        call MPI_Send(pa3d(:,:,:), size(pa3d(:,:,:), kind=4), MPI_DOUBLE_PRECISION, FIRST, ncg + cg_desc%tot_cg_n*i, MPI_COMM_WORLD, err_mpi)
                     enddo
                  endif
                  if (size(wr_lst) > 0) then
                     do i = lbound(wr_lst, dim=1, kind=4), ubound(wr_lst, dim=1, kind=4)
                        pa4d => cg%w(wr_lst(i))%span(pick_area(cg, wna%lst(wr_lst(i))%restart_mode))
                        call MPI_Send(pa4d(:,:,:,:), size(pa4d(:,:,:,:), kind=4), MPI_DOUBLE_PRECISION, FIRST, ncg + cg_desc%tot_cg_n*(size(qr_lst, kind=4)+i), MPI_COMM_WORLD, err_mpi)
                     enddo
                  endif
               endif
            endif
#ifdef NBODY
            call serial_nbody_datafields(cg_desc%pdset_id, gdf_translate(pvarn), ncg, cg_desc%cg_src_n(ncg), cg_desc%cg_src_p(ncg), cg_desc%tot_cg_n)
#endif /* NBODY */
            call ppp_main%stop(wrcg1s_label, PPP_IO + PPP_CG)
         enddo
      else ! perform parallel write
         ! This piece will be a generalization of the serial case. It should work correctly also for nproc_io == 1 so it should replace the serial code

         ! It seems that the fastest implementation would use _at most_ single thread per node.
         ! This way the MPI communication would be intra-node only and there would be not too many I/O-active threads.

         if (can_i_write) then
            ! write own
            n = 0
            cgl => leaves%first
            do while (associated(cgl))
               call ppp_main%start(wrcg1p_label, PPP_IO + PPP_CG)
               n = n + I_ONE
               ncg = cg_desc%offsets(proc) + n
               cg => cgl%cg
               ic = 0
               if (size(qr_lst) > 0) then
                  allocate(dims(rank3))
                  dims(:) = cg%n_b
                  do i = lbound(qr_lst, dim=1, kind=4), ubound(qr_lst, dim=1, kind=4)
                     ic = ic + 1
                     pa3d => cg%q(qr_lst(i))%span(pick_area(cg, qna%lst(qr_lst(i))%restart_mode))
                     dims(:) = pick_dims(cg, qna%lst(qr_lst(i))%restart_mode)
                     call h5dwrite_f(cg_desc%dset_id(ncg, ic), H5T_NATIVE_DOUBLE, pa3d, dims, error, xfer_prp = cg_desc%xfer_prp)
                  enddo
                  deallocate(dims)
               endif
               if (size(wr_lst) > 0) then
                  allocate(dims(rank4))
                  do i = lbound(wr_lst, dim=1, kind=4), ubound(wr_lst, dim=1, kind=4)
                     ic = ic + 1
                     pa4d => cg%w(wr_lst(i))%span(pick_area(cg, wna%lst(wr_lst(i))%restart_mode))
                     dims(:) = [ wna%lst(wr_lst(i))%dim4, pick_dims(cg, wna%lst(wr_lst(i))%restart_mode) ]
                     call h5dwrite_f(cg_desc%dset_id(ncg, ic), H5T_NATIVE_DOUBLE, pa4d, dims, error, xfer_prp = cg_desc%xfer_prp)
                  enddo
                  deallocate(dims)
               endif

#ifdef NBODY
               call parallel_nbody_datafields(cg_desc%pdset_id, gdf_translate(pvarn), ncg, cg)
#endif /* NBODY */

               cgl => cgl%nxt
               call ppp_main%stop(wrcg1p_label, PPP_IO + PPP_CG)
            enddo

            ! Parallel HDF calls have to be matched on each process, so we're doing some calls on non-existing data here.
            do ncg = I_ONE, maxval(cg_n)-n
               ic = 0
               if (size(qr_lst) > 0) then
                  allocate(dims(rank3))
                  dims(:) = 0
                  call h5screate_simple_f(rank3, dims, filespace_id, error)
                  call h5sselect_none_f(filespace_id, error)  ! empty filespace
                  call h5screate_simple_f(rank3, dims, memspace_id, error)
                  call h5sselect_none_f(memspace_id, error)   ! empty memoryscape

                  do i = lbound(qr_lst, dim=1, kind=4), ubound(qr_lst, dim=1, kind=4)
                     ic = ic + 1
                     pa3d => null_r3d
                     dims(:) = 0
                     call h5dwrite_f(cg_desc%dset_id(1, ic), H5T_NATIVE_DOUBLE, pa3d, dims, error, &
                        xfer_prp = cg_desc%xfer_prp, file_space_id = filespace_id, mem_space_id = memspace_id)
                  enddo
                  call h5sclose_f(memspace_id, error)
                  call h5sclose_f(filespace_id, error)
                  deallocate(dims)
               endif
               if (size(wr_lst) > 0) then
                  allocate(dims(rank4))
                  do i = lbound(wr_lst, dim=1, kind=4), ubound(wr_lst, dim=1, kind=4)
                     ic = ic + 1
                     dims(:) = 0
                     ! dims can change so h5s* is here as well, I don't know it
                     ! it's really necessary though,
                     ! \todo check if it can be done only once...
                     call h5screate_simple_f(rank4, dims, filespace_id, error)
                     call h5sselect_none_f(filespace_id, error)  ! empty filespace
                     call h5screate_simple_f(rank4, dims, memspace_id, error)
                     call h5sselect_none_f(memspace_id, error)   ! empty memoryscape

                     pa4d => null_r4d
                     call h5dwrite_f(cg_desc%dset_id(1, ic), H5T_NATIVE_DOUBLE, pa4d, dims, error, &
                        xfer_prp = cg_desc%xfer_prp, file_space_id = filespace_id, mem_space_id = memspace_id)
                     call h5sclose_f(memspace_id, error)
                     call h5sclose_f(filespace_id, error)
                  enddo
                  deallocate(dims)
               endif
            enddo
            ! receive (from whom?)
         else
            call die("[restart_hdf5_v2:write_cg_to_restart] nproc != nproc_io not implemented yet")
            ! send (where?)
         endif
      endif

      ! clean up
      call cg_desc%clean()
      deallocate(qr_lst, wr_lst, dsets)

      call ppp_main%stop(wrcg_label, PPP_IO)

   contains

      function pick_area(cg, mode) result(ijkse)

         use constants,  only: AT_OUT_B, AT_NO_B, AT_USER, ndims, LO, HI
         use dataio_pub, only: die

         implicit none

         type(grid_container), pointer, intent(in) :: cg
         integer(kind=4),               intent(in) :: mode

         integer(kind=4), dimension(ndims, LO:HI)  :: ijkse

         ijkse = 0 ! suppress compiler warning on maybe uninitialized variable
         select case (mode)
            case (AT_OUT_B)
               ijkse = cg%lh_out
            case (AT_NO_B)
               ijkse = cg%ijkse
            case (AT_USER)
               call die("[restart_hdf5_v2:write_cg_to_restart:pick_area] AT_USER not implemented")
            case default
               call die("[restart_hdf5_v2:write_cg_to_restart:pick_area] Non-recognized area_type.")
         end select
         return

      end function pick_area

      function pick_dims(cg, mode) result(n_b)

         use constants,  only: ndims, LO, HI

         implicit none

         type(grid_container), pointer, intent(in) :: cg
         integer(kind=4),               intent(in) :: mode

         integer(kind=4), dimension(ndims, LO:HI)  :: ijkse
         integer(kind=4), dimension(ndims) :: n_b

         ijkse = pick_area(cg, mode)
         n_b = ijkse(:,HI) - ijkse(:,LO) + I_ONE

      end function pick_dims

      function pick_size(ncg, mode) result(n_b)

         use constants,  only: AT_OUT_B, AT_NO_B, AT_USER, ndims, INVALID
         use dataio_pub, only: die

         implicit none

         integer(kind=4), intent(in) :: ncg
         integer(kind=4), intent(in) :: mode

         integer(kind=4), dimension(ndims) :: n_b

         n_b = INVALID
         select case (mode)
            case (AT_OUT_B)
               n_b = cg_all_n_o(:, ncg)
            case (AT_NO_B)
               n_b = cg_all_n_b(:, ncg)
            case (AT_USER)
               call die("[restart_hdf5_v2:write_cg_to_restart:pick_size] AT_USER not implemented")
            case default
               call die("[restart_hdf5_v2:write_cg_to_restart:pick_size] Non-recognized area_type.")
         end select

      end function pick_size

   end subroutine write_cg_to_restart

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
!! - Receive anything that had to be read on other processes
!! - Check for local completeness
!!
!! Optionally:
!! - When a single file is visible by multiple processes and the filesystem does not tolerate massively concurrent I/O load,
!!   choose one or few process to do the I/O and make the others to wait for the data.
!! - It should be possible to read runtime parameters from data stored in restart point rather than from actual problem.par file(s).
!!   A way to override some of them (in most cases  END_CONTROL::tend and nend) should be provided.
!<

   subroutine read_restart_hdf5_v2(status_v2)

      use cg_level_base,      only: base
      use cg_level_connected, only: cg_level_connected_t
      use cg_level_finest,    only: finest
      use cg_list,            only: cg_list_element
      use cmp_1D_mpi,         only: compare_array1D
      use common_hdf5,        only: d_gname, dir_pref, n_cg_name, d_size_aname, d_fc_aname, d_edge_apname, d_bnd_apname, &
           &                        cg_size_aname, cg_offset_aname, cg_lev_aname, base_d_gname, cg_cnt_aname, data_gname, &
           &                        output_fname, STAT_OK
      use constants,          only: cwdlen, dsetnamelen, cbuff_len, ndims, xdim, zdim, INVALID, RD, LO, HI, PPP_IO
      use dataio_pub,         only: die, warn, printio, msg, last_hdf_time, last_res_time, last_log_time, last_tsl_time, problem_name, res_id, domain_dump, &
           &                        require_problem_IC, piernik_hdf5_version2, nres, nhdf, fix_string
      use dataio_user,        only: user_reg_var_restart, user_attrs_rd
      use domain,             only: dom
      use func,               only: operator(.notequals.)
      use global,             only: t, dt, nstep
      use hdf5,               only: HID_T, H5F_ACC_RDONLY_F, h5open_f, h5close_f, h5fopen_f, h5fclose_f, h5gopen_f, h5gclose_f
      use mass_defect,        only: magic_mass
      use mpisetup,           only: master, piernik_MPI_Barrier
      use overlap,            only: is_overlap
      use ppp,                only: ppp_main
      use read_attr,          only: read_attribute
      use set_get_attributes, only: get_attr
      use timestep_pub,       only: c_all_old, cfl_c, stepcfl
#ifdef RANDOMIZE
      use randomization,      only: read_current_seed_from_restart
#endif /* RANDOMIZE */
#ifdef SN_SRC
      use snsources,          only: read_snsources_from_restart
#endif /* SN_SRC */
#if defined(MULTIGRID) && defined(SELF_GRAV)
      use multigrid_gravity,  only: read_oldsoln_from_restart
#endif /* MULTIGRID && SELF_GRAV */
#ifdef CRESP
      use initcrspectrum,     only: use_cresp
      use cresp_io,           only: cresp_read_smaps_from_hdf
#endif /* CRESP */

      implicit none

      integer, intent(out)                              :: status_v2

      integer(HID_T)                                    :: file_id              !< File identifier
      integer(HID_T)                                    :: doml_g_id, dom_g_id  !< domain list and domain group identifiers
      integer(HID_T)                                    :: cgl_g_id,  cg_g_id   !< cg list and cg group identifiers
      logical                                           :: file_exist, outside, overlapped
      integer                                           :: j
      integer(kind=4)                                   :: ia
      integer(kind=8)                                   :: tot_cells
      integer(kind=8), dimension(xdim:zdim, LO:HI)      :: my_box, other_box
      integer(kind=4)                                   :: error, nres_old  !< error perhaps should be of type integer(HID_T)

      ! buffers for reading attributes
      integer(kind=4),          dimension(:), allocatable :: ibuf
      real,                     dimension(:), allocatable :: rbuf
      character(len=cbuff_len), dimension(:), allocatable :: cbuf

      ! common attributes
      character(len=cbuff_len), dimension(10), parameter :: real_attrs = [ "time         ", &
           &                                                               "timestep     ", &
           &                                                               "last_hdf_time", &
           &                                                               "last_res_time", &
           &                                                               "last_log_time", &
           &                                                               "last_tsl_time", &
           &                                                               "c_all_old    ", &
           &                                                               "stepcfl      ", &
           &                                                               "cfl_c        ", &
           &                                                               "magic_mass   " ]
      character(len=cbuff_len), dimension(4), parameter :: int_attrs  = [ "nstep             ", &
           &                                                              "nres              ", &
           &                                                              "nhdf              ", &
           &                                                              "require_problem_IC" ]
      character(len=cbuff_len), dimension(3), parameter :: str_attrs  = [ "problem_name", &
           &                                                              "domain      ", &
           &                                                              "run_id      " ]
      !> \warning same strings are used independently in set_common_attributes*

      character(len=cwdlen)                             :: filename
      character(len=dsetnamelen)                        :: d_label
      type(cg_essentials), dimension(:), allocatable    :: cg_res
      type(cg_list_element), pointer                    :: cgl
      type(cg_level_connected_t), pointer               :: curl
      character(len=*), parameter :: rdr_label = "IO_read_restart_v2", rdrc_label = "IO_read_restart_v2_cg"

      call ppp_main%start(rdr_label, PPP_IO)
      if (master) call warn("[restart_hdf5_v2:read_restart_hdf5_v2] Experimental implementation")

      filename = output_fname(RD, '.res', nres, bcast=.true.)
      if (master) then
         write(msg, '(2a)') 'Reading restart file v2: ', trim(filename)
         call printio(msg)
      endif

      inquire(file = filename, exist = file_exist)
      if (.not. file_exist) then
         write(msg,'(3a)') '[restart_hdf5_v2:read_restart_hdf5_v2]: Restart file: ', trim(filename),' does not exist'
         call die(msg)
      endif

      status_v2 = INVALID

      call h5open_f(error)

      call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, error)

      ! Check file format version first

      call get_attr(file_id, "piernik", rbuf)

      if (size(rbuf) /= 1) call die("[restart_hdf5_v2:read_restart_hdf5_v2] Cannot read 'piernik' attribute from the restart file. The file may be either damaged or incompatible")
      if (rbuf(1) > piernik_hdf5_version2) then
         write(msg,'(2(a,f5.2))')"[restart_hdf5_v2:read_restart_hdf5_v2] Cannot read future versions of the restart file: ", rbuf(1)," > ", piernik_hdf5_version2
         call die(msg)
      else if (int(rbuf(1)) < int(piernik_hdf5_version2)) then
         write(msg,'(2(a,f5.2))')"[restart_hdf5_v2:read_restart_hdf5_v2] The restart file is too ancient. It is unlikely that it could work correctly: ", rbuf(1)," << ", piernik_hdf5_version2
         call warn(msg)
         call h5fclose_f(file_id, error)
         call h5close_f(error)
         return
      else if (rbuf(1) < piernik_hdf5_version2) then
         write(msg,'(2(a,f5.2))')"[restart_hdf5_v2:read_restart_hdf5_v2] Old versions of the restart file may not always work fully correctly: ", rbuf(1)," < ", piernik_hdf5_version2
         call warn(msg)
      endif
      call compare_array1D(rbuf(:))
      deallocate(rbuf)

      ! Compare attributes in the root of the restart point file with values read from problem.par
      !> \todo merge this code somehow with set_common_attributes_v2
      do ia = lbound(real_attrs, dim=1), ubound(real_attrs, dim=1)
         call get_attr(file_id, trim(real_attrs(ia)), rbuf)
         select case (real_attrs(ia))
            case ("time")
               t = rbuf(1)
            case ("timestep")
               dt = rbuf(1)
            case ("last_hdf_time")
               last_hdf_time = rbuf(1)
            case ("last_res_time")
               last_res_time = rbuf(1)
            case ("last_log_time")
               last_log_time = rbuf(1)
            case ("last_tsl_time")
               last_tsl_time = rbuf(1)
            case ("c_all_old")
               c_all_old = rbuf(1)
            case ("stepcfl")
               stepcfl = rbuf(1)
            case ("cfl_c")
               cfl_c = rbuf(1)
            case ("magic_mass")
               if (master) magic_mass(:) = rbuf(:)
            case default
               write(msg,'(3a,g15.5,a)')"[restart_hdf5_v2:read_restart_hdf5_v2] Real attribute '",trim(real_attrs(ia)),"' with value = ",rbuf(1)," was ignored"
               call warn(msg)
         end select
         deallocate(rbuf)
      enddo

      nres_old = nres
      do ia = lbound(int_attrs, dim=1), ubound(int_attrs, dim=1)
         call get_attr(file_id, trim(int_attrs(ia)), ibuf)
         select case (int_attrs(ia))
            case ("nstep")
               nstep = ibuf(1)
            case ("nres")
               nres = ibuf(1)
            case ("nhdf")
               nhdf = ibuf(1)
            case ("require_problem_IC")
               require_problem_IC = ibuf(1)
            case default
               write(msg,'(3a,i14,a)')"[restart_hdf5_v2:read_restart_hdf5_v2] Integer attribute '",trim(real_attrs(ia)),"' with value = ",ibuf(1)," was ignored"
               call warn(msg)
         end select
         deallocate(ibuf)
      enddo

#ifdef RANDOMIZE
      call read_current_seed_from_restart(file_id)
#endif /* RANDOMIZE */
#ifdef SN_SRC
      call read_snsources_from_restart(file_id)
#endif /* SN_SRC */
#if defined(MULTIGRID) && defined(SELF_GRAV)
      call read_oldsoln_from_restart(file_id)
#endif /* MULTIGRID && SELF_GRAV */

      do ia = lbound(str_attrs, dim=1), ubound(str_attrs, dim=1)
         call get_attr(file_id, trim(str_attrs(ia)), cbuf)
         select case (str_attrs(ia))
            case ("problem_name")
               problem_name = fix_string(trim(cbuf(1)))
            case ("domain")
               domain_dump = fix_string(trim(cbuf(1)(:len(domain_dump))))
            case ("run_id")
               res_id = fix_string(trim(cbuf(1)(:len(res_id))))
            case default
               write(msg,'(3a,i14,a)')"[restart_hdf5_v2:read_restart_hdf5_v2] String attribute '",trim(real_attrs(ia)),"' with value = ",cbuf," was ignored"
               call warn(msg)
         end select
         deallocate(cbuf)
      enddo

      !> \ToDo: try to simplify the code below and use get_attr in place of read_attribute

      ! Read domain description
      call h5gopen_f(file_id, d_gname, doml_g_id, error)       ! open "/domains"

      ! Read base domain
      call h5gopen_f(doml_g_id, base_d_gname, dom_g_id, error) ! open "/domains/base"

      allocate(ibuf(ndims))
      call read_attribute(dom_g_id, d_size_aname, ibuf)        ! read "/domains/base/n_d"
      if (any(ibuf(:) /= dom%n_d(:))) call die("[restart_hdf5_v2:read_restart_hdf5_v2] n_d doesn't match")
      deallocate(ibuf)

      allocate(ibuf(LO:HI), rbuf(LO:HI))
      do ia = xdim, zdim
         write(d_label, '(2a)') dir_pref(ia), d_edge_apname
         call read_attribute(dom_g_id, d_label, rbuf)          ! read "/domains/base/[xyz]-edge_position"
         if (any(rbuf(:).notequals.dom%edge(ia, :))) call die("[restart_hdf5_v2:read_restart_hdf5_v2] edge position doesn't match")
         write(d_label, '(2a)') dir_pref(ia), d_bnd_apname
         call read_attribute(dom_g_id, d_label, ibuf)          ! read "/domains/base/[xyz]-boundary_type"
         if (any(ibuf(:) /= dom%bnd(ia, :))) then
            write(msg,'(2a,2(a,2i3))')"[restart_hdf5_v2:read_restart_hdf5_v2] Boundary conditions in the ",dir_pref(ia),"-direction have changed. Saved type: ",ibuf(:), &
                 &                    ", new type:",dom%bnd(ia, :)
            call warn(msg)
         endif
      enddo
      deallocate(ibuf, rbuf)

      call h5gclose_f(dom_g_id, error)

      ! Read refined patches
      allocate(ibuf(1))
      call read_attribute(doml_g_id, d_fc_aname, ibuf)         ! read "/domains/fine_count"
      if (ibuf(1) > 0) call die("[restart_hdf5_v2:read_restart_hdf5_v2] /domains/fine_count /= 0 not implemented yet")
      deallocate(ibuf)

      call h5gclose_f(doml_g_id, error)

      ! Read available cg pieces
      call h5gopen_f(file_id, data_gname, cgl_g_id, error)         ! open "/data"

      allocate(ibuf(1))
      call read_attribute(cgl_g_id, cg_cnt_aname, ibuf)          ! open "/data/cg_count"
      if (ibuf(1) <= 0) call die("[restart_hdf5_v2:read_restart_hdf5_v2] Empty cg list")

      allocate(cg_res(ibuf(1)))
      deallocate(ibuf)

      do ia = lbound(cg_res, dim=1, kind=4), ubound(cg_res, dim=1, kind=4)
         call h5gopen_f(cgl_g_id, n_cg_name(ia), cg_g_id, error) ! open "/data/grid_%08d, ia"

         allocate(ibuf(1))
         call read_attribute(cg_g_id, cg_lev_aname, ibuf)        ! open "/data/grid_%08d/level"
         if (ibuf(1) < base%level%l%id) call die("[restart_hdf5_v2:read_restart_hdf5_v2] Grids coarser than base level are not supported")
         cg_res(ia)%level=ibuf(1)
         deallocate(ibuf)

         allocate(ibuf(ndims))
         call read_attribute(cg_g_id, cg_size_aname, ibuf)       ! open "/data/grid_%08d/n_b"
         if (any(ibuf(:)<=0)) call die("[restart_hdf5_v2:read_restart_hdf5_v2] Non-positive cg size detected")
         cg_res(ia)%n_b(:) = ibuf(:)
         call read_attribute(cg_g_id, cg_offset_aname, ibuf)     ! open "/data/grid_%08d/off"
         cg_res(ia)%off(:) = ibuf(:)
         deallocate(ibuf)

         call h5gclose_f(cg_g_id, error)
      enddo

      ! Check if whole base level is covered without any overlaps
      ! OPT: do overlap search on a separate list of base level blocks only because it is O(n^2) process
      tot_cells = 0
      outside = .false.
      overlapped = .false.
      do ia = lbound(cg_res, dim=1, kind=4), ubound(cg_res, dim=1, kind=4)
         if (cg_res(ia)%level == base%level%l%id) then
            tot_cells = tot_cells + product(cg_res(ia)%n_b(:))
            my_box(:,LO) = cg_res(ia)%off(:)
            my_box(:,HI) = cg_res(ia)%off(:) + cg_res(ia)%n_b(:) - 1
            outside = outside .or. any(my_box(:,LO) < 0) .or. any(my_box(:,HI) >= dom%n_d(:) .and. dom%has_dir(:))
            do j = ia+1, ubound(cg_res, dim=1)
               if (cg_res(j)%level == base%level%l%id) then
                  other_box(:,LO) = cg_res(j)%off(:)
                  other_box(:,HI) = cg_res(j)%off(:) + cg_res(j)%n_b(:) - 1
                  overlapped = overlapped .or. is_overlap(my_box, other_box)
               endif
            enddo
         endif
      enddo

      if (tot_cells /= product(dom%n_d(:)) .or. outside .or. overlapped) call die("[restart_hdf5_v2:read_restart_hdf5_v2] Improper coverage of base domain by available cg")

      call set_refinement(cg_res)

#ifdef CRESP
      if (use_cresp) call cresp_read_smaps_from_hdf(file_id)
#endif /* CRESP */

      ! set up things such as register user rank-3 and rank-4 arrays to be read by read_arr_from_restart. Read also anything that is not read by all read_cg_from_restart calls
      if (associated(user_attrs_rd)) call user_attrs_rd(file_id)
      if (associated(user_reg_var_restart)) call user_reg_var_restart

      ! Reset leafmap to allow checking if all existing grids were initialized
      ! This does not check whether all data on the grids in the restart file were used to initialize anything
      curl => base%level
      do while (associated(curl))
         cgl => curl%first
         do while (associated(cgl))
            cgl%cg%leafmap = .false.
            cgl => cgl%nxt
         enddo
         curl => curl%finer
      enddo

      call ppp_main%start(rdrc_label, PPP_IO)
      ! On each process determine which parts of the restart cg have to be read and where
      do ia = lbound(cg_res, dim=1, kind=4), ubound(cg_res, dim=1, kind=4)
         my_box(:,LO) = cg_res(ia)%off(:)
         my_box(:,HI) = cg_res(ia)%off(:) + cg_res(ia)%n_b(:) - 1
         curl => base%level
         do while (associated(curl))
            cgl => curl%first
            do while (associated(cgl))
               if (cg_res(ia)%level == curl%l%id) then
                  other_box(:, LO) = cgl%cg%my_se(:, LO) - curl%l%off(:)
                  other_box(:, HI) = cgl%cg%my_se(:, HI) - curl%l%off(:)
                  if (is_overlap(my_box, other_box)) call read_cg_from_restart(cgl%cg, cgl_g_id, ia, cg_res(ia))
               endif
               cgl => cgl%nxt
            enddo
            curl => curl%finer
         enddo
      enddo
      call ppp_main%stop(rdrc_label, PPP_IO)

      ! Use leafmap to check if all existing grids were initialized
      curl => base%level
      do while (associated(curl))
         cgl => curl%first
         do while (associated(cgl))
            if (.not. all(cgl%cg%leafmap)) call die("[restart_hdf5_v2:read_restart_hdf5_v2] Uninitialized grid found")
            cgl => cgl%nxt
         enddo
         call curl%sync_ru
         curl => curl%finer
      enddo

      call finest%level%restrict_to_base ! implies update of leafmap

      !> \todo update boundaries
      status_v2 = STAT_OK

      deallocate(cg_res)
      call h5gclose_f(cgl_g_id, error)

      call h5fclose_f(file_id, error)
      call h5close_f(error)

      if (status_v2 /= STAT_OK) nres = nres_old ! let's hope read_restart_hdf5_v1 will fix what could possibly got broken here
      call piernik_MPI_Barrier
      call ppp_main%stop(rdr_label, PPP_IO)

   end subroutine read_restart_hdf5_v2

!> \brief Find all non-base level grids in the restart file and distribute empty grids across the processes

   subroutine set_refinement(cg_res)

      use cg_leaves,          only: leaves
      use cg_level_base,      only: base
      use cg_level_connected, only: cg_level_connected_t
      use cg_level_finest,    only: finest
      use mpisetup,           only: master

      implicit none

      type(cg_essentials), dimension(:) :: cg_res

      integer :: lmax, i
      type(cg_level_connected_t), pointer :: curl

      lmax = base%level%l%id
      do i = lbound(cg_res(:), dim=1), ubound(cg_res(:), dim=1)
         if (cg_res(i)%level > lmax) lmax = cg_res(i)%level
      enddo

      if (lmax == base%level%l%id) return

      do while (lmax > finest%level%l%id)
         call finest%add_finer
      enddo

      ! Add all patches on master process and let the rebalancing routine do the work before the grids are created
      if (master) then
         curl => base%level%finer
         do while (associated(curl))
            do i = lbound(cg_res(:), dim=1), ubound(cg_res(:), dim=1)
               if (cg_res(i)%level == curl%l%id) call curl%add_patch(int(cg_res(i)%n_b, kind=8), cg_res(i)%off)
            enddo
            curl => curl%finer
         enddo
      endif

      curl => base%level%finer
      do while (associated(curl))
         call curl%init_all_new_cg
         curl => curl%finer
      enddo

      call leaves%update(" ( restart  ) ")

   end subroutine set_refinement

!> \brief Read as much as possible from stored cg to own cg
   subroutine read_cg_from_restart(cg, cgl_g_id, ncg, cg_r)

      use common_hdf5,      only: n_cg_name
      use constants,        only: xdim, ydim, zdim, ndims, LO, HI
      use dataio_pub,       only: die
      use domain,           only: dom
      use grid_cont,        only: grid_container
      use hdf5,             only: HID_T, HSIZE_T, H5S_SELECT_SET_F, H5T_NATIVE_DOUBLE, &
           &                      h5dopen_f, h5dclose_f, h5dget_space_f, h5dread_f, h5gopen_f, h5gclose_f, h5screate_simple_f, h5sselect_hyperslab_f
      use named_array_list, only: qna, wna
      use overlap,          only: is_overlap
#ifdef NBODY
      use common_hdf5,      only: part_gname, st_gname
      use data_hdf5,        only: gdf_translate
      use particles_io,     only: pvarn
      use particle_func,    only: particle_in_area
      use particle_utils,   only: is_part_in_cg
      use read_attr,        only: read_attribute
      use star_formation,   only: pid_gen
#endif /* NBODY */

      implicit none

      type(grid_container), pointer, intent(inout) :: cg        !< Own grid container
      integer(HID_T),                intent(in)    :: cgl_g_id  !< cg group identifier in the restart file
      integer(kind=4),               intent(in)    :: ncg       !< number of cg in the restart file
      type(cg_essentials),           intent(in)    :: cg_r      !< cg attributes that do not need to be reread

      integer(HID_T)                               :: cg_g_id !< cg group identifier
      integer(HID_T)                               :: dset_id
      integer(HID_T)                               :: filespace, memspace
      integer(HSIZE_T), dimension(:), allocatable  :: dims, off, cnt
      integer(kind=4)                              :: error                          !< error perhaps should be of type integer(HID_T)
      integer(kind=8), dimension(xdim:zdim)        :: own_off, restart_off, o_size   ! the position and size of the overlapped region
      integer(kind=8), dimension(xdim:zdim, LO:HI) :: own_box_nb, restart_box_nb              ! variants for AT_NO_B
      integer(kind=8), dimension(xdim:zdim)        :: own_off_nb, restart_off_nb, o_size_nb   !
      integer(kind=8), dimension(xdim:zdim, LO:HI) :: own_box_ob, restart_box_ob              ! variants for AT_OUT_B
      integer(kind=8), dimension(xdim:zdim)        :: own_off_ob, restart_off_ob, o_size_ob   ! as opposed to AT_NO_B
      integer, dimension(:), allocatable           :: qr_lst, wr_lst
      integer                                      :: i
      real, dimension(:,:,:),   allocatable        :: a3d
      real, dimension(:,:,:,:), allocatable        :: a4d
#ifdef NBODY
      integer(HID_T)                               :: part_g_id !< particles group identifier
      integer(HID_T)                               :: st_g_id   !< stars group identifier
      integer(HID_T)                               :: pdset_id
      integer(HSIZE_T)                             :: n_part
      integer(kind=8)                              :: j
      real,            dimension(:),   allocatable :: a1d
      real,            dimension(:),   allocatable :: mass, ener, tform, tdyn
      real,            dimension(:,:), allocatable :: pos, vel, acc
      integer(kind=4), dimension(:),   allocatable :: ibuf, pid
      logical                                      :: in, phy, out, fin
#endif /* NBODY */

      ! Find overlap between own cg and restart cg
      own_box_nb(:, :) = cg%my_se(:, :)
      restart_box_nb(:, LO) = cg%l%off(:) + cg_r%off(:)
      restart_box_nb(:, HI) = cg%l%off(:) + cg_r%off(:) + cg_r%n_b(:) - 1
      if (.not. is_overlap(own_box_nb, restart_box_nb)) call die("[restart_hdf5_v2:read_cg_from_restart] No overlap found") ! this condition should never happen

      call calc_off_and_size(restart_box_nb, own_box_nb, own_off_nb, restart_off_nb, o_size_nb)

      ! Find overlap between own cg and restart cg in case of AT_OUT_B
      own_box_ob(:, :) = cg%lh_out(:, :)
      restart_box_ob(:, LO) = cg%l%off(:) + cg_r%off(:)
      where (dom%has_dir(:) .and. (cg_r%off(:) <= 0)) restart_box_ob(:, LO) = restart_box_ob(:, LO) - dom%nb
      restart_box_ob(:, HI) = cg%l%off(:) + cg_r%off(:) + cg_r%n_b(:) - 1
      where (dom%has_dir(:) .and. (cg_r%off(:) + cg_r%n_b(:) >= cg%l%n_d)) restart_box_ob(:, HI) = restart_box_ob(:, HI) + dom%nb
      if (.not. is_overlap(own_box_ob, restart_box_ob)) call die("[restart_hdf5_v2:read_cg_from_restart] No overlap found (AT_OUT_B)") ! this condition should never happen

      call calc_off_and_size(restart_box_ob, own_box_ob, own_off_ob, restart_off_ob, o_size_ob)
      where (dom%has_dir(:) .and. (cg%lh_out(:, LO) < cg%my_se(:, LO))) own_off_ob(:) = own_off_ob(:) - dom%nb
      ! An extra correction for cg with low outer boundary. Blocks with high outer boundary will get it just through o_size_ob

      ! these conditions should never happen
      if (any(own_off_nb(:) > cg%n_b(:))) call die("[restart_hdf5_v2:read_cg_from_restart] own_off(:) > cg%n_b(:)")
      if (any(restart_off_nb(:) > cg_r%n_b(:))) call die("[restart_hdf5_v2:read_cg_from_restart] restart_off(:) > cg_r%n_b(:)")
      if (any(o_size_nb(:) > cg%n_b(:)) .or. any(o_size_nb(:) > cg_r%n_b(:))) call die("[restart_hdf5_v2:read_cg_from_restart] o_size(:) > cg%n_b(:) or o_size(:) > cg_r%n_b(:)")
      if (any(cg%leafmap(cg%is+own_off_nb(xdim):cg%is+own_off_nb(xdim)+o_size_nb(xdim)-1, &
           &             cg%js+own_off_nb(ydim):cg%js+own_off_nb(ydim)+o_size_nb(ydim)-1, &
           &             cg%ks+own_off_nb(zdim):cg%ks+own_off_nb(zdim)+o_size_nb(zdim)-1))) call die("[restart_hdf5_v2:read_cg_from_restart] Trying to initialize same area twice.")

      call qna%get_reslst(qr_lst)
      call wna%get_reslst(wr_lst)
      call h5gopen_f(cgl_g_id, n_cg_name(ncg), cg_g_id, error) ! open "/data/grid_%08d, ncg"

      if (size(qr_lst) > 0) then
         allocate(dims(ndims), off(ndims), cnt(ndims))
         do i = lbound(qr_lst, dim=1, kind=4), ubound(qr_lst, dim=1, kind=4)
            call pick_off_and_size(qna%lst(qr_lst(i))%restart_mode, o_size, restart_off, own_off)
            dims(:) = o_size(:)
            call h5dopen_f(cg_g_id, qna%lst(qr_lst(i))%name, dset_id, error) ! open "/data/grid_%08d/cg%q(qr_lst(i))%name"
            call h5dget_space_f(dset_id, filespace, error)
            off(:) = restart_off(:)
            cnt(:) = 1
            call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, off(:), cnt(:), error, block=dims(:))
            call h5screate_simple_f(size(dims, kind=4), dims(:), memspace, error)
            allocate(a3d(dims(xdim), dims(ydim), dims(zdim)))
            call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, a3d, dims(:), error, file_space_id = filespace, mem_space_id = memspace)
            cg%q(qr_lst(i))%arr(cg%is+own_off(xdim):cg%is+own_off(xdim)+o_size(xdim)-1, &
                 &              cg%js+own_off(ydim):cg%js+own_off(ydim)+o_size(ydim)-1, &
                 &              cg%ks+own_off(zdim):cg%ks+own_off(zdim)+o_size(zdim)-1) = a3d(:,:,:)
            deallocate(a3d)
            call h5dclose_f(dset_id, error)
         enddo
         deallocate(dims, off, cnt)
      endif

      if (size(wr_lst) > 0) then
         allocate(dims(ndims+1), off(ndims+1), cnt(ndims+1))
         do i = lbound(wr_lst, dim=1, kind=4), ubound(wr_lst, dim=1, kind=4)
            call pick_off_and_size(wna%lst(wr_lst(i))%restart_mode, o_size, restart_off, own_off)
            dims(:) = [ int(wna%lst(wr_lst(i))%dim4, kind=HSIZE_T), int(o_size(:), kind=HSIZE_T) ]
            call h5dopen_f(cg_g_id, wna%lst(wr_lst(i))%name, dset_id, error)
            call h5dget_space_f(dset_id, filespace, error)
            off(:) = [ 0_HSIZE_T, restart_off(:) ]
            cnt(:) = 1
            call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, off(:), cnt(:), error, block=dims(:))
            call h5screate_simple_f(size(dims, kind=4), dims(:), memspace, error)
            allocate(a4d(dims(1), dims(1+xdim), dims(1+ydim), dims(1+zdim)))
            call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, a4d, dims(:), error, file_space_id = filespace, mem_space_id = memspace)
            cg%w(wr_lst(i))%arr(:, cg%is+own_off(xdim):cg%is+own_off(xdim)+o_size(xdim)-1, &
                 &                 cg%js+own_off(ydim):cg%js+own_off(ydim)+o_size(ydim)-1, &
                 &                 cg%ks+own_off(zdim):cg%ks+own_off(zdim)+o_size(zdim)-1) = a4d(:,:,:,:)
            deallocate(a4d)
            call h5dclose_f(dset_id, error)
         enddo
         deallocate(dims, off, cnt)
      endif

#ifdef NBODY
      call h5gopen_f(cg_g_id,   part_gname, part_g_id, error)
      call h5gopen_f(part_g_id,   st_gname,   st_g_id, error)

      allocate(ibuf(1))
      call read_attribute(st_g_id, "n_part",  ibuf)
      n_part = ibuf(1)
      call read_attribute(st_g_id, "pid_max", ibuf)
      pid_gen = ibuf(1)
      deallocate(ibuf)

      allocate(pid(n_part), mass(n_part), ener(n_part), tform(n_part), tdyn(n_part))
      allocate(pos(n_part, ndims), vel(n_part, ndims), acc(n_part, ndims))
      allocate(a1d(n_part))
      do i = lbound(pvarn, dim=1), ubound(pvarn, dim=1)
         call h5dopen_f(st_g_id, gdf_translate(pvarn(i)), pdset_id, error)
         call h5dread_f(pdset_id, H5T_NATIVE_DOUBLE, a1d, [n_part], error)
         select case (pvarn(i))
            case ('ppid')
               pid = int(a1d, kind=4)
            case ('mass')
               mass = a1d
            case ('ener')
               ener = a1d
            case ('posx')
               pos(:, xdim) = a1d
            case ('posy')
               pos(:, ydim) = a1d
            case ('posz')
               pos(:, zdim) = a1d
            case ('velx')
               vel(:, xdim) = a1d
            case ('vely')
               vel(:, ydim) = a1d
            case ('velz')
               vel(:, zdim) = a1d
            case ('accx')
               acc(:, xdim) = a1d
            case ('accy')
               acc(:, ydim) = a1d
            case ('accz')
               acc(:, zdim) = a1d
            case ('tfor')
               tform = a1d
            case ('tdyn')
               tdyn = a1d
            case default
         end select
         call h5dclose_f(pdset_id, error)
      enddo
      deallocate(a1d)
      call h5gclose_f(st_g_id, error)
      call h5gclose_f(part_g_id, error)

      do j = 1, n_part
         call is_part_in_cg(cg, pos(j,:), particle_in_area(pos(j,:), dom%edge), in, phy, out, fin)
         call cg%pset%add(pid(j), mass(j), pos(j,:), vel(j,:), acc(j,:), ener(j), in, phy, out, .true., tform(j), tdyn(j))
      enddo
      deallocate(pid, mass, pos, vel, acc, ener, tform, tdyn)
#endif /* NBODY */

      call h5gclose_f(cg_g_id, error)
      deallocate(qr_lst, wr_lst)

      ! Mark the area as initialized
      cg%leafmap(cg%is+own_off_nb(xdim):cg%is+own_off_nb(xdim)+o_size_nb(xdim)-1, &
           &     cg%js+own_off_nb(ydim):cg%js+own_off_nb(ydim)+o_size_nb(ydim)-1, &
           &     cg%ks+own_off_nb(zdim):cg%ks+own_off_nb(zdim)+o_size_nb(zdim)-1) = .true.

   contains

      subroutine pick_off_and_size(mode, o_size, restart_off, own_off)

         use constants,  only: AT_OUT_B, AT_NO_B, AT_USER
         use dataio_pub, only: die

         implicit none

         integer(kind=4),                       intent(in)  :: mode
         integer(kind=8), dimension(xdim:zdim), intent(out) :: o_size
         integer(kind=8), dimension(xdim:zdim), intent(out) :: restart_off
         integer(kind=8), dimension(xdim:zdim), intent(out) :: own_off

         ! suppress compiler warnings on possibly use of uninitialized values
         o_size = 0
         restart_off = 0
         own_off = 0

         select case (mode)
            case (AT_OUT_B)
               o_size = o_size_ob
               restart_off = restart_off_ob
               own_off = own_off_ob
            case (AT_NO_B)
               o_size = o_size_nb
               restart_off = restart_off_nb
               own_off = own_off_nb
            case (AT_USER)
               call die("[restart_hdf5_v2:read_cg_from_restart:pick_off_and_size] AT_USER not implemented (w)")
            case default
               call die("[restart_hdf5_v2:read_cg_from_restart:pick_off_and_size] Non-recognized area_type. (w)")
         end select

      end subroutine pick_off_and_size

      subroutine calc_off_and_size(restart_box_my, own_box_my, own_off_my, restart_off_my, o_size_my)

         use constants, only: LONG, LO, HI, xdim, zdim

         implicit none

         integer(kind=8), dimension(xdim:zdim, LO:HI), intent(in)  :: restart_box_my
         integer(kind=8), dimension(xdim:zdim, LO:HI), intent(in)  :: own_box_my
         integer(kind=8), dimension(xdim:zdim),        intent(out) :: own_off_my
         integer(kind=8), dimension(xdim:zdim),        intent(out) :: restart_off_my
         integer(kind=8), dimension(xdim:zdim),        intent(out) :: o_size_my

         integer :: d

         own_off_my(:) = 0
         restart_off_my(:) = 0
         o_size_my(:) = 1
         do d = xdim, zdim
            if (dom%has_dir(d)) then
               own_off_my(d) = max(restart_box_my(d, LO) - own_box_my(d, LO), 0_LONG)
               restart_off_my(d) = max(own_box_my(d, LO) - restart_box_my(d, LO), 0_LONG)
               o_size_my(d) = min(restart_box_my(d, HI), own_box_my(d, HI)) - max(restart_box_my(d, LO), own_box_my(d, LO)) + 1
            endif
         enddo

      end subroutine calc_off_and_size

   end subroutine read_cg_from_restart

end module restart_hdf5_v2
