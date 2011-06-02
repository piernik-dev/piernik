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

module list_hdf5
! pulled by ANY

#ifdef NEW_HDF5
   use iso_c_binding, only: c_int
#endif /* NEW_HDF5 */

   implicit none

   private
   public :: write_arr, S_LEN, additional_attrs, problem_write_restart, problem_read_restart
#ifdef NEW_HDF5
   public :: add_lhdf5, lhdf5_info, iterate_lhdf5
   integer, parameter :: HID_T = c_int
#endif /* NEW_HDF5 */

   integer, parameter :: S_LEN = 30

   interface
      subroutine add_attr(file_id)
         use hdf5, only: HID_T
         implicit none
         integer(HID_T), intent(in) :: file_id
      end subroutine add_attr
   end interface

   procedure(add_attr), pointer :: additional_attrs
   procedure(add_attr), pointer :: problem_write_restart
   procedure(add_attr), pointer :: problem_read_restart

#ifdef NEW_HDF5
   type :: lhdf5_info
      character(len=S_LEN)                :: key
      real,    dimension(:), allocatable  :: rvec
      integer, dimension(:), allocatable  :: ivec
      procedure(get_tab) ,pointer, nopass :: p
   end type lhdf5_info

   abstract interface
      subroutine get_tab(ivec,rvec,tab)
         implicit none
         integer,  dimension(:), allocatable, intent(in)  :: ivec
         real,     dimension(:), allocatable, intent(in)  :: rvec
         real, dimension(:,:,:), allocatable, intent(out) :: tab
      end subroutine get_tab
   end interface

   type :: lhdf5_list
      type(lhdf5_node), pointer :: next => Null()
   end type lhdf5_list

   type :: lhdf5_node
      type(lhdf5_info) :: info
      type(lhdf5_list) :: node
   end type lhdf5_node

   type(lhdf5_list), target, save :: lhdf5_root
#endif /* NEW_HDF5 */

contains

#ifdef NEW_HDF5

   subroutine iterate_lhdf5(file_id)

      implicit none

      integer(HID_T)            :: file_id       !> File identifier
      type(lhdf5_list), pointer :: tp

      tp => lhdf5_root
      do while (associated(tp%next))
         call get_lhdf5(tp%next,file_id)
         tp => tp%next%node
      enddo

   contains

      subroutine get_lhdf5(tp,file_id)

         implicit none

         integer(HID_T)                      :: file_id       !> File identifier
         type(lhdf5_node), pointer           :: tp
         real, dimension(:,:,:), allocatable :: a

         call tp%info%p(tp%info%ivec,tp%info%rvec,a)
         call write_arr(real(a,4),tp%info%key,file_id)
         if (allocated(a)) deallocate(a)

      end subroutine get_lhdf5

   end subroutine iterate_lhdf5

   subroutine add_lhdf5(item)

      use dataio_pub,    only: msg, warn

      implicit none

      type(lhdf5_info), intent(in) :: item
      type(lhdf5_list), pointer :: tp

      tp => lhdf5_root
      do
         if ( associated(tp%next)) then
            if ( item%key == tp%next%info%key ) then
               write(msg,'(3a)') "[list_hdf5:add_lhdf5]: ",trim(item%key)," exists in the list"
               call warn(msg)
               return
            else if ( item%key < tp%next%info%key) then
               call insert_lhdf5(tp%next, item)
            else
               tp => tp%next%node
               cycle ! keep looking
            endif
         else
            call insert_lhdf5(tp%next,item)
         endif
         return
      enddo

   contains

      subroutine insert_lhdf5(tp, item)

         implicit none

         type(lhdf5_node), pointer :: tp
         type(lhdf5_info), intent(in) :: item
         type(lhdf5_node), pointer :: temp

         allocate(temp)
         temp%info%key  = item%key
         if (allocated(item%ivec)) then
            allocate(temp%info%ivec(size(item%ivec)))
            temp%info%ivec = item%ivec
         endif
         if (allocated(item%rvec)) then
            allocate(temp%info%rvec(size(item%rvec)))
            temp%info%rvec = item%rvec
         endif

         temp%info%p    => item%p
         temp%node%next => tp

         tp => temp

      end subroutine insert_lhdf5

   end subroutine add_lhdf5

#endif /* NEW_HDF5 */

   !> \todo this routine should be merged with dataio_hdf5:write_arr_to_restart

   subroutine write_arr(data, dsetname, file_id)

      use constants,  only: ndims
      use grid,       only: cga
      use grid_cont,  only: cg_list_element, grid_container
      use hdf5,       only: HID_T, HSIZE_T, H5FD_MPIO_COLLECTIVE_F, H5P_DATASET_CREATE_F, H5P_DATASET_XFER_F, H5S_SELECT_SET_F, &
           &                H5T_NATIVE_REAL, h5dwrite_f, h5screate_simple_f, h5pcreate_f, h5pset_chunk_f, h5dcreate_f, h5sclose_f, &
           &                h5dget_space_f, h5sselect_hyperslab_f, h5pset_dxpl_mpio_f, h5dclose_f
      use mpisetup,   only: dom, is_uneven

      implicit none

      real(kind=4), dimension(:,:,:) :: data    !< array of output values
      character(len=*) :: dsetname              !< Dataset name
      integer(HID_T) :: file_id                 !< File identifier

      integer, parameter :: rank = ndims        !< Dataset rank = 3
      integer(HID_T) :: dset_id                 !< Dataset identifier
      integer(HID_T) :: plist_id                !< Property list identifier
      integer(HID_T) :: filespace               !< Dataspace identifier in file
      integer(HID_T) :: memspace                !< Dataspace identifier in memory
      integer(HSIZE_T), dimension(rank) :: count, offset, stride, block, dimsf, chunk_dims
      integer :: error
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      cgl => cga%cg_leafs%cg_l(1)
      do while (associated(cgl))
         cg => cgl%cg

         chunk_dims = cg%n_b(:) ! Chunks dimensions
         dimsf  = dom%n_d(:)    ! Dataset dimensions
         !
         ! Create the data space for the  dataset.
         !
         call h5screate_simple_f(rank, dimsf, filespace, error)
         call h5screate_simple_f(rank, chunk_dims, memspace, error)

         !
         ! Create chunked dataset.
         !
         call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
         if (.not. is_uneven) call h5pset_chunk_f(plist_id, rank, chunk_dims, error) !> \todo check how much performance it gives (massively parallel I/O is required)
         call h5dcreate_f(file_id, dsetname, H5T_NATIVE_REAL, filespace, dset_id, error, plist_id)
         call h5sclose_f(filespace, error)

         !
         ! Each process defines dataset in memory and writes it to the hyperslab
         ! in the file.
         !
         stride(:) = 1
         count(:) =  1
         block(:) = chunk_dims(:)

         offset(:) = cg%off(:)
         !
         ! Select hyperslab in the file.
         !
         call h5dget_space_f(dset_id, filespace, error)
         call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error, stride, block)
         !
         ! Create property list for collective dataset write
         !
         call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
         call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

         !
         ! Write the dataset collectively.
         !
         call h5dwrite_f(dset_id, H5T_NATIVE_REAL, data, dimsf, error, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

         !
         ! Close dataspaces.
         !
         call h5sclose_f(filespace, error)
         call h5sclose_f(memspace, error)
         !
         ! Close the dataset.
         !
         call h5dclose_f(dset_id, error)

         cgl => cgl%nxt
      enddo

   end subroutine write_arr

end module list_hdf5
