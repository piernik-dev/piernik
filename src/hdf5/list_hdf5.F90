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
#ifdef NEW_HDF5
   use hdf5, only: HID_T
#endif /* NEW_HDF5 */
   private
   public :: write_arr, S_LEN, additional_attrs, problem_write_restart, problem_read_restart
#ifdef NEW_HDF5
   public :: add_lhdf5, lhdf5_info, iterate_lhdf5
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
      character(len=S_LEN)  :: key
      integer, dimension(6) :: ind
      integer, dimension(3) :: sz
      integer, dimension(5) :: opti
      procedure(func2) ,pointer, nopass :: p
   end type lhdf5_info

   abstract interface
      function func2(ind,sz,opti)
         implicit none
         integer, dimension(6), intent(in)    :: ind
         integer, dimension(3), intent(in)    :: sz
         integer, dimension(5), intent(in)    :: opti
         real, dimension(sz(1),sz(2),sz(3))   :: func2
      end function func2
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
         integer(HID_T) :: file_id       !> File identifier
         type(lhdf5_list), pointer :: tp
         tp => lhdf5_root
         do while (associated(tp%next))
            call get_lhdf5(tp%next,file_id)
            tp => tp%next%node
         enddo

         contains

            subroutine get_lhdf5(tp,file_id)
               implicit none
               integer(HID_T) :: file_id       !> File identifier
               type(lhdf5_node), pointer :: tp
               real, dimension(:,:,:),allocatable :: a
               allocate(a(tp%info%sz(1),tp%info%sz(2),tp%info%sz(3)))
               a = tp%info%p(tp%info%ind,tp%info%sz,tp%info%opti)
               call write_arr(real(a,4),tp%info%key,file_id)
               deallocate(a)
               return
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
            temp%info%sz   = item%sz
            temp%info%ind  = item%ind
            temp%info%opti = item%opti

            temp%info%p    => item%p
            temp%node%next => tp

            tp => temp

            return
         end subroutine insert_lhdf5

      end subroutine add_lhdf5

#endif /* NEW_HDF5 */

     subroutine write_arr(data,dsetname,file_id)
         use hdf5, only: HID_T, HSIZE_T, HSSIZE_T, H5FD_MPIO_INDEPENDENT_F, H5P_DATASET_CREATE_F, H5P_DATASET_XFER_F, H5S_SELECT_SET_F, &
            H5T_NATIVE_REAL, h5dwrite_f, h5screate_simple_f, h5pcreate_f, h5pset_chunk_f, h5dcreate_f, h5sclose_f, &
            h5dget_space_f, h5sselect_hyperslab_f, h5pset_dxpl_mpio_f, h5dclose_f
         use mpisetup, only: psize, pcoords

         implicit none
         real(kind=4), dimension(:,:,:) :: data

         integer, parameter :: rank = 3 ! Dataset rank

         CHARACTER(LEN=S_LEN) :: dsetname    !< Dataset name

         integer(HID_T) :: file_id       !< File identifier
         integer(HID_T) :: dset_id       !< Dataset identifier
         integer(HID_T) :: plist_id      !< Property list identifier
         integer(HID_T) :: filespace     !< Dataspace identifier in file
         integer(HID_T) :: memspace      !< Dataspace identifier in memory

         integer, parameter :: ndims = 3
         integer(HSIZE_T),  DIMENSION(ndims) :: count
         integer(HSSIZE_T), DIMENSION(ndims) :: offset
         integer(HSIZE_T),  DIMENSION(ndims) :: stride
         integer(HSIZE_T),  DIMENSION(ndims) :: block
         integer(HSIZE_T), DIMENSION(ndims) :: dimsf, dimsfi, chunk_dims
         integer :: error

         dimsf = (/size(data,1)*psize(1),size(data,2)*psize(2),size(data,3)*psize(3)/) ! Dataset dimensions
         dimsfi = dimsf
         chunk_dims = (/size(data,1),size(data,2),size(data,3)/) ! Chunks dimensions
         !
         ! Create the data space for the  dataset.
         !
         CALL h5screate_simple_f(rank, dimsf, filespace, error)
         CALL h5screate_simple_f(rank, chunk_dims, memspace, error)

         !
         ! Create chunked dataset.
         !
         CALL h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
         CALL h5pset_chunk_f(plist_id, rank, chunk_dims, error)
         CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_real, filespace, &
                      dset_id, error, plist_id)
         CALL h5sclose_f(filespace, error)

         !
         ! Each process defines dataset in memory and writes it to the hyperslab
         ! in the file.
         !
         stride(:) = 1
         count(:) =  1
         block(:) = chunk_dims(:)

         offset(:) = pcoords(:)*chunk_dims(:)
         !
         ! Select hyperslab in the file.
         !
         CALL h5dget_space_f(dset_id, filespace, error)
         CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error, &
                                 stride, block)
         !
         ! Create property list for collective dataset write
         !
         CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
         CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)

         !
         ! Write the dataset collectively.
         !
         CALL h5dwrite_f(dset_id, H5T_NATIVE_real, data, dimsfi, error, &
                     file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

         !
         ! Close dataspaces.
         !
         CALL h5sclose_f(filespace, error)
         CALL h5sclose_f(memspace, error)
         !
         ! Close the dataset.
         !
         CALL h5dclose_f(dset_id, error)

      end subroutine write_arr
end module list_hdf5
