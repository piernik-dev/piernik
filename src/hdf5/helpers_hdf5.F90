! $Id$
!
! PIERNIK Code Copyright (C) 2006-2012 Michal Hanasz
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
!! \brief Module that contains overloaded wrappers for HDF5 routines (create attribute, dataset etc)
!<
module helpers_hdf5

   implicit none

   private
   public :: create_attribute, create_dataset

!> \brief Add an attribute (1D array) to the given _id and initialize its value
   interface create_attribute
      module procedure create_int_attribute
      module procedure create_str_attribute
      module procedure create_int8_attribute
      module procedure create_real_attribute
   end interface

   interface create_dataset
      module procedure create_dataset_int4_dim2
      module procedure create_dataset_int8_dim2
      module procedure create_dataset_int8_dim1
   end interface

contains

!> \brief Create 32-bit integer dataset (rank-2 array) in the given place_id.
!
   subroutine create_dataset_int4_dim2(place, dname, ddata)

      use constants,     only: I_TWO
      use iso_c_binding, only: c_ptr, c_loc
      use hdf5,          only: HID_T, HSIZE_T, H5T_STD_I32LE, &
          &                    h5dcreate_f, h5dclose_f, h5screate_simple_f, h5sclose_f, h5dwrite_f, &
          &                    h5kind_to_type, H5_INTEGER_KIND

      implicit none

      integer(HID_T), intent(in)   :: place !> object id where dataset will be created
      character(len=*), intent(in) :: dname !> name of dataset
      integer(kind=4), dimension(:,:), pointer, intent(in) :: ddata !> data used to create dataset

      integer(HID_T)  :: dset, space, mem_type
      integer(kind=4) :: hdferr
      integer(HSIZE_T), dimension(2) :: dims
      type(c_ptr) :: f_ptr

      dims = shape(ddata)
      call h5screate_simple_f(I_TWO, dims, space, hdferr)
      call h5dcreate_f(place, dname, H5T_STD_I32LE, space, dset, hdferr)
      f_ptr = c_loc(ddata(1,1))
      mem_type = h5kind_to_type(int(KIND(ddata(1,1)), kind=4), H5_INTEGER_KIND)
      call h5dwrite_f(dset, mem_type, f_ptr, hdferr)
      call h5dclose_f(dset,  hdferr)
      call h5sclose_f(space, hdferr)

   end subroutine create_dataset_int4_dim2

!> \brief Create 64-bit integer dataset (rank-2 array) in the given place_id.
!
   subroutine create_dataset_int8_dim2(place, dname, ddata)

      use constants,     only: I_TWO
      use iso_c_binding, only: c_ptr, c_loc
      use hdf5,          only: HID_T, HSIZE_T, H5T_STD_I64LE, &
          &                    h5dcreate_f, h5dclose_f, h5screate_simple_f, h5sclose_f, h5dwrite_f, &
          &                    h5kind_to_type, H5_INTEGER_KIND

      implicit none

      integer(HID_T), intent(in)   :: place !> object id where dataset will be created
      character(len=*), intent(in) :: dname !> name of dataset
      integer(kind=8), dimension(:,:), pointer, intent(in) :: ddata !> data used to create dataset

      integer(HID_T)  :: dset, space, mem_type
      integer(kind=4) :: hdferr
      integer(HSIZE_T), dimension(2) :: dims
      type(c_ptr) :: f_ptr

      dims = shape(ddata)
      call h5screate_simple_f(I_TWO, dims, space, hdferr)
      call h5dcreate_f(place, dname, H5T_STD_I64LE, space, dset, hdferr)
      f_ptr = c_loc(ddata(1,1))
      mem_type = h5kind_to_type(int(KIND(ddata(1,1)), kind=4), H5_INTEGER_KIND)
      call h5dwrite_f(dset, mem_type, f_ptr, hdferr)
      call h5dclose_f(dset,  hdferr)
      call h5sclose_f(space, hdferr)

   end subroutine create_dataset_int8_dim2

!> \brief Create 64-bit integer dataset (rank-1 array) in the given place_id.
!
   subroutine create_dataset_int8_dim1(place, dname, ddata)

      use constants,     only: I_ONE
      use iso_c_binding, only: c_ptr, c_loc
      use hdf5,          only: HID_T, HSIZE_T, H5T_STD_I64LE, &
          &                    h5dcreate_f, h5dclose_f, h5screate_simple_f, h5sclose_f, h5dwrite_f, &
          &                    h5kind_to_type, H5_INTEGER_KIND

      implicit none

      integer(HID_T), intent(in)   :: place !> object id where dataset will be created
      character(len=*), intent(in) :: dname !> name of dataset
      integer(kind=8), dimension(:), pointer, intent(in) :: ddata !> data used to create dataset

      integer(HID_T)  :: dset, space, mem_type
      integer(kind=4) :: hdferr
      integer(HSIZE_T), dimension(1) :: dims
      type(c_ptr) :: f_ptr

      dims = shape(ddata)
      call h5screate_simple_f(I_ONE, dims, space, hdferr)
      call h5dcreate_f(place, dname, H5T_STD_I64LE, space, dset, hdferr)
      f_ptr = c_loc(ddata(1))
      mem_type = h5kind_to_type(int(KIND(ddata(1)), kind=4), H5_INTEGER_KIND)
      call h5dwrite_f(dset, mem_type, f_ptr, hdferr)
      call h5dclose_f(dset,  hdferr)
      call h5sclose_f(space, hdferr)

   end subroutine create_dataset_int8_dim1

!> \brief Attach an 32-bit integer attribute (scalar or rank-1 small array) to the given group.
!
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

!> \brief Attach an 64-bit integer attribute (scalar or rank-1 small array) to the given group.
!
   subroutine create_int8_attribute(g_id, name, int_array)

      use iso_c_binding, only: c_ptr, c_loc
      use constants, only: I_ONE
      use hdf5,      only: HID_T, HSIZE_T, h5kind_to_type, H5_INTEGER_KIND, H5T_STD_I64LE, &
          &               h5acreate_f, h5aclose_f, h5awrite_f, h5screate_simple_f, h5sclose_f

      implicit none

      integer(HID_T), intent(in)                         :: g_id      !< group id where to create the attribute
      character(len=*), intent(in)                       :: name      !< name
      integer(kind=8), dimension(:), pointer, intent(in) :: int_array !< the data

      integer(HID_T)  :: attr, space, mem_type
      integer(kind=4) :: hdferr
      integer(HSIZE_T), dimension(I_ONE) :: dims
      type(c_ptr) :: f_ptr

      dims = shape(int_array)
      call h5screate_simple_f(I_ONE, dims, space, hdferr)
      call h5acreate_f(g_id, name, H5T_STD_I64LE, space, attr, hdferr)
      f_ptr = c_loc(int_array(1))
      mem_type = h5kind_to_type(int(KIND(int_array(1)), kind=4), H5_INTEGER_KIND)
      call h5awrite_f(attr, mem_type, f_ptr, hdferr)
      call h5aclose_f(attr,  hdferr)
      call h5sclose_f(space, hdferr)

   end subroutine create_int8_attribute

!> \brief Attach an 64-bit real attribute (scalar or rank-1 small array) to the given group.
!
   subroutine create_real_attribute(g_id, name, real_array)

     use constants, only: I_ONE
     use hdf5,      only: H5T_NATIVE_DOUBLE, HID_T, HSIZE_T, &
          &               h5acreate_f, h5aclose_f, h5awrite_f, h5screate_simple_f, h5sclose_f

     implicit none

     integer(HID_T), intent(in)     :: g_id       !< group id where to create the attribute
     character(len=*), intent(in)   :: name       !< name
     real(kind=8), dimension(:), intent(in) :: real_array !< the data

     integer(HID_T)  :: aspace_id, attr_id
     integer(kind=4) :: error

     call h5screate_simple_f(I_ONE, [ size(real_array, kind=HSIZE_T) ], aspace_id, error)
     call h5acreate_f(g_id, name, H5T_NATIVE_DOUBLE, aspace_id, attr_id, error)
     call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, real_array, [ size(real_array, kind=HSIZE_T) ], error)
     call h5aclose_f(attr_id, error)
     call h5sclose_f(aspace_id, error)

   end subroutine create_real_attribute

!> \brief Attach an string attribute (scalar) to the given group.
!
   subroutine create_str_attribute(g_id, name, data)

      use iso_c_binding, only: c_loc, c_ptr
      use constants, only: I_ONE
      use hdf5,      only: HID_T, HSIZE_T, SIZE_T, &
           &               h5acreate_f, h5aclose_f, h5awrite_f, h5screate_simple_f, h5sclose_f, &
           &               H5Tcopy_f, H5T_C_S1, H5Tset_size_f, H5T_FORTRAN_S1, H5Tset_size_f, H5tclose_f

      implicit none

      integer(HID_T), intent(in)     :: g_id       !< group id where to create the attribute
      character(len=*), intent(in)   :: name       !< name
      character(len=*), intent(in), pointer :: data !< the data

      integer(HID_T)  :: space, attr_id, memtype, filetype
      integer(kind=4) :: error
      type(c_ptr) :: f_ptr

      call H5Tcopy_f(H5T_C_S1, filetype, error)
      call H5Tset_size_f(filetype, int(len(data)+1, SIZE_T), error)
      call H5Tcopy_f( H5T_FORTRAN_S1, memtype, error)
      call H5Tset_size_f(memtype, int(len(data), SIZE_T), error)

      call h5screate_simple_f(I_ONE, [ int(1, kind=HSIZE_T) ], space, error)
      call H5Acreate_f(g_id, name, filetype, space, attr_id, error)
      f_ptr = C_LOC(data(1:1))
      call H5Awrite_f(attr_id, memtype, f_ptr, error)
      call h5aclose_f(attr_id, error)
      call h5sclose_f(space, error)
      call H5tclose_f(filetype, error)
      call H5tclose_f(memtype, error)

   end subroutine create_str_attribute
end module helpers_hdf5
