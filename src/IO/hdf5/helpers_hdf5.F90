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

   use hdf5, only: SIZE_T

   implicit none

   private
   public :: create_attribute, create_dataset, create_corefile

   enum, bind(C)
      enumerator :: I_ONE = 1, I_TWO
   end enum

   integer(kind=SIZE_T), parameter :: default_increment = 1024**2 !< Specifies the increment by which allocated
                                                                  !< memory is to be increased each time more
                                                                  !< memory is required in core file.
   logical, parameter :: default_backing_store = .false.  !< Flag to indicate that entire file contents are flushed to
                                                          !< a file with the same name as the core file

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

!>
!! \brief Creates file in memory using "core" driver
!! \todo check if HDF5 library has been already initialized
!<
   subroutine create_corefile(fname, f_id, incr, bstore)
      use hdf5,          only: HID_T, H5P_FILE_ACCESS_F, H5F_ACC_TRUNC_F, H5P_DEFAULT_F, h5open_f, h5pcreate_f, &
         & h5pset_fapl_core_f, h5fcreate_f, SIZE_T
      implicit none
      character(len=*), intent(in)  :: fname   !< Filename
      integer(HID_T), intent(inout) :: f_id    !< File id
      integer(kind=SIZE_T), intent(in), optional :: incr    !< \copydoc helpers_hdf5::default_increment
      logical, intent(in), optional :: bstore  !< \copydoc helpers_hdf5::default_backing_store

      integer(hid_t) :: faplist_id
      integer(kind=SIZE_T) :: increment
      logical :: backing_store
      integer(kind=4) :: hdferr

      increment = default_increment
      backing_store = default_backing_store
      if (present(incr)) increment = incr
      if (present(bstore)) backing_store = bstore

      ! HDF5 library initialization
      call h5open_f(hdferr)

      ! Create a property list for file access
      call h5pcreate_f(H5P_FILE_ACCESS_F, faplist_id, hdferr)
!      if (hdferr /= 0) call die("[helpers_hdf5:create_corefile] Failed to create property list")

      ! Use magical "core"
      call h5pset_fapl_core_f(faplist_id, increment, backing_store, hdferr)
!      if (hdferr /= 0) call die("[helpers_hdf5:create_corefile] Failed to use core driver")

      ! Create the file with the property list
      call h5fcreate_f(fname, H5F_ACC_TRUNC_F, f_id, hdferr, H5P_DEFAULT_F, faplist_id)
!      if (hdferr /= 0) call die("[helpers_hdf5:create_corefile] Failed to create file in memory")
      return
   end subroutine create_corefile

!> \brief Create 32-bit integer dataset (rank-2 array) in the given place_id.
!
   subroutine create_dataset_int4_dim2(place, dname, ddata)

      use iso_c_binding, only: c_ptr, c_loc
      use hdf5,          only: HID_T, HSIZE_T, H5T_STD_I32LE, &
          &                    h5dcreate_f, h5dclose_f, h5screate_simple_f, h5sclose_f, h5dwrite_f, &
          &                    h5kind_to_type, H5_INTEGER_KIND

      implicit none

      integer(HID_T), intent(in)   :: place !< object id where dataset will be created
      character(len=*), intent(in) :: dname !< name of dataset
      integer(kind=4), dimension(:,:), pointer, intent(in) :: ddata !< data used to create dataset

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

      use iso_c_binding, only: c_ptr, c_loc
      use hdf5,          only: HID_T, HSIZE_T, H5T_STD_I64LE, &
          &                    h5dcreate_f, h5dclose_f, h5screate_simple_f, h5sclose_f, h5dwrite_f, &
          &                    h5kind_to_type, H5_INTEGER_KIND

      implicit none

      integer(HID_T), intent(in)   :: place !< object id where dataset will be created
      character(len=*), intent(in) :: dname !< name of dataset
      integer(kind=8), dimension(:,:), pointer, intent(in) :: ddata !< data used to create dataset

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

      use iso_c_binding, only: c_ptr, c_loc
      use hdf5,          only: HID_T, HSIZE_T, H5T_STD_I64LE, &
          &                    h5dcreate_f, h5dclose_f, h5screate_simple_f, h5sclose_f, h5dwrite_f, &
          &                    h5kind_to_type, H5_INTEGER_KIND

      implicit none

      integer(HID_T), intent(in)   :: place !< object id where dataset will be created
      character(len=*), intent(in) :: dname !< name of dataset
      integer(kind=8), dimension(:), pointer, intent(in) :: ddata !< data used to create dataset

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
