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
!! \brief Module that allows easily set and get attributes in dump and restart files.
!<

module set_get_attributes

! pulled by HDF5

   use constants,  only: I_ONE
   use dataio_pub, only: msglen

   implicit none

   private
   public :: set_attr, get_attr

   integer(kind=4) :: error
   character(len=I_ONE) :: root_path = "/"
   character(len=msglen) :: path_   !< a buffer for attribute path manipulation

!> \brief Set common types of attributes, autodetect data types and sizes
   interface set_attr
      module procedure set_attr_R
      module procedure set_attr_I
      module procedure set_attr_C
   end interface set_attr

!> \brief Get common types of attributes, autodetect data types and sizes
   interface get_attr
      module procedure get_attr_R
      module procedure get_attr_I
      module procedure get_attr_C
   end interface get_attr

!> \ToDo Add 64-bit integer attributes (can't use H5LT in fortran for this, need to follow {set,get}_attr_C)

contains

!> \brief Set a real attribute (by default in "/" if not specified otherwise), autodetect the size.

   subroutine set_attr_R(file_id, name, array, path)

      use dataio_pub, only: die
      use hdf5,       only: HID_T, SIZE_T
      use h5lt,       only: h5ltset_attribute_double_f

      implicit none

      integer(HID_T),             intent(in) :: file_id  !< File identifier
      character(len=*),           intent(in) :: name     !< Real attribute name
      real, dimension(:),         intent(in) :: array    !< Array of values (should contain at least one element)
      character(len=*), optional, intent(in) :: path     !< Path to attribute ("/" if not specified)

      if (size(array) < 1) call die("[set_get_attributes:set_attr_R] empty input array")

      path_ = checkpath(root_path)
      if (present(path)) path_ = checkpath(path)

      call h5ltset_attribute_double_f(file_id, path_, name, array, size(array, kind=SIZE_T), error)

   end subroutine set_attr_R

!> \brief Get a real attribute (by default in "/" if not specified otherwise), autodetect the size.

   subroutine get_attr_R(file_id, name, array, path)

      use cmp_1D_mpi, only: compare_array1D
      use dataio_pub, only: msg, warn
      use hdf5,       only: HID_T, HSIZE_T
      use h5lt,       only: h5ltget_attribute_double_f

      implicit none

      integer(HID_T),                  intent(in)  :: file_id  !< File identifier
      character(len=*),                intent(in)  :: name     !< Real attribute name
      real, allocatable, dimension(:), intent(out) :: array    !< Array of returned values (empty array marks an error)
      character(len=*), optional,      intent(in)  :: path     !< Path to attribute ("/" if not specified)

      integer(kind=4) :: rank
      integer(HSIZE_T), dimension(:), allocatable :: dims

      path_ = checkpath(root_path)
      if (present(path)) path_ = checkpath(path)

      call find_rank_dims(file_id, path_, name, rank, dims)
      if (rank /= 1 .or. .not. allocated(dims)) then
         allocate(array(0))
         write(msg, '(5a,i3,a,7i10)')"[set_get_attributes:get_attr_R] Failed to read attribute '", trim(name), " at '", trim(path_), "' rank = ", rank, " dims() = ", dims
         call warn(msg)
         return
      endif

      allocate(array(dims(1)))
      call h5ltget_attribute_double_f(file_id, path_, name, array, error)

      deallocate(dims)

      call compare_array1D(array)

   end subroutine get_attr_R

!---------------------------------------------------------------------------

!> \brief Set an int32 attribute (by default in "/" if not specified otherwise), autodetect the size.

! semi-spaghetti, but can't yet figure out how to really exploit polymorphism here without using class(*) and risking ICE and other kinds of compiler faults

   subroutine set_attr_I(file_id, name, array, path)

      use dataio_pub, only: die
      use hdf5,       only: HID_T, SIZE_T
      use h5lt,       only: h5ltset_attribute_int_f

      implicit none

      integer(HID_T),                intent(in) :: file_id  !< File identifier
      character(len=*),              intent(in) :: name     !< Integer attribute name
      integer(kind=4), dimension(:), intent(in) :: array    !< Array of values (should contain at least one element)
      character(len=*), optional,    intent(in) :: path     !< Path to attribute ("/" if not specified)

      if (size(array) < 1) call die("[set_get_attributes:set_attr_I] empty input array")

      path_ = checkpath(root_path)
      if (present(path)) path_ = checkpath(path)

      call h5ltset_attribute_int_f(file_id, path_, name, array, size(array, kind=SIZE_T), error)

   end subroutine set_attr_I

!> \brief Get an int32 attribute (by default in "/" if not specified otherwise), autodetect the size.

   subroutine get_attr_I(file_id, name, array, path)

      use cmp_1D_mpi, only: compare_array1D
      use dataio_pub, only: msg, warn
      use hdf5,       only: HID_T, HSIZE_T
      use h5lt,       only: h5ltget_attribute_int_f

      implicit none

      integer(HID_T),                             intent(in)  :: file_id  !< File identifier
      character(len=*),                           intent(in)  :: name     !< Integer attribute name
      integer(kind=4), allocatable, dimension(:), intent(out) :: array    !< Array of returned values (empty array marks an error)
      character(len=*), optional,                 intent(in)  :: path     !< Path to attribute ("/" if not specified)

      integer(kind=4) :: rank
      integer(HSIZE_T), dimension(:), allocatable :: dims

      path_ = checkpath(root_path)
      if (present(path)) path_ = checkpath(path)

      call find_rank_dims(file_id, path_, name, rank, dims)
      if (rank /= 1 .or. .not. allocated(dims)) then
         allocate(array(0))
         write(msg, '(5a,i3,a,7i10)')"[set_get_attributes:get_attr_I] Failed to read attribute '", trim(name), " at '", trim(path_), "' rank = ", rank, " dims() = ", dims
         call warn(msg)
         return
      endif

      allocate(array(dims(1)))
      call h5ltget_attribute_int_f(file_id, path_, name, array, error)

      deallocate(dims)

      call compare_array1D(array)

   end subroutine get_attr_I

!---------------------------------------------------------------------------
!> \brief Set a character attribute (by default in "/" if not specified otherwise), autodetect the size.

! semi-spaghetti, but can't yet figure out how to really exploit polymorphism here without using class(*) and risking ICE and other kinds of compiler faults

   subroutine set_attr_C(file_id, name, array, path)

      use constants,  only: I_ONE
      use dataio_pub, only: die
      use hdf5,       only: HID_T, HSIZE_T, SIZE_T, H5T_FORTRAN_S1, &
           &                h5acreate_f, h5aclose_f, h5awrite_f, &
           &                h5gopen_f, h5gclose_f, &
           &                h5screate_simple_f, h5sclose_f, &
           &                h5tcopy_f, h5tset_size_f, h5tclose_f

      implicit none

      integer(HID_T),                 intent(in) :: file_id  !< File identifier
      character(len=*),               intent(in) :: name     !< Character attribute name
      character(len=*), dimension(:), intent(in) :: array    !< Array of values (should contain at least one element)
      character(len=*), optional,     intent(in) :: path     !< Path to attribute ("/" if not specified)

      integer(kind=4), parameter         :: arank = I_ONE
      integer(HID_T)                     :: aspace_id, attr_id, g_id, type_id
      integer(kind=4)                    :: error
      integer(HSIZE_T), dimension(arank) :: dims

      if (size(array) < 1) call die("[set_get_attributes:set_attr_C] empty input array")

      path_ = checkpath(root_path)
      if (present(path)) path_ = checkpath(path)

      dims = size(array)
      call h5tcopy_f(H5T_FORTRAN_S1, type_id, error)
      call h5tset_size_f(type_id, int(len(array), kind=SIZE_T), error)
      call h5gopen_f(file_id, path_, g_id, error)
      call h5screate_simple_f(arank, dims, aspace_id, error)
      call h5acreate_f(g_id, name, type_id, aspace_id, attr_id, error)
      call h5awrite_f(attr_id, type_id, array, dims, error)
      call h5aclose_f(attr_id, error)
      call h5sclose_f(aspace_id, error)
      call h5gclose_f(g_id, error)
      call h5tclose_f(type_id, error)

   end subroutine set_attr_C

!> \brief Get an int32 attribute (by default in "/" if not specified otherwise), autodetect the size.

   subroutine get_attr_C(file_id, name, array, path)

      use constants,  only: cbuff_len
      use cmp_1D_mpi, only: compare_array1D
      use dataio_pub, only: msg, warn
      use hdf5,       only: HID_T, HSIZE_T, SIZE_T, H5T_FORTRAN_S1, H5S_SCALAR_F, &
           &                h5aopen_f, h5aclose_f, h5aread_f, &
           &                h5gopen_f, h5gclose_f, &
           &                h5tcopy_f, h5tset_size_f, h5tclose_f
      use h5lt,       only: h5ltget_attribute_string_f

      implicit none

      integer(HID_T),                                      intent(in)  :: file_id  !< File identifier
      character(len=*),                                    intent(in)  :: name     !< Character attribute name
      character(len=cbuff_len), allocatable, dimension(:), intent(out) :: array    !< Array of returned values (empty array marks an error)
      character(len=*), optional,                          intent(in)  :: path     !< Path to attribute ("/" if not specified)

      integer(kind=4)                             :: rank
      integer(HSIZE_T), dimension(:), allocatable :: dims
      integer(HID_T)                              :: attr_id, g_id, type_id
      integer(kind=4)                             :: error
      integer                                     :: i

      path_ = checkpath(root_path)
      if (present(path)) path_ = checkpath(path)

      call find_rank_dims(file_id, path_, name, rank, dims)
      if ((rank /= 1 .and. rank /= H5S_SCALAR_F) .or. .not. allocated(dims)) then
         allocate(array(0))
         write(msg, '(5a,i3,a,7i10)')"[set_get_attributes:get_attr_C] Failed to read attribute '", trim(name), " at '", trim(path_), "' rank = ", rank, " dims() = ", dims
         call warn(msg)
         return
      endif

      if (rank == H5S_SCALAR_F) then  ! keep compatibility with old restarts
         allocate(array(1))
         call h5ltget_attribute_string_f(file_id, path_, name, array(1), error)
      else
         allocate(array(dims(1)))
         call h5tcopy_f(H5T_FORTRAN_S1, type_id, error)
         call h5tset_size_f(type_id, int(len(array), kind=SIZE_T), error)
         call h5gopen_f(file_id, path_, g_id, error)
         call h5aopen_f(g_id, name, attr_id, error)
         call h5aread_f(attr_id, type_id, array, dims, error)
         call h5aclose_f(attr_id, error)
         call h5gclose_f(g_id, error)
         call h5tclose_f(type_id, error)
      endif

      deallocate(dims)

      do i = lbound(array, 1), ubound(array, 1)
         call compare_array1D(trim(array(i)))
      enddo

   end subroutine get_attr_C

!---------------------------------------------------------------------------

!> \brief Query attribute

   subroutine find_rank_dims(file_id, path, name, rank, dims)

      use constants,  only: INVALID, cbuff_len
      use dataio_pub, only: msg, warn
      use hdf5,       only: HID_T, HSIZE_T, SIZE_T, H5T_STRING_F, H5S_SCALAR_F
      use h5lt,       only: h5ltget_attribute_ndims_f, h5ltget_attribute_info_f

      implicit none

      integer(HID_T),                              intent(in)  :: file_id  !< File identifier
      character(len=*),                            intent(in)  :: path     !< Path to attribute
      character(len=*),                            intent(in)  :: name     !< Real attribute name
      integer(kind=4),                             intent(out) :: rank     !< Rank of attribute
      integer(HSIZE_T), dimension(:), allocatable, intent(out) :: dims     !< Dimensions of attribute

      integer(kind=4) :: tclass
      integer(SIZE_T) :: tsize

      rank = INVALID
      allocate(dims(0))

      call h5ltget_attribute_ndims_f(file_id, path, name, rank, error)
      if (error /= 0) then
         write(msg, '(4a)')"[set_get_attributes:find_rank_dims] cannot read rank of attribute '", trim(name), "' at ", trim(path)
         call warn(msg)  ! no need for if (master) here
         return
      endif
      if (rank > 1) then
         write(msg, '(a,i2,a)')"[set_get_attributes:find_rank_dims] rank-", rank, " attributes aren't supported, only rank-1 are allowed"
         call warn(msg)  ! no need for if (master) here
         return
      endif

      deallocate(dims)
      if (rank  == H5S_SCALAR_F) then
         allocate(dims(1))  ! promote scalar to an array for simplicity
      else
         allocate(dims(rank))
      endif
      call h5ltget_attribute_info_f(file_id, path, name, dims, tclass, tsize, error)

      if (error /= 0) then
         write(msg, '(4a)')"[set_get_attributes:find_rank_dims] cannot read info of attribute '", trim(name), "' at ", trim(path)
         call warn(msg)  ! no need for if (master) here
         deallocate(dims)
         return
      endif

      if (rank  == H5S_SCALAR_F) dims = 1   ! promote scalar to an array for simplicity

      if (tclass == H5T_STRING_F .and. tsize > cbuff_len) then
         write(msg, '(2(a,i3),a)')"[set_get_attributes:find_rank_dims] too long attribute: len(element): ", tsize, " > ", cbuff_len, ". Truncating."
         call warn(msg)
      endif

   end subroutine find_rank_dims

!> \brief Check if provided path is not too long

   function checkpath(path) result(path_)

      use dataio_pub, only: msg, die

      implicit none

      character(len=*), intent(in) :: path

      character(len=msglen) :: path_

      if (len(path) <= len(path_)) then
         path_ = path
      else
         write(msg, '(2(a,i4))')"[set_get_attributes:setpath] Attribute path too long (", len(path), " > ", len(path_)
         call die(msg)

      endif

   end function checkpath

end module set_get_attributes
