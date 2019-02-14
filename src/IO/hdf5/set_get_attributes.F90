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
!      module procedure set_attr_I
!      module procedure set_attr_C
   end interface set_attr

!> \brief Get common types of attributes, autodetect data types and sizes
   interface get_attr
      module procedure get_attr_R
!      module procedure get_attr_I
!      module procedure get_attr_C
   end interface get_attr

contains

!> \brief Set a real attribute (by default in "/" if not specified otherwise), autodetect the size.

   subroutine set_attr_R(file_id, name, arr_R, path)

      use dataio_pub, only: die
      use hdf5,       only: HID_T, SIZE_T
      use h5lt,       only: h5ltset_attribute_double_f

      implicit none

      integer(HID_T),             intent(in) :: file_id  !< File identifier
      character(len=*),           intent(in) :: name     !< Real attribute name
      real, dimension(:),         intent(in) :: arr_R    !< Array of values (should contain at least one element)
      character(len=*), optional, intent(in) :: path     !< Path to attribute ("/" if not specified)

      if (size(arr_R) < 1) call die("[set_get_attributes:set_attr_R] empty input array")

      path_ = checkpath(root_path)
      if (present(path)) path_ = checkpath(path)

      call h5ltset_attribute_double_f(file_id, path_, name, arr_R, size(arr_R, kind=SIZE_T), error)

   end subroutine set_attr_R

!> \brief Get a real attribute (by default in "/" if not specified otherwise), autodetect the size.

   subroutine get_attr_R(file_id, name, arr_R, path)

      use dataio_pub, only: msg, warn
      use hdf5,       only: HID_T, HSIZE_T, SIZE_T
      use h5lt,       only: h5ltget_attribute_ndims_f, h5ltget_attribute_info_f, h5ltget_attribute_double_f

      implicit none

      integer(HID_T),                  intent(in)  :: file_id  !< File identifier
      character(len=*),                intent(in)  :: name     !< Real attribute name
      real, allocatable, dimension(:), intent(out) :: arr_R    !< Array of returned values (empty array marks an error)
      character(len=*), optional,      intent(in)  :: path     !< Path to attribute ("/" if not specified)

      integer(kind=4) :: rank, tclass
      integer(SIZE_T) :: tsize
      integer(HSIZE_T), dimension(:), allocatable :: dims

      path_ = checkpath(root_path)
      if (present(path)) path_ = checkpath(path)

      call h5ltget_attribute_ndims_f(file_id, path_, name, rank, error)
      if (error /= 0) then
         write(msg, '(4a)')"[set_get_attributes:get_attr_R] cannot read rank of attribute '", trim(name), "' at ", trim(path_)
         call warn(msg)  ! no need for if (master) here
         allocate(arr_R(0))
         return
      endif
      if (rank > 1) then
         write(msg, '(a,i2,a)')"[set_get_attributes:get_attr_R] rank-", rank, " attributes aren't supported, only rank-1 are allowed"
         call warn(msg)  ! no need for if (master) here
         allocate(arr_R(0))
         return
      endif

      allocate(dims(rank))
      call h5ltget_attribute_info_f(file_id, path_, name, dims, tclass, tsize, error)
      if (error /= 0) then
         write(msg, '(4a)')"[set_get_attributes:get_attr_R] cannot read info of attribute '", trim(name), "' at ", trim(path_)
         call warn(msg)  ! no need for if (master) here
         allocate(arr_R(0))
         deallocate(dims)
        return
      endif

      allocate(arr_R(dims(1)))
      call h5ltget_attribute_double_f(file_id, path_, name, arr_R, error)

      deallocate(dims)

   end subroutine get_attr_R

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
