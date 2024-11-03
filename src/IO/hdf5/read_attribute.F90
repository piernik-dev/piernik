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

!> \brief Module that contains HDF5 I/O routines for reading attributes.

module read_attr
! pulled by HDF5

   implicit none

   private
   public :: read_attribute

!> \brief Read an attribute (1D array) from the given group
   interface read_attribute
      module procedure read_int_attribute
      module procedure read_real_attribute
   end interface

contains

!> \brief Read an integer attribute (1D array) from the given group

   subroutine read_int_attribute(g_id, name, int_array)

      use constants, only: I_ONE
      use hdf5,      only: HID_T, HSIZE_T, H5T_NATIVE_INTEGER, h5aopen_f, h5aclose_f, h5aread_f

      implicit none

      integer(HID_T),                intent(in)    :: g_id       !< group id where to create the attribute
      character(len=*),              intent(in)    :: name       !< name
      integer(kind=4), dimension(:), intent(inout) :: int_array  !< the data

      integer(HID_T)                     :: attr_id
      integer(kind=4)                    :: error    !< error perhaps should be of type integer(HID_T)
      integer(HSIZE_T), dimension(I_ONE) :: dims

      !> \todo Implement size checks
      call h5aopen_f(g_id, name, attr_id, error)
      dims = size(int_array, kind=HSIZE_T)
      call h5aread_f(attr_id, H5T_NATIVE_INTEGER, int_array, dims, error)
      call h5aclose_f(attr_id, error)

   end subroutine read_int_attribute

!> \brief Read a real attribute (1D array) from the given group

   subroutine read_real_attribute(g_id, name, real_array)

      use constants, only: I_ONE
      use hdf5,      only: HID_T, HSIZE_T, H5T_NATIVE_DOUBLE, h5aopen_f, h5aclose_f, h5aread_f

      implicit none

      integer(HID_T),     intent(in)    :: g_id        !< group id where to create the attribute
      character(len=*),   intent(in)    :: name        !< name
      real, dimension(:), intent(inout) :: real_array  !< the data

      integer(HID_T)                     :: attr_id
      integer(kind=4)                    :: error      !< error perhaps should be of type integer(HID_T)
      integer(HSIZE_T), dimension(I_ONE) :: dims

      !> \todo Implement size checks
      call h5aopen_f(g_id, name, attr_id, error)
      dims = size(real_array, kind=HSIZE_T)
      call h5aread_f(attr_id, H5T_NATIVE_DOUBLE, real_array, dims, error)
      call h5aclose_f(attr_id, error)

   end subroutine read_real_attribute

end module read_attr
