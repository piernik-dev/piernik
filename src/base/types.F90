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
!! \brief Definitions of some compound types used in several other modules
!<
module types

   use constants, only: ndims

   implicit none

   private
   public :: axes, value, real_vec_T

   type :: value
      real                      :: val
      real                      :: assoc
      real,    dimension(ndims) :: coords
      integer, dimension(ndims) :: loc
      integer                   :: proc
   end type value

   type :: axes
      real, allocatable, dimension(:) :: x      !< array of x-positions of %grid cells centers
      real, allocatable, dimension(:) :: y      !< array of y-positions of %grid cells centers
      real, allocatable, dimension(:) :: z      !< array of z-positions of %grid cells centers
   end type axes

   type :: real_vec_T
      real, dimension(:), pointer :: r => null()
   contains
      procedure :: associated => real_vec_T_associated
      procedure :: allocate => real_vec_T_allocate
      procedure :: deallocate => real_vec_T_deallocate
   end type real_vec_T

contains

   logical function real_vec_T_associated(this)

      implicit none
      class(real_vec_T), intent(in) :: this

      real_vec_T_associated = associated(this%r)
   end function real_vec_T_associated

   subroutine real_vec_T_allocate(this, n)

      implicit none
      class(real_vec_T), intent(inout) :: this
      integer(kind=4), intent(in) :: n

      allocate(this%r(n))

   end subroutine real_vec_T_allocate

   subroutine real_vec_T_deallocate(this)

      implicit none
      class(real_vec_T), intent(inout) :: this

      deallocate(this%r)
      nullify(this%r)

   end subroutine real_vec_T_deallocate

end module types
