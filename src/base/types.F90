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
!! \brief Definitions of some types used in several other modules
!!
!! \details The definitions here are either hard to categorize or too simple to go to separate file
!<
module types

   use constants, only: ndims

   implicit none

   private
   public :: value, ema_t

   type :: value
      real                      :: val
      real                      :: assoc
      real,    dimension(ndims) :: coords
      integer, dimension(ndims) :: loc
      integer(kind=4)           :: proc
   end type value

   ! exponential moving average and its variance
   type ema_t
      real          :: avg  !< exponential moving average
      real          :: var  !< exponential moving variance
      real, private :: fac  !< multiplier applied on each accumulation
   contains
      procedure     :: add  !< add next measurement
   end type ema_t

contains

!< \brief Add next measurement to the ema_t. Use factor to (re)set this object.

   subroutine add(this, x, factor)

      implicit none

      class(ema_t),   intent(inout) :: this    !< an object invoking the type-bound procedure
      real,           intent(in)    :: x       !< data
      real, optional, intent(in)    :: factor  !< moving average factor

      real :: d

      if (present(factor)) then  ! reset or initialize
         this%avg = x
         this%var = 0.
         this%fac = factor
      else
         d = x - this%avg
         this%avg = this%avg + this%fac * d
         this%var = (1. - this%fac) * (this%var + this%fac * d**2)
      endif

   end subroutine add

end module types
