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
!! \brief Some useful, generic functions. Hard to classify, because can be used in completely unrelated modules.
!!
!! \details Note that the functions here are elemental, so one can use them also on arrays to get an array of results.
!<
module func

   implicit none

   private
   public :: ekin, emag, L2norm, sq_sum3

contains

!> \brief Sum of squares of three arguments. Useful for calculating distance, energy, etc.
   elemental real function sq_sum3(x1, x2, x3)
      implicit none
      real, intent(in) :: x1, x2, x3

      sq_sum3 = x1**2 + x2**2 + x3**2

   end function sq_sum3

!> \brief Calculate L2-norm or cartesian distance
   elemental real function L2norm(x1, x2, x3, y1, y2, y3)
      implicit none
      real, intent(in) :: x1, x2, x3, y1, y2, y3

      L2norm = sqrt(sq_sum3(x1-y1, x2-y2, x3-y3))
   end function L2norm

!> \brief Calculate magnetic energy from magnetic field components
   elemental real function emag(bx, by, bz)
      use constants,  only: half
      implicit none
      real, intent(in) :: bx, by, bz

      emag = half*sq_sum3(bx, by, bz)

   end function emag

!> \brief Calculate kinetic energy from momenta and density
   elemental real function ekin(mx, my, mz, dn)
      implicit none
      real, intent(in) :: mx, my, mz, dn

      ekin = emag(mx, my, mz)/dn

   end function ekin

end module func
