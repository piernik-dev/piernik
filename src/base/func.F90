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
!! Run:
!!   grep -r "\*\*2[^\*]*\*\*2[^\*]*\*\*2" src {,../}problems/
!! in the trunk to see a list of potential places, where use of these function may simplify the code.
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

!> \brief Calculate mean value of 3D gaussian profile in n-times resampled cell
   elemental function resample_gauss(x, y, z, dx, dy, dz, sx, sy, sz, n) result(val)
      implicit none

      real, intent(in) :: x, y, z      !! coordinates of the cell center in the reference frame of gaussian, i.e. (x - x_0) ...
      real, intent(in) :: dx, dy, dz   !! cell sizes
      real, intent(in) :: sx, sy, sz   !! sigmas
      integer, intent(in) :: n
      real :: val

      real :: xi, yj, zk
      real :: odx, ody, odz

      integer :: i, j ,k

      odx = dx / n
      ody = dy / n
      odz = dz / n

      val = 0.0
      do k = 1, n
         zk = (z - 0.5*dz) + (k-0.5)*odz
         do j = 1, n
            yj = (y - 0.5*dy) + (j-0.5)*ody
            do i = 1, n
               xi = (x - 0.5*dx) + (i-0.5)*odx
               val = val + exp( -0.5*sq_sum3((xi/sx)**2, (yj/sy)**2, (zk/sz)**2) )
            enddo
         enddo
      enddo

      val = val / n**3

   end function resample_gauss

end module func
