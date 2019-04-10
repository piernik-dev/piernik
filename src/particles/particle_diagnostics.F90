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
!! \brief Diagnostics of particle module
!!
!! This module contains particle diagnostics related things that are needed by I/O and debug routines.
!<

module particle_diagnostics
! pulled by NBODY

   implicit none

   private
   public :: max_pvel_1d, max_pacc_3d

contains

   subroutine max_pvel_1d(max_v)

      use constants,      only: ndims, xdim, zdim
      use particle_types, only: pset

      implicit none

      real, dimension(ndims), intent(out) :: max_v
      integer(kind=4)                     :: cdim

      do cdim = xdim, zdim
         max_v(cdim)  = maxval(abs(pset%p(:)%vel(cdim)))
      enddo

   end subroutine max_pvel_1d

   subroutine max_pacc_3d(max_acc)

      use constants,      only: zero
      use particle_types, only: pset

      implicit none

      real, intent(out) :: max_acc
      real              :: acc2
      integer           :: i, n_part

      max_acc  = zero
      n_part = size(pset%p, dim=1)

      do i = 1, n_part
         acc2 = sum(pset%p(i)%acc(:)**2)
         max_acc = max(acc2, max_acc)
      enddo
      max_acc = sqrt(max_acc)

   end subroutine max_pacc_3d

end module particle_diagnostics
