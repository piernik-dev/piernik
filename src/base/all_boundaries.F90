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
!! /brief Global update of boundary routines
!!
!! /details moved to a separate file due to concerns about circular dependencies
!<
module all_boundaries

   implicit none

   private
   public :: all_bnd

contains

!>
!! Subroutine calling all type boundaries after initialization of new run or restart reading
!<

   subroutine all_bnd

      use fluidboundaries, only: all_fluid_boundaries
#ifdef MAGNETIC
      use magboundaries,   only: all_mag_boundaries
#endif /* MAGNETIC */

      implicit none

!      if (all(cg%bnd(:,:) /= BND_USER)) then
! \todo make sure that all_fluid_boundaries and all_mag_boundaries can handle BND_USER boundaries right now, or do the boundaries later
      call all_fluid_boundaries
#ifdef MAGNETIC
      call all_mag_boundaries
#endif /* MAGNETIC */
!     endif

   end subroutine all_bnd

end module all_boundaries
