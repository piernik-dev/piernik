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
#define RNG 2:n-1

!>
!! \brief Computation of advection %fluxes of Cosmic Rays
!!
!!
!<

module fluxcosmicrays
! pulled by COSM_RAYS
   implicit none

   private
   public :: flux_crs

contains

   subroutine flux_crs(fluxc, vion, uuc, n)

      use fluidindex, only: flind

      implicit none

      integer(kind=4), intent(in)                  :: n
      real, dimension(n),   intent(in)             :: vion
      real, dimension(:,:), intent(in),    pointer :: uuc
      real, dimension(:,:), intent(inout), pointer :: fluxc
      integer :: i

      fluxc   = 0.0

      do i = 1, flind%crs%all
         fluxc(RNG, i) = uuc(RNG, i)*vion(RNG)
      enddo

   end subroutine flux_crs

end module fluxcosmicrays
