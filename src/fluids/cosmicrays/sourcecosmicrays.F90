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
!! \brief (MH) Computation of advection fluxes of Cosmic Rays
!!
!!
!<

module sourcecosmicrays
! pulled by COSM_RAYS_SOURCES
   implicit none

   private
   public :: src_crn

contains

!==========================================================================================

   subroutine src_crn(uu, n, decrn)

      use fluids_pub,    only: has_ion, has_neu
      use fluidindex,    only: flind
      use cr_data, only: icr_Be10, icr_Be9, icr_C12, sigma_c12_be10, sigma_c12_be9, tau_Be10
      ! icr_Li7, icr_N14, icr_O16, sigma_c12_li7, sigma_n14_li7, sigma_o16_be10, sigma_o16_be9, sigma_o16_li7

      implicit none

      integer(kind=4), intent(in) :: n
      real, dimension(flind%all,n), intent(in)     :: uu
      real, dimension(flind%crn%all,n), intent(out) :: decrn

! locals
      real, allocatable               :: dgas(:)

      real, parameter  :: gamma_lor = 10.
      real(kind=8), parameter :: speed_of_light = 3e10*1e6*365.*24.*60.*60. !> cm/Myr \deprecated BEWARE: this line breaks unit consistency, move it to constants.F90 and use scaling
      real, parameter  :: ndim = 2.

      allocate(dgas(n))

      dgas(:) = 0.0

      if (has_ion) dgas(:) = dgas(:) + uu(flind%ion%idn,:)
      if (has_neu) dgas(:) = dgas(:) + uu(flind%neu%idn,:)

      decrn(:,:) = 0.0

      decrn(icr_C12,:)  = -1./ndim*dgas(:)*speed_of_light*sigma_C12_Be9*uu(flind%crn%beg-1+icr_C12,:) &
           &              -1./ndim*dgas(:)*speed_of_light*sigma_C12_Be10*uu(flind%crn%beg-1+icr_C12,:)

      decrn(icr_Be9,:)  =(  1./ndim*dgas(:)*speed_of_light*sigma_C12_Be9*uu(flind%crn%beg-1+icr_C12,:) )!&
!          &               +1./ndim*dgas(:)*speed_of_light*sigma_O16_Be9*uu(flind%crn%beg-1+icr_O16,:)  )

      decrn(icr_Be10,:) =(  1./ndim*dgas(:)*speed_of_light*sigma_C12_Be10*uu(flind%crn%beg-1+icr_C12,:) &
!          &               +1./ndim*dgas(:)*speed_of_light*sigma_O16_Be10*uu(flind%crn%beg-1+icr_O16,:) &
           &               -1./ndim*uu(flind%crn%beg-1+icr_Be10,:)/gamma_lor/tau_Be10 )

!                        -1./ndim*dgas(:)*speed_of_light*sigma_C12_Li7*uu(flind%crn%beg-1+icr_C12,:)
!    decrn(icr_Li7,:)  =( 1./ndim*dgas(:)*speed_of_light*sigma_C12_Li7*uu(flind%crn%beg-1+icr_C12,:) !&
!                        +1./ndim*dgas(:)*speed_of_light*sigma_N14_Li7*uu(flind%crn%beg-1+icr_N14,:) &
!                        +1./ndim*dgas(:)*speed_of_light*sigma_O16_Li7*uu(flind%crn%beg-1+icr_O16,:) )

!    decrn(icr_N14,:)  = -1./ndim*dgas(:)*speed_of_light*sigma_N14_Li7*uu(flind%crn%beg-1+icr_N14,:)

!    decrn(icr_O16,:)  = -1./ndim*dgas(:)*speed_of_light*sigma_O16_Li7*uu(flind%crn%beg-1+icr_O16,:) &
!                        -1./ndim*dgas(:)*speed_of_light*sigma_O16_Be9*uu(flind%crn%beg-1+icr_O16,:) &
!                        -1./ndim*dgas(:)*speed_of_light*sigma_O16_Be10*uu(flind%crn%beg-1+icr_O16,:)

   end subroutine src_crn

end module sourcecosmicrays
