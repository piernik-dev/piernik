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
!! \brief Computation of Cosmic Ray sources and pcr gradient pcr
!!
!!
!<
module sourcecosmicrays
! pulled by COSM_RAYS
   implicit none

   private
   public :: src_gpcr
#ifdef COSM_RAYS_SOURCES
   public :: src_crn
#endif /* COSM_RAYS_SOURCES */

contains

!>
!! \brief Computation of Cosmic ray pressure gradient and pcr div v
!<
   subroutine src_gpcr(uu, nn, dx, divv, decr, grad_pcr)

      use domain,         only: dom
      use fluidindex,     only: flind
      use initcosmicrays, only: iarr_crs, gamma_crs, cr_active

      implicit none

      integer(kind=4),                   intent(in)  :: nn                 !< array size
      real, dimension(flind%all,nn),     intent(in)  :: uu                 !< vector of conservative variables
      real, dimension(:), pointer,       intent(in)  :: divv               !< vector of velocity divergence used in cosmic ray advection
      real,                              intent(in)  :: dx                 !< cell length
      real, dimension(nn),               intent(out) :: grad_pcr
      real, dimension(flind%crs%all,nn), intent(out) :: decr
      integer                                        :: icr

      grad_pcr(:) = 0.0
      do icr = 1, flind%crs%all
         ! 1/eff_dim is because we compute the p_cr*dv in every sweep (3 times in 3D, twice in 2D and once in 1D experiments)
         decr(icr,:)      = -1./real(dom%eff_dim)*(gamma_crs(icr)-1.)*uu(iarr_crs(icr),:)*divv(:)
      enddo
         ! Only protons (p+) are dynamically important, we can neglect grad_pcr from heavier nuclei
         ! because of their lower abundancies: n(alpha) ~ 0.1 n(p+), other elements less abundant by orders of magnitude
      grad_pcr(2:nn-1) = grad_pcr(2:nn-1) + cr_active*(gamma_crs(1)-1.)*(uu(iarr_crs(1),1:nn-2)-uu(iarr_crs(icr),3:nn))/(2.*dx)
      grad_pcr(1:2) = 0.0 ; grad_pcr(nn-1:nn) = 0.0

   end subroutine src_gpcr

!==========================================================================================

#ifdef COSM_RAYS_SOURCES
!>
!! \brief Computation of Cosmic ray particles spallation and decay
!! \deprecated BEWARE: several lines in this routine break unit consistency, move it to units.F90 and use scaling
!<
   subroutine src_crn(uu, n, decrn, rk_coeff)

      use cr_data,       only: eCRSP, cr_table, cr_tau, cr_sigma, icr_Be10, icrH, icrL
      use fluids_pub,    only: has_ion, has_neu
      use fluidindex,    only: flind

      implicit none

      integer(kind=4),                  intent(in)  :: n
      real, dimension(flind%all,n),     intent(in)  :: uu
      real,                             intent(in)  :: rk_coeff   !< coeffecient used in RK step, while computing source term
      real, dimension(flind%crn%all,n), intent(out) :: decrn

! locals
      real, dimension(n)      :: dgas
      real, dimension(n)      :: dcr
      real,         parameter :: gamma_lor = 10.0
      real(kind=8), parameter :: speed_of_light = 3e10*1e6*365.*24.*60.*60. !< cm/Myr \deprecated BEWARE: this line breaks unit consistency, move it to units.F90 and use scaling
      real,         parameter :: ndim = 2.0
      real,         parameter :: c_n = speed_of_light / ndim
      real,         parameter :: gn = 1.0 / gamma_lor / ndim

      integer               :: i, j

      dgas = 0.0
      if (has_ion) dgas = dgas + uu(flind%ion%idn,:)
      if (has_neu) dgas = dgas + uu(flind%neu%idn,:)
      dgas = c_n*dgas

      decrn(:,:) = 0.0

      if (eCRSP(icr_Be10)) then
         decrn(cr_table(icr_Be10), :) = decrn(cr_table(icr_Be10), :) - &
            & gn * uu(flind%crn%beg - 1 + cr_table(icr_Be10), :) / cr_tau(cr_table(icr_Be10))
      endif

      do i = lbound(icrH, 1), ubound(icrH, 1)
         associate( Hi => icrH(i) )
            if (eCRSP(Hi)) then
               do j = lbound(icrL, 1), ubound(icrL, 1)
               associate( &
                  Lj => icrL(j), &
                  idx => flind%crn%beg - 1 + cr_table(Hi) &
               )
                  if (eCRSP(Lj)) then
                     dcr = cr_sigma(cr_table(Hi), cr_table(Lj)) * dgas * uu(idx, :)
                     dcr = min(uu(idx,:)/rk_coeff, dcr)  ! Don't decay more elements than available
                     decrn(cr_table(Hi), :) = decrn(cr_table(Hi), :) - dcr
                     decrn(cr_table(Lj), :) = decrn(cr_table(Lj), :) + dcr
                  endif
               end associate
               enddo
         endif
         end associate
      enddo

   end subroutine src_crn
#endif /* COSM_RAYS_SOURCES */
end module sourcecosmicrays
