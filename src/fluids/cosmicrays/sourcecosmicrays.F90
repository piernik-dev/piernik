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
!! \brief Computation of Cosmic Ray sources and pcr gradient pcr
!!
!!
!<
module sourcecosmicrays
! pulled by COSM_RAYS
   implicit none

   private
   public :: src_gpcr, src_cr_spallation_and_decay

contains


!>
!! \brief Computation of Cosmic ray pressure gradient and pcr div v
!<
   subroutine src_gpcr(uu, nn, usrc, sweep, i1, i2, cg, vx)

      use crhelpers,        only: set_div_v1d
      use domain,           only: dom
      use fluidindex,       only: flind, iarr_all_mx, iarr_all_en
      use grid_cont,        only: grid_container
      use initcosmicrays,   only: cr_active, gamma_cr_1, gpcr_ess_noncresp, iarr_crn
#ifdef CRESP
      use cr_data,          only: cr_gpess, cr_spectral, cr_table, icr_H1, iarr_spc
      use cresp_crspectrum, only: src_gpcresp
      use initcosmicrays,   only: iarr_crspc2_e
#endif /* CRESP */

      implicit none

      integer(kind=4),                intent(in)  :: nn                 !< array size
      real, dimension(nn, flind%all), intent(in)  :: uu                 !< vector of conservative variables
      integer(kind=4),                intent(in)  :: sweep              !< direction (x, y or z) we are doing calculations for
      integer,                        intent(in)  :: i1                 !< coordinate of sweep in the 1st remaining direction
      integer,                        intent(in)  :: i2                 !< coordinate of sweep in the 2nd remaining direction
      type(grid_container), pointer,  intent(in)  :: cg                 !< current grid piece
      real, dimension(:,:), pointer,  intent(in)  :: vx
      real, dimension(nn, flind%all), intent(out) :: usrc               !< u array update component for sources
!locals
      logical                                     :: full_dim
      integer                                     :: icr, jcr
      real, dimension(:), pointer                 :: divv               !< vector of velocity divergence used in cosmic ray advection
      real, dimension(nn)                         :: grad_pcr
#ifdef CRESP
      real, dimension(nn)                         :: grad_pcr_cresp
#endif /* CRESP */

      full_dim = nn > 1
      !print *, 'break 1'
      !print *, full_dim
      !print *, ' '
      !print *, flind%crspc%all

      usrc = 0.0
      if (.not. full_dim .or. flind%crspc%all < 1) return
      !print *, 'break 2'

      call set_div_v1d(divv, sweep, i1, i2, cg)
      do icr = 1, flind%crn%all
         ! 1/eff_dim is because we compute the p_cr*dv in every sweep (3 times in 3D, twice in 2D and once in 1D experiments)
         usrc(:, iarr_crn(icr)) = -1. / real(dom%eff_dim) * gamma_cr_1 * uu(:, iarr_crn(icr)) * divv(:)
      enddo

      !< gpcr_ess_noncresp includes species only for non-CRESP treatment and where cr_gpcr_ess = .true. for unnamed species and set in CR_SPECIES namelist for named species
      grad_pcr(:) = 0.0
      do icr = 1, size(gpcr_ess_noncresp)
         jcr = gpcr_ess_noncresp(icr)
         grad_pcr(2:nn-1) = grad_pcr(2:nn-1) + cr_active * gamma_cr_1 * (uu(1:nn-2, iarr_crn(jcr)) - uu(3:nn, iarr_crn(jcr))) / (2. * cg%dl(sweep))
      enddo
      grad_pcr(1:2) = 0.0 ; grad_pcr(nn-1:nn) = 0.0

      usrc(:, iarr_all_mx(flind%ion%pos)) = grad_pcr
      !print *, 'break 3'
#ifdef CRESP
      if (cr_spectral(cr_table(icr_H1)) .or. cr_gpess(cr_table(icr_H1))) then !< Primarily treat protons as the source of GPCR !TODO expand me for other CR spectral species (optional)
         !print *, 'uu array : ', uu(:,iarr_crspc2_e(iarr_spc(icr_H1), 4))
         call src_gpcresp(uu(:,iarr_crspc2_e(iarr_spc(icr_H1), :)), nn, cg%dl(sweep), grad_pcr_cresp, iarr_spc(icr_H1))         !< cg%dl(sweep) = dx, contribution due to pressure acted upon spectral components in CRESP via div_v
         !print *, 'break 4'
      endif
#endif /* CRESP */
#ifdef ISO
      return
#endif /* ISO */
      usrc(:, iarr_all_en(flind%ion%pos)) = vx(:, flind%ion%pos) * grad_pcr
      !print *, 'grad_pcr_cresp', grad_pcr_cresp
#ifdef CRESP
      usrc(:, iarr_all_mx(flind%ion%pos)) = grad_pcr_cresp
      usrc(:, iarr_all_en(flind%ion%pos)) = vx(:, flind%ion%pos) * grad_pcr_cresp !< BEWARE - check it
#endif /* CRESP */

   end subroutine src_gpcr

!==========================================================================================

!>
!! \brief Computation of Cosmic ray particles spallation and decay
!! \warning multiplying by cr_sigma may cause underflow errors in some unit sets
!<
   subroutine src_cr_spallation_and_decay(uu, n, usrc, rk_coeff)

      use cr_data,        only: eCRSP, cr_table, cr_tau, cr_sigma, icr_Be10, icr_prim, icr_sec
      use domain,         only: dom
      use fluids_pub,     only: has_ion, has_neu
      use fluidindex,     only: flind
      use initcosmicrays, only: iarr_crn
      use units,          only: clight, mH, mp

      implicit none

      integer(kind=4),               intent(in)  :: n
      real, dimension(n, flind%all), intent(in)  :: uu
      real,                          intent(in)  :: rk_coeff   !< coefficient used in RK step, while computing source term
      real, dimension(n, flind%all), intent(out) :: usrc       !< u array update component for sources

! locals
      real, dimension(n)                         :: dgas, dcr
      real, parameter                            :: gamma_lor = 10.0 !< \todo should taken from cosmic ray species settings
      real                                       :: gn
      integer                                    :: i, j

      gn = 1.0 / dom%eff_dim / gamma_lor
      dgas = 0.0
      if (has_ion) dgas = dgas + uu(:, flind%ion%idn) / mp
      if (has_neu) dgas = dgas + uu(:, flind%neu%idn) / mH
      dgas = dgas * clight / dom%eff_dim

      usrc(:,:) = 0.0

      i = cr_table(icr_Be10) ; j = iarr_crn(i)
      if (eCRSP(icr_Be10)) usrc(:, j) = usrc(:, j) - gn * uu(:, j) / cr_tau(i)

      do i = lbound(icr_prim, 1), ubound(icr_prim, 1)
         associate( cr_prim => cr_table(icr_prim(i)) )
            if (eCRSP(icr_prim(i))) then
               do j = lbound(icr_sec, 1), ubound(icr_sec, 1)
               associate( cr_sec => cr_table(icr_sec(j)) )
                  if (eCRSP(icr_sec(j))) then
                     dcr = cr_sigma(cr_prim, cr_sec) * dgas * uu(:, iarr_crn(cr_prim))
                     dcr = min(uu(:, iarr_crn(cr_prim))/rk_coeff, dcr)  ! Don't decay more elements than available
                     usrc(:, iarr_crn(cr_prim)) = usrc(:, iarr_crn(cr_prim)) - dcr
                     usrc(:, iarr_crn(cr_sec)) = usrc(:, iarr_crn(cr_sec)) + dcr
                  endif
               end associate
               enddo
         endif
         end associate
      enddo

   end subroutine src_cr_spallation_and_decay

end module sourcecosmicrays
