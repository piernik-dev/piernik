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
   public :: src_gpcr, src_cr_spallation_and_decay, cr_spallation_sources

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
         !print *, 'values :', grad_pcr_cresp
      endif
#endif /* CRESP */
      usrc(:, iarr_all_en(flind%ion%pos)) = vx(:, flind%ion%pos) * grad_pcr
      !print *, 'grad_pcr_cresp', grad_pcr_cresp
#ifdef CRESP
      usrc(:, iarr_all_mx(flind%ion%pos)) = grad_pcr_cresp
#ifdef ISO
      return
#endif /* ISO */
      !print *, 'grad pcr cresp : ', grad_pcr_cresp
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
      use dataio_pub,     only: die
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

      if (dom%eff_dim == 0) call die("[sourcecosmicrays:src_cr_spallation_and_decay] dom%eff_dim == 0 is not supported yet")

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

   subroutine cr_spallation_sources(i,j,k,u_cell,dt_doubled, q_spc_all)

      use all_boundaries,   only: all_fluid_boundaries
      use constants,        only: xdim, ydim, zdim, onet, one, zero
      use cresp_crspectrum, only: cresp_update_cell, q
      use cr_data,          only: eCRSP, ePRIM, ncrsp_prim, ncrsp_sec, cr_table, cr_tau, cr_sigma, icr_Be9, icr_prim, icr_sec, cr_tau, cr_mass, icr_C12, icr_N14, icr_O16, eC12, eO16, eN14, PRIM, cr_mass
      use dataio_pub,       only: msg, warn
      use func,             only: emag
      use global,           only: dt
      use grid_cont,        only: grid_container
      use initcosmicrays,   only: iarr_crspc2_e, iarr_crspc2_n, nspc, ncrb
      use initcrspectrum,   only: spec_mod_trms, synch_active, adiab_active, cresp, crel, dfpq, fsynchr, u_b_max, use_cresp_evol, bin_old, eps
      use initcrspectrum,   only: cresp_substep, n_substeps_max, p_fix
      use named_array_list, only: wna
      use ppp,              only: ppp_main
      use timestep_cresp,   only: cresp_timestep_cell
      use fluidindex,       only: flind
      use fluids_pub,       only: has_ion, has_neu
      use units,            only: clight, mH, mp
#ifdef DEBUG
      use cresp_crspectrum, only: cresp_detect_negative_content
#endif /* DEBUG */

      implicit none

      integer                        :: i,j,k,nssteps_max, i_prim, i_sec, i_sec_n, i_sec_e
      integer(kind=4)                :: nssteps, i_spc, i_bin
      type(spec_mod_trms)            :: sptab
      type(bin_old), dimension(:), allocatable :: crspc_bins_all
      real                           :: dt_crs_sstep, dt_cresp, dt_doubled
      logical                        :: inactive_cell, cfl_violation_step
      character(len=*), parameter    :: crug_label = "CRESP_upd_grid"
      real, dimension(flind%all)     :: u_cell
      real, dimension(flind%all)     :: usrc_cell
      real, dimension(1:ncrb)        :: dcr_n, dcr_e, Q_ratio_1, Q_ratio_2, S_ratio_1, S_ratio_2
      !real, dimension(1:ncrb)        :: Q_ratio_1, Q_ratio_2, S_ratio_1, S_ratio_2
      real                           :: dgas
      !real                           :: dt_doubled
      real, parameter                :: gamma_lor = 10.0
      real, dimension(ncrb,nspc)     :: q_spc_all

      inactive_cell       = .false.

      dgas = 0.0
      if (has_ion) dgas = dgas + u_cell(flind%ion%idn) / mp
      if (has_neu) dgas = dgas + u_cell(flind%neu%idn) / mH
      sptab%ud = 0.0 ; sptab%ub = 0.0 ; sptab%ucmb = 0.0

      usrc_cell = 0.0

      Q_ratio_1 = 0.0
      S_ratio_1 = 0.0
      Q_ratio_2 = 0.0
      S_ratio_2 = 0.0

      do i_prim = 1, ncrsp_prim
         associate( cr_prim => cr_table(icr_prim(i_prim)) )
         !print *, 'index prim : ', icr_prim(i_prim)
         !print *, 'i prim : ', i_prim

            if (eCRSP(icr_prim(i_prim))) then

               do i_sec = 1, ncrsp_sec
               associate( cr_sec => cr_table(icr_sec(i_sec)) )

                  if (eCRSP(icr_sec(i_sec))) then

                     do i_bin = 1, ncrb

                        if (p_fix(i_bin) > zero) then

                           if (p_fix(i_bin+1) > zero) then

                              Q_ratio_1(i_bin) = Q_ratio(cr_mass(icr_prim(i_prim)), cr_mass(icr_sec(i_sec)),q_spc_all(i_bin,icr_prim(i_prim)),p_fix(i_bin),p_fix(i_bin+1))
                              S_ratio_1(i_bin) = S_ratio(cr_mass(icr_prim(i_prim)), cr_mass(icr_sec(i_sec)),q_spc_all(i_bin,icr_prim(i_prim)),p_fix(i_bin),p_fix(i_bin+1)) !TODO print values of ratio in the last bin

                              Q_ratio_2(i_bin) = one - Q_ratio_1(i_bin)
                              !   Q_ratio_2(i_bin-1) = Q_ratio_2(i_bin)
                              !   !Q_ratio_2(i_bin) = 0.0

                              !
                              S_ratio_2(i_bin) = one - S_ratio_1(i_bin)
                              !   S_ratio_2(i_bin-1) = S_ratio_2(i_bin)
                              !   !S_ratio_2(i_bin) = 0.0


                                    !dcr_e(i_bin)/(dcr_n(i_bin)*p_fix(i_bin))
                              !if(Q_ratio_2(i_bin)*dcr_n(i_bin)*cr_mass(icr_sec(i_sec))/cr_mass(icr_prim(i_prim))*p_fix(i_bin-1) > zero) print *, 'Secondary Ratio ',i_bin,' : ', S_ratio_2(i_bin)*dcr_e(i_bin)/(Q_ratio_2(i_bin)*dcr_n(i_bin)*cr_mass(icr_sec(i_sec))/cr_mass(icr_prim(i_prim))*p_fix(i_bin-1))

                           endif

                        endif

                     enddo



                     !Q_ratio_2 = 1.0 - Q_ratio_1
                     !S_ratio_2 = 1.0 - S_ratio_1

                     dcr_n(:) = cr_sigma(cr_prim, cr_sec) * dgas * clight * u_cell(iarr_crspc2_n(cr_prim,:))
                     dcr_e(:) = cr_sigma(cr_prim, cr_sec) * dgas * clight * u_cell(iarr_crspc2_e(cr_prim,:))

                     dcr_n(:) = min(dcr_n,u_cell(iarr_crspc2_n(cr_prim,:)))
                     dcr_e(:) = min(dcr_e,u_cell(iarr_crspc2_e(cr_prim,:))) ! Don't decay more element than available

                     usrc_cell(iarr_crspc2_n(cr_prim,:)) = usrc_cell(iarr_crspc2_n(cr_prim,:)) - dcr_n(:)

                     usrc_cell(iarr_crspc2_e(cr_prim,:)) = usrc_cell(iarr_crspc2_e(cr_prim,:)) - dcr_e(:)

                        !if (i_bin == ncrb) then

                     !usrc_cell(iarr_crspc2_n(cr_sec,:)) = usrc_cell(iarr_crspc2_n(cr_sec,:)) + dcr_n(:)

                     !usrc_cell(iarr_crspc2_e(cr_sec,:)) = usrc_cell(iarr_crspc2_e(cr_sec,:)) + dcr_e(:)

                     do i_bin = 1, ncrb

                        if (i_bin == ncrb) then

                           usrc_cell(iarr_crspc2_n(cr_sec,i_bin)) = usrc_cell(iarr_crspc2_n(cr_sec,i_bin)) + Q_ratio_2(i_bin) * dcr_n(i_bin)
                         !
                           usrc_cell(iarr_crspc2_e(cr_sec,i_bin)) = usrc_cell(iarr_crspc2_e(cr_sec,i_bin)) + S_ratio_2(i_bin) * dcr_e(i_bin)*cr_mass(icr_sec(i_sec))/cr_mass(icr_prim(i_prim))

                        else

                           usrc_cell(iarr_crspc2_n(cr_sec,i_bin)) = usrc_cell(iarr_crspc2_n(cr_sec,i_bin)) + Q_ratio_2(i_bin) * dcr_n(i_bin) + Q_ratio_1(i_bin+1)*dcr_n(i_bin+1)
                         !
                           usrc_cell(iarr_crspc2_e(cr_sec,i_bin)) = usrc_cell(iarr_crspc2_e(cr_sec,i_bin)) + (S_ratio_2(i_bin) * dcr_e(i_bin) + S_ratio_1(i_bin+1)*dcr_e(i_bin+1))*cr_mass(icr_sec(i_sec))/cr_mass(icr_prim(i_prim))



                        endif

                        !if(icr_prim(i_prim)==icr_C12) then
                        !      if(icr_sec(i_sec)==icr_Be9) then
                        !         if (dcr_n(i_bin)*p_fix(i_bin) > zero) then
                        !
                        !            print *, 'Primary Ratio n',i_bin,' : ', Q_ratio_1(i_bin), Q_ratio_2(i_bin), 'de/dnp : ', dcr_e(i_bin)/(dcr_n(i_bin)*p_fix(i_bin))
                        !
                        !            print *, 'Primary Ratio e',i_bin,' : ', S_ratio_1(i_bin), S_ratio_2(i_bin), 'e_c/n_c p : ', u_cell(iarr_crspc2_e(cr_prim,i_bin))/(u_cell(iarr_crspc2_n(cr_prim,i_bin)) * p_fix(i_bin))
                        !
                        !         endif
                        !
                        !      endif
                        !endif







                        !usrc_cell(iarr_crspc2_n(cr_sec,i_bin)) = usrc_cell(iarr_crspc2_n(cr_sec,i_bin)) + Q_ratio_2(i_bin) * dcr_n(i_bin) + Q_ratio_1(i_bin)*dcr_n(i_bin)
                        !
                        !usrc_cell(iarr_crspc2_e(cr_sec,i_bin)) = usrc_cell(iarr_crspc2_e(cr_sec,i_bin)) + S_ratio_1(i_bin) * dcr_e(i_bin) + S_ratio_2(i_bin) * dcr_e(i_bin)

                        !endif

                     enddo

                     !if(icr_prim(i_prim)==icr_C12) then
                     !      if(icr_sec(i_sec)==icr_Be9) then
                     !         stop
                     !      endif
                     !endif

                     !if (icr_prim(i_prim)==icr_C12) then
                     !   if (icr_sec(i_sec)==icr_Be9) stop
                     !endif

                  endif
               end associate
               enddo

            endif
         end associate
      enddo

      do i_spc = 1, nspc

         u_cell(iarr_crspc2_n(i_spc,:)) = u_cell(iarr_crspc2_n(i_spc,:)) + dt_doubled*usrc_cell(iarr_crspc2_n(i_spc,:))
         u_cell(iarr_crspc2_e(i_spc,:)) = u_cell(iarr_crspc2_e(i_spc,:)) + dt_doubled*usrc_cell(iarr_crspc2_e(i_spc,:))

      enddo

   end subroutine cr_spallation_sources

   function Q_ratio(A_prim, A_sec, q_l, p_L, p_R)

      use constants,      only: zero, one, three
      use initcrspectrum, only: eps

      implicit none

      real :: A_prim, A_sec, q_l, p_L, p_R
      real :: Q_ratio

      Q_ratio = zero

      if(abs(q_l - three) > eps) then
         Q_ratio = ((A_prim/A_sec)**(three-q_l)-one)/((p_R/p_L)**(three-q_l)-one)
      else
         Q_ratio = log(A_prim/A_sec)/log(p_R/p_L)
      endif

   end function Q_ratio

   function S_ratio(A_prim, A_sec, q_l, p_L, p_R)

      use constants,      only: zero, one, four
      use initcrspectrum, only: eps

      implicit none

      real :: A_prim, A_sec, q_l, p_L, p_R
      real :: S_ratio

      S_ratio = zero

      if(abs(q_l - four) > eps) then
         S_ratio = ((A_prim/A_sec)**(four-q_l)-one)/((p_R/p_L)**(four-q_l)-one)
      else
         S_ratio = log(A_prim/A_sec)/log(p_R/p_L)
      endif

   end function S_ratio

end module sourcecosmicrays
