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
!!\brief This module serves purpose to update single cells using CRESP routines over the whole grid
!<

module cresp_grid
! pulled by CRESP

   implicit none

   public :: cresp_update_grid, cresp_init_grid, cfl_cresp_violation, cresp_clean_grid

   logical :: cfl_cresp_violation
   logical :: allow_loop_leave

contains

!----------------------------------------------------------------------------------------------------

   subroutine cresp_init_grid

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use cg_list_global,   only: all_cg
      use constants,        only: I_ONE, I_TWO
      use cresp_crspectrum, only: cresp_allocate_all, cresp_init_state, p_rch_init
      use cresp_NR_method,  only: cresp_initialize_guess_grids
      use dataio,           only: vars
      use dataio_pub,       only: printinfo, restarted_sim
      use global,           only: repetitive_steps, cflcontrol, disallow_CRnegatives
      use grid_cont,        only: grid_container
      use initcosmicrays,   only: iarr_crspc_n, iarr_crspc_e, ncrb, nspc
      use initcrspectrum,   only: norm_init_spectrum_n, norm_init_spectrum_e, dfpq, check_if_dump_fpq, use_cresp
      use mpisetup,         only: master
      use named_array_list, only: wna

      implicit none

      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
      integer                        :: icr

      if (.not. use_cresp) return

      call cresp_initialize_guess_grids
      call cresp_allocate_all

      call check_if_dump_fpq(vars)

      if (dfpq%any_dump) then
         if (dfpq%dump_f) call all_cg%reg_var(dfpq%f_nam, dim4=ncrb + I_ONE)
         if (dfpq%dump_p) call all_cg%reg_var(dfpq%p_nam, dim4=I_TWO)
         if (dfpq%dump_q) call all_cg%reg_var(dfpq%q_nam, dim4=ncrb)
      endif

      if (.not. restarted_sim) then
         cgl => leaves%first
         do while (associated(cgl))
            cg => cgl%cg

            cg%u(iarr_crspc_n,:,:,:)  = 0.0
            cg%u(iarr_crspc_e,:,:,:)  = 0.0

            if (dfpq%any_dump) then
               if (dfpq%dump_f) cg%w(wna%ind(dfpq%f_nam))%arr = 0.0
               if (dfpq%dump_p) cg%w(wna%ind(dfpq%p_nam))%arr = 0.0
               if (dfpq%dump_q) cg%w(wna%ind(dfpq%q_nam))%arr = 0.0
            endif

            cgl => cgl%nxt
         enddo

      endif

      call p_rch_init               !< sets the right pointer for p_rch function, based on used Taylor expansion coefficient

      call cresp_init_state(norm_init_spectrum_n, norm_init_spectrum_e)   !< initialize spectrum here, f_init should be 1.0

      allow_loop_leave = (disallow_CRnegatives .and. repetitive_steps .and. cflcontrol /= "flex" .and. cflcontrol /= "flexible")

      if (master) call printinfo(" [cresp_grid:cresp_init_grid] CRESP initialized")

   end subroutine cresp_init_grid

   subroutine cresp_update_grid

      use all_boundaries, only: all_fluid_boundaries
      use cg_cost_data,     only: I_MHD
      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: xdim, ydim, zdim, onet
      use cresp_crspectrum, only: cresp_update_cell, q
      use cr_data,          only: eCRSP, ePRIM, ncrsp_prim, ncrsp_sec, cr_table, cr_tau, cr_sigma, icr_Be10, icr_prim, icr_sec, cr_tau, cr_mass, icr_C12, icr_N14, icr_O16, eC12, eO16, eN14, PRIM
      use crhelpers,        only: divv_i
      use dataio_pub,       only: msg, warn
      use func,             only: emag
      use global,           only: dt
      use grid_cont,        only: grid_container
      use initcosmicrays,   only: iarr_crspc2_e, iarr_crspc2_n, nspc, ncrb
      use initcrspectrum,   only: spec_mod_trms, synch_active, adiab_active, cresp, crel, dfpq, fsynchr, u_b_max, use_cresp_evol, bin_old
      use initcrspectrum,   only: cresp_substep, n_substeps_max
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

      integer                        :: i, j, k, nssteps_max, i_prim, i_sec, i_sec_n, i_sec_e
      integer(kind=4)                :: nssteps, i_spc
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer  :: cg
      type(spec_mod_trms)            :: sptab
      type(bin_old), dimension(:), allocatable :: crspc_bins_all
      real                           :: dt_crs_sstep, dt_cresp, dt_doubled
      logical                        :: inactive_cell, cfl_violation_step
      character(len=*), parameter    :: crug_label = "CRESP_upd_grid"
      real, dimension(flind%all)     :: u_cell
      real, dimension(flind%all)     :: usrc_cell
      real                           :: dgas
      real, dimension(ncrb)          :: dcr_e, dcr_n
      real, parameter                :: gamma_lor = 10.0
      real, dimension(ncrb)          :: q_spc
      real, dimension(ncrb,nspc)     :: q_spc_all


      allocate(crspc_bins_all(2*nspc*ncrb))
  !    allocate(dcr_e(ncrb), dcr_n(ncrb))

      if (.not. use_cresp_evol) return

      call ppp_main%start(crug_label)

      cgl => leaves%first
      cfl_cresp_violation = .false.
      cfl_violation_step  = .false.
      inactive_cell       = .false.
      dt_doubled  = 2 * dt       !< used always when cresp_substep is not performed
      dt_cresp    = dt_doubled   !< computed for each cell if cresp_substep, using dt_doubled
      nssteps     = 1
      nssteps_max = 1
      q_spc = 0.
      q_spc_all = 0.
      usrc_cell = 0.0
      do while (associated(cgl))
         cg => cgl%cg
         call cg%costs%start

         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               do i = cg%is, cg%ie

                  u_cell = cg%u(:, i, j, k)
                  dgas = 0.0
                  if (has_ion) dgas = dgas + cg%u(flind%ion%idn, i, j, k) / mp
                  if (has_neu) dgas = dgas + cg%u(flind%neu%idn, i, j, k) / mH
                  sptab%ud = 0.0 ; sptab%ub = 0.0 ; sptab%ucmb = 0.0
                  dgas = dgas * clight

                  do i_spc = 1, nspc

                     cresp%n = cg%u(iarr_crspc2_n(i_spc,:), i, j, k)  !TODO OPTIMIZE ME PLEASE !!
                     cresp%e = cg%u(iarr_crspc2_e(i_spc,:), i, j, k)
                     if (synch_active(i_spc)) sptab%ub = min(emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k)) * fsynchr(i_spc), u_b_max)    !< WARNING assusmes that b is in mGs
                     if (adiab_active(i_spc)) sptab%ud = cg%q(divv_i)%point([i,j,k]) * onet

                     if (cresp_substep) then !< prepare substep timestep for each cell
                        call cresp_timestep_cell(cresp%n, cresp%e, sptab, dt_cresp, i_spc, inactive_cell)
                        call prepare_substep(dt_doubled, dt_cresp, dt_crs_sstep, nssteps)
                        dt_cresp = dt_crs_sstep    !< 2 * dt is equal to nssteps * dt_crs_sstep
                        nssteps_max = max(n_substeps_max, nssteps)
                     endif
#ifdef CRESP_VERBOSED
                     print *, 'Output of cosmic ray electrons module for grid cell with coordinates i,j,k:', i, j, k
#endif /* CRESP_VERBOSED */
                     if (.not. inactive_cell) call cresp_update_cell(dt_cresp, cresp%n, cresp%e, sptab, cfl_violation_step, q_spc, substeps = nssteps)
                    ! if (i_spc == icr_H1 .and. eH1(PRIM))  q_spc_all(:,i_spc) = q_spc
                     if (i_spc == icr_C12 .and. eC12(PRIM)) q_spc_all(:,i_spc) = q_spc
                     if (i_spc == icr_N14 .and. eN14(PRIM)) q_spc_all(:,i_spc) = q_spc
                     if (i_spc == icr_O16 .and. eO16(PRIM)) q_spc_all(:,i_spc) = q_spc
                     !print *, 'q spc : ', q_spc

                     !q_spc(:,i_spc) =
#ifdef DEBUG
                     call cresp_detect_negative_content(cfl_violation_step, [i, j, k])
#endif /* DEBUG */
                     if (cfl_violation_step) then
                        cfl_cresp_violation = cfl_violation_step
                        if (allow_loop_leave) then
                           call cg%costs%stop(I_MHD)
                           call ppp_main%stop(crug_label)
                           return ! nothing to do here!
                        endif
                     endif

                     !print *, 'i_spc : ', i_spc
                     usrc_cell = 0.0

                     !crspc_bins_all(i_spc) = crel
                     !print *, crspc_bins_all(i_spc)

                     !stop
 !                    print *, 'i_spc step ', i_spc, ' crel = ', crspc_bins_all%f
 ! spallation source terms
                     cg%u(iarr_crspc2_n(i_spc,:), i, j, k) = cresp%n
                     cg%u(iarr_crspc2_e(i_spc,:), i, j, k) = cresp%e


!                     call src_cr_spallation_and_decay_cresp(cresp%n, cresp%e)
!                      print *, 'cr_all ', cr_all(iarr_spc)


                     !print *, 'primaries : ', icr_prim, ' secondaries : ', icr_sec
                    ! print *,  ' lbound : ', lbound(icr_prim, 1), ' ubound : ', ubound(icr_prim, 1)


                     !stop
                     !print *, ' ncrsp_prim ', ncrsp_prim
                     !print *, ' ncrsp_sec ', ncrsp_sec
                     !stop


                     !print *, shape(usrc_cell(iarr_crspc2_n(:,:)))
                     !cresp%n = cresp%n + dt_doubled*usrc_cell(iarr_crspc2_n(:,:))
                     !cresp%e = cresp%e + dt_doubled*usrc_cell(iarr_crspc2_e(:,:))

                     !print *, 'crel q : ', crel%q




                     if (dfpq%any_dump) then
                     !print *, 'dump'
                     !stop
                        if (dfpq%dump_f) cg%w(wna%ind(dfpq%f_nam))%arr(:, i, j, k) = crel%f
                        if (dfpq%dump_p) cg%w(wna%ind(dfpq%p_nam))%arr(:, i, j, k) = crel%p(crel%i_cut)
                        if (dfpq%dump_q) cg%w(wna%ind(dfpq%q_nam))%arr(:, i, j, k) = crel%q
                     endif
                  enddo

                  !print *, 'q spc all : ', q_spc_all

                  do i_prim = 1, ncrsp_prim
                     associate( cr_prim => cr_table(icr_prim(i_prim)) )
                        !print *, 'cr_prim : ', cr_prim
                        !print *, 'i_prim : ', i_prim
                        if (eCRSP(icr_prim(i_prim))) then

                           do i_sec = 1, ncrsp_sec
                           associate( cr_sec => cr_table(icr_sec(i_sec)) )


                              if (eCRSP(icr_sec(i_sec))) then
                                 dcr_n = cr_sigma(cr_prim, cr_sec) * dgas * (cr_mass(icr_prim(i_prim))/cr_mass(icr_sec(i_sec)))**(3-q_spc_all(:,icr_prim(i_prim)))*u_cell(iarr_crspc2_n(cr_prim,:))
                                 dcr_n = min(u_cell(iarr_crspc2_n(cr_prim,:)), dcr_n)  ! Don't decay more elements than available
                                 !print *, 'i_sec : ', i_sec
                                 !print *, 'u_cell n sec :', u_cell(iarr_crspc2_n(cr_sec,:))
                                 !print *, 'u_cell e sec :', u_cell(iarr_crspc2_e(cr_sec,:))
                                 !print *, 'dgas : ', dgas
                                 !print *, 'u_cell : ', u_cell(iarr_crspc2_n(cr_prim,:))

                                 usrc_cell(iarr_crspc2_n(cr_prim,:)) = usrc_cell(iarr_crspc2_n(cr_prim,:)) - dcr_n
                                 usrc_cell(iarr_crspc2_n(cr_sec,:)) = usrc_cell(iarr_crspc2_n(cr_sec,:)) + dcr_n

                                 dcr_e = cr_sigma(cr_prim, cr_sec) * dgas *(cr_mass(icr_prim(i_prim))/cr_mass(icr_sec(i_sec)))**(4-q_spc_all(:,icr_prim(i_prim)))* u_cell(iarr_crspc2_e(cr_prim,:))
                                 dcr_e = min(u_cell(iarr_crspc2_e(cr_prim,:)), dcr_e)  ! Don't decay more elements than available

                                 usrc_cell(iarr_crspc2_e(cr_prim,:)) = usrc_cell(iarr_crspc2_e(cr_prim,:)) - dcr_e
                                 usrc_cell(iarr_crspc2_e(cr_sec,:)) = usrc_cell(iarr_crspc2_e(cr_sec,:)) + dcr_e
                                 !if (cr_sigma(cr_prim, cr_sec) .gt. 0.0) then
                                    !print *, i_prim, i_sec
                                    !print *, 'sigma : ', cr_sigma(cr_prim, cr_sec)
                                    !print *, 'dcr_n : ', dcr_n
                                    !print *, 'dcr_e : ', dcr_e
                                 !endif

                           !if (cr_sigma(cr_prim, cr_sec) .gt. 0.) then

                           !   print *, i_prim, icr_prim(i_prim)
                           !   print *, i_sec, icr_sec(i_sec)
                           !   print *, ' cr_prim : ', cr_prim
                           !   print *, 'cr_sec : ', cr_sec
                           !   print *, 'cr_table prim: ' ,cr_table(icr_prim(i_prim))
                           !   print *, 'cr_table sec : ',cr_table(icr_sec(i_sec))
                           !   print *, 'sigma : ', cr_sigma(cr_prim, cr_sec)
                           !   print *, 'prim mass : ', cr_mass(icr_prim(i_prim))
                           !   print *, 'sec mass : ', cr_mass(icr_sec(i_sec))
                           !print *, 'dcr_n : ', dcr_n
                           !print *, 'dcr_e : ', dcr_e
                           !endif
                              endif
                           end associate
                           enddo

                        endif
                     end associate
                  !print *, icr_prim(i_prim)
                  enddo

                  i_sec = cr_table(icr_Be10) !; i_sec_n = iarr_crspc2_n(i_sec,:) ; i_sec_e = iarr_crspc2_e(i_sec,:)
                  if (eCRSP(icr_Be10)) then
                     !print *, 'ok'
                     !print *, 'Be10 life time : ', cr_tau(i_sec), ' MyR'
                     !usrc(:, j) = usrc(:, j) - gn * uu(:, j) / cr_tau(i)
                     !print *, 'n before : ', usrc_cell(iarr_crspc2_e(i_sec,:))
                     !print *, 'e before : ', usrc_cell(iarr_crspc2_n(i_sec,:))
                     usrc_cell(iarr_crspc2_n(i_sec,:)) = usrc_cell(iarr_crspc2_n(i_sec,:)) -u_cell(iarr_crspc2_n(i_sec,:))/(gamma_lor*cr_tau(i_sec))
                     !print *, 'usrc_cell (n) : ', usrc_cell(iarr_crspc2_n(i_sec,:))
                     usrc_cell(iarr_crspc2_e(i_sec,:)) = usrc_cell(iarr_crspc2_e(i_sec,:)) -u_cell(iarr_crspc2_e(i_sec,:))/(gamma_lor*cr_tau(i_sec))
                     !print *, 'u_cell n :', u_cell(iarr_crspc2_n(i_sec,:))
                     !print *, 'u_cell e :', u_cell(iarr_crspc2_e(i_sec,:))
                  endif


                  do i_spc = 1, nspc
                     cg%u(iarr_crspc2_n(i_spc,:), i, j, k) = cg%u(iarr_crspc2_n(i_spc,:), i, j, k) + dt_doubled*usrc_cell(iarr_crspc2_n(i_spc,:))
                     cg%u(iarr_crspc2_e(i_spc,:), i, j, k) = cg%u(iarr_crspc2_e(i_spc,:), i, j, k) + dt_doubled*usrc_cell(iarr_crspc2_e(i_spc,:))
                    ! print *, i_spc
                    ! print *, 'usrc_cell n : ', usrc_cell(iarr_crspc2_n(i_spc,:))
                    ! print *, 'usrc_cell e : ', usrc_cell(iarr_crspc2_e(i_spc,:))
                    ! print *, 'dt_2 : ', dt_doubled
                     !if (maxval(usrc_cell(iarr_crspc2_n(i_spc,:))) .gt. 0.0) stop
                     !if (maxval(usrc_cell(iarr_crspc2_e(i_spc,:))) .gt. 0.0) stop
                  enddo

               enddo
            enddo
         enddo


         call cg%costs%stop(I_MHD)
         cgl=>cgl%nxt
      enddo
      if (nssteps_max > n_substeps_max) then
         write (msg,"(A50,I5, A27, I5)") "[cresp_grid:cresp_update_grid] max(n_substeps) = ", nssteps_max, " in this cg exceeds limit ", n_substeps_max
         call warn(msg)
      endif

      call all_fluid_boundaries

      call ppp_main%stop(crug_label)

   end subroutine cresp_update_grid


   !==========================================================================================

!>
!! \brief Computation of Cosmic ray particles spallation and decay
!! \warning multiplying by cr_sigma may cause underflow errors in some unit sets
!<
!   subroutine src_cr_spallation_and_decay_cresp(uu, )
!
!      use cr_data,        only: eCRSP, cr_table, cr_tau, cr_sigma, icr_Be10, icr_prim, icr_sec
!      use domain,         only: dom
!      use fluids_pub,     only: has_ion, has_neu
!      use fluidindex,     only: flind
!      use initcosmicrays, only: iarr_crspc2_e, iarr_crspc2_n, iarr_crn
!      use initcrspectrum, only: cresp
!      use units,          only: clight, mH, mp
!
!      implicit none
!
!      !integer(kind=4),               intent(in)  :: n
!      real, dimension(1, flind%all), intent(in)  :: uu
!      !real,                          intent(in)  :: rk_coeff   !< coefficient used in RK step, while computing source term
!      real, dimension(1, flind%all), intent(out) :: usrc       !< u array update component for sources
!
!! locals
!      real, dimension(1)                         :: dgas, dcr
!      real, parameter                            :: gamma_lor = 10.0 !< \todo should taken from cosmic ray species settings
!      real                                       :: gn
!      integer                                    :: i, j
!      !type(bin_old), dimension(:), allocatable   :: bins_all
!
!      gn = 1.0 / gamma_lor
!      dgas = 0.0
!      if (has_ion) dgas = dgas + uu(:, flind%ion%idn) / mp
!      if (has_neu) dgas = dgas + uu(:, flind%neu%idn) / mH
!      dgas = dgas * clight / dom%eff_dim
!
!      usrc(:,:) = 0.0
!
!      !i = cr_table(icr_Be10) ; j = iarr_crn(i)
!      !if (eCRSP(icr_Be10)) usrc(:, j) = usrc(:, j) - gn * uu(:, j) / cr_tau(i)
!
!      print *, 'primaries : ', icr_prim, ' secondaries : ', icr_sec
!      print *,  ' lbound : ', lbound(icr_prim, 1), ' ubound : ', ubound(icr_prim, 1)
!         do i = lbound(icr_prim, 1), ubound(icr_prim, 1)
!         associate( cr_prim => cr_table(icr_prim(i)) )
!            if (eCRSP(icr_prim(i))) then
!               do j = lbound(icr_sec, 1), ubound(icr_sec, 1)
!               associate( cr_sec => cr_table(icr_sec(j)) )
!                  if (eCRSP(icr_sec(j))) then
!                     dcr = cr_sigma(cr_prim, cr_sec) * dgas * uu(:, iarr_crn(cr_prim))
!                     dcr = min(uu(:, iarr_crn(cr_prim)), dcr)  ! Don't decay more elements than available
!                     usrc(:, iarr_crn(cr_prim)) = usrc(:, iarr_crn(cr_prim)) - dcr
!                     usrc(:, iarr_crn(cr_sec)) = usrc(:, iarr_crn(cr_sec)) + dcr
!                  endif
!               end associate
!               enddo
!         endif
!         end associate
!      enddo
!
!   end subroutine src_cr_spallation_and_decay_cresp





   subroutine prepare_substep(dt_simulation, dt_process_short, dt_substep, n_substeps)

      use initcrspectrum, only: n_substeps_max
#ifdef CRESP_VERBOSED
      use dataio_pub,     only: msg, warn
      use mpisetup,       only: master
#endif /* CRESP_VERBOSED */

      implicit none

      real,            intent(in)  :: dt_simulation, dt_process_short
      real,            intent(out) :: dt_substep
      integer(kind=4), intent(out) :: n_substeps

      n_substeps = ceiling(dt_simulation / dt_process_short, kind=4)  ! ceiling to assure resulting dt_substep .le. dt_process_short

#ifdef CRESP_VERBOSED
      if (n_substeps > n_substeps_max .and. master) then
         write (msg,"(A42,I5, A14, I5)") "[cresp_grid:prepare_substep] n_substeps = ", n_substeps, " exceeds limit ", n_substeps_max
         call warn(msg)
      endif
#endif /* CRESP_VERBOSED */

      dt_substep = dt_simulation / n_substeps

   end subroutine prepare_substep

!----------------------------------------------------------------------------------------------------

   subroutine cresp_clean_grid

      use cg_cost_data,     only: I_MHD
      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use cresp_crspectrum, only: detect_clean_spectrum
      use grid_cont,        only: grid_container
      use initcosmicrays,   only: iarr_crspc_e, iarr_crspc_n
      use initcrspectrum,   only: cresp, nullify_empty_bins
      use ppp,              only: ppp_main

      implicit none

      integer                        :: i, j, k
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
      character(len=*), parameter    :: crcg_label = "CRESP_clean_grid"

      if (.not.nullify_empty_bins) return

      call ppp_main%start(crcg_label)

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         call cg%costs%start

         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               do i = cg%is, cg%ie
                  cresp%n = cg%u(iarr_crspc_n, i, j, k)
                  cresp%e = cg%U(iarr_crspc_e, i, j, k)

                  call detect_clean_spectrum

                  cg%u(iarr_crspc_n, i, j, k) = cresp%n
                  cg%u(iarr_crspc_e, i, j, k) = cresp%e
               enddo
            enddo
         enddo

         call cg%costs%stop(I_MHD)
         cgl=>cgl%nxt
      enddo

      call ppp_main%stop(crcg_label)

   end subroutine cresp_clean_grid

end module cresp_grid
