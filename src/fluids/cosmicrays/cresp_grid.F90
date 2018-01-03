module cresp_grid
! pulled by COSM_RAY_ELECTRONS

! This module contains routines necessary to initialize, compute timestep for cre bins and to update spectrum in the whole domain
! as the crspectrum module operates on a single grid cell.
      use initcosmicrays, only: iarr_cre_e, iarr_cre_n
      use initcrspectrum, only: ncre
      use global,         only: dt, t

      private
      public        dt_cre, cresp_update_grid, cresp_init_grid, grid_cresp_timestep

      real(kind=8)                    :: dt_cre
      integer(kind=4), save           :: i_up_max_prev, i_mid=0, j_mid=0, k_mid=0
contains

 subroutine cresp_update_grid
  use cg_leaves,      only: leaves
  use cg_list,        only: cg_list_element
  use constants,      only: xdim, ydim, zdim, zero
  use grid_cont,      only: grid_container
  use cresp_crspectrum, only:cresp_update_cell, printer
  use initcrspectrum, only: spec_mod_trms, virtual_e, virtual_n, eps, prevent_neg_e
  use named_array_list, only: qna
  use crhelpers,      only: divv_n
  use func,           only: emag, ekin, operator(.equals.), operator(.notequals.)
  implicit none
    integer                         :: i, j, k
    type(cg_list_element),  pointer :: cgl
    type(grid_container),   pointer :: cg
    real(kind=8), dimension(1:ncre) :: n_cell, e_cell
    type(spec_mod_trms)  :: sptab
        i = 0; j = 0;  k = 0
        cgl => leaves%first
        n_cell = zero
        e_cell = zero
        do while (associated(cgl))
            cg => cgl%cg
            do k = cg%ks, cg%ke
                do j = cg%js, cg%je
                    do i = cg%is, cg%ie
                        if ( prevent_neg_e ) then ! e = eps where it dropped below zero due to diffusion algorithm - TEMP workaround
                            where (cg%u(iarr_cre_n, i, j, k) .lt. zero)
                                cg%u(iarr_cre_n, i, j, k) = eps
                            endwhere
                            where (cg%u(iarr_cre_e, i, j, k) .lt. zero)
                                cg%u(iarr_cre_e, i, j, k) = eps
                            endwhere
                        endif
                        sptab%ud = 0.0 ; sptab%ub = 0.0 ; sptab%ucmb = 0.0
                        n_cell    = cg%u(iarr_cre_n, i, j, k)
                        e_cell    = cg%u(iarr_cre_e, i, j, k)
                        sptab%ub = emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))
                        sptab%ud = cg%q(qna%ind(divv_n))%point([i,j,k])
#ifdef VERBOSE
                        print *, 'Output of cosmic ray electrons module for grid cell with coordinates i,j,k:', i, j, k
#endif /* VERBOSE */
                        call cresp_update_cell(2*dt, n_cell, e_cell, sptab, virtual_n(1:2,i,j,k), virtual_e(1:2,i,j,k))
                        cg%u(iarr_cre_n, i, j, k) = n_cell
                        cg%u(iarr_cre_e, i, j, k) = e_cell
                        if (i.eq.i_mid.and.j.eq.j_mid.and.k.eq.k_mid) then ! diagnostic:
                            call printer(t)
                        endif
                    enddo
                enddo
            enddo
            cgl=>cgl%nxt
        enddo
  end subroutine cresp_update_grid
!----------------------------------------------------------------------------------------------------
  subroutine cresp_init_grid
   use cg_leaves,      only: leaves
   use cg_list,        only: cg_list_element
   use constants,      only: I_ONE, I_FOUR, fpi, LO, HI, xdim, ydim, zdim, zero
   use grid_cont,      only: grid_container
   use initcrspectrum, only: ncre, f_init, p_up_init, p_lo_init, q_init, cre_eff, spec_mod_trms, initial_condition, bump_amp, &
                            virtual_e, virtual_n, e_small, e_small_approx_p_lo, e_small_approx_p_up
   use initcosmicrays, only: iarr_crn
   use cresp_crspectrum, only: cresp_init_state, cresp_allocate_all, printer, e_threshold_lo, e_threshold_up, &
                        & fail_count_interpol, fail_count_no_sol, fail_count_NR_2dim, fail_count_comp_q, second_fail
!    use units,          only: clight
   use cresp_variables, only: clight ! temp, TODO
   implicit none
    integer                         :: i, j, k
    type(cg_list_element),  pointer :: cgl
    type(grid_container),   pointer :: cg
    real(kind=8)                             :: max_amp_cr, f_amplitude, e_tot
    real(kind=8), dimension(1:ncre) :: n_cell, e_cell
    type(spec_mod_trms)  :: sptab
    logical, save :: first_run = .true.
      if (first_run .eqv. .true.) then
        call cresp_allocate_all
        
        fail_count_interpol = 0
        fail_count_no_sol   = 0
        fail_count_NR_2dim  = 0
        second_fail         = 0
        fail_count_comp_q   = 0

        e_threshold_lo = e_small * e_small_approx_p_lo
        e_threshold_up = e_small * e_small_approx_p_up
        
        i = 0; j = 0; k = 0
        cgl => leaves%first
        do while (associated(cgl))
            cg => cgl%cg
            cg%u(iarr_cre_e,:,:,:) = zero
            cg%u(iarr_cre_n,:,:,:) = zero
            max_amp_cr = maxval(cg%u(iarr_crn(1),:,:,:))

            i_mid = ((cg%lhn(xdim,LO)+cg%lhn(xdim,HI))/2)
            j_mid = ((cg%lhn(ydim,LO)+cg%lhn(ydim,HI))/2)
            k_mid = ((cg%lhn(zdim,LO)+cg%lhn(zdim,HI))/2)

            if (.not. allocated(virtual_e)) allocate(virtual_e(1:2, cg%lhn(xdim,LO):cg%lhn(xdim,HI), &
                cg%lhn(ydim,LO):cg%lhn(ydim,HI), cg%lhn(zdim,LO):cg%lhn(zdim,HI)))
            if (.not. allocated(virtual_n)) allocate(virtual_n(1:2, cg%lhn(xdim,LO):cg%lhn(xdim,HI), &
                cg%lhn(ydim,LO):cg%lhn(ydim,HI), cg%lhn(zdim,LO):cg%lhn(zdim,HI)))
            virtual_e = zero
            virtual_n = zero

            n_cell = zero
            e_cell = zero
            do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
                do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
                    do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
            ! Every initial condition should be normalized before initializing Cosmic Ray Electron SPectrum module
                        if (cg%u(iarr_crn(1),i,j,k)*cre_eff .gt. e_small) then ! early phase - fill cells only when total passed energy is greater than e_small
                            if ( initial_condition == "powl") then
            ! amplitude and distribution of electron energy density is inherited after those of nucleons, see crspectrum.pdf, eq. 29
                                e_tot = cg%u(iarr_crn(1),i,j,k) * cre_eff
                                f_amplitude = (e_tot / (fpi * clight * p_lo_init ** I_FOUR) ) * (I_FOUR - q_init) / &
                                                                 ((p_up_init/p_lo_init)**(I_FOUR - q_init) - I_ONE  )
                            endif
                            if (initial_condition == "bump") then
            ! for gaussian distribution & inheritance of spatial energy/number density after nucleons
                                f_amplitude = cre_eff * cg%u(iarr_crn(1),i,j,k)
!                               f_amplitude = bump_amp ! * clight
                            endif
                            call cresp_init_state(n_cell, e_cell, f_amplitude, sptab)
#ifdef VERBOSE
                            print *, 'Output of cosmic ray electrons module for grid cell with coordinates i,j,k:', i, j, k
#endif /* VERBOSE */
                            cg%u(iarr_cre_n, i, j, k) = n_cell
                            cg%u(iarr_cre_e, i, j, k) = e_cell
                        endif ! if total(cre_eff*e) less than e_small - nothing done, cell remains uninitialized
                        if (i.eq.i_mid.and.j.eq.j_mid.and.k.eq.k_mid) then ! diagnostics
                            call printer(t)
                        endif
                    enddo
                enddo
            enddo
            cgl=>cgl%nxt
            enddo
            i_up_max_prev = 0
            first_run = .false. ! FIXME uncommenting results inf SIGFPE for some reason; whole subroutine is called twice.
        endif
  end subroutine cresp_init_grid
!----------------------------------------------------------------------------------------------------
  subroutine grid_cresp_timestep
   use cg_leaves,        only: leaves
   use cg_list,          only: cg_list_element
   use crhelpers,        only: divv_n
   use func,             only: emag !, operator(.equals.), operator(.notequals.)
   use grid_cont,        only: grid_container
   use constants,        only: xdim, ydim, zdim
   use named_array_list, only: qna
   use constants,        only: one, half
   use initcrspectrum,   only: spec_mod_trms, cfl_cre
   use initcosmicrays,   only: K_cre_paral, K_cre_perp
   use timestep_cresp,   only: cresp_timestep
   implicit none
    integer(kind=4)                 :: i, j, k, i_up_max, i_up_max_tmp
    type(grid_container), pointer   :: cg
    type(cg_list_element), pointer  :: cgl
    real(kind=8)                    :: dt_cre_tmp, K_cre_max_sum
    real(kind=8),save               :: dt_cre_K
    type(spec_mod_trms)             :: sptab
        i_up_max     = 1
        i_up_max_tmp = 1

        dt_cre = huge(one)
        dt_cre_tmp = huge(one)

        cgl => leaves%first
        do while (associated(cgl))
            cg => cgl%cg
            do k = cg%ks, cg%ke
                do j = cg%js, cg%je
                    do i = cg%is, cg%ie
                        sptab%ud = 0.0 ; sptab%ub = 0.0 ; sptab%ucmb = 0.0
                        sptab%ub = emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))
                        sptab%ud = cg%q(qna%ind(divv_n))%point([i,j,k])
                        call cresp_timestep(dt_cre_tmp, sptab, cg%u(iarr_cre_n, i, j, k), cg%u(iarr_cre_e, i, j, k), i_up_max_tmp) ! gives dt_cre for the whole domain, but is unefficient
                        dt_cre = min(dt_cre, dt_cre_tmp)
                        i_up_max = max(i_up_max, i_up_max_tmp)
                    enddo
                enddo
            enddo
            cgl=>cgl%nxt
        enddo

        if ( i_up_max_prev .ne. i_up_max ) then ! dt_cre_K saved, computed again only if in the whole domain highest i_up changes.
            i_up_max_prev = i_up_max
            K_cre_max_sum = K_cre_paral(i_up_max) + K_cre_perp(i_up_max) ! assumes the same K for energy and number density
            if ( K_cre_max_sum <= 0) then                                ! K_cre dependent on momentum - maximal for highest bin number
                dt_cre_K = huge(one)
            else
                dt_cre_K = cfl_cre * half / K_cre_max_sum
                if (cg%dxmn < sqrt(huge(one))/dt_cre_K) then
                    dt_cre_K = dt_cre_K * cg%dxmn**2
                endif
            endif
        endif
        dt_cre = min(dt_cre, dt_cre_K)
        dt_cre = half * dt_cre ! dt comes in to cresp_crspectrum with factor * 2
        
  end subroutine grid_cresp_timestep
!----------------------------------------------------------------------------------------------------
 subroutine append_dissipative_terms(i,j,k) ! To be fixed
  use initcrspectrum,   only: spec_mod_trms
  use named_array_list, only: qna
  use crhelpers,        only: divv_n
  use func,             only: emag
  use grid_cont,        only: grid_container
  use cg_list,          only: cg_list_element
  use constants,        only: xdim, ydim, zdim
  use cg_leaves,        only: leaves
  implicit none
    type(spec_mod_trms)  :: sptab
    type(grid_container), pointer :: cg
    type(cg_list_element), pointer:: cgl
    integer :: i,j,k
!Below - magnetic energy density and velocity divergence values are passed to sptab
        cgl => leaves%first
        do while (associated(cgl))
            cg => cgl%cg
            sptab%ub = emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))
            sptab%ud = cg%q(qna%ind(divv_n))%point([i,j,k])
            sptab%ucmb = 0.0 ! Not included yet
            cgl =>cgl%nxt
        enddo
 end subroutine append_dissipative_terms

end module cresp_grid
