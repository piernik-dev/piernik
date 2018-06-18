module cresp_grid
! pulled by COSM_RAY_ELECTRONS

! This module contains routines necessary to initialize, compute timestep for cre bins and to update spectrum in the whole domain
! as the crspectrum module operates on a single grid cell.

   use global,          only: dt, t
   use initcosmicrays,  only: iarr_cre_e, iarr_cre_n
   use initcrspectrum,  only: ncre

   private
   public        dt_cre, cresp_update_grid, cresp_init_grid, grid_cresp_timestep, cfl_cresp_violation, cresp_clean_grid

   real(kind=8)                    :: dt_cre
   real(kind=8)                    :: bb_to_ub
   logical                         :: cfl_cresp_violation
   integer(kind=4), save           :: i_up_max_prev
! CRESP names
   character(len=*), parameter :: vn_n = "vn" !< helping array for CRESP number density
   character(len=*), parameter :: ve_n = "ve" !< helping array for CRESP energy density

   contains


   subroutine cresp_update_grid

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: xdim, ydim, zdim, onet
      use cresp_crspectrum, only: cresp_update_cell
      use crhelpers,        only: divv_n
      use func,             only: emag, ekin, operator(.equals.), operator(.notequals.)
      use grid_cont,        only: grid_container
      use initcrspectrum,   only: spec_mod_trms, synch_active, adiab_active, cresp
      use named_array,      only: p4
      use named_array_list, only: qna, wna

      implicit none

      real(kind=8), dimension(:), pointer :: virtual_e => null(), virtual_n => null()
      integer                             :: i, j, k
      type(cg_list_element), pointer      :: cgl
      type(grid_container), pointer       :: cg
      type(spec_mod_trms)                 :: sptab

      i = 0; j = 0;  k = 0
      cgl => leaves%first
      cfl_cresp_violation = .false.

      do while (associated(cgl))
         cg => cgl%cg
         p4 => cg%w(wna%fi)%arr
         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               do i = cg%is, cg%ie
                  sptab%ud = 0.0 ; sptab%ub = 0.0 ; sptab%ucmb = 0.0
                  cresp%n    = p4(iarr_cre_n, i, j, k)
                  cresp%e    = p4(iarr_cre_e, i, j, k)
                  virtual_n  => cg%w(wna%ind(vn_n))%point([i,j,k])
                  virtual_e  => cg%w(wna%ind(ve_n))%point([i,j,k])
                  if (synch_active) sptab%ub = emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k)) * bb_to_ub    !< WARNING assusmes that b is in mGs
                  if (adiab_active) sptab%ud = cg%q(qna%ind(divv_n))%point([i,j,k]) * onet
#ifdef CRESP_VERBOSED
                  print *, 'Output of cosmic ray electrons module for grid cell with coordinates i,j,k:', i, j, k
#endif /* CRESP_VERBOSED */
                  call cresp_update_cell(2 * dt, cresp%n, cresp%e, sptab, virtual_n, virtual_e, cfl_cresp_violation)
                  if ( cfl_cresp_violation ) return ! nothing to do here!
                  p4(iarr_cre_n, i, j, k) = cresp%n
                  p4(iarr_cre_e, i, j, k) = cresp%e
               enddo
            enddo
         enddo
         cg%u(iarr_cre_n, :,:,:) = p4(iarr_cre_n, :,:,:)
         cg%u(iarr_cre_e, :,:,:) = p4(iarr_cre_e, :,:,:)
         cgl=>cgl%nxt
      enddo

   end subroutine cresp_update_grid
!----------------------------------------------------------------------------------------------------
   subroutine cresp_clean_grid

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use cresp_crspectrum, only: detect_clean_spectrum
      use grid_cont,        only: grid_container
      use initcrspectrum,   only: cresp, nullify_empty_bins
      use named_array,      only: p4
      use named_array_list, only: wna

      implicit none

      integer                         :: i, j, k
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg
      logical                         :: empty_cell

      if (nullify_empty_bins) then ! else nothing is done here
         i = 0; j = 0;  k = 0
         cgl => leaves%first
         do while (associated(cgl))
            cg => cgl%cg
            p4 => cg%w(wna%fi)%arr
            do k = cg%ks, cg%ke
               do j = cg%js, cg%je
                  do i = cg%is, cg%ie
                     cresp%n    = p4(iarr_cre_n, i, j, k)
                     cresp%e    = p4(iarr_cre_e, i, j, k)
                     empty_cell = .true.

                     call detect_clean_spectrum(cresp%n, cresp%e, empty_cell)

                     p4(iarr_cre_n, i, j, k) = cresp%n
                     p4(iarr_cre_e, i, j, k) = cresp%e
                  enddo
               enddo
            enddo
            cg%u(iarr_cre_n, :,:,:) = p4(iarr_cre_n, :,:,:)
            cg%u(iarr_cre_e, :,:,:) = p4(iarr_cre_e, :,:,:)
            cgl=>cgl%nxt
         enddo
      endif

   end subroutine cresp_clean_grid
!----------------------------------------------------------------------------------------------------
   subroutine cresp_init_grid

      use cg_leaves,          only: leaves
      use cg_list,            only: cg_list_element
      use cg_list_global,     only: all_cg
      use constants,          only: pi
      use cresp_crspectrum,   only: cresp_allocate_all, e_threshold_lo, e_threshold_up, fail_count_interpol, fail_count_no_sol, &
                                    & fail_count_NR_2dim, fail_count_comp_q, second_fail, cresp_init_state
      use cresp_NR_method,    only: cresp_initialize_guess_grids
      use dataio_pub,         only: warn, printinfo, msg
      use grid_cont,          only: grid_container
      use initcosmicrays,     only: iarr_cre_n, iarr_cre_e
      use initcrspectrum,     only: e_small, e_small_approx_p_lo, e_small_approx_p_up, norm_init_spectrum, spec_mod_trms, f_init
      use mpisetup,           only: master
      use named_array_list,   only: wna
      use units,              only: cm, units_set

      implicit none

      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg
      logical, save                   :: first_run = .true., not_zeroed = .true.
      real(kind=8)                    :: sigma_T_cgs, me_cgs, myr_cgs, mGs_cgs, c_cgs, B_code_cgs_conversion

      if (first_run .eqv. .true.) then
         call cresp_initialize_guess_grids
         call cresp_allocate_all

         fail_count_interpol = 0
         fail_count_no_sol   = 0
         fail_count_NR_2dim  = 0
         second_fail         = 0
         fail_count_comp_q   = 0

         e_threshold_lo = e_small * e_small_approx_p_lo
         e_threshold_up = e_small * e_small_approx_p_up

         call all_cg%reg_var(vn_n, dim4=2) !< registering helper virtual arrays for CRESP number density
         call all_cg%reg_var(ve_n, dim4=2) !< registering helper virtual arrays for CRESP energy density


         sigma_T_cgs = 6.65245871571e-25 ! (cm ** 2)  ! < TODO: put this in the units module?
         me_cgs      = 9.1093835611e-28  ! g          ! TODO: "unitize" these quantities
         myr_cgs     = 3.1556952e+13     ! s          ! TODO: "unitize" these quantities
         mGs_cgs     = 1.0e-6            ! Gs         ! TODO: "unitize" these quantities
         c_cgs       = 29979245800.      ! cm/s       ! TODO: "unitize" these quantities
         B_code_cgs_conversion = 2.84

         if ( .not. ((trim(units_set) == "psm" ) .or. (trim(units_set) == "PSM")) ) then
            write(msg, *) "[cresp_grid:cresp_init_grid] units_set is not PSM. CRESP only works with PSM, other unit sets might cause crash."
            if (master) call warn(msg)
         endif

         bb_to_ub =  (4. / 3. ) * sigma_T_cgs / (me_cgs * c_cgs * 8. * pi) * (mGs_cgs)**2 * myr_cgs * B_code_cgs_conversion ** 2
         write (msg, *) "[cresp_grid:cresp_init_grid] 4/3 * sigma_T_cgs / ( me_cgs * c * 8 *  pi) * (mGs_cgs)**2  * myr_cgs = ", bb_to_ub         ! TODO: "unitize" these quantities
         if (master) call printinfo(msg)

         cgl => leaves%first
         do while (associated(cgl))
               cg => cgl%cg
                  cg%u(iarr_cre_n,:,:,:)  = 0.0
                  cg%u(iarr_cre_e,:,:,:)  = 0.0
                  cg%w(wna%ind(vn_n))%arr = 0.0
                  cg%w(wna%ind(ve_n))%arr = 0.0
                  not_zeroed = .false.
               cgl => cgl%nxt
            enddo

            call cresp_init_state(norm_init_spectrum%n, norm_init_spectrum%e, f_init)   !< initialize spectrum here, f_init should be 1.0

            if (master) call printinfo(" [cresp_grid:cresp_init_grid] CRESP initialized")
            first_run = .false.
      endif
      if (master) then
         if (first_run)  call warn("[cresp_grid:cresp_init_grid] CRESP might not be initialized!")
         if (not_zeroed) call warn("[cresp_grid:cresp_init_grid] CRESP virtual arrays might not be initialized properly!")
      endif

   end subroutine cresp_init_grid
!----------------------------------------------------------------------------------------------------
   subroutine grid_cresp_timestep

      use cg_leaves,          only: leaves
      use cg_list,            only: cg_list_element
      use constants,          only: xdim, ydim, zdim, one, half, onet
      use crhelpers,          only: div_v, divv_n
      use fluidindex,         only: flind
      use func,               only: emag !, operator(.equals.), operator(.notequals.)
      use grid_cont,          only: grid_container
      use initcosmicrays,     only: K_cre_paral, K_cre_perp
      use initcrspectrum,     only: spec_mod_trms, cfl_cre, synch_active, adiab_active, use_cresp
      use named_array_list,   only: qna
      use timestep_cresp,     only: cresp_timestep, dt_cre_min_ub, dt_cre_min_ud

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
      dt_cre_min_ub = huge(one)
      dt_cre_min_ud = huge(one)

      if (use_cresp) then
         cgl => leaves%first
         do while (associated(cgl))
            cg => cgl%cg
            if (adiab_active) call div_v(flind%ion%pos, cg)
            do k = cg%ks, cg%ke
               do j = cg%js, cg%je
                  do i = cg%is, cg%ie
                     sptab%ud = 0.0 ; sptab%ub = 0.0 ; sptab%ucmb = 0.0
                     if (synch_active) sptab%ub = emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k)) * bb_to_ub         !< WARNING assusmes that b is in mGs
                     if (adiab_active) sptab%ud = cg%q(qna%ind(divv_n))%point([i,j,k]) * onet

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
      endif
   end subroutine grid_cresp_timestep
!----------------------------------------------------------------------------------------------------
   subroutine append_dissipative_terms(i,j,k) ! To be fixed

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: xdim, ydim, zdim
      use crhelpers,        only: divv_n
      use func,             only: emag
      use grid_cont,        only: grid_container
      use initcrspectrum,   only: spec_mod_trms
      use named_array_list, only: qna

      implicit none

      type(spec_mod_trms)            :: sptab
      type(grid_container),  pointer :: cg
      type(cg_list_element), pointer :: cgl
      integer                        :: i,j,k
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
