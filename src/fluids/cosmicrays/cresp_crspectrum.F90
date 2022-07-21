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
!! \brief Core of Cosmic Ray Energy SPectrum (CRESP) algorithm
!<
module cresp_crspectrum
! pulled by CRESP

   use constants, only: LO, HI, I_ZERO, I_ONE

   implicit none

   private ! most of it
   public :: cresp_update_cell, cresp_init_state, cresp_get_scaled_init_spectrum, cleanup_cresp, cresp_allocate_all, &
      &      src_gpcresp, p_rch_init, detect_clean_spectrum, cresp_find_prepare_spectrum, cresp_detect_negative_content, fq_to_e, fq_to_n, q

   integer, dimension(1:2)            :: fail_count_NR_2dim, fail_count_interpol
   integer(kind=4), allocatable, dimension(:) :: fail_count_comp_q

! variables informing about change of bins
   integer, dimension(LO:HI)          :: del_i

! logical arrays / arrays determining use of p/n/e in some lines
   logical, allocatable, dimension(:) :: is_fixed_edge,   is_fixed_edge_next
   logical, allocatable, dimension(:) :: is_active_edge,  is_active_edge_next
   logical, allocatable, dimension(:) :: is_cooling_edge, is_cooling_edge_next
   logical, allocatable, dimension(:) :: is_heating_edge, is_heating_edge_next
   logical, allocatable, dimension(:) :: is_active_bin,   is_active_bin_next

! counters
   integer                            :: num_fixed_edges,   num_fixed_edges_next
   integer                            :: num_active_edges,  num_active_edges_next
   integer                            :: num_active_bins,   num_active_bins_next
   integer                            :: num_cooling_edges_next
   integer                            :: num_heating_edges_next

! dynamic arrays
   integer(kind=4), allocatable, dimension(:) :: fixed_edges,   fixed_edges_next
   integer(kind=4), allocatable, dimension(:) :: active_edges,  active_edges_next
   integer(kind=4), allocatable, dimension(:) :: active_bins,   active_bins_next
   integer(kind=4), allocatable, dimension(:) :: cooling_edges_next
   integer(kind=4), allocatable, dimension(:) :: heating_edges_next

   real, allocatable, dimension(:) :: r                                   !> r term for energy losses (Miniati 2001, eqn. 25)
   real, allocatable, dimension(:) :: q                                   !> power-law exponent array

! power-law
   real,    dimension(LO:HI)       :: p_cut_next, p_cut
   integer(kind=4), dimension(LO:HI) :: i_cut, i_cut_next
   real, allocatable, dimension(:) :: p                                   !> momentum table for piecewise power-law spectrum intervals
   real, allocatable, dimension(:) :: f                                   !> distribution function for piecewise power-law spectrum
   real, allocatable, dimension(:) :: p_next, p_upw , nflux, eflux        !> predicted and upwind momenta, number density / energy density fluxes
   real                            :: n_tot, n_tot0, e_tot, e_tot0        !> precision control for energy / number density transport and dissipation of energy
   real, allocatable, dimension(:) :: ndt, edt                            !> work array of number density and energy during algorithm execution
   real, allocatable, dimension(:) :: n, e                                !> in-algorithm energy & number density
   real, allocatable, dimension(:) :: e_amplitudes_l, e_amplitudes_r
   real,    dimension(LO:HI)       :: e_threshold                         !> lower / upper energy needed for bin activation
   integer(kind=4), dimension(LO:HI) :: approx_p                            !> if one bin, switch off cutoff p approximation
   integer(kind=4), dimension(LO:HI) :: max_ic                              !> maximum i_cut

   integer(kind=4), dimension(LO:HI), parameter:: oz = [I_ONE, I_ZERO], pm = [I_ONE, -I_ONE] !> auxiliary vectors /todo to be renamed

   abstract interface
      real function function_pointer_1D(x,y)
         real, intent(in) :: x, y
      end function function_pointer_1D
   end interface

   procedure (function_pointer_1D), pointer :: p_rch => null()

!-------------------------------------------------------------------------------------------------
!
contains
!
!-------------------------------------------------------------------------------------------------

!----- main subroutine -----

   subroutine cresp_update_cell(dt, n_inout, e_inout, sptab, cfl_cresp_violation, q1, substeps, p_out)

      use constants,      only: zero, one, I_ZERO, I_ONE
#ifdef CRESP_VERBOSED
      use dataio_pub,     only: msg, printinfo
#endif /* CRESP_VERBOSED */
      use diagnostics,    only: decr_vec
      use global,         only: disallow_CRnegatives
      use initcosmicrays, only: ncrb
      use initcrspectrum, only: allow_unnatural_transfer, crel, dfpq, e_small_approx_p, nullify_empty_bins, p_mid_fix, p_fix, spec_mod_trms

      implicit none

      real,                           intent(in)    :: dt
      real, dimension(1:ncrb),        intent(inout) :: n_inout, e_inout
      type(spec_mod_trms),            intent(in)    :: sptab
      logical,                        intent(inout) :: cfl_cresp_violation
      real, dimension(1:2), optional, intent(inout) :: p_out
      integer(kind=4), optional,      intent(in)    :: substeps

      logical                                       :: solve_fail_lo, solve_fail_up, empty_cell
      integer                                       :: i_sub, n_substep
      real, dimension(ncrb), intent(out)            :: q1

      e = zero; n = zero; edt = zero; ndt = zero
      solve_fail_lo = .false.
      solve_fail_up = .false.
      empty_cell    = .false.
      cfl_cresp_violation = .false.

      approx_p = e_small_approx_p

      p_cut_next = zero
      n_substep  = 1

      r = zero
      f = zero
      q = zero

      if (present(substeps)) then
         n_substep = substeps
      else
         n_substep = 1
      endif

      if (present(p_out)) then
         p_cut    = p_out
         p(i_cut) = p_cut
      else
         p_cut = zero
         if (dfpq%any_dump) then
            crel%p = zero
            crel%f = zero
            crel%q = zero
            crel%e = zero
            crel%n = zero
            crel%i_cut = max_ic
         endif
      endif

      call cresp_find_prepare_spectrum(n_inout, e_inout, empty_cell)

      if (empty_cell) then
         approx_p = e_small_approx_p         !< restore approximation before leaving
         if (nullify_empty_bins) then
            call nullify_all_bins(n_inout, e_inout)
         endif
#ifdef CRESP_VERBOSED
         write (msg, "(A)") "[cresp_crspectrum:cresp_update_cell] EMPTY CELL, returning"  ;  call printinfo(msg)
#endif /* CRESP_VERBOSED */
         return             ! if grid cell contains empty bins, no action is taken
      endif

! We pass values of external n_inout and e_inout to n and e after these've been preprocessed
      n = n_inout     ! number density of electrons passed to cresp module by the external module / grid
      e = e_inout     ! energy density of electrons passed to cresp module by the external module / grid

      if (approx_p(HI) > 0) then
         if (i_cut(HI) > 1) then
            call get_fqp_cutoff(HI, solve_fail_up)
         else                                                  !< spectrum cutoff beyond the fixed momentum grid
            p_cut(HI)     = p_fix(i_cut(HI))
            p(i_cut(HI))  = p_fix(i_cut(HI))
            solve_fail_up = .false.
         endif

         if (solve_fail_up) then                               !< exit_code support
            if (i_cut(HI) < ncrb) then
               if (allow_unnatural_transfer) call manually_deactivate_bin_via_transfer(i_cut(HI), -I_ONE, n, e)
               call decr_vec(active_bins, num_active_bins)
               call decr_vec(active_edges, num_active_edges)
               is_active_bin(i_cut(HI))  = .false.
               is_active_edge(i_cut(HI)) = .false.
               num_active_bins = num_active_bins - I_ONE
               f(i_cut(HI)-I_ONE:) = zero
               q(i_cut(HI))        = zero
               i_cut(HI) = i_cut(HI) - I_ONE
               p_cut(HI) = p_fix(i_cut(HI))
            else
               p_cut(HI) = p_mid_fix(i_cut(HI))
            endif
            p(i_cut(HI)) = p_cut(HI)
         endif
      endif

      if (approx_p(LO) > 0) then
         if (i_cut(LO) + 1 /= ncrb) then
            call get_fqp_cutoff(LO, solve_fail_lo)
         else                                                  !< spectrum cutoff beyond the fixed momentum grid
            p_cut(LO)     = p_fix(i_cut(LO))
            p(i_cut(LO))  = p_cut(LO)
            solve_fail_lo = .false.
         endif

         if (solve_fail_lo) then                               !< exit_code support
            if (i_cut(LO) > 0) then
               if (allow_unnatural_transfer)  call manually_deactivate_bin_via_transfer(i_cut(LO) + I_ONE, I_ONE, n, e)
               call decr_vec(active_bins, 1)
               call decr_vec(active_edges, 1)
               is_active_bin(i_cut(LO)+1) = .false.
               is_active_edge(i_cut(LO))  = .false.
               num_active_bins = num_active_bins - I_ONE
               f(i_cut(LO):)      = zero
               q(i_cut(LO)+I_ONE) = zero
               i_cut(LO) = i_cut(LO) + I_ONE
               p_cut(LO) = p_fix(i_cut(LO))
            else
               p_cut(LO) = p_mid_fix(1)
            endif
            p(i_cut(LO)) = p_cut(LO)
         endif
      endif

      if (num_active_bins < 1) then          !< if 2 active_bins and solution fails in both, return empty_cell
         approx_p = e_small_approx_p         !< restore approximation after momenta computed
         empty_cell = .true.
         return
      endif

      do i_sub = 1, n_substep                   !< if one substep, all is done the classic way
! Compute momentum changes in after time period [t,t+dt]
         call cresp_update_bin_index(sptab%ub*dt, sptab%ud*dt, p_cut, p_cut_next, cfl_cresp_violation)

         if (cfl_cresp_violation) then !< disallow_CRnegatives is not used here, as potential negatives do not appear in transfer of n,e but in p
            approx_p = e_small_approx_p         !< restore approximation after momenta computed
            call deallocate_active_arrays
#ifdef CRESP_VERBOSED
            write (msg, "(A)") "[cresp_crspectrum:cresp_update_cell] CFL violated, returning"   ;  call printinfo(msg)
#endif /* CRESP_VERBOSED */
            return                              !< returns with cfl_cresp_violation = T
         endif
! Compute fluxes through fixed edges in time period [t,t+dt], using f, q, p_cut(LO) and p_cut(HI) at [t]
! Note that new [t+dt] values of p_cut(LO) and p_cut(HI) in case new fixed edges appear or disappear.
! fill new bins
         call cresp_compute_fluxes(cooling_edges_next,heating_edges_next)

! Computing e and n at [t+dt]

         ndt(1:ncrb) = n(1:ncrb)  - (nflux(1:ncrb) - nflux(0:ncrb-1))
         edt(1:ncrb) = e(1:ncrb)  - (eflux(1:ncrb) - eflux(0:ncrb-1))

! edt(1:ncrb) = e(1:ncrb) *(one-0.5*dt*r(1:ncrb)) - (eflux(1:ncrb) - eflux(0:ncrb-1))/(one+0.5*dt*r(1:ncrb))   !!! oryginalnie u Miniatiego
! Compute coefficients R_i needed to find energy in [t,t+dt]
         call cresp_compute_r(sptab%ub, sptab%ud, p_next, active_bins_next)                 ! new active bins already received some particles, Ri is needed for those bins too

         edt(1:ncrb) = edt(1:ncrb) *(one-dt*r(1:ncrb))

         if (allow_unnatural_transfer) then
            if ((del_i(HI) == 0) .and. (approx_p(HI) > 0) .and. (i_cut_next(HI)-1 > 0)) then
               if (.not. assert_active_bin_via_nei(ndt(i_cut_next(HI)), edt(i_cut_next(HI)), i_cut_next(HI))) then
                  call manually_deactivate_bin_via_transfer(i_cut_next(HI), -I_ONE, ndt, edt)
               endif
            endif
            if ((del_i(LO) == 0) .and. (approx_p(LO) > 0) .and. (i_cut_next(LO)+2 <= ncrb)) then
               if (.not. assert_active_bin_via_nei(ndt(i_cut_next(LO)+1), edt(i_cut_next(LO)+1), i_cut_next(LO))) then
                  call manually_deactivate_bin_via_transfer(i_cut_next(LO) + I_ONE, I_ONE, ndt, edt)
               endif
            endif
         endif

         if (i_sub .ne. n_substep) then         !< proceed with end of substep
            p_cut = p_cut_next                  !< append changes for cutoffs
            i_cut = i_cut_next                  !< -//-
            p     = p_fix                       !< restore p
            p(i_cut_next(LO)) = p_cut_next(LO)  !< update cutoff p
            p(i_cut_next(HI)) = p_cut_next(HI)  !< update cutoff p

            e     = edt                         !< append changes for e
            n     = ndt                         !< append changes for e
            approx_p = [I_ZERO, I_ZERO]         !< switch off cutoff approximation

            deallocate(active_bins_next)        !< must be deallocated, reallocation in upcoming substep (update_bin_index)
            deallocate(cooling_edges_next)      !< -//-
            deallocate(heating_edges_next)      !< -//-

            if (disallow_CRnegatives) then  !< unless disallow_CRnegatives, negative n,e in CRESP are ignored
               call cresp_detect_negative_content(cfl_cresp_violation)  !< thus no point checking/returning cfl_cresp_violation here
               if (cfl_cresp_violation) then
                  call deallocate_active_arrays
#ifdef CRESP_VERBOSED
                  write (msg, "(A)") "[cresp_crspectrum:cresp_update_cell] CFL violated, returning"   ;  call printinfo(msg)
#endif /* CRESP_VERBOSED */
                  return
               endif
            endif

            call ne_to_q(n, e, q, active_bins)  !< begins new step
            f = nq_to_f(p(0:ncrb-1), p(1:ncrb), n(1:ncrb), q(1:ncrb), active_bins)  !< Compute values of distribution function in the new step
         endif

      enddo

      approx_p = e_small_approx_p         !< restore approximation after momenta computed

      p_cut = p_cut_next

#ifdef CRESP_VERBOSED
      write (msg, "(A)") "[cresp_crspectrum:cresp_update_cell] :"               ; call printinfo(msg)
      write (msg, "(A,2L2)") "[cresp_crspectrum:cresp_update_cell] solve_fail_lo, solve_fail_up:",solve_fail_lo, solve_fail_up   ; call printinfo(msg)
      write (msg, '(A5, 50E18.9)') "p_fix", p_fix      ; call printinfo(msg)
      write (msg, '(A5, 50E18.9)') "p_act", p          ; call printinfo(msg)
      write (msg, '(A5, 50E18.9)') "p_nex", p_next     ; call printinfo(msg)
      write (msg, '(A5, 50E18.9)') "p_upw", p_upw      ; call printinfo(msg)
      write (msg, '(A6, 1EN22.9, A9, 1EN22.9)') "p_lo ", p_cut(LO), ",  p_up ", p_cut(HI)  ; call printinfo(msg)

      write (msg, '(A5, 50E18.9)') "    n", n          ; call printinfo(msg)
      write (msg, '(A5, 50E18.9)') "nflux", nflux      ; call printinfo(msg)
      write (msg, '(A5, 50E18.9)') "  ndt", ndt        ; call printinfo(msg)

      write (msg, '(A5, 50E18.9)') "    e", e          ; call printinfo(msg)
      write (msg, '(A5, 50E18.9)') "eflux", eflux      ; call printinfo(msg)
      write (msg, '(A5, 50E18.9)') "  edt", edt        ; call printinfo(msg)

      write (msg, '(A5, 50E18.9)') "    r", r          ; call printinfo(msg)
      write (msg, '(A5, 50E18.9)') "    q", q          ; call printinfo(msg)
      write (msg, '(A5, 50E18.9)') "    f", f          ; call printinfo(msg)

      if (sum(approx_p) > 0) call print_failcounts

      call cresp_detect_negative_content(cfl_cresp_violation)
#endif /* CRESP_VERBOSED */

      if (disallow_CRnegatives) then   !< unless disallow_CRnegatives, negative n,e in CRESP are ignored
         call check_cutoff_ne(ndt(i_cut_next(LO) + I_ONE), edt(i_cut_next(LO) + I_ONE), i_cut_next(LO) + I_ONE, cfl_cresp_violation)
         call check_cutoff_ne(ndt(i_cut_next(HI)),         edt(i_cut_next(HI)),         i_cut_next(HI),         cfl_cresp_violation)
         if (cfl_cresp_violation) then
            call deallocate_active_arrays
#ifdef CRESP_VERBOSED
            write (msg, "(A)") "[cresp_crspectrum:cresp_update_cell] CFL violated, returning"   ;  call printinfo(msg)
#endif /* CRESP_VERBOSED */
            return
         endif
      endif

      if (nullify_empty_bins) call nullify_inactive_bins(ndt, edt)

      n = ndt
      e = edt
! --- for testing
      n_tot = sum(n)
      e_tot = sum(e)

      !print *, 'dfpq dump : ', dfpq%any_dump
      !stop
! --- saving the data to output arrays
      !if (dfpq%any_dump) then
      !  crel%p = p
      !crel%f = f
         crel%q = q ! Temporary fix for passing the q value to cresp_update_grid
      !crel%e = e
      !crel%n = n
      !crel%i_cut = i_cut
      !endif
      q1=q
      !print *, 'q : ', q1
      !stop
      n_inout = n  ! number density of electrons per bin passed back to the external module
      e_inout = e  ! energy density of electrons per bin passed back to the external module

      if (present(p_out)) then
         p_out = p_cut
      else ! if unapproximated case - clean momenta
         p_cut = zero
      endif

      call deallocate_active_arrays

   end subroutine cresp_update_cell

!----------------------------------------------------------------------------------------------------
   subroutine manually_deactivate_bin_via_transfer(ind_transfer_from, increment_transfer_to, n_array, e_array)

      use constants,       only: I_ONE
      use initcosmicrays,  only: ncrb

      implicit none

      real, dimension(I_ONE:ncrb), intent(inout) :: n_array, e_array
      integer(kind=4), intent(in)                :: ind_transfer_from, increment_transfer_to

      call transfer_quantities(n_array(ind_transfer_from), n_array(ind_transfer_from+increment_transfer_to))
      call transfer_quantities(e_array(ind_transfer_from), e_array(ind_transfer_from+increment_transfer_to))

   end subroutine manually_deactivate_bin_via_transfer

!----------------------------------------------------------------------------------------------------
   subroutine detect_clean_spectrum ! DEPRECATED

      use initcrspectrum, only: cresp

      implicit none

      logical :: empty_cell

      call find_i_bound(empty_cell)

      if (empty_cell) then
         call nullify_all_bins(cresp%n, cresp%e)
      else
         call nullify_inactive_bins(cresp%n, cresp%e)
      endif

   end subroutine detect_clean_spectrum

!----------------------------------------------------------------------------------------------------

   subroutine nullify_inactive_bins(ext_n, ext_e)

      use constants,      only: zero
      use initcosmicrays, only: ncrb

      implicit none

      real, dimension(ncrb), intent(inout) :: ext_n, ext_e

      ext_e(:i_cut(LO))   = zero
      ext_n(:i_cut(LO))   = zero
      ext_e(i_cut(HI)+1:) = zero
      ext_n(i_cut(HI)+1:) = zero

   end subroutine nullify_inactive_bins
!----------------------------------------------------------------------------------------------------
   subroutine nullify_all_bins(ext_n, ext_e)

      use constants,      only: zero
      use initcosmicrays, only: ncrb

      implicit none

      real, dimension(ncrb), intent(inout) :: ext_n, ext_e

      ext_e(:) = zero
      ext_n(:) = zero

   end subroutine nullify_all_bins

!-------------------------------------------------------------------------------------------------
! all the procedures below are called by cresp_update_cell subroutine or the driver
!-------------------------------------------------------------------------------------------------
   subroutine find_i_bound(empty_cell) ! DEPRECATED

      use constants,      only: zero, I_ONE
      use initcosmicrays, only: ncrb
      use initcrspectrum, only: cresp

      implicit none

      logical, intent(out) :: empty_cell
      integer(kind=4)      :: i

      empty_cell = .true.

      i_cut(LO) = 0
      do i = 1, ncrb                        ! if energy density is nonzero, so should be the number density
         i_cut(LO) = i - I_ONE
         if (cresp%e(i) > e_threshold(LO)) then
            if (cresp%n(i) > zero) then
               empty_cell = .false.
               exit
            endif
         endif
      enddo

      if (empty_cell) return   ! empty cell - nothing to do here!

      i_cut(HI) = ncrb
      do i = ncrb, 1, -1
         i_cut(HI) = i
         if (cresp%e(i) > e_threshold(HI)) then   ! if energy density is nonzero, so should be the number density
            if (cresp%n(i) > zero) exit
         endif
      enddo

   end subroutine find_i_bound
!-------------------------------------------------------------------------------------------------
   subroutine cresp_find_prepare_spectrum(n, e, empty_cell, i_up_out) ! EXPERIMENTAL

      use constants,      only: I_ZERO, zero, I_ONE
#ifdef CRESP_VERBOSED
      use dataio_pub,     only: msg, printinfo, warn
#endif /* CRESP_VERBOSED */
      use diagnostics,    only: incr_vec
      use initcosmicrays, only: ncrb
      use initcrspectrum, only: e_small, cresp_all_bins, p_fix, p_mid_fix

      implicit none

      real, dimension(ncrb),     intent(inout)   :: n, e
      logical,                   intent(out)     :: empty_cell
      integer(kind=4), optional, intent(out)     :: i_up_out
      integer(kind=8), dimension(:), allocatable :: nonempty_bins
      logical, dimension(ncrb)                   :: has_n_gt_zero, has_e_gt_zero
      integer(kind=4)                            :: i, num_has_gt_zero
      integer(kind=4), dimension(LO:HI)          :: approx_p_tmp, pre_i_cut

      has_n_gt_zero(:) = .false. ; has_e_gt_zero(:)  = .false.
      is_active_bin(:) = .false. ; is_active_edge(:) = .false.
      num_has_gt_zero  = I_ZERO  ; num_active_bins   = I_ZERO
      pre_i_cut = max_ic
      i_cut     = max_ic
      if (allocated(nonempty_bins)) deallocate(nonempty_bins)
      if (allocated(active_bins))   deallocate(active_bins)
      if (allocated(active_edges))  deallocate(active_edges)

! Detect where bins have nonzero values for both n and e; num_has_gt_zero stores preliminary active bins
      allocate(nonempty_bins(I_ZERO))
      do i = 1, ncrb
         has_n_gt_zero(i) = (n(i) > zero)
         has_e_gt_zero(i) = (e(i) > zero)
         if (has_n_gt_zero(i) .and. has_e_gt_zero(i)) then
            num_has_gt_zero = num_has_gt_zero + I_ONE
            call incr_vec(nonempty_bins, 1)
            nonempty_bins(num_has_gt_zero) = i
         endif
      enddo
! If cell is not empty, assume preliminary i_cut
      if (num_has_gt_zero == I_ZERO) then
         empty_cell = .true.
         return
      else
         pre_i_cut(LO) = max(int(nonempty_bins(I_ONE) - I_ONE, kind=4), I_ZERO)
         pre_i_cut(HI) = int(nonempty_bins(num_has_gt_zero), kind=4) !ubound(nonempty_bins,dim=1)
      endif

#ifdef CRESP_VERBOSED
      if (pre_i_cut(LO) == ncrb - I_ONE) then
         write(msg,*) "[cresp_crspectrum:cresp_find_prepare_spectrum] Whole spectrum moved beyond upper p_fix. Consider increasing p_up_init and restarting test."
         call warn(msg)
      endif
#endif /* CRESP_VERBOSED */

! Prepare p array
      p = zero
      p(pre_i_cut(LO)+I_ONE:pre_i_cut(HI)-I_ONE) = p_fix(pre_i_cut(LO)+I_ONE:pre_i_cut(HI)-I_ONE)
      p(pre_i_cut(LO)) = (I_ONE-approx_p(LO))*p_cut(LO) + approx_p(LO) * max(p_fix(pre_i_cut(LO)), p_mid_fix(I_ONE))  ! do not want to have zero here + p_out considered

      if (pre_i_cut(HI) < ncrb) then
         p(pre_i_cut(HI)) = (I_ONE-approx_p(HI))*p_cut(HI) + approx_p(HI) * p_fix(pre_i_cut(HI))                  ! do not want to have zero here + p_out considered
      else
         p(pre_i_cut(HI)) = (I_ONE-approx_p(HI))*p_cut(HI) + approx_p(HI) * p_mid_fix(pre_i_cut(HI))              ! do not want to have zero here + p_out considered
      endif

! preliminary allocation of active_bins
      allocate(active_bins(num_has_gt_zero))
      active_bins = int(nonempty_bins(:), kind=4)

      approx_p_tmp = approx_p                   !< Before computation of q and f for all bins approximation of cutoffs is disabled
      approx_p = I_ZERO
      i_cut = pre_i_cut                         !< make ne_to_q happy, FIXME - add cutoff indices to argument list

      call ne_to_q(n, e, q, active_bins)        !< Compute power indexes for each bin at [t] and f on left bin faces at [t]

      f = nq_to_f(p(I_ZERO:ncrb-I_ONE), p(I_ONE:ncrb), n(I_ONE:ncrb), q(I_ONE:ncrb), active_bins)  !< Compute values of distribution function f for active left edges at [t]

      approx_p = approx_p_tmp                   !< After computation of q and f for all bins approximation of cutoffs is reenabled (if was active)

      if (all(approx_p == I_ONE)) then
! compute energy density amplitudes
         e_amplitudes_l = zero   ;  e_amplitudes_r = zero
         do i = active_bins(I_ONE), active_bins(num_has_gt_zero)
            e_amplitudes_l(i) = fp_to_e_ampl(p(i-I_ONE), f(i-I_ONE))
            e_amplitudes_r(i) = fp_to_e_ampl(p(i), f(i-I_ONE) * (p(i) / p(i-I_ONE))**(-q(i)) )
         enddo

#ifdef CRESP_VERBOSED
         write (msg, "(A,50E12.4)") "[cresp_find_prepare_spectrum] e  :   ",e                   ; call printinfo(msg)
         write (msg, "(A,50E12.4)") "[cresp_find_prepare_spectrum] e_l:",      e_amplitudes_l   ; call printinfo(msg)
         write (msg, "(A,50E12.4)") "[cresp_find_prepare_spectrum] e_r:      ",e_amplitudes_r   ; call printinfo(msg)
#endif /* CRESP_VERBOSED */

! Find active bins
         is_active_bin = .false.

         do i = I_ONE, ncrb
            is_active_bin(i) = ((e_amplitudes_r(i) > e_small .or. e_amplitudes_l(i) > e_small ) .and. (e(i) > zero .and. n(i) > zero))
         enddo
         pre_i_cut = max_ic

         do i = I_ONE, ncrb
            pre_i_cut(LO) = i
            if (is_active_bin(i)) exit
         enddo
         if (pre_i_cut(LO) == ncrb) then
            empty_cell = .true.
            return
         endif
         pre_i_cut(LO) = pre_i_cut(LO) - I_ONE

         do i = ncrb, I_ONE, -I_ONE
            pre_i_cut(HI) = i
            if (is_active_bin(i)) exit
         enddo

      else
         is_active_bin(pre_i_cut(LO)+I_ONE:pre_i_cut(HI)) = .true.
      endif

      num_active_bins = count(is_active_bin)

      if (num_active_bins > I_ONE) then
         i_cut = pre_i_cut
      else if (num_active_bins == I_ONE) then
         if (i_cut(LO) > I_ZERO) then
            i_cut(HI) = active_bins(num_active_bins)
            i_cut(LO)    = i_cut(HI) -I_ONE
            p_cut(LO)    = (I_ONE-approx_p(LO))*p_cut(LO) + approx_p(LO) * p_fix(i_cut(LO));  p(i_cut(LO)) = p_cut(LO)
            approx_p(LO) = I_ZERO
         else
            i_cut(HI)    = active_bins(num_active_bins)
            p_cut(HI)    = (I_ONE-approx_p(HI))*p_cut(HI) + approx_p(HI) * p_fix(i_cut(HI));  p(i_cut(HI)) = p_cut(HI)
            approx_p(HI) = I_ZERO
         endif
      else
         empty_cell = .true.
         return
      endif

      is_active_bin(:i_cut(LO))       = .false.
      is_active_bin(i_cut(HI)+I_ONE:) = .false.

      num_active_bins = count(is_active_bin)

      if (present(i_up_out)) then
         approx_p = approx_p_tmp          !< restore approximation before leaving
         i_up_out = i_cut(HI)
         return
      endif

#ifdef CRESP_VERBOSED
      write (msg, *) "[cresp_crspectrum:cresp_find_prepare_spectrum] is_active_bin: ", is_active_bin, "|", count(is_active_bin), num_has_gt_zero      ; call printinfo(msg)
      write (msg, *) "[cresp_crspectrum:cresp_find_prepare_spectrum] i_cut:         ", i_cut
#endif /* CRESP_VERBOSED */
      if (allocated(active_bins))  deallocate(active_bins)
      if (allocated(active_edges)) deallocate(active_edges)
      if (allocated(fixed_edges))  deallocate(fixed_edges)

! allocate and prepare active bins for spectrum evolution
      if (num_active_bins > I_ZERO) then
         allocate(active_bins(num_active_bins))
         active_bins = I_ZERO
         active_bins = pack(cresp_all_bins, is_active_bin)

! Construct index arrays for fixed edges between p_cut(LO) and p_cut(HI), active edges
! before timestep

         call arrange_assoc_active_edge_arrays(fixed_edges,  is_fixed_edge,  num_fixed_edges,  i_cut+pm)
         call arrange_assoc_active_edge_arrays(active_edges, is_active_edge, num_active_edges, i_cut   )

#ifdef CRESP_VERBOSED
         write (msg, "(2(A9,i3))") "i_lo =", i_cut(LO), ", i_up = ", i_cut(HI)    ; call printinfo(msg)
#endif /* CRESP_VERBOSED */

         if (all(approx_p == I_ONE)) then
            p(:)   = p_fix(:)
            p(i_cut(LO)) = max(p_fix(i_cut(LO)), p_mid_fix(I_ONE))      ! do not want to have zero here
            p(i_cut(HI)) = max(p_fix(i_cut(HI)), p_fix(I_ONE))

            f(:i_cut(LO))  = zero
            f(i_cut(HI):)  = zero
         else
            return
         endif
      endif

   end subroutine cresp_find_prepare_spectrum
!----------------------------------------------------------------------------------------------------

   logical function assert_active_bin_via_nei(n_in, e_in, i_cutoff)

      use constants,       only: zero, fpi, one, three
      use cresp_variables, only: clight_cresp
      use cresp_NR_method, only: compute_q
      use initcosmicrays,  only: ncrb
      use initcrspectrum,  only: p_fix, e_small, eps
#ifdef CRESP_VERBOSED
      use dataio_pub,      only: printinfo, msg
#endif /* CRESP_VERBOSED */

      implicit none

      real,            intent(in) :: n_in, e_in
      integer(kind=4), intent(in) :: i_cutoff
      real                        :: f_one, q_one, q_1m3, p_l, p_r, e_amplitude_l, e_amplitude_r, alpha
      logical                     :: exit_code

      assert_active_bin_via_nei = .false.

      if (e_in > zero .and. n_in > zero .and. p_fix(max(i_cutoff-1,0)) > zero ) then

         if (i_cutoff <= 0 .or. i_cutoff >= ncrb) then
            assert_active_bin_via_nei = .true.
            return            ! WARN: returns .true. if bin of choice is the extreme one -- FIX_ME
         endif

         exit_code = .true.
         alpha = e_in/(n_in * clight_cresp * p_fix(i_cutoff-1))
         q_one = compute_q(alpha, exit_code)
         q_1m3 = three - q_one

         p_l = p_fix(i_cutoff-1)
         p_r = p_fix(i_cutoff)
         f_one = n_in / (fpi*p_l**3)

         if (abs(q_1m3) > eps) then
            f_one = f_one*(q_1m3) /((p_r/p_l)**(q_1m3) - one)
         else
            f_one = f_one/log(p_r/p_l)
         endif

         e_amplitude_l = fp_to_e_ampl(p_l, f_one)
         e_amplitude_r = fp_to_e_ampl(p_r, f_one)
         assert_active_bin_via_nei = ((e_amplitude_l > e_small) .and. (e_amplitude_r > e_small))

#ifdef CRESP_VERBOSED
         if (assert_active_bin_via_nei) then
            write(msg,*) "[cresp_crspectrum:verify_cutoff_i_next] No change to ", i_cutoff
         else
            write(msg,*) "[cresp_crspectrum:verify_cutoff_i_next] Cutoff index should change index no. ", i_cutoff
         endif
         call printinfo(msg)
#endif /* CRESP_VERBOSED */
      endif

   end function assert_active_bin_via_nei

!-----------------------------------------------------------------------

   subroutine cresp_detect_negative_content(negatives_found, location) ! Diagnostic measure - negative values should not show up:

      use constants,      only: zero, ndims
      use dataio_pub,     only: warn, msg
      use initcosmicrays, only: ncrb

      implicit none                       ! if they do, there's something wrong with last code modifications

      integer, dimension(ndims),optional :: location
      integer                            :: i
      logical, intent(inout)             :: negatives_found

      do i = 1, ncrb
         if (e(i) < zero .or. n(i) < zero .or. edt(i) < zero .or. ndt(i) < zero) then
            if (present(location)) then
               write(msg,'(A81,3I3,A7,I4,A9,E18.9,A9,E18.9)') '[cresp_crspectrum:cresp_detect_negative_content] Negative values @ (i j k ) = (', &
                           location, '): i=', i,': n(i)=', n(i), ', e(i)=',e(i)
            else
               write(msg,'(A66,A7,I4,A9,E18.9,A9,E18.9,A3,A9,I4,A9,E18.9,A9,E18.9,I2)') '[cresp_crspectrum:cresp_detect_negative_content] Negative values:',  &
                           'i=', i,': n(i)=', n(i), ', e(i)=',e(i), "|", 'i=', i,': ndt(i)=', ndt(i), ', edt(i)=',edt(i)
            endif
            negatives_found = .true. ! TODO usage with substep might prevent crash
            call warn(msg)
         endif
      enddo

   end subroutine cresp_detect_negative_content
!----------------------------------------------------------------------------------------------------
   subroutine check_cutoff_ne(n_bin, e_bin, i_bin, cfl_cresp_violated)     !< too high dt may result in negative edt, ndt (fluxes) -> cfl violation

      use constants,  only: zero
#ifdef CRESP_VERBOSED
      use dataio_pub, only: warn, msg
#endif /* CRESP_VERBOSED */

      implicit none

      real,            intent(in)  :: n_bin, e_bin
      integer(kind=4), intent(in)  :: i_bin
      logical,         intent(out) :: cfl_cresp_violated

      if (e_bin < zero .or. n_bin < zero) then
#ifdef CRESP_VERBOSED
         write(msg,'(A66,A5,E18.9,A6,E18.9,I4)')   '[cresp_crspectrum:check_cutoff_ne] Negative values:', ' n = ', n_bin, ', e = ', e_bin, i_bin
         call warn(msg)
#endif /* CRESP_VERBOSED */
         cfl_cresp_violated = .true.
      endif

      return
      if (.false.) cfl_cresp_violated = cfl_cresp_violated .or. abs(i_bin) < zero ! suppress compiler warnings

   end subroutine check_cutoff_ne

!----------------------------------------------------------------------------------------------------

   subroutine cresp_update_bin_index(ubdt, uddt, p_cut, p_cut_next, dt_too_high) ! evaluates only "next" momenta and is called after finding outer cutoff momenta

      use constants,      only: zero, I_ZERO, I_ONE, one
#ifdef CRESP_VERBOSED
      use dataio_pub,     only: msg, printinfo
#endif /* CRESP_VERBOSED */
      use initcosmicrays, only: ncrb
      use initcrspectrum, only: p_fix, cresp_all_bins, cresp_all_edges

      implicit none

      real,                   intent(in)  :: ubdt, uddt     !> in-algorithm energy dissipation terms
      real, dimension(LO:HI), intent(in)  :: p_cut
      real, dimension(LO:HI), intent(out) :: p_cut_next
      logical,                intent(out) :: dt_too_high
      integer                             :: i

      dt_too_high = .false.
! Compute p_cut at [t+dt] (update p_range)
      p_cut_next = p_cut * (one + [p_rch(uddt, ubdt*p_cut(LO)), p_rch(uddt, ubdt*p_cut(HI))]) ! changed from - to + for the sake of intuitiveness in p_rch subroutine
      p_cut_next = abs(p_cut_next)
! Compute likely cut-off indices after current timestep
      i_cut_next = get_i_cut(p_cut_next)
      i_cut_next(LO) = max(i_cut_next(LO), i_cut(LO)- I_ONE)
      i_cut_next(LO) = min(i_cut_next(LO), i_cut(LO)+ I_ONE)

      i_cut_next(HI) = max(i_cut_next(HI), i_cut(HI)- I_ONE)
      i_cut_next(HI) = min(i_cut_next(HI), i_cut(HI)+ I_ONE)

      if (p_cut_next(HI) < p_fix(i_cut_next(HI)- I_ONE)) then ! if no solution is found at the first try, approximation usually causes p_cut(HI) to jump
         dt_too_high = .true.                 ! towards higher values, which for sufficiently high dt can cause p_cut_next(HI) to even
         return                               ! become negative. As p_cut(HI) would propagate more than one bin this is clearly cfl violation.
      endif
! Detect changes in positions of lower an upper cut-ofs
      del_i = i_cut_next - i_cut

! Construct index arrays for fixed edges between p_cut(LO) and p_cut(HI), active edges after timestep
      call arrange_assoc_active_edge_arrays(fixed_edges_next,  is_fixed_edge_next,  num_fixed_edges_next,  [i_cut_next(LO) + I_ONE, i_cut_next(HI) - I_ONE])
      call arrange_assoc_active_edge_arrays(active_edges_next, is_active_edge_next, num_active_edges_next, [i_cut_next(LO),         i_cut_next(HI)])

! Active bins after timestep
      is_active_bin_next = .false.
      is_active_bin_next(i_cut_next(LO)+1:i_cut_next(HI)) = .true.
      num_active_bins_next = count(is_active_bin_next)
      allocate(active_bins_next(num_active_bins_next))
      active_bins_next = pack(cresp_all_bins, is_active_bin_next)

      p_next = zero
      p_next(fixed_edges_next) = p_fix(fixed_edges_next)
      p_next(i_cut_next) = p_cut_next

! Compute upwind momentum p_upw for all fixed edges
      p_upw = zero
      p_upw(1:ncrb) = [( p_fix(i)*(one - p_rch(uddt,ubdt*p_fix(i))), i=1,ncrb )] !< p_upw is computed with minus sign

#ifdef CRESP_VERBOSED
      write (msg, "(A, 2I3)") 'Change of  cut index lo,up:', del_i    ; call printinfo(msg)
#endif /* CRESP_VERBOSED */

! Detect cooling edges and heating edges
      is_cooling_edge_next = .false. ; num_cooling_edges_next = I_ZERO
      is_cooling_edge_next(fixed_edges_next)   = (p_upw(fixed_edges_next) &
                                                         > p_fix(fixed_edges_next))
      num_cooling_edges_next = count(is_cooling_edge_next)
      allocate(cooling_edges_next(num_cooling_edges_next))
      cooling_edges_next = pack(cresp_all_edges, is_cooling_edge_next)

      is_heating_edge_next = .false. ; num_heating_edges_next = I_ZERO
      is_heating_edge_next(fixed_edges_next) = (p_upw(fixed_edges_next) &
                                                         < p_fix(fixed_edges_next))
      num_heating_edges_next = count(is_heating_edge_next)
      allocate(heating_edges_next(num_heating_edges_next))
      heating_edges_next = pack(cresp_all_edges, is_heating_edge_next)
#ifdef CRESP_VERBOSED
      write (msg, "(A)")      'In update_bin_index'        ; call printinfo(msg)
      write (msg, "(A15,50L2, 50I3)") 'active_edges: ', is_active_edge, active_edges     ; call printinfo(msg)
      write (msg, "(A15,50L2, 50I3)") 'active edgesN:', is_active_edge_next, active_edges_next     ; call printinfo(msg)
      write (msg, "(A15,50L2, 50I3)") 'active bins  :', is_active_bin,       active_bins           ; call printinfo(msg)
      write (msg, "(A15,50L2, 50I3)") 'active binsN :', is_active_bin_next , active_bins_next      ; call printinfo(msg)
      write (msg, "(A15,50L2, 50I3)") 'fixed  edges: ', is_fixed_edge_next,  fixed_edges_next      ; call printinfo(msg)
      write (msg, "(A15,50L2, 50I3)") 'cooling edges:', is_cooling_edge_next,  cooling_edges_next  ; call printinfo(msg)
      write (msg, "(A15,50L2, 50I3)") 'heating edges:', is_heating_edge_next,  heating_edges_next  ; call printinfo(msg)
#endif /* CRESP_VERBOSED */

   end subroutine cresp_update_bin_index

!-------------------------------------------------------------------------------------------------
!
! arrays initialization | TODO: reorganize cresp_init_state
!
!-------------------------------------------------------------------------------------------------
   subroutine arrange_assoc_active_edge_arrays(index_array, where_true_array, count_is_true, range_in)

      use constants,          only: LO, HI, I_TWO
      use initcosmicrays,     only: ncrb
      use initcrspectrum,     only: cresp_all_edges

      implicit none

      integer(kind=4),          dimension(I_TWO), intent(in)   :: range_in
      integer,                                    intent(out)  :: count_is_true
      integer(kind=4), allocatable, dimension(:), intent(out)  :: index_array
      logical,                 dimension(0:ncrb), intent(out)  :: where_true_array

      where_true_array  = .false.
      where_true_array(range_in(LO):range_in(HI)) = .true.
      count_is_true  =  count(where_true_array)
      if (allocated(index_array)) deallocate(index_array)
      allocate(index_array(count_is_true))
      index_array    =  pack(cresp_all_edges, where_true_array)

   end subroutine arrange_assoc_active_edge_arrays
!-------------------------------------------------------------------------------------------------
   subroutine cresp_init_state(init_n, init_e)

      use constants,       only: zero, I_ZERO, I_ONE
      use cresp_NR_method, only: bound_name
      use dataio_pub,      only: warn, msg, die, printinfo
      use initcosmicrays,  only: ncrb, nspc
      use initcrspectrum,  only: q_init, p_init, initial_spectrum, eps, p_fix, f_init, dfpq, crel,   &
                              &  allow_source_spectrum_break, e_small_approx_init_cond, e_small_approx_p, total_init_cree, e_small, cresp_all_bins
      use mpisetup,        only: master

      implicit none

      real, dimension(I_ONE:nspc,I_ONE:ncrb), intent(out) :: init_n, init_e
      integer                                  :: i, k, i_spc
      integer(kind=4)                          :: co
      integer(kind=4), dimension(LO:HI)        :: i_ch
      real                                     :: c
      logical                                  :: exit_code

      max_ic = [I_ZERO, ncrb]
      fail_count_interpol = 0
      fail_count_NR_2dim  = 0
      fail_count_comp_q   = 0

      approx_p    = e_small_approx_p
      e_threshold = e_small_approx_p * e_small

      init_e = zero
      init_n = zero

      do i_spc = 1, nspc
        f = zero ; q = zero ; p = zero ; n = zero ; e = zero

! reading initial values of p_cut
        p_cut = p_init(:,i_spc)

        p         = p_fix       ! actual array of p including free edges, p_fix shared via initcrspectrum
        p(max_ic) = p_cut

! Sorting bin edges - arbitrary chosen p_cut may need to be sorted to appear in growing order
        do k = ncrb, 1, -1
            do i = 0, k-1
                if (p(i) > p(i+1)) then
                c = p(i)
                p(i) = p(i+1)
                p(i+1) = c
                endif
            enddo
        enddo

        i_cut = max_ic

! we only need cresp_init_state to derive (n, e) from initial (f, p_cut). For this purpose only 'active bins', i_cut are needed.

        i_cut = get_i_cut(p_cut)

        do co = LO, HI
            if (abs(p_init(co,i_spc) - p_fix(i_cut(co)+oz(co)-1)) <= eps ) then
                write(msg, *) "[cresp_crspectrum:cresp_init_state] p_init(",bound_name(co),") = p_fix(i_cut(",bound_name(co),"): (in/dec)crementing i_cut(",bound_name(co),") index to avoid FPE"
                if (master) call warn(msg)
                i_cut(co) = i_cut(co) + pm(co)
            endif
        enddo

        is_active_bin = .false.
        is_active_bin(i_cut(LO)+1:i_cut(HI)) = .true.
        num_active_bins = count(is_active_bin)
        allocate(active_bins(num_active_bins))
        active_bins = pack(cresp_all_bins, is_active_bin)

! Pure power law spectrum initial condition (default case)
        q = q_init(i_spc)
        f = f_init(i_spc) * (p/p_init(LO, i_spc))**(-q_init(i_spc))

        select case (initial_spectrum)
            case ('powl')
                call cresp_init_powl_spectrum(i_spc)
            case ('brpg')
                call cresp_init_brpg_spectrum(i_spc)
            case ('plpc')
                call cresp_init_plpc_spectrum(i_spc)
            case ('brpl')
                call cresp_init_brpl_spectrum(i_spc)
            case ('symf')
                call cresp_init_symf_spectrum(i_spc)
            case ('syme')
                call cresp_init_syme_spectrum(i_spc)
            case ('bump')
                call cresp_init_bump_spectrum(i_spc)
            case default
                write(msg,"(A,A,A)") "[cresp_crspectrum:cresp_init_state] Provided unrecognized initial_spectrum (",initial_spectrum,"). Make sure that value is correctly provided."
                call die(msg)
        end select

        if (e_small_approx_init_cond > 0) then
            do co = LO, HI
                call get_fqp_cutoff(co, exit_code)
                if (exit_code) then
                write(msg,"(a,a,a,i3)") "[cresp_crspectrum:cresp_init_state] e_small_approx_init_cond = 1, but failed to solve initial spectrum cutoff '",bound_name(co),"' for CR component #", i_spc
                call die(msg)
                endif
            enddo

            if (allow_source_spectrum_break) then

                i_ch = get_i_cut(p_cut)

                p(i_ch)     = p_cut
                q(i_ch+oz)  = q(i_cut+oz)
                f(i_ch(LO)) = e_small_to_f(p_cut(LO))
                f(i_ch(HI)) = e_small_to_f(p_cut(HI))

                do i = i_ch(LO)+1, i_cut(LO)
                p(i) = p_fix(i)
                f(i) = f(i_ch(LO)) * (p_fix(i)/p(i_ch(LO)))**(-q(i_ch(LO)+1))
                q(i+1) = q(i_ch(LO)+1)
                enddo

                do i = i_cut(HI), i_ch(HI)-1
                p(i) = p_fix(i)
                f(i) = f(i_cut(HI)-1)* (p_fix(i)/p_fix(i_cut(HI)-1))**(-q(i_cut(HI)))
                q(i) = q(i_cut(HI))
                enddo
#ifdef CRESP_VERBOSED
                write (msg,"(A,2I3,A,2I3)") "Boundary bins now (i_lo_new i_lo | i_up_new i_up)",  i_ch(LO), i_cut(LO), ' |', i_ch(HI), i_cut(HI)     ; call printinfo(msg)
#endif /* CRESP_VERBOSED */

                i_cut = i_ch
                p(i_cut(HI)) = p_cut(HI)

                is_active_bin = .false.
                is_active_bin(i_cut(LO)+1:i_cut(HI)) = .true.
                num_active_bins = count(is_active_bin) ! active arrays must be reevaluated - number of active bins and edges might have changed
                if (allocated(active_bins)) deallocate(active_bins)
                allocate(active_bins(num_active_bins)) ! active arrays must be reevaluated - number of active bins and edges might have changed
                active_bins = pack(cresp_all_bins, is_active_bin)

                e = fq_to_e(p(0:ncrb-1), p(1:ncrb), f(0:ncrb-1), q(1:ncrb), active_bins) ! once again we must count n and e
                n = fq_to_n(p(0:ncrb-1), p(1:ncrb), f(0:ncrb-1), q(1:ncrb), active_bins)
            endif
        endif

        if (dfpq%any_dump) then
            crel%p = p
            crel%f = f
            crel%q = q
            crel%e = e
            crel%n = n
            crel%i_cut = i_cut
        endif

        if (master) call check_init_spectrum(i_spc)

        !n_tot0 = sum(n)
        !e_tot0 = sum(e)

        init_n(i_spc,:) = n(:)
        init_e(i_spc,:) = e(:)

        total_init_cree(i_spc) = sum(e) !< total_init_cree value is used for initial spectrum scaling when spectrum is injected by source.
        call deallocate_active_arrays
    enddo

   end subroutine cresp_init_state

!>
!! \brief relaying e_small to f via its relation with momentum
!<
   real function e_small_to_f(p_outer)

      use constants,       only: three
      use cresp_variables, only: fpcc
      use initcrspectrum,  only: e_small

      implicit none

      real, intent(in) :: p_outer

      e_small_to_f = e_small / (fpcc * p_outer**three)

   end function e_small_to_f

!>
!! \brief Assumes power-law spectrum, without breaks. In principle the same thing is done in cresp_init_state, but init_state cannot be called from "outside".
!<
   subroutine cresp_init_powl_spectrum(i_spc)

      use constants,      only: zero
      use diagnostics,    only: my_deallocate
      use initcosmicrays, only: ncrb
      use initcrspectrum, only: cresp_all_bins, cresp_all_edges, f_init, p_fix, p_init, q_init

      implicit none

      real, dimension(0:ncrb)                    :: p_range_add
      integer(kind=4), allocatable, dimension(:) :: act_bins, act_edges
      integer(kind=4), dimension(LO:HI)          :: ic
      integer,                        intent(in) :: i_spc

      p_range_add = zero
      ic = get_i_cut(p_init(:,i_spc))

      p_range_add(ic(LO):ic(HI)) = p_fix(ic(LO):ic(HI))
      p_range_add(ic(LO)) = p_init(LO, i_spc)
      p_range_add(ic(HI)) = p_init(HI, i_spc)
      if (.not.allocated(act_edges)) allocate(act_edges(ic(HI) - ic(LO)  ))
      if (.not.allocated(act_bins )) allocate( act_bins(ic(HI) - ic(LO)+1))
      act_edges = cresp_all_edges(ic(LO)  :ic(HI))
      act_bins  = cresp_all_bins (ic(LO)+1:ic(HI))
      q(act_bins) = q_init(i_spc)

      f(act_edges) = f_init(i_spc) * (p_range_add(act_edges)/p_init(LO, i_spc))**(-q_init(i_spc))

      n = n + fq_to_n(p_range_add(0:ncrb-1), p_range_add(1:ncrb), f(0:ncrb-1), q(1:ncrb), act_bins)
      e = e + fq_to_e(p_range_add(0:ncrb-1), p_range_add(1:ncrb), f(0:ncrb-1), q(1:ncrb), act_bins)

      call my_deallocate(act_bins)
      call my_deallocate(act_edges)

   end subroutine cresp_init_powl_spectrum

!>
!! \brief "plpc": Power-law like spectrum parabolic (in log-log) cutoffs
!! \details In this case initial spectrum with a break at p_min_fix is assumed, the initial slope is parabolic
!! in ranges (p_init(LO); p_br_init(LO)) and (p_br_init(HI); p_init(HI)) and reaches e_small imposed value at cutoffs.
!<
   subroutine cresp_init_plpc_spectrum(i_spc)

      use constants,       only: zero, one, two, three, ten, I_ONE
      use cresp_variables, only: fpcc
      use diagnostics,     only: my_deallocate
      use initcosmicrays,  only: ncrb
      use initcrspectrum,  only: cresp_all_bins, e_small, f_init, p_br_init, p_fix, p_init, q_init

      implicit none

      real                                       :: c_1, c_2, c_3, lpb, lpu, lpl, a, b, lp_lpb
      real, dimension(0:ncrb)                    :: p_range_add
      integer(kind=4), allocatable, dimension(:) :: act_bins
      integer(kind=4), dimension(LO:HI)          :: ic
      integer(kind=4)                            :: i_br, i
      integer,                        intent(in) :: i_spc

      p_range_add(:) = zero
      i_br = int(minloc(abs(p_fix - p_br_init(LO, i_spc)), dim=1), kind=4) - I_ONE
      ic = get_i_cut(p_init(:,i_spc))

      p_range_add(ic(LO):ic(HI)) = p_fix(ic(LO):ic(HI))
      p_range_add(ic)            = p_init(:,i_spc)

      f(ic(LO):ic(HI)) = f_init(i_spc) * (p_range_add(ic(LO):ic(HI))/p_br_init(LO, i_spc))**(-q_init(i_spc))
      q(ic(LO):ic(HI)) = q_init(i_spc)

      if (.not.allocated(act_bins )) allocate( act_bins(ic(HI) - ic(LO)+1))
      act_bins = cresp_all_bins(ic(LO)+1:ic(HI))

      lpl = log10(p_init(LO, i_spc))
      lpb = log10(p_br_init(LO, i_spc))
      lp_lpb = lpl / lpb

      a = -q_init(i_spc)
      b = log10(f_init(i_spc) * p_br_init(LO, i_spc)**q_init(i_spc))

      c_3 =  (-three * lpl + log10(e_small / fpcc) + b * lp_lpb - a * lpl - two * b * lp_lpb ) / (lp_lpb - one)**two
      c_1 =  (c_3 - b) / lpb**two
      c_2 =  (a - two * c_1 * lpb)

      f(ic(LO):i_br-1) = ten**(c_1 * log10(p_range_add(ic(LO):i_br-1))**two + c_2 * log10(p_range_add(ic(LO):i_br-1)) + c_3)

      do i = ic(LO) + I_ONE, i_br
         q(i) = pf_to_q(p_range_add(i-1), p_range_add(i), f(i-1), f(i))
      enddo

! HIGH ENERGY CUTOFF; a and b remain unchanged
      i_br = int(minloc(abs(p_fix - p_br_init(HI, i_spc)), dim=1), kind=4)

      lpu = log10(p_init(HI, i_spc))
      lpb = log10(p_br_init(HI, i_spc))
      lp_lpb = lpu / lpb

      c_3 =  (-three * lpu + log10(e_small / fpcc) + b * lp_lpb - a * lpu - two * b * lp_lpb) / (lp_lpb - one)**two
      c_1 =  (c_3 - b) / lpb**two
      c_2 =  (a - two * c_1 * lpb)

      f(i_br:ic(HI)) = ten**(c_1 * log10(p_range_add(i_br:ic(HI)))**two + c_2 * log10(p_range_add(i_br:ic(HI))) + c_3)

      do i = i_br, ic(HI)
         q(i) = pf_to_q(p_range_add(i-1), p_range_add(i), f(i-1), f(i))
      enddo

      n = n + fq_to_n(p_range_add(0:ncrb-1), p_range_add(1:ncrb), f(0:ncrb-1), q(1:ncrb), act_bins)
      e = e + fq_to_e(p_range_add(0:ncrb-1), p_range_add(1:ncrb), f(0:ncrb-1), q(1:ncrb), act_bins)

      call my_deallocate(act_bins)
      !print *, 'n : ', n
      !print *,' e : ', e
   end subroutine cresp_init_plpc_spectrum

!>
!! \brief Power-law like spectrum with break at p_br_init_lo
!! \details In this case initial spectrum with a break at p_min_fix is assumed, the initial slope on the left side of the break is gaussian. q_br_init scales FWHM.
!<
   subroutine cresp_init_brpg_spectrum(i_spc)

      use constants,      only: I_ONE
      use initcosmicrays, only: ncrb
      use initcrspectrum, only: p_fix, p_br_init, p_init, q_br_init

      implicit none

      integer(kind=4)     :: i, i_br
      integer, intent(in) :: i_spc

      i_br = int(minloc(abs(p_fix - p_br_init(LO, i_spc)), dim=1), kind=4) - I_ONE
      f(i_cut(LO):i_br-1) = f(i_br-1) * exp(-(q_br_init(i_spc)*log(2.0) * log(p(i_cut(LO):i_br-1)/sqrt(p_init(LO, i_spc) * p(i_br)))**2))
      do i = 1, i_br
         q(i) = pf_to_q(p(i-1),p(i),f(i-1),f(i))
      enddo
      e = fq_to_e(p(0:ncrb-1), p(1:ncrb), f(0:ncrb-1), q(1:ncrb), active_bins)
      n = fq_to_n(p(0:ncrb-1), p(1:ncrb), f(0:ncrb-1), q(1:ncrb), active_bins)

   end subroutine cresp_init_brpg_spectrum

!>
!! \brief Power-law like spectrum with break at p_br_init_lo
!! \details In this case initial spectrum with a break at p_min_fix is assumed, the initial slope on the left side of the break is q_br_init.
!<
   subroutine cresp_init_brpl_spectrum(i_spc)

      use constants,      only: I_ONE
      use initcosmicrays, only: ncrb
      use initcrspectrum, only: p_fix, p_br_init, q_br_init, q_init

      implicit none

      integer(kind=4)     :: i_br
      integer, intent(in) :: i_spc

      i_br = int(minloc(abs(p_fix - p_br_init(LO, i_spc)), dim=1), kind=4) - I_ONE
      q(:i_br) = q_br_init(i_spc) ; q(i_br+1:) = q_init(i_spc)
      f(i_cut(LO):i_br-1) = f(i_br) * (p(i_cut(LO):i_br-1) / p(i_br))**(-q_br_init(i_spc))
      e = fq_to_e(p(0:ncrb-1), p(1:ncrb), f(0:ncrb-1), q(1:ncrb), active_bins)
      n = fq_to_n(p(0:ncrb-1), p(1:ncrb), f(0:ncrb-1), q(1:ncrb), active_bins)

   end subroutine cresp_init_brpl_spectrum

!>
!!\brief symf - similar to syme, but symmetric in distribution function
!<
   subroutine cresp_init_symf_spectrum(i_spc)

      use constants,      only: I_ONE
      use initcosmicrays, only: ncrb
      use initcrspectrum, only: p_fix, q_init

      implicit none

      integer(kind=4)     :: i, i_br
      integer, intent(in) :: i_spc

      i_br = int(sum(i_cut)/2, kind=4)
      q(i_cut(LO)+1:i_br) = -q_init(i_spc)
      f(i_br) = f(i_br+1)*(p(i_br+1)/p(i_br))**(-q(i_br))

      do i = 1, i_br-i_cut(LO)
         f(i_br-i) = f(i_br+i)
      enddo

      if ((i_cut(HI) - i_br /= i_br - i_cut(LO))) p_cut(HI) = p_cut(HI) - (p_cut(HI) - p_fix(i_cut(HI)-1))
      p(i_cut(HI)) = p_cut(HI) ; i_cut(HI) = i_cut(HI) - I_ONE
      e = fq_to_e(p(0:ncrb-1), p(1:ncrb), f(0:ncrb-1), q(1:ncrb), active_bins)
      n = fq_to_n(p(0:ncrb-1), p(1:ncrb), f(0:ncrb-1), q(1:ncrb), active_bins)

   end subroutine cresp_init_symf_spectrum

!>
!!\brief syme - symmetric energy distribution relative to the middle of the initial spectrum
!<
   subroutine cresp_init_syme_spectrum(i_spc)

      use constants,      only: I_ONE
      use initcosmicrays, only: ncrb
      use initcrspectrum, only: p_fix, q_init

      implicit none

      integer(kind=4)     :: i, i_br
      integer, intent(in) :: i_spc

      i_br = int(sum(i_cut)/2, kind=4)
      q(i_cut(LO)+1:i_br) = q_init(i_spc)-2.2
      do i = 1, i_br-i_cut(LO)
         f(i_br-i) = f(i_br)*(p(i_br)/p(i_br-i))**(q(i_br-i+1))
      enddo
      if ((i_cut(HI) - i_br /= i_br - i_cut(LO))) p_cut(HI) = p_cut(HI) - (p_cut(HI) - p_fix(i_cut(HI)-1))
      p(i_cut(HI)) = p_cut(HI) ; i_cut(HI) = i_cut(HI) -I_ONE
      e = fq_to_e(p(0:ncrb-1), p(1:ncrb), f(0:ncrb-1), q(1:ncrb), active_bins)
      n = fq_to_n(p(0:ncrb-1), p(1:ncrb), f(0:ncrb-1), q(1:ncrb), active_bins)

   end subroutine cresp_init_syme_spectrum

!>
!! \brief Gaussian bump-type initial condition for energy distribution
!! \todo @cresp_grid energy normalization and integral to scale cosmic ray electrons with nucleon energy density!
!<
   subroutine cresp_init_bump_spectrum(i_spc)

      use cresp_variables, only: fpcc
      use initcosmicrays,  only: ncrb
      use initcrspectrum,  only: f_init, p_init

      implicit none

      integer(kind=4)     :: i
      integer, intent(in) :: i_spc

      f = f_init(i_spc) * exp(-(4.*log(2.0)*log(p/sqrt(p_init(LO, i_spc)*p_init(HI, i_spc))/1.))**2) ! FWHM
      f(0:ncrb-1) = f(0:ncrb-1) / (fpcc * p(0:ncrb-1)**3) ! without this spectrum is gaussian for distribution function
      do i = 1, ncrb
         q(i) = pf_to_q(p(i-1),p(i),f(i-1),f(i)) !-log(f(i)/f(i-1))/log(p(i)/p(i-1))
      enddo
      e = fq_to_e(p(0:ncrb-1), p(1:ncrb), f(0:ncrb-1), q(1:ncrb), active_bins)
      n = fq_to_n(p(0:ncrb-1), p(1:ncrb), f(0:ncrb-1), q(1:ncrb), active_bins)

   end subroutine cresp_init_bump_spectrum

!-------------------------------------------------------------------------------------------------

   subroutine cresp_get_scaled_init_spectrum(n_inout, e_inout, e_in_total, i_spc) !< Using n,e spectrum obtained at initialization, obtain injected spectrum at given cell

      use initcosmicrays, only: ncrb
      use initcrspectrum, only: norm_init_spectrum_e, norm_init_spectrum_n, total_init_cree

      implicit none

      real, dimension(1:ncrb), intent(inout) :: n_inout, e_inout
      real,                    intent(in)    :: e_in_total
      integer,                 intent(in)    :: i_spc

      n_inout(:) = norm_init_spectrum_n(i_spc,:) * e_in_total / total_init_cree(i_spc)
      e_inout(:) = norm_init_spectrum_e(i_spc,:) * e_in_total / total_init_cree(i_spc)

   end subroutine cresp_get_scaled_init_spectrum

   function get_i_cut(pc) result(gic)

      use constants,       only: I_ONE
      use initcosmicrays,  only: ncrb
      use initcrspectrum,  only: p_fix, w

      implicit none

      real, dimension(LO:HI), intent(in) :: pc
      integer(kind=4), dimension(LO:HI)  :: gic
      integer(kind=4)                    :: side

      do side = LO, HI
         gic(side) = int(floor(log10(pc(side)/p_fix(1))/w), kind=4) + side
         gic(side) = max(gic(side), side - I_ONE   )
         gic(side) = min(gic(side), ncrb - oz(side))
      enddo

   end function get_i_cut

!-------------------------------------------------------------------------------------------------
!
! energy integral (eq. 21)
!
!-------------------------------------------------------------------------------------------------

   function fq_to_e(p_l, p_r, f_l, q, bins)

      use constants,       only: zero, one, four
      use cresp_variables, only: fpcc
      use initcosmicrays,  only: ncrb
      use initcrspectrum,  only: eps

      implicit none

      real,            dimension(:), intent(in) :: p_l, p_r, f_l, q
      integer(kind=4), dimension(:), intent(in) :: bins
      real,    dimension(size(bins))    :: e_bins
      real,    dimension(1:ncrb)        :: fq_to_e

      fq_to_e = zero
      e_bins = fpcc * f_l(bins) * p_l(bins)**4
      where (abs(q(bins) - four) > eps)
         e_bins = e_bins * ((p_r(bins)/p_l(bins))**(four-q(bins)) - one)/(four - q(bins))
      elsewhere
         e_bins = e_bins * log(p_r(bins)/p_l(bins))
      endwhere

      fq_to_e(bins) = e_bins

   end function fq_to_e
!-------------------------------------------------------------------------------------------------
!
! Compute edge values (amplitudes) of e
!
!-------------------------------------------------------------------------------------------------

   real function fp_to_e_ampl(p_1, f_1)

      use cresp_variables, only: fpcc2

      implicit none

      real, intent(in) :: p_1, f_1

      fp_to_e_ampl = fpcc2 * f_1 * p_1**3

   end function fp_to_e_ampl

!-------------------------------------------------------------------------------------------------
!
! density integral (eq. 9)
!
!-------------------------------------------------------------------------------------------------

   function fq_to_n(p_l, p_r, f_l, q, bins)

      use constants,      only: zero, one, three, fpi
      use initcosmicrays, only: ncrb
      use initcrspectrum, only: eps

      implicit none

      integer(kind=4), dimension(:), intent(in) :: bins
      real,            dimension(:), intent(in) :: p_l, p_r, f_l, q
      real,    dimension(size(bins))    :: n_bins
      real,    dimension(1:ncrb)        :: fq_to_n

      n_bins = zero

      n_bins = fpi*f_l(bins)*p_l(bins)**3
      where (abs(q(bins) - three) > eps)
         n_bins = n_bins*((p_r(bins)/p_l(bins))**(three-q(bins)) - one)/(three - q(bins))
      elsewhere
         n_bins = n_bins*log((p_r(bins)/p_l(bins)))
      endwhere

      fq_to_n = zero
      fq_to_n(bins) = n_bins

   end function fq_to_n

!---------------------------------------------------------------------------------------------------

   subroutine check_init_spectrum(i_spc)

      use constants,       only: one, I_ONE
      use cresp_NR_method, only: bound_name
      use dataio_pub,      only: msg, warn, printinfo
      use initcrspectrum,  only: e_small, e_small_approx_p, p_init

      implicit none

      real                   :: e_small_safe, rel_cut
      real, dimension(LO:HI) :: ec
      integer                :: co
      integer,    intent(in) :: i_spc

      e_small_safe = max(e_small, epsilon(e_small))

      do co = LO, HI
         ec = e_small_safe
         ec(co) = fp_to_e_ampl(p_init(co, i_spc), f(i_cut(co)-1+oz(co)))

         write(msg,*) "[cresp_crspectrum:check_init_spectrum] Amplitude of ", bound_name(co), " energy spectrum cutoff:", ec(co)
         call printinfo(msg)

         if (e_small_approx_p(co) == I_ONE) then
            rel_cut = (ec(LO) - ec(HI)) / ec(HI)
            write(msg,*) "[cresp_crspectrum:check_init_spectrum] Relative to e_small(", e_small_safe ,"): ", rel_cut
            if (abs(rel_cut) < one) then
               call printinfo(msg)
            else
               call warn(msg)
            endif
         endif
      enddo

   end subroutine check_init_spectrum

!-------------------------------------------------------------------------------------------------

   subroutine deallocate_active_arrays

      implicit none

      if (allocated(fixed_edges)) deallocate(fixed_edges)
      if (allocated(fixed_edges_next)) deallocate(fixed_edges_next)
      if (allocated(active_edges)) deallocate(active_edges)
      if (allocated(active_edges_next)) deallocate(active_edges_next)
      if (allocated(active_bins)) deallocate(active_bins)
      if (allocated(active_bins_next)) deallocate(active_bins_next)
      if (allocated(cooling_edges_next)) deallocate(cooling_edges_next)
      if (allocated(heating_edges_next)) deallocate(heating_edges_next)

   end subroutine deallocate_active_arrays
!-------------------------------------------------------------------------------------------------
!
! compute fluxes
!
!-------------------------------------------------------------------------------------------------
   subroutine cresp_compute_fluxes(ce,he)

      use constants,       only: zero, one, three, four, fpi
      use cresp_variables, only: fpcc
      use initcosmicrays,  only: ncrb
      use initcrspectrum,  only: eps, cresp_all_bins

      implicit none

      integer(kind=4), dimension(:), intent(in) :: ce, he    ! cooling edges, heating edges
      real,    dimension(1:ncrb-1)      :: pimh, pimth, fimh,fimth  ! *imh = i_minus_half, *imth = i_minus_third
      real,    dimension(1:ncrb-1)      :: dn_upw, de_upw, qi,qim1  ! *im1 = i_minus_one

      pimh(1:ncrb-1) = p(1:ncrb-1)
      pimth(1:ncrb-1) = p(0:ncrb-2)

      fimh(1:ncrb-1) = f(1:ncrb-1)
      fimth(1:ncrb-1) = f(0:ncrb-2)

      qi(1:ncrb-1)  = q(2:ncrb)
      qim1(1:ncrb-1) = q(1:ncrb-1)

      dn_upw = zero
      de_upw = zero
      nflux  = zero
      eflux  = zero

      dn_upw(ce) = fpi*fimh(ce)*pimh(ce)**3
      where (abs( qi(ce) - three ) > eps)
         dn_upw(ce) = dn_upw(ce)*((p_upw(ce)/pimh(ce))**(three-qi(ce)) - one)/(three - qi(ce))
      elsewhere
         dn_upw(ce) = dn_upw(ce)*log((p_upw(ce)/pimh(ce)))
      endwhere
      nflux(ce) = - dn_upw(ce)

      de_upw(ce) = fpcc * fimh(ce) * pimh(ce)**4
      where (abs(qi(ce) - four) > eps)
         de_upw(ce) = de_upw(ce)*((p_upw(ce)/pimh(ce))**(four-qi(ce)) - one)/(four - qi(ce))
      elsewhere
         de_upw(ce) = de_upw(ce)*log(p_upw(ce)/pimh(ce))
      endwhere
      eflux(ce) =  - de_upw(ce)

      if (del_i(HI) == -1) then
         nflux(i_cut(HI)-1) = -n(i_cut(HI))
         eflux(i_cut(HI)-1) = -e(i_cut(HI))
      endif

! filling empty empty bin - switch of upper boundary, condition is checked only once per flux computation and is very rarely satisfied.
      if (nflux(i_cut(HI)) > zero) then             ! If flux is greater than zero it will go through right edge, activating next bin in the next timestep.
         if (cresp_all_bins(i_cut(HI)+1) == i_cut(HI)+1) then  ! But it should only happen if there is bin with index i_up+1
            ndt(i_cut(HI)+1) = nflux(i_cut(HI))
            edt(i_cut(HI)+1) = eflux(i_cut(HI))
            del_i(HI) = +1
         endif
      endif

      if (nflux(i_cut(HI)-1)+n(i_cut(HI)) <= zero) then ! If flux is equal or greater than energy / density in a given bin,  these both shall migrate
         nflux(i_cut(HI)-1) =  -n(i_cut(HI))                   ! to an adjacent bin, thus making given bin detected as inactive (empty) in the next timestep
         eflux(i_cut(HI)-1) =  -e(i_cut(HI))
         del_i(HI) = -1
      endif

      dn_upw(he) = fpi*fimth(he)*p_upw(he)**3*(pimth(he)/p_upw(he))**qim1(he)
      where (abs(qim1(he) - three) > eps)
         dn_upw(he) = dn_upw(he)*((pimh(he)/p_upw(he))**(three-qim1(he)) - one)/(three - qim1(he))
      elsewhere
         dn_upw(he) = dn_upw(he)*log((pimh(he)/p_upw(he)))
      endwhere
      nflux(he) = dn_upw(he)

      de_upw(he) = fpcc * fimth(he) * p_upw(he)**4*(pimth(he)/p_upw(he))**qim1(he)
      where (abs(qi(he) - four) > eps)
         de_upw(he) = de_upw(he)*((pimh(he)/p_upw(he))**(four-qim1(he)) - one)/(four - qim1(he))
      elsewhere
         de_upw(he) = de_upw(he)*log(pimh(he)/p_upw(he))
      endwhere
      eflux(he) = de_upw(he)

      if (del_i(LO) == 1 .or. nflux(i_cut(LO)+1) >= n(i_cut(LO)+1)) then
         nflux(i_cut(LO)+1) = n(i_cut(LO)+1)
         eflux(i_cut(LO)+1) = e(i_cut(LO)+1)
! emptying lower boundary bins - in cases when flux gets greater than energy or number density
         del_i(LO) = 1   ! in case it hasn't yet been modified
      endif

   end subroutine cresp_compute_fluxes

!-------------------------------------------------------------------------------------------------
!
! compute R (eq. 25)
!
!-------------------------------------------------------------------------------------------------
   subroutine cresp_compute_r(u_b, u_d, p, bins)

      use constants,      only: zero, four, five
      use initcosmicrays, only: ncrb
      use initcrspectrum, only: eps

      implicit none

      integer(kind=4), dimension(:), intent(in) :: bins
      real,                          intent(in) :: u_b, u_d
      real, dimension(0:ncrb),       intent(in) :: p
      real, dimension(size(bins))               :: r_num, r_den

      r = zero

      ! Found here an FPE occurring in mcrwind/mcrwind_cresp
      ! bins = [ 11, 12, 13, 14, 15 ]
      ! q = [ 0, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, 4.922038984459582, -30 ]
      ! five-q(bins) = 35, that seems to be a bit high power to apply carelessly

      where (abs(q(bins) - five) > eps)
         r_num = (p(bins)**(five-q(bins)) - p(bins-1)**(five-q(bins)))/(five - q(bins))
      elsewhere
         r_num = log(p(bins)/p(bins-1))
      endwhere

      where (abs(q(bins) - four) > eps)
         r_den = (p(bins)**(four-q(bins)) - p(bins-1)**(four-q(bins)))/(four - q(bins))
      else where
         r_den = log(p(bins)/p(bins-1))
      endwhere

      where (abs(r_num) > zero .and. abs(r_den) > zero)                  !< BEWARE - regression: comparisons against
         r(bins) = u_d + u_b * r_num/r_den !all cooling effects will come here   !< eps and epsilon result in bad results;
      endwhere                                                                  !< range of values ofr_num and r_den is very wide

   end subroutine cresp_compute_r

!-------------------------------------------------------------------------------------------------
!
! find new q (eq. 29) and new f
!
!-------------------------------------------------------------------------------------------------

   subroutine ne_to_q(n, e, q, bins)

      use constants,       only: zero, I_ONE
      use dataio_pub,      only: warn
      use cresp_NR_method, only: compute_q
      use cresp_variables, only: clight_cresp
      use initcosmicrays,  only: ncrb
      use initcrspectrum,  only: e_small

      implicit none

      real, dimension(1:ncrb),       intent(in)  :: n, e
      real, dimension(1:ncrb),       intent(out) :: q
      integer(kind=4), dimension(:), intent(in)  :: bins
      integer                                    :: i, i_active
      real                                       :: alpha_in
      logical                                    :: exit_code

      q = zero

      do i_active = 1 + approx_p(LO), size(bins) - approx_p(HI)
         i = bins(i_active)
         if (e(i) > e_small .and. p(i-1) > zero) then
            exit_code = .true.
            if (abs(n(i)) < 1e-300) call warn("[cresp_crspectrum:ne_to_q] 1/|n(i)| > 1e300")
            ! n(i) of order 1e-100 does happen sometimes, but extreme values like 4.2346894890376292e-312 tend to create FPE in the line below
            ! these could be uninitialized values
            alpha_in = e(i)/(n(i)*p(i-1)*clight_cresp)
            if ((i == i_cut(LO)+1) .or. (i == i_cut(HI))) then ! for boundary case, when momenta are not approximated
               q(i) = compute_q(alpha_in, exit_code, p(i)/p(i-1))
            else
               q(i) = compute_q(alpha_in, exit_code)
            endif
         else
            q(i) = zero
         endif
         if (exit_code) fail_count_comp_q(i) = fail_count_comp_q(i) + I_ONE
      enddo

   end subroutine ne_to_q

!-------------------------------------------------------------------------------------------------
! Function used to obtain q for one cell out of f and p values - used to compute q after finding p_cut(HI)
!-------------------------------------------------------------------------------------------------
   real function pf_to_q(p_l, p_r, f_l, f_r)

      implicit none

      real, intent(in) :: p_l, p_r, f_l, f_r

      pf_to_q = -log(f_r/f_l)/log(p_r/p_l) ! append value of q for given p_cut(HI)

   end function pf_to_q

! -------------------------------------------------------------------------------------------------
!
! distribution function amplitudes (eq. 9)
!
! -------------------------------------------------------------------------------------------------
   function nq_to_f(p_l, p_r, n, q, bins)

      use constants,      only: zero, one, three, fpi
      use initcosmicrays, only: ncrb
      use initcrspectrum, only: eps

      implicit none

      integer(kind=4), dimension(:) :: bins
      real, dimension(1:ncrb)     :: p_l, p_r, n, q
      real, dimension(size(bins)) :: f_bins
      real, dimension(0:ncrb)     :: nq_to_f
      real, dimension(0:ncrb)     :: pr_by_pl   ! the array of values of p_r/p_l to avoid FPEs

      nq_to_f  = zero
      f_bins   = zero
      pr_by_pl = one

      where (p_r(bins) > zero .and. p_l(bins) > zero ) ! p(i) = 0 in inactive bins. This condition should be met by providing proper "bins" range - FIXME
         pr_by_pl(bins) = p_r(bins) / p_l(bins)                     ! + comparing reals with zero is still risky
         f_bins = n(bins) / (fpi*p_l(bins)**3)
         where (abs(q(bins)-three) > eps)
            f_bins = f_bins*(three - q(bins)) /((pr_by_pl(bins))**(three-q(bins)) - one)
         elsewhere
            f_bins = f_bins/log((p_r(bins)/p_l(bins)))
         endwhere
      endwhere

      nq_to_f(bins-1) = f_bins

   end function nq_to_f

!---------------------------------------------------------------------------------------------------
! Computing cosmic ray pressure (eq. 44)
!---------------------------------------------------------------------------------------------------
   real function get_pcresp(p_l, p_r, f_l, q, bins) ! computes cre pressure, not used currently

      use constants,       only: one, four
      use cresp_variables, only: fp3cc
      use initcrspectrum,  only: eps

      implicit none

      real,    dimension(:), intent(in) :: p_l, p_r, f_l, q
      integer, dimension(:), intent(in) :: bins
      real,    dimension(size(bins))    :: p_cresp

      p_cresp = fp3cc * f_l(bins)*p_l(bins)**4

      where (abs(q(bins) - four) > eps)
         p_cresp = p_cresp*((p_r(bins)/p_l(bins))**(four-q(bins)) - one)/(four - q(bins))
      elsewhere
         p_cresp = p_cresp*log(p_r(bins)/p_l(bins))
      endwhere

      get_pcresp = sum(p_cresp)

   end function get_pcresp
!---------------------------------------------------------------------------------------------------
! Computing cosmic ray spectrum component pressure gradient
!---------------------------------------------------------------------------------------------------
   subroutine src_gpcresp(u, n, dx, grad_pcresp, i_spc)

      use constants,      only: onet
      use initcosmicrays, only: ncrb
      use initcrspectrum, only: cre_active

      implicit none

      integer(kind=4),            intent(in)  :: n, i_spc
      real,                       intent(in)  :: dx
      real, dimension(n, 1:ncrb), intent(in)  :: u
      real, dimension(n),         intent(out) :: grad_pcresp
      real, dimension(n)                      :: P_cresp_r, P_cresp_l, P_cresp

      grad_pcresp = 0.0

      P_cresp_l = 0.0 ;  P_cresp_r = 0.0

!     (ultrarelativistic) !TODO Expand me to trans-relativistic
      !P_cresp_l(1:n-2) = onet * sum(u(1:n-2, :),dim=2)
      !P_cresp_r(3:n)   = onet * sum(u(3:n,   :),dim=2)
      !grad_pcresp(2:n-1) = cre_active(i_spc) * (P_cresp_l(1:n-2) - P_cresp_r(3:n) )/(2.*dx)

      P_cresp(:) = onet * sum(u(:, :),dim=2)
      grad_pcresp(2:n-1) = cre_active(i_spc) * (P_cresp(1:n-2) - P_cresp(3:n) )/(2.*dx)

      !print *, 'cre_active(i_spc) ', i_spc, ' ', cre_active(i_spc)

      !print *, 'Pressure gradient : ', grad_pcresp(2:n-1)
   end subroutine src_gpcresp
!---------------------------------------------------------------------------------------------------
! Preparation and computation of boundary momenta and and boundary
! distribution function amplitudes value on left bin edge, computing q for indicated cutoff bin (returns f, p, q)
!---------------------------------------------------------------------------------------------------
   subroutine get_fqp_cutoff(cutoff, exit_code)

      use constants,       only: zero, one, I_ONE, I_TWO, LO, HI
      use cresp_NR_method, only: intpol_pf_from_NR_grids, alpha, n_in, NR_algorithm, q_ratios, assoc_pointers
      use cresp_variables, only: clight_cresp
#ifdef CRESP_VERBOSED
      use cresp_NR_method, only: bound_name
      use dataio_pub,      only: msg, printinfo
#endif /* CRESP_VERBOSED */
      use initcrspectrum,  only: e_small, q_big, p_fix, NR_refine_pf

      implicit none

      integer(kind=4), intent(in)  :: cutoff
      logical,         intent(out) :: exit_code
      integer(kind=4)              :: ipfix, qi
      real, dimension(1:2)         :: x_NR, x_NR_init
      logical                      :: interpolated

      ipfix = i_cut(cutoff) + pm(cutoff)
      qi    = i_cut(cutoff) + oz(cutoff)
      call assoc_pointers(cutoff)

      alpha = e(qi)/(n(qi) * p_fix(ipfix) * clight_cresp)
      n_in  = n(qi)
      x_NR = intpol_pf_from_NR_grids(alpha, n_in, interpolated)
      if (.not. interpolated) then
         exit_code = .true.
         fail_count_interpol(cutoff) = fail_count_interpol(cutoff) + 1
         return
      else
         x_NR_init = x_NR
      endif

#ifdef CRESP_VERBOSED
      write (msg, "(A27,A2,A2,2E22.15)") "Input ratios(p, f) for NR (", bound_name(cutoff), "):", x_NR  ; call printinfo(msg)
#endif /* CRESP_VERBOSED */
      exit_code = .false.
      if (NR_refine_pf(cutoff)) then
         call NR_algorithm(x_NR, exit_code)
         if (exit_code) then ! some failures still take place
            fail_count_NR_2dim(cutoff) = fail_count_NR_2dim(cutoff) + I_ONE
            x_NR = x_NR_init
         endif
      endif

      x_NR = abs(x_NR) ! negative values cannot be allowed
      if (x_NR(I_ONE) < one) then
         exit_code = .true.
         return
      endif

      select case (cutoff)
         case (LO)
            p_cut(cutoff) = p_fix(ipfix) / x_NR(I_ONE)
            f(qi-1)       = e_small_to_f(p_cut(cutoff))
         case (HI)
            p_cut(cutoff) = p_fix(ipfix) * x_NR(I_ONE)
            f(qi-1)       = e_small_to_f(p_cut(cutoff)) / x_NR(I_TWO)
      end select
      p(i_cut(cutoff)) = p_cut(cutoff)
      q(qi)            = q_ratios(x_NR(I_TWO), x_NR(I_ONE))

      if (abs(q(qi)) > q_big) q(qi) = sign(one, q(qi)) * q_big
#ifdef CRESP_VERBOSED
      write (msg, "(A26,2E22.15)") " >>> Obtained (p, f):", p_cut(cutoff), f(qi-1) ; call printinfo(msg)
      write (msg, "(A26,2E22.15)") "     Corresponding ratios:", x_NR              ; call printinfo(msg)
      call printinfo(msg)
#endif /* CRESP_VERBOSED */
      alpha = zero ;  n_in = zero

   end subroutine get_fqp_cutoff

!>
!! \brief Relative change of momentum due to losses (u_b*p*dt) and compression u_d*dt (Taylor expansion up to 3rd order)
!<
   real function p_rch_ord_1(uddt, ubpdt)

      implicit none

      real, intent(in) :: uddt, ubpdt

      p_rch_ord_1 = -(uddt + ubpdt)

   end function p_rch_ord_1
!-------------------------------------------------------------------------------------------------
   real function p_rch_ord_2_1(uddt, ubpdt)     !< adds 2nd term and calls 1st order

      use constants, only: half

      implicit none

      real, intent(in) :: uddt, ubpdt

      p_rch_ord_2_1 = p_rch_ord_1(uddt, ubpdt) + ( half*uddt**2 + ubpdt**2)

   end function p_rch_ord_2_1
!-------------------------------------------------------------------------------------------------
   real function p_rch_ord_3_2_1(uddt, ubpdt)     !< adds 3rd term and calls 2nd and 1st order

      use constants, only: onesth

      implicit none

      real, intent(in) :: uddt, ubpdt

      p_rch_ord_3_2_1 = p_rch_ord_2_1(uddt, ubpdt) - onesth * uddt**3 - ubpdt**3

   end function p_rch_ord_3_2_1
!----------------------------------------------------------------------------------------------------
   subroutine p_rch_init

      use dataio_pub,     only: msg, die
      use initcrspectrum, only: expan_order

      implicit none

      select case (expan_order)
         case (1)
            p_rch => p_rch_ord_1
         case (2)
            p_rch => p_rch_ord_2_1
         case (3)
            p_rch => p_rch_ord_3_2_1
         case default
            write (msg,*) '[cresp_crspectrum:p_rch_init] expan_order =',expan_order,': value incorrect (accepted values [1;3]).'
            call die(msg)
      end select
   end subroutine p_rch_init
!----------------------------------------------------------------------------------------------------
   subroutine transfer_quantities(take_from, give_to)

      use constants, only: zero

      implicit none

      real, intent(inout) :: give_to, take_from

      give_to = give_to + take_from
      take_from = zero

   end subroutine transfer_quantities
!----------------------------------------------------------------------------------------------------
   subroutine cresp_allocate_all

      use constants,      only: I_ZERO, I_ONE
      use diagnostics,    only: my_allocate_with_index
      use initcosmicrays, only: ncrb

      implicit none

      integer(kind=4) :: ma1d

      ma1d = ncrb
      call my_allocate_with_index(fail_count_comp_q,ma1d, I_ONE)

      call my_allocate_with_index(n,ma1d, I_ONE)   !:: n, e, r
      call my_allocate_with_index(e,ma1d, I_ONE)
      call my_allocate_with_index(e_amplitudes_l,ma1d, I_ONE)   !:: n, e, r
      call my_allocate_with_index(e_amplitudes_r,ma1d, I_ONE)
      call my_allocate_with_index(r,ma1d, I_ONE)
      call my_allocate_with_index(q,ma1d, I_ONE)

      call my_allocate_with_index(f,ma1d, I_ZERO)
      call my_allocate_with_index(p,ma1d, I_ZERO)

      call my_allocate_with_index(edt,ma1d, I_ONE)
      call my_allocate_with_index(ndt,ma1d, I_ONE)

      call my_allocate_with_index(p_next,ma1d, I_ZERO)
      call my_allocate_with_index(p_upw,ma1d, I_ZERO)
      call my_allocate_with_index(nflux,ma1d, I_ZERO)
      call my_allocate_with_index(eflux,ma1d, I_ZERO)

      call my_allocate_with_index(is_fixed_edge,ma1d, I_ZERO)
      call my_allocate_with_index(is_fixed_edge_next,ma1d, I_ZERO)
      call my_allocate_with_index(is_active_edge,ma1d, I_ZERO)
      call my_allocate_with_index(is_active_edge_next,ma1d, I_ZERO)
      call my_allocate_with_index(is_cooling_edge,ma1d, I_ZERO)
      call my_allocate_with_index(is_cooling_edge_next,ma1d, I_ZERO)
      call my_allocate_with_index(is_heating_edge,ma1d, I_ZERO)
      call my_allocate_with_index(is_heating_edge_next,ma1d, I_ZERO)
      call my_allocate_with_index(is_active_bin,ma1d, I_ONE)
      call my_allocate_with_index(is_active_bin_next,ma1d, I_ONE)

   end subroutine cresp_allocate_all
!----------------------------------------------------------------------------------------------------
   subroutine cresp_deallocate_all

      use diagnostics, only: my_deallocate

      implicit none

      call my_deallocate(fail_count_comp_q)

      call my_deallocate(n)
      call my_deallocate(e)
      call my_deallocate(e_amplitudes_l)
      call my_deallocate(e_amplitudes_r)
      call my_deallocate(r)
      call my_deallocate(q)
      call my_deallocate(f)
      call my_deallocate(p)

      call my_deallocate(edt)
      call my_deallocate(ndt)

      call my_deallocate(p_next)
      call my_deallocate(p_upw)
      call my_deallocate(nflux)
      call my_deallocate(eflux)

      call my_deallocate(is_fixed_edge)
      call my_deallocate(is_fixed_edge_next)
      call my_deallocate(is_active_edge)
      call my_deallocate(is_active_edge_next)
      call my_deallocate(is_cooling_edge)
      call my_deallocate(is_cooling_edge_next)
      call my_deallocate(is_heating_edge)
      call my_deallocate(is_heating_edge_next)
      call my_deallocate(is_active_bin)
      call my_deallocate(is_active_bin_next)

      call deallocate_active_arrays

   end subroutine cresp_deallocate_all

!----------------------------------------------------------------------------------------------------

   subroutine cleanup_cresp

      implicit none

#ifdef CRESP_VERBOSED
      call print_failcounts   ! since active bin counting method changed some time ago, failcounts are overestimated, possibly DEPRECATED
#endif /* CRESP_VERBOSED */
      call cresp_deallocate_all

   end subroutine cleanup_cresp

!----------------------------------------------------------------------------------------------------

   subroutine print_failcounts

      use dataio_pub, only: msg, printinfo

      implicit none

      write(msg, '(A36,I8,A6,I8)') "NR_2dim:  convergence failure: p_lo", fail_count_NR_2dim(LO),  ", p_up", fail_count_NR_2dim(HI)  ; call printinfo(msg)
      write(msg, '(A36,I8,A6,I8)') "NR_2dim:interpolation failure: p_lo", fail_count_interpol(LO), ", p_up", fail_count_interpol(HI) ; call printinfo(msg)
      write(msg, '(A36,   100I8)') "NR_2dim:inpl/solve  q(bin) failure:", fail_count_comp_q                                          ; call printinfo(msg)

   end subroutine print_failcounts

end module cresp_crspectrum
