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
! pulled by COSM_RAY_ELECTRONS
   use dataio_pub,      only: msg, die, warn, printinfo

   implicit none

   private ! most of it
   public :: cresp_update_cell, cresp_init_state, printer, fail_count_interpol, fail_count_NR_2dim, cresp_get_scaled_init_spectrum, &
      &      cleanup_cresp, cresp_accuracy_test, b_losses, cresp_allocate_all, cresp_deallocate_all, e_threshold_lo, e_threshold_up, &
      &      fail_count_comp_q, src_gpcresp, cresp_init_powl_spectrum, p_rch_init,  &
      &      detect_clean_spectrum, cresp_find_prepare_spectrum, cresp_detect_negative_content

   integer, dimension(1:2), save      :: fail_count_NR_2dim, fail_count_interpol
   integer, allocatable, save         :: fail_count_comp_q(:)

! variables informing about change of bins
   integer                            :: del_i_lo, del_i_up

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
   integer, allocatable               :: fixed_edges(:),   fixed_edges_next(:)
   integer, allocatable               :: active_edges(:),  active_edges_next(:)
   integer, allocatable               :: active_bins(:),   active_bins_next(:)
   integer, allocatable               :: cooling_edges_next(:)
   integer, allocatable               :: heating_edges_next(:)

   real(kind=8), allocatable, dimension(:) :: r  ! r term for energy losses (Miniati 2001, eqn. 25)
   real(kind=8), allocatable, dimension(:) :: q  ! power-law exponent array

! power-law
   real(kind=8)                       :: p_lo_next, p_up_next, p_lo, p_up !, p_lo_bef, p_up_bef
   integer                            :: i_lo, i_up, i_lo_next, i_up_next
   real(kind=8), dimension(:), allocatable   :: p ! momentum table for piecewise power-law spectru intervals
   real(kind=8), dimension(:), allocatable   :: f ! distribution function for piecewise power-law spectrum

! predicted and upwind momenta, number density / energy density fluxes
   real(kind=8), dimension(:),allocatable   :: p_next, p_upw , nflux, eflux ! , p_fix

! precision control for energy / number density transport and dissipation of energy
   real(kind=8)                             :: n_tot, n_tot0, e_tot, e_tot0

! terms for energy dissipation tests
   real(kind=8)                 :: u_d_0
   real(kind=8)                 :: u_b_0

! in-algorithm energy dissipation terms
   real(kind=8)              :: u_b, u_d

! work array of number density and energy during algorithm execution
   real(kind=8),allocatable, dimension(:)    :: ndt, edt

! in-algorithm energy & number density
   real(kind=8), allocatable, dimension(:)  :: n, e ! dimension(1:ncre)
   real(kind=8), allocatable, dimension(:)  :: e_amplitudes_l, e_amplitudes_r
! lower / upper energy needed for bin activation
   real(kind=8), save :: e_threshold_lo, e_threshold_up
! if one bin, switch off cutoff p approximation
   integer  :: approx_p_lo, approx_p_up


   abstract interface
      function function_pointer_1D(x,y)
         real(kind=8)            :: function_pointer_1D
         real(kind=8),intent(in) :: x, y
      end function function_pointer_1D
   end interface

   procedure (function_pointer_1D), pointer :: p_rch => null()

!-------------------------------------------------------------------------------------------------
!
contains
!
!-------------------------------------------------------------------------------------------------

!----- main subroutine -----

   subroutine cresp_update_cell(dt, n_inout, e_inout, sptab, cfl_cresp_violation, p_out)

      use constants,      only: zero, one
      use diagnostics,    only: decr_vec
      use initcrspectrum, only: ncre, spec_mod_trms, e_small_approx_p_lo, e_small_approx_p_up, crel, p_mid_fix, nullify_empty_bins, p_fix

      implicit none

      real(kind=8),                           intent(in)    :: dt
      real(kind=8), dimension(1:ncre),        intent(inout) :: n_inout, e_inout
      type(spec_mod_trms),                    intent(in)    :: sptab
      logical,                                intent(inout) :: cfl_cresp_violation
      real(kind=8), dimension(1:2), optional, intent(inout) :: p_out
      logical                                               :: solve_fail_lo, solve_fail_up, empty_cell

      e = zero; n = zero; edt = zero; ndt = zero
      solve_fail_lo = .false.
      solve_fail_up = .false.
      empty_cell    = .false.
      cfl_cresp_violation = .false.

      approx_p_lo = e_small_approx_p_lo
      approx_p_up = e_small_approx_p_up

      p_lo_next = zero
      p_up_next = zero

      r = zero
      f = zero
      q = zero

      u_b = sptab%ub
      u_d = sptab%ud

      if (present(p_out)) then
         p_lo = p_out(1) ; p(i_lo) = p_lo
         p_up = p_out(2) ; p(i_up) = p_up
      else
         p_lo = zero
         p_up = zero
      endif

      call cresp_find_prepare_spectrum(n_inout, e_inout, empty_cell)

      if ( empty_cell ) then
         approx_p_lo = e_small_approx_p_lo         !< restore approximation before leaving
         approx_p_up = e_small_approx_p_up         !< restore approximation before leaving
         if (nullify_empty_bins) then
            call nullify_all_bins(n_inout, e_inout)
         endif
#ifdef CRESP_VERBOSED
         write (msg, "(A)") "[cresp_crspectrum:cresp_update_cell] EMPTY CELL, returning"   ;  call printinfo(msg)
#endif /* CRESP_VERBOSED */
         return             ! if grid cell contains empty bins, no action is taken
      endif

! We pass values of external n_inout and e_inout to n and e after these've been preprocessed
      n = n_inout     ! number density of electrons passed to cresp module by the external module / grid
      e = e_inout     ! energy density of electrons passed to cresp module by the external module / grid

      if (approx_p_up .gt. 0) then
         if (i_up .gt. 1) then
            call get_fqp_up(solve_fail_up)
         else                                                  !< spectrum cutoff beyond the fixed momentum grid
            p_up = p_fix(i_up)
            p(i_up) = p_fix(i_up)
            solve_fail_up = .false.
         endif

         if (solve_fail_up) then                               !< exit_code support
            if (i_up .lt. ncre) then
               call transfer_quantities(n(i_up-1),n(i_up))
               call transfer_quantities(e(i_up-1),e(i_up))
               call decr_vec(active_bins, num_active_bins)
               call decr_vec(active_edges, num_active_edges)
               is_active_bin(i_up) = .false.
               is_active_edge(i_up) = .false.
               num_active_bins = num_active_bins - 1
               i_up         = i_up - 1
               p_up         = p_fix(i_up) ; p(i_up)     = p_up
            else
               p_up         = p_mid_fix(i_up);  p(i_up)     = p_up
            endif
         endif
      endif

      if (approx_p_lo .gt. 0) then
         if ((i_lo+1) .ne. ncre) then
            call get_fqp_lo(solve_fail_lo)
         else                                                  !< spectrum cutoff beyond the fixed momentum grid
            p_lo           = p_fix(i_lo)
            p(i_lo)        = p_lo
            solve_fail_lo  = .false.
         endif

         if (solve_fail_lo) then                               !< exit_code support
            if (i_lo .gt. 0) then
               call transfer_quantities(n(i_lo+2),n(i_lo+1))
               call transfer_quantities(e(i_lo+2),e(i_lo+1))
               call decr_vec(active_bins, 1)
               call decr_vec(active_edges, 1)
               is_active_bin(i_lo+1) = .false.
               is_active_edge(i_lo) = .false.
               num_active_bins = num_active_bins - 1
               i_lo         = i_lo + 1
               p_lo         = p_fix(i_lo) ;     p(i_lo)     = p_lo
            else
               p_lo         = p_mid_fix(1);     p(i_lo)     = p_lo
            endif
         endif
      endif

      if (num_active_bins .lt. 1 ) then                        !< if 2 active_bins and solution fails in both, return empty_cell
         approx_p_lo = e_small_approx_p_lo         !< restore approximation after momenta computed
         approx_p_up = e_small_approx_p_up         !< restore approximation after momenta computed
         empty_cell = .true.
         return
      endif

      call cresp_update_bin_index(dt, p_lo, p_up, p_lo_next, p_up_next, cfl_cresp_violation)
      if ( cfl_cresp_violation ) then
         approx_p_lo = e_small_approx_p_lo         !< restore approximation after momenta computed
         approx_p_up = e_small_approx_p_up         !< restore approximation after momenta computed
         call deallocate_active_arrays
#ifdef CRESP_VERBOSED
         write (msg, "(A)") "[cresp_crspectrum:cresp_update_cell] EMPTY CELL, returning"   ;  call printinfo(msg)
#endif /* CRESP_VERBOSED */
         return
      endif
! Compute fluxes through fixed edges in time period [t,t+dt], using f, q, p_lo and p_up at [t]
! Note that new [t+dt] values of p_lo and p_up in case new fixed edges appear or disappear.
! fill new bins
      call cresp_compute_fluxes(cooling_edges_next,heating_edges_next)

! Computing e and n at [t+dt]

      ndt(1:ncre) = n(1:ncre)  - (nflux(1:ncre) - nflux(0:ncre-1))
      edt(1:ncre) = e(1:ncre)  - (eflux(1:ncre) - eflux(0:ncre-1))

! edt(1:ncre) = e(1:ncre) *(one-0.5*dt*r(1:ncre)) - (eflux(1:ncre) - eflux(0:ncre-1))/(one+0.5*dt*r(1:ncre))   !!! oryginalnie u Miniatiego
! Compute coefficients R_i needed to find energy in [t,t+dt]
      call cresp_compute_r(p_next, active_bins_next)                 ! new active bins already received some particles, Ri is needed for those bins too

      edt(1:ncre) = edt(1:ncre) *(one-dt*r(1:ncre))

      if ((del_i_up .eq. 0) .and. (approx_p_up .gt. 0)) then
         if (assert_active_bin_via_nei(ndt(i_up_next), edt(i_up_next), i_up_next) .eqv. .false.) then
            call transfer_quantities(ndt(i_up_next-1),ndt(i_up_next))
            call transfer_quantities(edt(i_up_next-1),edt(i_up_next))
         endif
      endif

      approx_p_lo = e_small_approx_p_lo         !< restore approximation after momenta computed
      approx_p_up = e_small_approx_p_up         !< restore approximation after momenta computed

      p_lo = p_lo_next
      p_up = p_up_next

#ifdef CRESP_VERBOSED
      write (msg, "(A)") "[cresp_crspectrum:cresp_update_cell] :"               ; call printinfo(msg)
      write (msg, "(A,2L2)") "[cresp_crspectrum:cresp_update_cell] solve_fail_lo, solve_fail_up:",solve_fail_lo, solve_fail_up   ; call printinfo(msg)
      write (msg, '(A5, 50E18.9)') "p_fix", p_fix      ; call printinfo(msg)
      write (msg, '(A5, 50E18.9)') "p_act", p          ; call printinfo(msg)
      write (msg, '(A5, 50E18.9)') "p_nex", p_next     ; call printinfo(msg)
      write (msg, '(A5, 50E18.9)') "p_upw", p_upw      ; call printinfo(msg)
      write (msg, '(A6, 1EN22.9, A9, 1EN22.9)') "p_lo ", p_lo, ",  p_up ", p_up  ; call printinfo(msg)

      write (msg, '(A5, 50E18.9)') "    n", n          ; call printinfo(msg)
      write (msg, '(A5, 50E18.9)') "nflux", nflux      ; call printinfo(msg)
      write (msg, '(A5, 50E18.9)') "  ndt", ndt        ; call printinfo(msg)

      write (msg, '(A5, 50E18.9)') "    e", e          ; call printinfo(msg)
      write (msg, '(A5, 50E18.9)') "eflux", eflux      ; call printinfo(msg)
      write (msg, '(A5, 50E18.9)') "  edt", edt        ; call printinfo(msg)

      write (msg, '(A5, 50E18.9)') "    r", r          ; call printinfo(msg)
      write (msg, '(A5, 50E18.9)') "    q", q          ; call printinfo(msg)
      write (msg, '(A5, 50E18.9)') "    f", f          ; call printinfo(msg)

      if ( (approx_p_lo+approx_p_up) .gt. 0 ) then
         write (msg, '(A36,I5,A6,I3)') "NR_2dim:  convergence failure: p_lo", fail_count_NR_2dim(1), ", p_up", fail_count_NR_2dim(2)      ; call printinfo(msg)
         write (msg, '(A36,I5,A6,I3)') "NR_2dim:interpolation failure: p_lo", fail_count_interpol(1), ", p_up", fail_count_interpol(2)    ; call printinfo(msg)
         write (msg, '(A36,   100I5)') "NR_2dim:inpl/solve  q(bin) failure:", fail_count_comp_q                                           ; call printinfo(msg)
      endif
      call cresp_detect_negative_content
#endif /* CRESP_VERBOSED */

      call check_cutoff_ne(ndt(i_lo_next+1), edt(i_lo_next+1), i_lo_next+1, cfl_cresp_violation)
      call check_cutoff_ne(ndt(i_up_next),   edt(i_up_next), i_up_next, cfl_cresp_violation)

      if ( cfl_cresp_violation ) then
         call deallocate_active_arrays
#ifdef CRESP_VERBOSED
         write (msg, "(A)") "[cresp_crspectrum:cresp_update_cell] EMPTY CELL, returning"
         call printinfo(msg)
#endif /* CRESP_VERBOSED */
         return
      endif

      if (nullify_empty_bins) call nullify_inactive_bins(ndt, edt)

      n = ndt
      e = edt
! --- for testing
      n_tot = sum(n)
      e_tot = sum(e)

      crel%p = p
      crel%f = f
      crel%q = q
      crel%e = e
      crel%n = n
      crel%i_lo = i_lo
      crel%i_up = i_up
! --- saving the data to output arrays
      n_inout  = n  ! number density of electrons per bin passed back to the external module
      e_inout  = e  ! energy density of electrons per bin passed back to the external module

      if (present(p_out)) then
         p_out(1) = p_lo
         p_out(2) = p_up
      endif

      call deallocate_active_arrays

      p_lo = zero
      p_up = zero

   end subroutine cresp_update_cell

!----------------------------------------------------------------------------------------------------
   subroutine detect_clean_spectrum(ext_n, ext_e, empty_cell) ! DEPRECATED

      use initcrspectrum,      only: ncre, nullify_empty_bins

      implicit none

      real(kind=8), dimension(ncre), intent(inout) :: ext_n, ext_e
      logical, intent(inout)  :: empty_cell

      call find_i_bound(ext_n, ext_e, empty_cell)

      if (empty_cell) then
         if (nullify_empty_bins) then
            call nullify_all_bins(ext_n, ext_e)
         endif
      endif
      if (nullify_empty_bins) then
         call nullify_inactive_bins(ext_n, ext_e)
      endif

   end subroutine detect_clean_spectrum

!----------------------------------------------------------------------------------------------------

   subroutine nullify_inactive_bins(ext_n, ext_e)

      use initcrspectrum,      only: ncre
      use constants,           only: zero

      implicit none

      real(kind=8), dimension(ncre), intent(inout) :: ext_n, ext_e

      ext_e(:i_lo)   = zero
      ext_n(:i_lo)   = zero
      ext_e(i_up+1:) = zero
      ext_n(i_up+1:) = zero

   end subroutine nullify_inactive_bins
!----------------------------------------------------------------------------------------------------
   subroutine nullify_all_bins(ext_n, ext_e)

      use initcrspectrum,      only: ncre
      use constants,           only: zero

      implicit none

      real(kind=8), dimension(ncre), intent(inout) :: ext_n, ext_e

      ext_e(:)   = zero
      ext_n(:)   = zero

   end subroutine nullify_all_bins

!-------------------------------------------------------------------------------------------------
! all the procedures below are called by cresp_update_cell subroutine or the driver
!-------------------------------------------------------------------------------------------------
   subroutine find_i_bound(ext_n, ext_e, empty_cell) ! DEPRECATED

      use initcrspectrum, only: ncre
      use constants, only: zero

      implicit none

      integer(kind=4)   :: i
      logical           :: empty_cell
      real(kind=8), dimension(ncre), intent(inout) :: ext_n, ext_e

      empty_cell    = .true.

      i_lo = 0
      do i = 1, ncre                        ! if energy density is nonzero, so should be the number density
         i_lo = i-1
         if ( ext_e(i) .gt. e_threshold_lo) then
           if (ext_n(i) .gt. zero) then
              empty_cell = .false.
              exit
           endif
        endif
     enddo

     if ( empty_cell .eqv. .true.) return   ! empty cell - nothing to do here!

     i_up = ncre
     do i = ncre, 1,-1
        i_up = i
        if (ext_e(i) .gt. e_threshold_up ) then   ! if energy density is nonzero, so should be the number density
           if ( ext_n(i) .gt. zero ) exit
        endif
     enddo

   end subroutine find_i_bound
!-------------------------------------------------------------------------------------------------
   subroutine cresp_find_prepare_spectrum(n, e, empty_cell, i_up_out) ! EXPERIMENTAL

      use constants,      only: I_ZERO, zero, I_ONE
      use diagnostics,    only: incr_vec
      use initcrspectrum, only: ncre, e_small, cresp_all_edges, cresp_all_bins, p_fix, p_mid_fix

      implicit none

      integer(kind=8), dimension(:), allocatable :: nonempty_bins
      logical, dimension(ncre) :: has_n_gt_zero, has_e_gt_zero
      logical                  :: empty_cell
      integer(kind=4)          :: i, pre_i_lo, pre_i_up, num_has_gt_zero, approx_p_lo_tmp, approx_p_up_tmp
      integer(kind=4), optional:: i_up_out
      real(kind=8), dimension(ncre), intent(inout)   :: n, e

      has_n_gt_zero(:) = .false. ; has_e_gt_zero(:)  = .false.
      is_active_bin(:) = .false. ; is_active_edge(:) = .false.
      num_has_gt_zero  = I_ZERO  ; num_active_bins   = I_ZERO
      pre_i_lo         = I_ZERO  ; pre_i_up          = ncre
      i_lo             = I_ZERO  ; i_up              = ncre
      if (allocated(nonempty_bins)) deallocate(nonempty_bins)
      if (allocated(active_bins))   deallocate(active_bins)
      if (allocated(active_edges))  deallocate(active_edges)

! Detect where bins have nonzero values for both n and e; num_has_gt_zero stores preliminary active bins
      allocate(nonempty_bins(I_ZERO))
      do i = 1, ncre
         has_n_gt_zero(i) = (n(i) .gt. zero)
         has_e_gt_zero(i) = (e(i) .gt. zero)
         if (has_n_gt_zero(i) .and. has_e_gt_zero(i)) then
            num_has_gt_zero = num_has_gt_zero + I_ONE
            call incr_vec(nonempty_bins, I_ONE)
            nonempty_bins(num_has_gt_zero) = i
         endif
      enddo
! If cell is not empty, assume preliminary i_lo and i_up
      if (num_has_gt_zero .eq. I_ZERO) then
         empty_cell = .true.
         return
      else
         pre_i_lo = max(int(nonempty_bins(I_ONE) - I_ONE,kind=4), I_ZERO)
         pre_i_up = int(nonempty_bins(num_has_gt_zero),kind=4) !ubound(nonempty_bins,dim=1)
      endif

      if (pre_i_lo .eq. (ncre - I_ONE)) then
         write(msg,*) "[cresp_crspectrum:cresp_find_prepare_spectrum] Whole spectrum moved beyond upper p_fix. Consider increasing p_up_init and restarting test."
         call warn(msg)
      endif

! Prepare p array
      p = zero
      p(pre_i_lo+I_ONE:pre_i_up-I_ONE) = p_fix(pre_i_lo+I_ONE:pre_i_up-I_ONE)
      p(pre_i_lo) = (I_ONE-approx_p_lo)*p_lo + approx_p_lo * max(p_fix(pre_i_lo), p_mid_fix(I_ONE))  ! do not want to have zero here + p_out considered

      if (pre_i_up .lt. ncre) then
         p(pre_i_up) = (I_ONE-approx_p_up)*p_up + approx_p_up * p_fix(pre_i_up)                  ! do not want to have zero here + p_out considered
      else
         p(pre_i_up) = (I_ONE-approx_p_up)*p_up + approx_p_up * p_mid_fix(pre_i_up)              ! do not want to have zero here + p_out considered
      endif

! preliminary allocation of active_bins
      allocate(active_bins(num_has_gt_zero))
      active_bins = int(nonempty_bins(:), kind=4)

      approx_p_lo_tmp = approx_p_lo                   !< Before computation of q and f for all bins approximation of cutoffs is disabled
      approx_p_up_tmp = approx_p_up
      approx_p_lo = I_ZERO
      approx_p_up = I_ZERO
      i_lo = pre_i_lo      ;  i_up = pre_i_up         !< make ne_to_q happy, FIXME - add cutoff indices to argument list

      call ne_to_q(n,e,q,active_bins)                                         !< Compute power indexes for each bin at [t] and f on left bin faces at [t]

      f = nq_to_f(p(I_ZERO:ncre-I_ONE), p(I_ONE:ncre), n(I_ONE:ncre), q(I_ONE:ncre), active_bins)  !< Compute values of distribution function f for active left edges at [t]

      approx_p_lo = approx_p_lo_tmp
      approx_p_up = approx_p_up_tmp                   !< After computation of q and f for all bins approximation of cutoffs is reenabled (if was active)

      if (approx_p_lo .eq. I_ONE .and. approx_p_up .eq. I_ONE) then
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

         do i = I_ONE, ncre
            is_active_bin(i) = ( (e_amplitudes_r(i) .gt. e_small .or. e_amplitudes_l(i) .gt. e_small ) .and. (e(i) > zero .and. n(i) > zero) )
         enddo
         pre_i_lo = I_ZERO ; pre_i_up = ncre

         do i = I_ONE, ncre
            pre_i_lo = i
            if (is_active_bin(i)) exit
         enddo
         if (pre_i_lo .eq. ncre) then
            empty_cell = .true.
            return
         endif
         pre_i_lo = pre_i_lo - I_ONE

         do i = ncre,I_ONE,-I_ONE
            pre_i_up = i
            if (is_active_bin(i)) exit
         enddo

      else
         is_active_bin(pre_i_lo+I_ONE:pre_i_up) = .true.
      endif

      num_active_bins = count(is_active_bin)

      if (num_active_bins .gt. I_ONE) then
         i_lo = pre_i_lo;   i_up = pre_i_up
      else if (num_active_bins .eq. I_ONE) then
         if (i_lo .gt. I_ZERO) then
            i_up = active_bins(num_active_bins)
            i_lo = i_up -I_ONE
            p_lo        = (I_ONE-approx_p_lo)*p_lo + approx_p_lo * p_fix(i_lo);  p(i_lo) = p_lo
            approx_p_lo = I_ZERO
         else
            i_up        = active_bins(num_active_bins)
            p_up        = (I_ONE-approx_p_up)*p_up + approx_p_up * p_fix(i_up);  p(i_up) = p_up
            approx_p_up = I_ZERO
         endif
      else
         empty_cell = .true.
         return
      endif

      is_active_bin(:i_lo)       = .false.
      is_active_bin(i_up+I_ONE:) = .false.

      num_active_bins = count(is_active_bin)

      if (present(i_up_out)) then
         approx_p_lo = approx_p_lo_tmp          !< restore approximation before leaving
         approx_p_up = approx_p_up_tmp          !< restore approximation before leaving
         i_up_out = i_up
         return
      endif

#ifdef CRESP_VERBOSED
      write (msg, *) "[cresp_crspectrum:cresp_find_prepare_spectrum] is_active_bin: ", is_active_bin, "|", count(is_active_bin), num_has_gt_zero      ; call printinfo(msg)
      write (msg, *) "[cresp_crspectrum:cresp_find_prepare_spectrum] i_lo, i_up:    ", i_lo, i_up
#endif /* CRESP_VERBOSED */
      if (allocated(active_bins))  deallocate(active_bins)
      if (allocated(active_edges)) deallocate(active_edges)
      if (allocated(fixed_edges))  deallocate(fixed_edges)

! allocate and prepare active bins for spectrum evolution
      if (num_active_bins .gt. I_ZERO) then
         allocate(active_bins(num_active_bins))
         active_bins = I_ZERO
         active_bins = pack(cresp_all_bins, is_active_bin)

! Construct index arrays for fixed edges betwen p_lo and p_up, active edges
! before timestep
         is_fixed_edge = .false.
         is_fixed_edge(i_lo+I_ONE:i_up-I_ONE) = .true.
         num_fixed_edges = count(is_fixed_edge)
         allocate(fixed_edges(num_fixed_edges))
         fixed_edges = pack(cresp_all_edges, is_fixed_edge)

         is_active_edge = .false.
         is_active_edge(i_lo:i_up) = .true.
         num_active_edges = count(is_active_edge)
         allocate(active_edges(i_lo:i_up))
         active_edges = pack(cresp_all_edges, is_active_edge)

#ifdef CRESP_VERBOSED
         write (msg, "(2(A9,i3))") "i_lo =", i_lo, ", i_up = ", i_up    ; call printinfo(msg)
#endif /* CRESP_VERBOSED */

         if (approx_p_lo .eq. I_ONE .and. approx_p_up .eq. I_ONE) then
            p(:)   = p_fix(:)
            p(i_lo) = max(p_fix(i_lo), p_mid_fix(I_ONE))      ! do not want to have zero here
            p(i_up) = max(p_fix(i_up), p_fix(I_ONE))

            f(:i_lo)  = zero
            f(i_up:)  = zero
         else
            return
         endif
      endif

   end subroutine cresp_find_prepare_spectrum
!----------------------------------------------------------------------------------------------------

   function assert_active_bin_via_nei(n_in, e_in, i_cutoff)
      use initcrspectrum,  only: p_fix, ncre, e_small, eps
      use constants,       only: zero, fpi, one, three
      use cresp_variables, only: clight ! use units,     only: clight
      use cresp_NR_method, only: compute_q
#ifdef CRESP_VERBOSED
      use dataio_pub,      only: printinfo, msg
#endif /* CRESP_VERBOSED */

      implicit none

      integer(kind=4),intent(in)    :: i_cutoff
      real(kind=8), dimension(1)    :: e_one, n_one, f_one, q_one, p_l, p_r
      real(kind=8)                  :: e_amplitude_l, e_amplitude_r, alpha, e_in, n_in
      logical                       :: exit_code, assert_active_bin_via_nei


      if (e_in .gt. zero .and. n_in .gt. zero .and. p_fix(i_cutoff-1) .gt. zero ) then
         assert_active_bin_via_nei = .true.

         e_one(1) = e_in   ;   n_one(1) = n_in
         exit_code = .true.

         if (i_cutoff .gt. 0 .and. i_cutoff .lt. ncre) then
            alpha =  e_one(1)/(n_one(1) * clight * p_fix(i_cutoff-1))
            q_one(1) =  compute_q(alpha, exit_code)
         else
            return            ! WARN: returns .true. if bin of choice is the extreme one -- FIX_ME
         endif

         p_l(1) = p_fix(i_cutoff-1)  ;   p_r(1) = p_fix(i_cutoff)
         f_one(1) = n_one(1) / (fpi*p_l(1)**3)

         if (abs(q_one(1)-three) .gt. eps) then
            f_one(1) = f_one(1)*(three - q_one(1)) /(( p_r(1)/p_l(1))**(three - q_one(1)) - one)
         else
            f_one(1) = f_one(1)/log(p_r(1)/p_l(1))
         endif

         e_amplitude_l = fp_to_e_ampl(p_l(1), f_one(1))
         e_amplitude_r = fp_to_e_ampl(p_r(1), f_one(1))

         if ( (e_amplitude_l .gt. e_small) .and.  (e_amplitude_r .gt. e_small) ) then
#ifdef CRESP_VERBOSED
            write(msg,*) "[cresp_crspectrum:verify_cutoff_i_next] No change to ", i_cutoff ; call printinfo(msg)
#endif /* CRESP_VERBOSED */
            assert_active_bin_via_nei = .true.
         else
#ifdef CRESP_VERBOSED
            write(msg,*) "[cresp_crspectrum:verify_cutoff_i_next] Cutoff index should change index no. ", i_cutoff; call printinfo(msg)
#endif /* CRESP_VERBOSED */
            assert_active_bin_via_nei = .false.
         endif
      else
         assert_active_bin_via_nei = .false.
         return
      endif

   end function assert_active_bin_via_nei

!-----------------------------------------------------------------------

   subroutine cresp_detect_negative_content(location) ! Diagnostic measure - negative values should not show up:
      use constants,       only: zero, ndims
      use dataio_pub,      only: warn, msg
      use initcrspectrum,  only: ncre

      implicit none                       ! if they do, there's something wrong with last code modifications

      integer :: i
      integer, dimension(ndims),optional      :: location

      do i = 1, ncre
         if (e(i) .lt. zero .or. n(i) .lt. zero .or. edt(i) .lt. zero .or. ndt(i) .lt. zero) then
            if (present(location)) then
               write(msg,'(A81,3I3,A7,I4,A9,E18.9,A9,E18.9)') '[cresp_crspectrum:cresp_detect_negative_content] Negative values @ (i j k ) = (', &
                           location, '): i=', i,': n(i)=', n(i), ', e(i)=',e(i)
            else
               write(msg,'(A66,A7,I4,A9,E18.9,A9,E18.9,A3,A9,I4,A9,E18.9,A9,E18.9)') '[cresp_crspectrum:cresp_detect_negative_content] Negative values:',  &
                           'i=', i,': n(i)=', n(i), ', e(i)=',e(i), "|", 'i=', i,': ndt(i)=', ndt(i), ', edt(i)=',edt(i)
            endif
            call warn(msg)
         endif
      enddo

   end subroutine cresp_detect_negative_content
!----------------------------------------------------------------------------------------------------
   subroutine check_cutoff_ne(n_bin, e_bin, i_bin, cfl_violated)     !< too high dt may result in negative edt, ndt (fluxes) -> cfl violation
      use constants,       only: zero
#ifdef CRESP_VERBOSED
      use dataio_pub,      only: warn, msg
#endif /* CRESP_VERBOSED */
      implicit none

      logical                    :: cfl_violated
      real(kind=8), intent(in)   :: n_bin, e_bin
      integer(kind=4),intent(in) :: i_bin

      if (e_bin .lt. zero .or. n_bin .lt. zero) then
#ifdef CRESP_VERBOSED
         write(msg,'(A66,A5,E18.9,A6,E18.9,I4)')   '[cresp_crspectrum:cresp_detect_negative_content] Negative values:', &
                                             &     ' n = ', n_bin, ', e = ', e_bin, i_bin
         call warn(msg)
#endif /* CRESP_VERBOSED */
         cfl_violated = .true.
      endif

   end subroutine check_cutoff_ne

!----------------------------------------------------------------------------------------------------

   subroutine cresp_update_bin_index(dt, p_lo, p_up, p_lo_next, p_up_next, dt_too_high) ! evaluates only "next" momenta and is called after finding outer cutoff momenta

      use constants, only: zero, I_ZERO, one
      use initcrspectrum, only: ncre, p_fix, w, cresp_all_bins, cresp_all_edges

      implicit none

      real(kind=8), intent(in)  :: dt
      real(kind=8), intent(in)  :: p_lo, p_up
      real(kind=8), intent(out) :: p_lo_next, p_up_next
      logical,      intent(out) :: dt_too_high
      integer                   :: i

      dt_too_high = .false.
! Compute p_lo and p_up at [t+dt]
      call p_update(dt, p_lo, p_lo_next)
      call p_update(dt, p_up, p_up_next)
      p_lo_next = abs(p_lo_next)
      p_up_next = abs(p_up_next)
! Compute likely cut-off indices after current timestep
      i_lo_next = int(floor(log10(p_lo_next/p_fix(1))/w)) + 1
      i_lo_next = max(0, i_lo_next, i_lo-1)
      i_lo_next = min(i_lo_next, ncre - 1, i_lo+1)

      i_up_next = int(floor(log10(p_up_next/p_fix(1))/w)) + 2
      i_up_next = max(1,i_up_next, i_up-1)
      i_up_next = min(i_up_next,ncre,i_up+1)

      if ( p_up_next .lt. p_fix(i_up_next-1) ) then ! if no solution is found at the first try, approximation usually causes p_up to jump
         dt_too_high = .true.                 ! towards higher values, which for sufficiently high dt can cause p_up_next to even
         return                               ! become negative. As p_up would propagate more than one bin this is clearly cfl violation.
      endif
! Detect changes in positions of lower an upper cut-ofs
      del_i_lo = i_lo_next - i_lo
      del_i_up = i_up_next - i_up

! Construct index arrays for fixed edges betwen p_lo and p_up, active edges
! after timestep
      is_fixed_edge_next = .false.
      is_fixed_edge_next(i_lo_next+1:i_up_next-1) = .true.
      num_fixed_edges_next = count(is_fixed_edge_next)
      allocate(fixed_edges_next(num_fixed_edges_next))
      fixed_edges_next = pack(cresp_all_edges, is_fixed_edge_next)

      is_active_edge_next = .false.
      is_active_edge_next(i_lo_next:i_up_next) = .true.
      num_active_edges_next = count(is_active_edge_next)
      allocate(active_edges_next(num_active_edges_next))
      active_edges_next = pack(cresp_all_edges, is_active_edge_next)

! Active bins after timestep
      is_active_bin_next = .false.
      is_active_bin_next(i_lo_next+1:i_up_next) = .true.
      num_active_bins_next = count(is_active_bin_next)
      allocate(active_bins_next(num_active_bins_next))
      active_bins_next = pack(cresp_all_bins, is_active_bin_next)

      p_next = zero
      p_next(fixed_edges_next) = p_fix(fixed_edges_next)
      p_next(i_lo_next) = p_lo_next
      p_next(i_up_next) = p_up_next

! Compute upwind momentum p_upw for all fixed edges
      p_upw = zero
      p_upw(1:ncre) = (/( p_fix(i)*(one - p_rch(dt,p_fix(i))), i=1,ncre )/) !< p_upw is computed with minus sign

#ifdef CRESP_VERBOSED
      write (msg, "(A, 2I3)") 'Change of  cut index lo,up:', del_i_lo, del_i_up    ; call printinfo(msg)
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
! update p_range
!
!-------------------------------------------------------------------------------------------------

   subroutine p_update(dt, p_old,  p_new)

      use constants, only: one

      implicit none

      real(kind=8), intent(in)  :: dt
      real(kind=8), intent(in)  :: p_old
      real(kind=8), intent(out) :: p_new

      p_new = p_old*(one + p_rch(dt, p_old)) ! changed from - to + for the sake of intuitiveness in p_rch subroutine

   end subroutine p_update
!-------------------------------------------------------------------------------------------------
!
! arrays initialization | TODO: reorganize cresp_init_state
!
!-------------------------------------------------------------------------------------------------
   subroutine cresp_init_state(init_n, init_e, f_amplitude, sptab)

      use constants, only: zero, I_ONE, fpi, three
      use cresp_NR_method, only: e_small_to_f
      use cresp_variables, only: clight ! use units, only: clight
      use dataio_pub,      only: warn, msg, die
      use initcrspectrum,  only: ncre, spec_mod_trms, q_init, p_lo_init, p_up_init, initial_condition, eps, p_fix, w,   &
                              &  allow_source_spectrum_break, e_small_approx_init_cond, e_small_approx_p_lo, crel,      &
                              &  e_small_approx_p_up, total_init_cree, e_small, cresp_all_bins,     &
                              &  q_br_init, p_br_init
      use mpisetup,        only: master

      implicit none

      integer                          :: i, k, i_lo_ch, i_up_ch, i_br
      real(kind=8)                     :: c
      real(kind=8), dimension(I_ONE:ncre)    :: init_n, init_e
      type (spec_mod_trms), intent(in), optional :: sptab
      real(kind=8), intent(in)         :: f_amplitude
      logical :: exit_code

      u_b = zero ; u_d = zero

      if (present(sptab)) u_b = sptab%ub
      if (present(sptab)) u_d = sptab%ud

      approx_p_lo = e_small_approx_p_lo
      approx_p_up = e_small_approx_p_up

      init_e = zero
      init_n = zero

      if (present(sptab)) u_b_0 = u_b
      if (present(sptab)) u_d_0 = u_d

      f = zero ; q = zero ; p = zero ; n = zero ; e = zero

      q = q_init
! reading initial values of p_lo and p_up
      p_lo = p_lo_init
      p_up = p_up_init

      p        = p_fix       ! actual array of p including free edges, p_fix shared via initcrspectrum
      p(0)     = p_lo
      p(ncre)  = p_up

! Sorting bin edges - arbitrary chosen p_lo and p_up may need to be sorted to appear in growing order
      do k = ncre, 1, -1
         do i = 0, k-1
            if (p(i)>p(i+1)) then
               c = p(i)
               p(i) = p(i+1)
               p(i+1) = c
            endif
         enddo
      enddo

      i_lo = 0
      i_up = ncre

! we only need cresp_init_state to derive (n, e) from initial (f, p_lo, p_up). For this purpose only 'active bins', i_lo & i_up are needed.

      i_lo = int(floor(log10(p_lo/p_fix(1))/w)) + 1
      i_lo = max(0, i_lo)
      i_lo = min(i_lo, ncre - 1)

      i_up = int(floor(log10(p_up/p_fix(1))/w)) + 2
      i_up = max(1,i_up)
      i_up = min(i_up,ncre)

      if (abs(p_lo_init - p_fix(i_lo)) .le. eps ) then
         write(msg, *) "[cresp_crspectrum:cresp_init_state] p_lo_init = p_fix(i_lo):  incrementing i_lo index to avoid FPE"
         if (master) call warn(msg)
         i_lo = i_lo+1
      endif

      if (abs(p_up_init - p_fix(i_up-1)) .le. eps ) then
         write(msg, *) "[cresp_crspectrum:cresp_init_state] p_up_init = p_fix(i_up-1): decrementing i_up index to avoid FPE"
         if (master) call warn(msg)
         i_up = i_up-1
      endif

      if (q_init .lt. three) then
         if (approx_p_lo .eq. I_ONE .or. approx_p_up .eq. I_ONE) then
            write(msg,*) "[cresp_crspectrum:cresp_init_state] Initial parameters: q_init < 3.0 and approximation of outer momenta is on,"
            call warn(msg)
            write(msg,*) "[cresp_crspectrum:cresp_init_state] approximation of outer momenta with hard energy spectrum might not work. You have been warned."
            call warn(msg)
         endif
      endif

      is_active_bin = .false.
      is_active_bin(i_lo+1:i_up) = .true.
      num_active_bins = count(is_active_bin)
      allocate(active_bins(num_active_bins))
      active_bins = pack(cresp_all_bins, is_active_bin)

! Pure power law spectrum initial condition (default case)
      q = q_init
      f = zero
      f = f_amplitude * (p/p_lo_init)**(-q_init)

      if (initial_condition == "powl") call cresp_init_powl_spectrum(n, e, f_amplitude, q_init, p_lo_init, p_up_init)

      if (initial_condition == 'brpl') then
!>
!!/brief Power-law like spectrum with break at p_br_init
!! In this case initial spectrum with a break at p_min_fix is assumed, the initial slope
!! on the left side of the break is q_br_init. If initial_condition = "brpl", but parameters
!! are not defined in problem.par, "powl" spectrum is initialized with a warning issued.
!<
         i_br = minloc(abs(p_fix - p_br_init),dim=1)-1
         q(:i_br) = q_br_init ; q(i_br+1:) = q_init
         f(i_lo:i_br-1) = f(i_br) * (p(i_lo:i_br-1) / p(i_br)) ** (-q_br_init)
         e = fq_to_e(p(0:ncre-1), p(1:ncre), f(0:ncre-1), q(1:ncre), active_bins)
         n = fq_to_n(p(0:ncre-1), p(1:ncre), f(0:ncre-1), q(1:ncre), active_bins)
      endif

      if (initial_condition == 'symf') then
         i_br = int((i_lo+i_up)/2)
         q(i_lo+1:i_br) = -q_init
         f(i_br) = f(i_br+1)*(p(i_br+1)/p(i_br))**(-q(i_br))

         do i=1,i_br-i_lo
            f(i_br-i) = f(i_br+i)
         enddo

         if ((i_up - i_br .ne. i_br - i_lo))  p_up = p_up - (p_up - p_fix(i_up-1))
         p(i_up) = p_up ; i_up = i_up -1
         e = fq_to_e(p(0:ncre-1), p(1:ncre), f(0:ncre-1), q(1:ncre), active_bins)
         n = fq_to_n(p(0:ncre-1), p(1:ncre), f(0:ncre-1), q(1:ncre), active_bins)
     endif

     if (initial_condition == 'syme' ) then
         i_br = int((i_lo+i_up)/2)
         q(i_lo+1:i_br) = q_init-2.2
         do i=1,i_br-i_lo
             f(i_br-i) = f(i_br)*(p(i_br)/p(i_br-i))**(q(i_br-i+1))
         enddo
         if ((i_up - i_br .ne. i_br - i_lo))  p_up = p_up - (p_up - p_fix(i_up-1))
         p(i_up) = p_up ; i_up = i_up -1
         e = fq_to_e(p(0:ncre-1), p(1:ncre), f(0:ncre-1), q(1:ncre), active_bins)
         n = fq_to_n(p(0:ncre-1), p(1:ncre), f(0:ncre-1), q(1:ncre), active_bins)
     endif

     if (initial_condition == 'bump') then  ! TODO - @cresp_grid energy normalization and integral to scale cosmic ray electrons with nucleon energy density!
! Gaussian bump-type initial condition for energy distribution
         f = f_amplitude * exp(-(4*log(2.0)*log(p/sqrt(p_lo_init*p_up_init/1.))**2)) ! FWHM
         f(0:ncre-1) = f(0:ncre-1) / (fpi * clight * p(0:ncre-1)**(3.0)) ! without this spectrum is gaussian for distribution function
         do i=1, ncre
            q(i) = pf_to_q(p(i-1),p(i),f(i-1),f(i)) !-log(f(i)/f(i-1))/log(p(i)/p(i-1))
         enddo
         e = fq_to_e(p(0:ncre-1), p(1:ncre), f(0:ncre-1), q(1:ncre), active_bins)
         n = fq_to_n(p(0:ncre-1), p(1:ncre), f(0:ncre-1), q(1:ncre), active_bins)
      endif

      crel%p = p
      crel%f = f
      crel%q = q
      crel%e = e
      crel%n = n
      crel%i_lo = i_lo
      crel%i_up = i_up

      if ( e_small_approx_init_cond .gt. 0) then
         if ( (approx_p_lo + e_small_approx_init_cond) .gt. 0 ) then
            call get_fqp_lo(exit_code)
            if (exit_code) then
               write(msg,*) "[cresp_crspectrum:cresp_init_state] approx_p_lo = 1, but solution for initial spectrum lower cutoff not found, exiting! "
               call die(msg)
            endif
         endif

         if ( (approx_p_up + e_small_approx_init_cond) .gt. 0 ) then
            call get_fqp_up(exit_code)
            if (exit_code) then
               write(msg,*) "[cresp_crspectrum:cresp_init_state] approx_p_up = 1, but solution for initial spectrum upper cutoff not found, exiting! "
               call die(msg)
            endif
         endif

         if (allow_source_spectrum_break) then

            i_lo_ch = int(floor(log10(p_lo/p_fix(1))/w)) + 1
            i_lo_ch = max(0, i_lo_ch)
            i_lo_ch = min(i_lo_ch, ncre - 1)

            i_up_ch = int(floor(log10(p_up/p_fix(1))/w)) + 2
            i_up_ch = max(1,i_up_ch)
            i_up_ch = min(i_up_ch,ncre)

            f(i_up_ch) = e_small_to_f(p_up)
            q(i_up_ch) = q(i_up)
            p(i_up_ch) = p_up

            p(i_lo_ch) = p_lo
            f(i_lo_ch) = e_small_to_f(p_lo)
            q(i_lo_ch+1) = q(i_lo+1)

            do i=i_lo_ch+1, i_lo
               p(i) = p_fix(i)
               f(i) = f(i_lo_ch) * (p_fix(i)/p(i_lo_ch))**(-q(i_lo_ch+1))
               q(i+1) = q(i_lo_ch+1)
#ifdef CRESP_VERBOSED
               write (msg, "(A)") 'Extending the range of lower boundary bin after NR_2dim momentum search'   ; call printinfo(msg)
#endif /* CRESP_VERBOSED */
            enddo

            do i=i_up, i_up_ch-1
               p(i) = p_fix(i)
               f(i) = f(i_up-1)* (p_fix(i)/p_fix(i_up-1))**(-q(i_up))
               q(i) = q(i_up)
#ifdef CRESP_VERBOSED
               write (msg, "(A)") 'Extending the range of upper boundary bin after NR_2dim momentum search'   ; call printinfo(msg)
#endif /* CRESP_VERBOSED */
            enddo
#ifdef CRESP_VERBOSED
            write (msg,"(A,2I3,A,2I3)") "Boundary bins now (i_lo_new i_lo | i_up_new i_up)",  i_lo_ch, i_lo, ' |', i_up_ch, i_up     ; call printinfo(msg)
#endif /* CRESP_VERBOSED */
            i_lo = i_lo_ch   ;   i_up = i_up_ch
            q(i_up_ch) = q(i_up)
            p(i_up) = p_fix(i_up);  p(i_up) = p_up

            is_active_bin = .false.
            is_active_bin(i_lo+1:i_up) = .true.
            num_active_bins = count(is_active_bin) ! active arrays must be reevaluated - number of active bins and edges might have changed
            if (allocated(active_bins)) deallocate(active_bins)
            allocate(active_bins(num_active_bins)) ! active arrays must be reevaluated - number of active bins and edges might have changed
            active_bins = pack(cresp_all_bins, is_active_bin)

            e = fq_to_e(p(0:ncre-1), p(1:ncre), f(0:ncre-1), q(1:ncre), active_bins) ! once again we must count n and e
            n = fq_to_n(p(0:ncre-1), p(1:ncre), f(0:ncre-1), q(1:ncre), active_bins)

            crel%p = p
            crel%f = f
            crel%q = q
            crel%e = e
            crel%n = n
            crel%i_lo = i_lo
            crel%i_up = i_up
         endif
      endif

      n_tot0 = sum(n)
      e_tot0 = sum(e)

      init_n = n
      init_e = e

      total_init_cree = sum(e) !< total_init_cree value is used for initial spectrum scaling when spectrum is injected by source.
      call deallocate_active_arrays

   end subroutine cresp_init_state

!-------------------------------------------------------------------------------------------------
! Assumes power-law spectrum, without breaks. In principle the same thing is done in cresp_init_state, but
! init_state cannot be called from "outside".
!-------------------------------------------------------------------------------------------------
   subroutine cresp_init_powl_spectrum(n_inout, e_inout, f_in, q_in, p_dist_lo, p_dist_up)

      use constants,      only: zero
      use diagnostics,    only: my_deallocate
      use initcrspectrum, only: ncre, p_fix, w, cresp_all_bins, cresp_all_edges

      implicit none

      real(kind=8), dimension(1:ncre), intent(inout) :: n_inout, e_inout
      real(kind=8), intent(in) ::     f_in, q_in, p_dist_lo, p_dist_up
      real(kind=8), dimension(1:ncre) :: n_add, e_add, q_add
      real(kind=8), dimension(0:ncre) :: p_range_add , f_add
      integer(kind=4), allocatable, dimension(:) :: act_bins, act_edges
      integer(kind=4) :: i_l, i_u !, n_bins

      n_add = zero  ; e_add = zero  ; q_add = zero  ; f_add = zero  ; p_range_add = zero

      i_l = int(floor(log10(p_dist_lo/p_fix(1))/w)) + 1
      i_l = max(0, i_l)
      i_l = min(i_l, ncre - 1)

      i_u = int(floor(log10(p_dist_up/p_fix(1))/w)) + 2
      i_u = max(1,i_u)
      i_u = min(i_u,ncre)

      p_range_add(i_l:i_u) = p_fix(i_l:i_u)
      p_range_add(i_l) = p_dist_lo
      p_range_add(i_u) = p_dist_up
      if (.not.allocated(act_edges)) allocate(act_edges(i_u - i_l  ))
      if (.not.allocated(act_bins )) allocate( act_bins(i_u - i_l+1))
      act_edges =  cresp_all_edges(i_l  :i_u)
      act_bins  =   cresp_all_bins(i_l+1:i_u)
      q_add(act_bins) = q_in

      f_add(act_edges) = f_in * (p_range_add(act_edges)/p_dist_lo)**(-q_in)

      n_add = fq_to_n(p_range_add(0:ncre-1), p_range_add(1:ncre), f_add(0:ncre-1), q_add(1:ncre), act_bins)
      e_add = fq_to_e(p_range_add(0:ncre-1), p_range_add(1:ncre), f_add(0:ncre-1), q_add(1:ncre), act_bins)

      n_inout = n_inout + n_add
      e_inout = e_inout + e_add

      call my_deallocate(act_bins)
      call my_deallocate(act_edges)

   end subroutine cresp_init_powl_spectrum
!-------------------------------------------------------------------------------------------------
   subroutine cresp_get_scaled_init_spectrum(n_inout, e_inout, e_in_total) !< Using n,e spectrum obtained at initialization, obtain injected spectrum at given cell
      use initcrspectrum, only: norm_init_spectrum, total_init_cree, ncre
      implicit none
      real(kind=8), dimension(1:ncre), intent(inout) :: n_inout, e_inout
      real(kind=8), intent(in)                       :: e_in_total

      n_inout = norm_init_spectrum%n * e_in_total / total_init_cree
      e_inout = norm_init_spectrum%e * e_in_total / total_init_cree

   end subroutine cresp_get_scaled_init_spectrum

!-------------------------------------------------------------------------------------------------
!
! energy integral (eq. 21)
!
!-------------------------------------------------------------------------------------------------

   function fq_to_e(p_l, p_r, f_l, q, bins)
      use constants, only: zero, one, four, fpi
      use cresp_variables, only: clight ! use units, only: clight
      use initcrspectrum, only: ncre, eps

      implicit none

      real(kind=8), dimension(:), intent (in)  :: p_l, p_r, f_l, q
      integer, dimension(:), intent(in)        :: bins
      real(kind=8), dimension(size(bins))  :: e_bins
      real(kind=8), dimension(1:ncre)      :: fq_to_e

      fq_to_e = zero
      e_bins = fpi*clight*f_l(bins)*p_l(bins)**4
      where (abs(q(bins) - four) .gt. eps)
         e_bins = e_bins*((p_r(bins)/p_l(bins))**(four-q(bins)) - one)/(four - q(bins))
      elsewhere
         e_bins = e_bins*log(p_r(bins)/p_l(bins))
      endwhere

      fq_to_e(bins) = e_bins

   end function fq_to_e
!-------------------------------------------------------------------------------------------------
!
! Compute edge values (amplitudes) of e
!
!-------------------------------------------------------------------------------------------------

   function fp_to_e_ampl(p_1, f_1)
      use constants,          only: fpi
      use cresp_variables,    only: clight

      implicit none

      real(kind=8), intent(in)    :: p_1, f_1
      real(kind=8)                :: fp_to_e_ampl

      fp_to_e_ampl = fpi * clight**2 * f_1 * p_1**3

   end function fp_to_e_ampl

!-------------------------------------------------------------------------------------------------
!
! density integral (eq. 9)
!
!-------------------------------------------------------------------------------------------------

   function fq_to_n(p_l, p_r, f_l, q, bins)
      use constants, only: zero, one, three, fpi
      use initcrspectrum, only: ncre, eps

      implicit none

      integer, dimension(:), intent(in)      :: bins
      real(kind=8), dimension(:), intent(in) :: p_l, p_r, f_l, q
      real(kind=8), dimension(size(bins))   :: n_bins
      real(kind=8), dimension(1:ncre)       :: fq_to_n

      n_bins = zero

      n_bins = fpi*f_l(bins)*p_l(bins)**3
      where (abs(q(bins) - three) .gt. eps)
         n_bins = n_bins*((p_r(bins)/p_l(bins))**(three-q(bins)) - one)/(three - q(bins))
      elsewhere
         n_bins = n_bins*log((p_r(bins)/p_l(bins)))
      endwhere

      fq_to_n = zero
      fq_to_n(bins) = n_bins

   end function fq_to_n

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
      use constants, only: zero, one, three, four, fpi
      use cresp_variables, only: clight ! use units, only: clight
      use initcrspectrum, only: ncre, eps, cresp_all_bins

      implicit none

      integer, dimension(:), intent(in) :: ce, he    ! cooling edges, heating edges
      real(kind=8), dimension(1:ncre-1) :: pimh, pimth, fimh,fimth  ! *imh = i_minus_half, *imth = i_minus_third
      real(kind=8), dimension(1:ncre-1) :: dn_upw, de_upw, qi,qim1  ! *im1 = i_minus_one

      pimh(1:ncre-1) = p(1:ncre-1)
      pimth(1:ncre-1) = p(0:ncre-2)

      fimh(1:ncre-1) = f(1:ncre-1)
      fimth(1:ncre-1) = f(0:ncre-2)

      qi(1:ncre-1)  = q(2:ncre)
      qim1(1:ncre-1) = q(1:ncre-1)

      dn_upw = zero
      de_upw = zero
      nflux  = zero
      eflux  = zero

      dn_upw(ce) = fpi*fimh(ce)*pimh(ce)**3
      where (abs( qi(ce) - three ) .gt. eps)
         dn_upw(ce) = dn_upw(ce)*((p_upw(ce)/pimh(ce))**(three-qi(ce)) - one)/(three - qi(ce))
      elsewhere
         dn_upw(ce) = dn_upw(ce)*log((p_upw(ce)/pimh(ce)))
      endwhere
      nflux(ce) = - dn_upw(ce)

      de_upw(ce) = fpi*clight*fimh(ce)*pimh(ce)**4
      where (abs(qi(ce) - four) .gt. eps)
         de_upw(ce) = de_upw(ce)*((p_upw(ce)/pimh(ce))**(four-qi(ce)) - one)/(four - qi(ce))
      elsewhere
         de_upw(ce) = de_upw(ce)*log(p_upw(ce)/pimh(ce))
      endwhere
      eflux(ce) =  - de_upw(ce)

      if (del_i_up == -1) then
         nflux(i_up-1) = -n(i_up)
         eflux(i_up-1) = -e(i_up)
      endif

! filling empty empty bin - switch of upper boundary, condition is checked only once per flux computation and is very rarely satisfied.
      if (nflux(i_up) .gt. zero) then             ! If flux is greater than zero it will go through right edge, activating next bin in the next timestep.
         if ( cresp_all_bins(i_up+1) .eq. i_up+1 ) then  ! But it shuld only happen if there is bin with index i_up+1
            ndt(i_up+1) = nflux(i_up)
            edt(i_up+1) = eflux(i_up)
            del_i_up = +1
         endif
      endif

      if ( nflux(i_up-1)+n(i_up) .le. zero) then ! If flux is equal or greater than energy / density in a given bin,  these both shall migrate
         nflux(i_up-1) =  -n(i_up)                   ! to an adjacent bin, thus making given bin detected as inactive (empty) in the next timestep
         eflux(i_up-1) =  -e(i_up)
         del_i_up = -1
      endif

      dn_upw(he) = fpi*fimth(he)*p_upw(he)**3*(pimth(he)/p_upw(he))**qim1(he)
      where (abs(qim1(he) - three) .gt. eps )
         dn_upw(he) = dn_upw(he)*((pimh(he)/p_upw(he))**(three-qim1(he)) - one)/(three - qim1(he))
      elsewhere
         dn_upw(he) = dn_upw(he)*log((pimh(he)/p_upw(he)))
      endwhere
      nflux(he) = dn_upw(he)

      de_upw(he) = fpi*clight*fimth(he)*p_upw(he)**4*(pimth(he)/p_upw(he))**qim1(he)
      where (abs(qi(he) - four) .gt. eps)
         de_upw(he) = de_upw(he)*((pimh(he)/p_upw(he))**(four-qim1(he)) - one)/(four - qim1(he))
      elsewhere
         de_upw(he) = de_upw(he)*log(pimh(he)/p_upw(he))
      endwhere
      eflux(he) = de_upw(he)

      if (del_i_lo == 1 .or. nflux(i_lo+1) .ge. n(i_lo+1) ) then
         nflux(i_lo+1) = n(i_lo+1)
         eflux(i_lo+1) = e(i_lo+1)
! emptying lower boundary bins - in cases when flux gets greater than energy or number density
         del_i_lo = 1   ! in case it hasn't yet been modified
      endif

   end subroutine cresp_compute_fluxes

!-------------------------------------------------------------------------------------------------
!
! compute R (eq. 25)
!
!-------------------------------------------------------------------------------------------------
   subroutine cresp_compute_r(p, bins)

      use constants, only: zero, four, five
      use initcrspectrum, only: ncre, eps

      implicit none

      integer, dimension(:), intent(in)            :: bins
      real(kind=8), dimension(0:ncre), intent(in)  :: p
      real(kind=8), dimension(size(bins))          :: r_num, r_den

      r = zero

      where (abs(q(bins) - five) .gt. eps)
         r_num = (p(bins)**(five-q(bins)) - p(bins-1)**(five-q(bins)))/(five - q(bins))
      elsewhere
         r_num = log(p(bins)/p(bins-1))
      end where

      where (abs(q(bins) - four) .gt. eps)
         r_den = (p(bins)**(four-q(bins)) - p(bins-1)**(four-q(bins)))/(four - q(bins))
      else where
         r_den = log(p(bins)/p(bins-1))
      end where

      where ( abs(r_num) .gt. zero .and. abs(r_den) .gt. zero )                  !< BEWARE - regression: comparisons against
         r(bins) = u_d + u_b * r_num/r_den !all cooling effects will come here   !< eps and epsilon result in bad results;
      end where                                                                  !< range of values ofr_num and r_den is very wide

   end subroutine cresp_compute_r

!-------------------------------------------------------------------------------------------------
!
! find new q (eq. 29) and new f
!
!-------------------------------------------------------------------------------------------------

   subroutine ne_to_q(n, e, q, bins)

      use constants, only: zero
      use cresp_NR_method, only: compute_q
      use cresp_variables, only: clight ! use units, only: clight
      use initcrspectrum, only: ncre, e_small

      implicit none

      integer(kind=4), dimension(:), intent(in)    :: bins
      integer          :: i, i_active
      real(kind=8)     :: alpha_in
      real(kind=8), dimension(1:ncre), intent(in)  :: n, e !
      real(kind=8), dimension(1:ncre), intent(out) :: q
      logical :: exit_code

      q = zero

      do i_active = 1 + approx_p_lo, size(bins) - approx_p_up
         i = bins(i_active)
         if (e(i) .gt. e_small .and. p(i-1) .gt. zero) then
            exit_code = .true.
            alpha_in = e(i)/(n(i)*p(i-1)*clight)
            if ((i .eq. i_lo+1) .or. (i .eq. i_up)) then ! for boudary case, when momenta are not approximated
               q(i) = compute_q(alpha_in, exit_code, p(i)/p(i-1))
            else
               q(i) = compute_q(alpha_in, exit_code)
            endif
         else
            q(i) = zero
        endif
        if ( exit_code .eqv. .true. ) fail_count_comp_q(i) = fail_count_comp_q(i) + 1
      enddo

   end subroutine ne_to_q

!-------------------------------------------------------------------------------------------------
! Function used to obtain q for one cell out of f and p values - used to compute q after finding p_up
!-------------------------------------------------------------------------------------------------
   function pf_to_q(p_l, p_r, f_l, f_r)

      implicit none
      real(kind=8), intent(in)   :: p_l, p_r, f_l, f_r
      real(kind=8)               :: pf_to_q

      pf_to_q = 0.0
      pf_to_q = -log(f_r/f_l)/log(p_r/p_l) ! append value of q for given p_up

   end function pf_to_q

! -------------------------------------------------------------------------------------------------
!
! distribution function amplitudes (eq. 9)
!
! -------------------------------------------------------------------------------------------------
   function nq_to_f(p_l, p_r, n, q, bins)

      use initcrspectrum, only: ncre, eps
      use constants, only: zero, one, three, fpi

      implicit none

      integer, dimension(:)                 :: bins
      real(kind=8), dimension(1:ncre)       :: p_l, p_r, n, q
      real(kind=8), dimension(size(bins))   :: f_bins
      real(kind=8), dimension(0:ncre)       :: nq_to_f
      real(kind=8), dimension(0:ncre)   :: pr_by_pl   ! the array of values of p_r/p_l to avoid FPEs

      nq_to_f= zero
      f_bins = zero
      pr_by_pl = one

      where (p_r(bins).gt. zero .and. p_l(bins) .gt. zero ) ! p(i) = 0 in inactive bins. This condition should be met by providing proper "bins" range - FIXME
         pr_by_pl(bins) = p_r(bins) / p_l(bins)                     ! + comparing reals with zero is still risky
         f_bins = n(bins) / (fpi*p_l(bins)**3)
         where (abs(q(bins)-three) .gt. eps)
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
   function get_pcresp(p_l, p_r, f_l, q, bins) ! computes cre pressure, not used currently

      use constants,       only: zero, one, three, four, fpi
      use cresp_variables, only: clight ! use units, only: clight
      use initcrspectrum,  only: eps

      implicit none

      real(kind=8), dimension(:), intent (in)  :: p_l, p_r, f_l, q
      integer, dimension(:), intent(in)        :: bins
      real(kind=8), dimension(size(bins))  :: p_cresp
      real(kind=8)              :: get_pcresp

      get_pcresp = zero
      p_cresp = (fpi/three) * clight*f_l(bins)*p_l(bins)**4

      where (abs(q(bins) - four) .gt. eps)
         p_cresp = p_cresp*((p_r(bins)/p_l(bins))**(four-q(bins)) - one)/(four - q(bins))
      elsewhere
         p_cresp = p_cresp*log(p_r(bins)/p_l(bins))
      endwhere

      get_pcresp = sum(p_cresp)

   end function get_pcresp
!---------------------------------------------------------------------------------------------------
! Computing cosmic ray spectrum component pressure gradient
!---------------------------------------------------------------------------------------------------
   subroutine src_gpcresp(u, n, dx, grad_pcresp)

      use constants, only: onet
      use initcrspectrum, only: ncre, cre_active, cre_gpcr_ess

      implicit none

      integer(kind=4), intent(in) :: n
      real(kind=8), intent(in)                       :: dx
      real(kind=8), dimension(n, 1:ncre), intent(in) :: u
      real(kind=8), dimension(n), intent(out)        :: grad_pcresp
      real(kind=8), dimension(n)                     :: P_cresp_r, P_cresp_l

      grad_pcresp = 0.0 ;  P_cresp_l = 0.0 ;  P_cresp_r = 0.0

      if (cre_gpcr_ess) then
!        (ultrarelativistic)
         P_cresp_l(1:n-2) = onet * sum(u(1:n-2, :),dim=2)
         P_cresp_r(3:n)   = onet * sum(u(3:n,   :),dim=2)
         grad_pcresp(2:n-1) = cre_active * (P_cresp_l(1:n-2) - P_cresp_r(3:n) )/(2.*dx)
      endif
   end subroutine src_gpcresp
!---------------------------------------------------------------------------------------------------
! Preparation and computatuon of upper boundary momentum "p_up" and and upper boundary
! distribution function value on left bin edge "f"
!---------------------------------------------------------------------------------------------------
   subroutine get_fqp_up(exit_code)
      use constants,       only: zero, one
      use cresp_variables, only: clight ! use units, only: clight
      use cresp_NR_method, only: intpol_pf_from_NR_grids, alpha, n_in, selected_function_2D, fvec_up, &
                           &     NR_algorithm, e_small_to_f, q_ratios, assoc_pointers_up
      use initcrspectrum,  only: e_small, q_big, p_fix, NR_refine_pf_up

      implicit none

      real(kind=8), dimension(1:2) :: x_NR, x_NR_init
      logical :: exit_code, interpolated

      x_NR = zero
      alpha = (e(i_up)/(n(i_up)*clight*p_fix(i_up-1)))
      n_in  = n(i_up)

      call assoc_pointers_up

      x_NR = intpol_pf_from_NR_grids(alpha, n_in, interpolated)
      if (interpolated .eqv. .false.) then
         exit_code = .true.
         fail_count_interpol(2) = fail_count_interpol(2) +1
         return
      else
         x_NR_init = x_NR
      endif

      exit_code = .not.interpolated
#ifdef CRESP_VERBOSED
      write(msg,"(A31,2E22.15)") "Input ratios(p, f) for NR (up):", x_NR      ; call printinfo(msg)
#endif /* CRESP_VERBOSED */

      if ( (NR_refine_pf_up .eqv. .true.) .or. (interpolated .eqv. .false.)) then
         call NR_algorithm(x_NR, exit_code)
         if (exit_code .eqv. .true.) then ! some failures still take place
            if (interpolated .eqv. .false.) then
               exit_code = .true.
#ifdef CRESP_VERBOSED
               write(msg,"(A,4E18.9)") " Interpolation AND NR failure (up)", alpha, n_in, x_NR_init      ; call printinfo(msg)
#endif /* CRESP_VERBOSED */
               return
            endif
            fail_count_NR_2dim(2) = fail_count_NR_2dim(2) +1
            x_NR = x_NR_init
         endif
      endif

      x_NR = abs(x_NR) ! negative values cannot be allowed

      if (x_NR(1) .lt. 1.0) then
         exit_code = .true.
         return
      endif

      p_up      = p_fix(i_up-1)*x_NR(1)
      p(i_up)   = p_up
      f(i_up-1) = e_small_to_f(p_up)/x_NR(2)
      q(i_up)   = q_ratios(x_NR(2), x_NR(1))

      if (abs(q(i_up)) .gt. q_big ) q(i_up) = sign(one, q(i_up)) * q_big
#ifdef CRESP_VERBOSED
      write(msg, "(A26,2E22.15)") " >>> Obtained (p_up, f_l):", p_up, f(i_up-1)     ; call printinfo(msg)
      write(msg, "(A26,2E22.15)") "     Corresponding ratios:", x_NR(1), x_NR(2)    ; call printinfo(msg)
#endif /* CRESP_VERBOSED */
      alpha = zero ;  n_in = zero

   end subroutine get_fqp_up
!--------------------------------------------------------------------------------------------------
! Preparation and computation of upper boundary momentum "p_lo" and and upper boundary
! distribution function value on the right bin edge "f"
!--------------------------------------------------------------------------------------------------
   subroutine get_fqp_lo(exit_code)

      use constants, only: zero, one
      use cresp_NR_method, only: intpol_pf_from_NR_grids, alpha, n_in, selected_function_2D, fvec_lo, &
                  NR_algorithm, e_small_to_f, q_ratios, assoc_pointers_lo
      use cresp_variables, only: clight ! use units, only: clight
      use initcrspectrum, only: e_small, q_big, p_fix, NR_refine_pf_lo

      implicit none

      real(kind=8), dimension(1:2) :: x_NR, x_NR_init
      logical :: exit_code, interpolated

      x_NR = zero
      alpha = (e(i_lo+1)/(n(i_lo+1)*clight*p_fix(i_lo+1)))
      n_in  = n(i_lo+1)

      call assoc_pointers_lo

      x_NR = intpol_pf_from_NR_grids(alpha, n_in, interpolated)
      if (interpolated .eqv. .false.) then
         exit_code = .true.
         fail_count_interpol(1) = fail_count_interpol(1) +1
         return
      else
         x_NR_init = x_NR
      endif

      exit_code = .not.interpolated !
#ifdef CRESP_VERBOSED
      write (msg, "(A31,2E22.15)") "Input ratios(p, f) for NR (lo):", x_NR    ; call printinfo(msg)
#endif /* CRESP_VERBOSED */
      if ( (NR_refine_pf_lo .eqv. .true.) .or. (interpolated .eqv. .false.)) then
         call NR_algorithm(x_NR, exit_code)
         if (exit_code .eqv. .true.) then ! some failures still take place
            if (interpolated .eqv. .false.) then
               exit_code = .true.
#ifdef CRESP_VERBOSED
               write (msg, "(A,3E18.9)") " Interpolation AND NR failure (lo)", alpha, n_in     ; call printinfo(msg)
#endif /* CRESP_VERBOSED */
               return
            endif
            fail_count_NR_2dim(1) = fail_count_NR_2dim(1) +1
            x_NR = x_NR_init
         endif
      endif

      x_NR = abs(x_NR) ! negative values cannot be allowed
      if (x_NR(1) .lt. 1.0) then
         exit_code = .true.
         return
      endif

      p_lo      = p_fix(i_lo+1)/ x_NR(1)
      p(i_lo)   = p_lo
      f(i_lo)   = e_small_to_f(p_lo)
      q(i_lo+1) = q_ratios(x_NR(2), x_NR(1))

      if (abs(q(i_lo+1)) .gt. q_big ) q(i_lo+1) = sign(one, q(i_lo+1)) * q_big
#ifdef CRESP_VERBOSED
      write (msg, "(A26,2E22.15)") " >>> Obtained (p_lo, f_r):", p_lo, x_NR(2)*f(i_lo)    ; call printinfo(msg)
      write (msg, "(A26,2E22.15)") "     Corresponding ratios:", x_NR(1), x_NR(2)         ; call printinfo(msg)
      call printinfo(msg)
#endif /* CRESP_VERBOSED */
      alpha = zero ;  n_in = zero

   end subroutine get_fqp_lo
!----------------------------------------------------------------------------------------------------
   function b_losses(p)
      implicit none

      real(kind=8), dimension(:), intent(in)  :: p
      real(kind=8), dimension(size(p)) :: b_losses

      b_losses = u_b*p**2  !!! b_sync_ic = 8.94e-25*(u_b+u_cmb)*gamma_l**2 ! erg/cm

   end function b_losses
!-------------------------------------------------------------------------------------------------
!>
!! \brief Relative change of momentum due to losses (u_b*p*dt) and compression u_d*dt (Taylor expansion up to 3rd order)
!<
!====================================================================================================
   function p_rch_ord_1(dt, p)
      use constants,    only: half, sixth

      implicit none

      real(kind=8), intent(in)  :: dt, p
      real(kind=8)              :: p_rch_ord_1

      p_rch_ord_1 = -( u_d + p * u_b ) *  dt

   end function p_rch_ord_1
!-------------------------------------------------------------------------------------------------
   function p_rch_ord_2_1(dt, p)     !< adds 2nd term and calls 1st order
      use constants,    only: half, sixth

      implicit none

      real(kind=8), intent(in)  :: dt, p
      real(kind=8)              :: p_rch_ord_2_1

      p_rch_ord_2_1 = p_rch_ord_1(dt, p)  + ( half*(u_d*dt)**2 + (u_b * p * dt)**2)

   end function p_rch_ord_2_1
!-------------------------------------------------------------------------------------------------
   function p_rch_ord_3_2_1(dt, p)     !< adds 3rd term and calls 2nd and 1st order
      use constants,    only: half, sixth

      implicit none

      real(kind=8), intent(in)  :: dt, p
      real(kind=8)              :: p_rch_ord_3_2_1

      p_rch_ord_3_2_1 = p_rch_ord_2_1(dt, p) - sixth*(u_d * dt)**3 - (u_b * p * dt)**3

   end function p_rch_ord_3_2_1
!----------------------------------------------------------------------------------------------------
   subroutine p_rch_init
      use initcrspectrum,  only: expan_order

      implicit none

      select case(expan_order)
         case(1)
            p_rch => p_rch_ord_1
         case(2)
            p_rch => p_rch_ord_2_1
         case(3)
            p_rch => p_rch_ord_3_2_1
         case default
            write (msg,*) '[cresp_crspectrum:p_rch_init] expan_order =',expan_order,': value incorrect (accepted values [1;3]).'
            call die(msg)
      end select
   end subroutine p_rch_init
!====================================================================================================
   subroutine transfer_quantities(give_to, take_from)
      use constants, only: zero

      implicit none

      real(kind=8), intent(inout) :: give_to, take_from

      give_to = give_to + take_from
      take_from = zero

   end subroutine transfer_quantities
!----------------------------------------------------------------------------------------------------
   subroutine cresp_allocate_all
      use diagnostics, only: my_allocate_with_index
      use initcrspectrum, only: ncre

      implicit none

      integer(kind = 4)          :: ma1d

      ma1d = ncre
      call my_allocate_with_index(fail_count_comp_q,ma1d,1)

      call my_allocate_with_index(n,ma1d,1)   !:: n, e, r
      call my_allocate_with_index(e,ma1d,1)
      call my_allocate_with_index(e_amplitudes_l,ma1d,1)   !:: n, e, r
      call my_allocate_with_index(e_amplitudes_r,ma1d,1)
      call my_allocate_with_index(r,ma1d,1)
      call my_allocate_with_index(q,ma1d,1)

      call my_allocate_with_index(f,ma1d,0)
      call my_allocate_with_index(p,ma1d,0)

      call my_allocate_with_index(edt,ma1d,1)
      call my_allocate_with_index(ndt,ma1d,1)

      call my_allocate_with_index(p_next,ma1d,0)
      call my_allocate_with_index(p_upw,ma1d,0)
      call my_allocate_with_index(nflux,ma1d,0)
      call my_allocate_with_index(eflux,ma1d,0)

      call my_allocate_with_index(is_fixed_edge,ma1d,0)
      call my_allocate_with_index(is_fixed_edge_next,ma1d,0)
      call my_allocate_with_index(is_active_edge,ma1d,0)
      call my_allocate_with_index(is_active_edge_next,ma1d,0)
      call my_allocate_with_index(is_cooling_edge,ma1d,0)
      call my_allocate_with_index(is_cooling_edge_next,ma1d,0)
      call my_allocate_with_index(is_heating_edge,ma1d,0)
      call my_allocate_with_index(is_heating_edge_next,ma1d,0)
      call my_allocate_with_index(is_active_bin,ma1d,1)
      call my_allocate_with_index(is_active_bin_next,ma1d,1)

   end subroutine cresp_allocate_all
!----------------------------------------------------------------------------------------------------
   subroutine cresp_deallocate_all
      use diagnostics, only: my_deallocate
      implicit none

      call my_deallocate(n)
      call my_deallocate(e)
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

!---------------------------------------------------------------------------------------------
   subroutine cresp_accuracy_test(t) ! Not for use with PIERNIK
      implicit none

      real(kind=8), intent(in)   :: t

      print*, 'Accuracy test for adabatic compression/expansion:'
      print*, 'n_tot = ', n_tot, 'n_tot0 = ', n_tot0, 'rel error = ', (n_tot - n_tot0)/n_tot0
      print*, 'e_tot = ', e_tot, 'e_anal = ', e_tot0*exp(-u_d_0*t), 'rel error = ',(e_tot-e_tot0*exp(-u_d_0*t))/(e_tot0*exp(-u_d_0*t))
      print*, 'e_tot0= ', e_tot0
      print *, '=================================='
      print*,  '! End of iteration               !'
      print *, '=================================='
      print*
      print*
      print*,'--------------------'
      print *,''

   end subroutine cresp_accuracy_test
!----------------------------------------------------------------------------------------------------
   subroutine cleanup_cresp
      implicit none

      write (msg, '(A36,I6,A6,I6)') "NR_2dim:  convergence failure: p_lo", fail_count_NR_2dim(1), ", p_up", fail_count_NR_2dim(2)   ; call printinfo(msg)
      write (msg, '(A36,I6,A6,I6)') "NR_2dim:interpolation failure: p_lo", fail_count_interpol(1), ", p_up", fail_count_interpol(2) ; call printinfo(msg)
      write (msg, '(A36,   100I5)') "NR_2dim:inpl/solve  q(bin) failure:", fail_count_comp_q    ; call printinfo(msg)
      call cresp_deallocate_all

   end subroutine cleanup_cresp
!----------------------------------------------------------------------------------------------------
   subroutine printer(t)
      use initcrspectrum, only: ncre, crel

      implicit none

      real(kind = 8)   :: t

      open(10, file="crs.dat", position='append')
      write(10, '(2e16.9, 3(1x,i8), 600(1x,ES18.9E3))') t, crel%dt, ncre, crel%i_lo, crel%i_up, crel%p, crel%f, crel%q
      close(10)

      open(11, file="crs_ne.dat", position='append')
      write(11, '(2I5,4x, e16.9, 600(1x,F18.9))') del_i_lo, del_i_up, t, crel%dt, crel%p(i_lo), crel%p(i_up), crel%n, crel%e
      close(11)

   end subroutine printer

!----------------------------------------------------------------------------------------------------

end module cresp_crspectrum
