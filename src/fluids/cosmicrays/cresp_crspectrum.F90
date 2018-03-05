module cresp_crspectrum
! pulled by COSM_RAY_ELECTRONS
 implicit none
  public :: cresp_update_cell, cresp_init_state, printer, fail_count_interpol, fail_count_no_sol, fail_count_NR_2dim, &
      &   cleanup_cresp, cresp_accuracy_test, b_losses, cresp_allocate_all, cresp_deallocate_all, e_threshold_lo, e_threshold_up, &
      &   fail_count_comp_q, second_fail, src_gpcresp, cresp_init_powl_spectrum, get_powl_f_ampl, e_tot_2_f_init_params, e_tot_2_en_powl_init_params
  private ! most of it
  real(kind=8)     , parameter      :: three   = 3.e0
  real(kind=8)     , parameter      :: four    = 4.e0
  real(kind=8)     , parameter      :: five    = 5.e0
  real(kind=8)     , parameter      :: ten     = 10.e0

  integer, dimension(1:2), save     :: fail_count_NR_2dim, fail_count_interpol, fail_count_no_sol, second_fail
  integer, allocatable, save   :: fail_count_comp_q(:)

! variables informing about change of bins
  integer                           :: del_i_lo, del_i_up

! logical arrays / arrays determining use of p/n/e in some lines
   logical, allocatable, dimension(:) :: is_fixed_edge,   is_fixed_edge_next 
   logical, allocatable, dimension(:) :: is_active_edge,  is_active_edge_next 
   logical, allocatable, dimension(:) :: is_cooling_edge, is_cooling_edge_next 
   logical, allocatable, dimension(:) :: is_heating_edge, is_heating_edge_next 
   logical, allocatable, dimension(:) :: is_active_bin,   is_active_bin_next
   logical, allocatable, dimension(:) :: not_spectrum_break

! counters
  integer                           :: num_fixed_edges,   num_fixed_edges_next
  integer                           :: num_active_edges,  num_active_edges_next
  integer                           :: num_active_bins,   num_active_bins_next
  integer                           :: num_cooling_edges_next
  integer                           :: num_heating_edges_next
!   integer                           :: num_not_spectrum_break

! dynamic arrays
  integer, allocatable              :: fixed_edges(:),   fixed_edges_next(:)
  integer, allocatable              :: active_edges(:),  active_edges_next(:)
  integer, allocatable              :: active_bins(:),   active_bins_next(:)
  integer, allocatable              :: cooling_edges_next(:)
  integer, allocatable              :: heating_edges_next(:)

   real(kind=8), allocatable, dimension(:) :: r  ! r term for energy losses (Miniati 2001, eqn. 25)
   real(kind=8), allocatable, dimension(:) :: q  ! power-law exponent array

! power-law 
   real(kind=8)                      :: p_lo_next, p_up_next, p_lo, p_up !, p_lo_bef, p_up_bef
   integer                           :: i_lo
   integer                           :: i_up
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
! virtual e,n arrays for cutoff, for cases when bins are only slightly filled and p ~ p_fix - it might not be possible to find solution via NR  
  real(kind=8), dimension(1:2) :: vrtl_n, vrtl_e
! lower / upper energy needed for bin activation
  real(kind=8), save :: e_threshold_lo, e_threshold_up
! if one bin, switch off cutoff p approximation
  integer  :: approx_p_lo, approx_p_up
!-------------------------------------------------------------------------------------------------
!
contains
!
!-------------------------------------------------------------------------------------------------

!----- main subroutine -----

  subroutine cresp_update_cell(dt, n_inout, e_inout, sptab, v_n, v_e, cfl_cresp_violation, p_out) !, p_lo_cell, p_up_cell)
   use initcrspectrum, only: ncre, spec_mod_trms, e_small_approx_p_lo, e_small_approx_p_up, e_small_approx_init_cond, crel, p_mid_fix
! #ifdef VERBOSE
   use initcrspectrum, only: p_fix
! #endif /* VERBOSE */
   use constants, only: zero, one
   implicit none
    real(kind=8), dimension(1:2), intent(inout), optional :: p_out
    real(kind=8), dimension(1:2), intent(inout) :: v_n, v_e
    logical, intent(inout)                      :: cfl_cresp_violation
    real(kind=8), intent(in)  :: dt
    real(kind=8), dimension(1:ncre), intent(inout)   :: n_inout, e_inout
    type(spec_mod_trms), intent(in)       :: sptab
    logical :: solve_fail_lo, solve_fail_up, index_changed, empty_cell
        e = zero; n = zero; edt = zero; ndt = zero
        solve_fail_lo = .false.
        solve_fail_up = .false.
        index_changed = .false.
        empty_cell    = .false.
        cfl_cresp_violation = .false.
        
        approx_p_lo = e_small_approx_p_lo
        approx_p_up = e_small_approx_p_up

        p_lo_next = zero
        p_up_next = zero
        if (present(p_out)) then
            p_lo = p_out(1)
            p_up = p_out(2)
        endif
        
        r = zero
        f = zero
        q = zero

        n = n_inout     ! number density of electrons passed to cresp module by the external module / grid
        e = e_inout     ! energy density of electrons passed to cresp module by the external module / grid
        u_b = sptab%ub
        u_d = sptab%ud

        vrtl_e = v_e
        vrtl_n = v_n

        call find_i_bound(empty_cell)
        if ( empty_cell ) return             ! if grid cell contains empty bins, no action is taken
        call cresp_find_active_bins
        call cresp_organize_p

! Compute power indexes for each bin at [t] and f on left bin faces at [t] 
        f = zero; q=zero
        call ne_to_q(n, e, q, active_bins)

! Here values of distribution function f for active left edges (excluding upper momentum boundary) are computed
        f = nq_to_f(p(0:ncre-1), p(1:ncre), n(1:ncre), q(1:ncre), active_bins)
        
        if (approx_p_up .gt. 0) then         ! momenta values stored only within module - for tests; will not work in PIERNIK
            if (i_up .gt. 1) then
                call get_fqp_up(solve_fail_up)
            else                                                   ! spectrum beyond the fixed momentum grid
                p_up = p_fix(i_up)
                p(i_up) = p_fix(i_up)
                solve_fail_up = .false.
            endif
        endif
        if (approx_p_lo .gt. 0) then  ! momenta values stored only within module - for tests; will not work in PIERNIK
            if ((i_lo+1) .ne. ncre) then
                call get_fqp_lo(solve_fail_lo)
            else                                                   ! spectrum beyond the fixed momentum grid
                p_lo = p_fix(i_lo)
                p(i_lo) = p_fix(i_lo)
                solve_fail_lo = .false.
            endif                                                  ! countermeasure against failure in finding boundary momentum
        endif
        if ((solve_fail_lo .eqv. .true.) .or. (solve_fail_up .eqv. .true.)) then
            call deallocate_active_arrays
            if (solve_fail_lo) then
                if (i_lo .gt. 0 .and. i_lo .lt. ncre-3) then
                    call transfer_quantities(e(i_lo+2),e(i_lo+1))
                    call transfer_quantities(n(i_lo+2),n(i_lo+1))
                    i_lo = i_lo + 1
                    call get_fqp_lo(solve_fail_lo)
                else
                    call transfer_quantities(vrtl_e(1),e(i_lo+1))
                    call transfer_quantities(vrtl_n(1),n(i_lo+1))
                    i_lo = i_lo+1
                    p_lo = p_fix(i_lo) ; p(i_lo) = p_fix(i_lo)
                endif
            endif
            if (solve_fail_up) then
                if (i_up .lt. ncre) then
                    call transfer_quantities(e(i_up-1),e(i_up)) ! instead of moving quantities to virtual
                    call transfer_quantities(n(i_up-1),n(i_up)) ! instead of moving quantities to virtual
                    i_up = i_up - 1
                    if (i_up .gt. 1) then
                        call get_fqp_up(solve_fail_up)
                    else
                        p_up = p_fix(i_up); p(i_up) = p_fix(i_up)
                    endif
                else
                    call transfer_quantities(vrtl_e(2),e(i_up))
                    call transfer_quantities(vrtl_n(2),n(i_up))
                endif
            endif
            if (i_up .eq. i_lo) then
                if (solve_fail_lo .and. solve_fail_up .eqv. .false.) then
                    i_lo = i_up - 1 ; p_lo = p_fix(i_up) ; p(i_lo) = p_fix(i_up)
                else if (solve_fail_up .and. solve_fail_lo .eqv. .false.) then
                    i_up = i_lo + 1 ; p_up = p_fix(i_lo) ; p(i_lo) = p_fix(i_lo)
                else
                    i_lo = i_up +1 ; p_lo = p_fix(i_lo) ; p(i_lo) = p_fix(i_lo) ; p_up = p_fix(i_up) ; p(i_up) = p_fix(i_up)
                endif
            endif
            if (p_lo .le. zero) then ! in case of second failure
                p_lo = p_mid_fix(i_lo)
                p(i_lo) = p_lo
                second_fail(1) = second_fail(1)+1
            endif
            if (p_up .le. zero) then ! in case of second failure
                p_up = p_mid_fix(i_up)
                p(i_up) = p_up
                second_fail(2) = second_fail(2)+1
            endif
            call cresp_find_active_bins
        endif

        call cresp_update_bin_index(dt, p_lo, p_up, p_lo_next, p_up_next, cfl_cresp_violation)      ! FIXME - must be modified in the future if this branch is connected to Piernik
        if ( cfl_cresp_violation ) then
            call deallocate_active_arrays
            return
        endif
! Compute fluxes through fixed edges in time period [t,t+dt], using f, q, p_lo and p_up at [t]
! Note that new [t+dt] values of p_lo and p_up in case new fixed edges appear or disappear.
! fill new bins
        call cresp_compute_fluxes(cooling_edges_next,heating_edges_next)

! Computing e and n at [t+dt]
        ndt(1:ncre) = n(1:ncre)  - (nflux(1:ncre) - nflux(0:ncre-1))
        edt(1:ncre) = e(1:ncre)  - (eflux(1:ncre) - eflux(0:ncre-1))

!         call boundary_flux_check ! If relative error between boundary momenta and p_fix is not tolerable, boundary fluxes are moved to "virtual" bins, preventing premature activation of new bins and associated numerical errors.
   
! edt(1:ncre) = e(1:ncre) *(one-0.5*dt*r(1:ncre)) - (eflux(1:ncre) - eflux(0:ncre-1))/(one+0.5*dt*r(1:ncre))   !!! oryginalnie u Miniatiego
! Compute coefficients R_i needed to find energy in [t,t+dt]
        call cresp_compute_r(p_next, active_bins_next)                 ! new active bins already received some particles, Ri is needed for those bins too

        edt(1:ncre) = edt(1:ncre) *(one-dt*r(1:ncre))

        p_lo = p_lo_next
        p_up = p_up_next
   
#ifdef VERBOSE
        write (*,'(A5, 50E18.9)') "p_fix", p_fix
        write (*,'(A5, 50E18.9)') "p_act", p
        write (*,'(A5, 50E18.9)') "p_nex", p_next
        write (*,'(A5, 50E18.9)') "p_upw", p_upw
        print '(A6, 1EN22.9, A9, 1EN22.9)', "p_lo ", p_lo, ",  p_up ", p_up

        print *, " "
        write (*,'(A5, 50E18.9)') "    n", n
        write (*,'(A5, 50E18.9)') "nflux", nflux
        write (*,'(A5, 50E18.9)') "  ndt", ndt

        print *, " "
        write (*,'(A5, 50E18.9)') "    e", e
        write (*,'(A5, 50E18.9)') "eflux", eflux
        write (*,'(A5, 50E18.9)') "  edt", edt

        print *, " "
        write (*,'(A5, 50E18.9)') "    r", r
        write (*,'(A5, 50E18.9)') "    q", q
        write (*,'(A5, 50E18.9)') "    f", f
        write (*,'(A15,2E18.9,A3,2E18.9)') "virtual e & n", vrtl_e, " | ", vrtl_n
        print *, " "
        if ( (approx_p_lo+approx_p_up) .gt. 0 ) then
            print '(A36,I5,A6,I3)', "NR_2dim:  convergence failure: p_lo", fail_count_NR_2dim(1), ", p_up", fail_count_NR_2dim(2)
            print '(A36,I5,A6,I3)', "NR_2dim:interpolation failure: p_lo", fail_count_interpol(1), ", p_up", fail_count_interpol(2)
            print '(A36,I5,A6,I3)', "NR_2dim:  no solution failure: p_lo", fail_count_no_sol(1), ", p_up", fail_count_no_sol(2)
            print '(A36,   100I5)', "NR_2dim:inpl/solve  q(bin) failure:", fail_count_comp_q
        endif
#endif /* VERBOSE */
        n = ndt
        e = edt
        call cresp_detect_negative_content ! for testing
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
        v_e = vrtl_e
        v_n = vrtl_n

        if (present(p_out)) then
            p_out(1) = p_lo
            p_out(2) = p_up
        endif
        call deallocate_active_arrays
  end subroutine cresp_update_cell
!-------------------------------------------------------------------------------------------------
! all the procedures below are called by cresp_update_cell subroutine or the driver
!-------------------------------------------------------------------------------------------------
  subroutine find_i_bound(exit_code)
  use initcrspectrum, only: ncre, e_small
  use constants, only: zero
  implicit none
    integer(kind=4) :: i
    logical :: i_lo_changed, i_up_changed, exit_code
        i_lo_changed = .false.
        i_up_changed = .false.
        exit_code    = .true.
! ! Locate cut-ofs before current timestep: indices are found without use of p_lo nor p_up and point to boundary edges    
        i_lo = 0
        do i = 1, ncre                        ! if energy density is nonzero, so should be the number density
            i_lo = i-1
            if ( e(i) .gt. e_threshold_lo ) then
                exit_code = .false.
                exit
            endif
        enddo

        if ( exit_code .eqv. .true.) return   ! empty cell - nothing to do here!

        i_up = ncre
        do i = ncre, 1,-1
            i_up = i
            if (e(i) .gt. e_threshold_up ) exit         ! if energy density is nonzero, so should be the number density
        enddo
        if ((e(i_lo+1)+vrtl_e(1)) .gt. e_small .and. vrtl_e(1) .gt. zero) then
            call transfer_quantities(e(i_lo+1), vrtl_e(1))
            call transfer_quantities(n(i_lo+1), vrtl_n(1))
        endif
        if ((e(i_up)+vrtl_e(2)) .gt. e_small .and. vrtl_e(2) .gt. zero) then
            call transfer_quantities(e(i_up), vrtl_e(2))
            call transfer_quantities(n(i_up), vrtl_n(2))
        endif
        if ( (approx_p_lo .eq. 1) .or. (approx_p_up .eq. 1) ) then ! TODO - this might need slight change of condition
            call threshold_energy_check_lo(e, n, i_lo_changed, .false.)
            call threshold_energy_check_up(e, n, i_up_changed, .false.)
            if (i_lo_changed) i_lo = i_lo + 1
            if (i_up_changed) i_up = i_up - 1
        endif
  end subroutine find_i_bound
!-------------------------------------------------------------------------------------------------
  subroutine cresp_find_active_bins
   use constants,      only: I_ZERO, zero
   use diagnostics,    only: my_allocate_with_index
   use initcrspectrum, only: ncre, e_small, cresp_all_edges, cresp_all_bins
   implicit none
      if(allocated(active_bins))  deallocate(active_bins)      
      if(allocated(active_edges)) deallocate(active_edges)

      is_active_bin = .false.
      is_active_bin(i_lo+1:i_up) = .true. ! unused if we use "not_spectrum_break"

      num_active_bins = count(is_active_bin)
      if (num_active_bins .gt. I_ZERO) then
        allocate(active_bins(num_active_bins))
        not_spectrum_break(:) = .false.
        where (e .gt. e_small .and. n .gt. zero)
            not_spectrum_break = .true.
        endwhere
        active_bins = I_ZERO
!         active_bins = pack(all_bins, is_active_bin)
        active_bins = pack(cresp_all_bins, not_spectrum_break)   ! not to iterate over spectrum break

! Construct index arrays for fixed edges betwen p_lo and p_up, active edges 
! before timestep  
        is_fixed_edge = .false.
        is_fixed_edge(i_lo+1:i_up-1) = .true.
        num_fixed_edges = count(is_fixed_edge)
        allocate(fixed_edges(num_fixed_edges))
        fixed_edges = pack(cresp_all_edges, is_fixed_edge)

        is_active_edge = .false.
        is_active_edge(i_lo:i_up) = .true.
        num_active_edges = count(is_active_edge)
        allocate(active_edges(i_lo:i_up))
        active_edges = pack(cresp_all_edges, is_active_edge)
#ifdef VERBOSE
        print "(2(A9,i3))", "i_lo =", i_lo, ", i_up = ", i_up
#endif /* VERBOSE */
! cleaning (TEST)
        e(:i_lo) = zero
        n(:i_lo) = zero
        e(i_up+1:) = zero
        n(i_up+1:) = zero

      endif
  end subroutine cresp_find_active_bins
!---------------! Compute p for all active edges !---------------------------------------------------  
  subroutine cresp_organize_p
  use initcrspectrum, only: p_fix, ncre
  use constants, only: zero
  implicit none
        p = zero
        p(fixed_edges) = p_fix(fixed_edges)
        p(i_lo) = p_lo * (1 - approx_p_lo ); p_lo = p(i_lo)
        p(i_up) = p_up * (1 - approx_p_up ); p_up = p(i_up) ! cutoff momenta are going to be evaluated with use of get_fqp_lo and get_fqp_up if e_small_approx_* is set
        if ( num_active_bins .eq. 1 .and. (i_lo .gt. 0 .and. i_up .lt. ncre) ) then
            approx_p_lo = 0
            approx_p_up = 0
            p_lo = p_fix(i_lo) ; p(i_lo) = p_fix(i_lo) ! one bin -> momentum fixed grid to boundary momenta
            p_up = p_fix(i_up) ; p(i_up) = p_fix(i_up) ! one bin -> momentum fixed grid to boundary momenta
        endif
  end subroutine cresp_organize_p
!-----------------------------------------------------------------------
  subroutine cresp_detect_negative_content ! Diagnostic measure - negative values should not show up:
   use initcrspectrum, only: ncre
   use constants, only: zero
   implicit none                       ! if they do, there's something wrong with last code modifications
    integer :: i
        do i = 1, ncre
            if (e(i) .lt. zero .or. n(i) .lt. zero .or. edt(i) .lt. zero .or. ndt(i) .lt. zero) then
!                 print *,'Negative value detected:',  'i=', i,'n(i)=', n(i), 'e(i)=',e(i)
!                 print *,'Negative value detected:',  'i=', i,'ndt(i)=', ndt(i), 'edt(i)=',edt(i)
!             call sleep(1)
            endif
        enddo
  end subroutine cresp_detect_negative_content
!----------------------------------------------------------------------------------------------------
  subroutine cresp_update_bin_index(dt, p_lo, p_up, p_lo_next, p_up_next, dt_too_high) ! evaluates only "next" momenta and is called after finding outer cutoff momenta
   use constants, only: zero, I_ZERO, one
   use initcrspectrum, only: ncre, p_fix, w, cresp_all_bins, cresp_all_edges
   implicit none
    real(kind=8), intent(in)  :: dt
    real(kind=8), intent(in)  :: p_lo, p_up
    real(kind=8), intent(out) :: p_lo_next, p_up_next
    logical,      intent(out) :: dt_too_high
    integer                   :: i_lo_next, i_up_next
        dt_too_high = .false.
! Compute p_lo and p_up at [t+dt]
        call p_update(dt, p_lo, p_lo_next)
        call p_update(dt, p_up, p_up_next)
        p_lo_next = abs(p_lo_next)
        p_up_next = abs(p_up_next)
! Locate cut-ofs after current timestep
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
        p_upw(1:ncre) = p_fix(1:ncre)*(one+p_upw_rch(dt,p_fix(1:ncre)))
            
#ifdef VERBOSE
        print*, 'Change of  cut index lo,up:', del_i_lo, del_i_up
#endif /* VERBOSE */      
      
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
#ifdef VERBOSE      
        print *, 'In update_bin_index'
        print *, 'p_lo_next, p_up_next:', p_lo_next, p_up_next
        write (*,"(A15,50L2, 50I3)") 'active_edges: ', is_active_edge, active_edges
        write (*,"(A15,50L2, 50I3)") 'active edgesN:', is_active_edge_next, active_edges_next
        write (*,"(A15,50L2, 50I3)") 'active binsN :', is_active_bin_next , active_bins_next
        write (*,"(A15,50L2, 50I3)") 'fixed  edges: ', is_fixed_edge_next,  fixed_edges_next
        write (*,"(A15,50L2, 50I3)") 'cooling edges:', is_cooling_edge_next,  cooling_edges_next
        write (*,"(A15,50L2, 50I3)") 'heating edges:', is_heating_edge_next,  heating_edges_next
#endif /* VERBOSE */

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
   use initcrspectrum, only: ncre, spec_mod_trms, q_init, p_lo_init, p_up_init, initial_condition, allow_init_extension, & ! f_init, bump_amp
                        e_small_approx_init_cond, e_small_approx_p_lo, e_small_approx_p_up, crel, p_fix, w,&
                        p_min_fix, p_max_fix, add_spectrum_base, e_small, test_spectrum_break, cresp_all_bins
   use cresp_NR_method,only: e_small_to_f
   use constants, only: zero, I_ONE, fpi
   use cresp_variables, only: clight ! use units, only: clight
   implicit none
    integer                          :: i, k, i_lo_ch, i_up_ch, i_br
    real(kind=8)                     :: c
    real(kind=8), dimension(I_ONE:ncre)    :: init_n, init_e
    type (spec_mod_trms), intent(in), optional :: sptab
    real(kind=8), intent(in)         :: f_amplitude
    logical :: exit_code
        u_b = zero ; u_d = zero
        if(present(sptab)) u_b = sptab%ub
        if(present(sptab)) u_d = sptab%ud
        
        approx_p_lo = e_small_approx_p_lo
        approx_p_up = e_small_approx_p_up
        
        init_e = zero
        init_n = zero

        if(present(sptab)) u_b_0 = u_b
        if(present(sptab)) u_d_0 = u_d
       
        f = zero
        q = zero
        p = zero
        n = zero
        e = zero

        q = q_init
! reading initial values of p_lo and p_up 
        p_lo = p_lo_init
        p_up = p_up_init
                  
        p          = p_fix       ! actual array of p including free edges, p_fix shared via initcrspectrum
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

        is_active_bin = .false.
        is_active_bin(i_lo+1:i_up) = .true.
        num_active_bins = count(is_active_bin)
        allocate(active_bins(num_active_bins))
        active_bins = pack(cresp_all_bins, is_active_bin)
!     TESTING ALGORITHM (finding p_lo and p_up using e_small)
#ifdef TEST_CRESP
!      call p_algorithm_accuracy_test  ! tests accuracy of algorithm which later seeks value of p_up using n, e, f and e_small
#endif /* TEST_CRESP */
! Pure power law spectrum initial condition (default case)
        q = q_init
        f = zero
        f = f_amplitude * (p/p_lo_init)**(-q_init)
        if (add_spectrum_base .gt. 0) then
            do i = 0, ncre-1
                if (f(i) .gt. zero ) f(i) = f(i) + e_small_to_f(p(i))
            enddo
        endif
        if (initial_condition == "powl") call cresp_init_powl_spectrum(n, e, f_amplitude, q_init, p_lo_init, p_up_init)
        if (initial_condition == 'brpl') then
! Power law with a break at p_lo_init initial condition
! In this case initial spectrum with a break at p_min_fix is assumed, the initial slope 
! on the left side of the break is just -q_init for simplicity.
            q(i_lo+1) = -q_init
            f(i_lo)   = f(i_lo+1) * (p(i_lo+1)/p_lo_init) ** q(i_lo+1)
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
        if(initial_condition == 'bump') then  ! TODO - @cresp_grid energy normalization and integral to scale cosmic ray electrons with nucleon energy density!
! Gaussian bump-type initial condition for energy distribution
#ifdef VERBOSE
            print *, 'init_state:',sqrt(p_lo_init*p_up_init/1.) !,sp_init_width
            print *, 'init_state:',log(p/sqrt(p_lo_init*p_up_init/1.))
#endif /* VERBOSE */
!             f = f_amplitude * exp(-(2.5*log(p/sqrt(p_lo_init*p_up_init/1.))**2))
            f = f_amplitude * exp(-(4*log(2.0)*log(p/sqrt(p_lo_init*p_up_init/1.))**2)) ! FWHM
            f(0:ncre-1) = f(0:ncre-1) / (fpi * clight * p(0:ncre-1)**(3.0)) ! without this spectrum is gaussian for distribution function
            if (add_spectrum_base .gt. 0) then
                do i = 0, ncre-1
                    if (f(i) .gt. zero )  f(i) = f(i) + e_small_to_f(p(i))
                enddo
            endif
#ifdef VERBOSE
            print *, "f", f
            print *, "clight =", clight
            print *, 'Initial distrib function: i, p, q'
#endif /* VERBOSE */
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
            if ( (approx_p_up + e_small_approx_init_cond ) .gt. 0 )  call get_fqp_up(exit_code)
            if ( (approx_p_lo + e_small_approx_init_cond) .gt. 0 )  call get_fqp_lo(exit_code)

            if (allow_init_extension) then
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
#ifdef VERBOSE
                  print *, 'Extending the range of lower boundary bin after NR_2dim momentum search'
#endif /* VERBOSE */
               enddo

               do i=i_up, i_up_ch-1
                  p(i) = p_fix(i)
                  f(i) = f(i_up-1)* (p_fix(i)/p_fix(i_up-1))**(-q(i_up))
                  q(i) = q(i_up)
#ifdef VERBOSE
                  print *, 'Extending the range of upper boundary bin after NR_2dim momentum search'
#endif /* VERBOSE */
               enddo
#ifdef VERBOSE
               print *, "Boundary bins now (i_lo_new i_lo | i_up_new i_up)",  i_lo_ch, i_lo, ' |', i_up_ch, i_up
#endif /* VERBOSE */
               i_lo = i_lo_ch   ;   i_up = i_up_ch
               q(i_up_ch) = q(i_up)
               p(i_up) = p_fix(i_up);  p(i_up) = p_up

               is_active_bin = .false.
               is_active_bin(i_lo+1:i_up) = .true.
               num_active_bins = count(is_active_bin) ! active arrays must be reevaluated - number of active bins and edges might have changed
               if(allocated(active_bins)) deallocate(active_bins)
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

! testing how the algorithm will handle discontinuity in the spectrum:
        if (test_spectrum_break) then
            e(int((i_lo+i_up)/2):int((i_lo+i_up)/2)+1) = 0.5*e_small ! some arbitrary values
            n(int((i_lo+i_up)/2):int((i_lo+i_up)/2)+1) = 0.1*e_small
        endif

        n_tot0 = sum(n)
        e_tot0 = sum(e)

        init_n = n
        init_e = e

#ifdef VERBOSE   
        print *, ''
        print *, 'n_tot0 =', n_tot0
        print *, 'e_tot0 =', e_tot0
        print *, 'Initialization finished'
#endif /* VERBOSE */
       call deallocate_active_arrays

  end subroutine cresp_init_state

!-------------------------------------------------------------------------------------------------
! Assumes power-law spectrum, without breaks. In principle the same thing is done in cresp_init_state, but
! init_state cannot be called from "outside".
!-------------------------------------------------------------------------------------------------
  subroutine cresp_init_powl_spectrum(n_inout, e_inout, f_in, q_in, p_dist_lo, p_dist_up)
   use constants,      only: zero
   use initcrspectrum, only: ncre, p_fix, w, cresp_all_bins, cresp_all_edges
   use diagnostics,    only: my_deallocate
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
  subroutine e_tot_2_en_powl_init_params(n_inout, e_inout, e_in_total)
   use initcrspectrum, only: ncre, p_lo_init, p_up_init, q_init
   use diagnostics,    only: my_deallocate
      real(kind=8), dimension(1:ncre), intent(inout):: n_inout, e_inout
      real(kind=8), intent(inout)     :: e_in_total
      real(kind=8) :: f_amplitude
         f_amplitude = get_powl_f_ampl(e_in_total, p_lo_init, p_up_init, q_init)
         call cresp_init_powl_spectrum(n_inout, e_inout, f_amplitude, q_init, p_lo_init, p_up_init)
  end subroutine e_tot_2_en_powl_init_params
!-------------------------------------------------------------------------------------------------
  function get_powl_f_ampl(e_tot, p_dist_lo, p_dist_up, q_dist)
  use constants,       only: zero, I_ONE, I_FOUR, fpi
  use cresp_variables, only: clight ! use units,    only: clight
    real(kind=8), intent(in) :: e_tot, p_dist_lo, p_dist_up, q_dist
    real(kind=8)             :: get_powl_f_ampl
        get_powl_f_ampl = zero
        get_powl_f_ampl = (e_tot / (fpi * clight * p_dist_lo ** I_FOUR) ) * ((I_FOUR - q_dist) / &
                          ((p_dist_up/p_dist_lo)**(I_FOUR - q_dist) - I_ONE  ))
  end function get_powl_f_ampl
!-------------------------------------------------------------------------------------------------
  function e_tot_2_f_init_params(e_in_total)
   use initcrspectrum, only: p_lo_init, p_up_init, q_init
   real(kind=8), intent(in) :: e_in_total
   real(kind=8) :: e_tot_2_f_init_params
      e_tot_2_f_init_params = get_powl_f_ampl(e_in_total, p_lo_init, p_up_init, q_init)
  end function e_tot_2_f_init_params
!-------------------------------------------------------------------------------------------------
! Testing p_up / p_lo finding algorithm, p_up -> p_u to avoid collisions with p_up values.
!-------------------------------------------------------------------------------------------------
! #ifdef TEST_CRESP
   subroutine p_algorithm_accuracy_test
   use initcrspectrum, only: ncre
   use cresp_NR_method, only : NR_get_solution_lo, NR_get_solution_up, e_small_to_f
    implicit none
     real(kind=8),dimension(1) :: q_up, q_lo   ! value of exponent to power law to be obtained using valus of f & p
     real(kind=8),dimension(1) :: f_l, f_r, p_l, p_r, q_new !, f_floor ! given "left" and "right" values of f & p
     real(kind=8), dimension(ncre) :: n, e, n_new, e_new   ! obtained from given, size necessary
     integer(kind=4), dimension(1) :: one_bin  ! 1-dim array, necessary to use existing functions & subroutines
     real(kind=8), dimension(2) :: sought_vars ! p_r and f_l vector
     real(kind=8), dimension(0:ncre)::  f_l_new
     logical :: exit_code
     
      one_bin = 1 ! simplest case
      e = 0.0; n = 0.0; e_new = 0.0; n_new = 0.0;  

      print *, '2D N-R p_u search test:'
      p_l = 1.0e2 
      p_r = 1.0e3  !
      f_l = 1.0 
      f_r = 0.1   ! Thus obtained q should be equal -1.0

      q_up(1) = pf_to_q(p_l(1), p_r(1), f_l(1), f_r(1))
      
      print *, 'Initial test parameters (p_u):'
      print '(A10, 1E18.9, A10, 1E18.9)','p_l:', p_l, ', p_r:', p_r
      print '(A10, 1E18.9, A10, 1E18.9)','f_l:', f_l, ', f_r:', f_r
      print '(A10, 1E18.9)', " q_up:", q_up
      ! computing energy & number density
      n = fq_to_n(p_l, p_r, f_l, q_up, one_bin)
      e = fq_to_e(p_l, p_r, f_l, q_up, one_bin)

      f_r = real(1.0e-2,kind=8) ! we change f_floor
      sought_vars(1) = p_l(1)*5
      sought_vars(2) = f_r(1)*50
      exit_code = .true.
      call NR_get_solution_up(sought_vars,p_l(1), e(1), n(1), exit_code)

      p_r(1) = sought_vars(2);   f_l_new(1) = sought_vars(2);   f_l = f_l_new(1) ;  f_r = e_small_to_f(p_r(1))
      print '(A10, 1E18.9, A10, 1E18.9)','NEW: p_r:', p_r(1), ', f_l_new:', f_l(1)
      q_new(1) = pf_to_q(p_l(1), p_r(1), f_l(1), f_r(1)); q_up(1) = q_new(1)
      
      n_new = fq_to_n(p_l, p_r, f_l, q_new, one_bin)
      e_new = fq_to_e(p_l, p_r, f_l, q_new, one_bin)

      print *, 'e_new=', e_new(1), '  n_new=', n_new(1) , ' e_new/n_new=',e_new(1)/n_new(1)
      print *, 'e    =', e(1),     '  n    =', n(1) ,     ' e/n        =',e(1)/n(1)
      print *, "Rel. error %(e):", (e(1) - e_new(1))/e(1)*100.0, ", (n): % ", (n(1) - n_new(1))/n(1) *100.0
      
      print *, '2D N-R p_lo search test:'
      p_l = 1.0
      p_r = 10.0  ! p_l was not altered
      f_l = 0.1 
      f_r = 1.0   ! Thus obtained q should be equal -1.0

      q_lo(1) = pf_to_q(p_l(1), p_r(1), f_l(1), f_r(1))

      print *, "Initial test parameters (p_lo):"
      print '(A10, 1E18.9, A10, 1E18.9)','p_l:', p_l, ', p_r:', p_r
      print '(A10, 1E18.9, A10, 1E18.9)','f_l:', f_l, ', f_r:', f_r
      print '(A10, 1E18.9)', " q_up:", q_up

      ! computing energy & number density
      n = fq_to_n(p_l, p_r, f_l, q_lo, one_bin)
      e = fq_to_e(p_l, p_r, f_l, q_lo, one_bin)
      
      f_l = real(1.0e-5,kind=8)
      
      sought_vars(1) = real(1.0e-0,kind=8)
      sought_vars(2) = real(f_l(1)*2d4,kind=8)
      print *, "init_nr_2dim params:", sought_vars

       call NR_get_solution_lo(sought_vars, p_r(1), e(1), n(1), exit_code)

      p_l = sought_vars(1);  f_r = sought_vars(2) ;  f_l = e_small_to_f(p_l(1))
      print *, p_l(1), p_r(1), f_l(1), f_r(1)

      q_new = pf_to_q(p_l(1), p_r(1), f_l(1), f_r(1))

      n_new = fq_to_n(p_l, p_r, f_l, q_new, one_bin)
      e_new = fq_to_e(p_l, p_r, f_l, q_new, one_bin)
      print *, 'e_new=', e_new(1), '  n_new=', n_new(1) , ' e_new/n_new=',e_new(1)/n_new(1)
      print *, 'e    =', e(1),     '  n    =', n(1) ,     ' e/n        =',e(1)/n(1)
      print *, "Rel. error %(e):", (e(1) - e_new(1))/e(1)*100.0, ", (n): % ", (n(1) - n_new(1))/n(1) *100.0
      
   end subroutine p_algorithm_accuracy_test

! #endif /* TEST_CRESP  */
!-------------------------------------------------------------------------------------------------
! 
! energy integral (eq. 21)
!
!-------------------------------------------------------------------------------------------------
 
   function fq_to_e(p_l, p_r, f_l, q, bins)
    use initcrspectrum, only: ncre, eps
    use constants, only: zero, one, fpi
    use cresp_variables, only: clight ! use units, only: clight
      implicit none
      real(kind=8), dimension(:), intent (in)  :: p_l, p_r, f_l, q
      integer, dimension(:), intent(in)        :: bins
      real(kind=8), dimension(size(bins))  :: e_bins
      real(kind=8), dimension(1:ncre)      :: fq_to_e

      fq_to_e = zero
      e_bins = fpi*clight*f_l(bins)*p_l(bins)**4
      where(abs(q(bins) - four) .gt. eps) 
         e_bins = e_bins*((p_r(bins)/p_l(bins))**(four-q(bins)) - one)/(four - q(bins))
      elsewhere
         e_bins = e_bins*log(p_r(bins)/p_l(bins))
      end where
      
      fq_to_e(bins) = e_bins
      
   end function fq_to_e
 
!-------------------------------------------------------------------------------------------------
! 
! density integral (eq. 9)
!
!-------------------------------------------------------------------------------------------------

   function fq_to_n(p_l, p_r, f_l, q, bins)
    use initcrspectrum, only: ncre, eps
    use constants, only: zero, one, fpi
      implicit none
      real(kind=8), dimension(:), intent(in) :: p_l, p_r, f_l, q
      integer, dimension(:), intent(in)     :: bins
      real(kind=8), dimension(size(bins))   :: n_bins
      real(kind=8), dimension(1:ncre)       :: fq_to_n
      
      n_bins = zero
      
      n_bins = fpi*f_l(bins)*p_l(bins)**3
      where(abs(q(bins) - three) .gt. eps) 
         n_bins = n_bins*((p_r(bins)/p_l(bins))**(three-q(bins)) - one)/(three - q(bins))
      elsewhere
         n_bins = n_bins*log((p_r(bins)/p_l(bins)))
      end where
      
      fq_to_n = zero
      fq_to_n(bins) = n_bins
      
   end function fq_to_n

!-------------------------------------------------------------------------------------------------
   
   subroutine deallocate_active_arrays
   
      implicit none
      if(allocated(fixed_edges)) deallocate(fixed_edges)
      if(allocated(fixed_edges_next)) deallocate(fixed_edges_next)
      if(allocated(active_edges)) deallocate(active_edges)
      if(allocated(active_edges_next)) deallocate(active_edges_next)
      if(allocated(active_bins)) deallocate(active_bins)
      if(allocated(active_bins_next)) deallocate(active_bins_next)
      if(allocated(cooling_edges_next)) deallocate(cooling_edges_next)
      if(allocated(heating_edges_next)) deallocate(heating_edges_next)

   end subroutine deallocate_active_arrays

!-------------------------------------------------------------------------------------------------
! 
! compute fluxes
!
!-------------------------------------------------------------------------------------------------
   subroutine cresp_compute_fluxes(ce,he)
    use initcrspectrum, only: ncre, eps, cresp_all_bins
    use constants, only: zero, one, fpi
    use cresp_variables, only: clight ! use units, only: clight
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
      where(abs( qi(ce) - three ) .gt. eps)
         dn_upw(ce) = dn_upw(ce)*((p_upw(ce)/pimh(ce))**(three-qi(ce)) - one)/(three - qi(ce))
      elsewhere
         dn_upw(ce) = dn_upw(ce)*log((p_upw(ce)/pimh(ce)))
      end where
      nflux(ce) = - dn_upw(ce)
            
      de_upw(ce) = fpi*clight*fimh(ce)*pimh(ce)**4
      where(abs(qi(ce) - four) .gt. eps) 
         de_upw(ce) = de_upw(ce)*((p_upw(ce)/pimh(ce))**(four-qi(ce)) - one)/(four - qi(ce))
      elsewhere
         de_upw(ce) = de_upw(ce)*log(p_upw(ce)/pimh(ce))
      end where
      eflux(ce) =  - de_upw(ce)
      
      if(del_i_up == -1) then 
         nflux(i_up-1) = -n(i_up)
         eflux(i_up-1) = -e(i_up  )
      endif
      
! filling empty empty bin - switch of upper boundary, condition is checked only once per flux computation and is very rarely satisfied.
      if (nflux(i_up) .gt. zero) then             ! If flux is greater than zero it will go through right edge, activating next bin in the next timestep.
        if ( cresp_all_bins(i_up+1) .eq. i_up+1 ) then  ! But it shuld only happen if there is bin with index i_up+1
         ndt(i_up+1) = nflux(i_up)
         edt(i_up+1) = eflux(i_up)
#ifdef VERBOSE
         print *, ' **** UPPER BOUND +1 ****'
#endif /* VERBOSE */
         del_i_up = +1
        endif
      endif
      
      if ( nflux(i_up-1)+n(i_up) .le. zero) then ! If flux is equal or greater than energy / density in a given bin,  these both shall migrate
         nflux(i_up-1) =  -n(i_up)                   ! to an adjacent bin, thus making given bin detected as inactive (empty) in the next timestep
         eflux(i_up-1) =  -e(i_up)
#ifdef VERBOSE
         print *, " **** UPPER BOUND -1 **** "
#endif /* VERBOSE */
         del_i_up = -1
      endif
      
      dn_upw(he) = fpi*fimth(he)*p_upw(he)**3*(pimth(he)/p_upw(he))**qim1(he)
      where(abs(qim1(he) - three) .gt. eps ) 
         dn_upw(he) = dn_upw(he)*((pimh(he)/p_upw(he))**(three-qim1(he)) - one)/(three - qim1(he))
      elsewhere
         dn_upw(he) = dn_upw(he)*log((pimh(he)/p_upw(he)))
      end where
      nflux(he) = dn_upw(he)
            
      de_upw(he) = fpi*clight*fimth(he)*p_upw(he)**4*(pimth(he)/p_upw(he))**qim1(he)
      where(abs(qi(he) - four) .gt. eps) 
         de_upw(he) = de_upw(he)*((pimh(he)/p_upw(he))**(four-qim1(he)) - one)/(four - qim1(he))
      elsewhere
         de_upw(he) = de_upw(he)*log(pimh(he)/p_upw(he))
      end where
      eflux(he) = de_upw(he)
      
      if(del_i_lo == 1 .or. nflux(i_lo+1) .ge. n(i_lo+1) ) then
         nflux(i_lo+1) = n(i_lo+1)
         eflux(i_lo+1) = e(i_lo+1)
! emptying lower boundary bins - in cases when flux gets greater than energy or number density 
#ifdef VERBOSE
         print *, ' **** LOWER BOUND +1 ****'
#endif /* VERBOSE */
         del_i_lo = 1   ! in case it hasn't yet been modified
      endif

   end subroutine cresp_compute_fluxes
      
!-------------------------------------------------------------------------------------------------
! 
! compute R (eq. 25)
!
!-------------------------------------------------------------------------------------------------

   subroutine cresp_compute_r(p, bins)
    use initcrspectrum, only: ncre, eps
    use constants, only: zero
    implicit none
      real(kind=8), dimension(0:ncre), intent(in) :: p
      integer, dimension(:), intent(in)     :: bins
      real(kind=8), dimension(size(bins)) :: r_num, r_den
      
      r = zero
      where(abs(q(bins) - five) .gt. eps) 
         r_num = (p(bins)**(five-q(bins)) - p(bins-1)**(five-q(bins)))/(five - q(bins))
      elsewhere
         r_num = log(p(bins)/p(bins-1))
      end where
      
      where(abs(q(bins) - four) .gt. eps) 
         r_den = (p(bins)**(four-q(bins)) - p(bins-1)**(four-q(bins)))/(four - q(bins))
      elsewhere
         r_den = log(p(bins)/p(bins-1))
      end where
      where (r_num .ne. zero .and. r_den .ne. zero)
      r(bins) = u_d + u_b * r_num/r_den   !!! all cooling effects will come here
      end where

   end subroutine cresp_compute_r

!-------------------------------------------------------------------------------------------------
! 
! find new q (eq. 29) and new f
!
!-------------------------------------------------------------------------------------------------

  subroutine ne_to_q(n, e, q, bins)
   use initcrspectrum, only: ncre, e_small
   use constants, only: zero
   use cresp_variables, only: clight ! use units, only: clight
   use cresp_NR_method, only: compute_q
   implicit none
    real(kind=8), dimension(1:ncre), intent(in)  :: n, e !
    real(kind=8), dimension(1:ncre), intent(out) :: q
    integer(kind=4), dimension(:), intent(in)    :: bins
    integer          :: i, i_active
    real(kind=8)     :: alpha_in
    logical :: exit_code

    q = zero
    do i_active = 1 + approx_p_lo, size(active_bins) - approx_p_up
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
   use constants, only: zero, one, fpi
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
            where(abs(q(bins)-three) .gt. eps) 
                f_bins = f_bins*(three - q(bins)) /((pr_by_pl(bins))**(three-q(bins)) - one)
            elsewhere
                f_bins = f_bins/log((p_r(bins)/p_l(bins)))
            end where
        end where
        nq_to_f(bins-1) = f_bins
   end function nq_to_f
!---------------------------------------------------------------------------------------------------
! Computing cosmic ray pressure (eq. 44)
!---------------------------------------------------------------------------------------------------
   function get_pcresp(p_l, p_r, f_l, q, bins)
    use initcrspectrum, only: ncre, eps
    use constants, only: zero, one, three, fpi
    use cresp_variables, only: clight ! use units, only: clight
      implicit none
      real(kind=8), dimension(:), intent (in)  :: p_l, p_r, f_l, q
      integer, dimension(:), intent(in)        :: bins
      real(kind=8), dimension(size(bins))  :: p_cresp
      real(kind=8)              :: get_pcresp

      get_pcresp = zero
      p_cresp = (fpi/three) * clight*f_l(bins)*p_l(bins)**4
      where(abs(q(bins) - four) .gt. eps)
         p_cresp = p_cresp*((p_r(bins)/p_l(bins))**(four-q(bins)) - one)/(four - q(bins))
      elsewhere
         p_cresp = p_cresp*log(p_r(bins)/p_l(bins))
      end where
    get_pcresp = sum(p_cresp)

   end function get_pcresp
!---------------------------------------------------------------------------------------------------
! Computing cosmic ray spectrum component pressure gradient
!---------------------------------------------------------------------------------------------------
  subroutine src_gpcresp(u, n, dx, grad_pcresp)
   use initcrspectrum, only: ncre, cre_active, cre_gpcr_ess
   use constants, only: onet
    implicit none
    integer(kind=4), intent(in) :: n
    real(kind=8), intent(in)                       :: dx
    real(kind=8), dimension(n, 1:ncre), intent(in) :: u
    real(kind=8), dimension(n), intent(out)        :: grad_pcresp
    real(kind=8), dimension(n)                     :: P_cresp_r, P_cresp_l
!     if (ultrarelativistic) then
    grad_pcresp = 0.0 ;  P_cresp_l = 0.0 ;  P_cresp_r = 0.0
    P_cresp_l(1:n-2) = onet * sum(u(1:n-2, :),dim=2)
    P_cresp_r(3:n)   = onet * sum(u(3:n,   :),dim=2)
    if (cre_gpcr_ess) grad_pcresp(2:n-1) = cre_active * (P_cresp_l(1:n-2) - P_cresp_r(3:n) )/(2.*dx)
!     endif
 end subroutine src_gpcresp
!---------------------------------------------------------------------------------------------------
! Preparation and computatuon of upper boundary momentum "p_up" and and upper boundary 
! distribution function value on left bin edge "f"
!---------------------------------------------------------------------------------------------------
  subroutine get_fqp_up(exit_code)
   use initcrspectrum, only: e_small, q_big, p_fix, NR_refine_solution_pf
   use constants, only: zero, one, two
   use cresp_NR_method,only: intpol_pf_from_NR_grids, alpha, n_in, selected_function_2D, fvec_up, &
            NR_algorithm, e_small_to_f, q_ratios
   use cresp_variables, only: clight ! use units, only: clight
   implicit none
    real(kind=8), dimension(1:2) :: x_NR, x_NR_init
    logical :: exit_code, interpolated, intpol_fail
    character(len=2) :: bound = "up"
        x_NR = zero
        alpha = (e(i_up)/(n(i_up)*clight*p_fix(i_up-1)))
        n_in  = n(i_up)
        x_NR = intpol_pf_from_NR_grids(bound, alpha, n_in, interpolated, intpol_fail)
        if (intpol_fail) then
            exit_code = .true.
            fail_count_no_sol(2) = fail_count_no_sol(2) + 1
            return
        endif
        x_NR_init = x_NR
        selected_function_2D => fvec_up
#ifdef VERBOSE
        write (*,"(A31,2E22.15)" ) "Input ratios(p, f) for NR (up):", x_NR
#endif /* VERBOSE */
        if ( (NR_refine_solution_pf .eqv. .true.) .or. (interpolated .eqv. .false.)) then
            if (interpolated .eqv. .false.) fail_count_interpol(2) = fail_count_interpol(2) +1
            call NR_algorithm(x_NR, exit_code)
            if (exit_code .eqv. .true.) then ! some failures still take place
                if (interpolated .eqv. .false.) then
                    exit_code = .true.
#ifdef VERBOSE
                    print *, " Interpolation AND NR failure (up)", alpha, n_in, x_NR_init
#endif /* VERBOSE */
                    return
                endif
                fail_count_NR_2dim(2) = fail_count_NR_2dim(2) +1
                x_NR = x_NR_init
                print *, "Interpolated?", interpolated, "NR_refine_solution_pf?", NR_refine_solution_pf,"solved?", exit_code
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
        if (abs(q(i_up)) .gt. two * q_big ) q(i_up) = sign(one, q(i_up)) * q_big
#ifdef VERBOSE
        write (*,"(A1)") " "
        write (*,"(A26,2E22.15)") " >>> Obtained (p_up, f_l):", p_up, f(i_up-1) &
                                 ,"     Corresponding ratios:", x_NR(1), x_NR(2)
#endif /* VERBOSE */
  end subroutine get_fqp_up
!--------------------------------------------------------------------------------------------------
! Preparation and computation of upper boundary momentum "p_lo" and and upper boundary 
! distribution function value on the right bin edge "f"
!--------------------------------------------------------------------------------------------------
  subroutine get_fqp_lo(exit_code)
   use initcrspectrum, only: e_small, q_big, p_fix, NR_refine_solution_pf
   use constants, only: zero, one
   use cresp_NR_method, only: intpol_pf_from_NR_grids, alpha, n_in, selected_function_2D, fvec_lo, &
                NR_algorithm, e_small_to_f, q_ratios
   use cresp_variables, only: clight ! use units, only: clight
   implicit none
    real(kind=8), dimension(1:2) :: x_NR, x_NR_init
    logical :: exit_code, interpolated, intpol_fail
    character(len=2) :: bound = "lo"
        x_NR = zero
        alpha = (e(i_lo+1)/(n(i_lo+1)*clight*p_fix(i_lo+1)))
        n_in  = n(i_lo+1)
        x_NR = intpol_pf_from_NR_grids(bound, alpha, n_in, interpolated, intpol_fail)
        if (intpol_fail) then
            exit_code = .true.
            fail_count_no_sol(1) = fail_count_no_sol(1) + 1
            return
        endif
        x_NR_init = x_NR
        selected_function_2D => fvec_lo 
#ifdef VERBOSE
        write (*,"(A31,2E22.15)" ) "Input ratios(p, f) for NR (lo):", x_NR
#endif /* VERBOSE */
        if ( (NR_refine_solution_pf .eqv. .true.) .or. (interpolated .eqv. .false.)) then
            if (interpolated .eqv. .false.) fail_count_interpol(1) = fail_count_interpol(1) +1
            call NR_algorithm(x_NR, exit_code)
            if (exit_code .eqv. .true.) then ! some failures still take place
                if (interpolated .eqv. .false.) then
                    exit_code = .true.
#ifdef VERBOSE
                    print *, " Interpolation AND NR failure (lo)", alpha, n_in
#endif /* VERBOSE */
                    return
                endif
                fail_count_NR_2dim(1) = fail_count_NR_2dim(1) +1
                x_NR = x_NR_init
                print *, "Interpolated?", interpolated, "NR_refine_solution_pf?", NR_refine_solution_pf,"solved?", exit_code
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
        if (abs(q(i_lo+1)) .gt. 2 * q_big ) q(i_lo+1) = sign(one, q(i_lo+1)) * q_big
#ifdef VERBOSE
        write (*,"(A1)") " "
        write (*,"(A26,2E22.15)") " >>> Obtained (p_lo, f_r):", p_lo, x_NR(2)*f(i_lo) &
                                 ,"     Corresponding ratios:", x_NR(1), x_NR(2)
#endif /* VERBOSE */
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
! 
! relative change of momentum due to losses (u_b*p*dt) and compression u_d*dt (Taylor expansion up to 3rd order)
!
!-------------------------------------------------------------------------------------------------
  function p_rch(dt, p)
   use initcrspectrum, only: taylor_coeff_2nd, taylor_coeff_3rd 
   use constants, only: half, sixth
   implicit none
    real(kind=8), intent(in)  :: dt
    real(kind=8), intent(in)  :: p
    real(kind=8)              :: p_rch
   !p_rch =  (u_b*p + u_d) * dt !!! b_sync_ic = 8.94e-25*(u_b+u_cmb)*gamma_l**2 ! erg/cm
        p_rch = (- u_d - p * u_b ) *  dt  + taylor_coeff_2nd *( half*u_d**2 + u_b**2 * p**2)*dt**2 &
                                     + taylor_coeff_3rd *(-sixth*u_d**3 - u_b**3 * p**3)*dt**3 ! analitycally correct solution
  end function p_rch
!-------------------------------------------------------------------------------------------------
  function p_upw_rch(dt, p)
   use initcrspectrum, only: taylor_coeff_2nd, taylor_coeff_3rd 
   use constants, only: half, sixth
   implicit none
    real(kind=8), intent(in)  :: dt
    real(kind=8), dimension(:), intent(in) :: p
    real(kind=8), dimension(size(p))  :: p_upw_rch
   !p_upw_rch =  (u_b*p + u_d) * dt      !!! b_sync_ic = 8.94e-25*(u_b+u_cmb)*gamma_l**2 ! erg/cm
        p_upw_rch = (u_d + p * u_b) * dt + taylor_coeff_2nd * (half*u_d**2 + u_b**2 * p**2)*dt**2 &
                                     + taylor_coeff_3rd *(sixth*u_d**3 + u_b**3 * p**3)*dt**3 ! analitycally correct
  end function p_upw_rch
  
!----------------------------------------------------------------------------------------------------
  subroutine boundary_flux_check
   use constants, only: I_ZERO, zero
   use initcrspectrum, only: p_fix
   implicit none
  
     if ( abs(p_up_next - p_fix(i_up))/p_fix(i_up) .lt. 3.0e-3 .and. del_i_up .ne. I_ZERO ) then
! If the difference of neighbouring momenta is too small, zero flux is assumed to avoid numerical errrors, and these quantities are moved to "virtual" array, postponing new bin activation, until
! the condition is satisfied; q can be therefore assumed zero as the bin not be subject to any energy losses.
       vrtl_n(2) = vrtl_n(2) + ndt(i_up+1) ; ndt(i_up+1) = zero
       vrtl_e(2) = vrtl_e(2) + edt(i_up+1) ; edt(i_up+1) = zero
!#ifdef VERBOSE
!        print *, "(p_up) relative momentum error exceeded, activating virtual n and e"
!        print *, "virtual n up", vrtl_n(2), "virtual n up", vrtl_e(2)
!#endif /* VERBOSE */
       q(i_up+1) = zero
    else
       ndt(i_up) = ndt(i_up) + vrtl_n(2) ; vrtl_n(2) = zero
       edt(i_up) = edt(i_up) + vrtl_e(2) ; vrtl_e(2) = zero
    endif

   if ( abs(p_lo_next - p_fix(i_lo))/p_fix(i_lo) .lt. 3.0e-3 .and. del_i_lo .ne. I_ZERO ) then
! If the difference of neighbouring momenta is too small, zero flux is assumed to avoid numerical errors, and these quantities are moved to "virtual" array, postponing new bins activation, until
! the condition is satisfied; q can be therefore assumed zero as the bin not be subject to any energy losses.
       vrtl_n(1) = vrtl_n(1) + ndt(i_lo) ; ndt(i_lo) = zero
       vrtl_e(1) = vrtl_e(1) + edt(i_lo) ; edt(i_lo) = zero
!#ifdef VERBOSE
!        print *, "(p_lo,p_fix) relative error exceeded, activating virtual n and e"
!        print *, "virtual n lo", vrtl_n(1), "virtual e lo", vrtl_e(1)
!#endif /* VERBOSE */
       q(i_lo+1) = zero
    else
       ndt(i_lo+1) = ndt(i_lo+1) + vrtl_n(1) ; vrtl_n(1) = zero
       edt(i_lo+1) = edt(i_lo+1) + vrtl_e(1) ; vrtl_e(1) = zero
    endif
  end subroutine boundary_flux_check
!----------------------------------------------------------------------------------------------------
  subroutine transfer_quantities(give_to, take_from)
  use constants, only: zero
  implicit none
    real(kind=8), intent(inout) :: give_to, take_from
        give_to = give_to + take_from
        take_from = zero
  end subroutine transfer_quantities
!----------------------------------------------------------------------------------------------------
  subroutine threshold_energy_check_lo(e_tab, n_tab, e_lo_lt_e_small, solution_failed)
  use initcrspectrum, only: e_small
  use constants, only: zero
  implicit none
    real(kind=8), dimension(:), intent(inout) :: e_tab, n_tab
    logical :: e_lo_lt_e_small, solution_failed
      if ( solution_failed .or. ((e_tab(i_lo+1) .le. e_small) .and. (approx_p_lo .eq. 1))) then
        vrtl_e(1) = vrtl_e(1) + e_tab(i_lo+1)  ; e_tab(i_lo+1) = zero
        vrtl_n(1) = vrtl_n(1) + n_tab(i_lo+1)  ; n_tab(i_lo+1) = zero
!         print *, "e_tab(i_lo+1) close to e_small, activating virtual n and e"
!         print *, "virtual n lo", vrtl_n(1), "virtual e lo", vrtl_e(1)
        e_lo_lt_e_small = .true.
    else
!         print *, "virtual n lo", vrtl_n(1), "virtual e lo", vrtl_e(1), "e, n: ", e_tab(i_lo+1), n_tab(i_lo+1)
        e_tab(i_lo+1) = e_tab(i_lo+1) + vrtl_e(1) ; vrtl_e(1) = zero
        n_tab(i_lo+1) = n_tab(i_lo+1) + vrtl_n(1) ; vrtl_n(1) = zero
        e_lo_lt_e_small = .false.
    endif
  end subroutine threshold_energy_check_lo
!----------------------------------------------------------------------------------------------------
  subroutine threshold_energy_check_up(e_tab, n_tab, e_up_lt_e_small, solution_failed)
  use initcrspectrum, only: e_small
  use constants, only: zero
  implicit none
    real(kind=8), dimension(:), intent(inout) :: e_tab, n_tab
    logical :: e_up_lt_e_small, solution_failed
    if ( solution_failed .or. ((e_tab(i_up) .le. e_small) .and. (approx_p_up .eq. 1))) then
        vrtl_e(2) = vrtl_e(2) + e_tab(i_up)  ; e_tab(i_up) = zero
        vrtl_n(2) = vrtl_n(2) + n_tab(i_up)  ; n_tab(i_up) = zero
!         print *, "e_tab(i_up) close to e_small, activating virtual n and e"
!         print *, "virtual n up", vrtl_n(2), "virtual e up", vrtl_e(2)
        e_up_lt_e_small = .true.
    else
!         print *, "virtual n up", vrtl_n(2), "virtual e up", vrtl_e(2), "e, n: ", e_tab(i_up), n_tab(i_up)
        e_tab(i_up) = e_tab(i_up) + vrtl_e(2) ; vrtl_e(2) = zero
        n_tab(i_up) = n_tab(i_up) + vrtl_n(2) ; vrtl_n(2) = zero
        e_up_lt_e_small = .false.
    endif
  end subroutine threshold_energy_check_up
!----------------------------------------------------------------------------------------------------
  subroutine cresp_allocate_all
   use initcrspectrum, only: ncre
   use diagnostics, only: my_allocate_with_index
   implicit none
     integer(kind = 4)          :: ma1d
   
   ma1d = ncre
   call my_allocate_with_index(fail_count_comp_q,ma1d,1)
   
   call my_allocate_with_index(n,ma1d,1)   !:: n, e, r
   call my_allocate_with_index(e,ma1d,1)
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
   call my_allocate_with_index(not_spectrum_break,ma1d,1)
   
  end subroutine cresp_allocate_all
  
  subroutine cresp_deallocate_all! called by driver
  use diagnostics, only: my_deallocate
  implicit none
  
   call my_deallocate(n)   !:: n, e, r
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
   
   call deallocate_active_arrays ! optional
   
  end subroutine cresp_deallocate_all

! !---------------------------------------------------------------------------------------------
   subroutine cresp_accuracy_test(t)
     implicit none
     real(kind=8), intent(in)   :: t 
   
   print *, " -------------------------- "

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
        print '(A36,I6,A6,I6)', "NR_2dim:  convergence failure: p_lo", fail_count_NR_2dim(1), ", p_up", fail_count_NR_2dim(2)
        print '(A36,I6,A6,I6)', "NR_2dim:interpolation failure: p_lo", fail_count_interpol(1), ", p_up", fail_count_interpol(2)
        print '(A36,I6,A6,I6)', "NR_2dim:  no solution failure: p_lo", fail_count_no_sol(1), ", p_up", fail_count_no_sol(2)
        print '(A36,I6,A6,I6)', "NR_2dim: second try failure  : p_lo", second_fail(1), ", p_up", second_fail(2)
        print '(A36,   100I5)', "NR_2dim:inpl/solve  q(bin) failure:", fail_count_comp_q
        call cresp_deallocate_all
 end subroutine cleanup_cresp
!----------------------------------------------------------------------------------------------------
   subroutine printer(t)
   use initcrspectrum, only: ncre, crel
   implicit none
     real(kind = 8)   :: t
      open(10, file="crs.dat", position='append')
      write(10, '(2e16.9, 3(1x,i8), 200(1x,ES18.9E3))') t, crel%dt, ncre, crel%i_lo, crel%i_up, crel%p, crel%f, crel%q
      close(10)
      
      open(11, file="crs_ne.dat", position='append')
      write(11, '(2I5,4x, e16.9, 100(1x,F18.9))') del_i_lo, del_i_up, t, crel%dt, crel%p(i_lo), crel%p(i_up), crel%n, crel%e
      close(11)

      end subroutine printer

! -------------------- 
  
end module cresp_crspectrum
