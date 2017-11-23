module initcrspectrum  
! pulled by COSM_RAY_ELECTRONS

 public ! QA_WARN no secrets are kept here
 private init_cresp_types

! contains routines reading namelist in problem.par file dedicated to cosmic ray electron spectrum and initializes types used.
! available via namelist COSMIC_RAY_SPECTRUM
 integer(kind=4)    :: ncre      = 0            ! < number of bins
 real(kind=8)       :: p_min_fix = 1.5e1        ! < fixed momentum grid lower cutoff
 real(kind=8)       :: p_max_fix = 1.65e4       ! < fixed momentum grid upper cutoff
 real(kind=8)       :: p_lo_init = 1.5e1        ! < initial lower cutoff momentum
 real(kind=8)       :: p_up_init = 7.5e2        ! < initial upper cutoff momentum
 character(len=4)   :: initial_condition = "bump"   ! < available types: bump, powl, brpl, symf, syme. Description below.
 real(kind=8)       :: f_init    = 1.0          ! < initial value of distr. func. for isolated case 
 real(kind=8)       :: q_init    = 2.5          ! < initial value of power law-like spectrum exponent
 real(kind=8)       :: bump_amp  = 0.5d0        ! < bump amplitude for gaussian-like energy spectrum
 real(kind=8)       :: q_big     = 30.0d0       ! < maximal amplitude of q
 real(kind=8)       :: cfl_cre   = 0.4          ! < CFL parameter  for cr electrons
 real(kind=8)       :: cre_eff   = 0.01         ! < fraction of energy passed to cr-electrons by nucleons (mainly protons)
 real(kind=8)       :: K_cre_paral_1 = 0        ! < parallell diffusion coefficient
 real(kind=8)       :: K_cre_perp_1  = 0        ! < perpendicular diffusion coefficient
 real(kind=8)       :: K_cre_pow     = 0        ! < exponent for power law-like diffusion-energy dependance
 integer(kind=4)    :: expan_order   = 1        ! < 1,2,3 order of Taylor expansion for p_update (cresp_crspectrum)
 real(kind=8)       :: e_small       = 1.0e-5   ! lower energy cutoff for energy-approximated cutoff momenta
 integer(kind=1)    :: e_small_approx_p_lo = 1  ! < 0,1 - turns off/on energy (e_small) approximated lower cutoff momentum in isolated case
 integer(kind=1)    :: e_small_approx_p_up = 1  ! < 0,1 - turns off/on energy (e_small) approximated upper cutoff momentum in isolated case
 integer(kind=1)    :: e_small_approx_init_cond = 1 ! < 0,1 - turns off/on energy (e_small) approximated momenta at initialization
 integer(kind=1)    :: add_spectrum_base = 1    ! < adds base to spectrum of any type of e_small value
 real(kind=8)       :: max_p_ratio = 2.5        ! < maximal ratio of momenta for solution grids resolved at initialization via cresp_NR_method
 integer(kind=2)    :: NR_iter_limit=100        ! < maximal number of iterations for NR algorithm
 logical            :: force_init_NR = .false.  ! < forces resolving new ratio solution grids at initialization
 logical            :: NR_run_refine_pf  = .false.      ! < enables "refine_grids" subroutines that fill empty spaces on the solution grid
 logical            :: NR_refine_solution_q = .false.    ! < enables NR_1D refinement for value of interpolated "q" value
 logical            :: NR_refine_solution_pf = .false.  ! < enables NR_2D refinement for interpolated values of "p" and "f". Note - algorithm tries to refine values if interpolation was unsuccessful.
 logical            :: prevent_neg_e = .true.   ! < forces e,n=eps where e or n drops below zero due to diffusion algorithm (TEMP workaround)
 logical            :: test_spectrum_break   = .false.  ! < introduce break in the middle of the spectrum (to see how algorithm handles it), TEMP
 
 real(kind=8)       :: tol_f = 1.0e-11          ! < tolerance for f abs. error in NR algorithm
 real(kind=8)       :: tol_x = 1.0e-11          ! < tolerance for x abs. error in NR algorithm
 real(kind=8)       :: tol_f_1D = 1.0e-14 ! < tolerance for f abs. error in NR algorithm (1D)
 real(kind=8)       :: tol_x_1D = 1.0e-14 ! < tolerance for x abs. error in NR algorithm (1D)
 integer(kind=4)    :: arr_dim = 200
!----------------------------------
 real(kind=8),parameter  :: eps   = 1.0e-15          ! < epsilon parameter for real number comparisons 
!----------------------------------
 real(kind=8),allocatable, dimension(:) :: p_fix, p_mid_fix
 real(kind=8)       :: w
! Types used in module: 
   type bin_old
    integer                           :: i_lo
    integer                           :: i_up
    real(kind=8),allocatable,dimension(:)   :: p
    real(kind=8),allocatable,dimension(:)   :: f
    real(kind=8),allocatable,dimension(:)   :: q
    real(kind=8),allocatable,dimension(:)   :: e
    real(kind=8),allocatable,dimension(:)   :: n
    real(kind=8) :: dt
  end type bin_old

  type(bin_old) crel
! For passing terms to compute energy sources / sinks
  type spec_mod_trms
   real(kind=8)    :: ub
   real(kind=8)    :: ud
   real(kind=8)    :: ucmb
  end type spec_mod_trms
  
  real(kind=8), allocatable, dimension(:,:,:,:) :: virtual_n, virtual_e ! arrays for storing n and e in bins that receive particles but are not yet activated, i.e. where the energy is less than e_small
      
  integer(kind=4)  :: taylor_coeff_2nd, taylor_coeff_3rd
  
  real(kind=8)     :: p_fix_ratio

  integer,allocatable, dimension(:) :: cresp_edges
!====================================================================================================
!  
 contains
!
!====================================================================================================

  subroutine init_cresp
   use constants,             only: I_ZERO, zero, ten
   use diagnostics, only: my_allocate_with_index
   implicit none
    integer                  :: i       ! enumerator
    logical, save            :: first_run = .true.

    call cresp_read_nml_module

    if (first_run .eqv. .true.) then
        if (ncre .ne. I_ZERO)  then
            print *,'[@init_cresp] Initial CRESP parameters read:'
            print '(A15, 1I3)'   ,'ncre        = ', ncre
            print '(A15, 1E15.7)','p_min_fix   = ', p_min_fix
            print '(A15, 1E15.7)','p_max_fix   = ', p_max_fix
            print '(A15, 1E15.7)','p_lo_init   = ', p_lo_init
            print '(A15, 1E15.7)','p_up_init   = ', p_up_init
            print '(A15, 1F15.7)','f_init      = ', f_init
            print '(A15, 1F15.7)','q_init      = ', q_init
            print '(A15, 1F15.7)','q_big       = ', q_big
            print '(A15, 1F15.7)','cfl_cre     = ', cfl_cre
            print '(A25, 1I3)'    ,'Taylor expansion order =', expan_order
            print '(A15, 10E15.7)','K_cre_paral_1=', K_cre_paral_1
            print '(A15, 10E15.7)','K_cre_perp_1 =', K_cre_perp_1
            print '(A15, 10E15.7)','K_cre_pow    =', K_cre_pow
            print '(A15, 1E15.7)', 'e_small      =', e_small
            print '(A20, A5)',    'initial_condition =' , initial_condition
            print '(A17, 1E15.7)','bump amplitude =', bump_amp
            print '(A27, I1)','e_small_approx_init_cond =', e_small_approx_init_cond
            print '(A22, I1)','e_small_approx_p_lo =', e_small_approx_p_lo
            print '(A22, I1)','e_small_approx_p_up =', e_small_approx_p_up
            print '(A22, I1)','add_spectrum_base   =', add_spectrum_base
            print '(A22, 1F15.7 )','max_p_ratio =', max_p_ratio
            print '(A22, L2 )' ,'force_init_NR   = ', force_init_NR
            print '(A22, I4)', 'NR_iter_limit  = ', NR_iter_limit
            print '(A22, 1E15.7)','epsilon(eps) = ', eps
   
            if ( (e_small_approx_p_lo+e_small_approx_p_up) .gt. 0 .and. e_small_approx_init_cond .lt. 1) then 
                e_small_approx_init_cond = 1  ! 
                print *, "[WARNING] @initcrspectrum: chosen approximation of boundary momenta"
                print *, " >>> Modifying e_small_approx_init_cond to 1 "
                call sleep(1)
            endif
! countermeasure - in case unrecognized or invalid parameters are provided
            if ( e_small_approx_p_lo .gt. 0 ) then ; e_small_approx_p_lo = 1 ; else ; e_small_approx_p_lo = 0 ; endif
            if ( e_small_approx_p_up .gt. 0 ) then ; e_small_approx_p_up = 1 ; else ; e_small_approx_p_up = 0 ; endif
            if ( e_small_approx_init_cond .gt. 0 ) then ; e_small_approx_init_cond = 1 ; else ; e_small_approx_init_cond = 0 ; endif
            
            if (e_small_approx_p_lo+e_small_approx_p_up .eq. 0) then
                NR_refine_solution_q = .true. ! for testing we leave precise solutions of q (especially for outer momenta)
                e_small = zero                ! no threshold energy for bin activation necessary
            endif
! arrays initialization
            call my_allocate_with_index(p_fix,ncre,0)
            call my_allocate_with_index(p_mid_fix,ncre,1)
            call my_allocate_with_index(cresp_edges,ncre,0)

            cresp_edges = (/ (i,i=0,ncre) /)
            p_fix = zero 
            w  = (log10(p_max_fix/p_min_fix))/real(ncre-2,kind=8)
            p_fix(1:ncre-1)  =  p_min_fix*ten**(w* real((cresp_edges(1:ncre-1)-1),kind=8) )
            p_fix(0)    = zero
            p_fix(ncre) = zero
            p_fix_ratio = ten**w
            write (*,'(A22, 50E15.7)') 'Fixed momentum grid: ', p_fix
            write (*,'(A22, 50E15.7)') 'Bin p-width (log10): ', w

            p_mid_fix = 0.0
            p_mid_fix(2:ncre-1) = sqrt(p_fix(1:ncre-2)*p_fix(2:ncre-1))
            p_mid_fix(1)    = p_mid_fix(2) / p_fix_ratio
            p_mid_fix(ncre) = p_mid_fix(ncre-1) * p_fix_ratio
            write (*,'(A22, 50E15.7)') 'Fixed mid momenta:   ',p_mid_fix(1:ncre)

! Input parameters check
            else
                write (*,"(A10,I4,A96)") 'ncre   = ', ncre, &
                      '; cr-electrons NOT initnialized. If COSM_RAY_ELECTRONS flag is on, please check your parameters.'
                stop
            endif
! Description of initial_condition - spectrum types: powl - pure power-law like, brpl - broken power-law like (with break in the first bin, making it easier for NR algorithm to find solution of lower cutoff momentum), bump - gaussian-like spectrum, syme - symmetric energy distribution relative to the middle of the initial spectrum, symf - similar, but symmetric in distribution function
            if (initial_condition .ne. 'powl' .and. initial_condition .ne. 'bump' .and. initial_condition .ne. 'brpl' &
                        .and. initial_condition .ne. 'symf' .and. initial_condition .ne. 'syme' ) then
                write(*,"(A41,A4,A47)") "Provided unrecognized initial_condition (",initial_condition,&
                                                "). Make sure that value is correctly provided."
                stop
            endif
            if ( f_init .lt. eps) then
                if (initial_condition == 'powl' .or. initial_condition == 'brpl') then
                write (*,"(A34,A4,A63)") "Provided power law type spectrum (",initial_condition &
                    ,") with initial amplitude f_init ~ zero. Check your parameters."
                stop
            endif
        endif
 
        if (initial_condition == 'bump' .and. bump_amp .lt. eps ) then
            write (*,"(A73,E16.8,A27)") "Provided gaussian type energy spectrum with initial amplitude bump_amp =",bump_amp, &
                                 "~ 0. Check your parameters."
            stop
        endif
        if (add_spectrum_base .gt. 0 ) then
            write (*,"(A73,E16.8,A27)") "add_spectrum_base is nonzero -> will assure energy .ge. e_small at initialization"
            if (add_spectrum_base .ne. 1) then
                add_spectrum_base = 1
            endif
        else
            write (*,"(A73,E16.8,A27)") "add_spectrum_base is zero -> will NOT assure energy.ge. e_small at initialization"
            if (add_spectrum_base .ne. 0) then
                add_spectrum_base = 0
            endif
        endif

        taylor_coeff_2nd = int(mod(2,expan_order) / 2 + mod(3,expan_order),kind=2 ) ! coefficient which is always equal 1 when order =2 or =3 and 0 if order = 1
        taylor_coeff_3rd = int((expan_order - 1)*(expan_order- 2) / 2,kind=2)        ! coefficient which is equal to 1 only when order = 3
        print '(A15, 2I3)', 'Taylor coeffs:', taylor_coeff_2nd, taylor_coeff_3rd
        call init_cresp_types
      endif
 end subroutine init_cresp
!---------------------------------------------------------------------------------------------------- 
 subroutine cresp_read_nml_module
 implicit none
    namelist /COSMIC_RAY_SPECTRUM/ cfl_cre, p_lo_init, p_up_init, f_init, q_init, q_big, ncre, initial_condition, &
    &                         p_min_fix, p_max_fix, cre_eff, K_cre_paral_1, K_cre_perp_1, &
    &                         K_cre_pow, expan_order, e_small, bump_amp, &
    &                         e_small_approx_init_cond, e_small_approx_p_lo, e_small_approx_p_up, force_init_NR,&
    &                         NR_iter_limit, max_p_ratio, add_spectrum_base !, arr_dim
           
        open(unit=101, file="problem.par", status="unknown")
        read(unit=101, nml=COSMIC_RAY_SPECTRUM)
        close(unit=101)
 end subroutine cresp_read_nml_module
!----------------------------------------------------------------------------------------------------
 subroutine init_cresp_types
  use diagnostics,   only: my_allocate_with_index
  use constants,               only: zero, I_ZERO
  implicit none
    call my_allocate_with_index(crel%p,ncre,0)
    call my_allocate_with_index(crel%f,ncre,0)
    call my_allocate_with_index(crel%q,ncre,1)
    call my_allocate_with_index(crel%e,ncre,1)
    call my_allocate_with_index(crel%n,ncre,1)
    crel%p = zero
    crel%q = zero
    crel%f = zero
    crel%e = zero
    crel%n = zero
    crel%i_lo = I_ZERO
    crel%i_up = I_ZERO
 end subroutine init_cresp_types
!----------------------------------------------------------------------------------------------------
 subroutine cleanup_cresp_virtual_en_arrays
  use diagnostics, only: my_deallocate
  implicit none
        if (allocated(virtual_e)) call my_deallocate(virtual_e)
        if (allocated(virtual_n)) call my_deallocate(virtual_n)

        if (allocated(p_fix)) call my_deallocate(p_fix)
        if (allocated(p_mid_fix)) call my_deallocate(p_mid_fix)
        if (allocated(cresp_edges)) call my_deallocate(cresp_edges)
        if (allocated(crel%p))  call my_deallocate(crel%p)
        if (allocated(crel%f)) call my_deallocate(crel%f)
        if (allocated(crel%q)) call my_deallocate(crel%q)
        if (allocated(crel%e)) call my_deallocate(crel%e)
        if (allocated(crel%n)) call my_deallocate(crel%n)
 end subroutine cleanup_cresp_virtual_en_arrays
!----------------------------------------------------------------------------------------------------
end module initcrspectrum
