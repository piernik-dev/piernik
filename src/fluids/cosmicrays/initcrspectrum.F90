module initcrspectrum  
! pulled by COSM_RAY_ELECTRONS

! contains routines reading namelist in problem.par file dedicated to cosmic ray electron spectrum and initializes types used.
! present in namelist COSMIC_RAY_SPECTRUM
 integer(kind=4)    :: ncre !       = 5
 real(kind=8)       :: p_min_fix ! = 1.0e2
 real(kind=8)       :: p_max_fix ! = 1.0e5
 real(kind=8)       :: f_init    ! = 1.0
 real(kind=8)       :: q_init    ! = 5.0
 real(kind=8)       :: q_big     ! maximal amplitude of q
 real(kind=8)       :: p_lo_init ! = 1.2e1
 real(kind=8)       :: p_up_init ! = 1.3e3
 real(kind=8)       :: cfl_cre   ! = 1.0
 real(kind=8)       :: cre_eff
 real(kind=8)       :: K_cre_paral_1 ! = 0
 real(kind=8)       :: K_cre_perp_1 ! = 0
 real(kind=8)       :: K_cre_pow    ! = 0 
 integer(kind=4)    :: expan_order  ! = 1
 real(kind=8)       :: e_small      = 1.0e-5
 real(kind=8)       :: bump_amp
 character(len=4)   :: initial_condition
 integer(kind=1)    :: e_small_approx_p_lo
 integer(kind=1)    :: e_small_approx_p_up
 integer(kind=1)    :: e_small_approx_init_cond
 integer(kind=1)    :: add_spectrum_base = 1
 real(kind=8)       :: max_p_ratio = 3.0
 integer(kind=2)    :: NR_iter_limit=300
!  integer(kind=2)    :: arr_dim
 logical            :: force_init_NR = .false.
 
 real(kind=8),allocatable, dimension(:) :: p_fix
 real(kind=8)       :: w
 
 real(kind=8)     , parameter      :: eps  = 1.0e-15  ! epsilon parameter for real comparisons
 real(kind=8)     , parameter      :: small_eps = eps*10e-10 ! a very small parameter, for comparisons with reals that should be greater than zero
!-----------
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

  type spec_mod_trms
   real(kind=8)    :: ub
   real(kind=8)    :: ud
   real(kind=8)    :: ucmb
  end type spec_mod_trms
  
  real(kind=8), allocatable, dimension(:,:,:,:) :: virtual_n, virtual_e ! arrays for storing n and e in bins that receive particles but are not yet activated
      
  integer(kind=4)  :: taylor_coeff_2nd, taylor_coeff_3rd

  integer,allocatable, dimension(:) :: cresp_edges
!----------------------------------------------  
!  
contains
!
!----------------------------------------------

  subroutine init_cresp
   use constants,             only: I_ZERO, zero, ten !, I_ONE
   use cresp_arrays_handling, only: allocate_with_index
   implicit none
    integer                  :: i       ! enumerator
    logical                  :: first_run = .true.
    
    namelist /COSMIC_RAY_SPECTRUM/ cfl_cre, p_lo_init, p_up_init, f_init, q_init, q_big, ncre, initial_condition, &
           &                         p_min_fix, p_max_fix, cre_eff, K_cre_paral_1, K_cre_perp_1, K_cre_pow, &
           &                         expan_order, e_small, bump_amp, &
           &                         e_small_approx_init_cond, e_small_approx_p_lo, e_small_approx_p_up, force_init_NR,&
           &                         NR_iter_limit, max_p_ratio, add_spectrum_base !, arr_dim
  

        open(unit=101, file="problem.par", status="unknown")
        read(unit=101, nml=COSMIC_RAY_SPECTRUM)
        close(unit=101)
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
            print '(A20, A5)',     'initial_condition =' , initial_condition
            print '(A17, 1E15.7)','bump amplitude =', bump_amp
            print '(A27, I1)','e_small_approx_init_cond =', e_small_approx_init_cond
            print '(A22, I1)','e_small_approx_p_lo =', e_small_approx_p_lo
            print '(A22, I1)','e_small_approx_p_up =', e_small_approx_p_up
            print '(A22, I1)','add_spectrum_base   =', add_spectrum_base
            print '(A22, 1F15.7 )','max_p_ratio =', max_p_ratio
            print '(A22, L2 )' ,'force_init_NR   = ', force_init_NR
            print '(A22, I4)', 'NR_iter_limit  = ', NR_iter_limit
        !    print '(A22,I4)',  'grids_dim      = ', arr_dim
   
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
! arrays initialization
            call allocate_with_index(p_fix,0,ncre)
            call allocate_with_index(cresp_edges ,0,ncre)

            cresp_edges = (/ (i,i=0,ncre) /)
            p_fix = zero 
            w  = (log10(p_max_fix/p_min_fix))/real(ncre-2,kind=8)
            p_fix(1:ncre-1)  =  p_min_fix*ten**(w* real((cresp_edges(1:ncre-1)-1),kind=8) )
            p_fix(0)    = zero
            p_fix(ncre) = zero
            write (*,'(A22, 50E15.7)') 'Fixed momentum grid: ', p_fix
            write (*,'(A22, 50E15.7)') 'Bin p-width (log10): ', w
! Input parameters check
            else
                write (*,"(A10,I4,A96)") 'ncre   = ', ncre, &
                      '; cr-electrons NOT initnialized. If COSM_RAY_ELECTRONS flag is on, please check your parameters.'
                stop
            endif
  
            if (initial_condition .ne. 'powl' .and. initial_condition .ne. 'bump' .and. initial_condition .ne. 'brpl' &
                        .and. initial_condition .ne. 'symf' .and. initial_condition .ne. 'syme' ) then
                write(*,"(A41,A4,A47)") "Provided unrecognized initial_condition (",initial_condition,&
                                                "). Make sure that value is correctly provided."
                stop
            endif
            if ( f_init .lt. small_eps) then
                if (initial_condition == 'powl' .or. initial_condition == 'brpl') then
                write (*,"(A34,A4,A63)") "Provided power law type spectrum (",initial_condition &
                    ,") with initial amplitude f_init ~ zero. Check your parameters."
                stop
            endif
        endif
 
        if (initial_condition == 'bump' .and. bump_amp .lt. small_eps ) then
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
 subroutine init_cresp_types
  use cresp_arrays_handling,   only: allocate_with_index
  use constants,               only: zero, I_ZERO !, I_ONE, zero, one
  implicit none
    call allocate_with_index(crel%p,0,ncre)
    call allocate_with_index(crel%f,0,ncre)
    call allocate_with_index(crel%q,1,ncre)
    call allocate_with_index(crel%e,1,ncre)
    call allocate_with_index(crel%n,1,ncre)
    crel%p = zero
    crel%q = zero
    crel%f = zero
    crel%e = zero
    crel%n = zero
    crel%i_lo = I_ZERO
    crel%i_up = I_ZERO
 end subroutine init_cresp_types
!----------------------------------------------------------------------------------------------------
 function compute_K(K0, alpha, alpha_n, length) ! /deprecated
 implicit none
    real(kind=8) :: K0, alpha, alpha_n
    integer :: length, i
    real(kind=8), dimension(length) :: compute_K
 
        do i = 1, length
            compute_K(i)   = K0*(p_min_fix*10.0**(w*real((cresp_edges(i-1)-1),kind=8))/p_max_fix) ** (alpha - alpha_n)
        enddo 
!  print *, '@compute_K = ', compute_K
 end function compute_K
!----------------------------------------------------------------------------------------------------
 subroutine cleanup_cresp_virtual_en_arrays
  use diagnostics, only: my_deallocate ! uncomment in PIERNIK
  use cresp_arrays_handling, only: deallocate_with_index
  implicit none
        if (allocated(virtual_e)) call my_deallocate(virtual_e)
        if (allocated(virtual_n)) call my_deallocate(virtual_n)
        
        if (allocated(p_fix)) call deallocate_with_index(p_fix)
        if (allocated(cresp_edges)) call my_deallocate(cresp_edges)
        if (allocated(crel%p))  call deallocate_with_index(crel%p)
        if (allocated(crel%f)) call deallocate_with_index(crel%f)
        if (allocated(crel%q)) call deallocate_with_index(crel%q)
        if (allocated(crel%e)) call deallocate_with_index(crel%e)
        if (allocated(crel%n)) call deallocate_with_index(crel%n)
 end subroutine cleanup_cresp_virtual_en_arrays

end module initcrspectrum
