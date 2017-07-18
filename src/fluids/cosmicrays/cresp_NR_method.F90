module cresp_NR_method
! pulled by COSM_RAY_ELECTRONS
 use initcrspectrum, only: max_p_ratio, eps !, arr_dim

 implicit none
  integer, parameter :: ndim = 2
  real(kind=8), dimension(:), allocatable :: p_space, q_space
  integer(kind=2), parameter :: ratio_ord_min = -6
  integer(kind=2), parameter :: ratio_ord_max = 1
  integer(kind=2), parameter :: p_p_dec       = 12
  integer(kind=2), parameter :: arr_dim = int(p_p_dec * (ratio_ord_max - ratio_ord_min))
  real(kind=8) :: alpha, n_in, e_in, p_im1, p_ip1
  real(kind=8), dimension(1:arr_dim) :: alpha_tab_lo, alpha_tab_up, n_tab_lo, n_tab_up
  real(kind=8), dimension(1:arr_dim,1:arr_dim) :: p_ratios_lo, f_ratios_lo, p_ratios_up, f_ratios_up
  integer(kind=2) :: ii
  real(kind=8) :: eps_det = eps * 1.0e-10
    
  abstract interface
    function arbitrary_function(z)
      real(kind=8),dimension(2) :: arbitrary_function
      real(kind=8),dimension(2) :: z
    end function arbitrary_function
  end interface
  procedure (arbitrary_function), pointer :: selected_function => null()
 
  contains !  -------*--------
  
  subroutine NR_algorithm(x,exit_code)
  use initcrspectrum, only: NR_iter_limit
  implicit none
    real(kind=8),dimension(ndim),intent(inout) :: x
    real(kind=8),dimension(ndim)   :: fun_vec_value
    real(kind=8),dimension(1:ndim) :: cor
    real(kind=8), dimension(size(x),size(x)) :: fun_vec_jac, fun_vec_inv_jac
    real(kind=8) :: err_f, err_x
    real(kind=8) :: det
    logical,intent(out)  :: exit_code
    integer(kind=2) :: i !, j
    
    err_f = 1.0e-15
    err_x = 1.0e-15
    exit_code=.true.

    fun_vec_value = selected_function(x)
    if (maxval(abs(fun_vec_value)) < 0.01 * err_f) then ! in case when f converges at initialization
            exit_code=.false.
!             write(*,"(A33,2E22.15)")"Convergence (f) at initialization", x
            return
    endif
!   do j = 1, 3
    do i = 1, NR_iter_limit                                 ! it is not possible to find solution with demanded precision.
        if (maxval(abs(fun_vec_value)) < err_f ) then    ! For convergence via value of f
            exit_code=.false. 
#ifdef VERBOSE
            write(*,"(A47,I4,A12)",advance="no") "Convergence via value of fun_vec_value after ",i, " iterations."!, x, fun_vec_value
#endif /* VERBOSE */
            ii = max(i, ii)
            return
        endif

        fun_vec_value = selected_function(x)
        fun_vec_jac = jac_fin_diff(x)                    ! function vector already explicitly provided to jac_fin_diff (finite difference method)

        det = determinant_2d_real(fun_vec_jac)           ! WARNING - algorithm is used for ndim = 2. For more dimensions LU or other methods should be implemented.
        if (abs(det) .lt. eps_det) then              ! Countermeasure against determinant = zero
!             write (*,"(A20)") "WARNING: det ~ 0.0"
            exit_code = .true. 
            return
        endif
        fun_vec_inv_jac = invert_2d_matrix(fun_vec_jac,det)

        cor(1) = fun_vec_inv_jac(1,1) *fun_vec_value(1) + fun_vec_inv_jac(1,2) * fun_vec_value(2)
        cor(2) = fun_vec_inv_jac(2,1) *fun_vec_value(1) + fun_vec_inv_jac(2,2) * fun_vec_value(2)
        x = x+cor
!  write(*,'(A20, 2E35.25, A5, 2E22.14)') "Obtained values (x): " , x,' | ', sum(abs(cor)), sum(abs(fun_vec_value))! ,maxval(abs(fun_vec_value)), maxval(abs(cor)),
        if (maxval(abs(cor)) < err_x) then                 ! For convergence via value of correction (cor) table.
#ifdef VERBOSE
           write(*,"(A47,I4,A12)",advance="no") "Convergence via value of cor array     after ",i," iterations."
#endif /* VERBOSE */
           exit_code = .false.
           ii = max(i, ii)
           return
        endif 
    enddo
!   err_f = 5.0*err_f
!   err_x = 5.0*err_x ! changing tolerance so that more solutions can be found
!   enddo
  
!     write(*,"(A45,I4,A24)") "  ... WARNING: Maximum number of iterations (",NR_iter_limit,") exceeded @global_newt!"
    exit_code = .true.

  end subroutine NR_algorithm

!----------------------------------------------------------------------------------------------------
  subroutine cresp_initialize_guess_grids
   use initcrspectrum, only: eps
   use constants,      only: zero
   implicit none
    logical   :: first_run = .true.
        
        if (first_run .eqv. .true. ) then
    
            if (.not. allocated(p_space))     allocate(p_space(1:arr_dim)) ! these will be deallocated once initialization is over
            if (.not. allocated(q_space))     allocate(q_space(1:arr_dim)) ! these will be deallocated once initialization is over

 
            p_ratios_up = zero ; f_ratios_up = zero
            p_ratios_lo = zero ; f_ratios_lo = zero
            
!             call sleep(1)
            call fill_guess_grids
        
            print *, "Are there zeros? (q_ratios)",    count(q_space.lt.eps)
            print *, "Are there zeros? (p_ratios_up)", count(p_ratios_up.lt.eps)
            print *, "Are there zeros? (f_ratios_up)", count(f_ratios_up.lt.eps)
            print *, "Are there zeros? (p_ratios_lo)", count(p_ratios_lo.lt.eps)
            print *, "Are there zeros? (f_ratios_lo)", count(f_ratios_lo.lt.eps)
            print *, "Count of array elements:", size(p_ratios_lo)
            print *,"----------"
!   
            if (allocated(p_space)) deallocate(p_space) ! only needed at initialization
            if (allocated(q_space)) deallocate(q_space)
            
            first_run = .false.
        endif
! 
  end subroutine cresp_initialize_guess_grids
!----------------------------------------------------------------------------------------------------
  subroutine fill_guess_grids
   use constants, only: zero, one, half, ten
   use initcrspectrum,  only: q_big, force_init_NR
   implicit none
    integer(kind=2)  :: i, j, int_logical_p, int_logical_f !, ierr
    logical  :: exit_code
    real(kind=8) :: a_min_lo=huge(one), a_max_lo=tiny(one), a_min_up=huge(one), a_max_up=tiny(one),& 
                    n_min_lo=huge(one), n_max_lo=tiny(one), n_min_up=huge(one), n_max_up=tiny(one)
        q_space = zero
        do i=1, int(half*arr_dim)
            q_space(i) = ln_eval_array_val(i, q_big, real(0.05,kind=8), int(1,kind=2), int(half*arr_dim,kind=2))
        enddo
        do i= 1, int(half*arr_dim)!, arr_dim
            q_space(int(half*arr_dim)+i) = -q_space(int(half*arr_dim)+1-i)
        enddo
        
! setting up a grids of ratios to be used as phase space for NR tabs, obtained later
        do i = 1, arr_dim
            p_space(i) =  max_p_ratio**(real(i)/real(arr_dim))
        enddo
        do i=1, arr_dim
            do j = 1,arr_dim
                a_min_lo = min(a_min_lo, abs(encp_func_2_zero_lo(p_space(i),zero, q_space(j))))
                n_min_lo = min(n_min_lo, abs(n_func_2_zero_lo(p_space(i), zero, q_space(j))))
                a_min_up = min(a_min_up, abs(encp_func_2_zero_up(p_space(i), zero ,q_space(j))))
                n_min_up = min(n_min_up, abs(n_func_2_zero_up(p_space(i),p_space(i)**(-q_space(j)), zero ,q_space(j))))

                a_max_lo = max(a_max_lo, abs(encp_func_2_zero_lo(p_space(i), zero, q_space(j))))
                n_max_lo = max(n_max_lo, abs(n_func_2_zero_lo(p_space(i), zero, q_space(j))))
                a_max_up = max(a_max_up, abs(encp_func_2_zero_up(p_space(i), zero ,q_space(j))))
                n_max_up = max(n_max_up, abs(n_func_2_zero_up(p_space(i),p_space(i)**(-q_space(j)), zero ,q_space(j))))
            enddo
        enddo
        
        a_min_lo = 0.8 * a_min_lo
        a_max_lo = 0.999999 !1 * a_max_lo
        a_min_up = 1.000005 ! 0.8 * a_min_up
        a_max_up = 1.1 * a_max_up
        n_min_lo = 0.001 * n_min_lo
        n_max_lo = 1.1 * n_max_lo
        n_min_up = 0.001 * n_min_up
        n_max_up = 1.1 * n_max_up
   
!             print *," min / max values"
!             print *, a_min_lo, a_max_lo
!             print *, n_min_lo, n_max_lo
!             print *, a_min_up, a_max_up
!             print *, n_min_up, n_max_up
        
        print *,"alpha_tab_lo(i),      alpha_tab_up(i),        n_tab_lo(i),        n_tab_up(i) |       p_space(i),     q_space(i)" 
        do i=1, arr_dim
            alpha_tab_lo(i) = ln_eval_array_val(i, a_min_lo, a_max_lo, int(1,kind=2), arr_dim)
            alpha_tab_up(i) = a_min_up * ten**((log10(a_max_up/a_min_up))/real(arr_dim-1,kind=8)*real((i-1),kind=8))
            n_tab_lo(i)     = n_min_lo * ten**((log10(n_max_lo/n_min_lo))/real(arr_dim-1,kind=8)*real((i-1),kind=8))
            n_tab_up(i)     = n_min_up * ten**((log10(n_max_up/n_min_up))/real(arr_dim-1,kind=8)*real((i-1),kind=8))
        enddo
! #ifdef VERBOSE
        do i = 1, arr_dim
            print *,i,"|", alpha_tab_lo(i), alpha_tab_up(i), n_tab_lo(i), n_tab_up(i),"|", p_space(i), q_space(i), &
                            p_space(i)**(-q_space(i))
        enddo
!         call sleep(1)
        print *, "-----------"
! #endif /* VERBOSE */

        call read_NR_guess_grid(p_ratios_up, "p_ratios_up", exit_code) ;  int_logical_p = logical_2_int(exit_code)
        call read_NR_guess_grid(f_ratios_up, "f_ratios_up", exit_code) ;  int_logical_f = logical_2_int(exit_code)
        if ( int_logical_f + int_logical_p .gt. 0 .or. force_init_NR .eqv. .true.) then
! Setting up the "guess grid" for p_up case
            call fill_boundary_grid("up", p_ratios_up, f_ratios_up)
            call save_NR_guess_grid(p_ratios_up,"p_ratios_up")
            call save_NR_guess_grid(f_ratios_up,"f_ratios_up")
        else
            print *," >> Will not solve ratios table (up), reading data from file instead."
        endif

        call read_NR_guess_grid(p_ratios_lo, "p_ratios_lo", exit_code) ;   int_logical_p = logical_2_int(exit_code)
        call read_NR_guess_grid(f_ratios_lo, "f_ratios_lo", exit_code) ;   int_logical_f = logical_2_int(exit_code)
        if ( int_logical_f + int_logical_p .gt. 0 .or. force_init_NR .eqv. .true.) then
! Setting up the "guess grid" for p_lo case
            call fill_boundary_grid("lo", p_ratios_lo, f_ratios_lo)
            call save_NR_guess_grid(p_ratios_lo,"p_ratios_lo")
            call save_NR_guess_grid(f_ratios_lo,"f_ratios_lo")
        else
            print *," >> Will not solve ratios table (lo), reading data from file instead."
        endif        
    
        print *, "max values:", maxval(p_ratios_up,MASK=p_ratios_up.gt.0.0), maxval(f_ratios_up,MASK=f_ratios_up.gt.0.0)
        print *, "min values:", minval(p_ratios_up,MASK=p_ratios_up.gt.0.0), minval(f_ratios_up,MASK=f_ratios_up.gt.0.0)
    
        print *, "max values:", maxval(p_ratios_lo,MASK=p_ratios_lo.gt.0.0), maxval(f_ratios_lo,MASK=f_ratios_lo.gt.0.0)
        print *, "min values:", minval(p_ratios_lo,MASK=p_ratios_lo.gt.0.0), minval(f_ratios_lo,MASK=f_ratios_lo.gt.0.0)
    
        print *, "maximal number of iterations with success:", ii
    
 end subroutine fill_guess_grids
!---------------------------------------------------------------------------------------------------- 
 function ln_eval_array_val(i, arr_min, arr_max, min_i, max_i)
 implicit none
  real(kind=8) :: b, arr_min, arr_max, ln_eval_array_val
  integer(kind=2) :: i, max_i, min_i
    b = (log(real(max_i)) -log(real(min_i)))/ (arr_max - arr_min)
    ln_eval_array_val = (arr_min-log(real(min_i))/b ) + log(real(i)) / b
 end function ln_eval_array_val
!---------------------------------------------------------------------------------------------------- 
 subroutine fill_boundary_grid(bound_case, fill_p, fill_f) ! to be paralelized
  use constants, only: zero
  implicit none
  real(kind=8), dimension(1:2) :: x_vec, prev_solution
  real(kind=8), dimension(:,:) :: fill_p, fill_f
  integer(kind=4) :: i, j, is, js
  logical         :: exit_code
  character(len=2) :: bound_case ! "up" or "lo"

    prev_solution(1) = p_space(1)
    prev_solution(2) = p_space(1)**q_space(1)
    call sleep(1)
    fill_p = zero ; fill_f = zero
        do i =1, arr_dim
            do j = 1, arr_dim
                exit_code = .true.
                if (bound_case == "lo") then
                    alpha = alpha_tab_lo(i)
                    n_in  = n_tab_lo(j)
                    selected_function => fvec_lo
                else if (bound_case == "up") then
                    alpha = alpha_tab_up(i)
                    n_in  = n_tab_up(j)
                    selected_function => fvec_up
                else
                    print *, "Wrong *bound_case* argument provided, stopping"
                    stop
                endif                    
                write(*,"(A12,2I4,A9,I4,A5,2E16.9)",advance="no") "Now solving",i,j,", sized ",arr_dim ,", (alpha,n): ",alpha,n_in
                x_vec(1) = prev_solution(1)
                x_vec(2) = prev_solution(2)
                call NR_algorithm(x_vec, exit_code)
                if (exit_code .eqv. .false.) then
                    fill_p(i,j) = x_vec(1)
                    fill_f(i,j) = x_vec(2)
                    write (*, "(A51, 2E22.15)",advance="no") " -> (prev) solution obtained, (p_ratio, f_ratio) = ", x_vec
                else
                    do is =1, arr_dim
                        do js = 1, arr_dim
                            x_vec(1) = p_space(is)
                            x_vec(2) = p_space(is)**(-q_space(js))
                            if (exit_code .eqv. .true.) then
                                call NR_algorithm(x_vec, exit_code)
                                if ( exit_code .eqv. .false. ) then
                                    fill_p(i,j) = x_vec(1) ! i index - alpha, j index - n_in
                                    fill_f(i,j) = x_vec(2)
                                    prev_solution(1) = x_vec(1)
                                    prev_solution(2) = x_vec(2)
                                    write (*, "(A45, 2E22.15)",advance="no") " -> solution obtained, (p_ratio, f_ratio) = ", x_vec
                                    exit
                                endif
                            endif
                        enddo
                    enddo
                endif
                print *,""
            enddo
        enddo
        print *,""
 
 
 end subroutine
 !----------------------------------------------------------------------------------------------------
  function q_ratios(f_ratio, p_ratio)
  use constants, only: zero
  implicit none
   real(kind=8), intent(in) :: f_ratio, p_ratio
   real(kind=8) :: q_ratios
        q_ratios = zero
        q_ratios = -log10(f_ratio) / log10(p_ratio)
  end function q_ratios
!---------------------------------------------------------------------------------------------------
! Function estimating values of jacobian via finite difference method
!---------------------------------------------------------------------------------------------------
  function jac_fin_diff(x) ! jacobian via finite difference method
  implicit none
   real(kind=8), dimension(ndim) :: x, xp, xm
   real(kind=8), dimension(size(x),size(x)) :: jac_fin_diff
   real(kind=8), dimension(size(x)) :: dx
   real(kind=8), parameter          :: dx_par = 1.0e-4
   integer(kind=2) :: j
        do j = 1,size(x)
            dx(:) = min(x(:)*dx_par,dx_par) ! the value of dx is scaled not to go over value of x
            xp = x ; xm = x
            xp(j) = x(j) - dx(j) ;  xm(j) = x(j) + dx(j)
            jac_fin_diff(:,j)  = ( selected_function(xp) - selected_function(xm)) / dx(j)
        enddo
        jac_fin_diff = 0.5 * jac_fin_diff
  end function jac_fin_diff
!----------------------------------------------------------------------------------------------------
  function determinant_2d_real(matrix_2d_real)
  implicit none
  real(kind=8), dimension(2,2), intent(in) :: matrix_2d_real
  real(kind=8) :: determinant_2d_real
        determinant_2d_real = matrix_2d_real(1,1) * matrix_2d_real (2,2) - ( matrix_2d_real(2,1) * matrix_2d_real(1,2) )
                ! in case determinant_2d_real = zero, algorithm will return exit_code=1 and new attempt to solve may take place
  end function determinant_2d_real
!----------------------------------------------------------------------------------------------------
  function get_cofactor_matrix_2d_real(matrix_2d_real)
  implicit none
  integer(kind=1),parameter :: m_dim = 2
  real(kind=8), dimension(m_dim,m_dim) :: matrix_2d_real
  real(kind=8), dimension(m_dim,m_dim) :: get_cofactor_matrix_2d_real
  integer(kind=1) :: i,j
        do i = 1,m_dim
            do j = 1,m_dim
                get_cofactor_matrix_2d_real(i,j) = ( (-1)**(i+j) * matrix_2d_real(m_dim+1-i, m_dim+1-j))
            enddo
        enddo
  end function get_cofactor_matrix_2d_real
!----------------------------------------------------------------------------------------------------
  function invert_2d_matrix(matrix,determinant)
  use constants, only: one
  implicit none
    integer(kind=1),parameter :: m_dim = 2
    real(kind=8), dimension(m_dim,m_dim), intent(in) :: matrix
    real(kind=8), dimension(m_dim,m_dim) :: invert_2d_matrix
    real(kind=8) :: determinant
        invert_2d_matrix = (one / determinant) * transpose( get_cofactor_matrix_2d_real(matrix) )
  end function invert_2d_matrix
!----------------------------------------------------------------------------------------------------
 function fvec_up(x)
 implicit none
    real(kind=8), dimension(ndim) :: x
    real(kind=8), dimension(ndim) :: fvec_up
    real(kind=8) :: q_in
        x = abs(x)
        q_in     = q_ratios(x(2),x(1))
        fvec_up(1) = encp_func_2_zero_up(x(1), alpha, q_in)
        fvec_up(2) = n_func_2_zero_up(x(1), x(2), n_in, q_in)
 end function fvec_up
 
!----------------------------------------------------------------------------------------------------
  function fvec_lo(x)
  implicit none
    real(kind=8), dimension(ndim) :: x
    real(kind=8), dimension(ndim) :: fvec_lo
    real(kind=8) :: q_in
        x = abs(x)    
        q_in     = q_ratios(x(2),x(1))
        fvec_lo(1) = encp_func_2_zero_lo(x(1), alpha, q_in)
        fvec_lo(2) = n_func_2_zero_lo(x(1), n_in, q_in)
 end function fvec_lo
!----------------------------------------------------------------------------------------------------
 function encp_func_2_zero_up(p_ratio, alpha_cnst, q_in) ! from eqn. 29
 use constants,     only: one, three, four
 implicit none
    real(kind=8) :: p_ratio, q_in, alpha_cnst
    real(kind=8) :: encp_func_2_zero_up
        if (abs(q_in - three) .lt. eps) then
            encp_func_2_zero_up = -alpha_cnst + (- one + p_ratio)/log(p_ratio) 
        else if (abs(q_in - four) .lt. eps) then
            encp_func_2_zero_up = -alpha_cnst + p_ratio*log(p_ratio)/(p_ratio - one)
        else
            encp_func_2_zero_up = -alpha_cnst + ((three - q_in)/(four - q_in))*((p_ratio **(four - q_in)- one)/ &
                                    (p_ratio **(three - q_in)- one))
        endif
  end function encp_func_2_zero_up
  
!---------------------------------------------------------------------------------------------------
 function encp_func_2_zero_lo(p_ratio, alpha_cnst, q_in) ! from eqn. 29
 use constants,     only: one, three, four
 implicit none
    real(kind=8) :: p_ratio, q_in, alpha_cnst
    real(kind=8) :: encp_func_2_zero_lo
        if (abs(q_in - three) .lt. eps) then
            encp_func_2_zero_lo = -alpha_cnst + ((- one + p_ratio)/log(p_ratio)) / p_ratio
        else if (abs(q_in - four) .lt. eps) then
            encp_func_2_zero_lo = -alpha_cnst + log(p_ratio)/(p_ratio - one) / p_ratio
        else
            encp_func_2_zero_lo = -alpha_cnst + ((three - q_in)/(four - q_in))*((p_ratio **(four - q_in) - one)/ &
                                (p_ratio **(three - q_in) - one)) / p_ratio
        endif
  end function encp_func_2_zero_lo

!----------------------------------------------------------------------------------------------------  
  function n_func_2_zero_up(p_ratio, f_ratio, n_cnst, q_in) ! from eqn. 9
  use constants,       only: one, two, three
  use initcrspectrum,  only: e_small
  use cresp_variables, only: clight ! use units, only: clight
  implicit none
    real(kind=8) :: p_ratio, f_ratio, q_in, n_cnst
    real(kind=8) :: n_func_2_zero_up
        if (abs(q_in - three) .lt. eps) then
            n_func_2_zero_up = - n_cnst + e_small / ((clight **two) * f_ratio * (p_ratio **three))* log(p_ratio)
        else
            n_func_2_zero_up = - n_cnst + e_small / ((clight **two) * f_ratio * (p_ratio **three)) &
                                            * ((p_ratio **(three-q_in) - one)/(three - q_in))
        endif
  end function n_func_2_zero_up
!----------------------------------------------------------------------------------------------------
  function n_func_2_zero_lo(p_ratio, n_cnst,q_in) ! from eqn. 9
  use constants,       only: one, two, three
  use initcrspectrum,  only: e_small
  use cresp_variables, only: clight ! use units, only: clight
  implicit none
    real(kind=8) :: p_ratio,  n_cnst, q_in
    real(kind=8) :: n_func_2_zero_lo
        if (abs(q_in - three) .lt. eps) then
            n_func_2_zero_lo = - n_cnst + (e_small / (clight **two)) * log(p_ratio)
        else
            n_func_2_zero_lo = - n_cnst + (e_small / (clight **two)) * ((p_ratio **(three-q_in) - one)/(three - q_in))
        endif
  end function n_func_2_zero_lo
!====================================================================================================
! Functions below are used to solve eqns 9 and 29 for p_up and f_l or p_lo and f_r 
! (value of f_r doesn't appear in the f array of cresp_crspectrum module, but is used to estimate q).
!----------------------------------------------------------------------------------------------------
 subroutine NR_get_solution_up(x, p_l, e_input, n_input, exit_code)
 use cresp_variables, only: clight ! use units, only: clight
 use constants,       only: zero
 implicit none
  real(kind=8) :: p_l, e_input, n_input
  real(kind=8), dimension(ndim) :: x, x_init
  logical :: exit_code
  integer :: k
        alpha = e_input / (n_input * clight * p_l)
        n_in  = n_input
        p_im1 = p_l
        x_init = x
        k = 1
        print *, " Determining p_up @ NR_get_solution_up"
        selected_function => func_val_vec_up_bnd

        do while(k .le. 5)
            call NR_algorithm(x, exit_code)
            if (exit_code .eqv. .false. ) then
                print *, "  NR root search done for p_up, exit_code 0"
                alpha = zero ;  p_im1 = zero; n_in = zero
                return
            else
                print *, "   CONVERGENCE FAILURE @ NR_get_solution_up"
                x = modify_params_up(x_init, k)
            endif
        enddo
! in cases when no solutions found:
        print *, "   CONVERGENCE FAILURE @ NR_get_solution_up, assuming initial parameters."
        x = x_init  ; call sleep(1)
 end subroutine NR_get_solution_up
!----------------------------------------------------------------------------------------------------
 function func_val_vec_up_bnd(x) ! called by cresp_crspectrum module via NR_get_solution_up
 implicit none
    real(kind=8), dimension(ndim) :: x
    real(kind=8), dimension(ndim) :: func_val_vec_up_bnd
        func_val_vec_up_bnd(1) = fun3_up(x(1), x(2), p_im1, alpha)
        func_val_vec_up_bnd(2) = fun5_up(x(1), x(2), p_im1, n_in)
 end function func_val_vec_up_bnd
 
!----------------------------------------------------------------------------------------------------
  function modify_params_up(x_init, k)
   use constants, only: zero
   implicit none
   real(kind=8), dimension(ndim), intent(in) :: x_init
   real(kind=8), dimension(ndim)             :: modify_params_up
   integer, intent(inout)   :: k
        modify_params_up = zero
        write(*,'(A48,I3)')    " Modifying parameters - det = 0, counter = ", k
        modify_params_up(1) = (x_init(1) * (1.0+k*0.005))
        modify_params_up(2) = (x_init(2) * (1.0 + 0.01*k) )
        write(*,'(A45,2E17.9)')" NR_algorithm new input (log10_p_r,log10_f_l) =", modify_params_up
        k=k+1
  end function modify_params_up

!----------------------------------------------------------------------------------------------------
 subroutine NR_get_solution_lo(x, p_r, e_input, n_input, exit_code)
 use cresp_variables, only: clight ! use units, only: clight
 use constants,       only: zero
 implicit none
  real(kind=8) :: p_r, e_input, n_input
  real(kind=8), dimension(ndim) :: x, x_init
  integer :: k
  logical :: exit_code
        alpha = e_input / (n_input * clight)
        n_in  = n_input
        p_ip1 = p_r
        x_init = x
        k = 1
        print *, " Determining p_lo @ NR_get_solution_lo"
        selected_function => func_val_vec_lo_bnd
        do while(k .le. 5)
            call NR_algorithm(x, exit_code)    
            if (exit_code .eqv. .false. ) then
                print *, " NR root search done for p_lo, exit_code 0"
                alpha = zero; p_ip1 = zero ; n_in = zero
                x = abs(x)
                return
            else
                print *, "   CONVERGENCE FAILURE @ NR_get_solution_lo"
                x = modify_params_lo(x_init, k)
            endif
        enddo
! in cases when no solutions found:
        print *, "   CONVERGENCE FAILURE @ NR_get_solution_lo, assuming initial parameters."
        x = x_init  ; call sleep(1)
 end subroutine NR_get_solution_lo
 
  !----------------------------------------------------------------------------------------------------
 function func_val_vec_lo_bnd(x) ! called by cresp_crspectrum module via NR_get_solution_lo
 implicit none
    real(kind=8), dimension(ndim) :: x
    real(kind=8), dimension(ndim) :: func_val_vec_lo_bnd
        func_val_vec_lo_bnd(1) = fun3_lo(x(1), x(2), p_ip1, alpha)
        func_val_vec_lo_bnd(2) = fun5_lo(x(1), x(2), p_ip1, n_in)
 end function func_val_vec_lo_bnd

!----------------------------------------------------------------------------------------------------
 function modify_params_lo(x_init, k)
   use constants, only: zero
   implicit none
   real(kind=8), dimension(ndim), intent(in) :: x_init
   real(kind=8), dimension(ndim)             :: modify_params_lo
   integer, intent(inout)   :: k
        modify_params_lo = zero
        write(*,'(A48,I3)')    "Modifying parameters - det = 0, counter = ", k
        modify_params_lo(1) = (x_init(1)/(1 + k*0.025)  )
        modify_params_lo(2) = (x_init(2) *( 1.0 + 0.01*k ))
        write(*,'(A45,2E17.9)')"NR_algorithm new input (log10_p_l,log10_f_r) =", modify_params_lo
        k=k+1
  end function modify_params_lo

!---------------------------------------------------------------------------------------------------    
! fun3_up - Alternative function introduced to find p_up using integrals of n, e and p_fix value. Does not seek 
!        q, therefore conditions for x = (3,4) are not necessary
!---------------------------------------------------------------------------------------------------
  function fun3_up(log10_p_r, log10_f_l, p_l, alpha ) ! used by func_val_vec_up_bnd to compute upper boundary p and f.
   use constants, only: three, four, one, ten
   implicit none
    real(kind=8),intent(in)   :: log10_p_r, log10_f_l
    real(kind=8) :: f_l, f_r, alpha, fun3_up, q_bin, p_l, p_r ! sought values will be x for single N-R and x, f_l_iup in NR-2dim
        p_r = ten **log10_p_r
        f_l = ten **log10_f_l
        f_r = e_small_to_f(p_r)
        q_bin = q_ratios(f_r/f_l,p_r/p_l)
    
        if      ( abs(q_bin - three) .lt. eps) then
            fun3_up = -alpha + (-one + p_r/p_l)/log(p_r/p_l) 
        else if ( abs(q_bin - four) .lt. eps) then
            fun3_up = -alpha + (p_r/p_l)*log(p_r/p_l)/(p_r/p_l - one)
        else
            fun3_up = -alpha + ((three - q_bin) / (four - q_bin)) * &
               (((p_r/p_l)**((four - q_bin)) - one ) / ((p_r/p_l)**((three - q_bin)) - one ))
        endif
   end function fun3_up
!---------------------------------------------------------------------------------------------------    
! fun5_up - Function similar to nq_to_f subroutine, used by compute_fp_NR_2dim to estimate value of p_up and f_l_iup
!---------------------------------------------------------------------------------------------------
  function fun5_up(log10_p_r, log10_f_l, p_l, n_in) ! used by func_val_vec_up_bnd to compute upper boundary p and f.
   use constants, only: three, ten, one, fpi
   implicit none
    real(kind=8),intent(in) :: log10_p_r, log10_f_l
    real(kind=8) :: f_l, f_r, n_in, q_bin, fun5_up, p_l, p_r
        p_r = ten **log10_p_r
        f_l = ten **log10_f_l
        f_r = e_small_to_f(p_r)
        q_bin = q_ratios(f_r/f_l,p_r/p_l)

        if (abs(q_bin - three) .lt. eps) then
            fun5_up = - f_l + n_in/((fpi * p_l **three) * log(p_r/p_l))
        else
            fun5_up = - f_l + ( n_in / (fpi * p_l **three) ) * ((three - q_bin) / ((p_r/p_l)**(three - q_bin) - one) )
        endif
  end function fun5_up
  
!---------------------------------------------------------------------------------------------------    
! fun3_lo - Alternative function introduced to find p_up using integrals of n, e and p_fix value. Does not seek 
!        q, therefore conditions for x = (3,4) are not necessary
!---------------------------------------------------------------------------------------------------
  function fun3_lo(log10_p_l, log10_f_r, p_r, alpha ) ! used by func_val_vec_lo_bnd to compute low boundary p and f.
   use constants, only: three, four, one, ten
   implicit none
    real(kind=8),intent(in)   :: log10_p_l, log10_f_r
    real(kind=8) :: f_l, f_r, alpha, fun3_lo, q_bin, p_r, p_l ! sought values will be x for single N-R and x, f_l_iup in NR-2dim
        p_l = ten **log10_p_l
        f_r = ten **log10_f_r
        f_l = e_small_to_f(p_l)
        q_bin = q_ratios(f_r/f_l,p_r/p_l)
    
        if   ( abs(q_bin - three) .lt. eps) then
            fun3_lo = -alpha/p_l + (-one + p_r/p_l)/log(p_r/p_l) 
        else if ( abs(q_bin - four) .lt. eps) then
            fun3_lo = -alpha/p_l + (p_r/p_l)*log(p_r/p_l)/(p_r/p_l - one)
        else
            fun3_lo = -alpha/p_l + ((three - q_bin) / (four - q_bin)) * &
               (((p_r/p_l)**((four - q_bin)) -one ) / ((p_r/p_l)**((three - q_bin)) -one ))
        endif
  end function fun3_lo

!---------------------------------------------------------------------------------------------------    
! fun5_lo - Function similar to nq_to_f subroutine, used by compute_fp_NR_2dim to estimate value of p_lo and f_r_lo
!---------------------------------------------------------------------------------------------------
  function fun5_lo(log10_p_l, log10_f_r, p_r, n_in) ! used by func_val_vec_lo_bnd to compute low boundary p and f.
   use constants, only: three, one, fpi, ten
   implicit none
    real(kind=8),intent(in) :: log10_p_l, log10_f_r
    real(kind=8) :: f_l, f_r, n_in, q_bin, fun5_lo, p_l, p_r
        p_l = ten **log10_p_l
        f_r = ten **log10_f_r
        f_l = e_small_to_f(p_l)
        q_bin = q_ratios(f_r/f_l,p_r/p_l)
        
        if (abs(q_bin - three) .lt. eps) then
            fun5_lo = - f_l + n_in/((fpi * p_l **three) * log(p_r/p_l))
        else
            fun5_lo = - f_l + ( n_in / (fpi * p_l **three) ) * ( (three - q_bin) / ((p_r/p_l)**(three - q_bin) - one) )
        endif
  end function fun5_lo
!----------------------------------------------------------------------------------------------------
! Here - relaying e_small to f via its relation with momentum
!----------------------------------------------------------------------------------------------------
  function e_small_to_f(p_outer) ! used by variety of procedures and functions
   use initcrspectrum,  only: e_small
   use cresp_variables, only: clight ! use units, only: clight
   use constants,       only: zero, three, two, fpi
   implicit none
    real (kind=8) :: e_small_to_f, p_outer
        e_small_to_f = zero
        e_small_to_f = e_small / (fpi * (clight **two)  * p_outer **three)

  end function e_small_to_f
!----------------------------------------------------------------------------------------------------
  subroutine find_both_indexes_lo(l_min1, l_min2, alpha_val, n_val, l_panic, exit_code)
  use initcrspectrum, only: eps
  implicit none
    integer(kind=2), dimension(1:2) :: l_min1, l_min2, l_panic
    integer(kind=2) :: i, j, i_prev
    real(kind=8)    :: min_alpha, min_n, min2_alpha, min2_n, d_alpha, d_n
    real(kind=8)    :: alpha_val, n_val
    logical :: exit_code
        min_alpha = huge(min_alpha); min_n = huge(min_n) ; min2_alpha = -huge(min2_alpha); min2_n = -huge(min2_n)
        i_prev = 0
        l_min1 = 1 ; l_min2 = 1
        do i = 1, arr_dim
            do j = 1, arr_dim
                if (p_ratios_lo(i,j) .gt. eps ) then
                    d_n     = n_tab_lo(j) - n_val
                    if ( i_prev .ne. i) then ! already checked, this omition should save up to 1/3 time needed for search each time this function is called
                        d_alpha = alpha_tab_lo(i) - alpha_val
                        call get_index_min_hi( d_alpha, min_alpha,  i, l_min1(1))
                        call get_index_min_lo( d_alpha, min2_alpha, i, l_min2(1))
                    endif
                    call get_index_min_hi(d_n, min_n,  j, l_min1(2))
                    call get_index_min_lo(d_n, min2_n, j, l_min2(2))
                    i_prev = i
                endif
            enddo
        enddo
        if (l_min1(1) .eq. l_min2(1) .or. l_min1(2) .eq. l_min2(2)) then
            exit_code = .true.
            l_panic(1) = find_indexes_panic(l_min1(1), l_min2(1), min_alpha, min2_alpha)
            l_panic(2) = find_indexes_panic(l_min1(2), l_min2(2), min_n, min2_n)
        endif
  end subroutine find_both_indexes_lo
!----------------------------------------------------------------------------------------------------
  subroutine find_both_indexes_up(l_min1, l_min2, alpha_val, n_val, l_panic, exit_code)
  use initcrspectrum, only: eps
  implicit none
    integer(kind=2), dimension(1:2) :: l_min1, l_min2, l_panic ! sets of indexes to point to p_ratios_.. and f_ratios_.. arrays via alpha and n arrays.
    integer(kind=2) :: i, j, i_prev
    real(kind=8)    :: min_alpha, min_n, min2_alpha, min2_n, d_alpha, d_n
    real(kind=8)    :: alpha_val, n_val
    logical :: exit_code
        min_alpha = huge(min_alpha); min_n = huge(min_n) ; min2_alpha = -huge(min2_alpha); min2_n = -huge(min2_n)
        i_prev = 0
        l_min1 = 1 ; l_min2 = 1
        do i = 1, arr_dim
            do j = 1, arr_dim
                if (p_ratios_up(i,j) .gt. eps ) then
                    d_n     = n_tab_up(j) - n_val
                    if ( i_prev .ne. i) then ! already checked, this omition should save up to 1/3 time needed for search each time this function is called
                    d_alpha = alpha_tab_up(i) - alpha_val
                        call get_index_min_hi( d_alpha, min_alpha,  i, l_min1(1))
                        call get_index_min_lo( d_alpha, min2_alpha, i, l_min2(1))
                    endif
                    call get_index_min_hi(d_n, min_n,  j, l_min1(2))
                    call get_index_min_lo(d_n, min2_n, j, l_min2(2))
                endif
                i_prev = i ! already got result for current i
            enddo
        enddo
        if (l_min1(1) .eq. l_min2(1) .or. l_min1(2) .eq. l_min2(2)) then
            exit_code = .true.
            l_panic(1) = find_indexes_panic(l_min1(1), l_min2(1), min_alpha, min2_alpha)
            l_panic(2) = find_indexes_panic(l_min1(2), l_min2(2), min_n, min2_n)
        endif
  end subroutine find_both_indexes_up
!----------------------------------------------------------------------------------------------------
  subroutine get_index_min_hi(delta, current_min_val, current_index, min_index)
  use constants, only: zero
  implicit none
    real(kind=8)    :: current_min_val, delta
    integer(kind=2) :: current_index, min_index
        min_index = min_index
        if (delta .lt. current_min_val .and. delta .gt. zero) then
            current_min_val = delta
            min_index = current_index
        endif
  end subroutine get_index_min_hi
!----------------------------------------------------------------------------------------------------
  subroutine get_index_min_lo(delta, current_min_val, current_index, min_index)
  use constants, only: zero
  implicit none
    real(kind=8)    :: current_min_val, delta
    integer(kind=2) :: current_index, min_index
        min_index = min_index
        if (delta .gt. current_min_val .and. delta .lt. zero) then
            current_min_val = delta
            min_index = current_index
        endif
  end subroutine get_index_min_lo            
!----------------------------------------------------------------------------------------------------
  function find_indexes_panic(index_hi, index_lo, min_pos, min_neg)
  implicit none
    integer(kind=2) :: index_hi, index_lo
    integer(kind=2) :: find_indexes_panic
    real(kind=8)    :: min_pos, min_neg
        if ( min_pos .lt. abs(min_neg) ) then
            find_indexes_panic = index_hi
        else
            find_indexes_panic = index_lo
        endif
  end function find_indexes_panic  
!----------------------------------------------------------------------------------------------------
  function bl_interpol(y11,y12,y21,y22,t,u) ! for details see paragraph "Bilinear interpolation" in Numerical Recipes for F77, page 117, eqn. 3.6.5
  use constants, only: one
   implicit none
    real(kind=8)  :: y11, y12, y21, y22, t, u ! y** - tabularized values of interpolated function, t, u - coefficients
    real(kind=8)  :: bl_interpol
        bl_interpol = (one - t)*(one - u) * y11 + t*(one - u)*y12 + (one - t)*u*y21 + t*u*y22
  end function bl_interpol
!----------------------------------------------------------------------------------------------------
  function bl_in_tu(val_left, val_mid, val_right) ! for details see paragraph "Bilinear interpolation" in Numerical Recipes for F77, page 117, eqn. 3.6.4
   implicit none
    real(kind=8)  :: val_left, val_right, val_mid
    real(kind=8)  :: bl_in_tu
        bl_in_tu = (val_mid - val_left) / (val_right - val_left)
  end function bl_in_tu
!----------------------------------------------------------------------------------------------------
  function intpol_pf_from_NR_grids(which_bound, a_val, n_val, a_arr, n_arr) ! for details see paragraph "Bilinear interpolation" in Numerical Recipes for F77, page 117, eqn. 3.6.4
   implicit none
    integer(kind=2), dimension(1:2) :: loc1, loc2, loc_no_interpol ! loc1, loc2 - indexes that points where alpha_tab_lo and up nad n_tab_lo and up are closest in value to a_val and n_val - indexes point to 
    real(kind=8) :: a_val, n_val                                   ! ratios arrays (p,f: lo and up), for which solutions have been obtained. loc_no_interpol - in case when interpolation is not possible, 
    real(kind=8), dimension(:) :: a_arr, n_arr                     ! indexes with best match and having solutions are chosen.
    real(kind=8), dimension(2) :: intpol_pf_from_NR_grids
    character(len=2) :: which_bound ! "lo" or "up"
    logical :: exit_code
        exit_code = .false.
        if (which_bound == "lo") then
            call find_both_indexes_lo(loc1, loc2, a_val, n_val, loc_no_interpol, exit_code)
        else
            call find_both_indexes_up(loc1, loc2, a_val, n_val, loc_no_interpol, exit_code)
        endif

        if (which_bound == "lo") then
            if (exit_code .eqv. .true. ) then ! interpolation won't work in this case, choosing closest values that have solutions.
                  intpol_pf_from_NR_grids(1) = p_ratios_lo(loc_no_interpol(1),loc_no_interpol(2))
                  intpol_pf_from_NR_grids(2) = f_ratios_lo(loc_no_interpol(1),loc_no_interpol(2))
                return
            endif
            intpol_pf_from_NR_grids(1) = bl_interpol(p_ratios_lo(loc1(1),loc1(2)),p_ratios_lo(loc1(1),loc2(2)), &
                p_ratios_lo(loc2(1),loc1(2)),p_ratios_lo(loc2(1),loc2(2)), bl_in_tu(a_arr(loc1(1)), a_val, a_arr(loc2(1))), &
                bl_in_tu(n_arr(loc1(2)), n_val, n_arr(loc2(2))) )
            intpol_pf_from_NR_grids(2) = bl_interpol(f_ratios_lo(loc1(1),loc1(2)),f_ratios_lo(loc1(1),loc2(2)), &
                f_ratios_lo(loc2(1),loc1(2)),f_ratios_lo(loc2(1),loc2(2)), bl_in_tu(a_arr(loc1(1)), a_val, a_arr(loc2(1))), &
                bl_in_tu(n_arr(loc1(2)), n_val, n_arr(loc2(2))) )
        else if (which_bound == "up") then
            if (exit_code .eqv. .true. ) then ! interpolation won't work in this case, choosing closest values that have solutions.
                intpol_pf_from_NR_grids(1) = p_ratios_up(loc_no_interpol(1),loc_no_interpol(2))
                intpol_pf_from_NR_grids(2) = f_ratios_up(loc_no_interpol(1),loc_no_interpol(2))
                return
            endif
            intpol_pf_from_NR_grids(1) = bl_interpol(p_ratios_up(loc1(1),loc1(2)), p_ratios_up(loc1(1),loc2(2)), &
                p_ratios_up(loc2(1),loc1(2)), p_ratios_up(loc2(1),loc2(2)),bl_in_tu(a_arr(loc1(1)), a_val, a_arr(loc2(1))), &
                bl_in_tu(n_arr(loc1(2)), n_val, n_arr(loc2(2))))
            intpol_pf_from_NR_grids(2) = bl_interpol(f_ratios_up(loc1(1),loc1(2)), f_ratios_up(loc1(1),loc2(2)),& 
                f_ratios_up(loc2(1),loc1(2)), f_ratios_up(loc2(1),loc2(2)),bl_in_tu(a_arr(loc1(1)), a_val, a_arr(loc2(1))), &
                bl_in_tu(n_arr(loc1(2)), n_val, n_arr(loc2(2))))
        else
            print *, "ERROR - provided wrong value of which_bound variable @intpol_pf_from_NR_grids"
            stop
        endif
  end function intpol_pf_from_NR_grids
!----------------------------------------------------------------------------------------------------
  subroutine save_NR_guess_grid(NR_guess_grid, var_name) 
  use initcrspectrum, only: e_small, q_big, max_p_ratio
  use cresp_variables, only: clight ! use units, only: clight
  implicit none
    real(kind=8),dimension(:,:) :: NR_guess_grid
    integer(kind=4) :: width, j
    character(len=11) :: var_name
    character(len=4)  :: extension
    character(len=15) :: f_name
        width     = size(NR_guess_grid, dim=1)
        extension =  ".dat"
        f_name = var_name // extension
        open(31, file=f_name, status="unknown", position="rewind")
            write(31,"(A57,A3,A95)") "This is a storage file for NR init grid, boundary case: ", var_name(9:), &
                    "Saved below: e_small, size(NR_guess_grid,dim=1), size(NR_guess_grid,dim=2), max_p_ratio, q_big,  clight. &
                    & Do not remove content from this file"
            write(31, "(1E15.8, 2I10,10E22.15)") e_small, size(NR_guess_grid,dim=1), size(NR_guess_grid, dim=2), &
                                                 max_p_ratio, q_big, clight
            write(31, "(A1)") " "                            ! Blank line for 
            do j=1, size(NR_guess_grid,dim=2)
                write(31, "(*(E22.15))") NR_guess_grid(:,j)  ! WARNING - THIS MIGHT NEED EXPLICIT INDICATION OF ELEMENTS COUNT IN LINE IN OLDER COMPILERS
            enddo
        print *, var_name,": saving done."
        close(31)

  end subroutine save_NR_guess_grid
!----------------------------------------------------------------------------------------------------  
  subroutine read_NR_guess_grid(NR_guess_grid, var_name, exit_code) ! must be improved, especially for cases when files do not exist
  use initcrspectrum, only: e_small, q_big, max_p_ratio
  use cresp_variables, only: clight ! use units, only: clight
  implicit none
    real(kind=8),dimension(:,:) :: NR_guess_grid
    real(kind=8)                :: svd_e_sm, svd_max_p_r, svd_q_big, svd_clight
    integer(kind=4) :: width, j, svd_cols, svd_rows
    character(len=11) :: var_name
    character(len=4)  :: extension
    character(len=15) :: f_name
    logical           :: exit_code
        width     = size(NR_guess_grid, dim=1)
        extension =  ".dat"
        f_name = var_name // extension
        open(31, file=f_name, status="unknown", position="rewind")
        
            read(31,*) ! Skipping comment line
            read(31,"(1E15.8,2I10,10E22.15)") svd_e_sm, svd_rows, svd_cols, svd_max_p_r, svd_q_big, svd_clight
            if ( abs(svd_e_sm-e_small)/e_small .le. 1.0e-6 .and. svd_rows .eq. size(NR_guess_grid,dim=1) &
                        .and. svd_cols .eq. size(NR_guess_grid,dim=2) .and. abs(max_p_ratio-svd_max_p_r)/max_p_ratio .le. 1.0e-6 &
                        .and. abs(q_big-svd_q_big)/q_big .le. 1.0e-6 .and. abs(clight-svd_clight)/clight .le. 1.0e-6 ) then ! Yeah, it's a long condition, but if met it saves a lot of time.
                read(31, *)
                do j=1, size(NR_guess_grid,dim=2)
                    read(31, "(*(E22.15))") NR_guess_grid(:,j)  ! WARNING - THIS MIGHT NEED EXPLICIT INDICATION OF ELEMENTS COUNT IN LINE IN OLDER COMPILERS
                enddo
                print *, var_name, ": reading done."
                exit_code = .false.
            else
                write(*,"(A70,A4,A6)") "Different initial parameters, proceeding to resolve ratio tables for",var_name(10:),"case."
                write(*,"(2I10,10E22.15)") svd_cols, svd_rows, svd_e_sm, svd_max_p_r, svd_q_big, svd_clight
                write(*,"(2I10,10E22.15)") size(NR_guess_grid, dim=1), size(NR_guess_grid,dim=2), e_small, max_p_ratio, q_big, clight
                exit_code = .true.
            endif
        close(31) 
  end subroutine read_NR_guess_grid
!----------------------------------------------------------------------------------------------------
  function logical_2_int(boolean_arg)
   implicit none
    logical :: boolean_arg
    integer(kind=2) :: logical_2_int
        if ( boolean_arg .eqv. .true.) then
            logical_2_int = 1
        else
            logical_2_int = 0
        endif
  end function logical_2_int
  
!====================================================================================================
! Solver test section
!----------------------------------------------------------------------------------------------------
 subroutine nr_test
 implicit none
   real(kind=8), dimension(ndim) :: x
   logical :: exit_code
        print *, "NR algorithm test, provided test functions are:"
        print *, "f1(x(:)) =  x(1)**2+x(2)**2+4.0d0"
        print *, "f2(x(:)) =  x(1) - x(2)"
        print *, "expected solution is x1 = x2 = sqrt(2)"
   
        selected_function => fvec_test
        x(1) = -5.0
        x(2) = 8.0
        print *, "starting with x =", x
        call NR_algorithm(x, exit_code)
        print *, "x**2 = ", x**2
  end subroutine nr_test
!----------------------------------------------------------------------------------------------------
  function fvec_test(x)
  use constants, only: four
  implicit none
    real(kind=8), dimension(ndim) :: x
    real(kind=8), dimension(ndim) :: fvec_test
    ! expected solution is sqrt(2), sqrt(2)
        fvec_test(1) = x(1)**2+x(2)**2 - four
        fvec_test(2) = x(1) - x(2)
  end function fvec_test

end module cresp_NR_method
