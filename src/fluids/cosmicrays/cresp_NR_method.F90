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
!! \brief This module contains Newton-Raphson 1,2-dim alhorithms, functions for use with these and used by CRESP.
!<

module cresp_NR_method
! pulled by COSM_RAY_ELECTRONS

   use constants,      only: LO, HI
   use initcrspectrum, only: eps

   implicit none


   private
   public :: alpha, n_in, NR_algorithm, NR_algorithm_1D, compute_q, intpol_pf_from_NR_grids, selected_function_1D, &
           &    selected_function_2D, selected_value_check_1D, initialize_arrays, e_small_to_f, q_ratios, fvec_lo, &
           &    fvec_up, fvec_test, cresp_initialize_guess_grids, assoc_pointers
   public :: e_in, nr_test, nr_test_1D, p_ip1, n_tab_up, alpha_tab_up, n_tab_lo, alpha_tab_lo, alpha_tab_q, q_control, & ! list for NR driver
           &    p_a, p_n, p_ratios_lo, f_ratios_lo, p_ratios_up, f_ratios_up, q_grid, lin_interpol_1D, alpha_to_q,   &   ! can be commented out for CRESP and PIERNIK
           &    lin_extrapol_1D, lin_interpolation_1D, nearest_solution

   integer, parameter                               :: ndim = 2
   real(kind=8), allocatable, dimension(:)          :: p_space, q_space
   real(kind=8)                                     :: alpha, p_ratio_4_q, n_in, e_in, p_im1, p_ip1
   real(kind=8), allocatable, dimension(:), target  :: alpha_tab_lo, alpha_tab_up, n_tab_lo, n_tab_up, alpha_tab_q, q_grid
   real(kind=8), allocatable, dimension(:,:),target :: p_ratios_lo, f_ratios_lo, p_ratios_up, f_ratios_up
   integer(kind=4)                                  :: helper_arr_dim
   real(kind=8)                                     :: eps_det = eps * 1.0e-15
   real(kind=8), pointer, dimension(:)              :: p_a => null(), p_n => null() ! pointers for alpha_tab_(lo,up) and n_tab_(lo,up) or optional - other 1-dim arrays
   real(kind=8), pointer, dimension(:,:)            :: p_p => null(), p_f => null() ! pointers for p_ratios_(lo,up) and f_ratios_(lo,up)
#ifdef CRESP_VERBOSED
   integer(kind=4)                                  :: current_bound
   integer, parameter                               :: blen = 2
   character(len=blen), dimension(LO:HI), parameter :: bound_name = ['lo', 'up']
#endif /* CRESP_VERBOSED */

   abstract interface
      function function_pointer_1D(z)
         real(kind=8) :: function_pointer_1D
         real(kind=8) :: z
      end function function_pointer_1D
      subroutine value_control_1D(z, exit_code)
         logical      :: exit_code
         real(kind=8) :: z
      end subroutine value_control_1D
      function function_pointer_2D(z)
         real(kind=8),dimension(2) :: function_pointer_2D
         real(kind=8),dimension(2) :: z
      end function function_pointer_2D
   end interface

   procedure (value_control_1D), pointer    :: selected_value_check_1D => null()
   procedure (function_pointer_1D), pointer :: selected_function_1D => null()
   procedure (function_pointer_2D), pointer :: selected_function_2D => null()

!----------------------------------------------------------------------------------------------------

contains

   subroutine NR_algorithm(x,exit_code)

      use initcrspectrum, only: NR_iter_limit, tol_f, tol_x

      implicit none

      real(kind=8), dimension(ndim), intent(inout) :: x
      logical,                       intent(out)   :: exit_code
      real(kind=8), dimension(ndim)                :: fun_vec_value
      real(kind=8), dimension(1:ndim)              :: cor
      real(kind=8), dimension(size(x),size(x))     :: fun_vec_jac, fun_vec_inv_jac
      real(kind=8)                                 :: err_f, err_x, det
      integer(kind=2)                              :: i

      err_f = tol_f
      err_x = tol_x
      exit_code = .true.

      fun_vec_value = selected_function_2D(x)
      if (maxval(abs(fun_vec_value)) < 0.01 * err_f) then ! in case when f converges at initialization
         exit_code=.false.
!         write(*,"(A33,2E22.15)")"Convergence (f) at initialization", x
         return
      endif

         do i = 1, NR_iter_limit
            if (maxval(abs(fun_vec_value)) < err_f ) then    ! For convergence via value of f
               exit_code=.false.
#ifdef CRESP_VERBOSED
               write(*,"(A47,I4,A12)",advance="no") "Convergence via value of fun_vec_value after ",i, " iterations."!, x, fun_vec_value
#endif /* CRESP_VERBOSED */
            return
         endif

         fun_vec_value = selected_function_2D(x)
         fun_vec_jac = jac_fin_diff(x)                    ! function vector already explicitly provided to jac_fin_diff (finite difference method)

         det = determinant_2d_real(fun_vec_jac)           ! WARNING - algorithm is used for ndim = 2. For more dimensions LU or other methods should be implemented.
         if (abs(det) .lt. eps_det) then              ! Countermeasure against determinant = zero
!      write (*,"(A20)") "WARNING: det ~ 0.0"
            exit_code = .true.
            return
         endif
         fun_vec_inv_jac = invert_2d_matrix(fun_vec_jac,det)

         cor(1) = fun_vec_inv_jac(1,1) *fun_vec_value(1) + fun_vec_inv_jac(1,2) * fun_vec_value(2)
         cor(2) = fun_vec_inv_jac(2,1) *fun_vec_value(1) + fun_vec_inv_jac(2,2) * fun_vec_value(2)
         x = x+cor
!          write(*,'(A20, 2E35.25, A5, 2E22.14)') "Obtained values (x): " , x,' | ', sum(abs(cor)), sum(abs(fun_vec_value))! ,maxval(abs(fun_vec_value)), maxval(abs(cor)),
         if (maxval(abs(cor)) < err_x) then                 ! For convergence via value of correction (cor) table.
#ifdef CRESP_VERBOSED
            write(*,"(A47,I4,A12)",advance="no") "Convergence via value of cor array     after ",i," iterations."
#endif /* CRESP_VERBOSED */
            exit_code = .false.
            return
         endif
      enddo
!       write(*,"(A45,I4,A24)") "  ... WARNING: Maximum number of iterations (",NR_iter_limit,") exceeded @global_newt!"
      exit_code = .true.

   end subroutine NR_algorithm

!----------------------------------------------------------------------------------------------------

   subroutine NR_algorithm_1D(x, exit_code)

      use constants,      only: zero
      use func,           only: operator(.equals.)
      use initcrspectrum, only: NR_iter_limit, tol_f_1D, tol_x_1D

      implicit none

      integer      :: i
      real(kind=8) :: x, delta, dfun_1D, fun1D_val
      logical      :: exit_code, func_check

      delta = huge(1.0)
      func_check = .false.

      do i = 1, NR_iter_limit
         fun1D_val = selected_function_1D(x)
         if ((abs(fun1D_val) .le. tol_f_1D) .and. (abs(delta) .le. tol_f_1D)) then ! delta .le. tol_f acceptable as in case of f convergence we must only check if the algorithm hasn't wandered astray
            exit_code = .false.
            return
         endif

         dfun_1D = derivative_1D(x)

         if ( abs(dfun_1D) .equals. zero ) then
            exit_code = .true.
            return
         endif

         delta   = fun1D_val / derivative_1D(x)

         x = x - delta
         if (abs(delta) .lt. tol_x_1D) then
            exit_code = .false.
            return
         endif
         call selected_value_check_1D(x, func_check)  ! necessary in some cases, when maximal value x can take is defined, for other cases dummy_check_1D function defined
         if ( func_check .eqv. .true. ) return
      enddo

      exit_code = .true. ! if all fails

   end subroutine NR_algorithm_1D

!----------------------------------------------------------------------------------------------------

   function derivative_1D(x) ! via finite difference method

      use constants,      only: half
      use initcrspectrum, only: eps

      implicit none

      real(kind=8),intent(in) :: x
      real(kind=8)            :: dx, derivative_1D
      real(kind=8)            :: dx_par = 1.0e-4

      dx = sign(1.0, x) * min(abs(x*dx_par), dx_par)
      dx = sign(1.0, x) * max(abs(dx), eps) ! dx = 0.0 must not be allowed
      derivative_1D = half * (selected_function_1D(x+dx) - selected_function_1D(x-dx))/dx

   end function derivative_1D

!----------------------------------------------------------------------------------------------------

   subroutine cresp_initialize_guess_grids

      use constants,       only: zero, I_FOUR
      use cresp_variables, only: clight ! use units, only: clight
      use initcrspectrum,  only: e_small, q_big, max_p_ratio, arr_dim, arr_dim_q
      use mpisetup,        only: master

      implicit none

      logical        :: first_run = .true. , save_to_log = .false.
      character(8)   :: date
      character(9)   :: time

      call initialize_arrays
      if ((master) .and. (first_run .eqv. .true. )) then
         helper_arr_dim = int(arr_dim/I_FOUR,kind=4)

         if (.not. allocated(p_space)) allocate(p_space(1:helper_arr_dim)) ! these will be deallocated once initialization is over
         if (.not. allocated(q_space)) allocate(q_space(1:helper_arr_dim)) ! these will be deallocated once initialization is over

         call date_and_time(date,time)
         call date_and_time(DATE=date)
         call date_and_time(TIME=time)

         p_ratios_up = zero ; f_ratios_up = zero
         p_ratios_lo = zero ; f_ratios_lo = zero

         q_grid      = q_big; q_grid(int(arr_dim_q/2):) = -q_big

         call fill_guess_grids

         print *, "Are there zeros? (q_ratios)",    count(abs(q_space).le.zero)
         print *, "Are there zeros? (p_ratios_up)", count(p_ratios_up.le.zero)
         print *, "Are there zeros? (f_ratios_up)", count(f_ratios_up.le.zero)
         print *, "Are there zeros? (p_ratios_lo)", count(p_ratios_lo.le.zero)
         print *, "Are there zeros? (f_ratios_lo)", count(f_ratios_lo.le.zero)
         print *, "Count of array elements:", size(p_ratios_lo)
         print *,"----------"
         if (save_to_log) then
            open(15, file="log_NR_solve",position="append")
            write (15,*) "------------------------------------------"
            write (15,"(A,2x,A,2x)") "Run on: ", date, "at: ", time
            write (15,*) "For set of parameters: e_small, size(NR_guess_grid,dim=1), size(NR_guess_grid, dim=2),", &
                                             "max_p_ratio, q_big, clight"
            write (15, "(1E15.8, 2I10,10E22.15)") e_small, size(p_space), size(q_space), &
                                             max_p_ratio, q_big, clight
            write (15,*) "Are there zeros? (q_ratios)",    count(abs(q_space).le.zero)
            write (15,*) "Are there zeros? (p_ratios_up)", count(p_ratios_up.le.zero), &
                              real(count(p_ratios_up.le.zero))/real(size(p_ratios_up)) * 100.0,"%"
            write (15,*) "Are there zeros? (f_ratios_up)", count(f_ratios_up.le.zero), &
                              real(count(f_ratios_up.le.zero))/real(size(f_ratios_up)) * 100.0,"%"
            write (15,*) "Are there zeros? (p_ratios_lo)", count(p_ratios_lo.le.zero), &
                              real(count(p_ratios_lo.le.zero))/real(size(p_ratios_lo)) * 100.0,"%"
            write (15,*) "Are there zeros? (f_ratios_lo)", count(f_ratios_lo.le.zero), &
                              real(count(f_ratios_lo.le.zero))/real(size(f_ratios_lo)) * 100.0,"%"
            write (15,*) "Count of array elements:", size(p_ratios_lo)
            write (15,*) "----------"
            close(15)
         else
            print *, "!!! Warning: save_to_log = F; result will not be registered in LOG file !!!"
         endif
!
         if (allocated(p_space)) deallocate(p_space) ! only needed at initialization
         if (allocated(q_space)) deallocate(q_space)

         first_run = .false.
      endif

      call cresp_NR_mpi_exchange

   end subroutine cresp_initialize_guess_grids

!----------------------------------------------------------------------------------------------------

   subroutine fill_guess_grids

      use constants,      only: zero, one, I_ONE, half
      use initcrspectrum, only: q_big, force_init_NR, NR_run_refine_pf, p_fix_ratio, e_small_approx_init_cond, arr_dim, arr_dim_q, max_p_ratio

      implicit none

      integer(kind=4)   :: i, j, int_logical_p, int_logical_f
      logical           :: exit_code
      real(kind=8) :: a_min_lo=huge(one), a_max_lo=tiny(one), a_min_up=huge(one), a_max_up=tiny(one),&
                    & n_min_lo=huge(one), n_max_lo=tiny(one), n_min_up=huge(one), n_max_up=tiny(one),&
                    & a_min_q=tiny(one),  a_max_q=tiny(one)

      q_space = zero
      do i=1, int(half*helper_arr_dim)
         q_space(i) = ln_eval_array_val(i, q_big, real(0.05,kind=8), int(1,kind=4), int(half*helper_arr_dim,kind=4)) ! BEWARE: magic number
      enddo

      do i= 1, int(half*helper_arr_dim)!, arr_dim
         q_space(int(half*helper_arr_dim,kind=4)+i) = -q_space(int(half*helper_arr_dim,kind=4)+1-i)
      enddo

! setting up a grids of ratios to be used as phase space for NR tabs, obtained later
      do i = 1, helper_arr_dim
         p_space(i) = max_p_ratio**(real(i)/real(helper_arr_dim)) ! ind_to_flog(i,1.000001 , max_p_ratio) ! max_p_ratio**(real(i)/real(arr_dim))
      enddo
      do i = 1, helper_arr_dim
         do j = 1, helper_arr_dim
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

      a_min_lo = 0.2       ! BEWARE: magic numbers!
      a_max_lo = 0.999999
      a_min_up = 1.000005
      a_max_up = 200.0
      n_min_lo = 1.0e-11
      n_max_lo = 5000.0
      n_min_up = 1.0e-12
      n_max_up = 1000.0

      do i=1, arr_dim
         alpha_tab_lo(i) = ind_to_flog(i, a_min_lo, a_max_lo, arr_dim) ! a_min_lo * ten**((log10(a_max_lo/a_min_lo))/real(arr_dim-1,kind=8)*real((i-1),kind=8))
         alpha_tab_up(i) = ind_to_flog(i, a_min_up, a_max_up, arr_dim) ! a_min_up * ten**((log10(a_max_up/a_min_up))/real(arr_dim-1,kind=8)*real((i-1),kind=8))
         n_tab_lo(i)     = ind_to_flog(i, n_min_lo, n_max_lo, arr_dim) ! n_min_lo * ten**((log10(n_max_lo/n_min_lo))/real(arr_dim-1,kind=8)*real((i-1),kind=8))
         n_tab_up(i)     = ind_to_flog(i, n_min_up, n_max_up, arr_dim) ! n_min_up * ten**((log10(n_max_up/n_min_up))/real(arr_dim-1,kind=8)*real((i-1),kind=8))
      enddo

      if (e_small_approx_init_cond .eq. 1) then
         write (*, "(A36)", advance="no") "Reading (up) boundary ratio files..."
         do j = 1,2
            call read_NR_guess_grid(p_ratios_up, "p_ratios_up", exit_code) ;  int_logical_p = logical_2_int(exit_code)
            call read_NR_guess_grid(f_ratios_up, "f_ratios_up", exit_code) ;  int_logical_f = logical_2_int(exit_code)

            if ( int_logical_f + int_logical_p .gt. 0 .or. force_init_NR .eqv. .true.) then
   ! Setting up the "guess grid" for p_up case
               call fill_boundary_grid(HI, p_ratios_up, f_ratios_up)
            else
               print *," >> Will not solve ratios table (up), reading data from file instead."
            endif

            if ( NR_run_refine_pf .eqv. .true.) then
               call assoc_pointers(HI)
               call refine_all_directions(HI)
            endif

            call save_NR_guess_grid(p_ratios_up,"p_ratios_up")
            call save_NR_guess_grid(f_ratios_up,"f_ratios_up")
         enddo

         write (*, "(A36)", advance="no") "Reading (lo) boundary ratio files"
         do j = 1,2
            call read_NR_guess_grid(p_ratios_lo, "p_ratios_lo", exit_code) ;   int_logical_p = logical_2_int(exit_code)
            call read_NR_guess_grid(f_ratios_lo, "f_ratios_lo", exit_code) ;   int_logical_f = logical_2_int(exit_code)

            if ( int_logical_f + int_logical_p .gt. 0 .or. force_init_NR .eqv. .true.) then
   ! Setting up the "guess grid" for p_lo case
               call fill_boundary_grid(LO, p_ratios_lo, f_ratios_lo)
            else
               print *," >> Will not solve ratios table (lo), reading data from file instead."
            endif

            if (NR_run_refine_pf .eqv. .true.) then
               call assoc_pointers(LO)
               call refine_all_directions(LO)
            endif

            call save_NR_guess_grid(p_ratios_lo,"p_ratios_lo")
            call save_NR_guess_grid(f_ratios_lo,"f_ratios_lo")
         enddo
      endif

      a_min_q = one  + epsilon(one)
      a_max_q = (one + epsilon(one)) * p_fix_ratio
      j = min(arr_dim_q - int(arr_dim_q/100 ,kind=4), arr_dim_q - I_ONE)               ! BEWARE: magic number

      do while (q_grid(j) .le. (-q_big) .and. (q_grid(arr_dim_q) .le. (-q_big)) )
         a_max_q = a_max_q - a_max_q*0.005                                             ! BEWARE: magic number
         do i = 1, arr_dim_q
            alpha_tab_q(i)  = ind_to_flog(i, a_min_q, a_max_q, arr_dim_q)
         enddo
         call fill_q_grid(i_incr=1) ! computing q_grid takes so little time, that saving the grid is not necessary.
      enddo

#ifdef CRESP_VERBOSED
      do i = 1, arr_dim_q
         print "(A1,I3,A7,2F18.12)", "[ ", i,"] a : q ", q_grid(i), alpha_tab_q(i)
      enddo

      print *,"alpha_tab_lo(i),      alpha_tab_up(i),        n_tab_lo(i),        n_tab_up(i)  |       p_space(i),     q_space(i)"
      do i = 1, arr_dim
         if (i .le. helper_arr_dim) then
            print *,i,"|",  alpha_tab_lo(i), alpha_tab_up(i), n_tab_lo(i), n_tab_up(i), alpha_tab_q(i), "| i = ", &
                          min(i,helper_arr_dim), p_space(min(i,helper_arr_dim)), q_space(min(i,helper_arr_dim)), &
                          p_space(min(i,helper_arr_dim))**(-q_space(min(i,helper_arr_dim)))
         else
            print *,i,"|",  alpha_tab_lo(i), alpha_tab_up(i), n_tab_lo(i), n_tab_up(i), alpha_tab_q(i)
        endif
      enddo
      print *, "-----------"
#endif /* CRESP_VERBOSED */

   end subroutine fill_guess_grids

!----------------------------------------------------------------------------------------------------
   subroutine refine_all_directions(bound_case)

      implicit none

      integer(kind=4), intent(in) :: bound_case

      print *,"Running refine for:", bound_name(bound_case), " boundary"

      call refine_ij(p_p, p_f,  1, -1)
      call refine_ji(p_p, p_f,  1, -1)
      call refine_ij(p_p, p_f, -1, -1)
      call refine_ji(p_p, p_f, -1, -1)
      call refine_ij(p_p, p_f,  1,  1)
      call refine_ji(p_p, p_f,  1,  1)
      call refine_ij(p_p, p_f, -1,  1)
      call refine_ji(p_p, p_f, -1,  1)

   end subroutine refine_all_directions

!----------------------------------------------------------------------------------------------------
   function ln_eval_array_val(i, arr_min, arr_max, min_i, max_i)

      implicit none

      real(kind=8)    :: b, arr_min, arr_max, ln_eval_array_val
      integer(kind=4) :: i, max_i, min_i

      b = (log(real(max_i)) -log(real(min_i)))/ (arr_max - arr_min)
      ln_eval_array_val = (arr_min-log(real(min_i))/b ) + log(real(i)) / b

   end function ln_eval_array_val

!----------------------------------------------------------------------------------------------------
   subroutine assoc_pointers(bound_case)

      implicit none

      integer(kind=4), intent(in) :: bound_case

      if (bound_case == LO) then
         p_a => alpha_tab_lo
         p_n => n_tab_lo
         p_p => p_ratios_lo
         p_f => f_ratios_lo
         selected_function_2D => fvec_lo
      endif
      if (bound_case == HI) then
         p_a => alpha_tab_up
         p_n => n_tab_up
         p_p => p_ratios_up
         p_f => f_ratios_up
         selected_function_2D => fvec_up
      endif

#ifdef CRESP_VERBOSED
      current_bound = bound_case
#endif /* CRESP_VERBOSED */

   end subroutine assoc_pointers

!----------------------------------------------------------------------------------------------------
   subroutine fill_boundary_grid(bound_case, fill_p, fill_f) ! to be paralelized

      use constants,      only: zero
      use initcrspectrum, only: arr_dim

      implicit none

      integer(kind=4), intent(in)  :: bound_case ! HI or LO
      real(kind=8), dimension(1:2) :: x_vec, prev_solution, prev_solution_1, x_step
#ifdef CRESP_VERBOSED
      real(kind=8), dimension(1:2) :: x_in
#endif /* CRESP_VERBOSED */
      real(kind=8), dimension(:,:) :: fill_p, fill_f
      integer(kind=4) :: i, j, is, js
      logical         :: exit_code, new_line
      character(len=6) :: nam = "Solve "

      prev_solution(1) = p_space(1)
      prev_solution(2) = p_space(1)**q_space(1)
      prev_solution_1 = prev_solution

      call assoc_pointers(bound_case)

      call sleep(1)

      fill_p = zero ; fill_f = zero
      x_step = zero

      do i =1, arr_dim
         new_line = .true.
         prev_solution = prev_solution_1 ! easier to find when not searching from the top
         do j = 1, arr_dim ! j_incr = 1
            exit_code = .true.
            alpha = p_a(i)
            n_in  = p_n(j)
#ifdef CRESP_VERBOSED
            write(*,"(A14,A2,A2,2I4,A9,I4,A1)",advance="no") "Now solving (",bound_name(bound_case),") ",i,j,", sized ",arr_dim," "
#endif /* CRESP_VERBOSED */

            call seek_solution_prev(fill_p(i,j), fill_f(i,j), prev_solution, nam, exit_code)

            if (exit_code .eqv. .false. .and. new_line .eqv. .true.) then
               prev_solution_1 = prev_solution
               new_line = .false.
            endif

            if ( exit_code .eqv. .true. ) then
               if (j-2 .ge. 1 .and. j-2 .le. arr_dim) then
                  call step_extr(fill_p(i,j-2:j),fill_f(i,j-2:j),p_n(j-2:j),nam,exit_code)
               endif
               if ((j-1) .ge. 1 ) then
                  if (fill_p(i,j-1).gt.zero) call seek_solution_step(fill_p(i,j),fill_f(i,j),prev_solution,i,j-1,nam,exit_code)
               endif
            endif
            if (exit_code .eqv. .true.) then !still...
               do is =1, helper_arr_dim
                  do js = 1, helper_arr_dim
                     x_vec(1) = p_space(is)
                     x_vec(2) = p_space(is)**(-q_space(js))
#ifdef CRESP_VERBOSED
                     x_in = x_vec
#endif /* CRESP_VERBOSED */
                     if (exit_code .eqv. .true.) then
                        call NR_algorithm(x_vec, exit_code)
                        if ( exit_code .eqv. .false. ) then
                           fill_p(i,j) = x_vec(1) ! i index - alpha, j index - n_in
                           fill_f(i,j) = x_vec(2)
                           prev_solution = x_vec
#ifdef CRESP_VERBOSED
                           call msg_success("    ",nam,x_in, x_vec)
#endif /* CRESP_VERBOSED */
                           exit
                        endif
                     endif
                  enddo
               enddo
            else
               if (prev_solution(1) .le. eps ) then
                  prev_solution(1) = prev_solution_1(1)
               else if (prev_solution(2) .le. eps) then
                  prev_solution(2) = prev_solution_1(2)
               else
                  prev_solution(1) = fill_p(i,j)
                  prev_solution(2) = fill_f(i,j)
               endif
            endif
#ifdef CRESP_VERBOSED
            if (exit_code .eqv. .true.) print *,""
#endif /* CRESP_VERBOSED */
         enddo
      enddo

      fill_p = abs(fill_p)
      fill_f = abs(fill_f)

#ifdef CRESP_VERBOSED
      print *,""
#endif /* CRESP_VERBOSED */

   end subroutine fill_boundary_grid
!----------------------------------------------------------------------------------------------------
   subroutine step_seek(x_step, prev_sol, ii, jj, i_sol, j_sol, exit_code, nstep)

      use constants, only: I_ONE

      implicit none

      real(kind=8), dimension(1:2), intent(in)    :: prev_sol
      real(kind=8), dimension(1:2), intent(out)   :: x_step
      integer(kind=4),              intent(in)    :: ii, jj, i_sol, j_sol, nstep
      logical,                      intent(inout) :: exit_code
      integer(kind=4)                             :: i, j

      do i = -1,1
         do j = -1,1
            if (exit_code .eqv. .true.) then
               x_step(1) = prev_sol(1) + (p_space(max(min(i_sol+i,helper_arr_dim),I_ONE)) - prev_sol(1)) / real(nstep - ii + 0.1)
               x_step(2) = prev_sol(2) + (p_space(max(min(i_sol+i,helper_arr_dim),I_ONE))**(-q_space(max(min(j_sol+j,  &
                        &  helper_arr_dim),1))) - prev_sol(2)) / real(nstep - jj + 0.1)
               call NR_algorithm(x_step, exit_code)
               if (exit_code .eqv. .false.) return
            endif
         enddo
      enddo

   end subroutine step_seek

!----------------------------------------------------------------------------------------------------
   subroutine refine_ji(ref_p, ref_f, i_incr, j_incr) ! ref_f and ref_p should already be partially filled with solutions

      use constants,      only: zero
      use initcrspectrum, only: arr_dim

      implicit none

      integer(kind=4),              intent(in)    :: i_incr, j_incr
      real(kind=8), dimension(:,:), intent(inout) :: ref_p, ref_f
      integer(kind=4)                             :: i, j, i_beg, i_end, j_beg, j_end
      real(kind=8), dimension(1:2)                :: prev_solution
      character(len=6)                            :: nam = "Refine"
      logical                                     :: exit_code, new_line

      if (allocated(p_space) .and. allocated(q_space)) then
         prev_solution(1) = p_space(1)              ! refine must be called before these are deallocated
         prev_solution(2) = p_space(1)**q_space(1)
         call prepare_indices(i_incr, i_beg, i_end)
         call prepare_indices(j_incr, j_beg, j_end)
         do j = j_beg, j_end, j_incr
            new_line = .true.
            do i = i_beg, i_end, i_incr
               alpha = p_a(i)
               n_in  = p_n(j)
               exit_code = .true.
               if (ref_p(i,j) .gt. zero .and. ref_f(i,j) .gt. zero) then
                  prev_solution(1) = ref_p(i,j)
                  prev_solution(2) = ref_f(i,j)
                  new_line = .false.
               else
                  call seek_solution_prev(ref_p(i,j), ref_f(i,j), prev_solution, nam, exit_code) ! works for most cases
                  if ( exit_code .eqv. .true. ) then
                     if (i-2*i_incr .ge. 1 .and. i-2*i_incr .le. arr_dim) then
                        call step_extr(ref_p(i-2*i_incr:i:i_incr,j),ref_f(i-2*i_incr:i:i_incr,j),&
                                    &  p_a(i-2*i_incr:i:i_incr),nam,exit_code)
                     endif
                     if (i-i_incr .ge. 1 .and. i+i_incr .ge. 1 .and. i-i_incr .le. arr_dim .and. i+i_incr .le. arr_dim) then
                        call step_inpl(ref_p(i-i_incr:i+i_incr:i_incr,j),ref_f(i-i_incr:i+i_incr:i_incr,j),&
                                    &  i_incr,p_a(i-i_incr:i+i_incr:i_incr),nam,exit_code)
                     endif
                  endif
               endif
            enddo
         enddo
      else
         print *, "@cresp_NR_method: refine_grids called after array deallocation, stopping"
         stop
      endif

   end subroutine refine_ji

!----------------------------------------------------------------------------------------------------
   subroutine refine_ij(ref_p, ref_f, i_incr, j_incr) ! ref_f and ref_p should already be partially filled with solutions

      use constants,      only: zero
      use initcrspectrum, only: arr_dim

      implicit none

      integer(kind=4),              intent(in)    :: i_incr, j_incr
      real(kind=8), dimension(:,:), intent(inout) :: ref_p, ref_f
      integer(kind=4)                             :: i, j, i_beg, i_end, j_beg, j_end
      real(kind=8), dimension(1:2)                :: prev_solution
      character(len=6)                            :: nam = "Refine"
      logical                                     :: exit_code, new_line, i_primary

      if (allocated(p_space) .and. allocated(q_space)) then
         prev_solution(1) = p_space(1)              ! refine must be called before these are deallocated
         prev_solution(2) = p_space(1)**q_space(1)
         i_primary = .true.
         call prepare_indices(i_incr, i_beg, i_end)
         call prepare_indices(j_incr, j_beg, j_end)
         do i = i_beg, i_end, i_incr
            new_line = .true.
            do j = j_beg, j_end, j_incr
               alpha = p_a(i)
               n_in  = p_n(j)
               exit_code = .true.
               if (ref_p(i,j) .gt. zero .and. ref_f(i,j) .gt. zero) then
                  prev_solution(1) = ref_p(i,j)
                  prev_solution(2) = ref_f(i,j)
                  new_line = .false.
               else
                  call seek_solution_prev(ref_p(i,j), ref_f(i,j), prev_solution, nam, exit_code) ! works for the most cases
                  if ( exit_code .eqv. .true. ) then
                     if (j-2*j_incr .ge. 1 .and. j-2*j_incr .le. arr_dim) then
                        call step_extr(ref_p(i,j-2*j_incr:j:j_incr), &
                                    &  ref_f(i,j-2*j_incr:j:j_incr),p_n(j-2*j_incr:j),nam,exit_code)
                     endif
                     if (j-j_incr .ge. 1 .and. j+j_incr .ge. 1 .and. j-j_incr .le. arr_dim &
                                    &  .and. j+j_incr .le. arr_dim) then
                        call step_inpl(ref_p(i,j-j_incr:j+j_incr:j_incr),ref_f(i,j-j_incr:j+j_incr:j_incr), &
                                          &  j_incr,p_n(j-j_incr:j+j_incr:j_incr),nam,exit_code)
                     endif
                  endif
               endif
            enddo
         enddo
      else
         print *, "@cresp_NR_method: refine_grids called after array deallocation, stopping"
         stop
      endif
   end subroutine refine_ij
 !----------------------------------------------------------------------------------------------------
   subroutine step_extr(p3, f3, arg, sought_by, exit_code) ! checked

      use constants, only: zero

      implicit none

      real(kind=8), dimension(1:3),intent(inout) :: f3, p3
      real(kind=8), dimension(1:3),intent(in)    :: arg
      character(len=6),            intent(inout) :: sought_by
      logical,                     intent(inout) :: exit_code
      real(kind=8), dimension(1:2)               :: x_vec_0, x_vec, delta, x_in
      integer(kind=4)                            :: nstep = 100, k
!         alpha and n are set !

      if ( minval(p3(1:2)) .gt. tiny(zero) .and. p3(3) .le. zero ) then ! sometimes NaNs and numbers of order e-317 appear; must be looked into
         x_vec_0 = (/ p3(2) , f3(2) /)
         delta(1) = lin_extrapol_1D(p3(1:2), arg(1:2), arg(3)) - p3(2) ! direction is not relevant in this case
         delta(2) = lin_extrapol_1D(f3(1:2), arg(1:2), arg(3)) - f3(2)
         delta = delta/nstep
         do k = 0, nstep
            x_vec = x_vec_0 + delta * k ! first iteration is a simple extrapolation
            x_in = x_vec
            call NR_algorithm(x_vec, exit_code)
            if (exit_code .eqv. .false.) then
               x_vec = abs(x_vec)
#ifdef CRESP_VERBOSED
               call msg_success("extr", sought_by,x_in,x_vec)
#endif /* CRESP_VERBOSED */
               p3(3) = x_vec(1)
               f3(3) = x_vec(2)
               return
            endif
         enddo
      endif

      return
      if (.false.) k = len(sought_by) ! suppress compiler warnings

   end subroutine step_extr

!----------------------------------------------------------------------------------------------------
   subroutine step_inpl(p3, f3, incr, args, sought_by, exit_code) ! checked

   use constants, only: zero

      implicit none

      real(kind=8), dimension(1:3), intent(inout) :: p3, f3
      integer(kind=4),              intent(in)    :: incr
      real(kind=8), dimension(1:3), intent(in)    :: args
      character(len=6),             intent(in)    :: sought_by
      logical,                      intent(inout) :: exit_code
      real(kind=8), dimension(1:2)                :: x_vec, x_vec_0, delta, x_in
      integer(kind=4)                             :: k, nstep = 100
!         alpha and n are set !

      if ( exit_code .eqv. .true.) then
         if ( min(p3(1),p3(3)) .gt. tiny(zero) .and. p3(2) .le. zero ) then ! sometimes NaNs and numbers of order e-317 appear; must be looked into
            x_vec_0 = (/ p3(1) , f3(1) /)
            delta(1) = lin_interpolation_1D( (/ p3(2-incr), p3(2+incr) /), (/ args(2-incr), args(2+incr) /), args(2) ) - p3(1)
            delta(2) = lin_interpolation_1D( (/ f3(2-incr), f3(2+incr) /), (/ args(2-incr), args(2+incr) /), args(2) ) - f3(1)
            x_in = x_vec_0 + delta ! gives the interpolated value as starting one
            delta = delta/nstep
            do k = 0, nstep
               x_vec = x_vec_0 + delta * k
               x_in = x_vec
               call NR_algorithm(x_vec, exit_code)
               if (exit_code .eqv. .false.) then ! first iteration is a simple extrapolation
                  x_vec = abs(x_vec)
#ifdef CRESP_VERBOSED
                  call msg_success("inpl", sought_by, x_in, x_vec)
#endif /* CRESP_VERBOSED */
                  p3(2) = x_vec(1)
                  f3(2) = x_vec(2)
                  return
               endif
               x_vec = x_vec + delta
            enddo
         endif
      endif

      return
      if (.false.) k = len(sought_by) ! suppress compiler warnings

   end subroutine step_inpl
!----------------------------------------------------------------------------------------------------
#ifdef CRESP_VERBOSED
   subroutine msg_success(met_name, sought_by, x_in, x_out)

      implicit none

      real(kind=8), dimension(1:), intent(in) :: x_in
      real(kind=8), dimension(1:), intent(in) :: x_out
      character(len=4),            intent(in) :: met_name
      character(len=6),            intent(in) :: sought_by

      write (*, "(A6,A13,2E16.9)",advance="no") sought_by," (alpha, n): ",alpha,n_in
      write (*, "(A5,A4,A42, 2E19.10e3)",advance="no") " -> (",met_name,") solution obtained, (p_ratio, f_ratio) = ", x_out
      write (*, "(A21, 2E17.10)",advance="no") ", provided input:", x_in ; print *,""

   end subroutine msg_success
#endif /* CRESP_VERBOSED */
!----------------------------------------------------------------------------------------------------
   function lin_interpolation_1D(fun, arg, arg_mid)

      use constants, only: one

      implicit none

      real(kind=8) :: weight, arg_mid, lin_interpolation_1D
      real(kind=8), dimension(1:2), intent(in):: arg, fun

      weight   = (arg_mid - arg(1)) / (arg(2) - arg(1))
      lin_interpolation_1D =  fun(1) * (one - weight) + fun(2) * (one - weight)

   end function lin_interpolation_1D
!----------------------------------------------------------------------------------------------------
   subroutine seek_solution_prev(p2ref, f2ref, prev_solution, sought_by, exit_code)

      implicit none

      real(kind=8),                 intent(inout) :: p2ref, f2ref
      real(kind=8), dimension(1:2), intent(inout) :: prev_solution
      character(len=6),             intent(in)    :: sought_by
      logical,                      intent(inout) :: exit_code
      real(kind=8), dimension(1:2)                :: x_vec

      if (exit_code .eqv. .true.) then
         x_vec = prev_solution
         call NR_algorithm(x_vec, exit_code)
         if (exit_code .eqv. .false.) then
            x_vec = abs(x_vec)
            p2ref = x_vec(1)
            f2ref = x_vec(2)
#ifdef CRESP_VERBOSED
            call msg_success("prev", sought_by, prev_solution, x_vec)
#endif /* CRESP_VERBOSED */
            prev_solution = x_vec
            return
         endif
      endif

      return
      if (.false.) x_vec(1) = float(len(sought_by)) ! suppress compiler warnings

   end subroutine seek_solution_prev
!----------------------------------------------------------------------------------------------------
   subroutine seek_solution_step(p2ref, f2ref, prev_solution, i_obt, j_obt, sought_by, exit_code)

      implicit none

      real(kind=8),                 intent(out)   :: p2ref, f2ref
      real(kind=8), dimension(1:2), intent(inout) :: prev_solution
      integer(kind=4),              intent(in)    :: i_obt, j_obt
      character(len=6),             intent(in)    :: sought_by
      logical,                      intent(inout) :: exit_code
      real(kind=8), dimension(1:2)                :: x_vec, x_step
      integer(kind=4)                             :: ii, jj, nstep = 3
!    alpha and n are set !

      if (exit_code .eqv. .true. ) then
         do ii = 0, nstep
            do jj = 0,nstep
               call step_seek(x_step, prev_solution, ii, jj, i_obt, j_obt, exit_code, nstep)
               x_vec = x_step
               if (exit_code .eqv. .false. ) then
                  x_step = abs(x_step)
                  p2ref = x_step(1)
                  f2ref = x_step(2)
                  prev_solution = x_step
#ifdef CRESP_VERBOSED
                  call msg_success("step", sought_by, x_step, x_vec)
#endif /* CRESP_VERBOSED */
                  return
               endif
            enddo
         enddo
      endif

      return
      if (.false.) ii = len(sought_by) ! suppress compiler warnings

   end subroutine seek_solution_step
!----------------------------------------------------------------------------------------------------
   subroutine fill_q_grid(i_incr)

      use initcrspectrum, only: p_fix_ratio, arr_dim_q

      implicit none

      integer(kind=4) :: i,j, i_beg, i_end, i_incr
      real(kind=8)    :: x, prev_solution
      logical         :: exit_code

      selected_function_1D => alpha_to_q
      selected_value_check_1D  => q_control

      i_beg = 1
      i_end = arr_dim_q

      p_ratio_4_q = p_fix_ratio

      prev_solution = q_space(int(helper_arr_dim/2))

      do i = i_beg, i_end, i_incr
         exit_code = .true.
         alpha = alpha_tab_q(i)
#ifdef CRESP_VERBOSED
         write(*,"(A25,1I4,A9,I4,A10,1E16.9)",advance="no") "Now solving (q_grid) no.",i,", sized ",arr_dim_q ,", (alpha): ",alpha
#endif /* CRESP_VERBOSED */
         x = prev_solution
         call NR_algorithm_1D(x, exit_code)
         if ( exit_code .eqv. .true. ) then
            do j = 1, helper_arr_dim
               if (exit_code .eqv. .true. ) then
                  x = q_space(j)
                  call NR_algorithm_1D(x,exit_code)
                  if ( exit_code .eqv. .false.) then
                     q_grid(i) = x
                     prev_solution = x
#ifdef CRESP_VERBOSED
                     write (*, "(A44, 2E22.15)",advance="no") " ->        solution obtained, q_grid = ", x
#endif /* CRESP_VERBOSED */
                  endif
               endif
         enddo
         else
#ifdef CRESP_VERBOSED
            write (*, "(A44, 1E22.15)",advance="no") " -> (prev) solution obtained, q_grid = ", x
#endif /* CRESP_VERBOSED */
            q_grid(i) = x
         endif
#ifdef CRESP_VERBOSED
         print *,""
#endif /* CRESP_VERBOSED */
      enddo

   end subroutine fill_q_grid
!----------------------------------------------------------------------------------------------------
   subroutine q_control(x, exit_code)

      use constants,      only: one
      use initcrspectrum, only: q_big

      implicit none

      real(kind=8) :: x
      logical      :: exit_code

      if ( abs(x) .ge. q_big ) then
         x = sign(one, x) * q_big
         exit_code = .true.
         return
      endif

   end subroutine q_control
 !----------------------------------------------------------------------------------------------------
   function alpha_to_q(x) ! this one (as of now) is only usable with fixed p_ratio_4_q bins (middle ones)

      use constants,      only: one, three, four
      use initcrspectrum, only: eps

      implicit none

      real(kind=8)  :: x, alpha_to_q

      if (abs(x - three) .lt. eps) then
         alpha_to_q = -alpha + (-one + p_ratio_4_q)/log(p_ratio_4_q)
      else if (abs(x - four) .lt. eps) then
         alpha_to_q = -alpha + p_ratio_4_q*log(p_ratio_4_q)/(p_ratio_4_q - one)
      else
         alpha_to_q = -alpha + ((three-x)/(four-x))*((p_ratio_4_q**(four-x)-one)/(p_ratio_4_q**(three-x)-one))
      endif

   end function alpha_to_q
!----------------------------------------------------------------------------------------------------
   subroutine prepare_indices(ind_incr, ind_beg, ind_end)

      use initcrspectrum, only: arr_dim

      implicit none

      integer(kind=4), intent(out) :: ind_beg, ind_end
      integer(kind=4), intent(in)  :: ind_incr

      if (ind_incr .eq. 1 ) then
         ind_beg = 1 ; ind_end = arr_dim
      else if (ind_incr .eq. -1) then
         ind_beg = arr_dim ; ind_end = 1
      endif

   end subroutine prepare_indices
 !----------------------------------------------------------------------------------------------------
   function q_ratios(f_ratio, p_ratio)

      use constants, only: zero

      implicit none

      real(kind=8), intent(in) :: f_ratio, p_ratio
      real(kind=8)             :: q_ratios

      q_ratios = zero
      q_ratios = -log10(f_ratio) / log10(p_ratio)

   end function q_ratios
!---------------------------------------------------------------------------------------------------
! Function estimating values of jacobian via finite difference method
!---------------------------------------------------------------------------------------------------
   function jac_fin_diff(x) ! jacobian via finite difference method

      use constants, only: half

      implicit none

      real(kind=8), dimension(ndim)            :: x, xp, xm
      real(kind=8), dimension(size(x),size(x)) :: jac_fin_diff
      real(kind=8), dimension(size(x))         :: dx
      real(kind=8), parameter                  :: dx_par = 1.0e-3, dx_min = epsilon(dx_par)
      integer(kind=2)                          :: j

      do j = 1,size(x)
         dx(j) = max(x(j), dx_min )          ! assure dx > zero
         dx(j) = min(dx(j)*dx_par, dx_par)   ! the value of dx is scaled not to go over value of x
         if (x(j) .eq. dx(j)) dx(j) = half * dx(j)
         xp = x ; xm = x
         xp(j) = x(j) - dx(j) ;  xm(j) = x(j) + dx(j)
         jac_fin_diff(:,j)  = half*( selected_function_2D(xp) - selected_function_2D(xm)) / dx(j)
      enddo

   end function jac_fin_diff
!----------------------------------------------------------------------------------------------------
  function determinant_2d_real(matrix_2d_real)

      implicit none

      real(kind=8), dimension(2,2), intent(in) :: matrix_2d_real
      real(kind=8)                             :: determinant_2d_real

      determinant_2d_real = matrix_2d_real(1,1) * matrix_2d_real (2,2) - ( matrix_2d_real(2,1) * matrix_2d_real(1,2) )

   end function determinant_2d_real
!----------------------------------------------------------------------------------------------------
   function get_cofactor_matrix_2d_real(matrix_2d_real)

      implicit none

      integer(kind=1),parameter            :: m_dim = 2
      real(kind=8), dimension(m_dim,m_dim) :: matrix_2d_real
      real(kind=8), dimension(m_dim,m_dim) :: get_cofactor_matrix_2d_real
      integer(kind=1)                      :: i,j

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
                            &  (p_ratio **(three - q_in)- one))
      endif

   end function encp_func_2_zero_up

!---------------------------------------------------------------------------------------------------
   function encp_func_2_zero_lo(p_ratio, alpha_cnst, q_in) ! from eqn. 29

      use constants,     only: one, three, four

      implicit none

      real(kind=8) :: p_ratio, q_in, alpha_cnst
      real(kind=8) :: encp_func_2_zero_lo

      if (abs(q_in - three) .lt. eps) then
         encp_func_2_zero_lo = -alpha_cnst + ((one - p_ratio)/log(p_ratio)) / p_ratio
      else if (abs(q_in - four) .lt. eps) then
         encp_func_2_zero_lo = -alpha_cnst + log(p_ratio)/(p_ratio - one)
      else
         encp_func_2_zero_lo = -alpha_cnst + ((three - q_in)/(four - q_in))*((p_ratio **(four - q_in) - one)/ &
                           &   (p_ratio **(three - q_in) - one)) / p_ratio
      endif

   end function encp_func_2_zero_lo

!----------------------------------------------------------------------------------------------------
   function n_func_2_zero_up(p_ratio, f_ratio, n_cnst, q_in) ! from eqn. 9

      use constants,       only: one, two, three
      use cresp_variables, only: clight ! use units, only: clight
      use initcrspectrum,  only: e_small

      implicit none

      real(kind=8) :: p_ratio, f_ratio, q_in, n_cnst
      real(kind=8) :: n_func_2_zero_up

      if (abs(q_in - three) .lt. eps) then
         n_func_2_zero_up = - n_cnst + e_small / ((clight **two) * f_ratio * (p_ratio **three))* log(p_ratio)
      else
         n_func_2_zero_up = - n_cnst + e_small / ((clight **two) * f_ratio * (p_ratio **three)) &
                          & * ((p_ratio **(three-q_in) - one)/(three - q_in))
      endif

   end function n_func_2_zero_up
!----------------------------------------------------------------------------------------------------
   function n_func_2_zero_lo(p_ratio, n_cnst,q_in) ! from eqn. 9

      use constants,       only: one, two, three
      use cresp_variables, only: clight ! use units, only: clight
      use initcrspectrum,  only: e_small

      implicit none

      real(kind=8) :: p_ratio,  n_cnst, q_in
      real(kind=8) :: n_func_2_zero_lo

      if (abs(q_in - three) .lt. eps) then
         n_func_2_zero_lo = - n_cnst + (e_small / (clight **two)) * log(p_ratio)
      else
         n_func_2_zero_lo = - n_cnst + (e_small / (clight **two)) * ((p_ratio **(three-q_in) - one)/(three - q_in))
      endif

   end function n_func_2_zero_lo
!----------------------------------------------------------------------------------------------------
   function func_val_vec_up_bnd(x) ! called by cresp_crspectrum

      implicit none

      real(kind=8), dimension(ndim) :: x
      real(kind=8), dimension(ndim) :: func_val_vec_up_bnd

      func_val_vec_up_bnd(1) = fun3_up(x(1), x(2), p_im1, alpha)
      func_val_vec_up_bnd(2) = fun5_up(x(1), x(2), p_im1, n_in)

   end function func_val_vec_up_bnd

!----------------------------------------------------------------------------------------------------
   function func_val_vec_lo_bnd(x) ! called by cresp_crspectrum module via NR_get_solution_lo

      implicit none

      real(kind=8), dimension(ndim) :: x
      real(kind=8), dimension(ndim) :: func_val_vec_lo_bnd

      func_val_vec_lo_bnd(1) = fun3_lo(x(1), x(2), p_ip1, alpha)
      func_val_vec_lo_bnd(2) = fun5_lo(x(1), x(2), p_ip1, n_in)

   end function func_val_vec_lo_bnd

!---------------------------------------------------------------------------------------------------
! fun3_up - Alternative function introduced to find p_up using integrals of n, e and p_fix value.
!---------------------------------------------------------------------------------------------------
   function fun3_up(log10_p_r, log10_f_l, p_l, alpha ) ! used by func_val_vec_up_bnd to compute upper boundary p and f.                                 ! DEPRECATED ?

      use constants, only: three, four, one, ten

      implicit none

      real(kind=8),intent(in)   :: log10_p_r, log10_f_l
      real(kind=8) :: f_l, f_r, alpha, fun3_up, q_bin, p_l, p_r ! sought values will be x for single N-R and x, f_l_iup in NR-2dim

      p_r = ten **log10_p_r
      f_l = ten **log10_f_l
      f_r = e_small_to_f(p_r)
      q_bin = q_ratios(f_r/f_l,p_r/p_l)

      if ( abs(q_bin - three) .lt. eps) then
         fun3_up = -alpha + (-one + p_r/p_l)/log(p_r/p_l)
      else if ( abs(q_bin - four) .lt. eps) then
         fun3_up = -alpha + (p_r/p_l)*log(p_r/p_l)/(p_r/p_l - one)
      else
         fun3_up = -alpha + ((three - q_bin) / (four - q_bin)) * &
            & (((p_r/p_l)**((four - q_bin)) - one ) / ((p_r/p_l)**((three - q_bin)) - one ))
      endif

   end function fun3_up
!---------------------------------------------------------------------------------------------------
! fun5_up - Function similar to nq_to_f subroutine, used by compute_fp_NR_2dim to estimate value of p_up and f_l_iup
!---------------------------------------------------------------------------------------------------
   function fun5_up(log10_p_r, log10_f_l, p_l, n_in) ! used by func_val_vec_up_bnd to compute upper boundary p and f.                                 ! DEPRECATED ?

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
! fun3_lo - Alternative function introduced to find p_up using integrals of n, e and p_fix value.
!---------------------------------------------------------------------------------------------------
   function fun3_lo(log10_p_l, log10_f_r, p_r, alpha ) ! used by func_val_vec_lo_bnd to compute low boundary p and f.                                 ! DEPRECATED ?

      use constants, only: three, four, one, ten

      implicit none

      real(kind=8),intent(in)   :: log10_p_l, log10_f_r
      real(kind=8) :: f_l, f_r, alpha, fun3_lo, q_bin, p_r, p_l ! sought values will be x for single N-R and x, f_l_iup in NR-2dim

      p_l = ten **log10_p_l
      f_r = ten **log10_f_r
      f_l = e_small_to_f(p_l)
      q_bin = q_ratios(f_r/f_l,p_r/p_l)

      if ( abs(q_bin - three) .lt. eps) then
         fun3_lo = -alpha/p_l + (-one + p_r/p_l)/log(p_r/p_l)
      else if ( abs(q_bin - four) .lt. eps) then
         fun3_lo = -alpha/p_l + (p_r/p_l)*log(p_r/p_l)/(p_r/p_l - one)
      else
         fun3_lo = -alpha/p_l + ((three - q_bin) / (four - q_bin)) * &
            &      (((p_r/p_l)**((four - q_bin)) -one ) / ((p_r/p_l)**((three - q_bin)) -one ))
      endif

  end function fun3_lo
!---------------------------------------------------------------------------------------------------
! fun5_lo - Function similar to nq_to_f subroutine, used by compute_fp_NR_2dim to estimate value of p_lo and f_r_lo
!---------------------------------------------------------------------------------------------------
   function fun5_lo(log10_p_l, log10_f_r, p_r, n_in) ! used by func_val_vec_lo_bnd to compute low boundary p and f.                                 ! DEPRECATED ?

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

      use constants,       only: zero, three, two, fpi
      use cresp_variables, only: clight ! use units, only: clight
      use initcrspectrum,  only: e_small

      implicit none

      real (kind=8) :: e_small_to_f, p_outer

      e_small_to_f = zero
      e_small_to_f = e_small / (fpi * (clight **two)  * p_outer **three)

   end function e_small_to_f
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
   function lin_interpol_1D(loc_1, loc_2, value)

      implicit none

      integer(kind=4), intent(in) :: loc_1, loc_2
      real(kind=8),    intent(in) :: value
      real(kind=8)                :: lin_interpol_1D

      lin_interpol_1D = p_n(loc_1) + (value - p_a(loc_1)) * ( p_n(loc_1) - p_n(loc_2) ) / (p_a(loc_1) - p_a(loc_2)) ! WARNING - uses p_a and p_n, that are usually used to point alpha and n arrays.

   end function lin_interpol_1D
!----------------------------------------------------------------------------------------------------
   function lin_extrapol_1D(fun, arg, arg_out)

      implicit none

      real(kind=8), dimension(1:2), intent(in) :: fun, arg
      real(kind=8),                 intent(in) :: arg_out
      real(kind=8)                             :: lin_extrapol_1D

      lin_extrapol_1D = fun(1) + (fun(2) - fun(1)) * (arg_out - arg(1))/(arg(2)-arg(1))

   end function lin_extrapol_1D
!----------------------------------------------------------------------------------------------------
  function intpol_pf_from_NR_grids(a_val, n_val, interpolation_successful) ! for details see paragraph "Bilinear interpolation" in Numerical Recipes for F77, page 117, eqn. 3.6.4

      implicit none                                              ! should return exit code as well

      real(kind=8),     intent(inout) :: a_val, n_val  ! ratios arrays (p,f: lo and up), for which solutions have been obtained. loc_no_ip - in case when interpolation is not possible,
      logical,          intent(out)   :: interpolation_successful
      real(kind=8), dimension(2)      :: intpol_pf_from_NR_grids ! indexes with best match and having solutions are chosen.
      integer(kind=4), dimension(1:2) :: loc1, loc2, loc_no_ip ! loc1, loc2 - indexes that points where alpha_tab_ and up nad n_tab_ and up are closest in value to a_val and n_val - indexes point to
      logical                         :: exit_code

      exit_code = .false.

#ifdef CRESP_VERBOSED
      write (*,"(A30,A2,A4)",advance="no") "Determining indices for case: ", bound_name(current_bound), "... "
#endif /* CRESP_VERBOSED */
      call determine_loc(a_val, n_val, loc1, loc2, loc_no_ip, exit_code)

#ifdef CRESP_VERBOSED
      call save_loc(current_bound, loc1(1), loc1(2))
      call save_loc(current_bound, loc2(1), loc2(2))
#endif /* CRESP_VERBOSED */

      if (exit_code .eqv. .true. ) then ! interpolation won't work in this case, choosing closest values that have solutions.
         intpol_pf_from_NR_grids(1) = p_p(loc_no_ip(1),loc_no_ip(2)) ! this countermeasure wont work = loc_no_ip not initialized !
         intpol_pf_from_NR_grids(2) = p_f(loc_no_ip(1),loc_no_ip(2))
         interpolation_successful = .false.
         return
      else
         intpol_pf_from_NR_grids(1) = bl_interpol(p_p(loc1(1),loc1(2)),p_p(loc1(1),loc2(2)), &
               &  p_p(loc2(1),loc1(2)),p_p(loc2(1),loc2(2)), bl_in_tu(p_a(loc1(1)), a_val, p_a(loc2(1))), &
               &  bl_in_tu(p_n(loc1(2)), n_val, p_n(loc2(2))) )
         intpol_pf_from_NR_grids(2) = bl_interpol(p_f(loc1(1),loc1(2)),p_f(loc1(1),loc2(2)), &
               &  p_f(loc2(1),loc1(2)),p_f(loc2(1),loc2(2)), bl_in_tu(p_a(loc1(1)), a_val, p_a(loc2(1))), &
               &  bl_in_tu(p_n(loc1(2)), n_val, p_n(loc2(2))) )
         interpolation_successful = .true.
         return
      endif

   end function intpol_pf_from_NR_grids
!----------------------------------------------------------------------------------------------------
   subroutine determine_loc(a_val, n_val, loc1, loc2, loc_panic, exit_code)

      use constants,      only: zero
      use initcrspectrum, only: arr_dim

      implicit none

      real(kind=8),                    intent(inout) :: a_val, n_val
      integer(kind=4), dimension(1:2), intent(inout) :: loc1, loc2, loc_panic
      logical,                         intent(inout) :: exit_code
      logical                                        :: hit_zero !, no_solution

      hit_zero = .false.
      loc1(1) = inverse_f_to_ind(a_val, p_a(1), p_a(arr_dim), arr_dim)
      loc1(2) = inverse_f_to_ind(n_val, p_n(1), p_n(arr_dim), arr_dim)

      if ( (minval(loc1) .ge. 1 .and. maxval(loc1) .le. arr_dim-1)) then ! only need to test loc1
         if (p_p(loc1(1),loc1(2)) .gt. zero ) then
            loc2 = loc1+1
#ifdef CRESP_VERBOSED
            write(*,"(A19, 2I8, A3, 2I8)") "Obtained indices:", loc1, " | ", loc2
#endif /* CRESP_VERBOSED */
            return        ! normal exit
         else
            exit_code = .true.
            hit_zero  = .true.
         endif
      else
         exit_code = .true.
      endif

      if ( exit_code ) then ! namely if ((minval(loc1) .le. 0 .or. maxval(loc1) .ge. arr_dim))
         loc1(1) = max(1,min(loc1(1),arr_dim))   ! Here we either give algorithm closest nonzero value relative to a row
         loc1(2) = max(1,min(loc1(2),arr_dim))   ! that was in the proper range or we just feed the algorithm ANY nonzero
         exit_code = .true.                      ! initial vector that will prevent it from crashing.
         loc2 = loc1

         if (loc1(1) .eq. arr_dim .or. hit_zero) then
            call nearest_solution(p_p(:,loc1(2)), loc1(1), 1, loc1(1), hit_zero)
         endif
         if (loc1(1) .le. 1 .or. hit_zero) then
            call nearest_solution(p_p(:,loc1(2)), max(1,loc1(1)), arr_dim, loc1(1), hit_zero)
         endif
         if (loc1(2) .eq. arr_dim.or. hit_zero) then
            call nearest_solution(p_p(loc1(1),:), loc1(2), 1, loc1(2), hit_zero)
         endif
         if (loc1(2) .le. 1 .or. hit_zero) then
            call nearest_solution(p_p(loc1(1),:), max(1,loc1(2)), arr_dim, loc1(2), hit_zero)
         endif

         loc_panic = loc1
      endif

      loc2 = loc1 + 1

   end subroutine determine_loc
!----------------------------------------------------------------------------------------------------
   subroutine nearest_solution(arr_lin, i_beg, i_end, i_solution, hit_zero)

      use constants, only: zero

      implicit none

      real(kind=8), dimension(:) :: arr_lin
      integer(kind=4) :: i, i_beg, i_end, i_incr, i_solution
      logical, intent(inout) :: hit_zero

      i_solution = i_beg
      i_incr = sign(1, i_end - i_beg)

      do  i = i_beg, i_end, i_incr
         if (arr_lin(i) .gt. zero) then
            i_solution = i ! next one is i_solution + i_incr
            hit_zero = .false.
            return
         endif
      enddo

   end subroutine nearest_solution
!----------------------------------------------------------------------------------------------------
   function compute_q(alpha_in, exit_code, outer_p_ratio)

      use constants,      only: zero, one, I_ZERO, I_ONE
      use initcrspectrum, only: NR_refine_solution_q, q_big, p_fix_ratio, arr_dim_q

      implicit none

      real(kind=8) :: compute_q
      real(kind=8), optional :: outer_p_ratio
      real(kind=8), intent(inout) :: alpha_in
      integer(kind=4) :: loc_1, loc_2
      logical, intent(inout) :: exit_code ! value should be .true. at input

      p_a => alpha_tab_q
      p_n => q_grid
      selected_function_1D => alpha_to_q
      selected_value_check_1D => q_control

      compute_q = zero
      if (present(outer_p_ratio)) then
         p_ratio_4_q = outer_p_ratio
      else
         p_ratio_4_q = p_fix_ratio
      endif

      loc_1 = inverse_f_to_ind(alpha_in, alpha_tab_q(1), alpha_tab_q(arr_dim_q), arr_dim_q)

      if ((loc_1 .le. I_ZERO) .or. (loc_1 .ge. arr_dim_q)) then
         if (loc_1 .le. I_ZERO) then
            compute_q = q_big          !< should be consistent with q_grid(1)
         else
            compute_q = -q_big         !< should be consistent with q_grid(arr_dim_q)
         endif
         return                        ! < returns compute_q withh exit_code = .true.
      endif

      loc_2 = loc_1 + I_ONE
      compute_q = lin_interpol_1D(loc_1, loc_2, alpha_in)

      if (NR_refine_solution_q) then
         alpha = alpha_in
         call q_control(compute_q,exit_code)
         call NR_algorithm_1D(compute_q, exit_code)
      endif

      if (abs(compute_q) .gt. q_big) compute_q = sign(one, compute_q) * q_big

   end function compute_q
!----------------------------------------------------------------------------------------------------
   subroutine save_NR_guess_grid(NR_guess_grid, var_name)

      use cresp_variables, only: clight ! use units, only: clight
      use initcrspectrum,  only: e_small, q_big, max_p_ratio

      implicit none

      real(kind=8),dimension(:,:) :: NR_guess_grid
      integer(kind=4) :: j
      character(len=11) :: var_name
      character(len=4)  :: extension
      character(len=15) :: f_name

      extension =  ".dat"
      f_name = var_name // extension
      open(31, file=f_name, status="unknown", position="rewind")
         write(31,"(A57,A3,A105)") "This is a storage file for NR init grid, boundary case: ", var_name(9:), &
            &     "Saved below: e_small, size(NR_guess_grid,dim=1), size(NR_guess_grid,dim=2), max_p_ratio, q_big,  clight. &
            &      Do not remove content from this file"
         write(31, "(1E15.8, 2I10,10E22.15)") e_small, size(NR_guess_grid,dim=1), size(NR_guess_grid, dim=2), &
            &     max_p_ratio, q_big, clight              ! TODO: remove max_p_ratio, swap cols, rows with just arr_dim
         write(31, "(A1)") " "                            ! Blank line for

         do j=1, size(NR_guess_grid,dim=2)
            write(31, "(*(E24.15E3))") NR_guess_grid(:,j)  ! WARNING - THIS MIGHT NEED EXPLICIT INDICATION OF ELEMENTS COUNT IN LINE IN OLDER COMPILERS
         enddo
         close(31)

   end subroutine save_NR_guess_grid
!----------------------------------------------------------------------------------------------------
#ifdef CRESP_VERBOSED
   subroutine save_loc(bound_case, loc1, loc2)

      implicit none

      integer(kind=4), intent(in) :: bound_case, loc1, loc2
      character(len=10)           :: f_name

      f_name = "loc_"//bound_name(bound_case)//".dat"
      open(32, file=f_name, status="unknown", position="append")
      write (32,"(2I5)") loc1, loc2
      close(32)

   end subroutine save_loc
#endif /* CRESP_VERBOSED */
!----------------------------------------------------------------------------------------------------
   subroutine read_NR_guess_grid(NR_guess_grid, var_name, exit_code) ! must be improved, especially for cases when files do not exist

      use cresp_variables, only: clight ! use units, only: clight
      use func,            only: operator(.equals.)
      use initcrspectrum,  only: e_small, q_big, max_p_ratio

      implicit none

      real(kind=8),dimension(:,:) :: NR_guess_grid
      real(kind=8)                :: svd_e_sm, svd_max_p_r, svd_q_big, svd_clight
      integer(kind=4) :: j, svd_cols, svd_rows, file_status = 0
      character(len=11) :: var_name
      character(len=4)  :: extension
      character(len=15) :: f_name
      logical           :: exit_code

      extension =  ".dat"
      f_name = var_name // extension
      open(31, file=f_name, status="old", position="rewind",IOSTAT=file_status)
      if (file_status .gt. 0) then
         write(*,"(A8,I4,A8,2A20)") "IOSTAT:", file_status, ": file ", var_name//".dat"," does not exist!"
         exit_code = .true.
         return
      else
         read(31,*) ! Skipping comment line
         read(31,"(1E15.8,2I10,10E22.15)") svd_e_sm, svd_rows, svd_cols, svd_max_p_r, svd_q_big, svd_clight
         if ( (svd_e_sm .equals. e_small) .and. (svd_rows .eq. size(NR_guess_grid, dim=1)) .and. &
           &  (svd_cols .eq. size(NR_guess_grid, dim=2)) .and. &       ! TODO: swap rows and cols with just arr_dim, drop max_p_ratio
           &  (svd_max_p_r .equals. max_p_ratio) .and. (svd_q_big .equals. q_big) .and. (svd_clight .equals. clight) ) then ! This saves a lot of time
            read(31, *)

            do j=1, size(NR_guess_grid,dim=2)
               read(31, "(*(E24.15E3))") NR_guess_grid(:,j)  ! WARNING - THIS MIGHT NEED EXPLICIT INDICATION OF ELEMENTS COUNT IN LINE IN OLDER COMPILERS
            enddo
            exit_code = .false.
         else
            write(*,"(A61,A4,A6)") "Different initial parameters: will resolve ratio tables for ",var_name(10:)," case."
            write(*,"(2I10,10E22.15)") svd_cols, svd_rows, svd_e_sm, svd_max_p_r, svd_q_big, svd_clight
            write(*,"(2I10,10E22.15)") size(NR_guess_grid,dim=1), size(NR_guess_grid,dim=2),e_small,max_p_ratio,q_big,clight
            exit_code = .true.
         endif
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

!----------------------------------------------------------------------------------------------------
   function ind_to_flog(ind,min_in, max_in, length)

      use constants, only: I_ONE, ten

      implicit none

      real(kind=8) :: min_in, max_in
      integer(kind=4) :: ind, length
      real(kind=8) :: ind_to_flog

      ind_to_flog = min_in * ten**(((log10(max_in/min_in))/real(length-I_ONE,kind=8))*real((ind-I_ONE),kind=8))

   end function ind_to_flog
!----------------------------------------------------------------------------------------------------
   function inverse_f_to_ind(value, min_in, max_in, length) ! returns lower index for a given value, will need limiter

      use constants, only: I_ONE

      implicit none

      integer(kind=4) :: inverse_f_to_ind, length
      real(kind=8)    :: value, min_in, max_in

      inverse_f_to_ind = int((log10(value/min_in)/log10(max_in/min_in)) * (length - I_ONE )) + I_ONE

   end function inverse_f_to_ind
!----------------------------------------------------------------------------------------------------
   subroutine cresp_NR_mpi_exchange

      use mpisetup, only: piernik_MPI_Bcast

      implicit none

      call piernik_MPI_Bcast(p_ratios_lo)
      call piernik_MPI_Bcast(p_ratios_up)
      call piernik_MPI_Bcast(f_ratios_lo)
      call piernik_MPI_Bcast(f_ratios_up)
      call piernik_MPI_Bcast(q_grid)
      call piernik_MPI_Bcast(n_tab_lo)
      call piernik_MPI_Bcast(n_tab_up)
      call piernik_MPI_Bcast(alpha_tab_up)
      call piernik_MPI_Bcast(alpha_tab_lo)
      call piernik_MPI_Bcast(alpha_tab_q)

   end subroutine cresp_NR_mpi_exchange
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

      selected_function_2D => fvec_test
      x(1) = -5.0
      x(2) = 8.0

      print *, "starting with x =", x
      call NR_algorithm(x, exit_code)
      print *, x
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
!----------------------------------------------------------------------------------------------------
   subroutine nr_test_1D

      implicit none

      real(kind=8) :: x
      logical      :: exit_code

      print *, "NR algorithm test, provided test functions are:"
      print *, "f(x) = 1.1*x + sin(x)"
      print *, "expected solution is x = 0.0 "

      selected_value_check_1D  => dummy_check_1D
      selected_function_1D => f_test_1D
      x = -5.0

      print *, "starting with x =", x, "f(x) = ", selected_function_1D(x)
      call NR_algorithm_1D(x, exit_code)
      print *, x

   end subroutine nr_test_1D
!----------------------------------------------------------------------------------------------------
   function f_test_1D(x)

      implicit none

      real(kind=8) :: f_test_1D
      real(kind=8) :: x

      f_test_1D = 1.1*x + sin(x)

   end function f_test_1D
 !----------------------------------------------------------------------------------------------------
   subroutine dummy_check_1D(x, exit_code)

      implicit none

      real(kind=8) :: x
      logical      :: exit_code

      x = x ! to turn off unused argument warnings
      exit_code = .false.

   end subroutine dummy_check_1D
 !----------------------------------------------------------------------------------------------------
   subroutine initialize_arrays

      use diagnostics,    only: my_allocate_with_index, my_allocate, ma1d, ma2d
      use initcrspectrum, only: arr_dim, arr_dim_q

      implicit none

      ma1d = arr_dim
      ma2d = [arr_dim, arr_dim]

      call my_allocate_with_index(alpha_tab_lo, arr_dim, 1)
      call my_allocate_with_index(alpha_tab_up, arr_dim, 1)
      call my_allocate_with_index(n_tab_lo, arr_dim, 1)
      call my_allocate_with_index(n_tab_up, arr_dim, 1)
      call my_allocate_with_index(alpha_tab_q, arr_dim_q, 1)
      call my_allocate_with_index(q_grid, arr_dim_q, 1)
      call my_allocate(p_ratios_lo, ma2d )
      call my_allocate(f_ratios_lo, ma2d )
      call my_allocate(p_ratios_up, ma2d )
      call my_allocate(f_ratios_up, ma2d )

   end subroutine initialize_arrays
!----------------------------------------------------------------------------------------------------

end module cresp_NR_method
