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

   use constants, only: LO, HI

   implicit none

   private
   public :: alpha, assoc_pointers, bound_name, cresp_initialize_guess_grids, compute_q, intpol_pf_from_NR_grids, n_in, NR_algorithm, q_ratios

   integer, parameter                        :: ndim = 2
   real, allocatable, dimension(:)           :: p_space, q_space
   real                                      :: alpha, p_ratio_4_q, n_in
   real, allocatable, dimension(:),   target :: alpha_tab_lo, alpha_tab_up, n_tab_lo, n_tab_up, alpha_tab_q, q_grid
   real, allocatable, dimension(:,:), target :: p_ratios_lo, f_ratios_lo, p_ratios_up, f_ratios_up
   integer(kind=4)                           :: helper_arr_dim
   real, pointer, dimension(:)               :: p_a => null(), p_n => null() ! pointers for alpha_tab_(lo,up) and n_tab_(lo,up) or optional - other 1-dim arrays
   real, pointer, dimension(:,:)             :: p_p => null(), p_f => null() ! pointers for p_ratios_(lo,up) and f_ratios_(lo,up)
#ifdef CRESP_VERBOSED
   integer(kind=4)                           :: current_bound, sought_by
   integer(kind=4), parameter                :: SLV = 1, RFN = 2
#endif /* CRESP_VERBOSED */
   integer, parameter                               :: blen = 2, extlen = 4, flen = 15
   character(len=blen), dimension(LO:HI), parameter :: bound_name = ['lo', 'up']
   character(len=extlen), parameter                 :: extension =  ".dat"

   abstract interface
      function function_pointer_2D(z)
         real, dimension(2), intent(inout) :: z
         real, dimension(2)                :: function_pointer_2D
      end function function_pointer_2D
   end interface

   procedure (function_pointer_2D), pointer :: selected_function_2D => null()

   type     map_header
      integer  :: s_dim1, s_dim2
      real     :: s_es
      real     :: s_pr
      real     :: s_qbig
      real     :: s_c
      real     :: s_amin, s_amax, s_nmin, s_nmax
   end type map_header

!----------------------------------------------------------------------------------------------------

contains

   subroutine NR_algorithm(x, exit_code)

      use initcrspectrum, only: NR_iter_limit, tol_f, tol_x, eps_det

      implicit none

      real, dimension(ndim), intent(inout) :: x
      logical,               intent(out)   :: exit_code
      real, dimension(ndim)                :: fun_vec_value
      real, dimension(1:ndim)              :: cor
      real, dimension(size(x),size(x))     :: fun_vec_jac, fun_vec_inv_jac
      real                                 :: err_f, err_x, det
      integer(kind=2)                      :: i

      err_f = tol_f
      err_x = tol_x
      exit_code = .true.

      fun_vec_value = selected_function_2D(x)
      if (maxval(abs(fun_vec_value)) < 0.01 * err_f) then ! in case when f converges at initialization
         exit_code = .false.
         return
      endif

      do i = 1, NR_iter_limit
         if (maxval(abs(fun_vec_value)) < err_f) then    ! For convergence via value of f
            exit_code = .false.
#ifdef CRESP_VERBOSED
            write(*,"(A47,I4,A12)",advance="no") "Convergence via value of fun_vec_value after ",i, " iterations."!, x, fun_vec_value
#endif /* CRESP_VERBOSED */
            return
         endif

         fun_vec_value = selected_function_2D(x)
         fun_vec_jac = jac_fin_diff(x)                    ! function vector already explicitly provided to jac_fin_diff (finite difference method)

         det = determinant_2d_real(fun_vec_jac)           ! WARNING - algorithm is used for ndim = 2. For more dimensions LU or other methods should be implemented.
         if (abs(det) < eps_det) then              ! Countermeasure against determinant = zero
            exit_code = .true.
            return
         endif
         fun_vec_inv_jac = invert_2d_matrix(fun_vec_jac,det)

         cor(1) = fun_vec_inv_jac(1,1) * fun_vec_value(1) + fun_vec_inv_jac(1,2) * fun_vec_value(2)
         cor(2) = fun_vec_inv_jac(2,1) * fun_vec_value(1) + fun_vec_inv_jac(2,2) * fun_vec_value(2)
         x = x+cor
         if (maxval(abs(cor)) < err_x) then                 ! For convergence via value of correction (cor) table.
#ifdef CRESP_VERBOSED
            write(*,"(A47,I4,A12)",advance="no") "Convergence via value of cor array     after ",i," iterations."
#endif /* CRESP_VERBOSED */
            exit_code = .false.
            return
         endif
      enddo
      ! "  ... WARNING: Maximum number of iterations (",NR_iter_limit,") exceeded @global_newt!"
      exit_code = .true.

   end subroutine NR_algorithm

!----------------------------------------------------------------------------------------------------

   subroutine NR_algorithm_1D(x, exit_code)

      use constants,      only: zero, big
      use func,           only: operator(.equals.)
      use initcrspectrum, only: NR_iter_limit, tol_f_1D, tol_x_1D

      implicit none

      integer :: i
      real    :: x, delta, dfun_1D, fun1D_val
      logical :: exit_code, func_check

      delta = big
      func_check = .false.

      do i = 1, NR_iter_limit
         fun1D_val = alpha_to_q(x)
         if ((abs(fun1D_val) <= tol_f_1D) .and. (abs(delta) <= tol_f_1D)) then ! delta <= tol_f acceptable as in case of f convergence we must only check if the algorithm hasn't wandered astray
            exit_code = .false.
            return
         endif

         dfun_1D = derivative_1D(x)

         if (abs(dfun_1D) .equals. zero) then
            exit_code = .true.
            return
         endif

         delta = fun1D_val / dfun_1D

         x = x - delta
         if (abs(delta) < tol_x_1D) then
            exit_code = .false.
            return
         endif
         call q_control(x, func_check)  ! necessary in some cases, when maximal value x can take is defined, for other cases dummy_check_1D function defined
         if (func_check) return
      enddo

      exit_code = .true. ! if all fails

   end subroutine NR_algorithm_1D

!----------------------------------------------------------------------------------------------------

   real function derivative_1D(x) ! via finite difference method

      use constants,      only: half
      use initcrspectrum, only: eps

      implicit none

      real, intent(in) :: x
      real             :: dx
      real, parameter  :: dx_par = 1.0e-4

      dx = sign(1.0, x) * min(abs(x*dx_par), dx_par)
      dx = sign(1.0, x) * max(abs(dx), eps) ! dx = 0.0 must not be allowed
      derivative_1D = half * (alpha_to_q(x+dx) - alpha_to_q(x-dx))/dx

   end function derivative_1D

!----------------------------------------------------------------------------------------------------

   subroutine cresp_initialize_guess_grids

      use constants,       only: zero, I_FOUR
      use cresp_variables, only: clight_cresp
      use initcrspectrum,  only: e_small, q_big, max_p_ratio, arr_dim, arr_dim_q
      use mpisetup,        only: master

      implicit none

      logical      :: first_run = .true., save_to_log = .false.
      character(8) :: date
      character(9) :: time

      call initialize_arrays
      if (master .and. first_run) then
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

         print *, "Are there zeros? (q_ratios)",    count(abs(q_space) <= zero)
         print *, "Are there zeros? (p_ratios_up)", count(p_ratios_up  <= zero)
         print *, "Are there zeros? (f_ratios_up)", count(f_ratios_up  <= zero)
         print *, "Are there zeros? (p_ratios_lo)", count(p_ratios_lo  <= zero)
         print *, "Are there zeros? (f_ratios_lo)", count(f_ratios_lo  <= zero)
         print *, "Count of array elements:", size(p_ratios_lo)
         print *,"----------"
         if (save_to_log) then
            open(15, file="log_NR_solve", position = "append")
            write (15,*) "------------------------------------------"
            write (15,"(A,2x,A,2x)") "Run on: ", date, "at: ", time
            write (15,*) "For set of parameters: e_small, size(NR_guess_grid,dim=1), size(NR_guess_grid, dim=2), max_p_ratio, q_big, clight_cresp"
            write (15, "(1E15.8, 2I10,10E22.15)") e_small, size(p_space), size(q_space), max_p_ratio, q_big, clight_cresp
            write (15,*) "Are there zeros? (q_ratios)",    count(abs(q_space) <= zero)
            write (15,*) "Are there zeros? (p_ratios_up)", count(p_ratios_up  <= zero), real(count(p_ratios_up <= zero))/real(size(p_ratios_up)) * 100.0,"%"
            write (15,*) "Are there zeros? (f_ratios_up)", count(f_ratios_up  <= zero), real(count(f_ratios_up <= zero))/real(size(f_ratios_up)) * 100.0,"%"
            write (15,*) "Are there zeros? (p_ratios_lo)", count(p_ratios_lo  <= zero), real(count(p_ratios_lo <= zero))/real(size(p_ratios_lo)) * 100.0,"%"
            write (15,*) "Are there zeros? (f_ratios_lo)", count(f_ratios_lo  <= zero), real(count(f_ratios_lo <= zero))/real(size(f_ratios_lo)) * 100.0,"%"
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

      use constants,      only: zero, half, one, three, I_ONE, big, small
      use dataio_pub,     only: die,msg, printinfo, warn
      use initcrspectrum, only: q_big, force_init_NR, NR_run_refine_pf, p_fix_ratio, e_small_approx_init_cond, arr_dim, arr_dim_q, max_p_ratio, e_small
      use cresp_variables,only: clight_cresp

      implicit none

      integer(kind=4)      :: i, j, ilim = 0, qmaxiter = 100
      logical              :: read_error, headers_match
      real                 :: a_min_q = big, a_max_q = small , q_in3, pq_cmplx
      real, dimension(2)   :: a_min = big, a_max = small, n_min = big, n_max = small
      type(map_header)     :: hdr_init, hdr_read

      q_space = zero
      do i = 1, int(half*helper_arr_dim)
         q_space(i) = ln_eval_array_val(i, q_big, real(0.05), int(1,kind=4), int(half*helper_arr_dim,kind=4)) ! BEWARE: magic number
      enddo

      do i = 1, int(half*helper_arr_dim)!, arr_dim
         q_space(int(half*helper_arr_dim,kind=4)+i) = -q_space(int(half*helper_arr_dim,kind=4)+1-i)
      enddo

! setting up a grids of ratios to be used as phase space for NR tabs, obtained later
      do i = 1, helper_arr_dim
         p_space(i) = max_p_ratio**(real(i)/real(helper_arr_dim)) ! ind_to_flog(i, 1.000001, max_p_ratio) ! max_p_ratio**(real(i)/real(arr_dim))
      enddo
      do i = 1, helper_arr_dim
         do j = 1, helper_arr_dim
            q_in3 = three - q_space(j)
            pq_cmplx = p_space(i)**q_in3

            a_min(LO) = min(a_min(LO), abs(encp_func_2_zero(LO, p_space(i),           zero, q_in3)))
            n_min(LO) = min(n_min(LO), abs(   n_func_2_zero(    p_space(i), one,      zero, q_in3)))
            a_min(HI) = min(a_min(HI), abs(encp_func_2_zero(HI, p_space(i),           zero, q_in3)))
            n_min(HI) = min(n_min(HI), abs(   n_func_2_zero(    p_space(i), pq_cmplx, zero, q_in3)))

            a_max(LO) = max(a_max(LO), abs(encp_func_2_zero(LO, p_space(i),           zero, q_in3)))
            n_max(LO) = max(n_max(LO), abs(   n_func_2_zero(    p_space(i), one,      zero, q_in3)))
            a_max(HI) = max(a_max(HI), abs(encp_func_2_zero(HI, p_space(i),           zero, q_in3)))
            n_max(HI) = max(n_max(HI), abs(   n_func_2_zero(    p_space(i), pq_cmplx, zero, q_in3)))
         enddo
      enddo

!       a_min(LO) = 0.2       ! BEWARE: magic numbers!
!       a_max(LO) = 0.999999
!       a_min(HI) = 1.000005
!       a_max(HI) = 200.0
!       n_min(LO) = 1.0e-11
!       n_max(LO) = 5000.0
!       n_min(HI) = 1.0e-12
!       n_max(HI) = 1000.0

      hdr_init%s_es     = e_small
      hdr_init%s_dim1   = arr_dim
      hdr_init%s_dim2   = arr_dim
      hdr_init%s_pr     = max_p_ratio
      hdr_init%s_c      = clight_cresp


      do i = 1, arr_dim
         alpha_tab_lo(i) = ind_to_flog(i, a_min(LO), a_max(LO), arr_dim) ! a_min_lo * ten**((log10(a_max_lo/a_min_lo))/real(arr_dim-1)*real(i-1))
         alpha_tab_up(i) = ind_to_flog(i, a_min(HI), a_max(HI), arr_dim) ! a_min_up * ten**((log10(a_max_up/a_min_up))/real(arr_dim-1)*real(i-1))
         n_tab_lo(i)     = ind_to_flog(i, n_min(LO), n_max(LO), arr_dim) ! n_min_lo * ten**((log10(n_max_lo/n_min_lo))/real(arr_dim-1)*real(i-1))
         n_tab_up(i)     = ind_to_flog(i, n_min(HI), n_max(HI), arr_dim) ! n_min_up * ten**((log10(n_max_up/n_min_up))/real(arr_dim-1)*real(i-1))
      enddo

      if (e_small_approx_init_cond == 1) then

         hdr_init%s_amin   = alpha_tab_up(1)
         hdr_init%s_amax   = alpha_tab_up(arr_dim)
         hdr_init%s_nmin   = n_tab_up(1)
         hdr_init%s_nmin   = n_tab_up(arr_dim)

         write (msg, "(A47,A2,A10)") "[cresp_NR_method] Preparing solution maps for (",bound_name(HI), ") boundary"
         call printinfo(msg)

         call read_NR_smap_header("p_ratios_"//bound_name(HI), hdr_read, read_error)
         if (.not. read_error) then
            call check_NR_smap_header(hdr_read, hdr_init, headers_match)
            if (headers_match .and. (.not. force_init_NR)) then
               call read_NR_smap(p_ratios_up, "p_ratios_", HI, read_error)
               call read_NR_smap(f_ratios_up, "f_ratios_", HI, read_error)
            endif
         endif
         if (read_error) then
            write(msg,"(A,A,A)") "[cresp_NR_method] Problem reading data for (",bound_name(HI), ") boundary"
            call warn(msg)
         endif
         if (force_init_NR .or. (read_error .or. .not. headers_match) ) then
            call fill_boundary_grid(HI, p_ratios_up, f_ratios_up)
         endif

         if (NR_run_refine_pf) then
            call assoc_pointers(HI)
            call refine_all_directions(HI)
         endif

         call save_NR_smap(p_ratios_up, hdr_init, "p_ratios_", HI)
         call save_NR_smap(f_ratios_up, hdr_init, "f_ratios_", HI)

!--------------------
         hdr_init%s_amin   = alpha_tab_lo(1)
         hdr_init%s_amax   = alpha_tab_lo(arr_dim)
         hdr_init%s_nmin   = n_tab_lo(1)
         hdr_init%s_nmin   = n_tab_lo(arr_dim)

         write (msg, "(A47,A2,A10)") "[cresp_NR_method] Preparing solution maps for (",bound_name(HI), ") boundary"
         call printinfo(msg)

         call read_NR_smap_header("p_ratios_"//bound_name(LO), hdr_read, read_error)
         if (.not. read_error) then
            call check_NR_smap_header(hdr_read, hdr_init, headers_match)
            if (headers_match .and. (.not. force_init_NR)) then
               call read_NR_smap(p_ratios_lo, "p_ratios_", LO, read_error)
               call read_NR_smap(f_ratios_lo, "f_ratios_", LO, read_error)
            endif
         endif
         if (read_error) then
            write(msg,"(A,A,A)") "[cresp_NR_method] Problem reading data for (",bound_name(HI), ") boundary"
            call warn(msg)
         endif
         if (force_init_NR .or. (read_error .or. .not. headers_match) ) then
            call fill_boundary_grid(LO, p_ratios_lo, f_ratios_lo)
         endif

         if (NR_run_refine_pf) then
            call assoc_pointers(LO)
            call refine_all_directions(LO)
         endif

         call save_NR_smap(p_ratios_lo, hdr_init, "p_ratios_", LO)
         call save_NR_smap(f_ratios_lo, hdr_init, "f_ratios_", LO)
      endif

      a_min_q = one  + epsilon(one)
      a_max_q = (one + epsilon(one)) * p_fix_ratio
      j = min(arr_dim_q - int(arr_dim_q/100, kind=4), arr_dim_q - I_ONE)               ! BEWARE: magic number

      do while ((q_grid(j) <= (-q_big) .and. (q_grid(arr_dim_q) <= (-q_big))) .and. (ilim .le. qmaxiter) )
         a_max_q = a_max_q - a_max_q*0.005                                             ! BEWARE: magic number
         do i = 1, arr_dim_q
            alpha_tab_q(i)  = ind_to_flog(i, a_min_q, a_max_q, arr_dim_q)
         enddo
         call fill_q_grid(i_incr=1) ! computing q_grid takes so little time, that saving the grid is not necessary.
         ilim = ilim + 1
      enddo
      if (ilim .ge. qmaxiter) call die ("[cresp_NR_method:fill_guess_grids] Maximal iteration limit exceeded, q_grid might not have converged!")
#ifdef CRESP_VERBOSED
      do i = 1, arr_dim_q
         print "(A1,I3,A7,2F18.12)", "[ ", i,"] a : q ", q_grid(i), alpha_tab_q(i)
      enddo

      print *,"alpha_tab_lo(i),      alpha_tab_up(i),        n_tab_lo(i),        n_tab_up(i)  |       p_space(i),     q_space(i)"
      do i = 1, arr_dim
         if (i <= helper_arr_dim) then
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

      use dataio_pub, only: die, msg, printinfo

      implicit none

      integer(kind=4), intent(in) :: bound_case

#ifdef CRESP_VERBOSED
      sought_by = RFN
#endif /* CRESP_VERBOSED */

      write(msg,'(3a)') "Running refine for:", bound_name(bound_case), " boundary"
      call printinfo(msg)
      if (.not. allocated(p_space) .or. .not. allocated(q_space)) call die("[cresp_NR_method:refine_all_directions] refine_grids called after array deallocation, stopping")

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
   real function ln_eval_array_val(i, arr_min, arr_max, min_i, max_i)

      implicit none

      real            :: b, arr_min, arr_max
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
   subroutine fill_boundary_grid(bound_case, fill_p, fill_f) ! TODO FIXME to be paralelized

      use constants,      only: zero
      use dataio_pub,     only: msg, printinfo
      use initcrspectrum, only: arr_dim, eps

      implicit none

      integer(kind=4), intent(in) :: bound_case ! HI or LO
      real, dimension(:,:)        :: fill_p, fill_f
      real, dimension(1:2)        :: x_vec, prev_solution, prev_solution_1, x_step
      integer(kind=4)             :: i, j, is, js, jm
      logical                     :: exit_code, new_line
#ifdef CRESP_VERBOSED
      real, dimension(1:2)        :: x_in

      sought_by = SLV
#endif /* CRESP_VERBOSED */

      prev_solution(1) = p_space(1)
      prev_solution(2) = p_space(1)**q_space(1)
      prev_solution_1 = prev_solution

      call assoc_pointers(bound_case)

      fill_p = zero ; fill_f = zero
      x_step = zero
      write(msg, "(A,A2,A,I3,A)") "[cresp_NR_method:fill_boundary_grid] Solving solution maps for cutoff case (",bound_name(bound_case),"): DIM=",arr_dim,"**2"
      call printinfo(msg)

      do i = 1, arr_dim
         call add_dot( i .eq. arr_dim )
         new_line = .true.
         prev_solution = prev_solution_1 ! easier to find when not searching from the top
         do j = 1, arr_dim
         ! j_incr = 1
            alpha = p_a(i)
            n_in  = p_n(j)
#ifdef CRESP_VERBOSED
            write(*,"(A14,A2,A2,2I4,A9,I4,A1)",advance="no") "Now solving (",bound_name(bound_case),") ",i,j,", sized ",arr_dim," "
#endif /* CRESP_VERBOSED */

            call seek_solution_prev(fill_p(i,j), fill_f(i,j), prev_solution, exit_code)

            if (.not. exit_code .and. new_line) then
               prev_solution_1 = prev_solution
               new_line = .false.
            endif

            if (exit_code) then
               jm = j-2
               if (check_dimm(jm)) call step_extr(fill_p(i,jm:j), fill_f(i,jm:j), p_n(jm:j), exit_code)
               if (j >= 2) then
                  jm = j-1
                  if (fill_p(i,jm) > zero) call seek_solution_step(fill_p(i,j), fill_f(i,j), prev_solution, i, jm, exit_code)
               endif
            endif
            if (exit_code) then !still...
               do is = 1, helper_arr_dim
                  do js = 1, helper_arr_dim
                     x_vec(1) = p_space(is)
                     x_vec(2) = p_space(is)**(-q_space(js))
#ifdef CRESP_VERBOSED
                     x_in = x_vec
#endif /* CRESP_VERBOSED */
                     if (exit_code) then
                        call NR_algorithm(x_vec, exit_code)
                        if (.not. exit_code) then
                           fill_p(i,j) = x_vec(1) ! i index - alpha, j index - n_in
                           fill_f(i,j) = x_vec(2)
                           prev_solution = x_vec
#ifdef CRESP_VERBOSED
                           call msg_success("    ", x_in, x_vec)
#endif /* CRESP_VERBOSED */
                           exit
                        endif
                     endif
                  enddo
               enddo
            else
               if (prev_solution(1) <= eps) then
                  prev_solution(1) = prev_solution_1(1)
               else if (prev_solution(2) <= eps) then
                  prev_solution(2) = prev_solution_1(2)
               else
                  prev_solution(1) = fill_p(i,j)
                  prev_solution(2) = fill_f(i,j)
               endif
            endif
#ifdef CRESP_VERBOSED
            if (exit_code) print *,""
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
   subroutine step_seek(x_step, prev_sol, ii, jj, i_sol, j_sol, exit_code, nssstep)

      use constants, only: I_ONE

      implicit none

      real, dimension(1:2), intent(in)    :: prev_sol
      real, dimension(1:2), intent(out)   :: x_step
      integer(kind=4),      intent(in)    :: ii, jj, i_sol, j_sol, nssstep
      logical,              intent(inout) :: exit_code
      integer(kind=4)                     :: i, j

      do i = -1, 1
         do j = -1, 1
            if (exit_code) then
               x_step(1) = prev_sol(1) + (p_space(max(min(i_sol+i,helper_arr_dim),I_ONE)) - prev_sol(1)) / real(nssstep - ii + 0.1)
               x_step(2) = prev_sol(2) + (p_space(max(min(i_sol+i,helper_arr_dim),I_ONE))**(-q_space(max(min(j_sol+j, helper_arr_dim), 1))) - prev_sol(2)) / real(nssstep - jj + 0.1)
               call NR_algorithm(x_step, exit_code)
               if (.not. exit_code) return
            endif
         enddo
      enddo

   end subroutine step_seek

!----------------------------------------------------------------------------------------------------

   subroutine refine_ji(ref_p, ref_f, i_incr, j_incr) ! ref_f and ref_p should already be partially filled with solutions

      use constants,      only: zero

      implicit none

      real, dimension(:,:), intent(inout) :: ref_p, ref_f
      integer(kind=4),      intent(in)    :: i_incr, j_incr
      integer(kind=4)                     :: i, j, i_beg, i_end, j_beg, j_end, i1m, i2m, i1p
      real, dimension(1:2)                :: prev_solution
      logical                             :: exit_code

      prev_solution(1) = p_space(1)              ! refine must be called before these are deallocated
      prev_solution(2) = p_space(1)**q_space(1)
      call prepare_indices(i_incr, i_beg, i_end)
      call prepare_indices(j_incr, j_beg, j_end)
      do j = j_beg, j_end, j_incr
         do i = i_beg, i_end, i_incr
            alpha = p_a(i)
            n_in  = p_n(j)
            if (ref_p(i,j) > zero .and. ref_f(i,j) > zero) then
               prev_solution(1) = ref_p(i,j)
               prev_solution(2) = ref_f(i,j)
            else
               call seek_solution_prev(ref_p(i,j), ref_f(i,j), prev_solution, exit_code) ! works for most cases
               if (exit_code) then
                  i1m = i-i_incr ; i2m = i-2*i_incr ; i1p = i+i_incr
                  if (check_dimm(i2m)                      ) call step_extr(ref_p(i2m:i  :i_incr,j), ref_f(i2m:i  :i_incr,j),         p_a(i2m:i  :i_incr), exit_code)
                  if (check_dimm(i1m) .and. check_dimm(i1p)) call step_inpl(ref_p(i1m:i1p:i_incr,j), ref_f(i1m:i1p:i_incr,j), i_incr, p_a(i1m:i1p:i_incr), exit_code)
               endif
            endif
         enddo
      enddo

   end subroutine refine_ji

!----------------------------------------------------------------------------------------------------

   subroutine refine_ij(ref_p, ref_f, i_incr, j_incr) ! ref_f and ref_p should already be partially filled with solutions

      use constants,      only: zero

      implicit none

      real, dimension(:,:), intent(inout) :: ref_p, ref_f
      integer(kind=4),      intent(in)    :: i_incr, j_incr
      integer(kind=4)                     :: i, j, i_beg, i_end, j_beg, j_end, j1m, j2m, j1p
      real, dimension(1:2)                :: prev_solution
      logical                             :: exit_code

      prev_solution(1) = p_space(1)              ! refine must be called before these are deallocated
      prev_solution(2) = p_space(1)**q_space(1)
      call prepare_indices(i_incr, i_beg, i_end)
      call prepare_indices(j_incr, j_beg, j_end)
      do i = i_beg, i_end, i_incr
         do j = j_beg, j_end, j_incr
            alpha = p_a(i)
            n_in  = p_n(j)
            if (ref_p(i,j) > zero .and. ref_f(i,j) > zero) then
               prev_solution(1) = ref_p(i,j)
               prev_solution(2) = ref_f(i,j)
            else
               call seek_solution_prev(ref_p(i,j), ref_f(i,j), prev_solution, exit_code) ! works for most cases
               if (exit_code) then
                  j1m = j-j_incr ; j2m = j-2*j_incr ; j1p = j+j_incr
                  if (check_dimm(j2m)                      ) call step_extr(ref_p(i,j2m:j  :j_incr), ref_f(i,j2m:j  :j_incr),         p_n(j2m:j  :j_incr), exit_code)
                  if (check_dimm(j1m) .and. check_dimm(j1p)) call step_inpl(ref_p(i,j1m:j1p:j_incr), ref_f(i,j1m:j1p:j_incr), j_incr, p_n(j1m:j1p:j_incr), exit_code)
               endif
            endif
         enddo
      enddo

   end subroutine refine_ij

 !----------------------------------------------------------------------------------------------------

   logical function check_dimm(ind)

      use initcrspectrum, only: arr_dim

      implicit none

      integer(kind=4), intent(in) :: ind

      check_dimm = (ind >= 1 .and. ind <= arr_dim)

   end function check_dimm

 !----------------------------------------------------------------------------------------------------

   subroutine step_extr(p3, f3, arg, exit_code)

      use constants, only: zero

      implicit none

      real, dimension(1:3), intent(inout) :: p3, f3
      real, dimension(1:3), intent(in)    :: arg
      logical,              intent(out)   :: exit_code
      real, dimension(1:2)                :: x_vec_0, x_vec, delta, x_in
      integer(kind=4)                     :: nsubstep = 100, k
!         alpha and n are set !

      if (minval(p3(1:2)) > tiny(zero) .and. p3(3) <= zero) then ! sometimes NaNs and numbers of order e-317 appear; must be looked into
         x_vec_0 = [p3(2), f3(2)]
         delta(1) = lin_extrapol_1D(p3(1:2), arg(1:2), arg(3)) - p3(2) ! direction is not relevant in this case
         delta(2) = lin_extrapol_1D(f3(1:2), arg(1:2), arg(3)) - f3(2)
         delta = delta/nsubstep
         do k = 0, nsubstep
            x_vec = x_vec_0 + delta * k ! first iteration is a simple extrapolation
            x_in = x_vec
            call NR_algorithm(x_vec, exit_code)
            if (.not. exit_code) then
               x_vec = abs(x_vec)
#ifdef CRESP_VERBOSED
               call msg_success("extr", x_in, x_vec)
#endif /* CRESP_VERBOSED */
               p3(3) = x_vec(1)
               f3(3) = x_vec(2)
               return
            endif
         enddo
      endif

   end subroutine step_extr

!----------------------------------------------------------------------------------------------------

   subroutine step_inpl(p3, f3, incr, args, exit_code)

      use constants, only: zero

      implicit none

      real, dimension(1:3), intent(inout) :: p3, f3
      integer(kind=4),      intent(in)    :: incr
      real, dimension(1:3), intent(in)    :: args
      logical,              intent(inout) :: exit_code
      real, dimension(1:2)                :: x_vec, x_vec_0, delta, x_in
      integer(kind=4)                     :: k, nsubstep = 100
!         alpha and n are set !

      if (exit_code) then
         if (min(p3(1),p3(3)) > tiny(zero) .and. p3(2) <= zero) then ! sometimes NaNs and numbers of order e-317 appear; must be looked into
            x_vec_0 = [p3(1), f3(1)]
            delta(1) = lin_interpolation_1D( [p3(2-incr), p3(2+incr)], [args(2-incr), args(2+incr)], args(2) ) - p3(1)
            delta(2) = lin_interpolation_1D( [f3(2-incr), f3(2+incr)], [args(2-incr), args(2+incr)], args(2) ) - f3(1)
            x_in = x_vec_0 + delta ! gives the interpolated value as starting one
            delta = delta/nsubstep
            do k = 0, nsubstep
               x_vec = x_vec_0 + delta * k
               x_in = x_vec
               call NR_algorithm(x_vec, exit_code)
               if (.not. exit_code) then ! first iteration is a simple extrapolation
                  x_vec = abs(x_vec)
#ifdef CRESP_VERBOSED
                  call msg_success("inpl", x_in, x_vec)
#endif /* CRESP_VERBOSED */
                  p3(2) = x_vec(1)
                  f3(2) = x_vec(2)
                  return
               endif
               x_vec = x_vec + delta
            enddo
         endif
      endif

   end subroutine step_inpl
!----------------------------------------------------------------------------------------------------
#ifdef CRESP_VERBOSED
   subroutine msg_success(met_name, x_in, x_out)

      implicit none

      real, dimension(1:), intent(in) :: x_in
      real, dimension(1:), intent(in) :: x_out
      character(len=*),    intent(in) :: met_name
      integer, parameter                      :: slen = 6
      character(len=slen), dimension(SLV:RFN) :: sought = ['Solve ', 'Refine']

      write (*, "(A6,A13,2E16.9)",advance="no") sought(sought_by)," (alpha, n): ",alpha,n_in
      write (*, "(A5,A4,A42, 2E19.10e3)",advance="no") " -> (",met_name,") solution obtained, (p_ratio, f_ratio) = ", x_out
      write (*, "(A21, 2E17.10)",advance="no") ", provided input:", x_in ; print *,""

   end subroutine msg_success
#endif /* CRESP_VERBOSED */
!----------------------------------------------------------------------------------------------------
   function lin_interpolation_1D(fun, arg, arg_mid)

      use constants, only: one

      implicit none

      real, dimension(1:2), intent(in) :: arg, fun
      real,                 intent(in) :: arg_mid
      real :: weight, lin_interpolation_1D

      weight   = (arg_mid - arg(1)) / (arg(2) - arg(1))
      lin_interpolation_1D =  fun(1) * (one - weight) + fun(2) * (one - weight)

   end function lin_interpolation_1D
!----------------------------------------------------------------------------------------------------
   subroutine seek_solution_prev(p2ref, f2ref, prev_solution, exit_code)

      implicit none

      real,                 intent(inout) :: p2ref, f2ref
      real, dimension(1:2), intent(inout) :: prev_solution
      logical,              intent(out)   :: exit_code
      real, dimension(1:2)                :: x_vec

      x_vec = prev_solution
      call NR_algorithm(x_vec, exit_code)
      if (.not. exit_code) then
         x_vec = abs(x_vec)
         p2ref = x_vec(1)
         f2ref = x_vec(2)
#ifdef CRESP_VERBOSED
         call msg_success("prev", prev_solution, x_vec)
#endif /* CRESP_VERBOSED */
         prev_solution = x_vec
         return
      endif

   end subroutine seek_solution_prev
!----------------------------------------------------------------------------------------------------
   subroutine seek_solution_step(p2ref, f2ref, prev_solution, i_obt, j_obt, exit_code)

      implicit none

      real,                 intent(out)   :: p2ref, f2ref
      real, dimension(1:2), intent(inout) :: prev_solution
      integer(kind=4),      intent(in)    :: i_obt, j_obt
      logical,              intent(inout) :: exit_code
      real, dimension(1:2)                :: x_vec, x_step
      integer(kind=4)                             :: ii, jj, nssstep = 3
!    alpha and n are set !

      if (exit_code) then
         do ii = 0, nssstep
            do jj = 0,nssstep
               call step_seek(x_step, prev_solution, ii, jj, i_obt, j_obt, exit_code, nssstep)
               x_vec = x_step
               if (.not. exit_code) then
                  x_step = abs(x_step)
                  p2ref = x_step(1)
                  f2ref = x_step(2)
                  prev_solution = x_step
#ifdef CRESP_VERBOSED
                  call msg_success("step", x_step, x_vec)
#endif /* CRESP_VERBOSED */
                  return
               endif
            enddo
         enddo
      endif

   end subroutine seek_solution_step
!----------------------------------------------------------------------------------------------------
   subroutine fill_q_grid(i_incr)

      use initcrspectrum, only: p_fix_ratio, arr_dim_q

      implicit none

      integer(kind=4), intent(in) :: i_incr
      integer(kind=4)             :: i, j, i_beg, i_end
      real                        :: x, prev_solution
      logical                     :: exit_code

      i_beg = 1
      i_end = arr_dim_q

      p_ratio_4_q = p_fix_ratio

      prev_solution = q_space(int(helper_arr_dim/2))

      do i = i_beg, i_end, i_incr
         exit_code = .true.
         alpha = alpha_tab_q(i)
#ifdef CRESP_VERBOSED
         write(*,"(A25,1I4,A9,I4,A10,1E16.9)",advance="no") "Now solving (q_grid) no.",i,", sized ",arr_dim_q, ", (alpha): ",alpha
#endif /* CRESP_VERBOSED */
         x = prev_solution
         call NR_algorithm_1D(x, exit_code)
         if (exit_code) then
            do j = 1, helper_arr_dim
               if (exit_code) then
                  x = q_space(j)
                  call NR_algorithm_1D(x,exit_code)
                  if (.not. exit_code) then
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

      real,    intent(inout) :: x
      logical, intent(out)   :: exit_code

      if (abs(x) >= q_big) then
         x = sign(one, x) * q_big
         exit_code = .true.
         return
      endif

   end subroutine q_control
 !----------------------------------------------------------------------------------------------------
   real function alpha_to_q(x) ! this one (as of now) is only usable with fixed p_ratio_4_q bins (middle ones)

      use constants,      only: one, three
      use initcrspectrum, only: q_eps

      implicit none

      real, intent(in) :: x
      real             :: q_in3, q_in4

      q_in3 = three - x
      q_in4 = one + q_in3
      if (abs(q_in3) < q_eps) then
         alpha_to_q = (p_ratio_4_q**q_in4 - one)/log(p_ratio_4_q)
      else if (abs(q_in4) < q_eps) then
         alpha_to_q = q_in3 * log(p_ratio_4_q)/(p_ratio_4_q**q_in3 - one)
      else
         alpha_to_q = q_in3/q_in4 * (p_ratio_4_q**q_in4 - one)/(p_ratio_4_q**q_in3 - one)
      endif
      alpha_to_q = alpha_to_q - alpha

   end function alpha_to_q
!----------------------------------------------------------------------------------------------------
   subroutine prepare_indices(ind_incr, ind_beg, ind_end)

      use initcrspectrum, only: arr_dim

      implicit none

      integer(kind=4), intent(in)  :: ind_incr
      integer(kind=4), intent(out) :: ind_beg, ind_end

      if (ind_incr == 1) then
         ind_beg = 1 ; ind_end = arr_dim
      else if (ind_incr == -1) then
         ind_beg = arr_dim ; ind_end = 1
      endif

   end subroutine prepare_indices
 !----------------------------------------------------------------------------------------------------
   real function q_ratios(f_ratio, p_ratio)

      implicit none

      real, intent(in) :: f_ratio, p_ratio

      q_ratios = -log10(f_ratio) / log10(p_ratio)

   end function q_ratios
!---------------------------------------------------------------------------------------------------
! Function estimating values of jacobian via finite difference method
!---------------------------------------------------------------------------------------------------
   function jac_fin_diff(x) ! jacobian via finite difference method

      use constants, only: half

      implicit none

      real, dimension(ndim)            :: x, xp, xm
      real, dimension(size(x),size(x)) :: jac_fin_diff
      real, dimension(size(x))         :: dx
      real, parameter                  :: dx_par = 1.0e-3, dx_min = epsilon(dx_par)
      integer(kind=2)                  :: j

      do j = 1, ndim
         dx(j) = max(x(j), dx_min )          ! assure dx > zero
         dx(j) = min(dx(j)*dx_par, dx_par)   ! the value of dx is scaled not to go over value of x
         if (x(j) == dx(j)) dx(j) = half * dx(j)
         xp = x ; xm = x
         xp(j) = x(j) - dx(j) ;  xm(j) = x(j) + dx(j)
         jac_fin_diff(:,j)  = half*( selected_function_2D(xp) - selected_function_2D(xm)) / dx(j)
      enddo

   end function jac_fin_diff
!----------------------------------------------------------------------------------------------------
   real function determinant_2d_real(matrix_2d_real)

      implicit none

      real, dimension(2,2), intent(in) :: matrix_2d_real

      determinant_2d_real = matrix_2d_real(1,1) * matrix_2d_real (2,2) - ( matrix_2d_real(2,1) * matrix_2d_real(1,2) )

   end function determinant_2d_real
!----------------------------------------------------------------------------------------------------
   function get_cofactor_matrix_2d_real(matrix_2d_real)

      implicit none

      real, dimension(ndim,ndim) :: matrix_2d_real
      real, dimension(ndim,ndim) :: get_cofactor_matrix_2d_real
      integer(kind=1)            :: i, j

      do i = 1, ndim
         do j = 1, ndim
            get_cofactor_matrix_2d_real(i,j) = ( (-1)**(i+j) * matrix_2d_real(ndim+1-i, ndim+1-j))
         enddo
      enddo

   end function get_cofactor_matrix_2d_real
!----------------------------------------------------------------------------------------------------
   function invert_2d_matrix(matrix, determinant)

      use constants, only: one

      implicit none

      real, dimension(ndim,ndim), intent(in) :: matrix
      real                                   :: determinant
      real, dimension(ndim,ndim)             :: invert_2d_matrix

      invert_2d_matrix = (one / determinant) * transpose( get_cofactor_matrix_2d_real(matrix) )

   end function invert_2d_matrix
!----------------------------------------------------------------------------------------------------
   function fvec_up(x)

      use constants, only: three

      implicit none

      real, dimension(ndim), intent(inout) :: x
      real, dimension(ndim)                :: fvec_up
      real                                 :: q_in3

      x = abs(x)
      q_in3      = three - q_ratios(x(2), x(1))
      fvec_up(1) = encp_func_2_zero(HI, x(1),                  alpha, q_in3)
      fvec_up(2) =    n_func_2_zero(    x(1), x(2)*x(1)**three, n_in, q_in3)

   end function fvec_up

!----------------------------------------------------------------------------------------------------
   function fvec_lo(x)

      use constants, only: one, three

      implicit none

      real, dimension(ndim), intent(inout) :: x
      real, dimension(ndim)                :: fvec_lo
      real                                 :: q_in3

      x = abs(x)
      q_in3      = three - q_ratios(x(2), x(1))
      fvec_lo(1) = encp_func_2_zero(LO, x(1),     alpha, q_in3)
      fvec_lo(2) =    n_func_2_zero(    x(1), one, n_in, q_in3)

   end function fvec_lo

!---------------------------------------------------------------------------------------------------
   real function encp_func_2_zero(side, p_ratio, alpha_cnst, q_in3) ! from eqn. 29

      use constants,      only: one
      use initcrspectrum, only: eps

      implicit none

      integer(kind=4), intent(in) :: side
      real,            intent(in) :: p_ratio, alpha_cnst, q_in3
      real                        :: q_in4

      q_in4 = one + q_in3
      if (abs(q_in3) < eps) then
         encp_func_2_zero = (p_ratio**q_in4)/(q_in4*log(p_ratio))  ! if q = 3
      else if (abs(q_in4) < eps) then
         encp_func_2_zero = q_in3 * log(p_ratio)/(p_ratio**q_in3 - one)
      else
         encp_func_2_zero = q_in3/q_in4*(p_ratio**q_in4 - one) / (p_ratio**q_in3 - one)
      endif
      if (side == LO) encp_func_2_zero = encp_func_2_zero / p_ratio
      encp_func_2_zero = encp_func_2_zero - alpha_cnst

   end function encp_func_2_zero

!----------------------------------------------------------------------------------------------------
   real function n_func_2_zero(p_ratio, fp_cmplx, n_cnst, q_in3) ! from eqn. 9

      use constants,       only: one
      use cresp_variables, only: clight_cresp
      use initcrspectrum,  only: e_small, eps

      implicit none

      real, intent(in) :: p_ratio, fp_cmplx, n_cnst, q_in3

      n_func_2_zero = e_small / (clight_cresp * fp_cmplx)
      if (abs(q_in3) < eps) then
         n_func_2_zero = n_func_2_zero * log(p_ratio)
      else
         n_func_2_zero = n_func_2_zero * (p_ratio**q_in3 - one)/q_in3
      endif
      n_func_2_zero = n_func_2_zero - n_cnst

   end function n_func_2_zero

!----------------------------------------------------------------------------------------------------
   real function bl_interpol(y11, y12, y21, y22, t, u) ! for details see paragraph "Bilinear interpolation" in Numerical Recipes for F77, page 117, eqn. 3.6.5

      use constants, only: one

      implicit none

      real, intent(in) :: y11, y12, y21, y22, t, u ! y** - tabularized values of interpolated function, t, u - coefficients

      bl_interpol = (one - t)*(one - u) * y11 + t*(one - u)*y12 + (one - t)*u*y21 + t*u*y22

   end function bl_interpol
!----------------------------------------------------------------------------------------------------
   real function bl_in_tu(val_left, val_mid, val_right) ! for details see paragraph "Bilinear interpolation" in Numerical Recipes for F77, page 117, eqn. 3.6.4

      implicit none

      real, intent(in) :: val_left, val_mid, val_right

      bl_in_tu = (val_mid - val_left) / (val_right - val_left)

   end function bl_in_tu
!----------------------------------------------------------------------------------------------------
   real function lin_interpol_1D(loc_1, loc_2, val)

      implicit none

      integer(kind=4), intent(in) :: loc_1, loc_2
      real,            intent(in) :: val

      lin_interpol_1D = p_n(loc_1) + (val - p_a(loc_1)) * ( p_n(loc_1) - p_n(loc_2) ) / (p_a(loc_1) - p_a(loc_2)) ! WARNING - uses p_a and p_n, that are usually used to point alpha and n arrays.

   end function lin_interpol_1D
!----------------------------------------------------------------------------------------------------
   real function lin_extrapol_1D(fun, arg, arg_out)

      implicit none

      real, dimension(1:2), intent(in) :: fun, arg
      real,                 intent(in) :: arg_out

      lin_extrapol_1D = fun(1) + (fun(2) - fun(1)) * (arg_out - arg(1))/(arg(2)-arg(1))

   end function lin_extrapol_1D
!----------------------------------------------------------------------------------------------------
  function intpol_pf_from_NR_grids(a_val, n_val, successful) ! for details see paragraph "Bilinear interpolation" in Numerical Recipes for F77, page 117, eqn. 3.6.4

      implicit none

      real,    intent(in)             :: a_val, n_val  ! ratios arrays (p,f: lo and up), for which solutions have been obtained. loc_no_ip (changed to l1) - in case when interpolation is not possible,
      logical, intent(out)            :: successful
      real, dimension(2)              :: intpol_pf_from_NR_grids ! indexes with best match and having solutions are chosen.
      real                            :: blin_a, blin_n
      integer(kind=4), dimension(1:2) :: l1, l2 ! indexes that points where alpha_tab_ and up nad n_tab_ and up are closest in value to a_val and n_val - indexes point to

#ifdef CRESP_VERBOSED
      write (*,"(A30,A2,A4)",advance="no") "Determining indices for case: ", bound_name(current_bound), "... "
#endif /* CRESP_VERBOSED */
      call determine_loc(a_val, n_val, l1, successful)
      l2 = l1 + 1

#ifdef CRESP_VERBOSED
      if (successful) write(*,"(A19, 2I8, A3, 2I8)") "Obtained indices:", l1, " | ", l2
      call save_loc(current_bound, l1, l2)
#endif /* CRESP_VERBOSED */

      if (successful) then
         blin_a = bl_in_tu(p_a(l1(1)), a_val, p_a(l2(1)))
         blin_n = bl_in_tu(p_n(l1(2)), n_val, p_n(l2(2)))
         intpol_pf_from_NR_grids(1) = bl_interpol(p_p(l1(1),l1(2)), p_p(l1(1),l2(2)), p_p(l2(1),l1(2)), p_p(l2(1),l2(2)), blin_a, blin_n)
         intpol_pf_from_NR_grids(2) = bl_interpol(p_f(l1(1),l1(2)), p_f(l1(1),l2(2)), p_f(l2(1),l1(2)), p_f(l2(1),l2(2)), blin_a, blin_n)
      else ! interpolation won't work in this case, choosing closest values that have solutions.
         intpol_pf_from_NR_grids(1) = p_p(l1(1), l1(2))
         intpol_pf_from_NR_grids(2) = p_f(l1(1), l1(2))
      endif

   end function intpol_pf_from_NR_grids
!----------------------------------------------------------------------------------------------------
   subroutine determine_loc(a_val, n_val, loc1, successful)

      use constants,      only: zero
      use initcrspectrum, only: arr_dim

      implicit none

      real,                            intent(in)  :: a_val, n_val
      integer(kind=4), dimension(1:2), intent(out) :: loc1
      logical,                         intent(out) :: successful
      logical                                      :: hit_zero

      hit_zero  = .false.
      loc1(1) = inverse_f_to_ind(a_val, p_a(1), p_a(arr_dim), arr_dim)
      loc1(2) = inverse_f_to_ind(n_val, p_n(1), p_n(arr_dim), arr_dim)

      if ((minval(loc1) >= 1 .and. maxval(loc1) <= arr_dim-1)) then ! only need to test loc1
         if (p_p(loc1(1), loc1(2)) > zero) then
            successful = .true.
            return        ! normal exit
         else
            hit_zero  = .true.
         endif
      endif
      successful = .false.  ! namely if ((minval(loc1) <= 0 .or. maxval(loc1) >= arr_dim))

      loc1(1) = max(1, min(loc1(1), arr_dim))   ! Here we either give algorithm closest nonzero value relative to a row that was in the proper range
      loc1(2) = max(1, min(loc1(2), arr_dim))   ! or we just feed the algorithm ANY nonzero initial vector that will prevent it from crashing.

      if (loc1(1) == arr_dim .or. hit_zero) call nearest_solution(p_p(:,loc1(2)), loc1(1),        1,       loc1(1), hit_zero)
      if (loc1(1) <= 1       .or. hit_zero) call nearest_solution(p_p(:,loc1(2)), max(1,loc1(1)), arr_dim, loc1(1), hit_zero)
      if (loc1(2) == arr_dim .or. hit_zero) call nearest_solution(p_p(loc1(1),:), loc1(2),        1,       loc1(2), hit_zero)
      if (loc1(2) <= 1       .or. hit_zero) call nearest_solution(p_p(loc1(1),:), max(1,loc1(2)), arr_dim, loc1(2), hit_zero)

   end subroutine determine_loc
!----------------------------------------------------------------------------------------------------
   subroutine nearest_solution(arr_lin, i_beg, i_end, i_solution, hit_zero)

      use constants, only: zero

      implicit none

      real, dimension(:), intent(in)    :: arr_lin
      integer(kind=4),    intent(in)    :: i_beg, i_end
      integer(kind=4),    intent(out)   :: i_solution
      logical,            intent(inout) :: hit_zero
      integer(kind=4)                   :: i, i_incr

      i_solution = i_beg
      i_incr = sign(1, i_end - i_beg)

      do  i = i_beg, i_end, i_incr
         if (arr_lin(i) > zero) then
            i_solution = i ! next one is i_solution + i_incr
            hit_zero = .false.
            return
         endif
      enddo

   end subroutine nearest_solution
!----------------------------------------------------------------------------------------------------
   real function compute_q(alpha_in, exit_code, outer_p_ratio)

      use constants,      only: zero, one, I_ZERO, I_ONE
      use initcrspectrum, only: NR_refine_solution_q, q_big, p_fix_ratio, arr_dim_q

      implicit none

      real,           intent(inout) :: alpha_in
      logical,        intent(inout) :: exit_code ! value should be .true. at input
      real, optional, intent(in)    :: outer_p_ratio
      integer(kind=4)               :: loc_1, loc_2

      p_a => alpha_tab_q
      p_n => q_grid

      compute_q = zero
      if (present(outer_p_ratio)) then
         p_ratio_4_q = outer_p_ratio
      else
         p_ratio_4_q = p_fix_ratio
      endif

      loc_1 = inverse_f_to_ind(alpha_in, alpha_tab_q(1), alpha_tab_q(arr_dim_q), arr_dim_q)

      if ((loc_1 <= I_ZERO) .or. (loc_1 >= arr_dim_q)) then
         if (loc_1 <= I_ZERO) then
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

      if (abs(compute_q) > q_big) compute_q = sign(one, compute_q) * q_big

   end function compute_q
!----------------------------------------------------------------------------------------------------
!    subroutine save_NR_guess_grid(NR_guess_grid, var_name, bc)
!
!       use cresp_variables, only: clight_cresp
!       use initcrspectrum,  only: e_small, q_big, max_p_ratio
!
!       implicit none
!
!       real, dimension(:,:), intent(in) :: NR_guess_grid
!       character(len=*),     intent(in) :: var_name
!       integer(kind=4),      intent(in) :: bc
!       integer(kind=4)                  :: j
!       character(len=flen)              :: f_name
!
!       f_name = var_name // bound_name(bc) // extension
!       open(31, file=f_name, status="unknown", position="rewind")
!          write(31,"(A56,A2,A110)") "This is a storage file for NR init grid, boundary case: ", bound_name(bc), &
!             &    " Saved below: e_small, size(NR_guess_grid,dim=1), size(NR_guess_grid,dim=2), max_p_ratio, q_big, clight_cresp. Do not remove content from this file"
!          write(31, "(1E15.8, 2I10,10E22.15)") e_small, size(NR_guess_grid,dim=1), size(NR_guess_grid, dim=2), max_p_ratio, q_big, clight_cresp ! TODO: remove max_p_ratio, swap cols, rows with just arr_dim
!          write(31, "(A1)") " "                            ! Blank line for
!
!          do j=1, size(NR_guess_grid,dim=2)
!             write(31, "(*(E24.15E3))") NR_guess_grid(:,j)  ! WARNING - MIGHT NEED EXPLICIT ELEMENT COUNT IN LINE IN OLDER COMPILERS
!          enddo
!          close(31)
!
!    end subroutine save_NR_guess_grid

   subroutine save_NR_smap(NR_smap, hdr, vname, bc)

      implicit none

      integer(kind=4),      intent(in) :: bc
      integer(kind=4), parameter       :: flun = 31
      character(len=flen)              :: fname
      type(map_header),     intent(in) :: hdr
      integer(kind=4)                  :: j
      real, dimension(:,:), intent(in) :: NR_smap
      character(len=*),     intent(in) :: vname

      fname = vname // bound_name(bc) // extension
      open(flun, file=fname, status="unknown", position="rewind")
         write(flun,"(A56,A2,A26)") "This is a storage file for NR init grid, boundary case: ", bound_name(bc), &
            & ". Do not append this file."

         write(flun, "(1E15.8, 2I10,10E22.15)") hdr%s_es, hdr%s_dim1, hdr%s_dim2, hdr%s_pr, hdr%s_qbig, hdr%s_c, hdr%s_amin, hdr%s_amax, hdr%s_nmin, hdr%s_nmax
         write(flun, "(A1)") " "                             ! Blank line
         do j=1, size(NR_smap,dim=2)
            write(flun, "(*(E24.15E3))") NR_smap(:,j)  ! WARNING - MIGHT NEED EXPLICIT ELEMENT COUNT IN LINE IN OLDER COMPILERS
         enddo
         close(flun)

   end subroutine save_NR_smap
!----------------------------------------------------------------------------------------------------
#ifdef CRESP_VERBOSED
   subroutine save_loc(bound_case, loc1, loc2)

      implicit none

      integer(kind=4),                 intent(in) :: bound_case
      integer(kind=4), dimension(1:2), intent(in) :: loc1, loc2
      integer(kind=4), parameter                  :: fnlen = 10, flun = 32
      character(len=fnlen)                        :: f_name

      f_name = "loc_"//bound_name(bound_case)//extension
      open(flun, file=f_name, status="unknown", position="append")
      write (flun,"(2I5)") loc1(1), loc1(2)
      write (flun,"(2I5)") loc2(1), loc2(2)
      close(flun)

   end subroutine save_loc
#endif /* CRESP_VERBOSED */
!----------------------------------------------------------------------------------------------------
   subroutine read_NR_smap(NR_smap, vname, bc, exit_code)

      use func,            only: operator(.equals.)
      use initcrspectrum,  only: e_small, q_big, max_p_ratio

      implicit none

      real, dimension(:,:), intent(inout) :: NR_smap
      character(len=*),     intent(in)    :: vname
      integer(kind=4),      intent(in)    :: bc
      logical,              intent(out)   :: exit_code
      integer(kind=4)                     :: j, rstat = 0, flun = 31
      character(len=flen)                 :: fname

      fname = vname // bound_name(bc) // extension
      open(flun, file=fname, status="old", position="rewind", IOSTAT=rstat)
      if (rstat > 0) then
         write(*,"(A8,I4,A8,2A20)") "IOSTAT:", rstat, ": file ", vname//bound_name(bc)//extension," does not exist!"
         exit_code = .true.
         return
      else
         read(flun, *) ! Skipping comment line
         read(flun, *) ! Skipping header
         read(flun, *) ! Skipping blank line
         do j=1, size(NR_smap, dim=2)
            read(flun, "(*(E24.15E3))", IOSTAT=rstat) NR_smap(:,j)  ! WARNING - THIS MIGHT NEED EXPLICIT INDICATION OF ELEMENTS COUNT IN LINE IN OLDER COMPILERS
         enddo
         exit_code = .false.
      endif
      if (rstat > 0) exit_code = .true.
      close(flun)

   end subroutine read_NR_smap

   subroutine read_NR_smap_header(var_name, hdr, exit_code)

      use dataio_pub,   only: msg, printinfo, warn
      use constants,    only: fmt_len

      implicit none

      logical                          :: exit_code
      integer(kind=4), parameter       :: flun = 31
      character(len=flen)              :: f_name
      type(map_header), intent(inout)  :: hdr
      character(len=fmt_len)           :: fmt
      character(len=*), intent(in)     :: var_name
      integer(kind=4)                  :: fstat, rstat

      fstat = 0
      rstat = 0
      f_name = var_name // extension
      fmt = "(1E15.8,2I10,10E22.15)"

      open(flun, file=f_name, status="old", position="rewind", IOSTAT=fstat)

      if (fstat > 0) then
         write(msg,"(A8,I4,A8,2A20)") "IOSTAT:", fstat, ": file ", f_name, " does not exist!"; call warn(msg)
         exit_code = .true.
         return
      endif

      read(flun, fmt, IOSTAT=rstat) hdr%s_es, hdr%s_dim1, hdr%s_dim2, hdr%s_pr, hdr%s_qbig, hdr%s_c, hdr%s_amin, hdr%s_amax, hdr%s_nmin, hdr%s_nmax
      if (rstat > 0 ) then  ! should work for older files using the same format
         read(flun, fmt, IOSTAT=rstat) hdr%s_es, hdr%s_dim1, hdr%s_dim2, hdr%s_pr, hdr%s_qbig, hdr%s_c
         hdr%s_amin = 0.
         hdr%s_amax = 0.
         hdr%s_nmin = 0.
         hdr%s_nmax = 0.
      endif

      exit_code = .false.
      close(flun)

   end subroutine read_NR_smap_header

   subroutine check_NR_smap_header(hdr, hdr_std, hdr_equal)

      use constants, only: zero
      use dataio_pub,only: msg, printinfo, warn
      use func,      only: operator(.equals.)

      implicit none

      type(map_header), intent(in)  :: hdr, hdr_std
      logical                       :: hdr_equal

      hdr_equal = .true.

      hdr_equal = hdr_equal .and. (hdr%s_es .equals. hdr_std%s_es)
      hdr_equal = hdr_equal .and. (hdr%s_dim1 .eq.   hdr_std%s_dim1)
      hdr_equal = hdr_equal .and. (hdr%s_dim2 .eq.   hdr_std%s_dim2)
      hdr_equal = hdr_equal .and. (hdr%s_pr .equals. hdr_std%s_pr)
      hdr_equal = hdr_equal .and. (hdr%s_c  .equals. hdr_std%s_c)

!  WARNING allowing to read old solution maps; without saved a_tab and n_tab limits
      hdr_equal = hdr_equal .and. ((hdr%s_amin .equals. hdr_std%s_amin) .or. (hdr%s_amin .equals. zero))
      hdr_equal = hdr_equal .and. ((hdr%s_amax .equals. hdr_std%s_amax) .or. (hdr%s_amax .equals. zero))
      hdr_equal = hdr_equal .and. ((hdr%s_nmin .equals. hdr_std%s_nmin) .or. (hdr%s_nmin .equals. zero))
      hdr_equal = hdr_equal .and. ((hdr%s_nmax .equals. hdr_std%s_nmax) .or. (hdr%s_nmax .equals. zero))

      if (.not. hdr_equal) then
         write(msg,"(A117)") "[cresp_NR_method:check_NR_smap_header] Headers differ (provided in ratios files vs. values resulting from parameters)"
         call warn(msg)
      else
         write(msg,"(A115)") "[cresp_NR_method:check_NR_smap_header] Headers match (provided in ratios files vs. values resulting from parameters)"
         call printinfo(msg)
      endif

   end subroutine check_NR_smap_header
!----------------------------------------------------------------------------------------------------
   real function ind_to_flog(ind, min_in, max_in, length)

      use constants, only: I_ONE, ten

      implicit none

      real,            intent(in) :: min_in, max_in
      integer(kind=4), intent(in) :: ind, length

      ind_to_flog = min_in * ten**(((log10(max_in/min_in))/real(length-I_ONE))*real(ind-I_ONE))

   end function ind_to_flog
!----------------------------------------------------------------------------------------------------
   integer(kind=4) function inverse_f_to_ind(val, min_in, max_in, length) ! returns lower index for a given value, will need limiter

      use constants, only: I_ONE

      implicit none

      real,            intent(in) :: val, min_in, max_in
      integer(kind=4), intent(in) :: length

      inverse_f_to_ind = int((log10(val/min_in)/log10(max_in/min_in)) * (length - I_ONE )) + I_ONE

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

   subroutine add_dot(is_finishing)

      use constants, only: stdout
      use mpisetup,  only: master

      implicit none

      logical   ::  is_finishing

      if (master) then
         if (.not. is_finishing) then
            write(stdout,'(a)',advance='no')"."
         else
            write(stdout,'(a)',advance='yes')"."
         endif
      endif

   end subroutine add_dot
end module cresp_NR_method
