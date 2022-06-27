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
!! \brief This module contains Newton-Raphson 1,2-dim algorithms, functions for use with these and used by CRESP.
!<

module cresp_NR_method
! pulled by CRESP

   use cresp_helpers, only: map_header

   implicit none

   private
   public :: alpha, assoc_pointers, cresp_initialize_guess_grids, compute_q, intpol_pf_from_NR_grids, n_in, NR_algorithm, q_ratios, cresp_write_smaps_to_hdf, cresp_read_smaps_from_hdf, deallocate_all_smaps

   integer, parameter                        :: ndim = 2
   real, allocatable, dimension(:)           :: p_space, q_space
   real                                      :: alpha, p_ratio_4_q, n_in
   real, allocatable, dimension(:)           :: alpha_tab_q, q_grid
   real, allocatable, dimension(:,:)         :: alpha_tab, n_tab
   integer(kind=4)                           :: helper_arr_dim
   real, pointer, dimension(:,:)             :: p_p => null(), p_f => null() ! pointers for p_ratios_(lo,up) and f_ratios_(lo,up)
#ifdef CRESP_VERBOSED
   integer(kind=4)                           :: sought_by
   integer(kind=4), parameter                :: SLV = 1, RFN = 2
#endif /* CRESP_VERBOSED */
   logical, save                             :: got_smaps_from_restart = .false.
   type(map_header), dimension(2)            :: hdr_res

   abstract interface
      function function_pointer_2D(z)
         real, dimension(2), intent(inout) :: z
         real, dimension(2)                :: function_pointer_2D
      end function function_pointer_2D
   end interface

   procedure (function_pointer_2D), pointer :: selected_function_2D => null()

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

      fun_vec_value = selected_function_2D(x) - [alpha, n_in]
      if (maxval(abs(fun_vec_value)) < 0.01 * err_f) then ! in case when f converges at initialization
         exit_code = .false.
         return
      endif

      do i = 1, NR_iter_limit
         if (maxval(abs(fun_vec_value)) < err_f) then    ! For convergence via value of f
            exit_code = .false.
#ifdef CRESP_VERBOSED
            write(*,"(A47,I4,A12)",advance="no") "Convergence via value of fun_vec_value after ",i, " iterations."!, x, fun_vec_value  ! QA_WARN debug
#endif /* CRESP_VERBOSED */
            return
         endif

         fun_vec_value = selected_function_2D(x) - [alpha, n_in]
         fun_vec_jac = jac_fin_diff(x)                    ! function vector already explicitly provided to jac_fin_diff (finite difference method)

         det = determinant_2d_real(fun_vec_jac)           ! WARNING - algorithm is used for ndim = 2. For more dimensions LU or other methods should be implemented.
         if (abs(det) < eps_det) then              ! Countermeasure against determinant = zero
            exit_code = .true.
            return
         endif
         fun_vec_inv_jac = invert_2d_matrix(fun_vec_jac,det)

         cor(1) = sum(fun_vec_inv_jac(1,:) * fun_vec_value(:))
         cor(2) = sum(fun_vec_inv_jac(2,:) * fun_vec_value(:))
         x = x + cor
         if (maxval(abs(cor)) < err_x) then                 ! For convergence via value of correction (cor) table.
#ifdef CRESP_VERBOSED
            write(*,"(A47,I4,A12)",advance="no") "Convergence via value of cor array     after ",i," iterations."  ! QA_WARN debug
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
         fun1D_val = alpha_to_q(x) - alpha
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

      use constants,      only: zero, I_FOUR, LO, HI
      use cresp_helpers,  only: map_header, hdr_io, p_ratios_lo, f_ratios_lo, p_ratios_up, f_ratios_up
      use cresp_io,       only: check_NR_smaps_headers, save_NR_smap
      use dataio_pub,     only: warn
      use initcrspectrum, only: arr_dim_a, force_init_NR, e_small_approx_init_cond, NR_allow_old_smaps, NR_smap_file
      use mpisetup,       only: master

      implicit none

      logical                        :: first_run = .true., solve_new_smaps
      logical, dimension(2)          :: hdr_match_res
      type(map_header), dimension(2) :: hdr

      call allocate_all_smap_arrays
      if (master .and. first_run) then
         helper_arr_dim = int(arr_dim_a/I_FOUR, kind=4)

         if (.not. allocated(p_space)) allocate(p_space(1:helper_arr_dim)) ! these will be deallocated once initialization is over
         if (.not. allocated(q_space)) allocate(q_space(1:helper_arr_dim)) ! these will be deallocated once initialization is over

         if (.not. got_smaps_from_restart) then
            p_ratios_up = zero ; f_ratios_up = zero
            p_ratios_lo = zero ; f_ratios_lo = zero
         endif

         call init_smap_array_values(hdr)
         ! now check if restart file was loaded
         solve_new_smaps = .not. got_smaps_from_restart
         if (got_smaps_from_restart) then
         ! if restart file loaded, check the header
               call check_NR_smaps_headers(hdr_res, hdr, hdr_match_res)
            solve_new_smaps = .not. all(hdr_match_res)
         endif

         if ( .not. got_smaps_from_restart .or. (hdr_match_res(LO) .eqv. .false. .or. hdr_match_res(HI) .eqv. .false.) ) then
            if (e_small_approx_init_cond == 1) then   ! implies solution maps will be needed at leas once. ! TODO make imports from initcrspectrum more readable
               call warn("[cresp_NR_method] Different smap parameters in restart or restart unavailable. Trying harder...")
               ! if header not in agreement, try to load backup, e.g., "NR_smap_file". If additionally user declared "NR_allow_old_smaps", read dat files
               call try_read_user_h5(NR_smap_file, hdr, solve_new_smaps)
               if (NR_allow_old_smaps) call try_read_old_smap_files(hdr, solve_new_smaps)
               if (solve_new_smaps .or. force_init_NR) call fill_refine_smap
            endif
         endif

         ! if "NR_allow_old_smaps" save to work dat files
         if (NR_allow_old_smaps) then
            call save_NR_smap(p_ratios_lo, hdr(LO), "pWratios_", LO)   ! save to work file
            call save_NR_smap(f_ratios_lo, hdr(LO), "fWratios_", LO)   ! save to work file
            call save_NR_smap(p_ratios_up, hdr(HI), "pWratios_", HI)   ! save to work file
            call save_NR_smap(f_ratios_up, hdr(HI), "fWratios_", HI)   ! save to work file
         endif

         ! if flag "NR_run_refine_pf", try to overwrite solution in "NR_smap_file". Die, if file has extension "res"
         ! finalize, broadcast, clean and exit.

         if (allocated(p_space)) deallocate(p_space) ! only needed at initialization
         if (allocated(q_space)) deallocate(q_space)

         first_run = .false.
      endif

      call cresp_NR_mpi_exchange(hdr)

      hdr_io = hdr   ! for consistency in hdr_io within all threads, though most important arrays were already broadcasted

   end subroutine cresp_initialize_guess_grids

!----------------------------------------------------------------------------------------------------
   subroutine init_smap_array_values(hdr_init)

      use constants,       only: zero, half, one, three, I_ONE, big, small, LO, HI
      use cresp_helpers,   only: map_header
      use cresp_variables, only: clight_cresp
      use dataio_pub,      only: die
      use initcrspectrum,  only: arr_dim_a, arr_dim_n, arr_dim_q, e_small, max_p_ratio, p_fix_ratio, q_big

      implicit none

      real               :: a_min_q = big, a_max_q = small , q_in3, pq_cmplx
      real, dimension(2) :: a_min   = big, a_max   = small, n_min = big, n_max = small
      integer(kind=4)    :: i, j, ilim = 0, qmaxiter = 100
      type(map_header), dimension(2), intent(inout) :: hdr_init


      q_space = zero
      do i = I_ONE, int(half*helper_arr_dim, kind=4)
         q_space(i) = ln_eval_array_val(i, q_big, real(0.05), int(1,kind=4), int(half*helper_arr_dim,kind=4)) ! BEWARE: magic number
      enddo

      do i = I_ONE, int(half*helper_arr_dim, kind=4)!, arr_dim_a
         q_space(int(half*helper_arr_dim,kind=4)+i) = -q_space(int(half*helper_arr_dim,kind=4)+1-i)
      enddo

! setting up a grids of ratios to be used as phase space for NR tabs, obtained later
      do i = 1, helper_arr_dim
         p_space(i) = max_p_ratio**(real(i)/real(helper_arr_dim)) ! ind_to_flog(i, 1.000001, max_p_ratio) ! max_p_ratio**(real(i)/real(arr_dim_a))
      enddo
      do i = 1, helper_arr_dim
         do j = 1, helper_arr_dim
            q_in3 = three - q_space(j)
            pq_cmplx = p_space(i)**q_in3

            a_min(LO) = min(a_min(LO), abs(encp_func_2_zero(LO, p_space(i),           q_in3)))
            n_min(LO) = min(n_min(LO), abs(   n_func_2_zero(    p_space(i), one,      q_in3)))
            a_min(HI) = min(a_min(HI), abs(encp_func_2_zero(HI, p_space(i),           q_in3)))
            n_min(HI) = min(n_min(HI), abs(   n_func_2_zero(    p_space(i), pq_cmplx, q_in3)))

            a_max(LO) = max(a_max(LO), abs(encp_func_2_zero(LO, p_space(i),           q_in3)))
            n_max(LO) = max(n_max(LO), abs(   n_func_2_zero(    p_space(i), one,      q_in3)))
            a_max(HI) = max(a_max(HI), abs(encp_func_2_zero(HI, p_space(i),           q_in3)))
            n_max(HI) = max(n_max(HI), abs(   n_func_2_zero(    p_space(i), pq_cmplx, q_in3)))
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

      do i = 1, arr_dim_a
         alpha_tab(LO, i) = ind_to_flog(i, a_min(LO), a_max(LO), arr_dim_a) ! a_min_lo * ten**((log10(a_max_lo/a_min_lo))/real(arr_dim_a-1)*real(i-1))
         alpha_tab(HI, i) = ind_to_flog(i, a_min(HI), a_max(HI), arr_dim_a) ! a_min_up * ten**((log10(a_max_up/a_min_up))/real(arr_dim_a-1)*real(i-1))
      enddo
      do i = 1, arr_dim_n
         n_tab(LO, i)     = ind_to_flog(i, n_min(LO), n_max(LO), arr_dim_n) ! n_min_lo * ten**((log10(n_max_lo/n_min_lo))/real(arr_dim_n-1)*real(i-1))
         n_tab(HI, i)     = ind_to_flog(i, n_min(HI), n_max(HI), arr_dim_n) ! n_min_up * ten**((log10(n_max_up/n_min_up))/real(arr_dim_n-1)*real(i-1))
      enddo


      q_grid      = q_big; q_grid(int(arr_dim_q/2):) = -q_big


      a_min_q = one  + epsilon(one)
      a_max_q = (one + epsilon(one)) * p_fix_ratio
      j = min(arr_dim_q - int(arr_dim_q/100, kind=4), arr_dim_q - I_ONE)               ! BEWARE: magic number

      do while ((q_grid(j) <= (-q_big) .and. (q_grid(arr_dim_q) <= (-q_big))) .and. (ilim .le. qmaxiter) )
         a_max_q = a_max_q - a_max_q*0.005                                             ! BEWARE: magic number
         do i = 1, arr_dim_q
            alpha_tab_q(i)  = ind_to_flog(i, a_min_q, a_max_q, arr_dim_q)
         enddo
         call fill_q_grid(i_incr=I_ONE) ! computing q_grid takes so little time, that saving the grid is not necessary.
         ilim = ilim + I_ONE
      enddo
      if (ilim .ge. qmaxiter) call die ("[cresp_NR_method:fill_guess_grids] Maximal iteration limit exceeded, q_grid might not have converged!")
#ifdef CRESP_VERBOSED
      do i = 1, arr_dim_q
         print "(A1,I3,A7,2F18.12)", "[ ", i,"] a : q ", q_grid(i), alpha_tab_q(i)
      enddo

      print *,"alpha_tab_lo(i),      alpha_tab_up(i),        n_tab_lo(i),        n_tab_up(i)  |       p_space(i),     q_space(i)"
      do i = 1, min(arr_dim_a, arr_dim_n)
         if (i <= helper_arr_dim) then
            print *,i,"|",  alpha_tab(LO,i), alpha_tab(HI,i), n_tab(LO,i), n_tab(HI,i), alpha_tab_q(i), "| i = ", &
                          min(i,helper_arr_dim), p_space(min(i,helper_arr_dim)), q_space(min(i,helper_arr_dim)), &
                          p_space(min(i,helper_arr_dim))**(-q_space(min(i,helper_arr_dim)))
         else
            print *,i,"|",  alpha_tab(LO,i), alpha_tab(HI,i), n_tab(LO,i), n_tab(HI,i), alpha_tab_q(i)
         endif
      enddo
      print *, "-----------"
#endif /* CRESP_VERBOSED */

      hdr_init(:)%s_es     = e_small
      hdr_init(:)%s_dim1   = arr_dim_a
      hdr_init(:)%s_dim2   = arr_dim_n
      hdr_init(:)%s_qbig   = q_big
      hdr_init(:)%s_pr     = max_p_ratio
      hdr_init(:)%s_c      = clight_cresp

      hdr_init(1)%s_amin = alpha_tab(LO,1)
      hdr_init(1)%s_amax = alpha_tab(LO,arr_dim_a)
      hdr_init(1)%s_nmin = n_tab(LO,1)
      hdr_init(1)%s_nmax = n_tab(LO,arr_dim_n)

      hdr_init(2)%s_amin = alpha_tab(HI,1)
      hdr_init(2)%s_amax = alpha_tab(HI,arr_dim_a)
      hdr_init(2)%s_nmin = n_tab(HI,1)
      hdr_init(2)%s_nmax = n_tab(HI,arr_dim_n)


   end subroutine init_smap_array_values

!----------------------------------------------------------------------------------------------------
   subroutine fill_refine_smap

      use constants,      only: LO, HI
      use initcrspectrum, only: NR_run_refine_pf

      implicit none

      integer(kind=4) :: i

      do i = LO, HI
         call assoc_pointers(i)
         call fill_boundary_grid(i)

         if (NR_run_refine_pf) call refine_all_directions(i)
      enddo

   end subroutine fill_refine_smap

!----------------------------------------------------------------------------------------------------
   subroutine try_read_old_smap_files(hdr_init, solve_new_smap)

      use constants,     only: LO, HI
      use cresp_helpers, only: extension, flen, map_header, p_ratios_lo, f_ratios_lo, p_ratios_up, f_ratios_up
      use cresp_io,      only: read_NR_smap, read_NR_smap_header, check_NR_smaps_headers

      implicit none

      character(len=flen-len(extension))         :: filename_read_lo, filename_read_up
      type(map_header), dimension(2), intent(in) :: hdr_init
      type(map_header), dimension(2)             :: hdr_read
      logical, dimension(2)                      :: read_error, hdr_match
      logical,                       intent(out) :: solve_new_smap

      call get_smap_filename("p_ratios_", LO, filename_read_lo)
      call read_NR_smap_header(filename_read_lo, hdr_read(LO), read_error(LO))
      call get_smap_filename("p_ratios_", HI, filename_read_up)
      call read_NR_smap_header(filename_read_up, hdr_read(HI), read_error(HI))

      solve_new_smap = read_error(LO) .and. read_error(HI)

      if (.not. (read_error(LO) .and. read_error(HI)) ) then

         call check_NR_smaps_headers(hdr_read, hdr_init, hdr_match)
         if (all(hdr_match)) then
            call read_NR_smap(p_ratios_lo, "p"//filename_read_lo(2:9), LO, read_error(LO))
            call read_NR_smap(f_ratios_lo, "f"//filename_read_lo(2:9), LO, read_error(HI))

            solve_new_smap = solve_new_smap .and. (read_error(LO) .and. read_error(HI))

            call read_NR_smap(p_ratios_up, "p"//filename_read_up(2:9), HI, read_error(LO))
            call read_NR_smap(f_ratios_up, "f"//filename_read_up(2:9), HI, read_error(HI))

            solve_new_smap = solve_new_smap .and. (read_error(LO) .and. read_error(HI))
         else
            solve_new_smap = .true.
         endif
      endif

   end subroutine try_read_old_smap_files
!----------------------------------------------------------------------------------------------------
   subroutine refine_all_directions(bound_case)

      use constants,     only: I_ONE
      use cresp_helpers, only: bound_name
      use dataio_pub,    only: die, msg, printinfo

      implicit none

      integer(kind=4), intent(in) :: bound_case

#ifdef CRESP_VERBOSED
      sought_by = RFN
#endif /* CRESP_VERBOSED */

      write(msg,'(3a)') "Running refine for:", bound_name(bound_case), " boundary"
      call printinfo(msg)
      if (.not. allocated(p_space) .or. .not. allocated(q_space)) call die("[cresp_NR_method:refine_all_directions] refine_grids called after array deallocation, stopping")

      call refine_ij(bound_case, p_p, p_f,  I_ONE, -I_ONE)
      call refine_ji(bound_case, p_p, p_f,  I_ONE, -I_ONE)
      call refine_ij(bound_case, p_p, p_f, -I_ONE, -I_ONE)
      call refine_ji(bound_case, p_p, p_f, -I_ONE, -I_ONE)
      call refine_ij(bound_case, p_p, p_f,  I_ONE,  I_ONE)
      call refine_ji(bound_case, p_p, p_f,  I_ONE,  I_ONE)
      call refine_ij(bound_case, p_p, p_f, -I_ONE,  I_ONE)
      call refine_ji(bound_case, p_p, p_f, -I_ONE,  I_ONE)

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

      use constants,     only: LO, HI
      use cresp_helpers, only: p_ratios_lo, f_ratios_lo, p_ratios_up, f_ratios_up

      implicit none

      integer(kind=4), intent(in) :: bound_case

      if (bound_case == LO) then
         p_p => p_ratios_lo
         p_f => f_ratios_lo
         selected_function_2D => fvec_lo
      endif
      if (bound_case == HI) then
         p_p => p_ratios_up
         p_f => f_ratios_up
         selected_function_2D => fvec_up
      endif

   end subroutine assoc_pointers

!----------------------------------------------------------------------------------------------------
! TODO FIXME to be parallelized
!
! The cost of this routine theoretically depends as O(arr_dim**4), which may be a severe limitation
! if we ever need finer arrays. Practically the exponent is a bit higher, around 4.5  to 5.0.
! Parallelization may help with doubling arr_dim but anything more than that require change
! of the algorithm to decrease exponent to arr_dim.

   subroutine fill_boundary_grid(bound_case)

      use constants,      only: zero, I_ONE, I_TWO
      use cresp_helpers,  only: bound_name
      use dataio_pub,     only: msg, printinfo
      use initcrspectrum, only: arr_dim_a, arr_dim_n, eps

      implicit none

      integer(kind=4), intent(in) :: bound_case
      real, dimension(:), pointer :: pfp, pff
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

      p_p = zero ; p_f = zero
      x_step = zero
      write(msg, "(A,A2,A,I3,A,I3)") "[cresp_NR_method:fill_boundary_grid] Solving solution maps for cutoff case (",bound_name(bound_case),"): DIM=",arr_dim_a, ' x ', arr_dim_n
      call printinfo(msg)

      do i = 1, arr_dim_a
         call add_dot( i .eq. arr_dim_a )
         new_line = .true.
         prev_solution = prev_solution_1 ! easier to find when not searching from the top
         do j = 1, arr_dim_n
            alpha = alpha_tab(bound_case, i)
            n_in  = n_tab(bound_case, j)
#ifdef CRESP_VERBOSED
            write(*,"(A14,A2,A2,2I4,A9,I4,I4,A1)",advance="no") "Now solving (",bound_name(bound_case),") ",i,j,", sized ",arr_dim_a, arr_dim_n," "  ! QA_WARN debug
#endif /* CRESP_VERBOSED */

            call seek_solution_prev(p_p(i,j), p_f(i,j), prev_solution, exit_code)

            if (.not. exit_code .and. new_line) then
               prev_solution_1 = prev_solution
               new_line = .false.
            endif

            if (exit_code) then
               jm = j - I_TWO
               if (check_dimm(jm, arr_dim_n)) then
                  pfp => p_p(i,jm:j)
                  pff => p_f(i,jm:j)
                  call step_extr(pfp, pff, n_tab(bound_case, jm:j), exit_code)
               endif
               if (j >= 2) then
                  jm = j - I_ONE
                  if (p_p(i,jm) > zero) call seek_solution_step(p_p(i,j), p_f(i,j), prev_solution, i, jm, exit_code)
               endif
            endif
            if (exit_code) then !still...
               do is = I_ONE, helper_arr_dim
                  do js = 1, helper_arr_dim  ! WARNING: here comes O(arr_dim**4)
                     x_vec(1) = p_space(is)
                     x_vec(2) = p_space(is)**(-q_space(js))
#ifdef CRESP_VERBOSED
                     x_in = x_vec
#endif /* CRESP_VERBOSED */
                     if (exit_code) then
                        call NR_algorithm(x_vec, exit_code)
                        if (.not. exit_code) then
                           p_p(i,j) = x_vec(1) ! i index - alpha, j index - n_in
                           p_f(i,j) = x_vec(2)
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
                  prev_solution(1) = p_p(i,j)
                  prev_solution(2) = p_f(i,j)
               endif
            endif
#ifdef CRESP_VERBOSED
            if (exit_code) print *,""
#endif /* CRESP_VERBOSED */
         enddo
      enddo

      p_p = abs(p_p)
      p_f = abs(p_f)

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
               x_step(2) = prev_sol(2) + (p_space(max(min(i_sol+i,helper_arr_dim),I_ONE))**(-q_space(max(min(j_sol+j, helper_arr_dim), I_ONE))) - prev_sol(2)) / real(nssstep - jj + 0.1)
               call NR_algorithm(x_step, exit_code)
               if (.not. exit_code) return
            endif
         enddo
      enddo

   end subroutine step_seek

!----------------------------------------------------------------------------------------------------

   subroutine refine_ji(co, ref_p, ref_f, i_incr, j_incr) ! ref_f and ref_p should already be partially filled with solutions

      use constants,      only: zero, I_TWO
      use initcrspectrum, only: arr_dim_a, arr_dim_n

      implicit none

      real, dimension(:,:), intent(inout) :: ref_p, ref_f
      integer(kind=4),      intent(in)    :: co, i_incr, j_incr
      integer(kind=4)                     :: i, j, i_beg, i_end, j_beg, j_end, i1m, i2m, i1p
      real, dimension(1:2)                :: prev_solution
      logical                             :: exit_code

      prev_solution(1) = p_space(1)              ! refine must be called before these are deallocated
      prev_solution(2) = p_space(1)**q_space(1)
      call prepare_indices(arr_dim_a, i_incr, i_beg, i_end)
      call prepare_indices(arr_dim_n, j_incr, j_beg, j_end)
      do j = j_beg, j_end, j_incr
         do i = i_beg, i_end, i_incr
            alpha = alpha_tab(co, i)
            n_in  =     n_tab(co, j)
            if (ref_p(i,j) > zero .and. ref_f(i,j) > zero) then
               prev_solution(1) = ref_p(i,j)
               prev_solution(2) = ref_f(i,j)
            else
               call seek_solution_prev(ref_p(i,j), ref_f(i,j), prev_solution, exit_code) ! works for most cases
               if (exit_code) then
                  i1m = i-i_incr ; i2m = i - I_TWO * i_incr ; i1p = i+i_incr
                  if (check_dimm(i2m, arr_dim_a)                                 ) call step_extr(ref_p(i2m:i  :i_incr,j), ref_f(i2m:i  :i_incr,j),         alpha_tab(co, i2m:i  :i_incr), exit_code)
                  if (check_dimm(i1m, arr_dim_a) .and. check_dimm(i1p, arr_dim_a)) call step_inpl(ref_p(i1m:i1p:i_incr,j), ref_f(i1m:i1p:i_incr,j), i_incr, alpha_tab(co, i1m:i1p:i_incr), exit_code)
               endif
            endif
         enddo
      enddo

   end subroutine refine_ji

!----------------------------------------------------------------------------------------------------

   subroutine refine_ij(co, ref_p, ref_f, i_incr, j_incr) ! ref_f and ref_p should already be partially filled with solutions

      use constants,      only: zero, I_TWO
      use initcrspectrum, only: arr_dim_a, arr_dim_n

      implicit none

      real, dimension(:,:), intent(inout) :: ref_p, ref_f
      integer(kind=4),      intent(in)    :: co, i_incr, j_incr
      integer(kind=4)                     :: i, j, i_beg, i_end, j_beg, j_end, j1m, j2m, j1p
      real, dimension(1:2)                :: prev_solution
      logical                             :: exit_code

      prev_solution(1) = p_space(1)              ! refine must be called before these are deallocated
      prev_solution(2) = p_space(1)**q_space(1)
      call prepare_indices(arr_dim_a, i_incr, i_beg, i_end)
      call prepare_indices(arr_dim_n, j_incr, j_beg, j_end)
      do i = i_beg, i_end, i_incr
         do j = j_beg, j_end, j_incr
            alpha = alpha_tab(co, i)
            n_in  =     n_tab(co, j)
            if (ref_p(i,j) > zero .and. ref_f(i,j) > zero) then
               prev_solution(1) = ref_p(i,j)
               prev_solution(2) = ref_f(i,j)
            else
               call seek_solution_prev(ref_p(i,j), ref_f(i,j), prev_solution, exit_code) ! works for most cases
               if (exit_code) then
                  j1m = j-j_incr ; j2m = j - I_TWO * j_incr ; j1p = j+j_incr
                  if (check_dimm(j2m, arr_dim_n)                                 ) call step_extr(ref_p(i,j2m:j  :j_incr), ref_f(i,j2m:j  :j_incr),         n_tab(co, j2m:j  :j_incr), exit_code)
                  if (check_dimm(j1m, arr_dim_n) .and. check_dimm(j1p, arr_dim_n)) call step_inpl(ref_p(i,j1m:j1p:j_incr), ref_f(i,j1m:j1p:j_incr), j_incr, n_tab(co, j1m:j1p:j_incr), exit_code)
               endif
            endif
         enddo
      enddo

   end subroutine refine_ij

!----------------------------------------------------------------------------------------------------

   subroutine prepare_indices(ind_max, ind_incr, ind_beg, ind_end)

      use constants,  only: INVALID
      use dataio_pub, only: die

      implicit none

      integer(kind=4), intent(in)  :: ind_max, ind_incr
      integer(kind=4), intent(out) :: ind_beg, ind_end

      if (ind_incr == 1) then
         ind_beg = 1 ; ind_end = ind_max
      else if (ind_incr == -1) then
         ind_beg = ind_max ; ind_end = 1
      else
         ind_beg = INVALID; ind_end = INVALID
         call die("[cresp_NR_method:prepare_indices] ind_incr /= +/-1")
      endif

   end subroutine prepare_indices

!----------------------------------------------------------------------------------------------------

   logical function check_dimm(ind, irange)

      implicit none

      integer(kind=4), intent(in) :: ind, irange

      check_dimm = (ind >= 1 .and. ind <= irange)

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

      ! alpha and n are set !

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

      write (*, "(A6,A13,2E16.9)",advance="no") sought(sought_by)," (alpha, n): ", alpha, n_in  ! QA_WARN debug
      write (*, "(A5,A4,A42, 2E19.10e3)",advance="no") " -> (",met_name,") solution obtained, (p_ratio, f_ratio) = ", x_out  ! QA_WARN debug
      write (*, "(A21, 2E17.10)",advance="no") ", provided input:", x_in ; print *,""  ! QA_WARN debug

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

      ! alpha and n are set !

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
         write(*,"(A25,1I4,A9,I4,A10,1E16.9)",advance="no") "Now solving (q_grid) no.",i,", sized ",arr_dim_q, ", (alpha): ",alpha  ! QA_WARN debug
#endif /* CRESP_VERBOSED */
         x = prev_solution
         call NR_algorithm_1D(x, exit_code)
         if (exit_code) then
            do j = 1, helper_arr_dim
               if (exit_code) then
                  x = q_space(j)
                  call NR_algorithm_1D(x, exit_code)
                  if (.not. exit_code) then
                     q_grid(i) = x
                     prev_solution = x
#ifdef CRESP_VERBOSED
                     write (*, "(A44, 2E22.15)",advance="no") " ->        solution obtained, q_grid = ", x  ! QA_WARN debug
#endif /* CRESP_VERBOSED */
                  endif
               endif
         enddo
         else
#ifdef CRESP_VERBOSED
            write (*, "(A44, 1E22.15)",advance="no") " -> (prev) solution obtained, q_grid = ", x  ! QA_WARN debug
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

   end function alpha_to_q
!----------------------------------------------------------------------------------------------------
   real function q_ratios(f_ratio, p_ratio)

      implicit none

      real, intent(in) :: f_ratio, p_ratio

      q_ratios = -log10(f_ratio) / log10(p_ratio)

   end function q_ratios

!>
!! \brief Function estimating values of Jacobian via finite difference method
!! \todo allow user to tune dx_par
!<
   function jac_fin_diff(x)

      use constants, only: half

      implicit none

      real, dimension(ndim)            :: x, xp, xm
      real, dimension(size(x),size(x)) :: jac_fin_diff
      real, dimension(size(x))         :: dx
      real, parameter                  :: dx_par = 1.0e-3
      integer(kind=2)                  :: j

      do j = 1, ndim
         dx(j) = min(x(j)*dx_par, dx_par) ! the value of dx is scaled not to go over value of x
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

      use constants, only: three, HI

      implicit none

      real, dimension(ndim), intent(inout) :: x
      real, dimension(ndim)                :: fvec_up
      real                                 :: q_in3

      x = abs(x)
      q_in3      = three - q_ratios(x(2), x(1))
      fvec_up(1) = encp_func_2_zero(HI, x(1),                   q_in3)
      fvec_up(2) =    n_func_2_zero(    x(1), x(2)*x(1)**three, q_in3)

   end function fvec_up

!----------------------------------------------------------------------------------------------------
   function fvec_lo(x)

      use constants, only: LO, one, three

      implicit none

      real, dimension(ndim), intent(inout) :: x
      real, dimension(ndim)                :: fvec_lo
      real                                 :: q_in3

      x = abs(x)
      q_in3      = three - q_ratios(x(2), x(1))
      fvec_lo(1) = encp_func_2_zero(LO, x(1),      q_in3)
      fvec_lo(2) =    n_func_2_zero(    x(1), one, q_in3)

   end function fvec_lo

!---------------------------------------------------------------------------------------------------
   real function encp_func_2_zero(side, p_ratio, q_in3) ! from eqn. 29

      use constants,      only: LO, one
      use initcrspectrum, only: eps

      implicit none

      integer(kind=4), intent(in) :: side
      real,            intent(in) :: p_ratio, q_in3
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

   end function encp_func_2_zero

!----------------------------------------------------------------------------------------------------
   real function n_func_2_zero(p_ratio, fp_cmplx, q_in3) ! from eqn. 9

      use constants,       only: one
      use cresp_variables, only: clight_cresp
      use initcrspectrum,  only: e_small, eps

      implicit none

      real, intent(in) :: p_ratio, fp_cmplx, q_in3

      n_func_2_zero = e_small / (clight_cresp * fp_cmplx)
      if (abs(q_in3) < eps) then
         n_func_2_zero = n_func_2_zero * log(p_ratio)
      else
         n_func_2_zero = n_func_2_zero * (p_ratio**q_in3 - one)/q_in3
      endif

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
   real function lin_interpol_1D(a1, a2, n1, n2, val)

      implicit none

      real, intent(in) :: a1, a2, n1, n2, val

      lin_interpol_1D = n1 + (val - a1) * ( n1 - n2 ) / (a1 - a2)

   end function lin_interpol_1D
!----------------------------------------------------------------------------------------------------
   real function lin_extrapol_1D(fun, arg, arg_out)

      implicit none

      real, dimension(1:2), intent(in) :: fun, arg
      real,                 intent(in) :: arg_out

      lin_extrapol_1D = fun(1) + (fun(2) - fun(1)) * (arg_out - arg(1))/(arg(2)-arg(1))

   end function lin_extrapol_1D
!----------------------------------------------------------------------------------------------------
   function intpol_pf_from_NR_grids(co, a_val, n_val, successful) ! for details see paragraph "Bilinear interpolation" in Numerical Recipes for F77, page 117, eqn. 3.6.4

      use constants,     only: I_ONE
#ifdef CRESP_VERBOSED
      use cresp_helpers, only: bound_name  ! QA_WARN debug
#endif /* CRESP_VERBOSED */

      implicit none

      integer(kind=4), intent(in)     :: co
      real,            intent(in)     :: a_val, n_val  ! ratios arrays (p,f: lo and up), for which solutions have been obtained. loc_no_ip (changed to l1) - in case when interpolation is not possible,
      logical,         intent(out)    :: successful
      real, dimension(2)              :: intpol_pf_from_NR_grids ! indexes with best match and having solutions are chosen.
      real                            :: blin_a, blin_n
      integer(kind=4), dimension(1:2) :: l1, l2 ! indexes that points where alpha_tab_ and up and n_tab_ and up are closest in value to a_val and n_val - indexes point to

#ifdef CRESP_VERBOSED
      write (*,"(A30,A2,A4)",advance="no") "Determining indices for case: ", bound_name(co), "... "  ! QA_WARN debug
#endif /* CRESP_VERBOSED */
      call determine_loc(co, a_val, n_val, l1, successful)
      l2 = l1 + I_ONE

#ifdef CRESP_VERBOSED
      if (successful) write(*,"(A19, 2I8, A3, 2I8)") "Obtained indices:", l1, " | ", l2  ! QA_WARN debug
      call save_loc(co, l1, l2)
#endif /* CRESP_VERBOSED */

      if (successful) then
         blin_a = bl_in_tu(alpha_tab(co, l1(1)), a_val, alpha_tab(co, l2(1)))
         blin_n = bl_in_tu(    n_tab(co, l1(2)), n_val,     n_tab(co, l2(2)))
         intpol_pf_from_NR_grids(1) = bl_interpol(p_p(l1(1),l1(2)), p_p(l1(1),l2(2)), p_p(l2(1),l1(2)), p_p(l2(1),l2(2)), blin_a, blin_n)
         intpol_pf_from_NR_grids(2) = bl_interpol(p_f(l1(1),l1(2)), p_f(l1(1),l2(2)), p_f(l2(1),l1(2)), p_f(l2(1),l2(2)), blin_a, blin_n)
      else ! interpolation won't work in this case, choosing closest values that have solutions.
         intpol_pf_from_NR_grids(1) = p_p(l1(1), l1(2))
         intpol_pf_from_NR_grids(2) = p_f(l1(1), l1(2))
      endif

   end function intpol_pf_from_NR_grids
!----------------------------------------------------------------------------------------------------
   subroutine determine_loc(co, a_val, n_val, loc1, successful)

      use constants,      only: zero, I_ONE
      use initcrspectrum, only: arr_dim_a, arr_dim_n

      implicit none

      integer(kind=4),                 intent(in)  :: co
      real,                            intent(in)  :: a_val, n_val
      integer(kind=4), dimension(1:2), intent(out) :: loc1
      logical,                         intent(out) :: successful
      logical                                      :: hit_zero

      hit_zero  = .false.
      loc1(1) = inverse_f_to_ind(a_val, alpha_tab(co, 1), alpha_tab(co, arr_dim_a), arr_dim_a)
      loc1(2) = inverse_f_to_ind(n_val,     n_tab(co, 1),     n_tab(co, arr_dim_n), arr_dim_n)

      if ((minval(loc1) >= 1 .and. maxval(loc1) <= arr_dim_a-1)) then ! only need to test loc1
         if (p_p(loc1(1), loc1(2)) > zero) then
            successful = .true.
            return        ! normal exit
         else
            hit_zero  = .true.
         endif
      endif
      successful = .false.  ! namely if ((minval(loc1) <= 0 .or. maxval(loc1) >= arr_dim_a))

      loc1(1) = max(I_ONE, min(loc1(1), arr_dim_a))   ! Here we either give algorithm closest nonzero value relative to a row that was in the proper range
      loc1(2) = max(I_ONE, min(loc1(2), arr_dim_n))   ! or we just feed the algorithm ANY nonzero initial vector that will prevent it from crashing.

      if (loc1(1) == arr_dim_a .or. hit_zero) call nearest_solution(p_p(:,loc1(2)), loc1(1),             I_ONE,     loc1(1), hit_zero)
      if (loc1(1) <= 1         .or. hit_zero) call nearest_solution(p_p(:,loc1(2)), max(I_ONE, loc1(1)), arr_dim_a, loc1(1), hit_zero)
      if (loc1(2) == arr_dim_n .or. hit_zero) call nearest_solution(p_p(loc1(1),:), loc1(2),             I_ONE,     loc1(2), hit_zero)
      if (loc1(2) <= 1         .or. hit_zero) call nearest_solution(p_p(loc1(1),:), max(I_ONE, loc1(2)), arr_dim_n, loc1(2), hit_zero)

   end subroutine determine_loc
!----------------------------------------------------------------------------------------------------
   subroutine nearest_solution(arr_lin, i_beg, i_end, i_solution, hit_zero)

      use constants, only: zero, I_ONE

      implicit none

      real, dimension(:), intent(in)    :: arr_lin
      integer(kind=4),    intent(in)    :: i_beg, i_end
      integer(kind=4),    intent(out)   :: i_solution
      logical,            intent(inout) :: hit_zero
      integer(kind=4)                   :: i, i_incr

      i_solution = i_beg
      i_incr = sign(I_ONE, i_end - i_beg)

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
         return                        ! < returns compute_q with exit_code = .true.
      endif

      loc_2 = loc_1 + I_ONE
      compute_q = lin_interpol_1D(alpha_tab_q(loc_1), alpha_tab_q(loc_2), q_grid(loc_1), q_grid(loc_2), alpha_in)

      if (NR_refine_solution_q) then
         alpha = alpha_in
         call q_control(compute_q,exit_code)
         call NR_algorithm_1D(compute_q, exit_code)
      endif

      if (abs(compute_q) > q_big) compute_q = sign(one, compute_q) * q_big

   end function compute_q
!----------------------------------------------------------------------------------------------------
#ifdef CRESP_VERBOSED
   subroutine save_loc(bound_case, loc1, loc2)

      use cresp_helpers, only: bound_name, extension

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

      inverse_f_to_ind = int((log10(val/min_in)/log10(max_in/min_in)) * (length - I_ONE ), kind=4) + I_ONE

   end function inverse_f_to_ind
!----------------------------------------------------------------------------------------------------
   subroutine cresp_NR_mpi_exchange(hdr_share)

      use constants,     only: LO, HI
      use cresp_helpers, only: p_ratios_lo, f_ratios_lo, p_ratios_up, f_ratios_up
      use mpisetup,      only: piernik_MPI_Bcast

      implicit none

      integer                        :: i
      type(map_header), dimension(2) :: hdr_share

      call piernik_MPI_Bcast(p_ratios_lo)
      call piernik_MPI_Bcast(p_ratios_up)
      call piernik_MPI_Bcast(f_ratios_lo)
      call piernik_MPI_Bcast(f_ratios_up)
      call piernik_MPI_Bcast(q_grid)
      call piernik_MPI_Bcast(n_tab)
      call piernik_MPI_Bcast(alpha_tab)
      call piernik_MPI_Bcast(alpha_tab_q)

      do i = LO, HI
         call piernik_MPI_Bcast(hdr_share(i)%s_dim1)
         call piernik_MPI_Bcast(hdr_share(i)%s_dim2)
         call piernik_MPI_Bcast(hdr_share(i)%s_es)
         call piernik_MPI_Bcast(hdr_share(i)%s_pr)
         call piernik_MPI_Bcast(hdr_share(i)%s_qbig)
         call piernik_MPI_Bcast(hdr_share(i)%s_c)
         call piernik_MPI_Bcast(hdr_share(i)%s_amin)
         call piernik_MPI_Bcast(hdr_share(i)%s_amax)
         call piernik_MPI_Bcast(hdr_share(i)%s_nmin)
         call piernik_MPI_Bcast(hdr_share(i)%s_nmax)
      enddo

   end subroutine cresp_NR_mpi_exchange
!----------------------------------------------------------------------------------------------------
!>
!! (Stub) If work (more recent) smap file exists in directory, the subroutine
!! returns its name, so that it will be used later for header check, etc.
!<
   subroutine get_smap_filename(var_name, bc, fname_no_ext)

      use cresp_helpers, only: bound_name, extension, flen
      use dataio_pub,    only: msg, printinfo

      implicit none

      character(len=*),                   intent(in)  :: var_name
      integer(kind=4),                    intent(in)  :: bc
      character(len=flen-len(extension)), intent(out) :: fname_no_ext
      character(len=flen-len(extension))              :: fname_no_ext_work
      logical                                         :: file_exists

      fname_no_ext_work = var_name // bound_name(bc)
      fname_no_ext_work(2:2) = "W"                                    ! WARNING magic numbers
      inquire(file=fname_no_ext_work // extension, exist=file_exists)

      if (file_exists) then
         fname_no_ext = fname_no_ext_work
         write(msg, "(A,A)") "[cresp_NR_method] More recent (work) solution map files found: {p,f}", fname_no_ext(2:)//extension
      else
         fname_no_ext = var_name // bound_name(bc)
         write(msg, "(A,A)") "[cresp_NR_method] Default solution map files will be used: {p,f}", fname_no_ext(2:)//extension
      endif

      call printinfo(msg)

   end subroutine get_smap_filename
!----------------------------------------------------------------------------------------------------
   subroutine allocate_all_smap_arrays

      use constants,      only: HI, I_ONE
      use diagnostics,    only: ma2d, my_allocate, my_allocate_with_index
      use initcrspectrum, only: arr_dim_a, arr_dim_n, arr_dim_q

      implicit none

      if (.not. allocated(alpha_tab_q)) call my_allocate_with_index(alpha_tab_q, arr_dim_q, I_ONE)
      if (.not. allocated(q_grid))      call my_allocate_with_index(q_grid,      arr_dim_q, I_ONE)

      ma2d = [HI, arr_dim_a]
      if (.not. allocated(alpha_tab))   call my_allocate(alpha_tab, ma2d)

      ma2d = [HI, arr_dim_n]
      if (.not. allocated(n_tab))       call my_allocate(n_tab, ma2d)

      call allocate_smaps(arr_dim_a, arr_dim_n)

   end subroutine allocate_all_smap_arrays
!----------------------------------------------------------------------------------------------------
   subroutine allocate_smaps(dim1, dim2)

      use diagnostics,   only: ma2d, my_allocate
      use cresp_helpers, only: p_ratios_lo, f_ratios_lo, p_ratios_up, f_ratios_up

      implicit none

      integer(kind=4), intent(in) :: dim1, dim2

      ma2d = [dim1, dim2]
      if (.not. allocated(p_ratios_lo)) call my_allocate(p_ratios_lo, ma2d )
      if (.not. allocated(f_ratios_lo)) call my_allocate(f_ratios_lo, ma2d )
      if (.not. allocated(p_ratios_up)) call my_allocate(p_ratios_up, ma2d )
      if (.not. allocated(f_ratios_up)) call my_allocate(f_ratios_up, ma2d )

   end subroutine allocate_smaps
!----------------------------------------------------------------------------------------------------
   subroutine deallocate_all_smaps

      use diagnostics, only: my_deallocate

      implicit none

      call my_deallocate(alpha_tab_q)
      call my_deallocate(q_grid)
      call my_deallocate(alpha_tab)
      call my_deallocate(n_tab)
      call deallocate_smaps

   end subroutine deallocate_all_smaps
!----------------------------------------------------------------------------------------------------
   subroutine deallocate_smaps

      use diagnostics,   only: my_deallocate
      use cresp_helpers, only: p_ratios_lo, f_ratios_lo, p_ratios_up, f_ratios_up

      implicit none

      call my_deallocate(p_ratios_lo)
      call my_deallocate(f_ratios_lo)
      call my_deallocate(p_ratios_up)
      call my_deallocate(f_ratios_up)

   end subroutine deallocate_smaps
!----------------------------------------------------------------------------------------------------
   subroutine cresp_write_smaps_to_hdf(file_id)

      use constants,     only: LO, HI
      use cresp_helpers, only: n_g_smaps, dset_attrs, p_ratios_lo, f_ratios_lo, p_ratios_up, f_ratios_up
      use cresp_io,      only: save_smap_to_open
      use hdf5,          only: HID_T

      implicit none

      integer(HID_T), intent(in) :: file_id

      call save_smap_to_open(file_id, n_g_smaps(LO), dset_attrs(1), p_ratios_lo)
      call save_smap_to_open(file_id, n_g_smaps(LO), dset_attrs(2), f_ratios_lo)
      call save_smap_to_open(file_id, n_g_smaps(HI), dset_attrs(1), p_ratios_up)
      call save_smap_to_open(file_id, n_g_smaps(HI), dset_attrs(2), f_ratios_up)

   end subroutine cresp_write_smaps_to_hdf
!----------------------------------------------------------------------------------------------------
   subroutine cresp_read_smaps_from_hdf(file_id)

      use constants,     only: LO, HI, I_ZERO
      use cresp_io,      only: read_real_arr2d_dset, read_smap_header_h5
      use cresp_helpers, only: map_header, dset_attrs, n_g_smaps, p_ratios_lo, f_ratios_lo, p_ratios_up, f_ratios_up
      use dataio_pub,    only: warn
      use hdf5,          only: HID_T

      implicit none

      integer(HID_T), intent(in) :: file_id

      call read_smap_header_h5(file_id, hdr_res)

      if (hdr_res(1)%s_dim1 .eq. I_ZERO .or. hdr_res(1)%s_dim2 .eq. I_ZERO) then
         call warn("[cresp_NR_method:cresp_read_smaps_from_hdf] Got solution map dimension = 0. Will solve for new maps.")
      else
         call deallocate_smaps ! TODO just in case. Reading should be called before "fill_guess_grids"

         call allocate_smaps(hdr_res(1)%s_dim1, hdr_res(1)%s_dim2) ! TODO decide whether the same dim is forced onto all maps (rather so)

         call read_real_arr2d_dset(file_id, n_g_smaps(LO)//"/"//dset_attrs(1), p_ratios_lo)
         call read_real_arr2d_dset(file_id, n_g_smaps(LO)//"/"//dset_attrs(2), f_ratios_lo)
         call read_real_arr2d_dset(file_id, n_g_smaps(HI)//"/"//dset_attrs(1), p_ratios_up)
         call read_real_arr2d_dset(file_id, n_g_smaps(HI)//"/"//dset_attrs(2), f_ratios_up)
         got_smaps_from_restart = .true.
      endif

   end subroutine cresp_read_smaps_from_hdf
!----------------------------------------------------------------------------------------------------
!> /brief Check if file of given name exists and contains readable solution maps. If describing parameters match, load data and proceed.
!
   subroutine try_read_user_h5(filename, hdr_init, unable_to_read)

      use constants,     only: I_ZERO, I_ONE
      use cresp_helpers, only: map_header, cresp_gname
      use cresp_io,      only: check_NR_smaps_headers
      use dataio_pub,    only: warn, printinfo
      use hdf5,          only: h5open_f, h5close_f, h5fopen_f, h5fclose_f, h5gopen_f, h5gclose_f, h5eset_auto_f, HID_T, H5F_ACC_RDONLY_F

      implicit none

      character(len=*),               intent(in)    :: filename
      type(map_header), dimension(2), intent(inout) :: hdr_init
      logical,                        intent(out)   :: unable_to_read
      integer, parameter                            :: min_fnamelen = 4 ! at least "?.h5"
      integer(kind=4)                               :: error
      integer(HID_T)                                :: file_id, group_id
      logical                                       :: file_exists, file_has_data
      logical, dimension(2)                         :: hdr_match

      unable_to_read = .true.
      if (len(filename) <= min_fnamelen) return

      inquire(file = trim(filename), exist = file_exists)    ! check if file "filename" exists

      if (file_exists) then
         file_has_data = .false.

         call h5open_f(error)
         call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, error)

         call h5eset_auto_f(I_ZERO, error)   ! Turn off printing messages
         call h5gopen_f(file_id, cresp_gname, group_id, error)

         unable_to_read = (error /= 0)

         call h5eset_auto_f(I_ONE, error)    ! Turn on  printing messages
         call h5gclose_f(group_id, error)

         if (.not. unable_to_read) then
            call cresp_read_smaps_from_hdf(file_id)   ! loads ratios and hdr_res
            call check_NR_smaps_headers(hdr_res, hdr_init, hdr_match)
            unable_to_read = .not. any(hdr_match)
            if (any(hdr_match)) call printinfo("[cresp_NR_method:try_read_user_h5] Successfully read data from provided file "//trim(filename))
         endif

         call h5fclose_f(file_id, error)
         call h5close_f(error)
      endif

      if (unable_to_read) call warn("[cresp_NR_method:try_read_user_h5] File provided as 'NR_smap_file' "//trim(filename)//" does not contain usable data.")

   end subroutine try_read_user_h5
!----------------------------------------------------------------------------------------------------
   subroutine add_dot(is_finishing)

      use constants, only: stdout
      use mpisetup,  only: master

      implicit none

      logical :: is_finishing

      if (master) then
         if (is_finishing) then
            write(stdout,'(a)', advance='yes')"."
         else
            write(stdout,'(a)', advance='no')"."
         endif
      endif

   end subroutine add_dot

end module cresp_NR_method
