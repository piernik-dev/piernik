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
!! \brief Preconditioned conjugate gradient solver.
!!
!! \details Preconditioned conjugate gradient solver. By default it uses single V-cycle as a preconditioner as it is most efficient.
!! Other possibilities are just smoothing and no preconditioning at all - both provided for academic interest only.
!!
!! Note that it is theoretically possible to use most of this module without Piernik's multigrid module.
!! We decided to leave dependencies on multigrid to keep the code simple.
!<

module pcg
! pulled by MULTIGRID && SELF_GRAV

   use constants, only: cbuff_len, dsetnamelen

   implicit none

   private
   public :: pcg_init, mgpcg, &
        &    use_CG, use_CG_outer, preconditioner, default_preconditioner, cg_corr

   logical                               :: use_CG                                !< .true. if we want to use multigrid-preconditioned conjugate gradient (MG-PCG) iterations
   logical                               :: use_CG_outer                          !< .true. if we want to use MG-PCG iterations for outer potential
   character(len=dsetnamelen), parameter :: cg_corr_n = "cg_correction"           !< correction vector for CG
   integer(kind=4)                       :: cg_corr                               !< index of the cg-correction vector
   character(len=cbuff_len)              :: preconditioner                        !< Multigrid (Huang-Greengard V-cycle) by default
   character(len=cbuff_len), parameter   :: default_preconditioner = "HG_V-cycle" !< default preconditioner

contains

!> \brief Register additional array required for PCG

   subroutine pcg_init

      use cg_list_global,   only: all_cg
      use constants,        only: base_level_id
      use dataio_pub,       only: warn
      use domain,           only: dom
      use mpisetup,         only: master
      use named_array_list, only: qna
      use refinement,       only: level_max, bsize

      implicit none

      if (use_CG .or. use_CG_outer) then
         call all_cg%reg_var(cg_corr_n)
         cg_corr = qna%ind(cg_corr_n)
         if (master) then
            call warn("[pcg:pcg_init] Multigrid-preconditioned conjugate gradient solver is experimental!")
            if (use_CG_outer) &
                 call warn("[pcg:pcg_init] Current implementation of Multigrid-Preconditioned Conjugate Gradient Method is known to converge poorly for outer potential.")
            if (level_max > base_level_id .and. any((bsize > 0) .and. dom%has_dir)) &
                 call warn("[pcg:pcg_init] Current implementation of Multigrid-Preconditioned Conjugate Gradient Method is known to converge poorly on refined grids.")
         endif
      endif

   end subroutine pcg_init

!>
!! \brief Preconditioned conjugate gradient solver
!!
!! \details This solver can converge faster than regular multigrid on domains that require high cell anisotropy e.g cg%dx >> cg%dy
!!
!! Note that forcing alpha = 1 and beta = 0 reduces this algorithm to bare multigrid
!!
!! Meaning of the symbols in comments:
!! {b}      - source (4 pi G * density)
!! {A}      - operator matrix (discrete Laplacian of chosen order and representation)
!! {x}_k    - solution (k-th approximation of the gravitational potential)
!! {r}_k    - residual (defect of the current solution)
!! {M}^{-1} - preconditioning matrix, approximation of A^{-1} (here, by default, we use multigrid with relaxation to cheaply obtain it)
!! {z}_k    - approximate correction for the solution obtained with the help of M^{-1} acting on current defect (preconditioner)
!! {d}_k    - the correction direction found by the conjugate gradient algorithm,
!!
!! For the derivation of the conjugate gradient method see:
!! Jonathan Richard Shewchuk, An Introduction to the Conjugate Gradient Method Without the Agonizing Pain
!! The algorithm of preconditioned conjugate gradient method is described in Ch. 12 of the above publication.
!!
!! \todo Find out why it converges so poorly for outer potential
!<

   subroutine mgpcg(max_cycles, norm_tol)

      use cg_leaves,          only: leaves
      use cg_list_dataop,     only: ind_val
      use constants,          only: tmr_mg, V_VERBOSE
      use dataio_pub,         only: msg, printinfo
      use mpisetup,           only: master
      use multigrid_Laplace,  only: residual, vT_A_v
      use multigridvars,      only: source, solution, defect, correction, tot_ts, ts
      use timer,              only: set_timer

      implicit none

      integer(kind=4),    intent(in)    :: max_cycles  !< Maximum allowed number of iterations
      real,               intent(in)    :: norm_tol    !< stop iterations when the ratio of norms ||residual||/||source|| is below this value

      real    :: alpha, beta, dc_k, dc_k1
      integer :: k
      real    :: norm_rhs, norm_lhs, norm_old
      real    :: loc_ts

      loc_ts = 0
      beta = 0.
      norm_rhs = leaves%norm_sq(source)
      norm_old = norm_rhs

      call residual(leaves, source, solution, defect)                                             ! {r}_0 := {b} - {A x}_0
      call precond(defect, correction)                                                            ! {z}_0 := {M}^{-1} {r}_0
      call leaves%q_copy(correction, cg_corr)                                                     ! {d}_0 := {z}_0
      dc_k = leaves%scalar_product(defect, correction)                                            ! dc_k := {r}_k^{T} {z}_k
      do k = 0, max_cycles

         alpha = dc_k/vT_A_v(cg_corr)                                                             ! \alpha_k := {{r}_k^{T} {z}_k}/{{d}_k^{T} {A d}_k}
         call leaves%q_lin_comb( [ ind_val(solution, 1.), ind_val(cg_corr, alpha) ], solution)    ! {x}_{k+1} := {x}_k + \alpha_k {d}_k
         call residual(leaves, source, solution, defect)                                          ! {r}_{k+1} := {r}_k - \alpha_k {A d}_k = {b} - {A x}_{k+1}
         !OPT: Use {A d}_k computed in vT_A_v_order to avoid communication required for residual (costs memory for storing one field, risk of drift due to boundary values)
         !OPT: Alternatively pass a flag to the residual routine that effectively switches off guardcell update as they're already have been updated by vT_A_v_order
         norm_lhs = leaves%norm_sq(defect)
         ts = set_timer(tmr_mg)
         tot_ts = tot_ts + ts
         loc_ts = loc_ts + ts
         if (norm_old/norm_lhs < 1e6) then
            write(msg,'(a,i3,a,f12.8,a,f10.2,a,f11.7,g14.6,a,f8.3)')"MG-PCG: ", k, " lhs/rhs= ",norm_lhs/norm_rhs, " improvement= ",norm_old/norm_lhs, &
                 &                                                 " a,b= ", alpha, beta, " time=", ts
         else
            write(msg,'(a,i3,a,f12.8,a,es10.3,a,f11.7,g14.6,a,f8.3)')"MG-PCG: ", k, " lhs/rhs= ",norm_lhs/norm_rhs, " improvement= ",norm_old/norm_lhs, &
                 &                                                 " a,b= ", alpha, beta, " time=", ts
         endif
         if (master) call printinfo(msg, V_VERBOSE)
         if (norm_lhs/norm_rhs <= norm_tol) exit                                                  ! if {r}_{k+1} is sufficiently small then exit loop
         norm_old = norm_lhs
         call precond(defect, correction)                                                         ! {z}_{k+1} := {M}^{-1} {r}_{k+1}
         dc_k1 = leaves%scalar_product(defect, correction)                                        ! dc_k1 := {z}_{k+1}^{T} {r}_{k+1}
         beta = dc_k1/dc_k                                                                        ! \beta_{k+1} := {{z}_{k+1}^{T} {r}_{k+1}}/{{z}_k^{T} {r}_k}
         call leaves%q_lin_comb( [ ind_val(correction, 1.), ind_val(cg_corr, beta) ], cg_corr )   ! {d}_{k+1} := {z}_{k+1} + \beta_{k+1} {d}_k
         dc_k = dc_k1

      enddo

      ts = set_timer(tmr_mg)
      tot_ts = tot_ts + ts
      loc_ts = loc_ts + ts
      write(msg, '(a,i4,a,f8.3)')"MG-PCG: ",k," iterations, total time spent =", loc_ts
      if (master) call printinfo(msg, V_VERBOSE)

   end subroutine mgpcg

!>
!! \brief Call selected preconditioner
!!
!! \details The V-cycle seems to be the most efficient preconditioner as yet.
!! It may be beneficial to increase nsmool from 4 to 8 or even 12 to get best convergence and minimize CPU time.
!!
!! Performing 2 V-cycles instead of 1 seems to cost more than it saves due to improvements to convergence rate.
!<

   subroutine precond(def, corr)

      use cg_leaves,  only: leaves
      use dataio_pub, only: die

      implicit none

      integer(kind=4), intent(in) :: def  !< Defect (source, residual)
      integer(kind=4), intent(in) :: corr !< Approximate solution for the defect (correction)

      select case (preconditioner)
         case (default_preconditioner, "MG", "V-cycle", "V")
            call single_v_cycle(def, corr)
         case ("smooth")
            call smoother(def, corr)
         case ("none")
            call leaves%q_copy(def, corr) ! non-preconditioned cg
         case default
            call die("[pcg:precond] Unknown preconditioner")
      end select

   end subroutine precond

!> \brief A single Huang-Greengard V-cycle

   subroutine single_v_cycle(def, corr)

      use cg_level_coarsest,        only: coarsest
      use cg_level_connected,       only: cg_level_connected_t
      use cg_level_finest,          only: finest
      use cg_list_global,           only: all_cg
      use constants,                only: dirtyH1
      use multigrid_gravity_helper, only: approximate_solution

      implicit none

      integer(kind=4), intent(in) :: def  !< Defect (source, residual)
      integer(kind=4), intent(in) :: corr !< Approximate solution for the defect (correction)

      type(cg_level_connected_t), pointer :: curl

      ! the Huang-Greengard V-cycle
      call finest%level%restrict_to_floor_q_1var(def)

      call all_cg%set_dirty(corr, 0.89*dirtyH1)

      curl => coarsest%level
      do while (associated(curl))
         call approximate_solution(curl, def, corr)
         call curl%check_dirty(corr, "Vup1 relax+")
         curl => curl%finer
      enddo

   end subroutine single_v_cycle

!> \brief Smoother as a preconditioner

   subroutine smoother(def, corr)

      use cg_level_finest,   only: finest
      use multigrid_Laplace, only: approximate_solution_relax
      use multigridvars,     only: nsmool

      implicit none

      integer(kind=4), intent(in) :: def  !< Defect (source, residual)
      integer(kind=4), intent(in) :: corr !< Approximate solution for the defect (correction)

      call finest%level%set_q_value(corr, 0.)
      call approximate_solution_relax(finest%level, def, corr, nsmool)

   end subroutine smoother

end module pcg
