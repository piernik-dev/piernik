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
!! \brief Here we perform pairs of timesteps
!!
!! \details Each pair of timesteps consists of two timesteps of equal length
!! in order to maintain second order accuracy in time for (M)HD solver.
!! In first timestep the sweeps are executed in order X->Y->Z, then in second
!! timestep the order is reversed (Z->Y->X). The state of the fluid between
!! sweeps X-Y and Y-Z has little physical sense because only selected terms of
!! Navier-Stokes equation were applied.
!!
!! Excessive accelerations may cause violation of the CFL timestep. When timestep
!! retry is enabled, the state of the simulation is rolled back to the beginning
!! of the pair of timesteps.
!!
!! In order to maintain symmetry of the sweep operators we don't perform any
!! refinement or derefinement in the middle of the pair of timesteps.
!<

module fluidupdate   ! SPLIT
! pulled by ANY

   implicit none

   private
   public :: fluid_update

contains

!> \brief choose between fully equipped solver and the simplified HLLC one

   subroutine fluid_update

      use constants,           only: RTVD_SPLIT, RIEMANN_SPLIT, HLLC_SPLIT, RIEMANN_UNSPLIT, I_ONE, I_TWO
      use dataio_pub,          only: die
      use domain,              only: dom, is_refined
      use global,              only: which_solver
      use fluidupdate_hllc,    only: fluid_update_simple
      use ppp,                 only: ppp_main
      use unsplit_fluidupdate, only: fluid_update_unsplit

      implicit none

      character(len=*), parameter :: fu_label = "fluid_update"

      call ppp_main%start(fu_label)

      if (is_refined .and. (mod(dom%nb, I_TWO) == I_ONE)) &
           call die("[fluidupdate:fluid_update] odd number of guardcells is known to cause inaccuracies in (M)HD and nonconvergence of V-cycles")

      select case (which_solver)
         case (HLLC_SPLIT)
            call fluid_update_simple
         case (RTVD_SPLIT)
            call fluid_update_full
         case (RIEMANN_SPLIT)
            call fluid_update_full
         case (RIEMANN_UNSPLIT)
            call fluid_update_unsplit
         case default
            call die("[fluidupdate:fluid_update] unknown solver")
      end select
      call ppp_main%stop(fu_label)

   end subroutine fluid_update

!>
!! \brief Advance the solution by two timesteps using directional splitting
!<

   subroutine fluid_update_full

      use dataio_pub,     only: halfstep
      use global,         only: dt, dtm, t
      use hdc,            only: update_chspeed
      use mass_defect,    only: update_magic_mass
      use timestep_retry, only: repeat_fluidstep
#ifdef CRESP
      use cresp_grid,     only: cresp_update_grid, cresp_clean_grid
#endif /* CRESP */

      implicit none

      call repeat_fluidstep
      call update_chspeed

      halfstep = .false.
      t = t + dt

      call make_3sweeps(.true.) ! X -> Y -> Z

! Sources should be hooked to problem_customize_solution with forward argument

#ifdef CRESP
      call cresp_update_grid     ! updating number density and energy density of cosmic ray electrons via CRESP module
#endif /* CRESP */

      halfstep = .true.
      t = t + dt
      dtm = dt

      call make_3sweeps(.false.) ! Z -> Y -> X
      call update_magic_mass
#ifdef CRESP
      call cresp_clean_grid ! BEWARE: due to diffusion some junk remains in the grid - this nullifies all inactive bins.
#endif /* CRESP */

   end subroutine fluid_update_full

!>
!! \brief Perform sweeps in all three directions plus sources that are calculated every timestep
!<
   subroutine make_3sweeps(forward)

      use cg_list_dataop,      only: expanded_domain
      use constants,           only: xdim, ydim, zdim, I_ONE
      use fargo,               only: make_fargosweep
      use global,              only: skip_sweep, use_fargo
      use hdc,                 only: glmdamping, eglm
      use ppp,                 only: ppp_main
      use sources,             only: external_sources
      use sweeps,              only: sweep
      use user_hooks,          only: problem_customize_solution
#ifdef GRAV
      use gravity,             only: source_terms_grav, compute_h_gpot, need_update
#ifdef NBODY
      use particle_solvers,    only: psolver
#endif /* NBODY */
#endif /* GRAV */
#ifdef COSM_RAYS
#ifdef MULTIGRID
      use multigrid_diffusion, only: inworth_mg_diff
#else /* !MULTIGRID */
      use initcosmicrays,      only: use_CRdiff
#endif /* !MULTIGRID */
      use crdiffusion,         only: make_diff_sweeps
#endif /* COSM_RAYS */
#ifdef SHEAR
      use shear,               only: shear_3sweeps
#endif /* SHEAR */

      implicit none

      logical, intent(in) :: forward  !< If .true. then do X->Y->Z sweeps, if .false. then reverse that order

      integer(kind=4) :: s, sFRST, sLAST, sCHNG
      character(len=*), parameter :: sw3_label = "sweeps"

      if (forward) then
         sFRST = xdim ; sLAST = zdim ; sCHNG = I_ONE
      else
         sFRST = zdim ; sLAST = xdim ; sCHNG = -I_ONE
      endif

#ifdef SHEAR
      call shear_3sweeps
#endif /* SHEAR */

#ifdef GRAV
      call compute_h_gpot
#endif /* GRAV */

#ifdef COSM_RAYS
#ifdef MULTIGRID
      if (inworth_mg_diff()) then
#else /* !MULTIGRID */
      if (use_CRdiff) then
#endif /* !MULTIGRID */
         call make_diff_sweeps(forward)
      endif
#endif /* COSM_RAYS */

      ! At this point everything should be initialized after domain expansion and we no longer need this list.
      call expanded_domain%delete

      ! The following block of code may be treated as a 3D (M)HD solver.
      ! Don't put anything inside unless you're sure it should belong to the (M)HD solver.
      call ppp_main%start(sw3_label)
      if (use_fargo) then
         if (.not.skip_sweep(zdim)) call make_adv_sweep(zdim, forward)
         if (.not.skip_sweep(xdim)) call make_adv_sweep(xdim, forward)
         if (.not.skip_sweep(ydim)) call make_fargosweep
      else
         do s = sFRST, sLAST, sCHNG
            if (.not.skip_sweep(s)) call make_adv_sweep(s, forward)
         enddo
      endif
      call ppp_main%stop(sw3_label)

#ifdef GRAV
      need_update = .true.
#ifdef NBODY
      if (associated(psolver)) call psolver(forward)  ! this will clear need_update it it would call source_terms_grav
#endif /* NBODY */
      if (need_update) call source_terms_grav
#endif /* GRAV */

      call external_sources(forward)
      if (associated(problem_customize_solution)) call problem_customize_solution(forward)

      call eglm
      call glmdamping

   end subroutine make_3sweeps

!>
!! \brief Perform single sweep in forward or backward direction.
!!
!! \details Effectively this is a 3D (M)HD solver that applies only terms related to the direction dir/
!<
   subroutine make_adv_sweep(dir, forward)

      use dataio_pub,     only: die
      use domain,         only: dom
      use global,         only: geometry25D
      use sweeps,         only: sweep
#ifdef MAGNETIC
      use constants,      only: DIVB_CT, RTVD_SPLIT
      use ct,             only: magfield
      use global,         only: divB_0_method, which_solver
#endif /* MAGNETIC */
#ifdef DEBUG
      use piernikiodebug, only: force_dumps
#endif /* DEBUG */

      implicit none

      integer(kind=4), intent(in) :: dir      !< direction, one of xdim, ydim, zdim
      logical,         intent(in) :: forward  !< if .false. then reverse operation order in the sweep

#ifdef MAGNETIC
      if ((which_solver == RTVD_SPLIT) .and. (divB_0_method /= DIVB_CT)) call die("[fluidupdate:make_sweep] only CT is implemented in RTVD")
#endif /* MAGNETIC */

      ! ToDo: check if changes of execution order here (block loop, direction loop, boundary update can change
      ! cost or allow for reduction of required guardcells

      if (dom%has_dir(dir)) then
         if (.not. forward) then
#ifdef MAGNETIC
            if (divB_0_method == DIVB_CT) call magfield(dir)
#endif /* MAGNETIC */
         endif

         call sweep(dir)

         if (forward) then
#ifdef MAGNETIC
            if (divB_0_method == DIVB_CT) call magfield(dir)
#endif /* MAGNETIC */
         endif
      else
         if (geometry25D) call sweep(dir)
      endif

#ifdef DEBUG
      call force_dumps
#endif /* DEBUG */

   end subroutine make_adv_sweep

end module fluidupdate
