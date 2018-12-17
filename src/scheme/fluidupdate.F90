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
!! \brief COMMENT ME
!<

module fluidupdate   ! SPLIT
! pulled by RTVD || RIEMANN

   implicit none

   private
   public :: fluid_update

contains

!>
!! \brief Advance the solution by two timesteps using directional splitting
!<

   subroutine fluid_update

      use dataio_pub,     only: halfstep
      use global,         only: dt, dtm, t
      use mass_defect,    only: update_magic_mass
      use timestep_retry, only: repeat_fluidstep
#ifdef RIEMANN
      use hdc,            only: update_chspeed
#endif /* RIEMANN */

      implicit none

      call repeat_fluidstep

#ifdef RIEMANN
      call update_chspeed
#endif /* RIEMANN */

      halfstep = .false.
      t = t + dt
      call make_3sweeps(.true.) ! X -> Y -> Z

! Sources should be hooked to problem_customize_solution with forward argument

      halfstep = .true.
      t = t + dt
      dtm = dt
      call make_3sweeps(.false.) ! Z -> Y -> X
      call update_magic_mass

   end subroutine fluid_update

!>
!! \brief Perform sweeps in all three directions plus sources that are calculated every timestep
!<
   subroutine make_3sweeps(forward)

      use cg_list_dataop,      only: expanded_domain
      use constants,           only: xdim, ydim, zdim, I_ONE
      use global,              only: skip_sweep, use_fargo
      use fargo,               only: make_fargosweep
      use sweeps,              only: sweep
      use user_hooks,          only: problem_customize_solution
#ifdef GRAV
      use global,              only: t, dt
      use gravity,             only: source_terms_grav
      use particle_pub,        only: psolver
      use particle_types,      only: pset
#endif /* GRAV */
#if defined(COSM_RAYS) && defined(MULTIGRID)
      use all_boundaries,      only: all_fluid_boundaries
      use initcosmicrays,      only: use_split
      use multigrid_diffusion, only: multigrid_solve_diff
#endif /* COSM_RAYS && MULTIGRID */
#ifdef SHEAR
      use shear,               only: shear_3sweeps
#endif /* SHEAR */
#ifdef RIEMANN
      use hdc,                 only: glmdamping, eglm
#endif /* RIEMANN */

      implicit none

      logical, intent(in) :: forward  !< If .true. then do X->Y->Z sweeps, if .false. then reverse that order

      integer(kind=4) :: s

#ifdef SHEAR
      call shear_3sweeps
#endif /* SHEAR */

#ifdef GRAV
      call source_terms_grav
#endif /* GRAV */

#if defined(COSM_RAYS) && defined(MULTIGRID)
      if (.not. use_split) then
         call multigrid_solve_diff
         call all_fluid_boundaries
      endif
#endif /* COSM_RAYS && MULTIGRID */

      call expanded_domain%delete ! at this point everything should be initialized after domain expansion and we no longer need this list
      if (use_fargo) then
         if (.not.skip_sweep(zdim)) call make_sweep(zdim, forward)
         if (.not.skip_sweep(xdim)) call make_sweep(xdim, forward)
         if (.not.skip_sweep(ydim)) call make_fargosweep
      else
         if (forward) then
            do s = xdim, zdim
               if (.not.skip_sweep(s)) call make_sweep(s, forward)
            enddo
         else
            do s = zdim, xdim, -I_ONE
               if (.not.skip_sweep(s)) call make_sweep(s, forward)
            enddo
         endif
      endif
#ifdef GRAV
      if (associated(psolver)) call pset%evolve(psolver, t-dt, dt, forward)
#endif /* GRAV */
      if (associated(problem_customize_solution)) call problem_customize_solution(forward)

#ifdef RIEMANN
      call eglm
      call glmdamping
#endif /* RIEMANN */

   end subroutine make_3sweeps

!>
!! \brief Perform single sweep in forward or backward direction
!<
   subroutine make_sweep(dir, forward)

      use constants,      only: DIVB_CT
      use dataio_pub,     only: die
      use domain,         only: dom
      use global,         only: geometry25D, divB_0_method
      use sweeps,         only: sweep
#ifdef COSM_RAYS
      use crdiffusion,    only: cr_diff
      use initcosmicrays, only: use_split
#endif /* COSM_RAYS */
#ifdef DEBUG
      use piernikiodebug, only: force_dumps
#endif /* DEBUG */
#ifdef MAGNETIC
      use ct,             only: magfield
#endif /* MAGNETIC */
#ifdef RIEMANN
      use bfc_bcc,        only: bfc2bcc
      use global,         only: force_cc_mag
#endif /* RIEMANN */

      implicit none

      integer(kind=4), intent(in) :: dir      !< direction, one of xdim, ydim, zdim
      logical,         intent(in) :: forward  !< if .false. then reverse operation order in the sweep

#ifdef RTVD
      if (divB_0_method /= DIVB_CT) call die("[fluidupdate:make_sweep] only CT is implemented in RTVD")
#endif /* RTVD */

      ! ToDo: check if changes of execution order here (block loop, direction loop, boundary update can change
      ! cost or allow for reduction of required guardcells
#ifdef RIEMANN
      if (.not. force_cc_mag) call bfc2bcc
#endif /* RIEMANN */

      if (dom%has_dir(dir)) then
         if (.not. forward) then
#ifdef COSM_RAYS
            if (use_split) call cr_diff(dir)
#endif /* COSM_RAYS */
#ifdef MAGNETIC
            if (divB_0_method == DIVB_CT) call magfield(dir)
#endif /* MAGNETIC */
         endif

         call sweep(dir)

         if (forward) then
#ifdef MAGNETIC
            if (divB_0_method == DIVB_CT) call magfield(dir)
#endif /* MAGNETIC */
#ifdef COSM_RAYS
            if (use_split) call cr_diff(dir)
#endif /* COSM_RAYS */
         endif
      else
         if (geometry25D) call sweep(dir)
      endif

#ifdef DEBUG
      call force_dumps
#endif /* DEBUG */

   end subroutine make_sweep

end module fluidupdate
