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
!! The job this module is to update boundaries. This was previously inside the module sweeps . Moved here so that this can be called inside
!! solver after each time step on the entire
!<

module unsplit_fluidupdate

! pulled by ANY

   implicit none

   private
   public :: fluid_update_unsplit

contains

    subroutine fluid_update_unsplit

      use cg_list_dataop,      only: expanded_domain
      use dataio_pub,          only: halfstep
      use global,              only: dt, t
      use hdc,                 only: update_chspeed,glmdamping, eglm
      use mass_defect,         only: update_magic_mass
      use timestep_retry,      only: repeat_fluidstep
      use unsplit_sweeps,      only: unsplit_sweep
      use sources,             only: external_sources

#ifdef GRAV
      use gravity,             only: source_terms_grav, compute_h_gpot, need_update
#ifdef NBODY
      use particle_solvers,    only: psolver
#endif /* NBODY */
#endif /* GRAV */
#ifdef CRESP
      use cresp_grid,          only: cresp_update_grid, cresp_clean_grid
#endif /* CRESP */
      use user_hooks,          only: problem_customize_solution
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

      call repeat_fluidstep
      halfstep = .true.
      call update_chspeed
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
         call make_diff_sweeps(.true.)
      endif
#endif /* COSM_RAYS */

      ! At this point everything should be initialized after domain expansion and we no longer need this list.
      call expanded_domain%delete

      call eglm
      call glmdamping(.true.)
      t = t + dt

      call unsplit_sweep
#ifdef GRAV
      need_update = .true.
#ifdef NBODY
      if (associated(psolver)) call psolver(.true.)  ! this will clear need_update it it would call source_terms_grav
#endif /* NBODY */
      if (need_update) call source_terms_grav
#endif /* GRAV */
      call external_sources(.true.)
      if (associated(problem_customize_solution)) call problem_customize_solution(.true.)
      call eglm
      call glmdamping(.true.)

#ifdef CRESP
      call cresp_update_grid     ! updating number density and energy density of cosmic ray electrons via CRESP module
#endif /* CRESP */
      call update_magic_mass
#ifdef CRESP
      call cresp_clean_grid ! BEWARE: due to diffusion some junk remains in the grid - this nullifies all inactive bins.
#endif /* CRESP */
    end subroutine fluid_update_unsplit

end module unsplit_fluidupdate
