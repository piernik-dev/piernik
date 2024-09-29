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
!! \brief Module for a general cleaning up routine after successful run
!<
module finalizepiernik

   implicit none
   private
   public :: cleanup_piernik

contains
!>
!! Meta subroutine responsible for cleaning after successful run
!<
   subroutine cleanup_piernik

      use dataio_pub,            only: flush_to_log
      use decomposition,         only: cleanup_decomposition
      use diagnostics,           only: cleanup_diagnostics
      use domain,                only: cleanup_domain
      use fargo,                 only: cleanup_fargo
      use fluidindex,            only: cleanup_fluidindex
      use global,                only: cleanup_global
      use grid,                  only: cleanup_grid
      use grid_container_ext,    only: cg_extptrs
      use initfluids,            only: cleanup_fluids
      use interactions,          only: cleanup_interactions
      use procnames,             only: pnames
      use tag_pool,              only: t_pool
      use timer,                 only: cleanup_timers
      use unified_ref_crit_list, only: urc_list
      use user_hooks,            only: cleanup_problem
#ifdef RANDOMIZE
      use randomization,         only: cleanup_randomization
#endif /* RANDOMIZE */
#ifdef MULTIGRID
      use multigrid,             only: cleanup_multigrid
#endif /* MULTIGRID */
#if defined(GRAV) && defined(NBODY)
      use particle_pub,          only: cleanup_particles
#endif /* GRAV && NBODY */
#ifdef PIERNIK_OPENCL
      use piernikcl,             only: cleanup_opencl
#endif /* PIERNIK_OPENCL */
#ifdef RESISTIVE
      use resistivity,           only: cleanup_resistivity
#endif /* RESISTIVE */
#ifdef COSM_RAYS
      use cr_data,               only: cleanup_cr_species
#endif /* COSM_RAYS */
#ifdef THERM
      use thermal,               only: cleanup_thermal
#endif /* THERM */

      implicit none

      call flush_to_log

      if (associated(cleanup_problem)) call cleanup_problem; call nextdot
      call urc_list%cleanup;       call nextdot
      call t_pool%cleanup;         call nextdot
      call cleanup_interactions;   call nextdot
      call cleanup_fargo;          call nextdot
#ifdef RESISTIVE
      call cleanup_resistivity;    call nextdot
#endif /* RESISTIVE */
      call cleanup_grid;           call nextdot
#ifdef MULTIGRID
      call cleanup_multigrid;      call nextdot
#endif /* MULTIGRID */
      call cleanup_fluids;         call nextdot
#if defined(GRAV) && defined(NBODY)
      call cleanup_particles;      call nextdot
#endif /* GRAV && NBODY */
      call cleanup_global;         call nextdot
      call cleanup_decomposition;  call nextdot
      call cleanup_domain;         call nextdot
      call cleanup_fluidindex;     call nextdot(print_t = .true.)
      call cleanup_timers;         call nextdot
      call cleanup_diagnostics;    call nextdot
      call cg_extptrs%epa_cleanup; call nextdot
      call pnames%cleanup;         call nextdot
#ifdef PIERNIK_OPENCL
      call cleanup_opencl;         call nextdot
#endif /* PIERNIK_OPENCL */
#ifdef RANDOMIZE
      call cleanup_randomization;  call nextdot
#endif /* RANDOMIZE */
#ifdef COSM_RAYS
      call cleanup_cr_species;     call nextdot
#endif /* COSM_RAYS */
#ifdef THERM
     call cleanup_thermal;         call nextdot
#endif /* THERM */

   end subroutine cleanup_piernik

!>
!! Just print a dot on the screen, do not put a newline unless asked to do so.
!<
   subroutine nextdot(print_t)

      use constants, only: stdout, tmr_fu
      use mpisetup,  only: master
      use timer,     only: set_timer

      implicit none

      logical, optional, intent(in) :: print_t

      if (master) then
         if (present(print_t)) then
            if (print_t) write(stdout,'(f7.2,a)',advance='no') set_timer(tmr_fu), " s "
         endif
         write(stdout,'(a)',advance='no')"."
      endif

   end subroutine nextdot

end module finalizepiernik
