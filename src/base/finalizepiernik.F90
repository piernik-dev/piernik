! $Id$
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

      use dataio,             only: cleanup_dataio
      use decomposition,      only: cleanup_decomposition
      use diagnostics,        only: cleanup_diagnostics
      use domain,             only: cleanup_domain
      use fluidindex,         only: cleanup_fluidindex
      use global,             only: cleanup_global
      use grid,               only: cleanup_grid
      use grid_container_ext, only: cg_extptrs
      use initfluids,         only: cleanup_fluids
      use interactions,       only: cleanup_interactions
      use mpisetup,           only: cleanup_mpi
      use timer,              only: cleanup_timers
      use user_hooks,         only: cleanup_problem
#ifdef MULTIGRID
      use multigrid,          only: cleanup_multigrid
#endif /* MULTIGRID */
#ifdef GRAV
      use particle_pub,       only: cleanup_particles
#endif /* GRAV */
#ifdef PIERNIK_OPENCL
      use piernikcl,          only: cleanup_opencl
#endif /* PIERNIK_OPENCL */
#ifdef RESISTIVE
      use resistivity,        only: cleanup_resistivity
#endif /* RESISTIVE */

      implicit none

      if (associated(cleanup_problem)) call cleanup_problem; call nextdot(.false.)
      call cleanup_interactions;   call nextdot(.false.)
      call cleanup_dataio;         call nextdot(.false.)
#ifdef RESISTIVE
      call cleanup_resistivity;    call nextdot(.false.)
#endif /* RESISTIVE */
#ifdef MULTIGRID
      call cleanup_multigrid;      call nextdot(.false.)
#endif /* MULTIGRID */
      call cleanup_grid;           call nextdot(.false.)
      call cleanup_fluids;         call nextdot(.false.)
#ifdef GRAV
      call cleanup_particles;      call nextdot(.false.)
#endif /* GRAV */
      call cleanup_global;         call nextdot(.false.)
      call cleanup_decomposition;  call nextdot(.false.)
      call cleanup_domain;         call nextdot(.false.)
      call cleanup_fluidindex;     call nextdot(.false., print_t = .true.)
      call cleanup_timers;         call nextdot(.false.)
      call cleanup_diagnostics;    call nextdot(.false.)
      call cg_extptrs%epa_cleanup; call nextdot(.false.)
#ifdef PIERNIK_OPENCL
      call cleanup_opencl;         call nextdot(.false.)
#endif /* PIERNIK_OPENCL */
      call cleanup_mpi;            call nextdot(.true.)

   end subroutine cleanup_piernik

!>
!! Just print a dot on the screen, do not put a newline unless asked to do so.
!<
   subroutine nextdot(advance, print_t)

      use mpisetup,  only: master
      use constants, only: stdout
      use timer,     only: set_timer, tmr_fu

      implicit none

      logical, intent(in) :: advance
      logical, optional, intent(in) :: print_t

      if (master) then
         if (present(print_t)) then
            if (print_t) &
               write(stdout,'(f7.2,a)',advance='no') set_timer(tmr_fu), " s "
         endif
         if (advance) then
            write(stdout,'(a)')"."
         else
            write(stdout,'(a)',advance='no')"."
         endif
      endif

   end subroutine nextdot

end module finalizepiernik
