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
!! Main program
!<
program piernik
! pulled by ANY
   use dataio,        only: write_data, user_msg_handler, check_log, check_tsl
   use dataio_pub,    only: nend, tend, msg, fplen, printinfo, warn, cwdlen, &
       &                    code_progress, PIERNIK_START, PIERNIK_INITIALIZED, PIERNIK_FINISHED, PIERNIK_CLEANUP
   use fluidupdate,   only: fluid_update
   use mpisetup,      only: comm, comm3d, ierr, proc, t, nstep, dt, dtm
   use timer,         only: time_left
   use types,         only: finalize_problem
   use timestep,      only: time_step
   use timer,         only: timer_
#ifdef PERFMON
   use timer,         only: timer_start, timer_stop
#endif /* PERFMON */

   implicit none

   logical              :: end_sim !< Used in main loop, to test whether to stop simulation or not
   character(len=fplen) :: nstr, tstr
   character(len=cwdlen), parameter :: fmt900 = "('   nstep = ',i7,'   dt = ',es22.16,'   t = ',es22.16,'   dWallClock = ',f7.2,' s')"
   logical, save                    :: first_step = .true.
   real                             :: ts   ! Timestep wallclock

   code_progress = PIERNIK_START

   call init_piernik

   call MPI_Barrier(comm3d,ierr)
!-------------------------------- MAIN LOOP ----------------------------------
#ifdef PERFMON
   call timer_start
#endif /* PERFMON */

   code_progress = PIERNIK_INITIALIZED

   end_sim = .false.

   if (proc == 0) then
      call printinfo("======================================================================================================", .false.)
      call printinfo("###############     Simulation     ###############", .false.)
   endif

   do while (t < tend .and. nstep < nend .and. .not.(end_sim) .and. time_left() ) ! main loop

      call time_step

      if (first_step) then
         ts  = 0.0
         dtm = 0.0
#ifdef RESISTIVE
         dt  = 0.0    ! BEWARE: smells like some dirty trick
#endif /* RESISTIVE */
      else
         dtm = dt
      endif
      if (proc == 0) then
         write(msg, fmt900) nstep, dt, t, ts
         call printinfo(msg, .true.)
      endif

      call check_log
      call check_tsl

      if (first_step) then
         ts=timer_("fluid_update",.true.)
         ts = 0.0
      else
         ts=timer_("fluid_update")
      endif
      call fluid_update
      ts=timer_("fluid_update")

      nstep=nstep+1

      call MPI_Barrier(comm3d,ierr)
      call write_data(output='all')

      call user_msg_handler(end_sim)

      first_step = .false.
   enddo ! main loop

   if (proc == 0) then
      write(msg, fmt900) nstep, dt, t, ts
      call printinfo(msg, .true.)
   endif

   code_progress = PIERNIK_FINISHED

   if (proc == 0) then
      write(tstr, '(g14.6)') t
      write(nstr, '(i7)') nstep
      tstr = adjustl(tstr)
      nstr = adjustl(nstr)
      call printinfo("======================================================================================================", .false.)
      call printinfo("###############     Finishing     ###############", .false.)
      if (t >= tend) then
         write(msg, '(2a)') "Simulation has reached final time t = ",trim(tstr)
         call printinfo(msg)
      endif
      if (nstep >= nend) then
         write(msg, '(4a)') "Maximum step count exceeded (",trim(nstr),") at t = ",trim(tstr)
         call warn(msg)
      endif
      if (end_sim) then
         write(msg, '(4a)') "Enforced stop at step ",trim(nstr),", t = ", trim(tstr)
         call warn(msg)
      endif
      if (.not.time_left()) then
         write(msg, '(4a)') "Wall time limit exceeded at step ",trim(nstr),", t = ", trim(tstr)
         call warn(msg)
      endif
   endif

   if (associated(finalize_problem)) call finalize_problem

#ifdef PERFMON
   call timer_stop
#endif /* PERFMON */
   call write_data(output='end')
!---------------------------- END OF MAIN LOOP ----------------------------------

   call MPI_Barrier(comm,ierr)

   code_progress = PIERNIK_CLEANUP

   if (proc == 0) write(*, '(a)', advance='no') "Finishing "       ! QA_WARN we allow non-advancing stdout on special occasions
   call cleanup_piernik

contains
!>
!! Meta subroutine responsible for initializing all piernik modules
!<
   subroutine init_piernik
      use arrays,                only: init_arrays
      use constants,             only: init_constants
      use dataio,                only: init_dataio, write_data
      use dataio_pub,            only: nrestart, cwd, par_file, tmp_log_file, msg, printio, die, warn, printinfo
      use diagnostics,           only: diagnose_arrays
      use fluidboundaries,       only: all_fluid_boundaries
      use fluidboundaries_pub,   only: init_fluidboundaries
      use fluidindex,            only: nvar
      use grid,                  only: nx, ny, nz, init_grid
      use initfluids,            only: init_fluids
      use initproblem,           only: init_prob, read_problem_par
      use problem_pub,           only: problem_name, run_id
      use mpiboundaries,         only: mpi_boundaries_prep
      use mpisetup,              only: init_mpi
      use types,                 only: grid_container
#ifdef MAGNETIC
      use magboundaries,         only: all_mag_boundaries
#ifdef RESISTIVE
      use resistivity,           only: init_resistivity
#endif /* RESISTIVE */
#endif /* MAGNETIC */
#ifdef SHEAR
      use shear,                 only: init_shear
#endif /* SHEAR */
#ifdef GRAV
      use gravity,               only: init_grav, grav_pot_3d
#endif /* GRAV */
#ifdef FLUID_INTERACTIONS
      use interactions,          only: init_interactions
#endif /* FLUID_INTERACTIONS */
#ifdef MULTIGRID
      use multigrid,             only: init_multigrid
#endif /* MULTIGRID */
#ifdef SN_SRC
      use snsources,             only: init_snsources
#endif /* SN_SRC */

      implicit none

      type(grid_container) :: cgrid
      logical :: tmp_log_exist
#ifdef GRAV
      logical              :: grav_pot_3d_called = .false.
#endif /* GRAV */

      call getarg(1, cwd)
      if (LEN_TRIM(cwd) == 0) cwd = '.'
      write(par_file,'(3a)') trim(cwd),'/problem.par'
      write(tmp_log_file,'(2a)') trim(cwd),'/tmp.log'
      inquire(file = tmp_log_file, exist = tmp_log_exist)
      if (tmp_log_exist) then
         open(3, file=tmp_log_file)
         close(3, status="delete")
      endif

      call init_mpi

      call init_constants

      call init_grid(cgrid)

      call init_fluidboundaries

#ifdef SHEAR
      call init_shear
#endif /* SHEAR */

      call read_problem_par
      if (proc == 0) then
         write(msg,'(4a)') "   Starting problem : ",trim(problem_name)," :: ",trim(run_id)
         call printinfo(msg, .true.)
         call printinfo("", .true.)
      endif

      call init_fluids(cgrid)

      call init_arrays(nx,ny,nz,nvar)

      call mpi_boundaries_prep

#ifdef RESISTIVE
      call init_resistivity
#endif /* RESISTIVE */

#ifdef GRAV
      call init_grav
! It is only temporary solution, but grav_pot_3d must be called after init_prob due to csim2,c_si,alpha clash!!!
      if (associated(grav_pot_3d)) then
         call grav_pot_3d
         grav_pot_3d_called = .true.
      else
#ifdef VERBOSE
         call warn("[piernik:init_piernik] grav_pot_3d is not associated! Will try to call it once more after init_problem.")
#endif /* VERBOSE */
      endif
#endif /* GRAV */

#ifdef FLUID_INTERACTIONS
      call init_interactions
#endif /* FLUID_INTERACTIONS */

#ifdef SN_SRC
      call init_snsources
#endif /* SN_SRC */

#ifdef MULTIGRID
      call init_multigrid(cgrid)
#endif /* MULTIGRID */

      call init_dataio

      ! It makes no sense to call (sometimes expensive) init_prob before reading restart file.
      ! BEWARE: Current change may break some problems that depend on other things set in init_prob.
      ! Move everything that is not regenerated by restart file to read_problem_par or create separate post-restart initialization
      if (nrestart>0) then
         if (proc == 0) then
            write(msg,'(a,i4,a)') "[piernik:init_piernik] Restart file #",nrestart," read. Skipping init_prob."
            call printio(msg)
         endif
      else
         call init_prob
         call all_fluid_boundaries ! Never assume that init_prob set guardcells correctly
#ifdef MAGNETIC
         call all_mag_boundaries
#endif /* MAGNETIC */
#ifdef GRAV
         if (.not.grav_pot_3d_called) then
            if (associated(grav_pot_3d)) then
               call grav_pot_3d
               grav_pot_3d_called = .true.
            else
               call die("[piernik:init_piernik] grav_pot_3d failed for the 2nd time!")
            endif
         endif
#endif /* GRAV */
         call write_data(output='all') ! moved from dataio::init_dataio
      endif
#ifdef VERBOSE
      call diagnose_arrays
#endif /* VERBOSE */

   end subroutine init_piernik

!>
!! Meta subroutine responsible for cleaning after successful run
!<

   subroutine cleanup_piernik

      use arrays,      only: cleanup_arrays
      use dataio,      only: cleanup_dataio
      use grid,        only: cleanup_grid
      use initfluids,  only: cleanup_fluids
      use fluidindex,  only: cleanup_fluidindex
      use mpisetup,    only: cleanup_mpi
      use timer,       only: cleanup_timers
#ifdef RESISTIVE
      use resistivity, only: cleanup_resistivity
#endif /* RESISTIVE */
#ifdef MULTIGRID
      use multigrid,   only: cleanup_multigrid
#endif /* MULTIGRID */
      implicit none

      call cleanup_grid;        call nextdot(.false.)
      call cleanup_dataio;      call nextdot(.false.)
#ifdef RESISTIVE
      call cleanup_resistivity; call nextdot(.false.)
#endif /* RESISTIVE */
#ifdef MULTIGRID
      call cleanup_multigrid;   call nextdot(.false.)
#endif /* MULTIGRID */
      call cleanup_arrays;      call nextdot(.false.)
      call cleanup_fluids;      call nextdot(.false.)
      call cleanup_fluidindex;  call nextdot(.false.)
      call cleanup_timers;      call nextdot(.false.)
      call cleanup_mpi;         call nextdot(.true.)

   end subroutine cleanup_piernik

!>
!! Just print a dot on the screen, do not put a newline unless asked to do so.
!<

   subroutine nextdot(advance)

      use mpisetup,      only: proc

      implicit none

      logical, intent(in) :: advance

      if (proc == 0) then
         if (advance) then
            write(*,'(a)')"."              ! QA_WARN we allow non-advancing stdout on special occasions
         else
            write(*,'(a)',advance='no')"." ! QA_WARN we allow non-advancing stdout on special occasions
         endif
      endif

   end subroutine nextdot

end program piernik
