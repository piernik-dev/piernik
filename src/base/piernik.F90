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

   use all_boundaries,    only: all_bnd
   use cg_leaves,         only: leaves
   use cg_list_global,    only: all_cg
   use constants,         only: PIERNIK_START, PIERNIK_INITIALIZED, PIERNIK_FINISHED, PIERNIK_CLEANUP, fplen, stdout, I_ONE, CHK, FINAL_DUMP, tmr_fu
   use dataio,            only: write_data, user_msg_handler, check_log, check_tsl, dump
   use dataio_pub,        only: nend, tend, msg, printinfo, warn, die, code_progress
   use finalizepiernik,   only: cleanup_piernik
   use fluidindex,        only: flind
   use fluidupdate,       only: fluid_update
   use func,              only: operator(.equals.)
   use global,            only: t, nstep, dt, dtm, cfl_violated
   use initpiernik,       only: init_piernik
   use list_of_cg_lists,  only: all_lists
   use mpisetup,          only: master, piernik_MPI_Barrier, piernik_MPI_Bcast
   use named_array_list,  only: qna, wna
   use refinement,        only: emergency_fix
   use refinement_update, only: update_refinement
   use timer,             only: walltime_end, set_timer
   use timestep,          only: time_step
   use user_hooks,        only: finalize_problem, problem_domain_update
#ifdef PERFMON
   use domain,            only: dom
   use timer,             only: timer_start, timer_stop
#endif /* PERFMON */
#if defined DEBUG && defined GRAV
   use particle_pub,      only: pset
#endif /* DEBUG && GRAV */

   implicit none

   logical              :: end_sim             !< Used in main loop, to test whether to stop simulation or not
   logical, save        :: tleft = .true.      !< Used in main loop, to test whether to stop simulation or not
   character(len=fplen) :: nstr, tstr
   logical, save        :: first_step = .true.
   real                 :: ts                  !< Timestep wallclock
   real                 :: tlast
   logical              :: try_rebalance

   try_rebalance = .false.
   ts=set_timer(tmr_fu,.true.)
   tlast = 0.0

   code_progress = PIERNIK_START

   call init_piernik

   call all_cg%check_na
   !call all_cg%check_for_dirt
#if defined DEBUG && defined GRAV
   call pset%print
#endif /* DEBUG && GRAV */

   call piernik_MPI_Barrier
!-------------------------------- MAIN LOOP ----------------------------------
#ifdef PERFMON
   call timer_start
#endif /* PERFMON */

   code_progress = PIERNIK_INITIALIZED

   end_sim = .false.

   if (master) then
      call printinfo("======================================================================================================", .false.)
      call printinfo("###############     Simulation     ###############", .false.)
      call printinfo("Named arrays present at start:", to_stdout=.false.)
      call qna%print_vars(to_stdout=.false.)
      call wna%print_vars(to_stdout=.false.)
      call printinfo("Grid lists present at start:", to_stdout=.false.)
   endif
   call all_lists%print(to_stdout = .false.)  ! needs all procs to participate
   if (master) then
      call printinfo("======================================================================================================", .false.)
   endif

   call print_progress(nstep)

   do while (t < tend .and. nstep < nend .and. .not.(end_sim)) ! main loop

      dump(:) = .false.
      if (associated(problem_domain_update)) then
         call problem_domain_update
         if (emergency_fix) try_rebalance = .true.
         call update_refinement(refinement_fixup_only=.true.)
         ! Full refinement here called rebalancing, which sometimes caused problems with initialization of expanded parts of the computational domain
      endif

      call all_cg%check_na
      !call all_cg%check_for_dirt

      call time_step(dt, flind)
      call grace_period

      if (first_step) then
         dtm = 0.0
      else
         if (.not.cfl_violated) dtm = dt
      endif

      if (.not.cfl_violated) then
        call check_log
        call check_tsl
      endif

      if (.not.cfl_violated) tlast = t
      call fluid_update
      nstep = nstep + I_ONE
      call print_progress(nstep)

      if ((t .equals. tlast) .and. .not. first_step .and. .not. cfl_violated) call die("[piernik] timestep is too small: t == t + 2 * dt")

      call piernik_MPI_Barrier

      call write_data(output=CHK)

      call user_msg_handler(end_sim)
      call update_refinement
      if (try_rebalance) then
         !> \todo try to rewrite this ugly chain of flags passed through global variables into something more fool-proof
         call leaves%balance_and_update(" (re-balance) ")
         call all_bnd ! For some strange reasons this call prevents MPI-deadlock
         try_rebalance = .false.
      endif

      if (master) tleft = walltime_end%time_left()
      call piernik_MPI_Bcast(tleft)

      if (.not.tleft) end_sim = .true.

      first_step = .false.
   enddo ! main loop

   code_progress = PIERNIK_FINISHED

   if (master) then
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
      if (.not.tleft) then
         write(msg, '(4a)') "Wall time limit exceeded at step ",trim(nstr),", t = ", trim(tstr)
         call warn(msg)
      endif

      call printinfo("Named arrays present at finish:", to_stdout=.false.)
      call qna%print_vars(to_stdout=.false.)
      call wna%print_vars(to_stdout=.false.)
      call printinfo("Grid lists present at finish:", to_stdout=.false.)
   endif
   call all_lists%print(to_stdout = .false.)  ! needs all procs to participate
   if (master) then
      call printinfo("======================================================================================================", .false.)
   endif

   if (associated(finalize_problem)) call finalize_problem

#ifdef PERFMON
   call timer_stop(nstep, dom%total_ncells)
#endif /* PERFMON */
   call write_data(output=FINAL_DUMP)
!---------------------------- END OF MAIN LOOP ----------------------------------

   call piernik_MPI_Barrier

   code_progress = PIERNIK_CLEANUP

   if (master) write(stdout, '(a)', advance='no') "Finishing "
   call cleanup_piernik

contains

   subroutine print_progress(nstep)

      use constants,  only: tmr_fu
      use dataio_pub, only: printinfo, msg
      use global,     only: dt, t
      use timer,      only: set_timer, get_timestamp

      implicit none

      integer(kind=4), intent(in) :: nstep
      character(len=*), parameter :: fmt900 = "('   nstep = ',i7,'   dt = ',es23.16,'   t = ',es23.16,'   dWallClock = ',f10.2,' s ',a)"
      logical, save :: first_run = .true.


      if (master) then
         if (first_run) then
            first_run = (set_timer(tmr_fu) < 0.0)
         endif
         write(msg, fmt900) nstep, dt, t, set_timer(tmr_fu), get_timestamp()
         call printinfo(msg, .true.)
      endif

   end subroutine print_progress
!>
!! Meta subroutine responsible for setting proper pointers or doing other magic
!! after relaxation/grace period has passed
!<
   subroutine grace_period

      use all_boundaries, only: all_fluid_boundaries
      use dataio_pub,     only: printinfo
      use global,         only: grace_period_passed, relax_time
      use interactions,   only: interactions_grace_passed
      use mpisetup,       only: master
      use user_hooks,     only: problem_grace_passed

      implicit none

      logical, save     :: runned = .false.

      if (runned) return
      if (grace_period_passed()) then
         if (relax_time > 0.0) then
            ! write info message only if relax_time was set
            write(msg,'(A,ES10.3)') "[piernik:grace_period] grace period has passed.", t
            if (master) call printinfo(msg)
         endif
         call interactions_grace_passed
         if (associated(problem_grace_passed)) call problem_grace_passed
         call all_fluid_boundaries
         runned = .true.
      endif
   end subroutine grace_period

end program piernik
