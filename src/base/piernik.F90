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

   use constants,   only: PIERNIK_START, PIERNIK_INITIALIZED, PIERNIK_FINISHED, PIERNIK_CLEANUP, fplen, stdout, I_ONE
   use dataio,      only: write_data, user_msg_handler, check_log, check_tsl
   use dataio_pub,  only: nend, tend, msg, printinfo, warn, die, code_progress
   use fluidupdate, only: fluid_update
   use global,      only: t, nstep, dt, dtm, cfl_violated
   use mpisetup,    only: comm, ierr, master
   use timer,       only: time_left, set_timer
   use timestep,    only: time_step
   use user_hooks,  only: finalize_problem
#ifdef PERFMON
   use timer,       only: timer_start, timer_stop
#endif /* PERFMON */

   implicit none

   logical              :: end_sim !< Used in main loop, to test whether to stop simulation or not
   character(len=fplen) :: nstr, tstr
   logical, save        :: first_step = .true.
   real                 :: ts    ! Timestep wallclock
   real                 :: tlast
   character(len=*),parameter :: tmr_fu = "fluid_update"

   ts=set_timer(tmr_fu,.true.)
   tlast = 0.0

   code_progress = PIERNIK_START

   call init_piernik

   call MPI_Barrier(comm,ierr)
!-------------------------------- MAIN LOOP ----------------------------------
#ifdef PERFMON
   call timer_start
#endif /* PERFMON */

   code_progress = PIERNIK_INITIALIZED

   end_sim = .false.

   if (master) then
      call printinfo("======================================================================================================", .false.)
      call printinfo("###############     Simulation     ###############", .false.)
   endif

   call print_progress(nstep)

   do while (t < tend .and. nstep < nend .and. .not.(end_sim) .and. time_left() ) ! main loop

      call time_step(dt)
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

      if (t == tlast .and. .not. first_step .and. .not. cfl_violated) call die("[piernik] timestep is too small: t == t + 2 * dt")

      call MPI_Barrier(comm,ierr)
      call write_data(output='all')

      call user_msg_handler(end_sim)

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

   if (master) write(stdout, '(a)', advance='no') "Finishing "
   call cleanup_piernik

contains

   subroutine print_progress(nstep)

      use constants,  only: cwdlen
      use dataio_pub, only: printinfo, msg
      use global,     only: dt, t

      implicit none

      integer(kind=4), intent(in) :: nstep

      character(len=cwdlen), parameter :: fmt900 = "('   nstep = ',i7,'   dt = ',es23.16,'   t = ',es23.16,'   dWallClock = ',f10.2,' s')"

      ts = set_timer(tmr_fu)
      if (master) then
         write(msg, fmt900) nstep, dt, t, ts
         call printinfo(msg, .true.)
      endif

   end subroutine print_progress

!>
!! Meta subroutine responsible for initializing all piernik modules
!<
   subroutine init_piernik

      use constants,             only: PIERNIK_INIT_MPI, PIERNIK_INIT_GLOBAL, PIERNIK_INIT_FLUIDS, PIERNIK_INIT_DOMAIN, PIERNIK_INIT_GRID, PIERNIK_INIT_IO_IC
      use dataio,                only: init_dataio, write_data
      use dataio_pub,            only: nrestart, wd_wr, wd_rd, par_file, tmp_log_file, msg, printio, die, warn, printinfo, require_init_prob, problem_name, run_id, code_progress
      use domain,                only: init_domain
      use diagnostics,           only: diagnose_arrays, check_environment
      use fluidboundaries,       only: all_fluid_boundaries
      use fluidboundaries_funcs, only: init_default_fluidboundaries
      use fluidindex,            only: flind
      use global,                only: init_global
      use grid,                  only: init_grid
      use gridgeometry,          only: init_geometry
      use initfluids,            only: init_fluids, sanitize_smallx_checks
      use interactions,          only: init_interactions
      use initproblem,           only: init_prob, read_problem_par, problem_pointers
      use mpisetup,              only: init_mpi
      use timestep,              only: init_time_step
      use units,                 only: init_units
#ifdef MAGNETIC
      use magboundaries,         only: all_mag_boundaries
#ifdef RESISTIVE
      use resistivity,           only: init_resistivity, compute_resist
#endif /* RESISTIVE */
#endif /* MAGNETIC */
#ifdef SHEAR
      use shear,                 only: init_shear
#endif /* SHEAR */
#ifdef GRAV
      use gravity,               only: init_grav, grav_pot_3d, grav_pot_3d_called, source_terms_grav, grav_pot_3d_bnd
      use hydrostatic,           only: cleanup_hydrostatic
#endif /* GRAV */
#ifdef MULTIGRID
      use multigrid,             only: init_multigrid
#endif /* MULTIGRID */
#ifdef SN_SRC
      use snsources,             only: init_snsources
#endif /* SN_SRC */
#ifdef DEBUG
      use piernikdebug,          only: init_piernikdebug
#endif /* DEBUG */
#ifdef CORIOLIS
      use coriolis,              only: init_coriolis
#endif /* CORIOLIS */
#ifdef COSM_RAYS
      use crdiffusion,           only: init_crdiffusion
#endif /* COSM_RAYS */

      implicit none

      call parse_cmdline
      write(par_file,'(2a)') trim(wd_rd),'problem.par'
      write(tmp_log_file,'(2a)') trim(wd_wr),'tmp.log'

      ! First, we must initialize the communication (and things that do not depend on init_mpi if there are any)
      call init_mpi
      code_progress = PIERNIK_INIT_MPI ! Now we can initialize grid and everything that depends at most on init_mpi. All calls prior to PIERNIK_INIT_GRID can be reshuffled when necessary
      call check_environment

#ifdef DEBUG
      call init_piernikdebug ! Make it available as early as possible - right after init_mpi
#endif /* DEBUG */

      call problem_pointers ! set up problem-specific pointers as early as possible to allow implementation of problem-specific hacks also during the initialization

      call init_units

      call init_default_fluidboundaries

      call init_global
      code_progress = PIERNIK_INIT_GLOBAL ! Global parameters are set up

      call init_time_step

      call init_fluids
      code_progress = PIERNIK_INIT_FLUIDS ! Fluid properties are set up

      call init_interactions ! requires flind and units

      call init_domain
      code_progress = PIERNIK_INIT_DOMAIN ! Base domain is known and initial domain decomposition is known

      call init_grid         ! Most of the cg's vars are now initialized, only arrays left
      code_progress = PIERNIK_INIT_GRID      ! Now we can initialize things that depend on all the above fundamental calls

      call init_geometry ! depends on domain; \todo move this call to init_domain or init_grid

#ifdef COSM_RAYS
      call init_crdiffusion(flind%crs%all) ! depends on grid and fluids
#endif /* COSM_RAYS */

#ifdef SHEAR
      call init_shear ! depends on fluids
#endif /* SHEAR */

#ifdef RESISTIVE
      call init_resistivity ! depends on grid
#endif /* RESISTIVE */

#ifdef GRAV
      call init_grav ! depends on units and arrays
!> \deprecated It is only temporary solution, but grav_pot_3d must be called after init_prob due to csim2,c_si,alpha clash!!!
      if (associated(grav_pot_3d)) then
         call grav_pot_3d ! depends on grav
         call grav_pot_3d_bnd
         grav_pot_3d_called = .true.
      else
#ifdef VERBOSE
         call warn("[piernik:init_piernik] grav_pot_3d is not associated! Will try to call it once more after init_problem.")
#endif /* VERBOSE */
      endif
#endif /* GRAV */

#ifdef CORIOLIS
      call init_coriolis ! depends on geometry
#endif /* CORIOLIS */

#ifdef SN_SRC
      call init_snsources ! depends on grid and fluids/cosmicrays
#endif /* SN_SRC */

#ifdef MULTIGRID
      call init_multigrid ! depends on grid, geometry, units and arrays
#endif /* MULTIGRID */

      code_progress = PIERNIK_INIT_IO_IC       ! Almost everything is initialized: do problem-related stuff here, set-up I/O and create or read the initial conditions.

      call read_problem_par ! may depend on anything but init_dataio, \todo add checks against PIERNIK_INIT_IO_IC to all initproblem::read_problem_par

      call init_dataio ! depends on units, fluids (through common_hdf5), fluidboundaries, arrays, grid and shear (through magboundaries::bnd_b or fluidboundaries::bnd_u) \todo split me

      if (master) then
         call printinfo("###############     Initial Conditions     ###############", .false.)
         write(msg,'(4a)') "   Starting problem : ",trim(problem_name)," :: ",trim(run_id)
         call printinfo(msg, .true.)
         call printinfo("", .true.)
      endif
      !>
      !! \deprecated It makes no sense to call (sometimes expensive) init_prob before reading restart file.
      !! BEWARE: If your problem requires to call init_prob add "require_init_prob = 1" to restart file
      !! Move everything that is not regenerated by restart file to read_problem_par or create separate post-restart initialization
      !<
      !> \warning Set initial conditions by hand when starting from scratch or read them from a restart file. Do not use both unless you REALLY need to do so.
      if (nrestart > 0 .and. require_init_prob /= 1) then
         if (master) then
            write(msg,'(a,i4,a)') "[piernik:init_piernik] Restart file #",nrestart," read. Skipping init_prob."
            call printio(msg)
         endif
      else
         call init_prob ! may depend on anything
         call all_fluid_boundaries !> \warning Never assume that init_prob set guardcells correctly
#ifdef MAGNETIC
         call all_mag_boundaries
#endif /* MAGNETIC */
#ifdef GRAV
         if (.not.grav_pot_3d_called) then
            if (associated(grav_pot_3d)) then
               call grav_pot_3d
               call grav_pot_3d_bnd
               grav_pot_3d_called = .true.
            else
               call die("[piernik:init_piernik] grav_pot_3d failed for the 2nd time!")
            endif
         endif
         call source_terms_grav ! make sure that all contributions to the gravitational potential are computed before first dump
         ! Possible side-effects: if variable_gp then grav_pot_3d may be called twice (second call from source_terms_grav)
         call cleanup_hydrostatic
#endif /* GRAV */
         call write_data(output='all') ! moved from dataio::init_dataio
      endif

#ifdef RESISTIVE
      call compute_resist  ! etamax%val is required by timestep_resist
#endif /* RESISTIVE */
#ifdef VERBOSE
      call diagnose_arrays ! may depend on everything
#endif /* VERBOSE */

      call sanitize_smallx_checks ! depends on init_prob || init_dataio/read_restart_hdf5

   end subroutine init_piernik

!>
!! Meta subroutine responsible for cleaning after successful run
!<
   subroutine cleanup_piernik

      use user_hooks,   only: cleanup_problem
      use dataio,       only: cleanup_dataio
      use diagnostics,  only: cleanup_diagnostics
      use domain,       only: cleanup_domain
      use global,       only: cleanup_global
      use grid,         only: cleanup_grid
      use initfluids,   only: cleanup_fluids
      use interactions, only: cleanup_interactions
      use fluidindex,   only: cleanup_fluidindex
      use fluxes,       only: cleanup_fluxes
      use mpisetup,     only: cleanup_mpi
      use timer,        only: cleanup_timers
#ifdef RESISTIVE
      use resistivity,  only: cleanup_resistivity
#endif /* RESISTIVE */
#ifdef MULTIGRID
      use multigrid,    only: cleanup_multigrid
#endif /* MULTIGRID */
      implicit none

      if (associated(cleanup_problem)) call cleanup_problem;        call  nextdot(.false.)
      call cleanup_interactions; call nextdot(.false.)
      call cleanup_fluxes;      call nextdot(.false.)
      call cleanup_dataio;      call nextdot(.false.)
#ifdef RESISTIVE
      call cleanup_resistivity; call nextdot(.false.)
#endif /* RESISTIVE */
#ifdef MULTIGRID
      call cleanup_multigrid;   call nextdot(.false.)
#endif /* MULTIGRID */
      call cleanup_grid;        call nextdot(.false.)
      call cleanup_fluids;      call nextdot(.false.)
      call cleanup_global;      call nextdot(.false.)
      call cleanup_domain;      call nextdot(.false.)
      call cleanup_fluidindex;  call nextdot(.false., print_t = .true.)
      call cleanup_timers;      call nextdot(.false.)
      call cleanup_diagnostics; call nextdot(.false.)
      call cleanup_mpi;         call nextdot(.true.)

   end subroutine cleanup_piernik

!>
!! Meta subroutine responsible for setting proper pointers or doing other magic
!! after relaxiation/grace period has passed
!<
   subroutine grace_period

      use dataio_pub,   only: printinfo
      use global,       only: grace_period_passed, relax_time
      use interactions, only: interactions_grace_passed
      use mpisetup,     only: master
      use user_hooks,   only: problem_grace_passed

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
         runned = .true.
      endif
   end subroutine grace_period
!>
!! Just print a dot on the screen, do not put a newline unless asked to do so.
!<
   subroutine nextdot(advance, print_t)

      use mpisetup,  only: master
      use constants, only: stdout

      implicit none

      logical, intent(in) :: advance
      logical, optional, intent(in) :: print_t

      if (master) then
         if (present(print_t)) then
            if (print_t) then
               ts = set_timer(tmr_fu)
               write(stdout,'(f7.2,a)',advance='no') ts," s "
            endif
         endif
         if (advance) then
            write(stdout,'(a)')"."
         else
            write(stdout,'(a)',advance='no')"."
         endif
      endif

   end subroutine nextdot

!-----------------------------------------------------------------------------

   subroutine parse_cmdline

      use constants,  only: stdout, cwdlen
      use dataio_pub, only: cmdl_nml, wd_rd, wd_wr, piernik_hdf5_version
      use version,    only: nenv,env, init_version

      implicit none

      integer :: i, j
      logical :: skip_next
      character(len=8)            :: date   ! QA_WARN len defined by ISO standard
      character(len=10)           :: time   ! QA_WARN len defined by ISO standard
      character(len=5)            :: zone   ! QA_WARN len defined by ISO standard
      character(len=cwdlen)       :: arg
!      character(len=*), parameter :: cmdlversion = '1.0'
      logical, save               :: do_time = .false.

      skip_next = .false.

      do i = 1, command_argument_count()
         if (skip_next) then
            skip_next = .false.
            cycle
         endif
         call get_command_argument(i, arg)

         select case (arg)
         case ('-v', '--version')
!            print '(2a)', 'cmdline version ', cmdlversion
            write(stdout, '(a,f5.2)') 'output version: ',piernik_hdf5_version
            call init_version
            write(stdout,'(a)') "###############     Source configuration     ###############"
            do j=1,nenv
               write(stdout,'(a)') env(j)
            enddo
            stop
         case ('-p', '--param')
            call get_command_argument(i+1,arg)
            write(wd_rd,'(a)') arg
            skip_next = .true.
         case ('-w', '--write')
            call get_command_argument(i+1,arg)
            write(wd_wr,'(a)') arg
            skip_next = .true.
         case ('-n', '--namelist')
            call get_command_argument(i+1,arg)
            write(cmdl_nml, '(3A)') cmdl_nml(1:len_trim(cmdl_nml)), " ", trim(arg)
            skip_next = .true.
         case ('-h', '--help')
            call print_help()
            stop
         case ('-t', '--time')
            do_time = .true.
         case default
            print '(a,a,/)', 'Unrecognized command-line option: ', arg
            call print_help()
            stop
         end select
      enddo

      if (wd_wr(len_trim(wd_wr):len_trim(wd_wr)) /= '/' ) write(wd_wr,'(a,a1)') trim(wd_wr),'/'
      if (wd_rd(len_trim(wd_rd):len_trim(wd_rd)) /= '/' ) write(wd_rd,'(a,a1)') trim(wd_rd),'/'

      ! Print the date and, optionally, the time
      call date_and_time(DATE=date, TIME=time, ZONE=zone)
      if (do_time) then
         write (stdout, '(a,"-",a,"-",a)', advance='no') date(1:4), date(5:6), date(7:8)
         write (stdout, '(1x,a,":",a,1x,a)') time(1:2), time(3:4), zone
         stop
      endif
   end subroutine parse_cmdline
!-----------------------------------------------------------------------------
   subroutine print_help

      implicit none

      print '(a)', 'usage: cmdline [OPTIONS]'
      print '(a)', ''
      print '(a)', 'Without further options, cmdline prints the date and exits.'
      print '(a)', ''
      print '(a)', 'cmdline options:'
      print '(a)', ''
      print '(a)', '  -v, --version     print version information and exit'
      print '(a)', '  -n, --namelist    read namelist from command line'
      print '(a)', '  -p, --param       path to the directory with problem.par and/or restarts (default: ".")'
      print '(a)', '  -w, --write       path to the directory where output will be written (default: ".")'
      print '(a)', '  -h, --help        print usage information and exit'
      print '(a)', '  -t, --time        print time and exit'

   end subroutine print_help

end program piernik
