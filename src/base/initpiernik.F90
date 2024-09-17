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

module initpiernik

   implicit none
   private
   public :: init_piernik

contains
!>
!! Meta subroutine responsible for initializing all piernik modules
!! \deprecated remove the "__INTEL_COMPILER" clauses as soon as Intel Compiler gets required features and/or bug fixes
!! \todo move init_geometry call to init_domain or init_grid
!! \todo add checks against PIERNIK_INIT_IO_IC to all initproblem::read_problem_par
!! \todo split init_dataio
!<

   subroutine init_piernik

      use all_boundaries,        only: all_bnd, all_bnd_vital_q
      use cg_level_finest,       only: finest
      use cg_list_global,        only: all_cg
      use constants,             only: PIERNIK_INIT_MPI, PIERNIK_INIT_GLOBAL, PIERNIK_INIT_FLUIDS, PIERNIK_INIT_DOMAIN, &
           &                           PIERNIK_INIT_GRID, PIERNIK_INIT_IO_IC, PIERNIK_POST_IC, &
           &                           INCEPTIVE, tmr_fu, cbuff_len, PPP_PROB, &
           &                           v_name, V_LOG, V_DEBUG, V_INFO, V_ESSENTIAL, V_LOWEST, V_HIGHEST
      use dataio,                only: init_dataio, init_dataio_parameters, write_data
      use dataio_pub,            only: nrestart, restarted_sim, wd_rd, par_file, tmp_log_file, msg, printio, printinfo, piernik_verbosity, &
           &                           warn, die, require_problem_IC, problem_name, run_id, code_progress, log_wr, set_colors
      use decomposition,         only: init_decomposition
      use domain,                only: init_domain
      use diagnostics,           only: diagnose_arrays, check_environment
      use fargo,                 only: init_fargo
      use fluidboundaries_funcs, only: init_default_fluidboundaries
      use global,                only: init_global
      use grid,                  only: init_grid
      use grid_container_ext,    only: cg_extptrs
      use gridgeometry,          only: init_geometry
      use hdc,                   only: init_psi
      use initfluids,            only: init_fluids, sanitize_smallx_checks
      use initproblem,           only: problem_initial_conditions, read_problem_par, problem_pointers
      use interpolations,        only: set_interpolations
      use lb_helpers,            only: costs_maintenance
      use load_balance,          only: init_load_balance
      use memory_usage,          only: init_memory
      use mpisetup,              only: init_mpi, master
      use ppp,                   only: init_profiling, ppp_main
      use procnames,             only: pnames
      use refinement,            only: init_refinement, level_max
      use refinement_update,     only: update_refinement
      use sources,               only: init_sources
      use timer,                 only: set_timer
      use unified_ref_crit_list, only: urc_list
      use units,                 only: init_units
      use user_hooks,            only: problem_post_restart, problem_post_IC
      use wrapper_stats,         only: init_wstats
#ifdef RESISTIVE
      use resistivity,           only: init_resistivity
#endif /* RESISTIVE */
#ifdef GRAV
      use gravity,               only: init_grav, init_terms_grav, source_terms_grav
      use hydrostatic,           only: init_hydrostatic, cleanup_hydrostatic
#endif /* GRAV */
#ifdef NBODY
      use particle_gravity,      only: update_particle_gravpot_and_acc
      use particle_pub,          only: init_particles
      use particle_solvers,      only: init_psolver, update_particle_kinetic_energy
      use particle_utils,        only: global_count_all_particles
#endif /* NBODY */
#ifdef MULTIGRID
      use multigrid,             only: init_multigrid, init_multigrid_ext, multigrid_par
#endif /* MULTIGRID */
#ifdef DEBUG
      use piernikdebug,          only: init_piernikdebug
      use piernikiodebug,        only: init_piernikiodebug
#endif /* DEBUG */
#ifdef COSM_RAYS
      use crdiffusion,           only: init_crdiffusion
#endif /* COSM_RAYS */
#ifdef CRESP
      use cresp_grid,            only: cresp_init_grid
#endif /* CRESP */
#ifdef RANDOMIZE
      use randomization,         only: init_randomization
#endif /* RANDOMIZE */
#ifdef PIERNIK_OPENCL
      use piernikcl,             only: init_opencl
#endif /* PIERNIK_OPENCL */
#if defined(__INTEL_COMPILER)
      !> \deprecated remove this clause as soon as Intel Compiler gets required features and/or bug fixes
      use timestep,              only: init_time_step
#ifdef COSM_RAYS
      use crhelpers,             only: init_div_v
#endif /* COSM_RAYS */
#endif /* __INTEL_COMPILER */

      implicit none

      integer :: nit, ac
      real    :: ts                  !< Timestep wallclock
      logical :: finished
      integer, parameter :: nit_over = 3 ! maximum number of auxiliary iterations after reaching level_max
      character(len=*), parameter :: ip_label = "init_piernik", ic_label = "IC_piernik", iter_label = "IC_iteration ", prob_label = "problem_IC"
      character(len=cbuff_len) :: label
#if defined(GRAV) && defined(NBODY)
      integer(kind=8) :: np
#endif /* GRAV && NBODY */

      call set_colors(.false.)               ! Make sure that we won't emit colorful messages before we are allowed to do so

      call parse_cmdline
      write(par_file,'(2a)') trim(wd_rd),'problem.par'
      write(tmp_log_file,'(2a)') trim(log_wr),'tmp.log'

      call init_wstats
      call init_mpi                          ! First, we must initialize the communication (and things that do not depend on init_mpi if there are any)
      code_progress = PIERNIK_INIT_MPI       ! Now we can initialize grid and everything that depends at most on init_mpi. All calls prior to PIERNIK_INIT_GRID can be reshuffled when necessary

      ! Timers should not be started before initializing MPI
      ts=set_timer(tmr_fu,.true.)

      call check_environment

#ifdef PIERNIK_OPENCL
      call init_opencl
#endif /* PIERNIK_OPENCL */

#ifdef DEBUG
      call init_piernikdebug                 ! Make it available as early as possible - right after init_mpi
      call init_piernikiodebug
#endif /* DEBUG */

      call cg_extptrs%epa_init

      call init_dataio_parameters            ! Required very early to call colormessage without side-effects
      ! Colors in stdout will be available from here if /OUTPUT_CONTROL/ colormode == .true.

      if (piernik_verbosity <= V_DEBUG) then
         do nit= V_LOWEST, V_HIGHEST
            write(msg, '(a,i1,a)')"Verbosity level ", nit, " (" // trim(v_name(nit)) // ")"
            if (master) call printinfo(msg, int(nit, kind=4))
         enddo
      endif

      call pnames%init ; call print_hostnames
      call init_load_balance
      call init_memory
      call init_profiling                    ! May require init_dataio_parameters and memory_usage set up
      call ppp_main%set_bb(ip_label)         ! can't call tst_cnt%start("init_piernik") before init_mpi

      call init_units

#ifdef RANDOMIZE
      call init_randomization
#endif /* RANDOMIZE */

      call init_default_fluidboundaries

      call problem_pointers                  ! set up problem-specific pointers as early as possible to allow implementation of problem-specific hacks also during the initialization
      call init_domain
      code_progress = PIERNIK_INIT_DOMAIN    ! Base domain is known and initial domain decomposition is known
      call init_geometry                     ! depends on domain

      call init_global
      code_progress = PIERNIK_INIT_GLOBAL    ! Global parameters are set up

      call set_interpolations

      call init_fluids
      code_progress = PIERNIK_INIT_FLUIDS    ! Fluid properties are set up

      call all_cg%init
      call all_cg%register_fluids            ! Register named fields for u, b and wa, depends on fluids and domain

#ifdef COSM_RAYS
#if defined(__INTEL_COMPILER)
      !> \deprecated remove this clause as soon as Intel Compiler gets required features and/or bug fixes
      call init_div_v
#endif /* __INTEL_COMPILER */
      call init_crdiffusion                  ! depends on fluids
#endif /* COSM_RAYS */

      call init_refinement
      call urc_list%init                     ! initialize unified refinement criteria

      call init_decomposition
#ifdef GRAV
      call init_grav                         ! Has to be called before init_grid
      call init_hydrostatic
#endif /* GRAV */
#ifdef NBODY
      call init_particles
      call init_psolver
#endif /* NBODY */
#ifdef MULTIGRID
      call init_multigrid_ext                ! Has to be called before init_grid
      call multigrid_par
#endif /* MULTIGRID */

      call init_grid                         ! Most of the cg's vars are now initialized, only arrays left
      code_progress = PIERNIK_INIT_GRID      ! Now we can initialize things that depend on all the above fundamental calls

#ifdef RESISTIVE
      call init_resistivity                  ! depends on grid
#endif /* RESISTIVE */

      call init_sources                      ! depends on: geometry, fluids, grid

#ifdef MULTIGRID
      call init_multigrid                    ! depends on grid, geometry, units and arrays
#endif /* MULTIGRID */

#if defined(__INTEL_COMPILER)
      !> \deprecated remove this clause as soon as Intel Compiler gets required features and/or bug fixes
      call init_time_step
#endif /* __INTEL_COMPILER */

      call init_fargo

      code_progress = PIERNIK_INIT_IO_IC     ! Almost everything is initialized: do problem-related stuff here, set-up I/O and create or read the initial conditions.

      call read_problem_par                  ! may depend on anything but init_dataio, \todo add checks against PIERNIK_INIT_IO_IC to all initproblem::read_problem_par

      call init_dataio                       ! depends on units, fluids (through common_hdf5), fluidboundaries, arrays, grid and shear (through magboundaries::bnd_b or fluidboundaries::bnd_u) \todo split me
      ! Initial conditions are read here from a restart file if possible

#ifdef CRESP
      call cresp_init_grid                   ! depends on cg
#endif /* CRESP */

#ifdef GRAV
      if (restarted_sim) call source_terms_grav
      call init_terms_grav
#endif /* GRAV */

      if (restarted_sim) then
         call all_bnd
         call all_bnd_vital_q
      endif

      call ppp_main%start(ic_label)
      if (master) then
         call printinfo("###############     Initial Conditions     ###############", V_LOG)
         write(msg,'(4a)') "   Starting problem : ",trim(problem_name)," :: ",trim(run_id)
         call printinfo("", V_ESSENTIAL, .true.)
         call printinfo(msg, V_ESSENTIAL, .true.)
         call printinfo("", V_ESSENTIAL, .true.)
      endif
      !>
      !! \deprecated It makes no sense to call (sometimes expensive) problem_initial_conditions before reading restart file.
      !! BEWARE: If your problem requires to call problem_initial_conditions add "require_problem_IC = 1" to restart file
      !! Move everything that is not regenerated by restart file to read_problem_par or create separate post-restart initialization
      !<
      !> \warning Set initial conditions by hand when starting from scratch or read them from a restart file. Do not use both unless you REALLY need to do so.
      if (restarted_sim .and. require_problem_IC /= 1) then
         if (master) then
            write(msg,'(a,i4,a)') "[initpiernik] Restart file #",nrestart," read. Skipping problem_initial_conditions."
            call printio(msg)
         endif
         if (associated(problem_post_restart)) then
            if (master) call printinfo("[initpiernik] Calling problem specific, post restart procedure", V_INFO)
            call problem_post_restart
         endif
      else

         call ppp_main%start(iter_label // "0", PPP_PROB)
         nit = 0
         finished = .false.

         call ppp_main%start(prob_label)
         call problem_initial_conditions ! may depend on anything
         call ppp_main%stop(prob_label)

         call init_psi ! initialize the auxiliary field for divergence cleaning when needed

         write(msg, '(a,f10.2)')"[initpiernik] IC on base level, time elapsed: ",set_timer(tmr_fu)
         if (master) call printinfo(msg, V_INFO)
         call ppp_main%stop(iter_label // "0", PPP_PROB)

         call costs_maintenance

         do while (.not. finished)
            write(label, '(i8)') nit + 1
            call ppp_main%start(iter_label // adjustl(label), PPP_PROB)

            call all_bnd !> \warning Never assume that problem_initial_conditions set guardcells correctly
#ifdef NBODY
            call update_particle_gravpot_and_acc  ! calls source_terms_grav
            call update_particle_kinetic_energy
#else /* !NBODY */
#ifdef GRAV
            call source_terms_grav
#endif /* GRAV */
#endif /* !NBODY */

            call update_refinement(act_count=ac)
            finished = (ac == 0) .or. (nit > level_max + nit_over) ! level_max iterations for creating refinement levels + level_max iterations for derefining excess of blocks

            call ppp_main%start(prob_label)
            call problem_initial_conditions ! reset initial conditions after possible changes of refinement structure
            call ppp_main%stop(prob_label)

            nit = nit + 1
            write(msg, '(2(a,i3),a,f10.2)')"[initpiernik] IC iteration: ",nit,", finest level:",finest%level%l%id,", time elapsed: ",set_timer(tmr_fu)
            if (master) call printinfo(msg, V_INFO)
            call ppp_main%stop(iter_label // adjustl(label), PPP_PROB)
            call costs_maintenance
         enddo
#ifdef GRAV
         call cleanup_hydrostatic
#endif /* GRAV */

         if (ac /= 0) then
            if (master) call warn("[initpiernik] The refinement structure does not seem to converge. Your refinement criteria may lead to oscillations of refinement structure. Bailing out.")
#ifdef GRAV
            call source_terms_grav  ! fix up gravitational potential when refiements did not converge
#endif /* GRAV */
         endif
#if defined(SELF_GRAV) && defined(NBODY)
         !  Do we need to do anything particle-related to be called here?
#endif /* SELF_GRAV && NBODY */
         if (associated(problem_post_IC)) call problem_post_IC
      endif
      call ppp_main%stop(ic_label)

      code_progress = PIERNIK_POST_IC

      write(msg, '(a,3i11,a,i2)')"[initpiernik] Effective resolution is [", finest%level%l%n_d(:), " ] at level ", finest%level%l%id
      !> \todo Do an MPI_Reduce in case the master process don't have any part of the globally finest level or ensure it is empty in such case
      if (master) call printinfo(msg, V_INFO)

#if defined(GRAV) && defined(NBODY)
      np = global_count_all_particles(message = "[initpiernik]")  ! we don't really need np, just call printinfo()
#endif /* GRAV && NBODY */

#ifdef VERBOSE
      call diagnose_arrays                   ! may depend on everything
#endif /* VERBOSE */
      call costs_maintenance

      call write_data(output=INCEPTIVE)

      call sanitize_smallx_checks            ! depends on problem_initial_conditions || init_dataio/read_restart_hdf5

      call ppp_main%stop(ip_label)

   contains

      !>
      !! \brief Print the list of hostnames in use and associated MPI ranks
      !!
      !! It was intentionally moved outside procnames module because of dataio_pub dependencies.
      !<

      subroutine print_hostnames

         use constants,  only: fmt_len, V_VERBOSE
         use dataio_pub, only: msg, printinfo
         use mpisetup,   only: slave, nproc
         use procnames,  only: pnames

         implicit none

         integer :: h, hl
         integer, parameter :: mpl = 16  ! maximum ranks per line to be printed (in non-consecutive case)
         character(len=*), parameter :: rah_o = "Ranks at host '", rah_c = "' :"
         character(len=fmt_len) :: fmtl, fmtr, fmt1, header
         logical :: successive, succ

         if (slave) return

         associate (intlen => int(log10(real(max(1, nproc-1)))) + 2)

            write(fmtl, *)"(a,", mpl, "i", intlen, ")"
            write(fmtr, *)"(a,i", intlen, ",' ..',i", intlen, ")"
            write(fmt1, *)"(a,i", intlen, ")"

         end associate

         successive = .true.
         do h = lbound(pnames%proc_on_node, 1), ubound(pnames%proc_on_node, 1)
            associate (p => pnames%proc_on_node(h))

               header = rah_o // p%nodename(:pnames%maxnamelen) // rah_c

               succ = .true.
               if (size(p%proc) > 1) succ = all(p%proc(:ubound(p%proc, 1)-1) + 1 == p%proc(lbound(p%proc, 1)+1:))
               successive = successive .and. succ

               if (succ) then

                  if (size(p%proc) > 1) then
                     write(msg, fmtr) trim(header), p%proc(lbound(p%proc, 1)), p%proc(ubound(p%proc, 1))
                  else
                     write(msg, fmt1) trim(header), p%proc(lbound(p%proc, 1))
                  endif
                  call printinfo(msg, V_VERBOSE)

               else

                  do hl = 0, int((ubound(p%proc, 1) - 1)/ mpl)
                     write(msg, fmtl) merge(repeat(" ", len_trim(header)), trim(header), hl>0), p%proc(hl*mpl+1:min((hl+1)*mpl, ubound(p%proc, 1)))
                     call printinfo(msg, V_VERBOSE)
                  enddo

               endif

            end associate
         enddo

         ! The load balance and cg distribution routines are designed to use
         ! space-filling curve that is placing spatially adjacent cgs closely
         ! on the block list in most cases.
         ! The use of mpirun options like "-map-by node" may defeat these efforts
         ! and significantly increase the amount of inter-node communication.
         if (.not. successive) call warn("[initpiernik:init_piernik] Non-successive MPI ranks on hosts detected. This may severely degrade the performance.")

         if (any(pnames%hostindex < 0)) call die("[initpiernik:init_piernik] pnames%hostindex contains invalid data")

      end subroutine print_hostnames

   end subroutine init_piernik
!-----------------------------------------------------------------------------
   subroutine parse_cmdline

      use constants,  only: stdout, cwdlen
      use dataio_pub, only: cmdl_nml, wd_rd, wd_wr, piernik_hdf5_version, piernik_hdf5_version2, log_wr
      use version,    only: nenv, env, init_version

      implicit none

      integer :: i, j
      logical :: skip_next
      character(len=8)            :: date   ! QA_WARN len defined by ISO standard
      character(len=10)           :: time   ! QA_WARN len defined by ISO standard
      character(len=5)            :: zone   ! QA_WARN len defined by ISO standard
      character(len=cwdlen)       :: arg
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
            call init_version
            write(stdout, '(a,f5.2)') 'GDF output version: ',piernik_hdf5_version2
            write(stdout, '(a,f5.2)') 'old output version: ',piernik_hdf5_version
            write(stdout,'(/,a)') "###############     Source configuration     ###############"
            do j = 1, nenv
               write(stdout,'(a)') env(j)
            enddo
            stop
         case ('-p', '--param')
            write(wd_rd,'(a)') get_next_arg(i+1, arg)
            skip_next = .true.
         case ('-w', '--write')
            write(wd_wr,'(a)') get_next_arg(i+1, arg)
            skip_next = .true.
         case ('-l', '--log')
            write(log_wr,'(a)') get_next_arg(i+1, arg)
            skip_next = .true.
         case ('-n', '--namelist')
            write(cmdl_nml(len_trim(cmdl_nml)+1:), '(2A)') " ", trim(get_next_arg(i+1, arg))
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

      if (wd_wr(len_trim(wd_wr):len_trim(wd_wr)) /= '/' ) write(wd_wr(len_trim(wd_wr)+1:),'(a1)') '/'
      if (wd_rd(len_trim(wd_rd):len_trim(wd_rd)) /= '/' ) write(wd_rd(len_trim(wd_rd)+1:),'(a1)') '/'

      ! Print the date and, optionally, the time
      call date_and_time(DATE=date, TIME=time, ZONE=zone)
      if (do_time) then
         write (stdout, '(a,"-",a,"-",a)', advance='no') date(1:4), date(5:6), date(7:8)
         write (stdout, '(1x,a,":",a,1x,a)') time(1:2), time(3:4), zone
         stop
      endif

   contains

      function get_next_arg(n, arg) result(param)

         use constants, only: stderr

         implicit none

         integer,               intent(in) :: n
         character(len=cwdlen), intent(in) :: arg

         character(len=cwdlen) :: param

         if (n > command_argument_count()) then
            write(stderr, '(2a)')"[initpiernik:parse_cmdline:get_next_arg] cannot find argument for option ", arg
            stop
         endif

         call get_command_argument(n, param)

      end function get_next_arg


   end subroutine parse_cmdline
!-----------------------------------------------------------------------------
   subroutine print_help

      use constants, only: cwdlen

      implicit none

      character(len=cwdlen) :: calledname

      call get_command_argument(0, calledname)

      print '(3a)','usage: ',trim(calledname),' [OPTIONS]'
      print '(a)', ''
      print '(a)', 'Recognized options:'
      print '(a)', ''
      print '(a)', '  -v, --version     print version information and exit'
      print '(6a)','  -n, --namelist    read namelist from command line, e.g. ',trim(calledname)," -n '",achar(38),'OUTPUT_CONTROL run_id="t22" /',"'"
      print '(a)', '  -p, --param       path to the directory with problem.par and/or restarts (default: ".")'
      print '(a)', '  -w, --write       path to the directory where output will be written (default: ".")'
      print '(a)', '  -l, --log         path to the directory where logs will be written (default: ".")'
      print '(a)', '  -h, --help        print usage information and exit'
      print '(a)', '  -t, --time        print time and exit'

   end subroutine print_help
!-----------------------------------------------------------------------------
end module initpiernik
