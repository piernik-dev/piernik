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
!    Initial implemetation of PIERNIK code was based on TVD split MHD code by
!    Ue-Li Pen
!        see: Pen, Arras & Wong (2003) for algorithm and
!             http://www.cita.utoronto.ca/~pen/MHD
!             for original source code "mhd.f90"
!
!    For full list of developers see $PIERNIK_HOME/license/pdt.txt
!
#include "piernik.def"
!>
!! Main program
!<
program piernik

  use mpisetup,      only : comm, comm3d, ierr, proc, t, dt, nstep, mpistop
  use dataio_public, only : nend, nstep_start, tend, log_file, log_lun, &
       &                    code_progress, PIERNIK_START, PIERNIK_INITIALIZED, PIERNIK_FINISHED, PIERNIK_CLEANUP
  use timer,         only : time_left
#ifdef PERFMON
  use timer,         only : timer_start, timer_stop
#endif
  use dataio,        only : write_data, user_msg_handler
  use fluidupdate,   only : fluid_update
  use types,         only : finalize_problem

  implicit none

  logical            :: end_sim !< Used in main loop, to test whether to stop simulation or not
  character(len=256) :: msg
  character(len=32)  :: nstr, tstr

  code_progress = PIERNIK_START

  call init_piernik

  call MPI_BARRIER(comm3d,ierr)
!-------------------------------- MAIN LOOP ----------------------------------
#ifdef PERFMON
  call timer_start
#endif

  code_progress = PIERNIK_INITIALIZED

  end_sim = .false.

  if (proc == 0) then
     write(*, '("======================================================================================================")')
     open(log_lun, file=log_file, position='append')
     write(log_lun, '(/,a,/)') "###############     Simulation     ###############"
     close(log_lun)
  end if

  do while(t < tend .and. nstep < nend .and. .not.(end_sim) .and. time_left() )

     nstep=nstep+1

     if (proc == 0) then
        open(log_lun, file=log_file, position='append')
        write(log_lun, '("   nstep = ",i7,"   dt = ",es22.16,"   t = ",es22.16)') nstep, dt, t
        close(log_lun)
     end if

     call fluid_update

     call MPI_BARRIER(comm3d,ierr)
     call write_data(output='all')

     call user_msg_handler(end_sim)

  end do ! main loop

  code_progress = PIERNIK_FINISHED

  if (proc == 0) then
     write(tstr, '(g14.6)') t
     write(nstr, '(i7)') nstep
     tstr = adjustl(tstr)
     nstr = adjustl(nstr)
     write(*, '("======================================================================================================")')
     open(log_lun, file=log_file, position='append')
     write(log_lun, '(/,a,/)') "###############     Finishing     ###############"
     if (t >= tend) then
        write(msg, '(2a)') "Simulation has reached final time t = ",trim(tstr)
        write(*, '(a)') trim(msg)
        write(log_lun, '(a)') trim(msg)
     end if
     if (nstep >= nend) then
        write(msg, '(4a)') "Maximum step count exceeded (",trim(nstr),") at t = ",trim(tstr)
        write(*, '(a)') trim(msg)
        write(log_lun, '(a)') trim(msg)
     end if
     if (end_sim) then
        write(msg, '(4a)') "Enforced stop at step ",trim(nstr),", t = ", trim(tstr)
        write(*, '(a)') trim(msg)
        write(log_lun, '(a)') trim(msg)
     end if
     if(.not.time_left()) then
        write(msg, '(4a)') "Wall time limit exceeded at step ",trim(nstr),", t = ", trim(tstr)
        write(*, '(a)') trim(msg)
        write(log_lun, '(a)') trim(msg)
     end if
     close(log_lun)
  end if

  if(associated(finalize_problem)) call finalize_problem

!  nstep=nstep+1
#ifdef PERFMON
  call timer_stop
#endif /* PERFMON */
  call write_data(output='end')
!---------------------------- END OF MAIN LOOP ----------------------------------

  call MPI_BARRIER(comm,ierr)

  code_progress = PIERNIK_CLEANUP

  if (proc == 0) write(*, '(a)', advance='no') "Finishing "
  call cleanup_piernik
  if (proc == 0) write(*, '("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")')

  call mpistop

contains
!>
!! Meta subroutine responsible for initializing all piernik modules
!<
   subroutine init_piernik
      use types,         only : grid_container
      use initfluids,    only : init_fluids
      use fluidindex,    only : nvar
      use arrays,        only : init_arrays
      use grid,          only : nx,ny,nz
      use grid,          only : init_grid,grid_xyz
      use initproblem,   only : init_prob, read_problem_par
      use problem_pub,   only : problem_name, run_id
      use dataio,        only : init_dataio, write_data
      use dataio_public, only : nrestart
      use mpisetup,      only : cwd, mpistart
      use mpiboundaries, only : mpi_boundaries_prep
      use fluidboundaries, only : all_fluid_boundaries
#ifdef MAGNETIC
      use magboundaries, only : all_mag_boundaries
#endif /* MAGNETIC */
#if defined MAGNETIC && defined RESISTIVE
      use resistivity,   only : init_resistivity
#endif /* MAGNETIC && RESISTIVE */
#ifdef SHEAR
      use shear,         only : init_shear
#endif /* SHEAR */
#ifdef GRAV
      use gravity,       only : init_grav, grav_pot_3d
#endif /* GRAV */
#ifdef FLUID_INTERACTIONS
      use interactions,  only : init_interactions
#endif /* FLUID_INTERACTIONS */
#ifdef SNE_DISTR
      use sndistr,       only : init_sndistr
#endif /* SNE_DISTR */
#ifdef MULTIGRID
      use multigrid,     only : init_multigrid
#endif /* MULTIGRID */

      implicit none

      type(grid_container) :: cgrid

      call getarg(1, cwd)
      if (LEN_TRIM(cwd) == 0) cwd = '.'

      call mpistart

      call init_grid(cgrid)

#ifdef SHEAR
      call init_shear
#endif /* SHEAR */

      call read_problem_par
      if (proc == 0) write(*,'(/,2a,/)')"   Starting problem : ",trim(problem_name)

      call init_fluids

      call init_arrays(nx,ny,nz,nvar)

      call mpi_boundaries_prep

#ifdef RESISTIVE
      call init_resistivity
#endif /* RESISTIVE */

#ifdef GRAV
      call init_grav
! It is only temporary solution, but grav_pot_3d must be called after init_prob due to csim2,c_si,alpha clash!!!
      call grav_pot_3d
#endif /* GRAV */

#ifdef FLUID_INTERACTIONS
      call init_interactions
#endif /* FLUID_INTERACTIONS */

#ifdef SNE_DISTR
      call init_sndistr
#endif /* SNE_DISTR */
#ifdef MULTIGRID
      call init_multigrid(cgrid)
#endif /* MULTIGRID */

      call init_dataio

      ! It makes no sense to call (sometimes expensive) init_prob before reading restart file.
      ! BEWARE: Current change may break some problems that depend on other things set in init_prob.
      ! Move everything that is not regenerated by restart file to read_problem_par or create separate post-restart initialization
      if (nrestart>0) then
         if (proc == 0) write(*,'(a,i4,a)')"[piernik:init_piernik] Restart file #",nrestart," read. Warning: skipping init_prob."
      else
         call init_prob
         call all_fluid_boundaries ! Never assume that init_prob set guardcells correctly
#ifdef MAGNETIC
         call all_mag_boundaries
#endif /* MAGNETIC */
         call write_data(output='all') ! moved from dataio::init_dataio
      end if

   end subroutine init_piernik

!>
!! Meta subroutine responsible for cleaning after successful run
!<

   subroutine cleanup_piernik

      use grid,        only : cleanup_grid
      use dataio,      only : cleanup_dataio
      use arrays,      only : cleanup_arrays
      use initfluids,  only : cleanup_fluids
      use timer,       only : cleanup_timers
#ifdef RESISTIVE
      use resistivity, only : cleanup_resistivity
#endif /* RESISTIVE */
#ifdef MULTIGRID
      use multigrid,   only : cleanup_multigrid
#endif /* MULTIGRID */

      call cleanup_grid;        if (proc == 0) write(*,'(a)',advance='no')"."
      call cleanup_dataio;      if (proc == 0) write(*,'(a)',advance='no')"."
#ifdef RESISTIVE
      call cleanup_resistivity; if (proc == 0) write(*,'(a)',advance='no')"."
#endif /* RESISTIVE */
#ifdef MULTIGRID
      call cleanup_multigrid;   if (proc == 0) write(*,'(a)',advance='no')"."
#endif /* MULTIGRID */
      call cleanup_arrays;      if (proc == 0) write(*,'(a)',advance='no')"."
      call cleanup_fluids;      if (proc == 0) write(*,'(a)',advance='no')"."
      call cleanup_timers;      if (proc == 0) write(*,'(a)')"."

   end subroutine cleanup_piernik

end program piernik
