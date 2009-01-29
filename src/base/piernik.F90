! $Id$
#include "piernik.def"

program piernik

! Written by: M. Hanasz, January 2006
! Uses:       Modified Pen, Arras & Wong (2003) algorithm + original
!             source code "mhd.f90" from the webpage:
!             http://www.cita.utoronto.ca/~pen/MHD
!             within the "mhdblock.f90"
! Modification history: see "changelog"

  use start, only  : t,dt,nstep, nend, nstep_start, tend
  use timer, only : timer_start, timer_stop
  use dataio, only : write_data, user_msg_handler
  use mhdstep, only  : fluid_update
  use mpisetup
    
  implicit none
  logical :: end_sim
 
  call init_piernik

  call MPI_BARRIER(comm3d,ierr)
!-------------------------------- MAIN LOOP ----------------------------------
  call timer_start

  end_sim = .false.
  do while(t <= tend .and. nstep < nend .and. .not.(end_sim) )
    nstep=nstep+1

      call fluid_update

      call MPI_BARRIER(comm3d,ierr)
      call write_data(output='all')

      call user_msg_handler(end_sim)

  end do ! main loop

  nstep=nstep-1
  call timer_stop

  call write_data(output='end')
!---------------------------- END OF MAIN LOOP ----------------------------------

  call MPI_BARRIER(comm,ierr)
  call cleanup_piernik

  call mpistop

contains

   subroutine init_piernik
      use start, only  : read_params
      use initfluids, only : init_fluids
      use fluidindex, only : fluid_index,nvar
      use arrays, only : arrays_allocate
      use grid, only : nx,ny,nz
      use grid, only : init_grid,grid_xyz
      use initproblem, only : init_prob, read_problem_par
      use dataio, only : init_dataio, write_data, user_msg_handler
      use mpisetup
      use mpiboundaries
#if defined MAGNETIC && defined RESISTIVE
      use resistivity, only : init_resistivity
#endif /* MAGNETIC && RESISTIVE */
#ifdef SHEAR
      use shear, only : init_shear
#endif /* SHEAR */
#ifdef SELF_GRAV
      use poissonsolver
#endif /* SELF_GRAV */
#ifdef GRAV
      use gravity, only : init_grav,grav_pot_3d
#endif /* GRAV */
#ifdef SNE_DISTR
      use sndistr, only  : init_sndistr
#endif /* SNE_DISTR */
      implicit none
      call getarg(1, cwd)
      if (LEN_TRIM(cwd) == 0) cwd = '.'
 
      call mpistart

      call read_params
 
      call init_grid
  
      call init_fluids
  
      call read_problem_par
  
      call fluid_index

      call arrays_allocate(nx,ny,nz,nvar)

      call grid_xyz

      call mpi_boundaries_prep

#ifdef RESISTIVE
      call init_resistivity
#endif /* RESISTIVE */

#ifdef SHEAR
      call init_shear
#endif /* SHEAR */

#ifdef GRAV
      call init_grav
      call grav_pot_3d
#endif /* GRAV */

#ifdef SNE_DISTR
      call init_sndistr
#endif /* SNE_DISTR */

      call init_prob

      call init_dataio

   end subroutine init_piernik

   subroutine cleanup_piernik
      use grid, only : cleanup_grid
      use arrays, only : arrays_deallocate
#ifdef RESISTIVE
      use resistivity, only : cleanup_resistivity
#endif RESISTIVE

      call cleanup_grid
#ifdef RESISTIVE
      call cleanup_resistivity
#endif /* RESISTIVE */
      call arrays_deallocate
   end subroutine cleanup_piernik
end program piernik

