! $Id$
#include "piernik.def"

program piernik

! Written by: M. Hanasz, January 2006
! Uses:       Modified Pen, Arras & Wong (2003) algorithm + original
!             source code "mhd.f90" from the webpage:
!             http://www.cita.utoronto.ca/~pen/MHD
!             within the "mhdblock.f90"
! Modification history: see "changelog"

  use start, only  : read_params
  use start, only  : t,dt,nstep, nend, nstep_start, tend

  use initfluids, only : init_fluids
  use fluidindex, only : fluid_index,nvar

  use arrays, only : arrays_deallocate, arrays_allocate
  use grid, only : nx,ny,nz
  use grid, only : init_grid,grid_xyz,cleanup_grid
  use init_problem, only : init_prob, read_problem_par
  use mhdstep, only  : mhd_step
  use dataio, only : init_dataio, write_data, user_msg_handler
  use timer, only : timer_start, timer_stop
  use mpi_setup
  use mpi_bnd
  use fluid_boundaries, only : compute_u_bnd
  
#ifdef MAGNETIC  
  use mag_boundaries, only   : compute_b_bnd
#ifdef RESISTIVE
  use resistivity, only : init_resistivity, cleanup_resistivity
#endif /* RESISTIVE */
#endif /* MAGNETIC */
#ifdef SHEAR
  use shear, only : init_shear
#endif /* SHEAR */
#ifdef SELF_GRAV
  use poisson_solver
#endif /* SELF_GRAV */
#ifdef GRAV
  use gravity, only : init_grav,grav_pot_3d
#endif /* GRAV */
#ifdef SNE_DISTR
    use sndistr, only  : init_sndistr
#endif /* SNE_DISTR */
    
  implicit none
  logical :: end_sim

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

  call mpi_bnd_prep

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

  call MPI_BARRIER(comm3d,ierr)
!-------------------------------- MAIN LOOP ----------------------------------
  call timer_start

  end_sim = .false.
  do
    nstep=nstep+1
    if (t>=tend .or. nstep>nend .or. end_sim ) exit

      call mhd_step

      call MPI_BARRIER(comm3d,ierr)
      call write_data(output='all')

      call user_msg_handler(end_sim)

   end do ! main loop

  nstep=nstep-1
  call timer_stop

  call write_data(output='end')
!---------------------------- END OF MAIN LOOP ----------------------------------


  call MPI_BARRIER(comm,ierr)

  call cleanup_grid
#ifdef RESISTIVE
  call cleanup_resistivity
#endif /* RESISTIVE */
  call arrays_deallocate

  call MPI_BARRIER(comm,ierr)
  call mpistop

end program piernik

