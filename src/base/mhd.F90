! $Id$
#include "piernik.def"

program mhd

! Written by: M. Hanasz, January 2006
! Uses:       Modified Pen, Arras & Wong (2003) algorithm + original 
!             source code "mhd.f90" from the webpage:
!             http://www.cita.utoronto.ca/~pen/MHD
!             within the "mhdblock.f90"
! Modification history: see "changelog"

  use start, only: read_params, new_id, restart, dt_hdf, dt_res,dt_log,dt_tsl
  use start, only: t,dt, tend, nstep, nend, nstep_start, nrestart
  use arrays, only : arrays_deallocate, arrays_allocate
  use grid, only : grid_xyz
  use init_problem, only : problem_name, run_id, init_prob, read_problem_par 
  use mod_mhdstep, only  : mhdstep
  use dataio, only : host, pid  
  use dataio, only : msg,step_res,step_hdf,log_file,nres,nhdf, &
      t_start, nres_start, nhdf_start, wait, msg_param
  use dataio, only : init_dataio,read_restart,write_data, check_disk, &
      read_file_msg, write_timeslice, write_log, write_hdf, write_restart, &
      find_last_restart, get_container, set_container
  use dataio_hdf5, only: read_restart_hdf5, write_plot
!  use diagnostics  
  use timer, only : timer_start, timer_stop
  use mpi_setup
  use mpi_bnd
  use fluid_boundaries, only : compute_u_bnd
  use mag_boundaries, only   : compute_b_bnd
#ifdef RESIST
  use resistivity
#endif /* RESIST */
#ifdef SELF_GRAV
  use poisson_solver
#endif /* SELF_GRAV */
#ifdef GRAV
  use gravity, only : grav_pot_3d
  use start, only : gpt_hdf
#endif /* GRAV */
#ifdef SNE_DISTR
    use sn_distr, only  : prepare_SNdistr
#endif /* SNE_DISTR */
#ifdef MASS_COMPENS
    use start,  only        :   init_mass, mass_loss, mass_loss_tot   
    use init_problem,  only :   save_init_dens, get_init_mass, mass_loss_compensate   
#endif  /* MASS_COMPENS */
    use types
    use comp_log, only : nenv,env

  implicit none
  character output*3
  integer system_status
  character system_command*160, cmd*256
  integer(kind=1)  :: system
  integer tsleep, i

  character tmp_log_file*(100)
  type(hdf) :: chdf

  call getarg(1, cwd)
  if (LEN_TRIM(cwd) == 0) cwd = '.'
  
  
  call mpistart

  call read_params
    
  call read_problem_par

  call arrays_allocate 
  
  call grid_xyz 

  call mpi_bnd_prep
                
  call init_dataio
#ifdef GRAV
  call grav_pot_3d
#endif /* GRAV */
  if(proc .eq. 0) then 
    if(restart .eq. 'last') then
      call find_last_restart(nrestart)
    endif
  endif
  call MPI_BARRIER(comm,ierr)
  call MPI_BCAST(nrestart, 1, MPI_INTEGER, 0, comm, ierr)

#ifdef SNE_DISTR
   call prepare_SNdistr
#endif /* SNE_DISTR */
      

  if(nrestart .eq. 0) then

    nstep_start = 0
    t_start     = 0.0

    call init_prob
#ifdef GRAV
!    if(gpt_hdf .eq. 'yes') call write_data('gpt')
#endif /* GRAV */

    call compute_u_bnd
    call compute_b_bnd

#ifdef MASS_COMPENS
    call save_init_dens
#endif /* MASS_COMPENS */

#ifdef SELF_GRAV
    call poisson
#endif /* SELF_GRAV */

    if (proc .eq.0) then
      write (log_file,'(a,a1,a3,a1,i3.3,a4)') &
              trim(problem_name),'_', run_id,'_',nrestart,'.log'
      tmp_log_file = trim(cwd)//'/tmp.log'
      log_file = trim(cwd)//'/'//trim(log_file)
      system_command = 'mv '//trim(tmp_log_file)//' '//trim(log_file)
      system_status = SYSTEM(system_command)
      open(3, file=log_file, position='append')
         do i=1,nenv
           write(3,*) trim(env(i))
         enddo
      close(3)
    endif

    call set_container(chdf); chdf%nres = nrestart
    call write_data(output='all')
   
  else  
    if (proc .eq. 0) then
      write (log_file,'(a,a1,a3,a1,i3.3,a4)') &
              trim(problem_name),'_', run_id,'_',nrestart,'.log'
      tmp_log_file = trim(cwd)//'/tmp.log'
      log_file = trim(cwd)//'/'//log_file
      system_command = 'mv '//trim(tmp_log_file)//' '//trim(log_file)
      system_status = SYSTEM(system_command)
    endif

#ifdef HDF5
    call set_container(chdf); chdf%nres = nrestart
    call read_restart_hdf5(chdf)
    call get_container(chdf); nstep = chdf%nstep
#else
    nres = nrestart
    call read_restart
#endif

    nstep_start = nstep
    t_start     = t
    nres_start  = nrestart
    nhdf_start  = nhdf-1

    if(new_id .ne. '') run_id=new_id

#ifdef MASS_COMPENS
      call get_init_mass
#endif /* MASS_COMPENS */     

  endif
  call MPI_BARRIER(comm3d,ierr)
!-------------------------------- MAIN LOOP ----------------------------------
  call timer_start
  
  do
    nstep=nstep+1
    if (t>=tend .or. nstep>nend ) exit
      call write_plot(chdf,"dens")

      call mhdstep

#ifdef MASS_COMPENS
      call mass_loss_compensate
#endif /* MASS_COMPENS */
      call MPI_BARRIER(comm3d,ierr)
      call write_data(output='all')

888   continue

!      call check_disk

!--- process 0 checks for messages

      if(proc .eq. 0) then
        call read_file_msg
      endif
      call MPI_BCAST(msg,       16, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(msg_param, 16, MPI_DOUBLE_PRECISION, 0, comm, ierr)
      
!---  if a user message is received then:
      if (len_trim(msg) .ne. 0) then
        if(trim(msg) .eq. 'res' .or. trim(msg) .eq. 'dump' ) then
          nres = max(nres,1)
          call write_restart
          nres = nres + 1
          step_res = nstep        
        endif  
        if(trim(msg) .eq. 'hdf') then
          call write_hdf        
          nhdf = nhdf + 1
          step_hdf = nstep        
        endif
        if(trim(msg) .eq. 'log')    call write_log      
        if(trim(msg) .eq. 'tsl')    call write_timeslice        
        if(trim(msg) .eq. 'tend')   tend   = msg_param  
        if(trim(msg) .eq. 'nend')   nend   = msg_param  
        if(trim(msg) .eq. 'dtres')  dt_res = msg_param  
        if(trim(msg) .eq. 'dthdf')  dt_hdf = msg_param  
        if(trim(msg) .eq. 'dtlog')  dt_log = msg_param  
        if(trim(msg) .eq. 'dttsl')  dt_tsl = msg_param  

        if(trim(msg) .eq. 'sleep') then
           tsleep = 60*msg_param
           call sleep(tsleep)   
        endif
 
        if(trim(msg) .eq. 'stop') goto 999
      endif
            
      if(wait .eqv. .true.) go to 888
    
   end do ! main loop

999 continue

  nstep=nstep-1
  call timer_stop

  call write_data(output='end')
!---------------------------- END OF MAIN LOOP ----------------------------------

 
  call MPI_BARRIER(comm,ierr)
  call arrays_deallocate 

  call MPI_BARRIER(comm,ierr)
  call mpistop 

end program mhd

