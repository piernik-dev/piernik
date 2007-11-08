#include "mhd.def"

program mhd

! Written by: M. Hanasz, January 2006
! Uses:       Modified Pen, Arras & Wong (2003) algorithm + original 
!             source code "mhd.f90" from the webpage:
!             http://www.cita.utoronto.ca/~pen/MHD
!             within the "mhdblock.f90"
! Modification history: see "changelog"

  use start, only: read_params
  use start, only: t,dt, tend, nstep, nend
  use grid
  use init_problem  
#ifdef RESIST
  use resistivity
#endif RESIST
  use mod_mhdstep
!  use poisson_solver
  use dataio
  use diagnostics  
  use timer
  use mpi_setup
#ifdef GRAV
  use gravity, only : grav_pot_3d
#endif 

  implicit none
  character output*3
  integer system_status
  character system_command*160
  real tsleep

  call mpistart

  call timer_start
  
  call read_params
    
  call arrays_allocate 
  
  call grid_xyz 
                
  call init_dataio
  
#ifdef GRAV
  call grav_pot_3d
#endif

  if(nrestart .eq. 0) then

    nstep_start = 0
    t_start     = 0.0

    call init_prob

    call compute_u_bnd
    call compute_b_bnd

!    if(selfgravity .ne. 'null') call poisson

    if (proc .eq.0) then
      write (log_file,'(a,a1,a3,a1,i3.3,a4)') &
              trim(problem_name),'_', run_id,'_',nrestart,'.log'

      system_command = 'mv tmp.log '//log_file
      system_status = SYSTEM(system_command)
    endif

    call write_data(output='all')
  else

    call read_problem_par

    if (proc .eq.0) then
      write (log_file,'(a,a1,a3,a1,i3.3,a4)') &
              trim(problem_name),'_', run_id,'_',nrestart,'.log'

      system_command = 'mv tmp.log '//log_file
      system_status = SYSTEM(system_command)
    endif

    nres = nrestart
    call read_restart
    nstep_start = nstep
    t_start     = t
    nres_start  = nrestart
    nhdf_start  = nhdf-1

    if(new_id .ne. '') run_id=new_id

  endif
  
!-------------------------------- MAIN LOOP ----------------------------------
  do
    nstep=nstep+1
    if (t>=tend .or. nstep>nend ) exit

      call mhdstep
      call write_data(output='all')

888   continue

      call check_disk

!--- process 0 checks for messages

      if(proc .eq. 0) then
        call read_file_msg
      endif
      call MPI_BCAST(msg,       16, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(msg_param, 16, MPI_DOUBLE_PRECISION, 0, comm, ierr)
      
!---  if a user message is received then:
      if (len_trim(msg) .ne. 0) then
!        write(*,*) proc, msg, msg_param
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
            
      if(wait .eq. .true.) go to 888
    
   end do ! main loop

999 continue

  nstep=nstep-1

  call write_data(output='end')
!---------------------------- END OF MAIN LOOP ----------------------------------
 
   call arrays_deallocate 
   
   call timer_stop

   call mpistop 
!   write(*,*)

end program mhd

