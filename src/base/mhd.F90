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
#endif 
  use mod_mhdstep
#ifdef SELF_GRAV
  use poisson_solver
#endif 
  use dataio
  use diagnostics  
  use timer
  use mpi_setup
  use fluid_boundaries
#ifdef GRAV
  use gravity, only : grav_pot_3d
#endif 
  use arrays

  implicit none
  character output*3
  integer system_status
  character system_command*160
  integer tsleep
!  real, dimension(:,:,:),allocatable :: ala
  
  call mpistart

  call timer_start
  
  call read_params
    
  call read_problem_par

  call arrays_allocate 
  
  call grid_xyz 
                
  call init_dataio

#ifdef GRAV
  call grav_pot_3d
#endif

  if(proc .eq. 0) then 
    if(restart .eq. 'last') then
      call find_last_restart(nrestart)
    endif
  endif
  call MPI_BCAST(nrestart, 1, MPI_INTEGER, 0, comm, ierr)
      

  if(nrestart .eq. 0) then

    nstep_start = 0
    t_start     = 0.0

    call init_prob

    call compute_u_bnd
    call compute_b_bnd

#ifdef MASS_COMPENS
    call save_init_dens
#endif        

#ifdef SELF_GRAV
    call poisson
#endif 

    if (proc .eq.0) then
      write (log_file,'(a,a1,a3,a1,i3.3,a4)') &
              trim(problem_name),'_', run_id,'_',nrestart,'.log'

      system_command = 'mv tmp.log '//log_file
      system_status = SYSTEM(system_command)
    endif

    call write_data(output='all')
   
  else  

    if (proc .eq. 0) then
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

#ifdef MASS_COMPENS
    call get_init_mass
#endif        

  endif
    
  
!-------------------------------- MAIN LOOP ----------------------------------
  do
    nstep=nstep+1
    if (t>=tend .or. nstep>nend ) exit

      call mhdstep

#ifdef MASS_COMPENS
      call mass_loss_compensate
#endif        
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
            
      if(wait .eqv. .true.) go to 888
    
   end do ! main loop

999 continue

  nstep=nstep-1

!  u(1,nb+1:nx-nb,nb+1:ny-nb,:) = & 
!     unshear_fft(u(1,nb+1:nx-nb,nb+1:ny-nb,:),x(nb+1:nx-nb),.true.)
!  allocate(ala(nxd,nyd,nz))
!  call unshear_fft_b(u(1,nb+1:nx-nb,nb+1:ny-nb,:),x(:),&
!       u(1,1:nb,nb+1:ny-nb,:),&
!       u(1,nxd+nb+1:nxd+2*nb,nb+1:ny-nb,:),&
!       ala(:,:,:),.true.)
!  u(1,nb+1:nx-nb,nb+1:ny-nb,:) = ala(:,:,:)
!  call write_hdf
!  nhdf = nhdf+1
!  u(1,nb+1:nx-nb,nb+1:ny-nb,:) = &
!     unshear_fft(u(1,nb+1:nx-nb,nb+1:ny-nb,:),x(nb+1:nx-nb))
!  call unshear_fft_b(u(1,nb+1:nx-nb,nb+1:ny-nb,:),x(:),&
!       u(1,1:nb,nb+1:ny-nb,:),&
!       u(1,nxd+nb+1:nxd+2*nb,nb+1:ny-nb,:),&
!       ala(:,:,:))
!  u(1,nb+1:nx-nb,nb+1:ny-nb,:) = ala(:,:,:)
!  call write_hdf
!  deallocate(ala)

!---------------------------- END OF MAIN LOOP ----------------------------------

 
   call arrays_deallocate 
   
   call timer_stop

   call mpistop 
!   write(*,*)

end program mhd

