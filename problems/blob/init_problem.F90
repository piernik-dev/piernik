#include "piernik.def"

module init_problem
  
! Initial condition for blob test
! Written by: D. Woltanski, March 2008

  use mpi_setup
  
  character problem_name*32,run_id*3
  real chi, rblob, blobxc, blobyc, blobzc, Mext, denv, tkh, vgal

  namelist /PROBLEM_CONTROL/  problem_name, run_id, &
                              chi, rblob, blobxc, blobyc, blobzc, Mext, denv, tkh, vgal


contains

!-----------------------------------------------------------------------------

  subroutine read_problem_par
    implicit none
  
    
    character par_file*(100), tmp_log_file*(100)

    par_file = trim(cwd)//'/problem.par'
    tmp_log_file = trim(cwd)//'/tmp.log'    


    problem_name = 'aaa'
    run_id  = 'aa'
    chi     = 10.0
    rblob   =  1.0
    blobxc  =  5.0
    blobyc  =  5.0
    blobzc  =  5.0
    Mext    =  2.7
    denv    =  1.0
    tkh     =  1.7
    vgal    =  0.0
    
    if(proc .eq. 0) then
      open(1,file=par_file)
        read(unit=1,nml=PROBLEM_CONTROL)
      close(1)
        write(*,nml=PROBLEM_CONTROL)
      open(3, file=tmp_log_file, position='append')
        write(3,nml=PROBLEM_CONTROL)
        write(3,*)
      close(3)
    endif

    if(proc .eq. 0) then


      cbuff(1) =  problem_name
      cbuff(2) =  run_id

      rbuff(1) = chi
      rbuff(2) = rblob
      rbuff(3) = blobxc
      rbuff(4) = blobyc
      rbuff(5) = blobzc
      rbuff(6) = Mext
      rbuff(7) = denv
      rbuff(8) = tkh
      rbuff(9) = vgal

    
      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

    else
    
      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)
      
      problem_name = cbuff(1)   
      run_id       = cbuff(2)   

      chi          = rbuff(1)  
      rblob        = rbuff(2)  
      blobxc	   = rbuff(3)
      blobyc	   = rbuff(4)
      blobzc	   = rbuff(5)
      Mext	   = rbuff(6)
      denv	   = rbuff(7)
      tkh          = rbuff(8)
      vgal	   = rbuff(9)
    
    endif

  end subroutine read_problem_par

!-----------------------------------------------------------------------------

  subroutine init_prob
    use arrays, only    :   u,x,y,z,nx,ny,nz
    use start, only     :   ymin,ymax,gamma,dimensions
    implicit none
    
    real penv, rcx, rcy, rcz, rrel
    integer i,j, k
 
    penv = 3.2*rblob*sqrt(chi)/tkh/(Mext*gamma/denv)
    
    do i = 1,nx
      rcx = x(i)
      do j = 1,ny
        rcy = y(j)
        do k = 1,nz
	  if(dimensions .eq. '3d') then
	    rcz = z(k)
	    rrel = sqrt((rcx-blobxc)**2+(rcy-blobyc)**2+(rcz-blobzc)**2)
	  else
	    rrel = sqrt((rcx-blobxc)**2+(rcy-blobyc)**2)
  	  endif
  	    if(rblob .ge. rrel) then
	      u(1,i,j,k) = chi*denv
	      u(2,i,j,k) = chi*denv*vgal
	      u(3,i,j,k) = 0.0
	      u(4,i,j,k) = 0.0
	    else
	      u(1,i,j,k) = denv
	      u(2,i,j,k) = denv*vgal
	      u(3,i,j,k) = Mext*gamma*penv
	      u(4,i,j,k) = 0.0
	    endif
            u(5,i,j,:) = penv/(gamma-1.0)
	  enddo
      enddo
    enddo
    if(proc .eq. 0) then
      open(5,file='cloudfrac.out',status='unknown')
      write(5,*) ' '
      close(5)
    endif

    return
  end subroutine init_prob  

!-----------------------------------------------------------------------------------------------------------------------------------

end module init_problem

