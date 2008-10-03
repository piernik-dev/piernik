#include "piernik.def"

module init_problem
  
! Initial condition for fluid flows for Kelvin-Helmholtz Instability
! Written by: D. Woltanski, February 2008

  use mpi_setup
  
  character problem_name*32,run_id*3
  real chi,dbot,lpert,Mtop,Mbot,dpert,tkh,vtransf

  namelist /PROBLEM_CONTROL/  problem_name, run_id, &
                              chi,dbot,lpert,Mtop,Mbot,dpert,tkh,vtransf


contains

!-----------------------------------------------------------------------------

  subroutine read_problem_par
    implicit none
  
    
    character par_file*(100), tmp_log_file*(100)
    integer :: cwd_status 

    par_file = trim(cwd)//'/problem.par'
    tmp_log_file = trim(cwd)//'/tmp.log'    


    problem_name = 'aaa'
    run_id  = 'aa'
    chi     = 8.0
    dbot    = 1.0
    lpert   = 0.05
    Mtop    = 0.11
    Mbot    = 0.34
    dpert   = 80.0
    tkh     = 1.70
    vtransf = 0.0
         
    
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
      rbuff(2) = dbot
      rbuff(3) = lpert
      rbuff(4) = Mtop
      rbuff(5) = Mbot
      rbuff(6) = dpert
      rbuff(7) = tkh
      rbuff(8) = vtransf

    
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
      dbot         = rbuff(2)  
      lpert	   = rbuff(3)
      Mtop	   = rbuff(4)
      Mbot	   = rbuff(5)
      dpert	   = rbuff(6)
      tkh	   = rbuff(7)
      vtransf      = rbuff(8)
    
    endif

  end subroutine read_problem_par

!-----------------------------------------------------------------------------

  subroutine init_prob
    use arrays, only    :   u,x,y,nx,ny
    use start, only     :   ymin,ymax,gamma,dimensions
    implicit none
    
    real dtop,lambda,boxlen,p0,vtop,vbot,k0,vp,rcx,rcy,rc
    integer i,j
 
dtop = dbot/chi
lambda = 1./6.
gamma = 5./3.
boxlen = ymax-ymin

       p0=lambda**2*(1.+chi)**2/(chi*tkh**2)*dbot/((Mtop*sqrt(chi)+Mbot)**2*gamma)
       vtop  =  1.*Mtop*sqrt(gamma*p0/dtop)
       vbot  = -1.*Mbot*sqrt(gamma*p0/dbot)
       k0    = 2.*3.141592652/lambda
       vp    = (Mtop*sqrt(chi)+Mbot)*sqrt(gamma*p0/dbot)/dpert

    do i = 1,nx
      rcx = x(i)
      do j = 1,ny
        rcy = y(j)
	rc=rcy-0.5*boxlen
	if(rc .gt. 0.0) then
	  u(1,i,j,:) = dtop
	  u(2,i,j,:) = vtop*dtop
	endif
	if(rc .le. 0.0) then
	  u(1,i,j,:) = dbot
	  u(2,i,j,:) = vbot*dbot
	endif
	if(abs(rc) .lt. lpert) then
	  u(3,i,j,:) = vp*sin(k0*rcx)*u(1,i,j,:)
	endif
	if(dimensions .eq. '3d') then
	  u(4,i,j,:) = vtransf*u(1,i,j,:)
	endif
          u(5,i,j,:) = p0/(gamma-1.0)
      enddo
    enddo


    return
  end subroutine init_prob  

!-----------------------------------------------------------------------------------------------------------------------------------

end module init_problem

