#include "piernik.def"


module init_problem
  
! Initial condition for Keplerian disk
! Written by: M. Hanasz, March 2006

  use arrays, only : x, u, b,x,y,z,nx,ny,nz, &
      idna,imxa,imya,imza
#ifndef ISO
  use arrays, only : iena
#endif /* ISO */
  use start,  only : proc, smallei, smalld, gamma, &
      rbuff, cbuff, ibuff
  use mpi_setup
  use shear, only: qshear, omega
!  use grid

  real :: d0,r0,bx0,by0,bz0
  character ::  problem_name*32,run_id*3,dir*1

  namelist /PROBLEM_CONTROL/  problem_name, run_id, d0, r0,bx0,by0,bz0


contains

!-----------------------------------------------------------------------------

  subroutine read_problem_par

    implicit none
  

    problem_name = 'slab'
    run_id  = 'tst'
    d0      = 1.0
    r0      = 0.25
    
    if(proc .eq. 0) then
      open(1,file='problem.par')
        read(unit=1,nml=PROBLEM_CONTROL)
        write(*,nml=PROBLEM_CONTROL)
      close(1)
      open(3, file='tmp.log', position='append')
        write(3,nml=PROBLEM_CONTROL)
        write(3,*)
      close(3)
    endif

    if(proc .eq. 0) then


      cbuff(1) =  problem_name
      cbuff(2) =  run_id

      rbuff(1) = d0
      rbuff(2) = r0
      rbuff(3) = bx0
      rbuff(4) = by0
      rbuff(5) = bz0
    
      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

    else
    
      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)
      
      problem_name = cbuff(1)   
      run_id       = cbuff(2)   

      d0           = rbuff(1)  
      r0           = rbuff(2)
      bx0          = rbuff(3)
      by0          = rbuff(4)
      bz0          = rbuff(5)
    
    endif

  end subroutine read_problem_par

!-----------------------------------------------------------------------------

  subroutine init_prob

    implicit none
 
    integer i,j,k
    real :: xi,yj,zk
    real :: vx,vy,vz
    real :: kn,Lx,kJ,Ly,Lz,Ln
    
    call read_problem_par

!   Secondary parameters


    do j = 1,ny
      yj = y(j)
      do i = 1,nx
        xi = x(i)
        do k = 1,nz
          zk = z(k)
          vx = 0.0
#ifdef SHEAR_MY
          vy = 0.0
#else
          vy = -qshear*omega*xi
#endif /* ~SHEAR_MY */
          vz = 0.0
          
          if(abs(yj) <= r0 ) then
            u(idna,i,j,k) = d0
          else
            u(idna,i,j,k) = 0.5*d0
          endif
                          
          u(imxa,i,j,k) = vx*u(idna,i,j,k)
          u(imya,i,j,k) = vy*u(idna,i,j,k)
          u(imza,i,j,k) = vz*u(idna,i,j,k)
#ifndef ISO	  
          u(iena,i,j,k) = 1.0/(gamma-1.0)!*u(idna,i,j,k)
          u(iena,i,j,k) = max(u(iena,i,j,k), smallei)
          u(iena,i,j,k) = u(iena,i,j,k) +0.5*(vx**2+vy**2+vz**2)*u(idna,i,j,k)
#endif /* ISO */
          b(1,i,j,k)   =  bx0
          b(2,i,j,k)   =  by0
          b(3,i,j,k)   =  bz0

#ifndef ISO	  
          u(iena,i,j,k)   = u(iena,i,j,k) +0.5*sum(b(:,i,j,k)**2,1)
#endif /* ISO */
        enddo
      enddo
    enddo
    return
  end subroutine init_prob  
  

end module init_problem

