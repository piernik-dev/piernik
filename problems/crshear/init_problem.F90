#include "piernik.def"

module init_problem
  
! Initial condition for Sedov-Taylor explosion
! Written by: M. Hanasz, March 2006

  use start
  use arrays
  use grid
  use mpi_setup

  real t_sn
  integer n_sn
  real d0,p0,bx0,by0,bz0, x0,y0,z0
  character problem_name*32,run_id*3

  namelist /PROBLEM_CONTROL/  problem_name, run_id, &
                              d0,p0, bx0,by0,bz0, x0,y0,z0

contains

!-----------------------------------------------------------------------------

  subroutine read_problem_par

    implicit none
    
      t_sn = 0.0
  
      problem_name = 'aaa'
      run_id  = 'aaa'
      d0      = 1.0
      p0      = 1.e-3
      bx0     =   0.
      by0     =   0.
      bz0     =   0.
      x0      = 0.0 
      y0      = 0.0 
      z0      = 0.0 
         
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
      rbuff(2) = p0
      rbuff(3) = bx0
      rbuff(4) = by0
      rbuff(5) = bz0
      rbuff(6) = x0
      rbuff(7) = y0
      rbuff(8) = z0
      
      ibuff(1) = n_sn
    
      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

    else
    
      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)
      
!  namelist /PROBLEM_CONTROL/  problem_name, run_id, &
!                              d0,p0, bx0,by0,bz0,  x0,y0,z0

      problem_name = cbuff(1)   
      run_id       = cbuff(2)   

      d0           = rbuff(1)  
      p0           = rbuff(2)  
      bx0          = rbuff(3)  
      by0          = rbuff(4)  
      bz0          = rbuff(5)  
      x0           = rbuff(6) 
      y0           = rbuff(7)
      z0           = rbuff(8)

    endif

  end subroutine read_problem_par

!-----------------------------------------------------------------------------

  subroutine init_prob

    implicit none

    integer i,j,k, n
    
    
    call read_problem_par
 
! Uniform equilibrium state

    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          u(idna,i,j,k)   = d0 
          u(imxa:imza,i,j,k) = 0.0
#ifdef SHEAR
          u(imya,i,j,k) = -qshear*omega*x(i)*u(idna,i,j,k)
#endif /* SHEAR */
#ifndef ISO	  	  	  
          u(iena,i,j,k)   = p0/(gamma-1.0)
	  u(iena,i,j,k)   = u(iena,i,j,k) &
	                  + 0.5*sum(u(imxa:imza,i,j,k)**2,1)/u(idna,i,j,k)
#endif /* ISO */

#ifdef COSM_RAYS
          u(iecr,i,j,k) =  beta_cr*c_si**2 * u(idna,i,j,k)/(gamma_cr-1.0)
#endif /* COSM_RAYS */

          b(ibx,i,j,k)   = bx0
          b(iby,i,j,k)   = by0
          b(ibz,i,j,k)   = bz0
#ifndef ISO	  	  
          u(iena,i,j,k)   = u(iena,i,j,k) + 0.5*sum(b(:,i,j,k)**2,1)
#endif /* ISO */
        enddo
      enddo
    enddo
    
! CR Explosions
      do k = ks,ke
        do j = nb+1,ny-nb
          do i = nb+1,nx-nb
#ifdef COSM_RAYS
            u(iecr,i,j,k)= u(iecr,i,j,k) &
	     + amp_cr*ethu*exp(-((x(i)-x0)**2+(y(j)-y0)**2+(z(k)-z0)**2)/r_sn**2) &
             + amp_cr*ethu*exp(-((x(i)-(x0+Lx))**2+(y(j)-y0)**2+(z(k)-z0)**2)/r_sn**2) &
             + amp_cr*ethu*exp(-((x(i)-x0)**2+(y(j)-(y0+Ly))**2+(z(k)-z0)**2)/r_sn**2) &
             + amp_cr*ethu*exp(-((x(i)-(x0+Lx))**2+(y(j)-(y0+Ly))**2+(z(k)-z0)**2)/r_sn**2)
#endif /* COSM_RAYS */
          enddo
        enddo
      enddo
    
    return
  end subroutine init_prob  

end module init_problem

