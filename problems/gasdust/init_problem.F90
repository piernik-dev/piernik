#include "piernik.def"

module init_problem
  
! Initial condition for dust fronts
! Written by: D. Wolt, June 2007

  use constants
  use arrays
  use start
  use grid

  real rhoa, rhol, rhor, distl, distr, vxl, vxr, collf
  integer mtr
  character problem_name*32,run_id*3

  namelist /PROBLEM_CONTROL/  problem_name, run_id, &
                              rhoa, rhol, rhor, distl, distr, vxl, vxr, collf, mtr


contains

!-----------------------------------------------------------------------------

  subroutine read_problem_par

    implicit none
  

    problem_name = 'aaa'
    run_id  = 'aa'
    rhoa      = 1.0
    rhol      = 1.0
    rhor      = 1.0
    distl     = 1.0
    distr     = 1.0
    vxl       = 1.0
    vxr       = 1.0
    collf     = 1.0
    mtr	      = 25
    
    if(proc .eq. 0) then
      open(1,file='problem.par')
        read(unit=1,nml=PROBLEM_CONTROL)
      close(1)
      open(99,file='simulation.out', position='append')
        write(99,nml=PROBLEM_CONTROL)
      close(99)
      open(3, file='tmp.log', position='append')
        write(3,nml=PROBLEM_CONTROL)
        write(3,*)
      close(3)
    endif

    if(proc .eq. 0) then


      cbuff(1) =  problem_name
      cbuff(2) =  run_id

      rbuff(1) = rhoa
      rbuff(2) = rhol
      rbuff(3) = rhor
      rbuff(4) = distl
      rbuff(5) = distr
      rbuff(6) = vxl
      rbuff(7) = vxr
      rbuff(8) = collf
      
      ibuff(1) = mtr
    
      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

    else
    
      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)
      
      problem_name = cbuff(1)   
      run_id       = cbuff(2)   

      rhoa         = rbuff(1)  
      rhol	   = rbuff(2)
      rhor	   = rbuff(3)
      distl	   = rbuff(4)
      distr	   = rbuff(5)
      vxl	   = rbuff(6)
      vxr	   = rbuff(7)
      collf	   = rbuff(8)
      
      mtr	   = ibuff(1)
    
    endif

  end subroutine read_problem_par

!-----------------------------------------------------------------------------

  subroutine init_prob

    use initdust,    only : idnd,imxd,imyd,imzd

#ifdef NEUTRAL
    use initneutral, only : gamma_neu
    use initneutral, only : idnn,imxn,imyn,imzn
#ifndef ISO    
    use initneutral, only : ienn   
#endif /* ISO */    
#endif /* NEUTRAL */ 

    implicit none
 
    integer i,j,k
    real xi,yj,zk, signmx, signpx
    
    call read_problem_par
!    collfaq = collf

    
    do i = 1,nx
      xi = x(i)
      signpx=0.0
      if(xi .ne. 0.0) signpx=abs(xi)/xi
      signmx=0.0
      if(xi .ne. 0.0) signmx=-abs(xi)/xi
      do j = 1,ny
        yj = y(j)
        do k = 1,nz

	  u(idnn,i,j,k)  = rhol
	  u(idnd,i,j,k)  = rhor

          u(imxn,i,j,k) =  vxl
	  u(imxd,i,j,k)  = -vxr
          u(imyn:imzn,i,j,k) = 0.0
          u(imyd:imzd,i,j,k) = 0.0
#ifndef ISO
	  u(ienn,i,j,k) = c_si**2*u(idna(1),i,j,k)/(gamma_neu-1.0)
	  u(ienn,i,j,k) = u(ienn,i,j,k) &
	                  +0.5*(u(imxn,i,j,k)**2+u(imyn,i,j,k)**2+u(imzn,i,j,k)**2)/u(idnn,i,j,k)
#endif /* ISO */


        enddo
      enddo
    enddo
    
    return
  end subroutine init_prob  
  

end module init_problem

