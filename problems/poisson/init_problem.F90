#include "mhd.def"

module init_problem
  
! Initial condition for Parker-Jeans instability (under construction)
!   in realistic galactic gravity and selfgravity 
! Written by: M. Hanasz, March 2006

  use arrays
  use start
  use grid
  use fluid_boundaries
  use hydrostatic

  character problem_name*32,run_id*3
  real coldens, d0, nbx0,nby0,nbz0, a_vp

  namelist /PROBLEM_CONTROL/  problem_name, run_id, &
                              d0, &
			      nbx0,nby0,nbz0, &
			      a_vp
    
contains


!-----------------------------------------------------------------------------

  subroutine read_problem_par

    implicit none
   
      problem_name = 'xxx'
      run_id  = 'aaa'
      d0      = 1.0
      nbx0    = 0.0
      nby0    = 1.0
      nbz0    = 0.0
      a_vp    = 0.0
    
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
      rbuff(2) = nbx0
      rbuff(3) = nby0
      rbuff(4) = nbz0
      rbuff(5) = a_vp

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
      nbx0         = rbuff(2)  
      nby0         = rbuff(3)  
      nbz0         = rbuff(4)  
      a_vp         = rbuff(5)  

    endif

  end subroutine read_problem_par

!-----------------------------------------------------------------------------

  subroutine init_prob

    use constants, only : pi
    implicit none

    integer i,j,k
    real b0, vz
    real, allocatable :: dprof(:)

    call read_problem_par

!   Secondary parameters

    b0 = sqrt(2.*alpha*d0*c_si**2) 
    
    allocate(dprof(nz))

    call hydrostatic_zeq(1, 1, d0, dprof)    

    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          u(1,i,j,k)   = dprof(k) 
          u(2:4,i,j,k) = 0.0

          vz = a_vp * cos(2.*pi*3*x(i)/Lx)*cos(2.*pi*y(j)/Ly)*cos(pi*z(k)/Lz)
          u(4,i,j,k)   = u(4,i,j,k) + u(1,i,j,k)*vz

          u(5,i,j,k)   = c_si**2/(gamma-1.0) * u(1,i,j,k) &
	                         +0.5*sum(u(2:4,i,j,k)**2,1)/u(1,i,j,k)
        enddo
      enddo
    enddo

  
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          b(1,i,j,k)   = b0*sqrt(u(1,i,j,k)/d0)* nbx0/sqrt(nbx0**2+nby0**2+nbz0**2)
          b(2,i,j,k)   = b0*sqrt(u(1,i,j,k)/d0)* nby0/sqrt(nbx0**2+nby0**2+nbz0**2)
          b(3,i,j,k)   = b0*sqrt(u(1,i,j,k)/d0)* nbz0/sqrt(nbx0**2+nby0**2+nbz0**2)
          u(5,i,j,k)   = u(5,i,j,k) +0.5*sum(b(:,i,j,k)**2,1)
        enddo
      enddo
    enddo

    deallocate(dprof)
    
    return
  end subroutine init_prob
  

end module init_problem

