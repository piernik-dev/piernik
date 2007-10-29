#include "mhd.def"

module init_problem
  
! Initial condition for Parker instability in realistic galactic gravity 
! Written by: M. Hanasz, February 2006

  use arrays
  use start
  use grid
  use fluid_boundaries
  use hydrostatic

  real ::  coldens, d0, nbx0,nby0,nbz0, a_vp, n_x, r0,x0,y0,z0
  character problem_name*32,run_id*3

  namelist /PROBLEM_CONTROL/  problem_name, run_id, &
                              d0, &
                              nbx0,nby0,nbz0, &
			      a_vp, n_x, r0,x0,y0,z0

contains

!-----------------------------------------------------------------------------

  subroutine read_problem_par

    implicit none
   
      problem_name = 'xxx'
      run_id  = 'aaa'
      d0      = 1.0
      x0      = 0.0
      y0      = 0.0
      z0      = 0.0
      r0      = 0.5
      nbx0    = 0.0
      nby0    = 1.0
      nbz0    = 0.0
      a_vp    = 0.0
      n_x     = 3.0
    
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
      rbuff(6) = omega
      rbuff(7) = n_x
      rbuff(8) = x0
      rbuff(9) = y0
      rbuff(10) = z0
      rbuff(11) = r0

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
      omega        = rbuff(6)
      n_x          = rbuff(7)
      x0           = rbuff(8)
      y0           = rbuff(9)
      z0           = rbuff(10)
      r0           = rbuff(11)


    endif

  end subroutine read_problem_par

!-----------------------------------------------------------------------------

  subroutine init_prob

    implicit none

    integer i,j,k
    real b0, vz
    real, allocatable :: dprof(:)
    real, dimension(3) :: rand


    call read_problem_par

!   Secondary parameters

    b0 = sqrt(2.*alpha*d0*c_si**2) 
    call random_seed()
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          call random_number(rand)
          u(imxa:imza,i,j,k) = 0.0
          u(imya,i,j,k) = -qshear*omega*x(i)*u(idna,i,j,k)
          u(imxa:imza,i,j,k) = u(imxa:imza,i,j,k) + 0.1*(rand(1)-0.5)

#ifndef ISO
          u(iena,i,j,k)   = 1.0/(gamma-1.0) +0.5*sum(u(imxa:imza,i,j,k)**2,1)
#endif
        enddo
      enddo
    enddo

  
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          b(ibx,i,j,k)   = nbx0 !b0*sqrt(u(idna,i,j,k)/d0)* nbx0/sqrt(nbx0**2+nby0**2+nbz0**2+small)
          b(iby,i,j,k)   = nby0 !b0*sqrt(u(idna,i,j,k)/d0)* nby0/sqrt(nbx0**2+nby0**2+nbz0**2+small)
          b(ibz,i,j,k)   = nbz0 !b0*sqrt(u(idna,i,j,k)/d0)* nbz0/sqrt(nbx0**2+nby0**2+nbz0**2+small)
#ifndef ISO
          u(iena,i,j,k)   = u(iena,i,j,k) +0.5*sum(b(:,i,j,k)**2,1)
#endif
        enddo
      enddo
    enddo

    return
  end subroutine init_prob
  

end module init_problem

