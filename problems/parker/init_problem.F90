#include "mhd.def"

module init_problem
  
! Initial condition for Parker instability in realistic galactic gravity 
! Written by: M. Hanasz, February 2006

  use arrays
  use start
  use grid
  use fluid_boundaries
  use hydrostatic

  real coldens, d0, nbx0,nby0,nbz0, a_vp, n_x
  character problem_name*32,run_id*3

  namelist /PROBLEM_CONTROL/  problem_name, run_id, &
                              d0, &
                              nbx0,nby0,nbz0, &
			      a_vp, n_x

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
      rbuff(7) = n_x

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
      n_x          = rbuff(7)

    endif

  end subroutine read_problem_par

!-----------------------------------------------------------------------------

  subroutine init_prob

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
          u(idna,i,j,k)   = dprof(k) 
          u(imxa:imza,i,j,k) = 0.0

          if ( abs(z(k)) .le. h_grav) then
            vz = a_vp * cos(2.*pi*n_x*x(i)/Lx)*cos(2.*pi*y(j)/Ly)*(cos(pi*z(k)/h_grav)+1.)/2.
          else
            vz = 0.0
          endif
!          vz = a_vp * cos(2.*pi*3*x(i)/Lx)*cos(2.*pi*y(j)/Ly)*cos(pi*z(k)/Lz)

          u(imza,i,j,k)   = u(imza,i,j,k) + u(idna,i,j,k)*vz
#ifndef ISO
          u(iena,i,j,k)   = c_si**2/(gamma-1.0) * u(idna,i,j,k) &
	                         +0.5*sum(u(imxa:imza,i,j,k)**2,1)
#endif
#ifdef COSM_RAYS
          u(iecr,i,j,k)   =  beta_cr*c_si**2 * u(idna,i,j,k)/(gamma_cr-1.0)
#endif COSM_RAYS
        enddo
      enddo
    enddo

  
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          b(ibx,i,j,k)   = b0*sqrt(u(idna,i,j,k)/d0)* nbx0/sqrt(nbx0**2+nby0**2+nbz0**2)
          b(iby,i,j,k)   = b0*sqrt(u(idna,i,j,k)/d0)* nby0/sqrt(nbx0**2+nby0**2+nbz0**2)
          b(ibz,i,j,k)   = b0*sqrt(u(idna,i,j,k)/d0)* nbz0/sqrt(nbx0**2+nby0**2+nbz0**2)
#ifndef ISO
          u(iena,i,j,k)   = u(iena,i,j,k) +0.5*sum(b(:,i,j,k)**2,1)
#endif
        enddo
      enddo
    enddo

    deallocate(dprof)
    
    return
  end subroutine init_prob
  

end module init_problem

