#include "mhd.def"

module init_problem
  
! Initial condition for Parker instability in realistic galactic gravity 
! Written by: M. Hanasz, February 2006

  use arrays
  use start
  use grid
  use fluid_boundaries
  use hydrostatic
  
  real d0, nbx0,nby0,nbz0, a_vp, n_x
  real x0,y0,z0,r0
  character problem_name*32,run_id*3

  namelist /PROBLEM_CONTROL/  problem_name, run_id, &
                              d0, &
                              nbx0,nby0,nbz0, &
			      a_vp, n_x, &
                              x0,y0,z0, r0 
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
      rbuff(2) = nbx0
      rbuff(3) = nby0
      rbuff(4) = nbz0
      rbuff(5) = a_vp
      rbuff(7) = n_x
      rbuff(8) = x0
      rbuff(9) = y0
      rbuff(10)= z0

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
      x0           = rbuff(8) 
      y0           = rbuff(9)
      z0           = rbuff(10)

    endif

  end subroutine read_problem_par

!-----------------------------------------------------------------------------

  subroutine init_prob
  
    use constants

    implicit none

    integer i,j,k
    real b0, vz
    real mass


    call read_problem_par

!   Secondary parameters

    b0 = sqrt(2.*alpha*d0*c_si**2) 
    
    call hydrostatic_zeq(1, 1, d0, dprof)    

    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          u(idna,i,j,k)   = max(smalld,dprof(k)) 
	  
          u(imxa:imza,i,j,k) = 0.0
#ifdef SHEAR
          u(imya,i,j,k) = -qshear*omega*x(i)*u(idna,i,j,k)
#endif SHEAR

          if ( abs(z(k)) .le. h_grav) then
            vz = a_vp * cos(2.*pi*n_x*x(i)/Lx)*cos(2.*pi*y(j)/Ly)*(cos(pi*z(k)/h_grav)+1.)/2.
          else
            vz = 0.0
          endif

          vz = a_vp * cos(2.*pi*3*x(i)/Lx)*cos(2.*pi*y(j)/Ly)*cos(pi*z(k)/Lz)
          u(imza,i,j,k)   = u(imza,i,j,k) + u(idna,i,j,k)*vz

#ifndef ISO
          u(iena,i,j,k)   = c_si**2/(gamma-1.0) * u(idna,i,j,k) &
	                         +0.5*sum(u(imxa:imza,i,j,k)**2,1)
#endif
#ifdef COSM_RAYS
          u(iecr,i,j,k)   =  beta_cr*c_si**2 * u(idna,i,j,k)/(gamma_cr-1.0)
#ifdef GALAXY
          u(iecr,i,j,k)= u(iecr,i,j,k) &
	     + amp_cr*ethu*exp(-((x(i)-x0)**2+(y(j)-y0)**2+(z(k)-z0)**2)/r_sn**2)  &
             + amp_cr*ethu*exp(-((x(i)-(x0+Lx))**2+(y(j)-y0)**2+(z(k)-z0)**2)/r_sn**2) &
             + amp_cr*ethu*exp(-((x(i)-x0)**2+(y(j)-(y0+Ly))**2+(z(k)-z0)**2)/r_sn**2) &
             + amp_cr*ethu*exp(-((x(i)-(x0+Lx))**2+(y(j)-(y0+Ly))**2+(z(k)-z0)**2)/r_sn**2)
#endif GALAXY
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
    
    mass = sum(u(idna,is:ie,js:je,ks:ke)) * dvol
    call mpi_allreduce(mass, init_mass, 1, mpi_real8, mpi_sum, comm3d, ierr)
      
    
    return
  end subroutine init_prob
  
!-----------------------------------------------------------------------------
  
  subroutine mass_loss_compensate
  
    implicit none
    
    real mass,tot_mass
    integer k
    
    mass = sum(u(idna,is:ie,js:je,ks:ke)) * dvol
    call mpi_allreduce(mass, tot_mass, 1, mpi_real8, mpi_sum, comm3d, ierr)

    mass_loss = max(0.0, init_mass - tot_mass )
    mass_loss_tot = mass_loss_tot + mass_loss

!     Correction for mass loss through open z-boundaries

    do k=1,nz
      u(idna,:,:,k) = u(idna,:,:,k) + mass_loss/init_mass * dprof(k)
#ifndef ISO
      u(iena,:,:,k) = u(iena,:,:,k) + mass_loss/init_mass * eprof(k)
#endif
    enddo
  
  end subroutine mass_loss_compensate

end module init_problem

