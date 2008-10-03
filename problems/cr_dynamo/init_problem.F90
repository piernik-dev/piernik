#include "piernik.def"

module init_problem
  
! Initial condition for the cosmic ray driven dynamo
! Based on Parker instability setup
! Written by: M. Hanasz, February 2006
! Modified by M.Hanasz for CR-driven dynamo

  use arrays
  use start
  use grid
  use fluid_boundaries
  use hydrostatic
  
  real d0, bxn,byn,bzn
  character problem_name*32,run_id*3

  namelist /PROBLEM_CONTROL/  problem_name, run_id, &
                              d0, &
                              bxn,byn,bzn, &
			      x0, y0, z0
contains

!-----------------------------------------------------------------------------

  subroutine read_problem_par

    implicit none
   
      problem_name = 'xxx'
      run_id  = 'aaa'
      d0      = 1.0
      bxn    = 0.0
      byn    = 1.0
      bzn    = 0.0
      x0     = 0.0
      y0     = 0.0
      z0     = 0.0
    
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
      rbuff(2) = bxn
      rbuff(3) = byn
      rbuff(4) = bzn
      rbuff(5) = x0
      rbuff(6) = y0
      rbuff(6) = z0

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
      bxn          = rbuff(2)  
      byn          = rbuff(3)  
      bzn          = rbuff(4)  
      x0           = rbuff(5) 
      y0           = rbuff(6) 
      z0           = rbuff(6) 

    endif

  end subroutine read_problem_par

!-----------------------------------------------------------------------------

  subroutine init_prob

    implicit none

    integer i,j,k
    real b0

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
#endif /* SHEAR */


#ifndef ISO
          u(iena,i,j,k)   = c_si**2/(gamma-1.0) * u(idna,i,j,k) &
	                         +0.5*sum(u(imxa:imza,i,j,k)**2,1)/u(idna,i,j,k)
#endif /* ISO */
#ifdef COSM_RAYS
          u(iecr,i,j,k)   =  beta_cr*c_si**2 * u(idna,i,j,k)/(gamma_cr-1.0)
#ifdef GALAXY
! Single SN explosion in x0,y0,z0 at t = 0 if amp_cr /= 0
          u(iecr,i,j,k)= u(iecr,i,j,k) &
	     + amp_cr*exp(-((x(i)-x0)**2+(y(j)-y0)**2+(z(k)-z0)**2)/r_sn**2)  &
             + amp_cr*exp(-((x(i)-(x0+Lx))**2+(y(j)-y0)**2+(z(k)-z0)**2)/r_sn**2) &
             + amp_cr*exp(-((x(i)-x0)**2+(y(j)-(y0+Ly))**2+(z(k)-z0)**2)/r_sn**2) &
             + amp_cr*exp(-((x(i)-(x0+Lx))**2+(y(j)-(y0+Ly))**2+(z(k)-z0)**2)/r_sn**2)
#endif /* GALAXY */
#endif /* COSM_RAYS */
        enddo
      enddo
    enddo
  
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          b(ibx,i,j,k)   = b0*sqrt(u(idna,i,j,k)/d0)* bxn/sqrt(bxn**2+byn**2+bzn**2)
          b(iby,i,j,k)   = b0*sqrt(u(idna,i,j,k)/d0)* byn/sqrt(bxn**2+byn**2+bzn**2)
          b(ibz,i,j,k)   = b0*sqrt(u(idna,i,j,k)/d0)* bzn/sqrt(bxn**2+byn**2+bzn**2)
#ifndef ISO
          u(iena,i,j,k)   = u(iena,i,j,k) +0.5*sum(b(:,i,j,k)**2,1)
#endif /* ISO */
        enddo
      enddo
    enddo
              
    return
  end subroutine init_prob
  
!-----------------------------------------------------------------------------

  subroutine mass_loss_compensate
  
    implicit none
    
    real tot_mass, dmass
    integer i
    
      call total_mass(tot_mass)

      mass_loss = max(0.0,init_mass - tot_mass) 

      dmass = mass_loss/init_mass

      do i=1,nx
        u(idna,:,:,:) = u(idna,:,:,:) + dmass * dinit(:,:,:)
#ifdef SHEAR
        u(imya,:,:,:) = u(imya,:,:,:) - dmass * qshear*omega*x(i) * dinit(:,:,:)
#endif /* SHEAR */
#ifndef ISO
        u(iena,:,:,:) = u(iena,:,:,:) + dmass * (c_si**2/(gamma-1.0) * dinit(:,:,:) 
#ifdef SHEAR
	u(iena,:,:,:) = u(iena,:,:,:) + dmass * 0.5 * (qshear*omega*x(i))**2) * dinit(:,:,:)
#endif /* SHEAR */
#endif /* ISO */
      enddo

  end subroutine mass_loss_compensate
  
!=============================================================================
! Te procedury powinny sie znalezc docelowo w jakims innym module. 
   
  subroutine save_init_dens
    implicit none
    real mass
    
      dinit(:,:,:) = u(idna,:,:,:)
      mass = sum(dinit(is:ie,js:je,ks:ke)) * dvol
      call MPI_ALLREDUCE(mass, init_mass, 1, mpi_real8, mpi_sum, comm3d, ierr)

  end subroutine save_init_dens
  
!-----------------------------------------------------------------------------
   
  subroutine get_init_mass
    implicit none
    real mass
    
      mass = sum(dinit(is:ie,js:je,ks:ke)) * dvol
      call MPI_ALLREDUCE(mass, init_mass, 1, mpi_real8, mpi_sum, comm3d, ierr)

  end subroutine get_init_mass
  
!-----------------------------------------------------------------------------
   
  
  subroutine total_mass(tmass)

    implicit none
    real mass, tmass
   
      mass = sum(u(idna,is:ie,js:je,ks:ke)) * dvol
      call MPI_ALLREDUCE(mass, tmass, 1, mpi_real8, mpi_sum, comm3d, ierr)

  end subroutine total_mass

end module init_problem

