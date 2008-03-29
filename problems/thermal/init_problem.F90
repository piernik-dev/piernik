#include "mhd.def"

module init_problem
  
! Initial condition for Parker instability in realistic galactic gravity 
! Written by: M. Hanasz, February 2006

  use arrays
  use constants
  use start
  use grid
  use fluid_boundaries
  use hydrostatic

  real coldens, d0, T0, nbx0,nby0,nbz0, a_pert, h_pert
  character problem_name*32,run_id*3


  namelist /PROBLEM_CONTROL/  problem_name, run_id, &
                              d0, T0, &
                              nbx0, nby0, nbz0, &
			      h_pert, a_pert

contains

!-----------------------------------------------------------------------------

  subroutine read_problem_par

    implicit none
    
      problem_name = 'xxx'
      run_id  = 'aaa'
      d0      = 1.0
      T0      = 7000.0
      nbx0    = 0.0
      nby0    = 1.0
      nbz0    = 0.0
      h_pert  = Lz
      a_pert  = 0.0
    
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
      rbuff(2) = T0
      rbuff(3) = nbx0
      rbuff(4) = nby0
      rbuff(5) = nbz0
      rbuff(6) = h_pert
      rbuff(7) = a_pert

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
      T0           = rbuff(2)  
      nbx0         = rbuff(3)  
      nby0         = rbuff(4)  
      nbz0         = rbuff(5)  
      h_pert       = rbuff(6)  
      a_pert       = rbuff(7)  

    endif

  end subroutine read_problem_par

!-----------------------------------------------------------------------------

  subroutine init_prob

    implicit none

    integer i,j,k
    real b0, vx,vy,vz,cs,va
    real, allocatable :: dprof(:), eprof(:), tprof(:), bprof(:), cfunc(:)



    call read_problem_par


    allocate(dprof(nz),eprof(nz),tprof(nz),bprof(nz),cfunc(nz))


#ifndef GRAV_NULL

      if (coolheat) then

        call hydro_thermal_zeq(0, 0, T0,dprof,eprof, tprof,bprof)    
	
      else

        call hydrostatic_zeq(0, 0, d0, dprof)    

      endif
      
#else /* GRAV_NULL */

      if (coolheat) then

        call hydro_thermal_zeq(0, 0, T0,dprof,eprof, tprof,bprof)    
        cs = sqrt(gamma*T0*k_B / hydro_mass)
	
	cs = sqrt(gamma*(gamma-1.)*eprof(1)/dprof(1))
        va = sqrt(bprof(1)**2/dprof(1))
	write(*,*) 'INIT:   cs=',cs, '    va=', va
	
      else
      
	dprof = d0
     
      endif

#endif /* GRAV_NULL */


    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          u(1,i,j,k)   = dprof(k) 
          u(2:4,i,j,k) = 0.0

          vx=0.0
          vy=0.0
          vz=0.0
#ifndef NO_GRAV
            if(abs(z(k)) .lt. h_pert) then
              vz = a_pert &
	           * cos(2.*pi*3*x(i)/Lx)*cos(2.*pi*y(j)/Ly)*cos(0.5*pi*z(k)/h_pert)
            endif
#else /* NO_GRAV */
!	    u(1,i,j,k) = u(1,i,j,k) + a_pert &
!	           * cos(2.*pi*x(i)/Lx)*cos(2.*pi*y(j)/Ly)*cos(2.*pi*z(k)/Lz)	  

            vx = -a_pert*cs* sin(2.*pi*x(i)/Lx)*cos(2.*pi*y(j)/Ly)*cos(2.*pi*z(k)/Lz)
            vy = -a_pert*cs* cos(2.*pi*x(i)/Lx)*sin(2.*pi*y(j)/Ly)*cos(2.*pi*z(k)/Lz)
            vz = -a_pert*cs* cos(2.*pi*x(i)/Lx)*cos(2.*pi*y(j)/Ly)*sin(2.*pi*z(k)/Lz)
#endif /* NO_GRAV */
	  
	  
          u(2,i,j,k)   = u(2,i,j,k) + u(1,i,j,k)*vx
          u(3,i,j,k)   = u(3,i,j,k) + u(1,i,j,k)*vy
          u(4,i,j,k)   = u(4,i,j,k) + u(1,i,j,k)*vz

          u(5,i,j,k)   = eprof(k) &
	                 +0.5*sum(u(2:4,i,j,k)**2,1)/u(1,i,j,k)
        enddo
      enddo
    enddo

  
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          b(1,i,j,k)   = bprof(k) * nbx0/sqrt(nbx0**2+nby0**2+nbz0**2)
          b(2,i,j,k)   = bprof(k) * nby0/sqrt(nbx0**2+nby0**2+nbz0**2)
          b(3,i,j,k)   = bprof(k) * nbz0/sqrt(nbx0**2+nby0**2+nbz0**2)
          u(5,i,j,k)   = u(5,i,j,k) +0.5*sum(b(:,i,j,k)**2,1)
        enddo
      enddo
    enddo

    deallocate(dprof)
    
    return
  end subroutine init_prob
  

end module init_problem

