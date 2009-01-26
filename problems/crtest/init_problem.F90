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
  real d0,p0,bx0,by0,bz0,x0,y0,z0,r0,beta_cr, amp_cr 
  character problem_name*32,run_id*3

  namelist /PROBLEM_CONTROL/  problem_name, run_id,      &
                              d0,p0, bx0,by0,bz0, &
                              x0,y0,z0, r0, &
                              beta_cr, amp_cr 

contains

!-----------------------------------------------------------------------------

  subroutine read_problem_par

    implicit none
    
      t_sn = 0.0
  
      problem_name = 'aaa'
      run_id  = 'aaa'
      d0      = 1.0
      p0      = 1.0
      bx0     =   0.
      by0     =   0.
      bz0     =   0.
      x0      = 0.0 
      y0      = 0.0 
      z0      = 0.0 
      r0      = dxmn/2.
      
      beta_cr    = 0.0
      amp_cr     = 1.0   
      
         
    if(proc .eq. 0) then
    
      open(1,file='problem.par')
        read(unit=1,nml=PROBLEM_CONTROL,end=1234,err=1234)

        goto 4567
 1234   write(*,*) 'Check namelist in "problem.par":'
        write(*,*) ''
        write(*,nml=PROBLEM_CONTROL) 
        write(*,*) ''
        stop
 4567   continue   
   
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
      rbuff(9) = r0
      
      rbuff(10)= beta_cr
      rbuff(11)= amp_cr
     
      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

    else
    
      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)
      
!  namelist /PROBLEM_CONTROL/  problem_name, run_id, &
!                              d0,p0, bx0,by0,bz0, x0,y0,z0, r0 &
!

      problem_name = cbuff(1)   
      run_id       = cbuff(2)   

      d0           = rbuff(1)  
      p0           = rbuff(2)  
      bx0          = rbuff(3)  
      by0          = rbuff(4)  
      bz0          = rbuff(5)  
      x0           = rbuff(7) 
      y0           = rbuff(8)
      z0           = rbuff(9)
      r0           = rbuff(10)
      
      beta_cr      = rbuff(10)
      amp_cr       = rbuff(11)

    endif

  end subroutine read_problem_par

!-----------------------------------------------------------------------------

  subroutine init_prob

    use fluidindex,     only : ibx,iby,ibz
    use initionized,    only : idni,imxi,imyi,imzi,ieni
    use initcosmicrays, only : gamma_cr,iecr
    use initionized,    only : gamma_ion   

    implicit none

    integer :: i, j, k, n
    real    :: cs_iso, r2
    
 
! Uniform equilibrium state

    cs_iso = sqrt(p0/d0)

    do k = 1,nz
      do j = 1,ny
        do i = 1,nx

          b(ibx,i,j,k)       = bx0
          b(iby,i,j,k)       = by0
          b(ibz,i,j,k)       = bz0

          u(idni,i,j,k)      = d0 
          u(imxi:imzi,i,j,k) = 0.0
#ifndef ISO	  	  	  
          u(ieni,i,j,k)      = p0/(gamma_ion-1.0)
          u(ieni,i,j,k)      = u(ieni,i,j,k) &
                    + 0.5*sum(u(imxi:imzi,i,j,k)**2,1)/u(idni,i,j,k)
          u(ieni,i,j,k)      = u(ieni,i,j,k) + 0.5*sum(b(:,i,j,k)**2,1)
#endif /* ISO */

#ifdef COSM_RAYS
          u(iecr,i,j,k)      =  beta_cr*cs_iso**2 * u(idni,i,j,k)/(gamma_cr-1.0)
#endif /* COSM_RAYS */

        enddo
      enddo
    enddo
    
! Explosions

#ifdef COSM_RAYS
    do k = ks,ke
      do j = nb+1,ny-nb
        do i = nb+1,nx-nb
            r2 = (x(i)-x0)**2+(y(j)-y0)**2+(z(k)-z0)**2
            u(iecr,i,j,k)= u(iecr,i,j,k) + amp_cr*exp(-r2/r0**2)
        enddo
      enddo
    enddo
#endif /* COSM_RAYS */
  
    write(*,*) 'maxecr =',maxval(u(iecr,:,:,:))
    write(*,*) amp_cr
    
    return
  end subroutine init_prob  

end module init_problem

