#include "piernik.def"

module init_problem
  
! Initial condition for Keplerian disk
! Written by: M. Hanasz, March 2006

  use constants
  use arrays
  use start
  use grid
  use hydrostatic

  real d0,v_zab
  character problem_name*32,run_id*3, mag_field_orient*32

  namelist /PROBLEM_CONTROL/  problem_name, run_id, &
                              d0,v_zab,mag_field_orient    


contains

!-----------------------------------------------------------------------------

  subroutine read_problem_par

    implicit none
  

    problem_name = 'aaa'
    run_id  = 'aa'
    d0      = 1.0
    v_zab   = 0.1
    mag_field_orient = 'toroidal'     
    
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
      cbuff(3) =  mag_field_orient

      rbuff(1) = d0
      rbuff(2) = v_zab
    
      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

    else
    
      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)
      
      problem_name = cbuff(1)   
      run_id       = cbuff(2)   
      mag_field_orient = cbuff(3)

      d0           = rbuff(1)  
      v_zab        = rbuff(2)
    
    endif

  end subroutine read_problem_par

!-----------------------------------------------------------------------------

  subroutine init_prob

    implicit none
 
    integer i,j,k,kmid
    real xi,yj,zk, rc, rs, vx, vy, vz, h2, dgdz, csim2, b0, sqr_gm, v_phi
    real vzab
    real, allocatable ::dprof(:), gpot(:)
    
    call read_problem_par

!   Secondary parameters

    sqr_gm = dsqrt(newtong*ptmass)

    b0 = dsqrt(2.*alpha*c_si**2)  ! multiplying by d0 is made below <it's in u(1,i,j,k)> )
    
    do j = 1,ny
      yj = y(j)
      do i = 1,nx
        xi = x(i)
        rc = sqrt(xi**2+yj**2)
        
        do k = 1,nz
          zk = z(k)

          if (rc .ge. r_smooth .and. rc .le. r_grav) then
            vzab = v_zab*c_si*dsin(pi*(rc-r_smooth)/(r_grav-r_smooth))*dsin(2.*pi*zk/Lz)
          else
            vzab = 0.0
          endif

          vx = sqr_gm * (-yj)/(rc**2+r_smooth**2)**0.75
          vy = sqr_gm * ( xi)/(rc**2+r_smooth**2)**0.75
          vz = 0.0
               
          u(1,i,j,k) = d0/cosh((rc/r_grav)**n_gravr)
          u(1,i,j,k) = max(u(1,i,j,k), smalld)
                          
          u(2,i,j,k) = vx*u(1,i,j,k)
          u(3,i,j,k) = vy*u(1,i,j,k)
          u(4,i,j,k) = vz*u(1,i,j,k)      
          u(5,i,j,k) = c_si**2/(gamma-1.0)*u(1,i,j,k)
          u(5,i,j,k) = max(u(5,i,j,k), smallei)
          
          u(5,i,j,k) = u(5,i,j,k) +0.5*(vx**2+vy**2+vz**2)*u(1,i,j,k)

          if(trim(mag_field_orient) .eq. 'toroidal') then
            b(1,i,j,k)   = -b0*dsqrt(u(1,i,j,k)/d0)*yj/rc
            b(2,i,j,k)   =  b0*dsqrt(u(1,i,j,k)/d0)*xi/rc
            b(3,i,j,k)   =  0.0
          else if(trim(mag_field_orient) .eq. 'vertical') then
            b(1,i,j,k)   =  0.0
            b(2,i,j,k)   =  0.0
            b(3,i,j,k)   =  b0*dsqrt(u(1,i,j,k))
          endif

          u(5,i,j,k)   = u(5,i,j,k) +0.5*sum(b(:,i,j,k)**2,1)
        enddo
      enddo
    enddo
    
    return
  end subroutine init_prob  
  

end module init_problem

