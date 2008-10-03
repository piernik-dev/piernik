#include "piernik.def"


module init_problem
  
! Initial condition for Keplerian disk
! Written by: M. Hanasz, March 2006

  use constants
  use arrays
  use start
  use grid
  use hydrostatic

  real d0
  character problem_name*32,run_id*3, mag_field_orient*32

  namelist /PROBLEM_CONTROL/  problem_name, run_id, &
                              d0,mag_field_orient    


contains

!-----------------------------------------------------------------------------

  subroutine read_problem_par

    implicit none
  

    problem_name = 'aaa'
    run_id  = 'aa'
    d0      = 1.0
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
    
    endif

  end subroutine read_problem_par

!-----------------------------------------------------------------------------

  subroutine init_prob

    implicit none
 
    integer i,j,k,kmid
    real xi,yj,zk, rc, rs, vx, vy, vz, h2, dgdz, csim2, b0, sqr_gm, v_phi
    real, allocatable ::dprof(:), gpot(:)
    
    call read_problem_par

!   Secondary parameters

    sqr_gm = sqrt(newtong*ptmass)

    b0 = sqrt(2.*alpha*d0*c_si**2) 
    write(*,*) b0 
    allocate(dprof(nz),gpot(nz))
    
    do k=1, nz
      if(z(k) .lt. 0.0) kmid = k       ! the midplane is in between 
    enddo                                  ! ksmid and ksmid+1
      
    
    do j = 1,ny
      yj = y(j)
      do i = 1,nx
        xi = x(i)
        rc = sqrt(xi**2+yj**2)
        
        if(dimensions .eq.'3d') then
          call hydrostatic_zeq(i, j, d0, dprof)    
        endif
                          
        do k = 1,nz
                  
          vx = sqr_gm * (-yj)/(rc**2+r_smooth**2)**0.75
          vy = sqr_gm * ( xi)/(rc**2+r_smooth**2)**0.75
          vz = 0.0
                  
          if(dimensions .eq.'3d') then
            u(idna,i,j,k) = min((rc/r_grav)**n_gravr,100.0)
            u(idna,i,j,k) = dprof(k)/cosh(u(idna,i,j,k))
          else
            u(idna,i,j,k) = min((rc/r_grav)**n_gravr,100.0)
            u(idna,i,j,k) = d0/cosh(u(idna,i,j,k))
          endif
          u(idna,i,j,k) = max(u(idna,i,j,k), smalld)
                          
          u(imxa,i,j,k) = vx*u(idna,i,j,k)
          u(imya,i,j,k) = vy*u(idna,i,j,k)
          u(imza,i,j,k) = vz*u(idna,i,j,k)      
#ifndef ISO	  
          u(iena,i,j,k) = c_si**2/(gamma-1.0)*u(idna,i,j,k)
          u(iena,i,j,k) = max(u(iena,i,j,k), smallei)
          u(iena,i,j,k) = u(iena,i,j,k) +0.5*(vx**2+vy**2+vz**2)*u(idna,i,j,k)
#endif /* ISO */
!          if(trim(mag_field_orient) .eq. 'toroidal') then
!            b(1,i,j,k)   = -b0*sqrt(u(idna,i,j,k)/d0)*yj/rc
!            b(2,i,j,k)   =  b0*sqrt(u(idna,i,j,k)/d0)*xi/rc
!            b(3,i,j,k)   =  0.0
!          else if(trim(mag_field_orient) .eq. 'vertical') then
            b(1,i,j,k)   =  1.0
            b(2,i,j,k)   =  1.0
            b(3,i,j,k)   =  0.0
!          endif

#ifndef ISO	  
          u(iena,i,j,k)   = u(iena,i,j,k) +0.5*sum(b(:,i,j,k)**2,1)
#endif /* ISO */
        enddo
      enddo
    enddo
    
    deallocate(dprof,gpot)
    
    return
  end subroutine init_prob  
  

end module init_problem

