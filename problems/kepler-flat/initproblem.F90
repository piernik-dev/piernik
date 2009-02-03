#include "piernik.def"

module initproblem
  
! Initial condition for Keplerian disk
! Written by: M. Hanasz, March 2006

  real :: d0, v_zab, alpha
  character(len=32) :: problem_name,mag_field_orient
  character(len=3)  :: run_id

  namelist /PROBLEM_CONTROL/  problem_name, run_id, alpha,&
                              d0, v_zab, mag_field_orient    

contains

!-----------------------------------------------------------------------------

  subroutine read_problem_par
    use mpisetup
    use func, only : namelist_errh
    implicit none
    integer :: errh

    problem_name = 'kepler'
    run_id  = 'flt'
    d0      = 1.0
    v_zab   = 0.1
    mag_field_orient = 'none'
    
    if(proc == 0) then
      open(1,file='problem.par')
        read(unit=1,nml=PROBLEM_CONTROL,iostat=errh)
        call namelist_errh(errh,'PROBLEM_CONTROL')
        write(*,nml=PROBLEM_CONTROL)
      close(1)
      open(3, file='tmp.log', position='append')
        write(3,nml=PROBLEM_CONTROL)
      close(3)
    endif

    if(proc == 0) then

      cbuff(1) =  problem_name
      cbuff(2) =  run_id
      cbuff(3) =  mag_field_orient

      rbuff(1) = d0
      rbuff(2) = v_zab
      rbuff(3) = alpha
    
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
      alpha        = rbuff(3)
 
    endif

  end subroutine read_problem_par

!-----------------------------------------------------------------------------

  subroutine init_prob
    use mpisetup, only : smalld, smallei
    use fluidindex, only : ibx, iby, ibz
    use gravity, only : ptmass, r_grav, r_smooth, n_gravr
    use arrays, only : u, b
    use grid, only : x, y, z, Lx, Ly, Lz, nx, ny, nz
    use constants, only   : newtong, pi, big
    use initionized, only : idni, imxi, imyi, imzi, ieni, gamma_ion, cs_ion
    implicit none
 
    integer :: i, j, k, kmid
    real :: xi,yj,zk, rc, rs, vx, vy, vz, h2, dgdz, csim2, b0, sqr_gm, v_phi
    real :: vzab
    
    call read_problem_par

!   Secondary parameters

    sqr_gm = dsqrt(newtong*ptmass)

    b0 = dsqrt(2.*alpha*cs_ion**2)  ! multiplying by d0 is made below <it's in u(1,i,j,k)> )
    
    do j = 1,ny
      yj = y(j)
      do i = 1,nx
        xi = x(i)
        rc = sqrt(xi**2+yj**2)
        
        do k = 1,nz
          zk = z(k)

          if (rc .ge. r_smooth .and. rc .le. r_grav) then
            vzab = v_zab*cs_ion*dsin(pi*(rc-r_smooth)/(r_grav-r_smooth))*dsin(2.*pi*zk/Lz)
          else
            vzab = 0.0
          endif

          vx = sqr_gm * (-yj)/(rc**2+r_smooth**2)**0.75
          vy = sqr_gm * ( xi)/(rc**2+r_smooth**2)**0.75
          vz = 0.0
          
          if( rc/r_grav > 1.5) then
             u(idni,i,j,k) = smalld
          else
             u(idni,i,j,k) = max(d0/dcosh((rc/r_grav)**n_gravr),smalld)
          endif 
          u(imxi,i,j,k) = vx*u(idni,i,j,k)
          u(imyi,i,j,k) = vy*u(idni,i,j,k)
          u(imzi,i,j,k) = vz*u(idni,i,j,k)      
          u(ieni,i,j,k) = cs_ion**2/(gamma_ion-1.0)*u(idni,i,j,k)
          u(ieni,i,j,k) = max(u(ieni,i,j,k), smallei)
          
          u(ieni,i,j,k) = u(ieni,i,j,k) +0.5*(vx**2+vy**2+vz**2)*u(idni,i,j,k)

          if(trim(mag_field_orient) .eq. 'toroidal') then
            b(ibx,i,j,k)   = -b0*dsqrt(u(idni,i,j,k)/d0)*yj/rc
            b(iby,i,j,k)   =  b0*dsqrt(u(idni,i,j,k)/d0)*xi/rc
            b(ibz,i,j,k)   =  0.0
          else if(trim(mag_field_orient) .eq. 'vertical') then
            b(ibx,i,j,k)   =  0.0
            b(iby,i,j,k)   =  0.0
            b(ibz,i,j,k)   =  b0*dsqrt(u(idni,i,j,k))
          endif

          u(ieni,i,j,k)   = u(ieni,i,j,k) +0.5*sum(b(:,i,j,k)**2,1)
        enddo
      enddo
    enddo
    
    return
  end subroutine init_prob  
  

end module initproblem

