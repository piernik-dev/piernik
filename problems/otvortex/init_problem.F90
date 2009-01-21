#include "piernik.def"


module init_problem

  use arrays, only : u,b
  use grid,   only : x,y,z,nx,ny,nz, xl, yl
  use initionized, only : idni,imxi,imyi,imzi
#ifndef ISO
  use initionized, only : ieni, gamma_ion
#endif /* ISO */
  use start,  only : proc, smallei, smalld, &
      rbuff, cbuff, ibuff
  use mpi_setup
  use grid, only : dx,dy
  use constants, only : pi,dpi,fpi

  real :: d0,r0,bx0,by0,bz0
  character ::  problem_name*32,run_id*3,dir*1

  namelist /PROBLEM_CONTROL/  problem_name, run_id, d0, r0,bx0,by0,bz0


contains

!-----------------------------------------------------------------------------

  subroutine read_problem_par

    implicit none
  

    problem_name = 'shock'
    run_id  = 'tst'
    d0      = 1.0
    r0      = 0.25
    
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
      rbuff(2) = r0
    
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
      r0           = rbuff(2)
    
    endif

  end subroutine read_problem_par

!-----------------------------------------------------------------------------

  subroutine init_prob

    implicit none
 
    integer i,j,k
    real :: xi,yj,zk
    real :: vx,vy,vz,rho,pre,bx,by,bz,b0
    real, dimension(:,:,:),allocatable :: A
    
    call read_problem_par

!   Secondary parameters

    if (.not.allocated(A)) allocate(A(nx,ny,1))

    rho = 25.0/(36.0*pi)
    pre =  5.0/(12.0*pi)
    b0  = 1./sqrt(fpi)
    vz = 0.0
    bz0 = 0.0

    do j=1,ny
       do i = 1,nx
         A(i,j,1) = b0*(dcos(fpi*xl(i))/fpi + dcos(dpi*yl(j))/dpi)
       enddo
    enddo

!    b(1,1:nx,1:ny-1,1)   = (A(1:nx,2:ny,1) - A(1:nx,1:ny-1,1))/dy
!    b(2,1:nx-1,1:ny,1)   = (A(2:nx,1:ny,1) - A(1:nx-1,1:ny,1))/dx
!    b(3,:,:,:)           =  0.0

    do j = 1,ny
      yj = y(j)
      do i = 1,nx
        xi = x(i)
        do k = 1,nz
          zk = z(k)

          vx  = -dsin(dpi*yj)
          vy  = dsin(dpi*xi)
          bx  = b0*vx
          by  = b0*dsin(fpi*xi)
          bz  = 0.0

          
          u(idni,i,j,k) = rho
                          
          u(imxi,i,j,k) = vx*u(idni,i,j,k)
          u(imyi,i,j,k) = vy*u(idni,i,j,k)
          u(imzi,i,j,k) = vz*u(idni,i,j,k)
#ifndef ISO	  
          u(ieni,i,j,k) = pre/(gamma_ion-1.0)
          u(ieni,i,j,k) = max(u(ieni,i,j,k), smallei)
          u(ieni,i,j,k) = u(ieni,i,j,k) +0.5*(vx**2+vy**2+vz**2)*u(idni,i,j,k)
#endif /* ISO */
          b(1,i,j,k)  = bx
          b(2,i,j,k)  = by
          b(3,i,j,k)  = bz

#ifndef ISO	  
          u(ieni,i,j,k)   = u(ieni,i,j,k) +0.5*sum(b(:,i,j,k)**2,1)
#endif /* ISO */
        enddo
      enddo
    enddo
    if (allocated(A)) deallocate(A)
    return
  end subroutine init_prob  
  

end module init_problem

