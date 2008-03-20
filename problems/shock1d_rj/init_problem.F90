#include "mhd.def"


module init_problem
! ----------------------------------------- ! 
! Initial condition for a 1D MHD shock tube !
! Written by: K. Kowalik, Mar 2008          !
! See: Ryu, Jones, ApJ 442:228-258, (1995)  !
!    and reference therein                  !
! ----------------------------------------- !

  use arrays, only : x, u, b,x,y,z,nx,ny,nz, &
      idna,imxa,imya,imza
#ifndef ISO
  use arrays, only : iena
#endif /* ISO */
  use start,  only : qshear, omega, proc, smallei, smalld, gamma, &
      rbuff, cbuff, ibuff
  use mpi_setup

  real :: dl,vxl,vyl,vzl,bxl,byl,bzl,el
  real :: dr,vxr,vyr,vzr,bxr,byr,bzr,er
  character ::  problem_name*32,run_id*3,dir*1

  namelist /PROBLEM_CONTROL/  problem_name, &
     run_id,dl,vxl,vyl,vzl,bxl,byl,bzl,el,  &
     dr,vxr,vyr,vzr,bxr,byr,bzr,er


contains

!-----------------------------------------------------------------------------

  subroutine read_problem_par

    implicit none
  

    problem_name = 'shock'
    run_id  = 'tst'
    
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

      rbuff(1) = dl
      rbuff(2) = vxl
      rbuff(3) = vyl
      rbuff(4) = vzl
      rbuff(5) = bxl
      rbuff(6) = byl
      rbuff(7) = bzl
      rbuff(8) = el
      rbuff(9) = dr
      rbuff(10) = vxr
      rbuff(11) = vyr
      rbuff(12) = vzr
      rbuff(13) = bxr
      rbuff(14) = byr
      rbuff(15) = bzr
      rbuff(16) = er
    
      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

    else
    
      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)
      
      problem_name = cbuff(1)   
      run_id       = cbuff(2)   

      dl  = rbuff(1)
      vxl = rbuff(2)
      vyl = rbuff(3)
      vzl = rbuff(4)
      bxl = rbuff(5)
      byl = rbuff(6)
      bzl = rbuff(7)
      el  = rbuff(8)
      dr  = rbuff(9)
      vxr = rbuff(10)
      vyr = rbuff(11)
      vzr = rbuff(12)
      bxr = rbuff(13)
      byr = rbuff(14)
      bzr = rbuff(15)
      er  = rbuff(16) 
    
    endif

  end subroutine read_problem_par

!-----------------------------------------------------------------------------

  subroutine init_prob

    implicit none
 
    integer i,j,k
    real :: xi,yj,zk
    real :: vx,vy,vz,rho,pre,bx,by,bz
    real :: kn,Lx,kJ,Ly,Lz,Ln
    
    call read_problem_par

!   Secondary parameters


    do j = 1,ny
      yj = y(j)
      do i = 1,nx
        xi = x(i)
        do k = 1,nz
          zk = z(k)

          if ((xi <= 0.5)) then
             rho = dl
             pre = el
             vx  = vxl
             vy  = vyl
             vz  = vzl
             bx  = bxl
             by  = byl
             bz  = bzl
          else
             rho = dr
             pre = er
             vx  = vxr
             vy  = vyr
             vz  = vzr
             bx  = bxr
             by  = byr
             bz  = bzr
          endif
          
        
          
          u(idna,i,j,k) = rho
                          
          u(imxa,i,j,k) = vx*u(idna,i,j,k)
          u(imya,i,j,k) = vy*u(idna,i,j,k)
          u(imza,i,j,k) = vz*u(idna,i,j,k)
#ifndef ISO	  
          u(iena,i,j,k) = pre ! pre here meand eint
          u(iena,i,j,k) = max(u(iena,i,j,k), smallei)
          u(iena,i,j,k) = u(iena,i,j,k) +0.5*(vx**2+vy**2+vz**2)*u(idna,i,j,k)
#endif /* ISO */
          b(1,i,j,k)   =  bx
          b(2,i,j,k)   =  by
          b(3,i,j,k)   =  bz

#ifndef ISO	  
          u(iena,i,j,k)   = u(iena,i,j,k) +0.5*sum(b(:,i,j,k)**2,1)
#endif /* ISO */
        enddo
      enddo
    enddo
    write(*,*) maxval(b(3,:,:,:)), minval(b(3,:,:,:))
    return
  end subroutine init_prob  
  

end module init_problem

