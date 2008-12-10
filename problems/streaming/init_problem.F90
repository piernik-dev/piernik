#include "piernik.def"


module init_problem
  
! Initial condition for Keplerian disk
! Written by: M. Hanasz, March 2006

  use arrays, only : x, u, b,x,y,z,nx,ny,nz, &
      idna,imxa,imya,imza
#ifndef ISO
  use arrays, only : iena,fmagn,fadiab
#endif /* ISO */
  use start,  only : qshear, omega, proc, smallei, smalld, gamma, &
      rbuff, cbuff, ibuff, c_si, collfaq
  use mpi_setup
!  use grid

  real :: dgas,taus,epsil,kmod,eta,xnfaq,ynfaq
  character ::  problem_name*32,run_id*3,dir*1

  namelist /PROBLEM_CONTROL/  problem_name, run_id, dgas, taus, epsil, kmod, eta, xnfaq, ynfaq


contains

!-----------------------------------------------------------------------------

  subroutine read_problem_par

    implicit none
  

    problem_name = 'streaming'
    run_id  = 'tst'
    dgas    = 1.0
    taus    = 0.1
    epsil   = 3.0
    kmod    = 30.
    eta     = 0.005
    xnfaq   = 1.0
    ynfaq   = 1.0
    
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

      rbuff(1) = dgas
      rbuff(2) = taus
      rbuff(3) = epsil
      rbuff(4) = kmod
      rbuff(5) = eta
      rbuff(6) = xnfaq
      rbuff(7) = ynfaq
    
      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

    else
    
      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)
      
      problem_name = cbuff(1)   
      run_id       = cbuff(2)   

      dgas         = rbuff(1)
      taus         = rbuff(2)
      epsil        = rbuff(3)
      kmod         = rbuff(4)  
      eta          = rbuff(5)
      xnfaq        = rbuff(6)
      ynfaq        = rbuff(7)
    
    endif

  end subroutine read_problem_par

!-----------------------------------------------------------------------------

  subroutine init_prob

    implicit none
 
    integer i,j,k
    real :: xi,yj,zk
    real :: vx,vy,vz
    real :: kn,Lx,kJ,Ly,Lz,Ln
    
    real :: comfrac,vk,tauf,ux,uy,uz,wx,wy,wz,som,eom,kx,kz,r0
    real,dimension(2,4,2) :: afluct
    real,dimension(2,4)   :: fluct
    
    call read_problem_par

!   Secondary parameters


    comfrac = taus/((1.+epsil)**2+taus**2)
    vk = c_si/0.1
    r0 = vk/omega
    tauf = taus/omega
    kx = kmod/eta/r0
    kz = kmod/eta/r0
#ifdef COLLS_PART
    collfaq = 1./tauf
#else /* COLLS_PART */
    collfaq = 1./tauf/dgas
#endif /* COLLS_PART */
    if(run_id .eq. 'lnA') then
    afluct(1,2,1) =-0.1691398*eta*vk
    afluct(1,2,2) = 0.0361553*eta*vk
    afluct(1,3,1) = 0.1336704*eta*vk
    afluct(1,3,2) = 0.0591695*eta*vk
    afluct(1,4,1) = 0.1691389*eta*vk
    afluct(1,4,2) =-0.0361555*eta*vk
    afluct(1,1,1) = 0.0000224*dgas
    afluct(1,1,2) = 0.0000212*dgas
    afluct(2,1,1) = 1.e-6*epsil*dgas
    afluct(2,1,2) = 1.e-6*epsil*dgas
    afluct(2,2,1) =-0.1398623*eta*vk
    afluct(2,2,2) = 0.0372951*eta*vk
    afluct(2,3,1) = 0.1305628*eta*vk
    afluct(2,3,2) = 0.0640574*eta*vk
    afluct(2,4,1) = 0.1639549*eta*vk
    afluct(2,4,2) =-0.0233277*eta*vk
    som =-0.3480127*omega
    eom = 0.4190204*omega
    endif
    if(run_id .eq. 'lnB') then
    afluct(1,2,1) =-0.0174121*eta*vk
    afluct(1,2,2) =-0.2770347*eta*vk
    afluct(1,3,1) = 0.2767976*eta*vk
    afluct(1,3,2) =-0.0187568*eta*vk
    afluct(1,4,1) = 0.0174130*eta*vk
    afluct(1,4,2) = 0.2770423*eta*vk
    afluct(1,1,1) =-0.0000067*dgas
    afluct(1,1,2) =-0.0000691*dgas
    afluct(2,1,1) = 1.e-6*epsil*dgas
    afluct(2,1,2) = 1.e-6*epsil*dgas
    afluct(2,2,1) = 0.0462916*eta*vk
    afluct(2,2,2) =-0.2743072*eta*vk
    afluct(2,3,1) = 0.2739304*eta*vk
    afluct(2,3,2) = 0.0039293*eta*vk
    afluct(2,4,1) = 0.0083263*eta*vk
    afluct(2,4,2) = 0.2768866*eta*vk
    som = 0.4998786*omega
    eom = 0.0154764*omega
    endif
    
    ! gas mean velocities
    ux = xnfaq * 2.*epsil*comfrac*eta*vk
    uy = ynfaq * (-(1.+epsil*taus*comfrac)*eta*vk/(1.+epsil))
    uz = 0.0
    ! dust mean velocities
    wx = xnfaq * (-2.*comfrac*eta*vk)
    wy = ynfaq * (-(1.-taus*comfrac)*eta*vk/(1.+epsil))
    wz = 0.0
    write(*,*) 'Mean gas x,y velocities:'
    write(*,*) ux,uy
    write(*,*) 'Mean dust x,y velocities:'
    write(*,*) wx,wy
	  
    do j = 1,ny
      yj = y(j)
      do i = 1,nx
        xi = x(i)+r0
        do k = 1,nz
          zk = z(k)
          vx = 0.0
          vy = 0.0 !-qshear*omega*xi
          vz = 0.0
	  
	  ! fluctuations
	  fluct(:,1:3) = (afluct(:,1:3,1)*cos(kx*xi)-afluct(:,1:3,2)*sin(kx*xi))*cos(kz*zk)
	  fluct(:,  4) = (afluct(:,  4,1)*sin(kx*xi)-afluct(:,  4,2)*cos(kx*xi))*sin(kz*zk)
          
	  u(idna(1),i,j,k) = dgas+fluct(1,1)
          u(imxa(1),i,j,k) = u(idna(1),i,j,k)*(ux+fluct(1,2))
          u(imya(1),i,j,k) = u(idna(1),i,j,k)*(uy+fluct(1,3))
          u(imza(1),i,j,k) = u(idna(1),i,j,k)*(uz+fluct(1,4))
	  u(idna(2),i,j,k) = epsil*u(idna(1),i,j,k)+fluct(2,1)
          u(imxa(2),i,j,k) = u(idna(2),i,j,k)*(wx+fluct(2,2))
          u(imya(2),i,j,k) = u(idna(2),i,j,k)*(wy+fluct(2,3))
          u(imza(2),i,j,k) = u(idna(2),i,j,k)*(wz+fluct(2,4))
#ifndef ISO	  
          u(iena,i,j,k) = c_si/(gamma(fadiab)-1.0)*u(idna(fadiab),i,j,k)
          u(iena,i,j,k) = max(u(iena,i,j,k), smallei)
          u(iena,i,j,k) = u(iena,i,j,k) +0.5*(u(imxa(fadiab),i,j,k)**2 &
	                                     +u(imya(fadiab),i,j,k)**2 &
					     +u(imza(fadiab),i,j,k)**2)/u(idna(fadiab),i,j,k)
#endif /* ISO */
          b(1,i,j,k)   =  0.0
          b(2,i,j,k)   =  0.0
          b(3,i,j,k)   =  0.0

#ifndef ISO	  
          u(iena(fmagn),i,j,k)   = u(iena(fmagn),i,j,k) +0.5*sum(b(:,i,j,k)**2,1)
#endif /* ISO */
        enddo
      enddo
    enddo
    return
  end subroutine init_prob  
  

end module init_problem

