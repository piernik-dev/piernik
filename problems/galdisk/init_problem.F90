#include "mhd.def"

module init_problem
  
! Initial condition for Galactic disk
! Written by: D. Wolt, June 2007

  use arrays
  use start
  use grid
  use hydrostatic
  use gravity
  use constants
#ifdef SNE_DISTR
  use sn_distr
#endif SNE_DISTR

  real d0, r_max, rhoa
  integer mtr
  character problem_name*32,run_id*3
  character mf_orient*32

  namelist /PROBLEM_CONTROL/  problem_name, run_id, &
                              rhoa,d0,r_max,mtr, mf_orient


contains

!-----------------------------------------------------------------------------

  subroutine read_problem_par

    implicit none
  

    problem_name = 'aaa'
    run_id  = 'aa'
    rhoa    = 1.0e-4
    d0      = 1.0
    r_max   = 0.8
    mtr      = 10
    mf_orient = 'null' ! 'toroidal', 'vertical'
         
    
    if(proc .eq. 0) then
      open(1,file='problem.par')
        read(unit=1,nml=PROBLEM_CONTROL)
      close(1)
        write(*,nml=PROBLEM_CONTROL)
      open(3, file='tmp.log', position='append')
        write(3,nml=PROBLEM_CONTROL)
        write(3,*)
      close(3)
    endif

    if(proc .eq. 0) then


      cbuff(1) =  problem_name
      cbuff(2) =  run_id
      cbuff(3) =  mf_orient

      rbuff(1) = rhoa
      rbuff(2) = d0
      rbuff(3) = r_max

      ibuff(1) = mtr
    
      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

    else
    
      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)
      
      problem_name = cbuff(1)   
      run_id       = cbuff(2)   
      mf_orient	   = cbuff(3)

      rhoa         = rbuff(1)  
      d0           = rbuff(2)  
      r_max	   = rbuff(3)
    
      mtr          = ibuff(1)
    
    endif

  end subroutine read_problem_par

!-----------------------------------------------------------------------------

  subroutine init_prob

    implicit none
 
    integer i,j,k
    real xi,yj,zk, rc, vx, vy, vz,rs
    real, allocatable ::dprof(:)
    real iOmega
    real alfar, dens0
    real dnmol, dncold, dnwarm, dnion, dnhot
    real dcolumn, dcolsmall, cd, cdprevious, d0previous, afactor, bfactor, dcolumnprevious
    integer iter, itermax, itermaxwrite, inzfac
    real xgradgp, xgradp
    character syscmd*36,syscmd2*40,syscmd3*37
    integer system, syslog

    call read_problem_par

    allocate(dprof(nz))

    if(proc .eq. 0) write(*,*) 'first loop: calculating dcolumn, dprof'

    if(proc .eq. 0) then
      write(syscmd,'(a33,i3.3)') "echo -n 'first loop: j = 001 of '",ny
      syslog=SYSTEM(syscmd)
    endif
    do j = 1,ny
      write(syscmd2,'(a30,i3.3,a4,i3.3)') "echo -n '\b\b\b\b\b\b\b\b\b\b'",j," of ",ny
      syslog=SYSTEM(syscmd2)
      yj = y(j)
      do i = 1,nx
        xi = x(i)
        rc = sqrt(xi**2+yj**2)
	dcolumn = 1.0e-6
        dcolumnprevious = 0.0
        inzfac = 0
        do while((dcolumn-dcolumnprevious).ge.(1.0e-6*dcolumn))
          inzfac = inzfac + 1
          dcolumnprevious = dcolumn
          dcolumn = 0.0
          do k = 1,2*inzfac*nz
            zk = inzfac*zmin + 0.5*dl(zdim) + k*dl(zdim)
            dnmol = 0.58/(cm**3) * exp(-((rc - 4.5*kpc)**2-(r_gc_sun - 4.5*kpc)**2)/(2.9*kpc)**2) &
         		& * (rc/r_gc_sun)**(-0.58) * exp(-(zk/(81.0*(rc/r_gc_sun)**(0.58)))**2)
	    if(rc.lt.r_gc_sun) then
	      alfar=1.0
	    else
	      alfar=rc/r_gc_sun
	    endif
	    dncold= 0.34/(cm**3)/alfar**2 * (0.859*exp(-(zk/(127.0*alfar))**2) &
	                                 & + 0.047*exp(-(zk/(318.0*alfar))**2) &
	    		                 & + 0.094*exp(-(abs(zk)/(403.0*alfar))))
	    dnwarm= 0.226/(cm**3)/alfar * ((1.745 - 1.289/alfar)*exp(-(zk/(127.0*alfar))**2) &
	    		               & + (0.473 -  0.07/alfar)*exp(-(zk/(318.0*alfar))**2) &
			               & + (0.283 - 0.142/alfar)*exp(-(abs(zk)/(403.0*alfar))))
	    dnion = 0.0237/(cm**3)*exp(-(rc**2 - r_gc_sun**2)/(37.0*kpc)**2)*exp(-abs(zk)/kpc) &
		& + 0.0013/(cm**3)*exp(-((rc - 4.0*kpc)**2 &
		                     & - (r_gc_sun - 4.0*kpc)**2)/(2.0*kpc)**2) &
			             & *exp(-abs(zk)/150.0/pc)
	    dnhot = 4.8e-4/(cm**3)*(0.12*exp(-(rc - r_gc_sun)/(4.9*kpc)) + 0.88*exp(-((rc - 4.5*kpc)**2 &
	    		& - (r_gc_sun - 4.5*kpc)**2)/(2.9*kpc)**2)) * (rc/r_gc_sun)**(-1.65)*exp(-abs(zk) &
			& /(1.5*kpc*(rc/r_gc_sun)**(1.65)))
	    dens0=1.36*mp*(dnmol+dncold+dnwarm+dnion+dnhot)
	    dcolumn=dcolumn+dens0*dl(zdim)
          enddo !k=1,2*iznfac*nz
        enddo !inzfac accuracy

        iter=0
        itermax=100
        dcolsmall=dcolumn*1.e-8
        d0 = dcolumn/dl(zdim)/nz
        call hydrostatic_zeq(i, j, d0, dprof)
        cd = 0.0
        do k=1,nz
	  cd = cd +dprof(k)*nz
        enddo
        do while((abs(cd - dcolumn) .gt. dcolsmall).and.(iter .le. itermax))
          if(iter .eq. 0) then
            d0previous = d0
            cdprevious = cd
            d0 = d0*2.
          else
            afactor = (cd - cdprevious)/(d0 - d0previous)
            bfactor = cd - afactor*d0
            d0previous = d0
            cdprevious = cd
            d0 = (dcolumn - bfactor)/afactor
          endif
          iter = iter+1
	  itermaxwrite = max(itermaxwrite, iter)
	  call hydrostatic_zeq(i, j, d0, dprof)
	  cd = 0.0
	  do k=1,nz
	    cd = cd + dprof(k)*nz
	  enddo
        enddo
	if(abs(cd-dcolumn).gt.dcolsmall) then
	    write(*,*) proc,i,j,'equatorial density accuracy different than required!'
	endif

        do k=1,nz
	  u(1,i,j,k)=rhoa + dprof(k)/cosh((rc/r_max)**mtr)
	  u(1,i,j,k)=max(u(1,i,j,k),smalld)
	enddo
      enddo
    enddo
    if(proc .eq. 0) then
      write(*,*)
      write(*,*) 'second loop: calculating iOmega, u(:,:,:)'
    endif
    if(proc .eq. 0) then
      write(syscmd3,'(a34,i3.3)') "echo -n 'second loop: j = 001 of '",ny
      syslog=SYSTEM(syscmd3)
    endif
    do j = 1,ny
      write(syscmd2,'(a30,i3.3,a4,i3.3)') "echo -n '\b\b\b\b\b\b\b\b\b\b'",j," of ",ny
      syslog=SYSTEM(syscmd2)
      yj = y(j)
      do i = 1,nx
	xi = x(i)
        rc = sqrt(xi**2+yj**2)
	vz = 0.0
	do k = 1,nz
	 zk=z(k)
	 rs = sqrt(xi**2+yj**2+zk**2)
!         u(1,i,j,k) = rhoa + dprof(k)/cosh((rc/r_max)**mtr)
!         u(1,i,j,k) = max(u(1,i,j,k), smalld)
	 u(2,i,j,k) = vx*u(1,i,j,k)
         u(3,i,j,k) = vy*u(1,i,j,k)
           if(i .ne. 1 .and. i .ne. nx) then
	     xgradgp=(gp(i+1,j,k)-gp(i-1,j,k))/2.0/dl(xdim)
	     xgradp =-0.5*c_si**2/gamma/u(1,i,j,k)*(u(1,i+1,j,k)-u(1,i-1,j,k))/dl(xdim)
           endif
           if(i .eq. 1) then
	     xgradgp=(gp(i+1,j,k)-gp(i,j,k))/dl(xdim)
	     xgradp =-c_si**2/gamma/u(1,i,j,k)*(u(1,i+1,j,k)-u(1,i,j,k))/dl(xdim)
           endif
           if(i .eq. nx) then
	     xgradgp=(gp(i,j,k)-gp(i-1,j,k))/dl(xdim)
	     xgradp =-c_si**2/gamma/u(1,i,j,k)*(u(1,i,j,k)-u(1,i-1,j,k))/dl(xdim)
           endif
           iOmega=sqrt(abs(xgradgp+xgradp)/abs(xi))
	     u(2,i,j,k)=-iOmega*yj*u(1,i,j,k)
             u(3,i,j,k)= iOmega*xi*u(1,i,j,k)
         u(4,i,j,k) = vz*u(1,i,j,k)
#ifndef ISO
         u(5,i,j,k) = c_si**2/(gamma-1.0)*u(1,i,j,k)
         u(5,i,j,k) = max(u(5,i,j,k), smallei)
	 u(5,i,j,k) = u(5,i,j,k) +0.5*(vx**2+vy**2+vz**2)*u(1,i,j,k)
#endif
         select case(mf_orient)
	 case('null')
         b(1,i,j,k)   = 0.0
         b(2,i,j,k)   = 0.0
         b(3,i,j,k)   = 0.0
         case('vertical')
         b(1,i,j,k)   = 0.0
         b(2,i,j,k)   = 0.0
         b(3,i,j,k)   = sqrt(2.*alpha*d0*c_si**2)
	 case('toroidal')
         b(1,i,j,k)   =-sqrt(2.*alpha*c_si**2*u(idna,i,j,k))*yj/rc
         b(2,i,j,k)   = sqrt(2.*alpha*c_si**2*u(idna,i,j,k))*xi/rc
         b(3,i,j,k)   = 0.0
	 end select
#ifndef ISO
         u(5,i,j,k)   = u(5,i,j,k) +0.5*sum(b(:,i,j,k)**2,1)
#endif
        enddo
      enddo
    enddo
      write(*,*) 'the longest iteration number of steps = ', itermaxwrite
    if(allocated(dprof)) deallocate(dprof)



#ifdef SNE_DISTR
   call compute_SNdistr
#endif SNE_DISTR

    return
  end subroutine init_prob  

!-----------------------------------------------------------------------------------------------------------------------------------

end module init_problem

