#include "mhd.def"

module init_problem
  
! Initial condition for Galactic disk with simple bulge
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
  integer mtr, przypadek
  character problem_name*32,run_id*3
  character bulge*4, rotsc*4, mf_orient*32

  namelist /PROBLEM_CONTROL/  problem_name, run_id, &
                              rhoa,d0,r_max,mtr, przypadek, bulge, rotsc, mf_orient


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
    przypadek= 2
    bulge    = 'null'  ! 'smpl'
    rotsc    = 'dcol'  ! 'dens', 'auto'
    mf_orient = 'null' ! 'toroidal', 'vertical'
         
    
    if(proc .eq. 0) then
      open(1,file='problem.par')
        read(unit=1,nml=PROBLEM_CONTROL)
      close(1)
#ifdef FILEINFO
      open(99, file='simulation.out', position='append')
        write(99,nml=PROBLEM_CONTROL)
      close(99)
#else FILEINFO
        write(*,nml=PROBLEM_CONTROL)
#endif FILEINFO
      open(3, file='tmp.log', position='append')
        write(3,nml=PROBLEM_CONTROL)
        write(3,*)
      close(3)
    endif

    if(proc .eq. 0) then


      cbuff(1) =  problem_name
      cbuff(2) =  run_id
      cbuff(3) =  bulge
      cbuff(4) =  rotsc
      cbuff(5) =  mf_orient

      rbuff(1) = rhoa
      rbuff(2) = d0
      rbuff(3) = r_max

      ibuff(1) = mtr
      ibuff(2) = przypadek
    
      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

    else
    
      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)
      
      problem_name = cbuff(1)   
      run_id       = cbuff(2)   
      bulge	   = cbuff(3)
      rotsc	   = cbuff(4)
      mf_orient	   = cbuff(5)

      rhoa         = rbuff(1)  
      d0           = rbuff(2)  
      r_max	   = rbuff(3)
    
      mtr          = ibuff(1)
      przypadek	   = ibuff(2)
    
    endif

  end subroutine read_problem_par

!-----------------------------------------------------------------------------

  subroutine init_prob

    implicit none
 
    integer i,j,k
    real xi,yj,zk, rc, vx, vy, vz,h2,dgdz,csim2,rs
    real, allocatable ::dprof(:),gravaccel(:,:,:,:)
    real ivsun, iOmega
    real alfar, dens0
    real dnmol, dncold, dnwarm, dnion, dnhot
    real densmax, densmin, densmaxall, densminall, maxdcolumn, mindcolumn, maxdcradius, mindcradius
    real dcolumn, dcolsmall, cd, cdprevious, d0previous, afactor, bfactor, dcolumnprevious
    integer iter, itermax, itermaxwrite, inzfac
    real densdiscmean, adensdiscmean,ddensdiscmean
    integer idensdiscmean, actproc
    real rotscalh, maxrotscalh, ivarotscalh
    real xgradgp, xgradp
    real, dimension(nx,ny) :: dcolumnarray
    real, allocatable :: rotscalharray(:,:)
    integer, allocatable :: arotscalh(:,:,:)
    integer ikk,jkk,kkk,ip
#ifdef FILEINFO
    character*14 simufile
    character syscmd*54,syscmd2*58,syscmd3*55
#else FILEINFO
    character syscmd*36,syscmd2*40,syscmd3*37
#endif FILEINFO
#ifdef ROTINFO
    character*14 omegatablproc
#endif ROTINFO
    character*44 system_command
    integer system_status, system, syslog

    call read_problem_par
#ifdef ROTINFO
    write(omegatablproc,'(a9,i1.1,a4)') 'omegatabl',proc,'.out'
    open(195,file=omegatablproc,status='unknown')
      write(195,*) ' '
    close(195)
#endif ROTINFO
#ifdef FILEINFO
    write(simufile,'(a9,i1.1,a4)') 'simulatio',proc,'.out'
    if(proc .eq. 0) simufile = 'simulation.out'
    write(system_command,'(a30,a14)') 'echo "Computing on " $HOST >> ',simufile
    system_status=SYSTEM(system_command)
#endif FILEINFO

    allocate(dprof(nz))
    if(rotsc .ne. 'auto') then  !========================================!
      allocate(arotscalh(psize(1),psize(2),psize(3)))                    !
      allocate(rotscalharray(nx,ny))                                     !
      ivsun=vsun                                                         !   
      maxrotscalh = 0.0                                                  !
    endif !==(rotsc .ne. 'auto')=========================================!

    !===diagnostyczne===!
    itermaxwrite = 1    !
    maxdcolumn = 0.0    !
    mindcolumn = 1.0e30 !
    densdiscmean = 0.0  !
    adensdiscmean = 0.0 !
    ddensdiscmean = 0.0 !
    idensdiscmean = 0   !
    !===================!

#ifdef FILEINFO
    open(99, file=simufile, position='append')  !=================================!
      write(99,*) 'The process number and block coordinates are:'		  !
      write(99,*) proc,pcoords                                                    !
      write(99,*) 'first loop: calculating dcolumn, dprof, rotscalh, maxrotscalh' !
    close(99) !===================================================================!
#else FILEINFO
    if(proc .eq. 0) write(*,*) 'first loop: calculating dcolumn, dprof, rotscalh, maxrotscalh'
#endif FILEINFO

!================poczatek petli I-go stopnia===========================================================================!
#ifdef FILEINFO
    write(syscmd,'(a33,i3.3,a4,a14)') "echo -n 'first loop: j = 001 of '",ny," >> ",simufile			       !
    syslog=SYSTEM(syscmd)											       !
#else FILEINFO
    if(proc .eq. 0) then
      write(syscmd,'(a33,i3.3)') "echo -n 'first loop: j = 001 of '",ny						       !
      syslog=SYSTEM(syscmd)											       !
    endif
#endif FILEINFO
    do j = 1,ny  !==========sumowanie dens i dcol po y-ku=============================================================!!
#ifdef FILEINFO
      write(syscmd2,'(a30,i3.3,a4,i3.3,a4,a14)') "echo -n '\b\b\b\b\b\b\b\b\b\b'",j," of ",ny," >> ",simufile	      !!
      syslog=SYSTEM(syscmd2)											      !!
#else FILEINFO
      write(syscmd2,'(a30,i3.3,a4,i3.3)') "echo -n '\b\b\b\b\b\b\b\b\b\b'",j," of ",ny				      !!
      syslog=SYSTEM(syscmd2)											      !!
#endif FILEINFO
      yj = y(j)													      !!
      do i = 1,nx  !=========sumowanie dens i dcol po x-ie===========================================================!!!
        xi = x(i)												     !!!
        rc = sqrt(xi**2+yj**2)											     !!!
        if(rotsc .ne. 'auto') then  !===========================================================!		     !!!
!	  if((grav_model .eq. 'galsmooth') .or. (grav_model .eq. 'gallksun')) then  !===!	!		     !!!
            iOmega=vsun*tanh(rc/3.0/kpc)/rc						!	!		     !!!
!         endif  !======================================================================!	!		     !!!
          vx = -iOmega * yj  ! * rc/rc								!		     !!!
          vy =  iOmega * xi  ! * rc/rc								!		     !!!
          vz = 0.0										!		     !!!
        endif !==(rotsc .ne. 'auto')============================================================!		     !!!
        !========starting column density calculation===============================================================! !!!
	dcolumn = 1.0e-6											   ! !!!
        dcolumnprevious = 0.0											   ! !!!
        inzfac = 0												   ! !!!
        do while((dcolumn-dcolumnprevious).ge.(1.0e-6*dcolumn))   !=====petla do precyzji sumowania===============!! !!!
          inzfac = inzfac + 1											  !! !!!
          dcolumnprevious = dcolumn										  !! !!!
          dcolumn = 0.0												  !! !!!
          do k = 1,2*inzfac*nz  !===sumowanie dens i dcol po z-cie===============================================!!! !!!
            zk = inzfac*zmin + 0.5*dl(zdim) + k*dl(zdim)							 !!! !!!
            !========================calculating of column density==============================================!!!! !!!
            dnmol = 0.58/(cm**3) * exp(-((rc - 4.5*kpc)**2-(r_gc_sun - 4.5*kpc)**2)/(2.9*kpc)**2) &		!!!! !!!
         		& * (rc/r_gc_sun)**(-0.58) * exp(-(zk/(81.0*(rc/r_gc_sun)**(0.58)))**2)			!!!! !!!
	    if(rc.lt.r_gc_sun) then  !==!									!!!! !!!
	      alfar=1.0			!									!!!! !!!
	    else			!									!!!! !!!
	      alfar=rc/r_gc_sun		!									!!!! !!!
	    endif  !====================!									!!!! !!!
	    dncold= 0.34/(cm**3)/alfar**2 * (0.859*exp(-(zk/(127.0*alfar))**2) &				!!!! !!!
	                                 & + 0.047*exp(-(zk/(318.0*alfar))**2) &				!!!! !!!
	    		                 & + 0.094*exp(-(abs(zk)/(403.0*alfar))))				!!!! !!!
	    dnwarm= 0.226/(cm**3)/alfar * ((1.745 - 1.289/alfar)*exp(-(zk/(127.0*alfar))**2) &			!!!! !!!
	    		               & + (0.473 -  0.07/alfar)*exp(-(zk/(318.0*alfar))**2) &			!!!! !!!
			               & + (0.283 - 0.142/alfar)*exp(-(abs(zk)/(403.0*alfar))))			!!!! !!!
	    dnion = 0.0237/(cm**3)*exp(-(rc**2 - r_gc_sun**2)/(37.0*kpc)**2)*exp(-abs(zk)/kpc) &		!!!! !!!
		& + 0.0013/(cm**3)*exp(-((rc - 4.0*kpc)**2 & 							!!!! !!!
		                     & - (r_gc_sun - 4.0*kpc)**2)/(2.0*kpc)**2) &				!!!! !!!
			             & *exp(-abs(zk)/150.0/pc)							!!!! !!!
	    dnhot = 4.8e-4/(cm**3)*(0.12*exp(-(rc - r_gc_sun)/(4.9*kpc)) + 0.88*exp(-((rc - 4.5*kpc)**2 &	!!!! !!!
	    		& - (r_gc_sun - 4.5*kpc)**2)/(2.9*kpc)**2)) * (rc/r_gc_sun)**(-1.65)*exp(-abs(zk) &	!!!! !!!
			& /(1.5*kpc*(rc/r_gc_sun)**(1.65)))							!!!! !!!
	    dens0=1.36*mp*(dnmol+dncold+dnwarm+dnion+dnhot)							!!!! !!!
	    dcolumn=dcolumn+dens0*dl(zdim)  !===================================================================!!!! !!!
          enddo !k=1,2*iznfac*nz   !=======sumowanie dens i dcol po z-cie========================================!!! !!!
        enddo !inzfac accuracy    !=====petla do w celu osiagniecia precyzji sumowania============================!! !!!
        !=============adding simple bulge=================================================!			   ! !!!
        if(bulge .eq. 'smpl') dcolumn = dcolumn + 30.0*Msun/pc**2/(cosh(rc/3.0/kpc))**10  !			   ! !!!
        !=================================================================================!    			   ! !!!
	dcolumnarray(i,j) = dcolumn										   ! !!!
														   ! !!!
        !===diagnostyczne================!									   ! !!!
        if(mindcolumn .gt. dcolumn) then !									   ! !!!
          mindcradius = rc		 !									   ! !!!
          mindcolumn = dcolumn		 !									   ! !!!
        endif  !=========================!									   ! !!!
        if(maxdcolumn .lt. dcolumn) then !									   ! !!!
          maxdcradius = rc		 !									   ! !!!
          maxdcolumn = dcolumn		 !									   ! !!!
        endif  !=========================!									   ! !!!
														   ! !!!
        !========================end of column density calculation=================================================! !!!
														     !!!
        !======================dprof(k) first calculation========================================================!   !!!
        iter=0													 !   !!!
        itermax=100												 !   !!!
        dcolsmall=dcolumn*1.e-8											 !   !!!
        d0 = dcolumn/dl(zdim)/nz											 !   !!!
        call hydrostatic_zeq(i, j, d0, dprof)									 !   !!!
        cd = 0.0												 !   !!!
        do k=1,nz  !============!										 !   !!!
	  cd = cd +dprof(k)*nz	!										 !   !!!
        enddo  !================!										 !   !!!
        do while((abs(cd - dcolumn) .gt. dcolsmall).and.(iter .le. itermax))  !==obliczanie iteracyjne dcolumn==!!   !!!
          if(iter .eq. 0) then  !===============================!						!!   !!!
            d0previous = d0					!						!!   !!!
            cdprevious = cd					!						!!   !!!
            d0 = d0*2.						!						!!   !!!
          else							!						!!   !!!
            afactor = (cd - cdprevious)/(d0 - d0previous)	!						!!   !!!
            bfactor = cd - afactor*d0				!						!!   !!!
            d0previous = d0					!						!!   !!!
            cdprevious = cd					!						!!   !!!
            d0 = (dcolumn - bfactor)/afactor			!						!!   !!!
          endif  !==============================================!						!!   !!!
          iter = iter+1												!!   !!!
	  itermaxwrite = max(itermaxwrite, iter)								!!   !!!
	  call hydrostatic_zeq(i, j, d0, dprof)									!!   !!!
	  cd = 0.0												!!   !!!
	  do k=1,nz  !==================!									!!   !!!
	    cd = cd + dprof(k)*nz	!									!!   !!!
	  enddo  !======================!									!!   !!!
        enddo  !==========obliczanie iteracyjne dcolumn=========================================================!!   !!!
	if(abs(cd-dcolumn).gt.dcolsmall) then  !========================================!			 !   !!!
#ifdef FILEINFO
	  open(99, file=simufile, position='append')					!			 !   !!!
	    write(99,*) i,j,'equatorial density accuracy different than required!'	!			 !   !!!
	  close(99)									!			 !   !!!
#else FILEINFO
	    write(*,*) proc,i,j,'equatorial density accuracy different than required!'	!			 !   !!!
#endif FILEINFO
	endif  !========================================================================!			 !   !!!
        !=======================end of the first dprof(k) calculation============================================!   !!!
        do k=1,nz  !============zapisanie gestosci w tablicy u==================!				     !!!
	  u(1,i,j,k)=rhoa + dprof(k)/cosh((rc/r_max)**mtr)			!				     !!!
	  u(1,i,j,k)=max(u(1,i,j,k),smalld)					!				     !!!
	enddo !=================================================================!				     !!!
      enddo  !==========sumowanie dens i dcol po x-ie================================================================!!!
    enddo  !==========sumowanie dens i dcol po y-ku===================================================================!!
!==============koniec petli I-go stopnia===============================================================================!
    if(rotsc .ne. 'auto') deallocate(arotscalh)
#ifdef FILEINFO
    open(99, file=simufile, position='append')  !===============!
      write(99,*)						!
      write(99,*) 'second loop: calculating iOmega, u(:,:,:)'	!
    close(99)  !================================================!
#else FILEINFO
    if(proc .eq. 0) then  !=====================================!
      write(*,*)						!
      write(*,*) 'second loop: calculating iOmega, u(:,:,:)'	!
    endif  !====================================================!
#endif FILEINFO
#ifdef ROTINFO
    if(rotsc .ne. 'auto') open(88, file='rotscal.out', status='unknown')
#endif ROTINFO
    if(rotsc .eq. 'dcol') then  !=======rotscalh(dcol) calculation======!
      rotscalharray=dcolumnarray					!
      maxrotscalh = maxval(dcolumnarray(:,:))				!
    endif  !=====================end of rotscalh(dcol) calculation======!
!==========poczatek petli I-go stopnia==============================================================================!
#ifdef FILEINFO
    write(syscmd3,'(a34,i3.3,a4,a14)') "echo -n 'second loop: j = 001 of '",ny," >> ",simufile			    !
    syslog=SYSTEM(syscmd3)											    !
#else FILEINFO
    if(proc .eq. 0) then
      write(syscmd3,'(a34,i3.3)') "echo -n 'second loop: j = 001 of '",ny					    !
      syslog=SYSTEM(syscmd3)											    !
    endif
#endif FILEINFO
    do j = 1,ny   !========przypisywanie wartosci tablicom u i b po y-ku===========================================!!
#ifdef FILEINFO
      write(syscmd2,'(a30,i3.3,a4,i3.3,a4,a14)') "echo -n '\b\b\b\b\b\b\b\b\b\b'",j," of ",ny," >> ",simufile	   !!
      syslog=SYSTEM(syscmd2)											   !!
#else FILEINFO
      write(syscmd2,'(a30,i3.3,a4,i3.3)') "echo -n '\b\b\b\b\b\b\b\b\b\b'",j," of ",ny				   !!
      syslog=SYSTEM(syscmd2)											   !!
#endif FILEINFO
      yj = y(j)													   !!
      do i = 1,nx   !===========przypisywanie wartosci tablicom u i b po x-ie=====================================!!!
	xi = x(i)												  !!!
        rc = sqrt(xi**2+yj**2)											  !!!
        if(rotsc .ne. 'auto') then  !====================================================!			  !!!
          if((grav_model .eq. 'galsmooth') .or. (grav_model .eq. 'gallksun')) then  !===!!			  !!!
            iOmega=vsun*tanh(rc/3.0/kpc)/rc						!!			  !!!
          endif  !======================================================================!!			  !!!
	  vx = -iOmega * yj  ! * rc/rc							 !			  !!!
	  vy =  iOmega * xi  ! * rc/rc							 !			  !!!
        endif !==(rotsc .ne. 'auto')=====================================================!			  !!!
	vz = 0.0												  !!!
#ifdef ROTINFO
        if(rotsc .ne. 'auto') write(88,*) rc,rotscalharray(i,j),maxrotscalh					  !!!
#endif ROTINFO
	do k = 1,nz    !===========przypisywanie wartosci tablicom u i b po z-cie================================!!!!
	 zk=z(k)												 !!!!
	 rs = sqrt(xi**2+yj**2+zk**2)										 !!!!
         if(grav_model .eq. 'galsmtdec') then  !================================================================!!!!!
           iOmega=vsun*tanh(rc/3.0/kpc)/rc*exp(-zk**2/2.0/(h_grav*rotscalharray(i,j)/maxrotscalh)**n_gravh)	!!!!!
	   vx = -iOmega * yj  ! * rc/rc										!!!!!
	   vy =  iOmega * xi  ! * rc/rc										!!!!!
	   vz = 0.0												!!!!!
         endif  !===============================================================================================!!!!!
!         u(1,i,j,k) = rhoa + dprof(k)/cosh((rc/r_max)**mtr)							 !!!!
!         u(1,i,j,k) = max(u(1,i,j,k), smalld)									 !!!!
	 u(2,i,j,k) = vx*u(1,i,j,k)										 !!!!
         u(3,i,j,k) = vy*u(1,i,j,k)										 !!!!
	 if(rotsc .eq. 'auto') then  !===================================================!			 !!!!
           if(i .ne. 1 .and. i .ne. nx) then  !=========================================!!			 !!!!
	     xgradgp=(gp(i+1,j,k)-gp(i-1,j,k))/2.0/dl(xdim)				!!			 !!!!
	     xgradp =-0.5*c_si**2/gamma/u(1,i,j,k)*(u(1,i+1,j,k)-u(1,i-1,j,k))/dl(xdim)	!!			 !!!!
           endif  !=====================================================================!!			 !!!!
           if(i .eq. 1) then  !=========================================================!!			 !!!!
	     xgradgp=(gp(i+1,j,k)-gp(i,j,k))/dl(xdim)					!!			 !!!!
	     xgradp =-c_si**2/gamma/u(1,i,j,k)*(u(1,i+1,j,k)-u(1,i,j,k))/dl(xdim)	!!			 !!!!
           endif  !=====================================================================!!			 !!!!
           if(i .eq. nx) then  !========================================================!!			 !!!!
	     xgradgp=(gp(i,j,k)-gp(i-1,j,k))/dl(xdim)					!!			 !!!!
	     xgradp =-c_si**2/gamma/u(1,i,j,k)*(u(1,i,j,k)-u(1,i-1,j,k))/dl(xdim)	!!			 !!!!
           endif  !=====================================================================!!			 !!!!
           iOmega=sqrt(abs(xgradgp+xgradp)/abs(xi))					 !			 !!!!
!==============================niepotrzebne=============================================!!			 !!!!
!           if(j .ne. 1 .and. j .ne. ny) then  !========================================!!			 !!!!
!	     ygradgp=(gp(i,j+1,k)-gp(i,j-1,k))/2.0/dl(ydim)				!!			 !!!!
!	     ygradp =-0.5*c_si**2/gamma/u(1,i,j,k)*(u(1,i,j+1,k)-u(1,i,j-1,k))/dl(ydim)	!!			 !!!!
!           endif									!!			 !!!!
!           if(j .eq. 1) then								!!			 !!!!
!	     ygradgp=(gp(i,j+1,k)-gp(i,j,k))/dl(ydim)					!!			 !!!!
!	     ygradp =-c_si**2/gamma/u(1,i,j,k)*(u(1,i,j+1,k)-u(1,i,j,k))/dl(ydim)	!!			 !!!!
!           endif									!!			 !!!!
!           if(j .eq. ny) then								!!			 !!!!
!	     ygradgp=(gp(i,j,k)-gp(i,j-1,k))/dl(ydim)					!!			 !!!!
!	     ygradp =-c_si**2/gamma/u(1,i,j,k)*(u(1,i,j,k)-u(1,i,j-1,k))/dl(ydim)	!!			 !!!!
!           endif									!!			 !!!!
!	   if(xgradp*xi .ge. 0.0) then							!!			 !!!!
!	     psign=1.0									!!			 !!!!
!	   else										!!			 !!!!
!	     psign=-1.0									!!			 !!!!
!	   endif									!!			 !!!!
!	   rgradp =sqrt(xgradp**2+ygrad**2)						!!			 !!!!
!	   rgradgp=-sqrt(xgradgp**2+ygradgp**2)						!!			 !!!!
!=======================================================================================!!			 !!!!
	     u(2,i,j,k)=-iOmega*yj*u(1,i,j,k)						 !			 !!!!
             u(3,i,j,k)= iOmega*xi*u(1,i,j,k)						 !			 !!!!
         endif !==(rotsc .eq. 'auto')====================================================!			 !!!!
#ifdef ROTINFO
         open(195,file=omegatablproc,position='append') !=======!			 			 !!!!
           write(195,*) rc,iOmega,iOmega*rc			!			 			 !!!!
	 close(195)	!=======================================!			 			 !!!!
#endif ROTINFO
         u(4,i,j,k) = vz*u(1,i,j,k)	  									 !!!!
#ifndef ISO
         u(5,i,j,k) = c_si**2/(gamma-1.0)*u(1,i,j,k)								 !!!!
         u(5,i,j,k) = max(u(5,i,j,k), smallei)									 !!!!
	 u(5,i,j,k) = u(5,i,j,k) +0.5*(vx**2+vy**2+vz**2)*u(1,i,j,k)						 !!!!
#endif
         select case(mf_orient)
	 case('null')
         b(1,i,j,k)   = 0.0											 !!!!
         b(2,i,j,k)   = 0.0											 !!!!
         b(3,i,j,k)   = 0.0											 !!!!
         case('vertical')
         b(1,i,j,k)   = 0.0											 !!!!
         b(2,i,j,k)   = 0.0											 !!!!
         b(3,i,j,k)   = sqrt(2.*alpha*d0*c_si**2)								 !!!!
	 case('toroidal')
         b(1,i,j,k)   =-sqrt(2.*alpha*c_si**2*u(idna,i,j,k))*yj/rc						 !!!!
         b(2,i,j,k)   = sqrt(2.*alpha*c_si**2*u(idna,i,j,k))*xi/rc						 !!!!
         b(3,i,j,k)   = 0.0											 !!!!
	 end select
#ifndef ISO
         u(5,i,j,k)   = u(5,i,j,k) +0.5*sum(b(:,i,j,k)**2,1)							 !!!!
#endif
        enddo  !============przypisywanie wartosci tablicom u i b po z-cie=======================================!!!!
	ddensdiscmean = maxval(u(1,i,j,:))									  !!!
	call MPI_ALLREDUCE(ddensdiscmean, adensdiscmean, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm, ierr)		  !!!
	if(adensdiscmean .gt. smalld) then  !===========!							  !!!
	  densdiscmean = densdiscmean + adensdiscmean	!							  !!!
	  idensdiscmean = idensdiscmean + 1		!							  !!!
	endif  !========================================!							  !!!
      enddo    !===========przypisywanie wartosci tablicom u i b po x-ie==========================================!!!
    enddo   !=========przypisywanie wartosci tablicom u i b po y-ku================================================!!
!==========koniec petli I-go stopnia================================================================================!
#ifdef ROTINFO
    if(rotsc .ne. 'auto') close(88)									!
#endif ROTINFO
    !=====diagnostyczne==========================================================================!	!
    densdiscmean = densdiscmean/real(idensdiscmean)						 !	!
    if(proc .eq. 0) then  !=====================================!				 !	!
#ifdef FILEINFO
      open(99, file=simufile, position='append') !----------------------------------------------!!	!
      write(99,*) ' '						!				!!	!
      write(99,*) 'Mean density of disk = ',densdiscmean	!				!!	!
#else FILEINFO
      write(*,*) ' '						!				!!	!
      write(*,*) 'Mean density of disk = ',densdiscmean		!				!!	!
#endif FILEINFO
    endif  !====================================================!				!!	!
    densmax=maxval(u(1,:,:,:))									!!	!
    densmin=minval(u(1,:,:,:))									!!	!
    call MPI_REDUCE(densmax, densmaxall, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierr)	!!	!
    call MPI_REDUCE(densmin, densminall, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, comm, ierr)	!!	!
    if(proc .eq. 0) then  !=====================================================!		!!	!
#ifdef FILEINFO
      write(99,*) 'maxdcolumn (radius = ',maxdcradius,' ) = ',maxdcolumn	!		!!	!
      write(99,*) 'mindcolumn (radius = ',mindcradius,' ) = ',mindcolumn	!		!!	!
      write(99,*) 'maxdens = ', densmaxall					!		!!	!
      write(99,*) 'mindens = ', densminall					!		!!	!
      close(99)	!-------------------------------------------------------------------------------!!	!
#else FILEINFO
      write(*,*) 'maxdcolumn (radius = ',maxdcradius,' ) = ',maxdcolumn	!			!!	!
      write(*,*) 'mindcolumn (radius = ',mindcradius,' ) = ',mindcolumn	!			!!	!
      write(*,*) 'maxdens = ', densmaxall					!		!!	!
      write(*,*) 'mindens = ', densminall					!		!!	!
#endif FILEINFO
    endif  !====================================================================!		 !	!
#ifdef FILEINFO
    open(99, file=simufile, position='append')  !===============================!		 !	!
      write(99,*) 'the longest iteration number of steps = ', itermaxwrite	!		 !	!
    close(99)  !================================================================!		 !	!
#else FILEINFO
      write(*,*) 'the longest iteration number of steps = ', itermaxwrite	!		 !	!
#endif FILEINFO
    !============================================================================================!	!
    if(allocated(dprof)) then										!
#ifdef FILEINFO
    open(99, file=simufile, position='append')  !===============================!		 	!
      write(99,*) 'dprof is allocated, deallocating...'				!		 	!
    close(99)  !================================================================!		 	!
#else FILEINFO
      write(*,*) 'dprof is allocated, deallocating...'				!		 	!
#endif FILEINFO
      deallocate(dprof)											!
#ifdef FILEINFO
    open(99, file=simufile, position='append')  !===============================!		 	!
      write(99,*) '...deallocated.'						!		 	!
    close(99)  !================================================================!		 	!
#else FILEINFO
      write(*,*) '...deallocated.'						!		 	!
#endif FILEINFO
    else												!
#ifdef FILEINFO
    open(99, file=simufile, position='append')  !===============================!		 	!
      write(99,*) 'dprof is not allocated, do nothing.'				!		 	!
    close(99)  !================================================================!		 	!
#else FILEINFO
      write(*,*) 'dprof is not allocated, do nothing.'				!		 	!
#endif FILEINFO
    endif												!
!==========koniec dla rotsc .eq. 'auto'=================================================================!

!==========poczatek warunku I-go stopnia=================================================================================================!
    if(rotsc .ne. 'auto') then  !=======================================================================================================!!
      if(allocated(gravaccel)) then													!!
#ifdef FILEINFO
    open(99, file=simufile, position='append')  !===============================!		 
      write(99,*) 'gravaccel is allocated, do nothing.'				!		 
    close(99)  !================================================================!		 
#else FILEINFO
      write(*,*) 'gravaccel is allocated, do nothing.'				
#endif FILEINFO
      else																!!
#ifdef FILEINFO
    open(99, file=simufile, position='append')  !===============================!		
      write(99,*) 'gravaccel is not allocated, allocating...'			!		
    close(99)  !================================================================!		
#else FILEINFO
      write(*,*) 'gravaccel is not allocated, allocating...'			
#endif FILEINFO
        allocate(gravaccel(2,nx,ny,nz))													!!
#ifdef FILEINFO
    open(99, file=simufile, position='append')  !===============================!		
      write(99,*) '...allocated.'						!		
    close(99)  !================================================================!		
#else FILEINFO
      write(*,*) '...allocated.'						
#endif FILEINFO
      endif														 		!!
#ifdef FILEINFO
      open(99, file=simufile, position='append')  !=====!										!!
        write(99,*)					!										!!
	write(99,*) 'third loop: upgrading gpot'	!										!!
      close(99)  !======================================!										!!
      write(syscmd,'(a33,i3.3,a4,a14)') "echo -n 'third loop: i = 001 of '",nx," >> ",simufile						!!
      syslog=SYSTEM(syscmd)														!!
#else FILEINFO
      if(proc .eq. 0) then  !===================!											!!
        write(*,*)				!											!!
	write(*,*) 'third loop: upgrading gpot'	!											!!
      endif  !==================================!											!!
      write(syscmd,'(a33,i3.3)') "echo -n 'third loop: i = 001 of '",nx									!!
      syslog=SYSTEM(syscmd)														!!
#endif FILEINFO
      do i=1,nx  !======pierwsza poprawka do potencjalu - obliczanie gravaccel - petla po x-ie=====================!			!!
#ifdef FILEINFO
      write(syscmd2,'(a30,i3.3,a4,i3.3,a4,a14)') "echo -n '\b\b\b\b\b\b\b\b\b\b'",i," of ",nx," >> ",simufile	   !			!!
      syslog=SYSTEM(syscmd2)											   !			!!
#else FILEINFO
      write(syscmd2,'(a30,i3.3,a4,i3.3)') "echo -n '\b\b\b\b\b\b\b\b\b\b'",i," of ",nx				   !			!!
      syslog=SYSTEM(syscmd2)											   !			!!
#endif FILEINFO
        xi=x(i)													   !			!!
        do j=1,ny   !=========pierwsza poprawka do potencjalu - obliczanie gravaccel - petla po y-ku==============!!			!!
          yj=y(j)												  !!			!!
          rc = sqrt(xi**2+yj**2)										  !!			!!
          if(grav_model .eq. 'galsmooth') then !!								  !!			!!
            iOmega=vsun*tanh(rc/3.0/kpc)/rc	!								  !!			!!
          endif  !==============================!								  !!			!!
          do k=1,nz    !=========pierwsza poprawka do potencjalu - obliczanie gravaccel - petla po z-cie=========!!!			!!
            zk=z(k)												 !!!			!!
            if(grav_model .eq. 'galsmtdec') then  !=============================================================!!!!			!!
              iOmega=vsun*tanh(rc/3.0/kpc)/rc*exp(-zk**2/2.0/(h_grav*rotscalharray(i,j)/maxrotscalh)**n_gravh)	!!!!			!!
            endif  !============================================================================================!!!!			!!
            if(i .ne. 1 .and. i .ne. nx) then  !================================================================!!!!			!!
              gravaccel(1,i,j,k)=-iOmega**2*xi-0.5*c_si**2/gamma/u(1,i,j,k)*(u(1,i+1,j,k)-u(1,i-1,j,k))/dl(xdim)!!!!			!!
            endif  !============================================================================================!!!!			!!
            if(i .eq. 1) then  !================================================================================!!!!			!!
              gravaccel(1,i,j,k)=-iOmega**2*xi-0.5*c_si**2/gamma/u(1,i,j,k)*(u(1,i+1,j,k)-u(1,i,j,k))/dl(xdim)	!!!!			!!
            endif  !============================================================================================!!!!			!!
            if(i .eq. nx) then  !===============================================================================!!!!			!!
              gravaccel(1,i,j,k)=-iOmega**2*xi-0.5*c_si**2/gamma/u(1,i,j,k)*(u(1,i,j,k)-u(1,i-1,j,k))/dl(xdim)	!!!!			!!
            endif  !============================================================================================!!!!			!!
            if(j .ne. 1 .and. j .ne. ny) then  !================================================================!!!!			!!
              gravaccel(2,i,j,k)=-iOmega**2*yj-0.5*c_si**2/gamma/u(1,i,j,k)*(u(1,i,j+1,k)-u(1,i,j-1,k))/dl(ydim)!!!!			!!
            endif  !============================================================================================!!!!			!!
            if(j .eq. 1) then  !================================================================================!!!!			!!
              gravaccel(2,i,j,k)=-iOmega**2*yj-0.5*c_si**2/gamma/u(1,i,j,k)*(u(1,i,j+1,k)-u(1,i,j,k))/dl(ydim)	!!!!			!!
            endif  !============================================================================================!!!!			!!
            if(j .eq. ny) then  !===============================================================================!!!!			!!
              gravaccel(2,i,j,k)=-iOmega**2*yj-0.5*c_si**2/gamma/u(1,i,j,k)*(u(1,i,j,k)-u(1,i,j-1,k))/dl(ydim)	!!!!			!!
            endif  !============================================================================================!!!!			!!
          enddo  !========pierwsza poprawka do potencjalu - obliczanie gravaccel - petla po z-cie================!!!			!!
        enddo   !=======pierwsza poprawka do potencjalu - obliczanie gravaccel - petla po y-ku====================!!			!!
      enddo   !=======pierwsza poprawka do potencjalu - obliczanie gravaccel - petla po x-ie=======================!			!!
      call grav_accel2pot_req(gravaccel)								!				!!
      deallocate(gravaccel)										!				!!
!=======koniec dla przypadku 1==========================================================================!				!!
																	!!
      !===========poczatek warunku II-go stopnia=================================================================================!	!!
      if((przypadek .eq. 2)) then  !===second calculation of dprof and gpot=====================================================!!	!!
        allocate(dprof(nz))													!!	!!
	!==diagnostyczne========!												!!	!!
        itermaxwrite = 1	!												!!	!!
        maxdcolumn = 0.0	!												!!	!!
        mindcolumn = 1.0e30	!												!!	!!
        densdiscmean = 0.0	!												!!	!!
        adensdiscmean = 0.0	!												!!	!!
        ddensdiscmean = 0.0	!												!!	!!
        idensdiscmean = 0	!												!!	!!
	!=======================!												!!	!!
#ifdef FILEINFO
	open(99, file=simufile, position='append')  !===========!								!!	!!
	  write(99,*)						!								!!	!!
	  write(99,*) 'fourth loop: upgrading dprof, u(:,:,:)'	!								!!	!!
	close(99)  !============================================!								!!	!!
    write(syscmd3,'(a34,i3.3,a4,a14)') "echo -n 'fourth loop: j = 001 of '",ny," >> ",simufile					!!	!!
    syslog=SYSTEM(syscmd3)													!!	!!
#else FILEINFO
	if(proc .eq. 0) then  !=================================!								!!	!!
	  write(*,*)						!								!!	!!
	  write(*,*) 'fourth loop: upgrading dprof, u(:,:,:)'	!								!!	!!
	endif  !================================================!								!!	!!
    write(syscmd3,'(a34,i3.3)') "echo -n 'fourth loop: j = 001 of '",ny								!!	!!
    syslog=SYSTEM(syscmd3)													!!	!!
#endif FILEINFO
        do j = 1,ny   !======poprawka do tablicy u - petla po y-ku=========================================================!	!!	!!
#ifdef FILEINFO
      write(syscmd2,'(a30,i3.3,a4,i3.3,a4,a14)') "echo -n '\b\b\b\b\b\b\b\b\b\b'",j," of ",ny," >> ",simufile		   !	!!	!!
      syslog=SYSTEM(syscmd2)											      	   !	!!	!!
#else FILEINFO
      write(syscmd2,'(a30,i3.3,a4,i3.3)') "echo -n '\b\b\b\b\b\b\b\b\b\b'",j," of ",ny					   !	!!	!!
      syslog=SYSTEM(syscmd2)											      	   !	!!	!!
#endif FILEINFO
          yj = y(j)													   !	!!	!!
          do i = 1,nx   !=========poprawka do tablicy u - petla po x-ie===================================================!!	!!	!!
	    xi = x(i)													  !!	!!	!!
            rc = sqrt(xi**2+yj**2)											  !!	!!	!!
            !============================dprof(k) second calculation=============================================!	  !!	!!	!!
	    dcolumn = dcolumnarray(i,j)										 !	  !!	!!	!!
	    iter=0												 !	  !!	!!	!!
	    itermax=100												 !	  !!	!!	!!
	    dcolsmall=dcolumn*1.e-8										 !	  !!	!!	!!
	    d0 = dcolumn/dl(zdim)/nz											 !	  !!	!!	!!
	    call hydrostatic_zeq(i, j, d0, dprof)								 !	  !!	!!	!!
	    cd = 0.0												 !	  !!	!!	!!
	    do k=1,nz  !================!									 !	  !!	!!	!!
	      cd = cd +dprof(k)*nz	!									 !	  !!	!!	!!
	    enddo  !====================!									 !	  !!	!!	!!
	    do while((abs(cd - dcolumn) .gt. dcolsmall).and.(iter .le. itermax))  !==iteracyjne szukanie dprof==!!	  !!	!!	!!
              if(iter .eq. 0) then  !===========================!						!!	  !!	!!	!!
                d0previous = d0					!						!!	  !!	!!	!!
                cdprevious = cd					!						!!	  !!	!!	!!
                d0 = d0*2.					!						!!	  !!	!!	!!
              else						!						!!	  !!	!!	!!
                afactor = (cd - cdprevious)/(d0 - d0previous)	!						!!	  !!	!!	!!
                bfactor = cd - afactor*d0			!						!!	  !!	!!	!!
                d0previous = d0					!						!!	  !!	!!	!!
                cdprevious = cd					!						!!	  !!	!!	!!
                d0 = (dcolumn - bfactor)/afactor		!						!!	  !!	!!	!!
              endif  !==========================================!						!!	  !!	!!	!!
              iter = iter+1											!!	  !!	!!	!!
	      itermaxwrite = max(itermaxwrite, iter)								!!	  !!	!!	!!
	      call hydrostatic_zeq(i, j, d0, dprof)								!!	  !!	!!	!!
	      cd = 0.0												!!	  !!	!!	!!
	      do k=1,nz  !==============!									!!	  !!	!!	!!
	        cd = cd + dprof(k)*nz	!									!!	  !!	!!	!!
	      enddo  !==================!									!!	  !!	!!	!!
            enddo  !========iteracyjne szukanie dprof===========================================================!!	  !!	!!	!!
	    if(abs(cd-dcolumn).gt.dcolsmall) then  !=====================================!			 !	  !!	!!	!!
#ifdef FILEINFO
	      open(99, file=simufile, position='append') !==============================!!			 !	  !!	!!	!!
	        write(99,*) i,j,'equatorial density accuracy different than required!'	!!			 !	  !!	!!	!!
	      close(99)  !==============================================================!!			 !	  !!	!!	!!
#else FILEINFO
	        write(*,*) proc,i,j,'equatorial density accuracy different than required!'			 !	  !!	!!	!!
#endif FILEINFO
	    endif  !=====================================================================!			 !	  !!	!!	!!
            !====================================end of the second dprof(k) calculation==========================!	  !!	!!	!!
															  !!	!!	!!
            !====rotscalh calculation===========!									  !!	!!	!!
	    rotscalh = rotscalharray(i,j)	!									  !!	!!	!!
            !===end of rotscalh calculation=====!									  !!	!!	!!
            if((grav_model .eq. 'galsmooth') .or. (grav_model .eq. 'gallksun')) then !==!				  !!	!!	!!
              iOmega=vsun*tanh(rc/3.0/kpc)/rc						!				  !!	!!	!!
            endif  !====================================================================!				  !!	!!	!!
	    vx = -iOmega * yj  ! * rc/rc										  !!	!!	!!
	    vy =  iOmega * xi  ! * rc/rc										  !!	!!	!!
	    vz = 0.0													  !!	!!	!!
	    do k = 1,nz  !===============================================================================================!!!	!!	!!
	      zk=z(k)													 !!!	!!	!!
              if(grav_model .eq. 'galsmtdec') then  !===================================================================!!!!	!!	!!
                iOmega=vsun*tanh(rc/3.0/kpc)/rc*exp(-zk**2/2.0/(h_grav*rotscalharray(i,j)/maxrotscalh)**n_gravh)	!!!!	!!	!!
	        vx = -iOmega * yj  ! * rc/rc										!!!!	!!	!!
	        vy =  iOmega * xi  ! * rc/rc										!!!!	!!	!!
	        vz = 0.0												!!!!	!!	!!
              endif  !==================================================================================================!!!!	!!	!!
              u(1,i,j,k) = rhoa + dprof(k)/cosh((rc/r_max)**mtr)							 !!!	!!	!!
              u(1,i,j,k) = max(u(1,i,j,k), smalld)									 !!!	!!	!!
              u(2,i,j,k) = vx*u(1,i,j,k)										 !!!	!!	!!
              u(3,i,j,k) = vy*u(1,i,j,k)										 !!!	!!	!!
              u(4,i,j,k) = vz*u(1,i,j,k)	  									 !!!	!!	!!
#ifndef ISO
              u(5,i,j,k) = c_si**2/(gamma-1.0)*u(1,i,j,k)/gamma								 !!!	!!	!!
              u(5,i,j,k) = max(u(5,i,j,k), smallei)									 !!!	!!	!!
	      u(5,i,j,k) = u(5,i,j,k) +0.5*(vx**2+vy**2+vz**2)*u(1,i,j,k)						 !!!	!!	!!
#endif
         select case(mf_orient)
	 case('null')
         b(1,i,j,k)   = 0.0												 !!!	!!	!!
         b(2,i,j,k)   = 0.0												 !!!	!!	!!
         b(3,i,j,k)   = 0.0												 !!!	!!	!!
         case('vertical')
         b(1,i,j,k)   = 0.0												 !!!	!!	!!
         b(2,i,j,k)   = 0.0												 !!!	!!	!!
         b(3,i,j,k)   = sqrt(2.*alpha*d0*c_si**2)									 !!!	!!	!!
	 case('toroidal')
         b(1,i,j,k)   =-sqrt(2.*alpha*c_si**2*u(idna,i,j,k))*yj/rc							 !!!	!!	!!
         b(2,i,j,k)   = sqrt(2.*alpha*c_si**2*u(idna,i,j,k))*xi/rc							 !!!	!!	!!
         b(3,i,j,k)   = 0.0												 !!!	!!	!!
	 end select
#ifndef ISO
              u(5,i,j,k)   = u(5,i,j,k) +0.5*sum(b(:,i,j,k)**2,1)							 !!!	!!	!!
#endif
            enddo  !=====================================================================================================!!!	!!	!!
															  !!	!!	!!
	    !==diagnostyczne====================================================================================!	  !!	!!	!!
	    ddensdiscmean = maxval(u(1,i,j,:))									!	  !!	!!	!!
	    call MPI_ALLREDUCE(ddensdiscmean, adensdiscmean, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm, ierr)	!	  !!	!!	!!
	    if(adensdiscmean .gt. smalld) then  !===============!						!	  !!	!!	!!
	      densdiscmean = densdiscmean + adensdiscmean	!						!	  !!	!!	!!
	      idensdiscmean = idensdiscmean + 1			!						!	  !!	!!	!!
	    endif  !============================================!						!	  !!	!!	!!
            !===================================================================================================!	  !!	!!	!!
	  enddo  !========================================================================================================!!	!!	!!
        enddo  !===========================================================================================================!	!!	!!
        !==diagnostyczne=========================================================================!				!!	!!
	densdiscmean = densdiscmean/real(idensdiscmean)						 !				!!	!!
        if(proc .eq. 0) then  !=================================!				 !				!!	!!
#ifdef FILEINFO
          open(99, file=simufile, position='append')  !-----------------------------------------!!				!!	!!
            write(99,*) ' '					!				!!				!!	!!
	    write(99,*) 'Mean density of disk = ',densdiscmean	!				!!				!!	!!
#else FILEINFO
            write(*,*) ' '					!				!!				!!	!!
	    write(*,*) 'Mean density of disk = ',densdiscmean	!				!!				!!	!!
#endif FILEINFO
        endif  !================================================!				!!				!!	!!
        densmax=maxval(u(1,:,:,:))								!!				!!	!!
        densmin=minval(u(1,:,:,:))								!!				!!	!!
        call MPI_REDUCE(densmax, densmaxall, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierr)	!!				!!	!!
        call MPI_REDUCE(densmin, densminall, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, comm, ierr)	!!				!!	!!
        if(proc .eq. 0) then  !=========================================================!	!!				!!	!!
#ifdef FILEINFO
            write(99,*) 'maxdcolumn (radius = ',maxdcradius,' ) = ',maxdcolumn		!	!!				!!	!!
            write(99,*) 'mindcolumn (radius = ',mindcradius,' ) = ',mindcolumn		!	!!				!!	!!
            write(99,*) 'maxdens = ', densmaxall					!	!!				!!	!!
            write(99,*) 'mindens = ', densminall					!	!!				!!	!!
          close(99)  !--------------------------------------------------------------------------!!				!!	!!
        endif  !=======================================================================!	 !				!!	!!
        open(99, file=simufile, position='append')  !=====================================!	 !				!!	!!
          write(99,*) 'the longest second iteration number of steps = ', itermaxwrite	  !	 !				!!	!!
        close(99)  !======================================================================!	 !				!!	!!
#else FILEINFO
            write(*,*) 'maxdcolumn (radius = ',maxdcradius,' ) = ',maxdcolumn		!	 !				!!	!!
            write(*,*) 'mindcolumn (radius = ',mindcradius,' ) = ',mindcolumn		!	 !				!!	!!
            write(*,*) 'maxdens = ', densmaxall						!	 !				!!	!!
            write(*,*) 'mindens = ', densminall						!	 !				!!	!!
          write(*,*) 'the longest second iteration number of steps = ', itermaxwrite	!	 !				!!	!!
        endif  !========================================================================!	 !				!!	!!
#endif FILEINFO
	!========================================================================================!			  	!!	!!
   															  	!!	!!
        deallocate(dprof)												   	!!	!!
        allocate(gravaccel(2,nx,ny,nz))											   	!!	!!
#ifdef FILEINFO
        open(99, file=simufile, position='append')  !===========!							   	!!	!!
	  write(99,*)						!							   	!!	!!
	  write(99,*) 'fifth loop: second upgrading gpot'	!							   	!!	!!
        close(99)  !============================================!							    	!!	!!
    write(syscmd,'(a33,i3.3,a4,a14)') "echo -n 'fifth loop: i = 001 of '",nx," >> ",simufile					!!	!!
    syslog=SYSTEM(syscmd)													!!	!!
#else FILEINFO
	  write(*,*)													   	!!	!!
	  write(*,*) 'fifth loop: second upgrading gpot'								   	!!	!!
    write(syscmd,'(a33,i3.3)') "echo -n 'fifth loop: i = 001 of '",nx								!!	!!
    syslog=SYSTEM(syscmd)													!!	!!
#endif FILEINFO
        do i=1,nx  !=======================================================================================================!	!!	!!
#ifdef FILEINFO
      write(syscmd2,'(a30,i3.3,a4,i3.3,a4,a14)') "echo -n '\b\b\b\b\b\b\b\b\b\b'",i," of ",nx," >> ",simufile		   !	!!	!!
      syslog=SYSTEM(syscmd2)											      	   !	!!	!!
#else FILEINFO
      write(syscmd2,'(a30,i3.3,a4,i3.3)') "echo -n '\b\b\b\b\b\b\b\b\b\b'",i," of ",nx					   !	!!	!!
      syslog=SYSTEM(syscmd2)											      	   !	!!	!!
#endif FILEINFO
          xi=x(i)													   !	!!	!!
          do j=1,ny  !====================================================================================================!!	!!	!!
            yj=y(j)													  !!	!!	!!
            rc = sqrt(xi**2+yj**2)											  !!	!!	!!
            if(grav_model .eq. 'galsmooth') then  !=====!								  !!	!!	!!
              iOmega=vsun*tanh(rc/3.0/kpc)/rc		!								  !!	!!	!!
            endif  !====================================!								  !!	!!	!!
            do k=1,nz  !=================================================================================================!!!	!!	!!
              zk=z(k)													 !!!	!!	!!
              if(grav_model .eq. 'galsmtdec') then  !===================================================================!!!!	!!	!!
                iOmega=vsun*tanh(rc/3.0/kpc)/rc*exp(-zk**2/2.0/(h_grav*rotscalharray(i,j)/maxrotscalh)**n_gravh)	!!!!	!!	!!
              endif  !==================================================================================================!!!!	!!	!!
															 !!!	!!	!!
              if(i .ne. 1 .and. i .ne. nx) then  !================================================================!	 !!!	!!	!!
                gravaccel(1,i,j,k)=-iOmega**2*xi-0.5*c_si**2/gamma/u(1,i,j,k)*(u(1,i+1,j,k)-u(1,i-1,j,k))/dl(xdim)!	 !!!	!!	!!
              endif  !============================================================================================!	 !!!	!!	!!
              if(i .eq. 1) then  !==============================================================================!	 !!!	!!	!!
                gravaccel(1,i,j,k)=-iOmega**2*xi-0.5*c_si**2/gamma/u(1,i,j,k)*(u(1,i+1,j,k)-u(1,i,j,k))/dl(xdim)!	 !!!	!!	!!
              endif  !==========================================================================================!	 !!!	!!	!!
              if(i .eq. nx) then  !=============================================================================!	 !!!	!!	!!
                gravaccel(1,i,j,k)=-iOmega**2*xi-0.5*c_si**2/gamma/u(1,i,j,k)*(u(1,i,j,k)-u(1,i-1,j,k))/dl(xdim)!	 !!!	!!	!!
              endif  !==========================================================================================!	 !!!	!!	!!
              if(j .ne. 1 .and. j .ne. ny) then  !================================================================!	 !!!	!!	!!
                gravaccel(2,i,j,k)=-iOmega**2*yj-0.5*c_si**2/gamma/u(1,i,j,k)*(u(1,i,j+1,k)-u(1,i,j-1,k))/dl(ydim)!	 !!!	!!	!!
              endif  !============================================================================================!	 !!!	!!	!!
              if(j .eq. 1) then  !==============================================================================!	 !!!	!!	!!
                gravaccel(2,i,j,k)=-iOmega**2*yj-0.5*c_si**2/gamma/u(1,i,j,k)*(u(1,i,j+1,k)-u(1,i,j,k))/dl(ydim)!	 !!!	!!	!!
              endif  !==========================================================================================!	 !!!	!!	!!
              if(j .eq. ny) then  !=============================================================================!	 !!!	!!	!!
                gravaccel(2,i,j,k)=-iOmega**2*yj-0.5*c_si**2/gamma/u(1,i,j,k)*(u(1,i,j,k)-u(1,i,j-1,k))/dl(ydim)!	 !!!	!!	!!
              endif  !==========================================================================================!	 !!!	!!	!!
            enddo  !=====================================================================================================!!!	!!	!!
          enddo  !========================================================================================================!!	!!	!!
        enddo  !===========================================================================================================!	!!	!!
        call grav_accel2pot_req(gravaccel)											!!	!!
        deallocate(gravaccel)													!!	!!
      endif  !=========================end of second calculation of dprof and gpot==============================================!!	!!
      !==============koniec warunku II-go stopnia================================================================================!	!!
      deallocate(rotscalharray)														!!
    endif !==(rotsc .ne. 'auto')========================================================================================================!!
!===================koniec warunku I-go stopnia==========================================================================================!


#ifdef SNE_DISTR
   call compute_SNdistr
#endif SNE_DISTR


#ifdef FILEINFO
    open(99, file=simufile, position='append')  !=======!
      write(99,*) 'end of init_prob'			!
    close(99)  !========================================!
#else FILEINFO
      write(*,*) 'end of init_prob'			
#endif FILEINFO
    
    return
  end subroutine init_prob  

!-----------------------------------------------------------------------------------------------------------------------------------

  subroutine grav_accel2pot_req(gravaccel)
  
    implicit none
    integer i, j, k, ip, pgpmax
    real gravaccel(2,nx,ny,nz)
    real gravrx(nx), gravry(ny), gravrz(nz)
    real gp_max, dgp(3), dgp0
    integer, dimension(3) :: loc_gp_max, gploc
    integer               :: proc_gp_max
    integer               :: px, py, pz, pc(3)
    real dgpx_proc, dgpx_all(0:nproc-1), &
         dgpy_proc, dgpy_all(0:nproc-1), &
	 dgpz_proc, dgpz_all(0:nproc-1), &
         dgpx(0:pxsize-1,0:pysize-1,0:pzsize-1), &
         dgpy(0:pxsize-1,0:pysize-1,0:pzsize-1), &
	 dgpz(0:pxsize-1,0:pysize-1,0:pzsize-1), &
	 ddgp(0:pxsize-1,0:pysize-1,0:pzsize-1)
    
    
! Gravitational potential gp(i,j,k) is defined in cell centers  
! First we find gp(:,:,:) in each block separately assuming:  
    
    gp(1,1,1) = 0.0

!    call grav_accel('xsweep',1,1, xr, nx, gravrx)    
    gravrx=gravaccel(1,:,1,1)
    do i = 1, nx-1      
      gp(i+1,1,1) = gp(i,1,1) - gravrx(i)*dl(xdim)
    enddo
    
    do i=1,nx
!      call grav_accel('ysweep',1, i, yr, ny, gravry)
      gravry=gravaccel(2,i,:,1)
      do j = 1, ny-1      
        gp(i,j+1,1) = gp(i,j,1) - gravry(j)*dl(ydim)
      enddo
    enddo

    do i=1,nx
      do j=1,ny
        call grav_accel('zsweep', i, j, zr, nz, gravrz, 'default')
        do k = 1, nz-1      
          gp(i,j,k+1) = gp(i,j,k) - gravrz(k)*dl(zdim)
        enddo
      enddo
    enddo
            
    dgpx_proc = gp(nxb+1,1,1)-gp(1,1,1)
    dgpy_proc = gp(1,nyb+1,1)-gp(1,1,1)
    dgpz_proc = gp(1,1,nzb+1)-gp(1,1,1)

    call MPI_GATHER ( dgpx_proc, 1, MPI_DOUBLE_PRECISION, &
                      dgpx_all,  1, MPI_DOUBLE_PRECISION, &
                      0, comm3d,err )
    call MPI_GATHER ( dgpy_proc, 1, MPI_DOUBLE_PRECISION, &
                      dgpy_all,  1, MPI_DOUBLE_PRECISION, &
                      0, comm3d,err )
    call MPI_GATHER ( dgpz_proc, 1, MPI_DOUBLE_PRECISION, &
                      dgpz_all,  1, MPI_DOUBLE_PRECISION, &
                      0, comm3d,err )
		      

    if(proc .eq. 0) then
    
      do ip = 0, nproc-1
        call MPI_CART_COORDS(comm3d, ip, ndims, pc, ierr)
       
        px = pc(1)
        py = pc(2)
        pz = pc(3)
      
        dgpx(px,py,pz) = dgpx_all(ip)
        dgpy(px,py,pz) = dgpy_all(ip)
        dgpz(px,py,pz) = dgpz_all(ip)

      enddo
    
      ddgp(0,0,0) = 0.0
      do i = 1, pxsize-1      
        ddgp(i,0,0) = ddgp(i-1,0,0) + dgpx(i-1,0,0)
      enddo
   
      do i=0, pxsize-1
        do j = 1, pysize-1      
          ddgp(i,j,0) = ddgp(i,j-1,0) + dgpy(i,j-1,0)
        enddo
      enddo

      do i=0, pxsize-1
        do j=0, pysize-1
          do k = 1, pzsize-1      
             ddgp(i,j,k)= ddgp(i,j,k-1) + dgpz(i,j,k-1)
          enddo
        enddo
      enddo
    
    endif

    call MPI_BCAST(ddgp, nproc, MPI_DOUBLE_PRECISION, 0, comm, ierr)       

    px = pcoords(1)
    py = pcoords(2)
    pz = pcoords(3)
    
    gp = gp + ddgp(px,py,pz)
      
    gp_max      = maxval(gp(nb+1:nb+nxb, nb+1:nb+nyb,nb+1:nb+nzb)) 
    loc_gp_max  = maxloc(gp(nb+1:nb+nxb, nb+1:nb+nyb,nb+1:nb+nzb)) &
                  + (/nb,nb,nb/)
		  
    call mpifind(gp_max, 'max', loc_gp_max, proc_gp_max)
    pgpmax = proc_gp_max
    
    call MPI_BCAST(gp_max, 1, MPI_DOUBLE_PRECISION, pgpmax, comm, ierr)       
    gp = gp - gp_max
   
  end subroutine grav_accel2pot_req
  

end module init_problem

