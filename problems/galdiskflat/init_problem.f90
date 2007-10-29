module init_problem
  
! Initial condition for Keplerian disk
! Written by: M. Hanasz, March 2006

  use arrays
  use start
  use grid
  use hydrostatic
  use gravity

  real d0, r_max
  integer mtr
  character problem_name*32,run_id*3

  namelist /PROBLEM_CONTROL/  problem_name, run_id, &
                              d0,r_max,mtr, nzfac    


contains

!-----------------------------------------------------------------------------

  subroutine read_problem_par

    implicit none
  

    problem_name = 'aaa'
    run_id  = 'aa'
    d0      = 1.0
    r_max   = 0.8
    mtr      = 10
         
    
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
      rbuff(2) = r_max

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

      d0           = rbuff(1)  
      r_max        = rbuff(2)  
    
      mtr          = ibuff(1)
    
    endif

  end subroutine read_problem_par

!-----------------------------------------------------------------------------

  subroutine init_prob

    implicit none
 
    integer i,j,k
    real xi,yj,zk, rc, vx, vy, vz,h2,dgdz,csim2
    real, allocatable ::dprof(:)
    real ivsun, iaconst, ibconst, icconst, iflat, iOmega
    real cm, mprot, alfar, dens0, skalamasy
    real dnmol, dncold, dnwarm, dnion, dnhot
    real densmax, densmin, densmaxall, densminall, maxdcolumn, mindcolumn, maxdcradius, mindcradius
    real dcolumn, dcolsmall, cd, cdprevious, d0previous, afactor, bfactor, dcolumnprevious
    integer iter, itermx, itermxwrite, inzfac
    
    call read_problem_par

    allocate(dprof(nz))

! ivsun=220.0*3600*24.0*365.2562e11/3.0856e18
! iaconst=(-1.0e-7)*ivsun
! ibconst=-3000.0*iaconst
! icconst=ivsun/5000.0-iaconst*5000.0-ibconst*log(5000.0)
! iflat=iaconst*3000.0+ibconst*log(3000.0)+icconst

 ivsun=220.0/3.0856e13 ![pc/s]
 iaconst=-ivsun/(5000.0*(5000.0-3000.0)) ![pc/s/pc2]
 ibconst=-3000.0*iaconst ![1/s]
 icconst=ivsun/5000.0-iaconst*5000.0-ibconst*log(5000.0) ![1/s]
 iflat=iaconst*3000.0+ibconst*log(3000.0)+icconst ![1/s]

itermxwrite = 1

    do j = 1,ny
      yj = y(j)
      do i = 1,nx
	xi = x(i)
        rc = sqrt(xi**2+yj**2)
            

      if(rc.le.3.0*kpc) then
      iOmega=iflat*3600.0*24.0*365.2562e6
      else
      if(rc.ge.5.0*kpc) then
      iOmega=ivsun/rc*3600.0*24.0*365.2562e6
      else
      iOmega=(iaconst*rc+ibconst*log(rc)+icconst)*3600.0*24.0*365.2562e6
      endif
      endif
      
!      write(*,*) 'nb = ',nb
!      write(*,*) 'dz0 = ',dz0
!      write(*,*) 'zmin = ',zmin
      
      maxdcolumn = 0.0
      mindcolumn = 1.0e30
      dcolumn = 1.0e-6
      dcolumnprevious = 0.0
      inzfac = 0
      do while((dcolumn-dcolumnprevious).ge.(1.0e-6*dcolumn))
      inzfac = inzfac + 1
      dcolumnprevious = dcolumn
      dcolumn = 0.0
	do k = 1,2*inzfac*nz
          zk = inzfac*zmin + 0.5*dz0 + k*dz0	  
	    
	    cm = 1./3.0856e18
	    skalamasy=1./1.989e33
	    mprot = 1.672614e-24 * skalamasy
	    dnmol = 0.58/(cm**3) * exp(-((rc - 4.5*kpc)**2-(r_gc_sun*kpc - 4.5*kpc)**2)/(2.9*kpc)**2) &
	    		& * (rc/r_gc_sun/kpc)**(-0.58) * exp(-(zk/(81.0*(rc/r_gc_sun/kpc)**(0.58)))**2)
	    
	    if(rc.lt.r_gc_sun*kpc) then
	    alfar=1.0
	    else
	    alfar=rc/r_gc_sun/kpc
	    endif
	    dncold= 0.34/(cm**3)/alfar**2 * (0.859*exp(-(zk/(127.0*alfar))**2) + 0.047*exp(-(zk/(318.0*alfar))**2) &
	    		& + 0.094*exp(-(zk/(403.0*alfar))**2))
	    dnwarm= 0.226/(cm**3)/alfar * ((1.745 - 1.289/alfar)*exp(-(zk/(127.0*alfar))**2) &
	    		& + (0.473 - 0.07/alfar)*exp(-(zk/(318.0*alfar))**2) + (0.283 - 0.142/alfar)*exp(-(zk/(403.0*alfar))**2))
	    dnion = 0.0237/(cm**3)*exp(-(rc**2 - (r_gc_sun*kpc)**2)/(37.0*kpc)**2)*exp(-abs(zk)/kpc) &
	    		& + 0.0013/(cm**3)*exp(-((rc - 4.0*kpc)**2 - (r_gc_sun*kpc - 4.0*kpc)**2)/(2.0*kpc)**2)*exp(-abs(zk)/150.0)
	    dnhot = 4.8e-4/(cm**3)*(0.12*exp(-(rc - r_gc_sun*kpc)/(4.9*kpc)) + 0.88*exp(-((rc - 4.5*kpc)**2 &
	    		& - (r_gc_sun*kpc - 4.5*kpc)**2)/(2.9*kpc)**2)) * (rc/r_gc_sun/kpc)**(-1.65)*exp(-abs(zk) &
			& /(1.5*kpc*(rc/r_gc_sun/kpc)**(1.65)))
	    dens0=1.36*mprot*(dnmol+dncold+dnwarm+dnion+dnhot)
	    dcolumn=dcolumn+dens0*dz0
	 enddo !k=1,2*iznfac*nz   
      enddo !inzfac accuracy
!      write(*,*) 'accuracy: inzfac = ',inzfac    
	    if(mindcolumn .gt. dcolumn) mindcradius = rc
	    mindcolumn=min(mindcolumn,dcolumn)
	    if(maxdcolumn .lt. dcolumn) maxdcradius = rc
	    maxdcolumn=max(maxdcolumn,dcolumn)
	
	iter=0
	itermx=100
	dcolsmall=dcolumn*1.e-8
	d0 = dcolumn/dz0/nz
	call hydrostatic_zeq(i, j, d0, dprof)
	cd = 0.0
	do k=1,nz
	cd = cd +dprof(k)*nz
	enddo
	do while((abs(cd - dcolumn) .gt. dcolsmall).and.(iter .le. itermx))
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
	itermxwrite = max(itermxwrite, iter)
	call hydrostatic_zeq(i, j, d0, dprof)
	cd = 0.0
	do k=1,nz
	cd = cd + dprof(k)*nz
	enddo
        enddo
	if(abs(cd-dcolumn).gt.dcolsmall) write(*,*) i,j,'equatorial density accuracy different than required!'
	

	 do k = 1,nz
	   zk=z(k)


	    vx = -iOmega * yj  ! * rc/rc
	    vy =  iOmega * xi  ! * rc/rc
	    vz = 0.0

          u(1,i,j,k) = maxval(dprof(:))/cosh((rc/r_max)**mtr)
!	  u(1,i,j,k) = dens0/cosh((rc/r_max)**mtr)
          u(1,i,j,k) = max(u(1,i,j,k), smalld)
	  	  	  
          u(2,i,j,k) = vx*u(1,i,j,k)
          u(3,i,j,k) = vy*u(1,i,j,k)
          u(4,i,j,k) = vz*u(1,i,j,k)	  
          u(5,i,j,k) = c_si**2/(gamma-1.0)*u(1,i,j,k)
          u(5,i,j,k) = max(u(5,i,j,k), smallei)
	  
	  u(5,i,j,k) = u(5,i,j,k) +0.5*(vx**2+vy**2+vz**2)*u(1,i,j,k)

          b(1,i,j,k)   = 0.0
          b(2,i,j,k)   = 0.0
          b(3,i,j,k)   = 0.0
          u(5,i,j,k)   = u(5,i,j,k) +0.5*sum(b(:,i,j,k)**2,1)
        enddo
      enddo
    enddo
    densmax=maxval(u(1,:,:,:))
    densmin=minval(u(1,:,:,:))
    call MPI_REDUCE(densmax, densmaxall, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierr)
    call MPI_REDUCE(densmin, densminall, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, comm, ierr)
    if(proc .eq. 0) then
    write(*,*) 'maxdcolumn (radius = ',maxdcradius,' ) = ',maxdcolumn
    write(*,*) 'mindcolumn (radius = ',mindcradius,' ) = ',mindcolumn
    write(*,*) 'maxdens = ', densmaxall
    write(*,*) 'mindens = ', densminall
    endif
    write(*,*) 'itermxwrite (prof = ',proc,') = ', itermxwrite
    
    deallocate(dprof)
    
    return
  end subroutine init_prob  
  

end module init_problem

