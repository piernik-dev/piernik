#include "mhd.def"

module init_problem
  
! Initial condition for Galactic disk
! Written by: D. Wolt, June 2007
! Restrictions: using only cor bnd_cond, disk radius must be smaller than xmax and ymax
  use mpi_setup

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
  
    
    character par_file*(100), tmp_log_file*(100)
    integer :: cwd_status 

    par_file = trim(cwd)//'/problem.par'
    tmp_log_file = trim(cwd)//'/tmp.log'    


    problem_name = 'aaa'
    run_id  = 'aa'
    rhoa    = 1.0e-4
    d0      = 1.0
    r_max   = 0.8
    mtr      = 10
    mf_orient = 'null' ! 'toroidal', 'vertical'
         
    
    if(proc .eq. 0) then
      open(1,file=par_file)
        read(unit=1,nml=PROBLEM_CONTROL)
      close(1)
        write(*,nml=PROBLEM_CONTROL)
      open(3, file=tmp_log_file, position='append')
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
    use arrays, only    :   u,b,x,y,z,nx,ny,nz,nxt,nxb,dl,gp,xdim,ydim,idna
    use constants, only :   cm,kpc,mp,r_gc_sun
    use grid, only      :   dx
    use hydrostatic, only : hydrostatic_zeq
    use start, only     :   xmin,nb,smalld,c_si,alpha,gamma
#ifndef ISO
    use start, only    :   smallei
#endif /* ISO */

! wolanie prepare_snedistr przeniesione do mhd - procedura musi byc wolana
! rowniez po restarcie 

#ifdef COSM_RAYS
    use arrays, only    :  iecr
    use start, only     : gamma_cr, beta_cr
#endif /* COSM_RAYS */
    implicit none
 
    integer i,j,k,iu,id,ju,jd
    real xi,yj,zk, rc, vx, vy, rs
    real, allocatable ::dprof(:)
    real iOmega, dens0
    real dcmol, dcneut, dcion, dchot
    real xgradgp, ygradgp, sfq
#ifdef PRESSURECORRECTION
    real xgradp, ygradp
#endif /* PRESSURECORRECTION */
    character syscmd*37,syscmd2*40,syscmd3*37
    integer system, syslog
    real,allocatable :: dxzprof(:,:),idxzprof(:,:),jdxzprof(:,:)
    integer,allocatable :: xproc(:), yproc(:)
    integer xtag, iproc, ilook, jproc, ytag

    allocate(dprof(nz))
    allocate(dxzprof(nxt,nz))
    allocate(idxzprof(nx,nz))
    allocate(jdxzprof(nx,nz))
    if(pcoords(2) .eq. 0) then
    do i = 1,nx
      rc = x(i)
      
      dcmol = 2.6*0.8*exp(-((rc - 4.5*kpc)**2-(r_gc_sun - 4.5*kpc)**2)/(2.9*kpc)**2)
      dcneut = 6.2*0.8
      dcion =  1.46*0.8*exp(-(rc**2 - r_gc_sun**2)/(37.0*kpc)**2) &
             + 1.20e-2*0.8*exp(-((rc - 4.0*kpc)**2 - (r_gc_sun - 4.0*kpc)**2)/(2.0*kpc)**2)
      dchot = 4.4e-2*0.8*(0.12*exp(-(rc-r_gc_sun)/4.9/kpc)+0.88*exp(-((rc-4.5*kpc)**2-(r_gc_sun-4.5*kpc)**2)/(2.9*kpc)**2))
      d0=(dcmol+dcneut+dcion+dchot) * 1.36 
      call hydrostatic_zeq(i, 0, d0, dprof)

      dxzprof(pcoords(1)*nxb+i,:) = dprof
      idxzprof(i,:) = dprof

    enddo
    
!    write(*,*) dprof
!    stop
    
    endif
    if(psize(2) .ne. 1) then
      allocate(yproc(psize(2)))
      do j = 1, psize(2)
        coords = (/pcoords(1),j-1,pcoords(3)/)
	call MPI_Cart_rank(comm3d, coords, jproc, ierr)
	yproc(j)=jproc
      enddo
      if(proc .eq. yproc(1)) then
        do j = 2,psize(2)
	  ytag=1000+100*pcoords(3)+10*j+pcoords(1)
	  call MPI_Send(idxzprof,nx*nz,MPI_DOUBLE_PRECISION,yproc(j),ytag,comm,ierr)
	enddo
      else
        ytag=1000+100*pcoords(3)+10*(pcoords(2)+1)+pcoords(1)
	call MPI_Recv(idxzprof,nx*nz,MPI_DOUBLE_PRECISION,yproc(1),ytag,comm,status,ierr)
	dxzprof(pcoords(1)*nxb+1:pcoords(1)*nxb+nx,:)=idxzprof
      endif
      deallocate(yproc)
    endif
    if(psize(1) .ne. 1) then
      allocate(xproc(psize(1)))
      do i = 1, psize(1)
        coords = (/i-1,pcoords(2),pcoords(3)/)
	call MPI_Cart_rank(comm3d, coords, iproc, ierr)
	xproc(i)=iproc
      enddo
      do i = 1, psize(1)
        xtag = 100*pcoords(3)+10*pcoords(2)+i
        if(proc .eq. xproc(i)) then
	  do j = 1, psize(1)
            if(proc .ne. xproc(j)) then
	      call MPI_Send(idxzprof,nx*nz,MPI_DOUBLE_PRECISION,xproc(j),xtag,comm,ierr)
	    endif
          enddo
	else
	  call MPI_Recv(jdxzprof,nx*nz, MPI_DOUBLE_PRECISION,xproc(i),xtag,comm,status,ierr)
	  dxzprof((i-1)*nxb+1:(i-1)*nxb+nx,:)=jdxzprof(:,:)
	endif
      enddo
      deallocate(xproc)
    endif
    deallocate(idxzprof,jdxzprof)

    if(proc .eq. 0) then
      write(syscmd,'(a34,i3.3)') "echo -n 'second loop: j = 001 of '",ny
      syslog=SYSTEM(syscmd)
    endif
    do j = 1,ny
!      write(syscmd2,'(a30,i3.3,a4,i3.3)') "echo -n '\b\b\b\b\b\b\b\b\b\b'",j," of ",ny
      if(proc .eq. 0) syslog=SYSTEM(syscmd2)
      yj = y(j)
      do i = 1,nx
        xi = x(i)
        rc = sqrt(xi**2+yj**2)

        ilook = (rc-xmin)/dx + 0.5 + nb
	if(int(ilook) .lt. nxt) then
	  dprof(:) = dxzprof(int(ilook),:)+(dxzprof(int(ilook)+1,:)-dxzprof(int(ilook),:))*(ilook-int(ilook))
	else
	  dprof(:) = dxzprof(nxt,:)
	endif
	
	do k=1,nz
	 zk=z(k)
	 rs = sqrt(xi**2+yj**2+zk**2)
         u(1,i,j,k) = rhoa + dprof(k)/cosh(min((rc/r_max)**mtr,100.0))
         u(1,i,j,k) = max(u(1,i,j,k), smalld)
	   if(i .ne. 1 .and. i .ne. nx) then
	     iu = i+1
	     id = i-1
	     sfq = 0.5
	   else
	     if(i .eq. 1) then
	       iu = i+1
	       id = i
	       sfq = 1.0
	     elseif(i .eq. nx) then
	       iu = i
	       id = i-1
	       sfq = 1.0
	     endif
           endif
	   xgradgp=(gp(iu,j,k)-gp(id,j,k))*sfq/dl(xdim)
#ifdef PRESSURECORRECTION
           xgradp =-sfq*c_si**2/gamma/u(1,i,j,k)*(u(1,iu,j,k)-u(1,id,j,k))/dl(xdim)
#endif /* PRESSURECORRECTION */
	   if(j .ne. 1 .and. j .ne. ny) then
	     ju = j+1
	     jd = j-1
	     sfq = 0.5
	   else
	     if(j .eq. 1) then
	       ju = j+1
	       jd = j
	       sfq = 1.0
	     elseif(j .eq. ny) then
	       ju = j
	       jd = j-1
	       sfq = 1.0
	     endif
           endif
	   ygradgp=(gp(i,ju,k)-gp(i,jd,k))*sfq/dl(ydim)
#ifdef PRESSURECORRECTION
           ygradp =-sfq*c_si**2/gamma/u(1,i,j,k)*(u(1,i,ju,k)-u(1,i,jd,k))/dl(ydim)
           iOmega=sqrt(abs(sqrt((xgradgp+xgradp)**2+(ygradgp+ygradp)**2))/rc)
#else /* PRESSURECORRECTION */
           iOmega=sqrt(abs(sqrt(xgradgp**2+ygradgp**2))/rc)
#endif /* PRESSURECORRECTION */
	     u(2,i,j,k)=-iOmega*yj*u(1,i,j,k)
             u(3,i,j,k)= iOmega*xi*u(1,i,j,k)
             u(4,i,j,k) = 0.0
#ifndef ISO
         u(5,i,j,k) = c_si**2/(gamma-1.0)*u(1,i,j,k)
         u(5,i,j,k) = max(u(5,i,j,k), smallei)
	 u(5,i,j,k) = u(5,i,j,k) +0.5*(u(2,i,j,k)**2+u(3,i,j,k)**2+u(4,i,j,k)**2)/u(1,i,j,k)
#endif /* ISO */

#ifdef COSM_RAYS
          u(iecr,i,j,k)   =  beta_cr*c_si**2 * u(idna,i,j,k)/(gamma_cr-1.0)
#endif /* COSM_RAYS */


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
#endif /* ISO */
        enddo
      enddo
    enddo
    write(*,*) ' '
    if(allocated(dprof)) deallocate(dprof)


    return
  end subroutine init_prob  

!-----------------------------------------------------------------------------------------------------------------------------------

end module init_problem

