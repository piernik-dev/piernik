#include "piernik.def"


module init_problem
  
! Initial condition for Keplerian disk
! Written by: M. Hanasz, March 2006

  use mpi_setup

  real d0, collf, r_max, alfasupp, rsup1, rsup2, rsup3, rsup4, dout
  character problem_name*32,run_id*3, mag_field_orient*32, suppkind*6

  namelist /PROBLEM_CONTROL/  problem_name, run_id, &
                              d0,dout,r_max,collf,alfasupp, &
			      rsup1,rsup2,rsup3,rsup4,suppkind,mag_field_orient    


contains

!-----------------------------------------------------------------------------

  subroutine read_problem_par

    implicit none

    problem_name = 'aaa'
    run_id  = 'aa'
    d0      = 1.0
    dout    = 1.0e-4
    r_max   = 1.0
    collf   = 0.0
    alfasupp= 0.0
    rsup1   = 0.0
    rsup2   = 1.0
    rsup3   = 2.0
    rsup4   = 3.0
    suppkind= 'full' !'smooth' !'sharp' !'smout' !'smpeak'
    mag_field_orient = 'none'
    
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
      cbuff(3) =  suppkind
      cbuff(4) =  mag_field_orient

      rbuff(1) = d0
      rbuff(2) = dout
      rbuff(3) = r_max
      rbuff(4) = collf
      rbuff(5) = alfasupp
      rbuff(6) = rsup1
      rbuff(7) = rsup2
      rbuff(8) = rsup3
      rbuff(9) = rsup4
    
      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

    else
    
      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)
      
      problem_name = cbuff(1)   
      run_id       = cbuff(2)   
      suppkind     = cbuff(3)
      mag_field_orient = cbuff(4)

      d0           = rbuff(1)
      dout         = rbuff(2)
      r_max	   = rbuff(3)
      collf	   = rbuff(4)
      alfasupp     = rbuff(5)
      rsup1        = rbuff(6)
      rsup2        = rbuff(7)
      rsup3        = rbuff(8)
      rsup4        = rbuff(9)
    
    endif

  end subroutine read_problem_par

!-----------------------------------------------------------------------------

  subroutine init_prob
    use arrays,      only : u,b,x,y,z,nx,ny,nz,nfluid,idna,imxa,imya,imza
    use arrays,      only : bndxrar,bndyrar,nu,fadiab,iena,fmagn,gp
    use constants,   only : newtong
    use grid,        only : dx,dy
    use hydrostatic, only : hydrostatic_zeq
    use start,       only : alpha,c_si,dimensions,r_smooth,r_grav,n_gravr,ptmass,smalld,collfaq
    use start,       only : smallei,gamma,nb,xmin
#ifdef KEPLER_SUPPRESSION
    use arrays,      only : alfsup,omx0,omy0
#endif /* KEPLER_SUPPRESSION */
    implicit none
 
    integer i,j,k,kmid,iu,id,ju,jd, ilook
    real xi,yj,zk, rc, rs, vx, vy, vz, h2, dgdz, csim2, b0, sqr_gm, v_phi
    real, allocatable ::dprof(:), gpot(:)
    real sfq, xgradgp, xgradp, ygradgp, ygradp, iOmega, gradp, gradgp, d00
    real, dimension(nx) :: omega
    
!    call read_problem_par

!   Secondary parameters

    sqr_gm = sqrt(newtong*ptmass)

    b0 = sqrt(2.*alpha*d0*c_si**2) 
    write(*,*) 'b0 = ',b0 
    allocate(dprof(nz),gpot(nz))
    
    do k=1, nz
      if(z(k) .lt. 0.0) kmid = k       ! the midplane is in between 
    enddo                                  ! ksmid and ksmid+1

    collfaq = collf/d0

#ifdef KEPLER_SUPPRESSION
    allocate(alfsup(nx,ny))
    allocate(omx0(nfluid,nx,ny,nz))
    allocate(omy0(nfluid,nx,ny,nz))
#endif /* KEPLER_SUPPRESSION */
      
    do j = 1,ny
      yj = y(j)
      do i = 1,nx
        xi = x(i)
        rc = sqrt(xi**2+yj**2)
        
        if(dimensions .eq.'3d') then
	  d00 = 1.0
          call hydrostatic_zeq(i, j, d00, dprof)    
        endif
                          
        do k = 1,nz
          if(dimensions .eq.'3d') then
            u(idna,i,j,k) = min((rc/r_grav)**n_gravr,100.0)
            u(idna,i,j,k) = dout + (dprof(k)-dout)/cosh(u(idna,i,j,k))
            u(idna,i,j,k) = 9.5*(rc/rsup1)**2*exp(-3.5*(rc/rsup1)**2)*d0*dprof(k) !/cosh((rc/r_max)**25.0)
            u(idna(1),i,j,k) = d0/cosh(min((rc/rsup2)**25.0,690.0))*dprof(k)
            u(idna(2),i,j,k) = d0/100.0/cosh(min((rc/rsup2)**25.0,690.0))*dprof(k)
          else
!            u(idna,i,j,k) = min((rc/r_grav)**n_gravr,100.0)
!            u(idna,i,j,k) = 3.5*rc/r_max*exp(-3.5*rc/r_max)*d0 !/cosh((rc/r_max)**25.0)
            u(idna,i,j,k) = 9.5*(rc/rsup1)**2*exp(-3.5*(rc/rsup1)**2)*d0 !/cosh((rc/r_max)**25.0)
!            u(idna(1),i,j,k) = d0/cosh(min((rc/rsup2)**25.0,690.0))
            u(idna(2),i,j,k) = d0/100.0/cosh(min((rc/rsup2)**25.0,690.0))
          endif
	  u(idna,i,j,k) = max(u(idna,i,j,k), smalld)
#ifdef KEPLER_SUPPRESSION
	  if(suppkind .eq. 'sharp') then
	    if(rc .le. rsup2) then
	      alfsup(i,j) = alfasupp*min((rc-rsup2)/(rsup1-rsup2),1.0)
	    endif
	    if(rc .ge. rsup3) then
	      alfsup(i,j) = alfasupp*min((rc-rsup3)/(rsup4-rsup3),1.0)
	    endif
	  endif
	  if(suppkind .eq. 'smooth') then
	    alfsup(i,j) = alfasupp*(tanh((rc-rsup4)*12) + tanh((-rc+rsup2)*20)+2.)/2.
	  endif
	  if(suppkind .eq. 'smout') then
	    alfsup(i,j) = alfasupp*(tanh((rc-rsup4)*12) + 1.)/2.
	  endif
	  if(suppkind .eq. 'smpeak') then
	    alfsup(i,j) = alfasupp*((tanh((rc-rsup2)*12) + 1.)/2.)*((tanh((-rc+rsup3)*12)+1.)/2.)
	  endif
	  if(suppkind .eq. 'full') then
	    alfsup = alfasupp
	  endif
#endif /* KEPLER_SUPPRESSION */
        enddo
      enddo
    enddo

    open(44,file='rotacja')
      sfq=0.5
      do i = 2,nx-1   ! 2d
        rc= x(i)*sqrt(2.0)
        xgradgp=(gp(i+1,i+1,max(nz/2,1))-gp(i-1,i-1,max(nz/2,1)))*sfq/dx/sqrt(2.)
        xgradp = -0.5*(u(idna(1),i+1,i+1,max(nz/2,1))-u(idna(1),i-1,i-1,max(nz/2,1)))/dx &
                           /sqrt(2.)*c_si**2/gamma(1)/u(idna(1),i,i,max(nz/2,1))
        omega(i)=sqrt(abs((xgradgp-xgradp)/rc))
        write(44,*) rc,omega(i),gp(i,i,max(nz/2,1))
      enddo
    close(44)
    do j = 1,ny
      yj = y(j)
      do i = 1,nx
        xi = x(i)
        rc = sqrt(xi**2+yj**2)
        
        if(dimensions .eq.'3d') then
          call hydrostatic_zeq(i, j, d0, dprof)    
        endif
                          
        do k = 1,nz
                  
          vz = 0.0
                  
!          if(dimensions .eq.'3d') then
!            u(idna,i,j,k) = min((rc/r_grav)**n_gravr,100.0)
!            u(idna,i,j,k) = dprof(k)/cosh(u(idna,i,j,k))
!          else
!            u(idna,i,j,k) = min((rc/r_grav)**n_gravr,100.0)
!            u(idna,i,j,k) = 3.5*rc/r_max*exp(-3.5*rc/r_max)*d0 !/cosh((rc/r_max)**25.0)
!          endif
!          u(idna(2),i,j,k) = u(idna(2),i,j,k)/100.0
!	   u(idna,i,j,k) = max(u(idna,i,j,k), smalld)

          call grad_factor(sfq,iu,id,i,nx)
          xgradgp=(gp(iu,j,k)-gp(id,j,k))*sfq/dx
#ifdef ISO
          xgradp = 0.0
#else /* ISO */
          xgradp = -sfq*(u(idna(1),iu,j,k)-u(idna(1),id,j,k))/dx !*c_si**2/gamma(1)/u(1,i,j,k)
#endif /* ISO */
          call grad_factor(sfq,ju,jd,j,ny)
          ygradgp=(gp(i,ju,k)-gp(i,jd,k))*sfq/dy
#ifdef ISO
          ygradp = 0.0
#else /* ISO */
          ygradp = -sfq*(u(1,i,ju,k)-u(1,i,jd,k))/dy !*c_si**2/gamma(1)/u(1,i,j,k)
#endif /* ISO */
          gradp = sqrt(xgradp**2+ygradp**2)
          gradgp = sqrt(xgradgp**2+ygradgp**2)
!          iOmega=sqrt(abs(sqrt((xgradgp+xgradp)**2+(ygradgp+ygradp)**2))/rc)
          iOmega=sqrt(abs(gradgp-gradp)/rc)
     
          ilook = (rc-xmin)/dx/sqrt(2.) + 0.5 + nb
!          if(int(ilook) .lt. nxt) then
!          iOmega = omega(int(ilook))+(omega(int(ilook)+1)-omega(int(ilook)))*(ilook-int(ilook))
          iOmega = omega(int(ilook))+(rc-x(int(ilook))*sqrt(2.))*(omega(int(ilook)+1)-omega(int(ilook))) &
	              /(x(int(ilook)+1)-x(int(ilook)))/sqrt(2.)
!	 else
!	   iOmega = omega(nxt)
!	 endif

	     
!           if((rc .le. r_max*0.6) .and. (rc .ge. r_max*0.4)) then
!	     collfaq = collf*iOmega/u(idna(1),i,j,k)
!	     write(*,*) 'iOmega = ',iOmega!*sek,' sek'
!	     write(*,*) 'iOmega = ',iOmega!*year,' yrs'
!	     write(*,*) 'dens = ',u(idna(1),i,j,k)!*cm**3/gram
!	     write(*,*) 'dcol = ',d0!*cm**2/gram
!	     write(*,*) 'rc   = ',rc!/cm,' cm'
!	     write(*,*) 'rc   = ',rc!/au,' AU'
!	   endif
	  u(imxa,i,j,k)=-iOmega*yj*u(idna,i,j,k)
          u(imya,i,j,k)= iOmega*xi*u(idna,i,j,k)
          u(imza,i,j,k) = 0.0
#ifdef KEPLER_SUPPRESSION
	  omx0(:,i,j,k)=-iOmega*yj
	  omy0(:,i,j,k)= iOmega*xi
#endif /* KEPLER_SUPPRESSION */
                  
          vx = sqr_gm * (-yj)/(rc**2+r_smooth**2)**0.75
          vy = sqr_gm * ( xi)/(rc**2+r_smooth**2)**0.75
                  
	  u(imxa(2),i,j,k) = vx*u(idna(2),i,j,k)
          u(imya(2),i,j,k) = vy*u(idna(2),i,j,k)
#ifdef KEPLER_SUPPRESSION
	  omx0(2,i,j,k) = vx
	  omy0(2,i,j,k) = vy
#endif /* KEPLER_SUPPRESSION */
!          u(imza(2),i,j,k) = vz*u(idna(2),i,j,k)      
#ifndef ISO	  
          u(iena(fadiab),i,j,k) = c_si**2/(gamma(fadiab)-1.0)*u(idna(fadiab),i,j,k)
          u(iena,i,j,k) = max(u(iena,i,j,k), smallei)
          u(iena(fadiab),i,j,k) = u(iena(fadiab),i,j,k) +0.5*(vx**2+vy**2+vz**2)*u(idna(fadiab),i,j,k)
#endif /* ISO */
          if(trim(mag_field_orient) .eq. 'toroidal') then
            b(1,i,j,k)   = -b0*sqrt(u(idna(1),i,j,k)/d0)*yj/rc
            b(2,i,j,k)   =  b0*sqrt(u(idna(1),i,j,k)/d0)*xi/rc
            b(3,i,j,k)   =  0.0
          else if(trim(mag_field_orient) .eq. 'vertical') then
            b(1,i,j,k)   =  0.0
            b(2,i,j,k)   =  0.0
            b(3,i,j,k)   =  b0
          else if(trim(mag_field_orient) .eq. 'none') then
            b(1,i,j,k)   =  0.0
            b(2,i,j,k)   =  0.0
            b(3,i,j,k)   =  0.0
          endif

#ifndef ISO	  
          u(iena(fmagn),i,j,k)   = u(iena(fmagn),i,j,k) +0.5*sum(b(:,i,j,k)**2,1)
#endif /* ISO */
        enddo
      enddo
    enddo
    
    
    allocate(bndxrar(nu,nb,ny,nz))
    allocate(bndyrar(nu,nx,nb,nz))
    do i=1,nb
    bndxrar(:,i,:,:)=u(:,nx-nb+i,:,:)
    bndyrar(:,:,i,:)=u(:,:,ny-nb+i,:)
    enddo

    
    deallocate(dprof,gpot)
    
    return
  end subroutine init_prob  

!-------------------------------------------------------
	   subroutine grad_factor(f,u,d,i,n)
	   implicit none
	   integer i,u,d,n
	   real f
	   
	   if(i .ne. 1 .and. i .ne. n) then
	     u = i+1
	     d = i-1
	     f = 0.5
	   else
	     if(i .eq. 1) then
	       u = i+1
	       d = i
	       f = 1.0
	     elseif(i .eq. n) then
	       u = i
	       d = i-1
	       f = 1.0
	     endif
           endif
	   return
	   end subroutine grad_factor
	   
  

end module init_problem

