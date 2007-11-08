#include "mhd.def" 
module tv ! split orig
  use start
  use arrays, only : idna,imxa,imya,imza,iecr,ibx,iby,ibz,nu,iena,wa
  contains
  subroutine tvdb(vibj,b,vg,n,dt,di)
    integer i,n,ip,ipp,im
    real :: dt,di,dti
    real, dimension(n) :: vibj,b,vg
! locals
    real, dimension(n) :: b1,vibj1,vh
    real w,wp,wm,dw,v

  ! unlike the B field, the vibj lives on the right cell boundary
    vh = 0.0
    vh(1:n-1) =(vg(1:n-1)+ vg(2:n))*0.5;     vh(n) = vh(n-1)

    dti = dt/di

    where(vh > 0.)
      vibj1=b*vg
    elsewhere
      vibj1=eoshift(b*vg,1,boundary=big)
    end where
    b1(2:n) = b(2:n) -(vibj1(2:n)-vibj1(1:n-1))*dti*0.5;    b1(1) = b(2)

    do i=3,n-3
      ip=i+1
      ipp=ip+1
      im=i-1
      v=vh(i)
      if (v .gt. 0.) then
        w=vg(i)*b1(i)
        wp=(vg(ip)*b1(ip)-w)*0.5
        wm=(w-vg(im)*b1(im))*0.5
      else
        w=vg(ip)*b1(ip)
        wp=(w-vg(ipp)*b1(ipp))*0.5
        wm=(vg(i)*b1(i)-w)*0.5
      end if
      dw=0.
      if(wm*wp .gt. 0.) dw=2.*wm*wp/(wm+wp)
      vibj(i)=(w+dw)*dt
    end do
  end subroutine tvdb

   subroutine relaxing_tvd(u,bb,sweep,i1,i2,dx,n,dt)
 ! Cooling and heating implemented following Rafal Kosinski
#ifdef GRAV
   use gravity, only :grav_pot2accel
#endif GRAV
   use grid, only : xr

    implicit none
    integer i1,i2, n
    real ::  dt,dx, dtx
    real, dimension(nu,n) :: u,cfr,ul,ur
    real, dimension(3,n)  :: bb
    real, dimension(n)    :: gravl,gravr,rotfl,rotfr,vx,vxr
    real, dimension(n)    :: dgrp,dgrm,dglp,dglm
           
    character sweep*6
          
!locals    
    real, dimension(nu,n) :: w,fr,fl,flux,dfrp,dfrm,dflm,dflp, &
                             dulf,durf
    real, dimension(nu,n) :: ul0,ur0,u1,ul1,ur1
    real, dimension(n)    :: ekin,eint,emag

#ifdef GRAV
    real, dimension(nu,n) :: duls,durs
#endif GRAV
#ifdef COSM_RAYS
    real, dimension(n)    :: divv,decr,gpcr,ecr,tmp
#endif COSM_RAYS
                            
    w         = 0.0
    cfr       = 0.0
    dtx       = dt / dx

#ifdef GRAV
    duls = 0.0 
    durs = 0.0
#endif GRAV
    
    u1 = u

  do istep=1,integration_order

    call mhdflux(w,cfr,u1,bb,n)

    fr = (u1*cfr+w)*0.5
    fl = (u1*cfr-w)*0.5
    if(istep == 1) then
      ur0 = fr/cfr
      ul0 = fl/cfr
    endif

    fl(:,1:n-1) = fl(:,2:n); fl(:,n) = fl(:,1)  
    !fl =  cshift(fl,shift=1,dim=2)

    if(istep == 2) then
       dfrp(:,1:n-1) = 0.5*(fr(:,2:n) - fr(:,1:n-1)); dfrp(:,n) = dfrp(:,1)       ! dfrp=(cshift(fr,shift=1,dim=2)-fr)*0.5      
       dfrm(:,2:n)   = dfrp(:,1:n-1);                 dfrm(:,1) = dfrm(:,n)       ! dfrm=cshift(dfrp,shift=-1,dim=2)            
      call flimiter(fr,dfrm,dfrp,nu,n)

       dflp(:,1:n-1) = 0.5*(fl(:,1:n-1) - fl(:,2:n)); dflp(:,n) = dflp(:,1)       ! dflp=(fl-cshift(fl,shift=1,dim=2))*0.5      
       dflm(:,2:n)   = dflp(:,1:n-1);                 dflm(:,1) = dflm(:,n)       ! dflm=cshift(dflp,shift=-1,dim=2)            
      call flimiter(fl,dflm,dflp,nu,n)
    endif

    durf(:,2:n) = dtx*(fr(:,2:n) - fr(:,1:n-1));     durf(:,1) = durf(:,n)        ! durf = (fr-cshift(fr,shift=-1,dim=2))/dx*dt 
    dulf(:,2:n) = dtx*(fl(:,2:n) - fl(:,1:n-1));     dulf(:,1) = durf(:,n)        ! dulf = (fl-cshift(fl,shift=-1,dim=2))/dx*dt
    
    ur1(:,:) = ur0 - cn(istep)*durf
    ul1(:,:) = ul0 + cn(istep)*dulf

    ur1(1,:) = max(ur1(1,:), smalld)
    ul1(1,:) = max(ul1(1,:), smalld)

#ifdef SHEAR
    vxr(1:n-1) = 0.5*(u1(imya,2:n)/u1(idna,2:n) + u1(imya,1:n-1)/u1(idna,1:n-1))
    if(sweep .eq. 'xsweep') then
      rotfr(:) =   2.0*omega*(vxr(:) + qshear*omega*xr(:))  !!! x(:) czy xr(:) ? 
    else if(sweep .eq. 'ysweep')  then
      rotfr(:) = - 2.0*omega*vxr(:)
    else
      rotfr(:) = 0.0
    endif
!    rotfl(2:n) = rotfr(1:n-1)                                                   ! gravl      = cshift(gravr,-1) 
!    rotfl(1)   = rotfl(2)
!    rotfr(n)   = rotfr(n-1)
#else  SHEAR
    rotfr(:) = 0.0
!    rotfl(:) = 0.0
#endif SHEAR

! Gravity source terms -------------------------------------
#ifdef GRAV
    call grav_pot2accel(sweep,i1,i2, n, gravr)
    gravr = gravr + rotfr
    gravl(2:n) = gravr(1:n-1)                                                   ! gravl      = cshift(gravr,-1) 
    gravl(1)   = gravl(2)
    gravr(n)   = gravr(n-1)
 
    if(istep == 2) then
      dgrp(1:n-1) = 0.5*(gravr(1:n-1) - gravr(2:n));  dgrp(n) = dgrp(1)         ! dgrp = (gravr-cshift(gravr,shift=1))*0.5
      dgrm(2:n) = dgrp(1:n-1);    dgrm(1) = dgrm(n)                             ! dgrm = (cshift(gravr,shift=-1)-gravr)*0.5
      call flimiter(gravr,dgrp,dgrm,1,n)

      dglp(1:n-1) = 0.5*(gravl(2:n) - gravl(1:n-1));   dglp(n) = dglp(1)        ! dglp = (cshift(gravl,shift=1)-gravl)*0.5
      dglm(2:n)   = dglp(1:n-1);  dglm(1) = dglm(n)                             ! dglm = (gravl-cshift(gravl,shift=-1))*0.5
      call flimiter(gravl,dglp,dglm,1,n)
    endif

!    gravr = gravr + rotfr
!    gravl = gravl + rotfl
     

#ifndef ISO
    duls(iena,:)  = gravr*ul0(imxa,:)*dt 
    durs(iena,:)  = gravl*ur0(imxa,:)*dt
#endif ISO
    duls(imxa,:)  = gravr*ul0(idna,:)*dt
    durs(imxa,:)  = gravl*ur0(idna,:)*dt
    ur1= ur1 + cn(istep)*durs
    ul1= ul1 + cn(istep)*duls
#endif GRAV
! ----------------------------------------------------------

    u1 = ul1 + ur1
    u1(1,:) = max(u1(1,:), smalld)

! ----------------------------------------------------------
          
#ifdef COSM_RAYS
    select case (sweep)
      case('xsweep') 
        divv = wa(:,i1,i2)
      case('ysweep') 
        divv = wa(i2,:,i1)
      case('zsweep') 
        divv = wa(i1,i2,:)
    end select

    decr(:)  = -(gamma_cr-1.)*u1(iecr,:)*divv(:)*dt 
    u1(iecr,:) = u1(iecr,:) + cn(istep)*decr(:)

    vx = u1(imxa,:)/u1(idna,:)
    ecr = u1(iecr,:)

    gpcr = (gamma_cr -1.)*(cshift(ecr,1)-cshift(ecr,-1))/(2.*dx)

!      where(vx > 0.)
!        gpcr = (gamma_cr -1.)*(ecr-cshift(ecr,-1))/dx
!      elsewhere
!        gpcr = (gamma_cr -1.)*(cshift(ecr,1)-ecr)/dx
!      end where
      
!    gpcr = 0.25*cshift(gpcr,-1)+0.5*gpcr+0.25*cshift(gpcr,1)

#ifndef ISO
    u1(iena,:) = u1(iena,:) - cn(istep)*u1(imxa,:)/u1(idna,:)*gpcr*dt
#endif ISO
    u1(imxa,:) = u1(imxa,:) - cn(istep)*gpcr*u1(idna,:)*dt

#endif COSM_RAYS
#ifndef ISO
    ekin = 0.5*(u1(imxa,:)*u1(imxa,:)+u1(imya,:)*u1(imya,:)+u1(imza,:)*u1(imza,:))/u1(idna,:)
    emag = 0.5*(bb(ibx,:)*bb(ibx,:) + bb(iby,:)*bb(iby,:) + bb(ibz,:)*bb(ibz,:)) 
    eint = u1(5,:)-ekin-emag
    eint = max(eint,smallei)
    u1(5,:) = eint+ekin+emag
#endif ISO
    
      u(:,:) = u1(:,:)
  enddo

    return

  end subroutine relaxing_tvd

!==========================================================================================
 
  subroutine mhdflux(flux,cfr,u,b,n)
    implicit none
    integer n
    real, dimension(nu,n)::flux,u,cfr
    real, dimension(3,n):: b
! locals
    real, dimension(n) :: vx,vy,vz,vt, ps,p,pmag

    flux   = 0.0
    vx  = 0.0
    vy  = 0.0
    vz  = 0.0
    
    cfr = 0.0 

    pmag(2:n-1)=(b(ibx,2:n-1)*b(ibx,2:n-1)+b(iby,2:n-1)*b(iby,2:n-1)+b(ibz,2:n-1)*b(ibz,2:n-1))*0.5
    vx(2:n-1)=u(imxa,2:n-1)/u(idna,2:n-1)

! UWAGA ZMIANA W OBLICZANIU LOKALNEJ PREDKOSCI NA PROBE
    vy(2:n-1)=u(imya,2:n-1)/u(idna,2:n-1)
    vz(2:n-1)=u(imza,2:n-1)/u(idna,2:n-1)
    vt = sqrt(vx*vx+vy*vy+vz*vz)
    
    
#ifdef ISO
    p(2:n-1) = csi2*u(idna,2:n-1)
    ps(2:n-1)= p(2:n-1) + pmag(2:n-1)
#else ISO    
    ps(2:n-1)=(u(iena,2:n-1)-(u(imxa,2:n-1)*u(imxa,2:n-1)+u(imya,2:n-1)*u(imya,2:n-1) & 
              +u(imza,2:n-1)*u(imza,2:n-1))/u(idna,2:n-1)*0.5)*(gamma-1)+(2-gamma)*pmag(2:n-1)
    p(2:n-1)=ps(2:n-1)-pmag(2:n-1)
#endif ISO
    flux(idna,2:n-1)=u(imxa,2:n-1)
    flux(imxa,2:n-1)=u(imxa,2:n-1)*vx(2:n-1)+ps(2:n-1)-b(ibx,2:n-1)*b(ibx,2:n-1)
    flux(imya,2:n-1)=u(imya,2:n-1)*vx(2:n-1)-b(iby,2:n-1)*b(ibx,2:n-1)
    flux(imza,2:n-1)=u(imza,2:n-1)*vx(2:n-1)-b(ibz,2:n-1)*b(ibx,2:n-1)
#ifndef ISO
    flux(iena,2:n-1)=(u(iena,2:n-1)+ps(2:n-1))*vx(2:n-1)-b(ibx,2:n-1)*(b(ibx,2:n-1)*u(imxa,2:n-1) &
                +b(iby,2:n-1)*u(imya,2:n-1)+b(ibz,2:n-1)*u(imza,2:n-1))/u(idna,2:n-1)
#endif ISO
#ifdef COSM_RAYS
    flux(iecr,2:n-1)= u(iecr,2:n-1)*vx(2:n-1)
#endif COSM_RAYS
#ifdef LOCAL_FR_SPEED

!       The freezing speed is now computed locally (in each cell) 
!       as in Trac & Pen (2003). This ensures much sharper shocks, 
!       but sometimes may lead to numerical instabilities   
#ifdef ISO
! UWAGA ZMIANA W OBLICZANIU LOKALNEJ PREDKOSCI NA PROBE
!        cfr(1,2:n-1) = abs(vx(2:n-1)) &
        cfr(1,2:n-1) = abs(vt(2:n-1)) &
                      +max(sqrt( abs(2*pmag(2:n-1) + p(2:n-1))/u(idna,2:n-1)),small)
#else ISO   
! UWAGA ZMIANA W OBLICZANIU LOKALNEJ PREDKOSCI NA PROBE
!        cfr(1,2:n-1) = abs(vx(2:n-1)) &
        cfr(1,2:n-1) = abs(vt(2:n-1)) &
                      +max(sqrt( abs(2*pmag(2:n-1) + gamma*p(2:n-1))/u(idna,2:n-1)),small)
#endif ISO
        cfr(1,1) = cfr(1,2)
        cfr(1,n) = cfr(1,n-1)   
        cfr = spread(cfr(1,:),1,nu)
#endif LOCAL_FR_SPEED
#ifdef GLOBAL_FR_SPEED
!       The freezing speed is now computed globally 
!       (c=const for the whole domain) in sobroutine 'timestep' 
!       Original computation of the freezing speed was done
!       for each sweep separately:
!       c=maxval(abs(vx(nb-2:n-1))+sqrt(abs(2*pmag(nb-2:n-1) &
!                            +gamma*p(nb-2:n-1))/u(idna,nb-2:n-1)))  
        cfr(:,:) = c
#endif GLOBAL_FR_SPEED

  end subroutine mhdflux

!==========================================================================================

  subroutine flimiter(f,a,b,m,n)
    implicit none
    integer m,n
    real, dimension(m,n) :: f,a,b,c
#ifdef VANLEER
!   'vanleer' flux limiter is used

      c = a*b
      where (c .gt. 0)                                        
        f=f+2*c/(a+b)
      endwhere      
#endif VANLEER
#ifdef MINMOD
!   'minmod' flux limiter is used

      f = f+(sign(1.,a)+sign(1.,b))*min(abs(a),abs(b))*0.5
#endif MINMOD
#ifdef SUPERBEE
!   'superbee' flux limiter is used

      where (abs(a) .gt. abs(b))
        f = f+(sign(1.,a)+sign(1.,b))*min(abs(a), abs(2*b))*0.5
      elsewhere
        f = f+(sign(1.,a)+sign(1.,b))*min(abs(2*a), abs(b))*0.5
      endwhere
#endif SUPERBEE
 
    return
  end subroutine flimiter
end module tv
