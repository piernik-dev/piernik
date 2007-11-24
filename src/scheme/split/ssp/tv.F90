#include "mhd.def"
module tv   ! split ssp
   use constants
   use start
   use arrays, only : idna,imxa,imya,imza,ibx,iby,ibz,nu,iena
   contains
   subroutine tvdb(vibj,b,vg,n,dt)
    use func, only : tvdb_emf
    integer i,n,ip,ipp,im
    real dt
    real, dimension(n) :: vibj,b,vg
! locals
    real, dimension(n) :: vh
    real w,wp,wm,dw,v

  ! unlike the B field, the vibj lives on the right cell boundary
    vh = 0.0
    vh(1:n-1) =(vg(1:n-1)+ vg(2:n))*0.5;     vh(n) = vh(n-1)
    
    vibj = tvdb_emf(vh,vg,b,dt)

  end subroutine tvdb

   subroutine relaxing_tvd(u,bb,sweep,i1,i2,dx,n,dt)

#ifdef GRAV
   use gravity, only :grav_pot2accel
#endif GRAV
   use grid, only : x

    implicit none
    integer i1,i2, n
    real :: dt,dx,dtx
    real, dimension(nu,n) :: u,cfr,ul,ur
    real, dimension(3,n)  :: bb
    real, dimension(n)    :: gravl,gravr,shear,vxr
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
    real :: q
    real, dimension(n)    :: divv,decr,tmp
#endif COSM_RAYS
                            
    w    = 0.0
    cfr  = 0.0
    dtx  = dt / dx

#ifdef GRAV
    duls = 0.0 
    durs = 0.0
#endif GRAV
    
    u1 = u

   ur1 = 0.0
   ul1 = 0.0

  do istep=1,integration_order

    call mhdflux(w,cfr,u1,bb,n)

    fr = (u1*cfr+w)*0.5
    fl = (u1*cfr-w)*0.5
    ur = fr/cfr
    ul = fl/cfr

    if(istep .eq. 1) then
      ur0 = ur
      ul0 = ul
    endif

    fl(:,1:n-1) = fl(:,2:n);                       fl(:,n) = fl(:,1)        ! fl=cshift(fl,shift=1,dim=2)                 

    dfrp(:,1:n-1) = 0.5*(fr(:,2:n) - fr(:,1:n-1)); dfrp(:,n) = dfrp(:,1)    ! dfrp=(cshift(fr,shift=1,dim=2)-fr)*0.5    
    dfrm(:,2:n)   = dfrp(:,1:n-1);                 dfrm(:,1) = dfrm(:,n)    ! dfrm=cshift(dfrp,shift=-1,dim=2)          
    call flimiter(fr,dfrm,dfrp,nu,n)
    
    dflp(:,1:n-1) = 0.5*(fl(:,1:n-1) - fl(:,2:n)); dflp(:,n) = dflp(:,1)    ! dflp=(fl-cshift(fl,shift=1,dim=2))*0.5    
    dflm(:,2:n)   = dflp(:,1:n-1);                 dflm(:,1) = dflm(:,n)    ! dflm=cshift(dflp,shift=-1,dim=2)          
    call flimiter(fl,dflm,dflp,nu,n)

    durf(:,2:n) = dtx*(fr(:,2:n) - fr(:,1:n-1));     durf(:,1) = durf(:,n)  ! durf = (fr-cshift(fr,shift=-1,dim=2))/dx*dt 
    dulf(:,2:n) = dtx*(fl(:,2:n) - fl(:,1:n-1));     dulf(:,1) = durf(:,n)  ! dulf = (fl-cshift(fl,shift=-1,dim=2))/dx*dt
    
    ur1= cn(1,istep)*ur0 + cn(2,istep)*(ur1 - durf)      
    ul1= cn(1,istep)*ul0 + cn(2,istep)*(ul1 + dulf)

    ur1(1,:) = max(ur1(1,:), smalld)
    ul1(1,:) = max(ul1(1,:), smalld)
    
#ifdef SHEAR
     vxr(1:n-1) = 0.5*(u(imya,2:n)/u(idna,2:n) + u(imya,1:n-1)/u(idna,1:n-1))
     if(sweep .eq. 'xsweep') then
       shear(:) =   2.0*omega*(vxr(:) + qshear*omega*x(:))
     else if(sweep .eq. 'ysweep')  then
       shear(:) = - 2.0*omega*vxr(:)
     else
       shear = 0.0
     endif
#else SHEAR
     shear = 0.0
#endif SHEAR
! Gravity source terms -------------------------------------
#ifdef GRAV
      call grav_pot2accel(sweep,i1,i2, n, gravr)
      gravr      = gravr + shear
      gravl(2:n) = gravr(1:n-1)                                             ! gravl      = cshift(gravr,-1) 
      gravl(1)   = gravl(2)
      gravr(n)   = gravr(n-1)

      dgrp(1:n-1) = 0.5*(gravr(1:n-1) - gravr(2:n));  dgrp(n) = dgrp(1)     ! dgrp = (gravr-cshift(gravr,shift=1))*0.5
      dgrm(2:n) = dgrp(1:n-1);                        dgrm(1) = dgrm(n)     ! dgrm = (cshift(gravr,shift=-1)-gravr)*0.5
      call flimiter(gravr,dgrp,dgrm,1,n)

      dglp(1:n-1) = 0.5*(gravl(2:n) - gravl(1:n-1));  dglp(n) = dglp(1)     ! dglp = (cshift(gravl,shift=1)-gravl)*0.5
      dglm(2:n)   = dglp(1:n-1);                      dglm(1) = dglm(n)     ! dglm = (gravl-cshift(gravl,shift=-1))*0.5
      call flimiter(gravl,dglp,dglm,1,n)
 
#ifndef ISO
      duls(iena,:)  = gravr*ul(imxa,:)*dt 
      durs(iena,:)  = gravl*ur(imxa,:)*dt
#endif ISO
      duls(imxa,:)  = gravr*ul(idna,:)*dt
      durs(imxa,:)  = gravl*ur(idna,:)*dt
      ur1= ur1 + cn(2,istep)*durs                 
      ul1= ul1 + cn(2,istep)*duls             
#endif GRAV
! ----------------------------------------------------------
    u1 = ul1 + ur1
    u1(1,:) = max(u1(1,:), smalld)
#ifdef SHEAR
     if(sweep .eq. 'xsweep') then
        u1(imxa,:) = u1(imxa,:) + 2.0*omega*(u(imya,:) + qshear*omega*x(:)*u(idna,:)) * dt! * cn(istep)
     else if(sweep .eq. 'ysweep')  then
        u1(imxa,:) = u1(imxa,:) - 2.0*omega*u(imya,:) * dt!* cn(istep)
    endif
#endif
          
#ifndef ISO
    ekin = 0.5*(u1(imxa,:)*u1(imxa,:)+u1(imya,:)*u1(imya,:)+u1(imza,:)*u1(imza,:))/u1(idna,:)
    emag = 0.5*(bb(ibx,:)*bb(ibx,:) + bb(iby,:)*bb(iby,:) + bb(ibz,:)*bb(ibz,:)) 
    eint = u1(iena,:)-ekin-emag
    eint = max(eint,smallei)
    u1(iena,:) = eint+ekin+emag
#endif ISO
    
    u = u1
  
  enddo

  return
  end subroutine relaxing_tvd

end module tv
