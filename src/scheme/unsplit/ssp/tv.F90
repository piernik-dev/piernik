#include "mhd.def"
module tv   ! unsplit ssp
  use start
  use arrays, only : idna,imxa,imya,imza,ibx,iby,ibz,nu,iena
  contains
  subroutine tvdb(vibj,b,vg,n,dt)
    integer i,n,ip,ipp,im
    real dt
    real, dimension(n) :: vibj,b,vg
! locals
    real, dimension(n) :: vh
    real w,wp,wm,dw,v

  ! unlike the B field, the vibj lives on the right cell boundary
    vh = 0.0
    vh =(vg+eoshift(vg,1,boundary=big))*0.5

    do i=3,n-3
      ip=i+1
      ipp=ip+1
      im=i-1
      v=vh(i)
      if (v .gt. 0.) then
        w=vg(i)*b(i)
        wp=(vg(ip)*b(ip)-w)*0.5
        wm=(w-vg(im)*b(im))*0.5
      else
        w=vg(ip)*b(ip)
        wp=(w-vg(ipp)*b(ipp))*0.5
        wm=(vg(i)*b(i)-w)*0.5
      end if
      dw=0.
      if(wm*wp .gt. 0.) dw=2.*wm*wp/(wm+wp)
      vibj(i)=(w+dw)*dt
    enddo
  end subroutine tvdb

   subroutine relaxing_tvd(u,bb,sweep,i1,i2,dx,n,dt)
 ! Cooling and heating implemented following Rafal Kosinski

   use gravity, only : grav_pot2accel

    implicit none
    integer i1,i2, n
    real :: dt,dx, dtx
    real, dimension(nu,n) :: u,cfr,ul,ur
    real, dimension(3,n)  :: bb
    real, dimension(n)    :: gravl,gravr
    real, dimension(n)    :: dgrp,dgrm,dglp,dglm
           
    character sweep*6
          
!locals    
    real, dimension(nu,n) :: w,fr,fl,flux,dfrp,dfrm,dflm,dflp, &
                             dulf,durf
#ifdef GRAV
    real, dimension(nu,n) :: duls,durs
#endif

#ifdef COSM_RAYS
    real q
    real, dimension(n)    :: divv,decr,tmp
#endif COSM_RAYS
                            
    w(:,:)    = 0.0
    cfr  = 0.0
    dtx = dt/dx
#ifdef GRAV
    duls = 0.0 
    durs = 0.0
#endif GRAC
    
    call mhdflux(w,cfr,u,bb,n)

    fr = (u*cfr+w)*0.5
    fl = (u*cfr-w)*0.5
    ur = fr/cfr
    ul = fl/cfr
    
    fl(:,1:n-1) = fl(:,2:n);                       fl(:,n) = fl(:,1)        ! fl=cshift(fl,shift=1,dim=2)                 

    dfrp(:,1:n-1) = 0.5*(fr(:,2:n) - fr(:,1:n-1)); dfrp(:,n) = dfrp(:,1)    ! dfrp=(cshift(fr,shift=1,dim=2)-fr)*0.5    
    dfrm(:,2:n)   = dfrp(:,1:n-1);                 dfrm(:,1) = dfrm(:,n)    ! dfrm=cshift(dfrp,shift=-1,dim=2)          
    call flimiter(fr,dfrm,dfrp,nu,n)
   
    dflp(:,1:n-1) = 0.5*(fl(:,1:n-1) - fl(:,2:n)); dflp(:,n) = dflp(:,1)    ! dflp=(fl-cshift(fl,shift=1,dim=2))*0.5    
    dflm(:,2:n)   = dflp(:,1:n-1);                 dflm(:,1) = dflm(:,n)    ! dflm=cshift(dflp,shift=-1,dim=2)
    call flimiter(fl,dflm,dflp,nu,n)

    durf(:,2:n) = dtx*(fr(:,2:n) - fr(:,1:n-1));     durf(:,1) = durf(:,n)  ! durf = (fr-cshift(fr,shift=-1,dim=2))/dx*dt 
    dulf(:,2:n) = dtx*(fl(:,2:n) - fl(:,1:n-1));     dulf(:,1) = durf(:,n)  ! dulf = (fl-cshift(fl,shift=-1,dim=2))/dx*dt
    
#ifdef GRAV

! Gravity source terms -------------------------------------

      call grav_pot2accel(sweep,i1,i2, n, gravr)
      gravl(2:n) = gravr(1:n-1)                                             ! gravl      = cshift(gravr,-1)
      gravl(1)  = gravl(2)
      gravr(n)  = gravr(n-1)
         
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

! ----------------------------------------------------------
#endif GRAV

   ur = -durf
   ul =  dulf

#ifdef GRAV
      ur = ur + durs
      ul = ul + duls
#endif

   u = ur + ul

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

  subroutine initials
    use arrays, only : ui,u,bi,b
    implicit none
    ui(:,:,:,:) = u(:,:,:,:)
    bi(:,:,:,:) = b(:,:,:,:)

  end subroutine initials

  subroutine integrate
    use fluid_boundaries, only : compute_u_bnd
    use mag_boundaries, only :   compute_b_bnd
    use arrays, only : ui,u,bi,b,Lu,Lb,nx,ny,nz
    implicit none
#ifndef ISO
  real, allocatable :: ekin(:,:,:), emag(:,:,:), eint(:,:,:)
  allocate(ekin(nx,ny,nz),emag(nx,ny,nz),eint(nx,ny,nz))
#endif ISO

    select case(istep)
     case(1)
        u(:,:,:,:) = ui(:,:,:,:) + Lu(:,:,:,:)
     case(2)
      if(integration_order .eq. 2) then
        u(:,:,:,:) = 0.5*ui(:,:,:,:) + 0.5*(u(:,:,:,:)+Lu(:,:,:,:))
      else
        u(:,:,:,:) = 0.75*ui(:,:,:,:) + 0.25*(u(:,:,:,:)+Lu(:,:,:,:))
      endif
     case(3)
        u(:,:,:,:) = 1./3.*ui(:,:,:,:) + 2./3.*(u(:,:,:,:)+Lu(:,:,:,:))
     case default
      write(*,*) 'That high integration order not implemented'
    end select

    u(idna,:,:,:) = max(u(idna,:,:,:),smalld)

#ifndef ISO
    ekin = 0.5*(u(imxa,:,:,:)*u(imxa,:,:,:)+u(imya,:,:,:)*u(imya,:,:,:) &
               +u(imza,:,:,:)*u(imza,:,:,:))/u(idna,:,:,:)
    emag = 0.125*((b(ibx,:,:,:)+cshift(b(ibx,:,:,:),shift=1,dim=1))**2&
                 +(b(iby,:,:,:)+cshift(b(iby,:,:,:),shift=1,dim=2))**2&
                 +(b(ibz,:,:,:)+cshift(b(ibz,:,:,:),shift=1,dim=3))**2)

    eint = u(iena,:,:,:)-ekin-emag
    eint = max(eint,smallei)
    u(iena,:,:,:) = eint+ekin+emag
#endif ISO
    call compute_u_bnd

    select case(istep)
     case(1)
        b(:,:,:,:) = bi(:,:,:,:) + Lb(:,:,:,:)
     case(2)
      if(integration_order .eq. 2) then
        b(:,:,:,:) = 0.5*bi(:,:,:,:) + 0.5*(b(:,:,:,:)+Lb(:,:,:,:))
      else
        b(:,:,:,:) = 0.75*bi(:,:,:,:) + 0.25*(b(:,:,:,:)+Lb(:,:,:,:))
      endif
     case(3)
        b(:,:,:,:) = 1./3.*bi(:,:,:,:) + 2./3.*(b(:,:,:,:)+Lb(:,:,:,:))
     case default
      write(*,*) 'That high integration order not implemented'
    end select

    call compute_b_bnd

#ifndef ISO
  deallocate(ekin,emag,eint)
#endif ISO

  end subroutine integrate  

end module tv
