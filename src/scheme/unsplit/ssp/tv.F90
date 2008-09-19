! $Id$
#include "mhd.def"
module tv   ! unsplit ssp
!  use constants
!  use start
  contains
  subroutine tvdb(vibj,b,vg,n,dt)
    use func, only: tvdb_emf
    integer, intent(in) :: n
    real, intent(in)    :: dt
    real, dimension(n)  :: vibj,b,vg
! locals
    real, dimension(n) :: vh

  ! unlike the B field, the vibj lives on the right cell boundary
    vh = 0.0
    vh(1:n-1) =(vg(1:n-1)+ vg(2:n))*0.5;     vh(n) = vh(n-1)

    vibj = tvdb_emf(vh,vg,b,dt)

  end subroutine tvdb

   subroutine relaxing_tvd(u,bb,sweep,i1,i2,dx,n,dt,fx,cx)
    use fluxes
    use arrays, only : nu,idna,imxa,imya,imza,iena
#ifdef GRAV
   use gravity, only : grav_pot2accel
#endif /* GRAV */
#ifdef GLOBAL_FR_SPEED
   use time, only : c
#endif /* GLOBAL_FR_SPEED */
#ifdef SHEAR
   use start, only : omega,qshear
   use arrays, only : xr
#endif /* SHEAR */

    implicit none
    integer i1,i2, n
    real :: dt,dx, dtx
    real, dimension(nu,n), optional :: fx,cx
    real, dimension(nu,n) :: u,cfr,ul,ur
    real, dimension(3,n)  :: bb
    real, dimension(n)    :: gravl,gravr,rotfr,vxr
    real, dimension(n)    :: dgrp,dgrm,dglp,dglm
           
    character sweep*6
          
!locals    
    real, dimension(nu,n) :: w,fr,fl,flux,dfrp,dfrm,dflm,dflp, &
                             dulf,durf
#ifdef GRAV
    real, dimension(nu,n) :: duls,durs
#endif /* GRAV */

#ifdef COSM_RAYS
    real q
    real, dimension(n)    :: divv,decr,tmp
#endif /* COSM_RAYS */
                            
    w(:,:)    = 0.0
    cfr  = 0.0
    dtx = dt/dx
#ifdef GRAV
    duls = 0.0 
    durs = 0.0
#endif /* GRAV */
    
    if (.not.present(fx)) then
       call mhdflux(w,cfr,u,bb,n)
    else
#ifndef GLOBAL_FR_SPEED
       cfr = cx
#else /* ~GLOBAL_FR_SPEED */
       cfr(:,:) = c
#endif /* ~GLOBAL_FR_SPEED */
       w   = fx
    endif

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
    dulf(:,2:n) = dtx*(fl(:,2:n) - fl(:,1:n-1));     dulf(:,1) = dulf(:,n)  ! dulf = (fl-cshift(fl,shift=-1,dim=2))/dx*dt

#ifdef SHEAR
    vxr(1:n-1) = 0.5*(u(imya,2:n)/u(idna,2:n) + u(imya,1:n-1)/u(idna,1:n-1))
    if(sweep .eq. 'xsweep') then
      rotfr(:) =   2.0*omega*(vxr(:) + qshear*omega*xr(:))
    else if(sweep .eq. 'ysweep')  then
      rotfr(:) = - 2.0*omega*vxr(:)
    else
      rotfr(:) = 0.0
    endif
#else /* ~SHEAR */ 
      rotfr(:) = 0.0
#endif /* ~SHEAR */ 

    
#ifdef GRAV

! Gravity source terms -------------------------------------

      call grav_pot2accel(sweep,i1,i2, n, gravr)
#ifdef SHEAR
      rotfr(n) = rotfr(n-1)
      gravr = gravr + rotfr
#endif /* SHEAR */
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
#endif /* ISO */
      duls(imxa,:)  = gravr*ul(idna,:)*dt
      durs(imxa,:)  = gravl*ur(idna,:)*dt

! ----------------------------------------------------------
#endif /* GRAV */

   ur = -durf
   ul =  dulf

#ifdef GRAV
      ur = ur + durs
      ul = ul + duls
#endif /* GRAV */

   u = ur + ul

    return

  end subroutine relaxing_tvd

!========================================================================

  subroutine initials
    use arrays, only : ui,u,bi,b
    implicit none
    ui(:,:,:,:) = u(:,:,:,:)
    bi(:,:,:,:) = b(:,:,:,:)

  end subroutine initials

!========================================================================

  subroutine integrate
    use start, only : istep, integration_order, smalld, smallei
    use fluid_boundaries, only : compute_u_bnd,bnd_u
    use mag_boundaries, only :   compute_b_bnd
    use arrays, only : ui,u,bi,b,Lu,Lb,nx,ny,nz
    use arrays, only : idna,imxa,imya,imza,ibx,iby,ibz,nu,iena
    implicit none
#ifndef ISO
  real, allocatable :: ekin(:,:,:), emag(:,:,:), eint(:,:,:)
  allocate(ekin(nx,ny,nz),emag(nx,ny,nz),eint(nx,ny,nz))
#endif /* ISO */
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
#endif /* ISO */
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
#endif /* ISO */

  end subroutine integrate  

end module tv
