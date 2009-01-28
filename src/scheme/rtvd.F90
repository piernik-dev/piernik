! $Id$
#include "piernik.def"
module rtvd ! split orig
  contains

  subroutine tvdb(vibj,b,vg,n,dt,di)
    use constants, only : big
    implicit none
    integer, intent(in) :: n
    real, intent(in)    :: dt,di
    real, dimension(n)  :: vibj,b,vg
! locals
    real, dimension(n)  :: b1,vibj1,vh
    real :: dti, v, w, dw, dwm, dwp
    integer :: i, ip, ipp, im

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

    do i = 3, n-3
       ip  = i  + 1
       ipp = ip + 1
       im  = i  - 1
       v   = vh(i)
       if (v > 0.0) then
         w=vg(i)*b1(i)
         dwp=(vg(ip)*b1(ip)-w)*0.5
         dwm=(w-vg(im)*b1(im))*0.5
       else
         w=vg(ip)*b1(ip)
         dwp=(w-vg(ipp)*b1(ipp))*0.5
         dwm=(vg(i)*b1(i)-w)*0.5
       end if
       dw=0.0
       if(dwm*dwp > 0.0) dw=2.0*dwm*dwp/(dwm+dwp)
       vibj(i)=(w+dw)*dt
    enddo

  end subroutine tvdb

  subroutine relaxing_tvd(u,bb,sweep,i1,i2,dx,n,dt)
  
#ifdef IONIZED
    use fluidindex, only : i_ion
#endif /* IONIZED */  
    
    
    use start,  only : smalld, integration_order  !!! ,cn
#ifndef ISO
    use start,  only : smallei
#endif /* ISO */
    
    use fluidindex, only : nvar,nmag,nfluid
    use fluidindex, only : iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz
#ifndef ISO
    use fluidindex, only : iarr_all_en
#endif /* ISO */
    use fluidindex, only : ibx,iby,ibz

    use fluxes,  only : flimiter,grav_limit,all_fluxes

#ifdef GRAV
    use gravity, only :grav_pot2accel
#endif /* GRAV */
#ifdef SHEAR
    use grid, only : xr
    use shear,  only : qshear, omega
#endif /* SHEAR */
#ifdef COSM_RAYS
    use initcosmicrays,  only : gamma_cr, cr_active, smallecr
    use initcosmicrays,  only : iecr
    use arrays,          only : divvel 
#endif /* COSM_RAYS */

    implicit none
    
    integer :: i1,i2, n, istep, i
    
    real    :: dt,dx,dtx
    real, dimension(nvar,n) :: u,cfr
    real, dimension(nmag,n)  :: bb
#ifdef COSM_RAYS
    real, dimension(n)    :: vx
#endif /* COSM_RAYS */
#ifdef GRAV
    real, dimension(n)    :: gravl,gravr
    real, dimension(n)    :: dgrp,dgrm,dglp,dglm
#endif /* GRAV */
#if defined GRAV || defined SHEAR
    real, dimension(n)    :: rotfr
#endif /* GRAV || SHEAR */
#ifdef SHEAR
    real, dimension(size(iarr_all_my),n)    :: vxr
#endif /* SHEAR */
    character sweep*6

!locals
    real, dimension(nvar,n) :: w,fr,fl,dfrp,dfrm,dflm,dflp,dulf,durf
    real, dimension(nvar,n) :: ul0,ur0,u1,ul1,ur1
#ifndef ISO
    real, dimension(nfluid,n)    :: ekin,eint
    real, dimension(n)           :: emag
#endif /* ISO */
#ifdef GRAV
    real, dimension(nvar,n) :: duls,durs
#endif /* GRAV */
#ifdef COSM_RAYS
    real, dimension(n)    :: divv,decr,grad_pcr,ecr
#endif /* COSM_RAYS */

    real, dimension(2)    :: cn 
    cn = [0.5 , 1.]

    w         = 0.0
    cfr       = 0.0
    dtx       = dt / dx

#ifdef GRAV
    duls = 0.0
    durs = 0.0
#endif /* GRAV */

    u1 = u

  do istep=1,integration_order

    call all_fluxes(w,cfr,u1,bb,n)

!   write(*,*) 'tv:'
!   write(*,*) w(1,:)
!   write(*,*) w(2,:)
!   write(*,*) w(5,:)
!   stop

    fr = (u1*cfr+w)*0.5
    fl = (u1*cfr-w)*0.5
    if(istep == 1) then
      ur0 = fr/cfr
      ul0 = fl/cfr
    endif

    fl(:,1:n-1) = fl(:,2:n)                         ; fl(:,n)   = fl(:,n-1)

    if(istep == 2) then
       dfrp(:,1:n-1) = 0.5*(fr(:,2:n) - fr(:,1:n-1)); dfrp(:,n) = dfrp(:,n-1)
       dfrm(:,2:n)   = dfrp(:,1:n-1);                 dfrm(:,1) = dfrm(:,2)
      call flimiter(fr,dfrm,dfrp,nvar,n)

       dflp(:,1:n-1) = 0.5*(fl(:,1:n-1) - fl(:,2:n)); dflp(:,n) = dflp(:,n-1)
       dflm(:,2:n)   = dflp(:,1:n-1);                 dflm(:,1) = dflm(:,2)
      call flimiter(fl,dflm,dflp,nvar,n)
    endif

    durf(:,2:n) = dtx*(fr(:,2:n) - fr(:,1:n-1));     durf(:,1) = durf(:,2)
    dulf(:,2:n) = dtx*(fl(:,2:n) - fl(:,1:n-1));     dulf(:,1) = dulf(:,2)

    ur1(:,:) = ur0 - cn(istep)*durf
    ul1(:,:) = ul0 + cn(istep)*dulf

    ur1(iarr_all_dn,:) = max(ur1(iarr_all_dn,:), 0.5*smalld)
    ul1(iarr_all_dn,:) = max(ul1(iarr_all_dn,:), 0.5*smalld)


#ifdef SHEAR
    vxr(:,1:n-1) = 0.5*(u1(iarr_all_my,2:n)/u1(iarr_all_dn,2:n) + u1(iarr_all_my,1:n-1)/u1(iarr_all_dn,1:n-1))
    if(sweep .eq. 'xsweep') then
      rotfr(:) =   2.0*omega*(vxr(1,:) + qshear*omega*xr(:))            ! it works only for ONE FLUID now
    else if(sweep .eq. 'ysweep')  then
      rotfr(:) = - 2.0*omega*vxr(1,:)
    else
      rotfr(:) = 0.0
    endif
#else /* SHEAR */
#ifdef GRAV
    rotfr(:) = 0.0
#endif /* GRAV */
#endif /* SHEAR */

! Gravity source terms -------------------------------------
#ifdef GRAV
    call grav_pot2accel(sweep,i1,i2, n, gravr)
    gravr      = gravr + rotfr                     ;  gravr(n)   = gravr(n-1)
    gravl(2:n) = gravr(1:n-1)                      ;  gravl(1)   = gravl(2)

    if(istep == 2) then
      dgrp(1:n-1) = 0.5*(gravr(1:n-1) - gravr(2:n));  dgrp(n) = dgrp(n-1)
      dgrm(2:n) = dgrp(1:n-1)                      ;  dgrm(1) = dgrm(2)
      call flimiter(gravr,dgrp,dgrm,1,n)

      dglp(1:n-1) = 0.5*(gravl(2:n) - gravl(1:n-1));  dglp(n) = dglp(n-1)
      dglm(2:n)   = dglp(1:n-1)                    ;  dglm(1) = dglm(2)
      call flimiter(gravl,dglp,dglm,1,n)
    endif

    do i=1, size(iarr_all_mx)
#ifndef ISO
       duls(iarr_all_en(i),:)  = gravr(:)*ul0(iarr_all_mx(i),:)*dt
       durs(iarr_all_en(i),:)  = gravl(:)*ur0(iarr_all_mx(i),:)*dt
#endif /* ISO */ 
       duls(iarr_all_mx(i),:)  = gravr(:)*ul0(iarr_all_dn(i),:)*dt
       durs(iarr_all_mx(i),:)  = gravl(:)*ur0(iarr_all_dn(i),:)*dt
    enddo
    ur1= ur1 + cn(istep)*durs
    ul1= ul1 + cn(istep)*duls
#endif /* GRAV */

    u1 = ul1 + ur1
    u1(iarr_all_dn,:) = max(u1(iarr_all_dn,:), smalld)

!   write(*,*) 'tv:'
!   write(*,*) dt, cn
!   write(*,*) u1(1,:)
!   write(*,*) u1(2,:)
!   write(*,*) u1(5,:)
!   stop




#ifdef VZ_LIMITS
    if(sweep .eq. 'zsweep') then
      where((z(:) .gt. 0.0) .and. (u1(iarr_all_mx,:) .lt.  floor_vz*u(iarr_all_dn,:)))
        u1(iarr_all_mx,:) =  floor_vz*u(iarr_all_dn,:)
      endwhere
      where((z(:) .lt. 0.0) .and. (u1(iarr_all_mx,:) .gt. -floor_vz*u(iarr_all_dn,:)))
        u1(iarr_all_mx,:) = -floor_vz*u(iarr_all_dn,:)
      endwhere
    endif
#endif /* VZ_LIMITS */

#if defined COSM_RAYS && defined IONIZED
    select case (sweep)
      case('xsweep')
        divv = divvel(:,i1,i2)
      case('ysweep')
        divv = divvel(i2,:,i1)
      case('zsweep')
        divv = divvel(i1,i2,:)
    end select

    decr(:)    = -(gamma_cr-1.)*u1(iecr,:)*divv(:)*dt
    u1(iecr,:) = u1(iecr,:) + cn(istep)*decr(:)
    u1(iecr,:) = max(smallecr,u1(iecr,:))

    vx  = u1(iarr_all_mx(i_ion),:)/u1(iarr_all_dn(i_ion),:)
    ecr = u1(iecr,:)

    grad_pcr(2:n-1) = cr_active*(gamma_cr -1.)*(ecr(3:n)-ecr(1:n-2))/(2.*dx) 
    grad_pcr(1:2)=0.0 ; grad_pcr(n-1:n) = 0.0

#ifndef ISO
    u1(iarr_all_en(i_ion),:) = u1(iarr_all_en(i_ion),:) &
                              - cn(istep)*u1(iarr_all_mx(i_ion),:)/u1(iarr_all_dn(i_ion),:)*grad_pcr*dt
#endif /* ISO */
    u1(iarr_all_mx(i_ion),:) = u1(iarr_all_mx(i_ion),:) - cn(istep)*grad_pcr*dt

#endif /* COSM_RAYS && IONIZED */

#ifndef ISO

    ekin = 0.5*( u1(iarr_all_mx,:)**2 + u1(iarr_all_my,:)**2 &
                +u1(iarr_all_mz,:)**2) /u1(iarr_all_dn,:)
    eint = u1(iarr_all_en,:)-ekin
#if defined IONIZED && defined MAGNETIC
    emag = 0.5*(bb(ibx,:)*bb(ibx,:) + bb(iby,:)*bb(iby,:) + bb(ibz,:)*bb(ibz,:))
    eint(i_ion,:) = eint(i_ion,:) - emag
#endif /* IONIZED && MAGNETIC */

    eint = max(eint,smallei)

    u1(iarr_all_en,:) = eint+ekin
#if defined IONIZED && defined MAGNETIC
    u1(iarr_all_en(i_ion),:) = u1(iarr_all_en(i_ion),:)+emag
#endif /* IONIZED && MAGNETIC */


#endif /* ISO */

      u(:,:) = u1(:,:)
  enddo

!   write(*,*) 'tv:',integration_order
!   write(*,*) u(1,:)
!   write(*,*) u(2,:)
!   write(*,*) u(5,:)
!   stop

    return

  end subroutine relaxing_tvd

!==========================================================================================
end module rtvd
