! $Id: tv.F90 438 2008-11-20 23:43:48Z wolt $
#include "piernik.def"
module rtvd ! split orig
  contains

  subroutine tvdb(vibj,b,vg,n,dt,di)
    use constants, only : big
    use func, only      : tvdb_emf
    implicit none
    integer, intent(in) :: n
    real, intent(in)    :: dt,di
    real, dimension(n)  :: vibj,b,vg
! locals
    real, dimension(n)  :: b1,vibj1,vh
    real :: dti

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

    vibj = tvdb_emf(vh,vg,b1,dt)

  end subroutine tvdb

  subroutine relaxing_tvd(u,bb,sweep,i1,i2,dx,n,dt)
  
    use fluidindex, only : i_ion
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
    use arrays, only : xr
    use start,  only : qshear, omega
#endif /* SHEAR */
#ifdef COSM_RAYS
    use start,  only : gamma_cr, cr_active, smallecr
    use arrays, only : wa, iecr
#endif /* COSM_RAYS */
#ifdef VZ_LIMITS
    use arrays, only : z
    use start, only : floor_vz, ceil_vz
#endif /* VZ_LIMITS */
#ifdef KEPLER_SUPPRESSION
    use rotsource, only : kepler_suppression
#endif /* KEPLER_SUPPRESSION */

    implicit none
    
    integer :: i1,i2, n, istep
    
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
    real, dimension(n)    :: vxr
#endif /* SHEAR */
#ifdef KEPLER_SUPPRESSION
    real, dimension(nvar,n) :: Duus
#endif /* KEPLER_SUPPRESSION */
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
    real, dimension(n)    :: divv,decr,gpcr,ecr
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
    vxr(1:n-1) = 0.5*(u1(iarr_all_my,2:n)/u1(iarr_all_dn,2:n) + u1(iarr_all_my,1:n-1)/u1(iarr_all_dn,1:n-1))
    if(sweep .eq. 'xsweep') then
      rotfr(:) =   2.0*omega*(vxr(:) + qshear*omega*xr(:))
    else if(sweep .eq. 'ysweep')  then
      rotfr(:) = - 2.0*omega*vxr(:)
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

#ifndef ISO
    duls(iarr_all_en,:)  = gravr*ul0(iarr_all_mx,:)*dt
    durs(iarr_all_en,:)  = gravl*ur0(iarr_all_mx,:)*dt
#endif /* ISO */
    duls(iarr_all_mx,:)  = gravr*ul0(iarr_all_dn,:)*dt
    durs(iarr_all_mx,:)  = gravl*ur0(iarr_all_dn,:)*dt
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
    u1(iecr,:) = max(smallecr,u1(iecr,:))

    vx = u1(iarr_all_mx,:)/u1(iarr_all_dn,:)
    ecr = u1(iecr,:)

    gpcr(2:n-1) = cr_active*(gamma_cr -1.)*(ecr(3:n)-ecr(1:n-2))/(2.*dx) ; gpcr(1:2)=0.0 ; gpcr(n-1:n) = 0.0

#ifndef ISO
    u1(iarr_all_en,:) = u1(iarr_all_en,:) - cn(istep)*u1(iarr_all_mx,:)/u1(iarr_all_dn,:)*gpcr*dt
#endif /* ISO */
    u1(iarr_all_mx,:) = u1(iarr_all_mx,:) - cn(istep)*gpcr*dt

#endif /* COSM_RAYS */
#ifndef ISO

    ekin = 0.5*(u1(iarr_all_mx,:)*u1(iarr_all_mx,:)+u1(iarr_all_my,:)*u1(iarr_all_my,:)+u1(iarr_all_mz,:)*u1(iarr_all_mz,:))/u1(iarr_all_dn,:)
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

#ifdef VZ_LIMITS
! Dwukrotne ograniczanie vz moze byc zbedne (patrz wyzej)
! Ten fragment bedzie mozna ewentualnie wykasowac po testach
    if(sweep .eq. 'zsweep') then
      where((z(:) .gt. 0.0) .and. (u1(iarr_all_mx,:) .lt.  floor_vz*u(iarr_all_dn,:)))
        u1(iarr_all_mx,:) =  floor_vz*u(iarr_all_dn,:)
      endwhere
      where((z(:) .lt. 0.0) .and. (u1(iarr_all_mx,:) .gt. -floor_vz*u(iarr_all_dn,:)))
        u1(iarr_all_mx,:) = -floor_vz*u(iarr_all_dn,:)
      endwhere
    endif
#endif /* VZ_LIMITS */

#ifdef KEPLER_SUPPRESSION
	call kepler_suppression(Duus,u1,sweep,i1,i2,n,dt)
#ifndef OVERLAP
	u1 = u1 + Duus
#endif /* OVERLAP */
#endif /* KEPLER_SUPPRESSION */

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
