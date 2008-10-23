! $Id$
#include "piernik.def"
module tv      ! unsplit org
  contains
  subroutine tvdb(vibj,b,vg,n,dt)
    use func, only      : tvdb_emf
    use constants, only : big
    use start, only     : istep
    implicit none
    integer, intent(in) :: n
    real, intent(in)    :: dt
    real, dimension(n)  :: vibj,b,vg
! locals
    real, dimension(n) :: vh

  ! unlike the B field, the vibj lives on the right cell boundary
    vh = 0.0
    vh(1:n-1) =(vg(1:n-1)+ vg(2:n))*0.5;     vh(n) = vh(n-1)
    if(istep .eq. 1) then
      where(vh > 0.)
        vibj=b*vg
      elsewhere
        vibj=eoshift(b*vg,1,boundary=big)
      end where
      vibj = vibj*dt

    else if(istep .eq. 2) then
      vibj = tvdb_emf(vh,vg,b,dt)
    else 
      stop
    endif

  end subroutine tvdb


  subroutine relaxing_tvd(ui,u,bi,bb,sweep,i1,i2,dx,n,dt)
 ! Cooling and heating implemented following Rafal Kosinski
   use arrays, only  : ibx,iby,ibz,nu,idna,imxa,imya,imza
   use start, only   : istep,smallei,smalld
#ifndef ISO
   use arrays, only  : iena
#endif
#ifdef GRAV
   use gravity, only : grav_pot2accel
#endif
#ifdef KEPLER_SUPPRESSION
    use rotsource, only : kepler_suppression
#endif /* KEPLER_SUPPRESSION */

    implicit none
    integer i1,i2, n
    real :: dt,dx,dtx
    real, dimension(nu,n) :: u,ui,cfr,ul,ur,ul0,ur0
    real, dimension(3,n)  :: bb,bi
    real, dimension(n)    :: gravl,gravr
    real, dimension(n)    :: dgrp,dgrm,dglp,dglm
#ifdef KEPLER_SUPPRESSION
    real, dimension(nu,n) :: Duus
#endif /* KEPLER_SUPPRESSION */
           
    character sweep*6
          
!locals    
    real, dimension(nu,n) :: w,fr,fl,flux,dfrp,dfrm,dflm,dflp, &
                             dulf,durf
#ifdef GRAV
    real, dimension(nu,n) :: durs,duls
#endif /* GRAV */
#ifdef COSM_RAYS
    real q
    real, dimension(n)    :: divv,decr,tmp
#endif /* COSM_RAYS */
                            
    w(:,:)    = 0.0
    cfr  = 0.0
    dtx  = dt / dx
#ifdef GRAV
    duls = 0.0 
    durs = 0.0
#endif /* GRAV */
    
    call mhdflux(w,cfr,ui,bi,n)

    fr = (ui*cfr+w)*0.5
    fl = (ui*cfr-w)*0.5
   ur0 = fr/cfr
   ul0 = fl/cfr
      
    fl(:,1:n-1) = fl(:,2:n);                       fl(:,n) = fl(:,1)        ! fl=cshift(fl,shift=1,dim=2)
    if(istep .eq. 2) then
      call mhdflux(w,cfr,u,bb,n)

      fr = (u*cfr+w)*0.5
      fl = (u*cfr-w)*0.5
      ur = fr/cfr
      ul = fl/cfr
    
      fl(:,1:n-1) = fl(:,2:n);                     fl(:,n) = fl(:,1)        ! fl=cshift(fl,shift=1,dim=2)

      dfrp(:,1:n-1) = 0.5*(fr(:,2:n) - fr(:,1:n-1)); dfrp(:,n) = dfrp(:,1)  ! dfrp=(cshift(fr,shift=1,dim=2)-fr)*0.5    
      dfrm(:,2:n)   = dfrp(:,1:n-1);                 dfrm(:,1) = dfrm(:,n)  ! dfrm=cshift(dfrp,shift=-1,dim=2)          

      call flimiter(fr,dfrm,dfrp,nu,n)
    
      dflp(:,1:n-1) = 0.5*(fl(:,1:n-1) - fl(:,2:n)); dflp(:,n) = dflp(:,1)  ! dflp=(fl-cshift(fl,shift=1,dim=2))*0.5    
      dflm(:,2:n)   = dflp(:,1:n-1);                 dflm(:,1) = dflm(:,n)  ! dflm=cshift(dflp,shift=-1,dim=2)
      call flimiter(fl,dflm,dflp,nu,n)
    endif

    durf(:,2:n) = dtx*(fr(:,2:n) - fr(:,1:n-1));     durf(:,1) = durf(:,n)  ! durf = (fr-cshift(fr,shift=-1,dim=2))/dx*dt 
    dulf(:,2:n) = dtx*(fl(:,2:n) - fl(:,1:n-1));     dulf(:,1) = dulf(:,n)  ! dulf = (fl-cshift(fl,shift=-1,dim=2))/dx*dt

    
#ifdef GRAV

! Gravity source terms -------------------------------------

      call grav_pot2accel(sweep,i1,i2, n, gravr)
      gravl(2:n) = gravr(1:n-1)
      gravl(1)  = gravl(2)
      gravr(n)  = gravr(n-1)
         
      if(istep .eq. 2) then
        dgrp(1:n-1) = 0.5*(gravr(1:n-1) - gravr(2:n));  dgrp(n) = dgrp(1)   ! dgrp = (gravr-cshift(gravr,shift=1))*0.5
        dgrm(2:n) = dgrp(1:n-1);                        dgrm(1) = dgrm(n)   ! dgrm = (cshift(gravr,shift=-1)-gravr)*0.5

        call flimiter(gravr,dgrp,dgrm,1,n)

        dglp(1:n-1) = 0.5*(gravl(2:n) - gravl(1:n-1));  dglp(n) = dglp(1)   ! dglp = (cshift(gravl,shift=1)-gravl)*0.5
        dglm(2:n)   = dglp(1:n-1);                      dglm(1) = dglm(n)   ! dglm = (gravl-cshift(gravl,shift=-1))*0.5
        call flimiter(gravl,dglp,dglm,1,n)
      endif

#ifndef /* ISO */
      duls(iena,:)  = gravr*ul0(imxa,:)*dt 
      durs(iena,:)  = gravl*ur0(imxa,:)*dt
#endif /* ISO */
      duls(imxa,:)  = gravr*ul0(idna,:)*dt
      durs(imxa,:)  = gravl*ur0(idna,:)*dt

! ----------------------------------------------------------
#endif /* GRAV */

   ur = -durf
   ul =  dulf

#ifdef GRAV
   ur = ur + durs
   ul = ul + duls
#endif /* GRAV */

   u = ur + ul
	
#ifdef KEPLER_SUPPRESSION
	call kepler_suppression(Duus,u,sweep,i1,i2,n,dt)
#ifndef OVERLAP
	u = u + Duus
#endif /* OVERLAP */
#endif /* KEPLER_SUPPRESSION */

    return

  end subroutine relaxing_tvd

!==========================================================================================
 
  subroutine initials
    use arrays, only : ui,bi,u,b
    implicit none

    ui(:,:,:,:) = u(:,:,:,:)
    bi(:,:,:,:) = b(:,:,:,:)

  end subroutine initials

  subroutine integrate

    use fluid_boundaries, only : compute_u_bnd
    use mag_boundaries, only :   compute_b_bnd
    use arrays, only : ui,u,bi,b,Lu,Lb,nx,ny,nz,idna,imxa,imya,imza,&
         ibx,iby,ibz
    use start, only  : istep,smalld
#ifndef ISO
    use arrays, only : iena
    use start, only  : smallei
#endif /* ~ISO */
    implicit none
#ifndef ISO
    real, allocatable :: ekin(:,:,:), emag(:,:,:), eint(:,:,:)

    allocate(ekin(nx,ny,nz),emag(nx,ny,nz),eint(nx,ny,nz))
#endif /* ISO */

    select case(istep)
     case(1)
        u(:,:,:,:) = ui(:,:,:,:) + 0.5*Lu(:,:,:,:)
     case(2)
        u(:,:,:,:) = ui(:,:,:,:) + Lu(:,:,:,:)
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
#endif /* ~ISO */
    call compute_u_bnd

    select case(istep)
     case(1)
        b(:,:,:,:) = bi(:,:,:,:) + 0.5*Lb(:,:,:,:)
     case(2)
        b(:,:,:,:) = b(:,:,:,:) + Lb(:,:,:,:)
     case default
      write(*,*) 'That high integration order not implemented'
    end select
    call compute_b_bnd

#ifndef ISO
  deallocate(ekin,emag,eint)
#endif /* ~ISO */

  end subroutine integrate  

end module tv
