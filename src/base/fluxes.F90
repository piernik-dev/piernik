! $Id$
#include "piernik.def"
#define SQR(var) ((var)*(var))
#define SUM_SQR(x,y,z) ( (x)**2+(y)**2+(z)**2 )
#define RNG 2:n-1
module fluxes
  implicit none

  contains
!==========================================================================================

  subroutine mhdflux(flux,cfr,uu,bb,n)
    use start, only : gamma !, gamma_cr
#ifdef ISO
    use start, only : csi2
#endif /* ISO */
    use constants
    use arrays, only : ibx,iby,ibz,idna,imxa,imya,imza, nu,magn,nfluid,nfmagn,fmagn,nadiab
#ifndef ISO
    use arrays, only : iena, fadiab
#endif /* ISO */
#ifdef COSM_RAYS
    use cosmic_rays, only : iecr, gamma_cr
#endif /* COSM_RAYS */
#ifdef DUST
    use arrays, only : fdust
#endif /* DUST */
    use allfluxes, only : all_flux
    use time, only : c
    implicit none
    integer :: n,i, indx
    real, dimension(nu,n)::flux,uu,cfr
    real, dimension(3,n):: bb
    real, dimension(nfmagn,3,n) :: bbn
! locals
    real, dimension(nfmagn,n) :: pmag
    real, dimension(nfluid,n) :: vx,vy,vz,vt,ps,p
!#ifdef COSM_RAYS
!    real, dimension(n) :: pcr
!#endif /* COSM_RAYS */

    flux   = 0.0
    vx  = 0.0
    vy  = 0.0
    vz  = 0.0

    cfr = 0.0

#ifdef IONIZED
    bbn    = spread(bb,1,nfmagn)
    pmag(:,RNG) = 0.5*SUM_SQR(bbn(:,ibx,RNG),bbn(:,iby,RNG),bbn(:,ibz,RNG))
#endif /* IONIZED */
    vx(:,RNG)=uu(imxa,RNG)/uu(idna,RNG)
#ifdef OLDLOCALCFR
    vt(:,RNG)=vx(:,RNG)
#else /* OLDLOCALCFR */
    vy(:,RNG)=uu(imya,RNG)/uu(idna,RNG)
    vz(:,RNG)=uu(imza,RNG)/uu(idna,RNG)
    vt = sqrt(SUM_SQR(vx,vy,vz))
#endif /* OLDLOCALCFR */

#ifdef ISO
    p(:,RNG) = csi2*uu(idna,RNG)
    ps(:,RNG)= p(:,RNG)
#ifdef IONIZED
    ps(fmagn,RNG)= p(fmagn,RNG) + pmag(fmagn,RNG)
#endif /* IONIZED */
#else /* ISO */
    ps(fadiab,RNG)=(uu(iena(fadiab),RNG) - &
       0.5*( uu(imxa(fadiab),RNG)**2 + &
             uu(imya(fadiab),RNG)**2 + &
             uu(imza(fadiab),RNG)**2) &
       / uu(idna(fadiab),RNG))
    do i=1,nadiab
      indx = fadiab(i)
      ps(indx,:) = ps(indx,:) * (gamma(indx)-1.0d0)
    enddo
    p(fadiab,RNG) = ps(fadiab,RNG)
#ifdef IONIZED
    do i=1,nfmagn
      indx = fmagn(i)
      ps(indx,RNG)= ps(indx,RNG) + (2.0d0 - gamma(indx))*pmag(indx,RNG)
    enddo
    p(fmagn,RNG) = ps(fmagn,RNG)- pmag(fmagn,RNG)
#endif /* IONIZED */
#endif /* ISO */
#ifdef DUST
    ps(fdust,RNG)=0.0
    p(fdust,RNG) =0.0
#endif /* DUST */
    flux(idna,RNG)=uu(imxa,RNG)
    flux(imxa,RNG)=uu(imxa,RNG)*vx(:,RNG)+ps(:,RNG)
    flux(imya,RNG)=uu(imya,RNG)*vx(:,RNG)
    flux(imza,RNG)=uu(imza,RNG)*vx(:,RNG)
#ifdef IONIZED
    flux(imxa(fmagn),RNG)=flux(imxa(fmagn),RNG)- SQR(bbn(fmagn,ibx,RNG))
    flux(imya(fmagn),RNG)=flux(imya(fmagn),RNG) - bbn(fmagn,iby,RNG)*bbn(fmagn,ibx,RNG)
    flux(imza(fmagn),RNG)=flux(imza(fmagn),RNG) - bbn(fmagn,ibz,RNG)*bbn(fmagn,ibx,RNG)
#endif /* IONIZED */
#ifndef ISO
    flux(iena,RNG)=(uu(iena,RNG)+ps(:,RNG))*vx(:,RNG)
#ifdef IONIZED
    flux(iena(fmagn),RNG)=flux(iena(fmagn),RNG) -bbn(fmagn,ibx,RNG) &
                           *(uu(imxa(fmagn),RNG)*bbn(fmagn,ibx,RNG) &
                            +uu(imya(fmagn),RNG)*bbn(fmagn,iby,RNG) &
                            +uu(imza(fmagn),RNG)*bbn(fmagn,ibz,RNG))/uu(idna(fmagn),RNG)
#endif /* IONIZED */
#endif /* ISO */
#ifdef LOCAL_FR_SPEED

!       The freezing speed is now computed locally (in each cell)
!       as in Trac & Pen (2003). This ensures much sharper shocks,
!       but sometimes may lead to numerical instabilities
#ifndef ISO
    do i=1,nadiab
      indx = fadiab(i)
      p(indx,RNG) = gamma(indx)*p(fadiab(indx),RNG) ! UWAGA : przy ewentualnym dalszym korzystaniu z p
    enddo
#endif /* ISO */
    p(fmagn,RNG)= p(fmagn,RNG)+ 2.0 * pmag(fmagn,RNG)  ! UWAGA : przy ewentualnym dalszym korzystaniu z p
    cfr(idna,RNG) = abs(vt(:,RNG))+max(sqrt(abs(p(:,RNG))/uu(idna,RNG)),small)
    cfr(idna,1) = cfr(idna,2)
    cfr(idna,n) = cfr(idna,n-1)
    cfr(imxa,:) = cfr(idna,:)
    cfr(imya,:) = cfr(idna,:)
    cfr(imza,:) = cfr(idna,:)
#ifndef ISO
    cfr(iena(fadiab),:)=cfr(idna(fadiab),:)
#endif /* ISO */
    call all_flux(flux,cfr,uu,bb,n)
#ifdef COMMONCFR
    cfr(1,:) = maxval(cfr(:,:),DIM=1)
    cfr(:,:) = spread(cfr(1,:),DIM=1,NCOPIES=nu)
#endif /* COMMONCFR */
#endif /* LOCAL_FR_SPEED */
#ifdef GLOBAL_FR_SPEED
!       The freezing speed is now computed globally
!       (c=const for the whole domain) in sobroutine 'timestep'
!       Original computation of the freezing speed was done
!       for each sweep separately:
!       c=maxval(abs(vx(nb-RNG))+sqrt(abs(2*pmag(nb-RNG) &
!                            +gamma*p(nb-RNG))/u(idna,nb-RNG)))
    cfr(:,:) = c
#endif /* GLOBAL_FR_SPEED */

  end subroutine mhdflux

!==========================================================================================


  subroutine flimiter(f,a,b,m,n)
    implicit none
    integer m,n
    real, dimension(m,n) :: f,a,b,c
#ifdef VANLEER
      c = a*b
      where (c .gt. 0.0)
        f = f+2.0*c/(a+b)
      endwhere
#endif /* VANLEER */
#ifdef MONCEN
        f = f+(sign(1.0,a)+sign(1.0,b))*min(2.*abs(a),2.*abs(b),0.5*abs(a+b))*0.5
#endif /* MONCEN */
#ifdef MINMOD
        f = f+(sign(1.0,a)+sign(1.0,b))*min(abs(a),abs(b))*0.5
#endif /* MINMOD */
#ifdef SUPERBEE
      where (abs(a) .gt. abs(b))
        f = f+(sign(1.0,a)+sign(1.0,b))*min(abs(a), abs(2.0*b))*0.5
      elsewhere
        f = f+(sign(1.0,a)+sign(1.0,b))*min(abs(2.0*a), abs(b))*0.5
      endwhere
#endif /* SUPERBEE */

    return
  end subroutine flimiter

  subroutine flux_limit(fr,fl,m,n)
    use start, only : istep

    implicit none
    integer             :: m,n
    real,dimension(m,n) :: fr,fl,dfrp,dfrm,dflp,dflm

#ifdef ORIG
    if(istep == 2) then
      dfrp(:,1:n-1) = 0.5*(fr(:,2:n) - fr(:,1:n-1)); dfrp(:,n) = dfrp(:,n-1)
      dfrm(:,2:n)   = dfrp(:,1:n-1);                 dfrm(:,1) = dfrm(:,2)
      call flimiter(fr,dfrm,dfrp,m,n)

      dflp(:,1:n-1) = 0.5*(fl(:,1:n-1) - fl(:,2:n)); dflp(:,n) = dflp(:,n-1)
      dflm(:,2:n)   = dflp(:,1:n-1);                 dflm(:,1) = dflm(:,2)
      call flimiter(fl,dflm,dflp,m,n)
    endif
#endif /* ORIG */
#ifdef SSP
      dfrp(:,1:n-1) = 0.5*(fr(:,2:n) - fr(:,1:n-1)); dfrp(:,n) = dfrp(:,n-1)
      dfrm(:,2:n)   = dfrp(:,1:n-1);                 dfrm(:,1) = dfrm(:,2)
      call flimiter(fr,dfrm,dfrp,m,n)

      dflp(:,1:n-1) = 0.5*(fl(:,1:n-1) - fl(:,2:n)); dflp(:,n) = dflp(:,n-1)
      dflm(:,2:n)   = dflp(:,1:n-1);                 dflm(:,1) = dflm(:,2)
      call flimiter(fl,dflm,dflp,m,n)
#endif /* SSP */

  end subroutine flux_limit

  subroutine grav_limit(gravr,gravl,n)
   use start, only : istep

   implicit none
   integer             :: n
   real,dimension(n)   :: gravr,gravl,dgrp,dgrm,dglp,dglm

#ifdef ORIG
    if(istep == 2) then
      dgrp(1:n-1) = 0.5*(gravr(1:n-1) - gravr(2:n));  dgrp(n) = dgrp(n-1)
      dgrm(2:n) = dgrp(1:n-1)                      ;  dgrm(1) = dgrm(2)
      call flimiter(gravr,dgrm,dgrp,1,n)

      dglp(1:n-1) = 0.5*(gravl(2:n) - gravl(1:n-1));  dglp(n) = dglp(n-1)
      dglm(2:n)   = dglp(1:n-1)                    ;  dglm(1) = dglm(2)
      call flimiter(gravl,dglm,dglp,1,n)
    endif
#endif /* ORIG */
#ifdef SSP
      dgrp(1:n-1) = 0.5*(gravr(1:n-1) - gravr(2:n));  dgrp(n) = dgrp(n-1)
      dgrm(2:n) = dgrp(1:n-1)                      ;  dgrm(1) = dgrm(2)
      call flimiter(gravr,dgrm,dgrp,1,n)

      dglp(1:n-1) = 0.5*(gravl(2:n) - gravl(1:n-1));  dglp(n) = dglp(n-1)
      dglm(2:n)   = dglp(1:n-1)                    ;  dglm(1) = dglm(2)
      call flimiter(gravl,dglm,dglp,1,n)
#endif /* SSP */

  end subroutine grav_limit

end module fluxes
