! $Id$
#include "piernik.def"
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
    use arrays, only : nu
#ifndef ISO
!    use arrays, only : ienn
#endif /* ISO */
    use time, only : c
    implicit none
    integer n,ifluid
    integer idnd,imxd,imyd,imzd,ienn
    real, dimension(nu,n)::flux,uu,cfr
    real, dimension(3,n):: bb
! locals
    real, dimension(n) :: vx,ps,p,pmag  

    ifluid=1
    idnd=1
    imxd=2
    imyd=3
    imzd=4

    flux   = 0.0
    vx  = 0.0

    cfr = 0.0

    vx(RNG)=uu(imxd,RNG)/uu(idnd,RNG)

    flux(idnd,RNG)=uu(imxd,RNG)
    flux(imxd,RNG)=uu(imxd,RNG)*vx(RNG)
    flux(imyd,RNG)=uu(imyd,RNG)*vx(RNG)
    flux(imzd,RNG)=uu(imzd,RNG)*vx(RNG)

#ifdef LOCAL_FR_SPEED

!       The freezing speed is now computed locally (in each cell)
!       as in Trac & Pen (2003). This ensures much sharper shocks,
!       but sometimes may lead to numerical instabilities

    cfr(1,RNG) = max(abs(vx(RNG)),small)

    cfr(1,1) = cfr(1,2)
    cfr(1,n) = cfr(1,n-1)
    cfr = spread(cfr(1,:),1,nu)
#endif /* LOCAL_FR_SPEED */
#ifdef GLOBAL_FR_SPEED
!       The freezing speed is now computed globally
!       (c=const for the whole domain) in sobroutine 'timestep'
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
