! $Id$
#include "piernik.def"
#define RNG 2:n-1

module fluxdust
  implicit none

  contains
!==========================================================================================

  subroutine flux_dst(fluxd,cfrd,uud,n)
    
    use constants,       only : small
    use fluidindex,      only : idn,imx,imy,imz,ien
    use fluidindex,      only : nvar_dst

    use timestepdust, only : c_dst

    implicit none
    integer n
    
! locals
    real, dimension(nvar_dst,n):: fluxd,uud,cfrd
    real, dimension(n) :: vx,p  
    
    fluxd   = 0.0
    cfrd    = 0.0
    vx      = 0.0

    vx(RNG)=uud(imx,RNG)/uud(idn,RNG)

    fluxd(idn,RNG)=uud(imx,RNG)
    fluxd(imx,RNG)=uud(imx,RNG)*vx(RNG)
    fluxd(imy,RNG)=uud(imy,RNG)*vx(RNG)
    fluxd(imz,RNG)=uud(imz,RNG)*vx(RNG)
    
#ifndef ISO
    fluxd(ien,RNG)=(uud(ien,RNG))*vx(RNG) 
#endif /* ISO */

#ifdef LOCAL_FR_SPEED
!       The freezing speed is now computed locally (in each cell)
!       as in Trac & Pen (2003). This ensures much sharper shocks,
!       but sometimes may lead to numerical instabilities

    cfrd(1,RNG) = max(abs(vx(RNG)),small) 

    cfrd(1,1) = cfrd(1,2)
    cfrd(1,n) = cfrd(1,n-1)
    cfrd = spread(cfrd(1,:),1,nvar_dst)
#endif /* LOCAL_FR_SPEED */

#ifdef GLOBAL_FR_SPEED
!       The freezing speed is now computed globally
!       (c=const for the whole domain) in sobroutine 'timestep'
    cfrd(:,:) = c_dst
#endif /* GLOBAL_FR_SPEED */

  end subroutine flux_dst


end module fluxdust
