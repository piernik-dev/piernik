! $Id: fluxes.F90 433 2008-11-20 22:27:57Z wolt $
#include "piernik.def"
#define RNG 2:n-1
module fluxneutral
  implicit none

  contains
!==========================================================================================

  subroutine flux_neu(fluxn,cfrn,uun,n)
    
    use constants,       only : small
    use fluidindex,      only : idn,imx,imy,imz,ien
    use fluidindex,      only : nvar_neu

    use initneutral,     only : gamma_neu, cs_iso_neu2
    use timestepneutral, only : c_neu

    implicit none
    integer n
    
! locals
    real, dimension(nvar_neu,n):: fluxn,uun,cfrn
    real, dimension(n) :: vx,p  
    
    fluxn   = 0.0
    cfrn    = 0.0
    vx      = 0.0

    vx(RNG)=uun(imx,RNG)/uun(idn,RNG)

#ifdef ISO
    p(RNG) = cs_iso_neu2*uun(idn,RNG)
#else /* ISO */
    p(RNG)=(uun(ien,RNG)  &
      - 0.5*( uun(imx,RNG)**2 + uun(imy,RNG)**2 + uun(imz,RNG)**2 ) &
          / uun(idn,RNG))*(gamma_neu-1.0) 
#endif /* ISO */

    fluxn(idn,RNG)=uun(imx,RNG)
    fluxn(imx,RNG)=uun(imx,RNG)*vx(RNG)+p(RNG)
    fluxn(imy,RNG)=uun(imy,RNG)*vx(RNG)
    fluxn(imz,RNG)=uun(imz,RNG)*vx(RNG)
#ifndef ISO
    fluxn(ien,RNG)=(uun(ien,RNG)+p(RNG))*vx(RNG) 
#endif /* ISO */

#ifdef LOCAL_FR_SPEED

!       The freezing speed is now computed locally (in each cell)
!       as in Trac & Pen (2003). This ensures much sharper shocks,
!       but sometimes may lead to numerical instabilities
#ifdef ISO
    cfrn(1,RNG) = abs(vx(RNG)) &
               + max(sqrt( abs(p(RNG))/uun(idn,RNG)),small)
#else /* ISO */
    cfrn(1,RNG) = abs(vx(RNG)) &
               + max(sqrt( abs(gamma_neu*p(RNG) &
                )/uun(idn,RNG)),small)
#endif /* ISO */
    cfrn(1,1) = cfrn(1,2)
    cfrn(1,n) = cfrn(1,n-1)
    cfrn = spread(cfrn(1,:),1,nvar_neu)
#endif /* LOCAL_FR_SPEED */

#ifdef GLOBAL_FR_SPEED
!       The freezing speed is now computed globally
!       (c=const for the whole domain) in sobroutine 'timestep'
    cfrn(:,:) = c_neu
#endif /* GLOBAL_FR_SPEED */

  end subroutine flux_neu


end module fluxneutral
