! $Id: fluxes.F90 433 2008-11-20 22:27:57Z wolt $
#include "piernik.def"
#define RNG 2:n-1
module fluxesneutral
  implicit none

  contains
!==========================================================================================

  subroutine flux_neu(flux,cfr,uu,bb,n)
    use initneutral, only : gamma_neu, cs_iso_neu2
    use constants
    use fluidindex,  only : nvar

    use initneutral, only : idnn,imxn,imyn,imzn
#ifndef ISO
    use initneutral, only : ienn
#endif /* ISO */

    use timestepneut, only : c_neu
    implicit none
    integer n
    real, dimension(nvar,n):: flux,uu,cfr
! locals
    real, dimension(n) :: vx,p  
    
    flux   = 0.0
    vx     = 0.0
    cfr    = 0.0

    vx(RNG)=uu(imxn,RNG)/uu(idnn,RNG)

#ifdef ISO
    p(RNG) = cs_iso_neu2*uu(idnn,RNG)
#else /* ISO */
    p(RNG)=(uu(ienn,RNG)  &
      - 0.5*( uu(imxn,RNG)**2 + uu(imyn,RNG)**2 + uu(imzn,RNG)**2 ) &
          / uu(idnn,RNG))*(gamma_neu-1.0) 
#endif /* ISO */

    flux(idnn,RNG)=uu(imxn,RNG)
    flux(imxn,RNG)=uu(imxn,RNG)*vx(RNG)+p(RNG)
    flux(imyn,RNG)=uu(imyn,RNG)*vx(RNG)
    flux(imzn,RNG)=uu(imzn,RNG)*vx(RNG)
#ifndef ISO
    flux(ienn,RNG)=(uu(ienn,RNG)+p(RNG))*vx(RNG) &
                )/uu(idnn,RNG)
#endif /* ISO */

#ifdef LOCAL_FR_SPEED

!       The freezing speed is now computed locally (in each cell)
!       as in Trac & Pen (2003). This ensures much sharper shocks,
!       but sometimes may lead to numerical instabilities
#ifdef ISO
    cfr(1,RNG) = abs(vx(RNG)) &
               + max(sqrt( abs(p(RNG))/uu(idnn,RNG)),small)
#else /* ISO */
    cfr(1,RNG) = abs(vx(RNG)) &
               + max(sqrt( abs(gamma_neu*p(RNG) &
                )/uu(idnn,RNG)),small)
#endif /* ISO */
    cfr(1,1) = cfr(1,2)
    cfr(1,n) = cfr(1,n-1)
    cfr = spread(cfr(1,:),1,nvar)
#endif /* LOCAL_FR_SPEED */
#ifdef GLOBAL_FR_SPEED
!       The freezing speed is now computed globally
!       (c=const for the whole domain) in sobroutine 'timestep'
    cfr(:,:) = c_neu
!    write(*,*) c_all
!    stop
#endif /* GLOBAL_FR_SPEED */

!   write(*,*) 'flux_neu:', idnn,imxn,imyn,imzn,ienn
!   write(*,*) flux(1,:)
!   write(*,*)
!   write(*,*) flux(2,:)
!   write(*,*)
!   write(*,*) flux(5,:)
!   write(*,*)
   

  end subroutine flux_neu


end module fluxesneutral
