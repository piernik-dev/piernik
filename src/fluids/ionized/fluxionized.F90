! $Id: fluxes.F90 433 2008-11-20 22:27:57Z wolt $
#include "piernik.def"
#define RNG 2:n-1


module fluxionized
  implicit none

  contains 
!==========================================================================================

  subroutine flux_ion(fluxi,cfri,uui,bb,n)
  
    use constants,       only : small
    use fluidindex,      only : nmag
    use fluidindex,      only : ibx,iby,ibz
    use fluidindex,      only : idn,imx,imy,imz,ien
    use fluidindex,      only : nvar_ion   
  
    use initionized,     only : gamma_ion, cs_iso_ion2
    use timestepionized, only : c_ion
    use constants,       only : small

!#ifdef COSM_RAYS
!    use arrays, only : iecr
!#endif /* COSM_RAYS */

    implicit none
    integer n

! locals
    real, dimension(nvar_ion,n):: fluxi,uui,cfri
    real, dimension(nmag,n):: bb
    real, dimension(n) :: vx,ps,p,pmag  
    
    fluxi   = 0.0
    cfri    = 0.0
    vx      = 0.0

    pmag(RNG)=0.5*( bb(ibx,RNG)**2 + bb(iby,RNG)**2 +bb(ibz,RNG)**2 )
    vx(RNG)=uui(imx,RNG)/uui(idn,RNG)

#ifdef ISO
    p(RNG) = cs_iso_ion2*uui(idn,RNG)
    ps(RNG)= p(RNG) + pmag(RNG)
#else /* ISO */
    ps(RNG)=(uui(ien,RNG) - &
      0.5*( uui(imx,RNG)**2 + uui(imy,RNG)**2 + uui(imz,RNG)**2 ) &
          / uui(idn,RNG))*(gamma_ion-1.0) + (2.0-gamma_ion)*pmag(RNG)
    p(RNG) = ps(RNG)- pmag(RNG)
#endif /* ISO */

    fluxi(idn,RNG)=uui(imx,RNG)
    fluxi(imx,RNG)=uui(imx,RNG)*vx(RNG)+ps(RNG) - bb(ibx,RNG)**2
    fluxi(imy,RNG)=uui(imy,RNG)*vx(RNG)-bb(iby,RNG)*bb(ibx,RNG)
    fluxi(imz,RNG)=uui(imz,RNG)*vx(RNG)-bb(ibz,RNG)*bb(ibx,RNG)
#ifndef ISO
    fluxi(ien,RNG)=(uui(ien,RNG)+ps(RNG))*vx(RNG)-bb(ibx,RNG)*(bb(ibx,RNG)*uui(imx,RNG) &
                +bb(iby,RNG)*uui(imy,RNG)+bb(ibz,RNG)*uui(imz,RNG))/uui(idn,RNG)
#endif /* ISO */

!#ifdef COSM_RAYS
!    fluxi(iecr,RNG)= uui(iecr,RNG)*vx(RNG)
!#endif /* COSM_RAYS */

#ifdef LOCAL_FR_SPEED

!       The freezing speed is now computed locally (in each cell)
!       as in Trac & Pen (2003). This ensures much sharper shocks,
!       but sometimes may lead to numerical instabilities
#ifdef ISO
    cfri(1,RNG) = abs(vx(RNG)) &
               + max(sqrt( abs(2.0*pmag(RNG) + p(RNG))/uui(idn,RNG)),small)
#else /* ISO */
    cfri(1,RNG) = abs(vx(RNG)) &
               + max(sqrt( abs(2.0*pmag(RNG) + gamma_ion*p(RNG) &
                )/uui(idn,RNG)),small)
#endif /* ISO */
    cfri(1,1) = cfri(1,2)
    cfri(1,n) = cfri(1,n-1)
    cfri = spread(cfri(1,:),1,nvar_ion)
#endif /* LOCAL_FR_SPEED */

#ifdef GLOBAL_FR_SPEED
!       The freezing speed is now computed globally
!       (c=const for the whole domain) in sobroutine 'timestep'
    cfri(:,:) = c_ion
#endif /* GLOBAL_FR_SPEED */

  end subroutine flux_ion


end module fluxionized
