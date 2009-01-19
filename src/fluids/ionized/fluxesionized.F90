! $Id: fluxes.F90 433 2008-11-20 22:27:57Z wolt $
#include "piernik.def"
#define RNG 2:n-1
module fluxesionized
  implicit none

  contains 
!==========================================================================================

  subroutine flux_ion(flux,cfr,uu,bb,n)
    use initionized, only : gamma_ion, cs_iso_ion2
    use constants
    use fluidindex, only  : nvar   
    use fluidindex, only  : nmag
    use initionized, only : ibx,iby,ibz
    use initionized, only : idni,imxi,imyi,imzi
#ifndef ISO
    use initionized, only : ieni
#endif /* ISO */
!#ifdef COSM_RAYS
!    use arrays, only : iecr
!#endif /* COSM_RAYS */
    use timestepion, only : c_ion
    implicit none
    integer n
    real, dimension(nvar,n):: flux,uu,cfr
    real, dimension(nmag,n):: bb
! locals
    real, dimension(n) :: vx,ps,p,pmag  
    
    flux   = 0.0
    vx     = 0.0
    cfr    = 0.0

    pmag(RNG)=0.5*( bb(ibx,RNG)**2 + bb(iby,RNG)**2 +bb(ibz,RNG)**2 )
    vx(RNG)=uu(imxi,RNG)/uu(idni,RNG)

#ifdef ISO
    p(RNG) = cs_iso_ion2*uu(idni,RNG)
    ps(RNG)= p(RNG) + pmag(RNG)
#else /* ISO */
    ps(RNG)=(uu(ieni,RNG) - &
      0.5*( uu(imxi,RNG)**2 + uu(imyi,RNG)**2 + uu(imzi,RNG)**2 ) &
          / uu(idni,RNG))*(gamma_ion-1.0) + (2.0-gamma_ion)*pmag(RNG)
    p(RNG) = ps(RNG)- pmag(RNG)
#endif /* ISO */

    flux(idni,RNG)=uu(imxi,RNG)
    flux(imxi,RNG)=uu(imxi,RNG)*vx(RNG)+ps(RNG) - bb(ibx,RNG)**2
    flux(imyi,RNG)=uu(imyi,RNG)*vx(RNG)-bb(iby,RNG)*bb(ibx,RNG)
    flux(imzi,RNG)=uu(imzi,RNG)*vx(RNG)-bb(ibz,RNG)*bb(ibx,RNG)
#ifndef ISO
    flux(ieni,RNG)=(uu(ieni,RNG)+ps(RNG))*vx(RNG)-bb(ibx,RNG)*(bb(ibx,RNG)*uu(imxi,RNG) &
                +bb(iby,RNG)*uu(imyi,RNG)+bb(ibz,RNG)*uu(imzi,RNG))/uu(idni,RNG)
#endif /* ISO */
!#ifdef COSM_RAYS
!    flux(iecr,RNG)= uu(iecr,RNG)*vx(RNG)
!#endif /* COSM_RAYS */
#ifdef LOCAL_FR_SPEED

!       The freezing speed is now computed locally (in each cell)
!       as in Trac & Pen (2003). This ensures much sharper shocks,
!       but sometimes may lead to numerical instabilities
#ifdef ISO
    cfr(1,RNG) = abs(vx(RNG)) &
               + max(sqrt( abs(2.0*pmag(RNG) + p(RNG))/uu(idni,RNG)),small)
#else /* ISO */
    cfr(1,RNG) = abs(vx(RNG)) &
               + max(sqrt( abs(2.0*pmag(RNG) + gamma_ion*p(RNG) &
                )/uu(idni,RNG)),small)
#endif /* ISO */
    cfr(1,1) = cfr(1,2)
    cfr(1,n) = cfr(1,n-1)
    cfr = spread(cfr(1,:),1,nvar)
#endif /* LOCAL_FR_SPEED */
#ifdef GLOBAL_FR_SPEED
!       The freezing speed is now computed globally
!       (c=const for the whole domain) in sobroutine 'timestep'
    cfr(:,:) = c_ion
!    write(*,*) c_all
!    stop
#endif /* GLOBAL_FR_SPEED */

!   write(*,*) 'flux_ion:', idni,imxi,imyi,imzi,ieni
!   write(*,*) flux(1,:)
!   write(*,*)
!   write(*,*) flux(2,:)
!   write(*,*)
!   write(*,*) flux(5,:)
!   write(*,*)
   

  end subroutine flux_ion


end module fluxesionized
