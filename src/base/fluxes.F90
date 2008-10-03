! $Id$
#include "piernik.def"
#define RNG 2:n-1
module fluxes
  implicit none

  contains
!==========================================================================================
 
  subroutine mhdflux(flux,cfr,uu,bb,n)
    use start, only : gamma, gamma_cr
#ifdef ISO
    use start, only : csi2
#endif /* ISO */
    use constants
    use arrays, only : ibx,iby,ibz,idna,imxa,imya,imza, nu
#ifndef ISO
    use arrays, only : iena
#endif /* ISO */
#ifdef COSM_RAYS
    use arrays, only :iecr
#endif /* COSM_RAYS */
    use time, only : c
    implicit none
    integer n
    real, dimension(nu,n)::flux,uu,cfr
    real, dimension(3,n):: bb
! locals
    real, dimension(n) :: vx,vy,vz,vt, ps,p,pmag
!#ifdef COSM_RAYS
!    real, dimension(n) :: pcr
!#endif /* COSM_RAYS */

    flux   = 0.0
    vx  = 0.0
    vy  = 0.0
    vz  = 0.0
    
    cfr = 0.0 

    pmag(RNG)=0.5*( bb(ibx,RNG)**2 + bb(iby,RNG)**2 +bb(ibz,RNG)**2 )
    vx(RNG)=uu(imxa,RNG)/uu(idna,RNG)

! UWAGA ZMIANA W OBLICZANIU LOKALNEJ PREDKOSCI NA PROBE
! mh 22-11-07 w problemach z udzialem promieniowania kosmicznego
! uzycie vx prowadzi do nieduzych oscylacji, a dla vt oscylacji nie widac,
! wobec tego po krotkim powrocie do vx (v.1.2) przywracam vt (v.1.3)

#ifndef OLDLOCALCFR
    vy(RNG)=uu(imya,RNG)/uu(idna,RNG)
    vz(RNG)=uu(imza,RNG)/uu(idna,RNG)
    vt = sqrt(vx**2 + vy**2 + vz**2)
#endif /* OLDLOCALCFR */    
    
#ifdef ISO
    p(RNG) = csi2*uu(idna,RNG)
    ps(RNG)= p(RNG) + pmag(RNG)
#else /* ISO */
    ps(RNG)=(uu(iena,RNG) - &
      0.5*( uu(imxa,RNG)**2 + uu(imya,RNG)**2 + uu(imza,RNG)**2 ) &
          / uu(idna,RNG))*(gamma-1.0) + (2.0-gamma)*pmag(RNG)
    p(RNG) = ps(RNG)- pmag(RNG)
#endif /* ISO */
!#ifdef COSM_RAYS
!    pcr = (gamma_cr-1)*uu(iecr,RNG)
!#endif COSM_RAYS
    flux(idna,RNG)=uu(imxa,RNG)
    flux(imxa,RNG)=uu(imxa,RNG)*vx(RNG)+ps(RNG) - bb(ibx,RNG)**2
    flux(imya,RNG)=uu(imya,RNG)*vx(RNG)-bb(iby,RNG)*bb(ibx,RNG)
    flux(imza,RNG)=uu(imza,RNG)*vx(RNG)-bb(ibz,RNG)*bb(ibx,RNG)
#ifndef ISO
    flux(iena,RNG)=(uu(iena,RNG)+ps(RNG))*vx(RNG)-bb(ibx,RNG)*(bb(ibx,RNG)*uu(imxa,RNG) &
                +bb(iby,RNG)*uu(imya,RNG)+bb(ibz,RNG)*uu(imza,RNG))/uu(idna,RNG)
#endif /* ISO */
#ifdef COSM_RAYS
    flux(iecr,RNG)= uu(iecr,RNG)*vx(RNG)
#endif /* COSM_RAYS */
#ifdef LOCAL_FR_SPEED

!       The freezing speed is now computed locally (in each cell) 
!       as in Trac & Pen (2003). This ensures much sharper shocks, 
!       but sometimes may lead to numerical instabilities   
#ifdef ISO
! UWAGA ZMIANA W OBLICZANIU LOKALNEJ PREDKOSCI NA PROBE
#ifdef OLDLOCALCFR
        cfr(1,RNG) = abs(vx(RNG)) &
#else /* OLDLOCALCFR */
        cfr(1,RNG) = abs(vt(RNG)) &
#endif /* OLDLOCALCFR */
                      +max(sqrt( abs(2.0*pmag(RNG) + p(RNG))/uu(idna,RNG)),small)
#else /* ISO */
! UWAGA ZMIANA W OBLICZANIU LOKALNEJ PREDKOSCI NA PROBE
#ifdef OLDLOCALCFR
        cfr(1,RNG) = abs(vx(RNG)) &
#else /* OLDLOCALCFR */
        cfr(1,RNG) = abs(vt(RNG)) &
#endif /* OLDLOCALCFR */
                      +max(sqrt( abs(2.0*pmag(RNG) + gamma*p(RNG) &
!#ifdef COSM_RAYS
!                                                   + gamma_cr*pcr(RNG) &
!#endif COSM_RAYS		         
                )/uu(idna,RNG)),small)
#endif /* ISO */
        cfr(1,1) = cfr(1,2)
        cfr(1,n) = cfr(1,n-1)   
        cfr = spread(cfr(1,:),1,nu)
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

end module fluxes
