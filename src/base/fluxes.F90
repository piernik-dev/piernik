#include "mhd.def"

module fluxes
  implicit none

  contains
!==========================================================================================
 
  subroutine mhdflux(flux,cfr,u,b,n)
    use constants
    use arrays
    use time, only : c
    implicit none
    integer n
    real, dimension(nu,n)::flux,u,cfr
    real, dimension(3,n):: b
! locals
    real, dimension(n) :: vx,vy,vz,vt, ps,p,pmag

    flux   = 0.0
    vx  = 0.0
    vy  = 0.0
    vz  = 0.0
    
    cfr = 0.0 

    pmag(2:n-1)=(b(ibx,2:n-1)*b(ibx,2:n-1)+b(iby,2:n-1)*b(iby,2:n-1)+b(ibz,2:n-1)*b(ibz,2:n-1))*0.5
    vx(2:n-1)=u(imxa,2:n-1)/u(idna,2:n-1)

! UWAGA ZMIANA W OBLICZANIU LOKALNEJ PREDKOSCI NA PROBE
! mh 22-11-07 w problemach z udzialem promieniowania kosmicznego
! uzycie vx prowadzi do nieduzych oscylacji, a dla vt oscylacji nie widac,
! wobec tego po krotkim powrocie do vx (v.1.2) przywracam vt (v.1.3)

    vy(2:n-1)=u(imya,2:n-1)/u(idna,2:n-1)
    vz(2:n-1)=u(imza,2:n-1)/u(idna,2:n-1)
    vt = sqrt(vx*vx+vy*vy+vz*vz)
    
    
#ifdef ISO
    p(2:n-1) = csi2*u(idna,2:n-1)
    ps(2:n-1)= p(2:n-1) + pmag(2:n-1)
#else ISO    
    ps(2:n-1)=(u(iena,2:n-1)-(u(imxa,2:n-1)*u(imxa,2:n-1)+u(imya,2:n-1)*u(imya,2:n-1) & 
              +u(imza,2:n-1)*u(imza,2:n-1))/u(idna,2:n-1)*0.5)*(gamma-1)+(2-gamma)*pmag(2:n-1)
    p(2:n-1)=ps(2:n-1)-pmag(2:n-1)
#endif ISO
    flux(idna,2:n-1)=u(imxa,2:n-1)
    flux(imxa,2:n-1)=u(imxa,2:n-1)*vx(2:n-1)+ps(2:n-1)-b(ibx,2:n-1)*b(ibx,2:n-1)
    flux(imya,2:n-1)=u(imya,2:n-1)*vx(2:n-1)-b(iby,2:n-1)*b(ibx,2:n-1)
    flux(imza,2:n-1)=u(imza,2:n-1)*vx(2:n-1)-b(ibz,2:n-1)*b(ibx,2:n-1)
#ifndef ISO
    flux(iena,2:n-1)=(u(iena,2:n-1)+ps(2:n-1))*vx(2:n-1)-b(ibx,2:n-1)*(b(ibx,2:n-1)*u(imxa,2:n-1) &
                +b(iby,2:n-1)*u(imya,2:n-1)+b(ibz,2:n-1)*u(imza,2:n-1))/u(idna,2:n-1)
#endif ISO
#ifdef COSM_RAYS
    flux(iecr,2:n-1)= u(iecr,2:n-1)*vx(2:n-1)
#endif COSM_RAYS
#ifdef LOCAL_FR_SPEED

!       The freezing speed is now computed locally (in each cell) 
!       as in Trac & Pen (2003). This ensures much sharper shocks, 
!       but sometimes may lead to numerical instabilities   
#ifdef ISO
! UWAGA ZMIANA W OBLICZANIU LOKALNEJ PREDKOSCI NA PROBE
!        cfr(1,2:n-1) = abs(vx(2:n-1)) &
        cfr(1,2:n-1) = abs(vt(2:n-1)) &
                      +max(sqrt( abs(2*pmag(2:n-1) + p(2:n-1))/u(idna,2:n-1)),small)
#else ISO   
! UWAGA ZMIANA W OBLICZANIU LOKALNEJ PREDKOSCI NA PROBE
!        cfr(1,2:n-1) = abs(vx(2:n-1)) &
        cfr(1,2:n-1) = abs(vt(2:n-1)) &
                      +max(sqrt( abs(2*pmag(2:n-1) + gamma*p(2:n-1))/u(idna,2:n-1)),small)
#endif ISO
        cfr(1,1) = cfr(1,2)
        cfr(1,n) = cfr(1,n-1)   
        cfr = spread(cfr(1,:),1,nu)
#endif LOCAL_FR_SPEED
#ifdef GLOBAL_FR_SPEED
!       The freezing speed is now computed globally 
!       (c=const for the whole domain) in sobroutine 'timestep' 
!       Original computation of the freezing speed was done
!       for each sweep separately:
!       c=maxval(abs(vx(nb-2:n-1))+sqrt(abs(2*pmag(nb-2:n-1) &
!                            +gamma*p(nb-2:n-1))/u(idna,nb-2:n-1)))  
        cfr(:,:) = c
#endif GLOBAL_FR_SPEED

  end subroutine mhdflux

!==========================================================================================
 

  subroutine flimiter(f,a,b,m,n)
    implicit none
    integer m,n
    real, dimension(m,n) :: f,a,b,c
#ifdef VANLEER
      c = a*b
      where (c .gt. 0)                                        
        f = f+2*c/(a+b)
      endwhere      
#endif VANLEER
#ifdef MINMOD
        f = f+(sign(1.,a)+sign(1.,b))*min(abs(a),abs(b))*0.5
#endif MINMOD
#ifdef SUPERBEE
      where (abs(a) .gt. abs(b))
        f = f+(sign(1.,a)+sign(1.,b))*min(abs(a), abs(2*b))*0.5
      elsewhere
        f = f+(sign(1.,a)+sign(1.,b))*min(abs(2*a), abs(b))*0.5
      endwhere
#endif SUPERBEE
 
    return
  end subroutine flimiter

end module fluxes
