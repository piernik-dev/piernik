! $Id$
#include "piernik.def"
#define SQR(var) var*var
#define SUM_SQR(x,y,z) ( SQR(x)+SQR(y)+SQR(z) )
#define RNG 2:n-1

module cosmic_rays

contains

subroutine flux_cr(flux,cfr,uu,bb,n)
    use arrays, only : idni,imxi,nfi,iecr,nu
    implicit none
    integer n
!    real, dimension(COSM_RAYS,n) :: flux
    real, dimension(nu,n)        :: uu,cfr,flux
    real, dimension(3,n)         :: bb
    real, dimension(nfi,n)       :: vx
!#ifdef COSM_RAYS
!    pcr = (gamma_cr-1)*uu(iecr,RNG)
!#endif COSM_RAYS
    vx(:,RNG)=uu(imxi(1:nfi),RNG)/uu(idni(1:nfi),RNG)
    flux(iecr,RNG)= uu(iecr,RNG)*spread(vx(1,RNG),1,COSM_RAYS)  ! niegotowe, vx to predkosc plynu zjonizowanego?
    cfr(iecr,:)=spread(cfr(idni(1),:),1,COSM_RAYS)

end subroutine flux_cr

subroutine src_cr
!    use start,  only : gamma_cr, cr_active, smallecr
!    use arrays, only : divvel
    implicit none
    real, dimension(nfluid,n) :: divv
    real, dimension(n)    :: decr,gpcr,ecr,tmp

    select case (sweep)
      case('xsweep')
        divv = divvel(:,:,i1,i2)
      case('ysweep')
        divv = divvel(:,i2,:,i1)
      case('zsweep')
        divv = divvel(:,i1,i2,:)
    end select

    Duu(iecr,:) = Duu(iecr,:) - (gamma_cr-1.)*uu(iecr,:)*divv(1,:)*dt

    vx  = uu(imxa,:)/uu(idna,:)
    ecr = uu(iecr,:)

    gpcr(2:n-1) = cr_active*(gamma_cr -1.)*(ecr(3:n)-ecr(1:n-2))/(2.*dx) ; gpcr(1:2)=0.0 ; gpcr(n-1:n) = 0.0

#ifndef ISO
    Duu(iena,:) = Duu(iena,:) - uu(imxa,:)/uu(idna,:)*gpcr*dt
#endif /* ISO */
    Duu(imxa,:) = Duu(imxa,:) - gpcr*dt

end subroutine src_cr

end module cosmic_rays
