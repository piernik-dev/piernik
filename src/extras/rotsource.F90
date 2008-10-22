!$Id$
#include "piernik.def"
module rotsource
#ifdef KEPLER_SUPPRESSION
contains
subroutine kepler_suppression(Duus,uu,sweep,i1,i2,n,dt)
    use arrays, only : alfsup,omx0,omy0,nu,x,y,z,nx,ny,idna,imxa,imya
    use start,  only : maxxyz
#ifndef ISO
  use arrays, only : nadiab,fadiab,iena
#endif /* !ISO */
  implicit none
  integer n,i,j,k,i1,i2,ifl,ni1,ni2
  character sweep*6
  real,dimension(nu,n)      :: Duus,uu
  real, dimension(n)        :: kplrsup,velcoor,vel0
  real, dimension(maxxyz)   :: r1,r2
  real :: a1,dt
  integer :: rend,ii

    Duus=0.0
    select case(sweep)
      case('xsweep')
        a1    = y(i1)
        rend  = size(x)
        r1(1:rend) = x(:)                            ! r1   max(size(x),size(y),size(z))
        r2(1:rend) = a1 / (r1(1:rend)*r1(1:rend) + a1 * a1)
	vel0(:)  = omx0(:,i1,i2)
      case('ysweep')
        a1    = x(i2)
        rend  = size(y)
        r1(1:rend) = y(1:rend)
        r2(1:rend) = a1 / (r1(1:rend)*r1(1:rend) + a1 * a1)
	vel0(:)  = omy0(i2,:,i1)
      case('zsweep')
        a1    = 1.0
        rend  = size(z)
        r1(1:rend) = 0.0
        r2(1:rend) = 1.0
    end select

    if(sweep .ne. 'zsweep') then
#ifdef KEPL_SUPP_SIMX
        velcoor(:)=uu(imxa,:)/uu(idna,:)
	kplrsup(:)=-alfsup(:,i1)*(velcoor(:)-vel0(:)) !*uu(idna,:)
#else /* KEPL_SUPP_SIMX */
        velcoor(:)=(uu(imxa,:)*a1-uu(imya,:)*r1)/uu(idna,:)*r2
        kplrsup(:)=-alfsup(:,i1)*(velcoor(:)-vel0(:))!*uu(idna,:)
#endif /* KEPL_SUPP_SIMX */
    else
      kplrsup(:)=0.0
    endif
#ifndef ISO
!    Duus(iena(fadiab),:)=(uu(imxa(fadiab),:)*kplrsup*uu(idna(fadiab),:)*dt &
!            +0.5*(kplrsup*uu(idna(fadiab),:)*dt)**2)/uu(idna(fadiab),:)
    Duus(iena,:)=kplrsup(:)*uu(imxa,:)*dt
#endif /* ISO */
    Duus(imxa,:)=kplrsup(:)*uu(idna,:)*dt

end subroutine kepler_suppression
#endif /* KEPLER_SUPPRESSION */
end module rotsource
