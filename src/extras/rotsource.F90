!$Id$
#include "piernik.def"
module rotsource
#ifdef KEPLER_SUPPRESSION
   contains
   subroutine kepler_suppression(Duus,uu,sweep,i1,i2,n,dt)
      use arrays, only : alfsup,omx0,omy0,nu,x,y,z,nx,ny,idna,imxa,imya,den0
      use start,  only : maxxyz
#ifndef ISO
      use arrays, only : nadiab,fadiab,iena
#endif /* !ISO */
      implicit none
      integer :: n, i1, i2
      character sweep*6
      real,dimension(nu,n)      :: Duus,uu
      real, dimension(n)        :: kplrsup, vel0, dns0
#ifdef KEPL_SUPP_SIMX
      real, dimension(nu,n)     :: velcoor
#endif /* KEPL_SUPP_SIMX */
      real, dimension(maxxyz)   :: r1,r2,alf
      real :: a1,dt
      integer :: rend

      Duus=0.0
      select case(sweep)
         case('xsweep')
            a1    = y(i1)
            rend  = size(x)
            r1(1:rend) = x(:)                            ! r1   max(size(x),size(y),size(z))
            r2(1:rend) = a1 / (r1(1:rend)*r1(1:rend) + a1 * a1)
            alf(1:rend) = alfsup(:,i1)
            vel0(:)  = omx0(:,i1,i2)
            dns0(:)  = den0(:,i1,i2)
         case('ysweep')
            a1    = x(i2)
            rend  = size(y)
            r1(1:rend) = y(1:rend)
            r2(1:rend) = a1 / (r1(1:rend)*r1(1:rend) + a1 * a1)
            alf(1:rend) = alfsup(i2,:)
            vel0(:)  = omy0(i2,:,i1)
            dns0(:)  = den0(i2,:,i1)
         case('zsweep')
            a1    = 1.0
            rend  = size(z)
            r1(1:rend) = 0.0
            r2(1:rend) = 1.0
      end select

      if(sweep .ne. 'zsweep') then
#ifdef OVERLAP
         uu(idna,:) = (1.-alf)*uu(idna,:) + alf*dns0
         uu(imxa,:) = (1.-alf)*uu(imxa,:) + alf*vel0*dns0
#else /* OVERLAP */
#ifdef KEPL_SUPP_SIMX
         velcoor(:)=uu(imxa,:)/uu(idna,:)
         kplrsup(:)=-alfsup(:,i1)*(velcoor(:)-vel0(:)) !*uu(idna,:)
#else /* KEPL_SUPP_SIMX */
         velcoor(:)=(uu(imxa,:)*a1-uu(imya,:)*r1)/uu(idna,:)*r2
         kplrsup(:)=-alfsup(:,i1)*(velcoor(:)-vel0(:))!*uu(idna,:)
#endif /* KEPL_SUPP_SIMX */
#endif /* OVERLAP */
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
