! $Id$
#include "piernik.def"
#ifdef SPLIT
#define UT ufc
#define UT_IMZA ufc(imxa,k)
#define UT_IDNA ufc(idna,k)
#define U_IECR uu(iecr,:)
#define VZ_LIM_DN uu(idna,:)
#define VZ_LIM_VZ uu(imxa,:)
#else /* SPLIT */
#define UT u
#define UT_IMZA u(imza,:,:,k)
#define UT_IDNA u(idna,:,:,k)
#define U_IECR u(iecr,:,:,:)
#define VZ_LIM_DN u(idna,:)
#define VZ_LIM_VZ u(imza,:)
#endif /* SPLIT */

module floor_ceil

contains

!================GENERAL SUBROUTINE OF FLOOR_CEIL MODULE=======================
subroutine floor_ceiling(callnumb,n)
#ifndef SPLIT
  use arrays, only : u
#endif /* SPLIT */
  use arrays, only : nu
  implicit none
  integer n,i1,i2,callnumb
  character sweep*6
  real,dimension(nu,n) :: ufc

#ifdef COSM_RAYS
!    ufc(iecr,:) = max(smallecr,ufc(iecr,:))
#endif /* COSM_RAYS */

#ifdef VZ_LIMITS
    if(sweep .eq. 'zsweep') then
      call vz_limits(ufc,n)
    endif
#endif /* VZ_LIMITS */

end subroutine floor_ceiling

!===============END OF GENERAL SUBROUTINE OF FLOOR_CEIL MODULE=====================

!================================ELEMENTARY SUBROUTINES===================================

#ifdef VZ_LIMITS
subroutine vz_limits(ufc,n)
    use arrays, only : z,idna,imxa,imza,nz,nu
    use start, only : floor_vz, ceil_vz
#ifndef SPLIT
    use arrays, only : u
#endif /* !SPLIT */
    implicit none
    integer :: n
    real,dimension(nu,n) :: ufc
    integer k

      do k=1,nz
        if(z(k) .gt. 0.0) then
          where(UT_IMZA .lt.  floor_vz*UT_IDNA)
            UT_IMZA =  floor_vz*UT_IDNA
          endwhere
        endif
        if(z(k) .lt. 0.0) then
          where(UT_IMZA .gt. -floor_vz*UT_IDNA)
            UT_IMZA = -floor_vz*UT_IDNA
          endwhere
        endif
      enddo
end subroutine vz_limits
#endif /* VZ_LIMITS */

!==========================END OF ELEMENTARY SUBROUTINES====================================

end module floor_ceil
