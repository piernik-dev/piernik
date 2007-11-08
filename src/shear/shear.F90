#include "mhd.def"
module shear
  real    :: ts, dely, eps
  integer :: delj

  contains

  subroutine yshift(ts)
    use start, only : qshear,omega,xmin,xmax,nyd
    use grid, only : dy
    implicit none
    real, intent(in) :: ts

    dely  = ts*qshear*omega*(xmax-xmin)
    delj  = mod(int(dely/dy),nyd)
    eps   = mod(dely,dy)/dy
!     write(*,*) ts,delj,eps
  end subroutine yshift
end module shear

