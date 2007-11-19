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

  function unshear(qty,x,ds,inv)
    use start, only  : nb
    use arrays, only : nyb
    logical, optional               :: inv
    real, dimension(:,:,:)          :: qty
    real, intent(in)                :: ds
    real, dimension(:), intent(in)  :: x
    real, dimension(size(qty,1),size(qty,2),size(qty,3))      :: unshear

    integer :: nx,i,sg
    real    :: lx,fx

    nx = size(qty,1)

    lx = maxval(x) - minval(x) + x(1) - x(0)
    fx = ds / lx
    sg = -1

    if(present(inv)) then
      fx = - fx
      sg = 1
    endif

    do i = 1, nx
      dl  = fx * x(i)
      ndl = int(dl)
      ddl = dl - ndl
      unshear(i,:,:) = cshift( qty(i,:,:),dim=2, shift=ndl)
      unshear(i,1:nb,:) = unshear(i,nyb+1:nyb+nb,:)
      unshear(i,nb+nyb+1:nyb+2*nb,:) = unshear(i,nb+1:2*nb,:)
      
      unshear(i,:,:) = (1.0+ddl)*(1.0-ddl) * unshear(i,:,:) &
            - 0.5*(ddl)*(1.0-ddl) * cshift(unshear(i,:,:),shift= sg,dim=3) &
            + 0.5*(ddl)*(1.0+ddl) * cshift(unshear(i,:,:),shift=-sg,dim=3) 

    enddo
    return
  end function
end module shear

