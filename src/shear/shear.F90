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

  function unshear(qty,x,inv)
    use start, only  : nb,xmax,xmin,nyd
    use grid, only   : dy
    
    logical, optional               :: inv
    real, dimension(:,:,:)          :: qty
    real, dimension(:), intent(in)  :: x
    real, dimension(size(qty,1),size(qty,2),size(qty,3)) :: unshear
    integer :: i,sg
    real    :: fx

    fx = dely / (xmax - xmin)
    sg = -1

    if(present(inv)) then
       fx = - fx
    endif
    do i = 1, size(qty,1)
      dl  = fx * x(i)
      ndl = mod(int(dl/dy),nyd)
      ddl = mod(dl,dy)/dy

      unshear(i,:,:) = cshift(qty(i,:,:),dim=1,shift=ndl)

      unshear(i,1:nb,:) = unshear(i,nyb+1:nyb+nb,:)
      unshear(i,nb+nyb+1:nyb+2*nb,:) = unshear(i,nb+1:2*nb,:)

     
      unshear(i,:,:) = (1.0+ddl)*(1.0-ddl) * unshear(i,:,:) &
            - 0.5*(ddl)*(1.0-ddl) * cshift(unshear(i,:,:),shift= sg,dim=1) &
            + 0.5*(ddl)*(1.0+ddl) * cshift(unshear(i,:,:),shift=-sg,dim=1) 

    enddo
    return
  end function
end module shear

