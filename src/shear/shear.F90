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
  end subroutine yshift

  function unshear(qty,x,inv)
    use start, only  : nb,xmax,xmin,nyd,smalld
    use grid, only   : dy
    
    logical, optional               :: inv
    real, dimension(:,:,:)          :: qty
    real, dimension(:), intent(in)  :: x
    real, dimension(size(qty,1),size(qty,2),size(qty,3)) :: unshear
    real, dimension(:,:), allocatable:: temp
    integer :: i,sg,my,nx,ny,nz
    real    :: fx

    nx = size(qty,1)
    ny = size(qty,2)
    nz = size(qty,3)

    my = 3*nyd+2*nb

    fx = dely / (xmax - xmin)
    sg = -1

    if(.not.allocated(temp)) allocate(temp(my,nz))

    unshear = 0.0

    if(present(inv)) then
       fx = - fx
    endif
    do i = 1, nx
      dl  = fx * x(i)
      ndl = mod(int(dl/dy),nyd)
      ddl = mod(dl,dy)/dy

      temp(         1:  nyd+nb,:)   = qty(i,   1:nyd+nb ,:)
      temp(  nyd+nb+1:2*nyd+nb,:)   = qty(i,nb+1:nyd+nb,:)
      temp(2*nyd+nb+1:3*nyd+2*nb,:) = qty(i,nb+1:ny    ,:)

      temp = cshift(temp,dim=1,shift=ndl)

!      temp(1:nb,:) = temp(nyb+1:nyb+nb,:)          ! not needed
!      temp(nb+nyb+1:nyb+2*nb,:) = temp(nb+1:2*nb,:)

      temp(:,:) = (1.0+ddl)*(1.0-ddl) * temp(:,:) &
            - 0.5*(ddl)*(1.0-ddl) * cshift(temp(:,:),shift= sg,dim=1) &
            + 0.5*(ddl)*(1.0+ddl) * cshift(temp(:,:),shift=-sg,dim=1) 

      unshear(i,nb+1:nb+nyd,:) = temp(nb+nyd+1:nb+2*nyd,:)

      unshear(i,1:nb,:)          = unshear(i,nyd+1:nyd+nb,:)
      unshear(i,nyd+nb+1:ny,:)   = unshear(i,nb+1 :2*nb,:)

      unshear(i,:,:) = max(unshear(i,:,:), smalld)
    enddo
    if (allocated(temp)) deallocate(temp)
    return
  end function
end module shear

