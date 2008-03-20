
#include "mhd.def"
module fluids     ! split fluids
  contains
  
  
  subroutine fluidx

    use start, only  : dimensions, magfield,dt,nb
    use arrays, only : u,b,nx,ny,nz,nu,ks,ke,ibx,iby,ibz,iuswpx
    use grid, only   : dx
    use tv, only     : relaxing_tvd
    use fluid_boundaries, only : compute_u_bnd 
#ifdef COSM_RAYS
    use cr_diffusion        
#endif /* COSM_RAYS */

    implicit none
    real, dimension(3,nx)  :: b_x
    real, dimension(nu,nx) :: u_x
    integer :: j,k,jp,kp
    real    :: yj,zk
  
    b_x = 0.0
    u_x = 0.0
    
#ifdef COSM_RAYS
    call div_v         
#endif /* COSM_RAYS */

    do k=ks,ke
      kp=k+1
      do j=1,ny-1
        if(magfield)then
          jp=j+1
          b_x=0.5*b(:,:,j,k)
          b_x(ibx,1:nx-1)=b_x(ibx,1:nx-1)+b_x(ibx,2:nx);       b_x(ibx,nx) = b_x(ibx,nx-1)
          b_x(iby,:)=b_x(iby,:)+0.5*b(iby,:,jp,k)
          if(dimensions .eq. '3d')  b_x(ibz,:)=b_x(ibz,:)+0.5*b(ibz,:,j,kp)
        else  ! tak "just in case"
          b_x = 0.0
        endif

        u_x(iuswpx,:)=u(:,:,j,k)
        call relaxing_tvd(u_x,b_x,'xsweep',j,k,dx,nx,dt)
        u(:,:,j,k)=u_x(iuswpx,:)
      end do
    end do

    call compute_u_bnd

  end subroutine fluidx

!------------------------------------------------------------------------------------------

  subroutine fluidy
    use start, only  : dimensions, magfield,dt,nb
    use arrays, only : u,b,nx,ny,nz,nu,ks,ke,ibx,iby,ibz,iuswpy
    use grid, only   : dy
    use tv, only     : relaxing_tvd
    use fluid_boundaries, only : compute_u_bnd 
#ifdef COSM_RAYS
    use cr_diffusion        
#endif /* COSM_RAYS */

    implicit none
    real, dimension(3,ny)  :: b_y
    real, dimension(nu,ny) :: u_y
    
    integer :: i,j,k,ip,jp,kp
    real    :: xi,zk
 
    b_y = 0.0
    u_y = 0.0
    
#ifdef COSM_RAYS
    call div_v         
#endif /* COSM_RAYS */

    do k=ks,ke
      kp=k+1
      do i=1,nx-1
        if(magfield)then
          ip=i+1
          b_y(:,:)=b(:,i,:,k)/2
          b_y(ibx,:)=b_y(ibx,:)+b(ibx,ip,:,k)/2
          b_y(iby,1:ny-1)=b_y(iby,1:ny-1)+b_y(iby,2:ny);       b_y(iby,ny) = b_y(iby,ny-1)
          if (dimensions .eq. '3d') b_y(ibz,:)=b_y(ibz,:)+b(ibz,i,:,kp)/2
          b_y((/iby,ibx,ibz/),:)=b_y(:,:)
        else
          b_y = 0.0
        endif

        u_y(iuswpy,:)=u(:,i,:,k) 

        call relaxing_tvd(u_y,b_y,'ysweep',k,i,dy,ny,dt)
        u(:,i,:,k)=u_y(iuswpy,:)

      end do
    end do

    call compute_u_bnd

  end subroutine fluidy

!------------------------------------------------------------------------------------------

  subroutine fluidz
    use start, only  : dimensions, magfield,dt,nb
    use arrays, only : u,b,nx,ny,nz,nu,ks,ke,ibx,iby,ibz,iuswpz
    use grid, only   : dz
    use tv, only     : relaxing_tvd
    use fluid_boundaries, only : compute_u_bnd 
#ifdef COSM_RAYS
    use cr_diffusion        
#endif /* COSM_RAYS */

    implicit none
    real, dimension(3,nz)  :: b_z
    real, dimension(nu,nz) :: u_z
    
    integer :: i,j,k,ip,jp,kp
    real    :: xi,yj

    b_z = 0.0
    u_z = 0.0
    
#ifdef COSM_RAYS
    call div_v         
#endif /* COSM_RAYS */

    do j=1,ny-1
      jp=j+1
      do i=1,nx-1
        if(magfield)then
          ip=i+1
          b_z(:,:)=b(:,i,j,:)/2
          b_z(ibx,:)=b_z(ibx,:)+b(ibx,ip,j,:)/2
          b_z(iby,:)=b_z(iby,:)+b(iby,i,jp,:)/2
          b_z(ibz,1:nz-1)=b_z(ibz,1:nz-1)+b_z(ibz,2:nz);       b_z(ibz,nz) = b_z(ibz,nz-1)
          b_z((/ibz,iby,ibx/),:)=b_z(:,:)
        else
          b_z = 0.0
        endif
        u_z(iuswpz,:)=u(:,i,j,:)

        call relaxing_tvd(u_z,b_z,'zsweep',i,j,dz,nz,dt)
        u(:,i,j,:)=u_z(iuswpz,:)
      end do
    end do
    
    call compute_u_bnd

  end subroutine fluidz
end module fluids
