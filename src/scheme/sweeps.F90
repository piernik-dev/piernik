! $Id: fluids.F90 362 2008-10-16 17:51:20Z wolt $
#include "piernik.def"
module sweeps     ! split sweeps
  contains
  
  
  subroutine sweepx

    use fluidindex,   only : nvar, iarr_all_swpx
#ifdef IONIZED
    use fluidindex,   only : nmag, ibx,iby,ibz
#endif /* IONIZED */
    use start, only  : dimensions, magfield,dt,nb
    use arrays, only : u,b,nx,ny,nz,ks,ke
    use grid, only   : dx
    use rtvd, only     : relaxing_tvd
    use fluid_boundaries, only : compute_u_bnd 
#ifdef COSM_RAYS
    use cr_diffusion        
#endif /* COSM_RAYS */

    implicit none
    real, dimension(nmag,nx)  :: b_x
    real, dimension(nvar,nx)  :: u_x
    integer :: j,k,jp,kp
  
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

        u_x(iarr_all_swpx,:)=u(:,:,j,k)

        call relaxing_tvd(u_x,b_x,'xsweep',j,k,dx,nx,dt)
        u(:,:,j,k)=u_x(iarr_all_swpx,:)
      end do
    end do

    call compute_u_bnd

  end subroutine sweepx

!------------------------------------------------------------------------------------------

  subroutine sweepy
    use fluidindex, only : nvar, iarr_all_swpy
#ifdef IONIZED
    use fluidindex,   only : nmag,ibx,iby,ibz
#endif /* IONIZED */
    use start, only  : dimensions, magfield,dt,nb
    use arrays, only : u,b,nx,ny,nz,ks,ke
    use grid, only   : dy
    use rtvd, only     : relaxing_tvd
    use fluid_boundaries, only : compute_u_bnd 
#ifdef COSM_RAYS
    use cr_diffusion        
#endif /* COSM_RAYS */

    implicit none
    real, dimension(nmag,ny)  :: b_y
    real, dimension(nvar,ny)  :: u_y
    
    integer :: i,k,ip,kp
 
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

        u_y(iarr_all_swpy,:)=u(:,i,:,k) 

        call relaxing_tvd(u_y,b_y,'ysweep',k,i,dy,ny,dt)
        u(:,i,:,k)=u_y(iarr_all_swpy,:)

      end do
    end do

    call compute_u_bnd

  end subroutine sweepy

!------------------------------------------------------------------------------------------

  subroutine sweepz
    use fluidindex, only   : nvar, iarr_all_swpz
#ifdef IONIZED
    use fluidindex, only   : nmag, ibx, iby, ibz
#endif /* IONIZED */
    use start, only  : dimensions, magfield,dt,nb
    use arrays, only : u,b,nx,ny,nz,ks,ke
    use grid, only   : dz
    use rtvd, only     : relaxing_tvd
    use fluid_boundaries, only : compute_u_bnd 
#ifdef COSM_RAYS
    use cr_diffusion        
#endif /* COSM_RAYS */

    implicit none
    real, dimension(nmag,nz)  :: b_z
    real, dimension(nvar,nz)  :: u_z
    
    integer :: i,j,ip,jp

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
        u_z(iarr_all_swpz,:)=u(:,i,j,:)

        call relaxing_tvd(u_z,b_z,'zsweep',i,j,dz,nz,dt)
        u(:,i,j,:)=u_z(iarr_all_swpz,:)
      end do
    end do
    
    call compute_u_bnd

  end subroutine sweepz
end module sweeps
