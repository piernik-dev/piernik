! $Id$
#include "piernik.def"
module sweeps     ! split sweeps
  contains
  
  
  subroutine sweepx

    use fluidindex,   only : nvar, iarr_all_swpx
    use fluidindex,   only : nmag
    use fluidindex,   only : ibx,iby,ibz
#ifdef IONIZED
    use fluidindex,   only : i_ion
#endif /* IONIZED */  

    use start, only  : dt
    use arrays, only : u,b
    use grid, only   : dx,nb,nx,ks,ke,js,je,nyd,nzd
    use rtvd, only   : relaxing_tvd
    use fluidboundaries, only : compute_u_bnd 

#ifdef COSM_RAYS
    use func,  only : div_v        
#endif /* COSM_RAYS */

    implicit none
    real, dimension(nmag,nx)  :: b_x
    real, dimension(nvar,nx)  :: u_x
    integer :: j,k,jp,kp
  
    b_x = 0.0
    u_x = 0.0
    
#ifdef COSM_RAYS
    call div_v(i_ion)         
#endif /* COSM_RAYS */

    do k=ks,ke
      kp=k+1
      do j=js,je
          jp=j+1
      
#ifdef MAGNETIC
          b_x=0.5*b(:,:,j,k)
          b_x(ibx,1:nx-1)=b_x(ibx,1:nx-1)+b_x(ibx,2:nx);       b_x(ibx,nx) = b_x(ibx,nx-1)
          if(nyd /= 1 .and. j < je)  b_x(iby,:)=b_x(iby,:)+0.5*b(iby,:,jp,k)
          if(nzd /= 1 .and. k < ke)  b_x(ibz,:)=b_x(ibz,:)+0.5*b(ibz,:,j,kp)
#endif /* MAGNETIC */

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
    use fluidindex,   only : nmag
    use fluidindex, only : ibx,iby,ibz
#ifdef IONIZED
    use fluidindex, only : i_ion
#endif /* IONIZED */  
       
    use start, only  : dt
    use arrays, only : u,b
    use grid, only   : dy,nb,ny,is,ie,ks,ke,nxd,nzd
    use rtvd, only     : relaxing_tvd
    use fluidboundaries, only : compute_u_bnd 
    
#ifdef COSM_RAYS
    use func, only   : div_v        
#endif /* COSM_RAYS */

    implicit none
    real, dimension(nmag,ny)  :: b_y
    real, dimension(nvar,ny)  :: u_y
    
    integer :: i,k,ip,kp
 
    b_y = 0.0
    u_y = 0.0
    
#ifdef COSM_RAYS
    call div_v(i_ion)         
#endif /* COSM_RAYS */

    do k=ks,ke
      kp=k+1
      do i=is,ie
          ip=i+1
        
#ifdef MAGNETIC	
          b_y(:,:) = 0.5*b(:,i,:,k)
          b_y(iby,1:ny-1)=b_y(iby,1:ny-1)+b_y(iby,2:ny);       b_y(iby,ny) = b_y(iby,ny-1)
          if(nxd /= 1 .and. i < ie) b_y(ibx,:)=b_y(ibx,:)+0.5*b(ibx,ip,:,k)
          if(nzd /= 1 .and. k < ke) b_y(ibz,:)=b_y(ibz,:)+0.5*b(ibz,i,:,kp)
          b_y((/iby,ibx,ibz/),:)=b_y(:,:)
#endif /* MAGNETIC */

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
    use fluidindex, only   : nmag
    use fluidindex, only : ibx,iby,ibz
#ifdef IONIZED
    use fluidindex, only : i_ion
#endif /* IONIZED */  
    

    use start, only  : dt
    use arrays, only : u,b
    use grid, only   : dz,nb,nz,is,ie,js,je,nxd,nyd
    use rtvd, only     : relaxing_tvd
    use fluidboundaries, only : compute_u_bnd 
    
#ifdef COSM_RAYS
    use func,  only : div_v       
#endif /* COSM_RAYS */

    implicit none
    real, dimension(nmag,nz)  :: b_z
    real, dimension(nvar,nz)  :: u_z
    
    integer :: i,j,ip,jp

    b_z = 0.0
    u_z = 0.0
    
#ifdef COSM_RAYS
    call div_v(i_ion)         
#endif /* COSM_RAYS */

    do j=js,je
      jp=j+1
      do i=is,ie
          ip=i+1
        
#ifdef MAGNETIC	
          b_z(:,:) = 0.5*b(:,i,j,:)
          b_z(ibz,1:nz-1) = b_z(ibz,1:nz-1) + b_z(ibz,2:nz);   b_z(ibz,nz) = b_z(ibz,nz-1)
          if(nxd /= 1 .and. i < ie) b_z(ibx,:) = b_z(ibx,:) + 0.5*b(ibx,ip,j,:)
          if(nyd /= 1 .and. j < je) b_z(iby,:) = b_z(iby,:) + 0.5*b(iby,i,jp,:)
          b_z((/ibz,iby,ibx/),:)=b_z(:,:)
#endif /* MAGNETIC */

        u_z(iarr_all_swpz,:)=u(:,i,j,:)

        call relaxing_tvd(u_z,b_z,'zsweep',i,j,dz,nz,dt)
        u(:,i,j,:)=u_z(iarr_all_swpz,:)
      end do
    end do
    
    call compute_u_bnd

  end subroutine sweepz
end module sweeps
