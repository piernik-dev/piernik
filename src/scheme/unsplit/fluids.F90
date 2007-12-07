
#include "mhd.def"
module fluids  ! unsplit fluids
  use grid
  use start
  use arrays
  use tv
  contains
!==========================================================================================

  subroutine fluidx
    implicit none
    real, dimension(3,nx)  :: b_x
    real, dimension(nu,nx) :: u_x
    
#ifdef ORIG
    real, dimension(3,nx)  :: b_i
    real, dimension(nu,nx) :: u_i
#endif /* ORIG */

    integer j,k,jp,kp
    real yj,zk

  
    b_x = 0.0
#ifdef ORIG
    b_i = 0.0
#endif /* ORIG */

#ifdef COSM_RAYS
   call div_v         
#endif /* COSM_RAYS */

    do k=ks,ke
      kp=k+1
      do j=1,ny-1
        if(magfield)then
          jp=j+1
          b_x=b(:,:,j,k)*0.5
          b_x(ibx,:)=b_x(ibx,:)+eoshift(b_x(ibx,:),shift=1,boundary=big)
          b_x(iby,:)=b_x(iby,:)+b(iby,:,jp,k)*0.5
          if(dimensions .eq. '3d')  b_x(ibz,:)=b_x(ibz,:)+b(ibz,:,j,kp)*0.5
#ifdef ORIG
          b_i=bi(:,:,j,k)*0.5
          b_i(ibx,:)=b_i(ibx,:)+eoshift(b_i(ibx,:),shift=1,boundary=big)
          b_i(iby,:)=b_i(iby,:)+bi(iby,:,jp,k)*0.5
          if(dimensions .eq. '3d')  b_i(ibz,:)=b_i(ibz,:)+bi(ibz,:,j,kp)*0.5
#endif /* ORIG */
        else
          b_x = 0.0
#ifdef ORIG
          b_i = 0.0
#endif /* ORIG */
        endif

        u_x(iuswpx,:)=u(:,:,j,k)
#ifdef ORIG
        u_i(iuswpx,:)=ui(:,:,j,k)
        call relaxing_tvd(u_i,u_x,b_i,b_x,'xsweep',j,k,dx,nx,dt)
#else /* ORIG */
        call relaxing_tvd(u_x,b_x,'xsweep',j,k,dx,nx,dt)
#endif /* ORIG */
        Lu(:,:,j,k)=Lu(:,:,j,k)+u_x(iuswpx,:)
      end do
    end do


  end subroutine fluidx

!------------------------------------------------------------------------------------------

  subroutine fluidy

    implicit none
    real, dimension(3,ny)  :: b_y
    real, dimension(nu,ny) :: u_y
    
#ifdef ORIG
    real, dimension(3,ny)  :: b_i
    real, dimension(nu,ny) :: u_i
#endif /* ORIG */
    integer i,j,k,ip,jp,kp
    real xi,zk

   
    b_y = 0.0
#ifdef ORIG
    b_i = 0.0
#endif /* ORIG */
#ifdef COSM_RAYS
   call div_v         
#endif /* COSM_RAYS */

    do k=ks,ke
      kp=k+1
      do i=1,nx-1
        if(magfield)then
          ip=i+1
          b_y(:,:)=b(:,i,:,k)*0.5
          b_y(ibx,:)=b_y(ibx,:)+b(ibx,ip,:,k)*0.5
          b_y(iby,:)=b_y(iby,:)+eoshift(b_y(iby,:),shift=1,boundary=big)
          if (dimensions .eq. '3d') b_y(ibz,:)=b_y(ibz,:)+b(ibz,i,:,kp)*0.5
          b_y((/iby,ibx,ibz/),:)=b_y(:,:)
#ifdef ORIG
          b_i(:,:)=bi(:,i,:,k)*0.5
          b_i(ibx,:)=b_i(ibx,:)+bi(ibx,ip,:,k)*0.5
          b_i(iby,:)=b_i(iby,:)+eoshift(b_i(iby,:),shift=1,boundary=big)
          if (dimensions .eq. '3d') b_i(ibz,:)=b_i(ibz,:)+bi(ibz,i,:,kp)*0.5
          b_i((/iby,ibx,ibz/),:)=b_i(:,:)
#endif /* ORIG */
        else
          b_y = 0.0
#ifdef ORIG 
          b_i = 0.0
#endif /* ORIG */
        endif

        u_y(iuswpy,:)=u(:,i,:,k) 

#ifdef ORIG
        u_i(iuswpy,:)=ui(:,i,:,k) 
        call relaxing_tvd(u_i,u_y,b_i,b_y,'ysweep',k,i,dy,ny,dt)
#else /* ORIG */
        call relaxing_tvd(u_y,b_y,'ysweep',k,i,dy,ny,dt)
#endif /* ORIG */
        Lu(:,i,:,k)=Lu(:,i,:,k)+u_y(iuswpy,:)

      end do
    end do

  end subroutine fluidy

!------------------------------------------------------------------------------------------

  subroutine fluidz

    implicit none
    real, dimension(3,nz)  :: b_z
    real, dimension(nu,nz) :: u_z
    
#ifdef ORIG
    real, dimension(3,nz)  :: b_i
    real, dimension(nu,nz) :: u_i
#endif /* ORIG */
    integer i,j,k,ip,jp,kp
    real xi,yj


    b_z = 0.0
#ifdef ORIG
    b_i = 0.0
#endif /* ORIG */
#ifdef COSM_RAYS
   call div_v         
#endif /* COSM_RAYS */

    do j=1,ny-1
      jp=j+1
      do i=1,nx-1
        if(magfield)then
          ip=i+1
          b_z(:,:)=b(:,i,j,:)*0.5
          b_z(ibx,:)=b_z(ibx,:)+b(ibx,ip,j,:)*0.5
          b_z(iby,:)=b_z(iby,:)+b(iby,i,jp,:)*0.5
          b_z(ibz,:)=b_z(ibz,:)+eoshift(b_z(ibz,:),shift=1,boundary=big)
          b_z((/ibz,iby,ibx/),:)=b_z(:,:)
#ifdef ORIG
          b_i(:,:)=bi(:,i,j,:)*0.5
          b_i(ibx,:)=b_i(ibx,:)+bi(ibx,ip,j,:)*0.5
          b_i(iby,:)=b_i(iby,:)+bi(iby,i,jp,:)*0.5
          b_i(ibz,:)=b_i(ibz,:)+eoshift(b_i(ibz,:),shift=1,boundary=big)
          b_i((/ibz,iby,ibx/),:)=b_i(:,:)
#endif /* ORIG */
        else
          b_z = 0.0
#ifdef ORIG 
          b_i = 0.0
#endif /* ORIG */
        endif

        u_z(iuswpz,:)=u(:,i,j,:)
#ifdef ORIG
        u_i(iuswpz,:)=ui(:,i,j,:) 
        call relaxing_tvd(u_i,u_z,b_i,b_z,'zsweep',i,j,dz,nz,dt)
#else /* ORIG */
        call relaxing_tvd(u_z,b_z,'zsweep',i,j,dz,nz,dt)
#endif /* ORIG */

        Lu(:,i,j,:)=Lu(:,i,j,:)+u_z(iuswpz,:)

      end do
    end do
    
  end subroutine fluidz

!==========================================================================================
end module fluids
