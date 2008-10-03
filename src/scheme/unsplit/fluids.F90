! $Id$
#include "piernik.def"
module fluids  ! unsplit fluids
!  use grid
!  use start
!  use arrays
!  use tv
  contains
!==========================================================================================

  subroutine fluidx
    use start, only     : dimensions, magfield,dt
    use arrays, only    : u,b,nx,ny,nz,nu,ks,ke,ibx,iby,ibz,iuswpx,Lu
    use grid, only      : dx
    use tv, only        : relaxing_tvd
    use constants, only : big
    use fluxes, only    : mhdflux

#ifdef FLX_BND
    use arrays, only  : cfr,flx
    use start, only   : nb,nxd,nyd,nzd
#ifdef SHEAR
    use start, only   : omega,qshear,xmax,xmin
    use arrays, only  : x,xl,xr
    use shear, only   : unshear
#endif /* SHEAR */
#endif /* FLX_BND */
    implicit none
    real, dimension(3,nx)  :: b_x
    real, dimension(nu,nx) :: u_x
#ifdef FLX_BND
    real, dimension(nu,nx) :: f_x,c_x
#endif /* FLX_BND */
    
#ifdef ORIG
    real, dimension(3,nx)  :: b_i
    real, dimension(nu,nx) :: u_i
#endif /* ORIG */

    integer :: j,k,jp,kp
    real    :: yj,zk,S

#ifdef SHEAR
!    S = qshear*omega*(xmax-xmin)
#endif /* SHEAR */

    b_x = 0.0
#ifdef ORIG
    b_i = 0.0
#endif /* ORIG */

#ifdef COSM_RAYS
   call div_v         
#endif /* COSM_RAYS */


#ifdef FLX_BND
   do k=ks,ke
     kp = k + 1
     do j = 1, ny-1
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

        u_x(iuswpx,:)= u(:,:,j,k)
        call mhdflux(f_x,c_x,u_x,b_x,nx)
        flx(:,:,j,k) = f_x
        cfr(:,j,k)   = c_x(1,:)
      enddo
    enddo

    do k=ks,ke
      do j=1,ny
        flx(3,:,j,k) = flx(3,:,j,k) + qshear*omega*x(:)*flx(1,:,j,k)
#ifndef ISO
        flx(5,:,j,k) = flx(5,:,j,k) + 0.5*(qshear*omega*x(:))**2*flx(1,:,j,k) &
            + flx(3,:,j,k)*qshear*omega*x(:) 
#endif /* ISO */
      enddo
    enddo

! UNSHEAR

    do j=LBOUND(flx,1),UBOUND(flx,1)
      flx(j,:,:,:) = unshear(flx(j,:,:,:), x(:),.true.)
    enddo

    flx(:,1:nb,:,:)              = flx(:,nxd+1:nxd+nb,:,:)
    flx(:,nxd+nb+1:nxd+2*nb,:,:) = flx(:,nb+1:2*nb,:,:)
    flx(:,:,1:nb,:)              = flx(:,:,nyd+1:nyd+nb,:)
    flx(:,:,nyd+nb+1:nyd+2*nb,:) = flx(:,:,nb+1:2*nb,:)
    if(dimensions == '3d') then
       flx(:,:,:,1:nb)              = flx(:,:,:,nzd+1:nzd+nb)
       flx(:,:,:,nzd+nb+1:nzd+2*nb) = flx(:,:,:,nb+1:2*nb)
    endif

    do j=LBOUND(flx,1),UBOUND(flx,1)
      flx(j,:,:,:) = unshear(flx(j,:,:,:), x(:))
    enddo
! SHEAR

    do k=ks,ke
      do j=1,ny
        flx(3,:,j,k) = flx(3,:,j,k) - qshear*omega*xl(:)*flx(1,:,j,k)
#ifndef ISO
        flx(5,:,j,k) = flx(5,:,j,k) + 0.5*(qshear*omega*xl(:))**2*flx(1,:,j,k) &
          - qshear*omega*xl(:)*flx(3,:,j,k)
#endif /* ISO */
      enddo
    enddo

#endif /* FLX_BND */

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

        u_x(iuswpx,:) = u(:,:,j,k)
#ifdef FLX_BND 
        c_x           = spread(cfr(:,j,k),1,5)
        f_x(:,:)      = flx(:,:,j,k)
#endif /* FLX_BND */
#ifdef ORIG
        u_i(iuswpx,:)=ui(:,:,j,k)
        call relaxing_tvd(u_i,u_x,b_i,b_x,'xsweep',j,k,dx,nx,dt)
#else /* ORIG */
#ifdef FLX_BND
        call relaxing_tvd(u_x,b_x,'xsweep',j,k,dx,nx,dt,f_x,c_x)
#else /* ~FLX_BND */
        call relaxing_tvd(u_x,b_x,'xsweep',j,k,dx,nx,dt)
#endif /* ~FLX_BND */
#endif /* ORIG */
        Lu(:,:,j,k)=Lu(:,:,j,k)+u_x(iuswpx,:)
      end do
    end do


  end subroutine fluidx

!------------------------------------------------------------------------------------------

  subroutine fluidy
    use start, only     : dimensions, magfield,dt
    use arrays, only    : u,b,nx,ny,nz,nu,ks,ke,ibx,iby,ibz,iuswpy,Lu
    use grid, only      : dy
    use tv, only        : relaxing_tvd
    use constants, only : big
    use fluxes, only    : mhdflux

#ifdef FLX_BND
    use arrays, only  : cfr,flx
#endif /* FLX_BND */
    use start, only   : nb,nxd,nyd,nzd
#ifdef SHEAR
    use start, only   : omega,qshear
    use arrays, only  : x
    use shear, only   : unshear
#endif /* SHEAR */

    implicit none
    real, dimension(3,ny)  :: b_y
    real, dimension(nu,ny) :: u_y,f_y,c_y
    
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

#undef FLX_BND_Y
#ifdef FLX_BND_Y
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
        call mhdflux(f_y,c_y,u_y,b_y,ny)
        flx(:,i,:,k) = f_y
        cfr(i,:,k)   = c_y(1,:)

      end do
    end do

    do k=ks,ke
      do j=1,ny
         flx(1,:,j,k) = flx(1,:,j,k) + qshear*omega*x(:)*u(1,:,j,k)
         flx(2,:,j,k) = flx(1,:,j,k)**2/u(1,:,j,k) + &
            2*flx(1,:,j,k)*qshear*omega*x(:) + &
            (qshear*omega*x(:))**2*u(1,:,j,k)
         flx(3,:,j,k) = flx(3,:,j,k) + flx(3,:,j,k)*(qshear*omega*x(:))
         flx(4,:,j,k) = flx(4,:,j,k) + flx(4,:,j,k)*(qshear*omega*x(:))
      enddo
    enddo

    do j=LBOUND(flx,1),UBOUND(flx,1)
      flx(j,:,:,:) = unshear(flx(j,:,:,:), x(:),.true.)
    enddo

    flx(:,1:nb,:,:)              = flx(:,nxd+1:nxd+nb,:,:)
    flx(:,nxd+nb+1:nxd+2*nb,:,:) = flx(:,nb+1:2*nb,:,:)
    flx(:,:,1:nb,:)              = flx(:,:,nyd+1:nyd+nb,:)
    flx(:,:,nyd+nb+1:nyd+2*nb,:) = flx(:,:,nb+1:2*nb,:)
    if(dimensions == '3d') then
       flx(:,:,:,1:nb)              = flx(:,:,:,nzd+1:nzd+nb)
       flx(:,:,:,nzd+nb+1:nzd+2*nb) = flx(:,:,:,nb+1:2*nb)
    endif

    do j=LBOUND(flx,1),UBOUND(flx,1)
      flx(j,:,:,:) = unshear(flx(j,:,:,:), x(:))
    enddo
    do k=ks,ke
      do j=1,ny
         flx(1,:,j,k) = flx(1,:,j,k) - qshear*omega*x(:)*u(1,:,j,k)
         flx(2,:,j,k) = flx(1,:,j,k)**2/u(1,:,j,k) - &
            2*flx(1,:,j,k)*qshear*omega*x(:) + &
            (qshear*omega*x(:))**2*u(1,:,j,k)
         flx(3,:,j,k) = flx(3,:,j,k) - flx(3,:,j,k)*(qshear*omega*x(:))
         flx(4,:,j,k) = flx(4,:,j,k) - flx(4,:,j,k)*(qshear*omega*x(:))
      enddo
    enddo
#endif /* FLX_BND_Y */
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
#ifdef FLX_BND_Y
        call relaxing_tvd(u_y,b_y,'ysweep',k,i,dy,ny,dt,f_y,c_y)
#else /* ~FLX_BND */
        call relaxing_tvd(u_y,b_y,'ysweep',k,i,dy,ny,dt)
#endif /* ~FLX_BND */
#endif /* ORIG */
        Lu(:,i,:,k)=Lu(:,i,:,k)+u_y(iuswpy,:)

      end do
    end do

  end subroutine fluidy

!------------------------------------------------------------------------------------------

  subroutine fluidz
    use start, only     : dimensions, magfield,dt
    use arrays, only    : u,b,nx,ny,nz,nu,ks,ke,ibx,iby,ibz,iuswpz,Lu
    use grid, only      : dz
    use tv, only        : relaxing_tvd
    use constants, only : big
    use fluxes, only    : mhdflux
#ifdef FLX_BND
    use arrays, only  : cfr,flx
#endif /* FLX_BND */
    use start, only   : nb,nxd,nyd,nzd

    implicit none
    real, dimension(3,nz)  :: b_z
    real, dimension(nu,nz) :: u_z,f_z,c_z
    
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

#ifdef FLX_BND
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
        call mhdflux(f_z,c_z,u_z,b_z,nz)
        flx(:,i,j,:) = f_z
        cfr(i,j,:)   = c_z(1,:)

      end do
    end do
    flx(:,1:nb,:,:)              = flx(:,nxd+1:nxd+nb,:,:)
    flx(:,nxd+nb+1:nxd+2*nb,:,:) = flx(:,nb+1:2*nb,:,:)
    flx(:,:,1:nb,:)              = flx(:,:,nyd+1:nyd+nb,:)
    flx(:,:,nyd+nb+1:nyd+2*nb,:) = flx(:,:,nb+1:2*nb,:)
    if(dimensions == '3d') then
       flx(:,:,:,1:nb)              = flx(:,:,:,nzd+1:nzd+nb)
       flx(:,:,:,nzd+nb+1:nzd+2*nb) = flx(:,:,:,nb+1:2*nb)
    endif
#endif /* FLX_BND */
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
