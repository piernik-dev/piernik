! $Id$
#include "piernik.def"

module fluids  ! unsplit fluids
  contains

  subroutine fluidx

    use start,  only : dimensions,magfield,dt
    use arrays, only : u,b,nx,ny,nz,nu,ks,ke,ibx,iby,ibz,iuswpx
    use grid,   only : dx
    use tv,     only : relaxing_tvd

#ifdef SPLIT
    use start, only  : istep,integration_order
    use tv,    only  : integrate
    use fluid_boundaries, only : compute_u_bnd
! below arrays are used in relaxing_tvd & integrate
    use arrays, only : ur0,ul0
#else /* SPLIT */
#ifdef ORIG
    use start,  only : istep
    use arrays, only : ui,bi
    use fluxes, only : mhdflux
! below arrays are used in relaxing_tvd & integrate
    use arrays, only : ur0,ul0
#endif /* ORIG */
    use arrays, only : Lu
#endif /* SPLIT */
#if defined COSM_RAYS || defined PRESS_GRAD_EXCH
    use func,  only : div_v
#ifndef SPLIT
   use func,   only : div_vx,divv
   use arrays, only : nfluid
#endif /* SPLIT */
#endif /* COSM_RAYS || PRESS_GRAD_EXCH */

    implicit none
    real, dimension(3,nx)  :: b_x
    real, dimension(nu,nx) :: u_x,L_x
#ifndef SPLIT
#ifdef ORIG
    real, dimension(nu,nx)  :: flx,cfr
#endif /* ORIG */
#endif /* SPLIT */

    integer :: j,k,jp,kp

    u_x = 0.0
    L_x = 0.0

#if defined SPLIT || defined ORIG
    allocate(ul0(nu,nx),ur0(nu,nx))
#endif /* SPLIT || ORIG */
#if defined COSM_RAYS || defined PRESS_GRAD_EXCH
#ifndef SPLIT
    allocate(divv(nfluid,nx))
#endif /* SPLIT */
    call div_v
#endif /* COSM_RAYS || PRESS_GRAD_EXCH */

    do k=ks,ke
      kp=k+1
      do j=1,ny-1
#if defined COSM_RAYS || defined PRESS_GRAD_EXCH
#ifndef SPLIT
        call div_vx(k,j)
#endif /* SPLIT */
#endif /* COSM_RAYS || PRESS_GRAD_EXCH */
#ifndef SPLIT
#ifdef ORIG
        if(istep /= 1) then
          if(magfield)then
            jp=j+1
            b_x=0.5*bi(:,:,j,k)
            b_x(ibx,1:nx-1)=b_x(ibx,1:nx-1)+b_x(ibx,2:nx);       b_x(ibx,nx) = b_x(ibx,nx-1)
            b_x(iby,:)=b_x(iby,:)+0.5*bi(iby,:,jp,k)
            if(dimensions .eq. '3d')  b_x(ibz,:)=b_x(ibz,:)+0.5*bi(ibz,:,j,kp)
          else
            b_x = 0.0
          endif
          u_x(iuswpx,:) = ui(:,:,j,k)
          call mhdflux(flx,cfr,u_x,b_x,nx)
          ul0 = (u_x-flx/cfr)*0.5
          ur0 = (u_x+flx/cfr)*0.5
        endif
#endif /* ORIG */
#endif /* SPLIT */
        if(magfield)then
          jp=j+1
          b_x=0.5*b(:,:,j,k)
          b_x(ibx,1:nx-1)=b_x(ibx,1:nx-1)+b_x(ibx,2:nx)
	  b_x(ibx,nx) = b_x(ibx,nx-1)
          b_x(iby,:)=b_x(iby,:)+0.5*b(iby,:,jp,k)
          if(dimensions .eq. '3d')  b_x(ibz,:)=b_x(ibz,:)+0.5*b(ibz,:,j,kp)
        else
          b_x = 0.0
        endif
        u_x(iuswpx,:) = u(:,:,j,k)
#ifdef SPLIT
        do istep=1,integration_order
          call relaxing_tvd(L_x,u_x,b_x,'xsweep',j,k,dx,nx,dt)
          call integrate(u_x,L_x,b_x,'xsweep',j,k,nx)
        enddo
        u(:,:,j,k)=u_x(iuswpx,:)
#else /* SPLIT */
        call relaxing_tvd(L_x,u_x,b_x,'xsweep',j,k,dx,nx,dt)
        Lu(:,:,j,k)=Lu(:,:,j,k)+L_x(iuswpx,:)
#endif /* SPLIT */
      end do
    end do

#if defined COSM_RAYS || defined PRESS_GRAD_EXCH
#ifndef SPLIT
    deallocate(divv)
#endif /* SPLIT */
#endif /* COSM_RAYS || PRESS_GRAD_EXCH */
#ifdef SPLIT
    deallocate(ul0,ur0)
    call compute_u_bnd
#else /* SPLIT */
#ifdef ORIG
    deallocate(ul0,ur0)
#endif /* ORIG */
#endif /* SPLIT */

  end subroutine fluidx

!------------------------------------------------------------------------------------------

  subroutine fluidy

    use start,  only : dimensions,magfield,dt
    use arrays, only : u,b,nx,ny,nz,nu,ks,ke,ibx,iby,ibz,iuswpy
    use grid,   only : dy
    use tv,     only : relaxing_tvd
#ifdef SPLIT
    use start,  only : istep,integration_order
    use tv,     only : integrate
    use fluid_boundaries, only : compute_u_bnd
! below arrays are used in relaxing_tvd & integrate
    use arrays, only : ur0,ul0
#else /* SPLIT */
#ifdef ORIG
    use start,  only : istep
    use arrays, only : ui,bi
    use fluxes, only : mhdflux
! below arrays are used in relaxing_tvd & integrate
    use arrays, only : ur0,ul0
#endif /* ORIG */
    use arrays, only : Lu
#endif /* SPLIT */
#if defined COSM_RAYS || defined PRESS_GRAD_EXCH
    use func,   only : div_v
#ifndef SPLIT
    use func,   only : div_vy,divv
    use arrays, only : nfluid
#endif /* SPLIT */
#endif /* COSM_RAYS || PRESS_GRAD_EXCH */

    implicit none
    real, dimension(3,ny)  :: b_y
    real, dimension(nu,ny) :: u_y,L_y
#ifndef SPLIT
#ifdef ORIG
    real, dimension(nu,ny)  :: flx,cfr
#endif /* ORIG */
#endif /* SPLIT */

    integer :: i,k,ip,kp

    u_y = 0.0
    L_y = 0.0

#if defined SPLIT || defined ORIG
    allocate(ul0(nu,ny),ur0(nu,ny))
#endif /* SPLIT || ORIG */
#if defined COSM_RAYS || defined PRESS_GRAD_EXCH
#ifndef SPLIT
   allocate(divv(nfluid,ny))
#endif /* SPLIT */
    call div_v
#endif /* COSM_RAYS || PRESS_GRAD_EXCH */

    do k=ks,ke
      kp=k+1
      do i=1,nx-1
#if defined COSM_RAYS || defined PRESS_GRAD_EXCH
#ifndef SPLIT
        call div_vy(k,i)
#endif /* SPLIT */
#endif /* COSM_RAYS || PRESS_GRAD_EXCH */
#ifndef SPLIT
#ifdef ORIG
        if(istep /= 1) then
          if(magfield)then
            ip=i+1
            b_y(:,:)=0.5*bi(:,i,:,k)
            b_y(ibx,:)=b_y(ibx,:)+0.5*bi(ibx,ip,:,k)
            b_y(iby,1:ny-1)=b_y(iby,1:ny-1)+b_y(iby,2:ny)
	    b_y(iby,ny) = b_y(iby,ny-1)
            if (dimensions .eq. '3d') b_y(ibz,:)=b_y(ibz,:)+0.5*bi(ibz,i,:,kp)
            b_y((/iby,ibx,ibz/),:)=b_y(:,:)
          else
            b_y = 0.0
          endif
          u_y(iuswpy,:) = ui(:,i,:,k)
          call mhdflux(flx,cfr,u_y,b_y,ny)
          ul0 = (u_y-flx/cfr)*0.5
          ur0 = (u_y+flx/cfr)*0.5
        endif
#endif /* ORIG */
#endif /* SPLIT */
        if(magfield)then
          ip=i+1
          b_y(:,:)=0.5*b(:,i,:,k)
          b_y(ibx,:)=b_y(ibx,:)+0.5*b(ibx,ip,:,k)
          b_y(iby,1:ny-1)=b_y(iby,1:ny-1)+b_y(iby,2:ny)
	  b_y(iby,ny) = b_y(iby,ny-1)
          if (dimensions .eq. '3d') b_y(ibz,:)=b_y(ibz,:)+0.5*b(ibz,i,:,kp)
          b_y((/iby,ibx,ibz/),:)=b_y(:,:)
        else
          b_y = 0.0
        endif
        u_y(iuswpy,:) = u(:,i,:,k)
#ifdef SPLIT
        do istep=1,integration_order
          call relaxing_tvd(L_y,u_y,b_y,'ysweep',k,i,dy,ny,dt)
          call integrate(u_y,L_y,b_y,'ysweep',k,i,ny)
        enddo
        u(:,i,:,k)=u_y(iuswpy,:)
#else /* SPLIT */
        call relaxing_tvd(L_y,u_y,b_y,'ysweep',k,i,dy,ny,dt)
        Lu(:,i,:,k)=Lu(:,i,:,k)+L_y(iuswpy,:)
#endif /* SPLIT */
      end do
    end do

#if defined COSM_RAYS || defined PRESS_GRAD_EXCH
#ifndef SPLIT
    deallocate(divv)
#endif /* SPLIT */
#endif /* COSM_RAYS || PRESS_GRAD_EXCH */
#ifdef SPLIT
    deallocate(ul0,ur0)
    call compute_u_bnd
#else /* SPLIT */
#ifdef ORIG
    deallocate(ul0,ur0)
#endif /* ORIG */
#endif /* SPLIT */

  end subroutine fluidy

!------------------------------------------------------------------------------------------

  subroutine fluidz

    use start,  only : dimensions,magfield,dt
    use arrays, only : u,b,nx,ny,nz,nu,ks,ke,ibx,iby,ibz,iuswpz
    use grid,   only : dz
    use tv,     only : relaxing_tvd
#ifdef SPLIT
    use start,  only : istep,integration_order
    use tv,     only : integrate
    use fluid_boundaries, only : compute_u_bnd
! below arrays are used in relaxing_tvd & integrate
    use arrays, only : ur0,ul0
#else /* SPLIT */
#ifdef ORIG
    use start,  only : istep
    use arrays, only : ui,bi
    use fluxes, only : mhdflux
! below arrays are used in relaxing_tvd & integrate
    use arrays, only : ur0,ul0
#endif /* ORIG */
    use arrays, only : Lu
#endif /* SPLIT */
#if defined COSM_RAYS || defined PRESS_GRAD_EXCH
    use func,   only : div_v
#ifndef SPLIT
    use func,   only : div_vz,divv
    use arrays, only : nfluid
#endif /* SPLIT */
#endif /* COSM_RAYS || PRESS_GRAD_EXCH */

    implicit none
    real, dimension(3,nz)  :: b_z
    real, dimension(nu,nz) :: u_z,L_z
#ifndef SPLIT
#ifdef ORIG
    real, dimension(nu,nz) :: flx,cfr
#endif /* ORIG */
#endif /* SPLIT */

    integer :: i,j,ip,jp

    u_z = 0.0
    L_z = 0.0

#if defined SPLIT || defined ORIG
    allocate(ul0(nu,nz),ur0(nu,nz))
#endif /* SPLIT || ORIG */
#if defined COSM_RAYS || defined PRESS_GRAD_EXCH
#ifndef SPLIT
    allocate(divv(nfluid,nz))
#endif /* SPLIT */
    call div_v
#endif /* COSM_RAYS || PRESS_GRAD_EXCH */

    do j=1,ny-1
      jp=j+1
      do i=1,nx-1
#if defined COSM_RAYS || defined PRESS_GRAD_EXCH
#ifndef SPLIT
        call div_vz(j,i)
#endif /* SPLIT */
#endif /* COSM_RAYS || PRESS_GRAD_EXCH */
#ifndef SPLIT
#ifdef ORIG
        if(istep /= 1) then
          if(magfield)then
            ip=i+1
            b_z(:,:)=0.5*bi(:,i,j,:)
            b_z(ibx,:)=b_z(ibx,:)+0.5*bi(ibx,ip,j,:)
            b_z(iby,:)=b_z(iby,:)+0.5*bi(iby,i,jp,:)
            b_z(ibz,1:nz-1)=b_z(ibz,1:nz-1)+b_z(ibz,2:nz);       b_z(ibz,nz) = b_z(ibz,nz-1)
            b_z((/ibz,iby,ibx/),:)=b_z(:,:)
          else
            b_z = 0.0
          endif
          u_z(iuswpz,:)=ui(:,i,j,:)
          call mhdflux(flx,cfr,u_z,b_z,nz)
          ul0 = (u_z-flx/cfr)*0.5
          ur0 = (u_z+flx/cfr)*0.5
        endif
#endif /* ORIG */
#endif /* SPLIT */
        if(magfield)then
          ip=i+1
          b_z(:,:)=0.5*b(:,i,j,:)
          b_z(ibx,:)=b_z(ibx,:)+0.5*b(ibx,ip,j,:)
          b_z(iby,:)=b_z(iby,:)+0.5*b(iby,i,jp,:)
          b_z(ibz,1:nz-1)=b_z(ibz,1:nz-1)+b_z(ibz,2:nz)
	  b_z(ibz,nz) = b_z(ibz,nz-1)
          b_z((/ibz,iby,ibx/),:)=b_z(:,:)
        else
          b_z = 0.0
        endif

        u_z(iuswpz,:)=u(:,i,j,:)
#ifdef SPLIT
        do istep=1,integration_order
          call relaxing_tvd(L_z,u_z,b_z,'zsweep',i,j,dz,nz,dt)
          call integrate(u_z,L_z,b_z,'zsweep',i,j,nz)
        enddo
        u(:,i,j,:)=u_z(iuswpz,:)
#else /* SPLIT */
        call relaxing_tvd(L_z,u_z,b_z,'zsweep',i,j,dz,nz,dt)
        Lu(:,i,j,:)=Lu(:,i,j,:)+L_z(iuswpz,:)
#endif /* SPLIT */
      end do
    end do

#if defined COSM_RAYS || defined PRESS_GRAD_EXCH
#ifndef SPLIT
    deallocate(divv)
#endif /* SPLIT */
#endif /* COSM_RAYS || PRESS_GRAD_EXCH */
#ifdef SPLIT
    deallocate(ul0,ur0)
    call compute_u_bnd
#else /* SPLIT */
#ifdef ORIG
    deallocate(ul0,ur0)
#endif /* ORIG */
#endif /* SPLIT */

  end subroutine fluidz
end module fluids
