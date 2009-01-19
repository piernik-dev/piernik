! $Id$
#include "piernik.def"

module advects    ! both advects
  contains

  subroutine advectby_x
    use start,  only : dimensions,dt
    use fluidindex,   only : nfluid
    use initionized, only : ibx,iby,ibz
    use initionized, only : idni,imxi
    use arrays, only : b,u,wa,nx,ny,nz
    use grid,   only : dx
#ifndef SPLIT
    use arrays, only : Lb,xdim,ydim
    use grid,   only : dy
#endif /* SPLIT */
    use rtvd,     only : tvdb
    use mag_boundaries, only : bnd_emf

    implicit none
    real, dimension(nx) :: vxby,by_x,vx
    integer                j,k,jm,ifluid

    vxby = 0.0

    do k=1,nz
      do j=2,ny
        jm=j-1
	vx=0.0
        vx=(u(imxi,:,jm,k)+u(imxi,:,j,k))/(u(idni,:,jm,k)+u(idni,:,j,k))
        vx(2:nx-1)=(vx(1:nx-2) + vx(3:nx) + 2.0*vx(2:nx-1))*0.25
        vx(1)  = vx(2)
        vx(nx) = vx(nx-1)
        by_x=b(iby,:,j,k)
#if defined ORIG && defined SPLIT
        call tvdb(vxby,by_x,vx,nx,dt,dx)
#else /* ORIG && SPLIT */
        call tvdb(vxby,by_x,vx,nx,dt)
#endif /* ORIG && SPLIT */
        wa(:,j,k) = vxby
      end do
    end do

    call bnd_emf(wa, 'vxby', 'xdim')
    call bnd_emf(wa, 'vxby', 'ydim')
    if(dimensions .eq. '3d') call bnd_emf(wa, 'vxby', 'zdim')

#ifndef SPLIT
    Lb(iby,:,:,:) = Lb(iby,:,:,:) - wa/dx
    wa = cshift(wa,shift=-1,dim=xdim)
    Lb(iby,:,:,:) = Lb(iby,:,:,:) + wa/dx
    Lb(ibx,:,:,:) = Lb(ibx,:,:,:) - wa/dy
    wa = cshift(wa,shift=1,dim=ydim)
    Lb(ibx,:,:,:) = Lb(ibx,:,:,:) + wa/dy

#endif /* SPLIT */
  end subroutine advectby_x

  subroutine advectbz_x
    use start,  only : dimensions,dt
    use fluidindex,   only : nfluid
    use initionized, only : ibx,iby,ibz
    use initionized, only : idni,imxi
    use arrays, only : b,u,wa,nx,ny,nz
    use grid,   only : dx
    use rtvd,   only : tvdb
#ifndef SPLIT
    use arrays, only : Lb,xdim,zdim
    use grid,   only : dz
#endif /* SPLIT */
    use mag_boundaries, only : bnd_emf

    implicit none
    real, dimension(nx) :: vxbz,bz_x,vx
    integer                j,k,km,ifluid

    vxbz = 0.0

    do k=2,nz
      km=k-1
      do j=1,ny
        vx=0.0
        vx=(u(imxi,:,j,km)+u(imxi,:,j,k))/(u(idni,:,j,km)+u(idni,:,j,k))
        vx(2:nx-1)=(vx(1:nx-2) + vx(3:nx) + 2.0*vx(2:nx-1))*0.25
        vx(1)  = vx(2)
        vx(nx) = vx(nx-1)
        bz_x=b(ibz,:,j,k)
#if defined ORIG && defined SPLIT
        call tvdb(vxbz,bz_x,vx,nx,dt,dx)
#else /* ORIG && SPLIT */
        call tvdb(vxbz,bz_x,vx,nx,dt)
#endif /* ORIG && SPLIT */
        wa(:,j,k) = vxbz
      end do
    end do

    call bnd_emf(wa, 'vxbz', 'xdim')
    call bnd_emf(wa, 'vxbz', 'ydim')
    if(dimensions .eq. '3d') call bnd_emf(wa, 'vxbz', 'zdim')

#ifndef SPLIT
    Lb(ibz,:,:,:) = Lb(ibz,:,:,:) - wa/dx
    wa = cshift(wa,shift=-1,dim=xdim)
    Lb(ibz,:,:,:) = Lb(ibz,:,:,:) + wa/dx
    Lb(ibx,:,:,:) = Lb(ibx,:,:,:) - wa/dz
    wa = cshift(wa,shift=1,dim=zdim)
    Lb(ibx,:,:,:) = Lb(ibx,:,:,:) + wa/dz

#endif /* SPLIT */
  end subroutine advectbz_x

  subroutine advectbz_y
    use start,  only : dimensions,dt
    use fluidindex,   only : nfluid
    use initionized, only : ibx,iby,ibz
    use initionized, only : idni,imyi
    use arrays, only : b,u,wa,nx,ny,nz
    use grid,   only : dy
    use rtvd,   only : tvdb
#ifndef SPLIT
    use arrays, only : Lb,ydim,zdim
    use grid,   only : dz
#endif /* SPLIT */
    use mag_boundaries, only : bnd_emf

    implicit none
    real, dimension(ny)   :: vybz,bz_y,vy
    integer                  i,k,km,ifluid

    vybz = 0.0

    do k=2,nz
      km=k-1
      do i=1,nx
        vy=0.0
        vy=vy+(u(imyi,i,:,km)+u(imyi,i,:,k))/(u(idni,i,:,km)+u(idni,i,:,k))
        vy(2:ny-1)=(vy(1:ny-2) + vy(3:ny) + 2.0*vy(2:ny-1))*0.25
        vy(1)  = vy(2)
        vy(ny) = vy(ny-1)
        bz_y=b(ibz,i,:,k)
#if defined ORIG && defined SPLIT
        call tvdb(vybz,bz_y,vy,ny,dt,dy)
#else /* ORIG && SPLIT */
        call tvdb(vybz,bz_y,vy,ny,dt)
#endif /* ORIG && SPLIT */
        wa(i,:,k) = vybz
      end do
    end do

    call bnd_emf(wa, 'vybz', 'ydim')
    if(dimensions .eq. '3d') call bnd_emf(wa, 'vybz', 'zdim')
    call bnd_emf(wa, 'vybz', 'xdim')

#ifndef SPLIT
    Lb(ibz,:,:,:) = Lb(ibz,:,:,:) - wa/dy
    wa = cshift(wa,shift=-1,dim=ydim)
    Lb(ibz,:,:,:) = Lb(ibz,:,:,:) + wa/dy
    Lb(iby,:,:,:) = Lb(iby,:,:,:) - wa/dz
    wa = cshift(wa,shift=1,dim=zdim)
    Lb(iby,:,:,:) = Lb(iby,:,:,:) + wa/dz

#endif /* SPLIT */
  end subroutine advectbz_y

  subroutine advectbx_y
    use start,  only : dimensions,dt
    use fluidindex,   only : nfluid
    use initionized, only : ibx,iby,ibz
    use initionized, only : idni,imyi
    use arrays, only : b,u,wa,nx,ny,nz
    use grid,   only : dy
    use rtvd,     only : tvdb
#ifndef SPLIT
    use arrays, only : Lb,ydim,xdim
    use grid,   only : dx
#endif /* SPLIT */
    use mag_boundaries, only : bnd_emf

    implicit none
    real, dimension(ny) :: vybx,bx_y,vy
    integer                k,i,im,ifluid

    vybx = 0.0

    do k=1,nz
      do i=2,nx
        im=i-1
        vy=0.0
	vy=(u(imyi,im,:,k)+u(imyi,i,:,k))/(u(idni,im,:,k)+u(idni,i,:,k))
        vy(2:ny-1)=(vy(1:ny-2) + vy(3:ny) + 2.0*vy(2:ny-1))*0.25
        vy(1)  = vy(2)
        vy(ny) = vy(ny-1)
        bx_y=b(ibx,i,:,k)
#if defined ORIG && defined SPLIT
        call tvdb(vybx,bx_y,vy,ny,dt,dy)
#else /* ORIG && SPLIT */
        call tvdb(vybx,bx_y,vy,ny,dt)
#endif /* ORIG && SPLIT */
        wa(i,:,k) = vybx
      end do
    end do

    call bnd_emf(wa, 'vybx', 'ydim')
    if(dimensions .eq. '3d') call bnd_emf(wa, 'vybx', 'zdim')
    call bnd_emf(wa, 'vybx', 'xdim')

#ifndef SPLIT
    Lb(ibx,:,:,:) = Lb(ibx,:,:,:) - wa/dy
    wa = cshift(wa,shift=-1,dim=ydim)
    Lb(ibx,:,:,:) = Lb(ibx,:,:,:) + wa/dy
    Lb(iby,:,:,:) = Lb(iby,:,:,:) - wa/dx
    wa = cshift(wa,shift=1,dim=xdim)
    Lb(iby,:,:,:) = Lb(iby,:,:,:) + wa/dx

#endif /* SPLIT */
  end subroutine advectbx_y

  subroutine advectbx_z
    use start,  only : dimensions,dt
    use fluidindex,   only : nfluid
    use initionized, only : ibx,iby,ibz
    use initionized, only : idni,imzi
    use arrays, only : b,u,wa,nx,ny,nz
    use grid,   only : dz
    use rtvd,     only : tvdb
#ifndef SPLIT
    use arrays, only : Lb,zdim,xdim
    use grid,   only : dx
#endif /* SPLIT */
    use mag_boundaries, only : bnd_emf

    implicit none
    real, dimension(nz)  :: vzbx,bx_z,vz
    integer                 j,i,im,ifluid

    vzbx = 0.0

    do j=1,ny
      do i=2,nx
        im=i-1
	vz=0.0
        vz=(u(imzi,im,j,:)+u(imzi,i,j,:))/(u(idni,im,j,:)+u(idni,i,j,:))
        vz(2:nz-1)=(vz(1:nz-2) + vz(3:nz) + 2.0*vz(2:nz-1))*0.25
        vz(1)  = vz(2)
        vz(nz) = vz(nz-1)
        bx_z=b(ibx,i,j,:)
#if defined ORIG && defined SPLIT
        call tvdb(vzbx,bx_z,vz,nz,dt,dz)
#else /* ORIG && SPLIT */
        call tvdb(vzbx,bx_z,vz,nz,dt)
#endif /* ORIG && SPLIT */
        wa(i,j,:) = vzbx
      end do
    end do

    if(dimensions .eq. '3d') call bnd_emf(wa, 'vzbx', 'zdim')
    call bnd_emf(wa, 'vzbx', 'xdim')
    call bnd_emf(wa, 'vzbx', 'ydim')

#ifndef SPLIT
    Lb(ibx,:,:,:) = Lb(ibx,:,:,:) - wa/dz
    wa = cshift(wa,shift=-1,dim=zdim)
    Lb(ibx,:,:,:) = Lb(ibx,:,:,:) + wa/dz
    Lb(ibz,:,:,:) = Lb(ibz,:,:,:) - wa/dx
    wa = cshift(wa,shift=1,dim=xdim)
    Lb(ibz,:,:,:) = Lb(ibz,:,:,:) + wa/dx

#endif /* SPLIT */
  end subroutine advectbx_z

  subroutine advectby_z
    use start,  only : dimensions,dt
    use fluidindex,   only : nfluid
    use initionized, only : ibx,iby,ibz
    use initionized, only : idni,imzi
    use arrays, only : b,u,wa,nx,ny,nz
    use grid,   only : dz
    use rtvd,     only : tvdb
#ifndef SPLIT
    use arrays, only : Lb,zdim,ydim
    use grid,   only : dy
#endif /* SPLIT */
    use mag_boundaries, only : bnd_emf

    implicit none
    real, dimension(nz)  :: vzby,by_z,vz
    integer                 i,j,jm,ifluid

    vzby = 0.0

    do j=2,ny
      jm=j-1
      do i=1,nx
        vz=0.0
	vz=(u(imzi,i,jm,:)+u(imzi,i,j,:))/(u(idni,i,jm,:)+u(idni,i,j,:))
        vz(2:nz-1)=(vz(1:nz-2) + vz(3:nz) + 2.0*vz(2:nz-1))*0.25
        vz(1)  = vz(2)
        vz(nz) = vz(nz-1)
        by_z=b(iby,i,j,:)
#if defined ORIG && defined SPLIT
        call tvdb(vzby,by_z,vz,nz,dt,dz)
#else /* ORIG && SPLIT */
        call tvdb(vzby,by_z,vz,nz,dt)
#endif /* ORIG && SPLIT */
        wa(i,j,:) = vzby
      end do
    end do

    if(dimensions .eq. '3d') call bnd_emf(wa, 'vzby', 'zdim')
    call bnd_emf(wa, 'vzby', 'xdim')
    call bnd_emf(wa, 'vzby', 'ydim')

#ifndef SPLIT
    Lb(iby,:,:,:) = Lb(iby,:,:,:) - wa/dz
    wa = cshift(wa,shift=-1,dim=zdim)
    Lb(iby,:,:,:) = Lb(iby,:,:,:) + wa/dz
    Lb(ibz,:,:,:) = Lb(ibz,:,:,:) - wa/dy
    wa = cshift(wa,shift=1,dim=ydim)
    Lb(ibz,:,:,:) = Lb(ibz,:,:,:) + wa/dy

#endif /* SPLIT */
  end subroutine advectby_z

end module advects
