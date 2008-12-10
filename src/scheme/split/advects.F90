! $Id$
#include "piernik.def"
module advects    ! split advects
  contains
  subroutine advectby_x
    use start, only  : dimensions,dt
    use arrays, only : b,u,wa,nx,ny,nz,idna,imxa,iby,magn,nfluid
    use grid, only   : dx
    use tv, only     : tvdb
    use mag_boundaries, only : bnd_emf
    implicit none

    real, dimension(nx) :: vxby,by_x,vx
    integer                j,k,jm,ifluid

    vxby = 0.0
    
    do k=1,nz
      do j=2,ny
        jm=j-1
	vx=0.0
	do ifluid=1,nfluid
        if(magn(ifluid) .eq. 1) vx=(u(imxa(ifluid),:,jm,k)+u(imxa(ifluid),:,j,k))/(u(idna(ifluid),:,jm,k)+u(idna(ifluid),:,j,k))
	enddo
        vx(2:nx-1)=(vx(1:nx-2) + vx(3:nx) + 2.0*vx(2:nx-1))*0.25
        vx(1)  = vx(2)
        vx(nx) = vx(nx-1)
        by_x=b(iby,:,j,k)
#ifdef ORIG
        call tvdb(vxby,by_x,vx,nx,dt,dx)
#else /* ORIG */
        call tvdb(vxby,by_x,vx,nx,dt)
#endif /* ORIG */
        wa(:,j,k) = vxby
      end do
    end do

    call bnd_emf(wa, 'vxby', 'xdim')
    call bnd_emf(wa, 'vxby', 'ydim')
    if(dimensions .eq. '3d') call bnd_emf(wa, 'vxby', 'zdim')

  end subroutine advectby_x

  subroutine advectbz_x
    use start, only  : dimensions,dt
    use arrays, only : b,u,wa,nx,ny,nz,idna,imxa,ibz,magn,nfluid
    use grid, only   : dx
    use tv, only     : tvdb
    use mag_boundaries, only : bnd_emf

    implicit none
    real, dimension(nx) :: vxbz,bz_x,vx
    integer                j,k,km,ifluid

    vxbz = 0.0

    do k=2,nz
      km=k-1
      do j=1,ny
        vx=0.0
	do ifluid=1,nfluid
        if(magn(ifluid) .eq. 1) vx=(u(imxa(ifluid),:,j,km)+u(imxa(ifluid),:,j,k))/(u(idna(ifluid),:,j,km)+u(idna(ifluid),:,j,k))
	enddo
        vx(2:nx-1)=(vx(1:nx-2) + vx(3:nx) + 2.0*vx(2:nx-1))*0.25
        vx(1)  = vx(2)
        vx(nx) = vx(nx-1)
        bz_x=b(ibz,:,j,k)
#ifdef ORIG
        call tvdb(vxbz,bz_x,vx,nx,dt,dx)
#else /* ORIG */
        call tvdb(vxbz,bz_x,vx,nx,dt)
#endif /* ORIG */
        wa(:,j,k) = vxbz
      end do
    end do

    call bnd_emf(wa, 'vxbz', 'xdim')
    call bnd_emf(wa, 'vxbz', 'ydim')
    if(dimensions .eq. '3d') call bnd_emf(wa, 'vxbz', 'zdim')

  end subroutine advectbz_x

  subroutine advectbz_y
    use start, only  : dimensions,dt
    use arrays, only : b,u,wa,nx,ny,nz,idna,imya,ibz,magn,nfluid
    use grid, only   : dy
    use tv, only     : tvdb
    use mag_boundaries, only : bnd_emf

    implicit none
    real, dimension(ny)   :: vybz,bz_y,vy
    integer                  i,k,km,ifluid

    vybz = 0.0

    do k=2,nz
      km=k-1
      do i=1,nx  
        vy=0.0
        do ifluid=1,nfluid
          if(magn(ifluid) .eq. 1)  then 
              vy=vy+(u(imya(ifluid),i,:,km)+u(imya(ifluid),i,:,k))/(u(idna(ifluid),i,:,km)+u(idna(ifluid),i,:,k))
          endif
        enddo
        vy(2:ny-1)=(vy(1:ny-2) + vy(3:ny) + 2.0*vy(2:ny-1))*0.25
        vy(1)  = vy(2)
        vy(ny) = vy(ny-1)
        bz_y=b(ibz,i,:,k)
#ifdef ORIG
        call tvdb(vybz,bz_y,vy,ny,dt,dy)
#else /* ORIG */
        call tvdb(vybz,bz_y,vy,ny,dt)
#endif /* ORIG */
        wa(i,:,k) = vybz
      end do
    end do

    call bnd_emf(wa, 'vybz', 'ydim')
    if(dimensions .eq. '3d') call bnd_emf(wa, 'vybz', 'zdim')
    call bnd_emf(wa, 'vybz', 'xdim')
  
  end subroutine advectbz_y

  subroutine advectbx_y
    use start, only  : dimensions,dt
    use arrays, only : b,u,wa,nx,ny,nz,idna,imya,ibx,magn,nfluid
    use grid, only   : dy
    use tv, only     : tvdb
    use mag_boundaries, only : bnd_emf

    implicit none
    real, dimension(ny) :: vybx,bx_y,vy
    integer                k,i,im,ifluid

    vybx = 0.0

    do k=1,nz
      do i=2,nx
        im=i-1
        vy=0.0
	do ifluid=1,nfluid
	if(magn(ifluid) .eq. 1) vy=(u(imya(ifluid),im,:,k)+u(imya(ifluid),i,:,k))/(u(idna(ifluid),im,:,k)+u(idna(ifluid),i,:,k))
	enddo
        vy(2:ny-1)=(vy(1:ny-2) + vy(3:ny) + 2.0*vy(2:ny-1))*0.25
        vy(1)  = vy(2)
        vy(ny) = vy(ny-1)
        bx_y=b(ibx,i,:,k)
#ifdef ORIG
        call tvdb(vybx,bx_y,vy,ny,dt,dy)
#else /* ORIG */
        call tvdb(vybx,bx_y,vy,ny,dt)
#endif /* ORIG */
        wa(i,:,k) = vybx
      end do
    end do

    call bnd_emf(wa, 'vybx', 'ydim')
    if(dimensions .eq. '3d') call bnd_emf(wa, 'vybx', 'zdim')
    call bnd_emf(wa, 'vybx', 'xdim')

  end subroutine advectbx_y

  subroutine advectbx_z
    use start, only  : dimensions,dt
    use arrays, only : b,u,wa,nx,ny,nz,idna,imza,ibx,magn,nfluid
    use grid, only   : dz
    use tv, only     : tvdb
    use mag_boundaries, only : bnd_emf

    implicit none
    real, dimension(nz)  :: vzbx,bx_z,vz
    integer                 j,i,im,ifluid

    vzbx = 0.0

    do j=1,ny
      do i=2,nx
        im=i-1
	vz=0.0
	do ifluid=1,nfluid
        if(magn(ifluid) .eq. 1) vz=(u(imza(ifluid),im,j,:)+u(imza(ifluid),i,j,:))/(u(idna(ifluid),im,j,:)+u(idna(ifluid),i,j,:))
	enddo
        vz(2:nz-1)=(vz(1:nz-2) + vz(3:nz) + 2.0*vz(2:nz-1))*0.25
        vz(1)  = vz(2)
        vz(nz) = vz(nz-1)
        bx_z=b(ibx,i,j,:)
#ifdef ORIG
        call tvdb(vzbx,bx_z,vz,nz,dt,dz)
#else /* ORIG */
        call tvdb(vzbx,bx_z,vz,nz,dt)
#endif /* ORIG */
        wa(i,j,:) = vzbx
      end do
    end do
    
    if(dimensions .eq. '3d') call bnd_emf(wa, 'vzbx', 'zdim')
    call bnd_emf(wa, 'vzbx', 'xdim')
    call bnd_emf(wa, 'vzbx', 'ydim')

  end subroutine advectbx_z

  subroutine advectby_z
    use start, only  : dimensions,dt
    use arrays, only : b,u,wa,nx,ny,nz,idna,imza,iby,magn,nfluid
    use grid, only   : dz
    use tv, only     : tvdb
    use mag_boundaries, only : bnd_emf
    
    implicit none
    real, dimension(nz)  :: vzby,by_z,vz
    integer                 i,j,jm,ifluid

    vzby = 0.0

    do j=2,ny
      jm=j-1
      do i=1,nx
        vz=0.0
	do ifluid=1,nfluid
	if(magn(ifluid) .eq. 1) vz=(u(imza(ifluid),i,jm,:)+u(imza(ifluid),i,j,:))/(u(idna(ifluid),i,jm,:)+u(idna(ifluid),i,j,:))
	enddo
        vz(2:nz-1)=(vz(1:nz-2) + vz(3:nz) + 2.0*vz(2:nz-1))*0.25
        vz(1)  = vz(2)
        vz(nz) = vz(nz-1)
        by_z=b(iby,i,j,:)
#ifdef ORIG
        call tvdb(vzby,by_z,vz,nz,dt,dz)
#else /* ORIG */
        call tvdb(vzby,by_z,vz,nz,dt)
#endif /* ORIG */
        wa(i,j,:) = vzby
      end do
    end do

    if(dimensions .eq. '3d') call bnd_emf(wa, 'vzby', 'zdim')
    call bnd_emf(wa, 'vzby', 'xdim')
    call bnd_emf(wa, 'vzby', 'ydim')

  end subroutine advectby_z

end module advects
