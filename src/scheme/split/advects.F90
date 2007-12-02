
#include "mhd.def"
module advects    ! split advects
  use start
  use arrays
  use grid
  use time, only : dt
  use tv
  use mag_boundaries
  contains
  subroutine advectby_x

    implicit none
    real, dimension(nx) :: vxby,by_x,vx
    integer                j,k,jm

    vxby = 0.0
    
    do k=1,nz
      do j=2,ny
        jm=j-1
        vx=(u(imxa,:,jm,k)+u(imxa,:,j,k))/(u(idna,:,jm,k)+u(idna,:,j,k))
        vx(2:nx-1)=(vx(1:nx-2) + vx(3:nx) + 2.0*vx(2:nx-1))*0.25
        vx(1)  = vx(2)
        vx(nx) = vx(nx-1)
        by_x=b(iby,:,j,k)
#ifdef ORIG
        call tvdb(vxby,by_x,vx,nx,dt,dx)
#else 
        call tvdb(vxby,by_x,vx,nx,dt)
#endif 
        wa(:,j,k) = vxby
      end do
    end do

    call bnd_emf(wa, 'vxby', 'xdim')
    call bnd_emf(wa, 'vxby', 'ydim')
    if(dimensions .eq. '3d') call bnd_emf(wa, 'vxby', 'zdim')

  end subroutine advectby_x

  subroutine advectbz_x

    implicit none
    real, dimension(nx) :: vxbz,bz_x,vx
    integer                j,k,km

    vxbz = 0.0

    do k=2,nz
      km=k-1
      do j=1,ny
        vx=(u(imxa,:,j,km)+u(imxa,:,j,k))/(u(idna,:,j,km)+u(idna,:,j,k))
        vx(2:nx-1)=(vx(1:nx-2) + vx(3:nx) + 2.0*vx(2:nx-1))*0.25
        vx(1)  = vx(2)
        vx(nx) = vx(nx-1)
        bz_x=b(ibz,:,j,k)
#ifdef ORIG
        call tvdb(vxbz,bz_x,vx,nx,dt,dx)
#else 
        call tvdb(vxbz,bz_x,vx,nx,dt)
#endif
        wa(:,j,k) = vxbz
      end do
    end do

    call bnd_emf(wa, 'vxbz', 'xdim')
    call bnd_emf(wa, 'vxbz', 'ydim')
    if(dimensions .eq. '3d') call bnd_emf(wa, 'vxbz', 'zdim')

  end subroutine advectbz_x

  subroutine advectbz_y

    implicit none
    real, dimension(ny)   :: vybz,bz_y,vy
    integer                  i,k,km

    vybz = 0.0

    do k=2,nz
      km=k-1
      do i=1,nx  
        vy=(u(imya,i,:,km)+u(imya,i,:,k))/(u(idna,i,:,km)+u(idna,i,:,k))
        vy(2:ny-1)=(vy(1:ny-2) + vy(3:ny) + 2.0*vy(2:ny-1))*0.25
        vy(1)  = vy(2)
        vy(ny) = vy(ny-1)
        bz_y=b(ibz,i,:,k)
#ifdef ORIG
        call tvdb(vybz,bz_y,vy,ny,dt,dy)
#else 
        call tvdb(vybz,bz_y,vy,ny,dt)
#endif 
        wa(i,:,k) = vybz
      end do
    end do

    call bnd_emf(wa, 'vybz', 'ydim')
    if(dimensions .eq. '3d') call bnd_emf(wa, 'vybz', 'zdim')
    call bnd_emf(wa, 'vybz', 'xdim')
  
  end subroutine advectbz_y

  subroutine advectbx_y

    implicit none
    real, dimension(ny) :: vybx,bx_y,vy
    integer                k,i,im

    vybx = 0.0

    do k=1,nz
      do i=2,nx
        im=i-1
        vy=(u(imya,im,:,k)+u(imya,i,:,k))/(u(idna,im,:,k)+u(idna,i,:,k))
        vy(2:ny-1)=(vy(1:ny-2) + vy(3:ny) + 2.0*vy(2:ny-1))*0.25
        vy(1)  = vy(2)
        vy(ny) = vy(ny-1)
        bx_y=b(ibx,i,:,k)
#ifdef ORIG
        call tvdb(vybx,bx_y,vy,ny,dt,dy)
#else 
        call tvdb(vybx,bx_y,vy,ny,dt)
#endif 
        wa(i,:,k) = vybx
      end do
    end do

    call bnd_emf(wa, 'vybx', 'ydim')
    if(dimensions .eq. '3d') call bnd_emf(wa, 'vybx', 'zdim')
    call bnd_emf(wa, 'vybx', 'xdim')

  end subroutine advectbx_y

  subroutine advectbx_z

    implicit none
    real, dimension(nz)  :: vzbx,bx_z,vz
    integer                 j,i,im

    vzbx = 0.0

    do j=1,ny
      do i=2,nx
        im=i-1
        vz=(u(imza,im,j,:)+u(imza,i,j,:))/(u(idna,im,j,:)+u(idna,i,j,:))
        vz(2:nz-1)=(vz(1:nz-2) + vz(3:nz) + 2.0*vz(2:nz-1))*0.25
        vz(1)  = vz(2)
        vz(nz) = vz(nz-1)
        bx_z=b(ibx,i,j,:)
#ifdef ORIG
        call tvdb(vzbx,bx_z,vz,nz,dt,dz)
#else
        call tvdb(vzbx,bx_z,vz,nz,dt)
#endif
        wa(i,j,:) = vzbx
      end do
    end do
    
    if(dimensions .eq. '3d') call bnd_emf(wa, 'vzbx', 'zdim')
    call bnd_emf(wa, 'vzbx', 'xdim')
    call bnd_emf(wa, 'vzbx', 'ydim')

  end subroutine advectbx_z

  subroutine advectby_z
    
    implicit none
    real, dimension(nz)  :: vzby,by_z,vz
    integer                 i,j,jm

    vzby = 0.0

    do j=2,ny
      jm=j-1
      do i=1,nx
        vz=(u(imza,i,jm,:)+u(imza,i,j,:))/(u(idna,i,jm,:)+u(idna,i,j,:))
        vz(2:nz-1)=(vz(1:nz-2) + vz(3:nz) + 2.0*vz(2:nz-1))*0.25
        vz(1)  = vz(2)
        vz(nz) = vz(nz-1)
        by_z=b(iby,i,j,:)
#ifdef ORIG
        call tvdb(vzby,by_z,vz,nz,dt,dz)
#else
        call tvdb(vzby,by_z,vz,nz,dt)
#endif
        wa(i,j,:) = vzby
      end do
    end do

    if(dimensions .eq. '3d') call bnd_emf(wa, 'vzby', 'zdim')
    call bnd_emf(wa, 'vzby', 'xdim')
    call bnd_emf(wa, 'vzby', 'ydim')

  end subroutine advectby_z

end module advects
