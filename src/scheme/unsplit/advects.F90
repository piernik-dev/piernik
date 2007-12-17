#include "mhd.def"
module advects    ! unsplit advects
  use start
  use arrays
  use grid
  use tv
  use mag_boundaries
  contains

   subroutine advectby_x
    use constants, only : big

    implicit none
    real, dimension(nx) :: vxby,by_x,vx
    integer                j,k,jm
    
    do k=1,nz
      do j=2,ny
        jm=j-1
        vx=(u(imxa,:,jm,k)+u(imxa,:,j,k))/(u(idna,:,jm,k)+u(idna,:,j,k))
        vx=(eoshift(vx,shift=-1,boundary=big) &
           +eoshift(vx,shift=1,boundary=big)+2.*vx)/4.
        by_x=b(iby,:,j,k)
        call tvdb(vxby,by_x,vx,nx,dt)
        wa(:,j,k) = vxby
      end do
    end do

    call bnd_emf(wa, 'vxby', 'xdim')
    call bnd_emf(wa, 'vxby', 'ydim')
    if(dimensions .eq. '3d') call bnd_emf(wa, 'vxby', 'zdim')

    Lb(iby,:,:,:) = Lb(iby,:,:,:) - wa/dx
    wa = cshift(wa,shift=-1,dim=xdim)
    Lb(iby,:,:,:) = Lb(iby,:,:,:) + wa/dx
    Lb(ibx,:,:,:) = Lb(ibx,:,:,:) - wa/dy
    wa = cshift(wa,shift=1,dim=ydim)
    Lb(ibx,:,:,:) = Lb(ibx,:,:,:) + wa/dy

  end subroutine advectby_x

  subroutine advectbz_x
    use constants, only : big

    implicit none
    real, dimension(nx) :: vxbz,bz_x,vx
    integer                j,k,km

    vxbz = 0.0

    do k=2,nz
      km=k-1
      do j=1,ny
        vx=(u(imxa,:,j,km)+u(imxa,:,j,k))/(u(idna,:,j,km)+u(idna,:,j,k))
        vx=(eoshift(vx,-1,boundary=big) &
           +eoshift(vx,1,boundary=big)+2*vx)/4
        bz_x=b(ibz,:,j,k)
        call tvdb(vxbz,bz_x,vx,nx,dt)
        wa(:,j,k) = vxbz
      end do
    end do

    call bnd_emf(wa, 'vxbz', 'xdim')
    call bnd_emf(wa, 'vxbz', 'ydim')
    if(dimensions .eq. '3d') call bnd_emf(wa, 'vxbz', 'zdim')

    Lb(ibz,:,:,:) = Lb(ibz,:,:,:) - wa/dx
    wa = cshift(wa,shift=-1,dim=xdim)
    Lb(ibz,:,:,:) = Lb(ibz,:,:,:) + wa/dx
    Lb(ibx,:,:,:) = Lb(ibx,:,:,:) - wa/dz
    wa = cshift(wa,shift=1,dim=zdim)
    Lb(ibx,:,:,:) = Lb(ibx,:,:,:) + wa/dz

  end subroutine advectbz_x

  subroutine advectbz_y
    use constants, only : big

    implicit none
    real, dimension(ny)   :: vybz,bz_y,vy
    integer                  i,k,km

    vybz = 0.0

    do k=2,nz
      km=k-1
      do i=1,nx  
        vy=(u(imya,i,:,km)+u(imya,i,:,k))/(u(idna,i,:,km)+u(idna,i,:,k))
        vy=(eoshift(vy,-1,boundary=big) &
           +eoshift(vy,1,boundary=big)+2*vy)/4
        bz_y=b(ibz,i,:,k)
        call tvdb(vybz,bz_y,vy,ny,dt)
        wa(i,:,k) = vybz
      end do
    end do

    call bnd_emf(wa, 'vybz', 'ydim')
    if(dimensions .eq. '3d') call bnd_emf(wa, 'vybz', 'zdim')
    call bnd_emf(wa, 'vybz', 'xdim')
  
    Lb(ibz,:,:,:) = Lb(ibz,:,:,:) - wa/dy
    wa = cshift(wa,shift=-1,dim=ydim)
    Lb(ibz,:,:,:) = Lb(ibz,:,:,:) + wa/dy
    Lb(iby,:,:,:) = Lb(iby,:,:,:) - wa/dz
    wa = cshift(wa,shift=1,dim=zdim)
    Lb(iby,:,:,:) = Lb(iby,:,:,:) + wa/dz

  end subroutine advectbz_y

  subroutine advectbx_y
    use constants, only : big

    implicit none
    real, dimension(ny) :: vybx,bx_y,vy
    integer                k,i,im

    vybx = 0.0

    do k=1,nz
      do i=2,nx
        im=i-1
        vy=(u(imya,im,:,k)+u(imya,i,:,k))/(u(idna,im,:,k)+u(idna,i,:,k))
        vy=(eoshift(vy,-1,boundary=big) &
           +eoshift(vy,1,boundary=big)+2*vy)/4
        bx_y=b(ibx,i,:,k)
        call tvdb(vybx,bx_y,vy,ny,dt)
        wa(i,:,k) = vybx
      end do
    end do

    call bnd_emf(wa, 'vybx', 'ydim')
    if(dimensions .eq. '3d') call bnd_emf(wa, 'vybx', 'zdim')
    call bnd_emf(wa, 'vybx', 'xdim')

    Lb(ibx,:,:,:) = Lb(ibx,:,:,:) - wa/dy
    wa = cshift(wa,shift=-1,dim=ydim)
    Lb(ibx,:,:,:) = Lb(ibx,:,:,:) + wa/dy
    Lb(iby,:,:,:) = Lb(iby,:,:,:) - wa/dx
    wa = cshift(wa,shift=1,dim=xdim)
    Lb(iby,:,:,:) = Lb(iby,:,:,:) + wa/dx

  end subroutine advectbx_y

  subroutine advectbx_z
    use constants, only : big

    implicit none
    real, dimension(nz)  :: vzbx,bx_z,vz
    integer                 j,i,im

    do j=1,ny
      do i=2,nx
        im=i-1
        vz=(u(imza,im,j,:)+u(imza,i,j,:))/(u(idna,im,j,:)+u(idna,i,j,:))
        vz=(eoshift(vz,-1,boundary=big) &
           +eoshift(vz,1,boundary=big)+2*vz)/4
        bx_z=b(ibx,i,j,:)
        call tvdb(vzbx,bx_z,vz,nz,dt)
        wa(i,j,:) = vzbx
      end do
    end do
    
    if(dimensions .eq. '3d') call bnd_emf(wa, 'vzbx', 'zdim')
    call bnd_emf(wa, 'vzbx', 'xdim')
    call bnd_emf(wa, 'vzbx', 'ydim')

    Lb(ibx,:,:,:) = Lb(ibx,:,:,:) - wa/dz
    wa = cshift(wa,shift=-1,dim=zdim)
    Lb(ibx,:,:,:) = Lb(ibx,:,:,:) + wa/dz
    Lb(ibz,:,:,:) = Lb(ibz,:,:,:) - wa/dx
    wa = cshift(wa,shift=1,dim=xdim)
    Lb(ibz,:,:,:) = Lb(ibz,:,:,:) + wa/dx

  end subroutine advectbx_z

  subroutine advectby_z
    use constants, only : big
    
    implicit none
    real, dimension(nz)  :: vzby,by_z,vz
    integer                 i,j,jm

    do j=2,ny
      jm=j-1
      do i=1,nx
        vz=(u(imza,i,jm,:)+u(imza,i,j,:))/(u(idna,i,jm,:)+u(idna,i,j,:))
        vz=(eoshift(vz,-1,boundary=big) &
           +eoshift(vz,1,boundary=big)+2*vz)/4
        by_z=b(iby,i,j,:)
        call tvdb(vzby,by_z,vz,nz,dt)
        wa(i,j,:) = vzby
      end do
    end do

    if(dimensions .eq. '3d') call bnd_emf(wa, 'vzby', 'zdim')
    call bnd_emf(wa, 'vzby', 'xdim')
    call bnd_emf(wa, 'vzby', 'ydim')

    Lb(iby,:,:,:) = Lb(iby,:,:,:) - wa/dz
    wa = cshift(wa,shift=-1,dim=zdim)
    Lb(iby,:,:,:) = Lb(iby,:,:,:) + wa/dz
    Lb(ibz,:,:,:) = Lb(ibz,:,:,:) - wa/dy
    wa = cshift(wa,shift=1,dim=ydim)
    Lb(ibz,:,:,:) = Lb(ibz,:,:,:) + wa/dy
 
  end subroutine advectby_z

end module advects
