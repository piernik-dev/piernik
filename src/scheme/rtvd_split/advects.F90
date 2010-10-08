! $Id$
!
! PIERNIK Code Copyright (C) 2006 Michal Hanasz
!
!    This file is part of PIERNIK code.
!
!    PIERNIK is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    PIERNIK is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with PIERNIK.  If not, see <http://www.gnu.org/licenses/>.
!
!    Initial implemetation of PIERNIK code was based on TVD split MHD code by
!    Ue-Li Pen
!        see: Pen, Arras & Wong (2003) for algorithm and
!             http://www.cita.utoronto.ca/~pen/MHD
!             for original source code "mhd.f90"
!
!    For full list of developers see $PIERNIK_HOME/license/pdt.txt
!
#include "piernik.def"

module advects
  contains

  subroutine advectby_x
    use mpisetup,  only: dt
    use fluidindex, only: ibx,iby,ibz
    use initionized, only: idni,imxi
    use arrays, only: b,u,wa
    use grid,   only: dx,nx,ny,nz,nxd,nyd,nzd,ks,ke
    use rtvd,     only: tvdb
    use magboundaries, only: bnd_emf

    implicit none
    real, dimension(nx) :: vxby,by_x,vx
    integer             :: j,k,jm, j_s,j_e

    vxby = 0.0

    if (nyd /= 1) then
        j_s = 2
        j_e = ny
    else
        j_s = 1
        j_e = 1
    endif

    do k=1,nz
      do j=j_s,j_e         ! nyd is /= 1 in by_x
        jm=j-1
        vx=0.0
        if (jm == 0) then
           vx = u(imxi,:,1,k) / u(idni,:,1,k)
        else
           vx=(u(imxi,:,jm,k)+u(imxi,:,j,k))/(u(idni,:,jm,k)+u(idni,:,j,k))
        endif
        vx(2:nx-1)=(vx(1:nx-2) + vx(3:nx) + 2.0*vx(2:nx-1))*0.25
        vx(1)  = vx(2)
        vx(nx) = vx(nx-1)
        by_x=b(iby,:,j,k)

        call tvdb(vxby,by_x,vx,nx,dt,dx)

        wa(:,j,k) = vxby
      enddo
    enddo

    if (nxd /= 1) call bnd_emf(wa, 'vxby', 'xdim')
    if (nyd /= 1) call bnd_emf(wa, 'vxby', 'ydim')
    if (nzd /= 1) call bnd_emf(wa, 'vxby', 'zdim')

  end subroutine advectby_x

  subroutine advectbz_x
    use mpisetup,  only: dt
    use fluidindex, only: ibx,iby,ibz
    use initionized, only: idni,imxi
    use arrays, only: b,u,wa
    use grid,   only: dx,nx,ny,nz,nxd,nyd,nzd,js,je
    use rtvd,   only: tvdb
    use magboundaries, only: bnd_emf

    implicit none
    real, dimension(nx) :: vxbz,bz_x,vx
    integer             :: j,k,km,k_s,k_e

    vxbz = 0.0

    if (nzd /= 1) then
        k_s = 2
        k_e = nz
    else
        k_s = 1
        k_e = 1
    endif

    do k=k_s,k_e
      km=k-1
      do j=js,je
        vx=0.0
        if (km == 0) then
           vx = u(imxi,:,j,1) / u(idni,:,j,1)
        else
           vx=(u(imxi,:,j,km)+u(imxi,:,j,k))/(u(idni,:,j,km)+u(idni,:,j,k))
        endif
        vx(2:nx-1)=(vx(1:nx-2) + vx(3:nx) + 2.0*vx(2:nx-1))*0.25
        vx(1)  = vx(2)
        vx(nx) = vx(nx-1)
        bz_x=b(ibz,:,j,k)

        call tvdb(vxbz,bz_x,vx,nx,dt,dx)

        wa(:,j,k) = vxbz
      enddo
    enddo

    if (nxd /= 1) call bnd_emf(wa, 'vxbz', 'xdim')
    if (nyd /= 1) call bnd_emf(wa, 'vxbz', 'ydim')
    if (nzd /= 1) call bnd_emf(wa, 'vxbz', 'zdim')

  end subroutine advectbz_x

  subroutine advectbz_y
    use mpisetup,  only: dt
    use fluidindex, only: ibx,iby,ibz
    use initionized, only: idni,imyi
    use arrays, only: b,u,wa
    use grid,   only: dy,nx,ny,nz,is,ie,nxd,nyd,nzd
    use rtvd,   only: tvdb
    use magboundaries, only: bnd_emf

    implicit none
    real, dimension(ny)   :: vybz,bz_y,vy
    integer               :: i,k,km,k_s,k_e

    vybz = 0.0
    if (nzd /= 1) then
        k_s = 2
        k_e = nz
    else
        k_s = 1
        k_e = 1
    endif

    do k=k_s,k_e
      km=k-1
      do i=is,ie
        vy=0.0
        if (km == 0) then
           vy = u(imyi,i,:,1) / u(idni,i,:,1)
        else
           vy=(u(imyi,i,:,km)+u(imyi,i,:,k))/(u(idni,i,:,km)+u(idni,i,:,k))
        endif
        vy(2:ny-1)=(vy(1:ny-2) + vy(3:ny) + 2.0*vy(2:ny-1))*0.25
        vy(1)  = vy(2)
        vy(ny) = vy(ny-1)
        bz_y=b(ibz,i,:,k)

        call tvdb(vybz,bz_y,vy,ny,dt,dy)

        wa(i,:,k) = vybz
      enddo
    enddo

    if (nyd /= 1) call bnd_emf(wa, 'vybz', 'ydim')
    if (nzd /= 1) call bnd_emf(wa, 'vybz', 'zdim')
    if (nxd /= 1) call bnd_emf(wa, 'vybz', 'xdim')

  end subroutine advectbz_y

  subroutine advectbx_y
    use mpisetup,  only: dt
    use fluidindex, only: ibx,iby,ibz
    use initionized, only: idni,imyi
    use arrays, only: b,u,wa
    use grid,   only: dy,nx,ny,nz,nxd,nzd,nyd,ks,ke
    use rtvd,     only: tvdb
    use magboundaries, only: bnd_emf

    implicit none
    real, dimension(ny) :: vybx,bx_y,vy
    integer             :: k,i,im,i_s,i_e

    vybx = 0.0

    if (nxd /= 1) then
        i_s = 2
        i_e = nx
    else
        i_s = 1
        i_e = 1
    endif

    do k=ks,ke
      do i=i_s,i_e
        im=i-1
        vy=0.0
        if (im == 0) then
           vy = u(imyi,1,:,k) / u(idni,1,:,k)
        else
           vy=(u(imyi,im,:,k)+u(imyi,i,:,k))/(u(idni,im,:,k)+u(idni,i,:,k))
        endif
        vy(2:ny-1)=(vy(1:ny-2) + vy(3:ny) + 2.0*vy(2:ny-1))*0.25
        vy(1)  = vy(2)
        vy(ny) = vy(ny-1)
        bx_y=b(ibx,i,:,k)

        call tvdb(vybx,bx_y,vy,ny,dt,dy)

        wa(i,:,k) = vybx
      enddo
    enddo

    if (nyd /= 1) call bnd_emf(wa, 'vybx', 'ydim')
    if (nzd /= 1) call bnd_emf(wa, 'vybx', 'zdim')
    if (nxd /= 1) call bnd_emf(wa, 'vybx', 'xdim')

  end subroutine advectbx_y

  subroutine advectbx_z
    use mpisetup,  only: dt
    use fluidindex, only: ibx,iby,ibz
    use initionized, only: idni,imzi
    use arrays, only: b,u,wa
    use grid,   only: dz,nx,ny,nz,nxd,nzd,nyd,js,je
    use rtvd,     only: tvdb
    use magboundaries, only: bnd_emf

    implicit none
    real, dimension(nz)  :: vzbx,bx_z,vz
    integer              :: j,i,im,i_s,i_e

    vzbx = 0.0

    if (nxd /= 1) then
        i_s = 2
        i_e = nx
    else
        i_s = 1
        i_e = 1
    endif
    do j=js,je
      do i=i_s,i_e
        im=i-1
        vz=0.0
        if (im == 0) then
           vz = u(imzi,1,j,:) / u(idni,1,j,:)
        else
           vz = (u(imzi,im,j,:)+u(imzi,i,j,:))/(u(idni,im,j,:)+u(idni,i,j,:))
        endif
        vz(2:nz-1)=(vz(1:nz-2) + vz(3:nz) + 2.0*vz(2:nz-1))*0.25
        vz(1)  = vz(2)
        vz(nz) = vz(nz-1)
        bx_z=b(ibx,i,j,:)

        call tvdb(vzbx,bx_z,vz,nz,dt,dz)

        wa(i,j,:) = vzbx
      enddo
    enddo

    if (nzd /= 1) call bnd_emf(wa, 'vzbx', 'zdim')
    if (nxd /= 1) call bnd_emf(wa, 'vzbx', 'xdim')
    if (nyd /= 1) call bnd_emf(wa, 'vzbx', 'ydim')

  end subroutine advectbx_z

  subroutine advectby_z
    use mpisetup,  only: dt
    use fluidindex, only: ibx,iby,ibz
    use initionized, only: idni,imzi
    use arrays, only: b,u,wa
    use grid,   only: dz,nx,ny,nz,nzd,nyd,nxd,ie,is
    use rtvd,     only: tvdb
    use magboundaries, only: bnd_emf

    implicit none
    real, dimension(nz)  :: vzby,by_z,vz
    integer              :: i,j,jm,j_s,j_e

    vzby = 0.0
    if (nyd /= 1) then
        j_s = 2
        j_e = ny
    else
        j_s = 1
        j_e = 1
    endif

    do j=j_s,j_e
      jm=j-1
      do i=is,ie
        vz=0.0
        if (jm == 0) then
           vz = u(imzi,i,1,:) / u(idni,i,1,:)
        else
           vz=(u(imzi,i,jm,:)+u(imzi,i,j,:))/(u(idni,i,jm,:)+u(idni,i,j,:))
        endif
        vz(2:nz-1)=(vz(1:nz-2) + vz(3:nz) + 2.0*vz(2:nz-1))*0.25
        vz(1)  = vz(2)
        vz(nz) = vz(nz-1)
        by_z=b(iby,i,j,:)

        call tvdb(vzby,by_z,vz,nz,dt,dz)

        wa(i,j,:) = vzby
      enddo
    enddo

    if (nzd /= 1) call bnd_emf(wa, 'vzby', 'zdim')
    if (nxd /= 1) call bnd_emf(wa, 'vzby', 'xdim')
    if (nyd /= 1) call bnd_emf(wa, 'vzby', 'ydim')

  end subroutine advectby_z

end module advects
