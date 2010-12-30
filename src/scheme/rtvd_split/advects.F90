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
!    Initial implementation of PIERNIK code was based on TVD split MHD code by
!    Ue-Li Pen
!        see: Pen, Arras & Wong (2003) for algorithm and
!             http://www.cita.utoronto.ca/~pen/MHD
!             for original source code "mhd.f90"
!
!    For full list of developers see $PIERNIK_HOME/license/pdt.txt
!
#include "piernik.h"

module advects
! pulled by MAGNETIC
   implicit none

   private
   public  :: advectby_x, advectbz_x, advectbz_y, advectbx_y, advectbx_z, advectby_z

contains

  subroutine advectby_x

    use arrays,        only: b, u, wa
    use fluidindex,    only: iby, nvar
    use grid,          only: cg
    use magboundaries, only: bnd_emf
    use mpisetup,      only: dt, xdim, ydim, zdim, has_dir
    use rtvd,          only: tvdb

    implicit none
    real, dimension(cg%nx) :: vxby, by_x, vx
    integer             :: j, k, jm, j_s, j_e

    vxby = 0.0

    if (has_dir(ydim)) then ! This can be converted by use grid, only: Dy; j_s = 1 + D_y; j_e = 1 + D_y * (cg%ny - 1) . How will it impact performance?
        j_s = 2
        j_e = cg%ny
    else
        j_s = 1
        j_e = 1
    endif

    do k = cg%ks, cg%ke
      do j = j_s, j_e         ! cg%nyb cg%is /= 1 in by_x
        jm=j-1
        vx=0.0
        if (jm == 0) then
           vx = u(nvar%ion%imx,:,1,k) / u(nvar%ion%idn,:,1,k)
        else
           vx=(u(nvar%ion%imx,:,jm,k)+u(nvar%ion%imx,:,j,k))/(u(nvar%ion%idn,:,jm,k)+u(nvar%ion%idn,:,j,k))
        endif
        vx(2:cg%nx-1)=(vx(1:cg%nx-2) + vx(3:cg%nx) + 2.0*vx(2:cg%nx-1))*0.25
        vx(1)  = vx(2)
        vx(cg%nx) = vx(cg%nx-1)
        by_x=b(iby,:,j,k)

        call tvdb(vxby,by_x,vx, cg%nx,dt, cg%idx)

        wa(:,j,k) = vxby
      enddo
    enddo

    if (has_dir(xdim)) call bnd_emf(wa, 'vxby', 'xdim')
    if (has_dir(ydim)) call bnd_emf(wa, 'vxby', 'ydim')
    if (has_dir(zdim)) call bnd_emf(wa, 'vxby', 'zdim')

  end subroutine advectby_x

  subroutine advectbz_x
    use arrays,        only: b, u, wa
    use fluidindex,    only: ibz, nvar
    use grid,          only: cg
    use magboundaries, only: bnd_emf
    use mpisetup,      only: dt, xdim, ydim, zdim, has_dir
    use rtvd,          only: tvdb

    implicit none
    real, dimension(cg%nx) :: vxbz, bz_x, vx
    integer             :: j, k, km, k_s, k_e

    vxbz = 0.0

    if (has_dir(zdim)) then
        k_s = 2
        k_e = cg%nz
    else
        k_s = 1
        k_e = 1
    endif

    do k=k_s,k_e
      km=k-1
      do j=cg%js, cg%je
        vx=0.0
        if (km == 0) then
           vx = u(nvar%ion%imx,:,j,1) / u(nvar%ion%idn,:,j,1)
        else
           vx=(u(nvar%ion%imx,:,j,km)+u(nvar%ion%imx,:,j,k))/(u(nvar%ion%idn,:,j,km)+u(nvar%ion%idn,:,j,k))
        endif
        vx(2:cg%nx-1)=(vx(1:cg%nx-2) + vx(3:cg%nx) + 2.0*vx(2:cg%nx-1))*0.25
        vx(1)  = vx(2)
        vx(cg%nx) = vx(cg%nx-1)
        bz_x=b(ibz,:,j,k)

        call tvdb(vxbz,bz_x,vx, cg%nx,dt, cg%idx)

        wa(:,j,k) = vxbz
      enddo
    enddo

    if (has_dir(xdim)) call bnd_emf(wa, 'vxbz', 'xdim')
    if (has_dir(ydim)) call bnd_emf(wa, 'vxbz', 'ydim')
    if (has_dir(zdim)) call bnd_emf(wa, 'vxbz', 'zdim')

  end subroutine advectbz_x

  subroutine advectbz_y
    use arrays,        only: b, u, wa
    use fluidindex,    only: ibz, nvar
    use grid,          only: cg
    use magboundaries, only: bnd_emf
    use mpisetup,      only: dt, xdim, ydim, zdim, has_dir
    use rtvd,          only: tvdb

    implicit none
    real, dimension(cg%ny) :: vybz, bz_y, vy
    integer             :: i, k, km, k_s, k_e

    vybz = 0.0
    if (has_dir(zdim)) then
        k_s = 2
        k_e = cg%nz
    else
        k_s = 1
        k_e = 1
    endif

    do k=k_s,k_e
      km=k-1
      do i=cg%is, cg%ie
        vy=0.0
        if (km == 0) then
           vy = u(nvar%ion%imy,i,:,1) / u(nvar%ion%idn,i,:,1)
        else
           vy=(u(nvar%ion%imy,i,:,km)+u(nvar%ion%imy,i,:,k))/(u(nvar%ion%idn,i,:,km)+u(nvar%ion%idn,i,:,k))
        endif
        vy(2:cg%ny-1)=(vy(1:cg%ny-2) + vy(3:cg%ny) + 2.0*vy(2:cg%ny-1))*0.25
        vy(1)  = vy(2)
        vy(cg%ny) = vy(cg%ny-1)
        bz_y=b(ibz,i,:,k)

        call tvdb(vybz,bz_y,vy, cg%ny,dt, cg%idy)

        wa(i,:,k) = vybz
      enddo
    enddo

    if (has_dir(ydim)) call bnd_emf(wa, 'vybz', 'ydim')
    if (has_dir(zdim)) call bnd_emf(wa, 'vybz', 'zdim')
    if (has_dir(xdim)) call bnd_emf(wa, 'vybz', 'xdim')

  end subroutine advectbz_y

  subroutine advectbx_y
    use arrays,        only: b, u, wa
    use fluidindex,    only: ibx, nvar
    use grid,          only: cg
    use magboundaries, only: bnd_emf
    use mpisetup,      only: dt, xdim, ydim, zdim, has_dir
    use rtvd,          only: tvdb

    implicit none
    real, dimension(cg%ny) :: vybx, bx_y, vy
    integer             :: k, i, im, i_s, i_e

    vybx = 0.0

    if (has_dir(xdim)) then
        i_s = 2
        i_e = cg%nx
    else
        i_s = 1
        i_e = 1
    endif

    do k=cg%ks, cg%ke
      do i=i_s,i_e
        im=i-1
        vy=0.0
        if (im == 0) then
           vy = u(nvar%ion%imy,1,:,k) / u(nvar%ion%idn,1,:,k)
        else
           vy=(u(nvar%ion%imy,im,:,k)+u(nvar%ion%imy,i,:,k))/(u(nvar%ion%idn,im,:,k)+u(nvar%ion%idn,i,:,k))
        endif
        vy(2:cg%ny-1)=(vy(1:cg%ny-2) + vy(3:cg%ny) + 2.0*vy(2:cg%ny-1))*0.25
        vy(1)  = vy(2)
        vy(cg%ny) = vy(cg%ny-1)
        bx_y=b(ibx,i,:,k)

        call tvdb(vybx,bx_y,vy, cg%ny,dt, cg%idy)

        wa(i,:,k) = vybx
      enddo
    enddo

    if (has_dir(ydim)) call bnd_emf(wa, 'vybx', 'ydim')
    if (has_dir(zdim)) call bnd_emf(wa, 'vybx', 'zdim')
    if (has_dir(xdim)) call bnd_emf(wa, 'vybx', 'xdim')

  end subroutine advectbx_y

  subroutine advectbx_z
    use arrays,        only: b, u, wa
    use fluidindex,    only: ibx, nvar
    use grid,          only: cg
    use magboundaries, only: bnd_emf
    use mpisetup,      only: dt, xdim, ydim, zdim, has_dir
    use rtvd,          only: tvdb

    implicit none
    real, dimension(cg%nz)  :: vzbx, bx_z, vz
    integer              :: j, i, im, i_s, i_e

    vzbx = 0.0

    if (has_dir(xdim)) then
        i_s = 2
        i_e = cg%nx
    else
        i_s = 1
        i_e = 1
    endif
    do j=cg%js, cg%je
      do i=i_s,i_e
        im=i-1
        vz=0.0
        if (im == 0) then
           vz = u(nvar%ion%imz,1,j,:) / u(nvar%ion%idn,1,j,:)
        else
           vz = (u(nvar%ion%imz,im,j,:)+u(nvar%ion%imz,i,j,:))/(u(nvar%ion%idn,im,j,:)+u(nvar%ion%idn,i,j,:))
        endif
        vz(2:cg%nz-1)=(vz(1:cg%nz-2) + vz(3:cg%nz) + 2.0*vz(2:cg%nz-1))*0.25
        vz(1)  = vz(2)
        vz(cg%nz) = vz(cg%nz-1)
        bx_z=b(ibx,i,j,:)

        call tvdb(vzbx,bx_z,vz, cg%nz,dt, cg%idz)

        wa(i,j,:) = vzbx
      enddo
    enddo

    if (has_dir(zdim)) call bnd_emf(wa, 'vzbx', 'zdim')
    if (has_dir(xdim)) call bnd_emf(wa, 'vzbx', 'xdim')
    if (has_dir(ydim)) call bnd_emf(wa, 'vzbx', 'ydim')

  end subroutine advectbx_z

  subroutine advectby_z
    use arrays,        only: b, u, wa
    use fluidindex,    only: iby, nvar
    use grid,          only: cg
    use magboundaries, only: bnd_emf
    use mpisetup,      only: dt, xdim, ydim, zdim, has_dir
    use rtvd,          only: tvdb

    implicit none
    real, dimension(cg%nz) :: vzby, by_z, vz
    integer             :: i, j, jm, j_s, j_e

    vzby = 0.0
    if (has_dir(ydim)) then
        j_s = 2
        j_e = cg%ny
    else
        j_s = 1
        j_e = 1
    endif

    do j=j_s,j_e
      jm=j-1
      do i=cg%is, cg%ie
        vz=0.0
        if (jm == 0) then
           vz = u(nvar%ion%imz,i,1,:) / u(nvar%ion%idn,i,1,:)
        else
           vz=(u(nvar%ion%imz,i,jm,:)+u(nvar%ion%imz,i,j,:))/(u(nvar%ion%idn,i,jm,:)+u(nvar%ion%idn,i,j,:))
        endif
        vz(2:cg%nz-1)=(vz(1:cg%nz-2) + vz(3:cg%nz) + 2.0*vz(2:cg%nz-1))*0.25
        vz(1)  = vz(2)
        vz(cg%nz) = vz(cg%nz-1)
        by_z=b(iby,i,j,:)

        call tvdb(vzby,by_z,vz, cg%nz,dt, cg%idz)

        wa(i,j,:) = vzby
      enddo
    enddo

    if (has_dir(zdim)) call bnd_emf(wa, 'vzby', 'zdim')
    if (has_dir(xdim)) call bnd_emf(wa, 'vzby', 'xdim')
    if (has_dir(ydim)) call bnd_emf(wa, 'vzby', 'ydim')

  end subroutine advectby_z

end module advects
