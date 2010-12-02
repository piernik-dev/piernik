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
module sweeps     ! split sweeps
! pulled by ANY
   implicit none
   private
   public  :: sweepx, sweepy, sweepz
#if defined SHEAR && defined FLUID_INTERACTIONS
   public  :: source_terms_y
#endif /* SHEAR && FLUID_INTERACTIONS */
   contains
#if defined SHEAR && defined FLUID_INTERACTIONS
   subroutine source_terms_y
      use arrays,          only: u
      use fluidboundaries, only: all_fluid_boundaries
      use fluidindex,      only: iarr_all_dn, iarr_all_mx, iarr_all_my, nvar
      use grid,            only: nx, nz
      use mpisetup,        only: dt
      use shear,           only: omega, qshear
#ifdef DUST
      use initdust,        only: dragc_gas_dust
#endif /* DUST */

    implicit none
    real, dimension(size(iarr_all_my),nx,nz) :: vxr, v_r, rotaccr
    real, dimension(nx,nz)      :: epsa
    real, dimension(nvar%all,nx,nz) :: u1
    integer :: ind,i
    real, parameter, dimension(2) :: fac = [0.5, 1.0]

    u1(:,:,:) = u(:,:,1,:)

    do i = 1,2
       where (u1(iarr_all_dn,:,:) > 0.0)
          vxr(:,:,:) = u1(iarr_all_mx,:,:) / u1(iarr_all_dn,:,:)
          v_r(:,:,:) = u1(iarr_all_my,:,:) / u1(iarr_all_dn,:,:)
       elsewhere
          vxr(:,:,:) = 0.0
          v_r(:,:,:) = 0.0
       endwhere
       epsa(:,:) =  u1(iarr_all_dn(2),:,:) / u1(iarr_all_dn(1),:,:)

       do ind=1,size(iarr_all_my)
          if (ind == 1) then
            rotaccr(ind,:,:) = - dragc_gas_dust * epsa(:,:) * (v_r(1,:,:) - v_r(2,:,:))
          else
            rotaccr(ind,:,:) = - dragc_gas_dust * 1.0    * (v_r(2,:,:) - v_r(1,:,:))
          endif
!         rotaccr(ind,:,:) = rotaccr(ind,:,:) - 2.0*omega*vxr(ind,:,:)
          rotaccr(ind,:,:) = rotaccr(ind,:,:) + (qshear-2.0)*omega*vxr(ind,:,:)
       enddo

       where (u1(iarr_all_dn,:,:) > 0.0)
          u1(iarr_all_my,:,:) = u1(iarr_all_my,:,:) + fac(i)*dt*rotaccr(:,:,:)*u(iarr_all_dn,:,1,:)
       endwhere
       u(iarr_all_my,:,1,:) = u1(iarr_all_my,:,:)
    enddo

    call all_fluid_boundaries

  end subroutine source_terms_y
#endif /* SHEAR && FLUID_INTERACTIONS */
  subroutine sweepx

    use arrays,          only: u, b
    use fluidboundaries, only: all_fluid_boundaries
    use fluidindex,      only: nvar, iarr_all_swpx, ibx, iby, ibz, nmag
    use grid,            only: dx, nx, ks, ke, js, je, ydim, zdim, has_dir
    use mpisetup,        only: dt
    use rtvd,            only: relaxing_tvd
#ifdef COSM_RAYS
    use crhelpers,       only: div_v
#endif /* COSM_RAYS */

    implicit none
    real, dimension(nmag,nx)     :: b_x
    real, dimension(nvar%all,nx) :: u_x
    integer :: j,k,jp,kp

    b_x = 0.0
    u_x = 0.0

#ifdef COSM_RAYS
    call div_v(nvar%ion%pos)
#endif /* COSM_RAYS */
    do k=ks,ke
      kp=k+1
      do j=js,je
          jp=j+1

#ifdef MAGNETIC
          b_x=0.5*b(:,:,j,k)
          b_x(ibx,1:nx-1)=b_x(ibx,1:nx-1)+b_x(ibx,2:nx);       b_x(ibx,nx) = b_x(ibx,nx-1)
          if (has_dir(ydim) .and. j <= je)  b_x(iby,:)=b_x(iby,:)+0.5*b(iby,:,jp,k)
          if (has_dir(zdim) .and. k <= ke)  b_x(ibz,:)=b_x(ibz,:)+0.5*b(ibz,:,j,kp)
#endif /* MAGNETIC */

        u_x(iarr_all_swpx,:)=u(:,:,j,k)

        call relaxing_tvd(nx, u_x, b_x, 'xsweep', j, k, dx, dt)
        u(:,:,j,k)=u_x(iarr_all_swpx,:)
      enddo
    enddo

    call all_fluid_boundaries


  end subroutine sweepx
!------------------------------------------------------------------------------------------
  subroutine sweepy
    use arrays,          only: u, b
    use fluidboundaries, only: all_fluid_boundaries
    use fluidindex,      only: nvar, iarr_all_swpy, ibx, iby, ibz, nmag
    use grid,            only: dy, ny, is, ie, ks, ke, xdim, zdim, has_dir
    use mpisetup,        only: dt
    use rtvd,            only: relaxing_tvd
#ifdef COSM_RAYS
    use crhelpers,       only: div_v
#endif /* COSM_RAYS */

    implicit none
    real, dimension(nmag,ny)     :: b_y
    real, dimension(nvar%all,ny) :: u_y
    integer :: i,k,ip,kp

    b_y = 0.0
    u_y = 0.0

#ifdef COSM_RAYS
    call div_v(nvar%ion%pos)
#endif /* COSM_RAYS */

    do k=ks,ke
      kp=k+1
      do i=is,ie
          ip=i+1

#ifdef MAGNETIC
          b_y(:,:) = 0.5*b(:,i,:,k)
          b_y(iby,1:ny-1)=b_y(iby,1:ny-1)+b_y(iby,2:ny);       b_y(iby,ny) = b_y(iby,ny-1)
          if (has_dir(xdim) .and. i <= ie) b_y(ibx,:)=b_y(ibx,:)+0.5*b(ibx,ip,:,k)
          if (has_dir(zdim) .and. k <= ke) b_y(ibz,:)=b_y(ibz,:)+0.5*b(ibz,i,:,kp)
          b_y((/iby,ibx,ibz/),:)=b_y(:,:)
#endif /* MAGNETIC */

        u_y(iarr_all_swpy,:)=u(:,i,:,k)

        call relaxing_tvd(ny, u_y, b_y, 'ysweep', k, i, dy, dt)
        u(:,i,:,k)=u_y(iarr_all_swpy,:)

      enddo
    enddo

    call all_fluid_boundaries

  end subroutine sweepy
!------------------------------------------------------------------------------------------
  subroutine sweepz
    use arrays,          only: u, b
    use fluidboundaries, only: all_fluid_boundaries
    use fluidindex,      only: nvar, iarr_all_swpz, ibx, iby, ibz, nmag
    use grid,            only: dz, nz, is, ie, js, je, xdim, ydim, has_dir
    use mpisetup,        only: dt
    use rtvd,            only: relaxing_tvd
#ifdef COSM_RAYS
    use crhelpers,       only: div_v
#endif /* COSM_RAYS */

    implicit none
    real, dimension(nmag,nz)     :: b_z
    real, dimension(nvar%all,nz) :: u_z

    integer :: i,j,ip,jp

    b_z = 0.0
    u_z = 0.0

#ifdef COSM_RAYS
    call div_v(nvar%ion%pos)
#endif /* COSM_RAYS */

    do j=js,je
      jp=j+1
      do i=is,ie
          ip=i+1

#ifdef MAGNETIC
          b_z(:,:) = 0.5*b(:,i,j,:)
          b_z(ibz,1:nz-1) = b_z(ibz,1:nz-1) + b_z(ibz,2:nz);   b_z(ibz,nz) = b_z(ibz,nz-1)
          if (has_dir(xdim) .and. i <= ie) b_z(ibx,:) = b_z(ibx,:) + 0.5*b(ibx,ip,j,:)
          if (has_dir(ydim) .and. j <= je) b_z(iby,:) = b_z(iby,:) + 0.5*b(iby,i,jp,:)
          b_z((/ibz,iby,ibx/),:)=b_z(:,:)
#endif /* MAGNETIC */

        u_z(iarr_all_swpz,:)=u(:,i,j,:)

        call relaxing_tvd(nz, u_z, b_z, 'zsweep', i, j, dz, dt)
        u(:,i,j,:)=u_z(iarr_all_swpz,:)
      enddo
    enddo

    call all_fluid_boundaries

  end subroutine sweepz
end module sweeps
