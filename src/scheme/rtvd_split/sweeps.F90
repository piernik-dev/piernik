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
      use fluidboundaries, only: all_fluid_boundaries
      use fluidindex,      only: iarr_all_dn, iarr_all_mx, iarr_all_my, flind
      use grid,            only: cg
      use mpisetup,        only: dt
      use shear,           only: omega, qshear
      use interactions,    only: dragc_gas_dust

      implicit none
      real, dimension(size(iarr_all_my), cg%nx, cg%nz) :: vxr, v_r, rotaccr
      real, dimension(cg%nx, cg%nz)      :: epsa
      real, dimension(flind%all, cg%nx, cg%nz) :: u1
      integer :: ind,i
      real, parameter, dimension(2) :: fac = [0.5, 1.0]

      u1(:,:,:) = cg%u%arr(:,:,1,:)

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
!           rotaccr(ind,:,:) = rotaccr(ind,:,:) - 2.0*omega*vxr(ind,:,:)
            rotaccr(ind,:,:) = rotaccr(ind,:,:) + (qshear-2.0)*omega*vxr(ind,:,:)
         enddo

         where (u1(iarr_all_dn,:,:) > 0.0)
            u1(iarr_all_my,:,:) = u1(iarr_all_my,:,:) + fac(i)*dt*rotaccr(:,:,:)*cg%u%arr(iarr_all_dn,:,1,:)
         endwhere
         cg%u%arr(iarr_all_my,:,1,:) = u1(iarr_all_my,:,:)
      enddo

      call all_fluid_boundaries

   end subroutine source_terms_y
#endif /* SHEAR && FLUID_INTERACTIONS */
!------------------------------------------------------------------------------------------
   subroutine sweepx

      use arrays,          only: uh
      use fluidboundaries, only: all_fluid_boundaries
      use fluidindex,      only: flind, iarr_all_swpx, ibx, iby, ibz, nmag
      use grid,            only: cg
      use gridgeometry,    only: set_geo_coeffs
      use mpisetup,        only: dt, has_dir, integration_order
      use constants,       only: xdim, ydim, zdim
      use rtvd,            only: relaxing_tvd
#ifdef COSM_RAYS
      use crhelpers,       only: div_v, set_div_v1d
#endif /* COSM_RAYS */

      implicit none
      real, dimension(nmag, cg%nx)      :: b_x
      real, dimension(flind%all, cg%nx) :: u_x, u0_x
      real, dimension(:), pointer       :: div_v1d => null()
      integer                           :: j, k, jp, kp, istep

      b_x = 0.0
      u_x = 0.0

#ifdef COSM_RAYS
      call div_v(flind%ion%pos)
#endif /* COSM_RAYS */
      uh = cg%u%arr
      do istep = 1, integration_order
         do k=cg%ks, cg%ke
            kp=k+1
            do j=cg%js, cg%je
               jp=j+1

#ifdef MAGNETIC
               b_x=0.5*cg%b%arr(:,:,j,k)
               b_x(ibx,1:cg%nx-1) = b_x(ibx,1:cg%nx-1)+b_x(ibx,2:cg%nx);       b_x(ibx, cg%nx) = b_x(ibx, cg%nx-1)
               if (has_dir(ydim) .and. j <= cg%je)  b_x(iby,:)=b_x(iby,:)+0.5*cg%b%arr(iby,:,jp,k)
               if (has_dir(zdim) .and. k <= cg%ke)  b_x(ibz,:)=b_x(ibz,:)+0.5*cg%b%arr(ibz,:,j,kp)
#endif /* MAGNETIC */

               call set_geo_coeffs(xdim,flind,j,k)
#ifdef COSM_RAYS
               call set_div_v1d(div_v1d,xdim,j,k)
#endif /* COSM_RAYS */

               u_x (iarr_all_swpx,:) = cg%u%arr(:,:,j,k)
               u0_x(iarr_all_swpx,:) = uh(:,:,j,k)

               call relaxing_tvd(cg%nx, u_x, u0_x, b_x, div_v1d, istep, xdim, j, k, cg%dx, dt)
               cg%u%arr(:,:,j,k)=u_x(iarr_all_swpx,:)
            enddo
         enddo
         call all_fluid_boundaries    ! \todo : call only x for istep=1, call all for istep=2
      enddo

   end subroutine sweepx
!------------------------------------------------------------------------------------------
   subroutine sweepy
      use arrays,          only: uh
      use fluidboundaries, only: all_fluid_boundaries
      use fluidindex,      only: flind, iarr_all_swpy, ibx, iby, ibz, nmag
      use grid,            only: cg
      use gridgeometry,    only: set_geo_coeffs
      use mpisetup,        only: dt, has_dir, integration_order
      use constants,       only: xdim, ydim, zdim
      use rtvd,            only: relaxing_tvd
#ifdef COSM_RAYS
      use crhelpers,       only: div_v, set_div_v1d
#endif /* COSM_RAYS */

      implicit none
      real, dimension(nmag, cg%ny)      :: b_y
      real, dimension(flind%all, cg%ny) :: u_y, u0_y
      real, dimension(:), pointer       :: div_v1d => null()
      integer                           :: i, k, ip, kp, istep
      b_y = 0.0
      u_y = 0.0

#ifdef COSM_RAYS
      call div_v(flind%ion%pos)
#endif /* COSM_RAYS */
      uh = cg%u%arr
      do istep = 1, integration_order
         do k=cg%ks, cg%ke
            kp=k+1
            do i=cg%is, cg%ie
               ip=i+1

#ifdef MAGNETIC
               b_y(:,:) = 0.5*cg%b%arr(:,i,:,k)
               b_y(iby,1:cg%ny-1)=b_y(iby,1:cg%ny-1)+b_y(iby,2:cg%ny);       b_y(iby, cg%ny) = b_y(iby, cg%ny-1)
               if (has_dir(xdim) .and. i <= cg%ie) b_y(ibx,:)=b_y(ibx,:)+0.5*cg%b%arr(ibx,ip,:,k)
               if (has_dir(zdim) .and. k <= cg%ke) b_y(ibz,:)=b_y(ibz,:)+0.5*cg%b%arr(ibz,i,:,kp)
               b_y((/iby,ibx,ibz/),:)=b_y(:,:)
#endif /* MAGNETIC */

               call set_geo_coeffs(ydim,flind,k,i)
#ifdef COSM_RAYS
               call set_div_v1d(div_v1d,ydim,k,i)
#endif /* COSM_RAYS */

               u_y (iarr_all_swpy,:) = cg%u%arr(:,i,:,k)
               u0_y(iarr_all_swpy,:) = uh(:,i,:,k)

               call relaxing_tvd(cg%ny, u_y, u0_y, b_y, div_v1d, istep, ydim, k, i, cg%dy, dt)
               cg%u%arr(:,i,:,k)=u_y(iarr_all_swpy,:)

            enddo
         enddo

         call all_fluid_boundaries
      enddo

   end subroutine sweepy
!------------------------------------------------------------------------------------------
   subroutine sweepz
      use arrays,          only: uh
      use fluidboundaries, only: all_fluid_boundaries
      use fluidindex,      only: flind, iarr_all_swpz, ibx, iby, ibz, nmag
      use grid,            only: cg
      use gridgeometry,    only: set_geo_coeffs
      use mpisetup,        only: dt, has_dir, integration_order
      use constants,       only: xdim, ydim, zdim
      use rtvd,            only: relaxing_tvd
#ifdef COSM_RAYS
      use crhelpers,       only: div_v, set_div_v1d
#endif /* COSM_RAYS */

      implicit none
      real, dimension(nmag, cg%nz)      :: b_z
      real, dimension(flind%all, cg%nz) :: u_z, u0_z
      real, dimension(:), pointer       :: div_v1d => null()
      integer                           :: i, j, ip, jp, istep

      b_z = 0.0
      u_z = 0.0

#ifdef COSM_RAYS
      call div_v(flind%ion%pos)
#endif /* COSM_RAYS */
      uh = cg%u%arr
      do istep = 1, integration_order
         do j=cg%js, cg%je
            jp=j+1
            do i=cg%is, cg%ie
               ip=i+1

#ifdef MAGNETIC
               b_z(:,:) = 0.5*cg%b%arr(:,i,j,:)
               b_z(ibz,1:cg%nz-1) = b_z(ibz,1:cg%nz-1) + b_z(ibz,2:cg%nz);   b_z(ibz, cg%nz) = b_z(ibz, cg%nz-1)
               if (has_dir(xdim) .and. i <= cg%ie) b_z(ibx,:) = b_z(ibx,:) + 0.5*cg%b%arr(ibx,ip,j,:)
               if (has_dir(ydim) .and. j <= cg%je) b_z(iby,:) = b_z(iby,:) + 0.5*cg%b%arr(iby,i,jp,:)
               b_z((/ibz,iby,ibx/),:)=b_z(:,:)
#endif /* MAGNETIC */

               call set_geo_coeffs(zdim,flind,i,j)
#ifdef COSM_RAYS
               call set_div_v1d(div_v1d,zdim,i,j)
#endif /* COSM_RAYS */

               !OPT: It looks that u_z gets re-assigned to something inside relaxing_tvd. \todo try to merge these assignments
               !OPT: 3% D1mr, 3% D2mr, 20% D1mw, Ir:Dr:Dw ~ 10:4:1
               !OPT: The same applies to sweepy and sweepz
               u_z (iarr_all_swpz,:) = cg%u%arr(:,i,j,:)
               u0_z(iarr_all_swpz,:) = uh(:,i,j,:)

               call relaxing_tvd(cg%nz, u_z, u0_z, b_z, div_v1d, istep, zdim, i, j, cg%dz, dt)
               cg%u%arr(:,i,j,:)=u_z(iarr_all_swpz,:)
            enddo
         enddo

         call all_fluid_boundaries
      enddo

   end subroutine sweepz
!------------------------------------------------------------------------------------------
end module sweeps
