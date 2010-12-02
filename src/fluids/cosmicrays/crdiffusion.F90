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

!>
!! \brief (MH) Numerical scheme for the diffusive transport of Cosmic Rays
!!
!!
!<

module crdiffusion
! pulled by COSM_RAYS
   implicit none

   private
   public :: cr_diff_x, cr_diff_y, cr_diff_z

contains

!****************************************************************************
!
   subroutine cr_diff_x
!
!  PURPOSE:  Diffusive transport of ecr in 1-direction
!
!-----------------------------------------------------------------------
      use arrays,         only: b, u, wcr
      use fluidindex,     only: ibx, iby, ibz, nvar
      use grid,           only: idx, idy, idz, xdim, ydim, zdim, has_dir, nx, js, je, ks, ke
      use initcosmicrays, only: iarr_crs, K_crs_paral, K_crs_perp
      use mpisetup,       only: dt

      implicit none

      integer :: i, j, k
      real    :: b1b, b2b, b3b, bb
      real, dimension(nvar%crs%all) :: decr1, decr2, decr3, fcrdif1
      real, dimension(nvar%crs%all) :: dqp, dqm

      !=======================================================================

      if (.not.has_dir(xdim)) return

      wcr(:, 1, :, :) = 0.
      wcr(:, :, :, :ks-1) = 0.
      wcr(:, :, :, ke+1:) = 0.
      wcr(:, :, :js-1, :) = 0.
      wcr(:, :, je+1:, :) = 0.

      do k=ks,ke
         do j=js,je
            do i=2,nx     ! if we are here this implies nxb /= 1

               decr1 =  (u(iarr_crs,i,  j,  k) - u(iarr_crs,i-1,j,  k)) * idx
               fcrdif1 = K_crs_perp * decr1

               b1b =  b(ibx,i,  j,  k)

               if (has_dir(ydim)) then
                  dqm = 0.5*((u(iarr_crs,i-1,j  ,k ) + u(iarr_crs,i ,j  ,k )) - (u(iarr_crs,i-1,j-1,k ) + u(iarr_crs,i ,j-1,k ))) * idy
                  dqp = 0.5*((u(iarr_crs,i-1,j+1,k ) + u(iarr_crs,i ,j+1,k )) - (u(iarr_crs,i-1,j  ,k ) + u(iarr_crs,i ,j  ,k ))) * idy
                  decr2 = (dqp+dqm)* (1.0 + sign(1.0, dqm*dqp))*0.25
                  b2b = sum(b(iby,i-1:i, j:j+1, k    ))*0.25
               else
                  decr2 = 0.
                  b2b = 0.
               endif

               if (has_dir(zdim)) then
                  dqm = 0.5*((u(iarr_crs,i-1,j ,k  ) + u(iarr_crs,i ,j ,k  )) - (u(iarr_crs,i-1,j ,k-1) + u(iarr_crs,i ,j ,k-1))) * idz
                  dqp = 0.5*((u(iarr_crs,i-1,j ,k+1) + u(iarr_crs,i ,j ,k+1)) - (u(iarr_crs,i-1,j ,k  ) + u(iarr_crs,i ,j ,k  ))) * idz
                  decr3 = (dqp+dqm)* (1.0 + sign(1.0, dqm*dqp))*0.25
                  b3b = sum(b(ibz,i-1:i, j,     k:k+1))*0.25
               else
                  decr3 = 0.
                  b3b = 0.
               endif

               bb = b1b**2 + b2b**2 + b3b**2
               if (bb /= 0.) fcrdif1 = fcrdif1 + K_crs_paral * b1b * (b1b*decr1 + b2b*decr2 + b3b*decr3) / bb

               wcr(:,i,j,k) = - fcrdif1 * dt * idx

            enddo
         enddo
      enddo

      u(iarr_crs,1:nx-1,:,:) = u(iarr_crs,1:nx-1,:,:) - ( wcr(:,2:nx,:,:) - wcr(:,1:nx-1,:,:) )
      u(iarr_crs,nx,:,:) = u(iarr_crs,nx-1,:,:) ! for sanity

   end subroutine cr_diff_x

!****************************************************************************
!
   subroutine cr_diff_y
!
!  PURPOSE:   Diffusive transport of ecr in 2-direction
!
!-----------------------------------------------------------------------
      use arrays,         only: b, u, wcr
      use fluidindex,     only: ibx, iby, ibz, nvar
      use grid,           only: idx, idy, idz, xdim, ydim, zdim, has_dir, ny, is, ie, ks, ke
      use initcosmicrays, only: iarr_crs, K_crs_paral, K_crs_perp
      use mpisetup,       only: dt

      implicit none

      integer :: i, j, k
      real    :: b1b, b2b, b3b, bb
      real, dimension(nvar%crs%all) :: decr1, decr2, decr3, fcrdif2
      real, dimension(nvar%crs%all) :: dqp, dqm

!=======================================================================

      if (.not.has_dir(ydim)) return

      wcr(:, :, 1, :) = 0.
      wcr(:, :is-1, :, :) = 0.
      wcr(:, ie+1:, :, :) = 0.
      wcr(:, :, :, :ks-1) = 0.
      wcr(:, :, :, ke+1:) = 0.

      do k=ks,ke
         do j=2,ny ! if we are here nyb /= 1
            do i=is,ie

               decr2 = (u(iarr_crs,i,j,k) - u(iarr_crs,i,j-1,k)) * idy
               fcrdif2 = K_crs_perp * decr2

               b2b =  b(iby,i,j,k)

               if (has_dir(xdim)) then
                  dqm = 0.5*((u(iarr_crs,i  ,j-1,k ) + u(iarr_crs,i  ,j  ,k )) - (u(iarr_crs,i-1,j-1,k ) + u(iarr_crs,i-1,j  ,k ))) * idx
                  dqp = 0.5*((u(iarr_crs,i+1,j-1,k ) + u(iarr_crs,i+1,j  ,k )) - (u(iarr_crs,i  ,j-1,k ) + u(iarr_crs,i  ,j  ,k ))) * idx
                  decr1 = (dqp+dqm)* (1.0 + sign(1.0, dqm*dqp))*0.25
                  b1b = sum(b(ibx,i:i+1, j-1:j, k))*0.25
               else
                  decr1 = 0.
                  b1b = 0.
               endif

               if (has_dir(zdim)) then
                  dqm = 0.5*((u(iarr_crs,i ,j-1,k  ) + u(iarr_crs,i ,j  ,k  )) - (u(iarr_crs,i ,j-1,k-1) + u(iarr_crs,i ,j  ,k-1))) * idz
                  dqp = 0.5*((u(iarr_crs,i ,j-1,k+1) + u(iarr_crs,i ,j  ,k+1)) - (u(iarr_crs,i ,j-1,k  ) + u(iarr_crs,i ,j  ,k  ))) * idz
                  decr3 = (dqp+dqm)* (1.0 + sign(1.0, dqm*dqp))*0.25
                  b3b = sum(b(ibz,i, j-1:j, k:k+1))*0.25
               else
                  decr3 = 0.
                  b3b = 0.
               endif

               bb = b1b**2 + b2b**2 + b3b**2
               if (bb /= 0.) fcrdif2 = fcrdif2 + K_crs_paral * b2b * (b1b*decr1 + b2b*decr2 + b3b*decr3) / bb

               wcr(:,i,j,k) = - fcrdif2 * dt * idy

            enddo
         enddo
      enddo

      u(iarr_crs,:,1:ny-1,:) = u(iarr_crs,:,1:ny-1,:) - ( wcr(:,:,2:ny,:) - wcr(:,:,1:ny-1,:) )
      u(iarr_crs,:,ny,:) = u(iarr_crs,:,ny-1,:) ! for sanity

   end subroutine cr_diff_y

!************************************************************************
!
   subroutine cr_diff_z
!
!  PURPOSE:   Diffusive transport of ecr in 3-direction
!
!-----------------------------------------------------------------------
      use arrays,         only: b, u, wcr
      use fluidindex,     only: ibx, iby, ibz, nvar
      use grid,           only: idx, idy, idz, xdim, ydim, zdim, has_dir, nz, is, ie, js, je
      use initcosmicrays, only: iarr_crs, K_crs_paral, K_crs_perp
      use mpisetup,       only: dt

      implicit none

      integer :: i, j, k
      real    :: b1b, b2b, b3b, bb
      real, dimension(nvar%crs%all) :: decr1, decr2, decr3, fcrdif3
      real, dimension(nvar%crs%all) :: dqp, dqm

!=======================================================================

      if (.not.has_dir(zdim)) return

      wcr(:, :, :, 1) = 0.
      wcr(:, :is-1, :, :) = 0.
      wcr(:, ie+1:, :, :) = 0.
      wcr(:, :, :js-1, :) = 0.
      wcr(:, :, je+1:, :) = 0.

      do k=2,nz      ! nzb /= 1
         do j=js,je
            do i=is,ie

               decr3 =  (u(iarr_crs,i,j,k    ) - u(iarr_crs,i,j,  k-1)) * idz
               fcrdif3 = K_crs_perp * decr3

               b3b =  b(ibz,i,j,k)

               if (has_dir(xdim)) then
                  dqm = 0.5*((u(iarr_crs,i  ,j ,k-1) + u(iarr_crs,i  ,j  ,k )) - (u(iarr_crs,i-1,j ,k-1) + u(iarr_crs,i-1,j  ,k ))) * idx
                  dqp = 0.5*((u(iarr_crs,i+1,j ,k-1) + u(iarr_crs,i+1,j  ,k )) - (u(iarr_crs,i  ,j ,k-1) + u(iarr_crs,i  ,j  ,k ))) * idx
                  decr1 = (dqp+dqm)* (1.0 + sign(1.0, dqm*dqp))*0.25
                  b1b = sum(b(ibx,i:i+1, j,     k-1:k))*0.25
               else
                  decr1 = 0.
                  b1b = 0.
               endif

               if (has_dir(ydim)) then
                  dqm = 0.5*((u(iarr_crs,i ,j  ,k-1) + u(iarr_crs,i  ,j  ,k )) - (u(iarr_crs,i ,j-1,k-1) + u(iarr_crs,i  ,j-1,k ))) * idy
                  dqp = 0.5*((u(iarr_crs,i ,j+1,k-1) + u(iarr_crs,i  ,j+1,k )) - (u(iarr_crs,i ,j  ,k-1) + u(iarr_crs,i  ,j  ,k ))) * idy
                  decr2 = (dqp+dqm)* (1.0 + sign(1.0, dqm*dqp))*0.25
                  b2b = sum(b(iby,i,     j:j+1, k-1:k))*0.25
               else
                  decr2 = 0.
                  b2b = 0.
               endif

               bb = b1b**2 + b2b**2 + b3b**2
               if (bb /= 0.) fcrdif3 = fcrdif3 + K_crs_paral * b3b * (b1b*decr1 + b2b*decr2 + b3b*decr3) / bb

               wcr(:,i,j,k) = - fcrdif3 * dt * idz

            enddo
         enddo
      enddo

      u(iarr_crs,:,:,1:nz-1) = u(iarr_crs,:,:,1:nz-1) - ( wcr(:,:,:,2:nz) - wcr(:,:,:,1:nz-1) )
      u(iarr_crs,:,:,nz) = u(iarr_crs,:,:,nz-1) ! for sanity

   end subroutine cr_diff_z

end module crdiffusion
