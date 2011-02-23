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
!<

module crdiffusion
! pulled by COSM_RAYS
   implicit none

   private
   public :: cr_diff_x, cr_diff_y, cr_diff_z

contains

!>
!! \brief Diffusive transport of ecr in 1-direction
!<
   subroutine cr_diff_x

      use arrays,         only: b, u, wcr
      use fluidindex,     only: ibx, iby, ibz, flind
      use grid,           only: cg
      use initcosmicrays, only: iarr_crs, K_crs_paral, K_crs_perp
      use mpisetup,       only: dt, has_dir
      use constants,      only: xdim, ydim, zdim

      implicit none

      integer :: i, j, k
      real    :: b1b, b2b, b3b, bb
      real, dimension(flind%crs%all) :: decr1, decr2, decr3, fcrdif1
      real, dimension(flind%crs%all) :: dqp, dqm

      !=======================================================================

      if (.not.has_dir(xdim)) return

      wcr(:, 1, :, :) = 0.
      wcr(:, :, :, :cg%ks-1) = 0.
      wcr(:, :, :, cg%ke+1:) = 0.
      wcr(:, :, :cg%js-1, :) = 0.
      wcr(:, :, cg%je+1:, :) = 0.

      do k=cg%ks, cg%ke
         do j=cg%js, cg%je
            do i=2, cg%nx     ! if we are here this implies nxb /= 1

               decr1 =  (u(iarr_crs,i,  j,  k) - u(iarr_crs,i-1,j,  k)) * cg%idx
               fcrdif1 = K_crs_perp * decr1

               b1b =  b(ibx,i,  j,  k)

               if (has_dir(ydim)) then
                  dqm = 0.5*((u(iarr_crs,i-1,j  ,k ) + u(iarr_crs,i ,j  ,k )) - (u(iarr_crs,i-1,j-1,k ) + u(iarr_crs,i ,j-1,k ))) * cg%idy
                  dqp = 0.5*((u(iarr_crs,i-1,j+1,k ) + u(iarr_crs,i ,j+1,k )) - (u(iarr_crs,i-1,j  ,k ) + u(iarr_crs,i ,j  ,k ))) * cg%idy
                  decr2 = (dqp+dqm)* (1.0 + sign(1.0, dqm*dqp))*0.25
                  b2b = sum(b(iby,i-1:i, j:j+1, k    ))*0.25
               else
                  decr2 = 0.
                  b2b = 0.
               endif

               if (has_dir(zdim)) then
                  dqm = 0.5*((u(iarr_crs,i-1,j ,k  ) + u(iarr_crs,i ,j ,k  )) - (u(iarr_crs,i-1,j ,k-1) + u(iarr_crs,i ,j ,k-1))) * cg%idz
                  dqp = 0.5*((u(iarr_crs,i-1,j ,k+1) + u(iarr_crs,i ,j ,k+1)) - (u(iarr_crs,i-1,j ,k  ) + u(iarr_crs,i ,j ,k  ))) * cg%idz
                  decr3 = (dqp+dqm)* (1.0 + sign(1.0, dqm*dqp))*0.25
                  b3b = sum(b(ibz,i-1:i, j,     k:k+1))*0.25
               else
                  decr3 = 0.
                  b3b = 0.
               endif

               bb = b1b**2 + b2b**2 + b3b**2
               if (bb /= 0.) fcrdif1 = fcrdif1 + K_crs_paral * b1b * (b1b*decr1 + b2b*decr2 + b3b*decr3) / bb

               wcr(:,i,j,k) = - fcrdif1 * dt * cg%idx

            enddo
         enddo
      enddo

      u(iarr_crs,1:cg%nx-1,:,:) = u(iarr_crs,1:cg%nx-1,:,:) - ( wcr(:,2:cg%nx,:,:) - wcr(:,1:cg%nx-1,:,:) )
      u(iarr_crs, cg%nx,:,:) = u(iarr_crs, cg%nx-1,:,:) ! for sanity

   end subroutine cr_diff_x

!>
!! \brief Diffusive transport of ecr in 2-direction
!<
   subroutine cr_diff_y

      use arrays,         only: b, u, wcr
      use fluidindex,     only: ibx, iby, ibz, flind
      use grid,           only: cg
      use initcosmicrays, only: iarr_crs, K_crs_paral, K_crs_perp
      use mpisetup,       only: dt, has_dir
      use constants,      only: xdim, ydim, zdim

      implicit none

      integer :: i, j, k
      real    :: b1b, b2b, b3b, bb
      real, dimension(flind%crs%all) :: decr1, decr2, decr3, fcrdif2
      real, dimension(flind%crs%all) :: dqp, dqm

!=======================================================================

      if (.not.has_dir(ydim)) return

      wcr(:, :, 1, :) = 0.
      wcr(:, :cg%is-1, :, :) = 0.
      wcr(:, cg%ie+1:, :, :) = 0.
      wcr(:, :, :, :cg%ks-1) = 0.
      wcr(:, :, :, cg%ke+1:) = 0.

      do k=cg%ks, cg%ke
         do j=2, cg%ny ! if we are here nyb /= 1
            do i=cg%is, cg%ie

               decr2 = (u(iarr_crs,i,j,k) - u(iarr_crs,i,j-1,k)) * cg%idy
               fcrdif2 = K_crs_perp * decr2

               b2b =  b(iby,i,j,k)

               if (has_dir(xdim)) then
                  dqm = 0.5*((u(iarr_crs,i  ,j-1,k ) + u(iarr_crs,i  ,j  ,k )) - (u(iarr_crs,i-1,j-1,k ) + u(iarr_crs,i-1,j  ,k ))) * cg%idx
                  dqp = 0.5*((u(iarr_crs,i+1,j-1,k ) + u(iarr_crs,i+1,j  ,k )) - (u(iarr_crs,i  ,j-1,k ) + u(iarr_crs,i  ,j  ,k ))) * cg%idx
                  decr1 = (dqp+dqm)* (1.0 + sign(1.0, dqm*dqp))*0.25
                  b1b = sum(b(ibx,i:i+1, j-1:j, k))*0.25
               else
                  decr1 = 0.
                  b1b = 0.
               endif

               if (has_dir(zdim)) then
                  dqm = 0.5*((u(iarr_crs,i ,j-1,k  ) + u(iarr_crs,i ,j  ,k  )) - (u(iarr_crs,i ,j-1,k-1) + u(iarr_crs,i ,j  ,k-1))) * cg%idz
                  dqp = 0.5*((u(iarr_crs,i ,j-1,k+1) + u(iarr_crs,i ,j  ,k+1)) - (u(iarr_crs,i ,j-1,k  ) + u(iarr_crs,i ,j  ,k  ))) * cg%idz
                  decr3 = (dqp+dqm)* (1.0 + sign(1.0, dqm*dqp))*0.25
                  b3b = sum(b(ibz,i, j-1:j, k:k+1))*0.25
               else
                  decr3 = 0.
                  b3b = 0.
               endif

               bb = b1b**2 + b2b**2 + b3b**2
               if (bb /= 0.) fcrdif2 = fcrdif2 + K_crs_paral * b2b * (b1b*decr1 + b2b*decr2 + b3b*decr3) / bb

               wcr(:,i,j,k) = - fcrdif2 * dt * cg%idy

            enddo
         enddo
      enddo

      u(iarr_crs,:,1:cg%ny-1,:) = u(iarr_crs,:,1:cg%ny-1,:) - ( wcr(:,:,2:cg%ny,:) - wcr(:,:,1:cg%ny-1,:) )
      u(iarr_crs,:, cg%ny,:) = u(iarr_crs,:, cg%ny-1,:) ! for sanity

   end subroutine cr_diff_y

!>
!! \brief Diffusive transport of ecr in 3-direction
!<
   subroutine cr_diff_z

      use arrays,         only: b, u, wcr
      use fluidindex,     only: ibx, iby, ibz, flind
      use grid,           only: cg
      use initcosmicrays, only: iarr_crs, K_crs_paral, K_crs_perp
      use mpisetup,       only: dt, has_dir
      use constants,      only: xdim, ydim, zdim

      implicit none

      integer :: i, j, k
      real    :: b1b, b2b, b3b, bb
      real, dimension(flind%crs%all) :: decr1, decr2, decr3, fcrdif3
      real, dimension(flind%crs%all) :: dqp, dqm

!=======================================================================

      if (.not.has_dir(zdim)) return

      wcr(:, :, :, 1) = 0.
      wcr(:, :cg%is-1, :, :) = 0.
      wcr(:, cg%ie+1:, :, :) = 0.
      wcr(:, :, :cg%js-1, :) = 0.
      wcr(:, :, cg%je+1:, :) = 0.

      do k=2, cg%nz      ! nzb /= 1
         do j=cg%js, cg%je
            do i=cg%is, cg%ie

               decr3 =  (u(iarr_crs,i,j,k    ) - u(iarr_crs,i,j,  k-1)) * cg%idz
               fcrdif3 = K_crs_perp * decr3

               b3b =  b(ibz,i,j,k)

               if (has_dir(xdim)) then
                  dqm = 0.5*((u(iarr_crs,i  ,j ,k-1) + u(iarr_crs,i  ,j  ,k )) - (u(iarr_crs,i-1,j ,k-1) + u(iarr_crs,i-1,j  ,k ))) * cg%idx
                  dqp = 0.5*((u(iarr_crs,i+1,j ,k-1) + u(iarr_crs,i+1,j  ,k )) - (u(iarr_crs,i  ,j ,k-1) + u(iarr_crs,i  ,j  ,k ))) * cg%idx
                  decr1 = (dqp+dqm)* (1.0 + sign(1.0, dqm*dqp))*0.25
                  b1b = sum(b(ibx,i:i+1, j,     k-1:k))*0.25
               else
                  decr1 = 0.
                  b1b = 0.
               endif

               if (has_dir(ydim)) then
                  dqm = 0.5*((u(iarr_crs,i ,j  ,k-1) + u(iarr_crs,i  ,j  ,k )) - (u(iarr_crs,i ,j-1,k-1) + u(iarr_crs,i  ,j-1,k ))) * cg%idy
                  dqp = 0.5*((u(iarr_crs,i ,j+1,k-1) + u(iarr_crs,i  ,j+1,k )) - (u(iarr_crs,i ,j  ,k-1) + u(iarr_crs,i  ,j  ,k ))) * cg%idy
                  decr2 = (dqp+dqm)* (1.0 + sign(1.0, dqm*dqp))*0.25
                  b2b = sum(b(iby,i,     j:j+1, k-1:k))*0.25
               else
                  decr2 = 0.
                  b2b = 0.
               endif

               bb = b1b**2 + b2b**2 + b3b**2
               if (bb /= 0.) fcrdif3 = fcrdif3 + K_crs_paral * b3b * (b1b*decr1 + b2b*decr2 + b3b*decr3) / bb

               wcr(:,i,j,k) = - fcrdif3 * dt * cg%idz

            enddo
         enddo
      enddo

      u(iarr_crs,:,:,1:cg%nz-1) = u(iarr_crs,:,:,1:cg%nz-1) - ( wcr(:,:,:,2:cg%nz) - wcr(:,:,:,1:cg%nz-1) )
      u(iarr_crs,:,:, cg%nz) = u(iarr_crs,:,:, cg%nz-1) ! for sanity

   end subroutine cr_diff_z

end module crdiffusion
