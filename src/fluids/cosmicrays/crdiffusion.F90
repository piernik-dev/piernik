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
#include "piernik.def"

!>
!! \brief (MH) Numerical scheme for the diffusive transport of Cosmic Rays
!!
!!
!<

module crdiffusion

  use arrays,         only: b, u, wcr
  use constants,      only: small
  use fluidindex,     only: ibx, iby, ibz, nvar
  use grid,           only: idx, idy, idz, nxd, nyd, nzd, nx, ny, nz, is, ie, js, je, ks, ke
  use initcosmicrays, only: iarr_crs, K_crs_paral, K_crs_perp
  use mpisetup,       only: dt

 contains

!
      subroutine cr_diff_x
!
!  PURPOSE:  Diffusive transport of ecr in 1-direction
!
      implicit none
      integer :: i, j, k
      real    :: b1b, b2b, b3b, n1b, n2b, n3b, bb
      real, dimension(nvar%crs%all) :: decr1, decr2, decr3, fcrdif1
      real, dimension(nvar%crs%all) :: dqp2, dqm2, dqp3, dqm3

!=======================================================================

      wcr(:, 1, :, :) = 0.
      wcr(:, :, :, :ks-1) = 0.
      wcr(:, :, :, ke+1:) = 0.
      wcr(:, :, :js-1, :) = 0.
      wcr(:, :, je+1:, :) = 0.

      do k=ks,ke
        do j=js,je
          do i=2,nx     ! if we are here this implies nxd /= 1

             b1b =  b(ibx,i,  j,  k)
             if (nyd /= 1) then
                b2b = (b(iby,i,  j,  k) + b(iby,i-1,j,  k)       &
                     + b(iby,i-1,j+1,k) + b(iby,i,  j+1,k))*0.25
             else
                b2b = 0.0
             endif
             if (nzd /= 1) then
                b3b = (b(ibz,i,  j,  k) + b(ibz,i-1,j,  k)       &
                     + b(ibz,i-1,j,k+1) + b(ibz,i,  j,k+1))*0.25
             else
                b3b = 0.0
             endif

             bb  = sqrt(b1b**2 + b2b**2 + b3b**2 + small)
             n1b = b1b/bb
             n2b = b2b/bb
             n3b = b3b/bb

             decr1 =  (u(iarr_crs,i,  j,  k) - u(iarr_crs,i-1,j,  k))*idx

             if (nyd /= 1) then
                dqm2  = 0.5*((u(iarr_crs,i-1,j  ,k ) + u(iarr_crs,i ,j  ,k ))    &
                     &      -(u(iarr_crs,i-1,j-1,k ) + u(iarr_crs,i ,j-1,k )))*idy
                dqp2  = 0.5*((u(iarr_crs,i-1,j+1,k ) + u(iarr_crs,i ,j+1,k ))    &
                     &      -(u(iarr_crs,i-1,j  ,k ) + u(iarr_crs,i ,j  ,k )))*idy

                decr2 = (dqp2+dqm2)* (1.0 + sign(1.0, dqm2*dqp2))*0.25
             else
                decr2 = 0.0
             endif

             if (nzd /= 1) then
                dqm3  = 0.5*((u(iarr_crs,i-1,j ,k  ) + u(iarr_crs,i ,j ,k  ))    &
                            -(u(iarr_crs,i-1,j ,k-1) + u(iarr_crs,i ,j ,k-1)))*idz
                dqp3  = 0.5*((u(iarr_crs,i-1,j ,k+1) + u(iarr_crs,i ,j ,k+1))    &
                            -(u(iarr_crs,i-1,j ,k  ) + u(iarr_crs,i ,j ,k  )))*idz

                decr3 = (dqp3+dqm3)* (1.0 + sign(1.0, dqm3*dqp3))*0.25
             else
                decr3 = 0.0
             endif

             fcrdif1 = K_crs_paral * n1b *   &
                   (n1b*decr1 + n2b*decr2 + n3b*decr3)  &
                   + K_crs_perp * decr1


             wcr(:,i,j,k) = - fcrdif1 * dt * idx

           enddo
         enddo
      enddo

!     do 60 k=ks,ke
!       do 50 j=js,je
!         do 40 i=is,ie
!            u(iarr_crs,i,j,k) = u(iarr_crs,i,j,k)  &
!                        -(wcr(:,i+1,j,k)-wcr(:,i,j,k))/dx
!40       continue
!50      continue
!60    continue
      u(iarr_crs,1:nx-1,:,:) = u(iarr_crs,1:nx-1,:,:) - ( wcr(:,2:nx,:,:) - wcr(:,1:nx-1,:,:) )
      u(iarr_crs,nx,:,:) = u(iarr_crs,nx-1,:,:) ! for sanity

      return
      end subroutine cr_diff_x


!****************************************************************************

!
      subroutine cr_diff_y
!
!  PURPOSE:   Diffusive transport of ecr in 2-direction
!
!-----------------------------------------------------------------------
      implicit none
      integer :: i, j, k
      real    :: b1b, b2b, b3b, n1b, n2b, n3b, bb
      real, dimension(nvar%crs%all) :: decr1, decr2, decr3, fcrdif2
      real, dimension(nvar%crs%all) :: dqp2, dqm1, dqm2, dqp1

!=======================================================================
!

      wcr(:, :, 1, :) = 0.
      wcr(:, :is-1, :, :) = 0.
      wcr(:, ie+1:, :, :) = 0.
      wcr(:, :, :, :ks-1) = 0.
      wcr(:, :, :, ke+1:) = 0.

      do k=ks,ke
        do j=2,ny ! if we are here nyd /= 1
          do i=is,ie

             if (nxd /= 1) then
                b1b = (b(ibx,i,j,k) + b(ibx,i,j-1,k)   &
                     + b(ibx,i+1,j-1,k) + b(ibx,i+1,j,k))*0.25
             else
                b1b = 0.0
             endif

             b2b =  b(iby,i,j,k)

             if (nzd /= 1) then
                b3b = (b(ibz,i,j,k) + b(ibz,i,j-1,k)   &
                     + b(ibz,i,j-1,k+1) + b(ibz,i,j,k+1))*0.25
             else
                b3b = 0.0
             endif

             bb  = sqrt(b1b**2 + b2b**2 + b3b**2 + small)
             n1b = b1b/bb
             n2b = b2b/bb
             n3b = b3b/bb

             if (nxd /= 1) then
                dqm1  = 0.5*((u(iarr_crs,i  ,j-1,k ) + u(iarr_crs,i  ,j  ,k ))   &
                            -(u(iarr_crs,i-1,j-1,k ) + u(iarr_crs,i-1,j  ,k )))*idx
                dqp1  = 0.5*((u(iarr_crs,i+1,j-1,k ) + u(iarr_crs,i+1,j  ,k ))   &
                            -(u(iarr_crs,i  ,j-1,k ) + u(iarr_crs,i  ,j  ,k )))*idx

                decr1 = (dqp1+dqm1)* (1.0 + sign(1.0, dqm1*dqp1))*0.25
             else
                decr1 = 0.0
             endif

             decr2 = (u(iarr_crs,i,j,k) - u(iarr_crs,i,j-1,k))*idy

             if (nzd /= 1) then
                dqm2  = 0.5*((u(iarr_crs,i ,j-1,k  ) + u(iarr_crs,i ,j  ,k  ))   &
                            -(u(iarr_crs,i ,j-1,k-1) + u(iarr_crs,i ,j  ,k-1)))*idz
                dqp2  = 0.5*((u(iarr_crs,i ,j-1,k+1) + u(iarr_crs,i ,j  ,k+1))   &
                            -(u(iarr_crs,i ,j-1,k  ) + u(iarr_crs,i ,j  ,k  )))*idz

                decr3 = (dqp2+dqm2)* (1.0 + sign(1.0, dqm2*dqp2))*0.25
             else
                decr3 = 0.0
             endif


             fcrdif2 = K_crs_paral * n2b * &
                      (n1b*decr1 + n2b*decr2 + n3b*decr3) &
                      + K_crs_perp * decr2

             wcr(:,i,j,k) = - fcrdif2 * dt * idy

            enddo
         enddo
      enddo

!      do 60 k=ks,ke
!        do 50 j=js,je
!          do 40 i=is,ie
!            u(iarr_crs,i,j,k) = u(iarr_crs,i,j,k)    &
!                      - (wcr(:,i,j+1,k)-wcr(:,i,j,k))/dy
!40        continue
!50      continue
!60    continue

      u(iarr_crs,:,1:ny-1,:) = u(iarr_crs,:,1:ny-1,:) - ( wcr(:,:,2:ny,:) - wcr(:,:,1:ny-1,:) )
      u(iarr_crs,:,ny,:) = u(iarr_crs,:,ny-1,:) ! for sanity

      return
      end subroutine cr_diff_y

!************************************************************************

      subroutine cr_diff_z
!
!  PURPOSE:   Diffusive transport of ecr in 3-direction
!
      implicit none
      integer :: i,j,k
      real    :: b1b, b2b, b3b, n1b, n2b, n3b, bb
      real, dimension(nvar%crs%all) :: decr1, decr2, decr3, fcrdif3
      real, dimension(nvar%crs%all) :: dqp2, dqm2, dqm1, dqp1

!=======================================================================
!

      wcr(:, :, :, 1) = 0.
      wcr(:, :is-1, :, :) = 0.
      wcr(:, ie+1:, :, :) = 0.
      wcr(:, :, :js-1, :) = 0.
      wcr(:, :, je+1:, :) = 0.

      do k=2,nz      ! nzd /= 1
        do j=js,je
          do i=is,ie

             if (nxd /= 1) then
                b1b = (b(ibx,i,  j,k  ) + b(ibx,i,  j,k-1)  &
                     + b(ibx,i+1,j,k-1) + b(ibx,i+1,j,k  ))*0.25
             else
                b1b = 0.0
             endif
             if (nyd /= 1) then
                b2b = (b(iby,i,j,  k  ) + b(iby,i,j,  k-1)  &
                     + b(iby,i,j+1,k-1) + b(iby,i,j+1,k  ))*0.25
             else
                b2b = 0.0
             endif
             b3b =  b(ibz,i,j,k)

             bb  = sqrt(b1b**2 + b2b**2 + b3b**2 + small)
             n1b = b1b/bb
             n2b = b2b/bb
             n3b = b3b/bb

             if (nxd /= 1) then
                dqm1  = 0.5*((u(iarr_crs,i  ,j ,k-1) + u(iarr_crs,i  ,j  ,k ))   &
                            -(u(iarr_crs,i-1,j ,k-1) + u(iarr_crs,i-1,j  ,k )))*idx
                dqp1  = 0.5*((u(iarr_crs,i+1,j ,k-1) + u(iarr_crs,i+1,j  ,k ))   &
                            -(u(iarr_crs,i  ,j ,k-1) + u(iarr_crs,i  ,j  ,k )))*idx

                decr1 = (dqp1+dqm1)* (1.0 + sign(1.0, dqm1*dqp1))*0.25
             else
                decr1 = 0.0
             endif

             if (nyd /= 1) then
                dqm2  = 0.5*((u(iarr_crs,i ,j  ,k-1) + u(iarr_crs,i  ,j  ,k ))   &
                            -(u(iarr_crs,i ,j-1,k-1) + u(iarr_crs,i  ,j-1,k )))*idy
                dqp2  = 0.5*((u(iarr_crs,i ,j+1,k-1) + u(iarr_crs,i  ,j+1,k ))   &
                            -(u(iarr_crs,i ,j  ,k-1) + u(iarr_crs,i  ,j  ,k )))*idy

                decr2 = (dqp2+dqm2)* (1.0 + sign(1.0, dqm2*dqp2))*0.25
             else
                decr2 = 0.0
             endif

             decr3 =  (u(iarr_crs,i,j,k    ) - u(iarr_crs,i,j,  k-1))*idz

             fcrdif3 = K_crs_paral * n3b * &
                  (n1b*decr1 + n2b*decr2 + n3b*decr3) &
                   + K_crs_perp * decr3

             wcr(:,i,j,k) = - fcrdif3 * dt * idz
          enddo
        enddo
      enddo

!      do 60 k=ks,ke
!        do 50 j=js,je
!          do 40 i=is,ie
!            u(iarr_crs,i,j,k) = u(iarr_crs,i,j,k)  &
!                        -(wcr(:,i,j,k+1)-wcr(:,i,j,k))/dz
!40        continue
!50      continue
!60    continue

      u(iarr_crs,:,:,1:nz-1) = u(iarr_crs,:,:,1:nz-1) - ( wcr(:,:,:,2:nz) - wcr(:,:,:,1:nz-1) )
      u(iarr_crs,:,:,nz) = u(iarr_crs,:,:,nz-1) ! for sanity

      return
      end subroutine cr_diff_z


end module crdiffusion
