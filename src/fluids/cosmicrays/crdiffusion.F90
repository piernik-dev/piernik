! $Id: cr_diffusion.F90 594 2009-01-21 12:33:44Z xarth $
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


!>
!! \brief (MH) Numerical scheme for the diffusive transport of Cosmic Rays
!!
!!
!<

module crdiffusion

  use initcosmicrays, only : iecr
  use initcosmicrays, only : K_cr_paral,K_cr_perp

  use fluidindex,     only : ibx,iby,ibz
  use arrays,         only : b,u,wa
  use constants,      only : small
  use grid,           only : dx,dy,dz,nxd,nyd,nzd,nx,ny,nz
  use grid,           only : is,ie,js,je,ks,ke
  use mpisetup,       only : dt


 contains


!
      subroutine cr_diff_x
!
!  PURPOSE:  Diffusive transport of ecr in 1-direction
!
      implicit none
      integer i,j,k, n
      real :: b1b, b2b, b3b, n1b, n2b, n3b, bb
      real :: decr1, decr2, decr3, fcrdif1
      real :: dqp2, dqm2, dqp3, dqm3

!=======================================================================

      do k=ks,ke
        do j=js,je
          do i=2,nx     ! if we are here this implies nxd /= 1

             b1b =  b(ibx,i,  j,  k)
             if(nyd /= 1) then
                b2b = (b(iby,i,  j,  k) + b(iby,i-1,j,  k)       &
                     + b(iby,i-1,j+1,k) + b(iby,i,  j+1,k))*0.25
             else
                b2b = 0.0
             endif
             if(nzd /= 1) then
                b3b = (b(ibz,i,  j,  k) + b(ibz,i-1,j,  k)       &
                     + b(ibz,i-1,j,k+1) + b(ibz,i,  j,k+1))*0.25
             else
                b3b = 0.0
             endif

             bb  = sqrt(b1b**2 + b2b**2 + b3b**2 + small)
             n1b = b1b/bb
             n2b = b2b/bb
             n3b = b3b/bb

             decr1 =  (u(iecr,i,  j,  k) - u(iecr,i-1,j,  k))/dx

             dqm2  = 0.5*((u(iecr,i-1,j  ,k ) + u(iecr,i ,j  ,k ))    &
                         -(u(iecr,i-1,j-1,k ) + u(iecr,i ,j-1,k )))/dy
             dqp2  = 0.5*((u(iecr,i-1,j+1,k ) + u(iecr,i ,j+1,k ))    &
                         -(u(iecr,i-1,j  ,k ) + u(iecr,i ,j  ,k )))/dy

             decr2 = (dqp2+dqm2)* (1.0 + sign(1.0, dqm2*dqp2))*0.25

             if(nzd /= 1) then
                dqm3  = 0.5*((u(iecr,i-1,j ,k  ) + u(iecr,i ,j ,k  ))    &
                            -(u(iecr,i-1,j ,k-1) + u(iecr,i ,j ,k-1)))/dz
                dqp3  = 0.5*((u(iecr,i-1,j ,k+1) + u(iecr,i ,j ,k+1))    &
                            -(u(iecr,i-1,j ,k  ) + u(iecr,i ,j ,k  )))/dz

                decr3 = (dqp3+dqm3)* (1.0 + sign(1.0, dqm3*dqp3))*0.25
             else
                decr3 = 0.0
             endif

             fcrdif1 = K_cr_paral * n1b *   &
                   (n1b*decr1 + n2b*decr2 + n3b*decr3)  &
                   + K_cr_perp * decr1


             wa(i,j,k) = - fcrdif1 * dt / dx

           enddo
         enddo
      enddo

!     do 60 k=ks,ke
!       do 50 j=js,je
!         do 40 i=is,ie
!            u(iecr,i,j,k) = u(iecr,i,j,k)  &
!                        -(wa(i+1,j,k)-wa(i,j,k))/dx
!40       continue
!50      continue
!60    continue
      u(iecr,1:nx-1,:,:) = u(iecr,1:nx-1,:,:) - ( wa(2:nx,:,:) - wa(1:nx-1,:,:) )
      u(iecr,nx,:,:) = u(iecr,nx-1,:,:) ! for sanity

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
      integer i,j,k,n
      real :: b1b, b2b, b3b, n1b, n2b, n3b, bb
      real :: decr1, decr2, decr3, fcrdif2
      real :: dqp2, dqm1, dqm2, dqp3, dqm3, dqp1

!=======================================================================
!

      do k=ks,ke
        do j=2,ny ! if we are here nyd /= 1
          do i=is,ie

             if(nxd /= 1) then
                b1b = (b(ibx,i,j,k) + b(ibx,i,j-1,k)   &
                     + b(ibx,i+1,j-1,k) + b(ibx,i+1,j,k))*0.25
             else
                b1b = 0.0
             endif

             b2b =  b(iby,i,j,k)

             if(nzd /= 1) then
                b3b = (b(ibz,i,j,k) + b(ibz,i,j-1,k)   &
                     + b(ibz,i,j-1,k+1) + b(ibz,i,j,k+1))*0.25
             else
                b3b = 0.0
             endif

             bb  = sqrt(b1b**2 + b2b**2 + b3b**2 + small)
             n1b = b1b/bb
             n2b = b2b/bb
             n3b = b3b/bb

             if(nxd /= 1) then
                dqm1  = 0.5*((u(iecr,i  ,j-1,k ) + u(iecr,i  ,j  ,k ))   &
                            -(u(iecr,i-1,j-1,k ) + u(iecr,i-1,j  ,k )))/dx
                dqp1  = 0.5*((u(iecr,i+1,j-1,k ) + u(iecr,i+1,j  ,k ))   &
                            -(u(iecr,i  ,j-1,k ) + u(iecr,i  ,j  ,k )))/dx

                decr1 = (dqp1+dqm1)* (1.0 + sign(1.0, dqm1*dqp1))*0.25
             else
                decr1 = 0.0
             endif

             decr2 = (u(iecr,i,j,k) - u(iecr,i,j-1,k))/dy

             if(nzd /= 1) then
                dqm2  = 0.5*((u(iecr,i ,j-1,k  ) + u(iecr,i ,j  ,k  ))   &
                            -(u(iecr,i ,j-1,k-1) + u(iecr,i ,j  ,k-1)))/dz
                dqp2  = 0.5*((u(iecr,i ,j-1,k+1) + u(iecr,i ,j  ,k+1))   &
                            -(u(iecr,i ,j-1,k  ) + u(iecr,i ,j  ,k  )))/dz

                decr3 = (dqp2+dqm2)* (1.0 + sign(1.0, dqm2*dqp2))*0.25
             else
                decr3 = 0.0
             endif


             fcrdif2 = K_cr_paral * n2b * &
                      (n1b*decr1 + n2b*decr2 + n3b*decr3) &
                      + K_cr_perp * decr2

             wa(i,j,k) = - fcrdif2 * dt / dy

            enddo
         enddo
      enddo

!      do 60 k=ks,ke
!        do 50 j=js,je
!          do 40 i=is,ie
!            u(iecr,i,j,k) = u(iecr,i,j,k)    &
!                      - (wa(i,j+1,k)-wa(i,j,k))/dy
!40        continue
!50      continue
!60    continue

      u(iecr,:,1:ny-1,:) = u(iecr,:,1:ny-1,:) - ( wa(:,2:ny,:) - wa(:,1:ny-1,:) )
      u(iecr,:,ny,:) = u(iecr,:,ny-1,:) ! for sanity

      return
      end subroutine cr_diff_y

!************************************************************************

      subroutine cr_diff_z
!
!  PURPOSE:   Diffusive transport of ecr in 3-direction
!
      implicit none
      integer i,j,k,n
      real :: b1b, b2b, b3b, n1b, n2b, n3b, bb
      real :: decr1, decr2, decr3, fcrdif3
      real :: dqp2, dqm2, dqp3, dqm3, dqm1, dqp1

!=======================================================================
!
      do k=2,nz      ! nzd /= 1
        do j=js,je
          do i=is,ie

             if(nxd /= 1) then
                b1b = (b(ibx,i,  j,k  ) + b(ibx,i,  j,k-1)  &
                     + b(ibx,i+1,j,k-1) + b(ibx,i+1,j,k  ))*0.25
             else
                b1b = 0.0
             endif
             if(nyd /= 1) then
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

             if(nxd /= 1) then
                dqm1  = 0.5*((u(iecr,i  ,j ,k-1) + u(iecr,i  ,j  ,k ))   &
                            -(u(iecr,i-1,j ,k-1) + u(iecr,i-1,j  ,k )))/dx
                dqp1  = 0.5*((u(iecr,i+1,j ,k-1) + u(iecr,i+1,j  ,k ))   &
                            -(u(iecr,i  ,j ,k-1) + u(iecr,i  ,j  ,k )))/dx

                decr1 = (dqp1+dqm1)* (1.0 + sign(1.0, dqm1*dqp1))*0.25
             else
                decr1 = 0.0
             endif

             if(nyd /= 1) then
                dqm2  = 0.5*((u(iecr,i ,j  ,k-1) + u(iecr,i  ,j  ,k ))   &
                            -(u(iecr,i ,j-1,k-1) + u(iecr,i  ,j-1,k )))/dy
                dqp2  = 0.5*((u(iecr,i ,j+1,k-1) + u(iecr,i  ,j+1,k ))   &
                            -(u(iecr,i ,j  ,k-1) + u(iecr,i  ,j  ,k )))/dy

                decr2 = (dqp2+dqm2)* (1.0 + sign(1.0, dqm2*dqp2))*0.25
             else
                decr2 = 0.0
             endif

             decr3 =  (u(iecr,i,j,k    ) - u(iecr,i,j,  k-1))/dz

             fcrdif3 = K_cr_paral * n3b * &
                  (n1b*decr1 + n2b*decr2 + n3b*decr3) &
                   + K_cr_perp * decr3

             wa(i,j,k) = - fcrdif3 * dt / dz
          enddo
        enddo
      enddo

!      do 60 k=ks,ke
!        do 50 j=js,je
!          do 40 i=is,ie
!            u(iecr,i,j,k) = u(iecr,i,j,k)  &
!                        -(wa(i,j,k+1)-wa(i,j,k))/dz
!40        continue
!50      continue
!60    continue

      u(iecr,:,:,1:nz-1) = u(iecr,:,:,1:nz-1) - ( wa(:,:,2:nz) - wa(:,:,1:nz-1) )
      u(iecr,:,:,nz) = u(iecr,:,:,nz-1) ! for sanity

      return
      end subroutine cr_diff_z


end module crdiffusion
