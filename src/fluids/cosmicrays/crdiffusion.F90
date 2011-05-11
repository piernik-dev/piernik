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
!! \brief boundaries for wcr
!! This procedure is a shameless copy of grid::arr3d_boundaries adapted for wcr
!! \todo REMOVE ME, FIX ME, MERGE ME with grid::arr3d_boundaries or at least
!!  have a decency to make me more general
!<
   subroutine all_wcr_boundaries
   use arrays,        only: wcr   !< \todo doesn't it belong only to this module? What's the point in keeping it in arrays...
      use dataio_pub,    only: die
      use mpi,           only: MPI_STATUS_SIZE, MPI_REQUEST_NULL
      use constants,     only: CR, xdim, ydim, zdim, LO, HI, BND, BLK, BND_PER, BND_MPI
      use mpisetup,      only: comm, ierr, has_dir, psize, procxl, procyl, proczl, procxr, procyr, proczr
      use grid,          only: cg
      implicit none
      integer, parameter                          :: nreq = 3 * 4
      integer, dimension(nreq)                    :: req3d
      integer, dimension(MPI_STATUS_SIZE, nreq)   :: status3d
      integer                                     :: i, d, lh, p

      req3d(:) = MPI_REQUEST_NULL

      do d = xdim, zdim
         if (has_dir(d)) then
            do lh = LO, HI
               select case (cg%bnd(d, lh))
                  case (BND_PER)
                     do i = 1, ceiling(cg%nb/real(cg%n_b(d))) ! Repeating is important for domains that are narrower than their guardcells (e.g. cg%n_b(d) = 2)
                        select case (2*d+lh)
                           case (2*xdim+LO)
                              wcr(:,1:cg%nb, :, :) = wcr(:,cg%ieb:cg%ie, :, :)
                           case (2*ydim+LO)
                              wcr(:,:, 1:cg%nb, :) = wcr(:,:, cg%jeb:cg%je, :)
                           case (2*zdim+LO)
                              wcr(:,:, :, 1:cg%nb) = wcr(:,:, :, cg%keb:cg%ke)
                           case (2*xdim+HI)
                              wcr(:,cg%ie+1:cg%nx, :, :) = wcr(:,cg%is:cg%isb, :, :)
                           case (2*ydim+HI)
                              wcr(:,:, cg%je+1:cg%ny, :) = wcr(:,:, cg%js:cg%jsb, :)
                           case (2*zdim+HI)
                              wcr(:,:, :, cg%ke+1:cg%nz) = wcr(:,:, :, cg%ks:cg%ksb)
                        end select
                     enddo
                  case (BND_MPI)
                     if (psize(d) > 1) then
                        select case (2*d+lh)
                           case (2*xdim+LO)
                              p = procxl
                           case (2*ydim+LO)
                              p = procyl
                           case (2*zdim+LO)
                              p = proczl
                           case (2*xdim+HI)
                              p = procxr
                           case (2*ydim+HI)
                              p = procyr
                           case (2*zdim+HI)
                              p = proczr
                        end select
                        call MPI_Isend(wcr(1,1,1,1), 1, cg%mbc(CR, d, lh, BLK),  p, 2*d+(LO+HI-lh), comm, req3d(4*(d-xdim)+1+2*(lh-LO)), ierr)
                        call MPI_Irecv(wcr(1,1,1,1), 1, cg%mbc(CR, d, lh, BND),  p, 2*d+       lh,  comm, req3d(4*(d-xdim)+2+2*(lh-LO)), ierr)
                    else
                        call die("[crdiffiusion:all_wcr_boundaries] bnd_[xyz][lr] == 'mpi' && psize([xyz]dim) <= 1")
                     endif
                  case default ! Set gradient == 0 on the external boundaries
                     do i = 1, cg%nb
                        select case (2*d+lh)
                           case (2*xdim+LO)
                              wcr(:, i, :, :) = wcr(:, cg%is, :, :)
                           case (2*ydim+LO)
                              wcr(:, :, i, :) = wcr(:, :, cg%js, :)
                           case (2*zdim+LO)
                              wcr(:, :, :, i) = wcr(:, :, :, cg%ks)
                           case (2*xdim+HI)
                              wcr(:, cg%ie+i, :, :) = wcr(:, cg%ie, :, :)
                           case (2*ydim+HI)
                              wcr(:, :, cg%je+i, :) = wcr(:, :, cg%je, :)
                           case (2*zdim+HI)
                              wcr(:, :, :, cg%ke+i) = wcr(:, :, :, cg%ke)
                        end select
                     enddo
               end select

            enddo
         endif
         !>
         !! \warning outside xdim-zdim loop MPI_Waitall may change the operations order (at least for openmpi-1.4.3)
         !! and as a result may leave mpi-corners uninitiallized
         !<
         call MPI_Waitall(nreq, req3d(:), status3d(:,:), ierr)
      enddo

   end subroutine all_wcr_boundaries

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
      real, dimension(flind%crs%all)  :: decr1, decr2, decr3, fcrdif1
      real, dimension(flind%crs%all)  :: dqp, dqm

      !=======================================================================

      if (.not.has_dir(xdim)) return

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

      call all_wcr_boundaries
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
      real, dimension(flind%crs%all)  :: decr1, decr2, decr3, fcrdif2
      real, dimension(flind%crs%all)  :: dqp, dqm

!=======================================================================

      if (.not.has_dir(ydim)) return

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

      call all_wcr_boundaries
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

      call all_wcr_boundaries
      u(iarr_crs,:,:,1:cg%nz-1) = u(iarr_crs,:,:,1:cg%nz-1) - ( wcr(:,:,:,2:cg%nz) - wcr(:,:,:,1:cg%nz-1) )
      u(iarr_crs,:,:, cg%nz) = u(iarr_crs,:,:, cg%nz-1) ! for sanity

   end subroutine cr_diff_z

end module crdiffusion
