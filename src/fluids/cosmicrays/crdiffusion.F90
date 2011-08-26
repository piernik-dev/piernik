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
   public :: cr_diff_x, cr_diff_y, cr_diff_z, init_crdiffusion, cleanup_crdiffusion

   real, allocatable, dimension(:,:,:,:), target :: wcr      !< Temporary array used in crdiffusion module

contains

   subroutine init_crdiffusion(crsall)

      use dataio_pub,  only: die
      use diagnostics, only: ma4d, my_allocate
      use grid,        only: all_cg
      use gc_list,     only: cg_list_element
      use grid_cont,   only: grid_container

      implicit none

      integer(kind=4), intent(in) :: crsall

      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      if (all_cg%cnt > 1) call die("[crdiffusion:init_crdiffusion] multiple grid pieces per procesor not implemented yet") !nontrivial crsall
      !> \todo provide hooks for rank-4 user/physics arrays in grid container

      cgl => all_cg%first
      do while (associated(cgl))
         cg => cgl%cg
         ma4d = [crsall, cg%n_(:) ]
         call my_allocate(wcr, ma4d, "wcr")
         cgl => cgl%nxt
      enddo

   end subroutine init_crdiffusion

   subroutine cleanup_crdiffusion
      use diagnostics, only: my_deallocate
      implicit none
      call my_deallocate(wcr)
   end subroutine cleanup_crdiffusion

!>
!! \brief boundaries for wcr
!! This procedure is a shameless copy of grid::arr3d_boundaries adapted for wcr
!! \todo REMOVE ME, FIX ME, MERGE ME with grid::arr3d_boundaries or at least
!!  have a decency to make me more general
!<
   subroutine all_wcr_boundaries

      use constants,  only: CR, xdim, ydim, zdim, LO, HI, BND, BLK, BND_PER, BND_MPI, I_ONE
      use dataio_pub, only: die
      use domain,     only: has_dir, cdd
      use grid,       only: all_cg
      use gc_list,    only: cg_list_element
      use grid_cont,  only: grid_container
      use mpi,        only: MPI_REQUEST_NULL, MPI_COMM_NULL
      use mpisetup,   only: comm, ierr, req, status

      implicit none

      integer :: i, d, lh
      real, dimension(:,:,:,:), pointer :: pwcr
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      if (all_cg%cnt > 1) call die("[crdiffusion:all_wcr_boundaries] multiple grid pieces per procesor not implemented yet") !nontrivial MPI_Waitall should be outside do while (associated(cgl)) loop
      cgl => all_cg%first
      do while (associated(cgl))
         cg => cgl%cg

         if (cdd%comm3d == MPI_COMM_NULL) then
            pwcr => wcr
            call cg%internal_boundaries(CR, pa4d=pwcr)
         endif

         req(:) = MPI_REQUEST_NULL

         do d = xdim, zdim
            if (has_dir(d)) then
               do lh = LO, HI
                  select case (cg%bnd(d, lh))
                     case (BND_PER)
                        if (cdd%comm3d /= MPI_COMM_NULL) then
                           do i = 1, ceiling(cg%nb/real(cg%n_b(d))) ! Repeating is important for domains that are narrower than their guardcells (e.g. cg%n_b(d) = 2)
                              select case (2*d+lh)
                                 case (2*xdim+LO)
                                    wcr(:,1:cg%nb, :, :) = wcr(:,cg%ieb:cg%ie, :, :)
                                 case (2*ydim+LO)
                                    wcr(:,:, 1:cg%nb, :) = wcr(:,:, cg%jeb:cg%je, :)
                                 case (2*zdim+LO)
                                    wcr(:,:, :, 1:cg%nb) = wcr(:,:, :, cg%keb:cg%ke)
                                 case (2*xdim+HI)
                                    wcr(:,cg%ie+1:cg%n_(xdim), :, :) = wcr(:,cg%is:cg%isb, :, :)
                                 case (2*ydim+HI)
                                    wcr(:,:, cg%je+1:cg%n_(ydim), :) = wcr(:,:, cg%js:cg%jsb, :)
                                 case (2*zdim+HI)
                                    wcr(:,:, :, cg%ke+1:cg%n_(zdim)) = wcr(:,:, :, cg%ks:cg%ksb)
                              end select
                           enddo
                        endif
                     case (BND_MPI)
                        if (cdd%comm3d /= MPI_COMM_NULL) then
                           if (cdd%psize(d) > 1) then
                              call MPI_Isend(wcr(1,1,1,1), I_ONE, cg%mbc(CR, d, lh, BLK), cdd%procn(d, lh), int(2*d+(LO+HI-lh), kind=4), cdd%comm3d, req(4*(d-xdim)+1+2*(lh-LO)), ierr)
                              call MPI_Irecv(wcr(1,1,1,1), I_ONE, cg%mbc(CR, d, lh, BND), cdd%procn(d, lh), int(2*d+       lh,  kind=4), cdd%comm3d, req(4*(d-xdim)+2+2*(lh-LO)), ierr)
                           else
                              call die("[crdiffiusion:all_wcr_boundaries] bnd_[xyz][lr] == 'mpi' && psize(:) <= 1")
                           endif
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
            if (cdd%comm3d /= MPI_COMM_NULL) call MPI_Waitall(size(req(:)), req(:), status(:,:), ierr)
         enddo

         cgl => cgl%nxt
      enddo

   end subroutine all_wcr_boundaries

!>
!! \brief Diffusive transport of ecr in 1-direction
!<
   subroutine cr_diff_x

      use constants,      only: xdim, ydim, zdim, half
      use dataio_pub,     only: die
      use domain,         only: has_dir
      use fluidindex,     only: ibx, iby, ibz, flind
      use global,         only: dt
      use grid,           only: all_cg
      use gc_list,        only: cg_list_element
      use grid_cont,      only: grid_container
      use initcosmicrays, only: iarr_crs, K_crs_paral, K_crs_perp

      implicit none

      integer :: i, j, k
      real    :: b1b, b2b, b3b, bb
      real, dimension(flind%crs%all)  :: decr1, decr2, decr3, fcrdif1
      real, dimension(flind%crs%all)  :: dqp, dqm
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      !=======================================================================

      if (.not.has_dir(xdim)) return

      if (all_cg%cnt > 1) call die("[crdiffusion:cr_diff_x] multiple grid pieces per procesor not implemented yet") !nontrivial wcr

      cgl => all_cg%first
      do while (associated(cgl))
         cg => cgl%cg

         do k=cg%ks, cg%ke
            do j=cg%js, cg%je
               do i=2, cg%n_(xdim)     ! if we are here this implies nxb /= 1

                  decr1 =  (cg%u%arr(iarr_crs,i,  j,  k) - cg%u%arr(iarr_crs,i-1,j,  k)) * cg%idx
                  fcrdif1 = K_crs_perp * decr1

                  b1b =  cg%b%arr(ibx,i,  j,  k)

                  if (has_dir(ydim)) then
                     dqm = half*((cg%u%arr(iarr_crs,i-1,j  ,k ) + cg%u%arr(iarr_crs,i ,j  ,k )) - (cg%u%arr(iarr_crs,i-1,j-1,k ) + cg%u%arr(iarr_crs,i ,j-1,k ))) * cg%idy
                     dqp = half*((cg%u%arr(iarr_crs,i-1,j+1,k ) + cg%u%arr(iarr_crs,i ,j+1,k )) - (cg%u%arr(iarr_crs,i-1,j  ,k ) + cg%u%arr(iarr_crs,i ,j  ,k ))) * cg%idy
                     decr2 = (dqp+dqm)* (1.0 + sign(1.0, dqm*dqp))*0.25
                     b2b = sum(cg%b%arr(iby,i-1:i, j:j+1, k    ))*0.25
                  else
                     decr2 = 0.
                     b2b = 0.
                  endif

                  if (has_dir(zdim)) then
                     dqm = half*((cg%u%arr(iarr_crs,i-1,j ,k  ) + cg%u%arr(iarr_crs,i ,j ,k  )) - (cg%u%arr(iarr_crs,i-1,j ,k-1) + cg%u%arr(iarr_crs,i ,j ,k-1))) * cg%idz
                     dqp = half*((cg%u%arr(iarr_crs,i-1,j ,k+1) + cg%u%arr(iarr_crs,i ,j ,k+1)) - (cg%u%arr(iarr_crs,i-1,j ,k  ) + cg%u%arr(iarr_crs,i ,j ,k  ))) * cg%idz
                     decr3 = (dqp+dqm)* (1.0 + sign(1.0, dqm*dqp))*0.25
                     b3b = sum(cg%b%arr(ibz,i-1:i, j,     k:k+1))*0.25
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
         cg%u%arr(iarr_crs,1:cg%n_(xdim)-1,:,:) = cg%u%arr(iarr_crs,1:cg%n_(xdim)-1,:,:) - ( wcr(:,2:cg%n_(xdim),:,:) - wcr(:,1:cg%n_(xdim)-1,:,:) )
         cg%u%arr(iarr_crs, cg%n_(xdim),:,:) = cg%u%arr(iarr_crs, cg%n_(xdim)-1,:,:) ! for sanity

         cgl => cgl%nxt
      enddo

   end subroutine cr_diff_x

!>
!! \brief Diffusive transport of ecr in 2-direction
!<
   subroutine cr_diff_y

      use constants,      only: xdim, ydim, zdim, half
      use dataio_pub,     only: die
      use domain,         only: has_dir
      use fluidindex,     only: ibx, iby, ibz, flind
      use global,         only: dt
      use grid,           only: all_cg
      use gc_list,        only: cg_list_element
      use grid_cont,      only: grid_container
      use initcosmicrays, only: iarr_crs, K_crs_paral, K_crs_perp

      implicit none

      integer :: i, j, k
      real    :: b1b, b2b, b3b, bb
      real, dimension(flind%crs%all)  :: decr1, decr2, decr3, fcrdif2
      real, dimension(flind%crs%all)  :: dqp, dqm
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

!=======================================================================

      if (.not.has_dir(ydim)) return

      if (all_cg%cnt > 1) call die("[crdiffusion:cr_diff_y] multiple grid pieces per procesor not implemented yet") !nontrivial wcr

      cgl => all_cg%first
      do while (associated(cgl))
         cg => cgl%cg

         do k=cg%ks, cg%ke
            do j=2, cg%n_(ydim) ! if we are here nyb /= 1
               do i=cg%is, cg%ie

                  decr2 = (cg%u%arr(iarr_crs,i,j,k) - cg%u%arr(iarr_crs,i,j-1,k)) * cg%idy
                  fcrdif2 = K_crs_perp * decr2

                  b2b =  cg%b%arr(iby,i,j,k)

                  if (has_dir(xdim)) then
                     dqm = half*((cg%u%arr(iarr_crs,i  ,j-1,k ) + cg%u%arr(iarr_crs,i  ,j  ,k )) - (cg%u%arr(iarr_crs,i-1,j-1,k ) + cg%u%arr(iarr_crs,i-1,j  ,k ))) * cg%idx
                     dqp = half*((cg%u%arr(iarr_crs,i+1,j-1,k ) + cg%u%arr(iarr_crs,i+1,j  ,k )) - (cg%u%arr(iarr_crs,i  ,j-1,k ) + cg%u%arr(iarr_crs,i  ,j  ,k ))) * cg%idx
                     decr1 = (dqp+dqm)* (1.0 + sign(1.0, dqm*dqp))*0.25
                     b1b = sum(cg%b%arr(ibx,i:i+1, j-1:j, k))*0.25
                  else
                     decr1 = 0.
                     b1b = 0.
                  endif

                  if (has_dir(zdim)) then
                     dqm = half*((cg%u%arr(iarr_crs,i ,j-1,k  ) + cg%u%arr(iarr_crs,i ,j  ,k  )) - (cg%u%arr(iarr_crs,i ,j-1,k-1) + cg%u%arr(iarr_crs,i ,j  ,k-1))) * cg%idz
                     dqp = half*((cg%u%arr(iarr_crs,i ,j-1,k+1) + cg%u%arr(iarr_crs,i ,j  ,k+1)) - (cg%u%arr(iarr_crs,i ,j-1,k  ) + cg%u%arr(iarr_crs,i ,j  ,k  ))) * cg%idz
                     decr3 = (dqp+dqm)* (1.0 + sign(1.0, dqm*dqp))*0.25
                     b3b = sum(cg%b%arr(ibz,i, j-1:j, k:k+1))*0.25
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
         cg%u%arr(iarr_crs,:,1:cg%n_(ydim)-1,:) = cg%u%arr(iarr_crs,:,1:cg%n_(ydim)-1,:) - ( wcr(:,:,2:cg%n_(ydim),:) - wcr(:,:,1:cg%n_(ydim)-1,:) )
         cg%u%arr(iarr_crs,:, cg%n_(ydim),:) = cg%u%arr(iarr_crs,:, cg%n_(ydim)-1,:) ! for sanity

         cgl => cgl%nxt
      enddo

   end subroutine cr_diff_y

!>
!! \brief Diffusive transport of ecr in 3-direction
!<
   subroutine cr_diff_z

      use constants,      only: xdim, ydim, zdim, half
      use dataio_pub,     only: die
      use domain,         only: has_dir
      use fluidindex,     only: ibx, iby, ibz, flind
      use global,         only: dt
      use grid,           only: all_cg
      use gc_list,        only: cg_list_element
      use grid_cont,      only: grid_container
      use initcosmicrays, only: iarr_crs, K_crs_paral, K_crs_perp

      implicit none

      integer :: i, j, k
      real    :: b1b, b2b, b3b, bb
      real, dimension(flind%crs%all) :: decr1, decr2, decr3, fcrdif3
      real, dimension(flind%crs%all) :: dqp, dqm
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

!=======================================================================

      if (.not.has_dir(zdim)) return

      if (all_cg%cnt > 1) call die("[crdiffusion:cr_diff_z] multiple grid pieces per procesor not implemented yet") !nontrivial wcr

      cgl => all_cg%first
      do while (associated(cgl))
         cg => cgl%cg

         do k=2, cg%n_(zdim)      ! nzb /= 1
            do j=cg%js, cg%je
               do i=cg%is, cg%ie

                  decr3 =  (cg%u%arr(iarr_crs,i,j,k    ) - cg%u%arr(iarr_crs,i,j,  k-1)) * cg%idz
                  fcrdif3 = K_crs_perp * decr3

                  b3b =  cg%b%arr(ibz,i,j,k)

                  if (has_dir(xdim)) then
                     dqm = half*((cg%u%arr(iarr_crs,i  ,j ,k-1) + cg%u%arr(iarr_crs,i  ,j  ,k )) - (cg%u%arr(iarr_crs,i-1,j ,k-1) + cg%u%arr(iarr_crs,i-1,j  ,k ))) * cg%idx
                     dqp = half*((cg%u%arr(iarr_crs,i+1,j ,k-1) + cg%u%arr(iarr_crs,i+1,j  ,k )) - (cg%u%arr(iarr_crs,i  ,j ,k-1) + cg%u%arr(iarr_crs,i  ,j  ,k ))) * cg%idx
                     decr1 = (dqp+dqm)* (1.0 + sign(1.0, dqm*dqp))*0.25
                     b1b = sum(cg%b%arr(ibx,i:i+1, j,     k-1:k))*0.25
                  else
                     decr1 = 0.
                     b1b = 0.
                  endif

                  if (has_dir(ydim)) then
                     dqm = half*((cg%u%arr(iarr_crs,i ,j  ,k-1) + cg%u%arr(iarr_crs,i  ,j  ,k )) - (cg%u%arr(iarr_crs,i ,j-1,k-1) + cg%u%arr(iarr_crs,i  ,j-1,k ))) * cg%idy
                     dqp = half*((cg%u%arr(iarr_crs,i ,j+1,k-1) + cg%u%arr(iarr_crs,i  ,j+1,k )) - (cg%u%arr(iarr_crs,i ,j  ,k-1) + cg%u%arr(iarr_crs,i  ,j  ,k ))) * cg%idy
                     decr2 = (dqp+dqm)* (1.0 + sign(1.0, dqm*dqp))*0.25
                     b2b = sum(cg%b%arr(iby,i,     j:j+1, k-1:k))*0.25
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
         cg%u%arr(iarr_crs,:,:,1:cg%n_(zdim)-1) = cg%u%arr(iarr_crs,:,:,1:cg%n_(zdim)-1) - ( wcr(:,:,:,2:cg%n_(zdim)) - wcr(:,:,:,1:cg%n_(zdim)-1) )
         cg%u%arr(iarr_crs,:,:, cg%n_(zdim)) = cg%u%arr(iarr_crs,:,:, cg%n_(zdim)-1) ! for sanity

         cgl => cgl%nxt
      enddo

   end subroutine cr_diff_z

end module crdiffusion
