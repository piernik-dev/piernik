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
   public :: cr_diff, init_crdiffusion, cleanup_crdiffusion

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

      use constants,    only: CR, xdim, ydim, zdim, LO, HI, BND, BLK, BND_PER, BND_MPI, I_ONE
      use dataio_pub,   only: die
      use domain,       only: has_dir, cdd
      use internal_bnd, only: internal_boundaries
      use grid,         only: all_cg
      use gc_list,      only: cg_list_element
      use grid_cont,    only: grid_container
      use mpi,          only: MPI_REQUEST_NULL, MPI_COMM_NULL
      use mpisetup,     only: comm, ierr, req, status

      implicit none

      integer :: i, d, lh
      real, dimension(:,:,:,:), pointer :: pwcr
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      if (all_cg%cnt > 1) call die("[crdiffusion:all_wcr_boundaries] multiple grid pieces per procesor not implemented yet") !nontrivial MPI_Waitall should be outside do while (associated(cgl)) loop

      if (cdd%comm3d == MPI_COMM_NULL) then
         pwcr => wcr
         call internal_boundaries(CR, pa4d=pwcr)
      endif

      cgl => all_cg%first
      do while (associated(cgl))
         cg => cgl%cg

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
!! \brief Diffusive transport of ecr in crdim/ibdir-direction
!!
!! cr_diff_x --> cr_diff(xdim,ibx)
!! cr_diff_y --> cr_diff(ydim,iby)
!! cr_diff_z --> cr_diff(zdim,ibz)
!<
   subroutine cr_diff(crdim, ibdir)

      use constants,      only: xdim, ydim, zdim, ndims, LO, HI, half
      use dataio_pub,     only: die
      use domain,         only: has_dir
      use fluidindex,     only: ibx, iby, ibz, flind
      use global,         only: dt
      use grid,           only: all_cg
      use gc_list,        only: cg_list_element
      use grid_cont,      only: grid_container
      use initcosmicrays, only: iarr_crs, K_crs_paral, K_crs_perp

      implicit none

      integer(kind=4), intent(in)          :: crdim, ibdir
      integer                              :: i, j, k, il, ih, jl, jh, kl, kh, ild, jld, kld
      integer, dimension(ndims)            :: idm, ndm, hdm, ldm
      real                                 :: bb
      real, dimension(ndims)               :: bcomp
      real, dimension(flind%crs%all)       :: fcrdif
      real, dimension(ndims,flind%crs%all) :: decr
      real, dimension(flind%crs%all)       :: dqp, dqm
      type(cg_list_element), pointer       :: cgl
      type(grid_container),  pointer       :: cg
      logical, dimension(ndims)            :: present_not_crdim

!=======================================================================

      if (.not.has_dir(crdim)) return

      if (all_cg%cnt > 1) call die("[crdiffusion:cr_diff] multiple grid pieces per procesor not implemented yet") !nontrivial wcr

      idm        = 0              ;      idm(crdim) = 1
      decr(:,:)  = 0.             ;      bcomp(:)   = 0.                 ! essential where ( .not.has_dir(dim) .and. (dim /= crdim) )
      present_not_crdim = has_dir .and. ( [ xdim,ydim,zdim ] /= crdim )

      cgl => all_cg%first
      do while (associated(cgl))
         cg => cgl%cg
                                                                         ! in case of integration with boundaries:
         ldm        = cg%ijkse(:,LO) ;      ldm(crdim) = 2               ! ldm =           1 + D_
         hdm        = cg%ijkse(:,HI) ;      hdm(crdim) = cg%n_(crdim)    ! hdm = cg%n_ + idm - D_

         do k = ldm(zdim), hdm(zdim)       ; kl = k-1 ; kh = k+1 ; kld = k-idm(zdim)
            do j = ldm(ydim), hdm(ydim)    ; jl = j-1 ; jh = j+1 ; jld = j-idm(ydim)
               do i = ldm(xdim), hdm(xdim) ; il = i-1 ; ih = i+1 ; ild = i-idm(xdim)

                  decr(crdim,:) = (cg%u%arr(iarr_crs,i,j,k) - cg%u%arr(iarr_crs,ild,jld,kld)) * cg%idl(crdim)
                  fcrdif = K_crs_perp * decr(crdim,:)

                  bcomp(ibdir) =  cg%b%arr(ibdir,i,j,k)

                  if (present_not_crdim(xdim)) then
                     dqm = half*((cg%u%arr(iarr_crs,i ,jld,kld) + cg%u%arr(iarr_crs,i ,j,k)) - (cg%u%arr(iarr_crs,il,jld,kld) + cg%u%arr(iarr_crs,il,j,k))) * cg%idx
                     dqp = half*((cg%u%arr(iarr_crs,ih,jld,kld) + cg%u%arr(iarr_crs,ih,j,k)) - (cg%u%arr(iarr_crs,i ,jld,kld) + cg%u%arr(iarr_crs,i ,j,k))) * cg%idx
                     decr(xdim,:) = (dqp+dqm)* (1.0 + sign(1.0, dqm*dqp))*0.25
                     bcomp(ibx)   = sum(cg%b%arr(ibx,i:ih, jld:j, kld:k))*0.25
                  endif

                  if (present_not_crdim(ydim)) then
                     dqm = half*((cg%u%arr(iarr_crs,ild,j ,kld) + cg%u%arr(iarr_crs,i,j ,k)) - (cg%u%arr(iarr_crs,ild,jl,kld) + cg%u%arr(iarr_crs,i,jl,k))) * cg%idy
                     dqp = half*((cg%u%arr(iarr_crs,ild,jh,kld) + cg%u%arr(iarr_crs,i,jh,k)) - (cg%u%arr(iarr_crs,ild,j ,kld) + cg%u%arr(iarr_crs,i,j ,k))) * cg%idy
                     decr(ydim,:) = (dqp+dqm)* (1.0 + sign(1.0, dqm*dqp))*0.25
                     bcomp(iby)   = sum(cg%b%arr(iby,ild:i, j:jh, kld:k))*0.25
                  endif

                  if (present_not_crdim(zdim)) then
                     dqm = half*((cg%u%arr(iarr_crs,ild,jld,k ) + cg%u%arr(iarr_crs,i,j,k )) - (cg%u%arr(iarr_crs,ild,jld,kl) + cg%u%arr(iarr_crs,i,j,kl))) * cg%idz
                     dqp = half*((cg%u%arr(iarr_crs,ild,jld,kh) + cg%u%arr(iarr_crs,i,j,kh)) - (cg%u%arr(iarr_crs,ild,jld,k ) + cg%u%arr(iarr_crs,i,j,k ))) * cg%idz
                     decr(zdim,:) = (dqp+dqm)* (1.0 + sign(1.0, dqm*dqp))*0.25
                     bcomp(ibz)   = sum(cg%b%arr(ibz,ild:i, jld:j, k:kh))*0.25
                  endif

                  bb = bcomp(ibx)**2 + bcomp(iby)**2 + bcomp(ibz)**2
                  if (bb > epsilon(0.d0)) fcrdif = fcrdif + K_crs_paral * bcomp(ibdir) * (bcomp(ibx)*decr(xdim,:) + bcomp(iby)*decr(ydim,:) + bcomp(ibz)*decr(zdim,:)) / bb

                  wcr(:,i,j,k) = - fcrdif * dt * cg%idl(crdim)

               enddo
            enddo
         enddo

         call all_wcr_boundaries

         ndm = cg%n_ - idm
         hdm = 1 + idm*ndm
         ldm = hdm - idm
         cg%u%arr(iarr_crs,:ndm(xdim), :ndm(ydim), :ndm(zdim)) =                         cg%u%arr(iarr_crs,:ndm(xdim),:ndm(ydim),:ndm(zdim)) &
           &    - ( wcr(:,1+idm(xdim):cg%n_(xdim),1+idm(ydim):cg%n_(ydim),1+idm(zdim):cg%n_(zdim)) - wcr(:,:ndm(xdim),:ndm(ydim),:ndm(zdim)) )
         cg%u%arr(iarr_crs,hdm(xdim):cg%n_(xdim),hdm(ydim):cg%n_(ydim),hdm(zdim):cg%n_(zdim)) = cg%u%arr(iarr_crs,ldm(xdim):ndm(xdim),ldm(ydim):ndm(ydim),ldm(zdim):ndm(zdim)) ! for sanity

         cgl => cgl%nxt
      enddo

   end subroutine cr_diff

end module crdiffusion
