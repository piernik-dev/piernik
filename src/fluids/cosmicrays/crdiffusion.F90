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
   public :: cr_diff, init_crdiffusion

   logical :: has_cr

contains

   subroutine init_crdiffusion(crsall)

      use constants,  only: wcr_n, AT_IGNORE
      use cr_data,    only: divv_n
      use dataio_pub, only: warn
      use grid,       only: all_cg

      implicit none

      integer(kind=4), intent(in) :: crsall

      has_cr = (crsall > 0)

      if (has_cr) then
         call all_cg%reg_var(wcr_n, AT_IGNORE, crsall)
      else
         call warn("[crdiffusion:init_crdiffusion] No CR species to diffuse")
      endif
      call all_cg%reg_var(divv_n, AT_IGNORE)

   end subroutine init_crdiffusion

!>
!! \brief boundaries for wcr
!! This procedure is a shameless copy of grid::arr3d_boundaries adapted for wcr
!! \todo REMOVE ME, FIX ME, MERGE ME with grid::arr3d_boundaries or at least
!!  have a decency to make me more general
!<
   subroutine all_wcr_boundaries

      use constants,    only: CR, xdim, ydim, zdim, LO, HI, BND, BLK, BND_PER, BND_MPI, I_ONE, wcr_n
      use dataio_pub,   only: die
      use domain,       only: dom
      use internal_bnd, only: internal_boundaries_4d
      use grid,         only: leaves, all_cg
      use gc_list,      only: cg_list_element
      use grid_cont,    only: grid_container
      use mpi,          only: MPI_REQUEST_NULL, MPI_COMM_NULL
      use mpisetup,     only: comm, ierr, req, status
      use types,        only: cdd

      implicit none

      integer :: i, d, lh
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg
      real, dimension(:,:,:,:), pointer :: wcr

      if (.not. has_cr) return

      if (cdd%comm3d == MPI_COMM_NULL) call internal_boundaries_4d(all_cg, all_cg%ind_4d(wcr_n)) ! should be more selective (modified leaves?)

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         wcr => cg%w(all_cg%ind_4d(wcr_n))%arr
         if (.not. associated(wcr)) call die("[crdiffusion:all_wcr_boundaries] cannot get wcr")

         req(:) = MPI_REQUEST_NULL

         do d = xdim, zdim
            if (dom%has_dir(d)) then
               do lh = LO, HI
                  select case (cg%bnd(d, lh))
                     case (BND_PER)
                        if (cdd%comm3d /= MPI_COMM_NULL) then
                           select case (2*d+lh)
                              case (2*xdim+LO)
                                 wcr(:,1:dom%nb, :, :) = wcr(:,cg%ieb:cg%ie, :, :)
                              case (2*ydim+LO)
                                 wcr(:,:, 1:dom%nb, :) = wcr(:,:, cg%jeb:cg%je, :)
                              case (2*zdim+LO)
                                 wcr(:,:, :, 1:dom%nb) = wcr(:,:, :, cg%keb:cg%ke)
                              case (2*xdim+HI)
                                 wcr(:,cg%ie+1:cg%n_(xdim), :, :) = wcr(:,cg%is:cg%isb, :, :)
                              case (2*ydim+HI)
                                 wcr(:,:, cg%je+1:cg%n_(ydim), :) = wcr(:,:, cg%js:cg%jsb, :)
                              case (2*zdim+HI)
                                 wcr(:,:, :, cg%ke+1:cg%n_(zdim)) = wcr(:,:, :, cg%ks:cg%ksb)
                           end select
                        endif
                     case (BND_MPI)
                        if (cdd%comm3d /= MPI_COMM_NULL) then
                           if (cdd%psize(d) > 1) then
                              call MPI_Isend(wcr(1,1,1,1), I_ONE, cg%mbc(CR, d, lh, BLK, dom%nb), cdd%procn(d, lh), int(2*d+(LO+HI-lh), kind=4), cdd%comm3d, req(4*(d-xdim)+1+2*(lh-LO)), ierr)
                              call MPI_Irecv(wcr(1,1,1,1), I_ONE, cg%mbc(CR, d, lh, BND, dom%nb), cdd%procn(d, lh), int(2*d+       lh,  kind=4), cdd%comm3d, req(4*(d-xdim)+2+2*(lh-LO)), ierr)
                           else
                              call die("[crdiffiusion:all_wcr_boundaries] bnd_[xyz][lr] == 'mpi' && psize(:) <= 1")
                           endif
                        endif
                     case default ! Set gradient == 0 on the external boundaries
                        do i = 1, dom%nb
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
!! cr_diff_x --> cr_diff(xdim)
!! cr_diff_y --> cr_diff(ydim)
!! cr_diff_z --> cr_diff(zdim)
!<
   subroutine cr_diff(crdim)

      use constants,      only: xdim, ydim, zdim, ndims, LO, HI, half, wcr_n
      use dataio_pub,     only: die
      use domain,         only: dom
      use fluidindex,     only: flind
      use global,         only: dt
      use grid,           only: leaves, all_cg
      use gc_list,        only: cg_list_element
      use grid_cont,      only: grid_container
      use initcosmicrays, only: iarr_crs, K_crs_paral, K_crs_perp

      implicit none

      integer(kind=4), intent(in)          :: crdim
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
      real, dimension(:,:,:,:), pointer    :: wcr
      integer                              :: wcri

      if (.not. has_cr) return
      if (.not.dom%has_dir(crdim)) return

      idm        = 0              ;      idm(crdim) = 1
      decr(:,:)  = 0.             ;      bcomp(:)   = 0.                 ! essential where ( .not.dom%has_dir(dim) .and. (dim /= crdim) )
      present_not_crdim = dom%has_dir .and. ( [ xdim,ydim,zdim ] /= crdim )
      wcri = all_cg%ind_4d(wcr_n)

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         wcr => cg%w(wcri)%arr
         if (.not. associated(wcr)) call die("[crdiffusion:cr_diff] cannot get wcr")
                                                                         ! in case of integration with boundaries:
         ldm        = cg%ijkse(:,LO) ;      ldm(crdim) = 2               ! ldm =           1 + D_
         hdm        = cg%ijkse(:,HI) ;      hdm(crdim) = cg%n_(crdim)    ! hdm = cg%n_ + idm - D_

         do k = ldm(zdim), hdm(zdim)       ; kl = k-1 ; kh = k+1 ; kld = k-idm(zdim)
            do j = ldm(ydim), hdm(ydim)    ; jl = j-1 ; jh = j+1 ; jld = j-idm(ydim)
               do i = ldm(xdim), hdm(xdim) ; il = i-1 ; ih = i+1 ; ild = i-idm(xdim)

                  decr(crdim,:) = (cg%u(iarr_crs,i,j,k) - cg%u(iarr_crs,ild,jld,kld)) * cg%idl(crdim)
                  fcrdif = K_crs_perp * decr(crdim,:)

                  bcomp(crdim) =  cg%b(crdim,i,j,k)

                  if (present_not_crdim(xdim)) then
                     dqm = half*((cg%u(iarr_crs,i ,jld,kld) + cg%u(iarr_crs,i ,j,k)) - (cg%u(iarr_crs,il,jld,kld) + cg%u(iarr_crs,il,j,k))) * cg%idx
                     dqp = half*((cg%u(iarr_crs,ih,jld,kld) + cg%u(iarr_crs,ih,j,k)) - (cg%u(iarr_crs,i ,jld,kld) + cg%u(iarr_crs,i ,j,k))) * cg%idx
                     decr(xdim,:) = (dqp+dqm)* (1.0 + sign(1.0, dqm*dqp))*0.25
                     bcomp(xdim)   = sum(cg%b(xdim,i:ih, jld:j, kld:k))*0.25
                  endif

                  if (present_not_crdim(ydim)) then
                     dqm = half*((cg%u(iarr_crs,ild,j ,kld) + cg%u(iarr_crs,i,j ,k)) - (cg%u(iarr_crs,ild,jl,kld) + cg%u(iarr_crs,i,jl,k))) * cg%idy
                     dqp = half*((cg%u(iarr_crs,ild,jh,kld) + cg%u(iarr_crs,i,jh,k)) - (cg%u(iarr_crs,ild,j ,kld) + cg%u(iarr_crs,i,j ,k))) * cg%idy
                     decr(ydim,:) = (dqp+dqm)* (1.0 + sign(1.0, dqm*dqp))*0.25
                     bcomp(ydim)   = sum(cg%b(ydim,ild:i, j:jh, kld:k))*0.25
                  endif

                  if (present_not_crdim(zdim)) then
                     dqm = half*((cg%u(iarr_crs,ild,jld,k ) + cg%u(iarr_crs,i,j,k )) - (cg%u(iarr_crs,ild,jld,kl) + cg%u(iarr_crs,i,j,kl))) * cg%idz
                     dqp = half*((cg%u(iarr_crs,ild,jld,kh) + cg%u(iarr_crs,i,j,kh)) - (cg%u(iarr_crs,ild,jld,k ) + cg%u(iarr_crs,i,j,k ))) * cg%idz
                     decr(zdim,:) = (dqp+dqm)* (1.0 + sign(1.0, dqm*dqp))*0.25
                     bcomp(zdim)   = sum(cg%b(zdim,ild:i, jld:j, k:kh))*0.25
                  endif

                  bb = sum(bcomp**2)
                  if (bb > epsilon(0.d0)) fcrdif = fcrdif + K_crs_paral * bcomp(crdim) * (bcomp(xdim)*decr(xdim,:) + bcomp(ydim)*decr(ydim,:) + bcomp(zdim)*decr(zdim,:)) / bb

                  wcr(:,i,j,k) = - fcrdif * dt * cg%idl(crdim)

               enddo
            enddo
         enddo
         cgl => cgl%nxt
      enddo

      call all_wcr_boundaries

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         wcr => cg%w(wcri)%arr

         ndm = cg%n_ - idm
         hdm = 1 + idm*ndm
         ldm = hdm - idm
         cg%u(iarr_crs,:ndm(xdim), :ndm(ydim), :ndm(zdim)) =                         cg%u(iarr_crs,:ndm(xdim),:ndm(ydim),:ndm(zdim)) &
           &    - ( wcr(:,1+idm(xdim):cg%n_(xdim),1+idm(ydim):cg%n_(ydim),1+idm(zdim):cg%n_(zdim)) - wcr(:,:ndm(xdim),:ndm(ydim),:ndm(zdim)) )
         cg%u(iarr_crs,hdm(xdim):cg%n_(xdim),hdm(ydim):cg%n_(ydim),hdm(zdim):cg%n_(zdim)) = cg%u(iarr_crs,ldm(xdim):ndm(xdim),ldm(ydim):ndm(ydim),ldm(zdim):ndm(zdim)) ! for sanity

         cgl => cgl%nxt
      enddo

   end subroutine cr_diff

end module crdiffusion
