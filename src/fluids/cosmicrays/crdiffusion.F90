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
!! \brief Numerical scheme for the diffusive transport of Cosmic Rays
!<

module crdiffusion

! pulled by COSM_RAYS

   implicit none

   private
   public :: cr_diff, init_crdiffusion

   logical :: has_cr

contains

   subroutine init_crdiffusion

      use cg_list_global,   only: all_cg
      use constants,        only: wcr_n
      use crhelpers,        only: divv_i, divv_n
      use dataio_pub,       only: warn
      use fluidindex,       only: flind
      use named_array_list, only: qna

      implicit none

      has_cr = (flind%crs%all > 0)

      if (has_cr) then
         call all_cg%reg_var(wcr_n, dim4 = flind%crs%all) !, ord_prolong = 2)  ! Smooth prolongation may help with interpolation from coarse grid to fine boundary
      else
         call warn("[crdiffusion:init_crdiffusion] No CR species to diffuse")
      endif
      call all_cg%reg_var(divv_n)
      divv_i = qna%ind(divv_n)

   end subroutine init_crdiffusion

!>
!! \brief boundaries for wcr
!! This procedure is a shameless copy of cg_list_bnd%arr3d_boundaries adapted for wcr
!! \todo REMOVE ME, FIX ME, MERGE ME with cg_list_bnd%arr3d_boundaries or at least have a decency to make me more general
!<
   subroutine all_wcr_boundaries

      use cg_cost_data,       only: I_DIFFUSE
      use cg_leaves,          only: leaves
      use cg_level_connected, only: cg_level_connected_t
      use cg_level_finest,    only: finest
      use cg_list,            only: cg_list_element
      use constants,          only: ndims, xdim, ydim, zdim, LO, HI, BND_PER, BND_MPI, BND_FC, BND_MPI_FC, I_TWO, I_THREE, wcr_n, PPP_CR
      use dataio_pub,         only: die
      use domain,             only: dom
      use grid_cont,          only: grid_container
      use initcosmicrays,     only: diff_max_lev, diff_prolong
      use named_array_list,   only: wna
      use ppp,                only: ppp_main

      implicit none

      integer(kind=4)                          :: i, d, lh, wcri
      integer(kind=4), dimension(ndims, LO:HI) :: l, r
      real, dimension(:,:,:,:),   pointer      :: wcr
      type(cg_list_element),      pointer      :: cgl
      type(grid_container),       pointer      :: cg
      type(cg_level_connected_t), pointer      :: diffl, curl
      character(len=*), parameter :: awb_label = "all_wcr_boundaries"

      if (.not. has_cr) return

      wcri = wna%ind(wcr_n)

      call ppp_main%start(awb_label, PPP_CR)
      if (diff_max_lev >= finest%level%l%id) then  ! do the old way
         call finest%level%restrict_to_base_w_1var(wcri)
         call leaves%leaf_arr4d_boundaries(wcri)
      else
         diffl => finest%level
         do while (associated(diffl))
            if (diffl%l%id == diff_max_lev) exit
            diffl => diffl%coarser
         enddo
         if (.not. associated(diffl)) call die("[crdiffusion:all_wcr_boundaries] .not. associated(diffl)")
         if (diffl%l%id /= diff_max_lev) call die("[crdiffusion:all_wcr_boundaries] diffl%l%id /= diff_max_lev")
         call diffl%restrict_to_base_w_1var(wcri)

         ! slightly modified copy of cg_leaves::leaf_arr4d_boundaries
         curl => diffl
         do while (associated(curl))
            call curl%level_4d_boundaries(wcri)
            curl => curl%coarser
         enddo

         curl => diffl
         do while (associated(curl))
            call curl%prolong_bnd_from_coarser(wcri, arr4d=.true., nocorners=.false.)
            ! corners are required on all levels except for finest anyway if prolongation order is greater than injection
            curl => curl%coarser
         enddo

      endif

      ! do the external boundaries
      cgl => leaves%up_to_level(min(diff_max_lev, finest%level%l%id))%p
      do while (associated(cgl))
         cg => cgl%cg
         call cg%costs%start

         wcr => cg%w(wcri)%arr
         if (.not. associated(wcr)) call die("[crdiffusion:all_wcr_boundaries] cannot get wcr")

         do d = xdim, zdim
            if (dom%has_dir(d)) then
               l = cg%lhn ; r = l
               do lh = LO, HI
                  select case (cg%bnd(d, lh))
                     case (BND_PER)
                     case (BND_MPI) !, BND_MPI_FC)
                     case (BND_FC, BND_MPI_FC)
                     case default ! Set gradient == 0 on the external boundaries
                        r(d,:) = cg%ijkse(d,lh)
                        do i = 1, dom%nb
                           l(d,:) = r(d,:) + (I_TWO*lh-I_THREE)*i
                           wcr(:,l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) = wcr(:,r(xdim,LO):r(xdim,HI),r(ydim,LO):r(ydim,HI),r(zdim,LO):r(zdim,HI))
                        enddo
                  end select

               enddo
            endif
         enddo

         call cg%costs%stop(I_DIFFUSE)
         cgl => cgl%nxt
      enddo

      call ppp_main%stop(awb_label, PPP_CR)

   end subroutine all_wcr_boundaries

!>
!! \brief Diffusive transport of ecr in crdim/ibdir-direction
!!
!! cr_diff_x --> cr_diff(xdim)
!! cr_diff_y --> cr_diff(ydim)
!! cr_diff_z --> cr_diff(zdim)
!<
   subroutine cr_diff(crdim)

      use all_boundaries,   only: all_bnd
      use cg_cost_data,     only: I_DIFFUSE
      use cg_leaves,        only: leaves
      use cg_level_finest,  only: finest
      use cg_list,          only: cg_list_element
      use constants,        only: xdim, ydim, zdim, ndims, LO, HI, oneeig, eight, wcr_n, GEO_XYZ, PPP_CR
      use dataio_pub,       only: die
      use domain,           only: dom
      use fluidindex,       only: flind
      use global,           only: dt
      use grid_cont,        only: grid_container
      use initcosmicrays,   only: iarr_crs, K_crs_paral, K_crs_perp, diff_max_lev
      use named_array,      only: p4
      use named_array_list, only: wna
      use ppp,              only: ppp_main
#ifdef MAGNETIC
      use constants,        only: four
#endif /* MAGNETIC */

      implicit none

      integer(kind=4), intent(in)          :: crdim
      integer                              :: i, j, k, il, ih, jl, jh, kl, kh, ild, jld, kld
      integer, dimension(ndims)            :: idm, ndm, hdm, ldm
      real                                 :: bb, f1, f2
      real, dimension(ndims)               :: bcomp
      real, dimension(flind%crs%all)       :: fcrdif
      real, dimension(ndims,flind%crs%all) :: decr
      real, dimension(flind%crs%all)       :: dqp, dqm
      type(cg_list_element), pointer       :: cgl
      type(grid_container),  pointer       :: cg
      logical, dimension(ndims)            :: present_not_crdim
      real, dimension(:,:,:,:), pointer    :: wcr
      integer                              :: wcri
      character(len=*), dimension(ndims), parameter :: crd_label = [ "cr_diff_X", "cr_diff_Y", "cr_diff_Z" ]

      if (.not. has_cr) return
      if (.not.dom%has_dir(crdim)) return

      call ppp_main%start(crd_label(crdim), PPP_CR)

      if (dom%geometry_type /= GEO_XYZ) call die("[crdiffusion:cr_diff] Unsupported geometry")

      idm        = 0              ;      idm(crdim) = 1
      decr(:,:)  = 0.             ;      bcomp(:)   = 0.                 ! essential where ( .not.dom%has_dir(dim) .and. (dim /= crdim) )
      present_not_crdim = dom%has_dir .and. ( [ xdim,ydim,zdim ] /= crdim )
      wcri = wna%ind(wcr_n)

      call finest%level%restrict_to_base ! overkill
      call all_bnd ! overkill

      cgl => leaves%up_to_level(min(diff_max_lev, finest%level%l%id))%p
      do while (associated(cgl))
         cg => cgl%cg
         call cg%costs%start

         wcr => cg%w(wcri)%arr
         if (.not. associated(wcr)) call die("[crdiffusion:cr_diff] cannot get wcr")

         f1 =    eight * cg%idl(crdim)
         f2 = - oneeig * cg%idl(crdim) * dt
                                                                                               ! in case of integration with boundaries:
         ldm        = cg%ijkse(:,LO) ;      ldm(crdim) = cg%lhn(crdim,LO) + dom%D_(crdim)      ! ldm =           1 + D_
         hdm        = cg%ijkse(:,HI) ;      hdm(crdim) = cg%lhn(crdim,HI)                      ! hdm = cg%n_ + idm - D_

         do k = ldm(zdim), hdm(zdim)       ; kl = k-1 ; kh = k+1 ; kld = k-idm(zdim)
            do j = ldm(ydim), hdm(ydim)    ; jl = j-1 ; jh = j+1 ; jld = j-idm(ydim)
               do i = ldm(xdim), hdm(xdim) ; il = i-1 ; ih = i+1 ; ild = i-idm(xdim)

                  decr(crdim,:) = (cg%u(iarr_crs,i,j,k) - cg%u(iarr_crs,ild,jld,kld)) * f1
                  fcrdif = K_crs_perp * decr(crdim,:)
#ifdef MAGNETIC
                  bcomp(crdim) =  cg%b(crdim,i,j,k) * four
#endif /* MAGNETIC */
                  if (present_not_crdim(xdim)) then
                     dqm = (cg%u(iarr_crs,i ,jld,kld) + cg%u(iarr_crs,i ,j,k)) - (cg%u(iarr_crs,il,jld,kld) + cg%u(iarr_crs,il,j,k))
                     dqp = (cg%u(iarr_crs,ih,jld,kld) + cg%u(iarr_crs,ih,j,k)) - (cg%u(iarr_crs,i ,jld,kld) + cg%u(iarr_crs,i ,j,k))
                     decr(xdim,:) = (dqp+dqm) * (1.0 + sign(1.0, dqm*dqp)) * cg%idx
#ifdef MAGNETIC
                     bcomp(xdim)  = sum(cg%b(xdim,i:ih, jld:j, kld:k))
#endif /* MAGNETIC */
                  endif

                  if (present_not_crdim(ydim)) then
                     dqm = (cg%u(iarr_crs,ild,j ,kld) + cg%u(iarr_crs,i,j ,k)) - (cg%u(iarr_crs,ild,jl,kld) + cg%u(iarr_crs,i,jl,k))
                     dqp = (cg%u(iarr_crs,ild,jh,kld) + cg%u(iarr_crs,i,jh,k)) - (cg%u(iarr_crs,ild,j ,kld) + cg%u(iarr_crs,i,j ,k))
                     decr(ydim,:) = (dqp+dqm) * (1.0 + sign(1.0, dqm*dqp)) * cg%idy
#ifdef MAGNETIC
                     bcomp(ydim)  = sum(cg%b(ydim,ild:i, j:jh, kld:k))
#endif /* MAGNETIC */
                  endif

                  if (present_not_crdim(zdim)) then
                     dqm = (cg%u(iarr_crs,ild,jld,k ) + cg%u(iarr_crs,i,j,k )) - (cg%u(iarr_crs,ild,jld,kl) + cg%u(iarr_crs,i,j,kl))
                     dqp = (cg%u(iarr_crs,ild,jld,kh) + cg%u(iarr_crs,i,j,kh)) - (cg%u(iarr_crs,ild,jld,k ) + cg%u(iarr_crs,i,j,k ))
                     decr(zdim,:) = (dqp+dqm) * (1.0 + sign(1.0, dqm*dqp)) * cg%idz
#ifdef MAGNETIC
                     bcomp(zdim)  = sum(cg%b(zdim,ild:i, jld:j, k:kh))
#endif /* MAGNETIC */
                  endif

                  bb = sum(bcomp**2)
                  if (bb > epsilon(0.d0)) fcrdif = fcrdif + K_crs_paral * bcomp(crdim) * (bcomp(xdim) * decr(xdim,:) + bcomp(ydim) * decr(ydim,:) + bcomp(zdim) * decr(zdim,:)) / bb

                  wcr(:,i,j,k) = fcrdif * f2

               enddo
            enddo
         enddo

         call cg%costs%stop(I_DIFFUSE)
         cgl => cgl%nxt
      enddo

      ! As one can see in crtest with refinements, the current implementation of f/c boundary is not too accurate for diffusion.
      ! The diffusion from fine to coarse goes too fast, leading to elevated accuracy error on f/c interface (on oth sides).
      ! ToDo: Improve this:
      ! * Check if improved interpolation levels and better fcoarse-to-fine reconstruction helps.
      ! * Check whether handling f/c fluxes in wcr can help
      ! * Other ideas?
      ! Current implementation may result in better accuracy when the diffusion is restricted to base or first than unrestricted config.

      call all_wcr_boundaries

      cgl => leaves%up_to_level(min(diff_max_lev, finest%level%l%id))%p
      do while (associated(cgl))
         cg => cgl%cg
         call cg%costs%start

         ndm = cg%lhn(:,HI) - idm
         hdm = cg%lhn(:,LO) ; hdm(crdim) = cg%lhn(crdim,HI)
         ldm = hdm - idm
         p4 => cg%w(wna%fi)%span(cg%lhn(:,LO), int(ndm, kind=4))
         p4(iarr_crs,:,:,:) = p4(iarr_crs,:,:,:) - (cg%w(wcri)%span(int(cg%lhn(:,LO)+idm, kind=4), cg%lhn(:,HI)) - cg%w(wcri)%span(cg%lhn(:,LO), int(ndm, kind=4)))
         cg%u(iarr_crs,hdm(xdim):cg%lhn(xdim,HI),hdm(ydim):cg%lhn(ydim,HI),hdm(zdim):cg%lhn(zdim,HI)) = cg%u(iarr_crs,ldm(xdim):ndm(xdim),ldm(ydim):ndm(ydim),ldm(zdim):ndm(zdim)) ! for sanity

         call cg%costs%stop(I_DIFFUSE)
         cgl => cgl%nxt
      enddo

      call prolong_crs

      call ppp_main%stop(crd_label(crdim), PPP_CR)

   end subroutine cr_diff

   subroutine prolong_crs

      use cg_level_connected, only: cg_level_connected_t
      use cg_level_finest,    only: finest
      use constants,          only: GEO_XYZ, O_INJ
      use dataio_pub,         only: die
      use domain,             only: dom
      use initcosmicrays,     only: diff_max_lev, diff_prolong, iarr_crs
      use named_array_list,   only: wna, qna

      implicit none

      type(cg_level_connected_t), pointer :: diffl, curl
      integer(kind=4) :: i

      if (diff_max_lev >= finest%level%l%id) return
      if (dom%geometry_type /= GEO_XYZ) call die("[crdiffusion:prolong_crs] Only cartesian geometry is supported by diff_max_lev right now")

      diffl => finest%level
      do while (associated(diffl))
         if (diffl%l%id == diff_max_lev) exit
         diffl => diffl%coarser
      enddo
      if (.not. associated(diffl)) call die("[crdiffusion:prolong_crs] .not. associated(diffl)")
      if (diffl%l%id /= diff_max_lev) call die("[crdiffusion:prolong_crs] diffl%l%id /= diff_max_lev")

      qna%lst(qna%wai)%ord_prolong = diff_prolong
      curl => diffl
      do while (associated(curl))
         if (associated(curl%finer)) then
            if (diff_prolong /= O_INJ) call curl%level_4d_boundaries(wna%fi)
            do i = lbound(iarr_crs, dim=1, kind=4), ubound(iarr_crs, dim=1, kind=4)
               call curl%wq_copy(wna%fi, iarr_crs(i), qna%wai)
               call curl%finer%wq_copy(wna%fi, iarr_crs(i), qna%wai)  ! a bit of overkill, but current implementation may need some of this data
               call curl%prolong_1var(qna%wai)
               call curl%finer%qw_copy(qna%wai, wna%fi, iarr_crs(i))
            enddo
            call curl%prolong_bnd_from_coarser(wna%fi, nocorners=.false.)
         endif
         curl => curl%finer
      enddo

   end subroutine prolong_crs

end module crdiffusion
