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
   public :: init_crdiffusion, make_diff_sweeps

   logical :: has_cr

contains

!> \brief Set up CR-specific auxiliary arrays

   subroutine init_crdiffusion

      use cg_list_global,   only: all_cg
      use constants,        only: wcr_n
      use crhelpers,        only: divv_i, divv_n
      use initcosmicrays,   only: ord_cr_prolong
      use dataio_pub,       only: warn
      use fluidindex,       only: flind
      use named_array_list, only: qna

      implicit none

      has_cr = (flind%crs%all > 0)

      if (has_cr) then
         call all_cg%reg_var(wcr_n, dim4 = flind%crs%all, ord_prolong = ord_cr_prolong)
         ! Dirty trick: Ensure that all prolongation buffers are ready for ord_cr_prolong,
         ! wcr_n itself doesn't need it.
      else
         call warn("[crdiffusion:init_crdiffusion] No CR species to diffuse")
      endif
      call all_cg%reg_var(divv_n)
      divv_i = qna%ind(divv_n)

   end subroutine init_crdiffusion

!> \brief Execute diffusion

   subroutine make_diff_sweeps(forward)

      use all_boundaries, only: all_fluid_boundaries, all_mag_boundaries
      use constants,      only: xdim, zdim, I_ONE, PPP_CR
      use global,         only: skip_sweep
      use ppp,            only: ppp_main

      implicit none

      logical, intent(in) :: forward  !< order of sweeps: XYZ or ZYX

      integer(kind=4) :: s
      character(len=*), parameter :: crdiff_label = "CR_diffusion"

      call ppp_main%start(crdiff_label, PPP_CR)

      call all_mag_boundaries

      do s = merge(xdim, zdim, forward), merge(zdim, xdim, forward), merge(I_ONE, -I_ONE, forward)
         if (.not. skip_sweep(s)) call make_diff_sweep(s)
      enddo

      ! This call prevents occurence of SIGFPE in the Riemann solver
      ! Strange thing is that it is not fully deterministic and sometimes the code may work well without this call
      call all_fluid_boundaries  ! overkill?

      call ppp_main%stop(crdiff_label, PPP_CR)

   contains

      !> \brief Perform single diffusion sweep in forward or backward direction

      subroutine make_diff_sweep(dir)

#ifdef DEBUG
         use piernikiodebug, only: force_dumps
#endif /* DEBUG */

         implicit none

         integer(kind=4), intent(in) :: dir !< direction, one of xdim, ydim, zdim

         call cr_diff(dir)

#ifdef DEBUG
         call force_dumps
#endif /* DEBUG */

      end subroutine make_diff_sweep

   end subroutine make_diff_sweeps

!>
!! \brief boundaries for wcr
!! This procedure is a shameless copy of cg_list_bnd%arr3d_boundaries adapted for wcr
!! \todo REMOVE ME, FIX ME, MERGE ME with cg_list_bnd%arr3d_boundaries or at least have a decency to make me more general
!<
   subroutine all_wcr_boundaries(crdim)

      use cg_cost_data,     only: I_DIFFUSE
      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: ndims, xdim, ydim, zdim, LO, HI, BND_PER, BND_MPI, BND_FC, BND_MPI_FC, I_TWO, I_THREE, wcr_n, PPP_CR
      use dataio_pub,       only: die
      use domain,           only: dom
      use grid_cont,        only: grid_container
      use named_array_list, only: wna
      use ppp,              only: ppp_main

      implicit none

      integer(kind=4), intent(in)             :: crdim

      integer(kind=4)                         :: i, d, lh
      integer(kind=4), dimension(ndims,LO:HI) :: l, r
      real, dimension(:,:,:,:), pointer       :: wcr
      type(cg_list_element),    pointer       :: cgl
      type(grid_container),     pointer       :: cg
      character(len=*), parameter :: awb_label = "all_wcr_boundaries"

      if (.not. has_cr) return

      call ppp_main%start(awb_label, PPP_CR)

      call leaves%leaf_arr4d_boundaries(wna%ind(wcr_n), no_fc = .true., dir=crdim)  ! skip coarse-to-fine prolongation as it doesn't work well for fluxes
      ! ToDo: override coarse fluxes with restricted fine fluxes

      ! do the external boundaries
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         call cg%costs%start

         wcr => cg%w(wna%ind(wcr_n))%arr
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

      use all_boundaries,   only: all_fluid_boundaries
      use cg_cost_data,     only: I_DIFFUSE
      use cg_leaves,        only: leaves
      use cg_level_finest,  only: finest
      use cg_list,          only: cg_list_element
      use constants,        only: xdim, ydim, zdim, ndims, LO, HI, oneeig, eight, wcr_n, GEO_XYZ, PPP_CR, half
      use dataio_pub,       only: die
      use domain,           only: dom
      use fc_fluxes,        only: compute_nr_recv, recv_cg_finebnd, send_cg_coarsebnd, finalize_fcflx
      use fluidindex,       only: flind
      use global,           only: dt
      use grid_cont,        only: grid_container
      use initcosmicrays,   only: iarr_crs, K_crs_paral, K_crs_perp, ord_cr_prolong
      use named_array,      only: p4
      use named_array_list, only: wna
      use ppp,              only: ppp_main
      use ppp_mpi,          only: piernik_Waitall
#ifdef MAGNETIC
      use constants,        only: four
      use global,           only: cc_mag
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
      integer(kind=4) :: nr, nr_recv
      logical :: all_received
      integer(kind=4) :: ord_save
      real, parameter :: flux_factor = half

      if (.not. has_cr) return
      if (.not.dom%has_dir(crdim)) return

      call ppp_main%start(crd_label(crdim), PPP_CR)

      if (dom%geometry_type /= GEO_XYZ) call die("[crdiffusion:cr_diff] Unsupported geometry")

      idm        = 0              ;      idm(crdim) = 1
      decr(:,:)  = 0.             ;      bcomp(:)   = 0.                 ! essential where ( .not.dom%has_dir(dim) .and. (dim /= crdim) )
      present_not_crdim = dom%has_dir .and. ( [ xdim,ydim,zdim ] /= crdim )
      wcri = wna%ind(wcr_n)

      ! Dirty trick: enforce prolongation order for CR in case someone wants smoother estimates for f/c ifnterpolation of cg%u(iarr_crs,:,:,:)
      ord_save = wna%lst(wna%fi)%ord_prolong
      wna%lst(wna%fi)%ord_prolong = ord_cr_prolong
      call finest%level%restrict_to_base  ! overkill?
      call all_fluid_boundaries  ! overkill?
      wna%lst(wna%fi)%ord_prolong = ord_save

      nr_recv = compute_nr_recv(crdim)
      nr = nr_recv

      cgl => leaves%first
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
                  if (cc_mag) then
                     bcomp(:) =  cg%b(:, i, j, k) + cg%b(:, ild, jld, kld)
                  else
                     bcomp(crdim) =  cg%b(crdim,i,j,k) * four
                  endif
#endif /* MAGNETIC */
                  if (present_not_crdim(xdim)) then
                     dqm = (cg%u(iarr_crs,i ,jld,kld) + cg%u(iarr_crs,i ,j,k)) - (cg%u(iarr_crs,il,jld,kld) + cg%u(iarr_crs,il,j,k))
                     dqp = (cg%u(iarr_crs,ih,jld,kld) + cg%u(iarr_crs,ih,j,k)) - (cg%u(iarr_crs,i ,jld,kld) + cg%u(iarr_crs,i ,j,k))
                     decr(xdim,:) = (dqp+dqm) * (1.0 + sign(1.0, dqm*dqp)) * cg%idx
#ifdef MAGNETIC
                     if (.not. cc_mag) bcomp(xdim)  = sum(cg%b(xdim,i:ih, jld:j, kld:k))
#endif /* MAGNETIC */
                  endif

                  if (present_not_crdim(ydim)) then
                     dqm = (cg%u(iarr_crs,ild,j ,kld) + cg%u(iarr_crs,i,j ,k)) - (cg%u(iarr_crs,ild,jl,kld) + cg%u(iarr_crs,i,jl,k))
                     dqp = (cg%u(iarr_crs,ild,jh,kld) + cg%u(iarr_crs,i,jh,k)) - (cg%u(iarr_crs,ild,j ,kld) + cg%u(iarr_crs,i,j ,k))
                     decr(ydim,:) = (dqp+dqm) * (1.0 + sign(1.0, dqm*dqp)) * cg%idy
#ifdef MAGNETIC
                     if (.not. cc_mag) bcomp(ydim)  = sum(cg%b(ydim,ild:i, j:jh, kld:k))
#endif /* MAGNETIC */
                  endif

                  if (present_not_crdim(zdim)) then
                     dqm = (cg%u(iarr_crs,ild,jld,k ) + cg%u(iarr_crs,i,j,k )) - (cg%u(iarr_crs,ild,jld,kl) + cg%u(iarr_crs,i,j,kl))
                     dqp = (cg%u(iarr_crs,ild,jld,kh) + cg%u(iarr_crs,i,j,kh)) - (cg%u(iarr_crs,ild,jld,k ) + cg%u(iarr_crs,i,j,k ))
                     decr(zdim,:) = (dqp+dqm) * (1.0 + sign(1.0, dqm*dqp)) * cg%idz
#ifdef MAGNETIC
                     if (.not. cc_mag) bcomp(zdim)  = sum(cg%b(zdim,ild:i, jld:j, k:kh))
#endif /* MAGNETIC */
                  endif

                  bb = sum(bcomp**2)
                  if (bb > epsilon(0.d0)) fcrdif = fcrdif + K_crs_paral * bcomp(crdim) * (bcomp(xdim) * decr(xdim,:) + bcomp(ydim) * decr(ydim,:) + bcomp(zdim) * decr(zdim,:)) / bb

                  wcr(:,i,j,k) = fcrdif * f2

               enddo
            enddo
         enddo

         ! Very general, perhaps not very efficient approach
         ! It can be optimized as there is either something to exchange or not (we're on fine side) and array assignment can be made
         do j = lbound(cg%coarsebnd(crdim, LO)%uflx, dim=2), ubound(cg%coarsebnd(crdim, LO)%uflx, dim=2)  ! bounds for (crdim, HI) are the same
            do k = lbound(cg%coarsebnd(crdim, LO)%uflx, dim=3), ubound(cg%coarsebnd(crdim, LO)%uflx, dim=3)
               select case (crdim)
                  case (xdim)
                     if (cg%coarsebnd(crdim, LO)%index(j, k) >= cg%is) &
                          cg%coarsebnd(crdim, LO)%uflx(:flind%crs%all, j, k) = wcr(:, cg%coarsebnd(crdim, LO)%index(j, k), j, k)
                     if (cg%coarsebnd(crdim, HI)%index(j, k) <= cg%ie) &
                          cg%coarsebnd(crdim, HI)%uflx(:flind%crs%all, j, k) = wcr(:, cg%coarsebnd(crdim, HI)%index(j, k)+1, j, k)
                  case (ydim)
                     if (cg%coarsebnd(crdim, LO)%index(j, k) >= cg%js) &
                          cg%coarsebnd(crdim, LO)%uflx(:flind%crs%all, j, k) = wcr(:, k, cg%coarsebnd(crdim, LO)%index(j, k), j)
                     if (cg%coarsebnd(crdim, HI)%index(j, k) <= cg%je) &
                          cg%coarsebnd(crdim, HI)%uflx(:flind%crs%all, j, k) = wcr(:, k, cg%coarsebnd(crdim, HI)%index(j, k)+1, j)
                  case (zdim)
                     if (cg%coarsebnd(crdim, LO)%index(j, k) >= cg%ks) &
                          cg%coarsebnd(crdim, LO)%uflx(:flind%crs%all, j, k) = wcr(:, j, k, cg%coarsebnd(crdim, LO)%index(j, k))
                     if (cg%coarsebnd(crdim, HI)%index(j, k) <= cg%ke) &
                          cg%coarsebnd(crdim, HI)%uflx(:flind%crs%all, j, k) = wcr(:, j, k, cg%coarsebnd(crdim, HI)%index(j, k)+1)
                  case default
                     call die("[crdiffusion:cr_diff] What to send?")
               end select
            enddo
         enddo

         call send_cg_coarsebnd(crdim, cg, nr)

         call cg%costs%stop(I_DIFFUSE)
         cgl => cgl%nxt
      enddo

      call piernik_Waitall(nr, "cr_diff")
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         call cg%costs%start

         call recv_cg_finebnd(crdim, cg, all_received)
         if (.not. all_received) call die("[crdiffusion:cr_diff] incomplete fluxes")
         call cg%costs%stop(I_DIFFUSE)
         cgl => cgl%nxt
      enddo

      call finalize_fcflx

      call all_wcr_boundaries(crdim)

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         call cg%costs%start

         ! Very general, perhaps not very efficient approach, again
         ! It can be optimized as there are at most four different regions to exchange (we're on coarse side) and array assignments can be made

         ! To Do: Figure out where the flux factor value come from. And why it does not depend on dimensionality.
         wcr => cg%w(wcri)%arr
         do j = lbound(cg%finebnd(crdim, LO)%uflx, dim=2), ubound(cg%finebnd(crdim, LO)%uflx, dim=2)  ! bounds for (crdim, HI) are the same
            do k = lbound(cg%finebnd(crdim, LO)%uflx, dim=3), ubound(cg%finebnd(crdim, LO)%uflx, dim=3)
               select case (crdim)
                  case (xdim)
                     if (cg%finebnd(crdim, LO)%index(j, k) > cg%is) &
                          wcr(:, cg%finebnd(crdim, LO)%index(j, k) + 1, j, k) = cg%finebnd(crdim, LO)%uflx(:flind%crs%all, j, k) * flux_factor
                     if (cg%finebnd(crdim, HI)%index(j, k) < cg%ie) &
                          wcr(:, cg%finebnd(crdim, HI)%index(j, k), j, k) = cg%finebnd(crdim, HI)%uflx(:flind%crs%all, j, k) * flux_factor
                  case (ydim)
                     if (cg%finebnd(crdim, LO)%index(j, k) >= cg%js) &
                          wcr(:, k, cg%finebnd(crdim, LO)%index(j, k) + 1, j) = cg%finebnd(crdim, LO)%uflx(:flind%crs%all, j, k) * flux_factor
                     if (cg%finebnd(crdim, HI)%index(j, k) <= cg%je) &
                          wcr(:, k, cg%finebnd(crdim, HI)%index(j, k), j) = cg%finebnd(crdim, HI)%uflx(:flind%crs%all, j, k) * flux_factor
                  case (zdim)
                     if (cg%finebnd(crdim, LO)%index(j, k) >= cg%ks) &
                          wcr(:, j, k, cg%finebnd(crdim, LO)%index(j, k) + 1) = cg%finebnd(crdim, LO)%uflx(:flind%crs%all, j, k) * flux_factor
                     if (cg%finebnd(crdim, HI)%index(j, k) <= cg%ke) &
                          wcr(:, j, k, cg%finebnd(crdim, HI)%index(j, k)) = cg%finebnd(crdim, HI)%uflx(:flind%crs%all, j, k) * flux_factor
                  case default
                     call die("[crdiffusion:cr_diff] What to receive?")
               end select
            enddo
         enddo

         ndm = cg%lhn(:,HI) - idm
         hdm = cg%lhn(:,LO) ; hdm(crdim) = cg%lhn(crdim,HI)
         ldm = hdm - idm
         p4 => cg%w(wna%fi)%span(cg%lhn(:,LO), int(ndm, kind=4))
         p4(iarr_crs,:,:,:) = p4(iarr_crs,:,:,:) - (cg%w(wcri)%span(int(cg%lhn(:,LO)+idm, kind=4), cg%lhn(:,HI)) - cg%w(wcri)%span(cg%lhn(:,LO), int(ndm, kind=4)))
         cg%u(iarr_crs,hdm(xdim):cg%lhn(xdim,HI),hdm(ydim):cg%lhn(ydim,HI),hdm(zdim):cg%lhn(zdim,HI)) = cg%u(iarr_crs,ldm(xdim):ndm(xdim),ldm(ydim):ndm(ydim),ldm(zdim):ndm(zdim)) ! for sanity

         call cg%costs%stop(I_DIFFUSE)
         cgl => cgl%nxt
      enddo

      call ppp_main%stop(crd_label(crdim), PPP_CR)

   end subroutine cr_diff

end module crdiffusion
