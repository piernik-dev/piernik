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
!! \brief Here we communicate fine-coarse fluxes obtained from a solver via cg%finebnd and cg%coarsebnd buffres
!!
!! These routines have to be called in particular order:
!!
!!     initiate_flx_recv(cdim) initializes the transfer for the specified direction and returns number of messages
!!     recv_cg_finebnd and send_cg_coarsebnd have to be called for each cg
!!     finalize_fcflx does a clean-up
!!
!! No other Piernik MPI communication should take place between calls to initiate_flx_recv and finalize_fcflx because of use the req%r array.
!!
!! See sweeps module for more hints about proper usage.
!!
!! OPT: some fluxes can be copied locally without invoking MPI.
!<

module fc_fluxes

! pulled by ANY

   implicit none

   private
   public :: initiate_flx_recv, recv_cg_finebnd, send_cg_coarsebnd

contains

!>
!! \brief Post a non-blocking MPI receives for all expected fluxes from fine grids.
!! Returns number of requests in `nr`
!<

   subroutine initiate_flx_recv(req, cdim, max_level)

      use cg_cost_data, only: I_MHD  ! ToDo: for explicit diffusion use I_DIFFUSE or I_REFINE
      use cg_leaves,    only: leaves
      use cg_list,      only: cg_list_element
      use constants,    only: LO, HI, base_level_id
      use pppmpi,       only: req_ppp

      implicit none

      type(req_ppp),             intent(inout) :: req
      integer(kind=4),           intent(in)    :: cdim
      integer(kind=4), optional, intent(in)    :: max_level

      type(cg_list_element), pointer :: cgl
      integer :: g

      call req%init(owncomm = .true., label = "fc_flx")

      nullify(cgl)
      if (present(max_level)) then  ! exclude some finest levels (useful in crdiffusion)
         if (max_level >= base_level_id) cgl => leaves%up_to_level(max_level)%p
      else  ! operate on the whole structure
         cgl => leaves%first
      endif

      do while (associated(cgl))
         call cgl%cg%costs%start

         cgl%cg%processed = .false.
         cgl%cg%finebnd(cdim, LO)%uflx(:, :, :) = 0. !> \warning overkill
         cgl%cg%finebnd(cdim, HI)%uflx(:, :, :) = 0.
         if (allocated(cgl%cg%finebnd(cdim, LO)%bflx)) cgl%cg%finebnd(cdim, LO)%bflx(:, :, :) = 0.
         if (allocated(cgl%cg%finebnd(cdim, HI)%bflx)) cgl%cg%finebnd(cdim, HI)%bflx(:, :, :) = 0.
         if (allocated(cgl%cg%rif_tgt%seg)) then
            associate ( seg => cgl%cg%rif_tgt%seg )
               do g = lbound(seg, dim=1), ubound(seg, dim=1)
                  if (seg(g)%se(cdim, LO) == seg(g)%se(cdim, HI)) call seg(g)%recv_buf(req)
               enddo
            end associate
         endif

         call cgl%cg%costs%stop(I_MHD)
         cgl => cgl%nxt
      enddo

   end subroutine initiate_flx_recv

!> \brief Test if expected fluxes from fine grids have already arrived.

   subroutine recv_cg_finebnd(req, cdim, cg, all_received)

      use constants,  only: LO, HI, INVALID, ORTHO1, ORTHO2, pdims, PPP_MPI
      use dataio_pub, only: die
      use fluidindex, only: flind
      use grid_cont,  only: grid_container
      use MPIF,       only: MPI_STATUS_IGNORE, MPI_Test, MPI_Wait
      use mpisetup,   only: err_mpi
      use ppp,        only: ppp_main
      use pppmpi,     only: req_ppp

      implicit none

      type(req_ppp),                 intent(inout) :: req
      integer(kind=4),               intent(in)    :: cdim
      type(grid_container), pointer, intent(inout) :: cg
      logical, optional,             intent(out)   :: all_received

      integer :: g, lh
      logical(kind=4) :: received
      integer(kind=8), dimension(LO:HI) :: j1, j2, jc
      character(len=*), parameter :: recv_label = "cg_recv_fine_bnd"

      call ppp_main%start(recv_label, PPP_MPI)

      if (present(all_received)) all_received = .true.
      if (allocated(cg%rif_tgt%seg)) then
         associate ( seg => cg%rif_tgt%seg )
         do g = lbound(seg, dim=1), ubound(seg, dim=1)
            jc = seg(g)%se(cdim, :)
            if (jc(LO) == jc(HI)) then
               if (present(all_received)) then
                  call MPI_Test(req%r(seg(g)%ireq), received, MPI_STATUS_IGNORE, err_mpi)
               else
                  call MPI_Wait(req%r(seg(g)%ireq), MPI_STATUS_IGNORE, err_mpi)
                  received = .true.
               endif
               if (received) then  !> \warning: partially duplicated code (see send_cg_coarsebnd())
                  j1 = seg(g)%se(pdims(cdim, ORTHO1), :)
                  j2 = seg(g)%se(pdims(cdim, ORTHO2), :)
                  if (all(cg%finebnd(cdim, LO)%index(j1(LO):j1(HI), j2(LO):j2(HI)) == jc(LO))) then
                     lh = LO
                  else if (all(cg%finebnd(cdim, HI)%index(j1(LO):j1(HI), j2(LO):j2(HI)) == jc(LO))) then
                     lh = HI
                  else
                     call die("[fc_fluxes:recv_cg_finebnd] Cannot determine side (Recv)")
                     lh = INVALID
                  endif
                  cg%finebnd(cdim, lh)%uflx(:, j1(LO):j1(HI), j2(LO):j2(HI)) = seg(g)%buf(:flind%all, :, :)
                  if (allocated(cg%finebnd(cdim, lh)%bflx)) cg%finebnd(cdim, lh)%bflx(:, j1(LO):j1(HI), j2(LO):j2(HI)) = seg(g)%buf(flind%all+1:, :, :)
               else
                  if (present(all_received)) all_received = .false.
               endif
            endif
         enddo
         end associate
      endif

      call ppp_main%stop(recv_label, PPP_MPI)

   end subroutine recv_cg_finebnd

!> \brief Do a non-blocking MPI Send of fluxes for coarse neighbors.

   subroutine send_cg_coarsebnd(req, cdim, cg)

      use constants,    only: pdims, LO, HI, ORTHO1, ORTHO2, INVALID, PPP_MPI
      use dataio_pub,   only: die
      use domain,       only: dom
      use grid_cont,    only: grid_container
      use grid_helpers, only: f2c_o
      use ppp,          only: ppp_main
      use pppmpi,       only: req_ppp

      implicit none

      type(req_ppp),                 intent(inout) :: req
      integer(kind=4),               intent(in)    :: cdim
      type(grid_container), pointer, intent(inout) :: cg

      integer :: g, lh
      integer(kind=8), dimension(LO:HI) :: j1, j2, jc
      integer(kind=8) :: j, k
      character(len=*), parameter :: send_label = "cg_send_coarse_bnd"

      call ppp_main%start(send_label, PPP_MPI)

      if (allocated(cg%rof_tgt%seg)) then
         associate ( seg => cg%rof_tgt%seg )
         do g = lbound(seg, dim=1), ubound(seg, dim=1)
            jc = seg(g)%se(cdim, :) !> \warning: partially duplicated code (see above)
            if (jc(LO) == jc(HI)) then
               j1 = seg(g)%se(pdims(cdim, ORTHO1), :)
               j2 = seg(g)%se(pdims(cdim, ORTHO2), :)
               if (all(cg%coarsebnd(cdim, LO)%index(j1(LO):j1(HI), j2(LO):j2(HI)) == jc(LO))) then
                  lh = LO
               else if (all(cg%coarsebnd(cdim, HI)%index(j1(LO):j1(HI), j2(LO):j2(HI)) == jc(LO))) then
                  lh = HI
               else
                  call die("[fc_fluxes:send_cg_coarsebnd] Cannot determine side (Send)")
                  lh = INVALID
               endif

               seg(g)%buf(:, :, :) = 0.
               do j = j1(LO), j1(HI)
                  do k = j2(LO), j2(HI)
                     if (allocated(cg%coarsebnd(cdim, lh)%bflx)) then
                        seg(g)%buf(:, f2c_o(j), f2c_o(k)) = seg(g)%buf(:, f2c_o(j), f2c_o(k)) + [ cg%coarsebnd(cdim, lh)%uflx(:, j, k), cg%coarsebnd(cdim, lh)%bflx(:, j, k) ]
                     else
                        seg(g)%buf(:, f2c_o(j), f2c_o(k)) = seg(g)%buf(:, f2c_o(j), f2c_o(k)) + cg%coarsebnd(cdim, lh)%uflx(:, j, k)
                     endif
                  enddo
               enddo
               seg(g)%buf = 1/2.**(dom%eff_dim-1) * seg(g)%buf
               call seg(g)%send_buf(req)
            endif
         enddo
         end associate
      endif

      call ppp_main%stop(send_label, PPP_MPI)

   end subroutine send_cg_coarsebnd

end module fc_fluxes
