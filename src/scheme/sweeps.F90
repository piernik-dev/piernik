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
!! \brief Module that implements a single sweep
!!
!! \details Here we communicate fine-coarse fluxes where needed and we perform
!! 1D-solves on all blocks.
!!
!! When a given block has any boundary with the coarse region, all flux data
!! that has to be sent right after calculation is finished. When a given block
!! has any boundary with the finer region, its calculation is delayed untill all
!! fine flux data is delivered. All communication is done quite asynchronously
!! in hope that it will nicely overlap with calculations. In some pessimistic cases
!! long stalls still may occur.
!!
!! This one is the most difficult to aggregate messages carrying flux data on fine/coarse interfaces as
!! there are complicated dependencies between grids. It is possible to calculate which fine grids should be
!! computed first in order to make critical fluxes available as early as possible.
!!
!! OPT: some fluxes can be copied locally without invoking MPI.
!<

module sweeps

! pulled by ANY

   implicit none

   private
   public  :: sweep

contains

!>
!! \brief Post a non-blocking MPI receives for all expected fluxes from fine grids.
!<

   integer(kind=4) function compute_nr_recv(cdim) result(nr)

      use constants, only: LO, HI, I_ONE
      use cg_leaves, only: leaves
      use cg_list,   only: cg_list_element
      use MPIF,      only: MPI_DOUBLE_PRECISION, MPI_Irecv
      use mpisetup,  only: comm, err_mpi, req, inflate_req

      implicit none

      integer(kind=4), intent(in)       :: cdim

      type(cg_list_element), pointer    :: cgl
      integer(kind=8), dimension(LO:HI) :: jc
      integer :: g

      nr = 0
      cgl => leaves%first
      do while (associated(cgl))
         cgl%cg%processed = .false.
         cgl%cg%finebnd(cdim, LO)%uflx(:, :, :) = 0. !> \warning overkill
         cgl%cg%finebnd(cdim, HI)%uflx(:, :, :) = 0.
         if (allocated(cgl%cg%finebnd(cdim, LO)%bflx)) cgl%cg%finebnd(cdim, LO)%bflx(:, :, :) = 0.
         if (allocated(cgl%cg%finebnd(cdim, HI)%bflx)) cgl%cg%finebnd(cdim, HI)%bflx(:, :, :) = 0.
         if (allocated(cgl%cg%rif_tgt%seg)) then
            associate ( seg => cgl%cg%rif_tgt%seg )
               do g = lbound(seg, dim=1), ubound(seg, dim=1)
                  jc = seg(g)%se(cdim, :)
                  if (jc(LO) == jc(HI)) then
                     nr = nr + I_ONE
                     if (nr > size(req, dim=1)) call inflate_req
                     call MPI_Irecv(seg(g)%buf, size(seg(g)%buf(:, :, :), kind=4), MPI_DOUBLE_PRECISION, seg(g)%proc, seg(g)%tag, comm, req(nr), err_mpi)
                     seg(g)%req => req(nr)
                  endif
               enddo
            end associate
         endif
         cgl => cgl%nxt
      enddo
   end function compute_nr_recv

!>
!! \brief Test if expected fluxes from fine grids have already arrived.
!<

   subroutine recv_cg_finebnd(cdim, cg, all_received)

      use constants,  only: LO, HI, INVALID, ORTHO1, ORTHO2, pdims, PPP_MPI
      use dataio_pub, only: die
      use fluidindex, only: flind
      use grid_cont,  only: grid_container
      use MPIF,       only: MPI_STATUS_IGNORE, MPI_Test
      use mpisetup,   only: err_mpi
      use ppp,        only: ppp_main

      implicit none

      integer(kind=4), intent(in)                  :: cdim
      type(grid_container), pointer, intent(inout) :: cg
      logical, intent(out)                         :: all_received

      integer :: g, lh
      logical(kind=4) :: received
      integer(kind=8), dimension(LO:HI) :: j1, j2, jc
      character(len=*), parameter :: recv_label = "cg_recv_fine_bnd"

      call ppp_main%start(recv_label, PPP_MPI)
      all_received = .true.
      if (allocated(cg%rif_tgt%seg)) then
         associate ( seg => cg%rif_tgt%seg )
         do g = lbound(seg, dim=1), ubound(seg, dim=1)
            jc = seg(g)%se(cdim, :)
            if (jc(LO) == jc(HI)) then
               call MPI_Test(seg(g)%req, received, MPI_STATUS_IGNORE, err_mpi)
               if (received) then
                  jc = seg(g)%se(cdim, :) !> \warning: partially duplicated code (see below)
                  j1 = seg(g)%se(pdims(cdim, ORTHO1), :)
                  j2 = seg(g)%se(pdims(cdim, ORTHO2), :)
                  if (jc(LO) /= jc(HI)) call die("[sweeps:sweep] layer too thick (Recv)")
                  if (all(cg%finebnd(cdim, LO)%index(j1(LO):j1(HI), j2(LO):j2(HI)) == jc(LO))) then
                     lh = LO
                  else if (all(cg%finebnd(cdim, HI)%index(j1(LO):j1(HI), j2(LO):j2(HI)) == jc(LO))) then
                     lh = HI
                  else
                     call die("[sweeps:sweep] Cannot determine side (Recv)")
                     lh = INVALID
                  endif
                  ! cg%finebnd(cdim, lh)%uflx(:, j1(LO):j1(HI), j2(LO):j2(HI)) = &
                  !     cg%finebnd(cdim, lh)%uflx(:, j1(LO):j1(HI), j2(LO):j2(HI)) + seg(g)%buf(:, :, :)
                  ! for more general decompositions with odd-offset patches it might be necessary to do sum, but it need to be debugged first
                  cg%finebnd(cdim, lh)%uflx(:, j1(LO):j1(HI), j2(LO):j2(HI)) = seg(g)%buf(:flind%all, :, :)
                  if (allocated(cg%finebnd(cdim, lh)%bflx)) cg%finebnd(cdim, lh)%bflx(:, j1(LO):j1(HI), j2(LO):j2(HI)) = seg(g)%buf(flind%all+1:, :, :)
               else
                  all_received = .false.
               endif
            endif
         enddo
         end associate
      endif

      call ppp_main%stop(recv_label, PPP_MPI)

   end subroutine recv_cg_finebnd

!>
!! \brief Do a non-blocking MPI Send of fluxes tfor coarse neighbors.
!<

   subroutine send_cg_coarsebnd(cdim, cg, nr)

      use constants,    only: pdims, LO, HI, ORTHO1, ORTHO2, I_ONE, INVALID, PPP_MPI
      use dataio_pub,   only: die
      use domain,       only: dom
      use grid_cont,    only: grid_container
      use grid_helpers, only: f2c_o
      use MPIF,         only: MPI_DOUBLE_PRECISION, MPI_Isend
      use mpisetup,     only: comm, err_mpi, req, inflate_req
      use ppp,          only: ppp_main

      implicit none

      integer(kind=4), intent(in)                  :: cdim
      type(grid_container), pointer, intent(inout) :: cg
      integer(kind=4), intent(inout)               :: nr

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
                  call die("[sweeps:sweep] Cannot determine side (Send)")
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

               nr = nr + I_ONE
               if (nr > size(req, dim=1)) call inflate_req
               call MPI_Isend(seg(g)%buf, size(seg(g)%buf(:, :, :), kind=4), MPI_DOUBLE_PRECISION, seg(g)%proc, seg(g)%tag, comm, req(nr), err_mpi)
               seg(g)%req => req(nr)
            endif
         enddo
         end associate
      endif

      call ppp_main%stop(send_label, PPP_MPI)

   end subroutine send_cg_coarsebnd

!>
!! \brief Call all boundaries, try to avoid unnecessary parts.
!!
!! For some reasons dir=cdim affect mcrwind tests if sweeps_mgu
!! \todo Find out why. Is it related to position of magnetic field components?
!!
!! \todo Once it gets simplified enough merge it back to sweep.
!<

   subroutine update_boundaries(cdim, istep)

      use all_boundaries, only: all_fluid_boundaries
!      use cg_leaves,      only: leaves
      use constants,      only: first_stage, DIVB_HDC
      use domain,         only: dom
      use global,         only: sweeps_mgu, integration_order, divB_0_method
#ifdef MAGNETIC
      use all_boundaries, only: all_mag_boundaries
#endif /* MAGNETIC */

      implicit none

      integer(kind=4), intent(in) :: cdim
      integer,         intent(in) :: istep

      if (divB_0_method == DIVB_HDC) then
#ifdef MAGNETIC
         call all_mag_boundaries ! ToDo: take care of psi boundaries
#endif /* MAGNETIC */
      endif

      if (dom%has_dir(cdim)) then
         if (sweeps_mgu) then
            if (istep == first_stage(integration_order)) then
               call all_fluid_boundaries(nocorners = .true., dir = cdim)
            else
               call all_fluid_boundaries(nocorners = .true.)
            endif
         else
            ! if (istep == first_stage(integration_order)) then
            !    call all_fluid_boundaries(nocorners = .true.) ! nocorners was doing something bad to periodic boundaries in Riemann with CT
            ! else
               call all_fluid_boundaries
            ! endif
         endif
      endif

   end subroutine update_boundaries

!>
!! \brief Perform 1D-solves on all blocks and send fine->coarse fluxes.
!! Update boundaries. Perform Runge-Kutta substeps.
!<

   subroutine sweep(cdim, fargo_vel)

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: ydim, ndims, first_stage, last_stage, uh_n, magh_n, psih_n, psi_n, INVALID, &
           &                      RTVD_SPLIT, RIEMANN_SPLIT, PPP_CG
      use dataio_pub,       only: die
      use global,           only: integration_order, use_fargo, which_solver
      use grid_cont,        only: grid_container
      use MPIF,             only: MPI_STATUSES_IGNORE, MPI_STATUS_IGNORE, MPI_Waitany, MPI_Waitall
      use mpisetup,         only: err_mpi, req
      use named_array_list, only: qna, wna
      use ppp,              only: ppp_main
      use solvecg_rtvd,     only: solve_cg_rtvd
      use solvecg_riemann,  only: solve_cg_riemann
      use sources,          only: prepare_sources

      implicit none

      interface

         subroutine solve_cg_sub(cg, ddim, istep, fargo_vel)

            use grid_cont, only: grid_container

            implicit none

            type(grid_container), pointer, intent(in) :: cg
            integer(kind=4),               intent(in) :: ddim
            integer,                       intent(in) :: istep     ! stage in the time integration scheme
            integer(kind=4), optional,     intent(in) :: fargo_vel

         end subroutine solve_cg_sub

      end interface

      integer(kind=4),           intent(in) :: cdim
      integer(kind=4), optional, intent(in) :: fargo_vel

      integer                        :: istep
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
      logical                        :: all_processed, all_received
      integer                        :: blocks_done
      integer(kind=4)                :: g, nr, nr_recv
      integer                        :: uhi, bhi, psii, psihi
      procedure(solve_cg_sub), pointer :: solve_cg => null()
      character(len=*), dimension(ndims), parameter :: sweep_label = [ "sweep_x", "sweep_y", "sweep_z" ]
      character(len=*), parameter :: solve_cgs_label = "solve_bunch_of_cg", cg_label = "solve_cg"

      call ppp_main%start(sweep_label(cdim))

      select case (which_solver)
         case (RTVD_SPLIT)
            solve_cg => solve_cg_rtvd
         case (RIEMANN_SPLIT)
            solve_cg => solve_cg_riemann
         case default
            call die("[sweeps:sweep] unsupported solver")
      end select

      if (use_fargo .and. cdim == ydim .and. .not. present(fargo_vel)) &
           call die("[sweeps:sweep] FARGO velocity keyword not present in y sweep")

      if (integration_order < lbound(first_stage, 1) .or. integration_order > ubound(first_stage, 1)) &
           call die("[sweeps:sweep] unknown integration_order")

      ! Despite of its name, cg%w(uhi) here contains unaltered fluid state right at
      ! the beginning of the timestep, not at half-step.
      ! For RK2, when istep==2, cg%u temporalily contains the state at half timestep.
      uhi = wna%ind(uh_n)
      bhi = INVALID
      if (wna%exists(magh_n)) bhi = wna%ind(magh_n)
      psii = INVALID
      psihi = INVALID
      if (qna%exists(psi_n)) then
         psii = qna%ind(psi_n)
         psihi = qna%ind(psih_n)
      endif

      ! for use with GLM divergence cleaning we also make a copy of b and psi fields
      cgl => leaves%first
      do while (associated(cgl))
         call prepare_sources(cgl%cg)
         cgl%cg%w(uhi)%arr = cgl%cg%u
         if (bhi  > INVALID) cgl%cg%w(bhi)%arr = cgl%cg%b
         if (psii > INVALID) cgl%cg%q(psihi)%arr = cgl%cg%q(psii)%arr
         cgl => cgl%nxt
      enddo

      ! This is the loop over Runge-Kutta stages
      do istep = first_stage(integration_order), last_stage(integration_order)
         nr_recv = compute_nr_recv(cdim)
         nr = nr_recv
         all_processed = .false.

         do while (.not. all_processed)
            all_processed = .true.
            blocks_done = 0
            ! OPT this loop should probably go from finest to coarsest for better compute-communicate overlap.
            cgl => leaves%first

            call ppp_main%start(solve_cgs_label)
            do while (associated(cgl))
               cg => cgl%cg

               if (.not. cg%processed) then
                  call recv_cg_finebnd(cdim, cg, all_received)

                  if (all_received) then
                     call ppp_main%start(cg_label, PPP_CG)
                     call solve_cg(cg, cdim, istep, fargo_vel)
                     call ppp_main%stop(cg_label, PPP_CG)

                     call send_cg_coarsebnd(cdim, cg, nr)
                     blocks_done = blocks_done + 1
                  else
                     all_processed = .false.
                  endif
               endif
               cgl => cgl%nxt
            enddo
            call ppp_main%stop(solve_cgs_label)

            if (.not. all_processed .and. blocks_done == 0) then
               if (nr_recv > 0) call MPI_Waitany(nr_recv, req(:nr_recv), g, MPI_STATUS_IGNORE, err_mpi)
               ! g is the number of completed operations
            endif
         enddo

         if (nr > 0) call MPI_Waitall(nr, req(:nr), MPI_STATUSES_IGNORE, err_mpi)

         call update_boundaries(cdim, istep)
      enddo
      call ppp_main%stop(sweep_label(cdim))

   end subroutine sweep

end module sweeps
