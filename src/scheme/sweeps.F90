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
!! \details Here we perform 1D-solves on all blocks including proper exchange of f/c fluxes.
!!
!! When a given block has any boundary with the coarse region, all flux data
!! that has to be sent right after calculation is finished. When a given block
!! has any boundary with the finer region, its calculation is delayed until all
!! fine flux data is delivered. All communication is done quite asynchronously
!! in hope that it will nicely overlap with calculations. In some pessimistic cases
!! long stalls still may occur.
!!
!! This one is the most difficult to aggregate messages carrying flux data on fine/coarse interfaces as
!! there are complicated dependencies between grids. It is possible to calculate which fine grids should be
!! computed first in order to make critical fluxes available as early as possible.
!<

module sweeps

! pulled by ANY

   implicit none

   private
   public :: sweep

contains

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

      if (dom%has_dir(cdim)) then
         if (sweeps_mgu) then
            if (istep == first_stage(integration_order)) then
               call all_fluid_boundaries(nocorners = .true., dir = cdim)
            else
               call all_fluid_boundaries(nocorners = .true.)
            endif
         else
            ! nocorners and dir = cdim can be used safely only when ord_fluid_prolong == 0 .and. cc_mag
            ! essential speedups here are possible but it requires c/f boundary prolongation that does not require corners

            ! if (istep == first_stage(integration_order)) then
            !    call all_fluid_boundaries(nocorners = .true.)
            ! else
               call all_fluid_boundaries !(nocorners = .true., dir = cdim)
            ! endif
         endif
      endif

      if (divB_0_method == DIVB_HDC) then
#ifdef MAGNETIC
         call all_mag_boundaries ! ToDo: take care of psi boundaries
#endif /* MAGNETIC */
      endif

   end subroutine update_boundaries

!>
!! \brief Perform 1D-solves on all blocks and send fine->coarse fluxes.
!! Update boundaries. Perform Runge-Kutta substeps.
!<

   subroutine sweep(cdim, fargo_vel)

      use cg_cost_data,     only: I_MHD, I_REFINE
      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use cg_list_dataop,   only: cg_list_dataop_t
      use constants,        only: ydim, ndims, first_stage, last_stage, uh_n, magh_n, psih_n, psi_n, INVALID, &
           &                      RTVD_SPLIT, RIEMANN_SPLIT, PPP_CG
      use dataio_pub,       only: die
      use fc_fluxes,        only: initiate_flx_recv, recv_cg_finebnd, send_cg_coarsebnd
      use global,           only: integration_order, use_fargo, which_solver
      use grid_cont,        only: grid_container
      use MPIF,             only: MPI_STATUS_IGNORE
      use MPIFUN,           only: MPI_Waitany
      use mpisetup,         only: err_mpi
      use named_array_list, only: qna, wna
      use ppp,              only: ppp_main
      use pppmpi,           only: req_ppp
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

      integer                          :: istep
      type(cg_list_element), pointer   :: cgl
      type(grid_container),  pointer   :: cg
      type(cg_list_dataop_t), pointer  :: sl
      type(req_ppp)                    :: req
      logical                          :: all_processed, all_received
      integer                          :: blocks_done
      integer(kind=4)                  :: n_recv, g
      integer                          :: uhi, bhi, psii, psihi
      procedure(solve_cg_sub), pointer :: solve_cg => null()
      character(len=*), dimension(ndims), parameter :: sweep_label = [ "sweep_x", "sweep_y", "sweep_z" ]
      character(len=*), parameter :: solve_cgs_label = "solve_bunch_of_cg", cg_label = "solve_cg", init_src_label = "init_src"

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
      ! For RK2, when istep==2, cg%u temporarily contains the state at half timestep.
      uhi = wna%ind(uh_n)
      bhi = INVALID
      if (wna%exists(magh_n)) bhi = wna%ind(magh_n)
      psii = INVALID
      psihi = INVALID
      if (qna%exists(psi_n)) then
         psii = qna%ind(psi_n)
         psihi = qna%ind(psih_n)
      endif

      sl => leaves%prioritized_cg(cdim, covered_too = .true.)
      ! We can't just skip the covered cg because it affects divvel (or
      ! other things that rely on data computed on coarse cells and are not
      ! restricted from fine blocks).
!      sl => leaves%leaf_only_cg()

      call ppp_main%start(init_src_label)
      ! for use with GLM divergence cleaning we also make a copy of b and psi fields
      cgl => leaves%first
      do while (associated(cgl))
         call prepare_sources(cgl%cg)
         cgl => cgl%nxt
      enddo
      call ppp_main%stop(init_src_label)

      ! This is the loop over Runge-Kutta stages
      do istep = first_stage(integration_order), last_stage(integration_order)

         call initiate_flx_recv(req, cdim)
         n_recv = req%n
         all_processed = .false.

         do while (.not. all_processed)
            all_processed = .true.
            blocks_done = 0
            ! OPT this loop should probably go from finest to coarsest for better compute-communicate overlap.
            cgl => sl%first

            call ppp_main%start(solve_cgs_label)
            do while (associated(cgl))
               cg => cgl%cg
               call cg%costs%start

               if (.not. cg%processed) then
                  call recv_cg_finebnd(req, cdim, cg, all_received)

                  if (all_received) then
                     call ppp_main%start(cg_label, PPP_CG)
                     call cg%costs%stop(I_REFINE)
                     ! The recv_cg_finebnd and send_cg_coarsebnd aren't MHD, so we should count them separately.
                     ! The tricky part is that we need to fit all the switching inside the conditional part
                     ! and don't mess pairing and don't let them to nest.

                     call cg%costs%start
                     call solve_cg(cg, cdim, istep, fargo_vel)
                     call cg%costs%stop(I_MHD)

                     call ppp_main%stop(cg_label, PPP_CG)

                     call cg%costs%start
                     call send_cg_coarsebnd(req, cdim, cg)
                     blocks_done = blocks_done + 1
                  else
                     all_processed = .false.
                  endif
               endif

               call cg%costs%stop(I_REFINE)
               cgl => cgl%nxt
            enddo
            call ppp_main%stop(solve_cgs_label)

            if (.not. all_processed .and. blocks_done == 0) then
               if (n_recv > 0) call MPI_Waitany(n_recv, req%r(:n_recv), g, MPI_STATUS_IGNORE, err_mpi)
               ! g is the number of completed operations
            endif
         enddo

         call req%waitall("sweeps")

         call update_boundaries(cdim, istep)
      enddo

      call sl%delete
      deallocate(sl)

      call ppp_main%stop(sweep_label(cdim))

   end subroutine sweep

end module sweeps
