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
!! \brief The job of this module is simple : Pass a  block of cg to solve to do a unsplit update of the state
!! Currently we dont not add AMR support or ppp monitoring. Sister module sweeps is used for directional sweep update and is called by
!! fluid update module. We will call this module from fluid_unsplit_update which is in turn mentioned in fluid_update to keep this line of
!! additions away from the main code and merger it later. We are not adding fargo support either. This will be the first update after this works
!<

module unsplit_sweeps

! pulled by ANY

   implicit none

   private
   public :: unsplit_sweep

contains

   subroutine update_boundaries(istep)

      use all_boundaries, only: all_fluid_boundaries
!      use cg_leaves,      only: leaves
      use constants,      only: first_stage, DIVB_HDC, xdim, zdim
      use domain,         only: dom
      use global,         only: sweeps_mgu, integration_order, divB_0_method
#ifdef MAGNETIC
      use all_boundaries, only: all_mag_boundaries
#endif /* MAGNETIC */

      implicit none

      integer, intent(in) :: istep

      integer(kind=4) :: ub_i

      if (sweeps_mgu) then
         if (istep == first_stage(integration_order)) then
            do ub_i = xdim, zdim
               if (.not. dom%has_dir(ub_i)) cycle
               call all_fluid_boundaries(nocorners = .true., dir = ub_i, istep = istep)
            enddo
         else
            call all_fluid_boundaries(nocorners = .true., istep = istep)
         endif
      else
         call all_fluid_boundaries(istep=istep)
      endif
      if (divB_0_method == DIVB_HDC) then
#ifdef MAGNETIC
         call all_mag_boundaries(istep) ! ToDo: take care of psi boundaries
#endif /* MAGNETIC */
      endif

   end subroutine update_boundaries

   subroutine unsplit_sweep()

      use cg_cost_data,      only: I_MHD, I_REFINE
      use cg_leaves,         only: leaves
      use cg_list,           only: cg_list_element
      use cg_list_dataop,    only: cg_list_dataop_t
      use constants,         only: first_stage, last_stage, INVALID, PPP_CG, RIEMANN_UNSPLIT
      use dataio_pub,        only: die
      use fc_fluxes_unsplit, only: initiate_flx_recv, recv_cg_finebnd, send_cg_coarsebnd
      use global,            only: integration_order, which_solver
      use grid_cont,         only: grid_container
      use MPIF,              only: MPI_STATUS_IGNORE
      use MPIFUN,            only: MPI_Waitany
      use mpisetup,          only: err_mpi
      use ppp,               only: ppp_main
      use pppmpi,            only: req_ppp
      use sources,           only: prepare_sources
      use solvecg_unsplit,   only: solve_cg_unsplit

      implicit none

      integer                          :: istep
      type(cg_list_element), pointer   :: cgl
      type(grid_container),  pointer   :: cg
      type(cg_list_dataop_t), pointer  :: sl
      type(req_ppp)                    :: req
      logical                          :: all_processed, all_received
      integer                          :: blocks_done
      integer(kind=4)                  :: n_recv, g
      character(len=*), parameter      :: solve_cgs_label = "solve_bunch_of_cg", cg_label = "solve_cg", init_src_label = "init_src"

      call ppp_main%start("unsplit_sweep")

      if (which_solver /= RIEMANN_UNSPLIT) call die("[unsplit_sweeps:unsplit_sweep] Only compatible with UNSPLIT solver")

      sl => leaves%prioritized_cg(INVALID, covered_too = .true.)

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

         call initiate_flx_recv(req, INVALID)
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
                  call recv_cg_finebnd(req, INVALID, cg, all_received)

                  if (all_received) then
                     call ppp_main%start(cg_label, PPP_CG)
                     call cg%costs%stop(I_REFINE)
                     ! The recv_cg_finebnd and send_cg_coarsebnd aren't MHD, so we should count them separately.
                     ! The tricky part is that we need to fit all the switching inside the conditional part
                     ! and don't mess pairing and don't let them to nest.

                     call cg%cleanup_flux()      ! Seems unnecessary.This just sets the flux array to 0.0

                     call cg%costs%start
                     call solve_cg_unsplit(cg, istep)
                     call cg%costs%stop(I_MHD)

                     call ppp_main%stop(cg_label, PPP_CG)

                     call cg%costs%start
                     call send_cg_coarsebnd(req, INVALID, cg)
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

         call update_boundaries(istep)
      enddo

      call sl%delete
      deallocate(sl)

      call ppp_main%stop("unsplit_sweep")

   end subroutine unsplit_sweep

end module unsplit_sweeps
