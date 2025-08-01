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
      use constants,      only: first_stage, DIVB_HDC,xdim,zdim, UNSPLIT
      use domain,         only: dom
      use global,         only: sweeps_mgu, integration_order, divB_0_method, which_solver
#ifdef MAGNETIC
      use all_boundaries, only: all_mag_boundaries
#endif /* MAGNETIC */

      implicit none

      integer,                  intent(in) :: istep

      integer(kind=4)                      :: ub_i

      do ub_i=xdim,zdim
            if (dom%has_dir(ub_i)) then
               if (sweeps_mgu) then
                  if (istep == first_stage(integration_order)) then
                     call all_fluid_boundaries(nocorners = .true., dir = ub_i)
                  else
                     call all_fluid_boundaries(nocorners = .true.)
                  endif
               else
                  ! nocorners and dir = cdim can be used safely only when ord_fluid_prolong == 0 .and. cc_mag
                  ! essential speedups here are possible but it requires c/f boundary prolongation that does not require corners

                  ! if (istep == first_stage(integration_order)) then
                  !    call all_fluid_boundaries(nocorners = .true.)
                  ! else
                     call all_fluid_boundaries(istep=istep) !(nocorners = .true., dir = cdim)
                  ! endif
               endif
            endif
         enddo
      if (divB_0_method == DIVB_HDC) then
#ifdef MAGNETIC
         if (which_solver==UNSPLIT) then
            call all_mag_boundaries(istep) ! ToDo: take care of psi boundaries
         else
            call all_mag_boundaries ! ToDo: take care of psi boundaries
         endif
#endif /* MAGNETIC */
      endif

   end subroutine update_boundaries

    subroutine unsplit_sweep()
        use cg_list,                            only: cg_list_element
        use cg_cost_data,                       only: I_MHD, I_REFINE
        use grid_cont,                          only: grid_container
        use cg_leaves,                          only: leaves
        use dataio_pub,                         only: die
        use MPIF,                               only: MPI_STATUS_IGNORE
        use MPIFUN,                             only: MPI_Waitany
        use mpisetup,                           only: err_mpi
        use solvecg_unsplit,                    only: solve_cg_unsplit
        use sources,                            only: prepare_sources
        use global,                             only: integration_order, which_solver 
        use constants,                          only: first_stage, last_stage, UNSPLIT, PPP_CG
        use cg_list_dataop,                     only: cg_list_dataop_t
        use pppmpi,                             only: req_ppp
        use MPIF,                               only: MPI_STATUS_IGNORE
        use MPIFUN,                             only: MPI_Waitany
        use fc_fluxes,                          only: initiate_flx_recv, recv_cg_finebnd, send_cg_coarsebnd


        implicit none


        type(cg_list_dataop_t), pointer  :: sl
        type(req_ppp)                    :: req
        logical                          :: all_processed, all_received
        integer                          :: blocks_done
        integer(kind=4)                  :: n_recv, g
        type(cg_list_element), pointer   :: cgl
        type(grid_container),  pointer   :: cg
        integer                          :: istep
        character(len=*), parameter :: solve_cgs_label = "solve_bunch_of_cg", cg_label = "solve_cg", init_src_label = "init_src"

        if (which_solver /= UNSPLIT) call die("[unsplit_sweeps:unsplit_sweep] Only compatible with UNSPLIT solver")
        sl => leaves%prioritized_cg(-1, covered_too = .true.)


        cgl => leaves%first
        do while (associated(cgl))
            !cgl%cg%f(:,:,:,:) = 0.0                         ! This looks unsafe. Might come back to bite your ass later. Careful 
            !cgl%cg%g(:,:,:,:) = 0.0
            !cgl%cg%h(:,:,:,:) = 0.0
            call prepare_sources(cgl%cg)
            cgl => cgl%nxt
        enddo



        !all_processed = .false.
        !blocks_done = 0
        
        do istep = first_stage(integration_order), last_stage(integration_order)
            
            call initiate_flx_recv(req, -1)
            n_recv = req%n
            all_processed = .false.
            do while (.not. all_processed)
                all_processed = .true.
                blocks_done = 0
                ! OPT this loop should probably go from finest to coarsest for better compute-communicate overlap.
                cgl => sl%first
            !cgl => leaves%first    ! This should be modified to point to the first leaf of prioritized leaves . Maybe important only in AMR ? 

                do while (associated(cgl))
                    cg => cgl%cg

                    if (.not. cg%processed) then
                        call recv_cg_finebnd(req,-1, cg, all_received)
                        if (all_received) then

                            call cg%cleanup_flux()

                            call solve_cg_unsplit(cg,istep)

                            call send_cg_coarsebnd(req, -1, cg)
                            blocks_done = blocks_done + 1
                        else
                            all_processed = .false.
                        endif
                    endif

                    cgl => cgl%nxt
                end do
                if (.not. all_processed .and. blocks_done == 0) then
                if (n_recv > 0) call MPI_Waitany(n_recv, req%r(:n_recv), g, MPI_STATUS_IGNORE, err_mpi)
                ! g is the number of completed operations
                endif
            enddo

            call req%waitall("sweeps")

            call update_boundaries(istep)
        end do
        call sl%delete
        deallocate(sl)


    end subroutine unsplit_sweep   

end module unsplit_sweeps