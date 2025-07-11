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

    subroutine unsplit_sweep()
        use cg_list,                            only: cg_list_element
        use grid_cont,                          only: grid_container
        use cg_leaves,                          only: leaves
        use dataio_pub,                         only: die
        use MPIF,                               only: MPI_STATUS_IGNORE
        use MPIFUN,                             only: MPI_Waitany
        use mpisetup,                           only: err_mpi
        use solvecg_unsplit,                    only: solve_cg_unsplit
        use sources,                            only: prepare_sources
        use global,                             only: integration_order, which_solver 
        use constants,                          only: first_stage, last_stage, UNSPLIT
        use unsplit_update_boundary,            only: update_boundaries
        implicit none



        type(cg_list_element), pointer   :: cgl
        type(grid_container),  pointer   :: cg
        integer                          :: istep

        if (which_solver /= UNSPLIT) call die("[unsplit_sweeps:unsplit_sweep] Only compatible with UNSPLIT solver")

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
            
            cgl => leaves%first    ! This should be modified to point to the first leaf of prioritized leaves . Maybe important only in AMR ? 

            do while (associated(cgl))
                cg => cgl%cg
                call cg%cleanup_flux()
                call solve_cg_unsplit(cg,istep)
                cgl => cgl%nxt
            end do

            call update_boundaries(istep)

        end do
    end subroutine unsplit_sweep   

end module unsplit_sweeps