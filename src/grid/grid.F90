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
#include "macros.h"
!>
!! \brief Module containing lists of grid containers for computational mesh and initialization and cleanup routines
!<
module grid

   implicit none

   private
   public :: init_grid, cleanup_grid

contains

!> \brief Routine that prepares base level and most importand cg lists

   subroutine init_grid

      use cg_list_global, only: all_cg
      use cg_list_level,  only: base_lev
      use constants,      only: PIERNIK_INIT_DOMAIN, I_ZERO, base_level_offset
      use dataio_pub,     only: printinfo, die, code_progress
      use domain,         only: dom

      implicit none

      if (code_progress < PIERNIK_INIT_DOMAIN) call die("[grid:init_grid] domain not initialized.")

#ifdef VERBOSE
      call printinfo("[grid:init_grid]: commencing...")
#endif /* VERBOSE */

      ! Create the empty main lists for base level only.
      ! Refinement lists will be added by iterating the initproblem::init_prob routine, in restart_hdf5::read_restart_hdf5 or in not_yet_implemented::refinement_update
      ! Underground levels will be added in multigrid::init_multigrid
      call all_cg%init
      all_cg%ord_prolong_nb = I_ZERO
      call base_lev%add_level(dom%n_d(:))

      call base_lev%add_patch(dom%n_d(:), base_level_offset)

      call base_lev%init_all_new_cg

      call all_cg%update_leaves

#ifdef VERBOSE
      call printinfo("[grid:init_grid]: cg finished. \o/")
#endif /* VERBOSE */

   end subroutine init_grid

! \brief Update the list of leaves
! subroutine update leaves
! end subroutine update leaves

!> \brief deallocate everything
   subroutine cleanup_grid

      use cg_list,        only: cg_list_element
      use cg_list_bnd,    only: leaves
      use cg_list_global, only: all_cg
      use cg_list_level,  only: cg_list_level_T, coarsest, finest, base_lev
      use named_array_list, only: qna, wna

      implicit none

      type(cg_list_element), pointer :: cgle
      type(cg_list_level_T), pointer :: curl
      integer :: p

      call leaves%clear
      curl => coarsest
      do while (associated(curl))
         call curl%delete
         if (allocated(curl%pse)) deallocate(curl%pse) ! curl%pse(:)%c should be deallocated automagically
         do p = lbound(curl%patches, dim=1), ubound(curl%patches, dim=1)
            deallocate(curl%patches(p)%pse) ! curl%patches(p)%pse(:)%c should be deallocated automagically
         enddo
         deallocate(curl%patches)

         curl => curl%finer
         if (associated(curl)) then
            if (associated(curl%coarser) .and. .not. associated(curl%coarser, base_lev)) deallocate(curl%coarser)
         endif
      enddo
      if (.not. associated(finest, base_lev)) deallocate(finest)

      ! manually deallocate all grid containers first
      cgle => all_cg%first
      do while (associated(cgle))
         deallocate(cgle%cg)
         cgle => cgle%nxt
      enddo
      if (allocated(qna%lst)) deallocate(qna%lst)
      if (allocated(wna%lst)) deallocate(wna%lst)
      call all_cg%delete

   end subroutine cleanup_grid

end module grid
