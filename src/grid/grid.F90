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
!! \brief Initializes the crucial grid lists and provides a cleanup routine.
!!
!! \todo For domains with non-0 offsets (e.g. crawling, expanding) provide some means to communicate offset requirements for multigrid
!<

module grid

   implicit none

   private
   public :: init_grid, cleanup_grid

contains

!> \brief Routine that prepares base level and most important cg lists

   subroutine init_grid

      use cg_leaves,          only: leaves
      use cg_level_base,      only: base
      use cg_level_coarsest,  only: coarsest
      use cg_level_finest,    only: finest
      use constants,          only: PIERNIK_INIT_DOMAIN, tmr_amr, V_DEBUG
      use dataio_pub,         only: printinfo, die, code_progress
      use domain,             only: dom
      use mpisetup,           only: master
      use timer,              only: set_timer

      implicit none

      real :: ts  !< time for runtime profiling

      ts =  set_timer(tmr_amr, .true.)  ! we need it here for call leaves%update()
      if (code_progress < PIERNIK_INIT_DOMAIN) call die("[grid:init_grid] domain not initialized.")

      if (master) call printinfo("[grid:init_grid]: commencing...", V_DEBUG)

      ! Create the empty main lists.with the base level

      allocate(base)
      call base%set(dom%n_d)
      finest%level => base%level
      coarsest%level => base%level
      if (master) call base%level%add_patch
      call base%level%init_all_new_cg

      ! Refinement lists will be added by iterating the initproblem::problem_initial_conditions routine, in restart_hdf5::read_restart_hdf5 or in refinement_update
      ! Underground levels will be added in multigrid::init_multigrid

      call leaves%update(" (base level) ")

      if (master) call printinfo("[grid:init_grid]: cg finished. \o/", V_DEBUG)

   end subroutine init_grid

!> \brief Deallocate everything

   subroutine cleanup_grid

      use cg_leaves,          only: leaves
      use cg_level_base,      only: base
      use cg_level_coarsest,  only: coarsest
      use cg_level_connected, only: cg_level_connected_t
      use cg_list_global,     only: all_cg
      use grid_cont_fcflx,    only: cleanup_flxp
      use list_of_cg_lists,   only: all_lists
      use named_array_list,   only: qna, wna

      implicit none

      type(cg_level_connected_t), pointer :: curl, aux

      call cleanup_flxp

      if (allocated(leaves%up_to_level)) deallocate(leaves%up_to_level)

      curl => coarsest%level
      do while (associated(curl))
         call curl%cleanup
         curl => curl%finer
      enddo

      call all_cg%delete_all
      call all_lists%delete

      curl => coarsest%level
      do while (associated(curl))
         aux => curl
         curl => curl%finer
         deallocate(aux)
      enddo
      deallocate(base)

      if (allocated(qna%lst)) deallocate(qna%lst)
      if (allocated(wna%lst)) deallocate(wna%lst)

   end subroutine cleanup_grid

end module grid
