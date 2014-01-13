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

!> \brief This module manages coarsest levels of refinement

module cg_level_coarsest

   use cg_level_connected, only: cg_level_connected_T

   implicit none

   private
   public :: coarsest

   !> \brief The pointer of the coarsest refinement level and a method to add a coarser one
   type :: cg_level_coarsest_T
      type(cg_level_connected_T), pointer :: level !< lowest refinement level
    contains
      procedure :: add_coarser                  !< add one level below current coarsest level
      procedure :: delete_coarsest              !< delete coarsest level
      procedure :: delete_coarser_than_base     !< delete all levels below base level (multigrid levels)
      !> \todo use delete_coarsest in cleanup
   end type cg_level_coarsest_T

   type(cg_level_coarsest_T) :: coarsest             !< coarsest level of refinement

contains

!> \brief Add a coarse level to the bottom of existing hierarchy

   subroutine add_coarser(this)

      use constants,        only: INVALID, I_ONE, refinement_factor
      use dataio_pub,       only: die, msg
      use domain,           only: dom
      use grid_helpers,     only: f2c_o, c2f_o
      use list_of_cg_lists, only: all_lists
      use mpisetup,         only: master

      implicit none

      class(cg_level_coarsest_T), intent(inout) :: this    !< object calling type-bound routine

      type(cg_level_connected_T), pointer       :: new_lev !< fresh refinement level to be added

      if (associated(this%level%coarser)) call die("[cg_level_coarsest:add_coarser] coarser level already exists")

      allocate(new_lev)
      call new_lev%init_level
      new_lev%n_d(:) = 1

      new_lev%level_id = this%level%level_id - I_ONE
      new_lev%off = f2c_o(this%level%off)
      if (any(c2f_o(new_lev%off) /= this%level%off)) then
         write(msg, '(a,3f10.1,a,i3)')"[cg_level_coarsest:add_coarser] Fractional offset: ", this%level%off(:)/real(refinement_factor), " at level ",new_lev%level_id
         call die(msg)
      endif
      where (dom%has_dir(:)) new_lev%n_d(:) = this%level%n_d(:) / refinement_factor
      if (master .and. any(new_lev%n_d(:)*refinement_factor /= this%level%n_d(:) .and. dom%has_dir(:))) then
         write(msg, '(a,3f10.1,a,i3)')"[cg_level_coarsest:add_coarser] Fractional number of domain cells: ", this%level%n_d(:)/real(refinement_factor), " at level ",new_lev%level_id
         call die(msg)
      endif

      !! make sure that vertical_prep will be called where necessary
      this%level%ord_prolong_set = INVALID
      write(msg, '(a,i3)')"level ",new_lev%level_id
      call all_lists%register(new_lev, msg)

      this%level%coarser => new_lev
      new_lev%finer => this%level
      this%level => new_lev

   end subroutine add_coarser

!> \brief delete coarsest level

   subroutine delete_coarsest(this)

      use cg_list,            only: cg_list_T
      use constants,          only: base_level_id
      use dataio_pub,         only: die
      use list_of_cg_lists,   only: all_lists

      implicit none

      class(cg_level_coarsest_T), intent(inout) :: this    !< object calling type-bound routine

      class(cg_list_T), pointer :: curl

      if (this%level%level_id >= base_level_id) call die("[cg_level_coarsest:delete_coarsest] Attempted to operate on base level or above")

      call this%level%free_all_cg
      curl => this%level
      call all_lists%unregister(curl)
      this%level => this%level%finer
      nullify(this%level%coarser)
      deallocate(curl)

   end subroutine delete_coarsest

!> \brief delete all levels below base level (multigrid levels)

   subroutine delete_coarser_than_base(this)

      use constants, only: base_level_id

      implicit none

      class(cg_level_coarsest_T), intent(inout) :: this    !< object calling type-bound routine

      do while (this%level%level_id < base_level_id)
         call this%delete_coarsest
      enddo

   end subroutine delete_coarser_than_base

end module cg_level_coarsest
