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

!> \brief This module manages finest levels of refinement

module cg_level_finest

   use cg_level_connected, only: cg_level_connected_t

   implicit none

   private
   public :: finest

   !> \brief The pointer of the finest refinement level and a method to add a finer one
   type :: cg_level_finest_t
      type(cg_level_connected_t), pointer :: level  !< highest refinement level
   contains
      procedure :: add_finer        !< Add a fine level to the top of existing hierarchy
      procedure :: equalize         !< Add an empty fine level when any of the other threads has one
      procedure :: find_finest_bnd  !< Find finest level with external boundary
      !> \todo Provide delete_finest and use it in cleanup
   end type cg_level_finest_t

   type(cg_level_finest_t) :: finest               !< finest level of refinement

contains

!>
!! \brief It is not allowed for different processes to have different height of level hierarchy, so let's add some empty levels, where necessary
!!
!! \details Multigrid levels are always created on all processes, so no need to fix anything on the bottom of the hierarchy
!!
!! \todo When global top level does not have any blocks then destroy it.
!<

   subroutine equalize(this)

      use allreduce,  only: piernik_MPI_Allreduce
      use constants,  only: pMAX
      use dataio_pub, only: die

      implicit none

      class(cg_level_finest_t), intent(inout) :: this    !< object calling type-bound routine

      integer(kind=4) :: g_finest_id

      g_finest_id = this%level%l%id
      call piernik_MPI_Allreduce(g_finest_id, pMAX)

      do while (g_finest_id > this%level%l%id)
         call this%add_finer
      enddo

      g_finest_id = this%level%l%id
      call piernik_MPI_Allreduce(g_finest_id, pMAX)

      if (g_finest_id /= this%level%l%id) call die("[cg_level_finest:equalize] failure")

   end subroutine equalize

!> \brief Add a fine level to the top of existing hierarchy

   subroutine add_finer(this)

      use constants,        only: INVALID, I_ONE, refinement_factor
      use dataio_pub,       only: die, msg
      use grid_helpers,     only: c2f_o
      use list_of_cg_lists, only: all_lists

      implicit none

      class(cg_level_finest_t), intent(inout) :: this    !< object calling type-bound routine

      type(cg_level_connected_t), pointer     :: new_lev !< fresh refinement level to be added

      allocate(new_lev)
      call new_lev%init_level

      if (associated(this%level%finer)) call die("[cg_level_finest:add_finer] finer level already exists")

      call new_lev%l%init(this%level%l%id+I_ONE, this%level%l%n_d*refinement_factor, c2f_o(this%level%l%off))

      !! make sure that vertical_prep will be called where necessary
      this%level%ord_prolong_set = INVALID
      write(msg, '(a,i3)')"level ",new_lev%l%id
      call all_lists%register(new_lev, msg)

      this%level%finer => new_lev
      new_lev%coarser => this%level
      this%level => new_lev

   end subroutine add_finer

!<
!! \brief Find finest level with external boundary
!!
!! Returns null() on periodic domains or pointer to highest level with any external boundary
!>

   function find_finest_bnd(this) result(level)

      use allreduce,  only: piernik_MPI_Allreduce
      use cg_list,    only: cg_list_element
      use constants,  only: pMAX
      use dataio_pub, only: die
      use domain,     only: dom

      implicit none

      class(cg_level_finest_t), intent(inout) :: this    !< object calling type-bound routine

      type(cg_level_connected_t), pointer :: level

      type(cg_level_connected_t), pointer :: l
      type(cg_list_element), pointer :: cgl
      integer(kind=4) :: il

      nullify(level)
      if (count(dom%periodic) == dom%eff_dim) return  ! fully periodic domain

      l => this%level
      do while (.not. associated(level) .and. associated(l))
         cgl => l%first
         do while (associated(cgl))
            if (any(cgl%cg%ext_bnd)) level => l
            cgl => cgl%nxt
         enddo
         l => l%coarser
      enddo

      if (associated(level)) then
         il = level%l%id
      else
         il = -1000  ! unrealistically small level
      endif
      call piernik_MPI_Allreduce(il, pMAX)

      level => this%level
      do while (level%l%id > il)
         level => level%coarser
         if (.not. associated(level)) call die("[cg_level_finest:find_finest_bnd] Failed to find level with finest boundary")
      enddo

   end function find_finest_bnd

end module cg_level_finest
