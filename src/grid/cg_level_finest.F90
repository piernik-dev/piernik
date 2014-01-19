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

   use cg_level_connected, only: cg_level_connected_T

   implicit none

   private
   public :: finest

   !> \brief The pointer of the finest refinement level and a method to add a finer one
   type :: cg_level_finest_T
      type(cg_level_connected_T), pointer :: level !< highest refinement level
    contains
      procedure          :: add_finer
      procedure          :: equalize
      !> \todo Provide delete_finest and use it in cleanup
   end type cg_level_finest_T

   type(cg_level_finest_T) :: finest               !< finest level of refinement

contains

!>
!! \brief It is not allowed for different processes to have different height of level hierarchy, so let's add some empty levels, where necessary
!!
!! \details Multigrid levels are always created on all processes, so no need to fix anything on the bottom of the hierarchy
!!
!! \todo When global top level do not have any blocks then destroy it.
!<

   subroutine equalize(this)

      use constants,  only: pMAX
      use dataio_pub, only: die
      use mpisetup,   only: piernik_MPI_Allreduce

      implicit none

      class(cg_level_finest_T), intent(inout) :: this    !< object calling type-bound routine

      integer(kind=4) :: g_finest_id

      g_finest_id = this%level%level_id
      call piernik_MPI_Allreduce(g_finest_id, pMAX)

      do while (g_finest_id > this%level%level_id)
         call this%add_finer
      enddo

      g_finest_id = this%level%level_id
      call piernik_MPI_Allreduce(g_finest_id, pMAX)

      if (g_finest_id /= this%level%level_id) call die("[cg_level_finest:equalize] failure")

   end subroutine equalize

!> \brief Add a fine level to the top of existing hierarchy

   subroutine add_finer(this)

      use constants,        only: INVALID, I_ONE, refinement_factor
      use dataio_pub,       only: die, msg
      use domain,           only: dom
      use grid_helpers,     only: c2f_o
      use list_of_cg_lists, only: all_lists

      implicit none

      class(cg_level_finest_T), intent(inout) :: this    !< object calling type-bound routine

      type(cg_level_connected_T), pointer     :: new_lev !< fresh refinement level to be added

      allocate(new_lev)
      call new_lev%init_level
      new_lev%n_d(:) = 1

      if (associated(this%level%finer)) call die("[cg_level_finest:add_finer] finer level already exists")

      new_lev%level_id = this%level%level_id + I_ONE
      new_lev%off = c2f_o(this%level%off)
      where (dom%has_dir(:)) new_lev%n_d(:) = this%level%n_d(:) * refinement_factor

      !! make sure that vertical_prep will be called where necessary
      this%level%ord_prolong_set = INVALID
      write(msg, '(a,i3)')"level ",new_lev%level_id
      call all_lists%register(new_lev, msg)

      this%level%finer => new_lev
      new_lev%coarser => this%level
      this%level => new_lev

   end subroutine add_finer

end module cg_level_finest
