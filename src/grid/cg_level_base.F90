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

!> \brief This module implements the base level and related methods

module cg_level_base

   use cg_level_connected, only: cg_level_connected_t

   implicit none

   private
   public :: base

   abstract interface
      subroutine no_args
         implicit none
      end subroutine no_args
   end interface

   !! \brief The pointer of the base level and a method to initialize it
   type :: cg_level_base_t
      type(cg_level_connected_t), pointer :: level            !< The base level
      procedure(no_args), nopass, pointer :: init_multigrid   !< a pointer to multigrid:init_multigrid or null()
   contains
      procedure          :: set                               !< initialize the base level
   end type cg_level_base_t

   type(cg_level_base_t), pointer :: base                     !< base level grid containers

contains

!> \brief Initialize the base level

   subroutine set(this, n_d)

      use constants,          only: base_level_id, ndims
      use dataio_pub,         only: die
      use domain,             only: dom
      use list_of_cg_lists,   only: all_lists

      implicit none

      class(cg_level_base_t),            intent(inout) :: this   !< object invoking type bound procedure
      integer(kind=4), dimension(ndims), intent(in)    :: n_d    !< size of global base grid in cells

      ! Multigrid and refinement work properly with non-0, even offset.
      ! Offset value equal to k*2**n, where k is odd will allow at most n levels of coarsening.
      ! Odd offsets or domain sizes prevent creation of coarse levels.
      !> \todo Find the limit that comes from multigrid: maximum refinement should not depend on base level offset
      ! Offset of the base domain may change after the domain gets expanded, shrinked or resized.

      if (any(n_d(:) < 1)) call die("[cg_level_base:set] non-positive base grid sizes")
      if (any(dom%has_dir(:) .neqv. (n_d(:) > 1))) call die("[cg_level_base:set] base grid size incompatible with has_dir masks")

      this%init_multigrid => null()  ! it will be set by init_multigrid, if it is included and called
      allocate(this%level)
      call this%level%init_level

      call this%level%l%init(base_level_id, int(n_d, kind=8), dom%off)

      call all_lists%register(this%level, "Base level")

   end subroutine set

end module cg_level_base
