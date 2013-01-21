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

!> \brief This module implements the base level and related methods

module cg_level_base

   use cg_level_connected, only: cg_level_connected_T

   implicit none

   private
   public :: base

   !! \brief The pointer of the base level and a method to initialize it
   !> \todo Domainshrinking, expanding and crawling should also be implemented here
   type, extends(cg_level_connected_T) :: cg_level_base_T
      type(cg_level_connected_T), pointer :: level            !< The base level
   contains
      procedure :: set                                        !< initialize the base level
   end type cg_level_base_T

   type(cg_level_base_T), pointer  :: base                    !< base level grid containers

contains

!> \brief Initialize the base level

   subroutine set(this, n_d, offset)

      use constants,        only: base_level_id, ndims
      use dataio_pub,       only: die
      use domain,           only: dom
      use list_of_cg_lists, only: all_lists

      implicit none

      class(cg_level_base_T),            intent(inout) :: this   !< object invoking type bound procedure
      integer(kind=4), dimension(ndims), intent(in)    :: n_d    !< size of global base grid in cells
      integer(kind=8), dimension(ndims), intent(in)    :: offset !< offset of global base grid in cells

      if (any(n_d(:) < 1)) call die("[cg_level_connected:set] non-positive base grid sizes")
      if (any(dom%has_dir(:) .neqv. (n_d(:) > 1))) call die("[cg_level_connected:set] base grid size incompatible with has_dir masks")

      allocate(this%level)
      call this%level%init_level
      this%level%level_id = base_level_id
      
      where (dom%has_dir(:))
         this%level%n_d(:) = n_d(:)
         this%level%off(:) = offset(:)
      elsewhere
         this%level%n_d(:) = 1
         this%level%off(:) = 0   
      endwhere
      call all_lists%register(this%level, "Base level")

   end subroutine set

end module cg_level_base
