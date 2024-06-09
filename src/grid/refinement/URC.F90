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
!! \brief Unified refinement criteria
!!
!! \details This module contains the abstract parts of the structures and some common stuff
!<

module unified_ref_crit

   use constants, only: INVALID

   implicit none

   private
   public :: urc

!>
!! \brief Things that should be common for all refinement criteria.
!!
!! \details The refinement criteria are meant to be ordered into a unidirectional list.
!! For each grid container the routines for marking refinements are called in the same order.
!! It is assumed that all vital fields will have valid boundaries for the routine mark.
!! Each criterion should be applicable on any cg and no global communication should take place in this routine.
!! Refinement criteria that require some global operations (like Toomre length in disks)
!! should be implemented separately and used with care.
!<

   type, abstract :: urc
      integer             :: iplot = INVALID  !< field index for storing the refinement criterion value
      class(urc), pointer :: next  => null()  !< next refinement criterion or null() (for unidirectional list)
   contains
      procedure(mark_urc), deferred :: mark   !< a routine that takes a cg and leaves suggestions on refining
   end type urc

   interface

!>
!! \brief Mark refinements on given grid container
!!
!! intent(inout) :: this allows for altering private variables (implicit initialisation)
!<

      subroutine mark_urc(this, cg)

         use grid_cont, only: grid_container
         import urc

         implicit none

         class(urc),                    intent(inout) :: this  !< an object invoking the type-bound procedure
         type(grid_container), pointer, intent(inout) :: cg    !< current grid piece

      end subroutine mark_urc

   end interface

end module unified_ref_crit
