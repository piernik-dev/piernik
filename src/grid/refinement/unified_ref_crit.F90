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
!    along with PIERNIK.  If not, see http://www.gnu.org/licenses/.
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

   implicit none

   private
   public :: urc, first_ref_crit

!>
!! \brief Things that should be common for all refinement criteria.
!!
!! \details The refinement criteria are meant to be ordered into a unidirectional list.
!! For each grid container the routines for marking refinements and derefinements are called in the same order.
!! It is assumed that all vital fields will have valid boundaries for the routine mark.
!! Each criterion should be applicable on any cg and no global communication should take place in this routine.
!! Refinement criteria that require some global operations should be implemented separately and used with care.
!<

   type, abstract :: urc
      integer                    :: iplot  !< field index for storing the refinement criterion value or INVALID
      class(urc), pointer        :: next   !< next refinement ctiterion or null() (for unidirectional list)
   contains
      procedure(p_urc), deferred :: mark   !< a routine that takes a cg and leaves suggestions on refining and derefining
   end type urc

!> \brief Mark refinements and derefinements on given grid container

   interface

      subroutine p_urc(this, cg)

         use grid_cont, only: grid_container
         import urc

         implicit none

         class(urc),                    intent(in)    :: this  !< an object invoking the type-bound procedure
         type(grid_container), pointer, intent(inout) :: cg    !< current grid piece

      end subroutine p_urc

   end interface

   class(urc), pointer :: first_ref_crit => null()  !< here the list should start
   ! A pointer to last refinement criterion would be of some use only during initialisation.

end module unified_ref_crit
