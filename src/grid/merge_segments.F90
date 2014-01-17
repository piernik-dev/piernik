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

!> \brief Module that merges individual segments for internal boundaries into large clumps to optimize MPI communication on AMR grids

module merge_segments

   use sort_segment_list, only: sort_segment_list_T

   implicit none

   private
   public :: merge_segments_T

   type :: merge_segments_T
      type(sort_segment_list_T), dimension(:), allocatable :: sl ! segment list (FIRST:LAST)
      logical :: valid
   contains
      procedure :: merge ! Merge segments
   end type merge_segments_T

contains

!> \brief Merge segments

   subroutine merge(this, list)

      use cg_list, only: cg_list_T

      implicit none

      class(merge_segments_T), intent(inout) :: this
      class(cg_list_T),        intent(in)    :: list

      !this%valid = .true.

   end subroutine merge

end module merge_segments
