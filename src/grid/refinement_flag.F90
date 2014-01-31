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

!> \brief This module contains variables and initialization routines related to refinement

module refinement_flag

   implicit none

   private
   public :: ref_flag, level_min, level_max

   ! A candidate for refinement
   type :: SFC_candidate
      integer(kind=8) :: level  ! level at which we want to put grid block
      integer(kind=8) :: SFC_id ! position at which we want to put grid block
   end type SFC_candidate

   type :: ref_flag
      logical :: refine   !> a request to refine
      logical :: derefine !> a request to derefine
      type(SFC_candidate), allocatable, dimension(:) :: SFC_refine_list
   contains
      procedure :: init      !> Initialize: (.false. , .false., allocate 0 elements)
      procedure :: add       !> Appends one element to SFC_refine_list
      procedure :: sanitize  !> Sanitize the refinement flags
   end type ref_flag

   integer(kind=4) :: level_min !< minimum allowed refinement
   integer(kind=4) :: level_max !< maximum allowed refinement (don't need to be reached if not necessary)

contains

!> \brief Initialize to (.false. , .false., allocate 0 elements)

   subroutine init(this)

      implicit none

      class(ref_flag), intent(inout) :: this     ! object invoking this procedure

      this%refine = .false.
      this%derefine = .false.
      if (allocated(this%SFC_refine_list)) deallocate(this%SFC_refine_list)
      allocate(this%SFC_refine_list(0))

   end subroutine init

!> \brief Sanitize the refinement flags

   subroutine sanitize(this, my_level)

      implicit none

      class(ref_flag), intent(inout) :: this     ! object invoking this procedure
      integer(kind=4), intent(in)    :: my_level ! refinement level at which the flag has to be sanitized

      if (size(this%SFC_refine_list) > 0) this%refine = .true.

      if (my_level >= level_max) this%refine   = .false.
      if (my_level <  level_min) this%refine   = .true.

      if (this%refine) this%derefine = .false.

      if (my_level >  level_max) this%derefine = .true.
      if (my_level <= level_min) this%derefine = .false.

   end subroutine sanitize

!> \brief Appends one element to SFC_refine_list

   subroutine add(this, level, SFC_id)

      implicit none

      class(ref_flag), intent(inout) :: this   ! object invoking this procedure
      integer(kind=4), intent(in)    :: level  ! level at which we want to put grid block
      integer(kind=8), intent(in)    :: SFC_id ! position at which we want to put grid block

      this%SFC_refine_list = [ this%SFC_refine_list, SFC_candidate(int(level, kind=8), SFC_id) ] !lhs reallocation

   end subroutine add

end module refinement_flag
