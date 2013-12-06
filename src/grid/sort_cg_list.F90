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

!>
!! \brief Stub module
!!
!! \details Don't forget to describe your stuff
!<

module sort_cg_list

   use cg_list,       only: cg_list_element
   use sortable_list, only: sortable_list_T

   implicit none

   private
   public :: sort_cg_list_T

   ! prepare pointer type to make an array of pointers
   type :: cg_piece
      type(cg_list_element), pointer :: cgl
   end type cg_piece

   type, extends(sortable_list_T) :: sort_cg_list_T
      type(cg_piece), dimension(:), allocatable :: list !< the list itself
      type(cg_piece) :: temp
   contains
      ! override abstract interface routines
      procedure :: init             !< Allocate the list
      procedure :: cleanup          !< Deallocate the list
      procedure :: l_bound          !< Get lower bound of the list
      procedure :: u_bound          !< Get upper bound of the list
      procedure :: assign_element   !< Make an assignment
      procedure :: compare_elements !< Make a comparision
   end type sort_cg_list_T

contains

!> \brief Allocate the list

   subroutine init(this, size)

      implicit none

      class(sort_cg_list_T), intent(inout) :: this
      integer(kind=4),        intent(in)    :: size

      allocate(this%list(size))

   end subroutine init

!> \brief Deallocate the list

   subroutine cleanup(this)

      implicit none

      class(sort_cg_list_T), intent(inout) :: this

      if (allocated(this%list)) deallocate(this%list)

   end subroutine cleanup

!>
!! \brief Tell if element at position a is greater than element at position b.
!! When the position equals temp_index, use temporary storage.
!<

   logical function compare_elements(this, a, b)

      use sortable_list, only: temp_index

      implicit none

      class(sort_cg_list_T), intent(inout) :: this
      integer,                intent(in)    :: a, b

      if (a == b) then
         compare_elements = .false.
         return
      endif

      if (b == temp_index) then ! this is the only case occuring in sortable_list::sort
         compare_elements = this%list(a)%cgl%cg%SFC_id > this%temp%cgl%cg%SFC_id
      else if (a == temp_index) then
         compare_elements = this%temp%cgl%cg%SFC_id > this%list(b)%cgl%cg%SFC_id
      else
         compare_elements = this%list(a)%cgl%cg%SFC_id > this%list(b)%cgl%cg%SFC_id
      endif

   end function compare_elements
!>
!! \brief Copy element at position b to position a.
!! When the position equals temp_index, use temporary storage.
!<

   subroutine assign_element(this, a, b)

      use sortable_list, only: temp_index

      implicit none

      class(sort_cg_list_T), intent(inout) :: this
      integer,                intent(in)    :: a, b

      if (a == b) return

      if (a == temp_index) then
         this%temp = this%list(b)
      else if (b == temp_index) then
         this%list(a) = this%temp
      else
         this%list(a) = this%list(b)
      endif

   end subroutine assign_element

!> \brief Get lower bound of the element list

   integer function l_bound(this)

      implicit none

      class(sort_cg_list_T), intent(in) :: this

      l_bound = lbound(this%list, dim=1)

   end function l_bound

!> \brief Get upper bound of the element list

   integer function u_bound(this)

      implicit none

      class(sort_cg_list_T), intent(in) :: this

      u_bound = ubound(this%list, dim=1)

   end function u_bound

end module sort_cg_list
