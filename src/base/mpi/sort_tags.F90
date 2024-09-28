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

!> \brief Module that contains a definition of sortable list of MPI tags

module sort_tags

   use sortable_list, only: sortable_list_t

   implicit none

   private
   public :: sort_tags_t, req_data

   ! basic info for requests stored here
   type :: req_data
      integer(kind=4) :: tag  !< the tag
      integer(kind=4) :: ir   !< request index
   end type req_data

   type, extends(sortable_list_t) :: sort_tags_t
      type(req_data), dimension(:), allocatable :: list  !< the list itself
      type(req_data)                            :: temp  !< the buffer for sorting
      integer                                   :: cnt   !< the number of stored elements
   contains
      procedure :: cleanup          !< Deallocate the list
      procedure :: dump             !< append an element to the list

      ! override abstract interface routines
      procedure :: l_bound          !< Get lower bound of the list
      procedure :: u_bound          !< Get upper bound of the list
      procedure :: assign_element   !< Make an assignment
      procedure :: compare_elements !< Make a comparison
   end type sort_tags_t

contains

!> \brief Append an element to the list. Allocate list first or expand it when necessary.

   subroutine dump(this, el)

      use dataio_pub, only: die

      implicit none

      class(sort_tags_t), intent(inout) :: this
      type(req_data),     intent(in)    :: el

      integer, parameter :: starting_size = 16
      type(req_data), dimension(:), allocatable :: new_list

      if (.not. allocated(this%list)) then  ! initialize
         allocate(this%list(starting_size))
         this%cnt = 0
      else if (this%cnt == ubound(this%list, dim=1)) then  !double the size
         allocate(new_list(2*this%cnt))
         new_list(1:this%cnt) = this%list(:)
         call move_alloc(from=new_list, to=this%list)
      endif

      if (.not. allocated(this%list)) call die("[sort_tags:dump] not allocated")  ! should never happen

      this%cnt = this%cnt + 1
      if (this%cnt < this%l_bound() .or. this%cnt > ubound(this%list, dim=1)) call die("[sort_tags:dump] run out of space")  ! perhaps can happen after undetected allocation error
      this%list(this%cnt) = el

   end subroutine dump

!> \brief Deallocate the list

   subroutine cleanup(this)

      use constants, only: INVALID

      implicit none

      class(sort_tags_t), intent(inout) :: this

      if (allocated(this%list)) deallocate(this%list)
      this%cnt = INVALID

   end subroutine cleanup

!>
!! \brief Tell if element at position a is greater than element at position b.
!! When the position equals temp_index, use temporary storage.
!<

   pure logical function compare_elements(this, a, b)

      use sortable_list, only: temp_index

      implicit none

      class(sort_tags_t), intent(in) :: this
      integer,               intent(in) :: a, b

      if (a == b) then
         compare_elements = .false.
         return
      endif

      if (b == temp_index) then ! this is the only case occurring in sortable_list::sort
         compare_elements = this%list(a)%tag > this%temp%tag
      else if (a == temp_index) then
         compare_elements = this%temp%tag > this%list(b)%tag
      else
         compare_elements = this%list(a)%tag > this%list(b)%tag
      endif

   end function compare_elements

!>
!! \brief Copy element at position b to position a.
!! When the position equals temp_index, use temporary storage.
!<

   subroutine assign_element(this, a, b)

      use sortable_list, only: temp_index

      implicit none

      class(sort_tags_t), intent(inout) :: this
      integer,            intent(in)    :: a, b

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

      class(sort_tags_t), intent(in) :: this

      l_bound = lbound(this%list, dim=1)

   end function l_bound

!> \brief Get upper bound of the element list

   integer function u_bound(this)

      implicit none

      class(sort_tags_t), intent(in) :: this

      u_bound = this%cnt

   end function u_bound

end module sort_tags
