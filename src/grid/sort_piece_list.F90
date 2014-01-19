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

!> \brief Module that contains a definition of sortable list of grid pieces, together with definition required by abstract type sortable_list_T

module sort_piece_list

   use constants, only: ndims
   use sortable_list, only: sortable_list_T

   implicit none

   private
   public :: grid_piece_list

   type :: grid_piece
      integer(kind=8)                   :: id        !< unique number used to sort
      integer(kind=8), dimension(ndims) :: off       !< offset
      integer(kind=4), dimension(ndims) :: n_b       !< size
      integer(kind=4)                   :: cur_gid   !< current grid_id
      integer(kind=4)                   :: cur_proc  !< current process number
      integer(kind=4)                   :: dest_proc !< process number according to ideal ordering
      real                              :: weight    !< number of cells relative to total number of cells on given level (unused)
      real                              :: cweight   !< cumulative weight for id <= own id (unused)
   contains
      procedure :: set_gp                            !< Set the primary properties, initialize derived properties with safe defaults
   end type grid_piece

   type, extends(sortable_list_T) :: grid_piece_list
      type(grid_piece), dimension(:), allocatable :: list !< the list itself
      type(grid_piece) :: temp
   contains
      ! override abstract interface routines
      procedure :: l_bound          !< Get lower bound of the list
      procedure :: u_bound          !< Get upper bound of the list
      procedure :: assign_element   !< Make an assignment
      procedure :: compare_elements !< Make a comparision

      ! own routines
      procedure :: init             !< Allocate the list
      procedure :: cleanup          !< Deallocate the list
      procedure :: set_id           !< Find grid id using a space-filling curve
      procedure :: set_weights      !< Find estimates of the cost of the grids
   end type grid_piece_list

contains

!> \brief initialize an element of the list to be sorted

   subroutine set_gp(this, off, n_b, gid, proc)

      use constants, only: ndims, INVALID

      implicit none

      class(grid_piece),                 intent(inout) :: this
      integer(kind=8), dimension(ndims), intent(in)    :: off   !< offset
      integer(kind=4), dimension(ndims), intent(in)    :: n_b   !< size
      integer(kind=4),                   intent(in)    :: gid   !< current grid_id
      integer(kind=4),                   intent(in)    :: proc  !< current process number

      this%off       = off
      this%n_b       = n_b
      this%cur_gid   = gid
      this%cur_proc  = proc
      this%dest_proc = INVALID
      this%id        = INVALID
      this%weight    = 0.
      this%cweight   = 0.

   end subroutine set_gp

!> \brief Allocate the list

   subroutine init(this, size)

      implicit none

      class(grid_piece_list), intent(inout) :: this
      integer(kind=4),        intent(in)    :: size

      allocate(this%list(size))

   end subroutine init

!> \brief Deallocate the list

   subroutine cleanup(this)

      implicit none

      class(grid_piece_list), intent(inout) :: this

      if (allocated(this%list)) deallocate(this%list)

   end subroutine cleanup

!> \brief Find grid id using a space-filling curve.

   subroutine set_id(this, off)

      use ordering, only: SFC_order

      implicit none

      class(grid_piece_list),            intent(inout) :: this
      integer(kind=8), dimension(ndims), intent(in)    :: off  !< offset of the level

      integer :: s

      do s = lbound(this%list, dim=1), ubound(this%list, dim=1)
         !> \todo use Hilbert ordering
         this%list(s)%id = SFC_order(this%list(s)%off-off)
      enddo

   end subroutine set_id

!> \brief Find estimates of the cost of the grids
!> \todo Do we want to include count of cg%leafmap, when available?

   subroutine set_weights(this)

      implicit none

      class(grid_piece_list), intent(inout) :: this

      integer :: s
      integer(kind=8) :: cc, c

      cc = 0
      do s = lbound(this%list, dim=1), ubound(this%list, dim=1)
         cc = cc + product(this%list(s)%n_b)
      enddo
      c = 0
      do s = lbound(this%list, dim=1), ubound(this%list, dim=1)
         c = c + product(this%list(s)%n_b)
         this%list(s)%weight = product(this%list(s)%n_b)/real(cc)
         this%list(s)%cweight = c/real(cc)
      enddo

   end subroutine set_weights

!>
!! \brief Tell if element at position a is greater than element at position b.
!! When the position equals temp_index, use temporary storage.
!<

   logical function compare_elements(this, a, b)

      use sortable_list, only: temp_index

      implicit none

      class(grid_piece_list), intent(inout) :: this
      integer,                intent(in)    :: a, b

      if (a == b) then
         compare_elements = .false.
         return
      endif

      if (b == temp_index) then ! this is the only case occuring in sortable_list::sort
         compare_elements = this%list(a)%id > this%temp%id
      else if (a == temp_index) then
         compare_elements = this%temp%id > this%list(b)%id
      else
         compare_elements = this%list(a)%id > this%list(b)%id
      endif

   end function compare_elements
!>
!! \brief Copy element at position b to position a.
!! When the position equals temp_index, use temporary storage.
!<

   subroutine assign_element(this, a, b)

      use sortable_list, only: temp_index

      implicit none

      class(grid_piece_list), intent(inout) :: this
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

      class(grid_piece_list), intent(in) :: this

      l_bound = lbound(this%list, dim=1)

   end function l_bound

!> \brief Get upper bound of the element list

   integer function u_bound(this)

      implicit none

      class(grid_piece_list), intent(in) :: this

      u_bound = ubound(this%list, dim=1)

   end function u_bound

end module sort_piece_list
