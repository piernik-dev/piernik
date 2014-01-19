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
!! \brief Module that contains a definition of sortable list of segments for internal boundary exchange
!!
!! \details When there are multiple segments to be exchanged between the same pair of processes (which occurs
!! quite often in AMR), we can merge these segments into one big chunk of data and make one large MPI transaction
!! instead of many small MPI ones. We can sort the segments on the sender and the recipient according to MPI tag
!! to avoid passing auxiliary description on how the actual chunk was composed.
!<

module sort_segment_list

   use constants,     only: xdim, zdim, LO, HI
   use grid_cont,     only: grid_container
   use sortable_list, only: sortable_list_T

   implicit none

   private
   public :: sort_segment_list_T

   ! prepare type to gather information on segments to be exchanged with one process
   type :: seg
      integer(kind=4) :: tag                              !< unique tag for data exchange, used for sorting
      integer(kind=8), dimension(xdim:zdim, LO:HI) :: se  !< range
      integer(kind=8) :: offset                           !< offset within the chunk
      type(grid_container), pointer :: cg                 !< source/target grid container
      integer :: dir                                      !< direction of the boundary segment
   end type seg

   type, extends(sortable_list_T) :: sort_segment_list_T
      type(seg), dimension(:), allocatable :: list !< the list itself
      type(seg) :: temp
      integer(kind=8) :: total_size                !< offset of the last segment data + its size
      real, allocatable, dimension(:) :: buf       !< buffer for the data to be sent or received (total_size)
   contains
      procedure :: add              !< Add a piece to the list
      procedure :: cleanup          !< Deallocate the list
      procedure :: find_offsets     !< Set up offsets

      ! override abstract interface routines
      procedure :: l_bound          !< Get lower bound of the list
      procedure :: u_bound          !< Get upper bound of the list
      procedure :: assign_element   !< Make an assignment
      procedure :: compare_elements !< Make a comparision
   end type sort_segment_list_T

contains

!> \brief Allocate the list

   subroutine add(this, tag, se, cg, dir)

      use grid_cont, only: grid_container

      implicit none

      class(sort_segment_list_T),                   intent(inout) :: this
      integer(kind=4),                              intent(in)    :: tag
      integer(kind=8), dimension(xdim:zdim, LO:HI), intent(in)    :: se
      type(grid_container), pointer,                intent(in)    :: cg
      integer,                                      intent(in)    :: dir

      type(seg), dimension(:), allocatable :: tmp

      !> \warning a lot of duplicated code with some othe rroutines such as named_array_list::add2lst
      if (.not. allocated(this%list)) then
         allocate(this%list(1))
      else
         allocate(tmp(lbound(this%list(:),dim=1):ubound(this%list(:), dim=1) + 1))
         tmp(:ubound(this%list(:), dim=1)) = this%list(:)
         call move_alloc(from=tmp, to=this%list)
      endif
      this%list(ubound(this%list(:), dim=1))%tag =  tag
      this%list(ubound(this%list(:), dim=1))%se  =  se
      this%list(ubound(this%list(:), dim=1))%cg  => cg
      this%list(ubound(this%list(:), dim=1))%dir =  dir

   end subroutine add

!< \brief Set up offsets

   subroutine find_offsets(this, dmask)

      use constants, only: xdim, cor_dim, LO, HI, I_ONE

      implicit none

      class(sort_segment_list_T),       intent(inout) :: this
      logical, dimension(xdim:cor_dim), intent(in)    :: dmask !< .true. for the directions we want to include

      integer :: i
      integer(kind=8) :: off

      if (allocated(this%list)) then
         off = 0
         do i = lbound(this%list, dim=1), ubound(this%list, dim=1)
            this%list(i)%offset = off
            if (dmask(this%list(i)%dir)) &
                 off = off + product(this%list(i)%se(:, HI) - this%list(i)%se(:, LO) + I_ONE)
         enddo
         this%total_size = off
      else
         this%total_size = 0
      endif

   end subroutine find_offsets

!> \brief Deallocate the list

   subroutine cleanup(this)

      implicit none

      class(sort_segment_list_T), intent(inout) :: this

      if (allocated(this%list)) deallocate(this%list)

   end subroutine cleanup

!>
!! \brief Tell if element at position a is greater than element at position b.
!! When the position equals temp_index, use temporary storage.
!<

   logical function compare_elements(this, a, b)

      use sortable_list, only: temp_index

      implicit none

      class(sort_segment_list_T), intent(inout) :: this
      integer,                    intent(in)    :: a, b

      if (a == b) then
         compare_elements = .false.
         return
      endif

      if (b == temp_index) then ! this is the only case occuring in sortable_list::sort
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

      class(sort_segment_list_T), intent(inout) :: this
      integer,                    intent(in)    :: a, b

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

      class(sort_segment_list_T), intent(in) :: this

      if (allocated(this%list)) then
         l_bound = lbound(this%list, dim=1)
      else
         l_bound = huge(1) ! safe default; should bail out from sorting
      endif

   end function l_bound

!> \brief Get upper bound of the element list

   integer function u_bound(this)

      implicit none

      class(sort_segment_list_T), intent(in) :: this

      if (allocated(this%list)) then
         u_bound = ubound(this%list, dim=1)
      else
         u_bound = -huge(1) ! safe default; should bail out from sorting
      endif

   end function u_bound

end module sort_segment_list
