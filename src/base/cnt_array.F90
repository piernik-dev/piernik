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

!> \brief Type for array of counters

module cnt_array

   use constants, only: LONG

   implicit none

   private
   public :: arrcnt

   type :: arrcnt
      integer(LONG),                   allocatable, private :: c0
      integer(LONG), dimension(:),     allocatable, private :: c1
      integer(LONG), dimension(:,:),   allocatable, private :: c2
      integer(LONG), dimension(:,:,:), allocatable, private :: c3
   contains
      procedure :: init     !< allocate the array of counters
      procedure :: cleanup  !< deallocate everything
      procedure :: incr     !< increase selected counter by 1
      procedure :: add      !< add some value to the selected counter
      procedure :: print    !< print the collected values
   end type arrcnt

contains

!> \brief Allocate the array of counters

   subroutine init(this, sizes)

      use dataio_pub, only: die

      implicit none

      class(arrcnt),         intent(inout) :: this   !< an object invoking the type-bound procedure
      integer, dimension(:), intent(in)    :: sizes  !< array of sizes for this%c (1-based)

      if (allocated(this%c0) .or. allocated(this%c1) .or. allocated(this%c2) .or. allocated(this%c3)) &
           call die("[cnt_array:init] already initialized")

      select case (size(sizes))
         case (0)
            allocate(this%c0)
            this%c0 = 0
         case (1)
            allocate(this%c1(sizes(1)))
            this%c1 = 0
         case (2)
            allocate(this%c2(sizes(1), sizes(2)))
            this%c2 = 0
         case (3)
            allocate(this%c3(sizes(1), sizes(2), sizes(3)))
            this%c3 = 0
         case default
            call die("[cnt_array:init] rank not implemented")
      end select

   end subroutine init

!> \brief Deallocate everything

   subroutine cleanup(this)

      implicit none

      class(arrcnt), intent(inout) :: this  !< an object invoking the type-bound procedure

      if (allocated(this%c0)) deallocate(this%c0)
      if (allocated(this%c1)) deallocate(this%c1)
      if (allocated(this%c2)) deallocate(this%c2)
      if (allocated(this%c3)) deallocate(this%c3)

   end subroutine cleanup

!> \brief Increase selected counter by 1

   subroutine incr(this, ind)

      implicit none

      class(arrcnt),         intent(inout) :: this  !< an object invoking the type-bound procedure
      integer, dimension(:), intent(in)    :: ind   !< position at which we need +1

      call this%add(ind, 1_LONG)

   end subroutine incr

!> \brief Add some value to the selected counter

   subroutine add(this, ind, val)

      use dataio_pub, only: die

      implicit none

      class(arrcnt),         intent(inout) :: this  !< an object invoking the type-bound procedure
      integer, dimension(:), intent(in)    :: ind   !< position at which we need +val
      integer(LONG),         intent(in)    :: val   !< value to add to selected counter

      select case (size(ind))
         case (0)
            if (.not. allocated(this%c0)) call die("[cnt_array:add] c0 not allocated")
            this%c0 = this%c0 + val
         case (1)
            if (.not. allocated(this%c1)) call die("[cnt_array:add] c1 not allocated")
            this%c1(ind(1)) = this%c1(ind(1)) + val
         case (2)
            if (.not. allocated(this%c2)) call die("[cnt_array:add] c2 not allocated")
            this%c2(ind(1), ind(2)) = this%c2(ind(1), ind(2)) + val
         case (3)
            if (.not. allocated(this%c3)) call die("[cnt_array:add] c3 not allocated")
            this%c3(ind(1), ind(2), ind(3)) = this%c3(ind(1), ind(2), ind(3)) + val
         case default
            call die("[cnt_array:add] rank not implemented")
      end select

   end subroutine add

!> \brief Print the collected values

   subroutine print(this, title, v)

      use dataio_pub, only: warn, printinfo, msg

      implicit none

      class(arrcnt),    intent(in) :: this   !< an object invoking the type-bound procedure
      character(len=*), intent(in) :: title  !< label for the title
      integer(kind=4),  intent(in) :: v      !< verbosity level

      integer :: i, j

      call printinfo(trim(title), v)

      if (allocated(this%c0)) then
         write(msg, '(a,i0)')"  ", this%c0
         call printinfo(trim(msg), v)
      else if (allocated(this%c1)) then
         write(msg, *)"  ", this%c1
         call printinfo(trim(msg), v)
      else if (allocated(this%c2)) then
         do i = lbound(this%c2, 1), ubound(this%c2, 1)
            write(msg, *)"  ", this%c2(i, :)
            call printinfo(trim(msg), v)
         enddo
      else if (allocated(this%c3)) then
         do i = lbound(this%c3, 1), ubound(this%c3, 1)
            do j = lbound(this%c3, 2), ubound(this%c3, 2)
               write(msg, *)"  ", this%c3(i, j, :)
               call printinfo(trim(msg), v)
            enddo
            call printinfo("", v)
         enddo
      else
         call warn("[cnt_array:print] nothing to print")
      endif

   end subroutine print

end module cnt_array
