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
   public :: arrcnt, arrsum

   type :: arrcnt  ! array of counters
      ! There's no polymorphic way to declare only one variable here in Fortran 2018.
      ! We're supposed to use only one counter at a time.
      integer(LONG),                   allocatable, private :: c0  ! single (scalar) counter
      integer(LONG), dimension(:),     allocatable, private :: c1  ! vector of counters
      integer(LONG), dimension(:,:),   allocatable, private :: c2  ! array (matrix) of counters
      integer(LONG), dimension(:,:,:), allocatable, private :: c3  ! rank-3 array of counters
      logical, private :: used  ! will remain .false. in there were no calls to incr or add
   contains
      procedure :: init     !< allocate the array of counters
      procedure :: cleanup  !< deallocate everything
      procedure :: incr     !< increase selected counter by 1
      procedure :: add      !< add some value to the selected counter
      procedure :: print    !< print the collected values
      procedure, private :: numdigits !< find minimum required field width
      procedure :: get_v1   !< get vector of counters
   end type arrcnt

   type, extends(arrcnt) :: arrsum  ! array of sumators + array of counters
      type(arrcnt) :: cnt
   end type arrsum

contains

!> \brief Allocate the array of counters

   recursive subroutine init(this, sizes)

      use dataio_pub, only: die

      implicit none

      class(arrcnt),                   intent(inout) :: this   !< an object invoking the type-bound procedure
      integer, dimension(:), optional, intent(in)    :: sizes  !< array of sizes for this%c (1-based)

      if (allocated(this%c0) .or. allocated(this%c1) .or. allocated(this%c2) .or. allocated(this%c3)) &
           call die("[cnt_array:init] already initialized")

      this%used = .false.

      if (present(sizes)) then
         select case (size(sizes))
            case (0)  ! the harder way for a scalar counter
               allocate(this%c0)
               this%c0 = 0
            case (1)
               allocate(this%c1(sizes(1)))
               this%c1 = 0
            case (2)
               allocate(this%c2(sizes(1), sizes(2)))
               this%c2 = 0
            case (3)  ! most likely never used
               allocate(this%c3(sizes(1), sizes(2), sizes(3)))
               this%c3 = 0
            case default
               call die("[cnt_array:init] rank not implemented")
         end select
      else  ! the easier way for a scalar counter
         allocate(this%c0)
         this%c0 = 0
      endif

      select type (this)
         type is (arrsum)
            call this%cnt%init(sizes)
      end select

   end subroutine init

!> \brief Deallocate everything

   recursive subroutine cleanup(this)

      implicit none

      class(arrcnt), intent(inout) :: this  !< an object invoking the type-bound procedure

      if (allocated(this%c0)) deallocate(this%c0)
      if (allocated(this%c1)) deallocate(this%c1)
      if (allocated(this%c2)) deallocate(this%c2)
      if (allocated(this%c3)) deallocate(this%c3)

      select type (this)
         type is (arrsum)
            call this%cnt%cleanup
      end select

   end subroutine cleanup

!> \brief Increase selected counter by 1

   subroutine incr(this, ind)

      implicit none

      class(arrcnt),                   intent(inout) :: this  !< an object invoking the type-bound procedure
      integer, dimension(:), optional, intent(in)    :: ind   !< position at which we need +1

      call this%add(ind, 1_LONG)

   end subroutine incr

!> \brief Add some value to the selected counter

   recursive subroutine add(this, ind, val)

      use dataio_pub, only: die

      implicit none

      class(arrcnt),                   intent(inout) :: this  !< an object invoking the type-bound procedure
      integer, dimension(:), optional, intent(in)    :: ind   !< position at which we need +val
      integer(LONG),                   intent(in)    :: val   !< value to add to selected counter

      if (present(ind)) then
         select case (size(ind))
            case (0)
               if (.not. allocated(this%c0)) call die("[cnt_array:add] c0 not allocated (1)")
               this%c0 = this%c0 + val
            case (1)
               if (.not. allocated(this%c1)) call die("[cnt_array:add] c1 not allocated")
               this%c1(ind(1)) = this%c1(ind(1)) + val
            case (2)
               if (.not. allocated(this%c2)) call die("[cnt_array:add] c2 not allocated")
               this%c2(ind(1), ind(2)) = this%c2(ind(1), ind(2)) + val
            case (3)  ! most likely never used
               if (.not. allocated(this%c3)) call die("[cnt_array:add] c3 not allocated")
               this%c3(ind(1), ind(2), ind(3)) = this%c3(ind(1), ind(2), ind(3)) + val
            case default
               call die("[cnt_array:add] rank not implemented")
         end select
      else
         if (.not. allocated(this%c0)) call die("[cnt_array:add] c0 not allocated (2)")
         this%c0 = this%c0 + val
      endif

      this%used = .true.

      select type (this)
         type is (arrsum)
            call this%cnt%incr(ind)
      end select

   end subroutine add

!> \brief Print the collected values

   subroutine print(this, title, v)

      use constants,  only: fmt_len
      use dataio_pub, only: warn, printinfo

      implicit none

      class(arrcnt),    intent(in) :: this   !< an object invoking the type-bound procedure
      character(len=*), intent(in) :: title  !< label for the title
      integer(kind=4),  intent(in) :: v      !< verbosity level

      integer :: i, j
      character(len=fmt_len) :: fmt  ! for formats like '(a,5(' ',i3,'(',i3,')'))'

      if (.not. this%used) then
         call printinfo(trim(title) // " (no counts)", v)
         return
      endif

      call printinfo(trim(title), v)

      ! find the number of entries per line
      if (allocated(this%c0)) then
         i = 1
      else if (allocated(this%c1)) then
         i = size(this%c1)
      else if (allocated(this%c2)) then
         i = size(this%c2, 2)
      else if (allocated(this%c3)) then  ! most likely never used
         i = size(this%c3, 3)
      else
         call warn("[cnt_array:print] not implemented")
      endif

      select type (this)
         type is (arrsum)  ! formats like '(a,5(' ',i3,'(',i3,')'))'
            write(fmt, '(a,3(i0,a))')"(a,", i, "(' ',i", this%numdigits(), ",'( ',i", this%cnt%numdigits(), ",')'))"
         type is (arrcnt)  ! formats like '(a,5(' ',i3))'
            write(fmt, '(a,2(i0,a))')"(a,", i, "(' ',i", this%numdigits(), "))"
      end select

      if (allocated(this%c0)) then
         select type (this)
            type is (arrsum)
               call data_line([this%c0], fmt, [this%cnt%c0])
            type is (arrcnt)
               call data_line([this%c0], fmt)
         end select
      else if (allocated(this%c1)) then
         select type (this)
            type is (arrsum)
               call data_line(this%c1, fmt, this%cnt%c1)
            type is (arrcnt)
               call data_line(this%c1, fmt)
         end select
      else if (allocated(this%c2)) then
         do i = lbound(this%c2, 1), ubound(this%c2, 1)
            select type (this)
               type is (arrsum)
                  call data_line(this%c2(i, :), fmt, this%cnt%c2(i, :))
               type is (arrcnt)
                  call data_line(this%c2(i, :), fmt)
            end select
         enddo
      else if (allocated(this%c3)) then  ! most likely never used
         do i = lbound(this%c3, 1), ubound(this%c3, 1)
            do j = lbound(this%c3, 2), ubound(this%c3, 2)
               select type (this)
                  type is (arrsum)
                     call data_line(this%c3(i, j, :), fmt, this%cnt%c3(i, j, :))
                  type is (arrcnt)
                     call data_line(this%c3(i, j, :), fmt)
               end select
            enddo
            call printinfo("", v)
         enddo
      else
         call warn("[cnt_array:print] nothing to print")
      endif

   contains

!> \brief print single line of counters

      subroutine data_line(d, f, c)

         use dataio_pub, only: msg

         implicit none

         integer(LONG), dimension(:),           intent(in) :: d  ! data line to print
         character(len=*),                      intent(in) :: f  ! format
         integer(LONG), dimension(:), optional, intent(in) :: c  ! extra data line (counters) to print

         integer :: i

         if (present(c)) then
            write(msg, f)"  ", (d(i), c(i), i = 1, ubound(d, 1))
         else
            write(msg, f)"  ", d
         endif
         call hide_zeroes(msg)
         call printinfo(trim(msg), v)

      end subroutine data_line

!> \brief replace "0( 0)" by "." in the log for better readability

      subroutine hide_zeroes(m)

         implicit none

         character(len=*), intent(inout) :: m

         integer :: i

         do while (index(m, " 0( ") /= 0)
            i = index(m, " 0( ") + 1
            m(i:i+1) = "  "
         enddo

         do while (index(m, " 0) ") /= 0)
            i = index(m, " 0) ") + 1
            m(i:i+1) = " ."
         enddo

      end subroutine hide_zeroes

   end subroutine print

!>
!! \brief compute the number of characters required for printing any element of integer array a
!!
!! Only "-0" won't fit, if occur
!<

   integer function numdigits(this)

      use constants,  only: INVALID
      use dataio_pub, only: warn

      implicit none

      class(arrcnt), intent(in) :: this   !< an object invoking the type-bound procedure

      if (allocated(this%c0)) then
         numdigits = int(log10(real(       abs(this%c0)+1 ))) + merge(2, 1,     this%c0 < 0 )
      else if (allocated(this%c1)) then
         numdigits = int(log10(real(maxval(abs(this%c1)+1)))) + merge(2, 1, any(this%c1 < 0))
      else if (allocated(this%c2)) then
         numdigits = int(log10(real(maxval(abs(this%c2)+1)))) + merge(2, 1, any(this%c2 < 0))
      else if (allocated(this%c3)) then  ! most likely never used
         numdigits = int(log10(real(maxval(abs(this%c3)+1)))) + merge(2, 1, any(this%c3 < 0))
      else
         call warn("[cnt_array:numdigits] not implemented")
         numdigits = INVALID
      endif

   end function numdigits

!>
!! \brief get vector of counters
!!
!! An ad-hoc implementation for log_bins. ToDo: generalize it.
!<

   function get_v1(this)

      use dataio_pub, only: die

      implicit none

      class(arrcnt), intent(in) :: this   !< an object invoking the type-bound procedure

      integer(LONG), dimension(size(this%c1)) :: get_v1

      if (.not. allocated(this%c1)) call die("[cnt_array:get_v1] not allocated")

      get_v1 = this%c1

   end function get_v1


end module cnt_array
