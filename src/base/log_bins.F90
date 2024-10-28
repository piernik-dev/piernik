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
!! \brief Type for array of logarithmically spaced bins
!!
!! \details The logarithmic bins for collecting e.g. Waitall execution time:
!!
!! [ < min,                                                   ! first bin
!!   min .. min * width,                                       ! first regular bin of declared width
!!   ...,                                                      ! regular bins of declared width
!!   min * width ** (nbins - 3) .. min * width ** (nbins - 2), ! last regular bin of declared width
!!   >= min * width ** (nbins - 2) ]                           ! last bin
!<

module log_bins

   use cnt_array, only: arrcnt

   implicit none

   private
   public :: lbins

   type :: lbins  ! array of logarithmically spaced bins
      real, private :: min     !< first bin will contain everything below this value
      real, private :: lwidth  !< each bin from 2 to (nbins - 1) has this width
      integer, private :: n    !< remember the size
      type(arrcnt) :: bins     !< the bins
      real, dimension(:), allocatable, private :: bnd  !< bin boundaries
   contains
      procedure :: init     !< allocate the bins
      procedure :: cleanup  !< deallocate everything
      procedure :: put      !< update bin counters according to given value
      procedure :: print    !< print the bins
   end type lbins

contains

!> \brief Allocate the array of bins

   recursive subroutine init(this, min, width, n)

      use dataio_pub, only: die

      implicit none

      class(lbins), intent(inout) :: this   !< an object invoking the type-bound procedure
      real,         intent(in)    :: min    !< set this%min
      real,         intent(in)    :: width  !< log(width)
      integer,      intent(in)    :: n      !< the number of bins

      integer :: i

      ! n=2 is an edge case: [ < min, >= min ]
      if (n < 2) call die("[log_bins:init] nonsense: this is a counter, not bins")

      if (width <= 1.) call die("[log_bins:init] nonsense: width of logarithmic bins has to be > 1.")

      this%min = min
      this%lwidth = log(width)
      this%n = n
      call this%bins%init([n])

      if (allocated(this%bnd)) call die("[log_bins:init] already initialized")
      allocate(this%bnd(n-1))

      this%bnd = [(this%min * width ** i, i = 0, n - 2)]

   end subroutine init

!> \brief Deallocate everything

   recursive subroutine cleanup(this)

      implicit none

      class(lbins), intent(inout) :: this  !< an object invoking the type-bound procedure

      call this%bins%cleanup
      if (allocated(this%bnd)) deallocate(this%bnd)

   end subroutine cleanup

!> \brief Add some value to the selected counter

   recursive subroutine put(this, val)

      implicit none

      class(lbins), intent(inout) :: this  !< an object invoking the type-bound procedure
      real,         intent(in)    :: val   !< value to be used to update bins

      integer :: ind

      ! ToDo: binary search would be faster for large arrays
      if (val >= this%bnd(ubound(this%bnd, 1))) then
         ind = this%n
      else
         do ind = lbound(this%bnd, 1), ubound(this%bnd, 1)
            if (val < this%bnd(ind)) exit
         enddo
      endif

      call this%bins%incr([ind])

   end subroutine put

!> \brief Print the collected values

   subroutine print(this, title, v)

      use constants,  only: fmt_len, LONG
      use dataio_pub, only: msg, warn, printinfo

      implicit none

      class(lbins),     intent(in) :: this   !< an object invoking the type-bound procedure
      character(len=*), intent(in) :: title  !< label for the title
      integer(kind=4),  intent(in) :: v      !< verbosity level

      integer :: maxn
      integer(LONG), dimension(this%n) :: bindata
      character(len=fmt_len) :: fmt  ! for formats like '(a,5(' ',i3,'(',i3,')'))'

      call printinfo(trim(title), v)

      ! Since msg has finite length we may run into errors here for large bin sets.
      maxn = min(this%n, int(len(msg) / 13.) - 1)  ! 13 is the current width occupied by one entry, with spacing
      if (maxn < this%n) then
         write(msg, '(a,2(i0,a))')"[log_bins:print] only first ", maxn, " out of ", this%n, " can be printed"
         call warn(msg)
      endif

      write(fmt, '(a,2(i0,a))')"(a,", maxn - 1, "('   ',es8.2))"
      write(msg, fmt)"    ", this%bnd(1:maxn-1)
      call printinfo(trim(msg), v)

      write(fmt, '(a,2(i0,a))')"(a,", maxn - 1, "(i8,' | '),i8)"
      bindata = this%bins%get_v1()
      write(msg, fmt)"  ", bindata(1:maxn)
      call printinfo(trim(msg), v)

   end subroutine print

end module log_bins
