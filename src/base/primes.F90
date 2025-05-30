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

!> \brief This module serves an array of prime numbers. Used primarily in decomposition module.

module primes_utils

   implicit none

   private
   public :: primes_t

   !> \brief The type that contain table of primes and initializer
   type :: primes_t
      integer(kind=4), allocatable, dimension(:) :: tab !< table of prime numbers
      integer, private :: max                           !< max number to which the search was performed
   contains
      procedure :: sieve                                !< routine used to initialize (and extend if necessary) the table of prime numbers
      procedure :: erase                                !< restore initial state
      procedure, private :: print                       !< print what was found
   end type primes_t

contains

!>
!! \brief Eratosthenes sieve
!!
!! \details On first call it creates a table of primes up to a given number (does an initialization).
!! On subsequent calls it may extend the table with primes, if required.
!<

   subroutine sieve(this, n)

      use constants,  only: I_ZERO, I_TWO
      use dataio_pub, only: die

      implicit none

      class(primes_t), intent(inout) :: this !< object invoking type-bound procedure
      integer(kind=4), intent(in)    :: n    !< max value of prime number

      integer(kind=4), dimension(n) :: numb
      integer(kind=4) :: i

      if (allocated(this%tab)) then ! we may expect this%max was already set
         if (n <= this%max) return
         deallocate(this%tab)
      endif

      if (n < 1) call die("[primes:sieve] Prime numbers are usually greater than 1 ...")

      ! The sieve
      numb = [I_ZERO, (i, i = I_TWO, n)]
      do i = 2, n
         if (numb(i) /= 0) numb( 2*i : n : i ) = 0
      enddo

      ! Gather the primes into a table
      allocate(this%tab(count(numb /= 0)))
      this%tab(:) = pack(numb, numb /= 0)
      this%max = n

#ifdef DEBUG
      call this%print
#endif /* DEBUG */

   end subroutine sieve

!> \brief Print the table of found primes

   subroutine print(this)

      use constants,  only: V_DEBUG, fmt_len
      use dataio_pub, only: msg, printinfo
      use mpisetup,   only: master

      implicit none

      class(primes_t), intent(inout) :: this !< object invoking type-bound procedure

      integer, parameter :: head_primenum = 10, tail_primenum = 5
      character(len=fmt_len) :: fmt

      if (master) then
         if (ubound(this%tab, 1) <= 0) then
            write(msg, '(a,i0)') "There are 0 prime numbers equal to or smaller than ", this%max
         else if (ubound(this%tab, 1) <= head_primenum + tail_primenum) then
            write(fmt, '(a,i0,a)')'(a,i0,a,i0,a,', ubound(this%tab, 1), '(" ",i0))'
            write(msg, fmt) "There are ", size(this%tab), " prime numbers equal to or smaller than ", this%max," :", this%tab
         else
            write(fmt, '(a,2(i0,a))')'(a,i0,a,i0,a,', head_primenum, '(" ",i0),a,', tail_primenum, '(" ",i0))'
            write(msg, fmt) "There are ", size(this%tab), " prime numbers equal to or smaller than ", this%max," :", this%tab(1:head_primenum), " .. ", this%tab(ubound(this%tab, 1)-tail_primenum+1:ubound(this%tab, 1))
         endif
         call printinfo(msg, V_DEBUG)
      endif

   end subroutine print

!> \brief Restore initial state

   subroutine erase(this)

      use constants, only: INVALID

      implicit none

      class(primes_t), intent(inout) :: this !< object invoking type-bound procedure

      deallocate(this%tab)
      this%max = INVALID

   end subroutine erase

end module primes_utils
