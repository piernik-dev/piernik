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

!> \brief Module that contains a definition of sortable list of grid pieces, a sorting routine and other useful small routines

module sort_piece_list

   use constants, only: ndims

   implicit none

   private
   public :: grid_piece_list, cleanup_piece_list

   type :: grid_piece
      integer(kind=8), dimension(ndims) :: off       !< offset
      integer(kind=4), dimension(ndims) :: n_b       !< size
      integer(kind=4)                   :: cur_gid   !< current grid_id
      integer(kind=4)                   :: cur_proc  !< current process number
      integer(kind=4)                   :: dest_proc !< process number according to ideal ordering
      integer(kind=8)                   :: id        !< unique number used to sort
      real                              :: weight    !< number of cells relative to total number of cells on given level (unused)
      real                              :: cweight   !< cumulative weight for id <= own id (unused)
   contains
      procedure :: set_gp                            !< Set the primary properties, initialize derived properties with safe defaults
   end type grid_piece

   type :: grid_piece_list
      type(grid_piece), dimension(:), allocatable :: list !< the list itself
   contains
      procedure :: init                                   !< Allocate the list
      procedure :: cleanup                                !< Deallocate the list
      procedure :: sort                                   !< Sorting routine (currently shellsort)
      procedure :: set_id                                 !< Find grid id using a space-filling curve
      procedure :: set_weights                            !< Find estimates of the cost of the grids
   end type grid_piece_list

   integer, dimension(:), allocatable :: gaps ! Auxiliary array for the sorting routine. Can be expanded and reused, so it is detached from the type grid_piece_list

contains

   subroutine set_gp(this, off, n_b, gid, proc)

      use constants, only: ndims, INVALID

      implicit none

      class(grid_piece),                 intent(inout) :: this
      integer(kind=8), dimension(ndims), intent(in)    :: off   !< offset
      integer(kind=4), dimension(ndims), intent(in)    :: n_b   !< size
      integer(kind=4),                   intent(in)    :: gid   !< current grid_id (unused)
      integer(kind=4),                   intent(in)    :: proc  !< current process number (unused)

      this%off       = off
      this%n_b       = n_b
      this%cur_gid   = gid
      this%cur_proc  = proc
      this%dest_proc = INVALID
      this%id        = INVALID
      this%weight    = 0.
      this%cweight   = 0.

   end subroutine set_gp

   !> \brief deallocate everything locally allocated

   subroutine cleanup_piece_list

      implicit none

      if (allocated(gaps)) deallocate(gaps)

   end subroutine cleanup_piece_list

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

!>
!! \brief Shell sort with nontrivial coefficients
!!
!! \details Gap sequence according to:
!! Tokuda, Naoyuki (1992). "An Improved Shellsort". In van Leeuven, Jan. Proceedings of the IFIP 12th World Computer Congress on Algorithms, Software, Architecture.
!! Alternatively one can use gaps provided by M. Ciura  "Best Increments for the Average Case of Shellsort".
!! Proceedings of the 13th International Symposium on Fundamentals of Computation Theory. London: pp. 106~117. ISBN 3-540-42487-3.
!! The Tokuda's gaps were chosen here due to their algebraic prescription.
!! Ciura's gaps are: [ 1, 4, 10, 23, 57, 132, 301, 701 ] and were obtained by numerical searching of optimal average sorting time of the algorithm
!! Larger Ciura's gaps can be obtained approximately with multiplier 2.25
!!
!! \todo Consider rewriting to mergesort if this routine consumes too much CPU power
!<

   subroutine sort(this)

#ifdef DEBUG
      use dataio_pub, only: msg, warn, die
#endif /* DEBUG */

      implicit none

      class(grid_piece_list), intent(inout) :: this

      integer :: g, i, j
      type(grid_piece) :: temp
#ifdef DEBUG
      logical :: fail
#endif /* DEBUG */

      if (.not. allocated(gaps)) then
         allocate(gaps(1))
         gaps(1) = 1
      endif

      do while (gaps(ubound(gaps, dim=1)) < ubound(this%list, dim=1))
         gaps = [ gaps, tokuda(ubound(gaps, dim=1)+1) ]
      enddo

      do g = ubound(gaps, dim=1), lbound(gaps, dim=1), -1
         do i = gaps(g)+1, ubound(this%list, dim=1)
            temp = this%list(i)
            j = i
            do while (j > gaps(g) .and. this%list(max(j - gaps(g), lbound(this%list, dim=1)))%id > temp%id)
               ! Either use max() here, or put this%list(j - gaps(g)) with an "if" inside the while loop
               this%list(j) = this%list(j - gaps(g))
               j = j - gaps(g)
            enddo
            this%list(j) = temp
         enddo
      enddo

#ifdef DEBUG
      fail = .false.
      do i = lbound(this%list, dim=1), ubound(this%list, dim=1) - 1
         if (this%list(i+1)%id < this%list(i)%id) then
            write(msg,*)"this%list(",i+1,")%id = ",this%list(i+1)%id ," < ", this%list(i)%id," = this%list(",i,")%id"
            call warn(msg)
            fail = .true.
         endif
      enddo
      if (fail) call die("[sort_piece_list:sort] failed")
#endif /* DEBUG */

   contains

      integer function tokuda(k)

         use dataio_pub, only: die

         implicit none

         integer, intent(in) :: k

         if (k < 1) call die("[sort_piece_list:sort:tokuda] k<1")

         tokuda = ceiling((9**k-4**k)/(5.*4**(k-1)))

      end function tokuda

   end subroutine sort

!> \brief Find grid id using a space-filling curve.

   subroutine set_id(this, off)

      use ordering, only: Morton_order

      implicit none

      class(grid_piece_list),            intent(inout) :: this
      integer(kind=8), dimension(ndims), intent(in)    :: off  !< offset of the level

      integer :: s

      do s = lbound(this%list, dim=1), ubound(this%list, dim=1)
         !> \todo use Hilbert ordering
         this%list(s)%id = Morton_order(this%list(s)%off-off)
      enddo

   end subroutine set_id

!> \brief Find estimates of the cost of the grids
!> \todo Do we want to include count of cg%leafmap, when available?

   subroutine set_weights(this)

      implicit none

      class(grid_piece_list),            intent(inout) :: this

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

end module sort_piece_list
