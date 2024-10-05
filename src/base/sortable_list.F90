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

!> \brief Module that contains a general sorting routine. To be extended for other more sophisticated structures

module sortable_list

   implicit none

   private
   public ::  sortable_list_t, cleanup_sortable_list, temp_index

   !>
   !! \brief Abstract array of sortable elements with sorting method
   !!
   !! \details This type does not contain any actual array since it is too hard to extend it later to something
   !! usable. Every type that extends this one has to provide own sortable array and methods to allocate it,
   !! deallocate it, get bounds, make assignments and comparisons based on integer indices.
   !! Note that we use special index (parameter temp_index) to denote temporary storage for swapping elements.
   !<
   type, abstract :: sortable_list_t
   contains
      procedure(lubound_list), deferred :: l_bound          !< Get lower bound of the list
      procedure(lubound_list), deferred :: u_bound          !< Get upper bound of the list
      procedure(assign_list),  deferred :: assign_element   !< Make an assignment
      procedure(compare_list), deferred :: compare_elements !< Make a comparison
      procedure                         :: sort             !< Sorting routine (currently shellsort)
   end type sortable_list_t

   interface

      subroutine assign_list(this, a, b)
         import sortable_list_t
         class(sortable_list_t), intent(inout) :: this
         integer,                intent(in)    :: a, b
      end subroutine assign_list

      pure logical function compare_list(this, a, b)
         import sortable_list_t
         class(sortable_list_t), intent(in) :: this
         integer,                intent(in) :: a, b
      end function compare_list

      integer function lubound_list(this)
         import sortable_list_t
         class(sortable_list_t), intent(in) :: this
      end function lubound_list

   end interface

   integer, parameter :: temp_index = -huge(1) !< Symbolic index indicating temporary storage

   integer, dimension(:), allocatable :: gaps  !< Auxiliary array for the sorting routine.
   ! It can be expanded and reused, so it is detached from the type sortable_list_t.
   ! allocate(gaps(ceiling(log(1.25*ub + 1) / log(2.25) )))
   ! gaps(:) = [ (tokuda(i), i = 1, size(gaps)) ]

contains

!> \brief deallocate everything locally allocated

   subroutine cleanup_sortable_list

      implicit none

      if (allocated(gaps)) deallocate(gaps)

   end subroutine cleanup_sortable_list

!>
!! \brief Shell sort with nontrivial coefficients
!!
!! \details Gap sequence according to:
!! Tokuda, Naoyuki (1992). "An Improved Shellsort". In van Leeuven, Jan. Proceedings of the IFIP 12th World Computer Congress on Algorithms, Software, Architecture.
!! Alternatively one can use gaps provided by M. Ciura or any other newer sequence but it seems that the gains are incremental at best.
!!
!! \todo Consider rewriting to mergesort if this routine consumes too much CPU power
!<

   subroutine sort(this)

      implicit none

      class(sortable_list_t), intent(inout) :: this

      integer :: g, i, j
      integer :: lb, ub

      lb = this%l_bound()
      ub = this%u_bound()

      if (.not. allocated(gaps)) then
         allocate(gaps(1))
         gaps(1) = 1
      endif

      do while (gaps(ubound(gaps, dim=1)) < ub)
         ! This is quite short list, so we can afford for slowest, simplest algorithm here
         gaps = [ gaps, tokuda(ubound(gaps, dim=1)+1) ]
      enddo

      do g = ubound(gaps, dim=1), lbound(gaps, dim=1), -1
         do i = gaps(g)+1, ub
            call this%assign_element(temp_index, i) !this%temp = this%list(i)
            j = i
!            do while (j > gaps(g) .and. this%list(max(j - gaps(g), lbound(this%list, dim=1)))%id > temp%id)
            do while (j > gaps(g) .and. this%compare_elements(max(j - gaps(g), lb), temp_index))
               ! Either use max() here, or put this%list(j - gaps(g)) with an "if" inside the while loop
               call this%assign_element(j, j-gaps(g)) !this%list(j) = this%list(j - gaps(g))
               j = j - gaps(g)
            enddo
            call this%assign_element(j, temp_index) !this%list(j) = this%temp
         enddo
      enddo

#ifdef DEBUG
      call check
#endif /* DEBUG */

   contains

      !> \brief Return k-th gap according to Tokuda prescription

      integer function tokuda(k)

         use dataio_pub, only: die

         implicit none

         integer, intent(in) :: k

         if (k < 1) call die("[sortable_list:sort:tokuda] k<1")

         tokuda = ceiling(0.8 * (2.25**k - 1.)) ! == ceiling((9**k-4**k)/(5.*4**(k-1)))

      end function tokuda

#ifdef DEBUG
      !> \brief Check if our sorting was performed correctly

      subroutine check

         use dataio_pub, only: msg, warn, die

         implicit none

         logical :: fail

         fail = .false.
         do i = this%l_bound(), this%u_bound() - 1
            if (this%compare_elements(i, i+1)) then
               write(msg,*)"this%list(", i+1, ") < this%list(", i, ")%id"
               call warn(msg)
               fail = .true.
            endif
         enddo
         if (fail) call die("[sortable_list:sort] failed")

      end subroutine check
#endif /* DEBUG */

   end subroutine sort

end module sortable_list
