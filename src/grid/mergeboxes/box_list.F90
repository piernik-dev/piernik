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

!> \brief A list of boxes (cuboids)

module box_list

   use constants, only: xdim, zdim, LO, HI

   implicit none

   private
   public :: box_list_T, box

   type :: box
      integer(kind=8), dimension(xdim:zdim, LO:HI) :: b !< Box description
   contains
      procedure :: sanitize                             !< Sort corners
   end type box

   type :: box_list_T
      type(box), allocatable, dimension(:) :: blist     !< List of boxes
   contains
      procedure :: add_b     !< Add a box
      procedure :: print     !< Print box
      procedure :: cleanup   !< Free memory
   end type box_list_T

contains

!> \brief Swap corners when LO has higher coordinate than HI

   subroutine sanitize(this)

      implicit none

      class(box), intent(inout) :: this

      integer :: d
      integer(kind=8) :: swap

      do d = xdim, zdim
         if (this%b(d, HI) < this%b(d, LO)) then
            swap = this%b(d, HI)
            this%b(d, HI) = this%b(d, LO)
            this%b(d, LO) = swap
         endif
      enddo

   end subroutine sanitize

!> \brief Add a box to the list (expand array by one entry)

   subroutine add_b(this, b)

      implicit none

      class(box_list_T), intent(inout) :: this
      type(box),         intent(in)    :: b

      if (.not. allocated(this%blist)) allocate(this%blist(0))
      this%blist = [ this%blist , b ] ! LHS realloc

      call this%blist(ubound(this%blist, dim=1))%sanitize

   end subroutine add_b

!> \brief Print boxes to stdout (useful only for debugging)

   subroutine print(this)

      use dataio_pub, only: msg, printinfo

      implicit none

      class(box_list_T), intent(inout) :: this

      integer :: i

      do i = lbound(this%blist, dim=1), ubound(this%blist, dim=1)
         write(msg,*)"[box_list:print] #",i," [ ",this%blist(i)%b(:, LO),"]..[",this%blist(i)%b(:, HI),"]",product(this%blist(i)%b(:, HI)-this%blist(i)%b(:, LO)+1)," cells"
         call printinfo(msg)
      enddo

   end subroutine print

!> \brief Free what can possibly be allocated

   subroutine cleanup(this)

      implicit none

      class(box_list_T), intent(inout) :: this

      if (allocated(this%blist)) deallocate(this%blist)

   end subroutine cleanup

end module box_list
