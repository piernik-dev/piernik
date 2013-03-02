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

!> \brief A list of convex corners

module corner_list

   use constants, only: ndims

   implicit none

   private
   public :: corner_list_T, corner

   type :: corner
      integer(kind=8), dimension(ndims) :: pos  !< Position of the convex corner
      integer(kind=8), dimension(ndims) :: ivec !< Direction to the interior of the set
      real :: dist                              !< Importance factor, e.g. distance from CoM
   end type corner

   type :: corner_list_T
      type(corner), allocatable, dimension(:) :: clist !< List of convex corners
   contains
      procedure :: add_c        !< Add a convex corner
      procedure :: print        !< Print box
      procedure :: cleanup      !< Free memory
      procedure :: update_dists !< update importance factors
   end type corner_list_T

contains

!> \brief Add a convex corner to the list (expand array by one entry)

   subroutine add_c(this, c)

      implicit none

      class(corner_list_T), intent(inout) :: this
      type(corner),         intent(in)    :: c

      if (.not. allocated(this%clist)) allocate(this%clist(0))
      this%clist = [ this%clist , c ] ! LHS realloc

   end subroutine add_c

!> \brief Print corners to stdout (useful only for debugging)

   subroutine print(this)

      use dataio_pub, only: warn, printinfo, msg

      implicit none

      class(corner_list_T), intent(inout) :: this

      integer :: i

      if (.not. allocated(this%clist)) then
         call warn("[corner_list:print] not allocated")
         return
      endif
      if (size(this%clist) < 1) then
         call warn("[corner_list:print] size == 0")
         return
      endif

      do i = lbound(this%clist, dim=1), ubound(this%clist, dim=1)
         write(msg,*)"[corner_list:print] #",i," @ ",this%clist(i)%pos," -> ",this%clist(i)%ivec, " dist(CoM) ", this%clist(i)%dist
         call printinfo(msg)
      enddo

   end subroutine print

!> \brief Free what can possibly be allocated

   subroutine cleanup(this)

      implicit none

      class(corner_list_T), intent(inout) :: this

      if (allocated(this%clist)) deallocate(this%clist)

   end subroutine cleanup

!> Recalculate importance factors: distance from Center of Mass

   subroutine update_dists(this, CoM)

      use constants, only: ndims

      implicit none

      class(corner_list_T),   intent(inout) :: this
      real, dimension(ndims), intent(in)    :: CoM

      integer :: i

      do i = lbound(this%clist, dim=1), ubound(this%clist, dim=1)
         this%clist(i)%dist = sqrt(sum((this%clist(i)%pos - CoM)**2)) ! sqrt can be dropped
      enddo

   end subroutine update_dists

end module corner_list
