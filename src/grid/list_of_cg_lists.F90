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

!> \brief This module contains list of grid container lists and related subroutines

module list_of_cg_lists

   use cg_list,   only: cg_list_T

   implicit none

   private
   public :: all_lists

   type cg_list_pointer
      class(cg_list_T), pointer  :: lp
   end type cg_list_pointer

   type all_cg_lists
      type(cg_list_pointer), dimension(:), allocatable :: entries
    contains
      procedure :: print    !< Print all cg lists for diagnostic purposes
      procedure :: register !< Reset (initialize) the given list and add it to the table if unique
      procedure :: forget   !< Erase given cg from all known lists
      procedure :: delete   !< Delete all lists
   end type all_cg_lists

   type(all_cg_lists) :: all_lists

contains

!> \brief Print all cg lists for diagnostic purposes

   subroutine print(this)

      use constants,  only: pSUM
      use dataio_pub, only: msg, printinfo, warn
      use mpisetup,   only: master, piernik_MPI_Allreduce

      implicit none

      class(all_cg_lists), intent(inout) :: this !< object invoking type-bound procedure

      integer :: i, g_cnt

      !> \todo use MPI_Gather and let the master process print everything

      if (.not. allocated(this%entries)) then
         call printinfo("[list_of_cg_lists:print] Unbelievable! No lists at all!")
      else
         if (master) call printinfo("[list_of_cg_lists:print] All known cg_lists:")
         do i = lbound(this%entries(:),dim=1), ubound(this%entries(:), dim=1)
            if (associated(this%entries(i)%lp)) then
               !> \todo Call MPI_Allgather and print detailed distribution of grid pieces across procesors
               g_cnt = this%entries(i)%lp%cnt
               call piernik_MPI_Allreduce(g_cnt, pSUM)
               write(msg, '(3a,i7,a)') "'", this%entries(i)%lp%label, "' : ", g_cnt, " element(s)"
               if (master) call printinfo(msg)
            else
               call warn("(null)")
            endif
         enddo
      endif

   end subroutine print

!> \brief Reset (initialize) the given list and add it to the table if unique

   subroutine register(this, cgl, label)

      implicit none

      class(all_cg_lists),      intent(inout) :: this  !< object invoking type-bound procedure
      class(cg_list_T), target, intent(inout) :: cgl   !< a cg list to be created or reset
      character(len=*),         intent(in)    :: label !< name of the list

      type(cg_list_pointer), dimension(:), allocatable :: new_lists
      integer :: i

      call cgl%init_new(label)

      ! update the list of lists
      if (.not. allocated(this%entries)) then
         allocate(this%entries(1))
      else
         do i = lbound(this%entries(:),dim=1), ubound(this%entries(:), dim=1)
            if (this%entries(i)%lp%label == label) return ! do not duplicate entries (e.g. leaves)
         enddo
         allocate(new_lists(lbound(this%entries(:),dim=1):ubound(this%entries(:), dim=1) + 1))
         new_lists(:ubound(this%entries(:), dim=1)) = this%entries(:)
         call move_alloc(from=new_lists, to=this%entries)
      endif

      this%entries(ubound(this%entries, dim=1))%lp => cgl

   end subroutine register

!> \brief Erase given cg from all known lists

   subroutine forget(this, cg)

      use cg_list,            only: cg_list_element
      use grid_cont,          only: grid_container
      use grid_container_ext, only: cg_extptrs

      implicit none

      class(all_cg_lists),           intent(inout) :: this  !< object invoking type-bound procedure
      type(grid_container), pointer, intent(inout) :: cg    !< grid piece deemed to be forgotten

      type(cg_list_element), pointer :: cgl, aux
      integer :: i, ep

      ! scan all lists and remove the element if found
      do i = lbound(this%entries(:),dim=1), ubound(this%entries(:), dim=1)
         cgl => this%entries(i)%lp%first
         if (this%entries(i)%lp%cnt > 0) then
            do while (associated(cgl))
               aux => cgl
               cgl => cgl%nxt
               if (associated(cg,aux%cg)) call this%entries(i)%lp%delete(aux)
               ! do not call exit here because it is safer to not assume single occurence on a list
            enddo
         endif
      enddo

      call cg%cleanup
      do ep = ubound(cg_extptrs%ext, dim=1), lbound(cg_extptrs%ext, dim=1), -1
         if (associated(cg_extptrs%ext(ep)%cleanup)) call cg_extptrs%ext(ep)%cleanup(cg)
      enddo
      deallocate(cg)

   end subroutine forget

!> \brief Delete all lists

   subroutine delete(this)

      implicit none

      class(all_cg_lists), intent(inout) :: this  !< object invoking type-bound procedure

      integer :: i

      do i = lbound(this%entries(:),dim=1), ubound(this%entries(:), dim=1)
         call this%entries(i)%lp%delete
      enddo
      if (allocated(all_lists%entries)) deallocate(all_lists%entries)

   end subroutine delete

end module list_of_cg_lists
