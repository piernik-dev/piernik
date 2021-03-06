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

!> \brief List of patches to be decomposed and then turned into grid containers

module patch_list

   use decomposition, only: box_t

   implicit none

   private
   public :: patch_list_t

   type :: patch_list_t
      type(box_t), dimension(:), allocatable :: patches  !< list of patches
   contains
      procedure :: p_deallocate !< Throw out patches list
      procedure :: expand       !< Expand the patch list by one
      procedure :: p_count      !< Count local patches
      procedure :: p2a          !< Copy patches to an array
   end type patch_list_t

contains

!> \brief Throw out patches list

   subroutine p_deallocate(this)

      implicit none

      class(patch_list_t), intent(inout) :: this !< object invoking type bound procedure

      if (allocated(this%patches)) deallocate(this%patches)
      ! this%patches(:)%pse should be deallocated automagically

   end subroutine p_deallocate

!> \brief Expand the patch list by one

   subroutine expand(this)

      use decomposition, only: box_t

      implicit none

      class(patch_list_t), target, intent(inout) :: this !< current level

      type(box_t), dimension(:), allocatable :: tmp
      integer :: i

      if (.not. allocated(this%patches)) then
         allocate(this%patches(1))
      else
         allocate(tmp(lbound(this%patches(:),dim=1):ubound(this%patches(:), dim=1) + 1))
         tmp(:ubound(this%patches(:), dim=1)) = this%patches(:)
         ! manually deallocate arrays inside user-types, as it seems that move_alloc is unable to do that
         do i = lbound(this%patches(:), dim=1), ubound(this%patches(:), dim=1)
            if (allocated(this%patches(i)%pse)) deallocate(this%patches(i)%pse)
         enddo
         call move_alloc(from=tmp, to=this%patches)
      endif

   end subroutine expand

!> \brief Count local patches

   integer function p_count(this)

      implicit none

      class(patch_list_t), intent(in) :: this

      integer :: p

      p_count = 0
      if (allocated(this%patches)) then
         do p = lbound(this%patches(:), dim=1), ubound(this%patches(:), dim=1)
            p_count = p_count + size(this%patches(p)%pse, dim=1)
         enddo
      endif

   end function p_count

!> Copy patches to an array

   subroutine p2a(this, gptemp)

      use constants, only: LO, HI

      implicit none

      class(patch_list_t),             intent(in)    :: this
      integer(kind=8), dimension(:,:), intent(inout) :: gptemp

      integer :: i, p, ss

      i = 0
      if (allocated(this%patches)) then
         do p = lbound(this%patches(:), dim=1), ubound(this%patches(:), dim=1)
            do ss = lbound(this%patches(p)%pse, dim=1), ubound(this%patches(p)%pse, dim=1)
               i = i + 1
               gptemp(:, i) = [ this%patches(p)%pse(ss)%se(:, LO), this%patches(p)%pse(ss)%se(:, HI) - this%patches(p)%pse(ss)%se(:, LO) + 1 ]
            enddo
         enddo
      endif

   end subroutine p2a

end module patch_list
