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

!>
!! \brief List of patches to be decomposed and then turned into grid containers
!!
!! \details Don't forget to describe your stuff
!<

module patch_list

   use decomposition, only: box_T

   implicit none

   private
   public :: patch_list_T

   type :: patch_list_T
      type(box_T), dimension(:), allocatable :: patches  !< list of patches
      !> \todo make it private
   contains
      procedure :: p_deallocate !< Throw out patches list
      procedure :: expand       !< Expand the patch list by one
   end type patch_list_T

contains

!> \brief Throw out patches list

   subroutine p_deallocate(this)

      implicit none

      class(patch_list_T), intent(inout) :: this !< object invoking type bound procedure

      if (allocated(this%patches)) deallocate(this%patches) 
      ! this%patches(:)%pse should be deallocated automagically

   end subroutine p_deallocate

!> \brief Expand the patch list by one

   subroutine expand(this)

      use decomposition, only: box_T

      implicit none

      class(patch_list_T), target, intent(inout) :: this !< current level

      type(box_T), dimension(:), allocatable :: tmp
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

end module patch_list
