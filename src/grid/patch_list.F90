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
   end type patch_list_T

contains

!> \brief Throw out patches list

   subroutine p_deallocate(this)

      implicit none

      class(patch_list_T), intent(inout) :: this !< object invoking type bound procedure

      if (allocated(this%patches)) deallocate(this%patches) 
      ! this%patches(:)%pse should be deallocated automagically

   end subroutine p_deallocate

end module patch_list
