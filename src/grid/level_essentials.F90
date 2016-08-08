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

!> \brief Module containing most basic properties ofgrid levels. Created to avoid circular dependencies between grid_container and cg_level_T

module level_essentials

   use constants, only: ndims

   implicit none
   private
   public :: level_T

   !> \brief the type that contains the most basic data to be linked, where appropriate

   type :: level_T
      integer(kind=8), dimension(ndims) :: off  !< offset of the level
      integer(kind=8), dimension(ndims) :: n_d  !< maximum number of grid cells in each direction (size of fully occupied level)
      integer(kind=4)                   :: id   !< level number (relative to base level). No arithmetic should depend on it.
      !type(level_T), pointer :: coarser
      !type(level_T), pointer :: finer
   contains
      procedure :: init  !< simple initialization
   end type level_T

contains

   subroutine init(this, id, n_d, off)

      use constants,  only: ndims
      use dataio_pub, only: msg, printinfo
      use domain,     only: dom

      implicit none

      class(level_T),                    intent(inout) :: this
      integer(kind=4),                   intent(in)    :: id
      integer(kind=8), dimension(ndims), intent(in)    :: n_d
      integer(kind=8), dimension(ndims), intent(in)    :: off

      this%id  = id
      where (dom%has_dir(:))
         this%n_d(:) = n_d(:)
         this%off(:) = off(:)
      elsewhere
         this%n_d(:) = 1
         this%off(:) = 0
      endwhere

      write(msg, '(a,i4,2(a,3i8),a)')"[level_essentials] Initializing level", this%id, ", size=[", this%n_d, "], offset=[", this%off, "]"
      call printinfo(msg, .false.)

   end subroutine init

end module level_essentials
