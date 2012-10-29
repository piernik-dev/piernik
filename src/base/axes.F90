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

!> \brief A type with 3 coordinate axes, used where axes contained in grid_container are not applicable

module axes_M

   implicit none

   private
   public :: axes

   type :: axes

      real, allocatable, dimension(:) :: x      !< array of x-positions of %grid cells centers
      real, allocatable, dimension(:) :: y      !< array of y-positions of %grid cells centers
      real, allocatable, dimension(:) :: z      !< array of z-positions of %grid cells centers

   contains

      procedure :: allocate_axes
      procedure :: deallocate_axes

   end type axes

contains

   subroutine allocate_axes(this, sizes)

      use constants, only: ndims, xdim, ydim, zdim

      implicit none

      class(axes),                       intent(inout) :: this
      integer(kind=4), dimension(ndims), intent(in)    :: sizes

      if (.not.allocated(this%x)) allocate(this%x(sizes(xdim)))
      if (.not.allocated(this%y)) allocate(this%y(sizes(ydim)))
      if (.not.allocated(this%z)) allocate(this%z(sizes(zdim)))

   end subroutine allocate_axes

   subroutine deallocate_axes(this)

      implicit none

      class(axes), intent(inout) :: this

      if (allocated(this%x)) deallocate(this%x)
      if (allocated(this%y)) deallocate(this%y)
      if (allocated(this%z)) deallocate(this%z)

   end subroutine deallocate_axes

end module axes_M
