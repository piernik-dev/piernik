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
!! \brief Definitions of some compound types used in several other modules
!<
module types

   use constants, only: ndims, LO, HI

   implicit none

   private
   public :: axes, value, cdd

   type :: value
      real                      :: val
      real                      :: assoc
      real,    dimension(ndims) :: coords
      integer, dimension(ndims) :: loc
      integer                   :: proc
   end type value

   type :: axes
      real, allocatable, dimension(:) :: x      !< array of x-positions of %grid cells centers
      real, allocatable, dimension(:) :: y      !< array of y-positions of %grid cells centers
      real, allocatable, dimension(:) :: z      !< array of z-positions of %grid cells centers
   end type axes

   type cart_decomposition
      integer(kind=4)                          :: comm3d  !< cartesian communicator
      integer(kind=4), dimension(ndims)        :: psize   !< number of divisions in each direction
      integer(kind=4), dimension(ndims)        :: pcoords !< own process coordinates within psize(:)-shaped array of processes
      integer(kind=4), dimension(ndims, LO:HI) :: procn   !< array of neighbours proc numbers
      integer(kind=4)                          :: procxyl !< neighbour in corner boundaries
      integer(kind=4)                          :: procyxl !< neighbour in corner boundaries
   end type cart_decomposition

   type(cart_decomposition) :: cdd !< Cartesian Domain Decomposition stuff \todo find a bertter place and protect it

end module types
