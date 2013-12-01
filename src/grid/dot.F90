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
!! \brief Depiction of Topology
!!
!! \details The structure that contains most important information of all blocks on all processes on a given level 
!! and information on decomposition as well.
!<

module dot

   use decomposition,   only: cuboid

   implicit none

   private
   public :: dot_T

   !> \brief A list of grid pieces (typically used as a list of all grids residing on a given process)
   type :: cuboids
      type(cuboid), allocatable, dimension(:) :: c !< an array of grid piece
   end type cuboids

   !> \brief Depiction of global Topology of a level. Use with care, because this is an antiparallel thing
   type :: dot_T
      type(cuboids), dimension(:), allocatable :: gse !< lists of grid chunks on each process (FIRST:LAST)
   contains
   end type dot_T

contains

end module dot
