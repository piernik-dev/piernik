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

module fluidboundaries_funcs
! pulled by ANY

   implicit none

   private
   public :: default_bnd, user_fluidbnd, init_default_fluidboundaries

   interface

      subroutine user_bnd(dir,side,cg)

         use grid_cont, only: grid_container

         implicit none

         integer(kind=4),               intent(in)    :: dir, side
         type(grid_container), pointer, intent(inout) :: cg

      end subroutine user_bnd

   end interface

   procedure(user_bnd), pointer :: user_fluidbnd

contains

!--------------------------------------------------------------------------------------------------
   subroutine init_default_fluidboundaries

      use constants,  only: PIERNIK_INIT_MPI
      use dataio_pub, only: code_progress, die
#ifdef VERBOSE
      use dataio_pub, only: printinfo
#endif /* VERBOSE */

      implicit none

      if (code_progress < PIERNIK_INIT_MPI) call die("[fluidboundaries_funcs:init_default_fluidboundaries] MPI not initialized.")

#ifdef VERBOSE
      call printinfo("[fluidboundaries_funcs:init_default_fluidboundaries]: commencing...")
#endif /* VERBOSE */

      user_fluidbnd => default_bnd

#ifdef VERBOSE
      call printinfo("[fluidboundaries_funcs:init_default_fluidboundaries]: finished. \o/")
#endif /* VERBOSE */

   end subroutine init_default_fluidboundaries

end module fluidboundaries_funcs
