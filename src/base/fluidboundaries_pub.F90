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
!    Initial implemetation of PIERNIK code was based on TVD split MHD code by
!    Ue-Li Pen
!        see: Pen, Arras & Wong (2003) for algorithm and
!             http://www.cita.utoronto.ca/~pen/MHD
!             for original source code "mhd.f90"
!
!    For full list of developers see $PIERNIK_HOME/license/pdt.txt
!
#include "piernik.def"

module fluidboundaries_pub
   implicit none
   private
   public :: user_bnd_xl, user_bnd_xr, user_bnd_yl, user_bnd_yr, user_bnd_zl, user_bnd_zr, &
      & init_fluidboundaries

   interface
      subroutine user_bnd
      end subroutine user_bnd
   end interface

   procedure(user_bnd), pointer :: user_bnd_xl, user_bnd_xr, user_bnd_yl, user_bnd_yr, user_bnd_zl, user_bnd_zr

contains

!--------------------------------------------------------------------------------------------------
   subroutine default_bnd
      use dataio_public, only: die
      implicit none
      call die("User boundaries are not defined")
   end subroutine default_bnd
!--------------------------------------------------------------------------------------------------
   subroutine init_fluidboundaries
#ifdef VERBOSE
      use dataio_public, only: printinfo
#endif /* VERBOSE */
      implicit none
#ifdef VERBOSE
      call printinfo("[fluidboundaries_pub:init_fluidboundaries]: commencing...")
#endif /* VERBOSE */
      user_bnd_xl => default_bnd
      user_bnd_xr => default_bnd
      user_bnd_yl => default_bnd
      user_bnd_yr => default_bnd
      user_bnd_zl => default_bnd
      user_bnd_zr => default_bnd
#ifdef VERBOSE
      call printinfo("[fluidboundaries_pub:init_fluidboundaries]: finished. \o/")
#endif /* VERBOSE */
   end subroutine init_fluidboundaries

end module fluidboundaries_pub
