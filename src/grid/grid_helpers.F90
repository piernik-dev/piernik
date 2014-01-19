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

!> \brief Module containing hepler functions that didn't fit elsewhere

module grid_helpers
   implicit none
   private
   public :: c2f, c2f_o, f2c, f2c_o

contains

   !> \brief Calculate minimal fine grid that completely covers given coarse grid

   pure function c2f(coarse) result (fine)

      use constants, only: xdim, zdim, LO, HI, refinement_factor
      use domain,    only: dom

      implicit none

      integer(kind=8), dimension(xdim:zdim, LO:HI), intent(in) :: coarse

      integer(kind=8), dimension(xdim:zdim, LO:HI) :: fine

      fine(:, LO) = c2f_o(coarse(:, LO))
      where (dom%has_dir(:))
         fine(:, HI) = c2f_o(coarse(:, HI)) + refinement_factor - 1
      elsewhere
         fine(:, HI) = c2f_o(coarse(:, HI))
      endwhere

   end function c2f

!> \brief Calculate minimal coarse grid that completely embeds given fine grid

   pure function f2c(fine) result (coarse)

      use constants, only: xdim, zdim, LO, HI

      implicit none

      integer(kind=8), dimension(xdim:zdim, LO:HI), intent(in) :: fine

      integer(kind=8), dimension(xdim:zdim, LO:HI) :: coarse

      coarse = f2c_o(fine)

   end function f2c

!> \brief Calculate refined offset

   elemental function c2f_o(o_coarse) result (o_fine)

      use constants, only: refinement_factor

      implicit none

      integer(kind=8), intent(in) :: o_coarse

      integer(kind=8) :: o_fine

      o_fine = o_coarse * refinement_factor

   end function c2f_o

!> \brief Calculate coarsened offset

   elemental function f2c_o(o_fine) result (o_coarse)

      use constants, only: refinement_factor

      implicit none

      integer(kind=8), intent(in) :: o_fine

      integer(kind=8) :: o_coarse

      o_coarse = floor(o_fine / real(refinement_factor))

   end function f2c_o

end module grid_helpers
