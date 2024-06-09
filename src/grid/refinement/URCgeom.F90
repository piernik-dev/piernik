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
!! \brief Unified refinement criteria for geometrical primitives
!!
!! \details Currently there are points, boxes and vertical cylinders implemented.
!! It should be relatively easy to add other shapes, like sphere or shell if needed.
!<

module unified_ref_crit_geometrical

   use unified_ref_crit, only: urc

   implicit none

   private
   public :: urc_geom

   !> \brief Things that should be common for all refinement criteria based on geometrical primitives.

   type, abstract, extends(urc) :: urc_geom
      integer :: level  !< desired level of refinement
   contains
      procedure :: enough_level       !< Level check for all dependent types.
      procedure, nopass :: coord2ind  !< Convert cordinates to indices.
   end type urc_geom

contains

   !> \brief Level check for all dependent types.

   pure logical function enough_level(this, lev)

      implicit none

      class(urc_geom), intent(in) :: this
      integer(kind=4), intent(in) :: lev

      enough_level = (lev >= this%level)

   end function enough_level

   !> \brief Convert coordinates to indices, perhaps it can go to some more general place

   pure function coord2ind(coords, l) result(ijk)

      use constants,        only: ndims, LO
      use domain,           only: dom
      use level_essentials, only: level_t

      implicit none

      real, dimension(ndims),  intent(in) :: coords
      class(level_t), pointer, intent(in) :: l

      integer(kind=8), dimension(ndims) :: ijk

      where (dom%has_dir)
         ijk(:) = l%off + floor((coords - dom%edge(:, LO))/dom%L_ * l%n_d)
         ! Excessively large this%coords will result in FPE exception.
         ! If FPE is not trapped, then huge() value will be assigned from uninit constant
         ! (checked on gfortran 7.3.1), which is safe.
         ! A wrapped value coming from integer overflow may be unsafe.
      elsewhere
         ijk(:) = l%off
      endwhere

   end function coord2ind

end module unified_ref_crit_geometrical
