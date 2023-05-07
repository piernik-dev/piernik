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

!> \brief Unified refinement criteria for filters that provide scalar indicator of the need for refinement.

module unified_ref_crit_filter

   use unified_ref_crit, only: urc

   implicit none

   private
   public :: urc_filter

!>
!! \brief Things that should be common for all refinement criteria based on filters.
!! The filters are supposed to decide whether local conditions deserve refinement or derefinement.
!<

   type, abstract, extends(urc) :: urc_filter
! These components should be private but then derived types wouldn't have access to them.
! Submodules doesn't seem to be a real solution for this.
!      private
      real    :: ref_thr    !< refinement threshold
      logical :: plotfield  !< create a 3D array to keep the value of refinement criterion
   end type urc_filter

end module unified_ref_crit_filter
