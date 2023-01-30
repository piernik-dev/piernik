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

!> \brief This module contains ext_fluxes type useful for flux enforcing and exchange in 1D (M)HD solvers

module fluxtypes

   use flx_cell, only: fluxpoint

   implicit none

   private
   public  :: ext_fluxes

   !> \brief Structure that may contain pointers to fluxes to be passed to or obtained from the RTVD/Riemann routine.
   !! Unassociated pointer means that no operation is required.

   type :: ext_fluxes
      type(fluxpoint), pointer :: li  !< incoming from the left (low index)
      type(fluxpoint), pointer :: lo  !< outgoing on the left
      type(fluxpoint), pointer :: ri  !< incoming from the right (high index)
      type(fluxpoint), pointer :: ro  !< outgoing on the right
   contains
      procedure :: init               !< Nullify point flux pointers
   end type ext_fluxes

contains

!> \brief Nullify point flux pointers

   subroutine init(this)

      implicit none

      class(ext_fluxes), intent(out) :: this  !< object invoking type bound procedure

      this%li => null()
      this%lo => null()
      this%ri => null()
      this%ro => null()

   end subroutine init

end module fluxtypes
