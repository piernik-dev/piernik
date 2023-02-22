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

!> \brief This module contains fluxpoint type, useful for flux enforcing and exchange in 1D solvers, like current directional-split (M)HD ones.

module flx_cell

   implicit none

   private
   public  :: fluxpoint

   !> \brief Structure that contains u-flux at a single cell index.

   type :: fluxpoint
      real, dimension(:), allocatable :: uflx   !< u-flux
      real, dimension(:), allocatable :: bflx   !< (b,psi)-flux
      integer                         :: index  !< Index where the flux has to be applied
   contains
      procedure :: init     !< Allocate flux vector
      procedure :: cleanup  !< Deallocate flux vector
   end type fluxpoint

contains

!> \brief Allocate flux vector

   subroutine init(this)

      use constants,  only: psidim, has_B
      use dataio_pub, only: die
      use fluidindex, only: flind

      implicit none

      class(fluxpoint), intent(inout) :: this  !< object invoking type bound procedure

      if (allocated(this%uflx)) call die("[flx_cell:init] uflx already allocated")
      allocate(this%uflx(flind%all))
      if (has_B) allocate(this%bflx(psidim))

   end subroutine init

!> \brief Deallocate flux vector

   subroutine cleanup(this)

      implicit none

      class(fluxpoint), intent(inout) :: this  !< object invoking type bound procedure

      if (allocated(this%uflx)) deallocate(this%uflx)
      if (allocated(this%bflx)) deallocate(this%bflx)

   end subroutine cleanup

end module flx_cell
