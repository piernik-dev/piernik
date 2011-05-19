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
!! Module storing all global %arrays used in simulation
!<
module arrays

   implicit none

   public  ! QA_WARN nothing to hide here

   real, allocatable, dimension(:,:,:,:), target :: u    !< Main array of all fluids' components
   real, allocatable, dimension(:,:,:,:), target :: uh   !< Main array of all fluids' components (for t += dt/2)
   real, allocatable, dimension(:,:,:,:), target :: b    !< Main array of magnetic field's components
   real, allocatable, dimension(:,:,:,:), target :: u0   !< Copy of main array of all fluids' components
   real, allocatable, dimension(:,:,:,:), target :: b0   !< Copy of main array of magnetic field's components
   real, allocatable, dimension(:,:,:),   target :: wa   !< Temporary array used for different purposes, usually has dimension (grid::nx, grid::ny, grid::nz)

#ifdef GRAV
   real, allocatable, dimension(:), target     :: dprof    !< Array used for storing density during calculation of hydrostatic equilibrium
   real, allocatable, dimension(:), target     :: eprof    !< Array used for storing energy during calculation of hydrostatic equilibrium
#endif /* GRAV */

contains

!>
!! Routine that allocates all arrays
!<

   subroutine init_arrays(flind)

      use constants,   only: PIERNIK_INIT_BASE, ndims
      use diagnostics, only: ma3d, ma4d, my_allocate
      use dataio_pub,  only: die, code_progress
      use fluidtypes,  only: var_numbers
      use grid,        only: cg
      use mpisetup,    only: repeat_step
#ifdef GRAV
      use diagnostics, only: ma1d
#endif /* GRAV */

      implicit none

      type(var_numbers), intent(in) :: flind !< fluid database; cannot use fluidindex::flind here due to circular dependencies in some setups

      if (code_progress < PIERNIK_INIT_BASE) call die("[arrays:init_arrays] grid or fluids not initialized.")

      ma4d = [flind%all, cg%nx, cg%ny, cg%nz]
      call my_allocate(u, ma4d, "u")
!     call cg%u%init(flind%all, cg%nx, cg%ny, cg%nz) ! \todo use after u -> cg%u transition
      call cg%u%init(u)

      if (repeat_step) call my_allocate(u0, ma4d, "u0")
      call my_allocate(uh, ma4d, "uh")

      ma4d = [ndims, cg%nx, cg%ny, cg%nz]
      call my_allocate(b, ma4d, "b")
!     call cg%b%init(ndims, cg%nx, cg%ny, cg%nz)  ! \todo use after b -> cg%b transition
      call cg%b%init(b)

      if (repeat_step) call my_allocate(b0, ma4d, "b0")

      ma3d = [cg%nx, cg%ny, cg%nz]
      call my_allocate(wa, ma3d, "wa")
!     call cg%wa%init(cg%nx, cg%ny, cg%nz)        ! \todo use after wa -> cg%wa transition
      call cg%wa%init(wa)

#ifdef GRAV
      ma1d = [cg%nz]
      call my_allocate(dprof, ma1d, "dprof")
      call my_allocate(eprof, ma1d, "eprof")
#endif /* GRAV */

   end subroutine init_arrays

!>
!! Routine that deallocates all arrays
!<

   subroutine cleanup_arrays
      use grid, only: cg
      implicit none

!      if (allocated(u))       deallocate(u)
      call cg%u%clean()
!      if (allocated(b))       deallocate(b)
      call cg%b%clean()
      if (allocated(u0))      deallocate(u0)
      if (allocated(uh))      deallocate(uh)
      if (allocated(b0))      deallocate(b0)
!      if (allocated(wa))      deallocate(wa)
      call cg%wa%clean()

#ifdef GRAV
      if (allocated(dprof))   deallocate(dprof)
      if (allocated(eprof))   deallocate(eprof)
#endif /* GRAV */

   end subroutine cleanup_arrays

end module arrays
