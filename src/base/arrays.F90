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
   real, allocatable, dimension(:,:,:,:), target :: b    !< Main array of magnetic field's components
   real, allocatable, dimension(:,:,:)           :: wa   !< Temporary array used for different purposes, usually has dimension (grid::nx, grid::ny, grid::nz)

#ifdef RESISTIVE
   real, allocatable, dimension(:,:,:)       :: wcu      !< Temporary array used in resistivity module
#endif /* RESISTIVE */

#ifdef GRAV
   real, allocatable, dimension(:,:,:)       :: gpot     !< Array for sum of gravitational potential at t += dt
   real, allocatable, dimension(:,:,:)       :: hgpot    !< Array for sum of gravitational potential at t += 0.5*dt
   real, allocatable, dimension(:,:,:)       :: gp       !< Array for gravitational potential from external fields
   real, allocatable, dimension(:)           :: dprof    !< Array used for storing density during calculation of hydrostatic equilibrium
   real, allocatable, dimension(:)           :: eprof    !< Array used for storing energy during calculation of hydrostatic equilibrium
#ifdef SELF_GRAV
   real, allocatable, dimension(:,:,:)       :: sgp      !< Array for gravitational potential from multigrid or FFT solver
   real, allocatable, dimension(:,:,:)       :: sgpm     !< Array for gravitational potential from multigrid or FFT solver at previous timestep saved by source_terms_grav.
#endif /* SELF_GRAV */
#endif /* GRAV */

#ifdef COSM_RAYS
   real, allocatable, dimension(:,:,:)       :: divvel   !< Array storing \f$\nabla\cdot\mathbf{v}\f$, needed in cosmic ray transport
   real, allocatable, dimension(:,:,:,:)     :: wcr      !< Temporary array used in crdiffusion module
#endif /* COSM_RAYS  */

#ifdef ISO_LOCAL
   real, allocatable, dimension(:,:,:),target:: cs_iso2_arr !< Array storing squared local isothermal sound speed
#endif /* ISO_LOCAL */

   contains

!>
!! Routine that allocates all arrays
!<

   subroutine init_arrays(nvar)

      use diagnostics, only: ma3d, ma4d, my_allocate
      use types,       only: var_numbers
      use grid,        only: cg
#ifdef GRAV
      use diagnostics, only: ma1d
#endif /* GRAV */

      implicit none

      type(var_numbers), intent(in) :: nvar !< fluid database; cannot use fluidindex::nvar here due to circular dependences in some setups

      ma4d = [nvar%all, cg%nx, cg%ny, cg%nz]
      call my_allocate(u, ma4d, "u")
      ma4d = [3, cg%nx, cg%ny, cg%nz]
      call my_allocate(b, ma4d, "b")
      ma3d = [cg%nx, cg%ny, cg%nz]
      call my_allocate(wa, ma3d, "wa")

#ifdef RESISTIVE
      call my_allocate(wcu, ma3d, "wcu")
#endif /* RESISTIVE */

#ifdef GRAV
      call my_allocate(gpot, ma3d, "gpot")
      call my_allocate(hgpot, ma3d, "hgpot")
      call my_allocate(gp, ma3d, "gp")
      ma1d = [cg%nz]
      call my_allocate(dprof, ma1d, "dprof")
      call my_allocate(eprof, ma1d, "eprof")
#ifdef SELF_GRAV
      call my_allocate(sgp, ma3d, "sgp")
      call my_allocate(sgpm, ma3d, "sgpm")
#endif /* SELF_GRAV */
#endif /* GRAV */

#ifdef COSM_RAYS
      call my_allocate(divvel, ma3d, "divvel")
      ma4d = [nvar%crs%all, cg%nx, cg%ny, cg%nz]
      call my_allocate(wcr, ma4d, "wcr")
#endif /* COSM_RAYS  */

#ifdef ISO_LOCAL
      call my_allocate(cs_iso2_arr, ma3d, "cs_iso2_arr")
#endif /* ISO_LOCAL */

   end subroutine init_arrays

!>
!! Routine that deallocates all arrays
!<

   subroutine cleanup_arrays
      implicit none

      if (allocated(u))       deallocate(u)
      if (allocated(b))       deallocate(b)
      if (allocated(wa))      deallocate(wa)

#ifdef RESISTIVE
      if (allocated(wcu))     deallocate(wcu)
#endif /* RESISTIVE */

#ifdef GRAV
      if (allocated(gpot))    deallocate(gpot)
      if (allocated(hgpot))   deallocate(hgpot)
      if (allocated(gp))      deallocate(gp)
      if (allocated(dprof))   deallocate(dprof)
      if (allocated(eprof))   deallocate(eprof)
#ifdef SELF_GRAV
      if (allocated(sgp))     deallocate(sgp)
      if (allocated(sgpm))    deallocate(sgpm)
#endif /* SELF_GRAV */
#endif /* GRAV */

#ifdef COSM_RAYS
      if (allocated(divvel))  deallocate(divvel)
      if (allocated(wcr))     deallocate(wcr)
#endif /* COSM_RAYS */

#ifdef ISO_LOCAL
      if (allocated(cs_iso2_arr)) deallocate(cs_iso2_arr)
#endif /* ISO_LOCAL */

   end subroutine cleanup_arrays

end module arrays
