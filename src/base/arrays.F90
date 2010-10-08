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
!>
!! Module storing all global %arrays used in simulation
!<
module arrays

   implicit none
   integer, parameter  :: nrlscal=100   !< Size of arrays::rlscal array
   integer, parameter  :: nintscal=100  !< Size of arrays::intscal array

   real, allocatable, dimension(:,:,:,:), target :: u    !< //Main array of all fluids' componets
   real, allocatable, dimension(:,:,:,:), target :: b    !< //Main array of magnetic field's components
   real, allocatable, dimension(:,:,:)           :: wa   !< Temporary array used for different purposes, usually has dimension (grid::nx, grid::ny, grid::nz)
   real, allocatable, dimension(:,:,:)       :: wcu      !< Temporary array used in resistivity module
#ifdef GRAV
   real, allocatable, dimension(:,:,:)       :: gpot     !< Array for sum of gravitational potential at t += dt
   real, allocatable, dimension(:,:,:)       :: hgpot    !< Array for sum of gravitational potential at t += 0.5*dt
   real, allocatable, dimension(:,:,:)       :: gp       !< Array for gravitational potential from external fields
   real, allocatable, dimension(:)           :: dprof    !< Array used for storing density during calculation of hydrostatic equilibrium
   real, allocatable, dimension(:)           :: eprof    !< Array used for storing energy during calculation of hydrostatic equilibrium
#if defined(MULTIGRID) || defined(POISSON_FFT)
   real, allocatable, dimension(:,:,:)       :: sgp      !< Array for gravitational potential from multigrid or FFT solver
   real, allocatable, dimension(:,:,:)       :: sgpm     !< Array for gravitational potential from multigrid or FFT solver at previous timestep saved by source_terms_grav.
#endif
#endif /* GRAV */

#ifdef COSM_RAYS
   real, allocatable, dimension(:,:,:)       :: divvel   !< Array storing \f$\nabla\cdot\mathbf{v}\f$, needed in cosmic ray transport
   real, allocatable, dimension(:,:,:,:)     :: wcr      !< Temporary array used in crdiffusion module
#endif /* COSM_RAYS  */
#ifdef ISO_LOCAL
   real, allocatable, dimension(:,:,:),target:: cs_iso2_arr !< Array storing squared local isothermal sound speed
#endif /* ISO_LOCAL */

   real,    allocatable, dimension(:)        :: rlscal   !< Arrays to store additional scalr quantities in hdf4 files (real)
   integer, allocatable, dimension(:)        :: intscal  !< Arrays to store additional scalr quantities in hdf4 files (int)

   contains
!>
!! Routine that allocates all arrays
!<

   subroutine init_arrays(nx,ny,nz,nvar)

      use diagnostics, only: my_allocate
      use fluidtypes,  only: var_numbers

      implicit none

      type(var_numbers), intent(in) :: nvar
      integer, intent(in) :: nx,ny,nz

      if(.not.allocated(rlscal))  allocate(rlscal(nrlscal))
      if(.not.allocated(intscal)) allocate(intscal(nintscal))

      call my_allocate(u, [nvar%all,nx,ny,nz],"u")
      call my_allocate(b, [3,nx,ny,nz], "b")

#ifdef GRAV
      if(.not.allocated(gp))      allocate(gp(nx,ny,nz))
      if(.not.allocated(gpot))    allocate(gpot(nx,ny,nz))
      if(.not.allocated(hgpot))   allocate(hgpot(nx,ny,nz))
#if defined(MULTIGRID) || defined(POISSON_FFT)
      if(.not.allocated(sgp))     allocate(sgp(nx,ny,nz))
      if(.not.allocated(sgpm))    allocate(sgpm(nx,ny,nz))
#endif
      if(.not.allocated(dprof))   allocate(dprof(nz))
      if(.not.allocated(eprof))   allocate(eprof(nz))
#endif /* GRAV */

      if(.not.allocated(wa))      allocate(wa(nx,ny,nz))
#ifdef RESISTIVE
      if(.not.allocated(wcu))     allocate(wcu(nx,ny,nz))
#endif /* RESISTIVE */
#ifdef COSM_RAYS
      if(.not.allocated(divvel)) allocate(divvel(nx,ny,nz))
      if(.not.allocated(wcr))    allocate(wcr(nvar%crs%all,nx,ny,nz))
#endif /* COSM_RAYS  */
#ifdef ISO_LOCAL
      if(.not.allocated(cs_iso2_arr)) allocate(cs_iso2_arr(nx,ny,nz))
#endif /* ISO_LOCAL */
   end subroutine init_arrays

!>
!! Routine that deallocates all arrays
!<
   subroutine cleanup_arrays

      if(allocated(rlscal))  deallocate(rlscal)
      if(allocated(intscal)) deallocate(intscal)
      if(allocated(u))       deallocate(u)
      if(allocated(b))       deallocate(b)

      if(allocated(wa))      deallocate(wa)
#ifdef RESISTIVE
      if(allocated(wcu))     deallocate(wcu)
#endif /* RESISTIVE */

#ifdef GRAV
      if(allocated(gp))      deallocate(gp)
      if(allocated(gpot))    deallocate(gpot)
      if(allocated(hgpot))   deallocate(hgpot)
#if defined(MULTIGRID) || defined(POISSON_FFT)
      if(allocated(sgp))     deallocate(sgp)
      if(allocated(sgpm))    deallocate(sgpm)
#endif
      if(allocated(dprof))   deallocate(dprof)
      if(allocated(eprof))   deallocate(eprof)
#endif /* GRAV */

#ifdef COSM_RAYS
      if(allocated(divvel))  deallocate(divvel)
      if(allocated(wcr))     deallocate(wcr)
#endif /* COSM_RAYS */
#ifdef ISO_LOCAL
      if(allocated(cs_iso2_arr)) deallocate(cs_iso2_arr)
#endif /* ISO_LOCAL */

   end subroutine cleanup_arrays

end module arrays
