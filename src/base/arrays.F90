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

   real, allocatable, dimension(:,:,:,:)     :: u        !< Main array of all fluids' componets
   real, allocatable, dimension(:,:,:,:)     :: b        !< Main array of magnetic field's components
   real, allocatable, dimension(:,:,:)       :: wa       !< Temporary array used for different purposes, usually has dimension (grid::nx, grid::ny, grid::nz)
   real, allocatable, dimension(:,:,:)       :: wcu      !< Temporary array used in resistivity module
#ifdef GRAV
   real, allocatable, dimension(:,:,:)       :: gp       !< Array for gravitational potential
   real, allocatable, dimension(:)           :: dprof    !< Array used for storing density during calculation of hydrostatic equilibrium
   real, allocatable, dimension(:)           :: eprof    !< Array used for storing energy during calculation of hydrostatic equilibrium
#endif /* GRAV */

#ifdef COSM_RAYS
   real, allocatable, dimension(:,:,:)       :: divvel   !< Array storing \f$\nabla\cdot\mathbf{v}\f$, needed in cosmic ray transport

#endif /* COSM_RAYS  */

   real,    allocatable, dimension(:)        :: rlscal   !< Arrays to store additional scalr quantities in hdf4 files (real)
   integer, allocatable, dimension(:)        :: intscal  !< Arrays to store additional scalr quantities in hdf4 files (int)

   contains
!>
!! Routine that allocates all arrays
!<
   subroutine arrays_allocate(nx,ny,nz,nvar)
      implicit none
      integer, intent(in) :: nx,ny,nz,nvar

      if(.not.allocated(rlscal))  allocate(rlscal(nrlscal))
      if(.not.allocated(intscal)) allocate(intscal(nintscal))

      if(.not.allocated(u))       allocate(u(nvar,nx,ny,nz))
      if(.not.allocated(b))       allocate(b(3,nx,ny,nz))

#ifdef GRAV
      if(.not.allocated(gp))      allocate(gp(nx,ny,nz))
      if(.not.allocated(dprof))   allocate(dprof(nz))
      if(.not.allocated(eprof))   allocate(eprof(nz))
#endif /* GRAV */

      if(.not.allocated(wa))      allocate(wa(nx,ny,nz))
#ifdef RESISTIVE
      if(.not.allocated(wcu))     allocate(wcu(nx,ny,nz))
#endif /* RESISTIVE */
#ifdef COSM_RAYS
      if(.not.allocated(divvel)) allocate(divvel(nx,ny,nz))
#endif /* COSM_RAYS  */

   end subroutine arrays_allocate

!>
!! Routine that deallocates all arrays
!<
   subroutine arrays_deallocate

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
      if(allocated(dprof))   deallocate(dprof)
      if(allocated(eprof))   deallocate(eprof)
#endif /* GRAV */

#ifdef COSM_RAYS
      if(allocated(divvel))  deallocate(divvel)
#endif /* COSM_RAYS */

   end subroutine arrays_deallocate

end module arrays
