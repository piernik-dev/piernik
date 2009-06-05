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
!! Module storing all global arrays used in simulation
!<
module arrays

   implicit none
   integer,parameter  :: nrlscal=100, nintscal=100

   real, allocatable, dimension(:,:,:,:)     :: u        !< Main array of all fluids' componets
   real, allocatable, dimension(:,:,:,:)     :: b        !< Main array of magnetic field's components
   real, allocatable, dimension(:,:,:)       :: wa       !< Temporary array used for different purposes, usually has dimension (nx, ny, nz)
   real, allocatable, dimension(:,:,:)       :: wcu      !< Temporary array used in resistivity module
#ifdef GRAV
   real, allocatable, dimension(:,:,:)       :: gp       !< Array for gravitational potential
   real, allocatable, dimension(:)           :: dprof    !< Array used for storing density during calculation of hydrostatic equilibrium
   real, allocatable, dimension(:)           :: eprof    !< Array used for storing energy during calculation of hydrostatic equilibrium
#endif /* GRAV */

#ifdef COSM_RAYS
   real, allocatable, dimension(:,:,:)       :: divvel   !< Array storing \f$\nabla\cdot\mathbf{v}\f$, needed in comic ray transport

#endif /* COSM_RAYS  */

   real,    allocatable, dimension(:)        :: rlscal   !< Arrays to store additional scalr quantities in hdf4 files (real)
   integer, allocatable, dimension(:)        :: intscal  !< Arrays to store additional scalr quantities in hdf4 files (int)

   contains

   subroutine arrays_allocate(nx,ny,nz,nvar)
      implicit none
      integer, intent(in) :: nx,ny,nz,nvar

      allocate(rlscal(nrlscal),intscal(nintscal)) ! arrays to store additional
                                                  ! scalar quantities in restart files
      allocate(u(nvar,nx,ny,nz),b(3,nx,ny,nz))

#ifdef GRAV
      allocate(gp(nx,ny,nz))
      allocate(dprof(nz),eprof(nz))
#endif /* GRAV */

      allocate(wa(nx,ny,nz),wcu(nx,ny,nz))

#ifdef COSM_RAYS
      allocate(divvel(nx,ny,nz))
#endif /* COSM_RAYS  */

   end subroutine arrays_allocate

   subroutine arrays_deallocate

      deallocate(rlscal,intscal)
      deallocate(u,b)
      deallocate(wa,wcu)

#ifdef GRAV
      deallocate(gp)
      deallocate(dprof,eprof)
#endif /* GRAV */

#ifdef COSM_RAYS
      deallocate(divvel)
#endif /* COSM_RAYS */

   end subroutine arrays_deallocate

end module arrays
