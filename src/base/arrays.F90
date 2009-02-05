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

module arrays

   implicit none
   integer,parameter  :: nrlscal=100, nintscal=100

   real, allocatable, dimension(:,:,:,:)     :: u, b
   real, allocatable, dimension(:,:,:)       :: wa, wcu
#ifdef GRAV
   real, allocatable, dimension(:,:,:)       :: gp
   real, allocatable, dimension(:)           :: dprof, eprof
#endif /* GRAV */

   real(kind=4), allocatable, dimension(:,:,:)  :: outwa, outwb, outwc

#ifdef COSM_RAYS
   real, allocatable, dimension(:,:,:)       :: divvel
#endif /* COSM_RAYS  */

   real,    allocatable, dimension(:)        :: rlscal
   integer, allocatable, dimension(:)        :: intscal
   real,    allocatable, dimension(:,:)      :: ul0, ur0
   real,    allocatable, dimension(:,:,:,:)  :: bndxrar, bndyrar

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
      allocate(outwa(nx,ny,nz),outwb(nx,ny,nz),outwc(nx,ny,nz))

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

      if(allocated(bndxrar)) deallocate(bndxrar)
      if(allocated(bndyrar)) deallocate(bndyrar)
      deallocate(outwa,outwb,outwc)

   end subroutine arrays_deallocate

end module arrays
