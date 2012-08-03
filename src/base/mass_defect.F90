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
#include "macros.h"

!> \brief This module tracks inaccuracies in total mass present on the grid

module mass_defect

   implicit none

   private
   public :: magic_mass, local_magic_mass, recent_magic_mass, init_magic_mass, update_magic_mass, cleanup_magic_mass

   real, dimension(:),   allocatable       :: magic_mass
   real, dimension(:),   allocatable, save :: local_magic_mass, recent_magic_mass

contains

   subroutine init_magic_mass

      use fluidindex, only: flind
      use mpisetup,   only: master

      implicit none

      allocate(local_magic_mass(flind%fluids), recent_magic_mass(flind%fluids))
      local_magic_mass  = 0.0
      recent_magic_mass = 0.0
      if (master) then
         allocate(magic_mass(flind%fluids))
         magic_mass = 0.0
         ! this variable should not be used on slaves
      endif

   end subroutine init_magic_mass

   subroutine update_magic_mass(tsl)

      use fluidindex,  only: flind
      use mpi,         only: MPI_DOUBLE_PRECISION, MPI_SUM
      use mpisetup,    only: comm, mpi_err, FIRST, master

      implicit none

      logical, optional             :: tsl
      real, dimension(flind%fluids) :: magic_mass0
      integer(kind=4)               :: ifl, pos

      call MPI_Reduce(local_magic_mass, magic_mass0, int(flind%fluids, kind=4), MPI_DOUBLE_PRECISION, MPI_SUM, FIRST, comm, mpi_err)
      local_magic_mass(:) = 0.0

      if (master) then
         magic_mass = magic_mass + magic_mass0

         if (present(tsl)) then
            do ifl = 1, flind%fluids
               pos = flind%all_fluids(ifl)%fl%pos
               flind%all_fluids(ifl)%fl%snap%mmass_cum = magic_mass(pos)
               flind%all_fluids(ifl)%fl%snap%mmass_cur = magic_mass(pos) - recent_magic_mass(pos)
            enddo
            recent_magic_mass(:) = magic_mass(:)
         endif
      endif

   end subroutine update_magic_mass

   subroutine cleanup_magic_mass

      implicit none

      if (allocated(magic_mass)) deallocate(magic_mass)
      deallocate(local_magic_mass, recent_magic_mass)

   end subroutine cleanup_magic_mass

end module mass_defect
