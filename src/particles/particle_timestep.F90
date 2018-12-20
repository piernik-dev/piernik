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
!! \brief Main particle module
!!
!! This module contains all particle related things that are needed by the rest
!! of the code. No other module should directly import anything from particle_*
!<

module particle_timestep
! pulled by NBODY

   implicit none

   private
   public :: timestep_nbody, dt_nbody

   real   :: dt_nbody           !< timestep depends on particles

contains

   subroutine timestep_nbody(dt, cg)

      use constants,      only: ndims, xdim, zdim, big, one, two, zero
      use dataio_pub,     only: msg, printinfo
      use func,           only: operator(.notequals.)
      use grid_cont,      only: grid_container
      use particle_pub,   only: lf_c
      use particle_types, only: pset

      implicit none

      real,                          intent(inout) :: dt
      type(grid_container), pointer, intent(in)    :: cg

      integer                                      :: n_part, i
      integer(kind=4)                              :: cdim
      real                                         :: acc2, max_acc, eta, eps, factor, dt_hydro
      real, dimension(ndims)                       :: max_v

#ifdef VERBOSE
      call printinfo('[particle_timestep:timestep_nbody] Commencing timestep_nbody')
#endif /* VERBOSE */

      n_part = size(pset%p, dim=1)

      eta      = one
      eps      = 1.0e-4
      factor   = one
      dt_nbody = big
      max_acc  = zero

      do i = 1, n_part
         acc2 = zero
         do cdim = xdim, ndims
            acc2 = acc2 + pset%p(i)%acc(cdim)**2
         enddo
         max_acc = max(acc2, max_acc)
      enddo
      max_acc = sqrt(max_acc)

      if(max_acc .notequals. zero) then
         dt_nbody = sqrt(two*eta*eps/max_acc)

         !> \todo following lines are related only to dust particles
         do cdim = xdim, zdim
            max_v(cdim)  = maxval(abs(pset%p(:)%vel(cdim)))
         enddo

         if (any(max_v*dt_nbody > cg%dl)) then
            factor = big
            do cdim = xdim, zdim
               if ((max_v(cdim) .notequals. zero)) then
                  factor = min(cg%dl(cdim)/max_v(cdim), factor)
               endif
            enddo
         endif

         dt_nbody  = lf_c * factor * dt_nbody

      endif

      dt_hydro = dt
      !> \todo verify this condition
      if (dt_nbody .notequals. 0.0) dt = min(dt, dt_nbody)
      write(msg,'(a,3f8.5)') '[particle_timestep:timestep_nbody] dt for hydro, nbody and both: ', dt_hydro, dt_nbody, dt
      call printinfo(msg)

#ifdef VERBOSE
      call printinfo('[particle_timestep:timestep_nbody] Finish timestep_nbody')
#endif /* VERBOSE */

   end subroutine timestep_nbody

end module particle_timestep
