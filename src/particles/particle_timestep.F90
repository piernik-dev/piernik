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
!! \brief Timestep particle module
!!
!! This module contains all particle related routines to determine dt limitations.
!<

module particle_timestep
! pulled by NBODY

   use types, only: value

   implicit none

   private
   public :: timestep_nbody, dt_nbody, pacc_max

   real        :: dt_nbody           !< timestep depends on particles
   type(value) :: pacc_max

contains

   real function max_pvel_dl(cg) result (factor)

      use constants,      only: ndims, zero
      use grid_cont,      only: grid_container
      use particle_types, only: particle

      implicit none

      type(grid_container), pointer, intent(in)  :: cg
      type(particle), pointer                    :: pset
      real, dimension(ndims)                     :: max_v, v_tmp

      max_v = zero
      pset => cg%pset%first
      do while (associated(pset))
         v_tmp = abs(pset%pdata%vel)
         where (v_tmp > max_v) max_v = v_tmp
         pset => pset%nxt
      enddo
      factor = maxval(max_v / cg%dl)

   end function max_pvel_dl

   subroutine max_pacc_3d(cg, max_pacc)

      use constants,      only: CENTER, half, LO, xdim, zdim, zero
      use grid_cont,      only: grid_container
      use mpisetup,       only: proc
      use types,          only: value
      use particle_types, only: particle

      implicit none

      type(grid_container), pointer, intent(in)  :: cg
      type(value),                   intent(out) :: max_pacc
      type(particle), pointer                    :: pset
      real                                       :: acc2, max_acc
      integer(kind=4)                            :: cdim

      max_pacc%assoc = huge(1.)

      max_acc  = zero
      pset => cg%pset%first
      do while (associated(pset))
         acc2 = sum(pset%pdata%acc(:)**2)
         if (acc2 > max_acc) then
            max_acc = acc2
            max_pacc%coords(:) = pset%pdata%pos(:)
            !max_pacc%proc = i !> \todo it might be an information about extremum particle, but the scheme of log file is to print the process number
         endif
         pset => pset%nxt
      enddo
      max_pacc%val = sqrt(max_acc)
      do cdim = xdim, zdim
         max_pacc%loc(cdim) = int( half + (max_pacc%coords(cdim) - cg%coord(CENTER,cdim)%r(cg%ijkse(cdim, LO))) * cg%idl(cdim) )
      enddo
      max_pacc%proc = proc

   end subroutine max_pacc_3d

   subroutine timestep_nbody(dt)

      use cg_leaves,    only: leaves
      use cg_list,      only: cg_list_element
      use constants,    only: zero, two, pMIN
      use dataio_pub,   only: msg, printinfo
      use global,       only: dt_max
      use grid_cont,    only: grid_container
      use mpisetup,     only: piernik_MPI_Allreduce, master
      use particle_pub, only: lf_c, ignore_dt_fluid
#ifdef DUST_PARTICLES
      use constants,    only: one
#endif /* DUST_PARTICLES */

      implicit none

      real,            intent(inout) :: dt
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
      real                           :: eta, eps, factor_a, factor_v, dt_hydro

#ifdef VERBOSE
      call printinfo('[particle_timestep:timestep_nbody] Commencing timestep_nbody')
#endif /* VERBOSE */

      eps      = 1.0e-1
      factor_a = zero
      factor_v = zero

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         eta = minval(cg%dl)   !scale timestep with cell size

         call max_pacc_3d(cg, pacc_max)
         factor_a = max(factor_a, pacc_max%val / eta)

#ifdef DUST_PARTICLES
         factor_v = max(factor_v, max_pvel_dl(cg))
#endif /* DUST_PARTICLES */

         cgl => cgl%nxt
      enddo

      dt_nbody = lf_c * sqrt(two * eps / factor_a)

#ifdef DUST_PARTICLES
      if (factor_v * dt_nbody > one)) dt_nbody = dt_nbody / factor_v
#endif /* DUST_PARTICLES */

      dt_nbody = min(dt_nbody, dt_max)
      call piernik_MPI_Allreduce(dt_nbody, pMIN)
      pacc_max%assoc = dt_nbody
      dt_hydro = dt

      if (ignore_dt_fluid) then
         dt = dt_nbody      !IGNORE HYDRO TIMESTEP
      else
         !> \todo verify this condition
         if (abs(dt_nbody) > tiny(1.)) dt = min(dt, dt_nbody)
      endif

      write(msg,'(a,3g12.5)') '[particle_timestep:timestep_nbody] dt for hydro, nbody and both: ', dt_hydro, dt_nbody, dt
      if (master) call printinfo(msg)

#ifdef VERBOSE
      call printinfo('[particle_timestep:timestep_nbody] Finish timestep_nbody')
#endif /* VERBOSE */

   end subroutine timestep_nbody

end module particle_timestep
