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

#ifdef DUST_PARTICLES
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
         if (pset%pdata%phy) then
            v_tmp = abs(pset%pdata%vel)
            where (v_tmp > max_v) max_v = v_tmp
         endif
         pset => pset%nxt
      enddo
      factor = maxval(max_v / cg%dl)

   end function max_pvel_dl
#endif /* DUST_PARTICLES */

   subroutine max_pacc_3d(cg, factor_a)

      use constants,      only: CENTER, half, LO, ndims, xdim, zdim, one, zero
      use grid_cont,      only: grid_container
      use mpisetup,       only: proc
      use particle_types, only: particle

      implicit none

      type(grid_container), pointer, intent(in)    :: cg
      real,                          intent(inout) :: factor_a
      type(particle), pointer                      :: pset
      real                                         :: acc2, eta, max_acc, new_fa
      real, dimension(ndims)                       :: max_crd
      integer(kind=4)                              :: cdim

      max_acc = -one
      max_crd = huge(1.)
      pset => cg%pset%first
      do while (associated(pset))
         if (pset%pdata%phy) then
            acc2 = sum(pset%pdata%acc(:)**2)
            if (acc2 > max_acc) then
               max_acc = acc2
               max_crd = pset%pdata%pos(:)
            endif
         endif
         pset => pset%nxt
      enddo

      if (max_acc < zero) return

      eta = minval(cg%dl)   !scale timestep with cell size
      max_acc = sqrt(max_acc)
      new_fa = max_acc / eta
      if (new_fa > factor_a) then
         factor_a        = new_fa
         pacc_max%val    = max_acc
         pacc_max%coords = max_crd
         pacc_max%proc   = proc
         do cdim = xdim, zdim
            pacc_max%loc(cdim) = int( half + (pacc_max%coords(cdim) - cg%coord(CENTER,cdim)%r(cg%ijkse(cdim, LO))) * cg%idl(cdim) )
         enddo
      endif

   end subroutine max_pacc_3d

!> \todo this is a part of cg_list_dataop::get_extremum; consider merge it together
   subroutine locate_extremum(prop)

      use constants,    only: I_ONE, ndims
      use MPIF,         only: MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_2DOUBLE_PRECISION, &
           &                  MPI_STATUS_IGNORE, MPI_COMM_WORLD, MPI_MAXLOC, MPI_IN_PLACE
      use MPIFUN,       only: MPI_Allreduce, MPI_Send, MPI_Recv
      use mpisetup,     only: err_mpi, master, proc, FIRST
      use types,        only: value

      implicit none

      type(value),   intent(out) :: prop    !< precise location of the extremum to be found
      integer(kind=4), parameter :: tag1 = 11, tag2 = tag1 + 1, tag3 = tag2 + 1
      enum, bind(C)
         enumerator :: I_V, I_P !< value and proc
      end enum
      real, dimension(I_V:I_P)   :: v_red

      v_red(I_V) = prop%val
      v_red(I_P) = real(proc)

      call MPI_Allreduce(MPI_IN_PLACE, v_red, I_ONE, MPI_2DOUBLE_PRECISION, MPI_MAXLOC, MPI_COMM_WORLD, err_mpi)

      prop%val = v_red(I_V)
      prop%proc = int(v_red(I_P), kind=4)

      if (prop%proc /= 0) then
         if (proc == prop%proc) then ! slave
            call MPI_Send (prop%loc,    ndims, MPI_INTEGER,          FIRST, tag1, MPI_COMM_WORLD, err_mpi)
            call MPI_Send (prop%coords, ndims, MPI_DOUBLE_PRECISION, FIRST, tag2, MPI_COMM_WORLD, err_mpi)
            call MPI_Send (prop%assoc,  I_ONE, MPI_DOUBLE_PRECISION, FIRST, tag3, MPI_COMM_WORLD, err_mpi)
         endif
         if (master) then
            call MPI_Recv (prop%loc,    ndims, MPI_INTEGER,          prop%proc, tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, err_mpi)
            call MPI_Recv (prop%coords, ndims, MPI_DOUBLE_PRECISION, prop%proc, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE, err_mpi)
            call MPI_Recv (prop%assoc,  I_ONE, MPI_DOUBLE_PRECISION, prop%proc, tag3, MPI_COMM_WORLD, MPI_STATUS_IGNORE, err_mpi)
         endif
      endif

   end subroutine locate_extremum

   subroutine timestep_nbody(dt)

      use cg_leaves,    only: leaves
      use cg_list,      only: cg_list_element
      use constants,    only: zero, one, two, pMIN, V_VERBOSE
      use dataio_pub,   only: msg, printinfo
      use grid_cont,    only: grid_container
      use mpisetup,     only: piernik_MPI_Allreduce, master
      use particle_pub, only: lf_c, eps, ignore_dt_fluid

      implicit none

      real,            intent(inout) :: dt
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
      real                           :: factor_a, factor_v, dt_hydro

#ifdef VERBOSE
      call printinfo('[particle_timestep:timestep_nbody] Commencing timestep_nbody')
#endif /* VERBOSE */

      factor_a = -one
      factor_v = -one

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         call max_pacc_3d(cg, factor_a)

#ifdef DUST_PARTICLES
         factor_v = max(factor_v, max_pvel_dl(cg))
#endif /* DUST_PARTICLES */

         cgl => cgl%nxt
      enddo

      if (factor_a > zero) then
         dt_nbody = lf_c * sqrt(two * eps / factor_a)
      else
         dt_nbody = huge(1.)
      endif
      pacc_max%assoc = dt_nbody

#ifdef DUST_PARTICLES
      if (factor_v * dt_nbody > one)) dt_nbody = dt_nbody / factor_v
#endif /* DUST_PARTICLES */

      call locate_extremum(pacc_max)
      call piernik_MPI_Allreduce(dt_nbody, pMIN)
      dt_hydro = dt

      if (ignore_dt_fluid) then
         dt = dt_nbody      !IGNORE HYDRO TIMESTEP
      else
         !> \todo verify this condition
         if (abs(dt_nbody) > tiny(1.)) dt = min(dt, dt_nbody)
      endif

      write(msg,'(a,3g12.5)') '[particle_timestep:timestep_nbody] dt for hydro, nbody and both: ', dt_hydro, dt_nbody, dt
      if (master) call printinfo(msg, V_VERBOSE)

#ifdef VERBOSE
      call printinfo('[particle_timestep:timestep_nbody] Finish timestep_nbody')
#endif /* VERBOSE */

   end subroutine timestep_nbody

end module particle_timestep
