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
!! \brief Diagnostics of particle module
!!
!! This module contains particle diagnostics related things that are needed by I/O and debug routines.
!<

module particle_utils
! pulled by NBODY

   implicit none

   private
   public :: max_pvel_1d, max_pacc_3d, particle_diagnostics, twodtscheme, dump_diagnose, tot_energy, d_energy, tot_angmom, d_angmom

   real    :: tot_angmom           !< angular momentum of set of the particles
   real    :: tot_energy           !< total energy of set of the particles
   real    :: d_energy             !< error of energy of set of the particles in succeeding timesteps
   real    :: d_angmom             !< error of angular momentum in succeeding timensteps
   logical :: twodtscheme
   logical :: dump_diagnose        !< dump diagnose for each particle to a seperate log file

contains

   subroutine max_pvel_1d(max_v)

      use constants,      only: ndims, xdim, zdim
      use particle_types, only: pset

      implicit none

      real, dimension(ndims), intent(out) :: max_v
      integer(kind=4)                     :: cdim

      do cdim = xdim, zdim
         max_v(cdim)  = maxval(abs(pset%p(:)%vel(cdim)))
      enddo

   end subroutine max_pvel_1d

   subroutine max_pacc_3d(cg, max_pacc)

      use constants,      only: big, CENTER, half, LO, xdim, zdim, zero
      use grid_cont,      only: grid_container
      use mpisetup,       only: proc
      use particle_types, only: pset
      use types,          only: value

      implicit none

      type(grid_container), pointer, intent(in)  :: cg
      type(value),                   intent(out) :: max_pacc
      real                                       :: acc2, max_acc
      integer                                    :: i, n_part
      integer(kind=4) :: cdim

      max_pacc%assoc = big

      max_acc  = zero
      n_part = size(pset%p, dim=1)

      do i = 1, n_part
         acc2 = sum(pset%p(i)%acc(:)**2)
         if (acc2 > max_acc) then
            max_acc = acc2
            max_pacc%coords(:) = pset%p(i)%pos(:)
            !max_pacc%proc = i !> \todo it might be an information about extremum particle, but the scheme of log file is to print the process number
         endif
      enddo
      max_pacc%val = sqrt(max_acc)
      do cdim = xdim, zdim
         max_pacc%loc(cdim) = int( half + (max_pacc%coords(cdim) - cg%coord(CENTER,cdim)%r(cg%ijkse(cdim, LO))) * cg%idl(cdim) )
      enddo
      max_pacc%proc = proc

   end subroutine max_pacc_3d

   subroutine particle_diagnostics(tg,kdt)

#ifdef VERBOSE
      use dataio_pub, only: msg, printinfo
#endif /* VERBOSE */

      implicit none

      real, intent(in) :: tg, kdt
      real, save       :: init_energy          !< total initial energy of set of the particles
      real, save       :: init_angmom          !< initial angular momentum of set of the particles
      logical, save    :: first_run_lf = .true.

      call get_angmom_totener(tot_angmom, tot_energy)

      if (first_run_lf) then
         init_energy = tot_energy
         init_angmom = tot_angmom
         first_run_lf = .false.
      endif
      d_energy = log_error(tot_energy, init_energy)
      d_angmom = log_error(tot_angmom, init_angmom)

#ifdef VERBOSE
      write(msg,'(a,3(1x,e12.5))') '[particle_utils:particle_diagnostics] Total energy: initial, current, error ', init_energy, tot_energy, d_energy
      call printinfo(msg)
      write(msg,'(a,3(1x,e12.5))') '[particle_utils:particle_diagnostics] ang_momentum: initial, current, error ', init_angmom, tot_angmom, d_angmom
      call printinfo(msg)
#endif /* VERBOSE */

      if (dump_diagnose) call dump_particles_to_textfile(tg,kdt)

   end subroutine particle_diagnostics

   real function log_error(vcurr,vinit)

      use constants, only: zero
      use func,      only: operator(.equals.)

      implicit none

      real, intent(in) :: vcurr, vinit

      if (vinit .equals. zero) then
         log_error = vcurr
         return
      endif
      log_error = (vcurr - vinit) / vinit
      if (log_error .equals. zero) then
         log_error = zero
      else
         log_error = log(abs(log_error))
      endif

   end function log_error

   subroutine get_angmom_totener(ang_momentum, total_energy)

      use constants,      only: xdim, ydim, zdim, zero
      use particle_types, only: pset

      implicit none

      real, intent(out) :: ang_momentum
      real, intent(out) :: total_energy !< total energy of set of particles
      integer           :: i
      real              :: L1, L2, L3

      ang_momentum = zero
      total_energy = zero

      do i = 1, size(pset%p, dim=1)
         L1 = pset%p(i)%pos(ydim) * pset%p(i)%vel(zdim) - pset%p(i)%pos(zdim) * pset%p(i)%vel(ydim)
         L2 = pset%p(i)%pos(zdim) * pset%p(i)%vel(xdim) - pset%p(i)%pos(xdim) * pset%p(i)%vel(zdim)
         L3 = pset%p(i)%pos(xdim) * pset%p(i)%vel(ydim) - pset%p(i)%pos(ydim) * pset%p(i)%vel(xdim)
         ang_momentum = ang_momentum + pset%p(i)%mass * sqrt(L1**2 + L2**2 + L3**2)
         total_energy = total_energy + pset%p(i)%energy
      enddo

   end subroutine get_angmom_totener

   subroutine dump_particles_to_textfile(tg,kdt)

      use constants,        only: ndims
      use particle_gravity, only: get_acc_model
      use particle_types,   only: pset

      implicit none

      real, intent(in)       :: tg, kdt
      real, dimension(ndims) :: acc2
      integer                :: i, lun_out

      open(newunit=lun_out, file='leapfrog_out.log', status='unknown',  position='append')
         do i = 1, size(pset%p, dim=1)
            call get_acc_model(i, 0.0, acc2)
            write(lun_out, '(a,I3.3,1X,19(E13.6,1X))') 'particle', i, tg+kdt, kdt, pset%p(i)%mass, pset%p(i)%pos, pset%p(i)%vel, pset%p(i)%acc, acc2(:), pset%p(i)%energy
         enddo
      close(lun_out)

   end subroutine dump_particles_to_textfile

end module particle_utils
