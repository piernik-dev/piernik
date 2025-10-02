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

module particle_diag
! pulled by NBODY

   implicit none

   private
   public :: particle_diagnostics, print_all_particles, tot_energy, d_energy, tot_angmom, d_angmom

   real :: tot_angmom           !< angular momentum of set of the particles
   real :: tot_energy           !< total energy of set of the particles
   real :: d_energy             !< error of energy of set of the particles in succeeding timesteps
   real :: d_angmom             !< error of angular momentum in succeeding timesteps

contains

   subroutine print_all_particles

      use cg_leaves,  only: leaves
      use cg_list,    only: cg_list_element
      use constants,  only: msg_len
      use dataio_pub, only: msg, printinfo
      use mpisetup,   only: master

      implicit none

      type(cg_list_element), pointer :: cgl

      character(len=msg_len) :: rem

      if (master) then
         call printinfo("[particle_diag:print_all_particles] Known particles:")
         write(msg, '(a,a12,2(a,a36),2a)')"# grid_id^ level@ proc # number  : ","mass"," [ ","position"," ] [ ","velocity"," ] outside, in, out, phy, fin"
         call printinfo(msg)
      endif

      cgl => leaves%first
      do while (associated(cgl))
         write(rem, '(a,i8,a,i6)')"#", cgl%cg%grid_id, "^", cgl%cg%l%id
         call cgl%cg%pset%print(rem)
         cgl => cgl%nxt
      enddo

   end subroutine print_all_particles

   subroutine particle_diagnostics(regular)

      use particle_pub, only: dump_diagnose
#ifdef VERBOSE
      use dataio_pub,   only: msg, printinfo
      use mpisetup,     only: master
#endif /* VERBOSE */

      implicit none

      logical, intent(in) :: regular              !< govern shorter diagnostics
      real, save          :: init_energy          !< total initial energy of set of the particles
      real, save          :: init_angmom          !< initial angular momentum of set of the particles
      logical, save       :: first_run_lf = .true.

      if (dump_diagnose .and. regular) call dump_particles_to_textfile

#ifndef VERBOSE
      if (regular) return
#endif /* !VERBOSE */

      call get_angmom_totener(tot_angmom, tot_energy)

      if (first_run_lf) then
         init_energy = tot_energy
         init_angmom = tot_angmom
         first_run_lf = .false.
      endif
      d_energy = log_error(tot_energy, init_energy)
      d_angmom = log_error(tot_angmom, init_angmom)

#ifdef VERBOSE
      if (master) then
         write(msg,'(a,3(1x,e12.5))') '[particle_diag:particle_diagnostics] Total energy: initial, current, error ', init_energy, tot_energy, d_energy
         call printinfo(msg)
         write(msg,'(a,3(1x,e12.5))') '[particle_diag:particle_diagnostics] ang_momentum: initial, current, error ', init_angmom, tot_angmom, d_angmom
         call printinfo(msg)
      endif
#endif /* VERBOSE */

   end subroutine particle_diagnostics

   real function log_error(vcurr, vinit)

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

      use allreduce,      only: piernik_MPI_Allreduce
      use cg_leaves,      only: leaves
      use cg_list,        only: cg_list_element
      use constants,      only: pSUM, xdim, ydim, zdim, zero
      use particle_types, only: particle

      implicit none

      real, intent(out)              :: ang_momentum
      real, intent(out)              :: total_energy !< total energy of set of particles
      real                           :: L1, L2, L3
      type(cg_list_element), pointer :: cgl
      type(particle), pointer        :: pset

      ang_momentum = zero
      total_energy = zero

      cgl => leaves%first
      do while (associated(cgl))

         pset => cgl%cg%pset%first
         do while (associated(pset))
            associate( part => pset%pdata )
               L1 = part%pos(ydim) * part%vel(zdim) - part%pos(zdim) * part%vel(ydim)
               L2 = part%pos(zdim) * part%vel(xdim) - part%pos(xdim) * part%vel(zdim)
               L3 = part%pos(xdim) * part%vel(ydim) - part%pos(ydim) * part%vel(xdim)
               ang_momentum = ang_momentum + part%mass * sqrt(L1**2 + L2**2 + L3**2)
               total_energy = total_energy + part%energy
            end associate
            pset => pset%nxt
         enddo

         cgl => cgl%nxt
      enddo

      call piernik_MPI_Allreduce(total_energy, pSUM)
      call piernik_MPI_Allreduce(ang_momentum, pSUM)

   end subroutine get_angmom_totener

   subroutine dump_particles_to_textfile

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: cwdlen, half, ndims, I_ZERO
      use dataio_pub,       only: log_wr, nrestart, problem_name, run_id
      use global,           only: t, dt
      use particle_gravity, only: get_acc_model
      use particle_pub,     only: twodtscheme
      use particle_types,   only: particle

      implicit none

      real                           :: kdt
      real, dimension(ndims)         :: acc2
      integer                        :: i, lun_out
      character(len=cwdlen)          :: plog_file
      type(cg_list_element), pointer :: cgl
      type(particle), pointer        :: pset

      kdt = dt ; if (.not.twodtscheme) kdt = half * dt

      write(plog_file,'(6a,i3.3,a)') trim(log_wr), '/', trim(problem_name), '_', trim(run_id), '_', max(I_ZERO, nrestart), '_out.log'
      open(newunit=lun_out, file=plog_file, status='unknown',  position='append')
      cgl => leaves%first
      do while (associated(cgl))

         pset => cgl%cg%pset%first
         do while (associated(pset))
            associate( part => pset%pdata )
               call get_acc_model(part%pos, part%mass, acc2)
               write(lun_out, '(a,I3.3,1X,19(E13.6,1X))') 'particle', i, t, kdt, part%mass, part%pos, part%vel, part%acc, acc2(:), part%energy
            end associate
            pset => pset%nxt
         enddo

         cgl => cgl%nxt
      enddo
      close(lun_out)

   end subroutine dump_particles_to_textfile

end module particle_diag
