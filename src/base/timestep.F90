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

!>
!! \brief This module gathers all applicable timestep limits and computes next timestep.
!<

module timestep

   implicit none

   private
   public :: init_time_step, time_step, c_all, cfl_manager

   real :: c_all, c_all_old
   procedure(), pointer :: cfl_manager => null()

contains

!>
!! \brief Initialization routine
!!
!! \details This routine sets cfl_manager according to mpisetup::cflcontrol parameter.
!<

   subroutine init_time_step

      use dataio_pub, only: msg, die, warn, code_progress, PIERNIK_INIT_MPI
      use mpisetup,   only: cflcontrol

      implicit none

      if (code_progress < PIERNIK_INIT_MPI) call die("[grid:init_grid] MPI not initialized.") ! cflcontrol

      if (associated(cfl_manager)) call die("[timestep:init_time_step] cfl_manager already associated.")
      select case (cflcontrol)
         case ('warn')
            cfl_manager => cfl_warn
         case ('auto', 'adaptive')
            cfl_manager => cfl_auto
         case ('none', '')
         case default
            write(msg, '(3a)')"[timestep:init_time_step] Unknown cfl_manager '",trim(cflcontrol),"'. Assuming 'none'."
            call warn(msg)
      end select

   end subroutine init_time_step

!>
!! \brief Timestep calculation
!!
!! \details This routine calls various routines associated with different modules.
!! These routines return limit for timestep due to various physical and numerical conditions.
!! At the end the timestep is checked against remaining simulation time, mimimum, and maximum allowed values etc.
!<
   subroutine time_step(dt)

      use dataio,               only: write_crashed
      use dataio_pub,           only: tend, msg, warn
      use mpisetup,             only: t, dt_old, dt_max_grow, dt_initial, dt_min, nstep, master, cflcontrol
#ifdef IONIZED
      use timestepionized,      only: timestep_ion, dt_ion, c_ion
#endif /* IONIZED */
#ifdef NEUTRAL
      use timestepneutral,      only: timestep_neu, dt_neu, c_neu
#endif /* NEUTRAL */
#ifdef DUST
      use timestepdust,         only: timestep_dst, dt_dst, c_dst
#endif /* DUST */
#ifdef COSM_RAYS
      use timestepcosmicrays,   only: timestep_crs, dt_crs
#endif /* COSM_RAYS */
#ifdef RESISTIVE
      use resistivity,          only: dt_resist, timestep_resist
#endif /* RESISTIVE */
#ifdef FLUID_INTERACTIONS
      use timestepinteractions, only: timestep_interactions
#endif /* FLUID_INTERACTIONS */
#ifdef DEBUG
      use piernikdebug,         only: has_const_dt, constant_dt
      use dataio_pub,           only: printinfo
#endif /* DEBUG */

      implicit none

      real, intent(inout) :: dt !< the timestep
! Timestep computation

      dt_old = dt

      c_all = 0.0
      dt = huge(1.0)

#ifdef IONIZED
      call timestep_ion
      dt=min(dt,dt_ion)
      c_all = max(c_all,c_ion)
#endif /* IONIZED */

#ifdef NEUTRAL
      call timestep_neu
      dt=min(dt,dt_neu)
      c_all = max(c_all,c_neu)
#endif /* NEUTRAL */

#ifdef DUST
      call timestep_dst
      dt=min(dt,dt_dst)
      c_all = max(c_all,c_dst)
#endif /* DUST */

#ifdef COSM_RAYS
      call timestep_crs
      dt=min(dt,dt_crs)
#endif /* COSM_RAYS */

#ifdef RESISTIVE
      call timestep_resist
      dt = min(dt,dt_resist)
#endif /* RESISTIVE */

#ifdef FLUID_INTERACTIONS
      dt = min(dt,timestep_interactions())
#endif /* FLUID_INTERACTIONS */

      ! finally apply some sanity factors
      if (nstep <=1) then
         if (dt_initial > 0.) dt = min(dt, dt_initial)
      else
         if (dt_old > 0.) dt = min(dt, dt_old*dt_max_grow)
      endif

      if (associated(cfl_manager)) call cfl_manager
      c_all_old = c_all

      if (dt < dt_min) then ! something nasty had happened
         if (master) then
            write(msg,'(2(a,es12.4))')"[timestep:time_step] dt = ",dt,", less than allowed minimum = ",dt_min
            call warn(msg)
         endif
         call write_crashed("[timestep:time_step] dt < dt_min")
      endif

      dt  = min(dt, (0.5*(tend-t)) + (2.*epsilon(1.)*((tend-t))))
#ifdef DEBUG
      ! We still need all above for c_all
      if (has_const_dt) then
         dt    = constant_dt
         write(msg,*) "[timestep:time_step]: (constant_dt) c_all = ", c_all
         call printinfo(msg)
      endif
#endif /* DEBUG */

   end subroutine time_step

!------------------------------------------------------------------------------------------
!>
!! \brief This routine detects sudden timestep changes due to strong velocity changes and interprets them as possible CFL criterion violations
!<

   subroutine cfl_warn

      use dataio_pub, only: msg, warn
      use mpisetup,   only: cfl, cfl_max, master

      implicit none

      real :: stepcfl

      stepcfl = cfl
      if (c_all_old > 0.) stepcfl = c_all/c_all_old*cfl

      if (master) then
         msg = ""
         if (stepcfl > cfl_max) then
            write(msg,'(a,g10.3)') "[timestep:cfl_warn] Possible violation of CFL: ",stepcfl
         else if (stepcfl < 2*cfl - cfl_max) then
            write(msg,'(2(a,g10.3))') "[timestep:cfl_warn] Low CFL: ", stepcfl, " << ", cfl
         endif
         if (len_trim(msg) > 0) call warn(msg)
      endif

   end subroutine cfl_warn

!------------------------------------------------------------------------------------------
!>
!! \brief This routine detects sudden timestep changes due to strong velocity changes and interprets them as possible CFL criterion violations
!! Timestep changes are used to estimate safe value of CFL for the next timestep (EXPERIMENTAL)
!<

   subroutine cfl_auto

      use dataio_pub, only: msg, warn
      use mpisetup,   only: cfl, cfl_max, master, dt, dt_old

      implicit none

      real, save :: stepcfl=0., cfl_c = 1.
      real       :: stepcfl_old

      stepcfl_old = stepcfl
      stepcfl = cfl
      if (c_all_old > 0.) then
         stepcfl = c_all/c_all_old*cfl
      else
         stepcfl_old = cfl
      endif

      if (stepcfl > 0. .and. dt_old > 0.) then
         cfl_c = min(1., 0.5 * (cfl_c + min(1., stepcfl_old/stepcfl *dt/dt_old)))
         dt = dt * cfl_c
      else
         cfl_c = 1.
      endif

      if (master) then
         msg = ""
         if (stepcfl > cfl_max) then
            write(msg,'(a,g10.3)') "[timestep:cfl_auto] Possible violation of CFL: ",stepcfl
         else if (stepcfl < 2*cfl - cfl_max) then
            write(msg,'(2(a,g10.3))') "[timestep:cfl_auto] Low CFL: ", stepcfl, " << ", cfl
         endif
         if (len_trim(msg) > 0) call warn(msg)
      endif

   end subroutine cfl_auto

end module timestep
