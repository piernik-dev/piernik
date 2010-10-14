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
#include "piernik.def"

module timestep

   implicit none
   real :: c_all ! BEWARE: assigned but unused

   contains

      subroutine time_step
         use mpisetup,      only: t, dt, dt_old, dt_max_grow, dt_initial, dt_min, nstep, proc
         use dataio_public, only: tend, msg, warn
         use constants,     only: small,big
         use dataio,        only: write_crashed

#ifdef IONIZED
         use timestepionized, only: timestep_ion
         use timestepionized, only: dt_ion,c_ion
#endif /* IONIZED */

#ifdef NEUTRAL
         use timestepneutral, only: timestep_neu
         use timestepneutral, only: dt_neu,c_neu
#endif /* NEUTRAL */

#ifdef DUST
         use timestepdust, only: timestep_dst
         use timestepdust, only: dt_dst,c_dst
#endif /* DUST */

#ifdef COSM_RAYS
         use timestepcosmicrays, only: timestep_crs
         use timestepcosmicrays, only: dt_crs
#endif /* COSM_RAYS */

#ifdef RESISTIVE
         use resistivity, only: dt_resist, timestep_resist
#endif /* RESISTIVE */

#ifdef FLUID_INTERACTIONS
         use timestepinteractions, only: timestep_interactions
         use timestepinteractions, only: dt_interact
#endif /* FLUID_INTERACTIONS */

         implicit none
! Timestep computation

         dt_old = dt

         c_all = 0.0
         dt    = (tend-t)/2.*(1+2.*epsilon(1.))

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
         call timestep_interactions
         dt = min(dt,dt_interact)
#endif /* FLUID_INTERACTIONS */

#ifdef SIMPLE_COOL
         dt = min(dt,0.01 * tauc)
#endif /* SIMPLE_COOL */

         ! finally apply some sanity factors
         if (nstep <=1) then
            if (dt_initial > 0.) dt = min(dt, dt_initial)
         else
            if (dt_old > 0.) dt = min(dt, dt_old*dt_max_grow)
         endif

         if (dt < dt_min) then ! something nasty had happened
            if (proc == 0) then
               write(msg,'(2(a,es12.4))')"[timestep:time_step] dt = ",dt,", less than allowed minimum = ",dt_min
               call warn(msg)
            endif
            call write_crashed("[timestep:time_step] dt < dt_min")
         endif

      end subroutine time_step
!------------------------------------------------------------------------------------------

end module timestep
