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

module timestep

! Based on: Pen Arras & Wong (2003)
! Modified extensively by M. Hanasz November 2005 - May 2006
! Optimalization of cpu-time efficiency of mhdflux, tvd1 by K. Kowalik
! Modification history:  see "changelog"  -- warning: out of date

   use mpisetup

   implicit none
   real :: c_all

   contains


      subroutine time_step
         use start, only : dt, tend, t
         use constants, only : small,big

#ifdef IONIZED
         use timestepionized, only : timestep_ion
         use timestepionized, only : dt_ion,c_ion
#endif /* IONIZED */

#ifdef NEUTRAL
         use timestepneutral, only : timestep_neu
         use timestepneutral, only : dt_neu,c_neu
#endif /* NEUTRAL */

#ifdef DUST
         use timestepdust, only : timestep_dst
         use timestepdust, only : dt_dst,c_dst
#endif /* DUST */

#ifdef COSM_RAYS
         use timestepcosmicrays, only : timestep_crs
         use timestepcosmicrays, only : dt_crs
#endif /* COSM_RAYS */

#ifdef SIMPLE_COOL
         use start, only : tauc
#endif /* SIMPLE_COOL */

#ifdef RESISTIVE
         use resistivity, only : dt_resist, timestep_resist
#endif /* RESISTIVE */

#ifdef ANY_INTERACTIONS
         use timestepinteractions, only : timestep_interactions
         use timestepinteractions, only : dt_interact
#endif /* ANY_INTERACTIONS */

         implicit none
! Timestep computation

         c_all = 0.0
         dt    = (tend-t)/2.

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

#ifdef ANY_INTERACTIONS
         call timestep_interactions
         dt = min(dt,dt_interact)
#endif /* ANY_INTERACTIONS */

#ifdef SIMPLE_COOL
         dt = min(dt,0.01 * tauc)
#endif /* SIMPLE_COOL */
      end subroutine time_step
!------------------------------------------------------------------------------------------

end module timestep
