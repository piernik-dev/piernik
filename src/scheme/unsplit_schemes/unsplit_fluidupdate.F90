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
!! The job this module is to update boundaries. This was previously inside the module sweeps . Moved here so that this can be called inside 
!! solver after each time step on the entire 
!<

module unsplit_fluidupdate

! pulled by ANY

   implicit none

   private
   public :: fluid_update_unsplit

contains

    subroutine fluid_update_unsplit

      use dataio_pub,     only: halfstep
      use global,         only: dt, dtm, t
      use hdc,            only: update_chspeed,glmdamping, eglm
      use mass_defect,    only: update_magic_mass
      use timestep_retry, only: repeat_fluidstep
      use unsplit_sweeps, only: unsplit_sweep

      implicit none

      call repeat_fluidstep
      call update_chspeed


      call eglm
      call glmdamping(.true.)
      t = t + dt

      call unsplit_sweep

      dtm = dt

      call eglm
      call glmdamping(.true.)
      call update_magic_mass

    end subroutine fluid_update_unsplit

end module unsplit_fluidupdate