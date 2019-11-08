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
!! \brief Module containing a routine to compute upper limit of %timestep due to fluids %interactions.
!<

module timestepthermal
! pulled by THERM

   implicit none

   private
   public :: timestep_thermal

contains

!>
!! \brief Routine that computes upper limit of %timestep due to cooling or heating processes.
!<
   real function timestep_thermal() result(dt_coolheat)

      use cg_leaves, only: leaves
      use cg_list,   only: cg_list_element
      use constants, only: big, pMIN, small
      use grid_cont, only: grid_container
      use mpisetup,  only: piernik_MPI_Allreduce
      use thermal,   only: maxdeint, cfl_coolheat, thermal_active

      implicit none

      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
      real                           :: mxdeint

      dt_coolheat = big
      if (.not.thermal_active) return

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         call maxdeint(cg, mxdeint)
         dt_coolheat = min(dt_coolheat, cfl_coolheat*abs(1./(mxdeint+small)))

         cgl => cgl%nxt
      enddo
      call piernik_MPI_Allreduce(dt_coolheat, pMIN)

   end function timestep_thermal

end module timestepthermal
