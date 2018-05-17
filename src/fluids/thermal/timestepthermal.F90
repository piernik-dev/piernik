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
   real function timestep_thermal(cg) result(dt)
      use grid_cont,    only: grid_container
      use constants,    only: small, I_ONE
      use grid_cont,    only: grid_container
      use mpi,          only: MPI_MIN, MPI_DOUBLE_PRECISION
      use mpisetup,     only: comm, mpi_err, FIRST, piernik_MPI_Bcast
      use thermal,      only: maxdeint, cfl_coolheat, thermal_active

      implicit none

      real :: dt_coolheat_proc        !< timestep due to cooling or heating processes for the current process (MPI block) only
      real :: dt_coolheat_all         !< timestep due to cooling or heating processes for all MPI blocks
      real :: mxdeint

      type(grid_container), pointer, intent(in) :: cg

      if (thermal_active) then
         call maxdeint(cg, mxdeint)
         dt_coolheat_proc = ABS(1./(mxdeint+small))
         call MPI_Reduce(dt_coolheat_proc, dt_coolheat_all, I_ONE, MPI_DOUBLE_PRECISION, MPI_MIN, FIRST, comm, mpi_err)
         call piernik_MPI_Bcast(dt_coolheat_all)
         dt = cfl_coolheat*dt_coolheat_all
      else
         dt = huge(1.)
      endif
    end function timestep_thermal

end module timestepthermal
