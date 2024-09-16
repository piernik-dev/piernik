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

!> \brief MPI wrappers

module mpi_wrappers

   implicit none

   private
   public :: piernik_MPI_Barrier, extra_barriers

   logical, save :: extra_barriers = .false.  !< when changed to .true. additional MPI_Barriers will be called.

   integer(kind=4) :: err_mpi  !< error status

contains

!> \brief Wrapper for MPI_Barrier

   subroutine piernik_MPI_Barrier

      use MPIF, only: MPI_COMM_WORLD, MPI_Barrier

      implicit none

      call MPI_Barrier(MPI_COMM_WORLD, err_mpi)

   end subroutine piernik_MPI_Barrier

end module mpi_wrappers
