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
!! Main program
!<
program hello_world
! pulled by ANY
! overrides src/base/piernik.F90

   use mpi, only: MPI_COMM_WORLD

   implicit none

   integer(kind=4) :: proc    !< rank of my process
   integer(kind=4) :: mpi_err !< error status

   call MPI_Init(mpi_err)
   call MPI_Comm_rank(MPI_COMM_WORLD, proc, mpi_err)

   if (proc == 0) then
      write(*,*)"This is a demonstration how to replace any of the standard Piernik files with your own version."                                                                                         ! QA_WARN intentional
      write(*,*)"It is not intended for regular and prolonged use as maintaining slightly changed files can be painful task."                                                                             ! QA_WARN intentional
      write(*,*)"This mechanism may be a bit useful for developing purposes, doing dirty experiments or for maintaining one-time problems that require special care without disturbing main source tree." ! QA_WARN intentional
   endif

   write(*,*)"Hello world from process #", proc ! QA_WARN intentional

   call MPI_Finalize(mpi_err)

end program hello_world
