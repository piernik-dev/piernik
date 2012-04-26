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
#include "macros.h"
!>
!! \brief (MH) Initialization of the dust fluid
!!
!!
!!
!! In this module following namelist of parameters is specified:
!! \copydetails initdust::init_dust
!<

module inittracer
! pulled by TRACER
   implicit none

   public ! QA_WARN no secrets are kept here

   integer               :: itrc, trace_fluid

contains

!>
!! \brief Routine to set parameter values from namelist FLUID_DUST
!!
!! \n \n
!! @b FLUID_DUST
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>selfgrav_trc  </td><td>.false.</td><td>logical   </td><td>\copydoc initdust::selfgrav_trc  </td></tr>
!! </table>
!! \n \n
!<
   subroutine init_tracer

      use dataio_pub,     only: par_file, ierrh, namelist_errh, compare_namelist, cmdl_nml, lun  ! QA_WARN required for diff_nml
      use mpisetup,       only: master, slave, ibuff, buffer_dim, comm, ierr
      use mpi,            only: MPI_INTEGER

      implicit none

      namelist /FLUID_TRACER/ trace_fluid

      trace_fluid = 1

      if (master) then
         diff_nml(FLUID_TRACER)

         ibuff(1)   = trace_fluid
      endif

      call MPI_Bcast(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)

      if (slave) then

         trace_fluid  = ibuff(1)

      endif

   end subroutine init_tracer

   subroutine tracer_index(flind)
      use fluidtypes,         only: var_numbers

      implicit none
      type(var_numbers), intent(inout) :: flind

      flind%trc%beg    = flind%all + 1

      itrc = flind%all + 1

      flind%trc%all  = 1
      flind%all      = itrc

      flind%trc%end    = flind%all
      flind%components = flind%components + 1
      flind%trc%pos    = flind%components

   end subroutine tracer_index

end module inittracer
