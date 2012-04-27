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

   private
   public :: init_tracer, tracer_index, iarr_trc, trace_fluid

   integer, dimension(:), allocatable :: iarr_trc, trace_fluid
   integer, private :: ntracers
   integer, dimension(10) :: tracers

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

      namelist /FLUID_TRACER/ tracers

      tracers = [1,0,0,0,0,0,0,0,0,0]

      if (master) then
         diff_nml(FLUID_TRACER)

         ibuff(1:10)   = tracers
      endif

      call MPI_Bcast(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)

      if (slave) then

         tracers  = ibuff(1:10)

      endif

      ntracers = count(tracers > 0)

      ! TODO: deallocate those arrays somewhere
      allocate(trace_fluid(ntracers))
      allocate(iarr_trc(ntracers))

      trace_fluid = pack(tracers, mask=(tracers > 0))

   end subroutine init_tracer

   subroutine tracer_index(flind)
      use constants,          only: I_ONE
      use fluidtypes,         only: var_numbers

      implicit none
      type(var_numbers), intent(inout) :: flind
      integer :: i

      flind%trc%beg    = flind%all + I_ONE
      flind%trc%all  = ntracers
      flind%all      = flind%all + flind%trc%all
      flind%trc%end    = flind%all

      iarr_trc = [(i, i = 0, ntracers-1)] + flind%trc%beg
!      flind%components = flind%components + 1
!      flind%trc%pos    = flind%components
      flind%trc%pos    = -1

   end subroutine tracer_index

end module inittracer
