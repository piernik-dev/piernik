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
!! \brief Initialization of the tracer fluid
!!
!!
!!
!! In this module following namelist of parameters is specified:
!! \copydetails inittracer::init_tracer
!<

module inittracer
! pulled by TRACER
   implicit none

   private
   public :: init_tracer, tracer_index, iarr_trc, trace_fluid

   integer(kind=4), dimension(:), allocatable :: iarr_trc, trace_fluid
   integer(kind=4) :: ntracers
   integer, parameter :: tracers_max = 10
   integer(kind=4), dimension(tracers_max) :: tracers

contains

!>
!! \brief Routine to set parameter values from namelist FLUID_DUST
!!
!! \n \n
!! @b FLUID_TRACER
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>tracers</td><td>0</td><td>10-element integer array</td><td>\copydoc inittracer::tracers</td></tr>
!! </table>
!! The list is active while \b "TRACER" is defined.
!! \n \n
!<
   subroutine init_tracer

      use constants,  only: INT4
      use dataio_pub, only: nh  ! QA_WARN required for diff_nml
      use dataio_pub, only: warn
      use mpisetup,   only: master, slave, ibuff, piernik_MPI_Bcast

      implicit none

      namelist /FLUID_TRACER/ tracers

      tracers(:) = 0_INT4; tracers(1) = 1_INT4 ! activate only first tracer fluid by default

      if (master) then
         diff_nml(FLUID_TRACER)

         ibuff(1:tracers_max)   = tracers
      endif

      call piernik_MPI_Bcast(ibuff)

      if (slave) then

         tracers  = int(ibuff(1:tracers_max), kind=4)

      endif

      ntracers = count(tracers > 0, kind=4)

      if (ntracers < 1) call warn("[inittracer:init_tracer] all tracer fluids are disabled")

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

      iarr_trc = int([(i, i = 0, ntracers-1)], kind=4) + flind%trc%beg
!      flind%components = flind%components + 1
!      flind%trc%pos    = flind%components
      flind%trc%pos    = -1

   end subroutine tracer_index

end module inittracer
