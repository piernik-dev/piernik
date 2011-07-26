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

module initdust
! pulled by DUST
   implicit none

   public ! QA_WARN no secrets are kept here

   integer               :: idnd, imxd, imyd, imzd
   logical               :: selfgrav_dst

contains

!>
!! \brief Routine to set parameter values from namelist FLUID_DUST
!!
!! \n \n
!! @b FLUID_DUST
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>selfgrav_dst  </td><td>.false.</td><td>logical   </td><td>\copydoc initdust::selfgrav_dst  </td></tr>
!! </table>
!! \n \n
!<
   subroutine init_dust

      use dataio_pub,     only: par_file, ierrh, namelist_errh, compare_namelist, cmdl_nml  ! QA_WARN required for diff_nml
      use mpisetup,       only: master, slave, lbuff, buffer_dim, comm, ierr
      use mpi,            only: MPI_LOGICAL

      implicit none

      namelist /FLUID_DUST/ selfgrav_dst

      selfgrav_dst = .false.

      if (master) then
         diff_nml(FLUID_DUST)

         lbuff(1)   = selfgrav_dst
      endif

      call MPI_Bcast(lbuff,    buffer_dim, MPI_LOGICAL,          0, comm, ierr)

      if (slave) then

         selfgrav_dst    = lbuff(1)

      endif

   end subroutine init_dust

   subroutine dust_index(flind)
      use diagnostics,   only: ma1d, my_allocate
      use fluidtypes,    only: var_numbers
      use constants,     only: DST

      implicit none
      type(var_numbers), intent(inout) :: flind

      flind%dst%beg    = flind%all + 1

      idnd = flind%all + 1
      imxd = flind%all + 2
      imyd = flind%all + 3
      imzd = flind%all + 4

      flind%dst%idn  = idnd
      flind%dst%imx  = imxd
      flind%dst%imy  = imyd
      flind%dst%imz  = imzd

      flind%dst%all  = 4
      flind%all      = imzd

      ma1d = [flind%dst%all]
      call my_allocate(flind%dst%iarr,      ma1d)
      call my_allocate(flind%dst%iarr_swpx, ma1d)
      call my_allocate(flind%dst%iarr_swpy, ma1d)
      call my_allocate(flind%dst%iarr_swpz, ma1d)

      flind%dst%iarr      = [idnd,imxd,imyd,imzd]
      flind%dst%iarr_swpx = [idnd,imxd,imyd,imzd]
      flind%dst%iarr_swpy = [idnd,imyd,imxd,imzd]
      flind%dst%iarr_swpz = [idnd,imzd,imyd,imxd]

      flind%dst%end    = flind%all
      flind%components = flind%components + 1
      flind%fluids     = flind%fluids + 1
      flind%dst%pos    = flind%components
      if (selfgrav_dst)  flind%fluids_sg = flind%fluids_sg + 1

      flind%dst%gam = -1.
      flind%dst%cs  = 0.0
      flind%dst%cs2 = 0.0
      flind%dst%tag = DST
      flind%dst%is_selfgrav   = selfgrav_dst
      flind%dst%is_magnetized = .false.
      flind%dst%has_energy    = .false.

   end subroutine dust_index

end module initdust
