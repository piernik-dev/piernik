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
!! \brief Constants used in communication with master process
!!
!! This module contains various constants used during communications with master
!! process, when Piernik is launched via MPI_Spawn.
!! After every change, please make sure that it\'s parseable by whatever external
!! scripts which uses it.
!<
module piernik_mpi_sig
! pulled by ANY
   implicit none
   public :: sig
   private

   enum, bind(C)
      enumerator :: MSG_CLEAN_EXIT = 0
      enumerator :: MSG_TSL_UPDATED
      enumerator :: MSG_LOG_UPDATED
      enumerator :: MSG_HDF_WRITTEN
      enumerator :: MSG_RES_WRITTEN
      enumerator :: MSG_DIED = -1
   end enum

   type :: piernik_mpi_signals_t
      integer(kind=4) :: clean_exit = MSG_CLEAN_EXIT
      integer(kind=4) :: tsl_updated = MSG_TSL_UPDATED
      integer(kind=4) :: log_updated = MSG_LOG_UPDATED
      integer(kind=4) :: hdf_written = MSG_HDF_WRITTEN
      integer(kind=4) :: res_written = MSG_RES_WRITTEN
      integer(kind=4) :: died = MSG_DIED
   end type piernik_mpi_signals_t

   type(piernik_mpi_signals_t), protected :: sig
end module piernik_mpi_sig
