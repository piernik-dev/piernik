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

module piernikdebug
! pulled by DEBUG

   implicit none

   private
   public :: init_piernikdebug, has_const_dt, constant_dt, force_hdf5_dump, force_log_dump

   real, protected    :: constant_dt               !< value of timestep regardless of fluid state
   logical, protected :: has_const_dt              !< true if piernikdebug::constant_dt > 0
   logical, protected :: force_hdf5_dump           !< dump hdf5 every sweep regardless of dataio_pub::dt_hdf
   logical, protected :: force_log_dump            !< dump log every sweep regardless of dataio_pub:dt_log

   namelist /PIERNIK_DEBUG/ constant_dt, force_hdf5_dump, force_log_dump

contains

   subroutine init_piernikdebug

      use mpisetup,              only: master, slave, comm, ierr, buffer_dim, rbuff, lbuff
      use mpi,                   only: MPI_DOUBLE_PRECISION, MPI_LOGICAL
      use dataio_pub,            only: par_file, ierrh, namelist_errh, compare_namelist, cmdl_nml  ! QA_WARN required for diff_nml
      use dataio_pub,            only: code_progress, PIERNIK_INIT_MPI

      implicit none

      if (code_progress < PIERNIK_INIT_MPI) call die("[debug:init_piernikdebug] MPI not initialized.")

      constant_dt = 0.0

      if (master) then
         diff_nml(PIERNIK_DEBUG)

         rbuff(1) = constant_dt

         lbuff(1) = force_hdf5_dump
         lbuff(2) = force_log_dump
      endif

      call MPI_Bcast(rbuff,           buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)
      call MPI_Bcast(lbuff,           buffer_dim, MPI_LOGICAL,          0, comm, ierr)

      if (slave) then
         constant_dt = rbuff(1)

         force_hdf5_dump = lbuff(1)
         force_log_dump  = lbuff(2)
      endif

      has_const_dt = (constant_dt > 0.0)

   end subroutine init_piernikdebug

end module piernikdebug
