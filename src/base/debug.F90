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

   use constants, only: cbuff_len

   implicit none

   private
   public :: init_piernikdebug, has_const_dt, constant_dt, force_dumps, aux_R, aux_I, aux_L, aux_S

   ! Auxiliary input parameters for debugging, quick tweaks and tests of new features.
   ! Their purpose is to avoid messing up existing namelists until it becomes clear that certain parameter is really useful.
   ! There is no reason to give them protected attribute.
   integer, parameter                        :: naux = 5 !< number of auxiliary variables of each kind
   real, dimension(naux)                     :: aux_R    !< real auxiliary parameter
   integer, dimension(naux)                  :: aux_I    !< integer auxiliary parameter
   logical, dimension(naux)                  :: aux_L    !< boolean auxiliary parameter
   character(len=cbuff_len), dimension(naux) :: aux_S    !< string auxiliary parameter

   real,    protected :: constant_dt               !< value of timestep regardless of fluid state
   logical, protected :: has_const_dt              !< true if piernikdebug::constant_dt > 0
   logical, protected :: force_hdf5_dump           !< dump hdf5 every sweep regardless of dataio_pub::dt_hdf
   logical, protected :: force_res_dump            !< dump restart every sweep regardless of dataio_pub::dt_res
   logical, protected :: force_allbnd_dump         !< dump restart with all boundaries every sweep regardless of dataio_pub::dt_res
   logical, protected :: force_log_dump            !< dump log every sweep regardless of dataio_pub:dt_log

   namelist /PIERNIK_DEBUG/ constant_dt, force_hdf5_dump, force_res_dump, force_allbnd_dump, force_log_dump, aux_R, aux_I, aux_L, aux_S

contains

   subroutine init_piernikdebug

      use constants,             only: PIERNIK_INIT_MPI
      use dataio_pub,            only: par_file, ierrh, namelist_errh, compare_namelist, cmdl_nml, lun  ! QA_WARN required for diff_nml
      use dataio_pub,            only: code_progress, die
      use mpi,                   only: MPI_DOUBLE_PRECISION, MPI_LOGICAL, MPI_CHARACTER, MPI_INTEGER
      use mpisetup,              only: master, slave, comm, ierr, rbuff, lbuff, cbuff, ibuff, buffer_dim, FIRST

      implicit none

      if (code_progress < PIERNIK_INIT_MPI) call die("[debug:init_piernikdebug] MPI not initialized.")

      constant_dt = 0.0
      aux_R(:) = 0.
      aux_I(:) = 0
      aux_L(:) = .false.
      aux_S(:) = ""

      if (master) then
         diff_nml(PIERNIK_DEBUG)

         rbuff(1) = constant_dt

         lbuff(1) = force_hdf5_dump
         lbuff(2) = force_log_dump
         lbuff(3) = force_res_dump
         lbuff(4) = force_allbnd_dump

         rbuff(buffer_dim-naux+1:buffer_dim) = aux_R(:)
         ibuff(buffer_dim-naux+1:buffer_dim) = aux_I(:)
         lbuff(buffer_dim-naux+1:buffer_dim) = aux_L(:)
         cbuff(buffer_dim-naux+1:buffer_dim) = aux_S(:)

      endif

      call MPI_Bcast(cbuff, cbuff_len*buffer_dim, MPI_CHARACTER,        FIRST, comm, ierr)
      call MPI_Bcast(ibuff,           buffer_dim, MPI_INTEGER,          FIRST, comm, ierr)
      call MPI_Bcast(rbuff,           buffer_dim, MPI_DOUBLE_PRECISION, FIRST, comm, ierr)
      call MPI_Bcast(lbuff,           buffer_dim, MPI_LOGICAL,          FIRST, comm, ierr)

      if (slave) then
         constant_dt       = rbuff(1)

         force_hdf5_dump   = lbuff(1)
         force_log_dump    = lbuff(2)
         force_res_dump    = lbuff(3)
         force_allbnd_dump = lbuff(4)

         aux_R(:) = rbuff(buffer_dim-naux+1:buffer_dim)
         aux_I(:) = ibuff(buffer_dim-naux+1:buffer_dim)
         aux_L(:) = lbuff(buffer_dim-naux+1:buffer_dim)
         aux_S(:) = cbuff(buffer_dim-naux+1:buffer_dim)

      endif

      has_const_dt = (constant_dt > 0.0)

   end subroutine init_piernikdebug

   subroutine force_dumps

      use common_hdf5,    only: set_container_chdf
      use dataio,         only: write_data
      use dataio_pub,     only: warn
      use data_hdf5,      only: write_hdf5
      use global,         only: nstep
      use restart_hdf5,   only: write_restart_hdf5

      implicit none

      call set_container_chdf(nstep)
      if (force_hdf5_dump)   call write_hdf5
      if (force_res_dump)    call write_restart_hdf5
      if (force_allbnd_dump) call warn("[fluidupdate:make_sweep] force_allbnd_dump has no effect for single-file HDF5 restart files")
      if (force_log_dump)    call write_data(output='log')

   end subroutine force_dumps

end module piernikdebug
