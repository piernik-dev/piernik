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

!> \brief Module that contains HDF5 I/O routines for reading and writing restart files.

module restart_hdf5
! pulled by HDF5

   implicit none

   private
   public :: read_restart_hdf5, write_restart_hdf5

contains

!>
!! \brief Wrapper routine that writes a version 2.x or version 1.x restart file, depending on use_v2_io switch
!<

   subroutine write_restart_hdf5(sequential)

      use common_hdf5,     only: dump_announcement, dump_announce_time, set_common_attributes
      use constants,       only: cwdlen, PPP_IO, RES
      use dataio_pub,      only: use_v2_io, nres, last_res_time
      use mpi_wrappers,    only: piernik_MPI_Barrier
      use ppp,             only: ppp_main
      use restart_hdf5_v1, only: write_restart_hdf5_v1
      use restart_hdf5_v2, only: write_restart_hdf5_v2
#if defined(MULTIGRID) && defined(SELF_GRAV)
      use multigrid_gravity, only: unmark_oldsoln
#endif /* MULTIGRID && SELF_GRAV */

      implicit none

      logical,         intent(in) :: sequential
      character(len=cwdlen)       :: filename  ! File name
      character(len=*), parameter :: wrr_label = "IO_write_restart"

      call ppp_main%start(wrr_label, PPP_IO)

      call dump_announcement(RES, nres, filename, last_res_time, sequential)

      call set_common_attributes(filename)

      if (use_v2_io) then
         call write_restart_hdf5_v2(filename)
      else
         call write_restart_hdf5_v1(filename)
      endif

#if defined(MULTIGRID) && defined(SELF_GRAV)
      call unmark_oldsoln
#endif /* MULTIGRID && SELF_GRAV */

      call piernik_MPI_Barrier

      call dump_announce_time

      call ppp_main%stop(wrr_label, PPP_IO)

   end subroutine write_restart_hdf5

!>
!! \brief Wrapper routine that reads a version 2.x or version 1.x restart file, depending on use_v2_io switch. On v2 failure it falls back to v1.
!<

   subroutine read_restart_hdf5

      use common_hdf5,     only: STAT_OK, STAT_INV
      use dataio_pub,      only: use_v2_io
      use restart_hdf5_v1, only: read_restart_hdf5_v1
      use restart_hdf5_v2, only: read_restart_hdf5_v2
#if defined(MULTIGRID) && defined(SELF_GRAV)
      use multigrid_gravity, only: unmark_oldsoln
#endif /* MULTIGRID && SELF_GRAV */

      implicit none

      integer :: status_v2

      status_v2 = STAT_INV
      if (use_v2_io) call read_restart_hdf5_v2(status_v2)
      if (status_v2 /= STAT_OK .or. .not. use_v2_io) call read_restart_hdf5_v1

#if defined(MULTIGRID) && defined(SELF_GRAV)
      call unmark_oldsoln
#endif /* MULTIGRID && SELF_GRAV */

   end subroutine read_restart_hdf5

end module restart_hdf5
