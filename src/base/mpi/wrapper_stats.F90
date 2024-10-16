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

!> \brief Call inits and cleanups for MPI wrappers

module wrapper_stats

   implicit none

   private
   public :: init_wstats, cleanup_wstats

contains

!> \brief Initialize MPI wrapper stats counters

   subroutine init_wstats

      use allreduce,    only: init_allreduce
      use barrier,      only: init_bar
      use bcast,        only: init_bcast
      use isend_irecv,  only: init_sr
      use req_array,    only: init_wall

      implicit none

      call init_allreduce
      call init_bcast
      call init_sr
      call init_wall
      call init_bar

   end subroutine init_wstats

!> \brief Clean up MPI wrapper stats counters

   subroutine cleanup_wstats

      use allreduce,    only: cleanup_allreduce
      use barrier,      only: cleanup_bar
      use bcast,        only: cleanup_bcast
      use constants,    only: stdout, V_DEBUG
      use dataio_pub,   only: piernik_verbosity
      use isend_irecv,  only: cleanup_sr
      use mpisetup,     only: master
      use req_array,    only: cleanup_wall

      implicit none

      ! Be nice to the stream of finalization dots on stdout
      if (master .and. (piernik_verbosity <= V_DEBUG)) write(stdout,'()')

      call cleanup_allreduce
      call cleanup_bcast
      call cleanup_sr
      call cleanup_wall
      call cleanup_bar

   end subroutine cleanup_wstats

end module wrapper_stats
