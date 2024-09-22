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

!> \brief Module for PPP-guarded MPI waitall calls

module ppp_mpi

   use req_array, only: req_arr

   implicit none

   private
   public :: init_wall, cleanup_wall, req_ppp

   type, extends(req_arr) :: req_ppp  ! array of requests with PPP capabilities
   contains
      procedure :: waitall_ppp      !< wait for completion of the requests on all procs
      generic, public :: waitall => waitall_ppp
   end type req_ppp

contains

!> \brief Initialize MPI_Waitall stat counters

   subroutine init_wall

      use mpi_wrappers, only: C_REQS, req_wall

      implicit none

      call req_wall%init(int([C_REQS]))

   end subroutine init_wall

!> \brief Print log and clean up MPI_Waitall stat counters

   subroutine cleanup_wall

      use constants,    only: V_DEBUG
      use mpi_wrappers, only: MPI_wrapper_stats, req_wall
      use mpisetup,     only: master

      implicit none

      if (master .and. MPI_wrapper_stats) &
           call req_wall%print("Waitall requests(calls). Columns: req_all%, req_some%.", V_DEBUG)

      call req_wall%cleanup

   end subroutine cleanup_wall

!>
!! \brief A req_arr PPP wrapper for MPI_Waitall
!!
!! \details This is a bit expanded version of req_arr::waitall_on_some that had to be
!! moved here to avoid circular dependencies at ppp_eventlist.
!!
!! BEWARE: Cannot use extra_barriers when this routine is called only by a subset of MPI ranks.
!<

   subroutine waitall_ppp(this, ppp_label, x_mask)

      use constants,    only: PPP_MPI, LONG
      use mpi_wrappers, only: piernik_MPI_Barrier, extra_barriers, C_REQA, req_wall
      use mpisetup,     only: err_mpi
      use MPIF,         only: MPI_STATUSES_IGNORE
      use MPIFUN,       only: MPI_Waitall, MPI_Comm_free
      use ppp,          only: ppp_main

      implicit none

      class(req_ppp),            intent(inout) :: this       !< an object invoking the type-bound procedure
      character(len=*),          intent(in)    :: ppp_label  !< identifier for PPP entry
      integer(kind=4), optional, intent(in)    :: x_mask     !< extra mask, if necessary

      character(len=*), parameter :: mpiw = "req:MPI_Waitall:"
      integer(kind=4) :: mask

      if (this%n > 0) then
         mask = PPP_MPI
         if (present(x_mask)) mask = mask + x_mask
         call ppp_main%start(mpiw // ppp_label, mask)

         call req_wall%add(int([C_REQA]), int(this%n, kind=LONG))

         call MPI_Waitall(this%n, this%r(:this%n), MPI_STATUSES_IGNORE, err_mpi)
         this%n = 0

         call ppp_main%stop(mpiw // ppp_label, mask)
      endif

      if (allocated(this%r)) deallocate(this%r)
      if (this%owncomm) then
         call MPI_Comm_free(this%comm, err_mpi)
         this%owncomm = .false.
      endif

      if (extra_barriers) call piernik_MPI_Barrier

   end subroutine waitall_ppp

end module ppp_mpi
