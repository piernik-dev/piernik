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

!> \brief Module for PPP-guarded MPI waitall call

module ppp_mpi

   use req_array, only: req_arr

   implicit none

   private
   public :: req_ppp

   type, extends(req_arr) :: req_ppp  ! array of requests with PPP capabilities
   contains
      procedure :: waitall_ppp  !< wait for completion of the requests on all procs
      generic, public :: waitall => waitall_ppp  ! overloading: the other method is waitall_on_some
   end type req_ppp

contains

!>
!! \brief A req_arr PPP wrapper for MPI_Waitall
!!
!! \details This is a bit expanded version of req_arr::waitall_on_some that had to be
!! moved here to avoid circular dependencies at ppp_eventlist.
!!
!! BEWARE: Cannot use extra_barriers when this routine is called only by a subset of MPI ranks.
!<

   subroutine waitall_ppp(this, ppp_label, x_mask)

      use barrier,   only: piernik_MPI_Barrier
      use constants, only: PPP_MPI
      use ppp,       only: ppp_main
      use req_array, only: C_REQA

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

         call this%waitall_wrapper(C_REQA, ppp_label)

         call ppp_main%stop(mpiw // ppp_label, mask)
      endif

      call this%cleanup

      call piernik_MPI_Barrier

   end subroutine waitall_ppp

end module ppp_mpi
