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

!> \brief Module for PPP-guarded MPI calls

module ppp_mpi

   implicit none

   private
   public :: piernik_Waitall

contains

!>
!! \brief a PPP wrapper for MPI_Waitall
!!
!! Cannot use extra_barriers when this routine is called only by a subset of MPI ranks.
!<

   subroutine piernik_Waitall(nr, ppp_label, x_mask, use_req2)

      use constants,    only: PPP_MPI
      use mpi_wrappers, only: piernik_MPI_Barrier, extra_barriers
      use mpisetup,     only: err_mpi, req, req2
      use MPIF,         only: MPI_STATUSES_IGNORE
      use MPIFUN,       only: MPI_Waitall
      use ppp,          only: ppp_main

      implicit none

      integer(kind=4),           intent(in) :: nr         !< number of requests in req(:) or req2(:)
      character(len=*),          intent(in) :: ppp_label  !< identifier for PPP entry
      integer(kind=4), optional, intent(in) :: x_mask     !< extra mask, if necessary
      logical,         optional, intent(in) :: use_req2   !< use req2 when .true.

      character(len=*), parameter :: mpiw = "MPI_Waitall:"
      integer(kind=4) :: mask
      logical :: r2

      if (nr > 0) then

         mask = PPP_MPI
         if (present(x_mask)) mask = mask + x_mask
         call ppp_main%start(mpiw // ppp_label, mask)

         r2 = .false.
         if (present(use_req2)) r2 = use_req2
         if (r2) then
            call MPI_Waitall(nr, req2(:nr), MPI_STATUSES_IGNORE, err_mpi)
         else
            call MPI_Waitall(nr, req(:nr), MPI_STATUSES_IGNORE, err_mpi)
         endif

         call ppp_main%stop(mpiw // ppp_label, mask)
      endif

      if (extra_barriers) call piernik_MPI_Barrier

   end subroutine piernik_Waitall

end module ppp_mpi
