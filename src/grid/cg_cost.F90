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

!> \brief This module provides a structure for accumulating computational costs of a cg

module cg_cost

   use cg_cost_data, only: cg_cost_data_t

   implicit none

   private
   public :: cg_cost_t

   type, extends(cg_cost_data_t) :: cg_cost_t
      real, private :: wstart  !< start value of the timer
   contains
      procedure :: reset  !< Set all counters to 0.
      procedure :: start  !< Remember start time
      procedure :: stop   !< Stop measuring the time and add to specified timer
   end type cg_cost_t

   real, parameter :: T_INVALID = -huge(1.)

contains

!> \brief Set all counters to 0.

   subroutine reset(this)

      implicit none

      class(cg_cost_t), intent(out) :: this

      call this%init
      this%wstart = T_INVALID

   end subroutine reset

!> \brief Remember start time

   subroutine start(this)

      use dataio_pub, only: msg, warn
      use MPIF,       only: MPI_Wtime

      implicit none

      class(cg_cost_t), intent(inout) :: this

      real :: t

      t = MPI_Wtime()
      if (this%wstart > T_INVALID) then
         write(msg, '(2(a,f18.6))')"[cg_cost:start] Some counting has already begin at ", this%wstart, ". Resetting to ", t
         call warn(msg)
      endif
      this%wstart = t

   end subroutine start

!> \brief Stop measuring the time and add to specified timer

   subroutine stop(this, ind, ppp_c)

      use cg_cost_data, only: cost_labels
      use dataio_pub,   only: msg, warn, die
      use MPIF,         only: MPI_Wtime
      use ppp,          only: ppp_main

      implicit none

      class(cg_cost_t),  intent(inout) :: this  !<
      integer(kind=4),   intent(in)    :: ind   !<
      logical, optional, intent(in)    :: ppp_c !< If present and .false. then don't add the interval to ppp

      real :: t
      logical :: ppp_call

      t = MPI_Wtime()
      if (this%wstart <= T_INVALID) then
         write(msg, '(a,i3,a,f18.6)')"[cg_cost:stop] Counting hasn't been started! index: ", ind, " time: ", t
         call warn(msg)
      endif

      if (ind >= lbound(this%wtime, 1) .and. ind <= ubound(this%wtime, 1)) then
         this%wtime(ind) = this%wtime(ind) + (t - this%wstart)
         ppp_call = .true.
         if (present(ppp_c)) ppp_call = ppp_c
         if (ppp_call) call ppp_main%single_cg_cost(this%wstart, t, "cg_cost:" // cost_labels(ind))
      else
         call die("[cg_cost:stop] invalid index")
      endif

      this%wstart = T_INVALID

   end subroutine stop

end module cg_cost
