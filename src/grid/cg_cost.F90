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

   implicit none

   private
   public :: cg_cost_data_t, cg_cost_t, cost_labels, I_MHD, I_MULTIGRID, I_MULTIPOLE, I_DIFFUSE, I_PARTICLE, I_REFINE, I_IC, I_OTHER

   enum, bind(C)
      enumerator :: I_MHD, I_MULTIGRID, I_MULTIPOLE, I_DIFFUSE, I_PARTICLE, I_REFINE, I_IC, I_OTHER
   end enum

   character(len=*), dimension(I_MHD:I_OTHER), parameter :: cost_labels = &
        [ "MHD       ", &  ! Riemann, RTVD, CT, divB cleaning
        & "multigrid ", &  ! self-gravity multigrid relaxation, residuals
        & "multipole ", &  ! multipole moments <=> potential conversions, not the costs of Q array manipulation
        & "diffusion ", &  ! explicit and multigrid diffussion costs
        & "particles ", &  ! particles in cg
        & "refines   ", &  ! prolongation, restriction, marking criteria
        & "init.cond.", &  ! for use only in the initproblems
        & "other     " ]   ! everything else that is related to cg but not tied to particular algorithm

   type :: cg_cost_data_t
      real, dimension(I_MHD:I_OTHER) :: wtime  ! walltime costs split into different categories
   end type cg_cost_data_t

   type, extends(cg_cost_data_t) :: cg_cost_t
      real, private :: wstart                           ! start value of the timer
   contains
      procedure :: reset  !< Set all counters to 0.
      procedure :: total  !< Return accumulated cost
      procedure :: start  !< Remember start time
      procedure :: stop   !< Stop measuring the time and add to specified timer
      procedure :: get    !< Read specified timer
   end type cg_cost_t

   real, parameter :: T_INVALID = -huge(1.)

contains

!> \brief Set all counters to 0.

   subroutine reset(this)

      implicit none

      class(cg_cost_t), intent(out) :: this

      this%wtime = 0.
      this%wstart = T_INVALID

   end subroutine reset

!> \brief Return accumulated cost

   real function total(this)

      implicit none

      class(cg_cost_t), intent(in) :: this

      total = sum(this%wtime(:))

   end function total

!> \brief Remember start time

   subroutine start(this)

      use dataio_pub, only: die
      use MPIF,       only: MPI_Wtime

      implicit none

      class(cg_cost_t), intent(inout) :: this

      if (this%wstart > T_INVALID) call die("[cg_cost:start] some counting has already begin")
      this%wstart =  MPI_Wtime()

   end subroutine start

!> \brief Stop measuring the time and add to specified timer

   subroutine stop(this, ind)

      use dataio_pub, only: die
      use MPIF,       only: MPI_Wtime

      implicit none

      class(cg_cost_t), intent(inout) :: this
      integer(kind=4),  intent(in)    :: ind

      if (this%wstart <= T_INVALID) call die("[cg_cost:start] counting hasn't been started")

      if (ind >= lbound(this%wtime, 1) .and. ind <= ubound(this%wtime, 1)) then
         this%wtime(ind) = this%wtime(ind) + MPI_Wtime() - this%wstart
      else
         call die("[cg_cost:start] invalid index")
      endif

      this%wstart = T_INVALID

   end subroutine stop

!> \brief Read specified timer

   real function get(this, ind)

      use dataio_pub, only: die

      implicit none

      class(cg_cost_t), intent(in) :: this
      integer(kind=4),  intent(in) :: ind

      if (ind >= lbound(this%wtime, 1) .and. ind <= ubound(this%wtime, 1)) then
         get = this%wtime(ind)
      else
         call die("[cg_cost:start] invalid index")
         get = 0.
      endif

   end function get

end module cg_cost
