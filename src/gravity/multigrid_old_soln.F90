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

!> \brief Multigrid historical solutions

module multigrid_old_soln
! pulled by MULTIGRID && GRAV

   implicit none

   private
   public :: nold_max, soln_history, ord_time_extrap

   ! solution recycling
   integer(kind=4), parameter :: nold_max=3                           !< maximum implemented extrapolation order

   type :: old_soln                                                   !< container for an old solution with its timestamp
      integer :: i_hist                                               !< index to the old solution
      real :: time                                                    !< time of the old solution
   end type old_soln

   type :: soln_history                                               !< container for a set of several old potential solutions
      type(old_soln), allocatable, dimension(:) :: old                !< indices and time points os stored solutions
      integer :: last                                                 !< index of the last stored potential
      logical :: valid                                                !< .true. when old(last) was properly initialized
    contains
      procedure :: init_history                                       !< Allocate arrays, register fields
      procedure :: cleanup_history                                    !< Deallocate arrays
      procedure :: init_solution                                      !< Construct first guess of potential based on previously obtained solution, if any.
      procedure :: store_solution                                     !< Manage old copies of potential for recycling.
   end type soln_history

   ! Namelist parameter
   integer(kind=4)    :: ord_time_extrap                              !< Order of temporal extrapolation for solution recycling; -1 means 0-guess, 2 does parabolic interpolation

contains

!> \brief Initialize structure for keeping historical potential fields

   subroutine init_history(this, nold, prefix)

      use cg_list_global,   only: all_cg
      use constants,        only: singlechar, dsetnamelen
      use dataio_pub,       only: die
      use multigridvars,    only: ord_prolong
      use named_array_list, only: qna

      implicit none

      class(soln_history),       intent(inout) :: this   !< Object invoking type-bound procedure
      integer,                   intent(in)    :: nold   !< how many historic points to store
      character(len=singlechar), intent(in)    :: prefix !< 'i' for inner, 'o' for outer potential

      integer :: i
      character(len=dsetnamelen) :: hname

      if (allocated(this%old)) call die("[multigrid_old_soln:init_history] Already allocated")
      allocate(this%old(nold))
      do i = lbound(this%old, dim=1), ubound(this%old, dim=1)
         write(hname,'(2a,i2.2)')prefix,"-h-",i
         call all_cg%reg_var(hname, vital = .true., ord_prolong = ord_prolong) ! no need for multigrid attribute here because history is defined only on leaves
         this%old(i) = old_soln(qna%ind(hname), -huge(1.0))
      enddo
      this%valid = .false.
      this%last  = 1

   end subroutine init_history

!> \brief Deallocate history arrays

   subroutine cleanup_history(this)

      use constants, only: INVALID

      implicit none

      class(soln_history), intent(inout) :: this   !< Object invoking type-bound procedure

      if (allocated(this%old)) deallocate(this%old)

      this%valid = .false.
      this%last  = INVALID

   end subroutine cleanup_history

!>
!! \brief This routine tries to construct first guess of potential based on previously obtained solution, if any.
!! If for some reason it is not desired to do the extrapolation (i.e. timestep varies too much) it would be good to have more control on this behavior.
!!
!! \details Quadratic extrapolation in time often gives better guess for smooth potentials, but is more risky for sharp-peaked potential fields (like moving self-bound clumps)
!! Rational extrapolation (1 + a t)/(b + c t) can give 2 times better and 10 times worse guess depending on timestep.
!!
!! Set history%valid to .false. to force start from scratch.
!!
!! \todo fix historical solutions after refinement update
!!
!! \todo Add an option to clear prolonged history to 0. and check it against high order prolongations
!<

   subroutine init_solution(this, prefix)

      use cg_list_dataop, only: ind_val
      use cg_leaves,      only: leaves
      use cg_list_global, only: all_cg
      use constants,      only: INVALID, O_INJ, O_LIN, O_I2
      use dataio_pub,     only: msg, die, printinfo
      use global,         only: t
      use mpisetup,       only: master
      use multigridvars,  only: stdout, solution

      implicit none

      class(soln_history), intent(inout) :: this    !< inner or outer potential history for recycling
      character(len=*),    intent(in)    :: prefix  !< informational string to be printed

      integer :: p0, p1, p2, ordt, h
      real, dimension(3) :: dt_fac

      call all_cg%set_dirty(solution)

      p0 = this%last

      associate (nold => size(this%old))
      if (nold > 0) then
         p1 = 1 + mod(p0 + nold - 2, nold) ! index of previous save
         p2 = 1 + mod(p1 + nold - 2, nold)
      else
         p1 = p0
         p2 = p0
      endif
      end associate

      ! selfgrav_clump/initproblem.F90 requires monotonic time sequence t > this%old(p0)%time > this%old(p1)%time > this%old(p2)%time
      ordt = ord_time_extrap
      if (this%valid) then
         if ( this%old(p2)%time < this%old(p1)%time .and. &        ! quadratic interpolation
              this%old(p1)%time < this%old(p0)%time .and. &
              this%old(p0)%time < t) then
            ordt = min(O_I2, ord_time_extrap)
         else if (this%old(p1)%time < this%old(p0)%time .and. &
              &   this%old(p0)%time < t) then                      ! linear extrapolation
            ordt = min(O_LIN, ord_time_extrap)
         else                                                      ! simple recycling
            ordt = min(O_INJ, ord_time_extrap)
         endif
      else                                                         ! coldstart
         ordt = min(INVALID, ord_time_extrap)
      endif

      select case (ordt)
         case (:INVALID)
            if (master .and. ord_time_extrap > ordt) then
               write(msg, '(3a)')"[multigrid_gravity:init_solution] Clearing ",trim(prefix),"solution."
               call printinfo(msg, stdout)
            endif
            call all_cg%set_q_value(solution, 0.)
            if (allocated(this%old)) this%old(:)%time = -huge(1.0)
         case (O_INJ)
            call leaves%check_dirty(this%old(p0)%i_hist, "history0")
            call leaves%q_copy(this%old(p0)%i_hist, solution)
            if (master .and. ord_time_extrap > ordt) then
               write(msg, '(3a)')"[multigrid_gravity:init_solution] No extrapolation of ",trim(prefix),"solution."
               call printinfo(msg, stdout)
            endif
         case (O_LIN)
            call leaves%check_dirty(this%old(p0)%i_hist, "history0")
            call leaves%check_dirty(this%old(p1)%i_hist, "history1")
            dt_fac(1) = (t - this%old(p0)%time) / (this%old(p0)%time - this%old(p1)%time)
            call leaves%q_lin_comb( [ ind_val(this%old(p0)%i_hist, (1.+dt_fac(1))), &
                 &                    ind_val(this%old(p1)%i_hist,    -dt_fac(1) ) ], solution )
            if (master .and. ord_time_extrap > ordt) then
               write(msg, '(3a)')"[multigrid_gravity:init_solution] Linear extrapolation of ",trim(prefix),"solution."
               call printinfo(msg, stdout)
            endif
         case (O_I2)
            call leaves%check_dirty(this%old(p0)%i_hist, "history0")
            call leaves%check_dirty(this%old(p1)%i_hist, "history1")
            call leaves%check_dirty(this%old(p2)%i_hist, "history2")
            dt_fac(:) = (t - this%old([ p0, p1, p2 ])%time) / (this%old([ p1, p2, p0 ])%time - this%old([ p2, p0, p1 ])%time)
            call leaves%q_lin_comb([ ind_val(this%old(p0)%i_hist, -dt_fac(2)*dt_fac(3)), &
                 &                   ind_val(this%old(p1)%i_hist, -dt_fac(1)*dt_fac(3)), &
                 &                   ind_val(this%old(p2)%i_hist, -dt_fac(1)*dt_fac(2)) ], solution )
         case default
            call die("[multigrid_gravity:init_solution] Extrapolation order not implemented")
      end select

      call leaves%check_dirty(solution, "init_soln")

      if (.not. this%valid) then
         do h = lbound(this%old, dim=1), ubound(this%old, dim=1)
            call leaves%set_q_value(this%old(h)%i_hist, 0.) ! set sane values for history to prevent spurious dirty exceptions in prolongation
         enddo
      endif

   end subroutine init_solution

!> \brief This routine manages old copies of potential for recycling.

   subroutine store_solution(this)

      use cg_leaves,     only: leaves
      use constants,     only: BND_XTRAP, BND_REF
      use domain,        only: dom
      use global,        only: t
      use multigridvars, only: solution, grav_bnd, bnd_isolated, bnd_givenval

      implicit none

      class(soln_history), intent(inout) :: this !< inner or outer potential history to store recent solution

      if (size(this%old) <= 0) return

      if (this%valid) then
         this%last = 1 + mod(this%last, size(this%old))
      else
         this%old(:)%time = t ! prevents extrapolation too early
      endif

      call leaves%q_copy(solution, this%old(this%last)%i_hist)
      this%old(this%last)%time = t
      this%valid = .true.

      ! Update guardcells of the solution before leaving. This can be done in higher-level routines that collect all the gravity contributions, but would be less safe.
      ! Extrapolate isolated boundaries, remember that grav_bnd is messed up by multigrid_solve_*
      if (grav_bnd == bnd_isolated .or. grav_bnd == bnd_givenval) then
         call leaves%arr3d_boundaries(solution, nb = dom%nb, bnd_type = BND_XTRAP)
      else
         call leaves%arr3d_boundaries(solution, nb = dom%nb, bnd_type = BND_REF)
      endif

   end subroutine store_solution

end module multigrid_old_soln
