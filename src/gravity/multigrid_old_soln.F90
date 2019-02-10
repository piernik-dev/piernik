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

!> \brief Multigrid historical solutions

module multigrid_old_soln
! pulled by MULTIGRID && GRAV

   use old_soln_list, only: old_soln, os_list_undef_T, os_list_T

   implicit none

   private
   public :: nold_max, soln_history, ord_time_extrap

   ! solution recycling
   integer(kind=4), parameter :: nold_max=3   !< maximum implemented extrapolation order

   type :: soln_history                       !< container for a set of several old potential solutions
      type(os_list_undef_T) :: invalid        !< a list of invalid slots ready to use
      type(os_list_T) :: old                  !< indices and time points of stored solutions
    contains
      procedure :: init_history               !< Allocate arrays, register fields
      procedure :: cleanup_history            !< Deallocate arrays
      procedure :: init_solution              !< Construct first guess of potential based on previously obtained solution, if any.
      procedure :: store_solution             !< Manage old copies of potential for recycling.
      procedure :: sanitize                   !< invalidate some stored solutions from the future i.e. when there was timestep retry
   end type soln_history

   ! Namelist parameter
   integer(kind=4)    :: ord_time_extrap      !< Order of temporal extrapolation for solution recycling; -1 means 0-guess, 2 does parabolic interpolation

contains

!> \brief Initialize structure for keeping historical potential fields

   subroutine init_history(this, nold, prefix)

      use cg_list_global,   only: all_cg
      use constants,        only: singlechar, dsetnamelen
      use dataio_pub,       only: die, msg
      use multigridvars,    only: ord_prolong
      use named_array_list, only: qna

      implicit none

      class(soln_history),       intent(inout) :: this   !< Object invoking type-bound procedure
      integer,                   intent(in)    :: nold   !< how many historic points to store
      character(len=singlechar), intent(in)    :: prefix !< 'i' for inner, 'o' for outer potential

      integer :: i
      character(len=dsetnamelen) :: hname

      if (associated(this%old%latest) .or. associated(this%invalid%latest)) then
         write(msg, '(3a)') "[multigrid_old_soln:init_history] ", prefix," already initialized."
         call die(msg)
      endif
      do i = 1, nold
         write(hname,'(2a,i2.2)')prefix,"-h-",i
         call all_cg%reg_var(hname, vital = .true., ord_prolong = ord_prolong) ! no need for multigrid attribute here because history is defined only on leaves
         call this%invalid%new(qna%ind(hname))
      enddo

      write(msg, '(2a)') prefix, "-invalid"
      this%invalid%label = trim(msg)

      write(msg, '(2a)') prefix, "-old"
      this%old%label = trim(msg)

   end subroutine init_history

!> \brief Deallocate history arrays

   subroutine cleanup_history(this)

      implicit none

      class(soln_history), intent(inout) :: this   !< Object invoking type-bound procedure

      call this%old%cleanup
      call this%invalid%cleanup

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
      use constants,      only: INVALID, O_INJ, O_LIN, O_I2, I_ONE
      use dataio_pub,     only: msg, die, printinfo
      use global,         only: t
      use mpisetup,       only: master
      use multigridvars,  only: stdout, solution

      implicit none

      class(soln_history), intent(inout) :: this    !< inner or outer potential history for recycling
      character(len=*),    intent(in)    :: prefix  !< informational string to be printed

      integer :: ordt
      real, dimension(3) :: dt_fac

      ! BEWARE selfgrav_clump/initproblem.F90 requires monotonic time sequence t > this%old(p0)%time > this%old(p1)%time > this%old(p2)%time

      call all_cg%set_dirty(solution)

      call this%sanitize

      ordt = min(this%old%cnt() - I_ONE, ord_time_extrap)

#ifdef DEBUG
      write(msg, '(a,g14.6,a,i2)')"[multigrid_old_soln:init_solution] init solution at time: ", t, " order: ", ordt
      if (master) call printinfo(msg)
      call this%old%print
      call this%invalid%print
#endif /* DEBUG */

      select case (ordt)
         case (:INVALID)
            if (master .and. ord_time_extrap > ordt) then
               write(msg, '(3a)')"[multigrid_old_soln:init_solution] Clearing ", trim(prefix), "solution."
               call printinfo(msg, stdout)
            endif
            call all_cg%set_q_value(solution, 0.)
            if (associated(this%old%latest)) call die("[multigrid_old_soln:init_solution] need to move %old to %invalid")
         case (O_INJ)
            call leaves%check_dirty(this%old%latest%i_hist, "history0")
            call leaves%q_copy(this%old%latest%i_hist, solution)
            if (master .and. ord_time_extrap > ordt) then
               write(msg, '(3a)')"[multigrid_old_soln:init_solution] No extrapolation of ",trim(prefix),"solution."
               call printinfo(msg, stdout)
            endif
         case (O_LIN)
            associate( p0 => this%old%latest, &
                 &     p1 => this%old%latest%earlier )
               call leaves%check_dirty(p0%i_hist, "history0")
               call leaves%check_dirty(p1%i_hist, "history1")
               dt_fac(1) = (t - p0%time) / (p0%time - p1%time)
               call leaves%q_lin_comb( [ ind_val(p0%i_hist, (1.+dt_fac(1))), &
                    &                    ind_val(p1%i_hist,    -dt_fac(1) ) ], solution )
               if (master .and. ord_time_extrap > ordt) then
                  write(msg, '(3a)')"[multigrid_old_soln:init_solution] Linear extrapolation of ",trim(prefix),"solution."
                  call printinfo(msg, stdout)
               endif
            end associate
         case (O_I2:)
            associate( p0 => this%old%latest, &
                 &     p1 => this%old%latest%earlier, &
                 &     p2 => this%old%latest%earlier%earlier )
               call leaves%check_dirty(p0%i_hist, "history0")
               call leaves%check_dirty(p1%i_hist, "history1")
               call leaves%check_dirty(p2%i_hist, "history2")
               dt_fac(:) = (t - [ p0%time, p1%time, p2%time ]) / ([ p1%time, p2%time, p0%time ] - [ p2%time, p0%time, p1%time ])
               call leaves%q_lin_comb([ ind_val(p0%i_hist, -dt_fac(2)*dt_fac(3)), &
                    &                   ind_val(p1%i_hist, -dt_fac(1)*dt_fac(3)), &
                    &                   ind_val(p2%i_hist, -dt_fac(1)*dt_fac(2)) ], solution )
            end associate
         case default
            call die("[multigrid_old_soln:init_solution] Extrapolation order not implemented")
      end select

      call leaves%check_dirty(solution, "init_soln")

   end subroutine init_solution

!> \brief This routine manages old copies of potential for recycling.

   subroutine store_solution(this)

      use cg_leaves,     only: leaves
      use dataio_pub,    only: die, msg
      use global,        only: t
      use multigridvars, only: solution
      use old_soln_list, only: old_soln
#ifdef DEBUG
      use dataio_pub,    only: printinfo
      use mpisetup,      only: master
#endif /* DEBUG */

      implicit none

      class(soln_history), intent(inout) :: this !< inner or outer potential history to store recent solution

      type(old_soln), pointer :: os

      if (.not. associated(this%old%latest) .and. .not. associated(this%invalid%latest)) return

      os => this%invalid%pick_head()
      if (.not. associated(os)) os => this%old%trim_tail()
      if (.not. associated(os)) then
         write(msg, '(4a)')"[multigrid_old_soln:store_solution] cannot get any slot from ", this%old%label, " or ", this%invalid%label
         call die(msg)
      endif

      os%time = t
      call this%old%new_head(os)
      call leaves%q_copy(solution, os%i_hist)

#ifdef DEBUG
      write(msg, '(a,g14.6)')"[multigrid_old_soln:store_solution] store solution at time ", t
      if (master) call printinfo(msg)
      call this%old%print
      call this%invalid%print
#endif /* DEBUG */

   end subroutine store_solution
!> \brief Invalidate some stored solutions from the future i.e. when there was timestep retry

   subroutine sanitize(this)

      use dataio_pub, only: printinfo, msg
      use global,     only: t
      use mpisetup,   only: master

      implicit none

      class(soln_history), intent(inout) :: this !< potential history to be sanitized

      type(old_soln), pointer :: os
      integer :: cnt

      cnt = this%old%cnt()
      do while (associated(this%old%latest))
         if (this%old%latest%time >= t) then
            os => this%old%pick_head()
            call this%invalid%new_head(os)
         else
            exit
         endif
      enddo

      if (cnt /= this%old%cnt()) then
         write(msg, '(a,g14.6,a,i2,a)')"[multigrid_old_soln:sanitize] sanitize solution history at time ", t, ", removed ", cnt - this%old%cnt(), " elements"
         if (master) call printinfo(msg)
#ifdef DEBUG
         call this%old%print
         call this%invalid%print
#endif /* DEBUG */
      endif

   end subroutine sanitize

end module multigrid_old_soln
