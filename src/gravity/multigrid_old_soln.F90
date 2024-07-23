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
!! Minimum restart version 2.02

module multigrid_old_soln
! pulled by MULTIGRID && SELF_GRAV

   use old_soln_list, only: old_soln, os_list_undef_t, os_list_t

   implicit none

   private
   public :: nold_max, soln_history, ord_time_extrap

   ! solution recycling
   integer(kind=4), parameter :: nold_max=3   !< maximum implemented extrapolation order

   type :: soln_history                       !< container for a set of several old potential solutions
      type(os_list_undef_t) :: invalid        !< a list of invalid slots ready to use
      type(os_list_t) :: old                  !< indices and time points of stored solutions
   contains
      procedure :: init_history               !< Allocate arrays, register fields
      procedure :: cleanup_history            !< Deallocate arrays
      procedure :: init_solution              !< Construct first guess of potential based on previously obtained solution, if any.
      procedure :: store_solution             !< Manage old copies of potential for recycling.
      procedure :: sanitize                   !< Invalidate some stored solutions from the future i.e. when there was timestep retry
      procedure :: sanitize_expanded          !< If the domain was recently expanded, initialize all history with zeroes
      procedure :: print                      !< Print the state of old solution list
      procedure :: unmark                     !< Reset restart flag of old soln
#ifdef HDF5
      procedure :: mark_and_create_attribute  !< Mark some old solutions for restarts and set up necessary attributes
      procedure :: read_os_attribute          !< Read old solutions identifiers, their times, and initialize history
#endif /* HDF5 */
   end type soln_history

   ! Namelist parameter
   integer(kind=4)    :: ord_time_extrap      !< Order of temporal extrapolation for solution recycling; -1 means 0-guess, 2 does parabolic interpolation

contains

!> \brief Initialize structure for keeping historical potential fields

   subroutine init_history(this, nold, prefix)

      use cg_list_global,   only: all_cg
      use constants,        only: singlechar, dsetnamelen, I_TWO
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
      do i = 1, nold + I_TWO  ! extra two slots for seamless timestep retries
         write(hname,'(2a,i2.2)')prefix,"-h-",i
         call all_cg%reg_var(hname, vital = .false., ord_prolong = ord_prolong) ! no need for multigrid attribute here because history is defined only on leaves
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
!! Rational extrapolation (1 + a t)/(b + c t) can give 2 times better or 10 times worse guess depending on timestep.
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
      use constants,      only: INVALID, O_INJ, O_LIN, O_I2, I_ONE, dirtyH1, V_DEBUG, V_VERBOSE
      use dataio_pub,     only: msg, die, printinfo
      use global,         only: t
      use mpisetup,       only: master
      use multigridvars,  only: solution

      implicit none

      class(soln_history), intent(inout) :: this    !< inner or outer potential history for recycling
      character(len=*),    intent(in)    :: prefix  !< informational string to be printed

      integer :: ordt
      real, dimension(3) :: dt_fac

      ! BEWARE selfgrav_clump/initproblem.F90 requires monotonic time sequence t > this%old(p0)%time > this%old(p1)%time > this%old(p2)%time

      call all_cg%set_dirty(solution, 0.98*dirtyH1)

      call this%sanitize

      ordt = min(this%old%cnt() - I_ONE, ord_time_extrap)

      write(msg, '(a,g14.6,a,i2)')"[multigrid_old_soln:init_solution] init solution at time: ", t, " order: ", ordt
      if (master) call printinfo(msg, V_DEBUG)
#ifdef DEBUG
      call this%print
#endif /* DEBUG */

      select case (ordt)
         case (:INVALID)
            if (master .and. ord_time_extrap > ordt) then
               write(msg, '(3a)')"[multigrid_old_soln:init_solution] Clearing ", trim(prefix), "solution."
               call printinfo(msg, V_VERBOSE)
            endif
            call all_cg%set_q_value(solution, 0.)
         case (O_INJ)
            call leaves%check_dirty(this%old%latest%i_hist, "history0")
            call leaves%q_copy(this%old%latest%i_hist, solution)
            if (master .and. ord_time_extrap > ordt) then
               write(msg, '(3a)')"[multigrid_old_soln:init_solution] No extrapolation of ",trim(prefix),"solution."
               call printinfo(msg, V_VERBOSE)
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
                  call printinfo(msg, V_VERBOSE)
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

      use cg_leaves,        only: leaves
      use constants,        only: V_DEBUG
      use dataio_pub,       only: die, msg, printinfo
      use global,           only: t
      use mpisetup,         only: master
      use multigridvars,    only: solution
      use named_array_list, only: qna
      use old_soln_list,    only: old_soln

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
      qna%lst(os%i_hist)%vital = .true.

      write(msg, '(a,g14.6)')"[multigrid_old_soln:store_solution] store solution at time ", t
      if (master) call printinfo(msg, V_DEBUG)
#ifdef DEBUG
      call this%print
#endif /* DEBUG */

   end subroutine store_solution
!> \brief Invalidate some stored solutions from the future i.e. when there was timestep retry

   subroutine sanitize(this)

      use constants,        only: V_VERBOSE
      use dataio_pub,       only: printinfo, msg
      use global,           only: t
      use mpisetup,         only: master
      use named_array_list, only: qna

      implicit none

      class(soln_history), intent(inout) :: this !< potential history to be sanitized

      type(old_soln), pointer :: os
      integer :: cnt

      cnt = this%old%cnt()
      do while (associated(this%old%latest))
         if (this%old%latest%time >= t) then
            os => this%old%pick_head()
            call this%invalid%new_head(os)
            qna%lst(os%i_hist)%vital = .false.
         else
            exit
         endif
      enddo

      if (cnt /= this%old%cnt()) then
         write(msg, '(a,g14.6,a,i2,a)')"[multigrid_old_soln:sanitize] sanitize solution history at time ", t, ", removed ", cnt - this%old%cnt(), " elements"
         if (master) call printinfo(msg, V_VERBOSE)
#ifdef DEBUG
         call this%print
#endif /* DEBUG */
      endif

   end subroutine sanitize

!> \brief If the domain was recently expanded, initialize history with zeroes

   subroutine sanitize_expanded(this)

      use cg_list_dataop, only: expanded_domain

      implicit none

      class(soln_history), intent(inout) :: this !< potential history to be sanitized after domain expansion

      type(old_soln), pointer :: os

      os => this%old%latest
      do while (associated(os))
         call expanded_domain%set_q_value(os%i_hist, 0.)
         os => os%earlier
      enddo

   end subroutine sanitize_expanded

!> \brief Print the state of old solution list

   subroutine print(this)

      implicit none

      class(soln_history), intent(in) :: this !< potential history to be printed

      call this%old%print
      call this%invalid%print

   end subroutine print
!>
!! \brief Reset restart flag of old soln
!!
!! This is required to avoid creating .retry field for copies in case of timestep retry
!!
!! This routine is safe to be called on uninitialized old_soln (.not. associated(this%old%latest))
!<

   subroutine unmark(this)

      use constants,        only: AT_IGNORE
      use named_array_list, only: qna

      implicit none

      class(soln_history), intent(inout) :: this !< potential history to be registered for restarts

      type(old_soln), pointer :: os

      os => this%old%latest
      do while (associated(os))
         qna%lst(os%i_hist)%restart_mode = AT_IGNORE
         os => os%earlier
      enddo

   end subroutine unmark

#ifdef HDF5
!>
!! \brief Mark some old solutions for restarts and set up necessary attributes
!!
!! This routine needs to be called before the datasets are written (before call write_restart_hdf5_v2).
!!
!! This routine is safe to be called on uninitialized old_soln (this%old%cnt() <= 0)
!<
   subroutine mark_and_create_attribute(this, file_id)

      use constants,          only: I_ONE, AT_IGNORE, AT_NO_B, cbuff_len, I_ONE, I_TWO
      use hdf5,               only: HID_T
      use named_array_list,   only: qna
      use mpisetup,           only: master
      use old_soln_list,      only: old_soln
      use set_get_attributes, only: set_attr

      implicit none

      class(soln_history), intent(in) :: this !< potential history to be registered for restarts
      integer(HID_T),      intent(in) :: file_id  !< File identifier

      integer(kind=4) :: n, i, b, found
      type(old_soln), pointer :: os
      character(len=cbuff_len), allocatable, dimension(:) :: namelist
      real, allocatable, dimension(:) :: timelist

      n = max(min(this%old%cnt(), ord_time_extrap + I_ONE), I_TWO)  ! try to save at least 2 points to recover also sgpm

      allocate(namelist(n), timelist(n))

      ! set the flags to mark which fields should go to the restart
      i = 1
      found = 0
      os => this%old%latest
      do while (associated(os))
         b = AT_IGNORE
         if (i <= n) then
            b = AT_NO_B
            namelist(i) = qna%lst(os%i_hist)%name
            timelist(i) = os%time
            found = found + I_ONE
         endif
         qna%lst(os%i_hist)%restart_mode = b
         i = i + I_ONE
         os => os%earlier
      enddo

      if (master .and. found > 0) then
         call set_attr(file_id, trim(this%old%label) // "_names", namelist(:found))
         call set_attr(file_id, trim(this%old%label) // "_times", timelist(:found))
      endif

      deallocate(namelist)
      deallocate(timelist)

   end subroutine mark_and_create_attribute

!>
!! \brief Read old solutions identifiers, their times, and initialize history
!!
!! This routine needs to be called before the datasets are read.
!!
!! Unlike mark_and_create_attribute and unmark this routine is NOT safe to be
!! called on uninitialized old_soln (non-fatal errors will occur).
!<

   subroutine read_os_attribute(this, file_id)

      use constants,          only: cbuff_len, AT_NO_B
      use dataio_pub,         only: msg, die, printio
      use hdf5,               only: HID_T, HSIZE_T, SIZE_T, &
           &                        h5aexists_f, h5gopen_f, h5gclose_f
      use h5lt,               only: h5ltget_attribute_ndims_f, h5ltget_attribute_info_f
      use mpisetup,           only: master
      use named_array_list,   only: qna
      use set_get_attributes, only: get_attr

      implicit none

      class(soln_history), intent(inout) :: this !< potential history to be registered for restarts
      integer(HID_T),      intent(in)    :: file_id  !< File identifier

      integer(kind=4) :: rank, error
      integer(HSIZE_T), dimension(1) :: dims
      integer(kind=4) :: tclass
      integer(SIZE_T) :: tsize
      character(len=cbuff_len), allocatable, dimension(:) :: namelist
      real, allocatable, dimension(:) :: timelist
      character(len=*), parameter, dimension(2) :: nt = [ "_names", "_times" ]
      logical(kind=4) :: a_exists
      integer(HID_T) :: g_id
      integer :: i
      type(old_soln), pointer :: os

      call h5gopen_f(file_id, "/", g_id, error)
      do i = lbound(nt, 1), ubound(nt, 1)
         call h5aexists_f(g_id, trim(this%old%label) // nt(i), a_exists, error)
         if (.not. a_exists) then
            if (master) call printio("[multigrid_old_soln:read_os_attribute] " // trim(this%old%label) // nt(i) // " does not exist. Coldboot")
            return
         endif
      enddo
      call h5gclose_f(g_id, error)

      call h5ltget_attribute_ndims_f(file_id, "/", trim(this%old%label) // "_names", rank, error)
      if (error /= 0 .or. rank /= 1) then
         write(msg, '(2(a,i4))')"[multigrid_old_soln:read_os_attribute] " // trim(this%old%label) // "_names: need rank=1, got ", rank, " , error = ", error
         call die(msg)
      endif
      call h5ltget_attribute_ndims_f(file_id, "/", trim(this%old%label) // "_times", rank, error)
      if (error /= 0 .or. rank /= 1) then
         write(msg, '(2(a,i4))')"[multigrid_old_soln:read_os_attribute] " // trim(this%old%label) // "_times: need rank=1, got ", rank, " , error = ", error
         call die(msg)
      endif
      call h5ltget_attribute_info_f(file_id, "/", trim(this%old%label) // "_names", dims, tclass, tsize, error)
      if (dims(1) <= 0) then
         if (master) call printio("[multigrid_old_soln:read_os_attribute] No " // trim(this%old%label) // "_names to read. Coldboot.")
         return
      endif
      allocate(namelist(dims(1)))
      call get_attr(file_id, trim(this%old%label) // "_names", namelist)
      call h5ltget_attribute_info_f(file_id, "/", trim(this%old%label) // "_times", dims, tclass, tsize, error)
      if (dims(1) /= size(namelist)) call die("[multigrid_old_soln:read_os_attribute] size("// trim(this%old%label) // "_names) /= size("// trim(this%old%label) // "_times")
      allocate(timelist(dims(1)))
      call get_attr(file_id, trim(this%old%label) // "_times", timelist)

      if (associated(this%old%latest)) call die("[multigrid_old_soln:read_os_attribute] " // trim(this%old%label) // "nonempty")
      do i = ubound(namelist, 1), lbound(namelist, 1), -1  ! it is stored with most recent entries first
         if (qna%exists(trim(namelist(i)))) then
            os => this%invalid%pick(qna%ind(namelist(i)))
            if (associated(os)) then
               call this%old%new_head(os)
               os%time = timelist(i)  ! new_head did put current time, need to enforce it
               qna%lst(os%i_hist)%restart_mode = AT_NO_B
               qna%lst(os%i_hist)%vital = .true.
            else
               call die("[multigrid_old_soln:read_os_attribute] Cannot find '" // trim(namelist(i)) // "' in free slots.")
            endif
         else
            call die("[multigrid_old_soln:read_os_attribute] Cannot find '" // trim(namelist(i)) // "' in qna.")
         endif

      enddo

   end subroutine read_os_attribute

#endif /* HDF5 */


end module multigrid_old_soln
