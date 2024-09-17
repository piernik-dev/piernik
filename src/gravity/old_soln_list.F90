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

!> \brief This module contains old solutions list and related methods

! pulled by MULTIGRID && SELF_GRAV

module old_soln_list

   use constants, only: dsetnamelen

   implicit none

   private
   public :: old_soln, os_list_undef_t, os_list_t

   real, parameter    :: invalid_time = -0.1 * huge(1.0)  !< don't trust too ancient solutions
   integer, parameter :: too_long = 100                   !< escape from loops

   !< container for an old solution with its timestamp
   type :: old_soln
      integer(kind=4)         :: i_hist   !< index to the old solution in named_array_list::qna
      real                    :: time     !< time of the old solution
      type(old_soln), pointer :: earlier  !< a pointer to earlier solution
      type(old_soln), pointer :: later    !< a pointer to later solution
   end type old_soln

   type, abstract :: os_list_at
      type(old_soln), pointer    :: latest    !< first (most recent) element of the chain of historical solutions
      character(len=dsetnamelen) :: label     !< name of the list for diagnostic and identification purposes
   contains
      procedure :: cleanup     !< deallocate memory
      procedure :: last        !< find the pointer to the last element
      procedure :: cnt         !< return length of the list
      procedure :: print       !< dump some info to stdout
      procedure :: pick_head   !< unlink head and return it to the caller
      procedure :: new_head    !< put an element to the front of the list
   end type os_list_at

   type, extends(os_list_at) :: os_list_undef_t
   contains
      procedure :: new         !< add a fresh element anywhere
      procedure :: pick        !< unlink an element that is matching i_hist
   end type os_list_undef_t

   type, extends(os_list_at) :: os_list_t
   contains
      procedure :: trim_tail   !< detach an element from the end of the list and return it to the caller
      procedure :: is_valid    !< can we trust this list?
   end type os_list_t

contains

!> \brief Deallocate memory

   subroutine cleanup(this)

      implicit none

      class(os_list_at), intent(inout) :: this

      type(old_soln), pointer :: os, prev

      if (.not. associated(this%latest)) return

      os => this%last()
      do while (associated(os))
         prev => os%later
         deallocate(os)
         os => prev
      enddo

   end subroutine cleanup

!> \brief Find the pointer to the last element or return null()

   function last(this) result(os)

      use constants,  only: I_ZERO, I_ONE
      use dataio_pub, only: die

      implicit none

      class(os_list_at), intent(inout) :: this

      type(old_soln), pointer :: os
      integer :: cnt

      if (.not. associated(this%latest)) then
         os => null()
         return
      endif

      cnt = I_ZERO
      os => this%latest
      do while (associated(os) .and. cnt <= too_long)
         cnt = cnt + I_ONE
         if (associated(os%earlier)) then
            os => os%earlier
         else
            exit  ! quit with os pointing to the last element
         endif
      enddo

      if (cnt > too_long) os => null()

      if (associated(os)) then
         if (associated(os%earlier)) call die("[old_soln_list:last] associated(os%earlier)")
      endif

   end function last

!> \brief Return length of the list

   function cnt(this)

      use constants,  only: I_ZERO, I_ONE, INVALID

      implicit none

      class(os_list_at), intent(in) :: this

      type(old_soln), pointer :: os
      integer(kind=4) :: cnt

      cnt = I_ZERO
      if (.not. associated(this%latest)) return

      os => this%latest
      do while (associated(os) .and. cnt <= too_long)
         cnt = cnt + I_ONE
         os => os%earlier
      enddo

      if (cnt > too_long) cnt = INVALID

   end function cnt

!> \brief Add a fresh element anywhere

   subroutine new(this, ind)

      implicit none

      class(os_list_undef_t), intent(inout) :: this
      integer(kind=4), intent(in) :: ind

      type(old_soln), pointer :: n

      allocate(n)
      n = old_soln(ind, invalid_time, this%latest, null())
      if (associated(this%latest)) this%latest%later => n
      this%latest => n

   end subroutine new

!> \brief Unlink head and return it to the caller

   function pick(this, ind) result(os)

      implicit none

      class(os_list_undef_t), intent(inout) :: this
      integer(kind=4),        intent(in)    :: ind

      type(old_soln), pointer :: os, ose, osl

      if (.not. associated(this%latest)) then
         os => null()
         return
      endif

      os => this%latest
      do while (associated(os))
         if (ind == os%i_hist) then
            ose => os%earlier
            osl => os%later
            if (associated(osl)) then
               osl%earlier => ose
            else
               this%latest => ose
            endif
            if (associated(ose)) ose%later => osl
            return
         endif
         os => os%earlier
      enddo

      os => null()

   end function pick

!> \brief Unlink head and return it to the caller

   function pick_head(this) result(os)

      implicit none

      class(os_list_at), intent(inout) :: this

      type(old_soln), pointer :: os

      if (.not. associated(this%latest)) then
         os => null()
         return
      endif

      os => this%latest
      this%latest => this%latest%earlier
      if (associated(this%latest)) this%latest%later => null()

   end function pick_head

!> \brief Put an element to the front of the list

   subroutine new_head(this, os)

      use dataio_pub, only: die
      use global,     only: t

      implicit none

      class(os_list_at),       intent(inout) :: this
      type(old_soln), pointer, intent(in)    :: os

      os%later => null()
      os%earlier => this%latest
      select type(this)
         type is (os_list_undef_t)
         type is (os_list_t)
            os%time = t
         class default
            call die("[old_soln_list:new_head] unknown type")
      end select
      this%latest => os
      if (associated(os%earlier)) os%earlier%later => os

   end subroutine new_head

!> \brief detach an element from the end of the list

   function trim_tail(this) result(os)

      use dataio_pub, only: warn, die
      use mpisetup,   only: master

      implicit none

      class(os_list_t), intent(inout) :: this

      type(old_soln), pointer :: os

      if (.not. associated(this%latest)) then
         if (master) call warn("[old_soln_list:trim_tail] no slot to pick")
         os => null()
         return
      endif

      os => this%last()

      if (associated(os%later)) os%later%earlier => null()
      if (associated(os%earlier)) call die("[old_soln_list:last] associated(last%earlier)")
      os%time = invalid_time

   end function trim_tail

!> \brief Check consistency of this list

   function is_valid(this)

      use constants,  only: I_ZERO, I_ONE
      use dataio_pub, only: warn, msglen
      use global,     only: t
      use mpisetup,   only: master

      implicit none

      class(os_list_t), intent(in) :: this

      logical :: is_valid

      type(old_soln), pointer :: os, prev
      real :: time
      integer :: cnt
      character(len=msglen) :: msg !< private buffer for messages, don't use dataio_pub::msg to allow is_valid() calls in other writes

      is_valid = .true.  ! treat empty list as consistent

      cnt = I_ZERO
      time = t
      os => this%latest

      if (associated(os)) then
         if (associated(os%later)) then
            is_valid = .false.
            write(msg, '(a)')"[old_soln_list:is_valid] associated(this%latest%later)"
            if (master) call warn(msg)
         endif
      endif

      do while (associated(os) .and. is_valid)
         cnt = cnt + I_ONE
         if ((.not. associated(os, this%latest) .and. os%time >= time) .or. &
              os%time > time .or. &
              os%time <= invalid_time) then
            is_valid = .false.
            write(msg, '(a,i3,3(a,g14.6))')"[old_soln_list:is_valid] cnt= ", cnt, " time= ", time, " os%time= ", os%time, " invalid_time= ", invalid_time
            if (master) call warn(msg)
         endif
         time = os%time
         prev => os
         os => os%earlier
         if (associated(os)) then
            if (.not. associated(os%later, prev)) then
               is_valid = .false.
               write(msg, '(a,i3,a,2i3)')"[old_soln_list:is_valid] cnt= ", cnt, " .not. associated(os%later, prev)", os%i_hist, prev%i_hist
               if (master) call warn(msg)
            endif
         endif
         if (cnt > too_long) then
            is_valid = .false.
            write(msg, '(2(a,i3))')"[old_soln_list:is_valid] cnt= ", cnt, " > too_long= ", too_long
            if (master) call warn(msg)
            os => null()
         endif
      enddo

   end function is_valid

!> \brief Dump some info to stdout

   subroutine print(this)

      use constants,        only: I_ZERO, I_ONE, V_VERBOSE
      use dataio_pub,       only: msg, printinfo
      use global,           only: t
      use mpisetup,         only: master
      use named_array_list, only: qna

      implicit none

      class(os_list_at), intent(in) :: this

      type(old_soln), pointer :: os
      integer :: cnt

      cnt = I_ZERO
      os => this%latest
      do while (associated(os))
         cnt = cnt + I_ONE
         select type(this)
            type is (os_list_undef_t)
               write(msg, '(2(a,i3),2a)') "(Undef) soln# ", cnt, " qna_index: ", os%i_hist, " qna_name: ", qna%lst(os%i_hist)%name
            type is (os_list_t)
               write(msg, '(a,i3,a,g14.6,a,i3,2a)') "(Old) soln# ", cnt, " time = ", os%time, " qna_index: ", os%i_hist, " qna_name: ", qna%lst(os%i_hist)%name
            class default
               write(msg, '(a,i3,a,g14.6,a,i3,2a)') "(Other ?) soln# ", cnt, " time = ", os%time, " qna_index: ", os%i_hist, " qna_name: ", qna%lst(os%i_hist)%name
         end select
         if (master) call printinfo(msg, V_VERBOSE)
         os => os%earlier
         if (cnt > too_long) os => null()
      enddo

      write(msg, '(a,g14.6,3a,i3,a)') "[old_soln_list] t= ", t, " name: '", trim(this%label), "' contains ", this%cnt(), " elements"
      select type(this)
         type is (os_list_t)
            write(msg(len_trim(msg)+1:), '(a,l2)') " is valid? ", this%is_valid()
         class default
      end select
      if (master) call printinfo(msg, V_VERBOSE)

   end subroutine print

end module old_soln_list
