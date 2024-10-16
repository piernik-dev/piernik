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

!> \brief Module for providing event list for precise parallel profiling data

module ppp_eventlist

   use constants,  only: cbuff_len, cwdlen
   use ppp_events, only: eventarray

   implicit none

   private
   public :: eventlist, use_profiling, disable_mask, profile_file, profile_lun

   ! namelist parameters
   logical :: use_profiling  !< control whether to do any PPProfiling or not

   character(len=cwdlen) :: profile_file  !< file name for the profile data
   integer :: profile_lun                 !< logical unit number for profile file
   integer(kind=4) :: disable_mask        !< logical mask for disabled events

   integer, parameter :: ev_arr_num = 10  !< number of event arrays allowed by default (bigger profiles may be difficult to handle)
   integer, parameter :: insane_arr_num = 2*ev_arr_num  !< as every next array is 2 times bigger, this allows for collecting 1024 times more events at the expense of RAM

   !> \brief list of events based on arrays of events, cheap to expand, avoid reallocation
   type eventlist
      private
      character(len=cbuff_len) :: label  !< label used to identify the event list
      type(eventarray), dimension(insane_arr_num) :: arrays  !< separate arrays to avoid lhs-reallocation
      integer :: arr_ind                !< currently used array
      integer :: ind                    !< first unused entry in currently used array
      logical :: overflown              !< emergency flag to turn off collecting more events
      logical :: xxl                    !< ignore ev_arr_num and use insane_arr_num instead (at your own risk)
   contains
      procedure :: init                 !< create new event list
      procedure :: cleanup              !< destroy this event list (typically called by publish)
      procedure :: start                !< add a beginning of an interval
      procedure :: stop                 !< add an end of an interval
      procedure :: single_cg_cost       !< add a cg-related interval
      procedure :: set_bb               !< add the initial event with bigbang time
      procedure, private :: next_event  !< for internal use in start, stop and put
      procedure, private :: expand      !< create next array for events
      procedure :: publish              !< write the collected data to a log file
   end type eventlist

contains

!> \brief Create new event list

   subroutine init(this, label, xxl)

      use dataio_pub, only: warn
      use mpisetup,   only: master

      implicit none

      class(eventlist), intent(inout) :: this   !< an object invoking the type-bound procedure
      character(len=*), intent(in)    :: label  !< event list label
      logical,          intent(in)    :: xxl    !< if .true. then ignore ev_arr_num until insane_arr_num is reached (risk of OOM)

      integer, parameter :: ev_arr_len = 1024  ! starting size of the array of events

      if (.not. use_profiling) return

      this%xxl = xxl
      this%overflown = .false.
      this%ind = 1
      this%arr_ind = 1
      this%label = label(1:min(cbuff_len, len_trim(label, kind=4)))
      call this%arrays(this%arr_ind)%init(ev_arr_len)

      if (this%xxl .and. master) call warn("[ppp_eventlist:init] XXL profiling enabled")
      ! To prevent running into OOM consider setting max_mem from MEMORY namelist.

   end subroutine init

!> \brief Destroy this event list

   subroutine cleanup(this)

      use constants, only: INVALID

      implicit none

      class(eventlist), intent(inout) :: this   !< an object invoking the type-bound procedure

      integer :: i

      do i = lbound(this%arrays, dim=1), ubound(this%arrays, dim=1)
         if (allocated(this%arrays(i)%ev_arr)) call this%arrays(i)%cleanup
      enddo

      this%ind = INVALID
      this%arr_ind = INVALID

   end subroutine cleanup

!> \brief Add a beginning of an interval

   subroutine start(this, label, mask)

      use MPIF,       only: MPI_Wtime
      use mpisetup,   only: bigbang_shift
      use ppp_events, only: event

      implicit none

      class(eventlist),          intent(inout) :: this   !< an object invoking the type-bound procedure
      character(len=*),          intent(in)    :: label  !< event label
      integer(kind=4), optional, intent(in)    :: mask   !< event category, if provided

      character(len=cbuff_len) :: l

      if (.not. use_profiling) return
      if (present(mask)) then
         if (iand(mask, disable_mask) /= 0) return
      endif

      l = label(1:min(cbuff_len, len_trim(label, kind=4)))
      call this%next_event(event(l, MPI_Wtime() + bigbang_shift))

   end subroutine start

!> \brief Add an end of an interval, use negative sign

   subroutine stop(this, label, mask)

      use MPIF,       only: MPI_Wtime
      use mpisetup,   only: bigbang_shift
      use ppp_events, only: event

      implicit none

      class(eventlist),          intent(inout) :: this   !< an object invoking the type-bound procedure
      character(len=*),          intent(in)    :: label  !< event label
      integer(kind=4), optional, intent(in)    :: mask   !< event category, if provided, should match the category provided in this%start call

      character(len=cbuff_len) :: l

      if (.not. use_profiling) return
      if (present(mask)) then
         if (iand(mask, disable_mask) /= 0) return
      endif

      l = label(1:min(cbuff_len, len_trim(label, kind=4)))
      call this%next_event(event(l, -MPI_Wtime() - bigbang_shift))

   end subroutine stop

!> \brief Add a cg-related interval

   subroutine single_cg_cost(this, t_start, t_stop, label)

      use constants,  only: PPP_CG
      use mpisetup,   only: bigbang_shift
      use ppp_events, only: event

      implicit none

      class(eventlist), intent(inout) :: this     !< an object invoking the type-bound procedure
      real,             intent(in)    :: t_start  !< start of the interval
      real,             intent(in)    :: t_stop   !< stop of the interval
      character(len=*), intent(in)    :: label    !< event label

      character(len=cbuff_len) :: l

      if (.not. use_profiling) return
      if (iand(PPP_CG, disable_mask) /= 0) return
      l = label(1:min(cbuff_len, len_trim(label, kind=4)))

      call this%next_event(event(l, t_start + bigbang_shift))
      call this%next_event(event(l, -t_stop - bigbang_shift))

   end subroutine single_cg_cost

!> \brief Add the initial event with bigbang time

   subroutine set_bb(this, label)

      use mpisetup,   only: bigbang, bigbang_shift
      use ppp_events, only: event

      implicit none

      class(eventlist), intent(inout) :: this   !< an object invoking the type-bound procedure
      character(len=*), intent(in)    :: label  !< event label

      character(len=cbuff_len) :: l

      if (.not. use_profiling) return

      l = label(1:min(cbuff_len, len_trim(label, kind=4)))
      call this%next_event(event(l, bigbang + bigbang_shift))

   end subroutine set_bb

!> \brief Start, stop and put

   subroutine next_event(this, ev)

      use dataio_pub, only: warn
      use ppp_events, only: event

      implicit none

      class(eventlist), intent(inout) :: this   !< an object invoking the type-bound procedure
      type(event),      intent(in)    :: ev     !< an event to be stored

      if (this%overflown) return

      this%arrays(this%arr_ind)%ev_arr(this%ind) = ev
      this%ind = this%ind + 1
      if (this%ind > ubound(this%arrays(this%arr_ind)%ev_arr, dim=1)) then
         if ((this%arr_ind >= ev_arr_num) .and. .not. this%xxl) then
            call warn("[ppp_eventlist:next_event] Run out of space allowed by ev_arr_num on '" // trim(ev%label) // "'")
            this%overflown = .true.
         else
            if (this%arr_ind >= insane_arr_num) then
               call warn("[ppp_eventlist:next_event] Run out of space allowed by insane_arr_num on '" // trim(ev%label) // "'")
               this%overflown = .true.
            else
               call this%expand
            endif
         endif
      else
         this%arrays(this%arr_ind)%ev_arr(this%ind)%wtime = 0.
      endif

   end subroutine next_event

!>
!! \brief Create next array for events
!!
!! \details Adding arrays should be cheaper than resizing existing ones.
!! The size of new array is double of the size of previous one to
!! prevent frequent, fragmented allocations on busy counters.
!<

   subroutine expand(this)

      use dataio_pub, only: warn

      implicit none

      class(eventlist), intent(inout) :: this   !< an object invoking the type-bound procedure

      integer, parameter :: grow_factor = 2  !< each next array should be bigger

      if (this%arr_ind >= ev_arr_num) &
           call warn("[ppp_eventlist:expand] Adding arrays beyond ev_arr_num limit for '" // &
           &         trim(this%arrays(this%arr_ind)%ev_arr(this%ind-1)%label) // "'")

      call this%arrays(this%arr_ind + 1)%init(grow_factor*size(this%arrays(this%arr_ind)%ev_arr))

      this%arr_ind = this%arr_ind + 1
      this%ind = 1
      this%arrays(this%arr_ind)%ev_arr(this%ind)%wtime = 0.

   end subroutine expand

!>
!! \brief Write the collected data to a log file and clear the log
!!
!! Perhaps MPI_TYPE_CREATE_STRUCT would simplify the code and improve the communication
!! but it requires some C-interoperability, which needs to be explored and tested first.
!!
!! \todo Consider HDF5, XML or JSON output but this would require adding proper support in ppp_plot.py
!<

   subroutine publish(this)

      use barrier,      only: piernik_MPI_Barrier
      use constants,    only: I_ZERO, I_ONE, V_INFO
      use dataio_pub,   only: warn, printinfo, msg
      use isend_irecv,  only: piernik_Isend
      use MPIF,         only: MPI_STATUS_IGNORE, MPI_CHARACTER, MPI_INTEGER, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD
      use MPIFUN,       only: MPI_Recv
      use mpisetup,     only: proc, master, slave, err_mpi, FIRST, LAST
      use req_array,    only: req_arr

      implicit none

      class(eventlist), intent(inout) :: this   !< an object invoking the type-bound procedure

      type(req_arr) :: req
      integer(kind=4) :: ne, p
      integer :: ia
      character(len=cbuff_len), dimension(:), allocatable  :: buflabel
      real(kind=8), dimension(:), allocatable :: buftime
      enum, bind(C)
         enumerator :: TAG_CNT = 1, TAG_ARR_L, TAG_ARR_T
      end enum

      if (.not. use_profiling) return

      if (this%overflown) then
         call warn("[ppp_eventlist:publish] Profile timings have overflown the allowed buffers and thus are partially broken. Skipping.")
         call this%cleanup
         return
      endif

      if (master) then
         if (disable_mask /= 0) then
            write(msg, '(a,b13)')"event disable mask = ", disable_mask
         else
            msg = "all events categories are enabled"
         endif

         call printinfo("[ppp_eventlist:publish] Profile timings will be written to '" // trim(profile_file) // "' file, " // trim(msg), V_INFO)
         open(newunit=profile_lun, file=profile_file)
      endif

      ! send
      ne = I_ZERO
      do ia = lbound(this%arrays, dim=1), ubound(this%arrays, dim=1)
         if (allocated(this%arrays(ia)%ev_arr)) then
            ne = ne + size(this%arrays(ia)%ev_arr, kind=4)
         else
            exit
         endif
      enddo
      if (slave) then
         call req%init(TAG_ARR_T, owncomm = .false., label = "ppp_ev")
         call piernik_Isend(ne, I_ONE, MPI_INTEGER, FIRST, TAG_CNT, req)
      endif

      if (ne > 0) then
         allocate(buflabel(ne), buftime(ne))
         p = I_ONE
         do ia = lbound(this%arrays, dim=1), ubound(this%arrays, dim=1)
            if (allocated(this%arrays(ia)%ev_arr)) then
               buflabel(p:p+size(this%arrays(ia)%ev_arr)-I_ONE) = this%arrays(ia)%ev_arr(:)%label
               buftime (p:p+size(this%arrays(ia)%ev_arr)-I_ONE) = this%arrays(ia)%ev_arr(:)%wtime
               p = p + size(this%arrays(ia)%ev_arr, kind=4)
            endif
         enddo
         if (master) then
            call publish_buffers(proc, buflabel, buftime)
            deallocate(buflabel, buftime)
         else
            call piernik_Isend(buflabel, size(buflabel, kind=4)*len(buflabel(1), kind=4), MPI_CHARACTER,        FIRST, TAG_ARR_L, req)
            call piernik_Isend(buftime,  size(buftime,  kind=4),                          MPI_DOUBLE_PRECISION, FIRST, TAG_ARR_T, req)
         endif
      endif

      if (master) then
         ! write(profile_lun, '(/,3a)') "#profile '", trim(this%label), "'"

         ! receive
         do p = FIRST + I_ONE, LAST
            call MPI_Recv(ne, I_ONE, MPI_INTEGER, p, TAG_CNT, MPI_COMM_WORLD, MPI_STATUS_IGNORE, err_mpi)
            if (ne /= 0) then
               allocate(buflabel(ne), buftime(ne))
               call MPI_Recv(buflabel, size(buflabel, kind=4)*len(buflabel(1), kind=4), MPI_CHARACTER,        p, TAG_ARR_L, MPI_COMM_WORLD, MPI_STATUS_IGNORE, err_mpi)
               call MPI_Recv(buftime,  size(buftime, kind=4),                           MPI_DOUBLE_PRECISION, p, TAG_ARR_T, MPI_COMM_WORLD, MPI_STATUS_IGNORE, err_mpi)
               call publish_buffers(p, buflabel, buftime)
               deallocate(buflabel, buftime)
            endif
         enddo
      else
         call req%waitall
         deallocate(buflabel, buftime)
      endif

      call piernik_MPI_Barrier

      call this%cleanup

   end subroutine publish

!> \brief Print the events from buffer arrays to the profile log

   subroutine publish_buffers(process, ev_label, ev_time)

      use dataio_pub, only: warn, die
      use func,       only: operator(.equals.)
      use mpisetup,   only: slave

      implicit none

      integer(kind=4),                        intent(in) :: process   !< origin of the event
      character(len=cbuff_len), dimension(:), intent(in) :: ev_label  !< array of event labels
      real(kind=8), dimension(:),             intent(in) :: ev_time   !< array of event times

      integer :: i, d

      if (slave) then
         call warn("[ppp_eventlist:publish_array] only master is supposed to write")
         return
      endif

      if (size(ev_label) /= size(ev_time)) call die("[ppp_eventlist:publish_array] arrays size mismatch")

      d = 1
      do i = lbound(ev_label, dim=1), ubound(ev_label, dim=1)
         if (ev_time(i) .equals. 0.) exit
         if (ev_time(i) < 0.) d = d - 1
         write(profile_lun, '(i5,a,f20.7,3a)') process, " ", ev_time(i), " ", repeat("  ", d), trim(ev_label(i))
         if (ev_time(i) > 0.) d = d + 1
      enddo
      flush(profile_lun)

   end subroutine publish_buffers

end module ppp_eventlist
