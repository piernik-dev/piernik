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
!>
!! \brief (KK)
!<
module timer

   implicit none

   integer, parameter, private :: S_LEN = 30

   private
   public :: cleanup_timers, time_left, set_timer, timer_start, timer_stop

   type, private :: timer_info
      character(len=S_LEN) :: key
      logical              :: reset
      real(kind=8)         :: time
   end type timer_info

   type, private :: timer_list
      type(timer_node), pointer :: next => Null()
   end type timer_list

   type, private :: timer_node
      type(timer_info) :: info
      type(timer_list) :: node
   end type timer_node

   type(timer_list), target, private, save :: timer_root

   integer :: cpuhours, cpumins, cpusecs, wchours, wcmins, wcsecs
   real    :: zcps, cputot, cpuallp, wctot, cpu_start, cpu_stop
   integer, dimension(3) :: iarray
   real(kind=4), dimension(2) :: tarray
   integer(kind=8) :: clock_start, clock_end

contains

   !>
   !! \brief Wrapper for timer handling.
   !! Only this routine should be used outside timer module.
   !! USAGE: call set_timer(name,reset)
   !!    1) if first called with "name", create timer "name" and set it with current MPI_Wtime(), no output
   !!    2) if timer "name" exists, print "name" - current time, set "name" with current MPI_Wtime()
   !!    (optional) if reset is true suppress output, set "name" with current MPI_Wtime
   !<
   real function set_timer(str,reset)

      implicit none

      character(len=*), intent(in) :: str    !< name of the timer
      logical, intent(in), optional :: reset !< if true all output is suppressed, use for resetting timers
      type(timer_info) :: temp

      temp%key = trim(str)
      temp%reset = present(reset)
      call search_timer(temp)
      set_timer = temp%time

   end function set_timer

   function delete_timer(tp) result (item)

      implicit none

      type(timer_node), pointer :: tp, temp
      type(timer_info) :: item

      temp => tp
      item = tp%info
      tp => tp%node%next
      deallocate(temp)

   end function delete_timer

   subroutine cleanup_timers

#ifdef VERBOSE
      use dataio_pub,    only: msg, warn
#endif /* VERBOSE */

      implicit none

      type(timer_list), pointer :: tp
      type(timer_info)          :: item

      tp => timer_root
#ifdef VERBOSE
      call warn("")
#endif /* VERBOSE */
      do while (associated(tp%next))
         item = delete_timer(tp%next)
#ifdef VERBOSE
         write(msg,'(3a)') "[timer:cleanup_timers]: Timer ",item%key," deleted"
         call warn(msg)
#endif /* VERBOSE */
      enddo

   end subroutine cleanup_timers

   subroutine search_timer(item)

      implicit none

      type(timer_info), intent(inout) :: item
      type(timer_list), pointer :: tp

      tp => timer_root
      do
         if ( associated(tp%next)) then
            if ( item%key == tp%next%info%key ) then
               item%time = modify_timer(tp%next, item%reset)
            else if ( item%key < tp%next%info%key) then
               call insert_timer(tp%next, item)
            else
               tp => tp%next%node
               cycle ! keep looking
            endif
         else
            call insert_timer(tp%next,item)
         endif
         return
      enddo

   contains

      real function modify_timer(tp,reset)

         use mpi, only: MPI_Wtime

         implicit none

         type(timer_node), pointer :: tp
         logical, intent(in) :: reset
         real :: time_old

         time_old = tp%info%time
         tp%info%time = MPI_Wtime()
         !if (.not.reset) write(*,'(3A,F7.3,A)') "Timer [",trim(tp%info%key),"] = ", tp%info%time - time_old, " s"  ! QA_WARN debug?
         if (.not.reset) then
            modify_timer = tp%info%time - time_old
         else
            modify_timer = tp%info%time
         endif

      end function modify_timer

      subroutine insert_timer(tp, item)

         use mpi, only: MPI_Wtime

         implicit none

         type(timer_node), pointer :: tp, temp
         type(timer_info), intent(in) :: item

         allocate(temp)
         temp%info = item
         temp%info%time = MPI_Wtime()
         temp%node%next => tp
         tp => temp

      end subroutine insert_timer

   end subroutine search_timer

!>
!! \brief  Initialize cpu and wall clocks.
!! \details "cputot" will be the total cpu time (in seconds) consumed by this job.
!! "wctot" will be the total elapsed wall clock time (in seconds) since the job began.
!<
   subroutine timer_start

      implicit none

      real(kind=4) :: dtime

#ifdef PERFMON
      call itime ( iarray )
#endif /* PERFMON */
      wctot  = iarray(1) * 3600. + iarray(2) * 60. + iarray(3)
      cputot  = dtime ( tarray )
      cpu_start = tarray(1) +tarray(2)

   end subroutine timer_start

!------------------------------------------------------------------------------------------

   function time_left(wend) result (tf)

      use dataio_pub, only: msg, printinfo, warn
      use mpisetup,   only: slave

      implicit none

      real(kind=8), intent(in), optional :: wend
      logical :: tf
      integer(kind=8) :: clock, cnt_rate, cnt_max
      real(kind=8)    :: r_clk_end

      integer(kind=8) :: hh, mm
      real(kind=8)    :: ss

      if (slave) then
         call warn("[timer:time_left] slaves are not supposed to call me. ABORT!")
      else

         if (present(wend)) then
            if (wend >= 0.0d0) then
               call system_clock(clock_start, cnt_rate, cnt_max)
               if (wend < 1d-4*huge(1.0d0)) then
                  r_clk_end = real(clock_start, kind=8) + wend*3600.0*real(cnt_rate, kind=8)
                  if (r_clk_end < cnt_max) then
                     clock_end = int(r_clk_end, kind=8)
                  else
                     clock_end = -cnt_max
                  endif
               else
                  clock_end = -cnt_max
               endif
            endif
         endif
!>
!! \deprecated BEWARE: gfortran gives 1ms resolution, but ifort can offer 0.1ms, which will result in an integer overflow in less than 5 days
!! Probably it is better to call date_and_time(VALUES) here
!! if clock, cnt_rate, cnt_max are of kind=8 gfortran put results in 1 ns
!<
         call system_clock(clock, cnt_rate, cnt_max)
         tf = .true.
         if (clock_end /= -cnt_max) then
            if ( clock_end - clock < 0 ) tf = .false.
         endif

         if (present(wend)) then
            if (wend < 0.0d0) then
               ss = (clock_end - clock)/(3600.0*cnt_rate)
               hh = int(ss,kind=8)
               ss = (ss - hh)*60.0
               mm = int(ss,kind=8)
               ss = (ss - mm)*60.0
               write(msg,'(A,I4,":",I2.2,":",F5.2)') "Walltime left until end of simulation: ", hh, mm, ss
               call printinfo(msg)
            endif
         endif
      endif

   end function time_left

!>
!! \brief Final wall clock time and cpu usage.
!! \details Final wall clock time, expressed in hours, minutes, and seconds.
!!          cpu usage, expressed in hours, minutes, and seconds.
!<
   subroutine timer_stop

      use constants,  only: I_ONE, half
      use dataio_pub, only: msg, printinfo
      use domain,     only: dom
      use global,     only: nstep
      use mpi,        only: MPI_DOUBLE_PRECISION, MPI_SUM
      use mpisetup,   only: comm, ierr, master, FIRST

      implicit none

      real(kind=4)       :: dtime

!      Final wall clock time, expressed in hours, minutes, and seconds.
!
#ifdef PERFMON
      call itime ( iarray )
#endif /* PERFMON */
      wctot  = iarray(1) * 3600. + iarray(2) * 60. + iarray(3) - wctot
      wchours  =  int ( wctot / 3600.0 )
      wcmins   =  int ( wctot / 60.0   ) - 60   * wchours
      wcsecs   =  int ( wctot + half   ) - 3600 * wchours &
                                         - 60   * wcmins
!
!      cpu usage, expressed in hours, minutes, and seconds.
!
      cputot  = dtime ( tarray )
      cpu_stop = tarray(1) +tarray(2)
      cputot = cpu_stop-cpu_start

      cpuhours =  int ( cputot / 3600.0 )
      cpumins  =  int ( cputot / 60.0   ) - 60   * cpuhours
      cpusecs  =  int ( cputot + half   ) - 3600 * cpuhours &
                                        - 60   * cpumins

      call MPI_Reduce(cputot, cpuallp, I_ONE, MPI_DOUBLE_PRECISION, MPI_SUM, FIRST, comm, ierr)

      if (master) then

         zcps  = real(nstep) * real(dom%total_ncells) / cpuallp

         call printinfo("", .true.)
         write(msg, "('CPU time        = ', f12.2,' s')") cpuallp
         call printinfo(msg)
         write(msg, "('Wall clock time = ', f12.2,' s')") wctot
         call printinfo(msg)
         write(msg, "('Zone-cycles / s = ',en14.5)") zcps
         call printinfo(msg)
         call printinfo("", .true.)

      endif

   end subroutine timer_stop

end module timer
