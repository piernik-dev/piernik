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
!    Initial implemetation of PIERNIK code was based on TVD split MHD code by
!    Ue-Li Pen
!        see: Pen, Arras & Wong (2003) for algorithm and
!             http://www.cita.utoronto.ca/~pen/MHD
!             for original source code "mhd.f90"
!
!    For full list of developers see $PIERNIK_HOME/license/pdt.txt
!
#include "piernik.def"
!>
!! \brief (KK)
!<
module timer

   implicit none
   integer, parameter, private :: S_LEN = 30

   private :: search_timer, delete_timer, clock_start, clock_end

   type, private :: timer_info
      character(len=S_LEN) :: key
      logical :: reset
      real :: time
   end type timer_info

   type, private :: timer_list
      type(timer_node), pointer :: next => Null()
   end type timer_list

   type, private :: timer_node
      type(timer_info) :: info
      type(timer_list) :: node
   end type timer_node

   type(timer_list), target, private, save :: timer_root

   integer :: nzones, cpuhours, cpumins, cpusecs , wchours , wcmins  , wcsecs
   real    :: zcps,  cputot, cpuallp, wctot, cpu_start, cpu_stop
   integer :: iarray(3)
   real(kind=4), dimension(2) :: tarray
   integer :: clock_start, clock_end

contains

   !>
   !! \brief Wrapper for timer handling.
   !! Only this routine should be used outside timer module.
   !! USAGE: call timer(name,reset)
   !!    1) if first called with "name", create timer "name" and set it with current cpu_time(), no output
   !!    2) if timer "name" exists, print "name" - current time, set "name" with current cpu_time()
   !!    (optional) if reset is true suppress output, set "name" with current cpu_time
   !<
   real function timer_(str,reset)
      implicit none
      character(len=*), intent(in) :: str    !< name of the timer
      logical, intent(in), optional :: reset !< if true all output is suppressed, use for resetting timers
      type(timer_info) :: temp
      temp%key = trim(str)
      if(present(reset)) then
         temp%reset = .true.
      else
         temp%reset = .false.
      endif
      call search_timer(temp)
      timer_ = temp%time

   end function timer_

   function delete_timer(tp) result (item)
      implicit none
      type(timer_node), pointer :: tp
      type(timer_info) :: item
      type(timer_node), pointer :: temp

      temp => tp
      item = tp%info
      tp => tp%node%next
      deallocate(temp)
      return
   end function delete_timer

   subroutine cleanup_timers

      implicit none

      type(timer_list), pointer :: tp
      type(timer_info) :: item

      tp => timer_root
#ifdef VERBOSE
      write(*,*)
#endif /* VERBOSE */
      do while (associated(tp%next))
         item = delete_timer(tp%next)
#ifdef VERBOSE
         write(*,*) "[timer:cleanup_timers]: Timer ",item%key," deleted"
#endif /* VERBOSE */
      enddo
      return

   end subroutine cleanup_timers

   subroutine search_timer(item)
      implicit none
      type(timer_info), intent(inout) :: item
      type(timer_list), pointer :: tp
      tp => timer_root

      do
         if( associated(tp%next)) then
            if( item%key == tp%next%info%key ) then
               item%time = modify_timer(tp%next, item%reset)
            else if( item%key < tp%next%info%key) then
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
         implicit none
         type(timer_node), pointer :: tp
         logical, intent(in) :: reset

         real :: time_old

         time_old = tp%info%time
         call cpu_time(tp%info%time)
         !if(.not.reset) write(*,'(A,F7.3,A)') "Timer ["//trim(tp%info%key)//"] = ", tp%info%time - time_old, " s"
         if(.not.reset) then
            modify_timer = tp%info%time - time_old
         else
            modify_timer = tp%info%time
         endif
      end function modify_timer

      subroutine insert_timer(tp, item)
         implicit none
         type(timer_node), pointer :: tp
         type(timer_info), intent(in) :: item
         type(timer_node), pointer :: temp

         allocate(temp)
         temp%info = item
         call cpu_time(temp%info%time)
         temp%node%next => tp
         tp => temp

         return
      end subroutine insert_timer

   end subroutine search_timer

   subroutine timer_start
      implicit none
      real(kind=4) :: dtime
!
!  Initialise cpu and wall clocks.  "cputot" will be the total cpu
!  time (in seconds) consumed by this job.  "wctot" will be the total
!  elapsed wall clock time (in seconds) since the job began.
!

#ifdef PERFMON
      call itime ( iarray )
#endif /* PERFMON */
      wctot  = iarray(1) * 3600. + iarray(2) * 60. + iarray(3)

      cputot  = dtime ( tarray )
      cpu_start = tarray(1) +tarray(2)

   end subroutine timer_start

!------------------------------------------------------------------------------------------

   function time_left(wend) result (tf)
      implicit none
      real, intent(in), optional :: wend
      logical :: tf
      integer :: clock, cnt_rate, cnt_max
      real    :: r_clk_end

      if(present(wend)) then
         call system_clock(clock_start, cnt_rate, cnt_max)
 !         clock_end = clock_start + int(wend*3600.*cnt_rate)
         r_clk_end = clock_start + wend*3600.*cnt_rate
         if (r_clk_end < cnt_max) then
            clock_end = int(r_clk_end)
         else
            clock_end = -cnt_max
         endif
      endif
 ! BEWARE: gfortran gives 1ms resolution, but ifort can offer 0.1ms, which will result in an integer overflow in less than 5 days
 ! Probably it is better to call date_and_time(VALUES) here
      call system_clock(clock, cnt_rate, cnt_max)
      tf = .true.
      if (clock_end /= -cnt_max) then
         if( clock_end - clock < 0 ) tf = .false.
      end if

   end function time_left

   subroutine timer_stop
      use mpisetup,      only : MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr, nstep, proc
      use dataio_public, only : log_file, log_lun
      use grid,          only : nxd,nyd,nzd

      implicit none
      real(kind=4) :: dtime

!      Final wall clock time, expressed in hours, minutes, and seconds.
!
#ifdef PERFMON
      call itime ( iarray )
#endif /* PERFMON */
      wctot  = iarray(1) * 3600. + iarray(2) * 60. + iarray(3) - wctot
      wchours  =  int ( wctot / 3600.0 )
      wcmins   =  int ( wctot / 60.0   ) - 60   * wchours
      wcsecs   =  int ( wctot + 0.5    ) - 3600 * wchours &
                                         - 60   * wcmins
!
!      cpu usage, expressed in hours, minutes, and seconds.
!
      cputot  = dtime ( tarray )
      cpu_stop = tarray(1) +tarray(2)
      cputot = cpu_stop-cpu_start

      cpuhours =  int ( cputot / 3600.0 )
      cpumins  =  int ( cputot / 60.0   ) - 60   * cpuhours
      cpusecs  =  int ( cputot + 0.5    ) - 3600 * cpuhours &
                                        - 60   * cpumins
      nzones = nxd * nyd * nzd


      call MPI_REDUCE(cputot, cpuallp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)

      if(proc == 0) then

         zcps  = real(nstep) * real(nzones) / cpuallp

         write(*,*)
         write (*, 10) cpuallp
         write (*, 20) wctot
         write (*, 30) zcps
         write(*,*)

         open(log_lun, file=log_file, position='append')
            write(log_lun,*)
            write (log_lun, 10) cpuallp
            write (log_lun, 20) wctot
            write (log_lun, 30) zcps
            write(log_lun,*)
         close(log_lun)
      endif

10    format('CPU time        = ', f12.2,' s')
20    format('Wall clock time = ', f12.2,' s')
30    format('Zone-cycles / s = ',es12.5)

   end subroutine timer_stop

end module timer
