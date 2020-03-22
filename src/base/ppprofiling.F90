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
!! \brief Module for obtaining precise parallel profiling data
!!
!! \details MPI_Wtime-based event log that will be useful for parallel profiling
!! of the Piernik code.
!!
!! General purpose profilers don't know what is important in Piernika and tend to
!! provide a lot of detailed information about irrelevant routines.
!! This sometimes adds too much overhear from instrumenting the code which can be
!! misleading. Here, the developer can choose which parts of the code need closer
!! look.
!!
!! Collecting the events should be as cheap as possible, thus we don't construct
!! any trees inside Piernik. All the collected data is meant to be used in
!! postprocessing. The developer can set up several logs, each focused on
!! different aspect of the code and not interfering with other logs.
!<

module ppp

   use constants, only: cbuff_len

   implicit none

   private
   public :: eventlist, init_profiling

   ! namelist parametrs
   logical               :: use_profiling  !< control whether to do any PPProfiling or not
   logical               :: profile_hdf5   !< Use HDF5 output when possible, fallback to ASCII

   integer, parameter :: ev_arr_num = 10    ! number of allowed event arrays

   !> \brief a cheap single-event entry
   type :: event
      character(len=cbuff_len) :: label  !< label used to identify the event
      real(kind=8)             :: wtime  !< output from MPI_Wtime(); positive for start, negative for end
   end type event

   !> \brief a cheap array of events
   type :: eventarray
      type(event), dimension(:), allocatable :: ev_arr
   contains
      procedure :: arr_init     !< allocate eventarray of given size
      procedure :: arr_cleanup  !< deallocate eventarray
   end type eventarray

   !> \briev list of events based on arrays of events, cheap to expand, avoid reallocation
   type eventlist
      private
      character(len=cbuff_len) :: label  !< label used to identify the event list
      type(eventarray), dimension(ev_arr_num) :: arrays  ! separate arrays to avoid lhs-reallocation
      integer :: arr_ind  ! currently used array
      integer :: ind      ! first unused entry in currently used array
   contains
      procedure :: init     !< create new event list
      procedure :: cleanup  !< destroy this event list
      procedure :: start    !< add a beginning of an interval
      procedure :: stop     !< add an end of an interval
      procedure :: put      !< add an event with time measured elsewhere
      procedure, private :: next_event  !< for internal use in start, stop and put
      procedure, private :: expand  !< create next array for events
      procedure :: publish  !< write the collected data to a log file
   end type eventlist

contains

!>
!! \brief initialize profiling output according to parameters from namelist PROFILING
!!
!! \n \n
!! @b PROFILING
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>use_profiling  </td><td>.false.  </td><td>logical value  </td><td>\copydoc ppp::use_profiling   </td></tr>
!! <tr><td>profile_hdf5   </td><td>.true.   </td><td>logical value  </td><td>\copydoc ppp::profile_hdf5    </td></tr>
!! </table>
!! \n \n
!<

   subroutine init_profiling

      use dataio_pub, only: nh, warn
      use mpisetup,   only: lbuff, master, slave, piernik_MPI_Bcast

      implicit none

      namelist /PROFILING/ use_profiling, profile_hdf5

      use_profiling = .false.
      profile_hdf5 = &
#ifdef HDF5
           .true.
#else /* !HDF5 */
           .false.
#endif /* HDF5 */

     if (master) then

         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=PROFILING)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=""
         read(unit=nh%lun, nml=PROFILING, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "PROFILING")
         read(nh%cmdl_nml,nml=PROFILING, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "PROFILING", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=PROFILING)
         close(nh%lun)
         call nh%compare_namelist()

         lbuff(1) = use_profiling
         lbuff(2) = profile_hdf5

      endif

      call piernik_MPI_Bcast(lbuff)

      if (slave) then

         use_profiling = lbuff(1)
         profile_hdf5  = lbuff(2)

      endif

      if (profile_hdf5) then
         profile_hdf5 = .false.
         if (master) call warn("[ppprofiling:init_profiling] profile_hdf5 not implemented yet")
      endif

   end subroutine init_profiling

!> \brief allocate eventarray of given size

   subroutine arr_init(this, asize)

      use constants,  only: PIERNIK_INIT_MPI
      use dataio_pub, only: die, code_progress
      use global,     only: check_mem_usage

      implicit none

      class(eventarray), intent(inout) :: this   !< an object invoking the type-bound procedure
      integer,           intent(in)    :: asize  !< size of the event array

      if (allocated(this%ev_arr)) call die("[ppprofiling:arr_init] already allocated")
      allocate(this%ev_arr(asize))
!      this%ev_arr(:)%wtime = 0.
      if (code_progress >= PIERNIK_INIT_MPI) call check_mem_usage

   end subroutine arr_init

!> \brief deallocate eventarray

   subroutine arr_cleanup(this)

      use dataio_pub, only: die

      implicit none

      class(eventarray), intent(inout) :: this   !< an object invoking the type-bound procedure

      if (.not. allocated(this%ev_arr)) call die("[ppprofiling:arr_init] not allocated")
      deallocate(this%ev_arr)

   end subroutine arr_cleanup

!> \brief Create new event list

   subroutine init(this, label)

      implicit none

      class(eventlist), intent(inout) :: this   !< an object invoking the type-bound procedure
      character(len=*), intent(in)    :: label  !< event list label

      integer, parameter :: ev_arr_len = 1024  ! starting size of the array of events

      if (.not. use_profiling) return

      this%ind = 1
      this%arr_ind = 1
      this%label = label(1:min(cbuff_len, len_trim(label)))
      call this%arrays(this%arr_ind)%arr_init(ev_arr_len)

   end subroutine init

!> \brief Destroy this event list

   subroutine cleanup(this)

      use constants, only: INVALID

      implicit none

      class(eventlist), intent(inout) :: this   !< an object invoking the type-bound procedure

      integer :: i

      do i = lbound(this%arrays, dim=1), ubound(this%arrays, dim=1)
         if (allocated(this%arrays(i)%ev_arr)) call this%arrays(i)%arr_cleanup
      enddo

      this%ind = INVALID
      this%arr_ind = INVALID

   end subroutine cleanup

!>
!! \brief Add a beginning of an interval
!!
!! Do not use inside die(), warn(), system_mem_usage() or check_mem_usage()
!<

   subroutine start(this, label)

      use mpi, only: MPI_Wtime

      implicit none

      class(eventlist), intent(inout) :: this   !< an object invoking the type-bound procedure
      character(len=*), intent(in)    :: label  !< event label

      character(len=cbuff_len) :: l

      if (.not. use_profiling) return

      l = label(1:min(cbuff_len, len_trim(label)))
      call this%next_event(event(l, MPI_Wtime()))

   end subroutine start

!>
!! \brief Add an end of an interval, use negative sign
!!
!! Do not use inside die(), warn(), system_mem_usage() or check_mem_usage()
!<

   subroutine stop(this, label)

      use mpi, only: MPI_Wtime

      implicit none

      class(eventlist), intent(inout) :: this   !< an object invoking the type-bound procedure
      character(len=*), intent(in)    :: label  !< event label

      character(len=cbuff_len) :: l

      if (.not. use_profiling) return

      l = label(1:min(cbuff_len, len_trim(label)))
      call this%next_event(event(l, -MPI_Wtime()))

   end subroutine stop

!>
!! \brief Add an event with time measured elsewhere
!!
!! Do not use inside die(), warn(), system_mem_usage() or check_mem_usage().
!!
!! Do not use at all, except for inserting event with mpisetup::bigbang
!<

   subroutine put(this, label, time)

      implicit none

      class(eventlist), intent(inout) :: this   !< an object invoking the type-bound procedure
      character(len=*), intent(in)    :: label  !< event label
      real(kind=8),     intent(in)    :: time   !< time of the event

      character(len=cbuff_len) :: l

      if (.not. use_profiling) return

      l = label(1:min(cbuff_len, len_trim(label)))
      call this%next_event(event(l, time))

   end subroutine put

!> \brief Start, stop and put

   subroutine next_event(this, ev)

      implicit none

      class(eventlist), intent(inout) :: this   !< an object invoking the type-bound procedure
      type(event),      intent(in)    :: ev     !< an event to be stored

      this%arrays(this%arr_ind)%ev_arr(this%ind) = ev
      this%ind = this%ind + 1
      if (this%ind > ubound(this%arrays(this%arr_ind)%ev_arr, dim=1)) then
         call this%expand
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
!!
!! \todo implement maximum of the log size and then politely give up instead of die()
!! when the size is exceeded.
!<

   subroutine expand(this)

      use dataio_pub, only: die

      implicit none

      class(eventlist), intent(inout) :: this   !< an object invoking the type-bound procedure

      integer, parameter :: grow_factor = 2  !< each next array should be bigger

      if (this%arr_ind >= ev_arr_num) call die("[ppprofiling:expand] Cannot add more arrays (ev_arr_num exceeded)")

      call this%arrays(this%arr_ind + 1)%arr_init(grow_factor*size(this%arrays(this%arr_ind)%ev_arr))

      this%arr_ind = this%arr_ind + 1
      this%ind = 1
      this%arrays(this%arr_ind)%ev_arr(this%ind)%wtime = 0.

   end subroutine expand

!>
!! \brief Write the collected data to a log file and clear the log
!!
!! \todo Use HDF5 format when available, stdout otherwise
!<

   subroutine publish(this)

      use dataio_pub, only: msg, printinfo
      use func,       only: operator(.equals.)
      use mpisetup,   only: proc

      implicit none

      class(eventlist), intent(inout) :: this   !< an object invoking the type-bound procedure

      integer :: i, ia

      if (.not. use_profiling) return

      do ia = lbound(this%arrays, dim=1), ubound(this%arrays, dim=1)
         if (allocated(this%arrays(ia)%ev_arr)) then
            do i = lbound(this%arrays(ia)%ev_arr, dim=1), ubound(this%arrays(ia)%ev_arr, dim=1)
               associate (ev => this%arrays(ia)%ev_arr(i))
                  if (ev%wtime .equals. 0.) exit
                  write(msg, '(2a,i4,2a,f20.7)') this%label, " @", proc, " ", ev%label, ev%wtime
                  call printinfo(msg)
               end associate
            enddo
         endif
      enddo

      call this%cleanup

   end subroutine publish

end module ppp
