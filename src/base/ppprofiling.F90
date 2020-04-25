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
!! General purpose profilers don't know what is important in Piernik and tend to
!! provide a lot of detailed information about irrelevant routines.
!! This sometimes adds too much overhead from instrumenting the code which can be
!! misleading. Here, the developer can choose which parts of the code need closer
!! look.
!!
!! Collecting the events should be as cheap as possible, thus we don't construct
!! any trees inside Piernik. All the collected data is meant to be used in
!! postprocessing. The developer can set up several logs, each focused on
!! different aspect of the code and not interfering with other logs.
!!
!! The profiling is turned off by default. One can turn it on by:
!! * enabling in problem.par or in CLI
!!     ./piernik -n "&PROFILING use_profiling = T/"
!!   to get the whole run profiled.
!! * enabling via user_message_file
!!      echo ppp > msg
!!   to enable profiling for just one step, or
!!      echo ppp 10 > msg
!!   to enable profiling for 10 steps.
!! Interacting with PPP via msg file will override the choice for use_profiling
!! namelist parameter and can be uset to stop profiling and flush the data to the
!! output file and continue the simulation without interrupting. The profiling may
!! be later turned on when needed. For example:
!!    ./piernik -n "&PROFILING use_profiling = T/" &\
!!        sleep 15 && echo ppp > msg &&\
!!        sleep 10 && echo ppp 2 > msg &&\
!!        sleep 1 && echo stop > msg
!! may give profiling data containing "init_piernik", first few steps, last step and
!! "finalize" (exact list depends on how long the steps are executed in the actual run).
!<

module ppp

   use constants, only: cbuff_len, cwdlen

   implicit none

   private
   public :: eventlist, init_profiling, cleanup_profiling, update_profiling, ppp_main, umsg_request

   ! namelist parametrs
   logical               :: use_profiling  !< control whether to do any PPProfiling or not
   logical               :: profile_hdf5   !< Use HDF5 output when possible, fallback to ASCII

   character(len=cwdlen) :: profile_file       !< file name for the profile data
   integer :: profile_lun                      !< logical unit number for profile file
   integer, save :: umsg_request = 0           !< turn on profiling for next umsg_request steps (read from msg file)
   logical, save :: profile_file_cr = .false.  !< .true after we open the profile file for writing
   enum, bind(C)
      enumerator :: TAG_CNT = 1, TAG_ARR_L, TAG_ARR_T
   end enum

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

   type(eventlist) :: ppp_main  ! main eventlist

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

      use dataio_pub, only: nh, warn, log_wr, problem_name, run_id, nrestart
      use mpisetup,   only: lbuff, master, slave, piernik_MPI_Bcast

      implicit none

      namelist /PROFILING/ use_profiling, profile_hdf5

      use_profiling = .false.
      profile_hdf5 = &
#ifdef HDF5
           .false. ! ToDo: change to .true. after it gets implemented.
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
         if (master) call warn("[ppprofiling:init_profiling] profile_hdf5 not implemented yet. Forcing ASCII output.")
      endif

      write(profile_file, '(6a,i3.3,a)') trim(log_wr), '/', trim(problem_name), '_', trim(run_id), '_', nrestart, '.ppprofile' // trim(merge(".h5   ", ".ascii", profile_hdf5))

      call ppp_main%init("main")

   end subroutine init_profiling

!> \brief Turn on and off profiling upon request

   subroutine update_profiling

      use dataio_pub, only: printinfo
      use global,     only: nstep
      use mpisetup,   only: master

      implicit none

      integer, save :: turn_off = huge(1)

      if (turn_off <= nstep) then  ! turn off ppp_main profiling
         call ppp_main%publish
         if (master) call printinfo("[ppprofiling:update_profiling] Stop PPP")
         use_profiling = .false.
         turn_off = huge(1)
      endif

      if (umsg_request < 1) return

      ! this allows for overriding previous umsg_request
      ! or to disable profiling enabled in problem.par and flush the data
      turn_off = nstep + umsg_request
      umsg_request = 0

      if (.not. use_profiling) then  ! turn on ppp_main profiling
         use_profiling = .true.
         call ppp_main%init("main")
         if (master) call printinfo("[ppprofiling:update_profiling] Start PPP")
      endif

   end subroutine update_profiling

!> \brief close the profile file

   subroutine cleanup_profiling

      use dataio_pub, only: close_txt_file, die
      use mpisetup,   only: master

      implicit none

      if (master) then
         if (profile_hdf5) then
            call die("[ppprofiling:cleanup_profiling] HDF5 profiles not implemented yet")
         else
            call close_txt_file(profile_file, profile_lun)
         endif
      endif

   end subroutine cleanup_profiling

!> \brief allocate eventarray of given size

   subroutine arr_init(this, asize)

      use constants,    only: PIERNIK_INIT_MPI
      use dataio_pub,   only: die, code_progress
      use memory_usage, only: check_mem_usage

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
!! \todo Use HDF5 format when available, ASCII otherwise
!!
!! Perhaps MPI_TYPE_CREATE_STRUCT would simplify the code and improve the communication
!! but it requires some C-interoperability, which needs to be explored and tested first.
!!
!! \todo XMLot JSON output?
!<

   subroutine publish(this)

      use constants,    only: I_ZERO, I_ONE
      use dataio_pub,   only: die, printinfo
      use memory_usage, only: check_mem_usage
      use mpi,          only: MPI_STATUS_IGNORE, MPI_CHARACTER, MPI_INTEGER, MPI_DOUBLE_PRECISION
      use mpisetup,     only: proc, master, slave, comm, mpi_err, FIRST, LAST, req, status, inflate_req

      implicit none

      class(eventlist), intent(inout) :: this   !< an object invoking the type-bound procedure

      integer(kind=4) :: t, ne, p
      integer :: ia
      character(len=cbuff_len), dimension(:), allocatable  :: buflabel
      real(kind=8), dimension(:), allocatable :: buftime
      integer(kind=4), dimension(:,:), pointer :: mpistatus

      if (.not. use_profiling) return

      if (master .and. .not. profile_file_cr) then
         if (profile_hdf5) call die("[ppprofiling:publish] HDF5 profiles not implemented yet")
         call printinfo("[ppprofiling:publish] Profile timings will be written to '" // trim(profile_file) // "' file.")
         open(newunit=profile_lun, file=profile_file)
         profile_file_cr = .true.
      endif

      ! send
      call inflate_req(int(TAG_ARR_T))
      t = TAG_CNT
      ne = I_ZERO
      do ia = lbound(this%arrays, dim=1), ubound(this%arrays, dim=1)
         if (allocated(this%arrays(ia)%ev_arr)) then
            ne = ne + size(this%arrays(ia)%ev_arr, kind=4)
         else
            exit
         endif
      enddo
      if (slave) call MPI_Isend(ne, I_ONE, MPI_INTEGER, FIRST, TAG_CNT, comm, req(TAG_CNT), mpi_err)

      if (ne > 0) then
         allocate(buflabel(ne), buftime(ne))
         call check_mem_usage
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
            call MPI_Isend(buflabel, size(buflabel)*len(buflabel(1)), MPI_CHARACTER,        FIRST, TAG_ARR_L, comm, req(TAG_ARR_L), mpi_err)
            call MPI_Isend(buftime,  size(buftime),                   MPI_DOUBLE_PRECISION, FIRST, TAG_ARR_T, comm, req(TAG_ARR_T), mpi_err)
            t = TAG_ARR_T
         endif
      endif

      if (master) then
         ! write(profile_lun, '(/,3a)') "#profile '", trim(this%label), "'"

         ! receive
         do p = FIRST + I_ONE, LAST
             call MPI_Recv(ne, I_ONE, MPI_INTEGER, p, TAG_CNT, comm, MPI_STATUS_IGNORE, mpi_err)
             if (ne /= 0) then
                allocate(buflabel(ne), buftime(ne))
                call check_mem_usage
                call MPI_Recv(buflabel, size(buflabel)*len(buflabel(1)), MPI_CHARACTER,        p, TAG_ARR_L, comm, MPI_STATUS_IGNORE, mpi_err)
                call MPI_Recv(buftime,  size(buftime),                   MPI_DOUBLE_PRECISION, p, TAG_ARR_T, comm, MPI_STATUS_IGNORE, mpi_err)
                call publish_buffers(p, buflabel, buftime)
                deallocate(buflabel, buftime)
             endif
         enddo
      else
         mpistatus => status(:, :t)
         call MPI_Waitall(t, req(:t), mpistatus, mpi_err)
         deallocate(buflabel, buftime)
      endif

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
         call warn("[ppprofiling:publish_array] only master is supposed to write")
         return
      endif

      if (profile_hdf5) call die("[ppprofiling:publish_array] HDF5 not implemented yet")

      if (size(ev_label) /= size(ev_time)) call die("[ppprofiling:publish_array] arrays size mismatch")

      d = 1
      do i = lbound(ev_label, dim=1), ubound(ev_label, dim=1)
         if (ev_time(i) .equals. 0.) exit
         if (ev_time(i) < 0.) d = d - 1
         write(profile_lun, '(i5,a,f20.7,3a)') process, " ", ev_time(i), " ", repeat("  ", d), trim(ev_label(i))
         if (ev_time(i) > 0.) d = d + 1
      enddo
      flush(profile_lun)

   end subroutine publish_buffers

end module ppp
