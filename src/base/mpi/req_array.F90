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

!> \brief Type for array of requests

module req_array

   use cnt_array,  only: arrsum
   use log_bins,   only: lbins
   use tag_arrays, only: tag_arrs
#ifdef MPIF08
   use MPIF,       only: MPI_Request, MPI_Comm
#endif /* MPIF08 */

   implicit none

   private
   public :: req_arr, req_wall, init_wall, cleanup_wall, C_REQA

   type :: req_arr  ! array of requests
      !< request array for MPI_Waitall
#ifdef MPIF08
      type(MPI_Request), allocatable, dimension(:) :: r     !< array of requests
      type(MPI_Comm)                               :: comm  !< own communicator
#else /* !MPIF08 */
      integer(kind=4),   allocatable, dimension(:) :: r
      integer(kind=4)                              :: comm
#endif /* !MPIF08 */
      integer(kind=4) :: n   !< number of requests
      logical :: owncomm     !< was the own communicator requested?
      type(tag_arrs) :: tags !< list of tags for each proc in current communication
   contains
      procedure :: nxt              !< give next unused request
      procedure :: cleanup          !< free the resources
      procedure :: store_tag        !< store tags for inspection upon Waitall
      procedure :: waitall_on_some  !< wait for completion of the requests on selected procs
      procedure :: waitall_wrapper  !< code common for waitall_on_some and req_ppp::waitall_ppp
      procedure :: setsize_req      !< set the storage for this%r
      procedure :: doublesize_req   !< double the storage size for this%r
      generic, public :: init => doublesize_req, setsize_req
      generic, public :: waitall => waitall_on_some
   end type req_arr

   ! logarithmic bins for collecting Waitall execution time:
   real, parameter :: min_bin = 1e-8  ! first bin will contain everything below this value
   real, parameter :: binwidth = 10.  ! each bin from 2 to (nbins - 1) has this width
   type(lbins) :: tbins_wall  ! bins for counting execution time of MPI_Waitall

   type(arrsum) :: req_wall  ! counter for ppp_mpi
   enum, bind(C)
      enumerator :: C_REQA = 1, C_REQS  ! for both req% waitall methods
   end enum

   real :: longest_wall  ! find the extreme Waital completion time

   integer(kind=4) :: err_mpi  !< error status

contains

!> \brief Initialize MPI_Waitall stat counters

   subroutine init_wall

      use global, only: waitall_timeout

      implicit none

      real, parameter :: def_max_bin = 100.  ! in seconds, this is a lot of time to complete waitalls

      integer :: nbins
      real :: mx_b

      longest_wall = 0.

      call req_wall%init(int([C_REQS]))

      ! assume sane values of min_bin, def_max_bin, binwidth and waitall_timeout
      mx_b = merge(waitall_timeout, def_max_bin, waitall_timeout > min_bin)
      nbins = int((log(mx_b / min_bin) + 0.1)/ log(binwidth)) + 3  ! with some safety factor
      call tbins_wall%init(min_bin, binwidth, nbins)

   end subroutine init_wall

!> \brief Print log and clean up MPI_Waitall stat counters

   subroutine cleanup_wall

      use constants,    only: V_DEBUG, I_ONE
      use dataio_pub,   only: printinfo, msg
      use mpi_wrappers, only: MPI_wrapper_stats
      use mpisetup,     only: master
      use MPIF,         only: MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD
      use MPIFUN,       only: MPI_Allreduce

      implicit none

      if (MPI_wrapper_stats) then
         call MPI_Allreduce(MPI_IN_PLACE, longest_wall, I_ONE, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, err_mpi)

         if (master) then
            call req_wall%print("Waitall requests(calls). Columns: req_all%, req_some%.", V_DEBUG)
            call tbins_wall%print("Waitall execution times (bin boundaries in seconds)", V_DEBUG)
            write(msg, '(a,g0.2,a)')"Longest completion of MPI_Waitall took ", longest_wall, " s (max over all processes)"
            call printinfo(msg, V_DEBUG)
         endif
      endif

      call req_wall%cleanup
      call tbins_wall%cleanup

   end subroutine cleanup_wall

!> \brief provide next free request

   function nxt(this)

      use constants, only: I_ONE

      implicit none

      class(req_arr), intent(inout), target :: this  !< an object invoking the type-bound procedure
#ifdef MPIF08
      type(MPI_Request), pointer :: nxt
#else /* !MPIF08 */
      integer(kind=4),   pointer :: nxt
#endif /* !MPIF08 */

      this%n = this%n + I_ONE
      if (this%n > ubound(this%r, 1)) call this%doublesize_req(owncomm = .false.)  ! the communicator was already duplicated, if required
      nxt => this%r(this%n)

   end function nxt

!> \brief clean up

   subroutine cleanup(this)

      use MPIFUN, only: MPI_Comm_free

      implicit none

      class(req_arr), intent(inout), target :: this  !< an object invoking the type-bound procedure

      if (allocated(this%r)) deallocate(this%r)
      if (this%owncomm) then
         call MPI_Comm_free(this%comm, err_mpi)
         this%owncomm = .false.
      endif
      call this%tags%cleanup

   end subroutine cleanup

!> brief Accumulate tags for inspection

   subroutine store_tag(this, tag, other_proc, recv)

      implicit none

      class(req_arr),  intent(inout) :: this        !< an object invoking the type-bound procedure
      integer(kind=4), intent(in)    :: tag         !< the tag
      integer(kind=4), intent(in)    :: other_proc  !< source or destination
      logical,         intent(in)    :: recv        !< send or receive?

      call this%tags%store_tag(tag, other_proc, this%n, recv)

   end subroutine store_tag

!>
!! \brief Simple wrapper for MPI_Waitall
!!
!! For PPP-instrumented variant see req_ppp::waitall_ppp (in ppp_mpi module)
!<

   subroutine waitall_on_some(this)

      implicit none

      class(req_arr), intent(inout) :: this  !< an object invoking the type-bound procedure

      call this%waitall_wrapper(C_REQS, "_" // trim(this%tags%label) // "_")
      call this%cleanup

   end subroutine waitall_on_some

!>
!! \brief Code common for waitall_on_some and req_ppp::waitall_ppp
!!
!! \details Use positive waitall_timeout if you experienced deadlocks in your runs.
!! With waitall_timeout > 0. this routine won't call MPI_Waitall, which is blocking
!! but will use nonblocking MPI_Test* calls instead. Chances are that some useful hints
!! will be printed on the stdout instead of going into silent deadlock.
!<

   subroutine waitall_wrapper(this, ind, label)

      use constants,  only: LONG, I_ONE, V_ESSENTIAL
      use dataio_pub, only: warn, msg, die, printinfo
      use global,     only: waitall_timeout
      use MPIF,       only: MPI_STATUSES_IGNORE, MPI_STATUS_IGNORE, MPI_INTEGER_KIND, MPI_UNDEFINED, MPI_Wtime
      use MPIFUN,     only: MPI_Waitall, MPI_Testsome, MPI_Test

      implicit none

      class(req_arr),   intent(inout) :: this  !< an object invoking the type-bound procedure
      integer(kind=4),  intent(in)    :: ind   !< the index for call counter
      character(len=*), intent(in)    :: label  !< identifier

      logical, save :: firstcall = .true.
      real :: t0, dtw
      integer(kind=MPI_INTEGER_KIND), dimension(:), allocatable :: inds
      integer(kind=MPI_INTEGER_KIND) :: tcnt, cnt_prev
      integer(kind=4) :: i
      logical(kind=4) :: flag, warn10, fs, fr
      logical, allocatable, dimension(:) :: tflag

      ! The initialization had to be delayed because waitall_timeout is initialized in module global, which is a bit late.
      if (firstcall) then
         call init_wall
         firstcall = .false.
      endif

      if (this%n > 0) then

         warn10 = .false.
         call req_wall%add(int([ind]), int(this%n, kind=LONG))
         t0 = MPI_Wtime()
         if (waitall_timeout > 0.) then  ! complete requests in a non-blocking way, print somediagnostics if completion is slow
            allocate(inds(this%n))
            inds(:) = -1

            ! First quickly check if there are already completed requests
            tcnt = 0
            do i = 1, this%n
               call MPI_Test(this%r(i), flag, MPI_STATUS_IGNORE, err_mpi)
               if (flag) then
                  tcnt = tcnt + I_ONE
                  inds(tcnt) = i
               endif
            enddo

            ! Look for remaining requests
            do while (tcnt < this%n)

               if (tcnt /= MPI_UNDEFINED) cnt_prev = tcnt
               tcnt = MPI_UNDEFINED

               ! If we're exceeding 0.1 * waitall_timeout then is starts to be fishy: either too short waitall_timeout was set or we're going to to get a deadlock
               dtw = MPI_Wtime() - t0
               if (dtw > 0.1 * waitall_timeout .and. .not. warn10) then
                  write(msg, '(2(a,i0),a)')"[req_array:waitall_wrapper] MPI_Testsome progress at 1/10 timeout: ", cnt_prev, " out of ", this%n, " '" // trim(label) // "'"
                  call printinfo(msg, V_ESSENTIAL)
                  warn10 = .true.
               endif

               ! Test remaining requests
               call MPI_Testsome(this%n, this%r(:this%n), tcnt, inds(cnt_prev+1:), MPI_STATUSES_IGNORE, err_mpi)
               if ((tcnt < 0 .and. tcnt /= MPI_UNDEFINED) .or. cnt_prev < 0. .or. err_mpi /= 0) then
                  write(msg, '(a,4(" ",i0),a)')"[req_array:waitall_wrapper] MPI_Testsome failing?", tcnt, cnt_prev, this%n, err_mpi, " '" // trim(label) // "'"
                  call die(msg)
               endif

               if (tcnt /= MPI_UNDEFINED) tcnt = cnt_prev + tcnt

               ! If waitall_timeout was exceeded then die with some meaningful messages
               if ((tcnt == MPI_UNDEFINED .or. tcnt < this%n) .and. dtw > waitall_timeout) then
                  call set_tflag
                  call sleep(1)
                  fs = one_more_test(this%tags%ts)
                  fr = one_more_test(this%tags%tr) ! make sure that both are executed
                  if (fs .or. fr) then
                     ! waitall_timeout was set too short
                     write(msg, '(a,f0.3,a)')"[req_array:waitall_wrapper] Late completions detected. Increase waitall to at least ", max(1. + waitall_timeout, 2. * waitall_timeout), " s"
                     call warn(msg)
                  endif
                  deallocate(tflag)

                  if ((tcnt == MPI_UNDEFINED .or. tcnt < this%n) .and. dtw > waitall_timeout) then
                     call set_tflag
                     call print_details(this%tags%ts, "ISend ")
                     call print_details(this%tags%tr, "IRecv ")
                     deallocate(tflag)
                     write(msg, '(2(a,i0),a,g0.2,a)')"[req_array:waitall_wrapper] only ", merge(tcnt, cnt_prev, tcnt /= MPI_UNDEFINED), " out of ", this%n, &
                          " requests completed in ", dtw, "s (timeout, '" // trim(label) // "')"
                     call sleep(3)  ! give the other processes a chance to print their state
                     call die(msg)
                  endif
               endif

            enddo
            deallocate(inds)
         else
            ! This is perhaps faster but gives no clue when deadlock occur
            call MPI_Waitall(this%n, this%r(:this%n), MPI_STATUSES_IGNORE, err_mpi)
         endif
         dtw = MPI_Wtime() - t0
         call tbins_wall%put(dtw)
         if (dtw > longest_wall) longest_wall = dtw

         this%n = 0

      endif

   contains

      subroutine set_tflag

         implicit none

         allocate(tflag(this%n))
         tflag = .false.
         do i = 1, this%n
            if (inds(i) > -1) then
               if (inds(i) >= lbound(tflag, 1) .and. inds(i) <= ubound(tflag, 1)) &
                    tflag(inds(i)) = .true.
            endif
         enddo

      end subroutine set_tflag

      subroutine print_details(ta, l)

         use tag_array, only: tag_arr

         implicit none

         type(tag_arr),    intent(in) :: ta
         character(len=*), intent(in) :: l

         integer :: p, i

         if (allocated(ta%t)) then
            write(msg, '(2a,2(i0,a))') trim(l) // " '" // trim(this%tags%label), "', procs: [ ", lbound(ta%t), " : ", ubound(ta%t), " ] "
            call printinfo(msg, V_ESSENTIAL)
            do p = lbound(ta%t, 1), ubound(ta%t, 1)
               if (allocated(ta%t(p)%list)) then
                  write(msg, '(a,3(a,i0))') trim(l) // " '" // trim(this%tags%label), "', checking requests [ ", ta%t(p)%l_bound(), " : ", ta%t(p)%u_bound(), " ] for proc# ", p
                  call printinfo(msg, V_ESSENTIAL)
                  do i = ta%t(p)%l_bound(), ta%t(p)%u_bound()

                     ! Print the unfinished requests
                     if (.not. tflag(ta%t(p)%list(i)%ir)) then
                        write(msg, '(2a,4(i0,a))') trim(l) // " '" // trim(this%tags%label), "', unfinished with proc ", p, " request ", i, " : ( tag: ", ta%t(p)%list(i)%tag, " , req# ", ta%t(p)%list(i)%ir, " )"
                        call printinfo(msg, V_ESSENTIAL)
                     endif

                  enddo
               endif
            enddo
         endif

      end subroutine print_details

      logical function one_more_test(ta)

         use tag_array, only: tag_arr

         implicit none

         type(tag_arr),    intent(in) :: ta

         integer :: p, i

         one_more_test = .false.

         if (allocated(ta%t)) then
            do p = lbound(ta%t, 1), ubound(ta%t, 1)
               if (allocated(ta%t(p)%list)) then
                  do i = ta%t(p)%l_bound(), ta%t(p)%u_bound()
                     ! Do one more test, just in case
                     call MPI_Test(this%r(ta%t(p)%list(i)%ir), flag, MPI_STATUS_IGNORE, err_mpi)
                     if (flag .neqv. tflag(ta%t(p)%list(i)%ir)) then
                        one_more_test = .true.
                        tcnt = tcnt + I_ONE
                        inds(tcnt) = ta%t(p)%list(i)%ir
                        ! we can return here, but let's go through remaining MPI_Test calls, just in case
                     endif
                  enddo
               endif
            enddo
         endif

      end function one_more_test

   end subroutine waitall_wrapper

!> \brief Initialize by setting the size of this%r(:) array for non-blocking communication on request.

   subroutine setsize_req(this, nreq, owncomm, label)

      use dataio_pub, only: die
      use MPIF,       only: MPI_COMM_WORLD
      use MPIFUN,     only: MPI_Comm_dup

      implicit none

      class(req_arr),             intent(inout) :: this     !< an object invoking the type-bound procedure
      integer(kind=4),            intent(in)    :: nreq     !< expected maximum number of concurrent MPI requests in non-blocking parts of the code
      logical,                    intent(in)    :: owncomm  !< operate on own communicator, duplicated from MPI_COMM_WORLD
      character(len=*), optional, intent(in)    :: label    !< identification

      integer :: sreq

      if (allocated(this%r)) then
         sreq = size(this%r)
         if (sreq < nreq) deallocate(this%r)
      else
         sreq = 0
         this%owncomm = .false.
         this%comm = MPI_COMM_WORLD  ! just copy it
         this%tags%label = ""
      endif

      if (present(label)) call this%tags%set_label(label)

      if (sreq < nreq) allocate(this%r(nreq))

      if (owncomm) then
         if (this%owncomm) then
            call die("[req_array:setsize_req] communicator already duplicated")
         else
            call MPI_Comm_dup(MPI_COMM_WORLD, this%comm, err_mpi)  ! create a separate duplicate
            this%owncomm = .true.
         endif
      endif
      this%n = 0

   end subroutine setsize_req

!>
!! \brief Double size of this%r(:) array for non-blocking communication on request.
!!
!! \details Perform a resize by a factor of 2. Save existing values stored in this%r(:).
!<

   subroutine doublesize_req(this, owncomm, label)

      use dataio_pub, only: die
#ifdef MPIF08
      use MPIF,       only: MPI_Request
#endif /* MPIF08 */

      implicit none

      class(req_arr),             intent(inout) :: this     !< an object invoking the type-bound procedure
      logical,                    intent(in)    :: owncomm  !< operate on own communicator, duplicated from MPI_COMM_WORLD
      character(len=*), optional, intent(in)    :: label    !< identification

      !< new request array for MPI_Waitall
#ifdef MPIF08
      type(MPI_Request), allocatable, dimension(:) :: new_req
#else /* !MPIF08 */
      integer(kind=4),   allocatable, dimension(:) :: new_req
#endif /* !MPIF08 */
      integer :: sreq
      integer(kind=4), parameter :: default_cnt = 16  ! with no better guess let's allocate some entries

      if (.not. allocated(this%r)) then
         call this%setsize_req(default_cnt, owncomm, label)
      else
         if (present(label)) call die("[req_array:doublesize_req] label for non-fresh req?")
         sreq = size(this%r)
         if (sreq <= 0) call die("[req_array:doublesize_req] size(req) <= 0")

         allocate(new_req(2*sreq))

         new_req(1:sreq) = this%r(:)
         call move_alloc(from=new_req, to=this%r)
      endif

   end subroutine doublesize_req

end module req_array
