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

#ifdef MPIF08
   use MPIF, only: MPI_Request, MPI_Comm
#endif /* MPIF08 */

   implicit none

   private
   public :: req_arr

   type :: req_arr  ! array of requests
      !< request array for MPI_Waitall
#ifdef MPIF08
      type(MPI_Request), allocatable, dimension(:) :: r     !< array of requests
      type(MPI_Comm)                               :: comm  !< own communicator
#else /* !MPIF08 */
      integer(kind=4),   allocatable, dimension(:) :: r
      integer(kind=4)                              :: comm
#endif /* !MPIF08 */
      integer(kind=4) :: n  !< number of requests
      logical :: owncomm    !< was the own communicator requested?
   contains
      procedure :: nxt              !< give next unused request
      procedure :: waitall_on_some  !< wait for completion of the requests on selected procs
      procedure :: setsize_req      !< set the storage for this%r
      procedure :: doublesize_req   !< double the storage size for this%r
      generic, public :: init => doublesize_req, setsize_req
      generic, public :: waitall => waitall_on_some
   end type req_arr

   integer(kind=4) :: err_mpi  !< error status

contains

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
      if (this%n > ubound(this%r, 1)) call this%doublesize_req
      nxt => this%r(this%n)

   end function nxt

!>
!! \brief Simple wrapper for MPI_Waitall
!!
!! For PPP-instrumented variant see req_ppp::waitall_ppp (in ppp_mpi module)
!<

   subroutine waitall_on_some(this)

      use constants,    only: LONG
      use mpi_wrappers, only: C_REQS, req_wall
      use MPIF,         only: MPI_STATUSES_IGNORE
      use MPIFUN,       only: MPI_Waitall, MPI_Comm_free

      implicit none

      class(req_arr), intent(inout) :: this  !< an object invoking the type-bound procedure

      if (this%n > 0) then

         call req_wall%add(int([C_REQS]), int(this%n, kind=LONG))
         call MPI_Waitall(this%n, this%r(:this%n), MPI_STATUSES_IGNORE, err_mpi)
         this%n = 0

      endif

      if (allocated(this%r)) deallocate(this%r)
      if (this%owncomm) then
         call MPI_Comm_free(this%comm, err_mpi)
         this%owncomm = .false.
      endif

   end subroutine waitall_on_some

!> \brief Initialize by setting the size of this%r(:) array for non-blocking communication on request.

   subroutine setsize_req(this, nreq, owncomm)

      use dataio_pub, only: die
      use MPIF,       only: MPI_COMM_WORLD
      use MPIFUN,     only: MPI_Comm_dup

      implicit none

      class(req_arr),    intent(inout) :: this     !< an object invoking the type-bound procedure
      integer(kind=4),   intent(in)    :: nreq     !< expected maximum number of concurrent MPI requests in non-blocking parts of the code
      logical, optional, intent(in)    :: owncomm  !< operate on own communicator, duplicated from MPI_COMM_WORLD

      integer :: sreq
      logical :: oc

      oc = .false.
      if (present(owncomm)) oc = owncomm

      if (allocated(this%r)) then
         sreq = size(this%r)
         if (sreq < nreq) deallocate(this%r)
      else
         sreq = 0
         this%owncomm = .false.
         this%comm = MPI_COMM_WORLD  ! just copy it
      endif

      if (sreq < nreq) allocate(this%r(nreq))

      if (oc) then
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

   subroutine doublesize_req(this, owncomm)

      use dataio_pub, only: die
#ifdef MPIF08
      use MPIF,       only: MPI_Request
#endif /* MPIF08 */

      implicit none

      class(req_arr),    intent(inout) :: this  !< an object invoking the type-bound procedure
      logical, optional, intent(in)    :: owncomm  !< operate on own communicator, duplicated from MPI_COMM_WORLD

      !< new request array for MPI_Waitall
#ifdef MPIF08
      type(MPI_Request), allocatable, dimension(:) :: new_req
#else /* !MPIF08 */
      integer(kind=4),   allocatable, dimension(:) :: new_req
#endif /* !MPIF08 */
      integer :: sreq
      integer(kind=4), parameter :: default_cnt = 16  ! with no better guess let's allocate some entries

      if (.not. allocated(this%r)) then
         call this%setsize_req(default_cnt, owncomm)
      else
         sreq = size(this%r)
         if (sreq <= 0) call die("[req_array:doublesize_req] size(req) <= 0")

         allocate(new_req(2*sreq))

         new_req(1:sreq) = this%r(:)
         call move_alloc(from=new_req, to=this%r)
      endif

   end subroutine doublesize_req

end module
