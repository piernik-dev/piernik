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
   use MPIF, only: MPI_Request
#endif /* MPIF08 */

   implicit none

   private
   public :: req_arr

   type :: req_arr  ! array of requests
      !< request array for MPI_Waitall
#ifdef MPIF08
      type(MPI_Request), allocatable, dimension(:) :: r
#else /* !MPIF08 */
      integer(kind=4),   allocatable, dimension(:) :: r
#endif /* !MPIF08 */
      integer(kind=4) :: n
   contains
      procedure :: nxt             !< give next unused request
      procedure :: waitall         !< wait for completion of the requests and do cleanup
      procedure :: setsize_req     !< set the storage for this%r
      procedure :: doublesize_req  !< double the storage size for this%r
      generic, public :: init => doublesize_req, setsize_req
   end type req_arr

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
!! \brief a PPP wrapper for MPI_Waitall
!!
!! BEWARE: Cannot use extra_barriers when this routine is called only by a subset of MPI ranks.
!<

   subroutine waitall(this, ppp_label, x_mask)

      use constants,    only: PPP_MPI, LONG
      use mpi_wrappers, only: piernik_MPI_Barrier, extra_barriers, C_REQA, req_wall
      use mpisetup,     only: err_mpi
      use MPIF,         only: MPI_STATUSES_IGNORE
      use MPIFUN,       only: MPI_Waitall
      use ppp,          only: ppp_main

      implicit none

      class(req_arr),            intent(inout) :: this       !< an object invoking the type-bound procedure
      character(len=*),          intent(in)    :: ppp_label  !< identifier for PPP entry
      integer(kind=4), optional, intent(in)    :: x_mask     !< extra mask, if necessary

      character(len=*), parameter :: mpiw = "req:MPI_Waitall:"
      integer(kind=4) :: mask

      if (this%n > 0) then
         mask = PPP_MPI
         if (present(x_mask)) mask = mask + x_mask
         call ppp_main%start(mpiw // ppp_label, mask)

         call req_wall%add(int([C_REQA]), int(this%n, kind=LONG))

         call MPI_Waitall(this%n, this%r(:this%n), MPI_STATUSES_IGNORE, err_mpi)
         this%n = 0

         call ppp_main%stop(mpiw // ppp_label, mask)
      endif

      if (allocated(this%r)) deallocate(this%r)

      if (extra_barriers) call piernik_MPI_Barrier

   end subroutine waitall

!> \brief Set size of this%r(:) array for non-blocking communication on request.

   subroutine setsize_req(this, nreq)

      implicit none

      class(req_arr),  intent(inout) :: this  !< an object invoking the type-bound procedure
      integer(kind=4), intent(in)    :: nreq  !< expected maximum number of concurrent MPI requests in non-blocking parts of the code

      integer :: sreq

      if (allocated(this%r)) then
         sreq = size(this%r)
         if (sreq < nreq) deallocate(this%r)
      else
         sreq = 0
      endif

      if (sreq < nreq) allocate(this%r(nreq))
      this%n = 0

   end subroutine setsize_req

!>
!! \brief Double size of this%r(:) array for non-blocking communication on request.
!!
!! \details Perform a resize by a factor of 2. Save existing values stored in this%r(:).
!<

   subroutine doublesize_req(this)

      use dataio_pub, only: die
#ifdef MPIF08
      use MPIF,       only: MPI_Request
#endif /* MPIF08 */

      implicit none

      class(req_arr), intent(inout) :: this  !< an object invoking the type-bound procedure

      !< new request array for MPI_Waitall
#ifdef MPIF08
      type(MPI_Request), allocatable, dimension(:) :: new_req
#else /* !MPIF08 */
      integer(kind=4),   allocatable, dimension(:) :: new_req
#endif /* !MPIF08 */
      integer :: sreq
      integer(kind=4), parameter :: default_cnt = 16  ! with no better guess let's allocate some entries

      if (.not. allocated(this%r)) then
         call this%setsize_req(default_cnt)
      else
         sreq = size(this%r)
         if (sreq <= 0) call die("[req_array:doublesize_req] size(req) <= 0")

         allocate(new_req(2*sreq))

         new_req(1:sreq) = this%r(:)
         call move_alloc(from=new_req, to=this%r)
      endif

   end subroutine doublesize_req

end module
