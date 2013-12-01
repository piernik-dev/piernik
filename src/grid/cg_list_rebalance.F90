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

!> \brief A module with an abstract type created to take out some load-balance related code from cg_level (old grids)

module cg_list_rebalance

   use cg_list_balance, only: cg_list_balance_T
   use constants,       only: ndims

   implicit none

   private
   public :: cg_list_rebalance_T

   !> \brief An abstract type created to take out some load-balance related code from cg_level (old grids)

   type, extends(cg_list_balance_T), abstract :: cg_list_rebalance_T
      integer(kind=4)                   :: level_id  !< level number (relative to base level). No arithmetic should depend on it.
      integer(kind=8), dimension(ndims) :: n_d       !< maximum number of grid cells in each direction (size of fully occupied level)
   contains
      procedure          :: rebalance_old  !< Routine for measuring disorder level in distribution of grids across processes
      procedure, private :: reshuffle      !< Routine for moving existing grids between processes
   end type cg_list_rebalance_T

contains

!> \brief Routine for measuring disorder level in distribution of grids across processes

!#define DEBUG
   subroutine rebalance_old(this)

      use cg_list,         only: cg_list_element
      use cg_list_balance, only: I_N_B, I_OFF
      use cg_list_dataop,  only: expanded_domain
      use constants,       only: LO, HI, I_ONE, pSUM, ndims
      use dataio_pub,      only: warn, msg, printinfo
      use mpisetup,        only: master, FIRST, LAST, nproc, piernik_MPI_Bcast, piernik_MPI_Allreduce
      use refinement,      only: oop_thr
      use sort_piece_list, only: grid_piece_list
#ifdef DEBUG
      use mpi,             only: MPI_INTEGER, MPI_INTEGER8, MPI_STATUS_IGNORE
      use mpisetup,        only: comm, mpi_err
#endif /* DEBUG */

      implicit none

      class(cg_list_rebalance_T), intent(inout) :: this

      integer(kind=4), dimension(FIRST:LAST) :: cnt_existing
      type(grid_piece_list) :: gp
      type(cg_list_element), pointer :: cgl
      integer :: s
      integer(kind=4) :: i, p
      integer(kind=8), dimension(:,:), allocatable :: gptemp
      enum, bind(C)
         enumerator :: I_GID = I_N_B + ndims
      end enum
#ifdef DEBUG
      integer(kind=4), parameter :: tag_gpt = 1
#else /* !DEBUG */
      integer(kind=4) :: ii
#endif /* DEBUG */

      this%recently_changed = .false.
      allocate(gptemp(I_OFF:I_GID, this%cnt))
      i = 0
      cgl => this%first
      do while (associated(cgl))
         i = i + I_ONE
         gptemp(:, i) = [ cgl%cg%my_se(:, LO), int(cgl%cg%n_b, kind=8), int(cgl%cg%grid_id, kind=8) ]
         cgl => cgl%nxt
      enddo
#ifdef DEBUG
      ! Gather complete grid list and compare with this%dot%gse
      call MPI_Gather(this%cnt, I_ONE, MPI_INTEGER, cnt_existing, I_ONE, MPI_INTEGER, FIRST, comm, mpi_err)
      if (master) then
         call gp%init(sum(cnt_existing))
         do i = I_ONE, this%cnt
            call gp%list(i)%set_gp(gptemp(I_OFF:I_OFF+ndims-1, i), int(gptemp(I_N_B:I_N_B+ndims-1, i), kind=4), int(gptemp(I_GID, i), kind=4), FIRST)
            if (any(this%dot%gse(FIRST)%c(i)%se(:, LO) /= gp%list(i)%off) .or. gp%list(i)%cur_gid /= i .or. &
                 any(this%dot%gse(FIRST)%c(i)%se(:, HI) - this%dot%gse(FIRST)%c(i)%se(:, LO) +1 /= gp%list(i)%n_b)) &
                 call warn("cl:bo this%dot%gse(FIRST) /= gptemp")
         enddo
         deallocate(gptemp)
         s = this%cnt
         do p = FIRST + 1, LAST
            if (cnt_existing(p) > 0) then
               allocate(gptemp(I_OFF:I_GID, cnt_existing(p)))
               call MPI_Recv(gptemp, size(gptemp), MPI_INTEGER8, p, tag_gpt, comm, MPI_STATUS_IGNORE, mpi_err)
               do i = I_ONE, cnt_existing(p)
                  call gp%list(i+s)%set_gp(gptemp(I_OFF:I_OFF+ndims-1, i), int(gptemp(I_N_B:I_N_B+ndims-1, i), kind=4), int(gptemp(I_GID, i), kind=4), p)
                  if (any(this%dot%gse(p)%c(i)%se(:, LO) /= gp%list(i+s)%off) .or. gp%list(i+s)%cur_gid /= i .or. &
                       any(this%dot%gse(p)%c(i)%se(:, HI) - this%dot%gse(p)%c(i)%se(:, LO) +1 /= gp%list(i+s)%n_b)) &
                       call warn("cl:bo this%dot%gse(p) /= gptemp")
               enddo
               s = s + cnt_existing(p)
               deallocate(gptemp)
            endif
         enddo
      else
         if (this%cnt > 0) call MPI_Send(gptemp, size(gptemp), MPI_INTEGER8, FIRST, tag_gpt, comm, mpi_err)
      endif
#else /* !DEBUG */
      ! Trust that this%dot%gse is updated
      if (master) then
         do p = FIRST, LAST
            cnt_existing(p) = size(this%dot%gse(p)%c, kind=4)
         enddo
         call gp%init(sum(cnt_existing))
         i = 0
         do p = FIRST, LAST
            ii = i
            if (ii /= sum(cnt_existing(:p-1))) call warn("cl:bo ii /= sum(cnt_existing(:p-1))")
            do s = lbound(this%dot%gse(p)%c, dim=1), ubound(this%dot%gse(p)%c, dim=1)
               i = i + I_ONE
               call gp%list(i)%set_gp(this%dot%gse(p)%c(s)%se(:, LO), int(this%dot%gse(p)%c(s)%se(:, HI) - this%dot%gse(p)%c(s)%se(:, LO) +1, kind=4), i - ii, p)
            enddo
         enddo
      endif
#endif /* DEBUG */
      if (allocated(gptemp)) deallocate(gptemp)

      if (master) then
         call gp%set_id(this%off)
         call gp%sort
         do p = FIRST, LAST
            gp%list(p*size(gp%list)/nproc+1 : ((p+1)*size(gp%list))/nproc)%dest_proc = p
         enddo
         s = 0
         if (size(gp%list) > 0) then
            s = count(gp%list(:)%cur_proc /= gp%list(:)%dest_proc)
            if (s/real(size(gp%list)) > oop_thr) then
               write(msg,'(a,i3,2(a,i6),a,f6.3,a)')"[cg_list_rebalance:balance_old] ^", this%level_id," Reshuffling OutOfPlace grids:",s, "/",size(gp%list)," (load balance: ",sum(cnt_existing)/real(maxval(cnt_existing)*size(cnt_existing)),")"
               call printinfo(msg)
            else
               s = 0
            endif
         endif
      endif

      call piernik_MPI_Bcast(s)
      if (s>0) then
         p = expanded_domain%cnt
         call piernik_MPI_Allreduce(p, pSUM)
         if (p /= 0) then
            write(msg,'(a,i5,a)')"[cg_list_rebalance:balance_old] Allreduce(expanded_domain%cnt) = ",p,", aborting reshuffling."
            if (master) call warn(msg)
         else
            call this%reshuffle(gp)
         endif
      endif

      if (master) call gp%cleanup

   end subroutine rebalance_old

!> \brief for moving existing grids between processes

   subroutine reshuffle(this, gp)

      use grid_container_ext, only: cg_extptrs
      use cg_list,            only: cg_list_element
      use cg_list_global,     only: all_cg
      use constants,          only: ndims, LO, HI, I_ONE, xdim, ydim, zdim, pMAX
      use dataio_pub,         only: die
      use grid_cont,          only: grid_container
      use list_of_cg_lists,   only: all_lists
      use mpi,                only: MPI_DOUBLE_PRECISION
      use mpisetup,           only: master, piernik_MPI_Bcast, piernik_MPI_Allreduce, proc, comm, mpi_err, req, status, inflate_req
      use named_array_list,   only: qna, wna
      use sort_piece_list,    only: grid_piece_list

      implicit none

      class(cg_list_rebalance_T), intent(inout) :: this
      type(grid_piece_list),      intent(in)    :: gp

      type(cg_list_element), pointer :: cgl
      integer :: s, n_gid, totfld, nr
      integer(kind=4) :: i, p
      integer(kind=8), dimension(:,:), allocatable :: gptemp
      integer(kind=8), dimension(ndims, LO:HI) :: se
      logical :: found
      type(grid_container),  pointer :: cg
      integer(kind=4), dimension(:,:), pointer :: mpistatus
      type :: cglep
         type(cg_list_element), pointer :: p
         real, dimension(:,:,:,:), allocatable :: tbuf
      end type cglep
      type(cglep), allocatable, dimension(:) :: cglepa
      logical, parameter :: only_vital = .false. ! set to true to minimize the amount of data to be transferred, may result in improper calculation of error in maclaurin test
      !> \todo measure how much it costs in reality
      enum, bind(C)
         enumerator :: I_OFF
         enumerator :: I_N_B = I_OFF + ndims
         enumerator :: I_GID = I_N_B + ndims
         enumerator :: I_C_P = I_GID + I_ONE
         enumerator :: I_D_P = I_C_P + I_ONE
      end enum

      totfld = 0
      cgl => this%first
      if (associated(cgl)) then
         do p = lbound(wna%lst, dim=1, kind=4), ubound(wna%lst, dim=1, kind=4)
            if ((.not. only_vital .or. wna%lst(p)%vital) .and. associated(cgl%cg%w(p)%arr)) then
               ! associated(cgl%cg%w(p)%arr)) .eqv. (this%level_id >= base_level_id) .or. wna%lst(p)%multigrid ?
               totfld = totfld + wna%lst(p)%dim4
            endif
         enddo
         do p = lbound(qna%lst, dim=1, kind=4), ubound(qna%lst, dim=1, kind=4)
            if ((.not. only_vital .or. qna%lst(p)%vital) .and. associated(cgl%cg%q(p)%arr)) then
               totfld = totfld + 1
            endif
         enddo
      endif
      call piernik_MPI_Allreduce(totfld, pMAX)
      ! communicate gp%list
      if (master) s = count(gp%list(:)%cur_proc /= gp%list(:)%dest_proc)
      call piernik_MPI_Bcast(s)
      allocate(gptemp(I_OFF:I_D_P, s))
      if (master) then
         p = 0
         do i = lbound(gp%list, dim=1, kind=4), ubound(gp%list, dim=1, kind=4)
            if (gp%list(i)%cur_proc /= gp%list(i)%dest_proc) then
               p = p + I_ONE
               gptemp(:, p) = [ gp%list(i)%off, int( [ gp%list(i)%n_b, gp%list(i)%cur_gid, gp%list(i)%cur_proc, gp%list(i)%dest_proc ], kind=8) ]
            endif
         enddo
      endif
      call piernik_MPI_Bcast(gptemp)

      ! Irecv & Isend
      nr = 0
      allocate(cglepa(size(gptemp)))
      do i = lbound(gptemp, dim=2, kind=4), ubound(gptemp, dim=2, kind=4)
         cglepa(i)%p => null()
         if (gptemp(I_C_P, i) == gptemp(I_D_P, i)) call die("[cg_list_rebalance:balance_old] can not send to self")
         if (gptemp(I_C_P, i) == proc) then ! send
            found = .false.
            cgl => this%first
            do while (associated(cgl))
               if (cgl%cg%grid_id == gptemp(I_GID,i)) then
                  found = .true.
                  cglepa(i)%p => cgl
                  allocate(cglepa(i)%tbuf(totfld, cgl%cg%n_(xdim), cgl%cg%n_(ydim), cgl%cg%n_(zdim)))
                  ! We communicate blocks with guardcells because current implementation of magnetic field evolves external guardcells in a way that makes it impossible
                  ! to reconstruct them from scratch. This results in much larger messages to communicate, but we won't need to call guardcell exchange afterwards.
                  s = lbound(cglepa(i)%tbuf, dim=1)
                  do p = lbound(wna%lst, dim=1, kind=4), ubound(wna%lst, dim=1, kind=4)
                     if ((.not. only_vital .or. wna%lst(p)%vital) .and. associated(cgl%cg%w(p)%arr)) then ! not associated for multigrid coarse levels
                        cglepa(i)%tbuf(s:s+wna%lst(p)%dim4-1, :, :, :) = cgl%cg%w(p)%arr(:, :, :, :)
                        s = s + wna%lst(p)%dim4
                     endif
                  enddo
                  do p = lbound(qna%lst, dim=1, kind=4), ubound(qna%lst, dim=1, kind=4)
                     if ((.not. only_vital .or. qna%lst(p)%vital) .and. associated(cgl%cg%q(p)%arr)) then
                        cglepa(i)%tbuf(s, :, :, :) = cgl%cg%q(p)%arr(:, :, :)
                        s = s + 1
                     endif
                  enddo
                  exit
               endif
               cgl => cgl%nxt
            enddo
            if (.not. found) call die("[cg_list_rebalance:balance_old] Grid id not found")
            nr = nr + 1
            if (nr > size(req, dim=1)) call inflate_req
            call MPI_Isend(cglepa(i)%tbuf, size(cglepa(i)%tbuf), MPI_DOUBLE_PRECISION, gptemp(I_D_P, i), i, comm, req(nr), mpi_err)
         endif
         if (gptemp(I_D_P, i) == proc) then ! receive
            n_gid = 1
            if (associated(this%last)) n_gid = this%last%cg%grid_id + 1
            se(:, LO) = gptemp(I_OFF:I_OFF+zdim-1, i)
            se(:, HI) = gptemp(I_OFF:I_OFF+zdim-1, i) + gptemp(I_N_B:I_N_B+zdim-1, i) - 1
            this%recently_changed = .true.
            call this%add
            cglepa(i)%p => this%last
            cgl => cglepa(i)%p
            call this%last%cg%init(this%n_d, this%off, se, n_gid, this%level_id) ! we cannot pass "this" as an argument because of circular dependencies
            allocate(cglepa(i)%tbuf(totfld, cgl%cg%n_(xdim), cgl%cg%n_(ydim), cgl%cg%n_(zdim)))
            do p = lbound(cg_extptrs%ext, dim=1, kind=4), ubound(cg_extptrs%ext, dim=1, kind=4)
               if (associated(cg_extptrs%ext(p)%init))  call cg_extptrs%ext(p)%init(this%last%cg)
            enddo
            call all_cg%add(this%last%cg)
            nr = nr + 1
            if (nr > size(req, dim=1)) call inflate_req
            call MPI_Irecv(cglepa(i)%tbuf, size(cglepa(i)%tbuf), MPI_DOUBLE_PRECISION, gptemp(I_C_P, i), i, comm, req(nr), mpi_err)
         endif
      enddo

      if (nr > 0) then
         mpistatus => status(:, :nr)
         call MPI_Waitall(nr, req(:nr), mpistatus, mpi_err)
      endif
      do i = lbound(gptemp, dim=2, kind=4), ubound(gptemp, dim=2, kind=4)
         cgl => cglepa(i)%p
         if (gptemp(I_C_P, i) == proc) then ! cleanup
            deallocate(cglepa(i)%tbuf)
            cg => cgl%cg
            call all_lists%forget(cg)
            this%recently_changed = .true.
         endif
         if (gptemp(I_D_P, i) == proc) then ! copy received
            s = lbound(cglepa(i)%tbuf, dim=1)
            do p = lbound(wna%lst, dim=1, kind=4), ubound(wna%lst, dim=1, kind=4)
               if ((.not. only_vital .or. wna%lst(p)%vital) .and. associated(cgl%cg%w(p)%arr)) then
                  cgl%cg%w(p)%arr(:, :, :, :) = cglepa(i)%tbuf(s:s+wna%lst(p)%dim4-1, :, :, :)
                  s = s + wna%lst(p)%dim4
               endif
            enddo
            do p = lbound(qna%lst, dim=1, kind=4), ubound(qna%lst, dim=1, kind=4)
               if ((.not. only_vital .or. qna%lst(p)%vital) .and. associated(cgl%cg%q(p)%arr)) then
                  cgl%cg%q(p)%arr(:, :, :) = cglepa(i)%tbuf(s, :, :, :)
                  s = s + 1
               endif
            enddo
            deallocate(cglepa(i)%tbuf)
         endif
      enddo
      deallocate(cglepa)

      s = 0
      cgl => this%first
      do while (associated(cgl))
         s = s + 1
         cgl%cg%grid_id = s
         cgl => cgl%nxt
      enddo

      deallocate(gptemp)

   end subroutine reshuffle

end module cg_list_rebalance
