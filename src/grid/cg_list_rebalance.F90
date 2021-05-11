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

!> \brief A module with routines for computing load balance and performing necessary cg moves.

module cg_list_rebalance

   implicit none

   private
   public :: rebalance_old

contains

   ! Collect cg costs over whole level
   !
   ! BEWARE: this is quite antiparallel approach
   ! Collecting everything on the master process is way easier than trying to perform global weighted load balance.
   ! In the future we may need to implement an alternative approach employing parallel sorting.
   !
   ! ToDo: reimplement it and base on all_cg instead of level-wise processing.

   subroutine collect_costs

      use cg_level_connected, only: cg_level_connected_t
      use cg_level_finest,    only: finest
      use cg_list,            only: cg_list_element
      use cg_list_balance,    only: I_N_B, I_OFF
      use cg_list_global,     only: all_cg
      use constants,          only: LO, I_ONE, ndims, PPP_AMR
      use dataio_pub,         only: warn
      use load_balance,       only: balance_cg, balance_host, balance_thread
      use mpisetup,           only: err_mpi, req, inflate_req, master, FIRST, LAST
      use MPIF,               only: MPI_INTEGER, MPI_INTEGER8, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, MPI_COMM_WORLD, MPI_Gather, MPI_Recv, MPI_Isend
      use procnames,          only: pnames
      use ppp,                only: ppp_main
      use ppp_mpi,            only: piernik_Waitall

      implicit none

      type(cg_level_connected_t), pointer :: curl

      integer(kind=8), dimension(:,:), allocatable :: gptemp
      real, allocatable, dimension(:) :: costs
      type(cg_list_element), pointer :: cgl
      integer(kind=4) :: p, i, ii
      real :: speed
!      logical :: invalid_speed
!   real, allocatable, dimension(:) :: glob_costs
      enum, bind(C)
         enumerator :: I_GID = I_N_B + ndims
      end enum
      integer(kind=4), parameter :: tag_gpt = 1, tag_cost = tag_gpt + 1  ! also used as counters for requests
      integer, parameter :: ALL_POS = 0
      character(len=*), parameter :: cc_label = "collect_costs"

      call ppp_main%start(cc_label, PPP_AMR)

      curl => finest%level
      do while (associated(curl))

         call inflate_req(tag_cost)

         ! invalid_speed = .false.
         curl%recently_changed = .false.
         allocate(gptemp(I_OFF:I_GID, curl%cnt), costs(ALL_POS:curl%cnt), curl%cnt_all(FIRST:LAST)) !, glob_costs(FIRST:LAST))

         ! We have to use cgl%cg%old_costs because cgl%cg%costs was reset just before the call to refinement update
         i = 0
         cgl => curl%first
         do while (associated(cgl))
            i = i + I_ONE
            gptemp(:, i) = [ cgl%cg%my_se(:, LO), int(cgl%cg%n_b, kind=8), int(cgl%cg%grid_id, kind=8) ]
            costs(i) = cgl%cg%old_costs%total()
            cgl => cgl%nxt
         enddo

         costs(ALL_POS) = 0.
         cgl => all_cg%first
         do while (associated(cgl))
            costs(ALL_POS) = costs(ALL_POS) + cgl%cg%old_costs%total()
            cgl => cgl%nxt
         enddo

         call MPI_Gather(curl%cnt, I_ONE, MPI_INTEGER, curl%cnt_all, I_ONE, MPI_INTEGER, FIRST, MPI_COMM_WORLD, err_mpi)
         if (master) then
            call curl%gp%init(sum(curl%cnt_all))
            ii = 0
            do p = FIRST, LAST
               if (curl%cnt_all(p) > 0) then
                  if (p /= FIRST) then
                     allocate(gptemp(I_OFF:I_GID, curl%cnt_all(p)), costs(ALL_POS:curl%cnt_all(p)))
                     call MPI_Recv(gptemp, size(gptemp, kind=4), MPI_INTEGER8,         p, tag_gpt,  MPI_COMM_WORLD, MPI_STATUS_IGNORE, err_mpi)
                     call MPI_Recv(costs,  size(costs, kind=4),  MPI_DOUBLE_PRECISION, p, tag_cost, MPI_COMM_WORLD, MPI_STATUS_IGNORE, err_mpi)
                  endif

                  ! Apply host or thread speed coefficients to costs(:), if applicable
                  if (balance_host > 0. .and. pnames%speed_avail) then
                     speed = merge(pnames%speed(p), pnames%proc_on_node(pnames%hostindex(p))%speed, balance_thread)
                     if (speed < 0.) call warn("[cg_list_rebalance:rebalance_old:normalize_costs] speed < 0.")
                     if (speed > 0.) then
                        costs(:) = costs(:) / speed
                        ! else
                        ! invalid_speed = .true.
                     endif
                  endif

                  do i = I_ONE, curl%cnt_all(p)
                     if (balance_cg > 0.) then
                        call curl%gp%list(i+ii)%set_gp(gptemp(I_OFF:I_OFF+ndims-1, i), int(gptemp(I_N_B:I_N_B+ndims-1, i), kind=4), int(gptemp(I_GID, i), kind=4), p, costs(i))
                     else
                        call curl%gp%list(i+ii)%set_gp(gptemp(I_OFF:I_OFF+ndims-1, i), int(gptemp(I_N_B:I_N_B+ndims-1, i), kind=4), int(gptemp(I_GID, i), kind=4), p)
                     endif
                  enddo
                  ii = ii + curl%cnt_all(p)
                  ! glob_costs(p) = costs(ALL_POS)
               else
                  ! glob_costs(p) = 0.
               endif
               if (allocated(gptemp)) deallocate(gptemp)
               if (allocated(costs))  deallocate(costs)
               ! if (p == FIRST .or. curl%cnt_all(p) > 0) deallocate(gptemp, costs)  ! should be correct but it is less defensive
            enddo
         else
            if (curl%cnt > 0) then
               call MPI_Isend(gptemp, size(gptemp, kind=4), MPI_INTEGER8,         FIRST, tag_gpt,  MPI_COMM_WORLD, req(tag_gpt),  err_mpi)
               call MPI_Isend(costs,  size(costs,  kind=4), MPI_DOUBLE_PRECISION, FIRST, tag_cost, MPI_COMM_WORLD, req(tag_cost), err_mpi)
               call piernik_Waitall(tag_cost, "rebalance")
            endif
            deallocate(gptemp, costs)
         endif

         ! deallocate(glob_costs)

         curl => curl%coarser
      enddo

      call ppp_main%stop(cc_label, PPP_AMR)

   end subroutine collect_costs

!> \brief Routine for measuring disorder level in distribution of grids across processes and restoring proper balance

   subroutine rebalance_old

      use cg_level_connected, only: cg_level_connected_t
      use cg_level_finest,    only: finest
      use cg_list_dataop,     only: expanded_domain
      use constants,          only: pSUM, PPP_AMR
      use dataio_pub,         only: warn, msg
      use mpisetup,           only: master, FIRST, LAST, piernik_MPI_Allreduce
      use ppp,                only: ppp_main

      implicit none

      type(cg_level_connected_t), pointer :: curl
      real, dimension(FIRST:LAST) :: speed, cumul  ! normalized thread speed and cumulative speed for all threads up to rank
      integer :: hmts
      integer(kind=4) :: edc
      character(len=*), parameter :: ro_label = "rebalance_old"

      call ppp_main%start(ro_label, PPP_AMR)

      call collect_costs

      edc = expanded_domain%cnt
      call piernik_MPI_Allreduce(edc, pSUM)

      hmts = how_many_to_shuffle()

      if (hmts > 0) then
         if (edc /= 0) then
            write(msg, '(a,i5,a)')"[cg_list_rebalance:balance_old] Allreduce(expanded_domain%cnt) = ", edc, ", aborting reshuffling."
            if (master) call warn(msg)
            ! Unfortunately this will bypass thread exclusion
         else
            call reshuffle
         endif
      endif

      ! Deallocate, whatever was allocated by local routines
      curl => finest%level
      do while (associated(curl))

         if (master) call curl%gp%cleanup
         deallocate(curl%cnt_all)

         curl => curl%coarser
      enddo

      call ppp_main%stop(ro_label, PPP_AMR)

   contains

      ! Count how many cg needs to be moved

      integer function how_many_to_shuffle() result(s)

         use constants,  only: INVALID, I_ONE
         use dataio_pub, only: printinfo, die
         use mpisetup,   only: piernik_MPI_Bcast
         use procnames,  only: pnames
         use refinement, only: oop_thr

         implicit none

         integer :: i
         integer(kind=4) :: p
         logical :: rebalance_necessary

         rebalance_necessary = .false.
         s = 0
         if (master) then

            curl => finest%level
            do while (associated(curl))

               if (size(curl%gp%list) > 0) then

                  if (all(pnames%exclude)) call die("[cg_list_rebalance:balance_old] all threads excluded")

                  call curl%gp%set_sort_weight(curl%l%off, .true.)
                  call compute_speed_cumul

                  p = FIRST
                  do i = lbound(curl%gp%list, 1), ubound(curl%gp%list, 1)
                     do while (curl%gp%list(i)%cweight - .5 * curl%gp%list(i)%weight > cumul(p))
                        p = p + I_ONE
                     enddo
                     if (p > LAST) call die("[cg_list_rebalance:balance_old] p > LAST")
                     curl%gp%list(i)%dest_proc = p
                  enddo

                  if (any(curl%gp%list(:)%dest_proc == INVALID)) call die("[cg_list_rebalance:balance_old] not all dest_proc have been set")

                  associate ( cnt_mov => count(curl%gp%list(:)%cur_proc /= curl%gp%list(:)%dest_proc))
                     s = s + cnt_mov
                     if (cnt_mov > 0) then
                        write(msg, '(a,i3,2(a,i6))')"[cg_list_rebalance:balance_old] ^", curl%l%id, " OutOfPlace grids:", cnt_mov, "/",size(curl%gp%list)
                        !a,f6.3,a , " (load balance: ", sum(curl%cnt_all) / real(maxval(curl%cnt_all) * size(curl%cnt_all)), ")"
                        call printinfo(msg)
                     endif
                     if (cnt_mov/real(size(curl%gp%list)) > oop_thr) rebalance_necessary = .true.
                  end associate

               endif

               if (sum(curl%cnt_all, mask = pnames%exclude) /= 0) rebalance_necessary = .true.

               curl => curl%coarser
            enddo

            ! if we exceed per-level threshold or any excluded thread has cg
            if (rebalance_necessary) then
               call printinfo("[cg_list_rebalance:balance_old] reshuffling OutOfPlace grids")
            else
               s = 0
            endif

         endif

         call piernik_MPI_Bcast(s)

      end function how_many_to_shuffle

      subroutine compute_speed_cumul

         use load_balance, only: balance_host, balance_thread
         use procnames,    only: pnames

         implicit none

         real :: cml
         integer(kind=4) :: p

         if (balance_host > 0. .and. pnames%speed_avail) then
            if (balance_thread) then
               speed(:) = pnames%speed(:)
            else
               speed(:) = pnames%proc_on_node(pnames%hostindex(:))%speed
            endif
            if (all(speed(:) > 0. .or. pnames%exclude(:))) then
               where (speed(:) > 0.)
                  speed(:) = 1. / speed(:)
               elsewhere
                  speed(:) = 0.
               endwhere
            else
               speed(:) = 1.
            endif
         else
            speed(:) = 1.
         endif

         ! ToDo here modify speed to account for over/underload on other levels

         where (pnames%exclude(:)) speed(:) = 0.

         cml = 0.
         do p = FIRST, LAST
            cml = cml + speed(p)
            cumul(p) = cml
         enddo

         speed(:) = speed(:) / cml
         cumul(:) = cumul(:) / cml

      end subroutine compute_speed_cumul

   end subroutine rebalance_old

!>
!! \brief Routine for moving existing grids between processes
!!
!! Warning: The allocates of cglepa(i)%tbuf can be pretty huge (totfld * product(cg%n_) * sum_od_incoming_and_outgoing_blocks * sizeof(double)).
!! Pessimistically, the amount of allocated memory can be doubled here afor a brief period of time.
!! ToDo: Make a check for available memory and implement also multi-pass variant of the exchange.
!<

   subroutine reshuffle

      use cg_level_connected, only: cg_level_connected_t
      use cg_level_finest,    only: finest
      use cg_list,            only: cg_list_element
      use cg_list_global,     only: all_cg
      use constants,          only: ndims, LO, HI, I_ZERO, I_ONE, xdim, ydim, zdim, pMAX, PPP_AMR
      use dataio_pub,         only: die
      use grid_cont,          only: grid_container
      use grid_container_ext, only: cg_extptrs
      use list_of_cg_lists,   only: all_lists
      use MPIF,               only: MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, MPI_Isend, MPI_Irecv
      use mpisetup,           only: master, piernik_MPI_Bcast, piernik_MPI_Allreduce, proc, err_mpi, req, inflate_req, tag_ub
      use named_array_list,   only: qna, wna
      use ppp,                only: ppp_main
      use ppp_mpi,            only: piernik_Waitall

      implicit none

      type(cg_level_connected_t), pointer :: curl
      type(cg_list_element), pointer :: cgl
      integer :: s, n_gid, totfld
      integer(kind=4) :: i, p, nr
      integer(kind=8), dimension(:,:), allocatable :: gptemp
      integer(kind=8), dimension(ndims, LO:HI) :: se
      logical :: found
      type(grid_container),  pointer :: cg
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
         enumerator :: I_LEV
         enumerator :: I_C_P
         enumerator :: I_D_P
      end enum
      character(len=*), parameter :: ISR_label = "reshuffle_Isend+Irecv", cp_label = "reshuffle_copy", gp_label = "reshuffle_gptemp"

#ifdef NBODY
      logical :: has_particles

      has_particles = .false.

      curl => finest%level
      do while (associated(curl))
         cgl => curl%first
         do while (associated(cgl))
            if (associated(cgl%cg%pset%first)) has_particles = .true.
            cgl => cgl%nxt
         enddo
         curl => curl%coarser
      enddo
      if (has_particles) call die("[cg_list_rebalance:reshuffle] Particles aren't supported yet")
#endif /* NBODY */

      ! Count the number of fields on any of the top level cg
      totfld = 0
      cgl => finest%level%first
      if (associated(cgl)) then
         do p = lbound(wna%lst, dim=1, kind=4), ubound(wna%lst, dim=1, kind=4)
            if ((.not. only_vital .or. wna%lst(p)%vital) .and. associated(cgl%cg%w(p)%arr)) totfld = totfld + wna%lst(p)%dim4
         enddo
         do p = lbound(qna%lst, dim=1, kind=4), ubound(qna%lst, dim=1, kind=4)
            if ((.not. only_vital .or. qna%lst(p)%vital) .and. associated(cgl%cg%q(p)%arr)) totfld = totfld + 1
         enddo
      endif
      call piernik_MPI_Allreduce(totfld, pMAX)

      s = 0
      if (master) then
         curl => finest%level
         do while (associated(curl))
            s = s + count(curl%gp%list(:)%cur_proc /= curl%gp%list(:)%dest_proc)
            curl => curl%coarser
         enddo
      endif
      call piernik_MPI_Bcast(s)

      call ppp_main%start(gp_label, PPP_AMR)
      allocate(gptemp(I_OFF:I_D_P, s), cglepa(s))
      if (master) then
         p = 0
         curl => finest%level
         do while (associated(curl))
            do i = lbound(curl%gp%list, dim=1, kind=4), ubound(curl%gp%list, dim=1, kind=4)
               if (curl%gp%list(i)%cur_proc /= curl%gp%list(i)%dest_proc) then
                  p = p + I_ONE
                  gptemp(:, p) = [ curl%gp%list(i)%off, int( [ curl%gp%list(i)%n_b, curl%gp%list(i)%cur_gid, curl%l%id, curl%gp%list(i)%cur_proc, curl%gp%list(i)%dest_proc ], kind=8) ]
               endif
            enddo
            curl => curl%coarser
         enddo
      endif
      call piernik_MPI_Bcast(gptemp)
      call ppp_main%stop(gp_label, PPP_AMR)

      ! Irecv & Isend
      call ppp_main%start(ISR_label, PPP_AMR)
      nr = I_ZERO
      curl => finest%level
      do while (associated(curl))

         if (ubound(gptemp, dim=2, kind=4) > tag_ub) call die("[cg_list_rebalance:balance_old] this MPI implementation has too low MPI_TAG_UB attribute")
         do i = lbound(gptemp, dim=2, kind=4), ubound(gptemp, dim=2, kind=4)
            if (gptemp(I_LEV, i) == curl%l%id) then

               cglepa(i)%p => null()
               if (gptemp(I_C_P, i) == gptemp(I_D_P, i)) call die("[cg_list_rebalance:balance_old] can not send to self")
               if (gptemp(I_C_P, i) == proc) then ! send
                  found = .false.
                  cgl => curl%first
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
                  nr = nr + I_ONE
                  if (nr > size(req, dim=1)) call inflate_req
                  call MPI_Isend(cglepa(i)%tbuf, size(cglepa(i)%tbuf, kind=4), MPI_DOUBLE_PRECISION, int(gptemp(I_D_P, i), kind=4), i, MPI_COMM_WORLD, req(nr), err_mpi)
               endif
               if (gptemp(I_D_P, i) == proc) then ! receive
                  n_gid = 1
                  if (associated(curl%last)) n_gid = curl%last%cg%grid_id + 1
                  se(:, LO) = gptemp(I_OFF:I_OFF+zdim-1, i)
                  se(:, HI) = gptemp(I_OFF:I_OFF+zdim-1, i) + gptemp(I_N_B:I_N_B+zdim-1, i) - 1
                  curl%recently_changed = .true.
                  call curl%add
                  cglepa(i)%p => curl%last
                  cgl => cglepa(i)%p
                  call curl%last%cg%init_gc(se, n_gid, curl%l)
                  allocate(cglepa(i)%tbuf(totfld, cgl%cg%n_(xdim), cgl%cg%n_(ydim), cgl%cg%n_(zdim)))
                  do p = lbound(cg_extptrs%ext, dim=1, kind=4), ubound(cg_extptrs%ext, dim=1, kind=4)
                     if (associated(cg_extptrs%ext(p)%init))  call cg_extptrs%ext(p)%init(curl%last%cg)
                  enddo
                  call all_cg%add(curl%last%cg)
                  nr = nr + I_ONE
                  if (nr > size(req, dim=1)) call inflate_req
                  call MPI_Irecv(cglepa(i)%tbuf, size(cglepa(i)%tbuf, kind=4), MPI_DOUBLE_PRECISION, int(gptemp(I_C_P, i), kind=4), i, MPI_COMM_WORLD, req(nr), err_mpi)
               endif
            endif
         enddo

         curl => curl%coarser
      enddo
      call ppp_main%stop(ISR_label, PPP_AMR)

      call piernik_Waitall(nr, "reshuffle", PPP_AMR)

      call ppp_main%start(cp_label, PPP_AMR)
      curl => finest%level
      do while (associated(curl))

         do i = lbound(gptemp, dim=2, kind=4), ubound(gptemp, dim=2, kind=4)
            if (gptemp(I_LEV, i) == curl%l%id) then
               cgl => cglepa(i)%p
               if (gptemp(I_C_P, i) == proc) then ! cleanup
                  deallocate(cglepa(i)%tbuf)
                  cg => cgl%cg
                  call all_lists%forget(cg)
                  curl%recently_changed = .true.
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
            endif
         enddo

         curl => curl%coarser
      enddo
      call ppp_main%stop(cp_label, PPP_AMR)

      deallocate(gptemp, cglepa)

      curl => finest%level
      do while (associated(curl))
         call curl%sort_SFC
         curl => curl%coarser
      enddo

   end subroutine reshuffle

end module cg_list_rebalance
