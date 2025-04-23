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

module rebalance

   implicit none

   private
   public :: rebalance_all

contains

   ! Collect cg costs over whole level
   !
   ! BEWARE: this is quite antiparallel approach
   ! Collecting everything on the master process is way easier than trying to perform global weighted load balance.
   ! In the future we may need to implement an alternative approach employing parallel sorting.
   !
   ! ToDo: reimplement it and base on all_cg instead of level-wise processing.

   subroutine collect_costs

      use barrier,            only: extra_barriers
      use cg_level_connected, only: cg_level_connected_t
      use cg_level_finest,    only: finest
      use cg_list,            only: cg_list_element
      use cg_list_balance,    only: I_N_B, I_OFF
      use constants,          only: LO, I_ONE, ndims, PPP_AMR
      use dataio_pub,         only: warn
      use isend_irecv,        only: piernik_Isend
      use load_balance,       only: balance_cg, balance_host, balance_thread, cost_mask
      use mpisetup,           only: err_mpi, master, FIRST, LAST
      use MPIF,               only: MPI_INTEGER, MPI_INTEGER8, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, MPI_COMM_WORLD
      use MPIFUN,             only: MPI_Gather, MPI_Recv
      use procnames,          only: pnames
      use ppp,                only: ppp_main
      use pppmpi,             only: req_ppp

      implicit none

      type(cg_level_connected_t), pointer :: curl
      type(req_ppp) :: req
      integer(kind=8), dimension(:,:), allocatable :: gptemp
      real, allocatable, dimension(:) :: costs
      type(cg_list_element), pointer :: cgl
      integer(kind=4) :: p, i, ii
      real :: wtime
!      logical :: invalid_speed
      enum, bind(C)
         enumerator :: I_GID = I_N_B + ndims
         enumerator :: I_NP = I_GID + I_ONE
      end enum
      integer(kind=4), parameter :: tag_gpt = 1, tag_cost = tag_gpt + 1  ! also used as counters for requests
      character(len=*), parameter :: cc_label = "collect_costs"

      call ppp_main%start(cc_label, PPP_AMR)

      curl => finest%level
      do while (associated(curl))

         ! invalid_speed = .false.
         curl%recently_changed = .false.
         allocate(gptemp(I_OFF:I_NP, curl%cnt), costs(curl%cnt), curl%cnt_all(FIRST:LAST))

         ! We have to use cgl%cg%old_costs because cgl%cg%costs was reset just before the call to refinement update
         i = 0
         cgl => curl%first
         do while (associated(cgl))
            i = i + I_ONE

            gptemp(:, i) = [ cgl%cg%my_se(:, LO), int(cgl%cg%n_b, kind=8), int(cgl%cg%grid_id, kind=8), int(cgl%cg%count_all_particles(), kind=8) ]
            costs(i) = sum(cgl%cg%old_costs%wtime, mask=cost_mask(:))
            cgl => cgl%nxt
         enddo

         call MPI_Gather(curl%cnt, I_ONE, MPI_INTEGER, curl%cnt_all, I_ONE, MPI_INTEGER, FIRST, MPI_COMM_WORLD, err_mpi)
         if (master) then
            call curl%gp%init(sum(curl%cnt_all))
            ii = 0
            do p = FIRST, LAST
               if (curl%cnt_all(p) > 0) then
                  if (p /= FIRST) then
                     allocate(gptemp(I_OFF:I_NP, curl%cnt_all(p)), costs(curl%cnt_all(p)))
                     call MPI_Recv(gptemp, size(gptemp, kind=4), MPI_INTEGER8,         p, tag_gpt,  MPI_COMM_WORLD, MPI_STATUS_IGNORE, err_mpi)
                     call MPI_Recv(costs,  size(costs, kind=4),  MPI_DOUBLE_PRECISION, p, tag_cost, MPI_COMM_WORLD, MPI_STATUS_IGNORE, err_mpi)
                  endif

                  ! Apply host or thread speed coefficients to costs(:), if applicable
                  if (balance_host > 0. .and. pnames%speed_avail) then
                     wtime = merge(pnames%wtime(p), pnames%proc_on_node(pnames%hostindex(p))%wtime, balance_thread)
                     if (wtime < 0.) call warn("[rebalance:collect_costs] wtime < 0.")
                     if (wtime > 0.) then
                        costs(:) = costs(:) / wtime
                        ! else
                        ! invalid_speed = .true.
                     endif
                  endif

                  do i = I_ONE, curl%cnt_all(p)
                     if (balance_cg > 0. .and. pnames%speed_avail) then
                        call curl%gp%list(i+ii)%set_gp(gptemp(I_OFF:I_OFF+ndims-1, i), int(gptemp(I_N_B:I_N_B+ndims-1, i), kind=4), int(gptemp(I_GID, i), kind=4), p, costs(i))
                     else
                        call curl%gp%list(i+ii)%set_gp(gptemp(I_OFF:I_OFF+ndims-1, i), int(gptemp(I_N_B:I_N_B+ndims-1, i), kind=4), int(gptemp(I_GID, i), kind=4), p)
                     endif
                     curl%gp%list(i+ii)%n_part = int(gptemp(I_NP, i), kind=4)  ! should go to set_gp
                  enddo
                  ii = ii + curl%cnt_all(p)
               endif
               if (allocated(gptemp)) deallocate(gptemp)
               if (allocated(costs))  deallocate(costs)
               ! if (p == FIRST .or. curl%cnt_all(p) > 0) deallocate(gptemp, costs)  ! should be correct but it is less defensive
            enddo
         else
            if (curl%cnt > 0) then
               call req%init(tag_cost, owncomm = .false., label = "rebalance:costs")
               call piernik_Isend(gptemp, size(gptemp, kind=4), MPI_INTEGER8,         FIRST, tag_gpt,  req)
               call piernik_Isend(costs,  size(costs,  kind=4), MPI_DOUBLE_PRECISION, FIRST, tag_cost, req)

               if (extra_barriers) then
                  call req%waitall  ! waitall_on_some, without MPI_Barrier call, without PPP
               else
                  call req%waitall("rebalance")  ! we can have PPP if there are no extra MPI_Barrier call
               endif
            endif
            deallocate(gptemp, costs)
         endif

         curl => curl%coarser
      enddo

      call ppp_main%stop(cc_label, PPP_AMR)

   end subroutine collect_costs

!> \brief Routine for measuring disorder level in distribution of grids across processes and restoring proper balance

   subroutine rebalance_all

      use allreduce,          only: piernik_MPI_Allreduce
      use bcast,              only: piernik_MPI_Bcast
      use cg_level_connected, only: cg_level_connected_t
      use cg_level_finest,    only: finest
      use cg_list_dataop,     only: expanded_domain
      use constants,          only: pSUM, PPP_AMR
      use dataio_pub,         only: warn, msg
      use mpisetup,           only: master, FIRST, LAST
      use ppp,                only: ppp_main
      use refinement,         only: is_blocky

      implicit none

      type(cg_level_connected_t), pointer :: curl
      real, dimension(FIRST:LAST) :: speed, cumul, load_fac  ! normalized thread speed and cumulative speed for all threads up to rank
      integer :: hmts
      integer(kind=4) :: edc
      character(len=*), parameter :: ro_label = "rebalance"

      if (.not. is_blocky()) return  !< Rebalancing non-blocky decompositions should be done in totally different way. Currently of no practical interest.

      call ppp_main%start(ro_label, PPP_AMR)

      call collect_costs

      edc = expanded_domain%cnt
      call piernik_MPI_Allreduce(edc, pSUM)

      hmts = 0
      if (master) hmts = how_many_to_shuffle()
      call piernik_MPI_Bcast(hmts)

      if (hmts > 0) then
         if (edc /= 0) then
            write(msg, '(a,i5,a)')"[rebalance:rebalance_all] Allreduce(expanded_domain%cnt) = ", edc, ", aborting reshuffling."
            if (master) call warn(msg)
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

         use cg_level_base,     only: base
         use cg_level_coarsest, only: coarsest
         use constants,         only: fmt_len, INVALID, I_ONE, base_level_id, V_VERBOSE
         use dataio_pub,        only: printinfo, die
         use load_balance,      only: balance_cg, balance_levels, oop_thr
         use mpisetup,          only: slave
         use procnames,         only: pnames

         implicit none

         integer :: i
         integer(kind=4) :: p
         logical :: rebalance_necessary, fine_processed
         integer, dimension(coarsest%level%l%id:finest%level%l%id) :: cnt_mv, cnt_gp
         real, dimension(FIRST:LAST) :: costs_above
         real :: global_costs_above
         character(len=fmt_len) :: fmt

         if (slave) call die("[rebalance:rebalance_all] slave in how_many_to_shuffle")
         if (all(pnames%exclude)) call die("[rebalance:rebalance_all] all threads excluded")

         call compute_speed
         costs_above(:) = 0.
         global_costs_above = 0.

         ! Use the following sequence of levels
         !     -1, ..., coarsest, finest, ..., base
         ! to be able to gather their costs while it is not allowed to  move them around.
         fine_processed = .false.
         curl => base%level%coarser
         if (.not. associated(curl)) curl => finest%level

         do while (associated(curl))

            if (size(curl%gp%list) > 0) then

               if (balance_cg > 0.) then
                  call curl%gp%set_sort_weight(curl%l%off, .true., balance_cg)
               else
                  call curl%gp%set_sort_weight(curl%l%off, .true.)
               endif
               global_costs_above = global_costs_above + curl%gp%w_norm
               where (speed(:) > 0.)
                  ! Increasing balance_levels towards 1. may improve the multigrid performance when few top consecutive levels have less blocks than there are MPI ranks.
                  ! An example of the problem is [8, 8] cg on 12 threads, where current algorithm will put two base-level cg on each thread unoccupied on finest level.
                  ! If this ever becomes an issue then add a postprocessing stage that will try to make some performance-neutral moves between consecutive threads.
                  load_fac(:) = max(0., 1. - (1. - max(0., balance_levels)) * costs_above(:) / (global_costs_above * speed(:)))
               elsewhere
                  load_fac(:) = 0.
               endwhere
               call compute_cumul

               if (curl%l%id >= base_level_id) then
                  fine_processed = .true.
                  p = FIRST
                  do i = lbound(curl%gp%list, 1), ubound(curl%gp%list, 1)
                     do while (curl%gp%list(i)%cweight - .5 * curl%gp%list(i)%weight > cumul(p))
                        p = p + I_ONE
                        if (p > LAST) call die("[rebalance:rebalance_all] p > LAST")
                     enddo
                     curl%gp%list(i)%dest_proc = p
                     costs_above(p) = costs_above(p) + curl%gp%list(i)%weight * curl%gp%w_norm
                  enddo
               else
                  ! Skip balancing the multigrid levels (levels below the base level)
                  ! because of unresolved problems with transfer of non-blocky levels.

                  ! This unfortunately means that these grids will not be
                  ! evacuated from excluded threads.
                  curl%gp%list(:)%dest_proc = curl%gp%list(:)%cur_proc
               endif

               if (any(curl%gp%list(:)%dest_proc == INVALID)) call die("[rebalance:rebalance_all] not all dest_proc have been set")
            endif

            curl => curl%coarser
            if (.not. associated(curl) .and. .not. fine_processed) curl => finest%level
         enddo

         rebalance_necessary = .false.
         cnt_mv(:) = 0
         curl => finest%level
         do while (associated(curl))
            cnt_gp(curl%l%id) = size(curl%gp%list)
            if (size(curl%gp%list) > 0) cnt_mv(curl%l%id) = count(curl%gp%list(:)%cur_proc /= curl%gp%list(:)%dest_proc)
            if (sum(curl%cnt_all, mask = pnames%exclude) /= 0) rebalance_necessary = .true.
            curl => curl%coarser
         enddo

         s = sum(cnt_mv)
         if (s / real(sum(cnt_gp)) > oop_thr) rebalance_necessary = .true.
         if (s > 0) then
            write(fmt, *)"(2(a,i2),a,", size(cnt_mv), "i", int(log10(real(maxval([cnt_mv, 1]))))+3, ",a,", &
                 &                      size(cnt_gp), "i", int(log10(real(maxval([cnt_gp, 1]))))+3, ",a,f6.3,a)"
            write(msg, fmt)"Rebalance: ^", lbound(cnt_gp, 1), " .. ", ubound(cnt_gp, 1), " OutOfPlace grids = [", cnt_mv, " ] / [ ", cnt_gp, &
                 " ] (", s/real(sum(cnt_gp)), " -> " // trim(merge("reshuffling)       ", "skipping reshuffle)", rebalance_necessary))
            call printinfo(msg, V_VERBOSE)
         endif

         if (.not. rebalance_necessary) s = 0.

      end function how_many_to_shuffle

      !> \brief Find the effective thread speed and calculate its cumulative distribution

      subroutine compute_speed

         use load_balance, only: balance_host, balance_thread
         use procnames,    only: pnames

         implicit none

         real :: cml

         ! Find the speed
         if (balance_host > 0. .and. pnames%speed_avail) then
            if (balance_thread) then
               speed(:) = pnames%wtime(:)
            else
               speed(:) = pnames%proc_on_node(pnames%hostindex(:))%wtime
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
         where (pnames%exclude(:)) speed(:) = 0.

         ! Normalize: speed(p) = 0.1 means that p-th MPI rank has capability to do 0.1 of sum of work of all ranks.
         cml = sum(speed(:))
         speed(:) = speed(:) / cml

         ! If balance_host is between 0. and 1. we will limit the feedback from the measured speed.
         if (balance_host > 0.) then
            where (speed(:) > 0.) speed(:) = balance_host * speed(:) + (1. - balance_host) / count(speed(:) > 0.)
         endif

      end subroutine compute_speed

      !>
      !! \brief Compute cumulative speed distribution: treat speed as a bin width.
      !! Faster threads/nodes will have wider bins, excluded threads will have 0-sized bins.
      !<

      subroutine compute_cumul

         implicit none

         integer(kind=4) :: p
         real :: cml

         cml = 0.
         do p = FIRST, LAST
            cml = cml + speed(p) * load_fac(p)
            cumul(p) = cml
         enddo

         ! Normalize
         cumul(:) = cumul(:) / cml

      end subroutine compute_cumul

   end subroutine rebalance_all

!>
!! \brief Routine for moving existing grids between processes
!!
!! Warning: The allocates of cglepa(i)%tbuf can be pretty huge (totfld * product(cg%n_) * sum_od_incoming_and_outgoing_blocks * sizeof(double)).
!! Pessimistically, the amount of allocated memory can be doubled here for a brief period of time.
!! ToDo: Make a check for available memory and implement also multi-pass variant of the exchange.
!!
!! The transfer of particles has to be implemented here alongside with the transfer of grid data.
!<

   subroutine reshuffle

      use allreduce,          only: piernik_MPI_Allreduce
      use bcast,              only: piernik_MPI_Bcast
      use cg_level_base,      only: base
      use cg_level_connected, only: cg_level_connected_t
      use cg_level_finest,    only: finest
      use cg_list,            only: cg_list_element
      use cg_list_global,     only: all_cg
      use constants,          only: ndims, LO, HI, I_ONE, xdim, ydim, zdim, pMAX, PPP_AMR
      use dataio_pub,         only: die
      use grid_cont,          only: grid_container
      use grid_container_ext, only: cg_extptrs
      use isend_irecv,        only: piernik_Isend, piernik_Irecv
      use list_of_cg_lists,   only: all_lists
      use MPIF,               only: MPI_DOUBLE_PRECISION
      use mpisetup,           only: master, proc, tag_ub
      use named_array_list,   only: qna, wna
      use ppp,                only: ppp_main
      use pppmpi,             only: req_ppp
#if defined(GRAV) && defined(NBODY)
      use domain,             only: dom
      use particle_func,      only: particle_in_area
      use particle_types,     only: particle, P_ID, P_MASS, P_POS_X, P_POS_Z, P_VEL_X, P_VEL_Z, P_ACC_X, P_ACC_Z, P_ENER, P_TFORM, P_TDYN, npf
      use particle_utils,     only: is_part_in_cg
#endif /* GRAV && NBODY */

      implicit none

      type(cg_level_connected_t), pointer :: curl
      type(cg_list_element), pointer :: cgl
      type(req_ppp) :: req
      integer :: s, n_gid, totfld
      integer(kind=4) :: i, p
      integer(kind=8), dimension(:,:), allocatable :: gptemp
      integer(kind=8), dimension(ndims, LO:HI) :: se
      logical :: found
      type(grid_container),  pointer :: cg
      type :: cglep
         type(cg_list_element), pointer :: p
         real, dimension(:,:,:,:), allocatable :: tbuf
#if defined(GRAV) && defined(NBODY)
         real, dimension(:,:), allocatable :: pbuf
#endif /* GRAV && NBODY */
      end type cglep
      type(cglep), allocatable, dimension(:) :: cglepa
      logical, parameter :: only_vital = .false. ! set to true to minimize the amount of data to be transferred, may result in improper calculation of error in maclaurin test
      !> \todo measure how much it costs in reality
      enum, bind(C)
         enumerator :: I_OFF                  ! cg offset
         enumerator :: I_N_B = I_OFF + ndims  ! cg size
         enumerator :: I_GID = I_N_B + ndims  ! cg%grid_id
         enumerator :: I_LEV                  ! cg level
         enumerator :: I_C_P                  ! cg process
         enumerator :: I_D_P                  ! cg destination process
         enumerator :: I_NP                   ! number of particles to transfer
      end enum
      character(len=*), parameter :: ISR_label = "reshuffle_Isend+Irecv", cp_label = "reshuffle_copy", gp_label = "reshuffle_gptemp"
#if defined(GRAV) && defined(NBODY)
      logical :: in, phy, out, fin, indomain
      type(particle), pointer :: part
      character(len=*), parameter :: cpp_label = "reshuffle_copy_part"
#endif /* GRAV && NBODY */

      ! Count the number of fields on any of the base level cg.
      ! Assume that base level and above have the same set of fields. Levels coarser than base may have different set of fields.
      totfld = 0
      cgl => base%level%first
      if (associated(cgl)) then
         do p = lbound(wna%lst, dim=1, kind=4), ubound(wna%lst, dim=1, kind=4)
            if ((.not. only_vital .or. wna%lst(p)%vital) .and. associated(cgl%cg%w(p)%arr)) totfld = totfld + wna%get_dim4(p)
         enddo
         do p = lbound(qna%lst, dim=1, kind=4), ubound(qna%lst, dim=1, kind=4)
            if ((.not. only_vital .or. qna%lst(p)%vital) .and. associated(cgl%cg%q(p)%arr)) totfld = totfld + 1
         enddo
      endif
      call piernik_MPI_Allreduce(totfld, pMAX)

      s = 0  ! The number of cg to send
      if (master) then
         curl => finest%level
         do while (associated(curl))
            s = s + count(curl%gp%list(:)%cur_proc /= curl%gp%list(:)%dest_proc)
            curl => curl%coarser
         enddo
      endif
      call piernik_MPI_Bcast(s)

      call ppp_main%start(gp_label, PPP_AMR)
      allocate(gptemp(I_OFF:I_NP, s), cglepa(s))
      if (master) then
         p = 0
         curl => finest%level
         do while (associated(curl))
            do i = lbound(curl%gp%list, dim=1, kind=4), ubound(curl%gp%list, dim=1, kind=4)
               if (curl%gp%list(i)%cur_proc /= curl%gp%list(i)%dest_proc) then
                  p = p + I_ONE
                  gptemp(:, p) = [ curl%gp%list(i)%off, int( [ curl%gp%list(i)%n_b, curl%gp%list(i)%cur_gid, curl%l%id, curl%gp%list(i)%cur_proc, curl%gp%list(i)%dest_proc, curl%gp%list(i)%n_part ], kind=8) ]
               endif
            enddo
            curl => curl%coarser
         enddo
      endif
      call piernik_MPI_Bcast(gptemp)
      call ppp_main%stop(gp_label, PPP_AMR)

      ! Irecv & Isend
      call ppp_main%start(ISR_label, PPP_AMR)
      call req%init(owncomm = .true., label = "reshuffle")
      curl => finest%level
      do while (associated(curl))

         if (ubound(gptemp, dim=2, kind=4) > tag_ub) call die("[rebalance:reshuffle] this MPI implementation has too low MPI_TAG_UB attribute")
         do i = lbound(gptemp, dim=2, kind=4), ubound(gptemp, dim=2, kind=4)
            associate (ip => i + size(gptemp, dim=2, kind=4))
            if (gptemp(I_LEV, i) == curl%l%id) then

               cglepa(i)%p => null()
               if (gptemp(I_C_P, i) == gptemp(I_D_P, i)) call die("[rebalance:reshuffle] can not send to self")
               if (gptemp(I_C_P, i) == proc) then ! send
                  found = .false.
                  cgl => curl%first
                  do while (associated(cgl))
                     if (cgl%cg%grid_id == gptemp(I_GID, i)) then
                        found = .true.
                        cglepa(i)%p => cgl
                        allocate(cglepa(i)%tbuf(totfld, cgl%cg%n_(xdim), cgl%cg%n_(ydim), cgl%cg%n_(zdim)))
                        ! We communicate blocks with guardcells because current implementation of magnetic field evolves external guardcells in a way that makes it impossible
                        ! to reconstruct them from scratch. This results in much larger messages to communicate, but we won't need to call guardcell exchange afterwards.
                        s = lbound(cglepa(i)%tbuf, dim=1)
                        do p = lbound(wna%lst, dim=1, kind=4), ubound(wna%lst, dim=1, kind=4)
                           if ((.not. only_vital .or. wna%lst(p)%vital) .and. associated(cgl%cg%w(p)%arr)) then ! not associated for multigrid coarse levels
                              cglepa(i)%tbuf(s:s+wna%get_dim4(p)-1, :, :, :) = cgl%cg%w(p)%arr(:, :, :, :)
                              s = s + wna%get_dim4(p)
                           endif
                        enddo
                        do p = lbound(qna%lst, dim=1, kind=4), ubound(qna%lst, dim=1, kind=4)
                           if ((.not. only_vital .or. qna%lst(p)%vital) .and. associated(cgl%cg%q(p)%arr)) then
                              cglepa(i)%tbuf(s, :, :, :) = cgl%cg%q(p)%arr(:, :, :)
                              s = s + 1
                           endif
                        enddo

#if defined(GRAV) && defined(NBODY)
                        ! Make copy of the particles here
                        allocate(cglepa(i)%pbuf(npf, gptemp(I_NP, i)))

                        part => cgl%cg%pset%first
                        p = 1
                        do while (associated(part))
                           cglepa(i)%pbuf(P_ID, p)            = part%pdata%pid
                           cglepa(i)%pbuf(P_MASS, p)          = part%pdata%mass
                           cglepa(i)%pbuf(P_POS_X:P_POS_Z, p) = part%pdata%pos
                           cglepa(i)%pbuf(P_VEL_X:P_VEL_Z, p) = part%pdata%vel
                           cglepa(i)%pbuf(P_ACC_X:P_ACC_Z, p) = part%pdata%acc
                           cglepa(i)%pbuf(P_ENER, p)          = part%pdata%energy
                           cglepa(i)%pbuf(P_TFORM, p)         = part%pdata%tform
                           cglepa(i)%pbuf(P_TDYN, p)          = part%pdata%tdyn

                           p = p + I_ONE
                           part => part%nxt
                        enddo
#endif /* GRAV && NBODY */

                        exit
                     endif
                     cgl => cgl%nxt
                  enddo
                  if (.not. found) call die("[rebalance:reshuffle] Grid id not found")
                  ! explicit buf(lbound(buf, ...), ...) needed to prevent valgrind complains on "Invalid read of size 8", at least with gfortran 12.3
                  call piernik_Isend(cglepa(i)%tbuf(lbound(cglepa(i)%tbuf, 1):, lbound(cglepa(i)%tbuf, 2):, lbound(cglepa(i)%tbuf, 3):, lbound(cglepa(i)%tbuf, 4):), size(cglepa(i)%tbuf, kind=4), MPI_DOUBLE_PRECISION, int(gptemp(I_D_P, i), kind=4), i,  req)

#if defined(GRAV) && defined(NBODY)
                  ! Isend for particles
                  call piernik_Isend(cglepa(i)%pbuf(lbound(cglepa(i)%pbuf, 1):, lbound(cglepa(i)%pbuf, 2):), size(cglepa(i)%pbuf, kind=4), MPI_DOUBLE_PRECISION, int(gptemp(I_D_P, i), kind=4), ip, req)
#endif /* GRAV && NBODY */

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
                  ! explicit buf(lbound(buf, ...), ...) needed to prevent valgrind complains on "Invalid read of size 8", at least with gfortran 12.3
                  call piernik_Irecv(cglepa(i)%tbuf(lbound(cglepa(i)%tbuf, 1):, lbound(cglepa(i)%tbuf, 2):, lbound(cglepa(i)%tbuf, 3):, lbound(cglepa(i)%tbuf, 4):), size(cglepa(i)%tbuf, kind=4), MPI_DOUBLE_PRECISION, int(gptemp(I_C_P, i), kind=4), i, req)

#if defined(GRAV) && defined(NBODY)
                  ! Irecv for particles
                  allocate(cglepa(i)%pbuf(npf, gptemp(I_NP, i)))
                  call piernik_Irecv(cglepa(i)%pbuf(lbound(cglepa(i)%pbuf, 1):, lbound(cglepa(i)%pbuf, 2):), size(cglepa(i)%pbuf, kind=4), MPI_DOUBLE_PRECISION, int(gptemp(I_C_P, i), kind=4), ip, req)
#endif /* GRAV && NBODY */

               endif
            endif
            end associate
         enddo

         curl => curl%coarser
      enddo
      call ppp_main%stop(ISR_label, PPP_AMR)

      call req%waitall("reshuffle", PPP_AMR)

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

#if defined(GRAV) && defined(NBODY)
                  deallocate(cglepa(i)%pbuf)
#endif /* GRAV && NBODY */

               endif
               if (gptemp(I_D_P, i) == proc) then ! copy received
                  s = lbound(cglepa(i)%tbuf, dim=1)
                  do p = lbound(wna%lst, dim=1, kind=4), ubound(wna%lst, dim=1, kind=4)
                     if ((.not. only_vital .or. wna%lst(p)%vital) .and. associated(cgl%cg%w(p)%arr)) then
                        cgl%cg%w(p)%arr(:, :, :, :) = cglepa(i)%tbuf(s:s+wna%get_dim4(p)-1, :, :, :)
                        s = s + wna%get_dim4(p)
                     endif
                  enddo
                  do p = lbound(qna%lst, dim=1, kind=4), ubound(qna%lst, dim=1, kind=4)
                     if ((.not. only_vital .or. qna%lst(p)%vital) .and. associated(cgl%cg%q(p)%arr)) then
                        cgl%cg%q(p)%arr(:, :, :) = cglepa(i)%tbuf(s, :, :, :)
                        s = s + 1
                     endif
                  enddo
                  deallocate(cglepa(i)%tbuf)

#if defined(GRAV) && defined(NBODY)
                  ! assign imported particles
                  call ppp_main%start(cpp_label, PPP_AMR)
                  do p = lbound(cglepa(i)%pbuf, 2, kind=4), ubound(cglepa(i)%pbuf, 2, kind=4)
                     indomain = particle_in_area(cglepa(i)%pbuf(P_POS_X:P_POS_Z, p), dom%edge)
                     call is_part_in_cg(cgl%cg, cglepa(i)%pbuf(P_POS_X:P_POS_Z, p), indomain, in, phy, out, fin)
                     call cgl%cg%pset%add(nint(cglepa(i)%pbuf(P_ID, p), kind=4), cglepa(i)%pbuf(P_MASS, p), &
                          cglepa(i)%pbuf(P_POS_X:P_POS_Z, p), cglepa(i)%pbuf(P_VEL_X:P_VEL_Z, p), &
                          cglepa(i)%pbuf(P_ACC_X:P_ACC_Z, p), cglepa(i)%pbuf(P_ENER, p), &
                          in, phy, out, fin, &
                          cglepa(i)%pbuf(P_TFORM, p), cglepa(i)%pbuf(P_TDYN, p))
                  enddo
                  deallocate(cglepa(i)%pbuf)
                  call ppp_main%stop(cpp_label, PPP_AMR)
#endif /* GRAV && NBODY */

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

end module rebalance
