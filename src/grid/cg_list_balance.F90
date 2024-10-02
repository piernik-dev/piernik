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

!> \brief A module with an abstract type created to take out some load-balance related code from cg_level (new grids)

module cg_list_balance

   use cg_list_bnd,      only: cg_list_bnd_t
   use constants,        only: ndims, I_ONE
   use dot,              only: dot_t
   use level_essentials, only: level_t
   use patch_list,       only: patch_list_t
   use sort_piece_list,  only: grid_piece_list

   implicit none

   private
   public :: cg_list_balance_t, I_OFF, I_N_B

   enum, bind(C)
      enumerator :: I_OFF
      enumerator :: I_N_B = I_OFF + ndims
      enumerator :: I_END = I_N_B + ndims - I_ONE
   end enum

   !> An abstract type created to take out some load-balance related code from cg_level (new grids)

   type, extends(cg_list_bnd_t), abstract :: cg_list_balance_t
      type(patch_list_t)                :: plist             !< list of patches that exist on the current level
      type(dot_t)                       :: dot               !< depiction of topology
      logical                           :: recently_changed  !< .true. when anything was added to or deleted from this level
      class(level_t), pointer           :: l                 !< single place to store off, n_d and id
      type(grid_piece_list)             :: gp                !< SFC-sortable structure for collecting cgs for rebalance
      integer(kind=4), allocatable, dimension(:) :: cnt_all  !< block count for all threads (used in rebalance, only on master)
   contains
      procedure          :: balance_new          !< Routine selector for moving proposed grids between processes
      procedure, private :: balance_fill_lowest  !< Routine for moving proposed grids between processes: add load to lightly-loaded processes
      procedure, private :: patches_to_list      !< Collect local proposed patches into an array on the master process
      procedure, private :: distribute_patches   !< Send balanced set patches from master to slaves and re-register them
      procedure, private :: add_patch_one_piece  !< Add a patch with only one grid piece
      procedure          :: sort_SFC             !< Sort list according to SFC id
   end type cg_list_balance_t

contains

!>
!! \brief Routine selector for moving proposed grids between processes
!!
!! There are several strategies that can be implemented:
!! * Local refinements go to local process. It is very simple, but for most
!!   simulations will build up load imbalance. Suitable for tests and global
!!   refinement. Removed.
!! * Local refinements can be assigned to remote processes, existing blocks
!!   stays in place. Should keep good load balance, but the amount of inter-process
!!   internal boundaries may grow significantly with time. Suitable for minor
!!   refinement updates and base level decomposition. This is the current
!!   implementation in balance_fill_lowest.
!! * All blocks (existing and new) have recalculated assignment and can be
!!   migrated to other processes. Most advanced. Should be used after reading
!!   restart data. ToDo.
!<

   subroutine balance_new(this)

      implicit none

      class(cg_list_balance_t), intent(inout) :: this

      call this%balance_fill_lowest ! never rebalances

   end subroutine balance_new

!>
!! \brief Routine for moving proposed grids between processes: add load to lightly-loaded processes
!!
!! \details Starts with list of already allocated blocks and list of patches which are not yet turned into blocks on a given level.
!! Move the patches between processes to maintain best possible work balance.
!!
!! First, all planned patches are gathered in an array on the master process, and deallocated locally.
!! Then, the patches are sorted according to Space-Filling Curve index and distributed among the processes.
!! The processes with least workload will get more patches.
!! After the distribution most processes should have roughly equal number of patches (+/- 1) with the possible exception
!! of few processes that were initially heavily loaded.
!! The load is counted only on current level, so the load imbalance may add up across several levels.
!!
!! As it is hard to predict future cg cost, we better refrain from sophisticated
!! load balancing here and wait until some performance data is collected.
!!
!! Note that this routine is not intended for moving existing blocks between processes.
!! A separate routine, called from cg_leaves::update will do that task when allowed and found worth the effort.
!!
!! Current implementation does all the work on master process which might be quite antiparallel.
!!
!! This is truly parallel-sorting problem. Note that the set of grid pieces has the following property:
!! * there is sorted or nearly-sorted list of existing grid pieces on each process, that means for most pieces on process p maximum id on process p-1 is less than own id
!!   and for process p+1 similarly
!! * there is chaotic (in practice not so much) set of grid pieces to be created
!! We may then sort iteratively:
!! * do long-range moves of chaotic pieces, based on distribution estimate
!! * iterate with short-range (+/-1 or at most +/-2 in process number) moves of all pieces until everything is sorted well enough
!<

   subroutine balance_fill_lowest(this)

      use allreduce,  only: piernik_MPI_Allreduce
      use constants,  only: pSUM, I_ZERO, I_ONE
      use dataio_pub, only: die
      use MPIF,       only: MPI_INTEGER, MPI_COMM_WORLD
      use MPIFUN,     only: MPI_Gather
      use mpisetup,   only: master, FIRST, LAST, err_mpi
      use procnames,  only: pnames

      use sort_piece_list, only: grid_piece_list

      implicit none

      class(cg_list_balance_t), intent(inout) :: this

      type(grid_piece_list) :: gp
      integer :: i
      integer(kind=4), dimension(FIRST:LAST+1) :: from
      integer(kind=4), dimension(FIRST:LAST) :: cnt_all
      integer(kind=4) :: ls, p, s

      ! count how many patches were requested on each process
      s = int(this%plist%p_count(), kind=4)
      ls = int(s, kind=4)
      call piernik_MPI_Allreduce(s, pSUM) !> \warning overkill: MPI_reduce is enough here
      if (s==0) return ! nihil novi

      if (master) call gp%init(s)
      call this%patches_to_list(gp, ls)

      call MPI_Gather(this%cnt, I_ONE, MPI_INTEGER, cnt_all, I_ONE, MPI_INTEGER, FIRST, MPI_COMM_WORLD, err_mpi)

      if (master) then !> \warning Antiparallel

         if (all(pnames%exclude)) call die("[cg_list_balance:balance_fill_lowest] all threads excluded")

         ! apply unique numbers to the grids and sort the list
         call gp%set_sort_weight(this%l%off, .false.)

         !> compute destination process corrected for current imbalance of existing grids as much as possible
         !> \todo replace counting of blocks with counting of weights - it will be required for merged blocks

         ! target number of cg per thread (ideal, level-balanced)
         s = int((size(gp%list) + sum(cnt_all))/real(count(.not. pnames%exclude)), kind=4)

         ! first estimate how to fill lowest-occupied threads
         i = (size(gp%list) + sum(cnt_all, mask=(cnt_all <= s .and. .not. pnames%exclude))) / &
              &               count(cnt_all <= s .and. .not. pnames%exclude)
         from(:) = [ lbound(gp%list, dim=1, kind=4), merge(I_ZERO, int(max(0, i - cnt_all(:)), kind=4), pnames%exclude(:)) ]
         i = size(gp%list) - sum(from(FIRST+1:LAST+1))

         if (i < 0) then
            p = FIRST
            do while (i < 0)
               if (from(p+1) > 0) then
                  from(p+1) = from(p+1) - I_ONE
                  i = i + 1
               endif
               p = p + I_ONE
               if (p > LAST) p = FIRST
            enddo
         else if (i > 0) then  ! this is the dominating case
            p = LAST
            do while (i > 0)
               ! here we tend to overload last processes on each level independently
               ! ToDo involve all_cg%cnt and try to improve global balancing
               if (.not. pnames%exclude(p)) then
                  from(p+1) = from(p+1) + I_ONE
                  i = i - 1
               endif
               p = p - I_ONE
               if (p < FIRST) p = LAST
            enddo
         endif

         if (size(gp%list) - sum(from(FIRST+1:LAST+1)) /= 0) call die("[cg_list_balance:balance_fill_lowest] i /= 0")
         do p = FIRST, LAST
            from(p+1) = from(p+1) + from(p)
         enddo
      endif

      call this%distribute_patches(gp, from)

      if (master) call gp%cleanup

   end subroutine balance_fill_lowest

!> \brief collect local proposed patches on a given level into an array on the master process

   subroutine patches_to_list(this, gp, ls)

      use constants,       only: ndims, INVALID, I_ONE
      use MPIF,            only: MPI_INTEGER, MPI_INTEGER8, MPI_STATUS_IGNORE, MPI_COMM_WORLD, MPI_Wait
      use MPIFUN,          only: MPI_Isend, MPI_Send, MPI_Recv
      use mpisetup,        only: master, slave, FIRST, LAST, req, err_mpi, inflate_req
      use sort_piece_list, only: grid_piece_list

      implicit none

      class(cg_list_balance_t), intent(inout) :: this !< object invoking type bound procedure
      type(grid_piece_list),    intent(inout) :: gp   !< list of grid pieces to be filled
      integer(kind=4),          intent(in)    :: ls   !< local number of patches

      integer(kind=8), dimension(:,:), allocatable :: gptemp
      integer :: i
      integer(kind=4) :: p, ss, rls
      integer(kind=4), parameter :: nreq = 1
      integer(kind=4), parameter :: tag_ls = 1, tag_gpt = tag_ls + 1

      call inflate_req(nreq)

      ! copy the patches data to a temporary array to be sent to the master
      if (slave) call MPI_Isend(ls, I_ONE, MPI_INTEGER, FIRST, tag_ls, MPI_COMM_WORLD, req(nreq), err_mpi)

      allocate(gptemp(I_OFF:I_END, ls))
      call this%plist%p2a(gptemp)

      if (master) then !> \warning Antiparallel

         ! put all the patches (own and obtained from slaves) on a list gp%list
         do ss = 1, ls
            call gp%list(ss)%set_gp(gptemp(I_OFF:I_OFF+ndims-1, ss), int(gptemp(I_N_B:I_N_B+ndims-1, ss), kind=4), INVALID, FIRST)
         enddo
         i = ls
         deallocate(gptemp)
         do p = FIRST + 1, LAST
            call MPI_Recv(rls, I_ONE, MPI_INTEGER, p, tag_ls, MPI_COMM_WORLD, MPI_STATUS_IGNORE, err_mpi)
            if (rls > 0) then
               allocate(gptemp(I_OFF:I_END, rls))
               call MPI_Recv(gptemp, size(gptemp, kind=4), MPI_INTEGER8, p, tag_gpt, MPI_COMM_WORLD, MPI_STATUS_IGNORE, err_mpi)
               do ss = 1, rls
                  call gp%list(i+ss)%set_gp(gptemp(I_OFF:I_OFF+ndims-1, ss), int(gptemp(I_N_B:I_N_B+ndims-1, ss), kind=4), INVALID, p)
               enddo
               i = i + rls
               deallocate(gptemp)
            endif
         enddo
      else

         ! send patches to master
         call MPI_Wait(req(nreq), MPI_STATUS_IGNORE, err_mpi)
         if (ls > 0) call MPI_Send(gptemp, size(gptemp, kind=4), MPI_INTEGER8, FIRST, tag_gpt, MPI_COMM_WORLD, err_mpi)
         deallocate(gptemp)

      endif

   end subroutine patches_to_list

!>
!! \brief Send balanced set patches from master to slaves and re-register them
!!
!! \todo Try to utilize gp%list(:)%dest_proc and not rely on from(:)
!<

   subroutine distribute_patches(this, gp, from)

      use constants,       only: ndims, I_ONE
      use MPIF,            only: MPI_INTEGER, MPI_INTEGER8, MPI_STATUS_IGNORE, MPI_COMM_WORLD
      use MPIFUN,          only: MPI_Send, MPI_Recv
      use mpisetup,        only: master, FIRST, LAST, err_mpi
      use sort_piece_list, only: grid_piece_list

      implicit none

      class(cg_list_balance_t),                 intent(inout) :: this !< object invoking type bound procedure
      type(grid_piece_list),                    intent(in)    :: gp   !< list of grid pieces to be filled
      integer(kind=4), dimension(FIRST:LAST+1), intent(in)    :: from !< indices that mark ranges in gp to be sent to processes

      integer(kind=8), dimension(:,:), allocatable :: gptemp
      integer(kind=4) :: ls, p, s
      integer(kind=4), parameter :: tag_lsR = 3, tag_gptR = tag_lsR+1

      call this%plist%p_deallocate
      if (master) then
         do p = from(FIRST), from(FIRST+1) - I_ONE
            call this%add_patch_one_piece(int(gp%list(p)%n_b, kind=8), gp%list(p)%off)
         enddo

         ! distribute proposed grids according to limits computed above
         do p = FIRST + I_ONE, LAST
            ls = int(from(p+1) - from(p), kind=4)
            call MPI_Send(ls, I_ONE, MPI_INTEGER, p, tag_lsR, MPI_COMM_WORLD, err_mpi)
            if (ls>0) then
               allocate(gptemp(I_OFF:I_END,from(p):from(p+1)-1))
               do s = lbound(gptemp, dim=2, kind=4), ubound(gptemp, dim=2, kind=4)
                  gptemp(:, s) = [ gp%list(s)%off, int(gp%list(s)%n_b, kind=8) ]
               enddo
               call MPI_Send(gptemp, size(gptemp, kind=4), MPI_INTEGER8, p, tag_gptR, MPI_COMM_WORLD, err_mpi)
               deallocate(gptemp)
            endif
         enddo
      else
         ! receive new, perhaps more balanced patches
         call MPI_Recv(ls, I_ONE, MPI_INTEGER, FIRST, tag_lsR, MPI_COMM_WORLD, MPI_STATUS_IGNORE, err_mpi)
         if (ls>0) then
            allocate(gptemp(I_OFF:I_END,ls))
            call MPI_Recv(gptemp, size(gptemp, kind=4), MPI_INTEGER8, FIRST, tag_gptR, MPI_COMM_WORLD, MPI_STATUS_IGNORE, err_mpi)
            do s = lbound(gptemp, dim=2, kind=4), ubound(gptemp, dim=2, kind=4)
               call this%add_patch_one_piece(gptemp(I_N_B:I_N_B+ndims-1, s), gptemp(I_OFF:I_OFF+ndims-1, s))
            enddo
            deallocate(gptemp)
         endif
      endif

   end subroutine distribute_patches

!> \brief Add a patch with only one grid piece

   subroutine add_patch_one_piece(this, n_d, off)

      use constants, only: ndims

      implicit none

      class(cg_list_balance_t), target,  intent(inout) :: this  !< current level
      integer(kind=8), dimension(ndims), intent(in)    :: n_d   !< number of grid cells
      integer(kind=8), dimension(ndims), intent(in)    :: off   !< offset (with respect to the base level, counted on own level)

      this%recently_changed = .true. ! assume that the new patches will change this level
      call this%plist%expand
      call this%plist%patches(ubound(this%plist%patches(:), dim=1))%one_piece_patch(n_d(:), off(:))

   end subroutine add_patch_one_piece

!> \brief Sort list according to SFC id

   subroutine sort_SFC(this)

      use cg_list,      only: cg_list_element
      use sort_cg_list, only: sort_cg_list_t

      implicit none

      class(cg_list_balance_t), intent(inout) :: this !< object invoking type-bound procedure

      type(cg_list_element), pointer :: cgl
      integer :: s
      type(sort_cg_list_t) :: l

      if (this%cnt <= 0) return ! nothing to sort or renumber

      ! Create auxiliary array of pointers and sort them
      call l%init(this%cnt)
      s = lbound(l%list, dim=1)
      cgl => this%first
      do while (associated(cgl))
         l%list(s)%cgl => cgl
         s = s + 1
         cgl => cgl%nxt
      enddo
      call l%sort

      ! relink all cg_list_elements and refresh cg%grid_id
      do s = lbound(l%list, dim=1), ubound(l%list, dim=1)
         if (s > lbound(l%list, dim=1)) then
            l%list(s)%cgl%prv => l%list(s-1)%cgl
         else
            l%list(s)%cgl%prv => null()
         endif
         if (s < ubound(l%list, dim=1)) then
            l%list(s)%cgl%nxt => l%list(s+1)%cgl
         else
            l%list(s)%cgl%nxt => null()
         endif
         l%list(s)%cgl%cg%grid_id = s
      enddo
      this%first => l%list(lbound(l%list, dim=1))%cgl
      this%last  => l%list(ubound(l%list, dim=1))%cgl

      call l%cleanup

   end subroutine sort_SFC

end module cg_list_balance
