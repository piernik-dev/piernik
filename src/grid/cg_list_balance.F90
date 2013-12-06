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

!> \brief A module with an abstract type created to take out some load-balance related code from cg_level (new grids)

module cg_list_balance

   use cg_list_bnd, only: cg_list_bnd_T
   use constants,   only: ndims, I_ONE
   use dot,         only: dot_T
   use patch_list,  only: patch_list_T

   implicit none

   private
   public :: cg_list_balance_T, I_OFF, I_N_B

   enum, bind(C)
      enumerator :: I_OFF
      enumerator :: I_N_B = I_OFF + ndims
      enumerator :: I_END = I_N_B + ndims - I_ONE
   end enum

   !> An abstract type created to take out some load-balance related code from cg_level (new grids)

   type, extends(cg_list_bnd_T), abstract :: cg_list_balance_T
      type(patch_list_T)                :: plist            !< list of patches that exist on the current level
      type(dot_T)                       :: dot              !< depiction of topology
      integer(kind=8), dimension(ndims) :: off              !< offset of the level
      logical                           :: recently_changed !< .true. when anything was added to or deleted from this level
   contains
      procedure          :: balance_new          !< Routine selector for moving proposed grids between processes
      procedure, private :: balance_strict_SFC   !< Routine for moving proposed grids between processes: keep strict SFC ordering
      procedure, private :: balance_fill_lowest  !< Routine for moving proposed grids between processes: add load to lightly-loaded processes
      procedure, private :: patches_to_list      !< Collect local proposed patches into an array on the master process
      procedure, private :: distribute_patches   !< Send balanced set patches from master to slaves and re-register them
      procedure, private :: add_patch_one_piece  !< Add a patch with only one grid piece
      procedure          :: sort_SFC             !< Sort list according to SFC id
   end type cg_list_balance_T

contains

!>
!! \brief Routine selector for moving proposed grids between processes
!!
!! There are several strategies that can be implemented:
!! * Local refinements go to local process. It is very simple, but for most simulations will build up load imbalance. Suitable for tests and global refinement.
!! * Local refinements can be assigned to remote processes, existing blocks stays in place. Should keep good load balance, but the amount of inter-process
!!   internal boundaries may grow significantly with time. Suitable for minor refinement updates and base level decomposition. This is the current implementation of balance_fill_lowest.
!! * All blocks (existing and new) have recalculated assignment and can be migrated to other processes. Most advanced. Should be used after reading restart data.
!<

   subroutine balance_new(this, prevent_rebalancing)

      use refinement, only: strict_SFC_ordering

      implicit none

      class(cg_list_balance_T), intent(inout) :: this
      logical, optional,        intent(in)    :: prevent_rebalancing !< if present and .true. then do not allow rebalancing during addition of new grids

      ! The only available strategy ATM
      if (strict_SFC_ordering) then
         call this%balance_strict_SFC(prevent_rebalancing)
      else
         call this%balance_fill_lowest ! never rebalances
      endif

   end subroutine balance_new

!>
!! \brief Routine for moving proposed grids between processes: keep strict SFC ordering
!!
!! \details Starts with list of already allocated blocks and list of patches which are not yet turned into blocks on a given level.
!! Assume that the existing blocks are distributed according to Space-Filling curve ordering: maximum id for process p is always lower than minimum id for process p+1.
!! Add new grids in a way that keeps the strict SFC ordering property even if it may introduce imbalance
!!
!! \todo do a global rebalance if it is allowed and worth the effort
!<

   subroutine balance_strict_SFC(this, prevent_rebalancing)

      use constants,       only: pSUM, LO, HI, I_ONE
      use dataio_pub,      only: warn
      use mpisetup,        only: piernik_MPI_Allreduce, master, FIRST, LAST, nproc
      use sort_piece_list, only: grid_piece_list

      implicit none

      class(cg_list_balance_T), intent(inout) :: this
      logical, optional,        intent(in)    :: prevent_rebalancing !< if present and .true. then do not allow rebalancing during addition of new grids

      logical :: rebalance
      type(grid_piece_list) :: gp
      integer(kind=4) :: ls, s, i, p
      integer(kind=4), dimension(FIRST:LAST+1) :: from

      rebalance = .true.
      if (present(prevent_rebalancing)) rebalance = .not. prevent_rebalancing

      if (.not. this%dot%check_SFC(this%off)) then
!!$         if (.not. rebalance) call die("[cg_list_balance:balance_strict_SFC] Cannot rebalence messy grid distribution.")
!!$         ! call reshuffle
!!$         call die("[cg_list_balance:balance_strict_SFC] reshuffling not implemented.")
         if (master) call warn("[cg_list_balance:balance_strict_SFC] non-SFC ordering!") ! May happen after resizing domain on the left sides
      endif

      ! gather patches id
      s = int(this%plist%p_count(), kind=4)
      ls = int(s, kind=4)
      call piernik_MPI_Allreduce(s, pSUM) !> \warning overkill: MPI_reduce is enough here
      if (s==0) return ! nihil novi

      if (master) call gp%init(s)
      call this%patches_to_list(gp, ls)

      ! if (rebalance) gather existing grids id

      ! sort id
      if (master) then !> \warning Antiparallel
         ! apply unique numbers to the grids and sort the list
         call gp%set_id(this%off)
         call gp%sort

         ! calculate patch distribution
         if (all(this%dot%SFC_id_range(:, HI) < this%dot%SFC_id_range(:, LO))) then !special case: empty level, huge values in this%dot%SFC_id_range(:,:)
            do p = FIRST, LAST
               from(p) = int(lbound(gp%list, dim=1) + (p*size(gp%list))/nproc, kind=4)
            enddo
            from(LAST+1) = ubound(gp%list, dim=1, kind=4) + I_ONE
!!$            do p = FIRST, LAST
!!$               gp%list(from(p):from(p+1)-1)%dest_proc = p
!!$            enddo
         else
            ! just keep SFC ordering, no attempts to balance things here
            !> \todo OPT try to do as much balance as possible
            p = FIRST
            do i = lbound(gp%list, dim=1, kind=4), ubound(gp%list, dim=1, kind=4)
               do while (this%dot%SFC_id_range(p, HI) < this%dot%SFC_id_range(p, LO)) ! skip processes with no grids for the sake of simplicity
                  p = p + I_ONE
                  if (p > LAST) exit
               enddo
               if (p > LAST) p = LAST
               do while (gp%list(i)%id > this%dot%SFC_id_range(p, HI) .and. p < LAST)
                  p = p + I_ONE
                  if (p > LAST) exit
               enddo
               gp%list(i)%dest_proc = p
            enddo

            p = FIRST
            from(FIRST)  = lbound(gp%list, dim=1, kind=4)
            from(FIRST+1:LAST+1) = ubound(gp%list, dim=1, kind=4) + I_ONE
            do i = lbound(gp%list, dim=1, kind=4), ubound(gp%list, dim=1, kind=4)
               if (gp%list(i)%dest_proc /= p) then
                  from(p+1:gp%list(i)%dest_proc) = i
                  p = gp%list(i)%dest_proc
               endif
            enddo
         endif
      endif
      ! if (rebalance) call reshuffle(distribution)

      ! send to slaves
      call this%distribute_patches(gp, from)

      if (master) call gp%cleanup

   end subroutine balance_strict_SFC

!>
!! \brief Routine for moving proposed grids between processes: add load to lightly-loaded processes
!!
!! \details Starts with list of already allocated blocks and list of patches which are not yet turned into blocks on a given level.
!! Move the patches between processes to maintain best possible work balance.
!!
!! First, all planned patches are gathered in an array on the master process, and deallocated locally.
!! Then, the patches are sorted according to Space-Filling Curve index and distributed among the prosesses.
!! The processes with least workload will get more patches.
!! After the distribution most processes should have roughly equal number of patches (+/- 1) with the possible exception
!! of few processes that were initially heavily loaded.
!! The load is counted only on current level, so the load imbalance may add up across several levels.
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

      use constants,       only: pSUM, I_ONE
      use dataio_pub,      only: die
      use mpi,             only: MPI_INTEGER
      use mpisetup,        only: piernik_MPI_Allreduce, master, FIRST, LAST, comm, mpi_err, nproc
      use sort_piece_list, only: grid_piece_list

      implicit none

      class(cg_list_balance_T), intent(inout) :: this

      type(grid_piece_list) :: gp
      integer :: i
      integer(kind=4), dimension(FIRST:LAST+1) :: from
      integer(kind=4), dimension(FIRST:LAST) :: cnt_existing
      integer(kind=4) :: ls, p, s

      ! count how many patches were requested on each process
      s = int(this%plist%p_count(), kind=4)
      ls = int(s, kind=4)
      call piernik_MPI_Allreduce(s, pSUM) !> \warning overkill: MPI_reduce is enough here
      if (s==0) return ! nihil novi

      if (master) call gp%init(s)
      call this%patches_to_list(gp, ls)

      call MPI_Gather(this%cnt, I_ONE, MPI_INTEGER, cnt_existing, I_ONE, MPI_INTEGER, FIRST, comm, mpi_err)

      if (master) then !> \warning Antiparallel
         ! apply unique numbers to the grids and sort the list
         call gp%set_id(this%off)
         call gp%sort

         ! measure their weight (unused yet)
         ! call gp%set_weights

         !> compute destination process corrected for current imbalance of existing grids as much as possible
         !> \todo replace counting of blocks with counting of weights - it will be required for merged blocks
         s = int((size(gp%list) + sum(cnt_existing))/real(nproc), kind=4)
         i = (size(gp%list) + sum(cnt_existing, mask=(cnt_existing <= s)))/(nproc - count(cnt_existing > s))
         from(FIRST) = lbound(gp%list, dim=1, kind=4)
         do p = FIRST, LAST
            from(p+1) = int(max(0, i - cnt_existing(p)), kind=4)
         enddo
         i = size(gp%list) - sum(from(FIRST+1:LAST+1))
         p = LAST
         do while (p >= FIRST .and. i /= 0)
            if (i<0) then
               if (from(p+1)>0) then
                  from(p+1) = from(p+1) - I_ONE
                  i = i + 1
               endif
            else if (i>0) then
               !> \deprecated this approach may result in building a small imbalance in favour of process with low id.
               if (cnt_existing(p) <= s) then
                  from(p+1) = from(p+1) + I_ONE
                  i = i - 1
               endif
            endif
            p = p - I_ONE
         enddo
         i = size(gp%list) - sum(from(FIRST+1:LAST+1))
         if (i /= 0) call die("[cg_list_balance:balance_fill_lowest] i /= 0")
         do p = FIRST, LAST
            from(p+1) = from(p+1) + from(p)
         enddo
      endif

      call this%distribute_patches(gp, from)

      if (master) call gp%cleanup

!!$      allocate(area(1:lmax))
!!$      area = 0
!!$      do i = lbound(cg_res(:), dim=1), ubound(cg_res(:), dim=1)
!!$         if (cg_res(i)%level >= 1) area(cg_res(i)%level) = area(cg_res(i)%level) + product(cg_res(i)%n_b)
!!$      enddo
!!$      deallocate(area)

   end subroutine balance_fill_lowest

!> \brief collect local proposed patches on a given level into an array on the master process

   subroutine patches_to_list(this, gp, ls)

      use constants,       only: ndims, INVALID, I_ONE
      use mpi,             only: MPI_INTEGER, MPI_INTEGER8, MPI_STATUS_IGNORE
      use mpisetup,        only: master, slave, FIRST, LAST, comm, req, mpi_err, status, inflate_req
      use sort_piece_list, only: grid_piece_list

      implicit none

      class(cg_list_balance_T), intent(inout) :: this !< object invoking type bound procedure
      type(grid_piece_list),    intent(inout) :: gp   !< list of grid pieces to be filled
      integer(kind=4),          intent(in)    :: ls   !< local number of patches

      integer(kind=8), dimension(:,:), allocatable :: gptemp
      integer :: i
      integer(kind=4) :: p, ss
      integer, parameter :: nreq = 1
      integer(kind=4), parameter :: tag_ls = 1, tag_gpt = tag_ls+1

      call inflate_req(nreq)

      ! copy the patches data to a temporary array to be sent to the master
      if (slave) call MPI_Isend(ls, I_ONE, MPI_INTEGER, FIRST, tag_ls, comm, req(nreq), mpi_err)

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
            call MPI_Recv(ls, I_ONE, MPI_INTEGER, p, tag_ls, comm, MPI_STATUS_IGNORE, mpi_err)
            if (ls > 0) then
               allocate(gptemp(I_OFF:I_END,ls))
               call MPI_Recv(gptemp, size(gptemp), MPI_INTEGER8, p, tag_gpt, comm, MPI_STATUS_IGNORE, mpi_err)
               do ss = 1, ls
                  call gp%list(i+ss)%set_gp(gptemp(I_OFF:I_OFF+ndims-1, ss), int(gptemp(I_N_B:I_N_B+ndims-1, ss), kind=4), INVALID, p)
               enddo
               i = i + ls
               deallocate(gptemp)
            endif
         enddo
      else

         ! send patches to master
         call MPI_Wait(req(nreq), status(:, nreq), mpi_err)
         if (ls > 0) call MPI_Send(gptemp, size(gptemp), MPI_INTEGER8, FIRST, tag_gpt, comm, mpi_err)
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
      use mpi,             only: MPI_INTEGER, MPI_INTEGER8, MPI_STATUS_IGNORE
      use mpisetup,        only: master, FIRST, LAST, comm, mpi_err
      use sort_piece_list, only: grid_piece_list

      implicit none

      class(cg_list_balance_T),                 intent(inout) :: this !< object invoking type bound procedure
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
            call MPI_Send(ls, I_ONE, MPI_INTEGER, p, tag_lsR, comm, mpi_err)
            if (ls>0) then
               allocate(gptemp(I_OFF:I_END,from(p):from(p+1)-1))
               do s = lbound(gptemp, dim=2, kind=4), ubound(gptemp, dim=2, kind=4)
                  gptemp(:, s) = [ gp%list(s)%off, int(gp%list(s)%n_b, kind=8) ]
               enddo
               call MPI_Send(gptemp, size(gptemp), MPI_INTEGER8, p, tag_gptR, comm, mpi_err)
               deallocate(gptemp)
            endif
         enddo
      else
         ! receive new, perhaps more balanced patches
         call MPI_Recv(ls, I_ONE, MPI_INTEGER, FIRST, tag_lsR, comm, MPI_STATUS_IGNORE, mpi_err)
         if (ls>0) then
            allocate(gptemp(I_OFF:I_END,ls))
            call MPI_Recv(gptemp, size(gptemp), MPI_INTEGER8, FIRST, tag_gptR, comm, MPI_STATUS_IGNORE, mpi_err)
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

      class(cg_list_balance_T), target,  intent(inout) :: this  !< current level
      integer(kind=8), dimension(ndims), intent(in)    :: n_d   !< number of grid cells
      integer(kind=8), dimension(ndims), intent(in)    :: off   !< offset (with respect to the base level, counted on own level)

      this%recently_changed = .true. ! assume that the new patches will change this level
      call this%plist%expand
      call this%plist%patches(ubound(this%plist%patches(:), dim=1))%one_piece_patch(n_d(:), off(:))

   end subroutine add_patch_one_piece

!> \brief Sort list according to SFC id

   subroutine sort_SFC(this)

      use cg_list,      only: cg_list_element
      use sort_cg_list, only: sort_cg_list_T

      implicit none

      class(cg_list_balance_T), intent(inout) :: this !< object invoking type-bound procedure

      type(cg_list_element), pointer :: cgl
      integer :: s
      type(sort_cg_list_T) :: l

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
         end if
         if (s < ubound(l%list, dim=1)) then
            l%list(s)%cgl%nxt => l%list(s+1)%cgl
         else
            l%list(s)%cgl%nxt => null()
         end if
         l%list(s)%cgl%cg%grid_id = s
      end do
      this%first => l%list(lbound(l%list, dim=1))%cgl
      this%last  => l%list(ubound(l%list, dim=1))%cgl

      call l%cleanup

   end subroutine sort_SFC

end module cg_list_balance
