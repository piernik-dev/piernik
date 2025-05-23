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
!! \brief This module contains list of grid containers that belong to a single level.
!! Everything that relies on connection between levels is placed in cg_level_connected module.
!<

module cg_level
#if defined(__INTEL_COMPILER)
   !! \deprecated remove this clause as soon as Intel Compiler gets required
   !! features and/or bug fixes, it's needed for 12.1, fixed in 13.0 but the
   !! latter is broken and we cannot use it yet
   use cg_list,           only: cg_list_t   ! QA_WARN intel
#endif /* __INTEL_COMPILER */
   use cg_list_neighbors, only: cg_list_neighbors_t

   implicit none

   private
   public :: cg_level_t

   !>
   !! \brief A list of all cg of the same resolution.
   !!
   !! \details For positive refinement levels the list may be composed of several disconnected subsets of cg
   !! (islands: made of one or more cg's).
   !! This type is not intended for direct use. It is extended in cg_level_connected into a functional object.
   !<
   type, extends(cg_list_neighbors_t), abstract :: cg_level_t

      integer                                    :: fft_type     !< type of FFT to employ in some multigrid solvers (depending on boundaries)

      ! FARGO
      real,    dimension(:, :), allocatable      :: omega_mean   !< mean angular velocity for each fluid
      real,    dimension(:, :), allocatable      :: omega_cr     !< constant residual angular velocity for each fluid
      integer, dimension(:, :), allocatable      :: nshift       !< number of cells that need to be shifted due to %omega_mean for each fluid
      real,    dimension(:, :), allocatable      :: local_omega  !< auxiliary array
      integer(kind=8), dimension(:), allocatable :: cell_count   !< auxiliary counter

   contains

      procedure          :: cleanup                                              !< deallocate arrays
      procedure          :: init_all_new_cg                                      !< initialize newest grid container
      procedure          :: print_segments                                       !< print detailed information about current level decomposition
      procedure, private :: update_decomposition_properties                      !< Update some flags in domain module
      procedure, private :: create                                               !< Get all decomposed patches and turn them into local grid containers
      generic,   public  :: add_patch => add_patch_fulllevel, add_patch_detailed !< Add a new piece of grid to the current level and decompose it
      procedure, private :: add_patch_fulllevel                                  !< Add a whole level to the list of patches
      procedure, private :: add_patch_detailed                                   !< Add a new piece of grid to the list of patches
      procedure          :: deallocate_patches                                   !< Throw out patches list
      procedure, private :: update_everything                                    !< Update all information on refinement structure and intra-level communication
      procedure          :: check_update_all                                     !< Check if it is necessary to call update_everything
      procedure          :: refresh_SFC_id                                       !< Recalculate SFC_id for grids, useful after domain expansion

   end type cg_level_t

contains

!> \brief deallocate arrays

   subroutine cleanup(this)

      implicit none

      class(cg_level_t), intent(inout) :: this !< object invoking type bound procedure

      call this%plist%p_deallocate
      call this%dot%cleanup
      if (allocated(this%omega_mean))   deallocate(this%omega_mean)
      if (allocated(this%omega_cr))     deallocate(this%omega_cr)
      if (allocated(this%nshift))       deallocate(this%nshift)
      if (allocated(this%local_omega))  deallocate(this%local_omega)
      call this%l%check
      deallocate(this%l)

   end subroutine cleanup

!>
!! \brief Throw out patches list
!!
!! \todo check if it is really necessary to call this from remote routines
!<

   subroutine deallocate_patches(this)

      implicit none

      class(cg_level_t), intent(inout) :: this !< object invoking type bound procedure

      call this%plist%p_deallocate

   end subroutine deallocate_patches

!> \brief Print detailed information about current level decomposition

   subroutine print_segments(this)

      use cg_list,    only: cg_list_element
      use constants,  only: LO, HI
      use dataio_pub, only: printinfo !, msg, warn
      use mpisetup,   only: FIRST, LAST, master!, nproc, proc
#ifdef VERBOSE
      use constants,  only: V_DEBUG
      use dataio_pub, only: msg
#endif /* VERBOSE */

      implicit none

      class(cg_level_t), intent(in)   :: this   !< object invoking type bound procedure

      integer                         :: p, i, hl, tot_cg
      integer(kind=8)                 :: ccnt
      real, allocatable, dimension(:) :: maxcnt
      type(cg_list_element), pointer  :: cgl
#ifdef VERBOSE
      character(len=len(msg))         :: header
#endif /* VERBOSE */

      i = 0
      cgl => this%first
      do while (associated(cgl))
         i = i + 1
         cgl => cgl%nxt
      enddo
!!$      if (i /= this%cnt .or. this%cnt /= size(this%dot%gse(proc)%c(:)) .or. size(this%dot%gse(proc)%c(:)) /= i) then
!!$         write(msg, '(2(a,i4),a,3i7)')"[cg_level:print_segments] Uncertain number of grid pieces @PE ",proc," on level ", this%l%id, &
!!$              &                       " : ",i,this%cnt,size(this%dot%gse(proc)%c(:))
!!$         call warn(msg)
!!$
!!$         cgl => this%first
!!$         do while (associated(cgl))
!!$            write(msg,'(2(a,i7),2(a,3i10),a)')" @",proc," #",cgl%cg%grid_id," : [", cgl%cg%my_se(:, LO), "] : [", cgl%cg%my_se(:, HI)," ]"
!!$            call printinfo(msg, V_DEBUG)
!!$            cgl => cgl%nxt
!!$         enddo
!!$      endif

      if (.not. master) return

      !call dom%print_me

      ! print segments according to list of patches
      allocate(maxcnt(FIRST:LAST))
      maxcnt(:) = 0
      tot_cg = 0
      do p = FIRST, LAST
         hl = 0
         tot_cg = tot_cg + size(this%dot%gse(p)%c(:))
         do i = lbound(this%dot%gse(p)%c(:), dim=1), ubound(this%dot%gse(p)%c(:), dim=1)
            ccnt = product(this%dot%gse(p)%c(i)%se(:, HI) - this%dot%gse(p)%c(i)%se(:, LO) + 1)
            maxcnt(p) = maxcnt(p) + ccnt
#ifdef VERBOSE
            if (i == 1) then
               write(header, '(a,i4,a,i3)')"[cg_level:print_segments] segment @", p, " ^", this%l%id
               hl = len_trim(header)
            else
               header = repeat(" ", hl)
            endif
            if (maxval(this%l%n_d(:)) < 1000000) then
               write(msg,'(2a,2(3i7,a),i8,a)') header(:hl), " : [", this%dot%gse(p)%c(i)%se(:, LO), "] : [", this%dot%gse(p)%c(i)%se(:, HI), "] #", ccnt, " cells"
            else if (maxval(this%l%n_d(:)) < 1000000000) then
               write(msg,'(2a,2(3i10,a),i8,a)') header(:hl), " : [", this%dot%gse(p)%c(i)%se(:, LO), "] : [", this%dot%gse(p)%c(i)%se(:, HI), "] #", ccnt, " cells"
            else
               write(msg,'(2a,2(3i18,a),i8,a)') header(:hl), " : [", this%dot%gse(p)%c(i)%se(:, LO), "] : [", this%dot%gse(p)%c(i)%se(:, HI), "] #", ccnt, " cells"
            endif
            call printinfo(msg, V_DEBUG)
#endif /* VERBOSE */
         enddo
      enddo

!      write(msg, '(a,i3,a,f5.1,a,i5,a,f8.5)')"[cg_level:print_segments] Level ", this%l%id, " filled in ",(100.*sum(maxcnt(:)))/product(real(this%n_d(:))), &
!           &                                 "%, ",tot_cg," grid(s), load balance : ", sum(maxcnt(:))/(nproc*maxval(maxcnt(:)))
      !> \todo add calculation of total internal boundary surface in cells
!      call printinfo(msg, V_DEBUG)
      deallocate(maxcnt)

   end subroutine print_segments

!>
!! \brief Initialize all grid containers on a new grid level
!!
!! \todo automagically rebalance existing grids unless it is explicitly forbidden
!<

   subroutine init_all_new_cg(this)

      implicit none

      class(cg_level_t), intent(inout) :: this !< object invoking type bound procedure

      ! First: do the balancing of new grids
      call this%balance_new

      ! Second: create new grids, invalidate most of this%dot%gse database
      call this%create

      ! Third: update all information on refinement structure and intra-level communication, update this%dot%gse database
      ! Remember that the communication between levels has to be updated as well, but we cannot do this here due to cyclic dependencies
      call this%update_everything

   end subroutine init_all_new_cg

!> \brief update all information on refinement structure and intra-level communication.

   subroutine update_everything(this)

      use constants, only: PPP_AMR
      use ppp,       only: ppp_main

      implicit none

      class(cg_level_t), intent(inout) :: this   !< object invoking type bound procedure

      character(len=*), parameter :: lue_label = "level_update_everything"

      call ppp_main%start(lue_label, PPP_AMR)

      call this%update_decomposition_properties
      call this%dot%update_global(this%first, this%cnt, this%l%off) ! communicate everything that was added before
      call this%dot%update_SFC_id_range(this%l%off)
      call this%find_neighbors ! requires access to whole this%dot%gse(:)%c(:)%se(:,:)
      call this%dot%update_tot_se
      call this%print_segments

      call ppp_main%stop(lue_label, PPP_AMR)

   end subroutine update_everything

!>
!! \brief Get all decomposed patches and turn them into local grid containers
!!
!! \details This routine starts with two lists:
!! * A list of blocks that survived derefinement attempts
!! * A list of blocks to be created due to requested refinement
!! The goal is to merge the two lists, create grid containers for the new grids and fix things such as cg%grid_id if anything was derefined
!<

   subroutine create(this)

      use cg_list_global,     only: all_cg
      use constants,          only: PPP_AMR
      use grid_cont,          only: grid_container
      use grid_container_ext, only: cg_extptrs
      use dataio_pub,         only: die
      use mpisetup,           only: proc
      use ppp,                only: ppp_main

      implicit none

      class(cg_level_t), intent(inout) :: this   !< object invoking type bound procedure

      integer                       :: i, p, ep
      integer(kind=8)               :: s
      type(grid_container), pointer :: cg
      character(len=*), parameter   :: gc_label = "init_gc"

      ! Find how many pieces are to be added and recreate local gse and make room for new pieces in the gse array
      call this%dot%update_local(this%first, int(this%cnt + this%plist%p_count(), kind=4))

      call this%dot%is_consitent(this%first) ! check local consistency

      call ppp_main%start(gc_label, PPP_AMR)
      ! create the new grid pieces
      i = this%cnt
      if (allocated(this%plist%patches)) then
         do p = lbound(this%plist%patches(:), dim=1), ubound(this%plist%patches(:), dim=1)
            do s = lbound(this%plist%patches(p)%pse, dim=1), ubound(this%plist%patches(p)%pse, dim=1)
               i = i + 1
               if (i > size(this%dot%gse(proc)%c(:))) call die("[cg_level:create] overflow")
               this%dot%gse(proc)%c(i)%se(:,:) = this%plist%patches(p)%pse(s)%se(:,:)
               call this%add
               cg => this%last%cg
               call cg%init_gc(this%dot%gse(proc)%c(i)%se(:, :), i, this%l)
               do ep = lbound(cg_extptrs%ext, dim=1), ubound(cg_extptrs%ext, dim=1)
                  if (associated(cg_extptrs%ext(ep)%init))  call cg_extptrs%ext(ep)%init(cg)
               enddo
               call all_cg%add(cg)
            enddo
         enddo
         call this%plist%p_deallocate
      endif
      call ppp_main%stop(gc_label, PPP_AMR)

      call this%sort_SFC

   end subroutine create

!> \brief Update some flags in domain module [ is_uneven, is_mpi_noncart, is_refined, is_multicg ]

   subroutine update_decomposition_properties(this)

      use allreduce, only: piernik_MPI_Allreduce
      use constants, only: base_level_id, pLOR
      use domain,    only: is_mpi_noncart, is_multicg, is_refined, is_uneven
      use mpisetup,  only: proc, proc

      implicit none

      class(cg_level_t), intent(inout) :: this   !< object invoking type bound procedure

      if (this%l%id > base_level_id) is_refined = .true.
      call piernik_MPI_Allreduce(is_refined, pLOR)
      if (is_refined) then
         is_mpi_noncart = .true.
         is_multicg = .true.
      endif
      if (is_mpi_noncart) is_uneven = .true.

      is_multicg = is_multicg .or. (ubound(this%dot%gse(proc)%c(:), dim=1) > 1)
      call piernik_MPI_Allreduce(is_multicg, pLOR)

      call this%dot%check_blocky ! check if all blocks in the domain have same size and shape

   end subroutine update_decomposition_properties

!> \brief Add a whole level to the list of patches on current refinement level and decompose it.

   subroutine add_patch_fulllevel(this, n_pieces)

      implicit none

      class(cg_level_t), target, intent(inout) :: this     !< current level
      integer(kind=4), optional, intent(in)    :: n_pieces !< how many pieces the patch should be divided to?

      call this%add_patch_detailed(this%l%n_d, this%l%off, n_pieces)

   end subroutine add_patch_fulllevel

!>
!! \brief Add a new piece of grid to the list of patches on current refinement level and decompose it
!!
!! \details Each patch should be added only once, i.e. the base level and multigrid coarse grids have to be added by the master and then redistributed in balance_new.
!<

   subroutine add_patch_detailed(this, n_d, off, n_pieces)

      use constants,  only: ndims
      use dataio_pub, only: msg, die

      implicit none

      class(cg_level_t), target,         intent(inout) :: this     !< current level
      integer(kind=8), dimension(ndims), intent(in)    :: n_d      !< number of grid cells
      integer(kind=8), dimension(ndims), intent(in)    :: off      !< offset (with respect to the base level, counted on own level)
      integer(kind=4), optional,         intent(in)    :: n_pieces !< how many pieces the patch should be divided to?

      this%recently_changed = .true. ! assume that the new patches will change this level
      call this%plist%expand
      if (.not. this%plist%patches(ubound(this%plist%patches(:), dim=1))%decompose_patch(n_d(:), off(:), this%l%id, n_pieces=n_pieces)) then
         write(msg,'(a,i4)')"[cg_level:add_patch_detailed] Decomposition failed at level ",this%l%id
         call die(msg)
      endif

   end subroutine add_patch_detailed

!> \brief Check if it is necessary to call update_everything

   subroutine check_update_all(this)

      use allreduce, only: piernik_MPI_Allreduce
      use constants, only: pLOR

      implicit none

      class(cg_level_t), intent(inout) :: this

      ! OPT: call this%update_gse inside this%update_everything can be quite long to complete
      call piernik_MPI_Allreduce(this%recently_changed, pLOR)
      if (this%recently_changed) call this%update_everything

   end subroutine check_update_all

!> \brief Recalculate SFC_id for grids, useful after domain expansion

   subroutine refresh_SFC_id(this)

      use cg_list,   only: cg_list_element
      use constants, only: LO
      use ordering,  only: SFC_order

      implicit none

      class(cg_level_t), intent(inout) :: this

      type(cg_list_element), pointer  :: cgl

      cgl => this%first
      do while (associated(cgl))
         cgl%cg%SFC_id = SFC_order(cgl%cg%my_se(:, LO) - this%l%off)
         cgl => cgl%nxt
      enddo

   end subroutine refresh_SFC_id

end module cg_level
