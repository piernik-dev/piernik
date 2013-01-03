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

!>
!! \brief This module contains list of grid containers that belong to a single level level.
!! Everything that relies on connection between levels is placed in cg_level_connected module.
!<

module cg_level
#if defined(__INTEL_COMPILER)
   !! \deprecated remove this clause as soon as Intel Compiler gets required
   !! features and/or bug fixes, it's needed for 12.1, fixed in 13.0 but the
   !! latter is broken and we cannot use it yet
   use cg_list,       only: cg_list_T   ! QA_WARN intel
#endif /* __INTEL_COMPILER */
   use cg_list_bnd,   only: cg_list_bnd_T
   use constants,     only: ndims, LO, HI
   use decomposition, only: box_T

   implicit none

   private
   public :: cg_level_T

   !> \brief A single grid piece plus auxiliary data
   !> \deprecated not to be confused with decomposition::cuboid
   type :: cuboid
      integer(kind=8), dimension(ndims, LO:HI) :: se     !< absolute index of grid segmenent on its refinement level wrt [0,0,0]
      logical                                  :: is_new !< a flag that marks newly added grid pieces
   end type cuboid

   !> \brief A list of grid pieces (typically used as a list of all grids residing on a given process)
   type :: cuboids
      type(cuboid), allocatable, dimension(:) :: c !< an array of grid piece
   end type cuboids

   !>
   !! \brief A list of all cg of the same resolution.
   !!
   !! \details For positive refinement levels the list may be composed of several disconnected subsets of cg ("islands: made of one or more cg: cg_list_patch).
   !! This type is not intended for direct use. It is extended in cg_level_connected into a functional object.
   !<
   type, extends(cg_list_bnd_T) :: cg_level_T

      integer(kind=4)                            :: level_id  !< level number (relative to base level). For printing, debug, and I/O use only. No arithmetic should depend on it.
      integer(kind=8), dimension(ndims)          :: n_d       !< maximum number of grid cells in each direction (size of fully occupied level)
      type(cuboids),   dimension(:), allocatable :: pse       !< lists of grid chunks on each process (FIRST:LAST); Use with care, because this is an antiparallel thing
      integer                                    :: tot_se    !< global number of segments on the level
      integer                                    :: fft_type  !< type of FFT to employ in some multigrid solvers (depending on boundaries)
      type(box_T),     dimension(:), allocatable :: patches   !< list of patches that exist on the current level
      integer(kind=8), dimension(ndims)          :: off       !< offset of the level

    contains

      procedure          :: init_all_new_cg                   !< initialize newest grid container
      procedure, private :: mpi_bnd_types                     !< create MPI types for boundary exchanges
      procedure          :: print_segments                    !< print detailed information about current level decomposition
      procedure, private :: update_decomposition_properties   !< Update some flags in domain module
      procedure, private :: distribute                        !< Get all decomposed patches and compute which pieces go to which process
      procedure, private :: calc_ord_range                    !< Compute which id\'s should belong to which process
      procedure, private :: simple_ordering                   !< This is just counting, not ordering
      procedure, private :: mark_new                          !< Detect which grid containers are new
      procedure          :: update_tot_se                     !< count all cg on current level for computing tags in vertical_prep
      generic,   public  :: add_patch => add_patch_fulllevel, add_patch_detailed !< Add a new piece of grid to the current level and decompose it
      procedure, private :: add_patch_fulllevel               !< Add a whole level to the list of patches
      procedure, private :: add_patch_detailed                !< Add a new piece of grid to the list of patches

   end type cg_level_T

contains

!> \brief Print detailed information about current level decomposition

   subroutine print_segments(this)

      use cg_list,    only: cg_list_element
      use constants,  only: LO, HI
      use dataio_pub, only: printinfo, msg, warn
      use mpisetup,   only: FIRST, LAST, master, nproc, proc

      implicit none

      class(cg_level_T), intent(in)   :: this   !< object invoking type bound procedure

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
      if (i /= this%cnt .or. this%cnt /= size(this%pse(proc)%c(:)) .or. size(this%pse(proc)%c(:)) /= i) then
         write(msg, '(2(a,i4),a,3i7)')"[cg_level:print_segments] Uncertain number of grid pieces @PE ",proc," on level ", this%level_id, &
              &                       " : ",i,this%cnt,size(this%pse(proc)%c(:))
         call warn(msg)

         cgl => this%first
         do while (associated(cgl))
            write(msg,'(2(a,i7),2(a,3i10),a)')" @",proc," #",cgl%cg%grid_id," : [", cgl%cg%my_se(:, LO), "] : [", cgl%cg%my_se(:, HI)," ]"
            call printinfo(msg)
            cgl => cgl%nxt
         enddo
      endif

      if (.not. master) return

      !call dom%print_me

      ! print segments according to list of patches
      allocate(maxcnt(FIRST:LAST))
      maxcnt(:) = 0
      tot_cg = 0
      do p = FIRST, LAST
         hl = 0
         tot_cg = tot_cg + size(this%pse(p)%c(:))
         do i = lbound(this%pse(p)%c(:), dim=1), ubound(this%pse(p)%c(:), dim=1)
            ccnt = product(this%pse(p)%c(i)%se(:, HI) - this%pse(p)%c(i)%se(:, LO) + 1)
            maxcnt(p) = maxcnt(p) + ccnt
#ifdef VERBOSE
            if (i == 1) then
               write(header, '(a,i4)')"[cg_level:print_segments] segment @", p
               hl = len_trim(header)
            else
               header = repeat(" ", hl)
            endif
            if (maxval(this%n_d(:)) < 1000000) then
               write(msg,'(2a,2(3i7,a),i8,a)') header(:hl), " : [", this%pse(p)%c(i)%se(:, LO), "] : [", this%pse(p)%c(i)%se(:, HI), "] #", ccnt, " cells"
            else if (maxval(this%n_d(:)) < 1000000000) then
               write(msg,'(2a,2(3i10,a),i8,a)') header(:hl), " : [", this%pse(p)%c(i)%se(:, LO), "] : [", this%pse(p)%c(i)%se(:, HI), "] #", ccnt, " cells"
            else
               write(msg,'(2a,2(3i18,a),i8,a)') header(:hl), " : [", this%pse(p)%c(i)%se(:, LO), "] : [", this%pse(p)%c(i)%se(:, HI), "] #", ccnt, " cells"
            endif
            call printinfo(msg)
#endif /* VERBOSE */
         enddo
      enddo

      write(msg, '(a,i3,a,f5.1,a,i5,a,f8.5)')"[cg_level:print_segments] Level ", this%level_id, " filled in ",(100.*sum(maxcnt(:)))/product(real(this%n_d(:))), &
           &                                 "%, ",tot_cg," grid(s), load balance : ", sum(maxcnt(:))/(nproc*maxval(maxcnt(:)))
      !> \todo add calculation of total internal boundary surface in cells
      call printinfo(msg)
      deallocate(maxcnt)

   end subroutine print_segments

!> \brief Initialize all grid containers on a new grid level

   subroutine init_all_new_cg(this)

      use cg_list_global,     only: all_cg
      use grid_cont,          only: grid_container
      use grid_container_ext, only: cg_extptrs
      use mpisetup,           only: proc

      implicit none

      class(cg_level_T), intent(inout) :: this   !< object invoking type bound procedure

      integer                          :: i, ep
      type(grid_container), pointer    :: cg

      call this%distribute
      call this%mark_new

      do i = lbound(this%pse(proc)%c(:), dim=1), ubound(this%pse(proc)%c(:), dim=1)
         if (this%pse(proc)%c(i)%is_new) then
            call this%add
            cg => this%last%cg
            call cg%init(this%n_d, this%off, this%pse(proc)%c(i)%se(:, :), i, this%level_id) ! we cannot pass "this" as an argument because of circular dependencies
            do ep = lbound(cg_extptrs%ext, dim=1), ubound(cg_extptrs%ext, dim=1)
               if (associated(cg_extptrs%ext(ep)%init))  call cg_extptrs%ext(ep)%init(cg)
            enddo
            call all_cg%add(cg)
         endif
      enddo

      call this%mpi_bnd_types ! require access to whole this%pse(:)%c(:)%se(:,:)

      call this%update_req    ! Perhaps this%mpi_bnd_types added some new entries
      call this%update_tot_se
      call this%print_segments

   end subroutine init_all_new_cg

!> \brief Detect which grid containers are new

   subroutine mark_new(this)

      use cg_list,    only: cg_list_element
      use dataio_pub, only: warn
      use mpisetup,   only: proc

      implicit none

      class(cg_level_T), intent(inout) :: this   !< object invoking type bound procedure

      type(cg_list_element), pointer   :: cgl
      integer                          :: i, occured

      do i = lbound(this%pse(proc)%c(:), dim=1), ubound(this%pse(proc)%c(:), dim=1)
         this%pse(proc)%c(i)%is_new = .true.
      enddo

      cgl => this%first
      do while (associated(cgl))

         occured = 0
         do i = lbound(this%pse(proc)%c(:), dim=1), ubound(this%pse(proc)%c(:), dim=1)
            if (all(this%pse(proc)%c(i)%se(:,:) == cgl%cg%my_se(:,:))) then
               this%pse(proc)%c(i)%is_new = .false.
               occured = occured + 1
            endif
         enddo

         if (occured <= 0) then
            call warn("[cg_level:mark_new] Existing cg not found in new list")
         else if (occured > 1) then
            call warn("[cg_level:mark_new] Existing cg found multiple times in new list")
         endif

         cgl => cgl%nxt
      enddo

   end subroutine mark_new

!>
!! \brief Get all decomposed patches and compute which pieces go to which process
!!
!! \details This routine starts with two lists:
!! * A list of blocks that survived derefinement attempts
!! * A list of blocks due to requested refinement
!! It has to decide if and how to do migration of grid pieces and to communicate this update global database of grid pieces.
!!
!! There are several strategies than can be implemented:
!! * Local refinements go to local process. It is very simple, but for most simulations will build up load imbalance. Suitable for tests and global refinement.
!! * Local refinements can be assigned to remote processes, existing blocks stays in place. Should keep good load balance, but the amount of inter-process
!!   internal boundariem may grow significantly with time. Suitable for minor refinement updates and base level decomposition.
!! * All blocks (existing and new) have recalculated assignment and can be migrated to other processes. Most advanced. Should be used after reading restart data.
!!
!! First startegy will be implemented first to get everything working. Second strategy will be used quite often. Third one do not need to be used on every refinement update.
!! It can be called when some benchmark of grid disorder exceeds particuklar threshold.
!<

   subroutine distribute(this)

      use dataio_pub, only: die
      use mpisetup,   only: FIRST, LAST

      implicit none

      class(cg_level_T), intent(inout)       :: this   !< object invoking type bound procedure

      integer                                :: i, p, s
      integer(kind=8), dimension(FIRST:LAST) :: min_id, max_id, pieces, filled

      call this%simple_ordering
      call this%calc_ord_range(min_id, max_id, pieces)
      if (.not. allocated(this%pse)) allocate(this%pse(FIRST:LAST))
      do i = FIRST, LAST
         if (allocated(this%pse(i)%c)) deallocate(this%pse(i)%c) !> \todo recycle previous list somehow?
         allocate(this%pse(i)%c(pieces(i)))
      enddo
      filled(:) = 0

      ! write the new grid pieces description to the pse array
      if (allocated(this%patches)) then
         do p = lbound(this%patches(:), dim=1), ubound(this%patches(:), dim=1)
            do s = lbound(this%patches(p)%pse, dim=1), ubound(this%patches(p)%pse, dim=1)
               do i = FIRST, LAST
                  if (this%patches(p)%pse(s)%id >= min_id(i) .and. this%patches(p)%pse(s)%id <= max_id(i)) then
                     filled(i) = filled(i) + 1
                     if (filled(i) > size(this%pse(i)%c(:))) call die("[cg_level:distribute] overflow")
                     this%pse(i)%c(filled(i)) = cuboid(this%patches(p)%pse(s)%se(:,:), .true.) ! The is_new flag will be fixed in mark_new subroutine
                     exit
                  endif
               enddo
            enddo
         enddo
      endif

      call this%update_decomposition_properties

   end subroutine distribute

!>
!! \brief Compute which id\'s should belong to which process
!!
!! \todo Reduce assumptions on the set of id only to uniqueness (i.e. some id might be absent, some might be <0 or >this%tot_se, perform partial sort (qsort? shell sort?))
!<

   subroutine calc_ord_range(this, min_id, max_id, pieces)

      use mpisetup, only: nproc, FIRST, LAST

      implicit none

      class(cg_level_T),                      intent(inout) :: this   !< object invoking type bound procedure
      integer(kind=8), dimension(FIRST:LAST), intent(out)   :: min_id !< \todo comment me
      integer(kind=8), dimension(FIRST:LAST), intent(out)   :: max_id !< \todo comment me
      integer(kind=8), dimension(FIRST:LAST), intent(out)   :: pieces !< \todo comment me

      integer :: p

      ! At the moment all id are in [0 .. this%tot_se-1] range. This will change in future
      !> \todo implement different weights of pieces

      do p = FIRST, LAST
         min_id(p) = ( p      * this%tot_se) / nproc
         max_id(p) = ((p + 1) * this%tot_se) / nproc - 1
      enddo
      pieces(:) = max_id(:) - min_id(:) + 1

   end subroutine calc_ord_range

!> \brief OMG! This is just counting, not ordering!
!> \todo Implement anything better. Peano-Hilbert ordering will appear here sooner or later :-)

   subroutine simple_ordering(this)

      implicit none

      class(cg_level_T), intent(inout) :: this   !< object invoking type bound procedure

      integer                          :: p, s
      integer(kind=8)                  :: id

      id = 0
      if (allocated(this%patches)) then
         do p = lbound(this%patches(:), dim=1), ubound(this%patches(:), dim=1)
            if (allocated(this%patches(p)%pse)) then
               do s = lbound(this%patches(p)%pse, dim=1), ubound(this%patches(p)%pse, dim=1)
                  this%patches(p)%pse(s)%id = id
                  id = id + 1
               enddo
            endif
         enddo
      endif

      this%tot_se = int(id, kind=4) ! obsoletes update_tot_se ?

   end subroutine simple_ordering

!> \brief Update some flags in domain module [ is_uneven, is_mpi_noncart, is_refined, is_multicg ]

   subroutine update_decomposition_properties(this)

      use constants,  only: base_level_id, pLOR
      use dataio_pub, only: warn
      use domain,     only: is_mpi_noncart, is_multicg, is_refined, is_uneven
      use mpisetup,   only: proc, master, piernik_MPI_Allreduce

      implicit none

      class(cg_level_T), intent(inout) :: this   !< object invoking type bound procedure

      if (this%level_id > base_level_id) is_refined = .true.
      call piernik_MPI_Allreduce(is_refined, pLOR)
      if (is_refined) then
         is_mpi_noncart = .true.
         is_multicg = .true.
         if (master) call warn("[cg_level:update_decomposition_properties] Refinements are experimental")
      endif
      if (is_mpi_noncart) is_uneven = .true.

      is_multicg = is_multicg .or. (ubound(this%pse(proc)%c(:), dim=1) > 1)
      call piernik_MPI_Allreduce(is_multicg, pLOR)

   end subroutine update_decomposition_properties

!> \brief Count all cg on current level. Useful for computing tags in vertical_prep

   subroutine update_tot_se(this)

      use mpisetup, only: FIRST, LAST

      implicit none

      class(cg_level_T), intent(inout) :: this   !< object invoking type bound procedure

      integer :: p

      this%tot_se = 0
      do p = FIRST, LAST
         if (allocated(this%pse)) this%tot_se = this%tot_se + ubound(this%pse(p)%c(:), dim=1)
      enddo

   end subroutine update_tot_se

!>
!! \brief Create MPI types for boundary exchanges
!!
!! \details this type can be a member of grid container type if we pass this%pse(:)%c(:) as an argument.
!! It would simplify dependencies and this%init_all_new_cg, but it could be quite a big object.
!! Assume that cuboids don't collide (no overlapping grid pieces on same refinement level are allowed)
!! \todo Put this%pse into a separate type and pass a pointer to it or even a pointer to pre-filtered segment list
!!
!! \todo Do not provide segments for each possible number of guardcells. Provide 1 layer, 1 layer with corners, all layers and all layers with cofrners instead.
!<

   subroutine mpi_bnd_types(this)

      use cg_list,          only: cg_list_element
      use constants,        only: xdim, ydim, zdim, ndims, LO, HI, I_ONE, BND_MPI_FC, BND_FC
      use domain,           only: dom
      use grid_cont,        only: grid_container, is_overlap
      use mpisetup,         only: FIRST, LAST, procmask

      implicit none

      class(cg_level_T), intent(inout)                 :: this    !< object invoking type bound procedure

      type(grid_container),  pointer                   :: cg      !< grid container that we are currently working on
      type(cg_list_element), pointer                   :: cgl
      integer                                          :: g, j, b
      integer(kind=8)                                  :: n_tot_face_cells
      integer(kind=8), dimension(LO:HI)                :: n_lbnd_face_cells
      integer(kind=4)                                  :: d, dd, hl, lh, ib
      integer(kind=8), dimension(xdim:zdim)            :: per
      integer(kind=8), dimension(xdim:zdim, LO:HI)     :: b_layer, bp_layer, poff
      logical,         dimension(:,:,:,:), allocatable :: facemap

      cgl => this%first
      do while (associated(cgl))
         cg => cgl%cg

         if (allocated(cg%i_bnd)) deallocate(cg%i_bnd)
         if (allocated(cg%o_bnd)) deallocate(cg%o_bnd)
         allocate(cg%i_bnd(xdim:zdim, dom%nb), cg%o_bnd(xdim:zdim, dom%nb))

         ! assume that cuboids fill the domain and don't collide
         per(:) = 0
         where (dom%periodic(:)) per(:) = this%n_d(:)

         do d = xdim, zdim
            if (dom%has_dir(d)) then

               ! identify processes with interesting neighbour data
               procmask(:) = 0
               do lh = LO, HI
                  hl = LO+HI-lh ! HI for LO, LO for HI
                  b_layer(:,:) = cg%my_se(:, :)
                  b_layer(d, lh) = b_layer(d, lh) + lh-hl ! -1 for LO, +1 for HI
                  b_layer(d, hl) = b_layer(d, lh) ! adjacent boundary layer, 1 cell wide, without corners
                  do j = FIRST, LAST
                     do b = lbound(this%pse(j)%c(:), dim=1), ubound(this%pse(j)%c(:), dim=1)
                        if (is_overlap(b_layer(:,:), this%pse(j)%c(b)%se(:,:), per(:))) procmask(j) = procmask(j) + 1 ! count how many boundaries we share with that process
                     enddo
                  enddo
               enddo
               do ib = 1, dom%nb
                  allocate(cg%i_bnd(d, ib)%seg(sum(procmask(:))), cg%o_bnd(d, ib)%seg(sum(procmask(:))))
               enddo

               ! set up segments to be sent or received
               g = 0
               do j = FIRST, LAST
                  if (procmask(j) /= 0) then
                     do lh = LO, HI
                        hl = LO+HI-lh
                        do b = lbound(this%pse(j)%c(:), dim=1), ubound(this%pse(j)%c(:), dim=1)
                           b_layer(:,:) = cg%my_se(:, :)
                           b_layer(d, lh) = b_layer(d, lh) + lh-hl
                           b_layer(d, hl) = b_layer(d, lh) ! adjacent boundary layer, 1 cell thick, without corners
                           bp_layer(:, :) = b_layer(:, :)
                           where (per(:) > 0)
                              where (bp_layer(:, LO) < this%off(:)) bp_layer(:, LO) = bp_layer(:, LO) + per(:)
                              where (bp_layer(:, HI) < this%off(:)) bp_layer(:, HI) = bp_layer(:, HI) + per(:)
                              where (bp_layer(:, LO) >= this%off(:)+this%n_d(:)) bp_layer(:, LO) = bp_layer(:, LO) - per(:)
                              where (bp_layer(:, HI) >= this%off(:)+this%n_d(:)) bp_layer(:, HI) = bp_layer(:, HI) - per(:)
                           endwhere
                           !> \todo save b_layer(:,:) and bp_layer(:,:) and move above calculations outside the b loop

                           if (is_overlap(bp_layer(:,:), this%pse(j)%c(b)%se(:,:))) then
                              poff(:,:) = bp_layer(:,:) - b_layer(:,:) ! displacement due to periodicity
                              bp_layer(:, LO) = max(bp_layer(:, LO), this%pse(j)%c(b)%se(:, LO))
                              bp_layer(:, HI) = min(bp_layer(:, HI), this%pse(j)%c(b)%se(:, HI))
                              b_layer(:,:) = bp_layer(:,:) - poff(:,:)
                              g = g + 1
                              do ib = 1, dom%nb
                                 cg%i_bnd(d, ib)%seg(g)%proc = j
                                 cg%i_bnd(d, ib)%seg(g)%se(:,LO) = b_layer(:, LO)
                                 cg%i_bnd(d, ib)%seg(g)%se(:,HI) = b_layer(:, HI)
                                 if (any(cg%i_bnd(d, ib)%seg(g)%se(d, :) < cg%lhn(d, LO))) &
                                      cg%i_bnd(d, ib)%seg(g)%se(d, :) = cg%i_bnd(d, ib)%seg(g)%se(d, :) + this%n_d(d)
                                 if (any(cg%i_bnd(d, ib)%seg(g)%se(d, :) > cg%lhn(d, HI))) &
                                      cg%i_bnd(d, ib)%seg(g)%se(d, :) = cg%i_bnd(d, ib)%seg(g)%se(d, :) - this%n_d(d)

                                 ! expand to cover corners (requires separate MPI_Waitall for each direction)
                                 !! \warning edges and corners will be filled multiple times
                                 do dd = xdim, zdim
                                    if (dd /= d .and. dom%has_dir(dd)) then
                                       cg%i_bnd(d, ib)%seg(g)%se(dd, LO) = cg%i_bnd(d, ib)%seg(g)%se(dd, LO) - ib
                                       cg%i_bnd(d, ib)%seg(g)%se(dd, HI) = cg%i_bnd(d, ib)%seg(g)%se(dd, HI) + ib
                                    endif
                                 enddo
                                 cg%o_bnd(d, ib)%seg(g) = cg%i_bnd(d, ib)%seg(g)
                                 cg%i_bnd(d, ib)%seg(g)%tag = int(HI*ndims*b          + (HI*d+lh-LO), kind=4) ! Assume that we won't mix communication with different ib
                                 cg%o_bnd(d, ib)%seg(g)%tag = int(HI*ndims*cg%grid_id + (HI*d+hl-LO), kind=4)
                                 select case (lh)
                                    case (LO)
                                       cg%i_bnd(d, ib)%seg(g)%se(d, LO) = cg%i_bnd(d, ib)%seg(g)%se(d, HI) - (ib - 1)
                                       cg%o_bnd(d, ib)%seg(g)%se(d, LO) = cg%i_bnd(d, ib)%seg(g)%se(d, HI) + 1
                                       cg%o_bnd(d, ib)%seg(g)%se(d, HI) = cg%o_bnd(d, ib)%seg(g)%se(d, LO) + (ib - 1)
                                    case (HI)
                                       cg%i_bnd(d, ib)%seg(g)%se(d, HI) = cg%i_bnd(d, ib)%seg(g)%se(d, LO) + (ib - 1)
                                       cg%o_bnd(d, ib)%seg(g)%se(d, HI) = cg%i_bnd(d, ib)%seg(g)%se(d, LO) - 1
                                       cg%o_bnd(d, ib)%seg(g)%se(d, LO) = cg%o_bnd(d, ib)%seg(g)%se(d, HI) - (ib - 1)
                                 end select

                              enddo
                           endif
                        enddo
                     enddo
                  endif
               enddo

               ! Detect fine-coarse boundaries and update boundary types

               b_layer(:,:) = cg%my_se(:, :)
               b_layer(d, HI) = b_layer(d, LO)
               n_tot_face_cells = product( b_layer(:, HI) - b_layer(:, LO) + 1 )
               allocate(facemap(b_layer(xdim,LO):b_layer(xdim,HI), b_layer(ydim,LO):b_layer(ydim,HI), b_layer(zdim,LO):b_layer(zdim,HI), LO:HI))

               facemap = .false.
               n_lbnd_face_cells = 0
               do g = lbound(cg%o_bnd(d, I_ONE)%seg, dim=1), ubound(cg%o_bnd(d, I_ONE)%seg, dim=1)
                  bp_layer(:, LO) = max(cg%o_bnd(d, I_ONE)%seg(g)%se(:, LO), cg%my_se(:, LO))
                  bp_layer(:, HI) = min(cg%o_bnd(d, I_ONE)%seg(g)%se(:, HI), cg%my_se(:, HI))
                  lh = LO
                  if (cg%o_bnd(d, I_ONE)%seg(g)%se(d, HI) > cg%ijkse(d, LO)) lh = HI
                  bp_layer(d, :) = b_layer(d, :)
                  facemap(bp_layer(xdim,LO):bp_layer(xdim,HI), bp_layer(ydim,LO):bp_layer(ydim,HI), bp_layer(zdim,LO):bp_layer(zdim,HI), lh) = .true.
               enddo
               do g = LO, HI
                  n_lbnd_face_cells(g) = count(facemap(:, :, :, g))
               enddo
               where (.not. cg%ext_bnd(d, :))
                  where (n_lbnd_face_cells(:) <  n_tot_face_cells) cg%bnd(d, :) = BND_MPI_FC
                  where (n_lbnd_face_cells(:) == 0)                cg%bnd(d, :) = BND_FC
               endwhere
               deallocate(facemap)
            endif
         enddo

         cgl => cgl%nxt
      enddo

   end subroutine mpi_bnd_types

!> \brief Add a whole level to the list of patches on current refinement level and decompose it

   subroutine add_patch_fulllevel(this, n_pieces)

      implicit none

      class(cg_level_T), target, intent(inout) :: this     !< current level
      integer(kind=4), optional, intent(in)    :: n_pieces !< how many pieces the patch should be divided to?

      call this%add_patch_detailed(this%n_d, this%off, n_pieces)

   end subroutine add_patch_fulllevel

!> \brief Add a new piece of grid to the list of patches on current refinement level and decompose it

   subroutine add_patch_detailed(this, n_d, off, n_pieces)

      use constants,     only: ndims
      use dataio_pub,    only: msg, die
      use decomposition, only: box_T

      implicit none

      class(cg_level_T), target,         intent(inout) :: this     !< current level
      integer(kind=8), dimension(ndims), intent(in)    :: n_d      !< number of grid cells
      integer(kind=8), dimension(ndims), intent(in)    :: off      !< offset (with respect to the base level, counted on own level)
      integer(kind=4), optional,         intent(in)    :: n_pieces !< how many pieces the patch should be divided to?

      type(box_T), dimension(:), allocatable :: tmp
      integer :: i

      if (.not. allocated(this%patches)) then
         allocate(this%patches(1))
      else
         allocate(tmp(lbound(this%patches(:),dim=1):ubound(this%patches(:), dim=1) + 1))
         tmp(:ubound(this%patches(:), dim=1)) = this%patches(:)
         ! manually deallocate arrays inside user-types, as it seems that move_alloc is unable to do that
         do i = lbound(this%patches(:), dim=1), ubound(this%patches(:), dim=1)
            if (allocated(this%patches(i)%pse)) deallocate(this%patches(i)%pse)
         enddo
         call move_alloc(from=tmp, to=this%patches)
      endif

      if (.not. this%patches(ubound(this%patches(:), dim=1))%decompose_patch(n_d(:), off(:), n_pieces)) then
         write(msg,'(a,i4)')"[cg_level:add_patch_detailed] Decomposition failed at level ",this%level_id
         call die(msg)
      endif

   end subroutine add_patch_detailed

end module cg_level
