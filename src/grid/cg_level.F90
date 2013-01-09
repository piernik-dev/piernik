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
!!
!! Current implementation (revision 7338) implies correct update of the corners, even on complicated refinement topologies (concave fine region - convect coarse region or
!! fine regions touching each other only by corners). Previous implementation could correctly fill the corners only on uniform grid and when it was called for
!! x, y and z dierctions separately. Warning: that change introduces measurable performance degradation! This is caused by the fact that in 3D it is required to make
!! 26 separate exchanges to fill all guardcells (in cg_list_bnd::internal_boundaries), while in previous approach only 6 exchanges were required.
!! Unfortunately the previous approach did not work properly for complicated refinements and I saw no easy solution for that.
!! Possible strategies to achieve previous (or better) performance
!! * provide separate communication pattern for calls where no corners are required,
!! * do local exchanges directly, without calling MPI.
!! * merge smaller blocks into larger ones,
!! * implement totally noncartesian domain decomposition for uniform grid (with 12 neighbours instead of 26),
!! * go back to directionally-split list and carefully supplement them by pieces that can not be correctly updated due to refinement topology.
!!
!! \todo Put this%pse into a separate type and pass a pointer to it or even a pointer to pre-filtered segment list
!!
!! \todo Do not provide segments for each possible number of guardcells. Provide 1 layer, 1 layer with corners, all layers and all layers with cofrners instead.
!<

   subroutine mpi_bnd_types(this)

      use cg_list,    only: cg_list_element
      use constants,  only: xdim, ydim, zdim, LO, HI, BND_MPI_FC, BND_FC
      use domain,     only: dom
      use grid_cont,  only: grid_container, is_overlap
      use mpisetup,   only: FIRST, LAST

      implicit none

      class(cg_level_T), intent(inout)                :: this    !< object invoking type bound procedure

      type(grid_container),  pointer                  :: cg      !< grid container that we are currently working on
      type(cg_list_element), pointer                  :: cgl
      integer                                         :: j, b, id, ix, iy, iz
      integer(kind=8)                                 :: n_lbnd_face_cells
      integer(kind=4)                                 :: d, hl, lh, tag
      integer(kind=8), dimension(xdim:zdim)           :: per
      integer(kind=8), dimension(xdim:zdim, LO:HI)    :: b_layer, poff
      type :: fmap
         logical, dimension(:,:,:), allocatable       :: map
         integer(kind=8), dimension(xdim:zdim, LO:HI) :: b_layer
      end type fmap
      type(fmap), dimension(xdim:zdim, LO:HI)         :: f
      integer(kind=8), dimension(ndims, LO:HI)        :: box_8   !< temporary storage

      cgl => this%first
      do while (associated(cgl))
         cg => cgl%cg

         if (allocated(cg%i_bnd)) deallocate(cg%i_bnd)
         if (allocated(cg%o_bnd)) deallocate(cg%o_bnd)
         allocate(cg%i_bnd(xdim:zdim), cg%o_bnd(xdim:zdim))

         per(:) = 0
         where (dom%periodic(:)) per(:) = this%n_d(:)

         ! Create maps to mark neighbouring face cells
         do d = xdim, zdim
            if (.not. allocated(cg%i_bnd(d)%seg)) allocate(cg%i_bnd(d)%seg(0))
            if (.not. allocated(cg%o_bnd(d)%seg)) allocate(cg%o_bnd(d)%seg(0))
            if (dom%has_dir(d)) then
               do lh = LO, HI
                  hl = LO+HI-lh ! HI for LO, LO for HI
                  f(d, lh)%b_layer(:,:) = cg%my_se(:, :)
                  f(d, lh)%b_layer(d, hl) = f(d, lh)%b_layer(d, lh)          ! interior cell layer, 1 cell thick, without corners
                  allocate(f(d, lh)%map(f(d, lh)%b_layer(xdim,LO):f(d, lh)%b_layer(xdim,HI), &
                       &                f(d, lh)%b_layer(ydim,LO):f(d, lh)%b_layer(ydim,HI), &
                       &                f(d, lh)%b_layer(zdim,LO):f(d, lh)%b_layer(zdim,HI)))
                  f(d, lh)%map = .false.
               enddo
            endif
         enddo

         do j = FIRST, LAST
            do b = lbound(this%pse(j)%c(:), dim=1), ubound(this%pse(j)%c(:), dim=1)
               box_8 = int(cg%lhn, kind=8)
               if (is_overlap(box_8, this%pse(j)%c(b)%se(:,:), per(:))) then                ! identify processes with interesting neighbour data

                  do d = xdim, zdim
                     if (dom%has_dir(d)) then
                        do lh = LO, HI
                           hl = LO+HI-lh ! HI for LO, LO for HI

                           ! First, update the map of faces with neighbours

                           ! create 1-layer thick map of neighbours
                           b_layer = this%pse(j)%c(b)%se
                           b_layer(d, hl) = b_layer(d, hl) - lh+hl  ! move the opposite boundary
                           b_layer(d, lh) = b_layer(d, hl)

                           do id = -1, 1 ! scan through periodic images of the domain
                              if (id == 0 .or. per(d)>0) then
                                 poff = b_layer
                                 poff(d, :) = poff(d, :) + id*per(d)
                                 poff(:, LO) = max(poff(:, LO), f(d, lh)%b_layer(:, LO))
                                 poff(:, HI) = min(poff(:, HI), f(d, lh)%b_layer(:, HI))
                                 ! construct the layer to be send to the _interior_ of neighbouring grid and set the flag map
                                 if (is_overlap(f(d, lh)%b_layer, poff)) &
                                      f(d, lh)%map(poff(xdim,LO):poff(xdim,HI), poff(ydim,LO):poff(ydim,HI), poff(zdim,LO):poff(zdim,HI)) = .true.
                              endif
                           enddo

                           ! Second, describe incoming data
                           b_layer = cg%my_se
                           b_layer(d, hl) = b_layer(d, lh) + (lh-hl)
                           b_layer(d, lh) = b_layer(d, lh) + (lh-hl)*dom%nb ! dom%nb thick layer without corners
                           b_layer(:d-1, LO) = b_layer(:d-1, LO) - dom%nb*dom%D_(:d-1)
                           b_layer(:d-1, HI) = b_layer(:d-1, HI) + dom%nb*dom%D_(:d-1) ! corners added in only one way
                           ! faces and corners are included in y and z direction to minimize number of pieces in non-cartesian grid decompositions

                           ! set up segments to be received
                           do iz = -1, 1 ! scan through all periodic possibilities
                              if (iz == 0 .or. per(zdim)>0) then
                                 do iy = -1, 1
                                    if (iy == 0 .or. per(ydim)>0) then
                                       do ix = -1, 1
                                          if (ix == 0 .or. per(xdim)>0) then
                                             poff = this%pse(j)%c(b)%se
                                             poff(:, LO) = poff(:, LO) + [ ix, iy, iz ] * per(:)
                                             poff(:, HI) = poff(:, HI) + [ ix, iy, iz ] * per(:)
                                             if (is_overlap(b_layer, poff)) then
                                                poff(:, LO) = max(b_layer(:, LO), poff(:, LO))
                                                poff(:, HI) = min(b_layer(:, HI), poff(:, HI))
                                                tag = int(27*(HI*ndims*b + (HI*d+lh-LO))+ixyz(ix, iy, iz), kind=4)
                                                call cg%i_bnd(d)%add_seg(j, poff, tag)
                                             endif
                                          endif
                                       enddo
                                    endif
                                 enddo
                              endif
                           enddo

                           ! Third, describe outgoing data
                           !> \warning replicated code, see above
                           b_layer = this%pse(j)%c(b)%se
                           b_layer(d, hl) = b_layer(d, lh) + (lh-hl)
                           b_layer(d, lh) = b_layer(d, lh) + (lh-hl)*dom%nb ! dom%nb thick layer without corners
                           b_layer(:d-1, LO) = b_layer(:d-1, LO) - dom%nb*dom%D_(:d-1)
                           b_layer(:d-1, HI) = b_layer(:d-1, HI) + dom%nb*dom%D_(:d-1) ! corners added in only one way
                           ! faces and corners are included in y and z direction to minimize number of pieces in non-cartesian grid decompositions

                           ! set up segments to be send
                           do iz = -1, 1 ! scan through all periodic possibilities
                              if (iz == 0 .or. per(zdim)>0) then
                                 do iy = -1, 1
                                    if (iy == 0 .or. per(ydim)>0) then
                                       do ix = -1, 1
                                          if (ix == 0 .or. per(xdim)>0) then
                                             poff = b_layer
                                             poff(:, LO) = poff(:, LO) + [ ix, iy, iz ] * per(:)
                                             poff(:, HI) = poff(:, HI) + [ ix, iy, iz ] * per(:)
                                             if (is_overlap(poff(:,:), cg%my_se)) then
                                                poff(:, LO) = max(cg%my_se(:, LO), poff(:, LO))
                                                poff(:, HI) = min(cg%my_se(:, HI), poff(:, HI))
                                                tag = int(27*(HI*ndims*cg%grid_id + (HI*d+lh-LO))+ixyz(-ix, -iy, -iz), kind=4)
                                                call cg%o_bnd(d)%add_seg(j, poff, tag)
                                             endif
                                          endif
                                       enddo
                                    endif
                                 enddo
                              endif
                           enddo

                        enddo
                     endif
                  enddo

               endif

            enddo
         enddo

         ! Detect fine-coarse boundaries and update boundary types. When not all mapped cells are facing neighbours, then we may deal with fine/coarse boundary (full or partial)
         do d = xdim, zdim
            if (dom%has_dir(d)) then
               do lh = LO, HI
                  n_lbnd_face_cells = count(f(d,lh)%map(:,:,:))
                  if (.not. cg%ext_bnd(d, lh)) then
                     if (n_lbnd_face_cells < size(f(d,lh)%map(:,:,:))) cg%bnd(d, lh) = BND_MPI_FC
                     if (n_lbnd_face_cells == 0)                       cg%bnd(d, lh) = BND_FC
                  endif
                  deallocate(f(d,lh)%map)
               enddo
            endif
         enddo

         cgl => cgl%nxt
      enddo

   end subroutine mpi_bnd_types

   ! Create unique number from periodic offsets
   ! Let's assume ix, iy and iz have values of -1, 0 or 1, then the output will belong to [0 .. 26] range
   pure function ixyz(i1, i2, i3)

      implicit none

      integer, intent(in) :: i1, i2, i3

      integer :: ixyz

      ixyz = (1+i1) + 3*(1+i2 + 3*(1+i3))

   end function ixyz

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
