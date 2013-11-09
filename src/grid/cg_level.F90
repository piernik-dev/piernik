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
   !! \details For positive refinement levels the list may be composed of several disconnected subsets of cg (islands: made of one or more cg's).
   !! This type is not intended for direct use. It is extended in cg_level_connected into a functional object.
   !!
   !! OPT: Searching through this%pse for neighbours, prolongation/restriction overlaps etc is quite costly - O(this%cnt^2)
   !! Provide a list, sorted according to Morton/Hilbert id's and do a bisection search instead of checking against all grids
   !! It will result in massive speedups on cg_level_T%mpi_bnd_types and cg_level_connected_T%{vertical_prep,vertical_b_prep).
   !! It may also simplify the process of fixing refinement structure in refinement_update::fix_refinement.
   !! Grids which are larger than AMR_bsize (merged grids, non-block decompositions, both not implemented yet) may be referred by several id's that correspond
   !! with AMR_bsize-d virtual grid pieces.
   !!
   !! Alternatively, construct a searchable binary tree or oct-tree and provide fast routines for searching grid pieces covering specified position.
   !!
   !! \todo Provide one of the structures described above
   !<
   type, extends(cg_list_bnd_T) :: cg_level_T

      integer(kind=4)                            :: level_id         !< level number (relative to base level). No arithmetic should depend on it.
      integer(kind=8), dimension(ndims)          :: n_d              !< maximum number of grid cells in each direction (size of fully occupied level)
      type(cuboids),   dimension(:), allocatable :: pse              !< lists of grid chunks on each process (FIRST:LAST); Use with care, because this is an antiparallel thing
      integer                                    :: tot_se           !< global number of segments on the level
      integer                                    :: fft_type         !< type of FFT to employ in some multigrid solvers (depending on boundaries)
      type(box_T),     dimension(:), allocatable :: patches          !< list of patches that exist on the current level
      integer(kind=8), dimension(ndims)          :: off              !< offset of the level
      logical                                    :: recently_changed !< .true. when anything was added to or deleted from this level

      ! FARGO
      real,    dimension(:, :), allocatable      :: omega_mean       !< mean angular velocity for each fluid
      real,    dimension(:, :), allocatable      :: omega_cr         !< constant residual angular velocity for each fluid
      integer, dimension(:, :), allocatable      :: nshift           !< number of cells that need to be shifted due to %omega_mean for each fluid
      real,    dimension(:, :), allocatable      :: local_omega      !< auxiliary array
      integer(kind=8), dimension(:), allocatable :: cell_count       !< auxiliary counter

    contains

      procedure          :: init_all_new_cg                                      !< initialize newest grid container
      procedure, private :: mpi_bnd_types                                        !< create MPI types for boundary exchanges
      procedure          :: print_segments                                       !< print detailed information about current level decomposition
      procedure, private :: update_decomposition_properties                      !< Update some flags in domain module
      procedure, private :: create                                               !< Get all decomposed patches and turn them into local grid containers
      procedure, private :: mark_new                                             !< Detect which grid containers are new
      procedure, private :: update_pse                                           !< Gather updated information about the level and overwrite it to this%pse
      procedure          :: update_tot_se                                        !< count all cg on current level for computing tags in vertical_prep
      generic,   public  :: add_patch => add_patch_fulllevel, add_patch_detailed !< Add a new piece of grid to the current level and decompose it
      procedure, private :: add_patch_fulllevel                                  !< Add a whole level to the list of patches
      procedure, private :: add_patch_detailed                                   !< Add a new piece of grid to the list of patches
      procedure, private :: add_patch_one_piece                                  !< Add a patch with only one grid piece
      procedure, private :: expand_list                                          !< Expand the patch list by one
      procedure, private :: balance_new                                          !< Routine for moving proposed grids between processes
      procedure          :: balance_old                                          !< Routine for measuring disorder level in distribution of grids across processes
      procedure, private :: reshuffle                                            !< Routine for moving existing grids between processes
   end type cg_level_T

contains

!> \brief Print detailed information about current level decomposition

   subroutine print_segments(this)

      use cg_list,    only: cg_list_element
      use constants,  only: LO, HI
      use dataio_pub, only: printinfo !, msg, warn
      use mpisetup,   only: FIRST, LAST, master!, nproc, proc
#ifdef VERBOSE
      use dataio_pub, only: msg
#endif /* VERBOSE */

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
!!$      if (i /= this%cnt .or. this%cnt /= size(this%pse(proc)%c(:)) .or. size(this%pse(proc)%c(:)) /= i) then
!!$         write(msg, '(2(a,i4),a,3i7)')"[cg_level:print_segments] Uncertain number of grid pieces @PE ",proc," on level ", this%level_id, &
!!$              &                       " : ",i,this%cnt,size(this%pse(proc)%c(:))
!!$         call warn(msg)
!!$
!!$         cgl => this%first
!!$         do while (associated(cgl))
!!$            write(msg,'(2(a,i7),2(a,3i10),a)')" @",proc," #",cgl%cg%grid_id," : [", cgl%cg%my_se(:, LO), "] : [", cgl%cg%my_se(:, HI)," ]"
!!$            call printinfo(msg)
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
         tot_cg = tot_cg + size(this%pse(p)%c(:))
         do i = lbound(this%pse(p)%c(:), dim=1), ubound(this%pse(p)%c(:), dim=1)
            ccnt = product(this%pse(p)%c(i)%se(:, HI) - this%pse(p)%c(i)%se(:, LO) + 1)
            maxcnt(p) = maxcnt(p) + ccnt
#ifdef VERBOSE
            if (i == 1) then
               write(header, '(a,i4,a,i3)')"[cg_level:print_segments] segment @", p, " ^", this%level_id
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

!      write(msg, '(a,i3,a,f5.1,a,i5,a,f8.5)')"[cg_level:print_segments] Level ", this%level_id, " filled in ",(100.*sum(maxcnt(:)))/product(real(this%n_d(:))), &
!           &                                 "%, ",tot_cg," grid(s), load balance : ", sum(maxcnt(:))/(nproc*maxval(maxcnt(:)))
      !> \todo add calculation of total internal boundary surface in cells
!      call printinfo(msg)
      deallocate(maxcnt)

   end subroutine print_segments

!>
!! \brief Initialize all grid containers on a new grid level
!!
!! \todo automagically rebalance existing grids unless it is explicitly forbidden
!<

   subroutine init_all_new_cg(this)

      implicit none

      class(cg_level_T), intent(inout) :: this   !< object invoking type bound procedure

      ! First: do the balancing of new grids, update this%pse database
      call this%balance_new

      ! Second: create new grids
      call this%create

      ! Third: update all information on refinement structure and intra-level communication.
      ! Remember that the communication between levels has to be updated as well, but we cannot do this here due to cyclic dependencies
      call this%update_decomposition_properties
      call this%update_pse    ! communicate everything that was added above
      call this%mpi_bnd_types ! require access to whole this%pse(:)%c(:)%se(:,:)
      call this%update_req    ! Perhaps this%mpi_bnd_types added some new entries
      call this%update_tot_se
      call this%print_segments

   end subroutine init_all_new_cg

!> \brief Gather information on cg's currently present on local level, and write new this%pse array

   subroutine update_pse(this)

      use cg_list,    only: cg_list_element
      use constants,  only: I_ZERO, I_ONE, ndims, LO, HI
      use dataio_pub, only: die
      use mpi,        only: MPI_IN_PLACE, MPI_DATATYPE_NULL, MPI_INTEGER
      use mpisetup,   only: FIRST, LAST, proc, comm, mpi_err !!$, master

      implicit none

      class(cg_level_T), intent(inout) :: this   !< object invoking type bound procedure

      integer(kind=4), dimension(FIRST:LAST) :: allcnt, alloff, ncub_allcnt, ncub_alloff
      integer(kind=4), allocatable, dimension(:) :: allse
      type(cg_list_element), pointer :: cgl
      integer :: i, p
      integer, parameter :: ncub = ndims*HI ! the number of integers in each cuboid

      ! get the count of grid pieces on each process
      ! Beware: int(this%cnt, kind=4) is not properly updated after calling this%distribute.
      ! Use size(this%pse(proc)%c) if you want to propagate pse before the grid containers are actually added to the level
      ! OPT: this call can be quite long to complete
      call MPI_Allgather(int(this%cnt, kind=4), I_ONE, MPI_INTEGER, allcnt, I_ONE, MPI_INTEGER, comm, mpi_err)

      ! compute offsets for  a composite table of all grid pieces
      alloff(FIRST) = I_ZERO
      do i = FIRST+I_ONE, LAST
         alloff(i) = alloff(i-1) + allcnt(i-1)
      enddo

      allocate(allse(ncub*sum(allcnt(:))))

      ! Collect definitions of own grid pieces
      allse = 0
      cgl => this%first
      do i = alloff(proc), alloff(proc) + allcnt(proc) - 1
         if (.not. associated(cgl)) call die("[cg_level:update_pse] Run out of cg.")
         allse(ncub*i      +1:ncub*i+   ndims) = int(cgl%cg%my_se(:, LO), kind=4)
         allse(ncub*i+ndims+1:ncub*i+HI*ndims) = int(cgl%cg%my_se(:, HI), kind=4) ! we do it in low-level way here. Is it worth using reshape() or something?
         if (any(cgl%cg%my_se > huge(allse(1)))) call die("[cg_level:update_pse] Implement 8-byter integers in MPI transactions for such huge refinements")
         cgl => cgl%nxt
      enddo
      if (associated(cgl)) call die("[cg_level:update_pse] Not all cg were read.")

      ! First use of MPI_Allgatherv in the Piernik Code!
      ncub_allcnt(:) = int(ncub * allcnt(:), kind=4)
      ncub_alloff(:) = int(ncub * alloff(:), kind=4)
      call MPI_Allgatherv(MPI_IN_PLACE, I_ZERO, MPI_DATATYPE_NULL, allse, ncub_allcnt, ncub_alloff, MPI_INTEGER, comm, mpi_err)

      ! Rewrite the pse array, forget about past.
      if (.not. allocated(this%pse)) allocate(this%pse(FIRST:LAST))
      do p = FIRST, LAST
         if (allocated(this%pse(p)%c)) deallocate(this%pse(p)%c)
         allocate(this%pse(p)%c(allcnt(p)))
         do i = alloff(p), alloff(p) + allcnt(p) - 1
            this%pse(p)%c(i-alloff(p)+1)%se(:, LO) = allse(ncub*i      +1:ncub*i+   ndims) ! we do it in low-level way here again.
            this%pse(p)%c(i-alloff(p)+1)%se(:, HI) = allse(ncub*i+ndims+1:ncub*i+HI*ndims)
         enddo
      enddo

   end subroutine update_pse

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
!! \brief Get all decomposed patches and turn them into local grid containers
!!
!! \details This routine starts with two lists:
!! * A list of blocks that survived derefinement attempts
!! * A list of blocks to be created due to requested refinement
!! It has to decide if and how to do migration of grid pieces and to communicate this update to global database of grid pieces.
!!
!! \deprecated: I have an impression that the most challenging work was moved to balance_new routine
!!
!!
!! First strategy will be implemented first to get everything working. Second strategy will be used quite often. Third one do not need to be used on every refinement update.
!! It can be called when some benchmark of grid disorder exceeds particular threshold.
!<

   subroutine create(this)

      use cg_list_global,     only: all_cg
      use cg_list,            only: cg_list_element
      use constants,          only: LONG
      use grid_cont,          only: grid_container
      use grid_container_ext, only: cg_extptrs
      use dataio_pub,         only: die
      use mpisetup,           only: FIRST, LAST, proc

      implicit none

      class(cg_level_T), intent(inout)        :: this   !< object invoking type bound procedure

      integer                                 :: i, p
      integer(kind=8)                         :: s
      integer(kind=8)                         :: pieces, filled
      type(cg_list_element), pointer          :: cgl
      type(cuboid), allocatable, dimension(:) :: new_c
      logical                                 :: found_id
      type(grid_container), pointer           :: cg

      call this%update_pse  ! required if anything was derefined

      ! Find how many pieces are to be added
      pieces = 0
      if (allocated(this%patches)) then
         do p = lbound(this%patches(:), dim=1), ubound(this%patches(:), dim=1)
            if (allocated(this%patches(p)%pse)) pieces = pieces + size(this%patches(p)%pse, dim=1)
         enddo
      endif

      ! make room for new pieces in the pse array
      filled = 0
      if (.not. allocated(this%pse)) allocate(this%pse(FIRST:LAST))
      if (allocated(this%pse(proc)%c)) then
         filled = size(this%pse(proc)%c)
         allocate(new_c(pieces + filled))
         new_c(:filled) = this%pse(proc)%c
         do s = filled+1_LONG, pieces + filled
            new_c(s)%se = -huge(1)
         enddo
         call move_alloc(from=new_c, to=this%pse(proc)%c)
      endif

      ! fix grid_id
      cgl => this%first
      do while (associated(cgl))
         found_id = .false.
         do i = lbound(this%pse(proc)%c, dim=1), ubound(this%pse(proc)%c, dim=1)
            if (all(this%pse(proc)%c(i)%se(:,:) == cgl%cg%my_se(:,:))) then
               if (found_id) call die("[cg_level:create] multiple occurrences")
               cgl%cg%grid_id = i
               found_id = .true.
            endif
         enddo
         if (.not. found_id) call die("[cg_level:create] no occurrences")
         cgl => cgl%nxt
      enddo

      ! write the new grid pieces description to the pse array
      if (allocated(this%patches)) then
         do p = lbound(this%patches(:), dim=1), ubound(this%patches(:), dim=1)
            do s = lbound(this%patches(p)%pse, dim=1), ubound(this%patches(p)%pse, dim=1)
               filled = filled + 1
               if (filled > size(this%pse(proc)%c(:))) call die("[cg_level:create] overflow")
               this%pse(proc)%c(filled)%se(:,:) = this%patches(p)%pse(s)%se(:,:)
            enddo
         enddo
         deallocate(this%patches)
      endif

      call this%mark_new

      do i = lbound(this%pse(proc)%c(:), dim=1), ubound(this%pse(proc)%c(:), dim=1)
         if (this%pse(proc)%c(i)%is_new) then
            this%pse(proc)%c(i)%is_new = .false.
            call this%add
            cg => this%last%cg
            call cg%init(this%n_d, this%off, this%pse(proc)%c(i)%se(:, :), i, this%level_id) ! we cannot pass "this" as an argument because of circular dependencies
            do p = lbound(cg_extptrs%ext, dim=1), ubound(cg_extptrs%ext, dim=1)
               if (associated(cg_extptrs%ext(p)%init))  call cg_extptrs%ext(p)%init(cg)
            enddo
            call all_cg%add(cg)
         endif
      enddo

   end subroutine create

!> \brief Update some flags in domain module [ is_uneven, is_mpi_noncart, is_refined, is_multicg ]

   subroutine update_decomposition_properties(this)

      use constants,  only: base_level_id, pLOR
      use dataio_pub, only: warn
      use domain,     only: is_mpi_noncart, is_multicg, is_refined, is_uneven
      use mpisetup,   only: proc, master, piernik_MPI_Allreduce

      implicit none

      logical, save :: warned = .false.

      class(cg_level_T), intent(inout) :: this   !< object invoking type bound procedure

      if (this%level_id > base_level_id) is_refined = .true.
      call piernik_MPI_Allreduce(is_refined, pLOR)
      if (is_refined) then
         is_mpi_noncart = .true.
         is_multicg = .true.
         if (master .and. .not. warned) then
            call warn("[cg_level:update_decomposition_properties] Refinements are experimental")
            warned = .true.
         endif
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
!! \details this type can be a member of grid container type if we pass this%pse(:) as an argument.
!! It would simplify dependencies and this%init_all_new_cg, but it could be quite a big object.
!! Assume that cuboids don't collide (no overlapping grid pieces on same refinement level are allowed)
!!
!! Current implementation (revision 7338) implies correct update of all corners, even on complicated refinement topologies (concave fine region - convect coarse region or
!! fine regions touching each other only by corners). Previous implementation could correctly fill the corners only on uniform grid and when it was called for
!! x, y and z directions separately. Warning: that change introduces measurable performance degradation! This is caused by the fact that in 3D it is required to make
!! 26 separate exchanges to fill all guardcells (in cg_list_bnd::internal_boundaries), while in previous approach only 6 exchanges were required.
!! Unfortunately the previous approach did not work properly for complicated refinements.
!!
!! Possible improvements of performance
!! * do local exchanges directly, without calling MPI.
!! * merge smaller blocks into larger ones,
!!
!! \todo Put this%pse into a separate type and pass a pointer to it or even a pointer to pre-filtered segment list
!!
!! \todo Rewrite this routine to achieve previous (pre-7338) performance and maintain correctness on corners on complicated topologies:
!! * Divide the descriptions of communicated regions into 4 categories: X-faces, Y-faces + XY-corners, Z-faces + [XY]Z-corners, other corners.
!!   The other corners would be non-empty only for some refinement local topologies, it would certainly be empty on an uniform grid.
!! * When no corners are required, perform simultaneous exchange described by the three directional categories. Some corners might be set up correctly by a chance,
!!   some might not.
!! * When corners are required, perform sequential exchange described by the three directional categories and supplement it with communication of "other corners".
!!   The sequence of Isend/Irecv should be as follows: Isend X-faces, Irecv X-faces, Waitall, Isend Y-faces, Irecv Y-faces, Waitall, Isend Z-faces, Irecv Z-faces, Waitall
!!   "Other corners" can be Isend at any time and must be Irecv after Z-faces are copied to the right place.
!<

   subroutine mpi_bnd_types(this)

      use cg_list,    only: cg_list_element
      use constants,  only: xdim, ydim, zdim, cor_dim, LO, HI, BND_MPI_FC, BND_FC
      use domain,     only: dom
      use grid_cont,  only: grid_container, is_overlap
      use mpisetup,   only: FIRST, LAST

      implicit none

      class(cg_level_T), intent(inout)                :: this    !< object invoking type bound procedure

      type(grid_container),  pointer                  :: cg      !< grid container that we are currently working on
      type(cg_list_element), pointer                  :: cgl
      integer                                         :: j, b, id, ix, iy, iz
      integer(kind=8)                                 :: n_lbnd_face_cells
      integer(kind=4)                                 :: d, dd, hl, lh, tag
      integer(kind=8), dimension(xdim:zdim)           :: per
      integer(kind=8), dimension(xdim:zdim, LO:HI)    :: b_layer, poff, aux
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
         allocate(cg%i_bnd(xdim:cor_dim), cg%o_bnd(xdim:cor_dim))

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
         allocate(cg%i_bnd(cor_dim)%seg(0), cg%o_bnd(cor_dim)%seg(0))

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
                                                aux = this%pse(j)%c(b)%se
                                                aux(:, LO) = aux(:, LO) + [ ix, iy, iz ] * per(:)
                                                aux(:, HI) = aux(:, HI) + [ ix, iy, iz ] * per(:)
                                                tag = uniq_tag(cg%my_se, aux, b)
                                                aux = poff
                                                aux(d, :) = aux(d, :) + [ -1, 1 ]
                                                if (is_overlap(cg%my_se, aux)) then
                                                   dd = d
                                                else
                                                   dd = cor_dim
                                                endif
                                                call cg%i_bnd(dd)%add_seg(j, poff, tag)
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
                                                aux = cg%my_se
                                                aux(:, LO) = aux(:, LO) - [ ix, iy, iz ] * per(:)
                                                aux(:, HI) = aux(:, HI) - [ ix, iy, iz ] * per(:)
                                                tag = uniq_tag(this%pse(j)%c(b)%se, aux, cg%grid_id)
                                                aux = poff
                                                aux(:, LO) = aux(:, LO) - [ ix, iy, iz ] * per(:)
                                                aux(:, HI) = aux(:, HI) - [ ix, iy, iz ] * per(:)
                                                aux(d, :) = aux(d, :) + [ -1, 1 ]
                                                if (is_overlap(this%pse(j)%c(b)%se, aux)) then
                                                   dd = d
                                                else
                                                   dd = cor_dim
                                                endif
                                                call cg%o_bnd(dd)%add_seg(j, poff, tag)
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

   contains

      !>
      !! \brief Create unique tag for cg - cg exchange
      !!
      !! \details If we put a constraint that a grid piece can not be smaller than dom%nb, then total number of neighbours that affect local guardcells is
      !! * 3 for cartesian decomposition (including AMR with equal-size blocks)
      !! * 2 to 4 for noncartesian decomposition
      !! * more for AMR with consolidated blocks (unimplemented yet, not compatible with current approach)
      !! Thus, in each direction we can describe realtive position as one of four cases, or a bit easier one of five cases:
      !! * FAR_LEFT, FAR_RIGHT - corner neighbours, either touching corner or a bit further away
      !! * LEFT, RIGHT - partially face, partially corner neighbours
      !! * FACE - face neighbour (may cover also some corners)
      !<
      pure function uniq_tag(se, nb_se, grid_id)

         use constants, only: LO, HI, xdim, ydim, zdim, INVALID

         implicit none

         integer(kind=8), dimension(xdim:zdim, LO:HI), intent(in) :: se       ! a grid piece
         integer(kind=8), dimension(xdim:zdim, LO:HI), intent(in) :: nb_se    ! neighboring grid piece
         integer,                                      intent(in) :: grid_id  ! grid piece id

         integer(kind=4) :: uniq_tag
         integer, dimension(xdim:zdim) :: r
         integer :: d
         enum, bind(C)
            enumerator :: FAR_LEFT=0, LEFT, FACE, RIGHT, FAR_RIGHT, N_POS
         end enum

         r = INVALID
         do d = xdim, zdim
            if (nb_se(d, LO) > se(d, HI)) then
               r(d) = FAR_RIGHT
            else if (nb_se(d, HI) < se(d, LO)) then
               r(d) = FAR_LEFT
            else if ((nb_se(d, LO) < se(d, LO)) .and. (nb_se(d, HI) < se(d, HI)) .and. (nb_se(d, HI) >= se(d, LO))) then
               r(d) = LEFT
            else if ((nb_se(d, HI) > se(d, HI)) .and. (nb_se(d, LO) > se(d, LO)) .and. (nb_se(d, LO) <= se(d, HI))) then
               r(d) = RIGHT
            else
               r(d) = FACE
            endif
         enddo
         uniq_tag = int(((grid_id*N_POS+r(zdim))*N_POS+r(ydim))*N_POS+r(xdim), kind=4)

      end function uniq_tag

   end subroutine mpi_bnd_types

!> \brief Add a whole level to the list of patches on current refinement level and decompose it.

   subroutine add_patch_fulllevel(this, n_pieces)

      implicit none

      class(cg_level_T), target, intent(inout) :: this     !< current level
      integer(kind=4), optional, intent(in)    :: n_pieces !< how many pieces the patch should be divided to?

      call this%add_patch_detailed(this%n_d, this%off, n_pieces)

   end subroutine add_patch_fulllevel

!>
!! \brief Add a new piece of grid to the list of patches on current refinement level and decompose it
!!
!! \details Each patch should be added only once, i.e. the base level and multigrid coarse grids have to be added by the master and then redistributed in balance_new.
!<

   subroutine add_patch_detailed(this, n_d, off, n_pieces)

      use constants,     only: ndims
      use dataio_pub,    only: msg, die

      implicit none

      class(cg_level_T), target,         intent(inout) :: this     !< current level
      integer(kind=8), dimension(ndims), intent(in)    :: n_d      !< number of grid cells
      integer(kind=8), dimension(ndims), intent(in)    :: off      !< offset (with respect to the base level, counted on own level)
      integer(kind=4), optional,         intent(in)    :: n_pieces !< how many pieces the patch should be divided to?

      this%recently_changed = .true. ! assume that the new patches will change this level
      call this%expand_list
      if (.not. this%patches(ubound(this%patches(:), dim=1))%decompose_patch(n_d(:), off(:), this%level_id, n_pieces=n_pieces)) then
         write(msg,'(a,i4)')"[cg_level:add_patch_detailed] Decomposition failed at level ",this%level_id
         call die(msg)
      endif

   end subroutine add_patch_detailed

!> \brief Add a patch with only one grid piece

   subroutine add_patch_one_piece(this, n_d, off)

      use constants,     only: ndims

      implicit none

      class(cg_level_T), target,         intent(inout) :: this     !< current level
      integer(kind=8), dimension(ndims), intent(in)    :: n_d      !< number of grid cells
      integer(kind=8), dimension(ndims), intent(in)    :: off      !< offset (with respect to the base level, counted on own level)

      this%recently_changed = .true. ! assume that the new patches will change this level
      call this%expand_list
      call this%patches(ubound(this%patches(:), dim=1))%one_piece_patch(n_d(:), off(:))

   end subroutine add_patch_one_piece

!> \brief Expand the patch list by one

   subroutine expand_list(this)

      use decomposition, only: box_T

      implicit none

      class(cg_level_T), target,         intent(inout) :: this     !< current level

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

   end subroutine expand_list

!>
!! \brief Routine for moving proposed grids between processes
!!
!! \details Starts with list of already allocated blocks and list of patches which are not yet turned into blocks on a given level
!! Move the patches between processes to maintain best possible work balance.
!!
!! First, all planned patches are gathered in an array on the master process, and deallocated locally.
!! Then, the patches are sorted according to Space-Filling Curve index and distributed among the prosesses.
!! The processes with least workload will get more patches.
!! After the distribution most processes should have roughly equal number of patches (+/- 1) with the possible exception
!! of few processes that were initially heavily loaded.
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
!!
!! There are several strategies than can be implemented:
!! * Local refinements go to local process. It is very simple, but for most simulations will build up load imbalance. Suitable for tests and global refinement.
!! * Local refinements can be assigned to remote processes, existing blocks stays in place. Should keep good load balance, but the amount of inter-process
!!   internal boundaries may grow significantly with time. Suitable for minor refinement updates and base level decomposition. This is the current implementation.
!! * All blocks (existing and new) have recalculated assignment and can be migrated to other processes. Most advanced. Should be used after reading restart data.
!<

   subroutine balance_new(this)

      use constants,       only: pSUM, ndims, INVALID, LO, HI, I_ONE
      use dataio_pub,      only: die
      use mpi,             only: MPI_INTEGER, MPI_INTEGER8, MPI_STATUS_IGNORE!, MPI_REQUEST_NULL
      use mpisetup,        only: piernik_MPI_Allreduce, master, FIRST, LAST, comm, req, mpi_err, status, nproc, inflate_req
      use sort_piece_list, only: grid_piece_list

      implicit none

      class(cg_level_T), intent(inout) :: this

      type(grid_piece_list) :: gp
      integer :: i
      integer(kind=4), dimension(FIRST:LAST+1) :: from
      integer(kind=4), dimension(FIRST:LAST) :: cnt_existing
      integer(kind=4) :: ls, p, s
      integer(kind=4), parameter :: tag_ls = 1, tag_gpt = tag_ls+1, tag_lsR = tag_gpt+1, tag_gptR = tag_lsR+1
      enum, bind(C)
         enumerator :: I_OFF
         enumerator :: I_N_B = I_OFF + ndims
         enumerator :: I_END = I_N_B + ndims - I_ONE
      end enum
      integer(kind=8), dimension(:,:), allocatable :: gptemp
      integer, parameter :: nreq = 1

      ! collect planned grids on a level

      call inflate_req(nreq)

      ! count how many patches were requested on each process
      s = 0
      if (allocated(this%patches)) then
         do p = lbound(this%patches(:), dim=1, kind=4), ubound(this%patches(:), dim=1, kind=4)
            s = s + size(this%patches(p)%pse, dim=1, kind=4)
         enddo
      endif
      ls = int(s, kind=4)
      call piernik_MPI_Allreduce(s, pSUM) !> \warning overkill: MPI_reduce is enough here

      if (s==0) return ! nihil novi

      ! copy the patches data to a temporary array to be sent to the master
      allocate(gptemp(I_OFF:I_END,ls))
      if (master) then
         call gp%init(s)
      else
         call MPI_Isend(ls, I_ONE, MPI_INTEGER, FIRST, tag_ls, comm, req(nreq), mpi_err)
      endif
      call MPI_Gather(this%cnt, I_ONE, MPI_INTEGER, cnt_existing, I_ONE, MPI_INTEGER, FIRST, comm, mpi_err)
      i = 0
      if (allocated(this%patches)) then
         do p = lbound(this%patches(:), dim=1, kind=4), ubound(this%patches(:), dim=1, kind=4)
            do s = lbound(this%patches(p)%pse, dim=1, kind=4), ubound(this%patches(p)%pse, dim=1, kind=4)
               i = i + 1
               gptemp(:, i) = [ this%patches(p)%pse(s)%se(:, LO), this%patches(p)%pse(s)%se(:, HI) - this%patches(p)%pse(s)%se(:, LO) + 1 ]
            enddo
         enddo
      endif
      if (allocated(this%patches)) deallocate(this%patches)

      if (master) then !> \warning Antiparallel

         ! put all the patches (own and obtained from slaves) on a list gp%list
         do s = 1, ls
            call gp%list(s)%set_gp(gptemp(I_OFF:I_OFF+ndims-1, s), int(gptemp(I_N_B:I_N_B+ndims-1, s), kind=4), INVALID, FIRST)
         enddo
         i = ls
         deallocate(gptemp)
         do p = FIRST + 1, LAST
            call MPI_Recv(ls, I_ONE, MPI_INTEGER, p, tag_ls, comm, MPI_STATUS_IGNORE, mpi_err)
            if (ls > 0) then
               allocate(gptemp(I_OFF:I_END,ls))
               call MPI_Recv(gptemp, size(gptemp), MPI_INTEGER8, p, tag_gpt, comm, MPI_STATUS_IGNORE, mpi_err)
               do s = 1, ls
                  call gp%list(i+s)%set_gp(gptemp(I_OFF:I_OFF+ndims-1, s), int(gptemp(I_N_B:I_N_B+ndims-1, s), kind=4), INVALID, p)
               enddo
               i = i + ls
               deallocate(gptemp)
            endif
         enddo

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
         if (i /= 0) call die("[cg_level:balance_new] i /= 0")
         do p = FIRST, LAST
            from(p+1) = from(p+1) + from(p)
         enddo
         do p = from(FIRST), from(FIRST+1) - I_ONE
            call this%add_patch_one_piece(int(gp%list(p)%n_b, kind=8), gp%list(p)%off)
         enddo

         ! distribute proposed grids according to limits computed above
         do p = FIRST + I_ONE, LAST
            ls = int(from(p+1) - from(p), kind=4)
            ! call MPI_Isend(ls, I_ONE, MPI_INTEGER, p, tag_lsR, comm, req(p), mpi_err) !can't reuse ls before MPI_Waitall
            call MPI_Send(ls, I_ONE, MPI_INTEGER, p, tag_lsR, comm, mpi_err)
            if (ls>0) then
               allocate(gptemp(I_OFF:I_END,from(p):from(p+1)-1))
               do s = lbound(gptemp, dim=2, kind=4), ubound(gptemp, dim=2, kind=4)
                  gptemp(:, s) = [ gp%list(s)%off, int(gp%list(s)%n_b, kind=8) ]
               enddo
               ! call MPI_Isend(gptemp, size(gptemp), MPI_INTEGER8, p, tag_gptR, comm, req(LAST+p), mpi_err) !can't deallocate gptemp before MPI_Waitall
               call MPI_Send(gptemp, size(gptemp), MPI_INTEGER8, p, tag_gptR, comm, mpi_err)
               deallocate(gptemp)
            else
               ! req(LAST+p) = MPI_REQUEST_NULL
            endif
         enddo
         ! call MPI_Waitall(2*LAST, req(:2*LAST), status(:,:2*LAST), mpi_err)
      else

         ! send patches to master
         call MPI_Wait(req(nreq), status(:, nreq), mpi_err)
         if (ls > 0) call MPI_Send(gptemp, size(gptemp), MPI_INTEGER8, FIRST, tag_gpt, comm, mpi_err)
         deallocate(gptemp)

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

      if (master) call gp%cleanup

!!$      allocate(area(1:lmax))
!!$      area = 0
!!$      do i = lbound(cg_res(:), dim=1), ubound(cg_res(:), dim=1)
!!$         if (cg_res(i)%level >= 1) area(cg_res(i)%level) = area(cg_res(i)%level) + product(cg_res(i)%n_b)
!!$      enddo
!!$      deallocate(area)

   end subroutine balance_new

!> \brief Routine for measuring disorder level in distribution of grids across processes

!#define DEBUG
   subroutine balance_old(this)

      use cg_list,         only: cg_list_element
      use cg_list_dataop,  only: expanded_domain
      use constants,       only: ndims, LO, HI, I_ONE, pSUM
      use dataio_pub,      only: warn, msg, printinfo
      use mpisetup,        only: master, FIRST, LAST, nproc, piernik_MPI_Bcast, piernik_MPI_Allreduce
      use refinement,      only: oop_thr
      use sort_piece_list, only: grid_piece_list
#ifdef DEBUG
      use mpi,             only: MPI_INTEGER, MPI_INTEGER8, MPI_STATUS_IGNORE
      use mpisetup,        only: comm, mpi_err
#endif /* DEBUG */

      implicit none

      class(cg_level_T), intent(inout) :: this

      integer(kind=4), dimension(FIRST:LAST) :: cnt_existing
      type(grid_piece_list) :: gp
      type(cg_list_element), pointer :: cgl
      integer :: s
      integer(kind=4) :: i, p
      integer(kind=8), dimension(:,:), allocatable :: gptemp
      enum, bind(C)
         enumerator :: I_OFF
         enumerator :: I_N_B = I_OFF + ndims
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
      ! Gather complete grid list and compare with this%pse
      call MPI_Gather(this%cnt, I_ONE, MPI_INTEGER, cnt_existing, I_ONE, MPI_INTEGER, FIRST, comm, mpi_err)
      if (master) then
         call gp%init(sum(cnt_existing))
         do i = I_ONE, this%cnt
            call gp%list(i)%set_gp(gptemp(I_OFF:I_OFF+ndims-1, i), int(gptemp(I_N_B:I_N_B+ndims-1, i), kind=4), int(gptemp(I_GID, i), kind=4), FIRST)
            if (any(this%pse(FIRST)%c(i)%se(:, LO) /= gp%list(i)%off) .or. gp%list(i)%cur_gid /= i .or. &
                 any(this%pse(FIRST)%c(i)%se(:, HI) - this%pse(FIRST)%c(i)%se(:, LO) +1 /= gp%list(i)%n_b)) &
                 call warn("cl:bo this%pse(FIRST) /= gptemp")
         enddo
         deallocate(gptemp)
         s = this%cnt
         do p = FIRST + 1, LAST
            if (cnt_existing(p) > 0) then
               allocate(gptemp(I_OFF:I_GID, cnt_existing(p)))
               call MPI_Recv(gptemp, size(gptemp), MPI_INTEGER8, p, tag_gpt, comm, MPI_STATUS_IGNORE, mpi_err)
               do i = I_ONE, cnt_existing(p)
                  call gp%list(i+s)%set_gp(gptemp(I_OFF:I_OFF+ndims-1, i), int(gptemp(I_N_B:I_N_B+ndims-1, i), kind=4), int(gptemp(I_GID, i), kind=4), p)
                  if (any(this%pse(p)%c(i)%se(:, LO) /= gp%list(i+s)%off) .or. gp%list(i+s)%cur_gid /= i .or. &
                       any(this%pse(p)%c(i)%se(:, HI) - this%pse(p)%c(i)%se(:, LO) +1 /= gp%list(i+s)%n_b)) &
                       call warn("cl:bo this%pse(p) /= gptemp")
               enddo
               s = s + cnt_existing(p)
               deallocate(gptemp)
            endif
         enddo
      else
         if (this%cnt > 0) call MPI_Send(gptemp, size(gptemp), MPI_INTEGER8, FIRST, tag_gpt, comm, mpi_err)
      endif
#else /* !DEBUG */
      ! Trust that this%pse is updated
      if (master) then
         do p = FIRST, LAST
            cnt_existing(p) = size(this%pse(p)%c, kind=4)
         enddo
         call gp%init(sum(cnt_existing))
         i = 0
         do p = FIRST, LAST
            ii = i
            if (ii /= sum(cnt_existing(:p-1))) call warn("cl:bo ii /= sum(cnt_existing(:p-1))")
            do s = lbound(this%pse(p)%c, dim=1), ubound(this%pse(p)%c, dim=1)
               i = i + I_ONE
               call gp%list(i)%set_gp(this%pse(p)%c(s)%se(:, LO), int(this%pse(p)%c(s)%se(:, HI) - this%pse(p)%c(s)%se(:, LO) +1, kind=4), i - ii, p)
            enddo
         enddo
      endif
#endif /* DEBUG */
      if (allocated(gptemp)) deallocate(gptemp)

      if (master) then
         call gp%set_id(this%off)
         call gp%sort
         do p = FIRST, LAST
            gp%list(p*size(gp%list)/nproc+1:(p+1)*size(gp%list)/nproc)%dest_proc = p
         enddo
         s = 0
         if (size(gp%list) > 0) then
            s = count(gp%list(:)%cur_proc /= gp%list(:)%dest_proc)
            if (s/real(size(gp%list)) > oop_thr) then
               write(msg,'(a,i3,2(a,i6),a,f6.3,a)')"[cg_level:balance_old] ^", this%level_id," Reshuffling OutOfPlace grids:",s, "/",size(gp%list)," (load balance: ",sum(cnt_existing)/real(maxval(cnt_existing)*size(cnt_existing)),")"
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
            write(msg,'(a,i5,a)')"[cg_level:balance_old] Allreduce(expanded_domain%cnt) = ",p,", aborting reshuffling."
            if (master) call warn(msg)
         else
            call this%reshuffle(gp)
         endif
      endif

      if (master) call gp%cleanup

   end subroutine balance_old

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

      class(cg_level_T),     intent(inout) :: this
      type(grid_piece_list), intent(in)    :: gp

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
         if (gptemp(I_C_P, i) == gptemp(I_D_P, i)) call die("[cg_level:balance_old] can not send to self")
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
            if (.not. found) call die("[cg_level:balance_old] Grid id not found")
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

      !> \deprecated partially copied code from init_all_new_cg
      ! OPT: call this%update_pse can be quite long to complete
      call this%update_pse
      call this%mpi_bnd_types ! require access to whole this%pse(:)%c(:)%se(:,:)
      call this%update_req    ! Perhaps this%mpi_bnd_types added some new entries
      call this%update_tot_se
      call this%print_segments

      deallocate(gptemp)

   end subroutine reshuffle

end module cg_level
