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
   use cg_list,           only: cg_list_T   ! QA_WARN intel
#endif /* __INTEL_COMPILER */
   use cg_list_neighbors, only: cg_list_neighbors_T

   implicit none

   private
   public :: cg_level_T

   !>
   !! \brief A list of all cg of the same resolution.
   !!
   !! \details For positive refinement levels the list may be composed of several disconnected subsets of cg
   !! (islands: made of one or more cg's).
   !! This type is not intended for direct use. It is extended in cg_level_connected into a functional object.
   !<
   type, extends(cg_list_neighbors_T), abstract :: cg_level_T

      integer                                    :: tot_se       !< global number of segments on the level
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
      procedure, private :: update_gse                                           !< Gather updated information about the level and overwrite it to this%dot%gse
      procedure          :: update_tot_se                                        !< count all cg on current level for computing tags in vertical_prep
      generic,   public  :: add_patch => add_patch_fulllevel, add_patch_detailed !< Add a new piece of grid to the current level and decompose it
      procedure, private :: add_patch_fulllevel                                  !< Add a whole level to the list of patches
      procedure, private :: add_patch_detailed                                   !< Add a new piece of grid to the list of patches
      procedure          :: deallocate_patches                                   !< Throw out patches list
      procedure, private :: update_everything                                    !< Update all information on refinement structure and intra-level communication
      procedure          :: balance_old                                          !< Wrapper for rebalance_old

   end type cg_level_T

contains

!> \brief deallocate arrays

   subroutine cleanup(this)

      implicit none

      class(cg_level_T), intent(inout) :: this !< object invoking type bound procedure

      call this%plist%p_deallocate
      call this%dot%cleanup
      if (allocated(this%omega_mean))   deallocate(this%omega_mean)
      if (allocated(this%omega_cr))     deallocate(this%omega_cr)
      if (allocated(this%nshift))       deallocate(this%nshift)
      if (allocated(this%local_omega))  deallocate(this%local_omega)
      if (allocated(this%SFC_id_range)) deallocate(this%SFC_id_range)

   end subroutine cleanup

!>
!! \brief Throw out patches list
!!
!! \todo check if it is really necessary to call this from remote routines
!<

   subroutine deallocate_patches(this)

      implicit none

      class(cg_level_T), intent(inout) :: this !< object invoking type bound procedure

      call this%plist%p_deallocate

   end subroutine deallocate_patches

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
!!$      if (i /= this%cnt .or. this%cnt /= size(this%dot%gse(proc)%c(:)) .or. size(this%dot%gse(proc)%c(:)) /= i) then
!!$         write(msg, '(2(a,i4),a,3i7)')"[cg_level:print_segments] Uncertain number of grid pieces @PE ",proc," on level ", this%level_id, &
!!$              &                       " : ",i,this%cnt,size(this%dot%gse(proc)%c(:))
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
         tot_cg = tot_cg + size(this%dot%gse(p)%c(:))
         do i = lbound(this%dot%gse(p)%c(:), dim=1), ubound(this%dot%gse(p)%c(:), dim=1)
            ccnt = product(this%dot%gse(p)%c(i)%se(:, HI) - this%dot%gse(p)%c(i)%se(:, LO) + 1)
            maxcnt(p) = maxcnt(p) + ccnt
#ifdef VERBOSE
            if (i == 1) then
               write(header, '(a,i4,a,i3)')"[cg_level:print_segments] segment @", p, " ^", this%level_id
               hl = len_trim(header)
            else
               header = repeat(" ", hl)
            endif
            if (maxval(this%n_d(:)) < 1000000) then
               write(msg,'(2a,2(3i7,a),i8,a)') header(:hl), " : [", this%dot%gse(p)%c(i)%se(:, LO), "] : [", this%dot%gse(p)%c(i)%se(:, HI), "] #", ccnt, " cells"
            else if (maxval(this%n_d(:)) < 1000000000) then
               write(msg,'(2a,2(3i10,a),i8,a)') header(:hl), " : [", this%dot%gse(p)%c(i)%se(:, LO), "] : [", this%dot%gse(p)%c(i)%se(:, HI), "] #", ccnt, " cells"
            else
               write(msg,'(2a,2(3i18,a),i8,a)') header(:hl), " : [", this%dot%gse(p)%c(i)%se(:, LO), "] : [", this%dot%gse(p)%c(i)%se(:, HI), "] #", ccnt, " cells"
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

   subroutine init_all_new_cg(this, prevent_rebalancing)

      implicit none

      class(cg_level_T), intent(inout) :: this                !< object invoking type bound procedure
      logical, optional, intent(in)    :: prevent_rebalancing !< if present and .true. then do not allow rebalancing during addition of new grids

      ! First: do the balancing of new grids
      call this%balance_new(prevent_rebalancing)

      ! Second: create new grids, invalidate most of this%dot%gse database
      call this%create

      ! Third: update all information on refinement structure and intra-level communication, update this%dot%gse database
      ! Remember that the communication between levels has to be updated as well, but we cannot do this here due to cyclic dependencies
      call this%update_everything

   end subroutine init_all_new_cg

!> \brief update all information on refinement structure and intra-level communication.

   subroutine update_everything(this)

      implicit none

      class(cg_level_T), intent(inout) :: this   !< object invoking type bound procedure

      call this%update_decomposition_properties
      call this%update_gse     ! communicate everything that was added before
      call this%find_neighbors ! requires access to whole this%dot%gse(:)%c(:)%se(:,:)
      call this%update_req     ! Perhaps this%find_neighbors added some new entries
      call this%update_tot_se
      call this%print_segments

   end subroutine update_everything

!>
!! \brief Gather information on cg's currently present on local level, and write new this%dot%gse array
!!
!! OPT: For strict SFC distribution it is possible to determine complete list of neighbors (on the same level and
!! also one level up and down) and exchange only that data. It might be a bit faster for massively parallel runs.
!<

   subroutine update_gse(this)

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
      ! Use size(this%dot%gse(proc)%c) if you want to propagate gse before the grid containers are actually added to the level
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
         if (.not. associated(cgl)) call die("[cg_level:update_gse] Run out of cg.")
         allse(ncub*i      +1:ncub*i+   ndims) = int(cgl%cg%my_se(:, LO), kind=4)
         allse(ncub*i+ndims+1:ncub*i+HI*ndims) = int(cgl%cg%my_se(:, HI), kind=4) ! we do it in low-level way here. Is it worth using reshape() or something?
         if (any(cgl%cg%my_se > huge(allse(1)))) call die("[cg_level:update_gse] Implement 8-byte integers in MPI transactions for such huge refinements")
         cgl => cgl%nxt
      enddo
      if (associated(cgl)) call die("[cg_level:update_gse] Not all cg were read.")

      ! First use of MPI_Allgatherv in the Piernik Code!
      ncub_allcnt(:) = int(ncub * allcnt(:), kind=4)
      ncub_alloff(:) = int(ncub * alloff(:), kind=4)
      call MPI_Allgatherv(MPI_IN_PLACE, I_ZERO, MPI_DATATYPE_NULL, allse, ncub_allcnt, ncub_alloff, MPI_INTEGER, comm, mpi_err)

      ! Rewrite the gse array, forget about past.
      if (.not. allocated(this%dot%gse)) allocate(this%dot%gse(FIRST:LAST))
      do p = FIRST, LAST
         if (allocated(this%dot%gse(p)%c)) deallocate(this%dot%gse(p)%c)
         allocate(this%dot%gse(p)%c(allcnt(p)))
         do i = alloff(p), alloff(p) + allcnt(p) - 1
            this%dot%gse(p)%c(i-alloff(p)+1)%se(:, LO) = allse(ncub*i      +1:ncub*i+   ndims) ! we do it in low-level way here again.
            this%dot%gse(p)%c(i-alloff(p)+1)%se(:, HI) = allse(ncub*i+ndims+1:ncub*i+HI*ndims)
         enddo
      enddo

   end subroutine update_gse

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
      use cg_list,            only: cg_list_element
      use grid_cont,          only: grid_container
      use grid_container_ext, only: cg_extptrs
      use dataio_pub,         only: die
      use mpisetup,           only: FIRST, LAST, proc

      implicit none

      class(cg_level_T), intent(inout)        :: this   !< object invoking type bound procedure

      integer                                 :: i, p, ep
      integer(kind=8)                         :: s
      integer(kind=8)                         :: pieces
      type(cg_list_element), pointer          :: cgl
      logical                                 :: found_id
      type(grid_container), pointer           :: cg

      ! Find how many pieces are to be added
      pieces = this%plist%p_count()

      ! recreate local gse in case anything was derefined, refresh grid_id and make room for new pieces in the gse array
      if (.not. allocated(this%dot%gse)) allocate(this%dot%gse(FIRST:LAST))
      if (allocated(this%dot%gse(proc)%c)) deallocate(this%dot%gse(proc)%c)
      allocate(this%dot%gse(proc)%c(this%cnt + pieces))
      i = 0
      cgl => this%first
      do while (associated(cgl))
         i = i + 1
         this%dot%gse(proc)%c(i)%se(:,:) = cgl%cg%my_se(:,:)
         cgl%cg%grid_id = i
         cgl => cgl%nxt
      enddo
      do while (i<ubound(this%dot%gse(proc)%c, dim=1))
         i = i + 1
         this%dot%gse(proc)%c(i)%se = -huge(1)
      enddo

      ! check local consistency
      cgl => this%first
      do while (associated(cgl))
         found_id = .false.
         do i = lbound(this%dot%gse(proc)%c, dim=1), lbound(this%dot%gse(proc)%c, dim=1) + this%cnt - 1
            if (all(this%dot%gse(proc)%c(i)%se(:,:) == cgl%cg%my_se(:,:))) then
               if (found_id) call die("[cg_level:create] multiple occurrences")
               found_id = .true.
            endif
         enddo
         if (.not. found_id) call die("[cg_level:create] no occurrences")
         cgl => cgl%nxt
      enddo

      ! write the new grid pieces description to the gse array
      i = this%cnt
      if (allocated(this%plist%patches)) then
         do p = lbound(this%plist%patches(:), dim=1), ubound(this%plist%patches(:), dim=1)
            do s = lbound(this%plist%patches(p)%pse, dim=1), ubound(this%plist%patches(p)%pse, dim=1)
               i = i + 1
               if (i > size(this%dot%gse(proc)%c(:))) call die("[cg_level:create] overflow")
               this%dot%gse(proc)%c(i)%se(:,:) = this%plist%patches(p)%pse(s)%se(:,:)
               call this%add
               cg => this%last%cg
               call cg%init(this%n_d, this%off, this%dot%gse(proc)%c(i)%se(:, :), i, this%level_id) ! we cannot pass "this" as an argument because of circular dependencies
               do ep = lbound(cg_extptrs%ext, dim=1), ubound(cg_extptrs%ext, dim=1)
                  if (associated(cg_extptrs%ext(ep)%init))  call cg_extptrs%ext(ep)%init(cg)
               enddo
               call all_cg%add(cg)
            enddo
         enddo
         call this%plist%p_deallocate
      endif

   end subroutine create

!> \brief Update some flags in domain module [ is_uneven, is_mpi_noncart, is_refined, is_multicg ]

   subroutine update_decomposition_properties(this)

      use cg_list,    only: cg_list_element
      use constants,  only: base_level_id, pLOR,  pLAND, ndims, I_ONE
      use dataio_pub, only: warn
      use domain,     only: is_mpi_noncart, is_multicg, is_refined, is_uneven
      use mpi,        only: MPI_INTEGER, MPI_REQUEST_NULL
      use mpisetup,   only: proc, master, slave, piernik_MPI_Allreduce, proc, req, status, comm, mpi_err, LAST, inflate_req

      implicit none

      logical, save :: warned = .false.

      class(cg_level_T), intent(inout) :: this   !< object invoking type bound procedure
      type(cg_list_element), pointer :: cgl
      integer(kind=4), dimension(ndims) :: shape, shape1
      integer(kind=4), parameter :: sh_tag = 7
      integer, parameter :: nr = 2

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

      is_multicg = is_multicg .or. (ubound(this%dot%gse(proc)%c(:), dim=1) > 1)
      call piernik_MPI_Allreduce(is_multicg, pLOR)

      ! check if all blocks in the domain have same size and shape
      call inflate_req(nr)
      this%is_blocky = .true.
      shape = 0
      shape1 = 0
      cgl => this%first
      if (associated(cgl)) then
         shape = cgl%cg%n_b
         cgl => cgl%nxt
         do while (associated(cgl))
            if (any(cgl%cg%n_b /= shape)) this%is_blocky = .false.
            cgl => cgl%nxt
         enddo
      endif
      req = MPI_REQUEST_NULL
      if (slave)     call MPI_Irecv(shape1, size(shape1), MPI_INTEGER, proc-I_ONE, sh_tag, comm, req(1), mpi_err)
      if (proc<LAST) call MPI_Isend(shape,  size(shape),  MPI_INTEGER, proc+I_ONE, sh_tag, comm, req(nr), mpi_err)
      call MPI_Waitall(nr, req(:nr), status(:, :nr), mpi_err)
      if (any(shape /= 0) .and. any(shape1 /= 0)) then
         if (any(shape /= shape1)) this%is_blocky = .false.
      endif
      call piernik_MPI_Allreduce(this%is_blocky, pLAND)

   end subroutine update_decomposition_properties

!> \brief Count all cg on current level. Useful for computing tags in vertical_prep

   subroutine update_tot_se(this)

      use mpisetup, only: FIRST, LAST

      implicit none

      class(cg_level_T), intent(inout) :: this   !< object invoking type bound procedure

      integer :: p

      this%tot_se = 0
      do p = FIRST, LAST
         if (allocated(this%dot%gse)) this%tot_se = this%tot_se + ubound(this%dot%gse(p)%c(:), dim=1)
      enddo

   end subroutine update_tot_se

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
      call this%plist%expand
      if (.not. this%plist%patches(ubound(this%plist%patches(:), dim=1))%decompose_patch(n_d(:), off(:), this%level_id, n_pieces=n_pieces)) then
         write(msg,'(a,i4)')"[cg_level:add_patch_detailed] Decomposition failed at level ",this%level_id
         call die(msg)
      endif

   end subroutine add_patch_detailed

!> \brief Wrapper for cg_list_rebalance%balance_old

   subroutine balance_old(this)

      use constants, only: pLOR
      use mpisetup,  only: piernik_MPI_Allreduce

      implicit none

      class(cg_level_T), intent(inout) :: this

      call this%rebalance_old
      ! OPT: call this%update_gse inside this%update_everything can be quite long to complete
      call piernik_MPI_Allreduce(this%recently_changed, pLOR)
      if (this%recently_changed) call this%update_everything

   end subroutine balance_old

end module cg_level
