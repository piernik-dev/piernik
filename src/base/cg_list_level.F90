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

!> \brief This module contains grid container level list and related methods

module cg_list_lev

   use constants, only: ndims
   use gc_list,   only: cg_list
   use domain,    only: cuboids

   implicit none

   private
   public :: cg_list_level, cg_list_patch

   !>
   !! \brief A list of all cg of the same resolution.
   !!
   !! \details For positive refinement levels the list may be composed of several disconnected subsets of cg ("islands: made of one or more cg: cg_list_patch).
   !<
   type, extends(cg_list) :: cg_list_level
      integer :: lev                           !< level number (relative to base level). For printing, debug, and I/O use only. No arithmetic should depend on it.
      integer(kind=8), dimension(ndims) :: n_d !< maximum number of grid cells in each direction (size of fully occupied level)
      type(cuboids), dimension(:), allocatable :: pse  !< lists of grid chunks on each process (FIRST:LAST); Use with care, because this is an antiparallel thing
      integer :: tot_se                        !< global number of segments on the level
      type(cg_list_level), pointer :: coarser  !< coarser level cg set or null()
      type(cg_list_level), pointer :: finer    !< finer level cg set or null()
    contains

      ! Level management
      procedure :: init_new_cg                 !< initialize newest grid container
      procedure :: mpi_bnd_types               !< create MPI types for boundary exchanges
      procedure :: print_segments              !< print detailed information about current level decomposition

      ! Prolongation and restriction
      procedure :: vertical_prep               !< initialize prolongation and restriction targets
      procedure :: update_tot_se               !< count all cg on current level for computing tags in vertical_prep
      procedure :: prolong0_q_1var              !< interpolate the grid data to this%finer level
      procedure :: restrict_q_1var              !< interpolate the grid data from this%coarser level
      procedure :: restrict_to_floor_q_1var                !< restrict as much as possible

      ! fine-coarse boundary exchanges may also belong to this type
   end type cg_list_level

   !>
   !! \brief A list of grid containers that cover single box (or rectangle) on a certain resolution level
   !!
   !! \details This set would be a result of base domain or patch decomposition
   !<
   type, extends(cg_list) :: cg_list_patch
      integer(kind=4), dimension(ndims) :: n_d                !< number of grid cells
      integer(kind=8), dimension(ndims) :: off                !< offset (with respect to the base level, counted on own level)
      type(cg_list_patch), pointer :: parent                  !< Parent patch (or null()). \todo Consider relaxing this restriction and allow multi-parent patches
      type(cg_list_patch), dimension(:), pointer :: children  !< refined patches
      type(cg_list_level), pointer :: list_level              !< all cg on the same level
      !> \todo consider creating neigbour list (or (ndims, LO:HI) lists)
   end type cg_list_patch

contains

!>
!! \brief Initialize prolongation and restriction targets. Called from init_multigrid.
!!
!! \details This is the most generic method of searching prolongation/restriction targets.
!! This routine does not assume any special sorting, so it may become costly at some point.
!!
!! \todo When we implement the Hilbert sorting, write a new routine that will take advantage of it and avoid using this one.
!! \todo Search parent/children patches only, should be faster than pairing between whole levels
!<

   subroutine vertical_prep(this)

      use constants,     only: xdim, ydim, zdim, LO, HI
      use dataio_pub,    only: die
      use gc_list,       only: cg_list_element
      use grid_cont,     only: pr_segment, grid_container, is_overlap
      use mpisetup,      only: proc, FIRST, LAST

      implicit none

      class(cg_list_level), intent(inout) :: this

      integer :: g, j, jf, fmax, tag
      integer(kind=8), dimension(xdim:zdim) :: ijks
      integer(kind=8), dimension(xdim:zdim, LO:HI) :: coarsened
      type(pr_segment), pointer :: seg
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg            !< current grid container
      type(cg_list_level),   pointer :: fine, coarse  !< shortcut
      type :: int_pair
         integer :: proc
         integer :: n_se
      end type int_pair
      type(int_pair), dimension(:), allocatable :: ps

      call this%update_tot_se

      fine => this%finer
      ! find fine target for receiving restricted data or sending data to be prolonged
      if (associated(fine)) then
         cgl => this%first
         do while (associated(cgl))
            cg => cgl%cg

            ijks(:) = cg%ijkse(:, LO) - cg%off(:)  ! add this to convert absolute cell coordinates to local indices. (+dom%nb - off(:))

            if (allocated(ps)) call die("cll:vp f a ps")
            fmax = 0
            do j = FIRST, LAST
               fmax = fmax + size(fine%pse(j)%sel(:, :, :), dim=1)
            enddo
            allocate(ps(fmax))

            if (allocated(cg%f_tgt%seg)) call die("cll:vp cg%f_tgt%seg a a")
            g = 0
            do j = FIRST, LAST
               do jf = lbound(fine%pse(j)%sel(:, :, :), dim=1), ubound(fine%pse(j)%sel(:, :, :), dim=1)
                  if (is_overlap(cg%my_se(:, :), fine%pse(j)%sel(jf, :, :)/2)) then
                     g = g + 1
                     ps(g) = int_pair(j, jf)
                  endif
               enddo
            enddo
            allocate(cg%f_tgt%seg(g))

            do g = lbound(cg%f_tgt%seg(:), dim=1), ubound(cg%f_tgt%seg(:), dim=1)
               seg => cg%f_tgt%seg(g)
               if (allocated(seg%buf)) call die("cll:vp f seg%buf a a")
               seg%proc = ps(g)%proc
               ! find cross-section of own segment with coarsened fine segment
               seg%se(:, LO) = max(cg%my_se(:, LO), fine%pse(seg%proc)%sel(ps(g)%n_se, :, LO)/2) + ijks(:)
               seg%se(:, HI) = min(cg%my_se(:, HI), fine%pse(seg%proc)%sel(ps(g)%n_se, :, HI)/2) + ijks(:)
               allocate(seg%buf(seg%se(xdim, HI)-seg%se(xdim, LO) + 1, &
                    &           seg%se(ydim, HI)-seg%se(ydim, LO) + 1, &
                    &           seg%se(zdim, HI)-seg%se(zdim, LO) + 1))
               tag = cg%grid_id + this%tot_se * ps(g)%n_se
               seg%tag = int(tag, kind=4) ! assumed that there is only one piece to be communicated from grid to grid (i.e. grids are not periodically wrapped around)
               if (tag /= int(seg%tag)) call die("[cg_list_level:vertical_prep] tag overflow (1)")
            enddo

            if (allocated(ps)) deallocate(ps)
            cgl => cgl%nxt
         enddo
      endif

      ! find coarse target for sending restricted data or receiving data to be prolonged
      !> \deprecated almost duplicated code
      coarse => this%coarser
      if (associated(coarse)) then

         call coarse%update_tot_se

         cgl => this%first
         do while (associated(cgl))
            cg => cgl%cg

            ijks(:) = cg%ijkse(:, LO) - cg%off(:)  ! add this to convert absolute cell coordinates to local indices. (+dom%nb - off(:))

            if (allocated(ps)) call die("cll:vp c a ps")
            fmax = 0
            do j = FIRST, LAST
               fmax = fmax + size(coarse%pse(j)%sel(:, :, :), dim=1)
            enddo
            allocate(ps(fmax))

            if (allocated(cg%c_tgt%seg)) call die("cll:vp cg%c_tgt%seg a a")
            g = 0
            coarsened(:,:) = cg%my_se(:, :)/2
            do j = FIRST, LAST
               do jf = lbound(coarse%pse(j)%sel(:, :, :), dim=1), ubound(coarse%pse(j)%sel(:, :, :), dim=1)
                  if (is_overlap(coarsened(:,:), coarse%pse(j)%sel(jf, :, :))) then
                     g = g + 1
                     ps(g) = int_pair(j, jf)
                  endif
               enddo
            enddo
            allocate(cg%c_tgt%seg(g))

            do g = lbound(cg%c_tgt%seg(:), dim=1), ubound(cg%c_tgt%seg(:), dim=1)
               seg => cg%c_tgt%seg(g)
               if (allocated(seg%buf)) call die("cll:vp c seg%buf a a")
               seg%proc = ps(g)%proc
               ! find cross-section of own segment with refined coarse segment
               seg%se(:, LO) = max(cg%my_se(:, LO), coarse%pse(seg%proc)%sel(ps(g)%n_se, :, LO)*2  )
               seg%se(:, HI) = min(cg%my_se(:, HI), coarse%pse(seg%proc)%sel(ps(g)%n_se, :, HI)*2+1)
               allocate(seg%buf(seg%se(xdim, HI)/2-seg%se(xdim, LO)/2 + 1, &
                    &           seg%se(ydim, HI)/2-seg%se(ydim, LO)/2 + 1, &
                    &           seg%se(zdim, HI)/2-seg%se(zdim, LO)/2 + 1))
               seg%se(:, LO) = seg%se(:, LO) + ijks(:)
               seg%se(:, HI) = seg%se(:, HI) + ijks(:)
               tag = ps(g)%n_se + coarse%tot_se * cg%grid_id
               seg%tag = int(tag, kind=4)
               if (tag /= int(seg%tag)) call die("[cg_list_level:vertical_prep] tag overflow (2)")
            enddo

            if (allocated(ps)) deallocate(ps)
            cgl => cgl%nxt
         enddo
      endif

   end subroutine vertical_prep

!> \brief Restrict as much as possible

   recursive subroutine restrict_to_floor_q_1var(this, iv)

      implicit none

      class(cg_list_level), target, intent(inout) :: this !< object invoking type-bound procedure
      integer,                      intent(in)    :: iv   !< variable to be restricted

      if (.not. associated(this%coarser)) return
      call this%restrict_q_1var(iv)
      call this%coarser%restrict_to_floor_q_1var(iv)

   end subroutine restrict_to_floor_q_1var

!>
!! \brief Simplest restriction (averaging).
!!
!! \todo implement high order restriction and test its influence on V-cycle convergence rate. Watch f/c boundaries.
!!
!! \details Some data can be locally copied without MPI, but this seems to have really little impact on the performance.
!! Some tests show that purely MPI code without local copies is marginally faster.
!<

   subroutine restrict_q_1var(this, iv)

      use constants,  only: xdim, ydim, zdim, LO, HI, LONG, I_ONE
      use dataio_pub, only: msg, warn!, die
      use domain,     only: dom
      use gc_list,    only: cg_list_element
      use grid_cont,  only: grid_container
      use mpisetup,   only: proc, comm, ierr, req, status
      use mpi,        only: MPI_DOUBLE_PRECISION

      implicit none

      class(cg_list_level), target, intent(inout) :: this !< object invoking type-bound procedure
      integer,                      intent(in)    :: iv   !< variable to be restricted

      type(cg_list_level), pointer :: coarse
      integer :: g
      integer(kind=8), dimension(:,:), pointer :: fse, cse ! shortcuts for fine segment and coarse segment
      integer(kind=8) :: i, j, k, ic, jc, kc
      integer(kind=8), dimension(xdim:zdim) :: off1
      real :: norm
      integer(kind=4) :: nr
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg            !< current grid container

      coarse => this%coarser
      if (.not. associated(coarse)) then ! can't restrict base level
         write(msg,'(a,i3)')"[cg_list_level:restrict_q_1var] no coarser level than ", this%lev
         call warn(msg)
         return
      endif

      ! be ready to receive everything into right buffers
      nr = 0
      cgl => coarse%first
      do while (associated(cgl))
         cg => cgl%cg
         if (allocated(cg%f_tgt%seg)) then
            do g = lbound(cg%f_tgt%seg(:), dim=1), ubound(cg%f_tgt%seg(:), dim=1)
               nr = nr + I_ONE
               call MPI_Irecv(cg%f_tgt%seg(g)%buf(1, 1, 1), size(cg%f_tgt%seg(g)%buf(:, :, :)), MPI_DOUBLE_PRECISION, cg%f_tgt%seg(g)%proc, cg%f_tgt%seg(g)%tag, comm, req(nr), ierr)
            enddo
         endif
         cgl => cgl%nxt
      enddo

      ! interpolate to coarse buffer and send it
      norm = 1./2**dom%eff_dim
      cgl => this%first
      do while (associated(cgl))
         cg => cgl%cg
         do g = lbound(cg%c_tgt%seg(:), dim=1), ubound(cg%c_tgt%seg(:), dim=1)

            fse => cg%c_tgt%seg(g)%se

            cg%c_tgt%seg(g)%buf(:, :, :) = 0.
            off1(:) = mod(cg%off(:)+fse(:, LO) - cg%ijkse(:, LO), 2_LONG)
            do k = fse(zdim, LO), fse(zdim, HI)
               kc = (k-fse(zdim, LO)+off1(zdim))/2 + 1
               do j = fse(ydim, LO), fse(ydim, HI)
                  jc = (j-fse(ydim, LO)+off1(ydim))/2 + 1
                  do i = fse(xdim, LO), fse(xdim, HI)
                     ic = (i-fse(xdim, LO)+off1(xdim))/2 + 1
                     cg%c_tgt%seg(g)%buf(ic, jc, kc) = cg%c_tgt%seg(g)%buf(ic, jc, kc) + cg%q(iv)%arr(i, j, k) * norm
                  enddo
               enddo
            enddo
            nr = nr + I_ONE
            call MPI_Isend(cg%c_tgt%seg(g)%buf(1, 1, 1), size(cg%c_tgt%seg(g)%buf(:, :, :)), MPI_DOUBLE_PRECISION, cg%c_tgt%seg(g)%proc, cg%c_tgt%seg(g)%tag, comm, req(nr), ierr)
         enddo
         cgl => cgl%nxt
      enddo

      if (nr>0) call MPI_Waitall(nr, req(:nr), status(:,:nr), ierr)

      ! copy the received buffers to the right places
      cgl => coarse%first
      do while (associated(cgl))
         cg => cgl%cg
         if (allocated(cg%f_tgt%seg)) then
            cg%q(iv)%arr(:,:,:) = 0. ! disables check_dirty
            do g = lbound(cg%f_tgt%seg(:), dim=1), ubound(cg%f_tgt%seg(:), dim=1)
               cse => cg%f_tgt%seg(g)%se
               cg%q(iv)%arr     (cse(xdim, LO):cse(xdim, HI), cse(ydim, LO):cse(ydim, HI), cse(zdim, LO):cse(zdim, HI)) = &
                    cg%q(iv)%arr(cse(xdim, LO):cse(xdim, HI), cse(ydim, LO):cse(ydim, HI), cse(zdim, LO):cse(zdim, HI)) + cg%f_tgt%seg(g)%buf(:, :, :)
            enddo
         endif
         cgl => cgl%nxt
      enddo

   end subroutine restrict_q_1var

!>
!! \brief 0th order prolongation : injection
!!
!! \todo implement high order prolongation. Watch f/c boundaries.
!<

   subroutine prolong0_q_1var(this, iv)

      use constants,  only: xdim, ydim, zdim, LO, HI, LONG, I_ONE
      use dataio_pub, only: msg, warn!, die
      use gc_list,    only: cg_list_element
      use grid_cont,  only: grid_container
      use mpisetup,   only: proc, comm, ierr, req, status
      use mpi,        only: MPI_DOUBLE_PRECISION

      implicit none

      class(cg_list_level), target, intent(inout) :: this !< object invoking type-bound procedure
      integer,                      intent(in)    :: iv   !< variable to be prolonged

      type(cg_list_level), pointer :: fine
      integer :: g
      integer(kind=8), dimension(:,:), pointer :: fse, cse ! shortcuts for fine segment and coarse segment
      integer(kind=8) :: i, j, k, ic, jc, kc
      integer(kind=8), dimension(xdim:zdim) :: off1
      integer(kind=4) :: nr
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg            !< current grid container

      fine => this%finer
      if (.not. associated(fine)) then ! can't prolong finest level
         write(msg,'(a,i3)')"[gc_list:restrict_q_1var] no finer level than: ", this%lev
         call warn(msg)
         return
      endif

      nr = 0
      ! be ready to receive everything into right buffers
      cgl => fine%first
      do while (associated(cgl))
         cg => cgl%cg
         if (allocated(cg%c_tgt%seg)) then
            do g = lbound(cg%c_tgt%seg(:), dim=1), ubound(cg%c_tgt%seg(:), dim=1)
               nr = nr + I_ONE
               call MPI_Irecv(cg%c_tgt%seg(g)%buf(1, 1, 1), size(cg%c_tgt%seg(g)%buf(:, :, :)), MPI_DOUBLE_PRECISION, cg%c_tgt%seg(g)%proc, cg%c_tgt%seg(g)%tag, comm, req(nr), ierr)
            enddo
         endif
         cgl => cgl%nxt
      enddo

      ! send coarse data
      cgl => this%first
      do while (associated(cgl))
         cg => cgl%cg
         do g = lbound(cg%f_tgt%seg(:), dim=1), ubound(cg%f_tgt%seg(:), dim=1)

            cse => cg%f_tgt%seg(g)%se

            nr = nr + I_ONE
            cg%f_tgt%seg(g)%buf(:, :, :) = cg%q(iv)%arr(cse(xdim, LO):cse(xdim, HI), cse(ydim, LO):cse(ydim, HI), cse(zdim, LO):cse(zdim, HI))
            call MPI_Isend(cg%f_tgt%seg(g)%buf(1, 1, 1), size(cg%f_tgt%seg(g)%buf(:, :, :)), MPI_DOUBLE_PRECISION, cg%f_tgt%seg(g)%proc, cg%f_tgt%seg(g)%tag, comm, req(nr), ierr)

         enddo
         cgl => cgl%nxt
      enddo

      if (nr>0) call MPI_Waitall(nr, req(:nr), status(:,:nr), ierr)

      ! interpolate received coarse data to the right place
      cgl => fine%first
      do while (associated(cgl))
         cg => cgl%cg
         if (allocated(cg%c_tgt%seg)) then

            cg%q(iv)%arr(:,:,:) = 0. ! disables check_dirty

            do g = lbound(cg%c_tgt%seg(:), dim=1), ubound(cg%c_tgt%seg(:), dim=1)
               fse => cg%c_tgt%seg(g)%se
               if (associated(this%first)) then
                  off1(:) = mod(cg%off(:)+fse(:, LO) - this%first%cg%ijkse(:, LO), 2_LONG)          !!!      this%first%cg%ijkse ???
               else
                  off1(:) = mod(cg%off(:)+fse(:, LO) - 1, 2_LONG) !!!
               endif
               do k = fse(zdim, LO), fse(zdim, HI)
                  kc = (k-fse(zdim, LO)+off1(zdim))/2 + 1
                  do j = fse(ydim, LO), fse(ydim, HI)
                     jc = (j-fse(ydim, LO)+off1(ydim))/2 + 1
                     do i = fse(xdim, LO), fse(xdim, HI)
                        ic = (i-fse(xdim, LO)+off1(xdim))/2 + 1
                        cg%q(iv)%arr(i, j, k) = cg%c_tgt%seg(g)%buf(ic, jc, kc)
                     enddo
                  enddo
               enddo
            enddo

         endif
         cgl => cgl%nxt
      enddo

   end subroutine prolong0_q_1var

!> \brief Print detailed information about current level decomposition

   subroutine print_segments(this)

      use constants,  only: LO, HI
      use dataio_pub, only: printinfo, msg
      use mpisetup,   only: proc, FIRST, LAST, master, nproc

      implicit none

      class(cg_list_level), intent(in) :: this

      integer :: p, i, hl
      integer(kind=8) :: ccnt
      real, allocatable, dimension(:) :: maxcnt
      character(len=len(msg)) :: header

      if (.not. master) return

      !call dom%print_me
      allocate(maxcnt(FIRST:LAST))
      maxcnt(:) = 0
      do p = FIRST, LAST
         hl = 0
         do i = lbound(this%pse(p)%sel(:, :, :), dim=1), ubound(this%pse(p)%sel(:, :, :), dim=1)
            ccnt = product(this%pse(p)%sel(i, :, HI) - this%pse(p)%sel(i, :, LO) + 1)
            maxcnt(p) = maxcnt(p) + ccnt
            if (i == 1) then
               write(header, '(a,i4)')"[cg_list_level:print_segments] segment @", p
               hl = len_trim(header)
            else
               header = repeat(" ", hl)
            endif
            write(msg,'(2a,2(3i6,a),i8,a)') header(:hl), " : [", this%pse(p)%sel(i, :, LO), "] : [", this%pse(p)%sel(i, :, HI), "] #", ccnt, " cells"
            call printinfo(msg)
         enddo
      enddo
      write(msg,'(a,i3,a,f8.5)')"[cg_list_level:print_segments] Load balance at level ", this%lev," : ",product(real(this%n_d(:)))/(nproc*maxval(maxcnt(:)))
      !> \todo add calculation of total internal boundary surface in cells
      call printinfo(msg)
      deallocate(maxcnt)

   end subroutine print_segments

!> \brief Initialize newest grid container

   subroutine init_new_cg(this, gr_id)

      use mpisetup, only: proc

      implicit none

      class(cg_list_level), intent(inout) :: this
      integer,              intent(in)    :: gr_id

      call this%add
      call this%last%cg%init(this%n_d, this%pse(proc)%sel(gr_id, :, :), gr_id, this%lev) ! we cannot pass "this" as an argument because of circular dependencies
      call mpi_bnd_types(this, this%last%cg)
      call this%last%cg%set_q_mbc

!! \todo add an optional argument, array of pointers to lists, where the cg should be added. Requires polymorphic array of pointers.

   end subroutine init_new_cg

!> \brief Count all cg on current level. Useful for computing tags in vertical_prep

   subroutine update_tot_se(this)

      use mpisetup, only: FIRST, LAST

      implicit none

      class(cg_list_level), intent(inout) :: this

      integer :: p

      this%tot_se = 0
      do p = FIRST, LAST
         this%tot_se = this%tot_se + ubound(this%pse(p)%sel(:,:,:), dim=1)
      enddo
   end subroutine update_tot_se

!> \brief Create MPI types for boundary exchanges

   subroutine mpi_bnd_types(this, cg)

      use constants,  only: FLUID, ARR, xdim, zdim, ndims, LO, HI, BND, BLK, I_ONE
      use dataio_pub, only: die
      use domain,     only: dom
      use fluidindex, only: flind
      use grid_cont,  only: grid_container, is_overlap
      use mpi,        only: MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, MPI_COMM_NULL
      use mpisetup,   only: ierr, proc, FIRST, LAST, procmask
      use types,      only: cdd

      implicit none

      class(cg_list_level), intent(in)    :: this
      type(grid_container), intent(inout) :: cg

      integer(kind=4), dimension(:), allocatable :: sizes, subsizes, starts
      integer :: t, g, j, b
      integer(kind=4) :: d, dd, hl, lh, ib
      integer(kind=4), parameter, dimension(FLUID:ARR) :: dims = [ I_ONE+ndims, I_ONE+ndims, I_ONE+ndims, ndims ] !< dimensionality of arrays
      integer(kind=4), dimension(FLUID:ARR) :: nc
      integer(kind=8), dimension(xdim:zdim) :: ijks, per
      integer(kind=8), dimension(xdim:zdim, LO:HI) :: b_layer, bp_layer, poff

      nc = [ flind%all, ndims, max(flind%crs%all,I_ONE), I_ONE ]      !< number of fluids, magnetic field components, CRs, and 1 for a rank-3 array

      if (cg%level_id /= this%lev) call die("[cg_list_level:mpi_bnd_types] Level mismatch")

      if (allocated(cg%i_bnd) .or. allocated(cg%o_bnd)) call die("[cg_list_level:mpi_bnd_types] cg%i_bnd or cg%o_bnd already allocated")
      allocate(cg%i_bnd(xdim:zdim, dom%nb), cg%o_bnd(xdim:zdim, dom%nb))

      ! assume that cuboids fill the domain and don't collide

      ijks(:) = cg%ijkse(:, LO) - cg%off(:)
      per(:) = 0
      where (dom%periodic(:)) per(:) = this%n_d(:)

      if (cdd%comm3d == MPI_COMM_NULL) then

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
                     do b = lbound(this%pse(j)%sel(:, :, :), dim=1), ubound(this%pse(j)%sel(:, :, :), dim=1)
                        if (is_overlap(b_layer(:,:), this%pse(j)%sel(b, :, :), per(:))) procmask(j) = procmask(j) + 1 ! count how many boundaries we share with that process
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
                        do b = lbound(this%pse(j)%sel(:, :, :), dim=1), ubound(this%pse(j)%sel(:, :, :), dim=1)
                           b_layer(:,:) = cg%my_se(:, :)
                           b_layer(d, lh) = b_layer(d, lh) + lh-hl
                           b_layer(d, hl) = b_layer(d, lh) ! adjacent boundary layer, 1 cell wide, without corners
                           bp_layer(:, :) = b_layer(:, :)
                           where (per(:) > 0)
                              bp_layer(:, LO) = mod(b_layer(:, LO) + per(:), per(:))
                              bp_layer(:, HI) = mod(b_layer(:, HI) + per(:), per(:)) ! adjacent boundary layer, 1 cell wide, without corners, corrected for periodicity
                           endwhere
                           !> \todo save b_layer(:,:) and bp_layer(:,:) and move above calculations outside the b loop

                           if (is_overlap(bp_layer(:,:), this%pse(j)%sel(b, :, :))) then
                              poff(:,:) = bp_layer(:,:) - b_layer(:,:) ! displacement due to periodicity
                              bp_layer(:, LO) = max(bp_layer(:, LO), this%pse(j)%sel(b, :, LO))
                              bp_layer(:, HI) = min(bp_layer(:, HI), this%pse(j)%sel(b, :, HI))
                              b_layer(:,:) = bp_layer(:,:) - poff(:,:)
                              g = g + 1
                              do ib = 1, dom%nb
                                 cg%i_bnd(d, ib)%seg(g)%proc = j
                                 cg%i_bnd(d, ib)%seg(g)%se(:,LO) = b_layer(:, LO) + ijks(:)
                                 cg%i_bnd(d, ib)%seg(g)%se(:,HI) = b_layer(:, HI) + ijks(:)
                                 if (any(cg%i_bnd(d, ib)%seg(g)%se(d, :) < 0)) &
                                      cg%i_bnd(d, ib)%seg(g)%se(d, :) = cg%i_bnd(d, ib)%seg(g)%se(d, :) + this%n_d(d)
                                 if (any(cg%i_bnd(d, ib)%seg(g)%se(d, :) > cg%n_b(d) + 2*dom%nb)) &
                                      cg%i_bnd(d, ib)%seg(g)%se(d, :) = cg%i_bnd(d, ib)%seg(g)%se(d, :) - this%n_d(d)

                                 ! expand to cover corners (requires separate MPI_Waitall for each direction)
                                 !! \todo create separate %mbc for corner-less exchange with one MPI_Waitall (can scale better)
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
            endif
         enddo

      else

         do d = xdim, zdim
            if (dom%has_dir(d)) then
               do t = FLUID, ARR  ! fluid, Bfield, wcr, grav

                  allocate(sizes(dims(t)), subsizes(dims(t)), starts(dims(t)))
                  if (dims(t) == 1+ndims) sizes(1) = nc(t)
                  sizes(dims(t)-zdim+xdim:dims(t)) = cg%n_(:)

                  do ib = 1, dom%nb

                     subsizes(:) = sizes(:)
                     subsizes(dims(t)-zdim+d) = ib
                     starts(:) = 0

                     starts(dims(t)-zdim+d) = dom%nb-ib
                     call MPI_Type_create_subarray(dims(t), sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, cg%mbc(t, d, LO, BND, ib), ierr)
                     call MPI_Type_commit(cg%mbc(t, d, LO, BND, ib), ierr)

                     starts(dims(t)-zdim+d) = dom%nb
                     call MPI_Type_create_subarray(dims(t), sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, cg%mbc(t, d, LO, BLK, ib), ierr)
                     call MPI_Type_commit(cg%mbc(t, d, LO, BLK, ib), ierr)

                     starts(dims(t)-zdim+d) = cg%n_b(d) + dom%nb - ib
                     call MPI_Type_create_subarray(dims(t), sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, cg%mbc(t, d, HI, BLK, ib), ierr)
                     call MPI_Type_commit(cg%mbc(t, d, HI, BLK, ib), ierr)

                     starts(dims(t)-zdim+d) = cg%ijkse(d, HI)
                     call MPI_Type_create_subarray(dims(t), sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, cg%mbc(t, d, HI, BND, ib), ierr)
                     call MPI_Type_commit(cg%mbc(t, d, HI, BND, ib), ierr)

                  enddo

                  deallocate(sizes, subsizes, starts)

               enddo
            endif
         enddo

      endif

   end subroutine mpi_bnd_types

end module cg_list_lev
