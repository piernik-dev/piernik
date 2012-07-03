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
      integer(kind=4) :: lev                   !< level number (relative to base level). For printing, debug, and I/O use only. No arithmetic should depend on it.
      integer(kind=8), dimension(ndims) :: n_d !< maximum number of grid cells in each direction (size of fully occupied level)
      type(cuboids), dimension(:), allocatable :: pse  !< lists of grid chunks on each process (FIRST:LAST); Use with care, because this is an antiparallel thing
      integer :: tot_se                        !< global number of segments on the level
      type(cg_list_level), pointer :: coarser  !< coarser level cg set or null()
      type(cg_list_level), pointer :: finer    !< finer level cg set or null()
      integer(kind=4) :: ord_prolong_set       !< Number of boundary cells for prolongation used in last update of cg_list_level%vertical_prep
      integer :: fft_type                      !< type of FFT to employ in some multigrid solvers (depending on boundaries)
    contains

      ! Level management
      procedure :: init_lev, init_lev_base
      generic, public :: init => init_lev, init_lev_base
      procedure :: init_all_new_cg             !< initialize newest grid container
      procedure, private :: mpi_bnd_types      !< create MPI types for boundary exchanges
      procedure :: print_segments              !< print detailed information about current level decomposition

      ! Prolongation and restriction
      procedure, private :: vertical_prep      !< initialize prolongation and restriction targets
      procedure, private :: update_tot_se      !< count all cg on current level for computing tags in vertical_prep
      procedure :: prolong                     !< interpolate the grid data which has the flag vital set to this%finer level
      procedure :: restrict                    !< interpolate the grid data which has the flag vital set from this%coarser level
      procedure :: prolong_q_1var              !< interpolate the grid data in specified q field to this%finer level
      procedure :: restrict_q_1var             !< interpolate the grid data in specified q field from this%coarser level
      procedure :: restrict_to_floor_q_1var    !< restrict specified q field as much as possible

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

!> \brief initialize the base level

   subroutine init_lev_base(this, n_d)

      use constants,  only: INVALID, base_level_id, fft_none
      use dataio_pub, only: die
      use domain,     only: dom

      implicit none

      class(cg_list_level),              intent(inout) :: this
      integer(kind=4), dimension(ndims), intent(in)    :: n_d    !< size of global base grid in cells

      if (any(n_d(:) < 1)) call die("[cg_list_level::init_lev_base] non-positive base grid sizes")
      if (any(dom%has_dir(:) .neqv. (n_d(:) > 1))) call die("[cg_list_level::init_lev_base] base grid size incompatible with has_dir masks")

      this%lev = base_level_id
      this%n_d(:) = n_d(:)
      this%tot_se = 0
      this%coarser => null()
      this%finer => null()
      this%ord_prolong_set = INVALID
      this%fft_type = fft_none
      call this%init

   end subroutine init_lev_base

!> \brief add a fine or coarse level to a existing one

   subroutine init_lev(this, link, coarse)

      use constants,  only: INVALID, I_ONE, refinement_factor, fft_none
      use dataio_pub, only: die, msg
      use domain,     only: dom
      use mpisetup,   only: master

      implicit none

      class(cg_list_level), target, intent(inout) :: this
      type(cg_list_level), pointer, intent(inout) :: link    !< lowest or highest refinement level (cannot refer to base_lev here due to cyclic deps)
      logical,                      intent(in)    :: coarse  !< if .true. then add a level below base level

      if (.not. associated(link)) call die("[cg_list_level::init_lev] cannot link to null")

      this%n_d(:) = 1
      this%tot_se = 0
      if (coarse) then
         if (associated(link%coarser)) call die("[cg_list_level::init_lev] coarser level already exists")
         this%lev = link%lev - I_ONE
         where (dom%has_dir(:))
            this%n_d(:) = link%n_d(:) / refinement_factor
         elsewhere
            this%n_d(:) = I_ONE
         endwhere
         if (master .and. any(this%n_d(:)*refinement_factor /= link%n_d(:) .and. dom%has_dir(:))) then
            write(msg, '(a,3f10.1,a,i3)')"[cg_list_level::init_lev] Fractional number of domain cells: ", link%n_d(:)/real(refinement_factor), " at level ",this%lev
            call die(msg)
         endif
         this%coarser => null()
         this%finer => link
         select type(this)
            type is (cg_list_level)
               link%coarser => this
            class default
               call die("[cg_list_level::init_lev] cannot call this routine for derivatives of cg_list_level (coarse)")
         end select
      else
         if (associated(link%finer)) call die("[cg_list_level::init_lev] finer level already exists")
         this%lev = link%lev + I_ONE
         this%n_d(:) = link%n_d(:) * refinement_factor
         this%coarser => link
         this%finer => null()
         select type(this)
            type is (cg_list_level)
               link%finer => this
            class default
               call die("[cg_list_level::init_lev] cannot call this routine for derivatives of cg_list_level (fine)")
         end select
      endif

      !! make sure that vertical_prep will be called where necessary
      this%ord_prolong_set = INVALID
      link%ord_prolong_set = INVALID
      this%fft_type = fft_none
      call this%init

   end subroutine init_lev

!> \brief Calculate minimal fine grid that completely covers given coarse grid

   pure function c2f(coarse) result (fine)

      use constants, only: xdim, zdim, LO, HI, refinement_factor
      use domain,    only: dom

      implicit none

      integer(kind=8), dimension(xdim:zdim, LO:HI), intent(in) :: coarse

      integer(kind=8), dimension(xdim:zdim, LO:HI) :: fine

      fine(:,:) = coarse(:,:) * refinement_factor
      where (dom%has_dir(:)) fine(:, HI) = fine(:, HI) + refinement_factor - 1

   end function c2f

!> \brief Calculate minimal coarse grid that completely embeds given fine grid

   pure function f2c(fine) result (coarse)

      use constants, only: xdim, zdim, LO, HI, refinement_factor, LONG

      implicit none

      integer(kind=8), dimension(xdim:zdim, LO:HI), intent(in) :: fine

      integer(kind=8), dimension(xdim:zdim, LO:HI) :: coarse

      integer(kind=8) :: bias !< prevents inconsistiences in arithmetic on integers due to rounding towards 0 (stencil_range should be more than enough)

      bias = max(0_LONG, -minval(fine(:,:)))
      coarse(:,:) =  (fine(:,:) + bias*refinement_factor) / refinement_factor - bias

   end function f2c

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

      use cg_list_global, only: all_cg
      use constants,      only: xdim, ydim, zdim, LO, HI
      use dataio_pub,     only: die
      use domain,         only: dom
      use gc_list,        only: cg_list_element
      use grid_cont,      only: pr_segment, grid_container, is_overlap
      use mpisetup,       only: FIRST, LAST

      implicit none

      class(cg_list_level), intent(inout) :: this

      integer :: g, j, jf, fmax, tag
      integer(kind=8), dimension(xdim:zdim, LO:HI) :: coarsened
      integer, dimension(xdim:zdim, LO:HI) :: enlargement
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

      if (all_cg%ord_prolong_nb == this%ord_prolong_set) return ! no need to update vertical communication on this level

      ! enlarge the fine blocks to allow for high orders of interpolation
      enlargement(:, LO) = -dom%D_(:) * all_cg%ord_prolong_nb
      enlargement(:, HI) =  dom%D_(:) * all_cg%ord_prolong_nb

      fine => this%finer
      ! find fine target for receiving restricted data or sending data to be prolonged
      if (associated(fine)) then
         cgl => this%first
         do while (associated(cgl))
            cg => cgl%cg

            if (allocated(ps)) call die("cll:vp f a ps")
            fmax = 0
            do j = FIRST, LAST
               fmax = fmax + size(fine%pse(j)%sel(:, :, :), dim=1)
            enddo
            allocate(ps(fmax))

            if (allocated(cg%ri_tgt%seg)) deallocate(cg%ri_tgt%seg) ! call warn("cll:vp cg%ri_tgt%seg a a")
            if (allocated(cg%po_tgt%seg)) deallocate(cg%po_tgt%seg)
            g = 0
            do j = FIRST, LAST
               do jf = lbound(fine%pse(j)%sel(:, :, :), dim=1), ubound(fine%pse(j)%sel(:, :, :), dim=1)
                  if (is_overlap(c2f(cg%my_se(:, :)), fine%pse(j)%sel(jf, :, :))) then
                     g = g + 1
                     ps(g) = int_pair(j, jf)
                  endif
               enddo
            enddo
            allocate(cg%ri_tgt%seg(g), cg%po_tgt%seg(g))

            !! \todo fuse the following two loops back into one
            do g = lbound(cg%ri_tgt%seg(:), dim=1), ubound(cg%ri_tgt%seg(:), dim=1)
               seg => cg%ri_tgt%seg(g)
               if (allocated(seg%buf)) call die("cll:vp fr seg%buf a a")
               seg%proc = ps(g)%proc
               ! find cross-section of own segment with coarsened fine segment
               coarsened(:,:) = f2c(fine%pse(seg%proc)%sel(ps(g)%n_se, :, :))
               seg%se(:, LO) = max(cg%my_se(:, LO), coarsened(:, LO))
               seg%se(:, HI) = min(cg%my_se(:, HI), coarsened(:, HI))
               allocate(seg%buf(seg%se(xdim, HI)-seg%se(xdim, LO) + 1, &
                    &           seg%se(ydim, HI)-seg%se(ydim, LO) + 1, &
                    &           seg%se(zdim, HI)-seg%se(zdim, LO) + 1))
               tag = cg%grid_id + this%tot_se * ps(g)%n_se
               seg%tag = int(tag, kind=4) ! assumed that there is only one piece to be communicated from grid to grid (i.e. grids are not periodically wrapped around)
               if (tag /= int(seg%tag)) call die("[cg_list_level:vertical_prep] tag overflow (ri)")
            enddo

            ! When fine and coarse pieces are within prolongation stencil length, but there is no direct overlap we rely on guardcell update on coarse side
            ! to provide valid prolongation data. An alternative approach would not require updating guardcells, but would usually result in setting up
            ! many small teransactions on the corners, especially on cartesian decomposition.
            !
            ! When the overlap is one fine cell (half coarse cell), not all communications are necessary, but we don't want to complicate the code too much at the moment
            do g = lbound(cg%po_tgt%seg(:), dim=1), ubound(cg%po_tgt%seg(:), dim=1)
               seg => cg%po_tgt%seg(g)
               if (allocated(seg%buf)) call die("cll:vp fp seg%buf a a")
               seg%proc = ps(g)%proc
               ! find cross-section of own segment with enlarged coarsened fine segment
               coarsened(:,:) = f2c(fine%pse(seg%proc)%sel(ps(g)%n_se, :, :)) + enlargement(:,:)
               seg%se(:, LO) = max(cg%my_se(:, LO) + enlargement(:, LO), coarsened(:, LO))
               seg%se(:, HI) = min(cg%my_se(:, HI) + enlargement(:, HI), coarsened(:, HI))
               allocate(seg%buf(seg%se(xdim, HI)-seg%se(xdim, LO) + 1, &
                    &           seg%se(ydim, HI)-seg%se(ydim, LO) + 1, &
                    &           seg%se(zdim, HI)-seg%se(zdim, LO) + 1))
               tag = cg%grid_id + this%tot_se * ps(g)%n_se
               seg%tag = int(tag, kind=4) ! assumed that there is only one piece to be communicated from grid to grid (i.e. grids are not periodically wrapped around)
               if (tag /= int(seg%tag)) call die("[cg_list_level:vertical_prep] tag overflow po)")
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

            if (allocated(ps)) call die("cll:vp c a ps")
            fmax = 0
            do j = FIRST, LAST
               fmax = fmax + size(coarse%pse(j)%sel(:, :, :), dim=1)
            enddo
            allocate(ps(fmax))

            if (allocated(cg%ro_tgt%seg)) deallocate(cg%ro_tgt%seg) ! call warn("cll:vp cg%ro_tgt%seg a a")
            if (allocated(cg%pi_tgt%seg)) deallocate(cg%pi_tgt%seg)
            g = 0
            do j = FIRST, LAST
               do jf = lbound(coarse%pse(j)%sel(:, :, :), dim=1), ubound(coarse%pse(j)%sel(:, :, :), dim=1)
                  if (is_overlap(cg%my_se(:, :), c2f(coarse%pse(j)%sel(jf, :, :)))) then
                     g = g + 1
                     ps(g) = int_pair(j, jf)
                  endif
               enddo
            enddo
            allocate(cg%ro_tgt%seg(g), cg%pi_tgt%seg(g))

            !! \todo fuse the following two loops back into one
            do g = lbound(cg%ro_tgt%seg(:), dim=1), ubound(cg%ro_tgt%seg(:), dim=1)
               seg => cg%ro_tgt%seg(g)
               if (allocated(seg%buf)) call die("cll:vp cr seg%buf a a")
               seg%proc = ps(g)%proc
               ! find cross-section of coarsened own segment with coarse segment
               coarsened(:,:) = f2c(cg%my_se(:,:))
               seg%se(:, LO) = max(coarsened(:, LO), coarse%pse(seg%proc)%sel(ps(g)%n_se, :, LO))
               seg%se(:, HI) = min(coarsened(:, HI), coarse%pse(seg%proc)%sel(ps(g)%n_se, :, HI))
               allocate(seg%buf(seg%se(xdim, HI)-seg%se(xdim, LO) + 1, &
                    &           seg%se(ydim, HI)-seg%se(ydim, LO) + 1, &
                    &           seg%se(zdim, HI)-seg%se(zdim, LO) + 1))
               coarsened(:,:) = c2f(coarse%pse(seg%proc)%sel(ps(g)%n_se, :, :)) ! should be renamed to refined(:,:)
               seg%se(:, LO) = max(cg%my_se(:, LO), coarsened(:, LO))
               seg%se(:, HI) = min(cg%my_se(:, HI), coarsened(:, HI))
               tag = ps(g)%n_se + coarse%tot_se * cg%grid_id
               seg%tag = int(tag, kind=4)
               if (tag /= int(seg%tag)) call die("[cg_list_level:vertical_prep] tag overflow (ro)")
            enddo

            do g = lbound(cg%pi_tgt%seg(:), dim=1), ubound(cg%pi_tgt%seg(:), dim=1)
               seg => cg%pi_tgt%seg(g)
               if (allocated(seg%buf)) call die("cll:vp cp seg%buf a a")
               seg%proc = ps(g)%proc
               ! find cross-section of coarsened own segment with enlarged coarse segment
               coarsened(:,:) = f2c(cg%my_se(:,:)) + enlargement(:,:)
               seg%se(:, LO) = max(coarsened(:, LO), coarse%pse(seg%proc)%sel(ps(g)%n_se, :, LO) + enlargement(:, LO))
               seg%se(:, HI) = min(coarsened(:, HI), coarse%pse(seg%proc)%sel(ps(g)%n_se, :, HI) + enlargement(:, HI))
               allocate(seg%buf(seg%se(xdim, HI)-seg%se(xdim, LO) + 1, &
                    &           seg%se(ydim, HI)-seg%se(ydim, LO) + 1, &
                    &           seg%se(zdim, HI)-seg%se(zdim, LO) + 1))
               tag = ps(g)%n_se + coarse%tot_se * cg%grid_id
               seg%tag = int(tag, kind=4)
               if (tag /= int(seg%tag)) call die("[cg_list_level:vertical_prep] tag overflow (pi)")
            enddo

            if (allocated(ps)) deallocate(ps)
            cgl => cgl%nxt
         enddo
      endif

      this%ord_prolong_set = all_cg%ord_prolong_nb

      !> \todo update mpisetup::req(:)

   end subroutine vertical_prep

!>
!! \brief interpolate the grid data which has the flag vital set to this%finer level
!!
!! \todo implement it in a more efficient way (reqiures a lot more temporary buffers)
!<
   subroutine prolong(this)

      use constants,   only: base_level_id, wa_n
      use dataio_pub,  only: warn
      use named_array, only: qna, wna

      implicit none

      class(cg_list_level), target, intent(inout) :: this !< object invoking type-bound procedure

      integer :: i, iw, iwa

      do i = lbound(qna%lst(:), dim=1), ubound(qna%lst(:), dim=1)
         if (qna%lst(i)%vital .and. (qna%lst(i)%multigrid .or. this%lev > base_level_id)) call this%prolong_q_1var(i)
      enddo

      iwa = qna%ind(wa_n)

      do i = lbound(wna%lst(:), dim=1), ubound(wna%lst(:), dim=1)
         if (wna%lst(i)%vital .and. (wna%lst(i)%multigrid .or. this%lev >= base_level_id)) then
            if (wna%lst(i)%multigrid) call warn("[cg_list_level:prolong] mg set for cg%w ???")
            do iw = 1, wna%lst(i)%dim4
               call this%wq_copy(i, iw, iwa)
               call this%prolong_q_1var(iwa)
               call this%finer%qw_copy(iwa, i, iw)
            enddo
         endif
      enddo

   end subroutine prolong

!>
!! \brief interpolate the grid data which has the flag vital set from this%coarser level
!!
!! \todo implement it in a more efficient way (reqiures a lot more temporary buffers)
!<
   subroutine restrict(this)

      use constants,   only: base_level_id, wa_n
      use dataio_pub,  only: warn
      use named_array, only: qna, wna

      implicit none

      class(cg_list_level), target, intent(inout) :: this !< object invoking type-bound procedure

      integer :: i, iw, iwa

      do i = lbound(qna%lst(:), dim=1), ubound(qna%lst(:), dim=1)
         if (qna%lst(i)%vital .and. (qna%lst(i)%multigrid .or. this%lev > base_level_id)) call this%restrict_q_1var(i)
      enddo

      iwa = qna%ind(wa_n)

      do i = lbound(wna%lst(:), dim=1), ubound(wna%lst(:), dim=1)
         if (wna%lst(i)%vital .and. (wna%lst(i)%multigrid .or. this%lev > base_level_id)) then
            if (wna%lst(i)%multigrid) call warn("[cg_list_level:restrict] mg set for cg%w ???")
            do iw = 1, wna%lst(i)%dim4
               call this%wq_copy(i, iw, iwa)
               call this%restrict_q_1var(iwa)
               call this%coarser%qw_copy(iwa, i, iw)
            enddo
         endif
      enddo

   end subroutine restrict

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

      use constants,   only: xdim, ydim, zdim, LO, HI, I_ONE, refinement_factor
      use dataio_pub,  only: msg, warn, die
      use domain,      only: dom
      use gc_list,     only: cg_list_element
      use grid_cont,   only: grid_container
      use mpisetup,    only: comm, mpi_err, req, status
      use mpi,         only: MPI_DOUBLE_PRECISION
      use named_array, only: p3

      implicit none

      class(cg_list_level), target, intent(inout) :: this !< object invoking type-bound procedure
      integer,                      intent(in)    :: iv   !< variable to be restricted

      type(cg_list_level), pointer :: coarse
      integer :: g
      integer(kind=8), dimension(xdim:zdim, LO:HI) :: fse, cse ! shortcuts for fine segment and coarse segment
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

      call this%vertical_prep
      call coarse%vertical_prep

      ! be ready to receive everything into right buffers
      nr = 0
      cgl => coarse%first
      do while (associated(cgl))
         cg => cgl%cg
         if (allocated(cg%ri_tgt%seg)) then
            do g = lbound(cg%ri_tgt%seg(:), dim=1), ubound(cg%ri_tgt%seg(:), dim=1)
               nr = nr + I_ONE
               if (nr > size(req, dim=1)) call die("[cg_list_level:restrict_q_1var] size(req) too small for Irecv")
               call MPI_Irecv(cg%ri_tgt%seg(g)%buf(1, 1, 1), size(cg%ri_tgt%seg(g)%buf(:, :, :)), MPI_DOUBLE_PRECISION, cg%ri_tgt%seg(g)%proc, cg%ri_tgt%seg(g)%tag, comm, req(nr), mpi_err)
            enddo
         endif
         cgl => cgl%nxt
      enddo

      ! interpolate to coarse buffer and send it
      norm = 1./refinement_factor**dom%eff_dim
      cgl => this%first
      do while (associated(cgl))
         cg => cgl%cg
         do g = lbound(cg%ro_tgt%seg(:), dim=1), ubound(cg%ro_tgt%seg(:), dim=1)

            fse(:,:) = cg%ro_tgt%seg(g)%se(:,:)
            fse(:, LO) = fse(:, LO) - cg%off(:)
            fse(:, HI) = fse(:, HI) - cg%off(:)

            cg%ro_tgt%seg(g)%buf(:, :, :) = 0.
            off1(:) = mod(cg%ro_tgt%seg(g)%se(:, LO), int(refinement_factor, kind=8))
            do k = fse(zdim, LO), fse(zdim, HI)
               kc = (k-fse(zdim, LO)+off1(zdim))/refinement_factor + 1
               do j = fse(ydim, LO), fse(ydim, HI)
                  jc = (j-fse(ydim, LO)+off1(ydim))/refinement_factor + 1
                  do i = fse(xdim, LO), fse(xdim, HI)
                     ic = (i-fse(xdim, LO)+off1(xdim))/refinement_factor + 1
                     cg%ro_tgt%seg(g)%buf(ic, jc, kc) = cg%ro_tgt%seg(g)%buf(ic, jc, kc) + cg%q(iv)%arr(i+cg%is, j+cg%js, k+cg%ks) * norm
                  enddo
               enddo
            enddo
            nr = nr + I_ONE
            if (nr > size(req, dim=1)) call die("[cg_list_level:restrict_q_1var] size(req) too small for Isend")
            call MPI_Isend(cg%ro_tgt%seg(g)%buf(1, 1, 1), size(cg%ro_tgt%seg(g)%buf(:, :, :)), MPI_DOUBLE_PRECISION, cg%ro_tgt%seg(g)%proc, cg%ro_tgt%seg(g)%tag, comm, req(nr), mpi_err)
         enddo
         cgl => cgl%nxt
      enddo

      if (nr>0) call MPI_Waitall(nr, req(:nr), status(:,:nr), mpi_err)

      ! copy the received buffers to the right places
      cgl => coarse%first
      do while (associated(cgl))
         cg => cgl%cg
         if (allocated(cg%ri_tgt%seg)) then
            cg%q(iv)%arr(:,:,:) = 0. ! disables check_dirty
            do g = lbound(cg%ri_tgt%seg(:), dim=1), ubound(cg%ri_tgt%seg(:), dim=1)
               cse(:,:) = cg%ri_tgt%seg(g)%se(:,:)
               cse(:, LO) = cse(:, LO) - cg%off(:) + cg%ijkse(:, LO)
               cse(:, HI) = cse(:, HI) - cg%off(:) + cg%ijkse(:, LO)

               p3 => cg%q(iv)%span(cse)
               p3 = p3 + cg%ri_tgt%seg(g)%buf(:, :, :) !errors on overlap?
            enddo
         endif
         cgl => cgl%nxt
      enddo

   end subroutine restrict_q_1var

!>
!! \brief high-order order prolongation interpolation
!!
!! \details
!! <table border="1" cellpadding="4" cellspacing="0">
!!   <tr><td> Cell-face prolongation stencils for fast convergence on uniform grid </td>
!!       <td> -1./12. </td><td> 7./12. </td><td> 7./12. </td><td> -1./12. </td><td> integral cubic </td></tr>
!!   <tr><td> Slightly slower convergence, less wide stencil  </td>
!!       <td>         </td><td> 1./2.  </td><td> 1./2.  </td><td>         </td><td> average; integral and direct linear </td></tr>
!! </table>
!!\n Prolongation of cell faces from cell centers are required for FFT local solver, red-black Gauss-Seidel relaxation don't use it.
!!
!!\n Cell-centered prolongation stencils, for odd fine cells, for even fine cells reverse the order.
!! <table border="1" cellpadding="4" cellspacing="0">
!!   <tr><td> 35./2048. </td><td> -252./2048. </td><td> 1890./2048. </td><td> 420./2048. </td><td> -45./2048. </td><td> direct quartic </td></tr>
!!   <tr><td>           </td><td>   -7./128.  </td><td>  105./128.  </td><td>  35./128.  </td><td>  -5./128.  </td><td> direct cubic </td></tr>
!!   <tr><td>           </td><td>   -3./32.   </td><td>   30./32.   </td><td>   5./32.   </td><td>            </td><td> direct quadratic </td></tr>
!!   <tr><td>           </td><td>             </td><td>    1.       </td><td>            </td><td>            </td><td> injection (0-th order), direct and integral approach </td></tr>
!!   <tr><td>           </td><td>             </td><td>    3./4.    </td><td>   1./4.    </td><td>            </td><td> linear, direct and integral approach </td></tr>
!!   <tr><td>           </td><td>   -1./8.    </td><td>    1.       </td><td>   1./8.    </td><td>            </td><td> integral quadratic </td></tr>
!!   <tr><td>           </td><td>   -5./64.   </td><td>   55./64.   </td><td>  17./64.   </td><td>  -3./64.   </td><td> integral cubic </td></tr>
!!   <tr><td>   3./128. </td><td>  -11./64.   </td><td>    1.       </td><td>  11./64.   </td><td>  -3./128.  </td><td> integral quartic </td></tr>
!! </table>
!!
!!\n General rule is that the second big positive coefficient should be assigned to closer neighbor of the coarse parent cell.
!!\n Thus a single coarse contributes to fine cells in the following way:
!! <table border="1" cellpadding="4" cellspacing="0">
!!   <tr><td> fine level   </td>
!!       <td> -3./128. </td><td> 3./128. </td><td> -11./64. </td><td>  11./64. </td><td> 1. </td><td> 1. </td>
!!       <td> 11./64. </td><td> -11./64. </td><td> 3./128.  </td><td> -3./128. </td><td> integral quartic coefficients </td></tr>
!!   <tr><td> coarse level </td>
!!       <td colspan="2">                </td><td colspan="2">                 </td><td colspan="2"> 1.  </td>
!!       <td colspan="2">                </td><td colspan="2">                 </td><td>                               </td></tr>
!! </table>
!!\n
!!\n The term "n-th order integral interpolation" here means that the prolonged values satisfy the following condition:
!!\n Integral over a cell of a n-th order polynomial fit to the nearest 5 points in each dimension on coarse level
!! is equal to the sum of similar integrals over fine cells covering the coarse cell.
!!\n
!!\n The term "n-th order direct interpolation" here means that the prolonged values are n-th order polynomial fit
!! to the nearest 5 points in each dimension on coarse level evaluated for fine cell centers.
!!\n
!!\n It seems that for 3D Cartesian grid with isolated boundaries and relaxation implemented in approximate_solution
!! direct quadratic and cubic interpolations give best norm reduction factors per V-cycle (maclaurin problem).
!!  For other boundary types, FFT implementation of approximate_solution, specific source distribution or
!!  some other multigrid scheme may give faster convergence rate.
!!\n
!!\n Estimated prolongation costs for integral quartic stencil:
!!\n  "gather" approach: loop over fine cells, each one collects weighted values from 5**3 coarse cells (125*n_fine multiplications
!!\n  "scatter" approach: loop over coarse cells, each one contributes weighted values to 10**3 fine cells (1000*n_coarse multiplications, roughly equal to cost of gather)
!!\n  "directionally split" approach: do the prolongation (either gather or scatter type) first in x direction (10*n_coarse multiplications -> 2*n_coarse intermediate cells
!!                                  result), then in y direction (10*2*n_coarse multiplications -> 4*n_coarse intermediate cells result), then in z direction
!!                                  (10*4*n_coarse multiplications -> 8*n_coarse = n_fine cells result). Looks like 70*n_coarse multiplications.
!!                                  Will require two additional arrays for intermediate results.
!!\n  "FFT" approach: do the convolution in Fourier space. Unfortunately it is not periodic box, so we cannot use power of 2 FFT sizes. No idea how fast or slow this can be.
!!\n
!!\n For AMR or nested grid low-order prolongation schemes (injection and linear interpolation at least) are known to produce artifacts
!! on fine-coarse boundaries. For uniform grid the simplest operators are probably the fastest and give best V-cycle convergence rates.
!!
!! \deprecated These routines assume simplest domain decomposition where fine grids cover exactly the same area as coarse grids
!!
!! \todo implement high order prolongation. Watch f/c boundaries.
!<

   subroutine prolong_q_1var(this, iv)

      use constants,      only: xdim, ydim, zdim, LO, HI, I_ZERO, I_ONE, I_TWO, BND_REF, O_INJ, O_LIN, O_D2, O_D3, O_D4, O_I2, O_I3, O_I4, refinement_factor
      use dataio_pub,     only: msg, warn, die
      use domain,         only: dom
      use external_bnd,   only: arr3d_boundaries
      use gc_list,        only: cg_list_element
      use grid_cont,      only: grid_container
      use mpisetup,       only: comm, mpi_err, req, status
      use mpi,            only: MPI_DOUBLE_PRECISION
      use named_array,    only: qna

      implicit none

      class(cg_list_level), target, intent(inout) :: this !< object invoking type-bound procedure
      integer,                      intent(in)    :: iv   !< variable to be prolonged

      type(cg_list_level), pointer :: fine
      integer :: g
      integer(kind=8), dimension(xdim:zdim, LO:HI) :: fse, cse ! shortcuts for fine segment and coarse segment
      integer(kind=8) :: iec, jec, kec
      integer(kind=8), dimension(xdim:zdim) :: off, odd, D
      integer(kind=4) :: nr
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg            !< current grid container
      real :: P_2, P_1, P0, P1, P2
      integer :: stencil_range

      fine => this%finer
      if (.not. associated(fine)) then ! can't prolong finest level
         write(msg,'(a,i3)')"[gc_list:prolong_q_1var] no finer level than: ", this%lev
         call warn(msg)
         return
      endif

      call this%vertical_prep
      call fine%vertical_prep

!> \todo something like this (connected with a todo below):      if (dirty_debug) fine%first%cg%q(iv)%arr(:, :, :) = dirtyH

      select case (qna%lst(iv)%ord_prolong)
         case (O_D4)
            P_2 = 35./2048.; P_1 = -252./2048.; P0 = 1890./2048.; P1 = 420./2048.; P2 = -45./2048.
         case (O_D3)
            P_2 = 0.;        P_1 = -7./128.;    P0 = 105./128.;   P1 = 35./128.;   P2 = -5./128.
         case (O_D2)
            P_2 = 0.;        P_1 = -3./32.;     P0 = 30./32.;     P1 = 5./32.;     P2 = 0.
         case (O_LIN)
            P_2 = 0.;        P_1 = 0.;          P0 = 3./4.;       P1 = 1./4.;      P2 = 0.
         case (O_INJ)
            P_2 = 0.;        P_1 = 0.;          P0 = 1.;          P1 = 0.;         P2 = 0.
         case (O_I2)
            P_2 = 0.;        P_1 = -1./8.;      P0 = 1.;          P1 = 1./8.;      P2 = 0.
         case (O_I3)
            P_2 = 0.;        P_1 = -5./64.;     P0 = 55./64;      P1 = 17./64.;    P2 = -3./64.
         case (O_I4)
            P_2 = 3./128.;   P_1 = -11./64.;    P0 = 1.;          P1 = 11./64.;    P2 = -3./128.
         case default
            call die("[cg_list_level:prolong_q_1var] Unsupported order")
            return
      end select

      ! this is just for optimization. Setting stencil_range = I_TWO should work correctly for all interpolations.
      stencil_range = I_ZERO
      if (P_1 /= 0. .or. P1 /= 0.) stencil_range = I_ONE
      if (P_2 /= 0. .or. P2 /= 0.) stencil_range = I_TWO

      if (qna%lst(iv)%ord_prolong /= O_INJ) then
         !> \todo some variables may need special care on external boundaries
         call arr3d_boundaries(this, iv, bnd_type = BND_REF, corners = .true.) ! nb =  int(stencil_range, kind=4) ! not needed for injection
      endif

      nr = 0
      ! be ready to receive everything into right buffers
      cgl => fine%first
      do while (associated(cgl))
         cg => cgl%cg
         if (allocated(cg%pi_tgt%seg)) then
            do g = lbound(cg%pi_tgt%seg(:), dim=1), ubound(cg%pi_tgt%seg(:), dim=1)
               nr = nr + I_ONE
               if (nr > size(req, dim=1)) call die("[cg_list_level:prolong_q_1var] size(req) too small for Irecv")
               call MPI_Irecv(cg%pi_tgt%seg(g)%buf(1, 1, 1), size(cg%pi_tgt%seg(g)%buf(:, :, :)), MPI_DOUBLE_PRECISION, cg%pi_tgt%seg(g)%proc, cg%pi_tgt%seg(g)%tag, comm, req(nr), mpi_err)
            enddo
         endif
         cgl => cgl%nxt
      enddo

      ! send coarse data
      cgl => this%first
      do while (associated(cgl))
         cg => cgl%cg
         do g = lbound(cg%po_tgt%seg(:), dim=1), ubound(cg%po_tgt%seg(:), dim=1)

            cse(:, LO) = cg%po_tgt%seg(g)%se(:,LO) - cg%off(:) + cg%ijkse(:, LO)
            cse(:, HI) = cg%po_tgt%seg(g)%se(:,HI) - cg%off(:) + cg%ijkse(:, LO)

            nr = nr + I_ONE
            if (nr > size(req, dim=1)) call die("[cg_list_level:prolong_q_1var] size(req) too small for Isend")
!            cg%po_tgt%seg(g)%buf(:, :, :) = cg%q(iv)%span(cse)
            cg%po_tgt%seg(g)%buf(:, :, :) = cg%q(iv)%arr(cse(xdim, LO):cse(xdim, HI), cse(ydim, LO):cse(ydim, HI), cse(zdim, LO):cse(zdim, HI))
            call MPI_Isend(cg%po_tgt%seg(g)%buf(1, 1, 1), size(cg%po_tgt%seg(g)%buf(:, :, :)), MPI_DOUBLE_PRECISION, cg%po_tgt%seg(g)%proc, cg%po_tgt%seg(g)%tag, comm, req(nr), mpi_err)
         enddo
         cgl => cgl%nxt
      enddo

      if (nr>0) call MPI_Waitall(nr, req(:nr), status(:,:nr), mpi_err)

      ! merge received coarse data into one array and interpolate it into the right place
      cgl => fine%first
      do while (associated(cgl))
         cg => cgl%cg
         if (allocated(cg%pi_tgt%seg)) then

            cg%q(iv)%arr(:,:,:) = 0. ! disables check_dirty
            cg%prolong_(:,:,:) = 0

            do g = lbound(cg%pi_tgt%seg(:), dim=1), ubound(cg%pi_tgt%seg(:), dim=1)
               fse(:,:) = cg%pi_tgt%seg(g)%se(:,:)

               off(:) = fse(:,LO) - cg%off(:)/refinement_factor + cg%ijkse(:,LO)
               cg%prolong_(off(xdim):off(xdim)+ubound(cg%pi_tgt%seg(g)%buf, dim=1)-1, &
                    &      off(ydim):off(ydim)+ubound(cg%pi_tgt%seg(g)%buf, dim=2)-1, &
                    &      off(zdim):off(zdim)+ubound(cg%pi_tgt%seg(g)%buf, dim=3)-1) = cg%pi_tgt%seg(g)%buf(:,:,:)

            enddo

            iec = cg%is + (cg%ie - cg%is - 1)/refinement_factor + dom%D_x
            jec = cg%js + (cg%je - cg%js - 1)/refinement_factor + dom%D_y
            kec = cg%ks + (cg%ke - cg%ks - 1)/refinement_factor + dom%D_z

            !! almost all occurences of number "2" are in fact connected to refinement_factor

            ! When the grid offset is odd, the coarse data is shifted by half coarse cell (or one fine cell)
            odd(:) = int(mod(cg%off(:), int(refinement_factor, kind=8)), kind=4)

            ! When the grid offset is odd we need to apply mirrored prolongation stencil (swap even and odd stencils)
            where (dom%has_dir(:))
               D(:) = 1 - 2 * odd(:)
            elsewhere
               D(:) = 0
            endwhere

            ! Perform directional-split interpolation
            select case (stencil_range*dom%D_x) ! stencil_range or I_ZERO if .not. dom%has_dir(xdim)
               case (I_ZERO)
                  cg%prolong_x(           cg%is          :cg%ie+dom%D_x:2,       cg%js:jec, cg%ks:kec) = &
                       &      cg%prolong_(cg%is          :iec,                   cg%js:jec, cg%ks:kec)
                  if (dom%has_dir(xdim)) &
                       cg%prolong_x(      cg%is+dom%D_x  :cg%ie        :2,       cg%js:jec, cg%ks:kec) = &
                       &      cg%prolong_(cg%is+odd(xdim):iec-dom%D_x+odd(xdim), cg%js:jec, cg%ks:kec)
               case (I_ONE)
                  cg%prolong_x(           cg%is                  :cg%ie+dom%D_x:2,               cg%js-dom%D_y:jec+dom%D_y, cg%ks-dom%D_z:kec+dom%D_z) = &
                       + P1 * cg%prolong_(cg%is-D(xdim)          :iec-D(xdim),                   cg%js-dom%D_y:jec+dom%D_y, cg%ks-dom%D_z:kec+dom%D_z)   &
                       + P0 * cg%prolong_(cg%is                  :iec,                           cg%js-dom%D_y:jec+dom%D_y, cg%ks-dom%D_z:kec+dom%D_z)   &
                       + P_1* cg%prolong_(cg%is+D(xdim)          :iec+D(xdim),                   cg%js-dom%D_y:jec+dom%D_y, cg%ks-dom%D_z:kec+dom%D_z)
                  cg%prolong_x(           cg%is+dom%D_x          :cg%ie        :2,               cg%js-dom%D_y:jec+dom%D_y, cg%ks-dom%D_z:kec+dom%D_z) = &
                       + P_1* cg%prolong_(cg%is+odd(xdim)-D(xdim):iec-dom%D_x+odd(xdim)-D(xdim), cg%js-dom%D_y:jec+dom%D_y, cg%ks-dom%D_z:kec+dom%D_z)   &
                       + P0 * cg%prolong_(cg%is+odd(xdim)        :iec-dom%D_x+odd(xdim),         cg%js-dom%D_y:jec+dom%D_y, cg%ks-dom%D_z:kec+dom%D_z)   &
                       + P1 * cg%prolong_(cg%is+odd(xdim)+D(xdim):iec-dom%D_x+odd(xdim)+D(xdim), cg%js-dom%D_y:jec+dom%D_y, cg%ks-dom%D_z:kec+dom%D_z)
               case (I_TWO)
                  cg%prolong_x(           cg%is                    :cg%ie+dom%D_x:2,                 cg%js-2*dom%D_y:jec+2*dom%D_y, cg%ks-2*dom%D_z:kec+2*dom%D_z) = &
                       + P2 * cg%prolong_(cg%is-2*D(xdim)          :iec-2*D(xdim),                   cg%js-2*dom%D_y:jec+2*dom%D_y, cg%ks-2*dom%D_z:kec+2*dom%D_z)   &
                       + P1 * cg%prolong_(cg%is-  D(xdim)          :iec-  D(xdim),                   cg%js-2*dom%D_y:jec+2*dom%D_y, cg%ks-2*dom%D_z:kec+2*dom%D_z)   &
                       + P0 * cg%prolong_(cg%is                    :iec,                             cg%js-2*dom%D_y:jec+2*dom%D_y, cg%ks-2*dom%D_z:kec+2*dom%D_z)   &
                       + P_1* cg%prolong_(cg%is+  D(xdim)          :iec+  D(xdim),                   cg%js-2*dom%D_y:jec+2*dom%D_y, cg%ks-2*dom%D_z:kec+2*dom%D_z)   &
                       + P_2* cg%prolong_(cg%is+2*D(xdim)          :iec+2*D(xdim),                   cg%js-2*dom%D_y:jec+2*dom%D_y, cg%ks-2*dom%D_z:kec+2*dom%D_z)
                  cg%prolong_x(           cg%is+dom%D_x            :cg%ie:2,                         cg%js-2*dom%D_y:jec+2*dom%D_y, cg%ks-2*dom%D_z:kec+2*dom%D_z) = &
                       + P_2* cg%prolong_(cg%is+odd(xdim)-2*D(xdim):iec-dom%D_x+odd(xdim)-2*D(xdim), cg%js-2*dom%D_y:jec+2*dom%D_y, cg%ks-2*dom%D_z:kec+2*dom%D_z)   &
                       + P_1* cg%prolong_(cg%is+odd(xdim)-  D(xdim):iec-dom%D_x+odd(xdim)-  D(xdim), cg%js-2*dom%D_y:jec+2*dom%D_y, cg%ks-2*dom%D_z:kec+2*dom%D_z)   &
                       + P0 * cg%prolong_(cg%is+odd(xdim)          :iec-dom%D_x+odd(xdim),           cg%js-2*dom%D_y:jec+2*dom%D_y, cg%ks-2*dom%D_z:kec+2*dom%D_z)   &
                       + P1 * cg%prolong_(cg%is+odd(xdim)+  D(xdim):iec-dom%D_x+odd(xdim)+  D(xdim), cg%js-2*dom%D_y:jec+2*dom%D_y, cg%ks-2*dom%D_z:kec+2*dom%D_z)   &
                       + P2 * cg%prolong_(cg%is+odd(xdim)+2*D(xdim):iec-dom%D_x+odd(xdim)+2*D(xdim), cg%js-2*dom%D_y:jec+2*dom%D_y, cg%ks-2*dom%D_z:kec+2*dom%D_z)
               case default
                  call die("[cg_list_level:prolong_q_1var] unsupported stencil size")
            end select

            select case (stencil_range*dom%D_y)
               case (I_ZERO)
                  cg%prolong_xy(           cg%is:cg%ie, cg%js          :cg%je+dom%D_y:2,       cg%ks:kec) = &
                       &      cg%prolong_x(cg%is:cg%ie, cg%js          :jec,                   cg%ks:kec)
                  if (dom%has_dir(ydim)) &
                       cg%prolong_xy(      cg%is:cg%ie, cg%js+dom%D_y  :cg%je        :2,       cg%ks:kec) = &
                       &      cg%prolong_x(cg%is:cg%ie, cg%js+odd(ydim):jec-dom%D_y+odd(ydim), cg%ks:kec)
               case (I_ONE)
                  cg%prolong_xy(           cg%is:cg%ie, cg%js                  :cg%je+dom%D_y:2,               cg%ks-dom%D_z:kec+dom%D_z) = &
                       + P1 * cg%prolong_x(cg%is:cg%ie, cg%js-D(ydim)          :jec-D(ydim),                   cg%ks-dom%D_z:kec+dom%D_z)   &
                       + P0 * cg%prolong_x(cg%is:cg%ie, cg%js                  :jec,                           cg%ks-dom%D_z:kec+dom%D_z)   &
                       + P_1* cg%prolong_x(cg%is:cg%ie, cg%js+D(ydim)          :jec+D(ydim),                   cg%ks-dom%D_z:kec+dom%D_z)
                  cg%prolong_xy(           cg%is:cg%ie, cg%js+dom%D_y          :cg%je        :2,               cg%ks-dom%D_z:kec+dom%D_z) = &
                       + P_1* cg%prolong_x(cg%is:cg%ie, cg%js+odd(ydim)-D(ydim):jec-dom%D_y+odd(ydim)-D(ydim), cg%ks-dom%D_z:kec+dom%D_z)   &
                       + P0 * cg%prolong_x(cg%is:cg%ie, cg%js+odd(ydim)        :jec-dom%D_y+odd(ydim),         cg%ks-dom%D_z:kec+dom%D_z)   &
                       + P1 * cg%prolong_x(cg%is:cg%ie, cg%js+odd(ydim)+D(ydim):jec-dom%D_y+odd(ydim)+D(ydim), cg%ks-dom%D_z:kec+dom%D_z)
               case (I_TWO)
                  cg%prolong_xy(           cg%is:cg%ie, cg%js                    :cg%je+dom%D_y:2,                 cg%ks-2*dom%D_z:kec+2*dom%D_z) = &
                       + P2 * cg%prolong_x(cg%is:cg%ie, cg%js-2*D(ydim)          :jec-2*D(ydim),                   cg%ks-2*dom%D_z:kec+2*dom%D_z)   &
                       + P1 * cg%prolong_x(cg%is:cg%ie, cg%js-  D(ydim)          :jec-  D(ydim),                   cg%ks-2*dom%D_z:kec+2*dom%D_z)   &
                       + P0 * cg%prolong_x(cg%is:cg%ie, cg%js                    :jec,                             cg%ks-2*dom%D_z:kec+2*dom%D_z)   &
                       + P_1* cg%prolong_x(cg%is:cg%ie, cg%js+  D(ydim)          :jec+  D(ydim),                   cg%ks-2*dom%D_z:kec+2*dom%D_z)   &
                       + P_2* cg%prolong_x(cg%is:cg%ie, cg%js+2*D(ydim)          :jec+2*D(ydim),                   cg%ks-2*dom%D_z:kec+2*dom%D_z)
                  cg%prolong_xy(           cg%is:cg%ie, cg%js+dom%D_y            :cg%je        :2,                 cg%ks-2*dom%D_z:kec+2*dom%D_z) = &
                       + P_2* cg%prolong_x(cg%is:cg%ie, cg%js+odd(ydim)-2*D(ydim):jec-dom%D_y+odd(ydim)-2*D(ydim), cg%ks-2*dom%D_z:kec+2*dom%D_z)   &
                       + P_1* cg%prolong_x(cg%is:cg%ie, cg%js+odd(ydim)-  D(ydim):jec-dom%D_y+odd(ydim)-  D(ydim), cg%ks-2*dom%D_z:kec+2*dom%D_z)   &
                       + P0 * cg%prolong_x(cg%is:cg%ie, cg%js+odd(ydim)          :jec-dom%D_y+odd(ydim),           cg%ks-2*dom%D_z:kec+2*dom%D_z)   &
                       + P1 * cg%prolong_x(cg%is:cg%ie, cg%js+odd(ydim)+  D(ydim):jec-dom%D_y+odd(ydim)+  D(ydim), cg%ks-2*dom%D_z:kec+2*dom%D_z)   &
                       + P2 * cg%prolong_x(cg%is:cg%ie, cg%js+odd(ydim)+2*D(ydim):jec-dom%D_y+odd(ydim)+2*D(ydim), cg%ks-2*dom%D_z:kec+2*dom%D_z)
               case default
                  call die("[cg_list_level:prolong_q_1var] unsupported stencil size")
            end select

            select case (stencil_range*dom%D_z)
               case (I_ZERO)
                  cg%q(iv)%arr (            cg%is:cg%ie, cg%js:cg%je, cg%ks          :cg%ke+dom%D_z:2      ) = &
                       &      cg%prolong_xy(cg%is:cg%ie, cg%js:cg%je, cg%ks          :kec                  )
                  if (dom%has_dir(zdim)) &
                       cg%q(iv)%arr (       cg%is:cg%ie, cg%js:cg%je, cg%ks+dom%D_z  :cg%ke        :2      ) = &
                       &      cg%prolong_xy(cg%is:cg%ie, cg%js:cg%je, cg%ks+odd(zdim):kec-dom%D_z+odd(zdim))
               case (I_ONE)
                  cg%q(iv)%arr(             cg%is:cg%ie, cg%js:cg%je, cg%ks                  :cg%ke+dom%D_z:2              ) = &
                       + P1 * cg%prolong_xy(cg%is:cg%ie, cg%js:cg%je, cg%ks-D(zdim)          :kec-D(zdim)                  )   &
                       + P0 * cg%prolong_xy(cg%is:cg%ie, cg%js:cg%je, cg%ks                  :kec                          )   &
                       + P_1* cg%prolong_xy(cg%is:cg%ie, cg%js:cg%je, cg%ks+D(zdim)          :kec+D(zdim)                  )
                  cg%q(iv)%arr(             cg%is:cg%ie, cg%js:cg%je, cg%ks+dom%D_z          :cg%ke        :2              ) = &
                       + P_1* cg%prolong_xy(cg%is:cg%ie, cg%js:cg%je, cg%ks+odd(zdim)-D(zdim):kec-dom%D_z+odd(zdim)-D(zdim))   &
                       + P0 * cg%prolong_xy(cg%is:cg%ie, cg%js:cg%je, cg%ks+odd(zdim)        :kec-dom%D_z+odd(zdim)        )   &
                       + P1 * cg%prolong_xy(cg%is:cg%ie, cg%js:cg%je, cg%ks+odd(zdim)+D(zdim):kec-dom%D_z+odd(zdim)+D(zdim))
               case (I_TWO)
                  cg%q(iv)%arr(             cg%is:cg%ie, cg%js:cg%je, cg%ks                    :cg%ke+dom%D_z:2) = &
                       + P2 * cg%prolong_xy(cg%is:cg%ie, cg%js:cg%je, cg%ks-2*D(zdim)          :kec-2*D(zdim)  )   &
                       + P1 * cg%prolong_xy(cg%is:cg%ie, cg%js:cg%je, cg%ks-  D(zdim)          :kec-  D(zdim)  )   &
                       + P0 * cg%prolong_xy(cg%is:cg%ie, cg%js:cg%je, cg%ks                    :kec            )   &
                       + P_1* cg%prolong_xy(cg%is:cg%ie, cg%js:cg%je, cg%ks+  D(zdim)          :kec+  D(zdim)  )   &
                       + P_2* cg%prolong_xy(cg%is:cg%ie, cg%js:cg%je, cg%ks+2*D(zdim)          :kec+2*D(zdim)  )
                  cg%q(iv)%arr(             cg%is:cg%ie, cg%js:cg%je, cg%ks+dom%D_z            :cg%ke        :2) = &
                       + P_2* cg%prolong_xy(cg%is:cg%ie, cg%js:cg%je, cg%ks+odd(zdim)-2*D(zdim):kec-dom%D_z+odd(zdim)-2*D(zdim))   &
                       + P_1* cg%prolong_xy(cg%is:cg%ie, cg%js:cg%je, cg%ks+odd(zdim)-  D(zdim):kec-dom%D_z+odd(zdim)-  D(zdim))   &
                       + P0 * cg%prolong_xy(cg%is:cg%ie, cg%js:cg%je, cg%ks+odd(zdim)          :kec-dom%D_z+odd(zdim)          )   &
                       + P1 * cg%prolong_xy(cg%is:cg%ie, cg%js:cg%je, cg%ks+odd(zdim)+  D(zdim):kec-dom%D_z+odd(zdim)+  D(zdim))   &
                       + P2 * cg%prolong_xy(cg%is:cg%ie, cg%js:cg%je, cg%ks+odd(zdim)+2*D(zdim):kec-dom%D_z+odd(zdim)+2*D(zdim))
               case default
                  call die("[cg_list_level:prolong_q_1var] unsupported stencil size")
            end select
            ! Alternatively, an FFT convolution may be employed after injection. No idea at what stencil size the FFT is faster. It is finite size for sure :-)
         endif
         cgl => cgl%nxt
      enddo

!> \todo (connected with todo above)      call check_dirty(fine, iv, "prolong+")

   end subroutine prolong_q_1var

!> \brief Print detailed information about current level decomposition

   subroutine print_segments(this)

      use constants,  only: LO, HI
      use dataio_pub, only: printinfo, msg
      use mpisetup,   only: FIRST, LAST, master, nproc

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

!> \brief Initialize all grid containers on a new grid level

   subroutine init_all_new_cg(this)

      use cg_list_global, only: all_cg
      use grid_cont,      only: grid_container
      use mpisetup,       only: proc
      use named_array,    only: qna, wna

      implicit none

      class(cg_list_level), intent(inout) :: this

      integer :: gr_id, i
      type(grid_container), pointer :: cg

      do gr_id = lbound(this%pse(proc)%sel(:,:,:), dim=1), ubound(this%pse(proc)%sel(:,:,:), dim=1)
         call this%add
         cg => this%last%cg

         !> \todo Integrate the following few calls into grid_container type, note that mpi_bnd_types will require access to whole this%pse(:)%sel(:,:,:)
         call cg%init(this%n_d, this%pse(proc)%sel(gr_id, :, :), gr_id, this%lev) ! we cannot pass "this" as an argument because of circular dependencies
         call this%mpi_bnd_types(cg)
         call cg%set_q_mbc
         ! register all known named arrays for this cg
         if (allocated(qna%lst)) then
            do i = lbound(qna%lst(:), dim=1), ubound(qna%lst(:), dim=1)
               call cg%add_na(qna%lst(i)%multigrid)
            enddo
         endif
         if (allocated(wna%lst)) then
            do i = lbound(wna%lst(:), dim=1), ubound(wna%lst(:), dim=1)
               call cg%add_na_4d(wna%lst(i)%dim4)
            enddo
         endif

         call all_cg%add(cg)
         !> \todo add an optional argument, array of pointers to lists, where the cg should be added. Requires polymorphic array of pointers.

      enddo

!#ifdef VERBOSE
      call this%print_segments
!#endif /* VERBOSE */

      call all_cg%update_req

   end subroutine init_all_new_cg

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

      use constants,      only: FLUID, MAG, CR, ARR, xdim, zdim, ndims, LO, HI, BND, BLK, I_ONE, wcr_n
      use cg_list_global, only: all_cg
      use dataio_pub,     only: die
      use domain,         only: dom
      use grid_cont,      only: grid_container, is_overlap
      use mpi,            only: MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, MPI_COMM_NULL
      use mpisetup,       only: mpi_err, FIRST, LAST, procmask
      use named_array,    only: wna
      use types,          only: cdd

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

      if (cg%level_id /= this%lev) call die("[cg_list_level:mpi_bnd_types] Level mismatch")

      if (allocated(cg%i_bnd) .or. allocated(cg%o_bnd)) call die("[cg_list_level:mpi_bnd_types] cg%i_bnd or cg%o_bnd already allocated")
      allocate(cg%i_bnd(xdim:zdim, dom%nb), cg%o_bnd(xdim:zdim, dom%nb))

      ! There are two completely different approaches: Very general and the old one. Can be put into separate routines.
      if (cdd%comm3d == MPI_COMM_NULL) then

         ! assume that cuboids fill the domain and don't collide

         ijks(:) = cg%ijkse(:, LO) - cg%off(:)
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

         !< number of fluids, magnetic field components, CRs, and 1 for a rank-3 array
         nc(:) = I_ONE ! set at least one component, even if there is none at all
         if (all_cg%fi     > 0) nc(FLUID) = wna%lst(all_cg%fi)%dim4
         if (all_cg%bi     > 0) nc(MAG)   = wna%lst(all_cg%bi)%dim4
         if (wna%exists(wcr_n)) nc(CR)    = wna%lst(wna%ind(wcr_n))%dim4

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
                     call MPI_Type_create_subarray(dims(t), sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, cg%mbc(t, d, LO, BND, ib), mpi_err)
                     call MPI_Type_commit(cg%mbc(t, d, LO, BND, ib), mpi_err)

                     starts(dims(t)-zdim+d) = dom%nb
                     call MPI_Type_create_subarray(dims(t), sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, cg%mbc(t, d, LO, BLK, ib), mpi_err)
                     call MPI_Type_commit(cg%mbc(t, d, LO, BLK, ib), mpi_err)

                     starts(dims(t)-zdim+d) = cg%n_b(d) + dom%nb - ib
                     call MPI_Type_create_subarray(dims(t), sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, cg%mbc(t, d, HI, BLK, ib), mpi_err)
                     call MPI_Type_commit(cg%mbc(t, d, HI, BLK, ib), mpi_err)

                     starts(dims(t)-zdim+d) = cg%ijkse(d, HI)
                     call MPI_Type_create_subarray(dims(t), sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, cg%mbc(t, d, HI, BND, ib), mpi_err)
                     call MPI_Type_commit(cg%mbc(t, d, HI, BND, ib), mpi_err)

                  enddo

                  deallocate(sizes, subsizes, starts)

               enddo
            endif
         enddo

      endif

   end subroutine mpi_bnd_types

end module cg_list_lev
