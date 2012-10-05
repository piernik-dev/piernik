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

!> \brief This module extends a refinement level with connections to other levels.

module cg_level_connected

   use cg_level,      only: cg_level_T

   implicit none

   private
   public :: cg_level_connected_T, base_lev, finest, coarsest, equalize_finest

   !! \brief A list of all cg of the same resolution with links to coarser and finer levels
   type, extends(cg_level_T) :: cg_level_connected_T

      type(cg_level_connected_T), pointer :: coarser    !< coarser level cg set or null()
      type(cg_level_connected_T), pointer :: finer      !< finer level cg set or null()
      integer(kind=4) :: ord_prolong_set                !< Number of boundary cells for prolongation used in last update of cg_level_connected_T%vertical_prep

    contains

      ! Level management
      procedure, private :: init_level                      !< common initialization for base level and other levels
      procedure :: add_lev                                  !< add a finer or coarser level
      procedure :: add_lev_base                             !< initialize the base level
      generic, public :: add_level => add_lev, add_lev_base

      ! Prolongation and restriction
      procedure, private :: vertical_prep                   !< initialize prolongation and restriction targets
      procedure :: prolong                                  !< interpolate the grid data which has the flag vital set to this%finer level
      procedure :: restrict                                 !< interpolate the grid data which has the flag vital set from this%coarser level
      procedure :: prolong_q_1var                           !< interpolate the grid data in specified q field to this%finer level
      procedure :: restrict_q_1var                          !< interpolate the grid data in specified q field from this%coarser level
      procedure :: restrict_to_floor_q_1var                 !< restrict specified q field as much as possible
      procedure :: restrict_to_base_q_1var                  !< restrict specified q field to the base level
      procedure :: arr3d_boundaries                         !< Set up all guardcells (internal, external and fine-coarse) for given rank-3 arrays.
      ! fine-coarse boundary exchanges may also belong to this type
   end type cg_level_connected_T

   type(cg_level_connected_T), pointer  :: base_lev             !< base level grid containers

   !! \todo finest and coarsest together with subroutine add_lev can go to separate module
   type(cg_level_connected_T), pointer :: finest               !< finest level of refinement
   type(cg_level_connected_T), pointer :: coarsest             !< coarsest level of refinement

contains

!>
!! \brief It is not allowed for different processes to have different height of level hierarchy, so let's add some empty levels, where necessary
!!
!! \details Multigrid levels are always created on all processes, so no need to fix anything on the bottom of the hierarchy
!!
!! \todo When global top level do not have any blocks then destroy it.
!<

   subroutine equalize_finest

      use constants, only: I_ONE
      use mpi,       only: MPI_INTEGER, MPI_MAX
      use mpisetup,  only: comm, mpi_err

      implicit none

      integer(kind=4) :: g_finest_id

      call MPI_Allreduce(finest%level_id, g_finest_id, I_ONE, MPI_INTEGER, MPI_MAX, comm, mpi_err)

      do while (g_finest_id > finest%level_id)
         call finest%add_level(coarse = .false.)
      enddo

   end subroutine equalize_finest

!> \brief Common initialization for base level and other levels

   subroutine init_level(this)

      use constants, only: INVALID, fft_none

      implicit none

      class(cg_level_connected_T), intent(inout) :: this   !< object invoking type bound procedure

      this%coarser => null()
      this%finer => null()
      this%tot_se = 0
      this%ord_prolong_set = INVALID
      this%fft_type = fft_none

   end subroutine init_level

!> \brief Initialize the base level

   subroutine add_lev_base(this, n_d)

      use constants,        only: base_level_id, ndims
      use dataio_pub,       only: die
      use domain,           only: dom
      use list_of_cg_lists, only: all_lists

      implicit none

      class(cg_level_connected_T),       intent(inout) :: this   !< object invoking type bound procedure
      integer(kind=4), dimension(ndims), intent(in)    :: n_d    !< size of global base grid in cells

      if (any(n_d(:) < 1)) call die("[cg_level_connected:add_lev_base] non-positive base grid sizes")
      if (any(dom%has_dir(:) .neqv. (n_d(:) > 1))) call die("[cg_level_connected:add_lev_base] base grid size incompatible with has_dir masks")

      call this%init_level
      this%level_id = base_level_id
      this%n_d(:) = n_d(:)
      call all_lists%register(this, "Base level")

   end subroutine add_lev_base

!> \brief Add a fine or coarse level to a existing one

   subroutine add_lev(this, coarse)

      use constants,        only: INVALID, I_ONE, refinement_factor
      use dataio_pub,       only: die, msg
      use domain,           only: dom
      use mpisetup,         only: master
      use list_of_cg_lists, only: all_lists

      implicit none

      class(cg_level_connected_T), target, intent(inout) :: this    !< lowest or highest refinement level
      logical,                             intent(in)    :: coarse  !< if .true. then add a level below base level

      type(cg_level_connected_T), pointer :: new_lev !< fresh refinement level to be added

      allocate(new_lev)
      call new_lev%init_level
      new_lev%n_d(:) = 1

      if (coarse) then
         if (associated(this%coarser)) call die("[cg_level_connected:add_lev] coarser level already exists")
         select type(this)
            type is (cg_level_connected_T)
               this%coarser => new_lev
               new_lev%finer => this
            class default
               call die("[cg_level_connected:add_lev] cannot call this routine for derivatives of cg_level_connected (coarse)")
         end select

         coarsest => new_lev
         new_lev%level_id = this%level_id - I_ONE
         where (dom%has_dir(:)) new_lev%n_d(:) = this%n_d(:) / refinement_factor
         if (master .and. any(new_lev%n_d(:)*refinement_factor /= this%n_d(:) .and. dom%has_dir(:))) then
            write(msg, '(a,3f10.1,a,i3)')"[cg_level_connected:add_lev] Fractional number of domain cells: ", this%n_d(:)/real(refinement_factor), " at level ",new_lev%level_id
            call die(msg)
         endif
      else
         if (associated(this%finer)) call die("[cg_level_connected:add_lev] finer level already exists")
         select type(this)
            type is (cg_level_connected_T)
               this%finer => new_lev
               new_lev%coarser => this
            class default
               call die("[cg_level_connected:add_lev] cannot call this routine for derivatives of cg_level_connected (fine)")
         end select

         finest => new_lev
         new_lev%level_id = this%level_id + I_ONE
         where (dom%has_dir(:)) new_lev%n_d(:) = this%n_d(:) * refinement_factor
      endif

      !! make sure that vertical_prep will be called where necessary
      this%ord_prolong_set = INVALID
      write(msg, '(a,i3)')"level ",new_lev%level_id
      call all_lists%register(new_lev, msg)

   end subroutine add_lev

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

      integer(kind=8) :: bias !< prevents inconsistencies in arithmetic on integers due to rounding towards 0 (stencil_range should be more than enough)

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
      use cg_list,        only: cg_list_element
      use grid_cont,      only: grid_container, is_overlap
      use mpisetup,       only: FIRST, LAST

      implicit none

      class(cg_level_connected_T), intent(inout) :: this   !< object invoking type bound procedure

      integer :: g, j, jf, fmax, tag
      integer(kind=8), dimension(xdim:zdim, LO:HI) :: coarsened
      integer, dimension(xdim:zdim, LO:HI) :: enlargement
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg            !< current grid container
      type(cg_level_connected_T), pointer :: fine, coarse  !< shortcut
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
               fmax = fmax + size(fine%pse(j)%c(:))
            enddo
            allocate(ps(fmax))

            if (allocated(cg%ri_tgt%seg)) deallocate(cg%ri_tgt%seg) ! call warn("cll:vp cg%ri_tgt%seg a a")
            if (allocated(cg%po_tgt%seg)) deallocate(cg%po_tgt%seg)
            g = 0
            do j = FIRST, LAST
               do jf = lbound(fine%pse(j)%c(:), dim=1), ubound(fine%pse(j)%c(:), dim=1)
                  if (is_overlap(c2f(cg%my_se(:, :)), fine%pse(j)%c(jf)%se(:,:))) then
                     g = g + 1
                     ps(g) = int_pair(j, jf)
                  endif
               enddo
            enddo
            allocate(cg%ri_tgt%seg(g), cg%po_tgt%seg(g))

            !! \todo fuse the following two loops back into one
            do g = lbound(cg%ri_tgt%seg(:), dim=1), ubound(cg%ri_tgt%seg(:), dim=1)
               associate ( seg => cg%ri_tgt%seg(g) )
               if (allocated(seg%buf)) call die("cll:vp fr seg%buf a a")
               seg%proc = ps(g)%proc
               ! find cross-section of own segment with coarsened fine segment
               coarsened(:,:) = f2c(fine%pse(seg%proc)%c(ps(g)%n_se)%se(:,:))
               seg%se(:, LO) = max(cg%my_se(:, LO), coarsened(:, LO))
               seg%se(:, HI) = min(cg%my_se(:, HI), coarsened(:, HI))
               allocate(seg%buf(seg%se(xdim, HI)-seg%se(xdim, LO) + 1, &
                    &           seg%se(ydim, HI)-seg%se(ydim, LO) + 1, &
                    &           seg%se(zdim, HI)-seg%se(zdim, LO) + 1))
               tag = cg%grid_id + this%tot_se * ps(g)%n_se
               seg%tag = int(tag, kind=4) ! assumed that there is only one piece to be communicated from grid to grid (i.e. grids are not periodically wrapped around)
               if (tag /= int(seg%tag)) call die("[cg_level_connected:vertical_prep] tag overflow (ri)")
               end associate
            enddo

            ! When fine and coarse pieces are within prolongation stencil length, but there is no direct overlap we rely on guardcell update on coarse side
            ! to provide valid prolongation data. An alternative approach would not require updating guardcells, but would usually result in setting up
            ! many small transactions on the corners, especially on Cartesian decomposition.
            !
            ! When the overlap is one fine cell (half coarse cell), not all communications are necessary, but we don't want to complicate the code too much at the moment
            do g = lbound(cg%po_tgt%seg(:), dim=1), ubound(cg%po_tgt%seg(:), dim=1)
               associate ( seg => cg%po_tgt%seg(g) )
               if (allocated(seg%buf)) call die("cll:vp fp seg%buf a a")
               seg%proc = ps(g)%proc
               ! find cross-section of own segment with enlarged coarsened fine segment
               coarsened(:,:) = f2c(fine%pse(seg%proc)%c(ps(g)%n_se)%se(:,:)) + enlargement(:,:)
               seg%se(:, LO) = max(cg%my_se(:, LO) + enlargement(:, LO), coarsened(:, LO))
               seg%se(:, HI) = min(cg%my_se(:, HI) + enlargement(:, HI), coarsened(:, HI))
               allocate(seg%buf(seg%se(xdim, HI)-seg%se(xdim, LO) + 1, &
                    &           seg%se(ydim, HI)-seg%se(ydim, LO) + 1, &
                    &           seg%se(zdim, HI)-seg%se(zdim, LO) + 1))
               tag = cg%grid_id + this%tot_se * ps(g)%n_se
               seg%tag = int(tag, kind=4) ! assumed that there is only one piece to be communicated from grid to grid (i.e. grids are not periodically wrapped around)
               if (tag /= int(seg%tag)) call die("[cg_level_connected:vertical_prep] tag overflow po)")
               end associate
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
               fmax = fmax + size(coarse%pse(j)%c(:))
            enddo
            allocate(ps(fmax))

            if (allocated(cg%ro_tgt%seg)) deallocate(cg%ro_tgt%seg) ! call warn("cll:vp cg%ro_tgt%seg a a")
            if (allocated(cg%pi_tgt%seg)) deallocate(cg%pi_tgt%seg)
            g = 0
            do j = FIRST, LAST
               do jf = lbound(coarse%pse(j)%c(:), dim=1), ubound(coarse%pse(j)%c(:), dim=1)
                  if (is_overlap(cg%my_se(:, :), c2f(coarse%pse(j)%c(jf)%se(:,:)))) then
                     g = g + 1
                     ps(g) = int_pair(j, jf)
                  endif
               enddo
            enddo
            allocate(cg%ro_tgt%seg(g), cg%pi_tgt%seg(g))

            !! \todo fuse the following two loops back into one
            do g = lbound(cg%ro_tgt%seg(:), dim=1), ubound(cg%ro_tgt%seg(:), dim=1)
               associate ( seg => cg%ro_tgt%seg(g) )
               if (allocated(seg%buf)) call die("cll:vp cr seg%buf a a")
               seg%proc = ps(g)%proc
               ! find cross-section of coarsened own segment with coarse segment
               coarsened(:,:) = f2c(cg%my_se(:,:))
               seg%se(:, LO) = max(coarsened(:, LO), coarse%pse(seg%proc)%c(ps(g)%n_se)%se(:, LO))
               seg%se(:, HI) = min(coarsened(:, HI), coarse%pse(seg%proc)%c(ps(g)%n_se)%se(:, HI))
               allocate(seg%buf(seg%se(xdim, HI)-seg%se(xdim, LO) + 1, &
                    &           seg%se(ydim, HI)-seg%se(ydim, LO) + 1, &
                    &           seg%se(zdim, HI)-seg%se(zdim, LO) + 1))
               coarsened(:,:) = c2f(coarse%pse(seg%proc)%c(ps(g)%n_se)%se(:,:)) ! should be renamed to refined(:,:)
               seg%se(:, LO) = max(cg%my_se(:, LO), coarsened(:, LO))
               seg%se(:, HI) = min(cg%my_se(:, HI), coarsened(:, HI))
               tag = ps(g)%n_se + coarse%tot_se * cg%grid_id
               seg%tag = int(tag, kind=4)
               if (tag /= int(seg%tag)) call die("[cg_level_connected:vertical_prep] tag overflow (ro)")
               end associate
            enddo

            do g = lbound(cg%pi_tgt%seg(:), dim=1), ubound(cg%pi_tgt%seg(:), dim=1)
               associate ( seg => cg%pi_tgt%seg(g) )
               if (allocated(seg%buf)) call die("cll:vp cp seg%buf a a")
               seg%proc = ps(g)%proc
               ! find cross-section of coarsened own segment with enlarged coarse segment
               coarsened(:,:) = f2c(cg%my_se(:,:)) + enlargement(:,:)
               seg%se(:, LO) = max(coarsened(:, LO), coarse%pse(seg%proc)%c(ps(g)%n_se)%se(:, LO) + enlargement(:, LO))
               seg%se(:, HI) = min(coarsened(:, HI), coarse%pse(seg%proc)%c(ps(g)%n_se)%se(:, HI) + enlargement(:, HI))
               allocate(seg%buf(seg%se(xdim, HI)-seg%se(xdim, LO) + 1, &
                    &           seg%se(ydim, HI)-seg%se(ydim, LO) + 1, &
                    &           seg%se(zdim, HI)-seg%se(zdim, LO) + 1))
               tag = ps(g)%n_se + coarse%tot_se * cg%grid_id
               seg%tag = int(tag, kind=4)
               if (tag /= int(seg%tag)) call die("[cg_level_connected:vertical_prep] tag overflow (pi)")
               end associate
            enddo

            if (allocated(ps)) deallocate(ps)
            cgl => cgl%nxt
         enddo
      endif

      this%ord_prolong_set = all_cg%ord_prolong_nb

      call all_cg%update_req

   end subroutine vertical_prep

!>
!! \brief interpolate the grid data which has the flag vital set to this%finer level
!!
!! \todo implement it in a more efficient way (requires a lot more temporary buffers)
!<
   subroutine prolong(this)

      use constants,   only: base_level_id, wa_n
      use dataio_pub,  only: warn
      use named_array_list, only: qna, wna

      implicit none

      class(cg_level_connected_T), target, intent(inout) :: this !< object invoking type-bound procedure

      integer :: i, iw, iwa

      do i = lbound(qna%lst(:), dim=1), ubound(qna%lst(:), dim=1)
         if (qna%lst(i)%vital .and. (qna%lst(i)%multigrid .or. this%level_id > base_level_id)) call this%prolong_q_1var(i)
      enddo

      iwa = qna%ind(wa_n)

      do i = lbound(wna%lst(:), dim=1), ubound(wna%lst(:), dim=1)
         if (wna%lst(i)%vital .and. (wna%lst(i)%multigrid .or. this%level_id >= base_level_id)) then
            if (wna%lst(i)%multigrid) call warn("[cg_level_connected:prolong] mg set for cg%w ???")
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
!! \todo implement it in a more efficient way (requires a lot more temporary buffers)
!<
   subroutine restrict(this)

      use constants,   only: base_level_id, wa_n
      use dataio_pub,  only: warn
      use named_array_list, only: qna, wna

      implicit none

      class(cg_level_connected_T), target, intent(inout) :: this !< object invoking type-bound procedure

      integer :: i, iw, iwa

      do i = lbound(qna%lst(:), dim=1), ubound(qna%lst(:), dim=1)
         if (qna%lst(i)%vital .and. (qna%lst(i)%multigrid .or. this%level_id > base_level_id)) call this%restrict_q_1var(i)
      enddo

      iwa = qna%ind(wa_n)

      do i = lbound(wna%lst(:), dim=1), ubound(wna%lst(:), dim=1)
         if (wna%lst(i)%vital .and. (wna%lst(i)%multigrid .or. this%level_id > base_level_id)) then
            if (wna%lst(i)%multigrid) call warn("[cg_level_connected:restrict] mg set for cg%w ???")
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

      class(cg_level_connected_T), intent(inout) :: this !< object invoking type-bound procedure
      integer,                     intent(in)    :: iv   !< variable to be restricted

      if (.not. associated(this%coarser)) return
      call this%restrict_q_1var(iv)
      call this%coarser%restrict_to_floor_q_1var(iv)

   end subroutine restrict_to_floor_q_1var

!> \brief Restrict to the base level

   recursive subroutine restrict_to_base_q_1var(this, iv)

      implicit none

      class(cg_level_connected_T), intent(inout) :: this !< object invoking type-bound procedure
      integer,                     intent(in)    :: iv   !< variable to be restricted

      if (this%level_id == base_lev%level_id) return
      call this%restrict_q_1var(iv)
      call this%coarser%restrict_to_base_q_1var(iv)

   end subroutine restrict_to_base_q_1var

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
      use dataio_pub,  only: msg, warn
      use domain,      only: dom
      use cg_list,     only: cg_list_element
      use grid_cont,   only: grid_container
      use mpisetup,    only: comm, mpi_err, req, status, inflate_req
      use mpi,         only: MPI_DOUBLE_PRECISION
      use named_array, only: p3

      implicit none

      class(cg_level_connected_T), target, intent(inout) :: this !< object invoking type-bound procedure
      integer,                        intent(in)    :: iv   !< variable to be restricted

      type(cg_level_connected_T), pointer :: coarse
      integer :: g
      integer(kind=8), dimension(xdim:zdim, LO:HI) :: fse, cse ! shortcuts for fine segment and coarse segment
      integer(kind=8) :: i, j, k, ic, jc, kc
      integer(kind=8), dimension(xdim:zdim) :: off1
      real :: norm
      integer(kind=4) :: nr
      integer(kind=4), dimension(:,:), pointer :: mpistatus
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg            !< current grid container

      coarse => this%coarser
      if (.not. associated(coarse)) then ! can't restrict base level
         write(msg,'(a,i3)')"[cg_level_connected:restrict_q_1var] no coarser level than ", this%level_id
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
               if (nr > size(req, dim=1)) call inflate_req
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
            if (nr > size(req, dim=1)) call inflate_req
            call MPI_Isend(cg%ro_tgt%seg(g)%buf(1, 1, 1), size(cg%ro_tgt%seg(g)%buf(:, :, :)), MPI_DOUBLE_PRECISION, cg%ro_tgt%seg(g)%proc, cg%ro_tgt%seg(g)%tag, comm, req(nr), mpi_err)
         enddo
         cgl => cgl%nxt
      enddo

      if (nr > 0) then
         mpistatus => status(:, :nr)
         call MPI_Waitall(nr, req(:nr), mpistatus, mpi_err)
      endif

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

      use constants,        only: xdim, ydim, zdim, LO, HI, I_ZERO, I_ONE, I_TWO, BND_REF, O_INJ, O_LIN, O_D2, O_D3, O_D4, O_I2, O_I3, O_I4, refinement_factor
      use dataio_pub,       only: msg, warn, die
      use domain,           only: dom
      use cg_list,          only: cg_list_element
      use grid_cont,        only: grid_container
      use mpisetup,         only: comm, mpi_err, req, status, inflate_req
      use mpi,              only: MPI_DOUBLE_PRECISION
      use named_array_list, only: qna

      implicit none

      class(cg_level_connected_T), target, intent(inout) :: this !< object invoking type-bound procedure
      integer,                        intent(in)    :: iv   !< variable to be prolonged

      type(cg_level_connected_T), pointer :: fine
      integer :: g
      integer(kind=8), dimension(xdim:zdim, LO:HI) :: cse ! shortcuts for fine segment and coarse segment
      integer(kind=8) :: iec, jec, kec
      integer(kind=8), dimension(xdim:zdim) :: off, odd, D
      integer(kind=4) :: nr
      integer(kind=4), dimension(:, :), pointer :: mpistatus
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg            !< current grid container
      real :: P_2, P_1, P0, P1, P2
      integer :: stencil_range

      fine => this%finer
      if (.not. associated(fine)) then ! can't prolong finest level
         write(msg,'(a,i3)')"[cg_level_connected:prolong_q_1var] no finer level than: ", this%level_id
         call warn(msg)
         return
      endif

      call this%vertical_prep
      call fine%vertical_prep

      call fine%set_dirty(iv)

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
            call die("[cg_level_connected:prolong_q_1var] Unsupported order")
            return
      end select

      ! this is just for optimization. Setting stencil_range = I_TWO should work correctly for all interpolations.
      stencil_range = I_ZERO
      if (P_1 /= 0. .or. P1 /= 0.) stencil_range = I_ONE
      if (P_2 /= 0. .or. P2 /= 0.) stencil_range = I_TWO

      if (this%ord_prolong_set /= O_INJ) then
         !> \todo some variables may need special care on external boundaries
         call this%arr3d_boundaries(iv, bnd_type = BND_REF, corners = .true.) ! nb =  int(stencil_range, kind=4) ! not needed for injection
      endif
      call this%check_dirty(iv, "prolong-")

      nr = 0
      ! be ready to receive everything into right buffers
      cgl => fine%first
      do while (associated(cgl))
         cg => cgl%cg
         associate( seg => cg%pi_tgt%seg )
         if (allocated(cg%pi_tgt%seg)) then
            do g = lbound(seg(:), dim=1), ubound(seg(:), dim=1)
               nr = nr + I_ONE
               if (nr > size(req, dim=1)) call inflate_req
               call MPI_Irecv(seg(g)%buf(1, 1, 1), size(seg(g)%buf(:, :, :)), MPI_DOUBLE_PRECISION, seg(g)%proc, seg(g)%tag, comm, req(nr), mpi_err)
            enddo
         endif
         end associate
         cgl => cgl%nxt
      enddo

      ! send coarse data
      cgl => this%first
      do while (associated(cgl))
         cg => cgl%cg
         associate( seg => cg%po_tgt%seg )
         do g = lbound(seg(:), dim=1), ubound(seg(:), dim=1)

            cse(:, LO) = seg(g)%se(:,LO) - cg%off(:) + cg%ijkse(:, LO)
            cse(:, HI) = seg(g)%se(:,HI) - cg%off(:) + cg%ijkse(:, LO)

            nr = nr + I_ONE
            if (nr > size(req, dim=1)) call inflate_req
!            seg(g)%buf(:, :, :) = cg%q(iv)%span(cse)
            seg(g)%buf(:, :, :) = cg%q(iv)%arr(cse(xdim, LO):cse(xdim, HI), cse(ydim, LO):cse(ydim, HI), cse(zdim, LO):cse(zdim, HI))
            call MPI_Isend(seg(g)%buf(1, 1, 1), size(seg(g)%buf(:, :, :)), MPI_DOUBLE_PRECISION, seg(g)%proc, seg(g)%tag, comm, req(nr), mpi_err)
         enddo
         end associate
         cgl => cgl%nxt
      enddo

      if (nr > 0) then
         mpistatus => status(:, :nr)
         call MPI_Waitall(nr, req(:nr), mpistatus, mpi_err)
      endif

      ! merge received coarse data into one array and interpolate it into the right place
      cgl => fine%first
      do while (associated(cgl))
         cg => cgl%cg
         if (allocated(cg%pi_tgt%seg)) then

            do g = lbound(cg%pi_tgt%seg(:), dim=1), ubound(cg%pi_tgt%seg(:), dim=1)

               off(:) = cg%pi_tgt%seg(g)%se(:,LO) - cg%off(:)/refinement_factor + cg%ijkse(:,LO)
               cg%prolong_(off(xdim):off(xdim)+ubound(cg%pi_tgt%seg(g)%buf, dim=1)-1, &
                    &      off(ydim):off(ydim)+ubound(cg%pi_tgt%seg(g)%buf, dim=2)-1, &
                    &      off(zdim):off(zdim)+ubound(cg%pi_tgt%seg(g)%buf, dim=3)-1) = cg%pi_tgt%seg(g)%buf(:,:,:)

               !> When this%ord_prolong_set /= O_INJ, the received cg%pi_tgt%seg(:)%buf(:,:,:) may overlap
               !! The incoming data thus must either contain valid guardcells (even if qna%lst(iv)%ord_prolong == O_INJ)
               !! or the guardcells must be zeroed before sending data and received buffer should be added to cg%prolong_(:,:,:) not just assigned

            enddo

            iec = cg%is + (cg%ie - cg%is - 1)/refinement_factor + dom%D_x
            jec = cg%js + (cg%je - cg%js - 1)/refinement_factor + dom%D_y
            kec = cg%ks + (cg%ke - cg%ks - 1)/refinement_factor + dom%D_z

            !! almost all occurrences of number "2" are in fact connected to refinement_factor

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
                  call die("[cg_level_connected:prolong_q_1var] unsupported stencil size")
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
                  call die("[cg_level_connected:prolong_q_1var] unsupported stencil size")
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
                  call die("[cg_level_connected:prolong_q_1var] unsupported stencil size")
            end select
            ! Alternatively, an FFT convolution may be employed after injection. No idea at what stencil size the FFT is faster. It is finite size for sure :-)
         endif
         cgl => cgl%nxt
      enddo

      call fine%check_dirty(iv, "prolong+")

   end subroutine prolong_q_1var

!> \brief This routine sets up all guardcells (internal, external and fine-coarse) for given rank-3 arrays.

   subroutine arr3d_boundaries(this, ind, nb, area_type, bnd_type, corners)

      implicit none

      class(cg_level_connected_T), intent(in) :: this       !< the list on which to perform the boundary exchange
      integer,                     intent(in) :: ind        !< Negative value: index of cg%q(:) 3d array
      integer(kind=4), optional,   intent(in) :: nb         !< number of grid cells to exchange (not implemented for comm3d)
      integer(kind=4), optional,   intent(in) :: area_type  !< defines how do we treat boundaries
      integer(kind=4), optional,   intent(in) :: bnd_type   !< Override default boundary type on external boundaries (useful in multigrid solver).
                                                            !< Note that BND_PER, BND_MPI, BND_SHE and BND_COR aren't external and cannot be overridden
      logical,         optional,   intent(in) :: corners    !< When present and .true. then call internal_boundaries_3d for each direction separately

      call this%clear_boundaries(ind) ! Apply BND_ZERO as long as there is no coarse-to-fine interpolation
      call this%level_3d_boundaries(ind, nb, area_type, bnd_type, corners)

   end subroutine arr3d_boundaries

end module cg_level_connected
