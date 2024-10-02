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

   use cg_level, only: cg_level_t

   implicit none

   private
   public :: cg_level_connected_t

   !! \brief A list of all cg of the same resolution with links to coarser and finer levels
   type, extends(cg_level_t) :: cg_level_connected_t

      type(cg_level_connected_t), pointer :: coarser          !< coarser level cg set or null()
      type(cg_level_connected_t), pointer :: finer            !< finer level cg set or null()
      integer(kind=4)                     :: ord_prolong_set  !< Number of boundary cells for prolongation used in last update of cg_level_connected_t%vertical_prep
      logical, private                    :: need_vb_update   !< If .true. then execute vertical_b_prep

   contains

      ! Level management
      procedure :: init_level                 !< common initialization for base level and other levels
      procedure :: sync_ru                    !< Synchronize this%recently_changed and set flags for update requests
      procedure :: free_all_cg                !< Erase all data on the level, leave it empty

      ! Prolongation and restriction
      procedure :: update_verticals           !< Go through all necessary vertical_prep and vertical_b_prep
      procedure, private :: vertical_prep     !< Initialize prolongation and restriction targets
      procedure, private :: vertical_b_prep   !< Initialize prolongation targets for fine-coarse boundary exchange
      procedure, private :: vertical_bf_prep  !< Initialize prolongation targets for fine->coarse flux exchange

      procedure :: prolong                    !< interpolate the grid data which has the flag vital set to this%finer level
      procedure :: prolong_1var               !< interpolate the grid data in specified q field to this%finer level
      procedure :: prolong_bnd_from_coarser   !< Interpolate boundaries from coarse level at fine-coarse interfaces

      procedure :: restrict                   !< interpolate the grid data which has the flag vital set from this%coarser level
      procedure :: restrict_to_base           !< restrict all variables to the base level
      procedure :: restrict_to_floor_q_1var   !< restrict specified q field as much as possible
      procedure :: restrict_to_base_q_1var    !< restrict specified q field to the base level
      procedure :: restrict_to_base_w_1var    !< restrict specified w field to the base level
      procedure :: restrict_1var              !< interpolate the grid data in specified q field to this%coarser level
      procedure :: restrict_w_1var            !< interpolate the grid data in specified w field to this%coarser level

      procedure :: arr3d_boundaries           !< Set up all guardcells (internal, external and fine-coarse) for given rank-3 arrays on a single level.
      procedure :: arr4d_boundaries           !< Set up all guardcells (internal, external and fine-coarse) for given rank-4 arrays on a single level.

   end type cg_level_connected_t

contains

!> \brief Common initialization for base level and other levels

   subroutine init_level(this)

      use constants, only: INVALID, fft_none

      implicit none

      class(cg_level_connected_t), intent(inout) :: this   !< object invoking type bound procedure

      this%coarser => null()
      this%finer => null()
      this%dot%tot_se = 0
      this%ord_prolong_set = INVALID
      this%fft_type = fft_none
      this%need_vb_update = .true.
      this%recently_changed = .false. ! this level was just created, but no grids were added yet, so no update is required.
      this%dot%is_blocky = .false.
      allocate(this%l)

   end subroutine init_level

!>
!! \brief Synchronize this%recently_changed and set flags for update requests
!!
!! \details Enforce update on all levels of refinement.
!! This is a bit overkill in most cases, but if refinement happens right after domain expansion, some grid_id may change due to radical change of SFC indices.
!! After change of grid_id, some tags for vertical exchanges change (in cg_level::create after call to dot%update local), so everything requires update.
!!
!! Possible solution: keep grid_id constant (but allow recycling old values)
!! Also important: implement SFC searching for vertical communication to cut costs.
!<

   subroutine sync_ru(this)

      use allreduce, only: piernik_MPI_Allreduce
      use constants, only: pLOR, INVALID

      implicit none

      class(cg_level_connected_t), target, intent(inout) :: this

      class(cg_level_connected_t), pointer :: curl

      call piernik_MPI_Allreduce(this%recently_changed, pLOR)
      if (this%recently_changed) then
         this%need_vb_update = .true.
         this%ord_prolong_set = INVALID
         curl => this
         do while (associated(curl%finer))
            curl%finer%need_vb_update = .true.
            curl%finer%ord_prolong_set = INVALID
            curl => curl%finer
         enddo
         curl => this
         do while (associated(curl%coarser))
            curl%coarser%need_vb_update = .true.
            curl%coarser%ord_prolong_set = INVALID
            curl => curl%coarser
         enddo
         this%recently_changed = .false.
      endif

   end subroutine sync_ru

!> \brief Erase all data on the level, leave it empty

   subroutine free_all_cg(this)

      use cg_list,          only: cg_list_element
      use grid_cont,        only: grid_container
      use list_of_cg_lists, only: all_lists

      implicit none

      class(cg_level_connected_t), intent(inout) :: this

      type(cg_list_element), pointer :: cgl, aux
      type(grid_container),  pointer :: cg

      call this%cleanup
      cgl => this%first
      do while (associated(cgl))
         aux => cgl
         cgl => cgl%nxt
         cg => aux%cg
         call all_lists%forget(cg)
      enddo
      this%recently_changed = .true.
      call this%sync_ru

   end subroutine free_all_cg

!>
!! \brief Go through all necessary vertical_prep and vertical_b_prep, from botom to top
!!
!! \details This is required to be called in leaves%balance_and_update to prevent deadlocks later in fluidupdate.
!! The horizontal exchanges are usually updated via update_everything called from check_update_all.
!>

   recursive subroutine update_verticals(this)

      implicit none

      class(cg_level_connected_t), intent(inout) :: this   !< object invoking type bound procedure

      call this%vertical_prep  ! Is it necessary just here? Perhaps it has to be called somewhere anyway
      call this%vertical_b_prep

      if (associated(this%finer)) call this%finer%update_verticals

   end subroutine update_verticals

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

      use cg_cost_data,   only: I_REFINE
      use cg_list,        only: cg_list_element
      use cg_list_global, only: all_cg
      use constants,      only: xdim, ydim, zdim, LO, HI, PPP_AMR
      use dataio_pub,     only: die
      use domain,         only: dom
      use grid_cont,      only: grid_container
      use grid_helpers,   only: f2c, c2f
      use mpisetup,       only: FIRST, LAST, tag_ub
      use overlap,        only: is_overlap
      use ppp,            only: ppp_main

      implicit none

      class(cg_level_connected_t), intent(inout)   :: this   !< object invoking type bound procedure

      integer                                      :: g, jf, fmax
      integer(kind=4)                              :: j
      integer(kind=8)                              :: tag
      integer(kind=8), dimension(xdim:zdim, LO:HI) :: coarsened
      integer,         dimension(xdim:zdim, LO:HI) :: enlargement
      type(cg_list_element),      pointer          :: cgl
      type(grid_container),       pointer          :: cg            !< current grid container
      type(cg_level_connected_t), pointer          :: fine, coarse  !< shortcut
      type :: int_pair
         integer(kind=4) :: proc
         integer :: n_se
      end type int_pair
      type(int_pair), dimension(:), allocatable    :: ps
      character(len=*), parameter :: vp_label = "vertical_prep"

      call this%dot%update_tot_se

      if (all_cg%ord_prolong_nb == this%ord_prolong_set) return ! no need to update vertical communication on this level
      call ppp_main%start(vp_label, PPP_AMR)

      ! enlarge the fine blocks to allow for high orders of interpolation
      enlargement(:, LO) = -dom%D_(:) * all_cg%ord_prolong_nb
      enlargement(:, HI) =  dom%D_(:) * all_cg%ord_prolong_nb

      fine => this%finer
      ! find fine target for receiving restricted data or sending data to be prolonged
      if (associated(fine)) then
         cgl => this%first
         do while (associated(cgl))
            cg => cgl%cg
            call cg%costs%start

            if (allocated(ps)) call die("cll:vp f a ps")
            fmax = 0
            do j = FIRST, LAST
               fmax = fmax + size(fine%dot%gse(j)%c(:))
            enddo
            allocate(ps(fmax))

            if (allocated(cg%ri_tgt%seg)) deallocate(cg%ri_tgt%seg) ! call warn("cll:vp cg%ri_tgt%seg a a")
            if (allocated(cg%po_tgt%seg)) deallocate(cg%po_tgt%seg)
            g = 0
            do j = FIRST, LAST
               do jf = lbound(fine%dot%gse(j)%c(:), dim=1), ubound(fine%dot%gse(j)%c(:), dim=1)
                  if (is_overlap(c2f(cg%my_se(:, :)), fine%dot%gse(j)%c(jf)%se(:,:))) then
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
               coarsened(:,:) = f2c(fine%dot%gse(seg%proc)%c(ps(g)%n_se)%se(:,:))
               seg%se(:, LO) = max(cg%my_se(:, LO), coarsened(:, LO))
               seg%se(:, HI) = min(cg%my_se(:, HI), coarsened(:, HI))
               allocate(seg%buf(seg%se(xdim, HI)-seg%se(xdim, LO) + 1, &
                    &           seg%se(ydim, HI)-seg%se(ydim, LO) + 1, &
                    &           seg%se(zdim, HI)-seg%se(zdim, LO) + 1))
               tag = cg%grid_id + this%dot%tot_se * ps(g)%n_se
               seg%tag = int(tag, kind=4) ! assumed that there is only one piece to be communicated from grid to grid (i.e. grids are not periodically wrapped around)
               if (tag /= seg%tag .or. tag < 0 .or. tag > tag_ub) call die("[cg_level_connected:vertical_prep] tag overflow (ri)")
               nullify(seg%local)
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
               coarsened(:,:) = f2c(fine%dot%gse(seg%proc)%c(ps(g)%n_se)%se(:,:)) + enlargement(:,:)
               seg%se(:, LO) = max(cg%my_se(:, LO) + enlargement(:, LO), coarsened(:, LO))
               seg%se(:, HI) = min(cg%my_se(:, HI) + enlargement(:, HI), coarsened(:, HI))
               allocate(seg%buf(seg%se(xdim, HI)-seg%se(xdim, LO) + 1, &
                    &           seg%se(ydim, HI)-seg%se(ydim, LO) + 1, &
                    &           seg%se(zdim, HI)-seg%se(zdim, LO) + 1))
               tag = cg%grid_id + this%dot%tot_se * ps(g)%n_se
               seg%tag = int(tag, kind=4) ! assumed that there is only one piece to be communicated from grid to grid (i.e. grids are not periodically wrapped around)
               if (tag /= seg%tag .or. tag < 0 .or. tag > tag_ub) call die("[cg_level_connected:vertical_prep] tag overflow (po))")
               nullify(seg%local)
               end associate
            enddo

            if (allocated(ps)) deallocate(ps)

            call cg%update_leafmap

            call cg%costs%stop(I_REFINE)
            cgl => cgl%nxt
         enddo
      endif

      ! find coarse target for sending restricted data or receiving data to be prolonged
      !> \deprecated almost duplicated code
      coarse => this%coarser
      if (associated(coarse)) then

         call coarse%dot%update_tot_se

         cgl => this%first
         do while (associated(cgl))
            cg => cgl%cg
            call cg%costs%start

            if (allocated(ps)) call die("cll:vp c a ps")
            fmax = 0
            do j = FIRST, LAST
               fmax = fmax + size(coarse%dot%gse(j)%c(:))
            enddo
            allocate(ps(fmax))

            if (allocated(cg%ro_tgt%seg)) deallocate(cg%ro_tgt%seg) ! call warn("cll:vp cg%ro_tgt%seg a a")
            if (allocated(cg%pi_tgt%seg)) deallocate(cg%pi_tgt%seg)
            g = 0
            do j = FIRST, LAST
               do jf = lbound(coarse%dot%gse(j)%c(:), dim=1), ubound(coarse%dot%gse(j)%c(:), dim=1)
                  if (is_overlap(cg%my_se(:, :), c2f(coarse%dot%gse(j)%c(jf)%se(:,:)))) then
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
               seg%se(:, LO) = max(coarsened(:, LO), coarse%dot%gse(seg%proc)%c(ps(g)%n_se)%se(:, LO))
               seg%se(:, HI) = min(coarsened(:, HI), coarse%dot%gse(seg%proc)%c(ps(g)%n_se)%se(:, HI))
               allocate(seg%buf(seg%se(xdim, HI)-seg%se(xdim, LO) + 1, &
                    &           seg%se(ydim, HI)-seg%se(ydim, LO) + 1, &
                    &           seg%se(zdim, HI)-seg%se(zdim, LO) + 1))
               coarsened(:,:) = c2f(coarse%dot%gse(seg%proc)%c(ps(g)%n_se)%se(:,:)) ! should be renamed to refined(:,:)
               seg%se(:, LO) = max(cg%my_se(:, LO), coarsened(:, LO))
               seg%se(:, HI) = min(cg%my_se(:, HI), coarsened(:, HI))
               tag = ps(g)%n_se + coarse%dot%tot_se * cg%grid_id
               seg%tag = int(tag, kind=4)
               if (tag /= seg%tag .or. tag < 0 .or. tag > tag_ub) call die("[cg_level_connected:vertical_prep] tag overflow (ro)")
               nullify(seg%local)
               end associate
            enddo

            do g = lbound(cg%pi_tgt%seg(:), dim=1), ubound(cg%pi_tgt%seg(:), dim=1)
               associate ( seg => cg%pi_tgt%seg(g) )
               if (allocated(seg%buf)) call die("cll:vp cp seg%buf a a")
               seg%proc = ps(g)%proc
               ! find cross-section of coarsened own segment with enlarged coarse segment
               coarsened(:,:) = f2c(cg%my_se(:,:)) + enlargement(:,:)
               seg%se(:, LO) = max(coarsened(:, LO), coarse%dot%gse(seg%proc)%c(ps(g)%n_se)%se(:, LO) + enlargement(:, LO))
               seg%se(:, HI) = min(coarsened(:, HI), coarse%dot%gse(seg%proc)%c(ps(g)%n_se)%se(:, HI) + enlargement(:, HI))
               allocate(seg%buf(seg%se(xdim, HI)-seg%se(xdim, LO) + 1, &
                    &           seg%se(ydim, HI)-seg%se(ydim, LO) + 1, &
                    &           seg%se(zdim, HI)-seg%se(zdim, LO) + 1))
               tag = ps(g)%n_se + coarse%dot%tot_se * cg%grid_id
               seg%tag = int(tag, kind=4)
               if (tag /= seg%tag .or. tag < 0 .or. tag > tag_ub) call die("[cg_level_connected:vertical_prep] tag overflow (pi)")
               nullify(seg%local)
               end associate
            enddo

            if (allocated(ps)) deallocate(ps)

            call cg%costs%stop(I_REFINE)
            cgl => cgl%nxt
         enddo
      endif

      this%ord_prolong_set = all_cg%ord_prolong_nb

      call all_cg%update_req
      call ppp_main%stop(vp_label, PPP_AMR)

   end subroutine vertical_prep

!>
!! \brief Initialize prolongation targets for fine-coarse boundary exchange
!!
!! \details This routine allows prolongation of guardcells only from interior cells of coarser level.
!! It means that no prolongation will happen on a 2-level step, even if it would be possible to interpolate on coarser level guardcells.
!! While it is possible to implement such data transport, we don't want to rely on double- or multiple-times interpolated data.
!!
!! OPT: As prolong_bnd_from_coarser calls coarse%level_3d_boundaries in some cases (except for all_cg%ord_prolong_nb == O_INJ) it should be possible to reduce number
!! of MPI messages by relying on corner data on the coarse level (corners can be obtained along with face values). Be careful as this may break things in vertical_bf_prep.
!!
!! \todo implement local copies without MPI (provide appropriate pointers)
!<

   subroutine vertical_b_prep(this)

      use cg_cost_data,   only: I_REFINE
      use cg_list,        only: cg_list_element
      use cg_list_global, only: all_cg
      use constants,      only: xdim, ydim, zdim, LO, HI, I_ZERO, I_ONE, ndims, INVALID, base_level_id
      use dataio_pub,     only: warn, msg, die
      use domain,         only: dom
      use grid_cont,      only: grid_container
      use grid_helpers,   only: f2c
      use mergebox,       only: wmap  ! this is the last place that uses this module
      use MPIF,           only: MPI_INTEGER, MPI_INTEGER8, MPI_COMM_WORLD
      use MPIFUN,         only: MPI_Alltoall, MPI_Alltoallv
      use mpisetup,       only: FIRST, LAST, err_mpi, proc
      use overlap,        only: is_overlap
      use tag_pool,       only: t_pool

      implicit none

      class(cg_level_connected_t), intent(inout), target :: this !< the list on which to update connectivity data for fine-coarse boundary exchange

      type(cg_level_connected_t), pointer :: coarse
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg            !< current grid container
      integer :: d, b, rp, dd, ix, iy, iz
      integer(kind=4) :: j, ls
      integer(kind=8), dimension(xdim:zdim, LO:HI) :: seg, segp, seg2, segp2, segf
      integer(kind=8), dimension(xdim:zdim) :: per, ext_buf
      integer :: mpifc_cnt
      integer(kind=4) :: tag, tag_min, tag_max
      type :: fc_seg !< the absolutely minimal set of data that defines the communication consists of [ grid_id, tag, and, seg ]. The proc numbers are for convenience only.
         integer(kind=4) :: proc     ! it can be rewritten in a way that does not need this number to be explicitly stored, but it is easier to have it
         integer :: grid_id
         integer(kind=4) :: tag
         integer(kind=4) :: src_proc ! this can be computed on destination process, but it is easier to have it
         integer(kind=8), dimension(xdim:zdim, LO:HI) :: seg
         integer(kind=8), dimension(xdim:zdim, LO:HI) :: seg2
         integer(kind=8), dimension(xdim:zdim, LO:HI) :: fse
         integer(kind=8), dimension(xdim:zdim, LO:HI) :: fse2
      end type fc_seg
      enum, bind(C) !< set of constants for converting between fc_seg type and crude integer(kind=8) tables
         enumerator :: I_PROC = 1
         enumerator :: I_GRID
         enumerator :: I_TAG
         enumerator :: I_SRC_PROC
         enumerator :: I_SEG
         enumerator :: I_SEG2 = I_SEG  + HI*zdim
         enumerator :: I_LAST = I_SEG2 + HI*zdim - I_ONE
         ! no need to create I_FSE at the moment
      end enum
      type(fc_seg), allocatable, dimension(:) :: seglist, tmp
      integer(kind=4), dimension(FIRST:LAST) :: pscnt, prcnt, psdispl, prdispl, psind
      integer(kind=4), dimension(FIRST:LAST) :: sendcounts, sdispls, recvcounts, rdispls
      integer(kind=8), allocatable, dimension(:) :: sseg, rseg
      type(wmap) :: cgmap
      integer(kind=8), dimension(ndims, LO:HI)  :: box_8   !< temporary storage
      logical :: found_flux
      integer, parameter :: initial_size = 16 ! for seglist
      real, parameter :: grow_ratio = 2.      ! for seglist
      integer(kind=4) :: isl                          ! current position in seglist

      if (.not. this%need_vb_update) return

      coarse => this%coarser
      if (.not. associated(coarse)) then
         this%need_vb_update = .false.
         return
      endif
      if (coarse%l%id < base_level_id) then  ! no need to calculate f/c below base level
         this%need_vb_update = .false.
         return
      endif

      ext_buf = dom%D_ * all_cg%ord_prolong_nb ! extension of the buffers due to stencil range

      ! define areas on the fine side at BND_FC and BND_MPI_FC faces that require coarse data
      per(:) = 0
      where (dom%periodic(:)) per(:) = coarse%l%n_d(:)

      call t_pool%release(this%l%id)
      call t_pool%get(this%l%id, tag_min, tag_max)
      tag = tag_min
      mpifc_cnt = 0
      allocate(seglist(initial_size))
      isl = I_ZERO
      cgl => this%first
      do while (associated(cgl))
         cg => cgl%cg
         call cgl%cg%costs%start

         ls = isl

         box_8 = int(cg%lhn, kind=8)
         call cgmap%init(box_8)
         call cgmap%set
         box_8 = int(cg%ijkse, kind=8)
         call cgmap%clear(box_8)
         do dd = lbound(cg%i_bnd, dim=1), ubound(cg%i_bnd, dim=1)
            if (allocated(cg%i_bnd(dd)%seg)) then
               do b = lbound(cg%i_bnd(dd)%seg, dim=1), ubound(cg%i_bnd(dd)%seg, dim=1)
                  call cgmap%clear(cg%i_bnd(dd)%seg(b)%se)
               enddo
            endif
         enddo
         call cgmap%find_boxes
         do j = FIRST, LAST !> \warning Antiparallel
            do b = lbound(coarse%dot%gse(j)%c(:), dim=1), ubound(coarse%dot%gse(j)%c(:), dim=1)
               do d = lbound(cgmap%blist%blist, dim=1), ubound(cgmap%blist%blist, dim=1)
                  if (is_overlap(f2c(cgmap%blist%blist(d)%b), coarse%dot%gse(j)%c(b)%se(:,:), per(:))) then
                     do iz = -1, 1 ! scan through all periodic possibilities
                        if (iz == 0 .or. per(zdim)>0) then
                           do iy = -1, 1
                              if (iy == 0 .or. per(ydim)>0) then
                                 do ix = -1, 1
                                    if (ix == 0 .or. per(xdim)>0) then
                                       seg = f2c(cgmap%blist%blist(d)%b)
                                       seg(:, LO) = seg(:, LO) + [ ix, iy, iz ] * per(:)
                                       seg(:, HI) = seg(:, HI) + [ ix, iy, iz ] * per(:)  ! try variants corrected for periodicity
                                       if (is_overlap(coarse%dot%gse(j)%c(b)%se(:,:), seg)) then
                                          seg(:, LO) = max( seg(:, LO), coarse%dot%gse(j)%c(b)%se(:, LO))
                                          seg(:, HI) = min( seg(:, HI), coarse%dot%gse(j)%c(b)%se(:, HI)) ! this is what we want
                                          tag = tag + I_ONE
                                          if (tag > tag_max) then
                                             call t_pool%get(this%l%id, tag_min, tag_max)
                                             tag = tag_min
                                          endif
                                          if (tag<0) call die("[cg_level_connected:vertical_b_prep] tag overflow")
                                          segp (:, LO) = seg  (:, LO) - [ ix, iy, iz ] * per(:)
                                          segp (:, HI) = seg  (:, HI) - [ ix, iy, iz ] * per(:)
                                          ! Find 1-layer thick areas which will be involved in fine->coarse flux exchanges
                                          segp2(:, LO) = seg  (:, LO) - [ ix, iy, iz ] * per(:)
                                          segp2(:, HI) = seg  (:, HI) - [ ix, iy, iz ] * per(:)
                                          found_flux = .false.
                                          do dd = xdim, zdim
                                             if (dom%has_dir(dd)) then
                                                segf = cg%my_se
                                                segf(dd, :) = segf(dd, :) + [ -1, 1 ]
                                                segf = f2c(segf)
                                                if (is_overlap(segf,segp2)) then
                                                   if (found_flux) call die("[cg_level_connected:vertical_b_prep] matches multiple directions")
                                                   segf(:, LO) = max(segf(:, LO), segp2(:, LO))
                                                   segf(:, HI) = min(segf(:, HI), segp2(:, HI))
                                                   segp2 = segf
                                                   found_flux = .true.
                                                endif
                                             endif
                                          enddo
                                          if (.not. found_flux) then
                                             segp2(:, LO) = this%l%off(:)
                                             segp2(:, HI) = this%l%off(:) - dom%D_(:)
                                          endif
                                          seg2 (:, LO) = segp2(:, LO) + [ ix, iy, iz ] * per(:)
                                          seg2 (:, HI) = segp2(:, HI) + [ ix, iy, iz ] * per(:)
                                          isl = isl + I_ONE
                                          if (isl > ubound(seglist, dim=1)) then
                                             allocate(tmp(lbound(seglist(:),dim=1):int(abs(grow_ratio*ubound(seglist(:), dim=1)))))
                                             tmp(:ubound(seglist(:), dim=1)) = seglist(:)
                                             tmp(ubound(seglist(:), dim=1)+1:)%proc = INVALID
                                             call move_alloc(from=tmp, to=seglist)
                                          endif
                                          seglist(isl) = fc_seg(j, b, tag, proc, seg, seg2, segp, segp2)
                                       endif
                                    endif
                                 enddo
                              endif
                           enddo
                        endif
                     enddo
                  endif
               enddo
            enddo
         enddo

         if (allocated(cg%pib_tgt%seg)) deallocate(cg%pib_tgt%seg)
         allocate(cg%pib_tgt%seg(isl-ls))
         do j = ls + I_ONE, isl
            associate ( se => cg%pib_tgt%seg(j-ls) )
            se%proc = seglist(j)%proc
            se%tag  = seglist(j)%tag
            se%se   = seglist(j)%fse
            se%se2  = seglist(j)%fse2
            allocate(se%buf(seglist(j)%seg(xdim, HI)-seglist(j)%seg(xdim, LO) + 1 + 2*ext_buf(xdim), &
                    &       seglist(j)%seg(ydim, HI)-seglist(j)%seg(ydim, LO) + 1 + 2*ext_buf(ydim), &
                    &       seglist(j)%seg(zdim, HI)-seglist(j)%seg(zdim, LO) + 1 + 2*ext_buf(zdim)))
            end associate
         enddo

         call cgmap%cleanup

         call cgl%cg%costs%stop(I_REFINE)
         cgl => cgl%nxt
      enddo
      if (mpifc_cnt > 0) call warn("[cg_level_connected:vertical_b_prep] mixed MPI and fine-coarse boundary will be treated as FC here and then fixed in intenal_boundaries")
      !> \warning OPT Try to exclude parts that will be overwritten by guardcell exchange on the same level

      pscnt = 0
      if (isl > 0) then
         do j = lbound(seglist, dim=1, kind=4), isl
            pscnt(seglist(j)%proc) = pscnt(seglist(j)%proc) + I_ONE
         enddo
      endif

      ! communicate to the processes with coarse data how many segments are required
      call MPI_Alltoall(pscnt, I_ONE, MPI_INTEGER, prcnt, I_ONE, MPI_INTEGER, MPI_COMM_WORLD, err_mpi)

      psdispl(FIRST) = 0; prdispl(FIRST) = 0
      do j = FIRST+1, LAST
         psdispl(j) = psdispl(j-1) + pscnt(j-1)
         prdispl(j) = prdispl(j-1) + prcnt(j-1)
      enddo

      allocate(sseg(I_LAST*sum(pscnt)))
      allocate(rseg(I_LAST*sum(prcnt)))
      psind = 0
      do j = lbound(seglist, dim=1, kind=4), isl
         b = (psdispl(seglist(j)%proc) + psind(seglist(j)%proc)) * I_LAST
         sseg(b+I_PROC:b+I_LAST) = [ int(seglist(j)%proc, kind=8), int(seglist(j)%grid_id, kind=8), int(seglist(j)%tag, kind=8), int(seglist(j)%src_proc, kind=8), seglist(j)%seg, seglist(j)%seg2 ]
         psind(seglist(j)%proc) = psind(seglist(j)%proc) + I_ONE
      enddo
      ! communicate to the processes with coarse data the segments that are required
      sendcounts = I_LAST * pscnt
      sdispls = I_LAST * psdispl
      recvcounts = I_LAST * prcnt
      rdispls = I_LAST * prdispl
      ! OPT: this call can be quite long to complete
      call MPI_Alltoallv(sseg, sendcounts, sdispls, MPI_INTEGER8, rseg, recvcounts, rdispls, MPI_INTEGER8, MPI_COMM_WORLD, err_mpi)

      ! define areas on the coarse side at fine BND_FC and BND_MPI_FC faces that have to be sent
      cgl => coarse%first
      do while (associated(cgl))
         cg => cgl%cg
         rp = FIRST
         b = 0
         do j = 1, sum(prcnt)
            if (rp < LAST) then
               do while (prdispl(rp+1) < j)  ! update remote process number
                  rp = rp + 1
                  if (rp == LAST) exit
               enddo
            endif
            if (rseg(I_GRID + (j-1)*I_LAST) == cg%grid_id) then
               if (rseg(I_PROC + (j-1)*I_LAST) /= proc) call warn("[cg_level_connected:vertical_b_prep] I_PROC /= proc")
               if (rseg(I_SRC_PROC + (j-1)*I_LAST) /= rp) then
                  write(msg,*)"[cg_level_connected:vertical_b_prep] I_SRC_PROC /= rp @",proc," @@",rp, rseg(I_SRC_PROC + (j-1)*I_LAST)
                  call warn(msg)
               endif
               ! call warn("[cg_level_connected:vertical_b_prep] I_TAG: visited twice")
               b = b + 1
            endif
         enddo
         if (allocated(cg%pob_tgt%seg)) deallocate(cg%pob_tgt%seg)
         allocate(cg%pob_tgt%seg(b))
         b = 0
         do j = 1, sum(prcnt)
            if (rseg(I_GRID + (j-1)*I_LAST) == cg%grid_id) then
               b = b + 1
               associate ( se => cg%pob_tgt%seg(b) )
               se%proc = int(rseg(I_SRC_PROC      + (j-1)*I_LAST), kind=4)
               se%tag  = int(rseg(I_TAG           + (j-1)*I_LAST), kind=4)
               se%se(:, LO) = rseg(I_SEG          + (j-1)*I_LAST:I_SEG  +   zdim - xdim + (j-1)*I_LAST)
               se%se(:, HI) = rseg(I_SEG + zdim   + (j-1)*I_LAST:I_SEG  + 2*zdim - xdim + (j-1)*I_LAST)
               se%se2(:, LO) = rseg(I_SEG2        + (j-1)*I_LAST:I_SEG2 +   zdim - xdim + (j-1)*I_LAST)
               se%se2(:, HI) = rseg(I_SEG2 + zdim + (j-1)*I_LAST:I_SEG2 + 2*zdim - xdim + (j-1)*I_LAST)
               allocate(se%buf(se%se(xdim, HI)-se%se(xdim, LO) + 1 + 2*ext_buf(xdim), &
                    &          se%se(ydim, HI)-se%se(ydim, LO) + 1 + 2*ext_buf(ydim), &
                    &          se%se(zdim, HI)-se%se(zdim, LO) + 1 + 2*ext_buf(zdim)))
               end associate
            endif
         enddo

         cgl => cgl%nxt
      enddo

      deallocate(sseg)
      deallocate(rseg)

      deallocate(seglist)

      call this%vertical_bf_prep

      this%need_vb_update = .false.

   end subroutine vertical_b_prep

!>
!! \brief Initialize targets for fine->coarse flux exchange
!!
!! \details This routine relies on the data generated by vertical_b_prep. Rejects corner locations, leaves only face layer.
!!
!! The problem with dom%nb <= 2 can be solved but practically we won't ever use such thin buffer of guardcells.
!<

   subroutine vertical_bf_prep(this)

      use cg_cost_data,     only: I_REFINE
      use cg_list,          only: cg_list_element
      use constants,        only: LO, HI, pdims, ORTHO1, ORTHO2, xdim, zdim, psidim, mag_n
      use dataio_pub,       only: die
      use domain,           only: dom
      use fluidindex,       only: flind
      use grid_helpers,     only: c2f
      use named_array_list, only: wna

      implicit none

      class(cg_level_connected_t), intent(inout) :: this !< the list on which to update connectivity data for fine->coarse flux exchange

      type(cg_level_connected_t), pointer :: coarse
      type(cg_list_element), pointer :: cgl
      integer :: g, d, dd, lh, fc_fluxes

      if (dom%nb <= 2) call die("[cg_level_connected:vertical_bf_prep] dom%nb <= 2 leads to incorrect calculations")

      fc_fluxes = flind%all
      if (wna%exists(mag_n)) fc_fluxes = fc_fluxes + psidim

      cgl => this%first
      do while (associated(cgl))
         call cgl%cg%costs%start

         associate ( seg => cgl%cg%pib_tgt%seg )
         if (allocated(cgl%cg%rof_tgt%seg)) deallocate(cgl%cg%rof_tgt%seg)
         ! outside means invalid (no f/c boundary to exchange fluxes)
         do dd = xdim, zdim
            cgl%cg%coarsebnd(dd, LO)%index(:, :) = cgl%cg%ijkse(dd, LO) - 1
            cgl%cg%coarsebnd(dd, HI)%index(:, :) = cgl%cg%ijkse(dd, HI) + 1
         enddo
         if (allocated(cgl%cg%pib_tgt%seg)) then
            d = 0
            do g = lbound(seg(:), dim=1), ubound(seg(:), dim=1)
               if (all(seg(g)%se2(:, HI)>=seg(g)%se2(:, LO))) d = d + 1
            enddo
            if (allocated(cgl%cg%rof_tgt%seg)) deallocate(cgl%cg%rof_tgt%seg)
            allocate(cgl%cg%rof_tgt%seg(d))
            d = 0
            do g = lbound(seg(:), dim=1), ubound(seg(:), dim=1)
               if (all(seg(g)%se2(:, HI)>=seg(g)%se2(:, LO))) then
                  dd = guess_dir(seg(g)%se2)
                  d = d + 1
                  associate( seg2 => cgl%cg%rof_tgt%seg )
                  seg2(d)%proc = seg(g)%proc
                  seg2(d)%tag  = seg(g)%tag
                  seg2(d)%se   = seg(g)%se2
                  allocate(seg2(d)%buf(fc_fluxes, seg2(d)%se(pdims(dd, ORTHO1), LO):seg2(d)%se(pdims(dd, ORTHO1), HI), &
                       &                          seg2(d)%se(pdims(dd, ORTHO2), LO):seg2(d)%se(pdims(dd, ORTHO2), HI)))
                  if (seg(g)%se(dd, LO) == seg2(d)%se(dd, LO)) then  ! This doesn't work correctly when dom%nb == 2 or less
                     lh = HI
                  else
                     lh = LO
                  endif
                  seg2(d)%se(dd, lh) = seg2(d)%se(dd, lh) + LO + HI - 2*lh
                  seg2(d)%se = c2f(seg2(d)%se)
                  seg2(d)%se(dd, HI + LO - lh) = seg2(d)%se(dd, lh)
                  cgl%cg%coarsebnd(dd, lh)%index(seg2(d)%se(pdims(dd, ORTHO1), LO):seg2(d)%se(pdims(dd, ORTHO1), HI), &
                       &                         seg2(d)%se(pdims(dd, ORTHO2), LO):seg2(d)%se(pdims(dd, ORTHO2), HI)) = int(seg2(d)%se(dd, lh))
                  end associate

               endif
            enddo
         endif
         end associate

         call cgl%cg%costs%stop(I_REFINE)
         cgl => cgl%nxt
      enddo

      coarse => this%coarser
      if (associated(coarse)) then
         cgl => coarse%first
         do while (associated(cgl))
            call cgl%cg%costs%start

            if (allocated(cgl%cg%rif_tgt%seg)) deallocate(cgl%cg%rif_tgt%seg)
            do dd = xdim, zdim
               ! outside means invalid (no f/c boundary to exchange fluxes)
               cgl%cg%finebnd(dd, LO)%index(:, :) = cgl%cg%ijkse(dd, LO) - 1
               cgl%cg%finebnd(dd, HI)%index(:, :) = cgl%cg%ijkse(dd, HI) + 1
            enddo
            associate( seg => cgl%cg%pob_tgt%seg )
            if (allocated(cgl%cg%pob_tgt%seg)) then
               d = 0
               do g = lbound(seg(:), dim=1), ubound(seg(:), dim=1)
                  if (all(seg(g)%se2(:, HI)>=seg(g)%se2(:, LO))) d = d + 1
               enddo
               if (allocated(cgl%cg%rif_tgt%seg)) deallocate(cgl%cg%rif_tgt%seg)
               allocate(cgl%cg%rif_tgt%seg(d))
               d = 0
               do g = lbound(seg(:), dim=1), ubound(seg(:), dim=1)
                  if (all(seg(g)%se2(:, HI)>=seg(g)%se2(:, LO))) then
                     dd = guess_dir(seg(g)%se2)
                     d = d + 1
                     associate( seg2 => cgl%cg%rif_tgt%seg )
                     seg2(d)%proc = seg(g)%proc
                     seg2(d)%tag  = seg(g)%tag
                     seg2(d)%se   = seg(g)%se2
                     allocate(seg2(d)%buf(fc_fluxes, seg2(d)%se(pdims(dd, ORTHO1), LO):seg2(d)%se(pdims(dd, ORTHO1), HI), &
                          &                          seg2(d)%se(pdims(dd, ORTHO2), LO):seg2(d)%se(pdims(dd, ORTHO2), HI)))
                     if (seg(g)%se(dd, LO) == seg2(d)%se(dd, LO)) then
                        lh = HI
                     else
                        lh = LO
                     endif
                     cgl%cg%finebnd(dd, lh)%index(seg2(d)%se(pdims(dd, ORTHO1), LO):seg2(d)%se(pdims(dd, ORTHO1), HI), &
                          &                       seg2(d)%se(pdims(dd, ORTHO2), LO):seg2(d)%se(pdims(dd, ORTHO2), HI)) = int(seg2(d)%se(dd, LO)) ! ? +1 for HI
                     end associate
                  endif
               enddo
            endif
            end associate

            call cgl%cg%costs%stop(I_REFINE)
            cgl => cgl%nxt
         enddo
      endif

   contains

      integer function guess_dir(se) result(dir)

         use constants,  only: LO, HI, xdim, zdim, INVALID
         use dataio_pub, only: die
         use domain,     only: dom

         implicit none

         integer(kind=8), dimension(xdim:zdim, LO:HI) :: se

         integer :: d

         dir = INVALID
         do d = xdim, zdim
            if (dom%has_dir(d) .and. se(d, HI) == se(d, LO)) then
               if (dir /= INVALID) call die("[cg_level_connected:vertical_bf_prep:guess_dir] point-like?")
               dir = d
            endif
         enddo
         if (dir == INVALID) call die("[cg_level_connected:vertical_bf_prep:guess_dir] undefined direction?")

      end function guess_dir

   end subroutine vertical_bf_prep

!>
!! \brief Interpolate the grid data which has the flag vital set to this%finer level
!!
!! Here we can look for particles that became covered by finer levels and send them to the child cgs.
!<

   subroutine prolong(this, bnd_type)

      use constants,        only: base_level_id, GEO_RPZ, PPP_AMR
      use dataio_pub,       only: warn
      use domain,           only: dom
      use fluidindex,       only: iarr_all_my
      use named_array_list, only: qna, wna
      use ppp,              only: ppp_main

      implicit none

      class(cg_level_connected_t), target, intent(inout) :: this       !< object invoking type-bound procedure
      integer(kind=4), optional,           intent(in)    :: bnd_type   !< Override default boundary type on external boundaries (useful in multigrid solver).

      integer(kind=4) :: i, iw
      character(len=*), parameter :: proq_label = "prolong_qna", prow_label = "prolong_wna"

      call ppp_main%start(proq_label, PPP_AMR)
      do i = lbound(qna%lst(:), dim=1, kind=4), ubound(qna%lst(:), dim=1, kind=4)
         if (qna%lst(i)%vital .and. (qna%lst(i)%multigrid .or. this%l%id >= base_level_id)) call this%prolong_1var(i, bnd_type = bnd_type)
         ! Although we aren't worried too much by nonconservation of psi or multigrid fields here
         ! but it will be worth checking if conservative high-order prolongation can help.
      enddo
      call ppp_main%stop(proq_label, PPP_AMR)

      call ppp_main%start(prow_label, PPP_AMR)
      do i = lbound(wna%lst(:), dim=1, kind=4), ubound(wna%lst(:), dim=1, kind=4)
         if (wna%lst(i)%vital .and. (wna%lst(i)%multigrid .or. this%l%id >= base_level_id)) then

            if (dom%geometry_type == GEO_RPZ .and. i == wna%fi) then  ! take the slow way
               qna%lst(qna%wai)%ord_prolong = 0 !> \todo implement high order conservative prolongation and use wna%lst(i)%ord_prolong here
               if (wna%lst(i)%multigrid) call warn("[cg_level_connected:prolong] mg set for cg%w ???")
               do iw = 1, wna%lst(i)%dim4
                  call this%wq_copy(i, iw, qna%wai)
                  if (dom%geometry_type == GEO_RPZ .and. i == wna%fi .and. any(iw == iarr_all_my)) call this%mul_by_r(qna%wai) ! angular momentum conservation
                  if (.true.) then  !> Quick and dirty fix for cases when cg%ignore_prolongation == .true.
                     call this%finer%wq_copy(i, iw, qna%wai)
                     if (dom%geometry_type == GEO_RPZ .and. i == wna%fi .and. any(iw == iarr_all_my)) call this%finer%mul_by_r(qna%wai)
                  endif
                  call this%prolong_1var(qna%wai, wna%lst(i)%position(iw), bnd_type = bnd_type)
                  if (dom%geometry_type == GEO_RPZ .and. i == wna%fi .and. any(iw == iarr_all_my)) call this%finer%div_by_r(qna%wai) ! angular momentum conservation
                  call this%finer%qw_copy(qna%wai, i, iw) !> \todo filter this through cg%ignore_prolongation
               enddo
            else
               call this%prolong_1var(i, dim4 = .true., bnd_type = bnd_type)
            endif
         endif
      enddo
      call ppp_main%stop(prow_label, PPP_AMR)

   end subroutine prolong

!>
!! \brief Perform prolongation of one rank-3 or rank-4 array.
!!
!! \details This routine communicates selected named array from coarse to fine grid.
!! The prolonged data is then copied to the destination if the cg%ignore_prolongation allows it.
!!
!! OPT: Find a way to prolong only what is really needed (some communication can be skipped).
!!
!! OPT Usually there are many messages that are sent between the same pairs of processes
!! \todo Sort all messages according to e.g. tag and send/receive aggregated message with everything
!! \todo implement local copies without MPI
!<

   subroutine prolong_1var(this, iv, pos, bnd_type, dim4)

      use cg_cost_data,     only: I_REFINE
      use cg_list,          only: cg_list_element
      use constants,        only: xdim, ydim, zdim, LO, HI, I_ONE, I_ZERO, VAR_CENTER, ndims, PPP_AMR  !, dirtyH1
      use dataio_pub,       only: msg, warn
      use grid_cont,        only: grid_container
      use grid_helpers,     only: f2c, c2f
      use mpisetup,         only: err_mpi, req, inflate_req, master
      use MPIF,             only: MPI_DOUBLE_PRECISION, MPI_COMM_WORLD
      use MPIFUN,           only: MPI_Irecv, MPI_Isend, MPI_Comm_dup, MPI_Comm_free
      use named_array_list, only: qna, wna
      use ppp,              only: ppp_main
      use ppp_mpi,          only: piernik_Waitall
#ifdef MPIF08
      use MPIF,             only: MPI_Comm
#endif /* MPIF08 */

      implicit none

      class(cg_level_connected_t), target, intent(inout) :: this     !< object invoking type-bound procedure
      integer(kind=4),                     intent(in)    :: iv       !< variable to be prolonged
      integer(kind=4), optional,           intent(in)    :: pos      !< position of the variable within cell
      integer(kind=4), optional,           intent(in)    :: bnd_type !< Override default boundary type on external boundaries (useful in multigrid solver).
      logical, optional,                   intent(in)    :: dim4     !< operate on wna instead

      type(cg_level_connected_t), pointer                :: fine
      integer                                            :: g
      integer(kind=8), dimension(xdim:zdim, LO:HI)       :: cse, fse         !< shortcuts for coarse and fine segments
      integer(kind=4)                                    :: nr, iw
      type(cg_list_element),            pointer          :: cgl
      type(grid_container),             pointer          :: cg               !< current grid container
      real, dimension(:,:,:),           pointer          :: p3d
      real, dimension(:,:,:,:),         pointer          :: p4d
      logical, save                                      :: warned = .false.
      integer                                            :: position
      integer(kind=8), dimension(ndims, LO:HI)           :: box_8            !< temporary storage
      character(len=*), parameter                        :: pq1_label = "prolong_1v"
      logical                                            :: d4
#ifdef MPIF08
      type(MPI_Comm)  :: p1v_comm
#else /* !MPIF08 */
      integer(kind=4) :: p1v_comm
#endif /* !MPIF08 */

      d4 = .false.
      if (present(dim4)) d4 = dim4

      call ppp_main%start(pq1_label, PPP_AMR)

      if (d4) then
         position = wna%lst(iv)%position(I_ONE)
      else
         position = qna%lst(iv)%position(I_ONE)
      endif
      if (present(pos)) position = pos
      if (position /= VAR_CENTER .and. .not. warned) then
         if (master) call warn("[cg_level_connected:prolong_1var] Only cell-centered interpolation scheme is implemented. Expect inaccurate results for variables that are placed on faces or corners")
         warned = .true.
      endif

      fine => this%finer
      if (.not. associated(fine)) then ! can't prolong finest level
         write(msg,'(a,i3)')"[cg_level_connected:prolong_1var] no finer level than: ", this%l%id
         call warn(msg)
         return
      endif

      call this%vertical_prep
      call fine%vertical_prep

!      call fine%set_dirty(iv, (0.895+0.0001*fine%l%id)*dirtyH1) !> \todo filter this through cg%ignore_prolongation

      if (this%ord_prolong_set /= I_ZERO) then
         !> \todo some variables may need special care on external boundaries
         if (d4) then
            call this%arr4d_boundaries(iv) !, bnd_type = bnd_type)
         else
            call this%arr3d_boundaries(iv, bnd_type = bnd_type)
         endif
      endif
      if (d4) then
         do iw = 1, wna%lst(iv)%dim4
            call this%check_dirty(iv, "prolong-", subfield=iw)
         enddo
      else
         call this%check_dirty(iv, "prolong-")
      endif

      call MPI_Comm_dup(MPI_COMM_WORLD, p1v_comm, err_mpi)
      nr = 0
      ! be ready to receive everything into right buffers
      cgl => fine%first
      do while (associated(cgl))
         cg => cgl%cg
         call cg%costs%start

         associate( seg => cg%pi_tgt%seg )
         if (allocated(cg%pi_tgt%seg)) then
            do g = lbound(seg(:), dim=1), ubound(seg(:), dim=1)
               nr = nr + I_ONE
               if (nr > size(req, dim=1)) call inflate_req
               if (d4) then
                  allocate(seg(g)%buf4(wna%lst(iv)%dim4, size(seg(g)%buf, dim=1), size(seg(g)%buf, dim=2), size(seg(g)%buf, dim=3)))
                  call MPI_Irecv(seg(g)%buf4(1, 1, 1, 1), size(seg(g)%buf4(:, :, :, :), kind=4), MPI_DOUBLE_PRECISION, seg(g)%proc, seg(g)%tag, p1v_comm, req(nr), err_mpi)
               else
                  call MPI_Irecv(seg(g)%buf(1, 1, 1), size(seg(g)%buf(:, :, :), kind=4), MPI_DOUBLE_PRECISION, seg(g)%proc, seg(g)%tag, p1v_comm, req(nr), err_mpi)
               endif
            enddo
         endif
         end associate

         call cg%costs%stop(I_REFINE)
         cgl => cgl%nxt
      enddo

      ! send coarse data
      cgl => this%first
      do while (associated(cgl))
         cg => cgl%cg
         call cg%costs%start

         associate( seg => cg%po_tgt%seg )
         do g = lbound(seg(:), dim=1), ubound(seg(:), dim=1)

            cse = seg(g)%se
            nr = nr + I_ONE
            if (nr > size(req, dim=1)) call inflate_req
            if (d4) then
               allocate(seg(g)%buf4(wna%lst(iv)%dim4, size(seg(g)%buf, dim=1), size(seg(g)%buf, dim=2), size(seg(g)%buf, dim=3)))
               p4d => cg%w(iv)%span(cse)
               seg(g)%buf4(:, :, :, :) = p4d
               call MPI_Isend(seg(g)%buf4(1, 1, 1, 1), size(seg(g)%buf4(:, :, :, :), kind=4), MPI_DOUBLE_PRECISION, seg(g)%proc, seg(g)%tag, p1v_comm, req(nr), err_mpi)
            else
               p3d => cg%q(iv)%span(cse)
               seg(g)%buf(:, :, :) = p3d
               call MPI_Isend(seg(g)%buf(1, 1, 1), size(seg(g)%buf(:, :, :), kind=4), MPI_DOUBLE_PRECISION, seg(g)%proc, seg(g)%tag, p1v_comm, req(nr), err_mpi)
            endif
         enddo
         end associate

         call cg%costs%stop(I_REFINE)
         cgl => cgl%nxt
      enddo

      call piernik_Waitall(nr, "prolong_1v")
      call MPI_Comm_free(p1v_comm, err_mpi)

      ! merge received coarse data into one array and interpolate it into the right place
      cgl => fine%first
      do while (associated(cgl))
         cg => cgl%cg
         call cg%costs%start

         if (allocated(cg%pi_tgt%seg) .and. .not. cg%ignore_prolongation) then

            box_8 = int(cg%ijkse, kind=8)
            cse = f2c(box_8)
            fse = c2f(cse)  ! what about odd-sized or odd-offset cg?

            if (d4) then
               qna%lst(qna%wai)%ord_prolong = 0  !> QUIRKY \todo implement high order conservative prolongation and use wna%lst(i)%ord_prolong here
               do iw = 1, wna%lst(iv)%dim4
                  do g = lbound(cg%pi_tgt%seg(:), dim=1), ubound(cg%pi_tgt%seg(:), dim=1)

                     associate (csep => cg%pi_tgt%seg(g)%se)
                        cg%prolong_(csep(xdim, LO):csep(xdim, HI), csep(ydim, LO):csep(ydim, HI), csep(zdim, LO):csep(zdim, HI)) = cg%pi_tgt%seg(g)%buf4(iw, :, :, :)
                     end associate
                  enddo
                  call cg%prolong(qna%wai, cse, p_xyz = .true.) ! prolong to auxiliary array cg%prolong_xyz
                  cg%w(iv)%arr(iw,       fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), fse(zdim, LO):fse(zdim, HI)) = &
                       &  cg%prolong_xyz(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), fse(zdim, LO):fse(zdim, HI))

               enddo
            else
               do g = lbound(cg%pi_tgt%seg(:), dim=1), ubound(cg%pi_tgt%seg(:), dim=1)

                  associate (csep => cg%pi_tgt%seg(g)%se)
                     cg%prolong_(csep(xdim, LO):csep(xdim, HI), csep(ydim, LO):csep(ydim, HI), csep(zdim, LO):csep(zdim, HI)) = cg%pi_tgt%seg(g)%buf(:,:,:)
                  end associate

                  !> When this%ord_prolong_set /= I_ZERO, the received cg%pi_tgt%seg(:)%buf(:,:,:) may overlap
                  !! The incoming data thus must either contain valid guardcells (even if qna%lst(iv)%ord_prolong == O_INJ)
                  !! or the guardcells must be zeroed before sending data and received buffer should be added to cg%prolong_(:,:,:) not just assigned

               enddo

               call cg%prolong(iv, cse, p_xyz = .false.) ! prolong directly to cg%q(iv)%arr
            endif

         endif

         call cg%costs%stop(I_REFINE)
         cgl => cgl%nxt
      enddo
      call ppp_main%stop(pq1_label, PPP_AMR)

      if (d4) then
         do iw = 1, wna%lst(iv)%dim4
            call fine%check_dirty(iv, "prolong_w+", subfield=iw)
         enddo
      else
         call fine%check_dirty(iv, "prolong_q+")
      endif

      if (d4) then
         cgl => fine%first
         do while (associated(cgl))
            associate( seg => cgl%cg%pi_tgt%seg )
               if (allocated(cgl%cg%pi_tgt%seg)) then
                  do g = lbound(seg(:), dim=1), ubound(seg(:), dim=1)
                     deallocate(seg(g)%buf4)
                  enddo
               endif
            end associate
            cgl => cgl%nxt
         enddo

         cgl => this%first
         do while (associated(cgl))
            associate( seg => cgl%cg%po_tgt%seg )
               do g = lbound(seg(:), dim=1), ubound(seg(:), dim=1)
                  deallocate(seg(g)%buf4)
               enddo
            end associate
            cgl => cgl%nxt
         enddo

      endif

   end subroutine prolong_1var

!>
!! \brief Interpolate boundaries from coarse level at fine-coarse interfaces
!!
!! \details There are two possible approaches to the problem of prolongation of the data from coarse level to the fine guardcells on the fine-coarse interface
!! * Interpolate the coarse data only
!! * interpolate the coarse and fine data
!! When the order of interpolation is 0 (injection) both methods degenerate into one.
!! Both methods have their area of applicability and both should be implemented.
!! \warning This routine does only the first approach.
!!
!! \todo implement local copies without MPI
!<

   subroutine prolong_bnd_from_coarser(this, ind, arr4d, dir, nocorners)

      use cg_cost_data,     only: I_REFINE
      use cg_list,          only: cg_list_element
      use cg_list_global,   only: all_cg
      use constants,        only: I_ONE, xdim, ydim, zdim, LO, HI, base_level_id, PPP_AMR  !, dirtyH1
      use dataio_pub,       only: warn, die
      use domain,           only: dom
      use grid_cont,        only: grid_container
      use grid_helpers,     only: c2f
      use MPIF,             only: MPI_DOUBLE_PRECISION, MPI_COMM_WORLD
      use MPIFUN,           only: MPI_Irecv, MPI_Isend, MPI_Comm_dup, MPI_Comm_free
      use mpisetup,         only: err_mpi, req, inflate_req, master
      use named_array_list, only: qna, wna
      use ppp,              only: ppp_main
      use ppp_mpi,          only: piernik_Waitall
#ifdef MPIF08
      use MPIF,             only: MPI_Comm
#endif /* MPIF08 */

      implicit none

      class(cg_level_connected_t), intent(inout) :: this      !< the list on which to perform the boundary exchange
      integer(kind=4),             intent(in)    :: ind       !< index of the prolonged variable
      logical,         optional,   intent(in)    :: arr4d     !< present only when ind refers to rank-4 arrays
      integer(kind=4), optional,   intent(in)    :: dir       !< select only this direction
      logical,         optional,   intent(in)    :: nocorners !< when .true. then don't care about proper edge and corner update

      type(cg_level_connected_t), pointer :: coarse
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg            !< current grid container
      integer(kind=8), dimension(xdim:zdim, LO:HI) :: cse, fse ! shortcuts for fine segment and coarse segment
      integer(kind=8), dimension(xdim:zdim) :: per, ext_buf
      integer(kind=4) :: nr, iw
      integer :: g
      logical, save :: firstcall = .true.
      character(len=*), parameter :: pbc_label = "prolong_bnd_from_coarser" , pbcv_label = "prolong_bnd_from_coarser:vbp"
#ifdef MPIF08
      type(MPI_Comm)  :: pbfc_comm
#else /* !MPIF08 */
      integer(kind=4) :: pbfc_comm
#endif /* !MPIF08 */

      if (present(dir)) then
         if (firstcall .and. master) call warn("[cg_level_connected:prolong_bnd_from_coarser] dir present but not implemented yet")
      endif
      if (present(nocorners)) then
         if (firstcall .and. master) call warn("[cg_level_connected:prolong_bnd_from_coarser] nocorners present but not implemented yet")
      endif
      if (present(arr4d)) then
         if (.not. arr4d) call die("[cg_level_connected:prolong_bnd_from_coarser] nonsense, fix it in the caller")
      endif

      firstcall = .false.

      coarse => this%coarser

      if (.not. associated(coarse)) return
      if (this%l%id == base_level_id) return ! There are no fine/coarse boundaries on the base level by definition

      call ppp_main%start(pbc_label, PPP_AMR)
      call ppp_main%start(pbcv_label, PPP_AMR)
      call this%vertical_b_prep
      call coarse%vertical_b_prep
      call ppp_main%stop(pbcv_label, PPP_AMR)

      !call this%clear_boundaries(ind, (0.885+0.0001*this%l%id)*dirtyH1) ! not implemented yet
      ext_buf = dom%D_ * all_cg%ord_prolong_nb ! extension of the buffers due to stencil range
      ! OPT: actual stencil range should be used instead

      call MPI_Comm_dup(MPI_COMM_WORLD, pbfc_comm, err_mpi)
      nr = 0
      ! be ready to receive everything into right buffers
      cgl => this%first
      do while (associated(cgl))
         call cgl%cg%costs%start

         associate ( seg => cgl%cg%pib_tgt%seg )
         if (allocated(cgl%cg%pib_tgt%seg)) then
            do g = lbound(seg(:), dim=1), ubound(seg(:), dim=1)
               nr = nr + I_ONE
               if (nr > size(req, dim=1)) call inflate_req
               if (present(arr4d)) then
                  if (allocated(seg(g)%buf4)) call die("[cg_level_connected:prolong_bnd_from_coarser] allocated pib buf4")
                  allocate(seg(g)%buf4(wna%lst(ind)%dim4, size(seg(g)%buf, dim=1), size(seg(g)%buf, dim=2), size(seg(g)%buf, dim=3)))
                  call MPI_Irecv(seg(g)%buf4(1, 1, 1, 1), size(seg(g)%buf4(:, :, :, :), kind=4), MPI_DOUBLE_PRECISION, seg(g)%proc, seg(g)%tag, pbfc_comm, req(nr), err_mpi)
               else
                  call MPI_Irecv(seg(g)%buf(1, 1, 1), size(seg(g)%buf(:, :, :), kind=4), MPI_DOUBLE_PRECISION, seg(g)%proc, seg(g)%tag, pbfc_comm, req(nr), err_mpi)
               endif
            enddo
         endif
         end associate

         call cgl%cg%costs%stop(I_REFINE)
         cgl => cgl%nxt
      enddo

      ! send coarse data
      cgl => coarse%first
      do while (associated(cgl))
         call cgl%cg%costs%start

         associate( seg => cgl%cg%pob_tgt%seg )
         if (allocated(cgl%cg%pob_tgt%seg)) then
            do g = lbound(seg(:), dim=1), ubound(seg(:), dim=1)

               cse = seg(g)%se
               cse(:, LO) = cse(:, LO) - ext_buf
               cse(:, HI) = cse(:, HI) + ext_buf

               nr = nr + I_ONE
               if (nr > size(req, dim=1)) call inflate_req
               if (present(arr4d)) then
                  if (allocated(seg(g)%buf4)) call die("[cg_level_connected:prolong_bnd_from_coarser] allocated pob buf4")
                  allocate(seg(g)%buf4(wna%lst(ind)%dim4, size(seg(g)%buf, dim=1), size(seg(g)%buf, dim=2), size(seg(g)%buf, dim=3)))
                  seg(g)%buf4(:, :, :, :) = cgl%cg%w(ind)%arr(:, cse(xdim, LO):cse(xdim, HI), cse(ydim, LO):cse(ydim, HI), cse(zdim, LO):cse(zdim, HI))
                  call MPI_Isend(seg(g)%buf4(1, 1, 1, 1), size(seg(g)%buf4(:, :, :, :), kind=4), MPI_DOUBLE_PRECISION, seg(g)%proc, seg(g)%tag, pbfc_comm, req(nr), err_mpi)
               else
                  seg(g)%buf(:, :, :)     = cgl%cg%q(ind)%arr(   cse(xdim, LO):cse(xdim, HI), cse(ydim, LO):cse(ydim, HI), cse(zdim, LO):cse(zdim, HI))
                  call MPI_Isend(seg(g)%buf(1, 1, 1), size(seg(g)%buf(:, :, :), kind=4), MPI_DOUBLE_PRECISION, seg(g)%proc, seg(g)%tag, pbfc_comm, req(nr), err_mpi)
               endif
            enddo
         endif
         end associate

         call cgl%cg%costs%stop(I_REFINE)
         cgl => cgl%nxt
      enddo

      call piernik_Waitall(nr, "prolong_bnd_from_coarser", PPP_AMR)
      call MPI_Comm_free(pbfc_comm, err_mpi)

      ! merge received coarse data into one array and interpolate it into the right place
      per(:) = 0
      where (dom%periodic(:)) per(:) = coarse%l%n_d(:)

      cgl => this%first
      do while (associated(cgl))
         cg => cgl%cg
         call cg%costs%start

         associate ( seg => cg%pib_tgt%seg )
         if (allocated(cg%pib_tgt%seg)) then

            !> \todo restore dirty checks when the implementation will be complete
!!$            cg%prolong_ = 3. !dirtyH
!!$            cg%prolong_x = 7.
!!$            cg%prolong_xy = 15.
!!$            cg%prolong_xyz = 31.

            if (size(seg) > 0) then

               do g = lbound(seg(:), dim=1), ubound(seg(:), dim=1)

                  cse = seg(g)%se
                  cse(:, LO) = cse(:, LO) - ext_buf
                  cse(:, HI) = cse(:, HI) + ext_buf

                  fse = c2f(seg(g)%se)
                  fse(:, LO) = max(fse(:, LO), int(cg%lhn(:, LO), kind=8))
                  fse(:, HI) = min(fse(:, HI), int(cg%lhn(:, HI), kind=8))
                  !> When this%ord_prolong_set /= I_ZERO, the incoming data thus must contain valid guardcells

                  if (present(arr4d)) then
                     qna%lst(qna%wai)%ord_prolong = wna%lst(ind)%ord_prolong  ! QUIRKY
                     do iw = 1, wna%lst(ind)%dim4
                        cg%prolong_(cse(xdim, LO):cse(xdim, HI), cse(ydim, LO):cse(ydim, HI), cse(zdim, LO):cse(zdim, HI)) = seg(g)%buf4(iw, :, :, :)
                        call cg%prolong(qna%wai, seg(g)%se, p_xyz=.true.)  ! prolong rank-4 to auxiliary array cg%prolong_xyz.
                        ! qna%wai is required only for indirect determination of prolongation order (TOO QUIRKY)
                        cg%w(ind)%arr(iw,               fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), fse(zdim, LO):fse(zdim, HI)) = &
                             &           cg%prolong_xyz(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), fse(zdim, LO):fse(zdim, HI))
                     enddo
                  else
                     cg%prolong_(cse(xdim, LO):cse(xdim, HI), cse(ydim, LO):cse(ydim, HI), cse(zdim, LO):cse(zdim, HI)) = seg(g)%buf(:,:,:)
                     call cg%prolong(ind, seg(g)%se, p_xyz=.false.)
                  endif

                  ! The cg%prolong above consumes about half of the prolong_bnd_from_coarser execution time
                  ! ToDo: implement special cases, especially O_INJ and O_LIN, perhaps O_I2 too,
                  ! at least for for 3D (precompute full stencils), watch for even/odd issues
                  ! For higher orders should issue a one-time warning

                  ! OPT: in many cases (like in (M)HD directionally split solver we don't need all f/c guardcells.
                  ! Implement efficient directional selection and avoid updating edge/corner f/c when not strictly necessary.
                  ! Full corner update is important in multigrid and may be important for divB cleaning.
                  ! Full corner update is required for prolongation orders better than injection

               enddo

            endif
         endif
         end associate

         call cg%costs%stop(I_REFINE)
         cgl => cgl%nxt
      enddo

      if (present(arr4d)) then

         cgl => this%first
         do while (associated(cgl))
            if (allocated(cgl%cg%pib_tgt%seg)) then
               do g = lbound(cgl%cg%pib_tgt%seg(:), dim=1), ubound(cgl%cg%pib_tgt%seg(:), dim=1)
                  if (allocated(cgl%cg%pib_tgt%seg(g)%buf4)) deallocate(cgl%cg%pib_tgt%seg(g)%buf4)
               enddo
            endif
            cgl => cgl%nxt
         enddo

         ! send coarse data
         cgl => coarse%first
         do while (associated(cgl))
            if (allocated(cgl%cg%pob_tgt%seg)) then
               do g = lbound(cgl%cg%pob_tgt%seg(:), dim=1), ubound(cgl%cg%pob_tgt%seg(:), dim=1)
                  if (allocated(cgl%cg%pob_tgt%seg(g)%buf4)) deallocate(cgl%cg%pob_tgt%seg(g)%buf4)
               enddo
            endif
            cgl => cgl%nxt
         enddo

      endif

      call ppp_main%stop(pbc_label, PPP_AMR)

   end subroutine prolong_bnd_from_coarser

!> \brief interpolate the grid data which has the flag vital set from this%coarser level

   subroutine restrict(this)

      use constants,        only: base_level_id, PPP_AMR
      use dataio_pub,       only: warn
      use named_array_list, only: qna, wna
      use ppp,              only: ppp_main

      implicit none

      class(cg_level_connected_t), target, intent(inout) :: this !< object invoking type-bound procedure

      integer(kind=4) :: i
      character(len=*), parameter :: resq_label = "restrict_qna", resw_label = "restrict_wna"

      call ppp_main%start(resq_label, PPP_AMR)
      do i = lbound(qna%lst(:), dim=1, kind=4), ubound(qna%lst(:), dim=1, kind=4)
         if (qna%lst(i)%vital .and. (qna%lst(i)%multigrid .or. this%l%id > base_level_id)) call this%restrict_1var(i)
      enddo
      call ppp_main%stop(resq_label, PPP_AMR)

      call ppp_main%start(resw_label, PPP_AMR)
      do i = lbound(wna%lst(:), dim=1, kind=4), ubound(wna%lst(:), dim=1, kind=4)
         if (wna%lst(i)%vital .and. (wna%lst(i)%multigrid .or. this%l%id > base_level_id)) then
            if (wna%lst(i)%multigrid) call warn("[cg_level_connected:restrict] mg set for cg%w ???")
            call this%restrict_w_1var(i)
         endif
      enddo

      call ppp_main%stop(resw_label, PPP_AMR)

   end subroutine restrict

!> \brief Restrict all variables to the base level

   recursive subroutine restrict_to_base(this)

      use constants, only: base_level_id

      implicit none

      class(cg_level_connected_t), intent(inout) :: this !< object invoking type-bound procedure

      if (this%l%id <= base_level_id) return
      call this%restrict
      call this%coarser%restrict_to_base

   end subroutine restrict_to_base

!> \brief Restrict as much as possible

   recursive subroutine restrict_to_floor_q_1var(this, iv)

      use dataio_pub,       only: warn
      use named_array_list, only: qna

      implicit none

      class(cg_level_connected_t), intent(inout) :: this !< object invoking type-bound procedure
      integer(kind=4),             intent(in)    :: iv   !< variable to be restricted

      if (.not. associated(this%coarser)) return
      if (.not. qna%lst(iv)%multigrid) then
         call warn("[cg_level_connected:restrict_to_floor_q_1var] '" // trim(qna%lst(iv)%name) // "' is a non-multigrid field, redirecting to restrict_to_base_q_1var")
         call this%restrict_to_base_q_1var(iv)
      else
         call this%restrict_1var(iv)
         call this%coarser%restrict_to_floor_q_1var(iv)
      endif

   end subroutine restrict_to_floor_q_1var

!> \brief Restrict to the base level

   recursive subroutine restrict_to_base_q_1var(this, iv)

      use constants, only: base_level_id

      implicit none

      class(cg_level_connected_t), intent(inout) :: this !< object invoking type-bound procedure
      integer(kind=4),             intent(in)    :: iv   !< variable to be restricted

      if (this%l%id <= base_level_id) return
      call this%restrict_1var(iv)
      call this%coarser%restrict_to_base_q_1var(iv)

   end subroutine restrict_to_base_q_1var

   recursive subroutine restrict_to_base_w_1var(this, iv)

      use constants, only: base_level_id

      implicit none

      class(cg_level_connected_t), intent(inout) :: this !< object invoking type-bound procedure
      integer(kind=4),             intent(in)    :: iv   !< variable to be restricted

      if (this%l%id <= base_level_id) return
      call this%restrict_w_1var(iv)
      call this%coarser%restrict_to_base_w_1var(iv)

   end subroutine restrict_to_base_w_1var

!>
!! \brief Simplest restriction (averaging).
!!
!! \todo implement high order restriction and test its influence on V-cycle convergence rate. Watch f/c boundaries.
!!
!! \details Some data can be locally copied without MPI, but this seems to have really little impact on the performance.
!! Some tests show that purely MPI code without local copies is marginally faster.
!!
!! OPT Usually there are many messages that are sent between the same pairs of processes
!! \todo Sort all messages according to e.g. tag and send/receive aggregated message with everything
!!
!! \todo implement local copies without MPI anyway
!<

   subroutine restrict_1var(this, iv, pos, dim4)

      use cg_cost_data,     only: I_REFINE
      use constants,        only: xdim, ydim, zdim, ndims, LO, HI, I_ONE, refinement_factor, VAR_CENTER, GEO_XYZ, GEO_RPZ
      use dataio_pub,       only: msg, warn, die
      use domain,           only: dom
      use cg_list,          only: cg_list_element
      use grid_cont,        only: grid_container
      use mpisetup,         only: err_mpi, req, inflate_req, master
      use MPIF,             only: MPI_DOUBLE_PRECISION, MPI_COMM_WORLD
      use MPIFUN,           only: MPI_Irecv, MPI_Isend, MPI_Comm_dup, MPI_Comm_free
      use named_array,      only: p3, p4
      use named_array_list, only: qna, wna
      use ppp_mpi,          only: piernik_Waitall
#ifdef MPIF08
      use MPIF,             only: MPI_Comm
#endif /* MPIF08 */

      implicit none

      class(cg_level_connected_t), target, intent(inout) :: this  !< object invoking type-bound procedure
      integer(kind=4),                     intent(in)    :: iv    !< variable to be restricted
      integer(kind=4), optional,           intent(in)    :: pos   !< position of the variable within cell
      logical, optional,                   intent(in)    :: dim4  !< operate on wna instead

      type(cg_level_connected_t), pointer                :: coarse
      integer                                            :: g
      integer(kind=8), dimension(xdim:zdim, LO:HI)       :: fse, cse              !< shortcuts for fine segment and coarse segment
      integer(kind=8)                                    :: i, j, k, ic, jc, kc
      integer(kind=8), dimension(xdim:zdim)              :: off1
      real                                               :: norm
      integer(kind=4)                                    :: nr
      type(cg_list_element), pointer                     :: cgl
      type(grid_container),  pointer                     :: cg                    !< current grid container
      logical, save                                      :: warned = .false.
      integer                                            :: position
      logical                                            :: d4
#ifdef MPIF08
      type(MPI_Comm)  :: r1v_comm
#else /* !MPIF08 */
      integer(kind=4) :: r1v_comm
#endif /* !MPIF08 */

      d4 = .false.
      if (present(dim4)) d4 = dim4

      position = qna%lst(iv)%position(I_ONE)
      if (present(pos)) position = pos
      if (position /= VAR_CENTER .and. .not. warned) then
         if (master) call warn("[cg_level_connected:restrict_1var] Only cell-centered interpolation scheme is implemented. Expect inaccurate results for variables that are placed on faces or corners")
         warned = .true.
      endif

      ! ToDo warn about positions when d4

      coarse => this%coarser
      if (.not. associated(coarse)) then ! can't restrict base level
         write(msg,'(a,i3)')"[cg_level_connected:restrict_1var] no coarser level than ", this%l%id
         call warn(msg)
         return
      endif

      call this%vertical_prep
      call coarse%vertical_prep

      ! be ready to receive everything into right buffers
      call MPI_Comm_dup(MPI_COMM_WORLD, r1v_comm, err_mpi)
      nr = 0
      cgl => coarse%first
      do while (associated(cgl))
         cg => cgl%cg
         call cg%costs%start

         if (allocated(cg%ri_tgt%seg)) then
            do g = lbound(cg%ri_tgt%seg(:), dim=1), ubound(cg%ri_tgt%seg(:), dim=1)
               nr = nr + I_ONE
               if (nr > size(req, dim=1)) call inflate_req
               associate (seg => cg%ri_tgt%seg(g))
                  if (d4) then
                     allocate(seg%buf4(wna%lst(iv)%dim4, size(seg%buf, dim=1), size(seg%buf, dim=2), size(seg%buf, dim=3)))
                     call MPI_Irecv(seg%buf4(1, 1, 1, 1), size(seg%buf4(:, :, :, :), kind=4), MPI_DOUBLE_PRECISION, seg%proc, seg%tag, r1v_comm, req(nr), err_mpi)
                  else
                     call MPI_Irecv(seg%buf(1, 1, 1), size(seg%buf(:, :, :), kind=4), MPI_DOUBLE_PRECISION, seg%proc, seg%tag, r1v_comm, req(nr), err_mpi)
                  endif
               end associate
            enddo
         endif

         call cg%costs%stop(I_REFINE)
         cgl => cgl%nxt
      enddo

      ! interpolate to coarse buffer and send it
      norm = 1./refinement_factor**dom%eff_dim
      cgl => this%first
      do while (associated(cgl))
         cg => cgl%cg
         call cg%costs%start

         do g = lbound(cg%ro_tgt%seg(:), dim=1), ubound(cg%ro_tgt%seg(:), dim=1)

            associate (seg => cg%ro_tgt%seg(g))
               if (d4) then
                  allocate(seg%buf4(wna%lst(iv)%dim4, size(seg%buf, dim=1), size(seg%buf, dim=2), size(seg%buf, dim=3)))
                  seg%buf4(:, :, :, :) = 0.
               else
                  seg%buf(:, :, :) = 0.
               endif

               fse(:,:) = seg%se(:,:)
               off1(:) = mod(seg%se(:, LO), int(refinement_factor, kind=8))
               if (all(off1 == 0) .and. all(mod(fse(:, HI)-fse(:, LO), int(refinement_factor, kind=8)) == 1) .and. dom%eff_dim == ndims) then
                  ! This is the easy, even offset/even size case. Happens in AMR and when UG has regular cartesian decomposition.
                  ! It is few times faster than the code for odd cases below
                  select case (dom%geometry_type)
                     case (GEO_XYZ)
                        if (d4) then
                           seg%buf4(:, 1:1+(fse(xdim, HI)-fse(xdim, LO))/refinement_factor, &
                                &      1:1+(fse(ydim, HI)-fse(ydim, LO))/refinement_factor, &
                                &      1:1+(fse(zdim, HI)-fse(zdim, LO))/refinement_factor) = ( &
                                cg%w(iv)%arr(:, fse(xdim, LO):fse(xdim, HI)-1:2, fse(ydim, LO):fse(ydim, HI)-1:2, fse(zdim, LO):fse(zdim, HI)-1:2) + &
                                cg%w(iv)%arr(:, fse(xdim, LO)+1:fse(xdim, HI):2, fse(ydim, LO):fse(ydim, HI)-1:2, fse(zdim, LO):fse(zdim, HI)-1:2) + &
                                cg%w(iv)%arr(:, fse(xdim, LO):fse(xdim, HI)-1:2, fse(ydim, LO)+1:fse(ydim, HI):2, fse(zdim, LO):fse(zdim, HI)-1:2) + &
                                cg%w(iv)%arr(:, fse(xdim, LO)+1:fse(xdim, HI):2, fse(ydim, LO)+1:fse(ydim, HI):2, fse(zdim, LO):fse(zdim, HI)-1:2) + &
                                cg%w(iv)%arr(:, fse(xdim, LO):fse(xdim, HI)-1:2, fse(ydim, LO):fse(ydim, HI)-1:2, fse(zdim, LO)+1:fse(zdim, HI):2) + &
                                cg%w(iv)%arr(:, fse(xdim, LO)+1:fse(xdim, HI):2, fse(ydim, LO):fse(ydim, HI)-1:2, fse(zdim, LO)+1:fse(zdim, HI):2) + &
                                cg%w(iv)%arr(:, fse(xdim, LO):fse(xdim, HI)-1:2, fse(ydim, LO)+1:fse(ydim, HI):2, fse(zdim, LO)+1:fse(zdim, HI):2) + &
                                cg%w(iv)%arr(:, fse(xdim, LO)+1:fse(xdim, HI):2, fse(ydim, LO)+1:fse(ydim, HI):2, fse(zdim, LO)+1:fse(zdim, HI):2)) * norm
                        else
                           seg%buf(1:1+(fse(xdim, HI)-fse(xdim, LO))/refinement_factor, &
                                &  1:1+(fse(ydim, HI)-fse(ydim, LO))/refinement_factor, &
                                &  1:1+(fse(zdim, HI)-fse(zdim, LO))/refinement_factor) = ( &
                                cg%q(iv)%arr(fse(xdim, LO):fse(xdim, HI)-1:2, fse(ydim, LO):fse(ydim, HI)-1:2, fse(zdim, LO):fse(zdim, HI)-1:2) + &
                                cg%q(iv)%arr(fse(xdim, LO)+1:fse(xdim, HI):2, fse(ydim, LO):fse(ydim, HI)-1:2, fse(zdim, LO):fse(zdim, HI)-1:2) + &
                                cg%q(iv)%arr(fse(xdim, LO):fse(xdim, HI)-1:2, fse(ydim, LO)+1:fse(ydim, HI):2, fse(zdim, LO):fse(zdim, HI)-1:2) + &
                                cg%q(iv)%arr(fse(xdim, LO)+1:fse(xdim, HI):2, fse(ydim, LO)+1:fse(ydim, HI):2, fse(zdim, LO):fse(zdim, HI)-1:2) + &
                                cg%q(iv)%arr(fse(xdim, LO):fse(xdim, HI)-1:2, fse(ydim, LO):fse(ydim, HI)-1:2, fse(zdim, LO)+1:fse(zdim, HI):2) + &
                                cg%q(iv)%arr(fse(xdim, LO)+1:fse(xdim, HI):2, fse(ydim, LO):fse(ydim, HI)-1:2, fse(zdim, LO)+1:fse(zdim, HI):2) + &
                                cg%q(iv)%arr(fse(xdim, LO):fse(xdim, HI)-1:2, fse(ydim, LO)+1:fse(ydim, HI):2, fse(zdim, LO)+1:fse(zdim, HI):2) + &
                                cg%q(iv)%arr(fse(xdim, LO)+1:fse(xdim, HI):2, fse(ydim, LO)+1:fse(ydim, HI):2, fse(zdim, LO)+1:fse(zdim, HI):2)) * norm
                        endif
                     case (GEO_RPZ)
                        if (d4) call die("clc:r1v GEO_RPZ not implemented (restrict)")
                        do i = fse(xdim, LO), fse(xdim, HI)
                           seg%buf     (  1+(i            -fse(xdim, LO))/refinement_factor, &
                                &                    1:1+(fse(ydim, HI)-fse(ydim, LO))/refinement_factor, &
                                &                    1:1+(fse(zdim, HI)-fse(zdim, LO))/refinement_factor) = &
                                seg%buf(  1+(i            -fse(xdim, LO))/refinement_factor, &
                                &                    1:1+(fse(ydim, HI)-fse(ydim, LO))/refinement_factor, &
                                &                    1:1+(fse(zdim, HI)-fse(zdim, LO))/refinement_factor) + ( &
                                cg%q(iv)%arr(i, fse(ydim, LO):fse(ydim, HI)-1:2, fse(zdim, LO):fse(zdim, HI)-1:2) + &
                                cg%q(iv)%arr(i, fse(ydim, LO)+1:fse(ydim, HI):2, fse(zdim, LO):fse(zdim, HI)-1:2) + &
                                cg%q(iv)%arr(i, fse(ydim, LO):fse(ydim, HI)-1:2, fse(zdim, LO)+1:fse(zdim, HI):2) + &
                                cg%q(iv)%arr(i, fse(ydim, LO)+1:fse(ydim, HI):2, fse(zdim, LO)+1:fse(zdim, HI):2)) * norm * cg%x(i)
                        enddo
                     case default
                        call die("[cg_level_connected:restrict_1var] Unknown geometry (1)")
                  end select
               else
                  ! OPT: Split the problem into the core that can be done by array arithmetic and finish the edges where necessary
                  do k = fse(zdim, LO), fse(zdim, HI)
                     kc = (k-fse(zdim, LO)+off1(zdim))/refinement_factor + 1
                     do j = fse(ydim, LO), fse(ydim, HI)
                        jc = (j-fse(ydim, LO)+off1(ydim))/refinement_factor + 1
                        do i = fse(xdim, LO), fse(xdim, HI)
                           ic = (i-fse(xdim, LO)+off1(xdim))/refinement_factor + 1
                           select case (dom%geometry_type)
                              case (GEO_XYZ)
                                 if (d4) then
                                    seg%buf4(:, ic, jc, kc) = seg%buf4(:, ic, jc, kc) + cg%w(iv)%arr(:, i, j, k) * norm
                                 else
                                    seg%buf(ic, jc, kc) = seg%buf(ic, jc, kc) + cg%q(iv)%arr(i, j, k) * norm
                                 endif
                              case (GEO_RPZ)
                                 if (d4) call die("clc:r1v odd not implemented for GEO_RPZ")
                                 seg%buf(ic, jc, kc) = seg%buf(ic, jc, kc) + cg%q(iv)%arr(i, j, k) * norm * cg%x(i)
                              case default
                                 call die("[cg_level_connected:restrict_1var] Unknown geometry (2)")
                           end select
                        enddo
                     enddo
                  enddo
               endif
               nr = nr + I_ONE
               if (nr > size(req, dim=1)) call inflate_req
               if (d4) then
                  call MPI_Isend(seg%buf4(1, 1, 1, 1), size(seg%buf4(:, :, :, :), kind=4), MPI_DOUBLE_PRECISION, seg%proc, seg%tag, r1v_comm, req(nr), err_mpi)
               else
                  call MPI_Isend(seg%buf(1, 1, 1), size(seg%buf(:, :, :), kind=4), MPI_DOUBLE_PRECISION, seg%proc, seg%tag, r1v_comm, req(nr), err_mpi)
               endif
            end associate
         enddo

         call cg%costs%stop(I_REFINE)
         cgl => cgl%nxt
      enddo

      call piernik_Waitall(nr, "restrict_1v")
      call MPI_Comm_free(r1v_comm, err_mpi)

      ! copy the received buffers to the right places
      cgl => coarse%first
      do while (associated(cgl))
         cg => cgl%cg
         call cg%costs%start

         if (allocated(cg%ri_tgt%seg)) then

            ! disables check_dirty
            if (d4) then
               do g = 1, wna%lst(iv)%dim4
                  where (.not. cg%leafmap(:,:,:)) cg%w(iv)%arr(g, RNG) = 0.
               enddo
            else
               where (.not. cg%leafmap(:,:,:)) cg%q(iv)%arr(RNG) = 0.
            endif

            do g = lbound(cg%ri_tgt%seg(:), dim=1), ubound(cg%ri_tgt%seg(:), dim=1)
               cse(:,:) = cg%ri_tgt%seg(g)%se(:,:)
               select case (dom%geometry_type)
                  case (GEO_XYZ) ! do nothing
                  case (GEO_RPZ)
                     if (d4) call die("clc:r1v GEO_RPZ not implemented (copy)")
                     do i = lbound(cg%ri_tgt%seg(g)%buf, dim=1), ubound(cg%ri_tgt%seg(g)%buf, dim=1)
                        ic = cse(xdim, LO) +i - lbound(cg%ri_tgt%seg(g)%buf, dim=1)
                        cg%ri_tgt%seg(g)%buf(i, :, :) = cg%ri_tgt%seg(g)%buf(i, :, :) / cg%x(ic)
                     enddo
                  case default
                     call die("[cg_level_connected:restrict_1var] Unknown geometry")
               end select
               if (d4) then
                  p4 => cg%w(iv)%span(cse)
                  p4 = p4 + cg%ri_tgt%seg(g)%buf4(:, :, :, :) !errors on overlap?
               else
                  p3 => cg%q(iv)%span(cse)
                  p3 = p3 + cg%ri_tgt%seg(g)%buf(:, :, :) !errors on overlap?
               endif
            enddo
         endif

         call cg%costs%stop(I_REFINE)
         cgl => cgl%nxt
      enddo

      if (d4) then

         cgl => this%first
         do while (associated(cgl))
            cg => cgl%cg
            if (allocated(cg%ro_tgt%seg)) then
               do g = lbound(cg%ro_tgt%seg(:), dim=1), ubound(cg%ro_tgt%seg(:), dim=1)
                  if (allocated(cg%ro_tgt%seg(g)%buf4)) deallocate(cg%ro_tgt%seg(g)%buf4)
               enddo
            endif
            cgl => cgl%nxt
         enddo

         cgl => coarse%first
         do while (associated(cgl))
            cg => cgl%cg
            if (allocated(cg%ri_tgt%seg)) then
               do g = lbound(cg%ri_tgt%seg(:), dim=1), ubound(cg%ri_tgt%seg(:), dim=1)
                  if (allocated(cg%ri_tgt%seg(g)%buf4)) deallocate(cg%ri_tgt%seg(g)%buf4)
               enddo
            endif
            cgl => cgl%nxt
         enddo

      endif

   end subroutine restrict_1var

!> \brief Quick and dirty restriction of 4D arrays. OPTIMIZE ME!

   subroutine restrict_w_1var(this, i)

      use constants,        only: GEO_RPZ
      use dataio_pub,       only: warn
      use domain,           only: dom
      use fluidindex,       only: iarr_all_my
      use named_array_list, only: qna, wna

      implicit none

      class(cg_level_connected_t), target, intent(inout) :: this  !< object invoking type-bound procedure
      integer(kind=4),                     intent(in)    :: i     !< variable to be restricted

      integer(kind=4) :: iw
      logical, save :: warned = .false.

      if (dom%geometry_type == GEO_RPZ .and. i == wna%fi) then ! take the slow way
         if (.not. warned) then
            call warn("[cg_level_connected:restrict_w_1var] Using the slow w-q-w copy algorithm")
            warned = .true.
         endif
         do iw = 1, wna%lst(i)%dim4
            call this%wq_copy(i, iw, qna%wai)
            if (any(iw == iarr_all_my)) call this%mul_by_r(qna%wai) ! angular momentum conservation
            ! this is required because we don't use (.not. cg%leafmap) mask in the this%coarser%qw_copy call below
            call this%coarser%wq_copy(i, iw, qna%wai)

            if (any(iw == iarr_all_my)) call this%coarser%mul_by_r(qna%wai) ! angular momentum conservation
            call this%restrict_1var(qna%wai, wna%lst(i)%position(iw))

            if (any(iw == iarr_all_my)) call this%coarser%div_by_r(qna%wai) ! angular momentum conservation
            call this%coarser%qw_copy(qna%wai, i, iw)
         enddo
      else
         call this%restrict_1var(i, dim4 = .true.)
      endif

   end subroutine restrict_w_1var

!>
!! \brief This routine sets up all guardcells (internal, external and fine-coarse)
!! for given rank-3 arrays on a single level.
!!
!! Use of this routine for whole stack of levels (leaves) will result in double update
!! of internal boundaries on all levels below the finest%level when prolongation
!! order is higher than injection.
!!
!! Currently only multigrid solver use level-wise updates on rank-3 array,
!! and leaves%leaf_arr3d_boundaries can't be used instead there.
!<

   subroutine arr3d_boundaries(this, ind, area_type, bnd_type, dir, nocorners)

      use constants,        only: base_level_id, PPP_AMR, O_INJ
      use global,           only: dirty_debug
      use named_array_list, only: qna
      use ppp,              only: ppp_main

      implicit none

      class(cg_level_connected_t), intent(inout) :: this      !< the list on which to perform the boundary exchange
      integer(kind=4),             intent(in)    :: ind       !< index of cg%q(:) 3d array
      integer(kind=4), optional,   intent(in)    :: area_type !< defines how do we treat boundaries
      integer(kind=4), optional,   intent(in)    :: bnd_type  !< Override default boundary type on external boundaries (useful in multigrid solver).
                                                              !< Note that BND_PER, BND_MPI, BND_SHE and BND_COR aren't external and cannot be overridden
      integer(kind=4), optional,   intent(in)    :: dir       !< select only this direction
      logical,         optional,   intent(in)    :: nocorners !< .when .true. then don't care about proper edge and corner update

      character(len=*), parameter :: a3b_label = "level:arr3d_boundaries", a3bp_label = "level:arr3d_boundaries:prolong"

      call ppp_main%start(a3b_label)

      if (dirty_debug) call this%dirty_boundaries(ind)

      call ppp_main%start(a3bp_label, PPP_AMR)
      if (associated(this%coarser) .and. this%l%id > base_level_id .and. qna%lst(ind)%ord_prolong /= O_INJ) &
           call this%coarser%level_3d_boundaries(ind, bnd_type = bnd_type)
      ! bnd_type = BND_NEGREF above is critical for convergence of multigrid with isolated boundaries.
      call this%prolong_bnd_from_coarser(ind, dir=dir, nocorners=nocorners)
      call ppp_main%stop(a3bp_label, PPP_AMR)
      call this%level_3d_boundaries(ind, area_type=area_type, bnd_type=bnd_type, dir=dir, nocorners=nocorners)
      ! The correctness of the sequence of calls above may depend on the implementation of internal boundary exchange

      call ppp_main%stop(a3b_label)

   end subroutine arr3d_boundaries

!>
!! \brief This routine sets up all guardcells (internal, external and fine-coarse)
!! for given rank-4 arrays on a single level.
!!
!! Use of this routine for whole stack of levels (leaves) will result in double update
!! of internal boundaries on all levels below the finest%level when prolongation
!! order is higher than injection.
!!
!! This routine should be used only for level-wise updates of rank-4 arrays (like in prolongation called from refinement update).
!! In most other cases use leaves%leaf_arr4d_boundaries instead.
!<

   subroutine arr4d_boundaries(this, ind, area_type, dir, nocorners)

      use constants,        only: base_level_id, O_INJ
      use named_array_list, only: wna
      use ppp,              only: ppp_main

      implicit none

      class(cg_level_connected_t), intent(inout) :: this      !< the list on which to perform the boundary exchange
      integer(kind=4),             intent(in)    :: ind       !< index of cg%w(:) 4d array
      integer(kind=4), optional,   intent(in)    :: area_type !< defines how do we treat boundaries
      integer(kind=4), optional,   intent(in)    :: dir       !< select only this direction
      logical,         optional,   intent(in)    :: nocorners !< .when .true. then don't care about proper edge and corner update

      character(len=*), parameter :: a4b_label = "lev:a4d_bnd"

      call ppp_main%start(a4b_label)

!      call this%dirty_boundaries(ind)
!      call this%clear_boundaries(ind, value=10.)

      if (associated(this%coarser) .and. this%l%id > base_level_id) then
         if (wna%lst(ind)%ord_prolong /= O_INJ) call this%coarser%level_4d_boundaries(ind)  ! overkill in most places
         ! here we can use any high order prolongation without destroying conservation
         call this%prolong_bnd_from_coarser(ind, arr4d=.true., dir=dir, nocorners=nocorners)
      endif

      call this%level_4d_boundaries(ind, area_type=area_type, dir=dir, nocorners=nocorners)

      call ppp_main%stop(a4b_label)

   end subroutine arr4d_boundaries

end module cg_level_connected
