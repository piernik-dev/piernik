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

!> \brief Module providing the full, usable grid container type and its methods that don't fit to any abstract subtypes of grid container

module grid_cont

   use constants,         only: LO, HI
   use grid_cont_bnd,     only: segment
   use grid_cont_prolong, only: grid_container_prolong_t
   use refinement_flag,   only: ref_flag

   implicit none

   private
   public :: grid_container

   !< \brief target list container for prolongations, restrictions and boundary exchanges
   type :: tgt_list
      type(segment), dimension(:), allocatable :: seg              !< a segment of data to be received or sent
   end type tgt_list

   !> \brief Everything required for autonomous computation of a single sweep on a portion of the domain on a single process
   type, extends(grid_container_prolong_t) :: grid_container

      ! Prolongation and restriction

      ! these lists are initialized and maintained in cg_level_connected
      ! perhaps parts of vertical_bf_prep can be moved here
      ! vertical_b_prep has to stay in cg_level_connected because we don't have dot here
      type(tgt_list) :: ri_tgt                                    !< description of incoming restriction data (this should be a linked list)
      type(tgt_list) :: ro_tgt                                    !< description of outgoing restriction data
      type(tgt_list) :: pi_tgt                                    !< description of incoming prolongation data
      type(tgt_list) :: po_tgt                                    !< description of outgoing prolongation data
      type(tgt_list) :: pib_tgt                                   !< description of incoming boundary prolongation data
      type(tgt_list) :: pob_tgt                                   !< description of outgoing boundary prolongation data
      type(tgt_list) :: rif_tgt                                   !< description of fluxes incoming from fine grid
      type(tgt_list) :: rof_tgt                                   !< description of fluxes outgoing to coarse grid

      ! Refinements

      logical, allocatable, dimension(:,:,:) :: leafmap           !< .true. when a cell is not covered by finer cells, .false. otherwise
      logical, allocatable, dimension(:,:,:) :: refinemap         !< .true. when a cell triggers refinement criteria, .false. otherwise

      ! Misc
      integer(kind=8) :: SFC_id                                  !< position of the grid on space-filling curve
      type(ref_flag) :: refine_flags                             !< refine or derefine this grid container?
      integer :: membership                                      !< How many cg lists use this grid piece?
      logical :: ignore_prolongation                             !< When .true. do not upgrade interior with incoming prolonged values
      logical :: is_old                                          !< .true. if a given grid existed prior to  upgrade_refinement call
      logical :: processed                                       !< for use in sweeps.F90

   contains

      procedure          :: init_gc                              !< Initialization
      procedure          :: cleanup                              !< Deallocate all internals
      procedure          :: update_leafmap                       !< Check if the grid container has any parts covered by finer grids and update appropriate map
      procedure          :: refinemap2SFC_list                   !< create list of SFC indices to be created from refine flags

   end type grid_container

contains

!> \brief This method sets up remaining grid container variables and arrays.

   subroutine init_gc(this, my_se, grid_id, l)

      use constants,        only: PIERNIK_INIT_DOMAIN, xdim, ydim, zdim, LO, HI
      use dataio_pub,       only: die, code_progress
      use level_essentials, only: level_t
      use ordering,         only: SFC_order

      implicit none

      class(grid_container), target,   intent(inout) :: this  ! intent(out) would silently clear everything, that was already set
                                                              ! (also the fields in types derived from grid_container)
      integer(kind=8), dimension(:,:), intent(in)    :: my_se    !< my segment
      integer,                         intent(in)    :: grid_id  !< ID which should be unique across level
      class(level_t), pointer,         intent(in)    :: l        !< level essential data

      if (code_progress < PIERNIK_INIT_DOMAIN) call die("[grid_container:init_gc] MPI not initialized.")

      call this%init_gc_base(my_se, grid_id, l)
      call this%init_gc_bnd
      call this%add_all_na
      call this%init_gc_prolong

      this%membership = 1
      this%SFC_id     = SFC_order(this%my_se(:, LO) - l%off)

      allocate(this%leafmap  (this%ijkse(xdim, LO):this%ijkse(xdim, HI), this%ijkse(ydim, LO):this%ijkse(ydim, HI), this%ijkse(zdim, LO):this%ijkse(zdim, HI)), &
           &   this%refinemap(this%lhn(xdim, LO):this%lhn(xdim, HI), this%lhn(ydim, LO):this%lhn(ydim, HI), this%lhn(zdim, LO):this%lhn(zdim, HI)))

      this%leafmap    (:, :, :) = .true.
      this%refinemap  (:, :, :) = .false.
      call this%refine_flags%init
      this%ignore_prolongation = .false.
      this%is_old = .false.
      this%has_previous_timestep = .false.

   end subroutine init_gc

!> \brief Routines that deallocates all internals of the grid container

   subroutine cleanup(this)

      implicit none

      class(grid_container), intent(inout) :: this !< object invoking type-bound procedure

      integer :: b, g
      integer, parameter :: nseg = 4*2
      type(tgt_list), dimension(nseg) :: rpio_tgt

      call this%cleanup_base
      call this%cleanup_na
      call this%cleanup_bnd
      call this%cleanup_prolong

      rpio_tgt(1:nseg) = [ this%ri_tgt,  this%ro_tgt,  this%pi_tgt,  this%po_tgt, &
           &               this%pib_tgt, this%pob_tgt, this%rif_tgt, this%rof_tgt ]
      do b = 1, nseg
         if (allocated(rpio_tgt(b)%seg)) then
            do g = lbound(rpio_tgt(b)%seg, dim=1), ubound(rpio_tgt(b)%seg, dim=1)
               if (allocated(rpio_tgt(b)%seg(g)%buf)) deallocate(rpio_tgt(b)%seg(g)%buf)
            enddo
            deallocate(rpio_tgt(b)%seg)
         endif
      enddo

      ! arrays not handled through named_array feature
      if (allocated(this%leafmap))   deallocate(this%leafmap)
      if (allocated(this%refinemap)) deallocate(this%refinemap)

   end subroutine cleanup

!> \brief Check if the grid container has any parts covered by finer grids and update appropriate map

   subroutine update_leafmap(this)

      use constants, only: xdim, ydim, zdim, LO, HI

      implicit none

      class(grid_container), intent(inout) :: this !< object invoking type-bound procedure

      integer(kind=8), dimension(xdim:zdim, LO:HI) :: se
      integer :: g

      this%leafmap = .true.
      if (allocated(this%ri_tgt%seg)) then
         do g = lbound(this%ri_tgt%seg(:), dim=1), ubound(this%ri_tgt%seg(:), dim=1)
            se(:, :) = this%ri_tgt%seg(g)%se(:, :)
            this%leafmap(se(xdim, LO):se(xdim, HI), se(ydim, LO):se(ydim, HI), se(zdim, LO):se(zdim, HI)) = .false.
         enddo
      endif

   end subroutine update_leafmap

!< \brief Create list of SFC indices to be created from refine flags

   subroutine refinemap2SFC_list(this)

      use constants,  only: refinement_factor, xdim, ydim, zdim, I_ONE
      use dataio_pub, only: die, warn
      use domain,     only: dom
      use refinement, only: bsize

      implicit none

      class(grid_container), intent(inout) :: this !< object invoking type-bound procedure

      integer :: i, j, k, ifs, ife, jfs, jfe, kfs, kfe
      enum, bind(C)
         enumerator :: NONE, REFINE, LEAF
      end enum
      integer :: type
      logical, save :: warned = .false.

      this%refinemap(this%is:this%ie, this%js:this%je, this%ks:this%ke) = &
           this%refinemap(this%is:this%ie, this%js:this%je, this%ks:this%ke) .and. this%leafmap
      type = NONE
      if (any(this%refinemap)) then
         type = REFINE
      else if (this%refine_flags%refine) then
         type = LEAF
         if (.not. warned) then
            warned = .true.
            call warn("[grid_container:refinemap2SFC_list] direct use of cg%refine_flags%refine is deprecated")
         endif
      endif

      if (type == NONE) return

      if (any((bsize <= 0) .and. dom%has_dir)) return ! this routine works only with blocky AMR

      !! ToDo: precompute refinement decomposition in this%init_gc and simplify the code below.
      !! It should also simplify decomposition management and make it more flexible in case we decide to work on uneven AMR blocks

      !! beware: consider dropping this%l%off feature for simplicity. It will require handling the shift due to domain expansion (some increase CPU cost)

      associate( b_size => merge(bsize, huge(I_ONE), dom%has_dir))
         do i = int(((this%is - this%l%off(xdim))*refinement_factor) / b_size(xdim)), int(((this%ie - this%l%off(xdim))*refinement_factor + I_ONE) / b_size(xdim))
            ifs = max(int(this%is), int(this%l%off(xdim)) + (i*b_size(xdim))/refinement_factor)
            ife = min(int(this%ie), int(this%l%off(xdim)) + ((i+I_ONE)*b_size(xdim)-I_ONE)/refinement_factor)

            do j = int(((this%js - this%l%off(ydim))*refinement_factor) / b_size(ydim)), int(((this%je - this%l%off(ydim))*refinement_factor + I_ONE) / b_size(ydim))
               jfs = max(int(this%js), int(this%l%off(ydim)) + (j*b_size(ydim))/refinement_factor)
               jfe = min(int(this%je), int(this%l%off(ydim)) + ((j+I_ONE)*b_size(ydim)-I_ONE)/refinement_factor)

               do k = int(((this%ks - this%l%off(zdim))*refinement_factor) / b_size(zdim)), int(((this%ke - this%l%off(zdim))*refinement_factor + I_ONE) / b_size(zdim))
                  kfs = max(int(this%ks), int(this%l%off(zdim)) + (k*b_size(zdim))/refinement_factor)
                  kfe = min(int(this%ke), int(this%l%off(zdim)) + ((k+I_ONE)*b_size(zdim)-I_ONE)/refinement_factor)
                  select case (type)
                     case (REFINE)
                        if (any(this%refinemap(ifs:ife, jfs:jfe, kfs:kfe))) call this%refine_flags%add(this%l%id+I_ONE, int([i, j, k]*b_size, kind=8)+refinement_factor*this%l%off, refinement_factor*this%l%off)
                     case (LEAF)
                        if (all(this%leafmap(ifs:ife, jfs:jfe, kfs:kfe))) then
                           call this%refine_flags%add(this%l%id+I_ONE, int([i, j, k]*b_size, kind=8)+refinement_factor*this%l%off, refinement_factor*this%l%off)
                        else if (any(this%leafmap(ifs:ife, jfs:jfe, kfs:kfe))) then
                           call die("[grid_container:refinemap2SFC_list] cannot refine partially leaf parf of the grid")
                        endif
                     case default
                        call die("[grid_container:refinemap2SFC_list] invalid type")
                  end select
               enddo
            enddo
         enddo
      end associate
      this%refinemap = .false.

   end subroutine refinemap2SFC_list

end module grid_cont
