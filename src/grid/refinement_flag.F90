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

!> \brief This module contains variables and routines related to refinement flags and derefinement requests

module refinement_flag

   use constants, only: ndims

   implicit none

   private
   public :: ref_flag_t

   !> A candidate for refinement
   type :: SFC_candidate_t
      integer(kind=8) :: level                  ! level at which we want to put grid block
      integer(kind=8) :: SFC_id                 ! position at which we want to put grid block
      integer(kind=8), dimension(ndims) :: off  ! offset of grid block (formally can be obtained from SFC_id)
   end type SFC_candidate_t

   !> Refinement map, derefinement request and a list of new grids to create
   type :: ref_flag_t
      logical :: derefine                                                  !> a request to derefine
      type(SFC_candidate_t), allocatable, dimension(:) :: SFC_refine_list  !> A list of children to create; \todo make it private
      logical, private, allocatable, dimension(:,:,:) :: map               !> .true. when a cell triggers refinement criteria, .false. otherwise
   contains
      procedure          :: init           !> Initialize: (.false., allocate 0 elements, clear the map)
      procedure          :: add            !> Appends one element to SFC_refine_list
      procedure          :: initmap        !> Allocate map
      procedure          :: pending_blocks !> A shortcut for size(SFC_refine_list) > 0
      procedure          :: reset_blocks   !> make SFC_refine_list empty

      !> Set the refinement flag in the map. This is the only call supposed to be used in user routines.
      generic,   public  :: set => set_cell, set_arrcell, set_arrng, set_all, set_mask  !, set_range
      procedure, private :: set_cell       !> Mark a cell for refinement (separate indices)
      procedure, private :: set_arrcell    !> Mark a cell for refinement (array indices)
!      procedure, private :: set_range      !> Mark a box of cells for refinement (list of indices)
      procedure, private :: set_arrng      !> Mark a box of cells for refinement (se-type array of indices)
      procedure, private :: set_all        !> Mark all cells for refinement
      procedure, private :: set_mask       !> Mark cells according to mask

      !> Check if any refinement flag is set
      generic,   public  :: get => get_cell, get_any_range, get_any_arrng, get_any
      procedure, private :: get_cell       !> Read refinement status of a cell
      procedure, private :: get_any_range  !> Read refinement status of a box of cells (list of indices)
      procedure, private :: get_any_arrng  !> Read refinement status of a box of cells (se-type array of indices)
      procedure, private :: get_any        !> Read refinement status of all cells

      !> Clear the refinement flag in the map
      generic,   public  :: clear => clear_mask, clear_all
      procedure, private :: clear_mask     !> Unmark cells according to mask
      procedure, private :: clear_all      !> Unmark all cells
   end type ref_flag_t

contains

!> \brief Initialize to (.false., allocate 0 elements, clear the map)

   subroutine init(this)

      use dataio_pub, only: die

      implicit none

      class(ref_flag_t), intent(inout) :: this  !> object invoking this procedure

      this%derefine = .false.
      call this%reset_blocks
      if (.not. allocated(this%map)) call die("[refinement_flag:init] map not allocated")
      call this%clear

   end subroutine init

!> \brief Allocate refinement map

   subroutine initmap(this, ijkse)

      use constants,  only: xdim, ydim, zdim, LO, HI
      use dataio_pub, only: die

      implicit none

      class(ref_flag_t),                            intent(inout) :: this   !> object invoking this procedure
      integer(kind=4), dimension(xdim:zdim, LO:HI), intent(in)    :: ijkse  !> grid bounds, including guardcells

      if (allocated(this%map)) call die("[refinement_flag:initmap] already allocated")
      allocate(this%map(ijkse(xdim, LO):ijkse(xdim, HI), ijkse(ydim, LO):ijkse(ydim, HI), ijkse(zdim, LO):ijkse(zdim, HI)))

   end subroutine initmap

!> \brief Tell whether we have nonempty list of blocks to be created

   logical function pending_blocks(this)

      implicit none

      class(ref_flag_t), intent(inout) :: this  !> object invoking this procedure

      pending_blocks = (size(this%SFC_refine_list) > 0)

   end function pending_blocks

!> \brief Make SFC_refine_list empty

   subroutine reset_blocks(this)

      implicit none

      class(ref_flag_t), intent(inout) :: this  !> object invoking this procedure

      ! This is perhaps a bit suboptimal: todo let cgl%cg%flag%SFC_refine_list(:) grow and implement INVALID entries (will require to change this%pending_blocks() as well).
      if (allocated(this%SFC_refine_list)) deallocate(this%SFC_refine_list)
      allocate(this%SFC_refine_list(0))

   end subroutine reset_blocks

!> \brief  Read refinement status of a cell

   logical function get_cell(this, i, j, k) result (refine)

#ifdef DEBUG
      use dataio_pub, only: die
#endif /* DEBUG */

      implicit none

      class(ref_flag_t), intent(inout) :: this     !> object invoking this procedure
      integer(kind=4),   intent(in)    :: i, j, k  !> cell coordinates

#ifdef DEBUG
      if (any([i, j, k] < lbound(this%map)) .or. any([i, j, k] > ubound(this%map))) call die("[refinement_flag:set_cell] out of range")  ! this can be costly check
#endif /* DEBUG */

      refine = this%map(i, j, k)

   end function get_cell

!> \brief  Read refinement status of all cells

   logical function get_any(this) result (refine)

      implicit none

      class(ref_flag_t), intent(inout) :: this  !> object invoking this procedure

      refine = any(this%map)

   end function get_any

!> \brief  Read refinement status of a box of cells (list of indices)

   logical function get_any_range(this, il, ih, jl, jh, kl, kh, mask) result (any_refine)

#ifdef DEBUG
      use dataio_pub, only: die
#endif /* DEBUG */

      implicit none

      class(ref_flag_t), intent(inout) :: this     !> object invoking this procedure
      integer(kind=8), intent(in)      :: il, jl, kl, ih, jh, kh  !> box coordinates
      logical, dimension(il:ih, jl:jh, kl:kh), optional, intent(in) :: mask  !> the mask to filter out something

#ifdef DEBUG
      if (any([il, jl, kl] < lbound(this%map)) .or. any([il, jl, kl] > ubound(this%map))) call die("[refinement_flag:get_any_range] out of range (l)")  ! this can be costly check
      if (any([ih, jh, kh] < lbound(this%map)) .or. any([ih, jh, kh] > ubound(this%map))) call die("[refinement_flag:get_any_range] out of range (h)")  ! this can be costly check
#endif /* DEBUG */

      if (present(mask)) then
         any_refine = any(this%map(il:ih, jl:jh, kl:kh) .and. mask)
      else
         any_refine = any(this%map(il:ih, jl:jh, kl:kh))
      endif

   end function get_any_range

!> \brief  Read refinement status of a box of cells (se-type array of indices)

   logical function get_any_arrng(this, se, mask) result (any_refine)

      use constants,  only: xdim, ydim, zdim, LO, HI
#ifdef DEBUG
      use dataio_pub, only: die
#endif /* DEBUG */

      implicit none

      class(ref_flag_t), intent(inout) :: this  !> object invoking this procedure
      integer(kind=8), dimension(xdim:zdim, LO:HI), intent(in) :: se  !> box coordinates
      logical, dimension(se(xdim,LO):se(xdim,HI), se(ydim,LO):se(ydim,HI), se(zdim,LO):se(zdim,HI)), optional, intent(in) :: mask  !> the mask to filter out something

#ifdef DEBUG
      if ( any(se(:, LO) > se(:, HI)) .or. &
           any(se(:, LO) < lbound(this%map)) .or. &
           any(se(:, HI) > ubound(this%map))) call die("[refinement_flag:get_any_arrng] out of range")  ! this can be costly check
#endif /* DEBUG */

      if (present(mask)) then
         any_refine = any(this%map(se(xdim,LO):se(xdim,HI), se(ydim,LO):se(ydim,HI), se(zdim,LO):se(zdim,HI)) .and. mask)
      else
         any_refine = any(this%map(se(xdim,LO):se(xdim,HI), se(ydim,LO):se(ydim,HI), se(zdim,LO):se(zdim,HI)))
      endif

   end function get_any_arrng

!> \brief Mark a cell for refinement (separate indices)

   subroutine set_cell(this, i, j, k)

#ifdef DEBUG
      use dataio_pub, only: die
#endif /* DEBUG */

      implicit none

      class(ref_flag_t), intent(inout) :: this     !> object invoking this procedure
      integer(kind=4),   intent(in)    :: i, j, k  !> cell coordinates

#ifdef DEBUG
      if (any([i, j, k] < lbound(this%map)) .or. any([i, j, k] > ubound(this%map))) call die("[refinement_flag:set_cell] out of range")  ! this can be costly check
#endif /* DEBUG */

      this%map(i, j, k) = .true.

   end subroutine set_cell

!> \brief Mark a cell for refinement (array indices)

   subroutine set_arrcell(this, ijk)

      use constants, only: xdim, ydim, zdim

#ifdef DEBUG
      use dataio_pub, only: die
#endif /* DEBUG */

      implicit none

      class(ref_flag_t),                     intent(inout) :: this  !> object invoking this procedure
      integer(kind=8), dimension(xdim:zdim), intent(in)    :: ijk   !> cell coordinates

#ifdef DEBUG
      if (any(ijk < lbound(this%map)) .or. any(ijk > ubound(this%map))) call die("[refinement_flag:set_arrcell] out of range")  ! this can be costly check
#endif /* DEBUG */

      this%map(ijk(xdim), ijk(ydim), ijk(zdim)) = .true.

   end subroutine set_arrcell

!> \brief Mark a box of cells for refinement

!!$   subroutine set_range(this, il, ih, jl, jh, kl, kh)
!!$
!!$#ifdef DEBUG
!!$      use dataio_pub, only: die
!!$#endif /* DEBUG */
!!$
!!$      implicit none
!!$
!!$      class(ref_flag_t), intent(inout) :: this                    !> object invoking this procedure
!!$      integer(kind=4),   intent(in)    :: il, jl, kl, ih, jh, kh  !> box coordinates
!!$
!!$#ifdef DEBUG
!!$      if (any([il, jl, kl] < lbound(this%map)) .or. any([il, jl, kl] > ubound(this%map))) call die("[refinement_flag:set_range] out of range (l)")  ! this can be costly check
!!$      if (any([ih, jh, kh] < lbound(this%map)) .or. any([ih, jh, kh] > ubound(this%map))) call die("[refinement_flag:set_range] out of range (h)")  ! this can be costly check
!!$#endif /* DEBUG */
!!$
!!$      this%map(il:ih, jl:jh, kl:kh) = .true.
!!$
!!$   end subroutine set_range

!> \brief Mark a box of cells for refinement (se-type array of indices)

   subroutine set_arrng(this, se)

      use constants,  only: xdim, ydim, zdim, LO, HI
#ifdef DEBUG
      use dataio_pub, only: die
#endif /* DEBUG */

      implicit none

      class(ref_flag_t), intent(inout) :: this  !> object invoking this procedure
      integer(kind=8), dimension(xdim:zdim, LO:HI), intent(in) :: se  !> box coordinates

#ifdef DEBUG
      if ( any(se(:, LO) > se(:, HI)) .or. &
           any(se(:, LO) < lbound(this%map)) .or. &
           any(se(:, HI) > ubound(this%map))) call die("[refinement_flag:set_arrng] out of range")  ! this can be costly check
#endif /* DEBUG */

      this%map(se(xdim,LO):se(xdim,HI), se(ydim,LO):se(ydim,HI), se(zdim,LO):se(zdim,HI)) = .true.

   end subroutine set_arrng

!> \brief Mark all cells for refinement

   subroutine set_all(this)

      implicit none

      class(ref_flag_t), intent(inout) :: this   !> object invoking this procedure

      this%map = .true.

   end subroutine set_all

!> \brief Mark cells according to mask

   subroutine set_mask(this, mask)

      use dataio_pub, only: die

      implicit none

      class(ref_flag_t),         intent(inout) :: this  !> object invoking this procedure
      logical, dimension(:,:,:), intent(in)    :: mask  !> set flags according to this mask, preserve flags already set

      if (any(shape(mask) /= shape(this%map))) call die("[refinement_flag::set_mask] wrong mask shape")

      this%map = this%map .or. mask

   end subroutine set_mask

!> \brief Unmark cells according to mask

   subroutine clear_mask(this, mask)

      implicit none

      class(ref_flag_t),                      intent(inout) :: this  !> object invoking this procedure
      logical, dimension(:,:,:), allocatable, intent(in)    :: mask  !> set flags according to this mask, preserve flags already set
      ! allocatable is required to get proper bounds, as they were allocated

      if (any(lbound(mask) /= lbound(this%map)) .or. any(ubound(mask) /= ubound(this%map))) then
         ! we expect a leafmap to be passed here, hence we don't calculate overlap (it would be easy, just more lenghty)
         this%map     (lbound(mask, dim=1):ubound(mask, dim=1), lbound(mask, dim=2):ubound(mask, dim=2), lbound(mask, dim=3):ubound(mask, dim=3)) = &
              this%map(lbound(mask, dim=1):ubound(mask, dim=1), lbound(mask, dim=2):ubound(mask, dim=2), lbound(mask, dim=3):ubound(mask, dim=3)) .and. mask
      else
         this%map = this%map .and. mask
      endif

   end subroutine clear_mask

!> \brief Unmark cells according to mask

   subroutine clear_all(this)

      implicit none

      class(ref_flag_t), intent(inout) :: this  !> object invoking this procedure

      this%map = .false.

   end subroutine clear_all

!> \brief Appends one element to SFC_refine_list

   subroutine add(this, level, off, l_off)

      use constants, only: ndims
      use ordering,  only: SFC_order

      implicit none

      class(ref_flag_t),                 intent(inout) :: this   !> object invoking this procedure
      integer(kind=4),                   intent(in)    :: level  !> level at which we want to put grid block
      integer(kind=8), dimension(ndims), intent(in)    :: off    !> offset of grid block
      integer(kind=8), dimension(ndims), intent(in)    :: l_off  !> offset of the level

      this%SFC_refine_list = [ this%SFC_refine_list, SFC_candidate_t(int(level, kind=8), SFC_order(off-l_off), off) ] !lhs reallocation

   end subroutine add

end module refinement_flag
