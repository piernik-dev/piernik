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

!> \brief Module extending grid container type by structures required for boundary exchanges (same level and f/c)

module grid_cont_ref

   use grid_cont_bseg,  only: grid_container_bseg_t
   use refinement_flag, only: ref_flag_t

   implicit none

   private
   public :: grid_container_ref_t

   !> \brief Arrays in grid container that are required for refinement
   type, extends(grid_container_bseg_t), abstract :: grid_container_ref_t

      ! Refinements
      logical, allocatable, dimension(:,:,:) :: leafmap  !< .true. when a cell is not covered by finer cells, .false. otherwise
      type(ref_flag_t) :: flag                           !< refine or derefine this grid container?

      ! Due to the memory cost, this%leafmap should be used only as a mask in array operations. It is safe but a bit too general performance-wise.
      ! \todo Implement a smaller array for oct-tree indicators of children coverage.

   contains

      procedure :: init_gc_ref         !< Initialization
      procedure :: cleanup_ref         !< Deallocate all internals
      procedure :: refinemap2SFC_list  !< Create list of SFC indices to be created from refine flags
      procedure :: has_leaves          !< Returns .true. if there are any non-covered cells on this

   end type grid_container_ref_t

contains

!> \brief Initialization of f/c fluxes for the grid container

   subroutine init_gc_ref(this)

      use constants, only: xdim, ydim, zdim, LO, HI

      implicit none

      class(grid_container_ref_t), target, intent(inout) :: this  !< object invoking type-bound procedure

      allocate(this%leafmap(this%ijkse(xdim, LO):this%ijkse(xdim, HI), &
           &                this%ijkse(ydim, LO):this%ijkse(ydim, HI), &
           &                this%ijkse(zdim, LO):this%ijkse(zdim, HI)))

      this%leafmap(:, :, :) = .true.
      call this%flag%initmap(this%lhn)

   end subroutine init_gc_ref

!> \brief Routines that deallocates all internals of the grid container

   subroutine cleanup_ref(this)

      implicit none

      class(grid_container_ref_t), intent(inout) :: this  !< object invoking type-bound procedure

      ! arrays not handled through named_array feature
      if (allocated(this%leafmap)) deallocate(this%leafmap)

   end subroutine cleanup_ref

!< \brief Create list of SFC indices to be created from refine flags

   subroutine refinemap2SFC_list(this)

      use constants,  only: refinement_factor, xdim, ydim, zdim, I_ONE
      use domain,     only: dom
      use refinement, only: bsize

      implicit none

      class(grid_container_ref_t), intent(inout) :: this !< object invoking type-bound procedure

      integer(kind=8) :: i, j, k, ifs, ife, jfs, jfe, kfs, kfe

      call this%flag%clear(this%leafmap) ! no parent correction possible beyond this point

      if (.not. this%flag%get()) return

      if (any((bsize <= 0) .and. dom%has_dir)) return ! this routine works only with blocky AMR

      !! ToDo: precompute refinement decomposition in this%init_gc and simplify the code below.
      !! It should also simplify decomposition management and make it more flexible in case we decide to work on uneven AMR blocks

      !! beware: consider dropping this%l%off feature for simplicity. It will require handling the shift due to domain expansion (some increase CPU cost)

      associate( b_size => merge(bsize, huge(I_ONE), dom%has_dir))
         do i = int(((this%is - this%l%off(xdim))*refinement_factor) / b_size(xdim)), int(((this%ie - this%l%off(xdim))*refinement_factor + I_ONE) / b_size(xdim))
            ifs = max(int(this%is, kind=8), int(this%l%off(xdim)) + (i*b_size(xdim))/refinement_factor)
            ife = min(int(this%ie, kind=8), int(this%l%off(xdim)) + ((i+I_ONE)*b_size(xdim)-I_ONE)/refinement_factor)

            do j = int(((this%js - this%l%off(ydim))*refinement_factor) / b_size(ydim)), int(((this%je - this%l%off(ydim))*refinement_factor + I_ONE) / b_size(ydim))
               jfs = max(int(this%js, kind=8), int(this%l%off(ydim)) + (j*b_size(ydim))/refinement_factor)
               jfe = min(int(this%je, kind=8), int(this%l%off(ydim)) + ((j+I_ONE)*b_size(ydim)-I_ONE)/refinement_factor)

               do k = int(((this%ks - this%l%off(zdim))*refinement_factor) / b_size(zdim)), int(((this%ke - this%l%off(zdim))*refinement_factor + I_ONE) / b_size(zdim))
                  kfs = max(int(this%ks, kind=8), int(this%l%off(zdim)) + (k*b_size(zdim))/refinement_factor)
                  kfe = min(int(this%ke, kind=8), int(this%l%off(zdim)) + ((k+I_ONE)*b_size(zdim)-I_ONE)/refinement_factor)
                  if (this%flag%get(ifs, ife, jfs, jfe, kfs, kfe)) call this%flag%add(this%l%id+I_ONE, [i, j, k]*b_size + refinement_factor*this%l%off, refinement_factor*this%l%off)
               enddo
            enddo
         enddo
      end associate
      call this%flag%clear

   end subroutine refinemap2SFC_list

!>
!! \brief Returns .true. if there are any non-covered cells on this
!!
!! \details This is quite general but expensive test. For an oct-tree
!! (complete or incomplete) it would be enough to check only 8 values.
!!
!! \todo Implement cheaper variant together with more metadata in %dot or similar structures.
!<

   pure logical function has_leaves(this)

      implicit none

      class(grid_container_ref_t), intent(in) :: this

      has_leaves = any(this%leafmap)

   end function has_leaves

end module grid_cont_ref
