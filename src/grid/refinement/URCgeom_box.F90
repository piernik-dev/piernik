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

!> \brief Unified refinement criterion for a single box.

module unified_ref_crit_geometrical_box

   use constants,                    only: ndims, LO, HI
   use unified_ref_crit_geometrical, only: urc_geom

   implicit none

   private
   public :: urc_box

!> \brief A type for box refinement

   type, extends(urc_geom) :: urc_box
      real, dimension(ndims, LO:HI)                          :: coords  !< coordinates, where to refine
      integer(kind=8), allocatable, dimension(:, :), private :: ijk_lo  !< integer coordinates of "bottom left corner" at allowed levels; shape: [ base_level_id:this%level-1, ndims ]
      integer(kind=8), allocatable, dimension(:, :), private :: ijk_hi  !< integer coordinates of "top right corner" at allowed levels; shape: [ base_level_id:this%level-1, ndims ]
   contains
      procedure          :: mark => mark_box
      procedure, private :: init_lev
   end type urc_box

   interface urc_box
      procedure :: init
   end interface urc_box

   integer(kind=8), parameter :: uninit = huge(1_8)

contains

!> \brief A simple constructor fed by parameters read from problem.par.

   function init(rp) result(this)

      use constants,  only: base_level_id, ndims, LO, HI, V_VERBOSE
      use dataio_pub, only: printinfo, msg
      use mpisetup,   only: master
      use refinement, only: ref_box

      implicit none

      type(ref_box), intent(in) :: rp  !< the data read from problem.par

      type(urc_box) :: this  !< an object to be constructed

      this%level  = rp%level
      this%coords = rp%coords

      allocate(this%ijk_lo(base_level_id:this%level, ndims), &
           &   this%ijk_hi(base_level_id:this%level, ndims))
      this%ijk_lo = uninit
      this%ijk_hi = uninit
      if (master) then
         write(msg, '(a,3g13.5,a,3g13.5,a,i3)')"[URC box]   Initializing refinement at box: [ ", this%coords(:, LO), " ]..[ ", this%coords(:, HI), " ] at level ", this%level
         call printinfo(msg, V_VERBOSE)
      endif

   end function init

!>
!! \brief Mark a single box in the domain for refinement.
!!
!! \details this%iplot is ignored here because it is not very interesting.
!<

   subroutine mark_box(this, cg)

      use constants,  only: xdim, zdim, LO, HI
      use dataio_pub, only: die
      use grid_cont,  only: grid_container

      implicit none

      class(urc_box),                intent(inout) :: this  !< an object invoking the type-bound procedure
      type(grid_container), pointer, intent(inout) :: cg    !< current grid piece

      integer(kind=8), dimension(xdim:zdim, LO:HI) :: ijk

      if (this%enough_level(cg%l%id)) return

      if (allocated(this%ijk_lo) .neqv. allocated(this%ijk_hi)) call die("[URCgeom_box:mark_box] inconsistent alloc")
      if (.not. allocated(this%ijk_lo)) call die("[URCgeom_box:mark_box] ijk_{lo,hi} not allocated")

      ! Have some new levels of refinement appeared in the meantime?
      if ( any(this%ijk_lo(cg%l%id, :) == uninit) .or. &
           any(this%ijk_hi(cg%l%id, :) == uninit)) call this%init_lev

      if (all(this%ijk_hi(cg%l%id, :) >= cg%ijkse(:, LO)) .and. all(this%ijk_lo(cg%l%id, :) <= cg%ijkse(:, HI))) then
         ijk(:, LO) = min(max(int(this%ijk_lo(cg%l%id, :), kind=4), cg%ijkse(:, LO)), cg%ijkse(:, HI))
         ijk(:, HI) = min(max(int(this%ijk_hi(cg%l%id, :), kind=4), cg%ijkse(:, LO)), cg%ijkse(:, HI))

         call cg%flag%set(ijk)
      endif

   end subroutine mark_box

!>
!! \brief Initialize this%ijk.
!!
!! \details Initialize using available levels. If more refinements appear then call this again to reinitialize.
!>

   subroutine init_lev(this)

      use cg_level_base,      only: base
      use cg_level_connected, only: cg_level_connected_t
      use constants,          only: LO, HI

      implicit none

      class(urc_box), intent(inout)  :: this  !< an object invoking the type-bound procedure

      type(cg_level_connected_t), pointer :: l

      l => base%level
      do while (associated(l))
         if (l%l%id <= ubound(this%ijk_lo, dim=1)) then
            if ( any(this%ijk_lo(l%l%id, :) == uninit) .or. &
                 any(this%ijk_hi(l%l%id, :) == uninit)) then

               this%ijk_lo(l%l%id, :) = this%coord2ind(this%coords(:, LO), l%l)
               this%ijk_hi(l%l%id, :) = this%coord2ind(this%coords(:, HI), l%l)

            endif
         endif
         l => l%finer
      enddo

   end subroutine init_lev

end module unified_ref_crit_geometrical_box
