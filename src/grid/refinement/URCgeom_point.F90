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

!> \brief Unified refinement criterion for a single point.

module unified_ref_crit_geometrical_point

   use constants,                    only: ndims
   use unified_ref_crit_geometrical, only: urc_geom

   implicit none

   private
   public :: urc_point

!> \brief A type for point refinement

   type, extends(urc_geom) :: urc_point
      real, dimension(ndims)                                 :: coords  !< coordinates, where to refine
      integer(kind=8), allocatable, dimension(:, :), private :: ijk     !< integer coordinates at allowed levels; shape: [ base_level_id:this%level-1, ndims ]
   contains
      procedure          :: mark => mark_point
      procedure, private :: init_lev
   end type urc_point

   interface urc_point
      procedure :: init
   end interface urc_point

   integer(kind=4), parameter :: uninit = huge(1_4)

contains

!> \brief A simple constructor fed by parameters read from problem.par.

   function init(rp) result(this)

      use constants,  only: base_level_id, ndims
      use dataio_pub, only: printinfo, msg
      use mpisetup,   only: master
      use refinement, only: ref_point

      implicit none

      type(ref_point), intent(in) :: rp  !< the data read from problem.par

      type(urc_point) :: this  !< an object to be constructed

      this%level  = rp%level
      this%coords = rp%coords

      allocate(this%ijk(base_level_id:this%level, ndims))
      this%ijk = uninit
      if (master) then
         write(msg, '(a,3g13.5,a,i3)')"[URC point] Initializing refinement at point: [ ", this%coords, " ] at level ", this%level
         call printinfo(msg)
      endif

   end function init

!>
!! \brief Mark a single point in the domain for refinement.
!!
!! \details this%iplot is ignored here because it is not very interesting.
!<

   subroutine mark_point(this, cg)

      use constants,  only: LO, HI
      use dataio_pub, only: die
      use grid_cont,  only: grid_container

      implicit none

      class(urc_point),              intent(inout) :: this  !< an object invoking the type-bound procedure
      type(grid_container), pointer, intent(inout) :: cg    !< current grid piece

      if (this%enough_level(cg%l%id)) return

      if (.not. allocated(this%ijk)) call die("[URCgeom_point:mark_point] ijk not allocated")

      ! Have some new levels of refinement appeared in the meantime?
      if (any(this%ijk(cg%l%id, :) == uninit)) call this%init_lev

      if (all(this%ijk(cg%l%id, :) >= cg%ijkse(:, LO)) .and. all(this%ijk(cg%l%id, :) <= cg%ijkse(:, HI))) &
           call cg%flag%set(this%ijk(cg%l%id, :))

   end subroutine mark_point

!>
!! \brief Initialize this%ijk.
!!
!! \details Initialize using available levels. If more refinements appear then call this again to reinitialize.
!>

   subroutine init_lev(this)

      use cg_level_base,      only: base
      use cg_level_connected, only: cg_level_connected_t

      implicit none

      class(urc_point), intent(inout)  :: this  !< an object invoking the type-bound procedure

      type(cg_level_connected_t), pointer :: l

      l => base%level
      do while (associated(l))
         if (l%l%id <= ubound(this%ijk, dim=1)) then
            if (any(this%ijk(l%l%id, :) == uninit)) &
                 this%ijk(l%l%id, :) = this%coord2ind(this%coords, l%l)
         endif
         l => l%finer
      enddo

   end subroutine init_lev

end module unified_ref_crit_geometrical_point
