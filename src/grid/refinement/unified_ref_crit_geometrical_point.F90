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
!    along with PIERNIK.  If not, see http://www.gnu.org/licenses/.
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

!> \brief Unified refinement criterium for a single point

module unified_ref_crit_geometrical_point

   use constants,                    only: ndims
   use unified_ref_crit_geometrical, only: urc_geom

   implicit none

   private
   public :: urc_point

!> \brief A type for point refinement

   type, extends(urc_geom) :: urc_point
      real, dimension(ndims) :: coords  !< coordinates, where to refine
      integer(kind=8), allocatable, dimension(:, :), private :: ijk  !< integer coordinates at allowed levels; shape: [ base_level_id:this%level-1, ndims ]
   contains
      procedure :: mark    => mark_point
      procedure :: cleanup => cleanup_point
      procedure :: init    => init_point
   end type urc_point

   integer(kind=8), parameter :: uninit = huge(1_8)

contains

!>
!! \brief Mark a single point in the domain for refinement
!!
!! \details this%iplot is ignored here as not very interesting
!<

   subroutine mark_point(this, cg)

      use constants, only: xdim, ydim, zdim, LO, HI
      use grid_cont, only: grid_container

      implicit none

      class(urc_point),              intent(inout) :: this  !< an object invoking the type-bound procedure
      type(grid_container), pointer, intent(inout) :: cg    !< current grid piece

      if (cg%l%id >= this%level) return

      if (.not. allocated(this%ijk)) call this%init
      if (any(this%ijk(cg%l%id, :) == uninit)) call this%init  ! new levels of refinement have appeares in the meantime

      if (all(this%ijk(cg%l%id, :) >= cg%ijkse(:, LO)) .and. all(this%ijk(cg%l%id, :) <= cg%ijkse(:, HI))) then
         cg%refinemap(this%ijk(cg%l%id, xdim), this%ijk(cg%l%id, ydim), this%ijk(cg%l%id, zdim)) = .true.
         cg%refine_flags%derefine = .false.  ! this should go one level up (sanitizing)
      endif

   end subroutine mark_point

!>
!! \brief Initialize ths%ijk
!!
!! \details Initialize uning available levels. If more refinements appear then call this again to reinitialize.
!>

   subroutine init_point(this)

      use cg_level_base,      only: base
      use cg_level_connected, only: cg_level_connected_T
      use constants,          only: base_level_id, ndims, LO
      use dataio_pub,         only: printinfo, msg
      use domain,             only: dom
      use mpisetup,           only: master

      implicit none

      class(urc_point), intent(inout)  :: this  !< an object invoking the type-bound procedure

      type(cg_level_connected_T), pointer :: l
      logical, parameter :: verbose = .false.

      if (allocated(this%ijk)) then
         if (verbose .and. master) call printinfo("[unified_ref_crit_geometrical_point:init_point] Re-initializing")
      else
         allocate(this%ijk(base_level_id:this%level-1, ndims))
         this%ijk = uninit
         if (master) then
            write(msg, '(a,3g13.5,a)')"[URC point] Initializing refinement at point: [ ", this%coords, " ]"
            call printinfo(msg)
         endif
      endif

      l => base%level
      do while (associated(l))
         if (l%l%id < this%level) then
            if (any(this%ijk(l%l%id, :) == uninit)) then
               where (dom%has_dir)
                  this%ijk(l%l%id, :) = l%l%off + floor((this%coords - dom%edge(:, LO))/dom%L_*l%l%n_d)
               elsewhere
                  this%ijk(l%l%id, :) = l%l%off
               endwhere
               if (verbose .and. master) then
                  write(msg, '(a,i3,a,3i8,a)')"[URC point]   point coordinates at level ", l%l%id, " are: [ ", this%ijk(l%l%id, :), " ]"
                  call printinfo(msg)
               endif
            endif
         endif
         l => l%finer
      enddo

   end subroutine init_point

!> \brief Deallocate integer coordinates of the point

   subroutine cleanup_point(this)

      implicit none

      class(urc_point), intent(inout)  :: this  !< an object invoking the type-bound procedure

      if (allocated(this%ijk)) deallocate(this%ijk)

   end subroutine cleanup_point

end module unified_ref_crit_geometrical_point
