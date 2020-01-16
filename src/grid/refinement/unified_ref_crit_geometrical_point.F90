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

!> \brief Unified refinement criterion for a single point

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
      procedure          :: mark => mark_point
      procedure, private :: init_lev
   end type urc_point

   interface urc_point
      procedure :: init
   end interface urc_point

   integer(kind=8), parameter :: uninit = huge(1_8)

contains

!> \brief A simple constructor fed by parameters read from problem.par

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
         write(msg, '(a,3g13.5,a)')"[URC point] Initializing refinement at point: [ ", this%coords, " ]"
         call printinfo(msg)
      endif

   end function init

!>
!! \brief Mark a single point in the domain for refinement
!!
!! \details this%iplot is ignored here as not very interesting
!<

   subroutine mark_point(this, cg)

      use constants,  only: xdim, ydim, zdim, LO, HI
      use dataio_pub, only: die
      use grid_cont,  only: grid_container

      implicit none

      class(urc_point),              intent(inout) :: this  !< an object invoking the type-bound procedure
      type(grid_container), pointer, intent(inout) :: cg    !< current grid piece

      if (cg%l%id > this%level) return

      if (.not. allocated(this%ijk)) call die("[unified_ref_crit_geometrical_point:mark_point] ijk not allocated")
      if (any(this%ijk(cg%l%id, :) == uninit)) call this%init_lev  ! new levels of refinement have appears in the meantime

      if (all(this%ijk(cg%l%id, :) >= cg%ijkse(:, LO)) .and. all(this%ijk(cg%l%id, :) <= cg%ijkse(:, HI))) then
         if (cg%l%id < this%level) cg%refinemap(this%ijk(cg%l%id, xdim), this%ijk(cg%l%id, ydim), this%ijk(cg%l%id, zdim)) = .true.
         cg%refine_flags%derefine = .false.  ! this should go one level up (sanitizing)
      endif

   end subroutine mark_point

!>
!! \brief Initialize ths%ijk
!!
!! \details Initialize uning available levels. If more refinements appear then call this again to reinitialize.
!>

   subroutine init_lev(this)

      use cg_level_base,      only: base
      use cg_level_connected, only: cg_level_connected_t
      use constants,          only: LO
      use dataio_pub,         only: printinfo, msg
      use domain,             only: dom
      use mpisetup,           only: master

      implicit none

      class(urc_point), intent(inout)  :: this  !< an object invoking the type-bound procedure

      type(cg_level_connected_t), pointer :: l
      logical, parameter :: verbose = .false.  ! for debugging only

      l => base%level
      do while (associated(l))
         if (l%l%id <= ubound(this%ijk, dim=1)) then
            if (any(this%ijk(l%l%id, :) == uninit)) then
               where (dom%has_dir)
                  this%ijk(l%l%id, :) = l%l%off + floor((this%coords - dom%edge(:, LO))/dom%L_*l%l%n_d, kind=8)
                  ! Excessively large this%coords will result in FPE exception.
                  ! If not trapped then huge() value will be assigned (checked on gfortran 7.3.1), which is safe.
                  ! A wrapped value coming from integer overflow may be unsafe.
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

   end subroutine init_lev

end module unified_ref_crit_geometrical_point
