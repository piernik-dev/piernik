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

!> \brief Unified refinement criterion for a cylinder with axis in the z-direction that fits to a simple box

module unified_ref_crit_geometrical_zcyl

   use constants,                    only: ndims, LO, HI
   use unified_ref_crit_geometrical, only: urc_geom

   implicit none

   private
   public :: urc_zcyl

!> \brief A type for box refinement

   type, extends(urc_geom) :: urc_zcyl
      real, dimension(ndims, LO:HI)                          :: coords  !< coordinates of the box that contains the cylinder, where to refine
      real, dimension(ndims), private                        :: center  !< coordinates of the center of the cylinder
      real, dimension(ndims), private                        :: size    !< x- and y-radius (elliptical cylinders are allowed) and half-height
      integer(kind=8), allocatable, dimension(:, :), private :: ijk_lo  !< integer coordinates of "bottom left corner" at allowed levels; shape: [ base_level_id:this%level-1, ndims ]
      integer(kind=8), allocatable, dimension(:, :), private :: ijk_hi  !< integer coordinates of "top right corner" at allowed levels; shape: [ base_level_id:this%level-1, ndims ]
   contains
      procedure          :: mark => mark_zcyl
      procedure, private :: init_lev
   end type urc_zcyl

   interface urc_zcyl
      procedure :: init
   end interface urc_zcyl

   integer(kind=8), parameter :: uninit = huge(1_8)

contains

!> \brief A simple constructor fed by parameters read from problem.par

   function init(rp) result(this)

      use constants,  only: base_level_id, ndims, LO, HI
      use dataio_pub, only: printinfo, msg
      use domain,     only: dom
      use mpisetup,   only: master
      use refinement, only: ref_box

      implicit none

      type(ref_box), intent(in) :: rp  !< the data read from problem.par

      type(urc_zcyl) :: this  !< an object to be constructed

      this%level  = rp%level
      this%coords = rp%coords

      where (dom%has_dir)
         this%center = (this%coords(:, HI) + this%coords(:, LO)) / 2.
         this%size   = (this%coords(:, HI) - this%coords(:, LO)) / 2.
      elsewhere
         this%center = dom%C_
         this%size   = dom%L_
      endwhere

      allocate(this%ijk_lo(base_level_id:this%level, ndims), this%ijk_hi(base_level_id:this%level, ndims))
      this%ijk_lo = uninit
      this%ijk_hi = uninit
      if (master) then
         write(msg, '(a,3g13.5,a,3g13.5,a)')"[URC zcyl]   Initializing cylindrical refinement in a box: [ ", this%coords(:, LO), " ]..[ ", this%coords(:, HI), " ]"
         call printinfo(msg)
      endif

   end function init

!>
!! \brief Mark a z-oriented cylinted inside a specified box in the domain for refinement
!!
!! Please note that it uses grid topology, so in GEO_RPZ may result in funny shapes, especially for "cylinders" extended in angular direction.
!!
!! \details this%iplot is ignored here as not very interesting
!<

   subroutine mark_zcyl(this, cg)

      use constants,  only: xdim, ydim, zdim, LO, HI
      use dataio_pub, only: die
      use grid_cont,  only: grid_container

      implicit none

      class(urc_zcyl),               intent(inout) :: this  !< an object invoking the type-bound procedure
      type(grid_container), pointer, intent(inout) :: cg    !< current grid piece

      real :: dist2_y, dist2_x
      integer :: i, j, k

      if (cg%l%id > this%level) return

      if (allocated(this%ijk_lo) .neqv. allocated(this%ijk_hi)) call die("[unified_ref_crit_geometrical_box:mark_zcyl] inconsistent alloc")

      if (.not. allocated(this%ijk_lo)) call die("[unified_ref_crit_geometrical_box:mark_zcyl] ijk_{lo,hi} not allocated")
      if (any(this%ijk_lo(cg%l%id, :) == uninit) .or. any(this%ijk_hi(cg%l%id, :) == uninit)) call this%init_lev  ! new levels of refinement have appeared in the meantime

      if (all(this%ijk_hi(cg%l%id, :) >= cg%ijkse(:, LO)) .and. all(this%ijk_lo(cg%l%id, :) <= cg%ijkse(:, HI))) then
         do k = cg%ks, cg%ke
            if ((cg%z(k) >= this%coords(zdim, LO)) .and. (cg%z(k) <= this%coords(zdim, HI))) then
               do j = cg%js, cg%je
                  dist2_y = ((cg%y(j) - this%center(ydim)) / this%size(ydim))**2
                  do i = cg%is, cg%ie
                     dist2_x = ((cg%x(i) - this%center(xdim)) / this%size(xdim))**2
                     if (dist2_x + dist2_y <= 1.) call cg%flag%set(i, j, k)
                  enddo
               enddo
            endif
         enddo
      endif

   end subroutine mark_zcyl

!>
!! \brief Initialize ths%ijk
!!
!! \details Initialize using available levels. If more refinements appear then call this again to reinitialize.
!!
!! Spaghetti warning: practically identical to unified_ref_crit_geometrical_box:init_lev
!>

   subroutine init_lev(this)

      use cg_level_base,      only: base
      use cg_level_connected, only: cg_level_connected_t
      use constants,          only: LO, HI
      use dataio_pub,         only: printinfo, msg
      use domain,             only: dom
      use mpisetup,           only: master

      implicit none

      class(urc_zcyl), intent(inout)  :: this  !< an object invoking the type-bound procedure

      type(cg_level_connected_t), pointer :: l
      logical, parameter :: verbose = .false.  ! for debugging only

      l => base%level
      do while (associated(l))
         if (l%l%id <= ubound(this%ijk_lo, dim=1)) then
            if (any(this%ijk_lo(l%l%id, :) == uninit) .or. any(this%ijk_hi(l%l%id, :) == uninit)) then
               where (dom%has_dir)
                  this%ijk_lo(l%l%id, :) = l%l%off + floor((this%coords(:, LO) - dom%edge(:, LO))/dom%L_*l%l%n_d, kind=8)
                  this%ijk_hi(l%l%id, :) = l%l%off + floor((this%coords(:, HI) - dom%edge(:, LO))/dom%L_*l%l%n_d, kind=8)
                  ! Excessively large this%coords will result in FPE exception.
                  ! If not trapped then huge() value will be assigned (checked on gfortran 7.3.1), which is safe.
                  ! A wrapped value coming from integer overflow may be unsafe.
               elsewhere
                  this%ijk_lo(l%l%id, :) = l%l%off
                  this%ijk_hi(l%l%id, :) = l%l%off
               endwhere
               if (verbose .and. master) then
                  write(msg, '(a,i3,a,3i8,a,3i8,a)')"[URC zcyl]  z-cylinder coordinates at level ", l%l%id, " are: [ ", this%ijk_lo(l%l%id, :), " ]..[ ", this%ijk_hi(l%l%id, :), " ]"
                  call printinfo(msg)
               endif
            endif
         endif
         l => l%finer
      enddo

   end subroutine init_lev

end module unified_ref_crit_geometrical_zcyl
