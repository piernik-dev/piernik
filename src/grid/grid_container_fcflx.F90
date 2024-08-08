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

module grid_cont_fcflx

   use constants,       only: xdim, zdim, LO, HI
   use grid_cont_ref,   only: grid_container_ref_t
   use flx_cell,        only: fluxpoint
   use flx_arr,         only: fluxarray

   implicit none

   private
   public :: grid_container_fcflx_t, cleanup_flxp

   type(fluxpoint), target :: fpl, fpr, cpl, cpr  !< Auxiliary buffers for cell fluxes.
   ! Use with care! Each call to save_outfluxes should be preceded by call to set_fluxpointers.
   ! Do not make another call to set_fluxpointers before completion of save_outfluxes.
   ! The routines set_fluxpointers/save_outfluxes were designed to wrap a 1D solver (line by line execution).

   !> \brief Buffers for f/c flux exchange
   type, extends(grid_container_ref_t), abstract :: grid_container_fcflx_t

      ! f/c boundaries
      type(fluxarray), dimension(xdim:zdim, LO:HI) :: finebnd    !< indices and flux arrays for fine/coarse flux updates on coarse side
      type(fluxarray), dimension(xdim:zdim, LO:HI) :: coarsebnd  !< indices and flux arrays for fine/coarse flux updates on fine side

   contains

      procedure :: init_gc_fcflx       !< Initialization
      procedure :: cleanup_fcflx       !< Deallocate all internals
      procedure :: set_fluxpointers    !< Calculate fluxes incoming from fine grid for 1D solver
      procedure :: save_outfluxes      !< Collect outgoing fine fluxes, do curvilinear scaling and store in appropriate array

   end type grid_container_fcflx_t

contains

!> \brief Initialization of f/c fluxes for the grid container

   subroutine init_gc_fcflx(this)

      use constants, only: xdim, ydim, zdim, LO, HI

      implicit none

      class(grid_container_fcflx_t), target, intent(inout) :: this  !< object invoking type-bound procedure

      integer :: i
      logical, save :: firstcall = .true.

      ! For simplicity we allocate all possible buffers for f/c fluxes.
      ! Some of this memory remains unused and may get swapped out.
      ! Index values outside of the block mark no flux to process.
      do i = LO, HI
         call this%finebnd  (xdim, i)%init([ this%js, this%je ], [ this%ks, this%ke ])
         call this%finebnd  (ydim, i)%init([ this%ks, this%ke ], [ this%is, this%ie ])
         call this%finebnd  (zdim, i)%init([ this%is, this%ie ], [ this%js, this%je ])
         call this%coarsebnd(xdim, i)%init([ this%js, this%je ], [ this%ks, this%ke ])
         call this%coarsebnd(ydim, i)%init([ this%ks, this%ke ], [ this%is, this%ie ])
         call this%coarsebnd(zdim, i)%init([ this%is, this%ie ], [ this%js, this%je ])
      enddo
      do i = xdim, zdim
         this%finebnd  (i, LO)%index = this%lhn(i, LO) - 1
         this%finebnd  (i, HI)%index = this%lhn(i, HI) + 1
         this%coarsebnd(i, LO)%index = this%lhn(i, LO) - 1
         this%coarsebnd(i, HI)%index = this%lhn(i, HI) + 1
      enddo

      ! These pointers aren't thread-safe
      if (firstcall) then
         call fpl%init
         call fpr%init
         call cpl%init
         call cpr%init
         firstcall = .false.
      endif

   end subroutine init_gc_fcflx

!> \brief Routines that deallocates all internals of the grid container

   subroutine cleanup_fcflx(this)

      use constants, only: xdim, zdim, LO, HI

      implicit none

      class(grid_container_fcflx_t), intent(inout) :: this  !< object invoking type-bound procedure

      integer :: d, g

      do d = xdim, zdim
         do g = LO, HI
            call this%finebnd  (d, g)%cleanup
            call this%coarsebnd(d, g)%cleanup
         enddo
      enddo

   end subroutine cleanup_fcflx

   subroutine cleanup_flxp

      implicit none

      call fpl%cleanup
      call fpr%cleanup
      call cpl%cleanup
      call cpr%cleanup

   end subroutine cleanup_flxp

!>
!! \brief Calculate fluxes incoming from fine grid for 1D solver
!!
!! Somewhat spaghetti style
!!
!! The index has to be corrected to match the way of allocation of flux arrays in the RTVD and Riemann solver.
!<

   subroutine set_fluxpointers(this, cdim, i1, i2, eflx)

      use constants,  only: LO, HI, ydim, zdim, GEO_RPZ
      use domain,     only: dom
      use fluidindex, only: iarr_all_mx, iarr_all_my
      use fluxtypes,  only: ext_fluxes

      implicit none

      class(grid_container_fcflx_t), intent(in)    :: this  !< object invoking type-bound procedure
      integer(kind=4),               intent(in)    :: cdim  !< direction of the flux
      integer,                       intent(in)    :: i1    !< coordinate_1, perpendicular to the f/c face
      integer,                       intent(in)    :: i2    !< coordinate_2, perpendicular to the f/c face
      type(ext_fluxes),              intent(inout) :: eflx  !< fluxes stored for the selected cell

      call eflx%init

      if (this%finebnd(cdim, LO)%index(i1, i2) >= this%ijkse(cdim, LO)) then
         call this%finebnd(cdim, LO)%fa2fp(fpl, i1, i2)
         eflx%li => fpl
         eflx%li%index = eflx%li%index - this%lhn(cdim, LO) + 1
      endif
      if (this%finebnd(cdim, HI)%index(i1, i2) <= this%ijkse(cdim, HI)) then
         call this%finebnd(cdim, HI)%fa2fp(fpr, i1, i2)
         eflx%ri => fpr
         eflx%ri%index = eflx%ri%index - this%lhn(cdim, LO)
      endif

      if (dom%geometry_type == GEO_RPZ) then
         if (cdim == ydim) then
            !> BEWARE: iarr_all_mx points to the y-momentum in y-sweep
            if (associated(eflx%li)) eflx%li%uflx(iarr_all_mx) = eflx%li%uflx(iarr_all_mx) / this%x(i2)
            if (associated(eflx%ri)) eflx%ri%uflx(iarr_all_mx) = eflx%ri%uflx(iarr_all_mx) / this%x(i2)
         else if (cdim == zdim) then
            if (associated(eflx%li)) then
               eflx%li%uflx = eflx%li%uflx / this%x(i1)
               eflx%li%uflx(iarr_all_my) = eflx%li%uflx(iarr_all_my) / this%x(i1) ! that makes this%x(i1)**2
            endif
            if (associated(eflx%ri)) then
               eflx%ri%uflx = eflx%ri%uflx / this%x(i1)
               eflx%ri%uflx(iarr_all_my) = eflx%ri%uflx(iarr_all_my) / this%x(i1) ! that makes this%x(i1)**2
            endif
         endif
      endif

      if (this%coarsebnd(cdim, LO)%index(i1, i2) >= this%ijkse(cdim, LO)) then
         cpl%index = this%coarsebnd(cdim, LO)%index(i1, i2) - this%lhn(cdim, LO)
         eflx%lo => cpl
      endif
      if (this%coarsebnd(cdim, HI)%index(i1, i2) <= this%ijkse(cdim, HI)) then
         cpr%index = this%coarsebnd(cdim, HI)%index(i1, i2) - this%lhn(cdim, LO) + 1
         eflx%ro => cpr
      endif

   end subroutine set_fluxpointers

!> \brief Collect outgoing fine fluxes from 1D solver, do curvilinear scaling and store in appropriate array

   subroutine save_outfluxes(this, cdim, i1, i2, eflx)

      use constants, only: LO, HI
      use fluxtypes, only: ext_fluxes

      implicit none

      class(grid_container_fcflx_t), intent(inout) :: this  !< object invoking type-bound procedure
      integer(kind=4),               intent(in)    :: cdim  !< direction of the flux
      integer,                       intent(in)    :: i1    !< coordinate_1, perpendicular to the f/c face
      integer,                       intent(in)    :: i2    !< coordinate_2, perpendicular to the f/c face
      type(ext_fluxes),              intent(inout) :: eflx  !< fluxes stored for the selected cell

      if (associated(eflx%lo)) then
         eflx%lo%index = eflx%lo%index + this%lhn(cdim, LO)
         call cyl_scale(eflx%lo)
         call this%coarsebnd(cdim, LO)%fp2fa(eflx%lo, i1, i2)
      endif

      if (associated(eflx%ro)) then
         eflx%ro%index = eflx%ro%index + this%lhn(cdim, LO) - 1
         call cyl_scale(eflx%ro)
         call this%coarsebnd(cdim, HI)%fp2fa(eflx%ro, i1, i2)
      endif

   contains

      ! It would be nice to have this routine as fluxpoint::cyl_scale() or inside fp2fa()
      ! Unfortunately it requires cdim and this%x(i1) or this%x(i2)

      subroutine cyl_scale(fp)

         use constants,  only: ydim, zdim, GEO_RPZ
         use domain,     only: dom
         use fluidindex, only: iarr_all_mx, iarr_all_my
         use flx_cell,   only: fluxpoint

         implicit none

         type(fluxpoint), intent(inout) :: fp

         if (dom%geometry_type == GEO_RPZ) then
            if (cdim == ydim) then
               !> BEWARE: iarr_all_mx points to the y-momentum in the y-sweep
               fp%uflx(iarr_all_mx) = fp%uflx(iarr_all_mx) * this%x(i2)
            else if (cdim == zdim) then
               fp%uflx = fp%uflx * this%x(i1)
               fp%uflx(iarr_all_my) = fp%uflx(iarr_all_my) * this%x(i1) ! that makes this%x(i1)**2
            endif
         endif

      end subroutine cyl_scale

   end subroutine save_outfluxes

end module grid_cont_fcflx
