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

module grid_cont_bnd

   use constants,    only: xdim, zdim, LO, HI
   use grid_cont_na, only: grid_container_na_t
   use fluxtypes,    only: fluxarray, fluxpoint

   implicit none

   private
   public :: grid_container_bnd_t, segment

   type(fluxpoint), target :: fpl, fpr, cpl, cpr

   !> \brief Specification of segment of data for boundary exchange
   type :: segment
      integer :: proc                                      !< target process
      integer(kind=8), dimension(xdim:zdim, LO:HI) :: se   !< range
      integer(kind=4) :: tag                               !< unique tag for data exchange
      real, allocatable, dimension(:,:,:)   :: buf         !< buffer for the 3D (scalar) data to be sent or received
      real, allocatable, dimension(:,:,:,:) :: buf4        !< buffer for the 4D (vector) data to be sent or received
      integer(kind=4), pointer :: req                      !< request ID, used for most asynchronous communication, such as fine-coarse flux exchanges
      integer(kind=8), dimension(xdim:zdim, LO:HI) :: se2  !< auxiliary range, used in cg_level_connected:vertical_bf_prep
      class(grid_container_bnd_t), pointer :: local        !< set this pointer to non-null when the exchange is local
   end type segment

   !> \brief Array of boundary segments to exchange
   type :: bnd_list
      type(segment), dimension(:), allocatable :: seg !< segments
   contains
      procedure :: add_seg !< Add an new segment, reallocate if necessary
   end type bnd_list

   !> \brief Everything required for autonomous computation of a single sweep on a portion of the domain on a single process
   type, extends(grid_container_na_t), abstract :: grid_container_bnd_t

      ! External boundary conditions and internal boundaries

      type(fluxarray), dimension(xdim:zdim, LO:HI) :: finebnd     !< indices and flux arrays for fine/coarse flux updates on coarse side
      type(fluxarray), dimension(xdim:zdim, LO:HI) :: coarsebnd   !< indices and flux arrays for fine/coarse flux updates on fine side

      ! initialization of i_bnd and o_bnd are done in cg_list_neighbors because we don't have access to cg_level%dot here
      type(bnd_list),  dimension(:), allocatable   :: i_bnd       !< description of incoming boundary data
      type(bnd_list),  dimension(:), allocatable   :: o_bnd       !< description of outgoing boundary data

   contains

      procedure          :: init_gc_bnd       !< Initialization
      procedure          :: cleanup_bnd       !< Deallocate all internals
      procedure          :: set_fluxpointers  !< Calculate fluxes incoming from fine grid for 1D solver
      procedure          :: save_outfluxes    !< Collect outgoing fine fluxes, do curvilinear scaling and store in appropriate array

   end type grid_container_bnd_t

contains

!> \brief Initialization of f/c fluxes for the grid container

   subroutine init_gc_bnd(this)

      use constants, only: xdim, ydim, zdim, LO, HI

      implicit none

      class(grid_container_bnd_t), target, intent(inout) :: this  !< object invoking type-bound procedure

      integer :: i

      do i = LO, HI
         call this%finebnd  (xdim, i)%fainit([ this%js, this%je ], [ this%ks, this%ke ])
         call this%finebnd  (ydim, i)%fainit([ this%ks, this%ke ], [ this%is, this%ie ])
         call this%finebnd  (zdim, i)%fainit([ this%is, this%ie ], [ this%js, this%je ])
         call this%coarsebnd(xdim, i)%fainit([ this%js, this%je ], [ this%ks, this%ke ])
         call this%coarsebnd(ydim, i)%fainit([ this%ks, this%ke ], [ this%is, this%ie ])
         call this%coarsebnd(zdim, i)%fainit([ this%is, this%ie ], [ this%js, this%je ])
      enddo
      do i = xdim, zdim
         this%finebnd  (i, LO)%index = this%lhn(i, LO) - 1
         this%finebnd  (i, HI)%index = this%lhn(i, HI) + 1
         this%coarsebnd(i, LO)%index = this%lhn(i, LO) - 1
         this%coarsebnd(i, HI)%index = this%lhn(i, HI) + 1
      enddo

   end subroutine init_gc_bnd

!> \brief Routines that deallocates all internals of the grid container

   subroutine cleanup_bnd(this)

      use constants, only: xdim, zdim, LO, HI

      implicit none

      class(grid_container_bnd_t), intent(inout) :: this  !< object invoking type-bound procedure

      integer :: d, g

      if (allocated(this%i_bnd)) then
         do d = lbound(this%i_bnd, dim=1), ubound(this%i_bnd, dim=1)
            if (allocated(this%i_bnd(d)%seg)) deallocate(this%i_bnd(d)%seg)
         enddo
         deallocate(this%i_bnd)
      endif
      if (allocated(this%o_bnd)) then
         do d = lbound(this%o_bnd, dim=1), ubound(this%o_bnd, dim=1)
            if (allocated(this%o_bnd(d)%seg)) deallocate(this%o_bnd(d)%seg)
         enddo
         deallocate(this%o_bnd)
      endif

      do d = xdim, zdim
         do g = LO, HI
            call this%finebnd  (d, g)%facleanup
            call this%coarsebnd(d, g)%facleanup
         enddo
      enddo

      call fpl%fpcleanup
      call fpr%fpcleanup
      call cpl%fpcleanup
      call cpr%fpcleanup

   end subroutine cleanup_bnd

!> \brief Add a new segment, reallocate if necessary

   subroutine add_seg(this, proc, se, tag)

      use dataio_pub, only: die

      implicit none

      class(bnd_list),                              intent(inout) :: this !< object invoking type-bound procedure
      integer,                                      intent(in)    :: proc !< process to be communicated
      integer(kind=8), dimension(xdim:zdim, LO:HI), intent(in)    :: se   !< segment definition
      integer(kind=4),                              intent(in)    :: tag  !< tag for MPI calls

      type(segment), dimension(:), allocatable :: tmp
      integer :: g

      if (tag <0) call die("[grid_container_bnd:add_seg] tag<0")

      if (allocated(this%seg)) then
         allocate(tmp(lbound(this%seg, dim=1):ubound(this%seg, dim=1)+1))
         tmp(:ubound(this%seg, dim=1)) = this%seg
         call move_alloc(from=tmp, to=this%seg)
      else
         allocate(this%seg(1))
      endif

      g = ubound(this%seg, dim=1)

      this%seg(g)%proc = proc
      this%seg(g)%se = se
      this%seg(g)%tag = tag
      nullify(this%seg(g)%local)

   end subroutine add_seg

!> \brief Calculate fluxes incoming from fine grid for 1D solver

   subroutine set_fluxpointers(this, cdim, i1, i2, eflx)

      use constants,  only: LO, HI, ydim, zdim, GEO_RPZ
      use domain,     only: dom
      use fluidindex, only: iarr_all_mx, iarr_all_my
      use fluxtypes,  only: ext_fluxes

      implicit none

      class(grid_container_bnd_t), intent(in)    :: this    !< object invoking type-bound procedure
      integer(kind=4),             intent(in)    :: cdim    !< direction of the flux
      integer,                     intent(in)    :: i1      !< coordinate
      integer,                     intent(in)    :: i2      !< coordinate
      type(ext_fluxes),            intent(inout) :: eflx

      if (this%finebnd(cdim, LO)%index(i1, i2) >= this%ijkse(cdim, LO)) then
         fpl = this%finebnd(cdim, LO)%fa2fp(i1, i2)
         if (.not. allocated(fpl%uflx)) call fpl%fpinit
         eflx%li => fpl
         eflx%li%index = eflx%li%index - this%lhn(cdim, LO) + 1
      else
         nullify(eflx%li)
      endif
      if (this%finebnd(cdim, HI)%index(i1, i2) <= this%ijkse(cdim, HI)) then
         fpr = this%finebnd(cdim, HI)%fa2fp(i1, i2)
         if (.not. allocated(fpr%uflx)) call fpr%fpinit
         eflx%ri => fpr
         eflx%ri%index = eflx%ri%index - this%lhn(cdim, LO)
      else
         nullify(eflx%ri)
      endif

      if (dom%geometry_type == GEO_RPZ) then
         if (cdim == ydim) then
            !> BEWARE: iarr_all_mx points to the y-momentum in y-sweep
            if (associated(eflx%li)) eflx%li%uflx(iarr_all_mx) = eflx%li%uflx(iarr_all_mx) / this%x(i2)
            if (associated(eflx%ri)) eflx%ri%uflx(iarr_all_mx) = eflx%ri%uflx(iarr_all_mx) / this%x(i2)
         else if (cdim == zdim) then
            if (associated(eflx%li)) eflx%li%uflx = eflx%li%uflx / this%x(i1)
            if (associated(eflx%li)) eflx%li%uflx(iarr_all_my) = eflx%li%uflx(iarr_all_my) / this%x(i1) ! that makes this%x(i1)**2
            if (associated(eflx%ri)) eflx%ri%uflx = eflx%ri%uflx / this%x(i1)
            if (associated(eflx%ri)) eflx%ri%uflx(iarr_all_my) = eflx%ri%uflx(iarr_all_my) / this%x(i1) ! that makes this%x(i1)**2
         endif
      endif

      if (this%coarsebnd(cdim, LO)%index(i1, i2) >= this%ijkse(cdim, LO)) then
         cpl%index = this%coarsebnd(cdim, LO)%index(i1, i2)
         if (.not. allocated(cpl%uflx)) call cpl%fpinit
         eflx%lo => cpl
         eflx%lo%index = eflx%lo%index - this%lhn(cdim, LO)
      else
         nullify(eflx%lo)
      endif
      if (this%coarsebnd(cdim, HI)%index(i1, i2) <= this%ijkse(cdim, HI)) then
         cpr%index = this%coarsebnd(cdim, HI)%index(i1, i2)
         if (.not. allocated(cpr%uflx)) call cpr%fpinit
         eflx%ro => cpr
         eflx%ro%index = eflx%ro%index - this%lhn(cdim, LO) + 1
      else
         nullify(eflx%ro)
      endif

   end subroutine set_fluxpointers

!> \brief Collect outgoing fine fluxes from 1D solver, do curvilinear scaling and store in appropriate array

   subroutine save_outfluxes(this, cdim, i1, i2, eflx)

      use constants, only: LO, HI
      use fluxtypes, only: ext_fluxes

      implicit none

      class(grid_container_bnd_t), intent(inout) :: this    !< object invoking type-bound procedure
      integer(kind=4),             intent(in)    :: cdim
      integer,                     intent(in)    :: i1
      integer,                     intent(in)    :: i2
      type(ext_fluxes),            intent(inout) :: eflx

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

      subroutine cyl_scale(fp)

         use constants,  only: ydim, zdim, GEO_RPZ
         use domain,     only: dom
         use fluidindex, only: iarr_all_mx, iarr_all_my
         use fluxtypes,  only: fluxpoint

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

end module grid_cont_bnd
