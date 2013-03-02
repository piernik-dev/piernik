! $Id$
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

!> \brief This module contains types useful for flux enforcing and exchange in (M)HD solvers

module fluxtypes

   implicit none

   private
   public  :: fluxpoint, ext_fluxes, fluxarray

   !> \brief Structure that contains u-flux at a single cell index.

   type :: fluxpoint
      real, dimension(:), allocatable :: uflx  !< u-flux
      integer                         :: index !< Index where the flux has to be applied
   contains
      procedure :: fpinit                      !< Allocate flux vector
      procedure :: fpcleanup                   !< Deallocate flux vector
   end type fluxpoint

   !> \brief Structure that may contain pointers to fluxes to be passed to or obtained from the rtvd routine. Unassociated pointer means that no aperation is required.

   type :: ext_fluxes
      type(fluxpoint), pointer :: li  !< incoming from the left (low index)
      type(fluxpoint), pointer :: lo  !< outgoing on the left
      type(fluxpoint), pointer :: ri  !< incoming from the right (high index)
      type(fluxpoint), pointer :: ro  !< outgoing on the right
   contains
      procedure :: init               !< Nullify point flux pointers
   end type ext_fluxes

   !>
   !! \brief Structure that contains u-flux at a single face of a grid container
   !!
   !! \todo Consider moving this to a separate file and s/fainit/nit/, s/facleanup/cleanup/
   !<

   type :: fluxarray
      real,    dimension(:,:,:), allocatable :: uflx  !< u-fluxes, shape (flind%all, n_b(dir1), n_b(dir2))
      integer, dimension(:,:),   allocatable :: index !< Index where the flux has to be applied, shape (n_b(dir1), n_b(dir2))
   contains
      procedure :: fainit                      !< Allocate flux array
      procedure :: facleanup                   !< Deallocate flux array
      procedure :: fa2fp                       !< Pick a point flux
      procedure :: fp2fa                       !< Store a point flux
   end type fluxarray

contains

!> \brief Allocate flux vector

   subroutine fpinit(this)

      use dataio_pub, only: die
      use fluidindex, only: flind

      implicit none

      class(fluxpoint), intent(inout) :: this

      if (allocated(this%uflx)) call die("[fluxtypes:fpinit] uflx already allocated")
      allocate(this%uflx(flind%all))

   end subroutine fpinit

!> \brief Deallocate flux vector

   subroutine fpcleanup(this)

      implicit none

      class(fluxpoint), intent(inout) :: this

      if (allocated(this%uflx)) deallocate(this%uflx)

   end subroutine fpcleanup

!> \brief Nullify point flux pointers

   subroutine init(this)

      implicit none

      class(ext_fluxes), intent(out) :: this

      this%li => null()
      this%lo => null()
      this%ri => null()
      this%ro => null()

   end subroutine init

!> \brief Allocate flux array

   subroutine fainit(this, i1, i2)

      use constants,  only: LO, HI
      use dataio_pub, only: die
      use fluidindex, only: flind

      implicit none

      class(fluxarray),          intent(inout) :: this
      integer, dimension(LO:HI), intent(in)    :: i1 !< 1st range
      integer, dimension(LO:HI), intent(in)    :: i2 !< 2nd range

      if (allocated(this%index) .or. allocated(this%uflx)) call die("[fluxtypes:fainit] already allocated")
      allocate(this%index(i1(LO):i1(HI), i2(LO):i2(HI)), this%uflx(flind%all, i1(LO):i1(HI), i2(LO):i2(HI)))

   end subroutine fainit

!> \brief Deallocate flux array

   subroutine facleanup(this)

      implicit none

      class(fluxarray), intent(inout) :: this

      if (allocated(this%index)) deallocate(this%index)
      if (allocated(this%uflx))  deallocate(this%uflx)

   end subroutine facleanup

!> \brief Pick a point flux

   function fa2fp(this, i1, i2) result(fp)

      implicit none

      class(fluxarray), intent(in) :: this
      integer,          intent(in) :: i1 !< 1st index
      integer,          intent(in) :: i2 !< 2nd index

      type(fluxpoint) :: fp

      fp%index = this%index(   i1, i2)
      fp%uflx  = this%uflx (:, i1, i2)

   end function fa2fp

!> \brief Store a point flux

   subroutine fp2fa(this, fp, i1, i2)

      use dataio_pub, only: die

      implicit none

      class(fluxarray), intent(inout) :: this
      type(fluxpoint),  intent(in)    :: fp !< point flux
      integer,          intent(in)    :: i1 !< 1st index
      integer,          intent(in)    :: i2 !< 2nd index

      if (this%index(i1, i2) /= fp%index) call die("[fluxtypes:fp2fa] inconsistent index")
      this%uflx (:, i1, i2) = fp%uflx

   end subroutine fp2fa

end module fluxtypes
