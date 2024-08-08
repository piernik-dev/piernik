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

!> \brief This module contains fluxarray type useful for flux enforcing and exchange in (M)HD solvers at f/c boundary.

module flx_arr

   implicit none

   private
   public  :: fluxarray

   !> \brief Structure that contains u- or b-flux at a single face of a grid container

   type :: fluxarray
      real,    dimension(:,:,:), allocatable :: uflx   !< u-fluxes, shape (flind%all, n_b(dir1), n_b(dir2))
      real,    dimension(:,:,:), allocatable :: bflx   !< b-fluxes, shape (psidim,    n_b(dir1), n_b(dir2)) (magnetic field components + psi)
      integer, dimension(:,:),   allocatable :: index  !< Index where the flux has to be applied, shape (n_b(dir1), n_b(dir2))
   contains
      procedure :: init     !< Allocate flux array
      procedure :: cleanup  !< Deallocate flux array
      procedure :: fa2fp    !< Pick a point flux
      procedure :: fp2fa    !< Store a point flux
   end type fluxarray

contains

!> \brief Allocate the flux array for a f/c face

   subroutine init(this, i1, i2)

      use constants,  only: LO, HI, psidim, has_B
      use dataio_pub, only: die
      use fluidindex, only: flind

      implicit none

      class(fluxarray),                  intent(inout) :: this  !< object invoking type bound procedure
      integer(kind=4), dimension(LO:HI), intent(in)    :: i1    !< 1st range
      integer(kind=4), dimension(LO:HI), intent(in)    :: i2    !< 2nd range

      if (allocated(this%index) .or. allocated(this%uflx)) call die("[flx_arr:init] already allocated")
      allocate(           this%index(           i1(LO):i1(HI), i2(LO):i2(HI)), &
           &              this%uflx (flind%all, i1(LO):i1(HI), i2(LO):i2(HI)))
      if (has_B) allocate(this%bflx (psidim,    i1(LO):i1(HI), i2(LO):i2(HI)))

   end subroutine init

!> \brief Deallocate flux arrays

   subroutine cleanup(this)

      implicit none

      class(fluxarray), intent(inout) :: this  !< object invoking type bound procedure

      if (allocated(this%index)) deallocate(this%index)
      if (allocated(this%uflx))  deallocate(this%uflx)
      if (allocated(this%bflx))  deallocate(this%bflx)

   end subroutine cleanup

!> \brief Pick a point flux

   subroutine fa2fp(this, fp, i1, i2)

      use constants,  only: has_B
      use dataio_pub, only: die
      use flx_cell,   only: fluxpoint

      implicit none

      class(fluxarray), intent(in)    :: this  !< object invoking type bound procedure
      integer,          intent(in)    :: i1    !< 1st index
      integer,          intent(in)    :: i2    !< 2nd index
      type(fluxpoint),  intent(inout) :: fp    !< buffer for the data

      fp%index = this%index(i1, i2)

      if (.not. allocated(this%uflx)) call die("[flx_arr:fa2fp] .not. allocated(this%uflx)")
      if (.not. allocated(fp%uflx)) call die("[flx_arr:fa2fp] .not. allocated(fp%uflx)")
      fp%uflx = this%uflx(:, i1, i2)

      if (has_B) then
         if (.not. allocated(this%bflx)) call die("[flx_arr:fa2fp] .not. allocated(this%bflx)")
         if (.not. allocated(fp%bflx)) call die("[flx_arr:fa2fp] .not. allocated(fp%bflx)")
         fp%bflx = this%bflx(:, i1, i2)
      endif

   end subroutine fa2fp

!> \brief Store a point flux

   subroutine fp2fa(this, fp, i1, i2)

      use constants,  only: has_B
      use dataio_pub, only: die
      use flx_cell,   only: fluxpoint

      implicit none

      class(fluxarray), intent(inout) :: this  !< object invoking type bound procedure
      type(fluxpoint),  intent(in)    :: fp    !< point flux
      integer,          intent(in)    :: i1    !< 1st index
      integer,          intent(in)    :: i2    !< 2nd index

      if (this%index(i1, i2) /= fp%index) call die("[flx_arr:fp2fa] inconsistent index")
      this%uflx(:, i1, i2) = fp%uflx
      if (has_B) this%bflx(:, i1, i2) = fp%bflx

   end subroutine fp2fa

end module flx_arr
