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
!>
!! \brief Definitions of compound types and subroutine templates for user hooks
!<
module types

   use constants, only: ndims

   implicit none

   private
   public :: axes, tsl_container, value, array4d, array3d

   type :: array4d
      real, dimension(:,:,:,:), pointer :: arr => null()
      contains
         procedure :: array4d_init           ! \todo check why private here does not  work as expected
         procedure :: array4d_associate
         procedure :: clean => array4d_clean
         procedure :: array4d_get_sweep
         procedure :: array4d_get_sweep_one_var
         procedure :: check => array4d_check_if_dirty
         procedure :: lb => array4d_lbound
         procedure :: ub => array4d_ubound
         generic, public :: init => array4d_init, array4d_associate
         generic, public :: get_sweep => array4d_get_sweep_one_var, array4d_get_sweep
   end type array4d

   type :: array3d
      real, dimension(:,:,:), pointer :: arr => null()
      contains
         procedure :: array3d_init
         procedure :: array3d_associate
         procedure :: clean => array3d_clean
         procedure :: check => array3d_check_if_dirty
         procedure :: get_sweep => array3d_get_sweep
         procedure :: lb => array3d_lbound
         procedure :: ub => array3d_ubound
         generic, public :: init => array3d_init, array3d_associate
   end type array3d

   type :: value
      real                      :: val
      real,    dimension(ndims) :: coords
      integer, dimension(ndims) :: loc
      integer                   :: proc
   end type value

   type :: axes
      real, allocatable, dimension(:) :: x      !< array of x-positions of %grid cells centers
      real, allocatable, dimension(:) :: y      !< array of y-positions of %grid cells centers
      real, allocatable, dimension(:) :: z      !< array of z-positions of %grid cells centers
   end type axes

   type :: tsl_container
      logical :: dummy
#ifdef COSM_RAYS
      real :: encr_min, encr_max
#endif /* COSM_RAYS */

#ifdef RESISTIVE
      real :: etamax
#endif /* RESISTIVE */

#ifdef MAGNETIC
      real :: b_min, b_max, divb_max
#endif /* MAGNETIC */

#ifdef IONIZED
      real :: vai_max
#endif /* IONIZED */

#ifdef VARIABLE_GP
      real :: gpxmax, gpymax, gpzmax
#endif /* VARIABLE_GP */

   end type tsl_container

contains

   subroutine array3d_init(this, n3)

      use constants, only: big_float

      implicit none

      class(array3d), intent(inout) :: this
      integer(kind=4), dimension(3), intent(in) :: n3

      if (.not.associated(this%arr)) allocate(this%arr(n3(1), n3(2), n3(3)))
      this%arr = big_float
!     if (.not.associated(this%arr)) this%arr = reshape( [ ( big_float, i=1, product(n3(:)) ) ], [ n3(1), n3(2), n3(3) ] ) ! lhs realloc

   end subroutine array3d_init

   subroutine array3d_clean(this)
      implicit none
      class(array3d), intent(inout) :: this

      if (associated(this%arr)) deallocate(this%arr)
      return
   end subroutine array3d_clean

   logical function array3d_check_if_dirty(this)
      use constants, only: big_float
      implicit none
      class(array3d), intent(inout) :: this                  !! \todo i want to become polymorphic class(*) :/

      array3d_check_if_dirty = any( this%arr >= big_float )

   end function array3d_check_if_dirty

   subroutine array3d_associate(this,other)
      implicit none
      class(array3d), intent(inout) :: this
      real, allocatable, dimension(:,:,:), target :: other

      if (.not.associated(this%arr)) this%arr => other

   end subroutine array3d_associate

   subroutine array4d_init(this, n4)

      use constants, only: big_float

      implicit none

      class(array4d), intent(inout) :: this
      integer(kind=4), dimension(4), intent(in) :: n4

      if (.not.associated(this%arr)) allocate(this%arr(n4(1), n4(2), n4(3), n4(4)))
      this%arr = big_float
!     if (.not.associated(this%arr)) this%arr = reshape( [ ( big_float, i=1, product(n4(:)) ) ], [ n4(1), n4(2), n4(3), n4(4) ] ) ! lhs realloc

   end subroutine array4d_init

   subroutine array4d_associate(this,other)
      implicit none
      class(array4d), intent(inout) :: this
      real, allocatable, dimension(:,:,:,:), target :: other

      if (.not.associated(this%arr)) this%arr => other

   end subroutine array4d_associate

   subroutine array4d_clean(this)
      implicit none
      class(array4d), intent(inout) :: this                  !! Unlimited polymorphism at (1) not yet supported

      if (associated(this%arr)) deallocate(this%arr)

   end subroutine array4d_clean

   logical function array4d_check_if_dirty(this)
      use constants, only: big_float
      implicit none
      class(array4d), intent(inout) :: this                  !! \todo i want to become polymorphic class(*) when I grow older

      array4d_check_if_dirty = any( this%arr >= big_float )

   end function array4d_check_if_dirty

   function array3d_get_sweep(this,ndim,i1,i2) result(p1d)
      use constants, only: xdim, ydim, zdim
      implicit none
      class(array3d), intent(inout) :: this
      real, dimension(:),  pointer  :: p1d
      integer(kind=4), intent(in)   :: ndim
      integer, intent(in)           :: i1, i2

      if (.not.associated(this%arr)) then
         p1d => null()
      else
         select case (ndim)
            case (xdim)
               p1d => this%arr(:,i1,i2)
            case (ydim)
               p1d => this%arr(i2,:,i1)
            case (zdim)
               p1d => this%arr(i1,i2,:)
         end select
      endif
   end function array3d_get_sweep

   function array4d_get_sweep_one_var(this,ndim,nn,i1,i2) result(p1d)
      use constants, only: xdim, ydim, zdim
      implicit none
      class(array4d), intent(inout)      :: this

      integer(kind=4), intent(in)        :: ndim, nn
      integer, intent(in)                :: i1, i2

      real, dimension(:),  pointer       :: p1d

      select case (ndim)
         case (xdim)
            p1d => this%arr(nn,:,i1,i2)
         case (ydim)
            p1d => this%arr(nn,i2,:,i1)
         case (zdim)
            p1d => this%arr(nn,i1,i2,:)
      end select
   end function array4d_get_sweep_one_var

   function array4d_get_sweep(this,ndim,i1,i2) result(p1d)

      use constants, only: xdim, ydim, zdim

      implicit none

      class(array4d), intent(inout)      :: this
      integer(kind=4), intent(in)        :: ndim
      integer, intent(in)                :: i1, i2

      real, dimension(:,:),  pointer     :: p1d

      select case (ndim)
         case (xdim)
            p1d => this%arr(:,:,i1,i2)
         case (ydim)
            p1d => this%arr(:,i2,:,i1)
         case (zdim)
            p1d => this%arr(:,i1,i2,:)
      end select
   end function array4d_get_sweep

   function array4d_ubound(this,dim_) result(n)
      implicit none
      class(array4d), intent(in) :: this
      integer(kind=4), intent(in) :: dim_
      integer(kind=8) :: n

      n = ubound(this%arr, dim=dim_, kind=8)
   end function array4d_ubound

   function array4d_lbound(this,dim_) result(n)

      implicit none
      class(array4d), intent(in) :: this
      integer(kind=4), intent(in) :: dim_
      integer(kind=8) :: n

      n = lbound(this%arr,dim=dim_,kind=8)
   end function array4d_lbound

   function array3d_ubound(this,dim_) result(n)
      implicit none
      class(array3d), intent(in) :: this
      integer(kind=4), intent(in) :: dim_
      integer(kind=8) :: n

      n = ubound(this%arr, dim=dim_, kind=8)
   end function array3d_ubound

   function array3d_lbound(this,dim_) result(n)

      implicit none
      class(array3d), intent(in) :: this
      integer(kind=4), intent(in) :: dim_
      integer(kind=8) :: n

      n = lbound(this%arr,dim=dim_,kind=8)
   end function array3d_lbound

end module types
