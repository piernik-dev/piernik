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
!! \brief Definitions of 3D and 4D data storage for named arrays
!!
!! \details The named arrays are used exclusively through grid container, where anyone may register a rank-3 (cg%q(:)) or rank-4 (cg%w(:)) array.
!!
!! The maintenance of named arrays(initialization, cleanup, boundary cell exchange, I/O) is unified as much as possible,
!! which saves us from writing a lot of partially duplicated code.
!<
module named_array

   implicit none

   private
   public :: named_array4d, named_array3d, p3, p4

   real, dimension(:,:,:),   pointer :: p3   !< auxiliary pointer to 3D named_arrays
   real, dimension(:,:,:,:), pointer :: p4   !< auxiliary pointer to 4D named_arrays

   !> \brief Common methods for 3D and 4D named arrays
   type, abstract :: generic_na
    contains
      procedure :: lb
      procedure :: ub
      procedure :: named_array_init
      procedure :: clean
      procedure :: check
      generic, public :: init => named_array_init
      !> \todo add also span and get_sweep
   end type generic_na

   !> \brief A named array for multi-scalar and vector fields
   type, extends(generic_na) :: named_array4d
      real, dimension(:,:,:,:), pointer :: arr => null()
      contains
         procedure :: array4d_associate
         procedure :: array4d_get_sweep
         procedure :: array4d_get_sweep_one_var
         procedure :: array4d_point
         procedure :: array4d_point_one_var
         procedure :: array4d_span
         procedure :: array4d_span_one_var
         procedure :: array4d_span_ijkse
         procedure :: array4d_span_ijkse8
         procedure :: array4d_span_one_var_ijkse
         generic, public :: init      => array4d_associate
         generic, public :: span      => array4d_span_one_var, array4d_span, array4d_span_one_var_ijkse, array4d_span_ijkse, array4d_span_ijkse8
         generic, public :: get_sweep => array4d_get_sweep_one_var, array4d_get_sweep
         generic, public :: point     => array4d_point, array4d_point_one_var
   end type named_array4d

   !> \brief A named array for scalar fields
   type, extends(generic_na) :: named_array3d
      real, dimension(:,:,:), pointer :: arr => null()
      contains
         procedure :: array3d_associate
         procedure :: array3d_span
         procedure :: array3d_span_ijkse
         procedure :: array3d_span_ijkse8
         generic, public :: init      => array3d_associate
         generic, public :: span      => array3d_span, array3d_span_ijkse, array3d_span_ijkse8
         procedure       :: get_sweep => array3d_get_sweep
         procedure       :: point     => array3d_point
   end type named_array3d

contains

!>
!! \brief Initialize a 4d named array
!!
!! \details The mbc component is initialized separately
!!
!! \warning Please note that maxloc and minloc return positions as all the declared lower bounds of array were 1, so whenever you plan to use these functions
!! on this%arr remember to add lbound(this%arr) - 1 to the result
!<
   subroutine named_array_init(this, n1, n2)

      use constants,  only: big_float, ndims, xdim, ydim, zdim, I_ONE
      use dataio_pub, only: die

      implicit none

      class(generic_na),             intent(inout) :: this
      integer(kind=4), dimension(:), intent(in)    :: n1
      integer(kind=4), dimension(:), intent(in)    :: n2

      if (size(n1) /= size(n2)) call die("[named_array:array_init] sizes differ")
      select type(this)
         type is (named_array3d)
            if (size(n1) /= ndims) call die("[named_array:array_init] expected 3d shape")
            if (.not.associated(this%arr)) allocate(this%arr(n1(xdim):n2(xdim), n1(ydim):n2(ydim), n1(zdim):n2(zdim)))
            this%arr = big_float
         type is (named_array4d)
            if (size(n1) /= I_ONE + ndims) call die("[named_array:array_init] expected 4d shape")
            if (.not.associated(this%arr)) allocate(this%arr(n1(I_ONE):n2(I_ONE), n1(I_ONE+xdim):n2(I_ONE+xdim), n1(I_ONE+ydim):n2(I_ONE+ydim), n1(I_ONE+zdim):n2(I_ONE+zdim)))
            this%arr = big_float
         class default
            call die("[named_array:named_array_init] No initialization for generic named array")
      end select

   end subroutine named_array_init

!> \brief deallocate named array

   subroutine clean(this)

      use dataio_pub, only: die

      implicit none

      class(generic_na), intent(inout) :: this

      select type(this)
         type is (named_array3d)
            if (associated(this%arr)) deallocate(this%arr)
         type is (named_array4d)
            if (associated(this%arr)) deallocate(this%arr)
         class default
            call die("[named_array:clean] No cleanup for generic named array")
      end select
   end subroutine clean

!> \brief check if the array was initialized with sane values

   logical function check(this)

      use constants,  only: big_float
      use dataio_pub, only: warn, die

      implicit none

      class(generic_na), intent(inout) :: this                  !! \todo i want to become polymorphic class(*) :/

      check = .false.
      select type(this)
         type is (named_array3d)
            if (associated(this%arr)) then
               check = any( this%arr >= big_float )
            else
               call warn("[named_array:check] Array not allocated!")
            endif
         type is (named_array4d)
            if (associated(this%arr)) then
               check = any( this%arr >= big_float )
            else
               call warn("[named_array:check] Array not allocated!")
            endif
         class default
            call die("[named_array:ccheck] No check for generic named array")
      end select

   end function check

!> \brief Get the upper bound of a named array

   function ub(this, dim_) result(n)

      use constants,  only: INVALID
      use dataio_pub, only: die

      implicit none

      class(generic_na), intent(in) :: this
      integer(kind=4),      intent(in) :: dim_
      integer(kind=4) :: n

      n = INVALID
      select type(this)
         type is (named_array3d)
            n = ubound(this%arr, dim=dim_, kind=4)
         type is (named_array4d)
            n = ubound(this%arr, dim=dim_, kind=4)
         class default
            call die("[named_array:ub] No upper bound for generic named array")
      end select

   end function ub

!> \brief Get the lower bound of a named array

   function lb(this, dim_) result(n)

      use constants,  only: INVALID
      use dataio_pub, only: die

      implicit none

      class(generic_na), intent(in) :: this
      integer(kind=4),      intent(in) :: dim_
      integer(kind=4) :: n

      n = INVALID
      select type(this)
         type is (named_array3d)
            n = lbound(this%arr,dim=dim_,kind=4)
         type is (named_array4d)
            n = lbound(this%arr,dim=dim_,kind=4)
         class default
            call die("[named_array:ub] No lower bound for generic named array")
      end select

   end function lb

!> \brief Initialize named array with a predefined simple array

   subroutine array3d_associate(this,other)

      implicit none

      class(named_array3d),           intent(inout) :: this
      real, dimension(:,:,:), target, intent(in)    :: other

      if (.not.associated(this%arr)) this%arr => other

   end subroutine array3d_associate

!> \brief Initialize named array with a predefined simple array

   subroutine array4d_associate(this,other)

      implicit none

      class(named_array4d),             intent(inout) :: this
      real, dimension(:,:,:,:), target, intent(in)    :: other

      if (.not.associated(this%arr)) this%arr => other

   end subroutine array4d_associate

!> \brief Get a selected line of values

   function array3d_get_sweep(this,ndim,i1,i2) result(p1d)

      use constants, only: xdim, ydim, zdim

      implicit none

      class(named_array3d), intent(inout) :: this
      integer(kind=4),      intent(in)    :: ndim
      integer,              intent(in)    :: i1, i2

      real, dimension(:), pointer         :: p1d

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

!> \brief Get a selected line of values of one variable

   function array4d_get_sweep_one_var(this,ndim,nn,i1,i2) result(p1d)

      use constants, only: xdim, ydim, zdim

      implicit none

      class(named_array4d), intent(inout) :: this
      integer(kind=4),      intent(in)    :: ndim, nn
      integer,              intent(in)    :: i1, i2

      real, dimension(:),  pointer        :: p1d

      select case (ndim)
         case (xdim)
            p1d => this%arr(nn,:,i1,i2)
         case (ydim)
            p1d => this%arr(nn,i2,:,i1)
         case (zdim)
            p1d => this%arr(nn,i1,i2,:)
      end select

   end function array4d_get_sweep_one_var

!> \brief Get a selected line of values of all variables

   function array4d_get_sweep(this,ndim,i1,i2) result(p1d)

      use constants, only: xdim, ydim, zdim

      implicit none

      class(named_array4d), intent(inout) :: this
      integer(kind=4),      intent(in)    :: ndim
      integer,              intent(in)    :: i1, i2

      real, dimension(:,:), pointer       :: p1d

      select case (ndim)
         case (xdim)
            p1d => this%arr(:,:,i1,i2)
         case (ydim)
            p1d => this%arr(:,i2,:,i1)
         case (zdim)
            p1d => this%arr(:,i1,i2,:)
      end select

   end function array4d_get_sweep

!> \brief Get a selected value from the rank-3 array

   function array3d_point(this, v) result(p)

      use constants,  only: xdim, ydim, zdim
      use dataio_pub, only: die

      implicit none

      class(named_array3d),  intent(inout) :: this
      integer, dimension(:), intent(in)    :: v

      real                                     :: p

      if (associated(this%arr)) then
         p = this%arr(v(xdim),v(ydim),v(zdim))
      else
         call die("[named_array:array3d_point] this%arr not associated")
         p = -huge(1.) ! suppress use of uninitialized variable warning
      endif

   end function array3d_point

!> \brief Get a selected value from the rank-4 array

   function array4d_point(this, v) result(p1d)

      use constants,  only: xdim, ydim, zdim
      use dataio_pub, only: die

      implicit none

      class(named_array4d),          intent(inout) :: this
      integer(kind=4), dimension(:), intent(in)    :: v

      real,    dimension(:),     pointer       :: p1d

      if (associated(this%arr)) then
         p1d => this%arr(:,v(xdim),v(ydim),v(zdim))
      else
         call die("[named_array:array4d_point] this%arr not associated")
         p1d => null() ! suppress use of uninitialized variable warning
      endif

   end function array4d_point

!> \brief Get a selected vector from the rank-4 array

   function array4d_point_one_var(this, nn, v) result(p)

      use constants,  only: xdim, ydim, zdim
      use dataio_pub, only: die

      implicit none

      class(named_array4d),          intent(inout) :: this
      integer(kind=4),               intent(in)    :: nn
      integer(kind=4), dimension(:), intent(in)    :: v

      real                                     :: p

      if (associated(this%arr)) then
         p = this%arr(nn,v(xdim),v(ydim),v(zdim))
      else
         call die("[named_array:array4d_point_one_var] this%arr not associated")
         p = -huge(1.) ! suppress use of uninitialized variable warning
      endif

   end function array4d_point_one_var

!> \brief Get a selected range of values

   function array3d_span(this,v1,v2) result(p3d)

      use constants, only: xdim, ydim, zdim

      implicit none

      class(named_array3d),          intent(inout) :: this
      integer(kind=4), dimension(:), intent(in)    :: v1, v2

      real,    dimension(:,:,:), pointer       :: p3d

      if (.not.associated(this%arr)) then
         p3d => null()
      else
         p3d => this%arr(v1(xdim):v2(xdim),v1(ydim):v2(ydim),v1(zdim):v2(zdim))
      endif

   end function array3d_span

!> \brief Get a selected range of values

   function array3d_span_ijkse(this,v) result(p3d)

      use constants, only: xdim, ydim, zdim, LO, HI

      implicit none

      class(named_array3d),            intent(inout) :: this
      integer(kind=4), dimension(:,:), intent(in)    :: v

      real,    dimension(:,:,:), pointer             :: p3d

      if (.not.associated(this%arr)) then
         p3d => null()
      else
         p3d => this%arr(v(xdim,LO):v(xdim,HI),v(ydim,LO):v(ydim,HI),v(zdim,LO):v(zdim,HI))
      endif

   end function array3d_span_ijkse

!> \brief Get a selected range of values

   function array3d_span_ijkse8(this,v) result(p3d)

      use constants, only: xdim, ydim, zdim, LO, HI

      implicit none

      class(named_array3d),            intent(inout) :: this
      integer(kind=8), dimension(:,:), intent(in)    :: v

      real,    dimension(:,:,:), pointer             :: p3d

      if (.not.associated(this%arr)) then
         p3d => null()
      else
         p3d => this%arr(v(xdim,LO):v(xdim,HI),v(ydim,LO):v(ydim,HI),v(zdim,LO):v(zdim,HI))
      endif

   end function array3d_span_ijkse8

!> \brief Get a selected range of values of one variable

   function array4d_span_one_var(this,nn,v1,v2) result(p3d)

      use constants, only: xdim, ydim, zdim

      implicit none

      class(named_array4d),          intent(inout) :: this
      integer(kind=4),               intent(in)    :: nn
      integer(kind=4), dimension(:), intent(in)    :: v1, v2

      real,    dimension(:,:,:),  pointer      :: p3d

      if (.not.associated(this%arr)) then
         p3d => null()
      else
         p3d => this%arr(nn,v1(xdim):v2(xdim),v1(ydim):v2(ydim),v1(zdim):v2(zdim))
      endif

   end function array4d_span_one_var

!> \brief Get a selected range of values of one variable

   function array4d_span_one_var_ijkse(this,nn,v) result(p3d)

      use constants, only: xdim, ydim, zdim, LO, HI

      implicit none

      class(named_array4d),             intent(inout) :: this
      integer(kind=4),                  intent(in)    :: nn
      integer(kind=4), dimension(:, :), intent(in)    :: v

      real,    dimension(:,:,:),  pointer            :: p3d

      if (.not.associated(this%arr)) then
         p3d => null()
      else
         p3d => this%arr(nn,v(xdim,LO):v(xdim,HI),v(ydim,LO):v(ydim,HI),v(zdim,LO):v(zdim,HI))
      endif

   end function array4d_span_one_var_ijkse

!> \brief Get a selected line of values of all variables

   function array4d_span(this,v1,v2) result(p3d)

      use constants, only: xdim, ydim, zdim

      implicit none

      class(named_array4d),          intent(inout) :: this
      integer(kind=4), dimension(:), intent(in)    :: v1, v2

      real,    dimension(:,:,:,:), pointer     :: p3d

      if (.not.associated(this%arr)) then
         p3d => null()
      else
         p3d => this%arr(:,v1(xdim):v2(xdim),v1(ydim):v2(ydim),v1(zdim):v2(zdim))
      endif

   end function array4d_span

!> \brief Get a selected line of values of all variables

   function array4d_span_ijkse(this,v) result(p3d)

      use constants, only: xdim, ydim, zdim, LO, HI

      implicit none

      class(named_array4d),             intent(inout) :: this
      integer(kind=4), dimension(:, :), intent(in)    :: v

      real,    dimension(:,:,:,:), pointer           :: p3d

      if (.not.associated(this%arr)) then
         p3d => null()
      else
         p3d => this%arr(:,v(xdim,LO):v(xdim,HI),v(ydim,LO):v(ydim,HI),v(zdim,LO):v(zdim,HI))
      endif

   end function array4d_span_ijkse

   function array4d_span_ijkse8(this,v) result(p3d)

      use constants,  only: xdim, ydim, zdim, LO, HI
      use dataio_pub, only: msg, die

      implicit none

      class(named_array4d),             intent(inout) :: this
      integer(kind=8), dimension(:, :), intent(in)    :: v

      real,    dimension(:,:,:,:), pointer           :: p3d

      integer :: d

      if (.not.associated(this%arr)) then
         p3d => null()
      else
         do d = xdim, zdim
            if (v(d, LO) < lbound(this%arr, dim=1+d) .or. v(d, HI) > ubound(this%arr, dim=1+d)) then
               write(msg,'(2(a,4i6),2(a,3i6),a)')"[named_array:array4d_span_ijkse8] this%arr[",lbound(this%arr),"]..[",ubound(this%arr),"], not embedding v[",v(:,LO),"]..[",v(:,HI),"]"
               call die(msg) ! we want backtrace
            endif
         enddo
         p3d => this%arr(:,v(xdim,LO):v(xdim,HI),v(ydim,LO):v(ydim,HI),v(zdim,LO):v(zdim,HI))
      endif

   end function array4d_span_ijkse8

end module named_array
