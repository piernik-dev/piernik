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
!! \brief Definitions of 3D and 4D named arrays
!!
!! \details The named arrays are used almost exclusively in grid container, where anyone may register a rank-3 (cg%q(:)) or rank-4 (cg%w) array.
!! The mainteance of named arrays(initialization, cleanup, boundary cell exchange, I/O) is unified as much as possible,
!! which saves us from writing a lot of partially duplicated code.
!<
module named_array

   use constants, only: dsetnamelen

   implicit none

   private
   public :: named_array4d, named_array3d, mbc_list, p3, p4, na_var, qna, wna

   !> \brief Common properties of 3D and 4D named arrays
   type :: na_var
      character(len=dsetnamelen)                 :: name          !< a user-provided id for the array
      logical                                    :: vital         !< fields that are subject of automatic prolongation and restriction (e.g. state variables)
      integer(kind=4)                            :: restart_mode  !< AT_IGNORE: do not write to restart, AT_NO_B write without ext. boundaries, AT_OUT_B write with ext. boundaries
      integer(kind=4)                            :: ord_prolong   !< Prolongation order for the variable
      integer(kind=4), allocatable, dimension(:) :: position      !< VAR_CENTER by default, also possible VAR_CORNER and VAR_[XYZ]FACE
      integer(kind=4)                            :: dim4          !< <=0 for 3D arrays, >0 for 4D arrays
      logical                                    :: multigrid     !< .true. for variables that may exist below base level (e.g. work fields for multigrid solver)
   end type na_var

   type :: na_var_list
      type(na_var), dimension(:), allocatable :: lst
    contains
      procedure :: ind                                           !< Get the index of a named array of given name.
      procedure :: exists                                        !< Check if a named array of given name is already registered
      procedure :: print_vars                                    !< Write a summary on registered fields. Can be useful for debugging
   end type na_var_list

   type(na_var_list) :: qna !< list of registered 3D named arrays
   type(na_var_list) :: wna !< list of registered 4D named arrays

   real, dimension(:,:,:),   pointer :: p3   !< auxiliary pointer to 3D named_arrays
   real, dimension(:,:,:,:), pointer :: p4   !< auxiliary pointer to 4D named_arrays

   !> \brief Common methods for 3D and 4D named arrays
   type, abstract :: generic_na
    contains
      procedure(g_na_clean), deferred, pass(this) :: clean
      procedure(g_na_check), deferred, pass(this) :: check
      procedure(g_na_b),     deferred, pass(this) :: lb
      procedure(g_na_b),     deferred, pass(this) :: ub
      !> \todo add also init, span and get_sweep
   end type generic_na

   interface
      subroutine g_na_clean(this)
         import generic_na
         implicit none
         class(generic_na), intent(inout) :: this
      end subroutine g_na_clean

      logical function g_na_check(this)
         import generic_na
         implicit none
         class(generic_na), intent(inout) :: this
      end function g_na_check

      function g_na_b(this,dim_) result(n)
         import generic_na
         implicit none
         class(generic_na), intent(in) :: this
         integer(kind=4), intent(in) :: dim_
         integer(kind=4) :: n
      end function g_na_b
   end interface

   !> \brief List of MPI Boundary conditions Containers for boundary exchanges
   type :: mbc_list
      integer(kind=4), dimension(:), allocatable :: mbc  !< MPI Boundary conditions Container for each segment
   end type mbc_list

   !> \brief A named array for multi-scalar and vector fields
   type, extends(generic_na) :: named_array4d
      real, dimension(:,:,:,:), pointer :: arr => null()
      type(mbc_list), dimension(:,:), allocatable :: w_i_mbc  !< MPI Boundary conditions Containers for incoming guardcell updates on the w arrays
      type(mbc_list), dimension(:,:), allocatable :: w_o_mbc  !< MPI Boundary conditions Containers for outgoing guardcell updates on the w arrays
      contains
         procedure :: array4d_init
         procedure :: array4d_associate
         procedure :: clean => array4d_clean
         procedure :: array4d_get_sweep
         procedure :: array4d_get_sweep_one_var
         procedure :: array4d_point
         procedure :: array4d_point_one_var
         procedure :: array4d_span
         procedure :: array4d_span_one_var
         procedure :: array4d_span_ijkse
         procedure :: array4d_span_one_var_ijkse
         procedure :: check => array4d_check_if_dirty
         procedure :: lb => array4d_lbound
         procedure :: ub => array4d_ubound
         generic, public :: init => array4d_init, array4d_associate
         generic, public :: get_sweep => array4d_get_sweep_one_var, array4d_get_sweep
         generic, public :: point => array4d_point, array4d_point_one_var
         generic, public :: span => array4d_span_one_var, array4d_span, array4d_span_one_var_ijkse, array4d_span_ijkse
   end type named_array4d

   !> \brief A named array for scalar fields
   type, extends(generic_na) :: named_array3d
      real, dimension(:,:,:), pointer :: arr => null()
      contains
         procedure :: array3d_init
         procedure :: array3d_associate
         procedure :: clean => array3d_clean
         procedure :: check => array3d_check_if_dirty
         procedure :: get_sweep => array3d_get_sweep
         procedure :: point => array3d_point
         procedure :: array3d_span
         procedure :: array3d_span_ijkse
         procedure :: array3d_span_ijkse8
         procedure :: lb => array3d_lbound
         procedure :: ub => array3d_ubound
         generic, public :: init => array3d_init, array3d_associate
         generic, public :: span => array3d_span, array3d_span_ijkse, array3d_span_ijkse8
   end type named_array3d

contains

!>
!! \brief Get the index of a named array of given name.
!!
!! \warning OPT The indices aren't updated so cache them, whenever possible
!<
   function ind(this, name) result(rind)

      use dataio_pub,  only: die, msg, warn

      implicit none

      class(na_var_list), intent(inout) :: this
      character(len=*),   intent(in)    :: name

      integer :: rind, i

      rind = 0

      do i = lbound(this%lst, dim=1, kind=4), ubound(this%lst, dim=1, kind=4)
         if (trim(name) ==  this%lst(i)%name) then
            if (rind /= 0) then
               write(msg, '(2a)') "[named_array:ind] multiple entries with the same name: ", trim(name)
               call die(msg)
            endif
            rind = i
         endif
      enddo

      if (rind == 0) then
         write(msg, '(2a)') "[named_array:ind] requested entry not found: ", trim(name)
         call warn(msg)
      endif

   end function ind

!> \brief Check if a named array of given name is already registered

   function exists(this, name)

      use dataio_pub,  only: die, msg

      implicit none

      class(na_var_list), intent(inout) :: this
      character(len=*),   intent(in)    :: name

      logical :: exists
      integer :: i

      exists = .false.

      if (allocated(this%lst)) then
         do i = lbound(this%lst(:), dim=1), ubound(this%lst(:), dim=1)
            if (trim(name) ==  this%lst(i)%name) then
               if (exists) then
                  write(msg, '(2a)') "[named_array:exists] multiple entries with the same name: ", trim(name)
                  call die(msg)
               endif
               exists = .true.
            endif
         enddo
      endif

   end function exists

!> \brief Summarize all registered fields and their properties

   subroutine print_vars(this)

      use constants,  only: INVALID
      use dataio_pub, only: printinfo, warn, msg
      use mpisetup,   only: slave

      implicit none

      class(na_var_list), intent(inout) :: this

      integer :: i, d3

      if (slave) return

      d3 = count(this%lst(:)%dim4 == INVALID)

      if (d3 /= 0) then
         write(msg,'(a,i2,a)')"[named_array:print_vars] Found ",size(this%lst(:))," rank-3 arrays:"
         call printinfo(msg)
      endif
      if (count(this%lst(:)%dim4 /= INVALID) /= 0) then
         write(msg,'(a,i2,a)')"[named_array:print_vars] Found ",size(this%lst(:))," rank-4 arrays:"
         call printinfo(msg)
         if (d3 /=0) call warn("[named_array:print_vars] Both rank-3 and rank-4 named arrays are present in the same list!")
      endif

      do i = lbound(this%lst(:), dim=1), ubound(this%lst(:), dim=1)
         if (this%lst(i)%dim4 == INVALID) then
            write(msg,'(3a,l2,a,i2,a,l2,2(a,i2))')"'", this%lst(i)%name, "', vital=", this%lst(i)%vital, ", restart_mode=", this%lst(i)%restart_mode, &
                 &                                ", multigrid=", this%lst(i)%multigrid, ", ord_prolong=", this%lst(i)%ord_prolong, ", position=", this%lst(i)%position(:)
         else
            write(msg,'(3a,l2,a,i2,a,l2,2(a,i2),a,100i2)')"'", this%lst(i)%name, "', vital=", this%lst(i)%vital, ", restart_mode=", this%lst(i)%restart_mode, &
                 &                                        ", multigrid=", this%lst(i)%multigrid, ", ord_prolong=", this%lst(i)%ord_prolong, &
                 &                                        ", components=", this%lst(i)%dim4, ", position=", this%lst(i)%position(:)
         endif
         call printinfo(msg)
      enddo

   end subroutine print_vars

!>
!! \brief Initialize a 3d named array
!!
!! \details The mbc component is initialized separately. Note that mbc is common for all 3d named arrays and is a member of the grid container type.
!<

   subroutine array3d_init(this, n)

      use constants,  only: big_float, ndims, xdim, ydim, zdim
      use dataio_pub, only: die

      implicit none

      class(named_array3d),          intent(inout) :: this
      integer(kind=4), dimension(:), intent(in)    :: n

      if (size(n) /= ndims) call die("[named_array:array_init] expected 3d shape")
      if (.not.associated(this%arr)) allocate(this%arr(n(xdim), n(ydim), n(zdim)))
      this%arr = big_float
      ! if (.not.associated(this%arr)) this%arr = reshape( [ ( big_float, i=1, product(n(:)) ) ], [ n(1), n(2), n(3) ] ) ! lhs realloc

   end subroutine array3d_init

!>
!! \brief Initialize a 4d named array
!!
!! \details The mbc component is initialized separately
!<

   subroutine array4d_init(this, n)

      use constants,  only: big_float, ndims, xdim, ydim, zdim, I_ONE
      use dataio_pub, only: die

      implicit none

      class(named_array4d),          intent(inout) :: this
      integer(kind=4), dimension(:), intent(in)    :: n

      if (size(n) /= I_ONE + ndims) call die("[named_array:array_init] expected 4d shape")
      if (.not.associated(this%arr)) allocate(this%arr(n(I_ONE), n(I_ONE+xdim), n(I_ONE+ydim), n(I_ONE+zdim)))
      this%arr = big_float
      ! if (.not.associated(this%arr)) this%arr = reshape( [ ( big_float, i=1, product(n(:)) ) ], [ n(1), n(2), n(3), n(4) ] ) ! lhs realloc

   end subroutine array4d_init

!>
!! \brief deallocate array
!!
!! \details The mbc component is a member of the grid container type and this is cleaned up elsewhere
!<

   subroutine array3d_clean(this)

      implicit none

      class(named_array3d), intent(inout) :: this

      if (associated(this%arr)) deallocate(this%arr)

   end subroutine array3d_clean

!> \brief check if the array was initialized with sane values

   logical function array3d_check_if_dirty(this)

      use constants,  only: big_float
      use dataio_pub, only: warn

      implicit none

      class(named_array3d), intent(inout) :: this                  !! \todo i want to become polymorphic class(*) :/

      if (associated(this%arr)) then
         array3d_check_if_dirty = any( this%arr >= big_float )
      else
         call warn("[named_array:array3d_check_if_dirty] Array not allocated!")
         array3d_check_if_dirty = .false.
      endif

   end function array3d_check_if_dirty

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

!> \brief deallocate array and mbc

   subroutine array4d_clean(this)

      use constants, only: xdim, zdim, INVALID
      use mpisetup,  only: mpi_err

      implicit none

      class(named_array4d), intent(inout) :: this                  !! Unlimited polymorphism at (1) not yet supported

      integer :: d, b, g

      if (associated(this%arr)) deallocate(this%arr)

      if (allocated(this%w_i_mbc) .and. allocated(this%w_o_mbc)) then
         do d = xdim, zdim
            do b = lbound(this%w_o_mbc, dim=2), ubound(this%w_o_mbc, dim=2)
               if (allocated(this%w_i_mbc(d, b)%mbc)) then
                  do g = lbound(this%w_i_mbc(d, b)%mbc, dim=1), ubound(this%w_i_mbc(d, b)%mbc, dim=1)
                     if (this%w_i_mbc(d, b)%mbc(g) /= INVALID) call MPI_Type_free(this%w_i_mbc(d, b)%mbc(g), mpi_err)
                  enddo
                  deallocate(this%w_i_mbc(d, b)%mbc)
               endif
               if (allocated(this%w_o_mbc(d, b)%mbc)) then
                  do g = lbound(this%w_o_mbc(d, b)%mbc, dim=1), ubound(this%w_o_mbc(d, b)%mbc, dim=1)
                     if (this%w_o_mbc(d, b)%mbc(g) /= INVALID) call MPI_Type_free(this%w_o_mbc(d, b)%mbc(g), mpi_err)
                  enddo
                  deallocate(this%w_o_mbc(d, b)%mbc)
               endif
            enddo
         enddo
         deallocate(this%w_i_mbc, this%w_o_mbc)
      endif

   end subroutine array4d_clean

!> \brief check if the array was initialized with sane values

   logical function array4d_check_if_dirty(this)

      use constants,  only: big_float
      use dataio_pub, only: warn

      implicit none

      class(named_array4d), intent(inout) :: this                  !> \todo i want to become polymorphic class(*) when I grow older

      if (associated(this%arr)) then
         array4d_check_if_dirty = any( this%arr >= big_float )
      else
         call warn("[named_array:array4d_check_if_dirty] Array not allocated!")
         array4d_check_if_dirty = .false.
      endif

   end function array4d_check_if_dirty

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

      use constants,  only: xdim, ydim, zdim, ndims
      use dataio_pub, only: die

      implicit none

      class(named_array3d),      intent(inout) :: this
      integer, dimension(ndims), intent(in)    :: v

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

      use constants,  only: xdim, ydim, zdim, ndims
      use dataio_pub, only: die

      implicit none

      class(named_array4d),              intent(inout) :: this
      integer(kind=4), dimension(ndims), intent(in)    :: v

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

      use constants,  only: xdim, ydim, zdim, ndims
      use dataio_pub, only: die

      implicit none

      class(named_array4d),              intent(inout) :: this
      integer(kind=4),                   intent(in)    :: nn
      integer(kind=4), dimension(ndims), intent(in)    :: v

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

      use constants, only: xdim, ydim, zdim, ndims

      implicit none

      class(named_array3d),              intent(inout) :: this
      integer(kind=4), dimension(ndims), intent(in)    :: v1, v2

      real,    dimension(:,:,:), pointer       :: p3d

      if (.not.associated(this%arr)) then
         p3d => null()
      else
         p3d => this%arr(v1(xdim):v2(xdim),v1(ydim):v2(ydim),v1(zdim):v2(zdim))
      endif

   end function array3d_span

!> \brief Get a selected range of values

   function array3d_span_ijkse(this,v) result(p3d)

      use constants, only: xdim, ydim, zdim, ndims, LO, HI

      implicit none

      class(named_array3d),                    intent(inout) :: this
      integer(kind=4), dimension(ndims,LO:HI), intent(in)    :: v

      real,    dimension(:,:,:), pointer             :: p3d

      if (.not.associated(this%arr)) then
         p3d => null()
      else
         p3d => this%arr(v(xdim,LO):v(xdim,HI),v(ydim,LO):v(ydim,HI),v(zdim,LO):v(zdim,HI))
      endif

   end function array3d_span_ijkse

!> \brief Get a selected range of values

   function array3d_span_ijkse8(this,v) result(p3d)

      use constants, only: xdim, ydim, zdim, ndims, LO, HI

      implicit none

      class(named_array3d),                    intent(inout) :: this
      integer(kind=8), dimension(ndims,LO:HI), intent(in)    :: v

      real,    dimension(:,:,:), pointer             :: p3d

      if (.not.associated(this%arr)) then
         p3d => null()
      else
         p3d => this%arr(v(xdim,LO):v(xdim,HI),v(ydim,LO):v(ydim,HI),v(zdim,LO):v(zdim,HI))
      endif

   end function array3d_span_ijkse8

!> \brief Get a selected range of values of one variable

   function array4d_span_one_var(this,nn,v1,v2) result(p3d)

      use constants, only: xdim, ydim, zdim, ndims

      implicit none

      class(named_array4d),              intent(inout) :: this
      integer(kind=4),                   intent(in)    :: nn
      integer(kind=4), dimension(ndims), intent(in)    :: v1, v2

      real,    dimension(:,:,:),  pointer      :: p3d

      if (.not.associated(this%arr)) then
         p3d => null()
      else
         p3d => this%arr(nn,v1(xdim):v2(xdim),v1(ydim):v2(ydim),v1(zdim):v2(zdim))
      endif

   end function array4d_span_one_var

!> \brief Get a selected range of values of one variable

   function array4d_span_one_var_ijkse(this,nn,v) result(p3d)

      use constants, only: xdim, ydim, zdim, ndims, LO, HI

      implicit none

      class(named_array4d),                    intent(inout) :: this
      integer(kind=4),                         intent(in)    :: nn
      integer(kind=4), dimension(ndims,LO:HI), intent(in)    :: v

      real,    dimension(:,:,:),  pointer            :: p3d

      if (.not.associated(this%arr)) then
         p3d => null()
      else
         p3d => this%arr(nn,v(xdim,LO):v(xdim,HI),v(ydim,LO):v(ydim,HI),v(zdim,LO):v(zdim,HI))
      endif

   end function array4d_span_one_var_ijkse

!> \brief Get a selected line of values of all variables

   function array4d_span(this,v1,v2) result(p3d)

      use constants, only: xdim, ydim, zdim, ndims

      implicit none

      class(named_array4d),              intent(inout) :: this
      integer(kind=4), dimension(ndims), intent(in)    :: v1, v2

      real,    dimension(:,:,:,:), pointer     :: p3d

      if (.not.associated(this%arr)) then
         p3d => null()
      else
         p3d => this%arr(:,v1(xdim):v2(xdim),v1(ydim):v2(ydim),v1(zdim):v2(zdim))
      endif

   end function array4d_span

!> \brief Get a selected line of values of all variables

   function array4d_span_ijkse(this,v) result(p3d)

      use constants, only: xdim, ydim, zdim, ndims, LO, HI

      implicit none

      class(named_array4d),                    intent(inout) :: this
      integer(kind=4), dimension(ndims,LO:HI), intent(in)    :: v

      real,    dimension(:,:,:,:), pointer           :: p3d

      if (.not.associated(this%arr)) then
         p3d => null()
      else
         p3d => this%arr(:,v(xdim,LO):v(xdim,HI),v(ydim,LO):v(ydim,HI),v(zdim,LO):v(zdim,HI))
      endif

   end function array4d_span_ijkse

!> \brief Get the upper bound of a 4D named array
   function array4d_ubound(this,dim_) result(n)

      implicit none

      class(named_array4d), intent(in) :: this
      integer(kind=4),      intent(in) :: dim_
      integer(kind=4) :: n

      n = ubound(this%arr, dim=dim_, kind=4)

   end function array4d_ubound

!> \brief Get the lower bound of a 4D named array
   function array4d_lbound(this,dim_) result(n)

      implicit none

      class(named_array4d), intent(in) :: this
      integer(kind=4),      intent(in) :: dim_
      integer(kind=4) :: n

      n = lbound(this%arr,dim=dim_,kind=4)

   end function array4d_lbound

!> \brief Get the upper bound of a 3D named array
   function array3d_ubound(this,dim_) result(n)

      implicit none

      class(named_array3d), intent(in) :: this
      integer(kind=4),      intent(in) :: dim_
      integer(kind=4) :: n

      n = ubound(this%arr, dim=dim_, kind=4)

   end function array3d_ubound

!> \brief Get the lower bound of a 3D named array
   function array3d_lbound(this,dim_) result(n)

      implicit none

      class(named_array3d), intent(in) :: this
      integer(kind=4),      intent(in) :: dim_
      integer(kind=4) :: n

      n = lbound(this%arr,dim=dim_,kind=4)

   end function array3d_lbound

end module named_array
