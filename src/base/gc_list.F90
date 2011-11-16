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

!> \brief This module contains grid container list and related methods

module gc_list

   use grid_cont, only: grid_container

   implicit none

   private
   public :: cg_list, cg_sublist, cg_list_element

   !> \brief A grid container with two links to other cg_list_elements
   !> \details the prv and nxt pointers are not elements of the grid_container type to allow membership in several lists simultaneously
   type cg_list_element
      type(grid_container),  pointer :: cg       !< the current grid container
      type(cg_list_element), pointer :: prv, nxt !< pointers to previous and next grid container or null() at the end of the list
   end type cg_list_element

   !> \brief A list of grid containers
   type cg_sublist

      type(cg_list_element), pointer :: first !< first element of the chain of grid containers, the most important one
      type(cg_list_element), pointer :: last  !< last element of the chain - useful for expanding and merging lists
      integer :: cnt                          !< number of chain links

    contains

!      procedure :: init_el
      procedure :: init_new
      generic, public :: init => init_new!, init_el

      procedure :: add_new
!      procedure :: add_el
      generic, public :: add => add_new !, add_el

      procedure :: un_link
      procedure :: get_extremum

   end type cg_sublist

   !>
   !! \brief A list of grid containers that are supposed to have the same variables registered
   !!
   !! \details The main purpose of this type is to provide a type for the set of all grid containers with methods and properties
   !! that should not be available for any arbitrarily composed subset of grid containers. Typically there will be only one variable
   !! of this type available in the code: grid::all_cg. Perhaps the multigrid solver will use another such variable for the levels
   !! below the base level, where it is not necessary to keep all the fluid and magnetic field data.
   !!
   !! It should be possible to use this type for more fancy, multi-domain  grid configurations such as:
   !! - Yin-Yang grid (covering sphere with two domains shaped as parts of the sphere to avoid  polar singulaties),
   !! - Cylindrical grid with cartesian core covering singularity at the axis.
   !! - Simulations with mixed dimensionality (e.g. 2d grid for dust particles and 3d grid for gas) should probably also use separate cg_list
   !! for their data (and additional routine for coupling the two grid sets).
   !<
   type, extends(cg_sublist) :: cg_list
   !> \todo store the information about registered entries also here (name, restart_mode, index, dim4)
    contains
      procedure :: reg_var
   end type cg_list

contains
!>
!! \brief sets counter and pointer to the last element, updates pointer to the first element if necessary.
!<

! unused
!!$   subroutine init_el(this, cgle)
!!$
!!$      use dataio_pub, only: warn, die
!!$
!!$      implicit none
!!$
!!$      class(cg_sublist), intent(inout) :: this
!!$      type(cg_list_element), pointer, intent(inout) :: cgle
!!$
!!$      type(cg_list_element), pointer :: cur, prv
!!$
!!$      if (.not. associated(cgle)) then
!!$         call warn("[gc_list:init_el] tried to initialize with null() element")
!!$         return
!!$      endif
!!$
!!$      cur => cgle
!!$      if (associated(cur%prv)) call warn("[gc_list:init_el] this is not the first element of the chain")
!!$
!!$      do while (associated(cur%prv))
!!$         prv => cur
!!$         cur => cur%prv
!!$         if (.not. associated(cur%nxt, prv)) call die("[gc_list:init_el] this is not a straight list (rev)")
!!$         if (associated(cur, cgle)) call die("[gc_list:init_el] loops are not allowed (rev)")
!!$      enddo
!!$
!!$      this%first => cur
!!$      this%cnt = 1
!!$
!!$      do while (associated(cur%nxt))
!!$         prv => cur
!!$         cur => cur%nxt
!!$         this%cnt = this%cnt + 1
!!$         if (.not. associated(cur%prv, prv)) call die("[gc_list:init_el] this is not a straight list (fwd)")
!!$         ! we don't need second loop check here
!!$      enddo
!!$
!!$      this%last => cur
!!$
!!$   end subroutine init_el

!>
!! \brief a constructor for an empty list
!<

   subroutine init_new(this)

      implicit none

      class(cg_sublist), intent(out) :: this

      this%first => null()
      this%last => null()
      this%cnt = 0

   end subroutine init_new

!>
!! \brief add new element to the list
!<

   subroutine add_new(this)

      use dataio_pub, only: die

      implicit none

      class(cg_sublist), intent(inout) :: this

      type(cg_list_element), pointer :: new

      allocate(new)
      allocate(new%cg)
      new%nxt => null()

      if (.not. associated(this%first)) then ! the list was empty
         if (associated(this%last)) call die("[gc_list:add_new] last without first")
         this%first => new
         new%prv => null()
      else
         if (.not. associated(this%last)) call die("[gc_list:add_new] first without last")
         this%last%nxt => new
         new%prv => this%last
      endif

      this%last => new
      this%cnt = this%cnt + 1

      !> \todo update the lists of registered variables

   end subroutine add_new

! unused
!!$   subroutine add_el(this, cgle)
!!$
!!$      use dataio_pub, only: die, warn
!!$
!!$      implicit none
!!$
!!$      class(cg_sublist), intent(inout) :: this
!!$      type(cg_list_element), pointer, intent(inout) :: cgle
!!$
!!$      if (.not. associated(cgle)) then
!!$         call warn("[gc_list:add_el] tried to add null() element")
!!$         return
!!$      endif
!!$
!!$      if (associated(cgle%prv)) call die("[gc_list:add_el] this is not a first link")
!!$      if (associated(cgle%nxt)) call die("[gc_list:add_el] this is a chain") !> todo convert to warn and count the links properly
!!$
!!$      cgle%prv => this%last
!!$      this%last%nxt => cgle
!!$      this%cnt = this%cnt + 1
!!$
!!$   end subroutine add_el

!>
!! \brief destroy the element
!<

   subroutine un_link(this, cgle)

      use dataio_pub, only: warn, die

      implicit none

      class(cg_sublist), intent(inout) :: this
      type(cg_list_element), pointer, intent(inout) :: cgle

      type(cg_list_element), pointer :: cur
      integer :: cnt

      if (.not. associated(cgle)) then
         call warn("[gc_list:un_link] tried to remove null() element")
         return
      endif

      cur => this%first
      cnt = this%cnt

      if (.not. associated(this%first)) call die("[gc_list:un_link] Cannot remove anything from an empty list")
      if (cnt <= 0) call die("[gc_list:un_link] this%cnt <=0 .and. associated(this%first)")

      do while (associated(cur))
         if (associated(cur, cgle)) then
            if (associated(this%first, cgle)) this%first => this%first%nxt
            if (associated(this%last,  cgle)) this%last  => this%last%prv
            if (associated(cur%prv)) cur%prv%nxt => cur%nxt
            if (associated(cur%nxt)) cur%nxt%prv => cur%prv
            nullify(cur%nxt, cur%prv)
            this%cnt = this%cnt - 1
            exit
         endif

         cur => cur%nxt
      enddo

      if (this%cnt == cnt) call warn("[gc_list:un_link] element not found on the list")

   end subroutine un_link

!> \todo merge lists

!>
!! \brief Use this routine to add a variable (cg%q or cg%w) to all grid containers.
!!
!! \details Register a rank-3 array of given name in each grid container (in cg%q) and decide what to do with it on restart.
!! When dim4 is present then create a rank-4 array instead.(in cg%w)
!<

   subroutine reg_var(this, name, restart_mode, dim4)

      use dataio_pub, only: die

      implicit none

      class(cg_list), intent(in)    :: this
      character(len=*), intent(in)  :: name
      integer(kind=4), intent(in)   :: restart_mode
      integer(kind=4), optional, intent(in) :: dim4

      type(cg_list_element), pointer :: cgl

      cgl => this%first
      do while (associated(cgl))
         if (present(dim4)) then
            if (dim4<=0) call die("[gc_list:reg_var] dim4<=0")
            call cgl%cg%add_na_4d(name, restart_mode, dim4)
         else
            call cgl%cg%add_na(name, restart_mode)
         endif
         cgl => cgl%nxt
      enddo

   end subroutine reg_var

!>
!! \brief Find munimum or maximum value over a specified list of grid containers
!!
!! \details It should be possible to find an extremum over a given level or leaf blocks or something
!<
   subroutine get_extremum(this, ind, minmax, prop)

      use constants,  only: MINL, MAXL, I_ONE, ndims, xdim, ydim, zdim
      use dataio_pub, only: msg, warn, die
      use domain,     only: dom
      use grid_cont,  only: grid_container
      use mpi,        only: MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_STATUS_IGNORE, MPI_2DOUBLE_PRECISION, MPI_MINLOC, MPI_MAXLOC, MPI_IN_PLACE
      use mpisetup,   only: comm, ierr, master, proc, FIRST
      use types,      only: value

      implicit none

      class(cg_sublist), intent(in)  :: this
      integer(kind=4),   intent(in)  :: ind  !< Index in cg%q(:)
      integer(kind=4),   intent(in)  :: minmax
      type(value),       intent(out) :: prop


      type(grid_container), pointer :: cg
      type(cg_list_element), pointer :: cgl
      integer, parameter :: tag1 = 11, tag2 = tag1 + 1
      integer, dimension(MINL:MAXL), parameter :: op = [ MPI_MINLOC, MPI_MAXLOC ]
      real, dimension(:,:,:), pointer :: tab
      enum, bind(C)
         enumerator :: I_V, I_P
      end enum
      real, dimension(I_V:I_P)  :: v_red

      prop%loc(:) = 0
      select case (minmax)
         case (MINL)
            prop%val = huge(1.)
         case (MAXL)
            prop%val = -huge(1.)
         case default
            write(msg,*) "[gc_list:get_extremum]: I don't know what to do with minmax = ", minmax
            call warn(msg)
      end select

      cgl => this%first
      if (ind > ubound(cgl%cg%q(:), dim=1) .or. ind < lbound(cgl%cg%q(:), dim=1)) call die("[gc_list:get_extremum] Wrong index")
      do while (associated(cgl))
         cg => cgl%cg

         tab => cg%q(ind)%arr(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke)
         select case (minmax)
            case (MINL)
               if (minval(tab) < prop%val) then
                  prop%val = minval(tab)
                  prop%loc = minloc(tab) + cg%nb
               endif
            case (MAXL)
               if (maxval(tab) > prop%val) then
                  prop%val = maxval(tab)
                  prop%loc = maxloc(tab) + cg%nb
               endif
         end select
         cgl => cgl%nxt
      enddo

      v_red(I_V) = prop%val; v_red(I_P) = real(proc)

      call MPI_Allreduce(MPI_IN_PLACE, v_red, I_ONE, MPI_2DOUBLE_PRECISION, op(minmax), comm, ierr)

      prop%val = v_red(I_V)
      prop%proc = int(v_red(I_P))

      if (proc == prop%proc) then
         where (.not. dom%has_dir(:)) prop%coords(:) = 0.
         if (dom%has_dir(xdim)) prop%coords(xdim) = cg%x(prop%loc(xdim))
         if (dom%has_dir(ydim)) prop%coords(ydim) = cg%y(prop%loc(ydim))
         if (dom%has_dir(zdim)) prop%coords(zdim) = cg%z(prop%loc(zdim))
      endif

      if (prop%proc /= 0) then
         if (proc == prop%proc) then ! slave
            call MPI_Send (prop%loc,    ndims, MPI_INTEGER,          FIRST, tag1, comm, ierr)
            call MPI_Send (prop%coords, ndims, MPI_DOUBLE_PRECISION, FIRST, tag2, comm, ierr)
         endif
         if (master) then
            call MPI_Recv (prop%loc,    ndims, MPI_INTEGER,          prop%proc, tag1, comm, MPI_STATUS_IGNORE, ierr)
            call MPI_Recv (prop%coords, ndims, MPI_DOUBLE_PRECISION, prop%proc, tag2, comm, MPI_STATUS_IGNORE, ierr)
         endif
      endif

   end subroutine get_extremum

end module gc_list
