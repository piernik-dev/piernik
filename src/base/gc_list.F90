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
!! \brief This module contains grid container list and related methods
!<

module gc_list

   use grid_cont, only: grid_container

   implicit none

   private
   public :: cg_list, cg_list_element

   ! the prv and nxt pointers are not elements of the grid_container type to allow membership in several lists simultaneously
   type cg_list_element
      type(grid_container),  pointer :: cg       !< the current grid container
      type(cg_list_element), pointer :: prv, nxt !< pointers to previous and next grid container or null() at the end of the list
   end type cg_list_element

   type cg_list

      type(cg_list_element), pointer :: first !< first element of the chain of grid containers, the most important one
      type(cg_list_element), pointer :: last  !< last element of the chain - useful for expanding and merging lists
      integer :: cnt                          !< number of chain links

    contains

      procedure :: init_el
      procedure :: init_new
      generic, public :: init => init_new, init_el

      procedure :: add_new
      procedure :: add_el
      generic, public :: add => add_new, add_el

      procedure :: un_link

      procedure :: reg_var

   end type cg_list

contains
!>
!! \brief sets counter and pointer to the last element, updates pointer to the first element if necessary.
!<

   subroutine init_el(this, cgle)

      use dataio_pub, only: warn, die

      implicit none

      class(cg_list), intent(inout) :: this
      type(cg_list_element), pointer, intent(inout) :: cgle

      type(cg_list_element), pointer :: cur, prv

      if (.not. associated(cgle)) then
         call warn("[gc_list:init_el] tried to initialize with null() element")
         return
      endif

      cur => cgle
      if (associated(cur%prv)) call warn("[gc_list:init_el] this is not the first element of the chain")

      do while (associated(cur%prv))
         prv => cur
         cur => cur%prv
         if (.not. associated(cur%nxt, prv)) call die("[gc_list:init_el] this is not a straight list (rev)")
         if (associated(cur, cgle)) call die("[gc_list:init_el] loops are not allowed (rev)")
      enddo

      this%first => cur
      this%cnt = 1

      do while (associated(cur%nxt))
         prv => cur
         cur => cur%nxt
         this%cnt = this%cnt + 1
         if (.not. associated(cur%prv, prv)) call die("[gc_list:init_el] this is not a straight list (fwd)")
         ! we don't need second loop check here
      enddo

      this%last => cur

   end subroutine init_el

!>
!! \brief a constructor for an empty list
!<

   subroutine init_new(this)

      implicit none

      class(cg_list), intent(out) :: this

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

      class(cg_list), intent(inout) :: this

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

   end subroutine add_new

   subroutine add_el(this, cgle)

      use dataio_pub, only: die, warn

      implicit none

      class(cg_list), intent(inout) :: this
      type(cg_list_element), pointer, intent(inout) :: cgle

      if (.not. associated(cgle)) then
         call warn("[gc_list:add_el] tried to add null() element")
         return
      endif

      if (associated(cgle%prv)) call die("[gc_list:add_el] this is not a first link")
      if (associated(cgle%nxt)) call die("[gc_list:add_el] this is a chain") !> todo convert to warn and count the links properly

      cgle%prv => this%last
      this%last%nxt => cgle
      this%cnt = this%cnt + 1

   end subroutine add_el

!>
!! \brief destroy the element
!<

   subroutine un_link(this, cgle)

      use dataio_pub, only: warn, die

      implicit none

      class(cg_list), intent(inout) :: this
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

end module gc_list
