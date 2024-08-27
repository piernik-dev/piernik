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

module cg_list

   use constants, only: dsetnamelen
   use grid_cont, only: grid_container

   implicit none

   private
   public :: cg_list_t, cg_list_element, cg_list_element_ptr

   !>
   !! \brief A grid container with two links to other cg_list_elements
   !!
   !! \details the prv and nxt pointers are not elements of the grid_container type to allow membership in several lists simultaneously
   !<
   type :: cg_list_element
      type(grid_container),  pointer :: cg       !< the current grid container
      type(cg_list_element), pointer :: prv, nxt !< pointers to previous and next grid container or null() at the end of the list
   end type cg_list_element

   !> \brief Just a wrapper around limitation of Fortran syntax
   type :: cg_list_element_ptr
      type(cg_list_element), pointer :: p
   end type cg_list_element_ptr

   !> \brief Arbitrary list of grid containers, not for direct use.
   type, abstract :: cg_list_t
      type(cg_list_element), pointer :: first !< first element of the chain of grid containers, the most important one
      type(cg_list_element), pointer :: last  !< last element of the chain - useful for quick expanding and merging lists
      integer(kind=4) :: cnt                  !< number of chain links
      character(len=dsetnamelen) :: label     !< name of the list for diagnostic and identification purposes
   contains
      ! List management
!      procedure       :: init_el
      procedure       :: init_new                          !< A constructor for an empty list
      procedure       :: add_new                           !< Add new element to the list
      generic, public :: add => add_new                    !< All methods of adding a new element
      procedure       :: del_lst                           !< Destroy the list
      procedure       :: un_link                           !< Un-link the element
      generic, public :: delete => un_link, del_lst        !< All methods of destroying
      procedure       :: print_list                        !< Print the list and associated cg ID (for debugging)
!> \todo merge lists
   end type cg_list_t

contains

!> \brief a constructor for an empty list
   subroutine init_new(this, label)

      implicit none

      class(cg_list_t), intent(inout) :: this  !< object invoking type-bound procedure
      character(len=*), intent(in)    :: label !< name of the list

      this%first => null()
      this%last  => null()
      this%cnt   =  0
      this%label =  trim(label)

   end subroutine init_new

!> \brief add new element to the list
   subroutine add_new(this, cg)

      use constants,  only: I_ONE
      use dataio_pub, only: die
      use grid_cont,  only: grid_container

      implicit none

      class(cg_list_t), intent(inout) :: this !< object invoking type-bound procedure
      type(grid_container), optional, pointer, intent(in) :: cg !< new grid container that will be added to add_new::this

      type(cg_list_element), pointer :: new

      allocate(new)
      if (present(cg)) then
         if (associated(cg)) then
            new%cg => cg
         else
            call die("[cg_list:add_new] tried to add null() element")
         endif
      else
         allocate(new%cg)
      endif
      new%nxt => null()

      if (.not. associated(this%first)) then ! the list was empty
         if (associated(this%last)) call die("[cg_list:add_new] last without first")
         this%first => new
         new%prv => null()
      else
         if (.not. associated(this%last)) call die("[cg_list:add_new] first without last")
         this%last%nxt => new
         new%prv => this%last
      endif

      this%last => new
      this%cnt = this%cnt + I_ONE

   end subroutine add_new

!> \brief destroy the list
   subroutine del_lst(this)

      use dataio_pub, only: die

      implicit none

      class(cg_list_t), intent(inout) :: this !< object invoking type-bound procedure

      type(cg_list_element), pointer  :: cgl

      do while (associated(this%first))
         cgl => this%last
         call this%delete(cgl) ! cannot just pass this%last because it will change after un_link and wrong element will be deallocated
      enddo

      if (this%cnt > 0) call die("[cg_list:del_lst] The list is still not empty")
      nullify(this%first)
      nullify(this%last)
      this%cnt = 0

   end subroutine del_lst

!> \brief Remove the element from the list, keep its contents

   subroutine un_link(this, cgle)

      use constants,  only: I_ONE
      use dataio_pub, only: die

      implicit none

      class(cg_list_t),               intent(inout) :: this !< object invoking type-bound procedure
      type(cg_list_element), pointer, intent(inout) :: cgle !< the element to be unlinked

      if (.not. associated(cgle)) call die("[cg_list:un_link] tried to remove null() element")
      if (.not. associated(this%first)) call die("[cg_list:un_link] Cannot remove anything from an empty list")
      if (this%cnt <= 0) call die("[cg_list:un_link] this%cnt <=0 .and. associated(this%first)")

      if (associated(this%first, cgle)) this%first => this%first%nxt
      if (associated(this%last,  cgle)) this%last  => this%last%prv
      if (associated(cgle%prv)) cgle%prv%nxt => cgle%nxt
      if (associated(cgle%nxt)) cgle%nxt%prv => cgle%prv
      deallocate(cgle)
      this%cnt = this%cnt - I_ONE

   end subroutine un_link

!> \brief print the list and associated cg ID for debugging purposes
   subroutine print_list(this)

      use constants,  only: V_DEBUG
      use dataio_pub, only: warn, printinfo, msg

      implicit none

      class(cg_list_t), intent(inout) :: this !< object invoking type-bound procedure

      type(cg_list_element), pointer :: cur
      integer :: cnt

      if (.not. associated(this%first)) then
         call warn("[cg_list:print_list] Empty list")
         if (this%cnt /= 0) then
            write(msg,'(a,i6)')"[cg_list:print_list] Empty list length /= 0 : ",this%cnt
            call warn(msg)
         endif
         if (associated(this%last)) then
            call warn("[cg_list:print_list] Tail without head")
            if (associated(this%last%cg)) then
               write(msg,'(a,i7)')"Last element #",this%last%cg%grid_id
               call warn(msg)
            endif
         endif
         return
      endif

      if (.not. associated(this%last)) call warn("[cg_list:print_list] Head without tail")

      if (associated(this%first%cg)) then
         write(msg,'(a,i7)')"First element #",this%first%cg%grid_id
         call printinfo(msg, V_DEBUG)
      endif
      if (associated(this%last%cg)) then
         write(msg,'(a,i7)')"Last element #",this%last%cg%grid_id
         call printinfo(msg, V_DEBUG)
      endif

      cnt = 0
      cur => this%first
      do while (associated(cur))
         cnt = cnt + 1
         if (associated(cur%cg)) then
            write(msg,'(i5,a,i7)')cnt,"-th element #",cur%cg%grid_id
            call printinfo(msg, V_DEBUG)
         endif
         cur => cur%nxt
      enddo

      if (cnt /= this%cnt) then
         write(msg, '(2(a,i5))')"[cg_list:print_list] this%cnt = ",this%cnt," /= ",cnt
         call warn(msg)
      endif

      cur => this%last
      do while (associated(cur))
         if (associated(cur%cg)) then
            write(msg,'(i5,a,i7)')cnt,"-th element #",cur%cg%grid_id
            call printinfo(msg, V_DEBUG)
         endif
         cnt = cnt - 1
         cur => cur%prv
      enddo

   end subroutine print_list

! unused
!!$!>
!!$!! \brief sets counter and pointer to the last element, updates pointer to the first element if necessary.
!!$!<
!!$   subroutine init_el(this, cgle)
!!$
!!$      use dataio_pub, only: warn, die
!!$
!!$      implicit none
!!$
!!$      class(cg_list_t), intent(inout) :: this !< object invoking type-bound procedure
!!$      type(cg_list_element), pointer, intent(inout) :: cgle
!!$
!!$      type(cg_list_element), pointer :: cur, prv
!!$
!!$      if (.not. associated(cgle)) then
!!$         call warn("[cg_list:init_el] tried to initialize with null() element")
!!$         return
!!$      endif
!!$
!!$      cur => cgle
!!$      if (associated(cur%prv)) call warn("[cg_list:init_el] this is not the first element of the chain")
!!$
!!$      do while (associated(cur%prv))
!!$         prv => cur
!!$         cur => cur%prv
!!$         if (.not. associated(cur%nxt, prv)) call die("[cg_list:init_el] this is not a straight list (rev)")
!!$         if (associated(cur, cgle)) call die("[cg_list:init_el] loops are not allowed (rev)")
!!$      enddo
!!$
!!$      this%first => cur
!!$      this%cnt = 1
!!$
!!$      do while (associated(cur%nxt))
!!$         prv => cur
!!$         cur => cur%nxt
!!$         this%cnt = this%cnt + 1
!!$         if (.not. associated(cur%prv, prv)) call die("[cg_list:init_el] this is not a straight list (fwd)")
!!$         ! we don't need second loop check here
!!$      enddo
!!$
!!$      this%last => cur
!!$
!!$   end subroutine init_el

end module cg_list
