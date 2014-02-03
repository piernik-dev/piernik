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
   public :: cg_list_T, cg_list_element

   !>
   !! \brief A grid container with two links to other cg_list_elements
   !!
   !! \details the prv and nxt pointers are not elements of the grid_container type to allow membership in several lists simultaneously
   !<
   type :: cg_list_element
      type(grid_container),  pointer :: cg       !< the current grid container
      type(cg_list_element), pointer :: prv, nxt !< pointers to previous and next grid container or null() at the end of the list
   end type cg_list_element

   !> \brief Arbitrary list of grid containers, not for direct use.
   type, abstract :: cg_list_T
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

      ! Misc
      procedure       :: print_list                        !< Print the list and associated cg ID
      procedure       :: numbered_ascii_dump               !< Construct name of emergency ASCII dump
      procedure       :: ascii_dump                        !< Emergency routine for quick ASCII dumps
      procedure       :: update_req                        !< Update mpisetup::req(:)
      procedure       :: prevent_prolong                   !< Mark grids as untouchable for prolongation
      procedure       :: enable_prolong                    !< Mark grids eligible for prolongation
      procedure       :: set_is_old                        !< Mark grids as existing in the previous timestep
      procedure       :: clear_ref_flags                   !< Clear refinement flags everywhere
      procedure       :: count_ref_flags                   !< Count refinement flags

!> \todo merge lists

   end type cg_list_T

contains

!> \brief a constructor for an empty list
   subroutine init_new(this, label)

      implicit none

      class(cg_list_T), intent(inout) :: this  !< object invoking type-bound procedure
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

      class(cg_list_T), intent(inout) :: this !< object invoking type-bound procedure
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

      class(cg_list_T), intent(inout) :: this !< object invoking type-bound procedure

      type(cg_list_element), pointer  :: cgl

      do while (associated(this%first))
         cgl => this%last
         call this%delete(cgl) ! cannot just pass this%last because if will change after un_link and wrong element will be deallocated
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

      class(cg_list_T),               intent(inout) :: this !< object invoking type-bound procedure
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

      use dataio_pub, only: warn, printinfo, msg

      implicit none

      class(cg_list_T), intent(inout) :: this !< object invoking type-bound procedure

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
         call printinfo(msg)
      endif
      if (associated(this%last%cg)) then
         write(msg,'(a,i7)')"Last element #",this%last%cg%grid_id
         call printinfo(msg)
      endif

      cnt = 0
      cur => this%first
      do while (associated(cur))
         cnt = cnt + 1
         if (associated(cur%cg)) then
            write(msg,'(i5,a,i7)')cnt,"-th element #",cur%cg%grid_id
            call printinfo(msg)
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
            call printinfo(msg)
         endif
         cnt = cnt - 1
         cur => cur%prv
      enddo
    end subroutine print_list

!> \brief Construct name of emergency ASCII dump

   subroutine numbered_ascii_dump(this, qlst, basename, a)

      use dataio_pub, only: halfstep, msg
      use global,     only: nstep, do_ascii_dump
      use mpisetup,   only: proc

      implicit none

      class(cg_list_T),              intent(inout) :: this     !< list for which do the dump (usually all_cg)
      integer(kind=4), dimension(:), intent(in)    :: qlst     !< list of scalar fields to be printed
      character(len=*),              intent(in)    :: basename !< first part of the filename
      integer, optional,             intent(in)    :: a        !< additional number

      integer                              :: l, n

      if (.not. do_ascii_dump) return

      n = 2 * nstep
      if (halfstep) n = n + 1

      if (present(a)) then
         write(msg, '(a,i4,i6,i3)') trim(basename), proc, n, a
      else
         write(msg, '(a,i4,i6)')    trim(basename), proc, n
      endif
      do l = 1, len_trim(msg)
         if (msg(l:l) == " ") msg(l:l) = "_"
      enddo
      call this%ascii_dump(trim(msg), qlst)

   end subroutine numbered_ascii_dump

!>
!! \brief Emergency routine for quick ASCII dumps
!!
!! \details Absolute integer coordinates also allow seamless concatenation of dumps made by all PEs.
!!
!! \warning This routine is intended only for debugging. It is strongly discouraged to use it for data dumps in production runs.
!<

   subroutine ascii_dump(this, filename, qlst)

      use dataio_pub,       only: msg, printio
      use named_array_list, only: qna

      implicit none

      class(cg_list_T),              intent(inout) :: this     !< list for which do the dump (usually all_cg)
      character(len=*),              intent(in)    :: filename !< name to write the emergency dump (should be different on each process)
      integer(kind=4), dimension(:), intent(in)    :: qlst     !< list of scalar fields to be printed

      integer, parameter                   :: fu=30
      integer                              :: i, j, k, q
      type(cg_list_element), pointer       :: cgl

      open(fu, file=filename, status="unknown")
      write(fu, '("#",a3,2a4,a6,3a20)', advance='no')"i", "j", "k", "level", "x(i)", "y(j)", "z(k)"
      do q = lbound(qlst(:), dim=1), ubound(qlst(:), dim=1)
         write(fu, '(a20)', advance='no') trim(qna%lst(qlst(q))%name)
      enddo
      write(fu, '(/)')

      cgl => this%first
      do while (associated(cgl))
         do i = cgl%cg%is, cgl%cg%ie
            do j = cgl%cg%js, cgl%cg%je
               do k = cgl%cg%ks, cgl%cg%ke
                  write(fu, '(3i4,i6,3es20.11e3)', advance='no') i, j, k, cgl%cg%level_id, cgl%cg%x(i), cgl%cg%y(j), cgl%cg%z(k)
                  do q = lbound(qlst(:), dim=1), ubound(qlst(:), dim=1)
                     write(fu, '(es20.11e3)', advance='no') cgl%cg%q(qlst(q))%arr(i, j, k)
                  enddo
                  write(fu, '()')
               enddo
               write(fu, '()')
            enddo
            write(fu, '()')
         enddo
         write(fu, '()')
         cgl => cgl%nxt
      enddo

      close(fu)

      write(msg,'(3a)') "[cg_list:ascii_dump] Wrote dump '",filename,"'"
      call printio(msg)

   end subroutine ascii_dump

!> \brief Update mpisetup::req(:)

   subroutine update_req(this)

      use mpisetup,  only: inflate_req

      implicit none

      class(cg_list_T), intent(in)   :: this

      integer                        :: nrq, d, dr, dp
      type(cg_list_element), pointer :: cgl

      ! calculate number of boundaries to communicate
      nrq = 0
      cgl => this%first
      do while (associated(cgl))

         do d = lbound(cgl%cg%i_bnd, dim=1), ubound(cgl%cg%i_bnd, dim=1)
            if (allocated(cgl%cg%i_bnd(d)%seg)) nrq = nrq + 2 * size(cgl%cg%i_bnd(d)%seg)
         enddo

         cgl => cgl%nxt
      enddo
      call inflate_req(nrq)

      ! calculate number of prolongation-restriction pairs
      nrq = 0
      cgl => this%first
      do while (associated(cgl))
         dr = 0
         if (allocated(cgl%cg%ri_tgt%seg)) dr =      size(cgl%cg%ri_tgt%seg(:), dim=1)
         if (allocated(cgl%cg%ro_tgt%seg)) dr = dr + size(cgl%cg%ro_tgt%seg(:), dim=1)

         dp = 0
         if (allocated(cgl%cg%pi_tgt%seg)) dp =      size(cgl%cg%pi_tgt%seg(:), dim=1)
         if (allocated(cgl%cg%po_tgt%seg)) dp = dp + size(cgl%cg%po_tgt%seg(:), dim=1)

         nrq = nrq + max(dr, dp)

         cgl => cgl%nxt
      enddo
      call inflate_req(nrq)

   end subroutine update_req

!> \brief Mark grids as untouchable for prolongation

   subroutine prevent_prolong(this)

      implicit none

      class(cg_list_T), intent(in)   :: this

      type(cg_list_element), pointer :: cgl

      cgl => this%first
      do while (associated(cgl))
         cgl%cg%ignore_prolongation = .true.
         cgl => cgl%nxt
      enddo

   end subroutine prevent_prolong

!> \brief Mark grids as eligible for prolongation

   subroutine enable_prolong(this)

      implicit none

      class(cg_list_T), intent(in)   :: this

      type(cg_list_element), pointer :: cgl

      cgl => this%first
      do while (associated(cgl))
         cgl%cg%ignore_prolongation = .false.
         cgl => cgl%nxt
      enddo

   end subroutine enable_prolong

!> \brief Mark grids as existing in the previous timestep

   subroutine set_is_old(this)

      implicit none

      class(cg_list_T), intent(in)   :: this

      type(cg_list_element), pointer :: cgl

      cgl => this%first
      do while (associated(cgl))
         cgl%cg%is_old = .true.
         cgl => cgl%nxt
      enddo

   end subroutine set_is_old

!> \brief Clear refinement flags everywhere

   subroutine clear_ref_flags(this)

      implicit none

      class(cg_list_T), intent(in) :: this !< object invoking type-bound procedure

      type(cg_list_element), pointer :: cgl

      cgl => this%first
      do while (associated(cgl))
         cgl%cg%refine_flags%refine   = .false.
         cgl%cg%refine_flags%derefine = .false.
         cgl => cgl%nxt
      enddo

   end subroutine clear_ref_flags

!> \brief Count refinement flags everywhere

   function count_ref_flags(this) result(cnt)

      implicit none

      class(cg_list_T), intent(in) :: this !< object invoking type-bound procedure
      integer :: cnt                       !< returned counter

      type(cg_list_element), pointer :: cgl

      cnt = 0
      cgl => this%first
      do while (associated(cgl))
         if ( cgl%cg%refine_flags%refine .or. &
              any(cgl%cg%refinemap .and. cgl%cg%leafmap) .or. &
              size(cgl%cg%refine_flags%SFC_refine_list) > 0) cnt = cnt + 1
         cgl => cgl%nxt
      enddo

   end function count_ref_flags

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
!!$      class(cg_list_T), intent(inout) :: this !< object invoking type-bound procedure
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
