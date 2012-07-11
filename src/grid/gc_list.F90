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
   public :: cg_list, cg_list_element, ind_val, dirty_label

   !>
   !! \brief A grid container with two links to other cg_list_elements
   !!
   !! \details the prv and nxt pointers are not elements of the grid_container type to allow membership in several lists simultaneously
   !<
   type cg_list_element
      type(grid_container),  pointer :: cg       !< the current grid container
      type(cg_list_element), pointer :: prv, nxt !< pointers to previous and next grid container or null() at the end of the list
   end type cg_list_element

   !> \brief Abitrary list of grid containers
   type cg_list

      type(cg_list_element), pointer :: first !< first element of the chain of grid containers, the most important one
      type(cg_list_element), pointer :: last  !< last element of the chain - useful for quick expanding and merging lists
      integer :: cnt                          !< number of chain links

    contains

      ! List management
!      procedure :: init_el
      procedure :: init_new                          !< A constructor for an empty list
      generic, public :: init => init_new!, init_el  !< All constructors of a new list
      procedure :: add_new                           !< Add new element to the list
      generic, public :: add => add_new              !< All methods of adding a new element
      procedure :: del_lnk                           !< Destroy the element
      procedure :: del_lst                           !< Destroy the list
      generic, public :: delete => del_lnk, del_lst  !< All methods of destroying

      procedure :: un_link                           !< Un-link the element

      ! Misc
      procedure :: get_extremum                      !< Find munimum or maximum value over a s list
      procedure :: print_list                        !< Print the list and associated cg ID
      procedure :: numbered_ascii_dump               !< Construct name of emergency ASCII dump
      procedure :: ascii_dump                        !< Emergency routine for quick ASCII dumps
      procedure :: set_dirty                         !< Pollute selected array with an insane value dirtyH.
      procedure :: check_dirty                       !< Check for detectable traces of set_dirty calls.
      procedure :: check_for_dirt                    !< Check all named arrays for constants:big_float
      procedure :: update_req                        !< Update mpisetup::req(:)

      ! Arithmetic on the fields
      procedure :: set_q_value                       !< reset given field to the value
      procedure :: q_copy                            !< copy a given field to another
      procedure :: qw_copy                           !< copy a given rank-3 field to a component of rank-4 field
      procedure :: wq_copy                           !< copy a component of rank-4 field to a given rank-3 field
      procedure :: q_add                             !< add a field to another
      procedure :: q_add_val                         !< add a value to a field
      procedure :: q_lin_comb                        !< assign linear combination of q fields
      procedure :: subtract_average                  !< subtract average value from the list
      procedure :: norm_sq                           !< calculate L2 norm

      ! Multigrid
      procedure :: zero_boundaries                   !< Clear boundary values
!> \todo merge lists

   end type cg_list

   !> \brief Index - value pairs for calling arithmetics on the grids with q_lin_comb
   type ind_val
      integer :: ind  !< index in cg%q
      real    :: val  !< value for multiplication
   end type ind_val

   integer, parameter    :: dl_len = 64 !< length of label buffer
   character(len=dl_len) :: dirty_label !< buffer for label for check_dirty subroutine

contains

!> \brief a constructor for an empty list
   subroutine init_new(this)

      implicit none

      class(cg_list), intent(inout) :: this !< object invoking type-bound procedure

      this%first => null()
      this%last => null()
      this%cnt = 0

   end subroutine init_new

!> \brief add new element to the list
   subroutine add_new(this, cg)

      use dataio_pub, only: die
      use grid_cont,  only: grid_container

      implicit none

      class(cg_list), intent(inout) :: this !< object invoking type-bound procedure
      type(grid_container), optional, pointer, intent(in) :: cg !< new grid container that will be added to add_new::this

      type(cg_list_element), pointer :: new

      allocate(new)
      if (present(cg)) then
         if (associated(cg)) then
            new%cg => cg
         else
            call die("[gc_list:add_new] tried to add null() element")
         endif
      else
         allocate(new%cg)
      endif
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

!> \brief destroy the element
   subroutine del_lnk(this, cgle)

      use constants,  only: INVALID
      use dataio_pub, only: warn

      implicit none

      class(cg_list), intent(inout) :: this !< object invoking type-bound procedure
      type(cg_list_element), pointer, intent(inout) :: cgle !< the element to be eradicated

      if (.not. associated(cgle)) then
         call warn("[gc_list:del_lnk] tried to remove null() element")
         return
      endif

      call this%un_link(cgle)
      if (associated(cgle%cg)) then
         if (cgle%cg%grid_id > INVALID) then ! dirty trick
            call cgle%cg%cleanup
!            deallocate(cgle%cg) ! if we deallocate now, we won't be able to determine this in the other lists
         endif
      endif
      deallocate(cgle)

   end subroutine del_lnk

!> \brief destroy the list
   subroutine del_lst(this)

     implicit none

     class(cg_list), intent(inout) :: this !< object invoking type-bound procedure
     type(cg_list_element), pointer :: cgl

     do while (associated(this%first))
        cgl => this%last
        call this%delete(cgl) ! cannot juss pass this%last because if will change after un_link and wrong element will be deallocated
     enddo

   end subroutine del_lst

!> \brief remove the element from the list, keep its contents
   subroutine un_link(this, cgle)

      use dataio_pub, only: warn, die

      implicit none

      class(cg_list), intent(inout) :: this !< object invoking type-bound procedure
      type(cg_list_element), pointer, intent(in) :: cgle !< the element to be unlinked

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

!> \brief print the list and associated cg ID for debugging purposes
   subroutine print_list(this)

      use dataio_pub, only: warn, printinfo, msg

      implicit none

      class(cg_list), intent(inout) :: this !< object invoking type-bound procedure

      type(cg_list_element), pointer :: cur
      integer :: cnt

      if (.not. associated(this%first)) then
         call warn("[gc_list:print_list] Empty list")
         if (this%cnt /= 0) then
            write(msg,'(a,i6)')"[gc_list:print_list] Empty list length /= 0 : ",this%cnt
            call warn(msg)
         endif
         if (associated(this%last)) then
            call warn("[gc_list:print_list] Tail without head")
            if (associated(this%last%cg)) then
               write(msg,'(a,i7)')"Last element #",this%last%cg%grid_id
               call warn(msg)
            endif
         endif
         return
      endif

      if (.not. associated(this%last)) call warn("[gc_list:print_list] Head without tail")

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
         write(msg, '(2(a,i5))')"[gc_list:print_list] this%cnt = ",this%cnt," /= ",cnt
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

!>
!! \brief Find munimum or maximum value over a specified list of grid containers
!!
!! \details It should be possible to find an extremum over a given level or leaf blocks or something
!<
   subroutine get_extremum(this, ind, minmax, prop, dir)

      use constants,  only: MINL, MAXL, I_ONE, ndims, xdim, ydim, zdim, big_float
      use dataio_pub, only: msg, warn, die
      use domain,     only: dom
      use grid_cont,  only: grid_container
      use mpi,        only: MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_STATUS_IGNORE, MPI_2DOUBLE_PRECISION, MPI_MINLOC, MPI_MAXLOC, MPI_IN_PLACE
      use mpisetup,   only: comm, mpi_err, master, proc, FIRST
      use types,      only: value

      implicit none

      class(cg_list),            intent(in)  :: this    !< object invoking type-bound procedure
      integer,                   intent(in)  :: ind     !< Index in cg%q(:)
      integer(kind=4),           intent(in)  :: minmax  !< minimum or maximum ?
      type(value),               intent(out) :: prop    !< precise location of the extremum to be found
      integer(kind=4), optional, intent(in)  :: dir     !< order the cell size in dir direction


      type(grid_container),   pointer :: cg, cg_x
      type(cg_list_element),  pointer :: cgl
      real, dimension(:,:,:), pointer :: tab
      integer,                       parameter :: tag1 = 11, tag2 = tag1 + 1, tag3 = tag2 + 1
      integer, dimension(MINL:MAXL), parameter :: op = [ MPI_MINLOC, MPI_MAXLOC ]
      enum, bind(C)
         enumerator :: I_V, I_P !< value and proc
      end enum
      real, dimension(I_V:I_P)  :: v_red

      if (.not. present(dir)) prop%assoc = - big_float ! means uninitialized

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

      nullify(cg_x)
      cgl => this%first
      if (ind > ubound(cgl%cg%q(:), dim=1) .or. ind < lbound(cgl%cg%q(:), dim=1)) call die("[gc_list:get_extremum] Wrong index")
      do while (associated(cgl))
         cg => cgl%cg

         tab => cg%q(ind)%span(cg%ijkse)
         select case (minmax)
            case (MINL)
               if (minval(tab) < prop%val) then
                  prop%val = minval(tab)
                  prop%loc = minloc(tab) + dom%nb
                  cg_x => cg
               endif
            case (MAXL)
               if (maxval(tab) > prop%val) then
                  prop%val = maxval(tab)
                  prop%loc = maxloc(tab) + dom%nb
                  cg_x => cg
               endif
         end select
         cgl => cgl%nxt
      enddo

      v_red(I_V) = prop%val; v_red(I_P) = real(proc)

      call MPI_Allreduce(MPI_IN_PLACE, v_red, I_ONE, MPI_2DOUBLE_PRECISION, op(minmax), comm, mpi_err)

      prop%val = v_red(I_V)
      prop%proc = int(v_red(I_P))

      if (proc == prop%proc) then
         where (.not. dom%has_dir(:)) prop%coords(:) = 0.
         if (associated(cg_x)) then
            if (dom%has_dir(xdim)) prop%coords(xdim) = cg_x%x(prop%loc(xdim))
            if (dom%has_dir(ydim)) prop%coords(ydim) = cg_x%y(prop%loc(ydim))
            if (dom%has_dir(zdim)) prop%coords(zdim) = cg_x%z(prop%loc(zdim))
            if (present(dir))      prop%assoc        = cg_x%dl(dir)
            ! else prop%assoc = minval(cg_x%dl(:)) ???
         else
            call die("[gc_list:get_extremum] not associated(cg_x)")
         endif
      endif

      if (prop%proc /= 0) then
         if (proc == prop%proc) then ! slave
            call MPI_Send (prop%loc,    ndims, MPI_INTEGER,          FIRST, tag1, comm, mpi_err)
            call MPI_Send (prop%coords, ndims, MPI_DOUBLE_PRECISION, FIRST, tag2, comm, mpi_err)
            if (present(dir)) call MPI_Send (prop%assoc, I_ONE, MPI_DOUBLE_PRECISION, FIRST, tag3, comm, mpi_err)
         endif
         if (master) then
            call MPI_Recv (prop%loc,    ndims, MPI_INTEGER,          prop%proc, tag1, comm, MPI_STATUS_IGNORE, mpi_err)
            call MPI_Recv (prop%coords, ndims, MPI_DOUBLE_PRECISION, prop%proc, tag2, comm, MPI_STATUS_IGNORE, mpi_err)
            if (present(dir)) call MPI_Recv (prop%assoc, I_ONE, MPI_DOUBLE_PRECISION, prop%proc, tag3, comm, MPI_STATUS_IGNORE, mpi_err)
         endif
      endif

   end subroutine get_extremum

!> \brief reset given field to the value (usually 0. or dirty)

   subroutine set_q_value(this, ind, val)

      implicit none

      class(cg_list), intent(in) :: this    !< object invoking type-bound procedure
      integer,        intent(in) :: ind     !< Index in cg%q(:)
      real,           intent(in) :: val     !< value to put

      type(cg_list_element), pointer :: cgl

      cgl => this%first
      do while (associated(cgl))
         cgl%cg%q(ind)%arr(:, :, :) = val
         cgl => cgl%nxt
      enddo

   end subroutine set_q_value

!> \brief copy a given field to another

   subroutine q_copy(this, i_from, i_to)

      implicit none

      class(cg_list), intent(in) :: this    !< object invoking type-bound procedure
      integer,        intent(in) :: i_from  !< Index of source in cg%q(:)
      integer,        intent(in) :: i_to    !< Index of destination in cg%q(:)

      type(cg_list_element), pointer :: cgl

      cgl => this%first
      do while (associated(cgl))
         cgl%cg%q(i_to)%arr(:, :, :) = cgl%cg%q(i_from)%arr(:, :, :)
         cgl => cgl%nxt
      enddo

   end subroutine q_copy

!> \brief copy a given rank-3 field to a component of a rank-4 field (cg%w(w_to)%arr(w_ind,:,:,:) = cg%q(q_from)%arr(:, :, :))

   subroutine qw_copy(this, q_from, w_to, w_ind)

      implicit none

      class(cg_list), intent(in) :: this    !< object invoking type-bound procedure
      integer,        intent(in) :: q_from  !< Index of source in cg%q(:)
      integer,        intent(in) :: w_to    !< Index of destination in cg%w(:)
      integer,        intent(in) :: w_ind   !< First index of destination in cg%w(w_to)%arr(:,:,:,:)

      type(cg_list_element), pointer :: cgl

      cgl => this%first
      do while (associated(cgl))
         cgl%cg%w(w_to)%arr(w_ind, :, :, :) = cgl%cg%q(q_from)%arr(:, :, :)
         cgl => cgl%nxt
      enddo

   end subroutine qw_copy

!> \brief copy a component of rank-4 field to a given rank-3 field (cg%q(q_to)%arr(:, :, :) = cg%w(w_from)%arr(w_ind,:,:,:))

   subroutine wq_copy(this, w_from, w_ind, q_to)

      implicit none

      class(cg_list), intent(in) :: this    !< object invoking type-bound procedure
      integer,        intent(in) :: w_from  !< Index of source in cg%w(:)
      integer,        intent(in) :: w_ind   !< First index of source in cg%w(w_from)%arr(:,:,:,:)
      integer,        intent(in) :: q_to    !< Index of destination in cg%q(:)

      type(cg_list_element), pointer :: cgl

      cgl => this%first
      do while (associated(cgl))
         cgl%cg%q(q_to)%arr(:, :, :) = cgl%cg%w(w_from)%arr(w_ind, :, :, :)
         cgl => cgl%nxt
      enddo

   end subroutine wq_copy

!> \brief Add a field to another (e.g. apply a correction)

   subroutine q_add(this, i_add, i_to)

      implicit none

      class(cg_list), intent(in) :: this    !< object invoking type-bound procedure
      integer,        intent(in) :: i_add   !< Index of field to be aded in cg%q(:)
      integer,        intent(in) :: i_to    !< Index of field to be modified in cg%q(:)


      type(cg_list_element), pointer :: cgl

      cgl => this%first
      do while (associated(cgl))
         cgl%cg%q(i_to)%arr(:, :, :) = cgl%cg%q(i_to)%arr(:, :, :) + cgl%cg%q(i_add)%arr(:, :, :)
         cgl => cgl%nxt
      enddo

   end subroutine q_add

!> \brief Add a value to a field (e.g. correct for average value)

   subroutine q_add_val(this, i_add, val)

      implicit none

      class(cg_list), intent(in) :: this    !< object invoking type-bound procedure
      integer,        intent(in) :: i_add   !< Index of field to be modified in cg%q(:)
      real,           intent(in) :: val     !< Value to be added


      type(cg_list_element), pointer :: cgl

      cgl => this%first
      do while (associated(cgl))
         cgl%cg%q(i_add)%arr(:, :, :) = cgl%cg%q(i_add)%arr(:, :, :) + val
         cgl => cgl%nxt
      enddo

   end subroutine q_add_val

!>
!! \brief Assign linear combination of q fields on the whole list
!!
!! \details On the whole list cg%q(ind) is assigned to sum of iv%val * cg%q(iv%ind)
!!
!! \todo add an option to select region of cg
!<

   subroutine q_lin_comb(this, iv, ind)

      use dataio_pub, only: die, warn

      implicit none

      class(cg_list),              intent(in) :: this    !< object invoking type-bound procedure
      type(ind_val), dimension(:), intent(in) :: iv      !< list of (coefficient, index) pairs
      integer,                     intent(in) :: ind     !< Index in cg%q(:)

      integer :: i
      type(ind_val), dimension(size(iv)) :: iv_safe !< sanitized copy of iv
      logical :: swapped
      type(cg_list_element), pointer :: cgl

      if (size(iv) <= 0) then
         call warn("[gc_list::q_lin_comb] Nothing to do")
         return
      endif

      iv_safe(:) = iv(:)
      ! if own field (ind) is involved then move it to the first position to avoid side effects and allow in-place operation
      swapped =.false.
      do i = lbound(iv, dim=1)+1, ubound(iv, dim=1)
         if (ind == iv(i)%ind) then
            if (swapped) call die("[gc_list::q_lin_comb] Cannot use own field twice due to side effects")
            iv_safe(lbound(iv, dim=1)) = iv(i)
            iv_safe(i) = iv(lbound(iv, dim=1))
            swapped = .true.
         endif
      enddo

      cgl => this%first
      do while (associated(cgl))
         cgl%cg%q(ind)%arr(:, :, :) = iv_safe(lbound(iv_safe, dim=1))%val * cgl%cg%q(iv_safe(lbound(iv_safe, dim=1))%ind)%arr(:, :, :)
         do i = lbound(iv_safe, dim=1)+1, ubound(iv_safe, dim=1)
            cgl%cg%q(ind)%arr(:, :, :) = cgl%cg%q(ind)%arr(:, :, :) + iv_safe(i)%val * cgl%cg%q(iv_safe(i)%ind)%arr(:, :, :)
         enddo
         cgl => cgl%nxt
      enddo

    end subroutine q_lin_comb

!>
!! \brief Compute the average value over a list and subtract it
!!
!! \details Typically it is used on a list of leaves or on a single level
!<

   subroutine subtract_average(this, iv)

      use constants,  only: GEO_XYZ, GEO_RPZ, I_ONE
      use dataio_pub, only: die
      use domain,     only: dom
      use grid_cont,  only: grid_container
      use mpi,        only: MPI_DOUBLE_PRECISION, MPI_SUM, MPI_IN_PLACE
      use mpisetup,   only: comm, mpi_err

      implicit none

      class(cg_list), intent(in) :: this !< list for which we want to subtract its average from
      integer,        intent(in) :: iv   !< index of variable in cg%q(:) which we want to have zero average

      real :: avg, vol
      integer :: i
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg

      avg = 0.
      vol = 0.
      cgl => this%first
      do while (associated(cgl))
         cg => cgl%cg
         select case (dom%geometry_type)
            case (GEO_XYZ)
               avg = avg + sum(cg%q(iv)%span(cg%ijkse)) * cg%dvol
            case (GEO_RPZ)
               do i = cg%is, cg%ie
                  avg = avg + sum(cg%q(iv)%arr(i, cg%js:cg%je, cg%ks:cg%ke)) * cg%dvol * cg%x(i)
               enddo
            case default
               call die("[gc_list:subtract_average] Unsupported geometry.")
         end select
         vol = vol + cg%vol
         cgl => cgl%nxt
      enddo
      call MPI_Allreduce(MPI_IN_PLACE, avg, I_ONE, MPI_DOUBLE_PRECISION, MPI_SUM, comm, mpi_err)
      call MPI_Allreduce(MPI_IN_PLACE, vol, I_ONE, MPI_DOUBLE_PRECISION, MPI_SUM, comm, mpi_err) !! \todo calculate this in some init routine
      avg = avg / vol

      call this%q_add_val(iv, -avg)

   end subroutine subtract_average

!>
!! \brief Calculate L2 norm
!!
!! \todo modify the code for reusing in subtract_average?
!<

   real function norm_sq(this, iv) result (norm)

      use constants,  only: GEO_XYZ, GEO_RPZ, I_ONE
      use dataio_pub, only: die
      use domain,     only: dom
      use grid_cont,  only: grid_container
      use mpi,        only: MPI_DOUBLE_PRECISION, MPI_SUM, MPI_IN_PLACE
      use mpisetup,   only: comm, mpi_err

      implicit none

      class(cg_list), intent(in) :: this !< list for which we want to calculate the L2 norm
      integer, intent(in)  :: iv   !< index of variable in cg%q(:) for which we want to find the norm

      integer :: i
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg

      norm = 0.
      cgl => this%first
      do while (associated(cgl))
         cg => cgl%cg
         select case (dom%geometry_type)
            case (GEO_XYZ)
               norm = norm + sum(cg%q(iv)%span(cg%ijkse)**2) * cg%dvol
            case (GEO_RPZ)
               do i = cg%is, cg%ie
                  norm = norm + sum(cg%q(iv)%arr(i, cg%js:cg%je, cg%ks:cg%ke)**2) * cg%dvol * cg%x(i)
               enddo
            case default
               call die("[gc_list:norm_sq] Unsupported geometry.")
         end select
         cgl => cgl%nxt
      enddo
      call MPI_Allreduce(MPI_IN_PLACE, norm, I_ONE, MPI_DOUBLE_PRECISION, MPI_SUM, comm, mpi_err)
      norm = sqrt(norm)

   end function norm_sq

!> \brief Clear boundary values

   subroutine zero_boundaries(this)

      implicit none

      class(cg_list), intent(inout) :: this  !< list for which clear the boundary values (typically a single level)

      type(cg_list_element), pointer :: cgl

      cgl => this%first
      do while (associated(cgl))
         cgl%cg%mg%bnd_x(:,:,:) = 0.
         cgl%cg%mg%bnd_y(:,:,:) = 0.
         cgl%cg%mg%bnd_z(:,:,:) = 0.
         cgl => cgl%nxt
      enddo

   end subroutine zero_boundaries

!> \brief Construct name of emergency ASCII dump

   subroutine numbered_ascii_dump(this, qlst, basename, a)

      use dataio_pub, only: halfstep, msg
      use global,     only: nstep, do_ascii_dump
      use mpisetup,   only: proc

      implicit none

      class(cg_list),        intent(inout) :: this     !< list for which do the dump (usually all_cg)
      integer, dimension(:), intent(in)    :: qlst     !< list of scalar fields to be printed
      character(len=*),      intent(in)    :: basename !< first part of the filename
      integer, optional,     intent(in)    :: a        !< additional number

      integer             :: l, n

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

      use constants,   only: LO
      use dataio_pub,  only: msg, printio
      use named_array, only: qna

      implicit none

      class(cg_list),        intent(inout) :: this     !< list for which do the dump (usually all_cg)
      character(len=*),      intent(in)    :: filename !< name to write the emergency dump (should be different on each process)
      integer, dimension(:), intent(in)    :: qlst     !< list of scalar fields to be printed

      integer, parameter :: fu=30
      integer            :: i, j, k, q
      type(cg_list_element), pointer :: cgl

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
                  write(fu, '(3i4,i6,3es20.11e3)', advance='no') [ i, j, k ] - cgl%cg%ijkse(:, LO) + cgl%cg%off(:), cgl%cg%level_id, cgl%cg%x(i), cgl%cg%y(j), cgl%cg%z(k)
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

      write(msg,'(3a)') "[gc_list:ascii_dump] Wrote dump '",filename,"'"
      call printio(msg)

   end subroutine ascii_dump

!>
!! \brief This routine pollutes selected array with an insane value dirtyH.
!!
!! \details If anything in the multigrid works by accident, through compiler-dependent initialization or unintentional relying on outdated values,
!! the insane value should pollute the solution in an easily visible way.
!<

   subroutine set_dirty(this, iv)

      use constants, only: dirtyH
      use global,    only: dirty_debug

      implicit none

      class(cg_list), intent(inout) :: this !< list for which we want to apply pollution
      integer,        intent(in)    :: iv   !< index of variable in cg%q(:) which we want to pollute

      if (.not. dirty_debug) return

      call this%set_q_value(iv, dirtyH)

   end subroutine set_dirty

!> \brief This routine checks for detectable traces of set_dirty calls.

   subroutine check_dirty(this, iv, label, expand)

      use constants,   only: dirtyL, LO
      use dataio_pub,  only: warn, msg, die
      use domain,      only: dom
      use global,      only: dirty_debug
      use mpisetup,    only: proc
      use named_array, only: qna

      implicit none

      class(cg_list),            intent(inout) :: this   !< level which we are checking
      integer,                   intent(in)    :: iv     !< index of variable in cg%q(:) which we want to pollute
      character(len=*),          intent(in)    :: label  !< label to indicate the origin of call
      integer(kind=4), optional, intent(in)    :: expand !< also check guardcells

      integer :: i, j, k, ng, cnt
      type(cg_list_element), pointer :: cgl

      if (.not. dirty_debug) return
      if (iv < lbound(qna%lst, dim=1) .or. iv > ubound(qna%lst, dim=1)) call die("[gc_list:check_dirty] Invalid variable index.")

      ng = 0
      if (present(expand)) ng = min(dom%nb, expand)

      cnt = 0
      cgl => this%first
      do while (associated(cgl))
         do k = cgl%cg%ks-ng*dom%D_z, cgl%cg%ke+ng*dom%D_z
            do j = cgl%cg%js-ng*dom%D_y, cgl%cg%je+ng*dom%D_y
               do i = cgl%cg%is-ng*dom%D_x, cgl%cg%ie+ng*dom%D_x
                  if (abs(cgl%cg%q(iv)%arr(i, j, k)) > dirtyL) then
                     ! if (count([i<cgl%cg%is .or. i>cgl%cg%ie, j<cgl%cg%js .or. j>cgl%cg%je, k<cgl%cg%ks .or. k>cgl%cg%ke]) <=1) then ! excludes corners
                     write(msg, '(3a,i4,a,i3,a,i5,3a,3i6,a,g20.12)') "[gc_list:check_dirty] ", trim(label), "@", proc, " lvl^", cgl%cg%level_id, &
                          &                                          " cg#", cgl%cg%grid_id, " '", trim(qna%lst(iv)%name), "'(", &
                          &                                          [ i, j, k ] - cgl%cg%ijkse(:, LO) + cgl%cg%off(:), ") = ", cgl%cg%q(iv)%arr(i, j, k)
                     call warn(msg)
                     cnt = cnt + 1
                     ! endif
                  endif
               enddo
            enddo
         enddo
         cgl => cgl%nxt
      enddo

      if (cnt /= 0) then
         write(msg,'(a,i8,a,i5)')"[gc_list:check_dirty] Found ", cnt, " dirty value @ process ", proc
         call die(msg)
      endif

   end subroutine check_dirty

!> \brief Check values of all named arrays for big_float

   subroutine check_for_dirt(this)

      use constants,   only: big_float
      use dataio_pub,  only: warn, msg
      use named_array, only: qna, wna

      implicit none

      class(cg_list), intent(in) :: this          !< object invoking type-bound procedure

      integer :: i
      type(cg_list_element), pointer :: cgl

      cgl => this%first
      do while (associated(cgl))
         do i = lbound(qna%lst(:), dim=1), ubound(qna%lst(:), dim=1)
            if (cgl%cg%q(i)%check()) then
               write(msg,'(3a,I12,a)') "[cg_list_global:check_for_dirt] Array ", trim(qna%lst(i)%name), " has ", &
                  & count(cgl%cg%q(i)%arr >= big_float), " wrong values."
               call warn(msg)
            endif
         enddo
         do i = lbound(wna%lst(:), dim=1), ubound(wna%lst(:), dim=1)
            if (cgl%cg%w(i)%check()) then
               write(msg,'(3a,I12,a)') "[cg_list_global:check_for_dirt] Array ", trim(wna%lst(i)%name), " has ", &
                  & count(cgl%cg%w(i)%arr >= big_float), " wrong values."
               call warn(msg)
            endif
         enddo
         cgl => cgl%nxt
      enddo

   end subroutine check_for_dirt

!> \brief Update mpisetup::req(:)

   subroutine update_req(this)

      use constants, only: INVALID, xdim, zdim
      use domain,    only: dom
      use mpisetup,  only: inflate_req

      implicit none

      class(cg_list), intent(in) :: this

      integer :: nrq, d
      type(cg_list_element), pointer :: cgl

      nrq = 0
      cgl => this%first
      do while (associated(cgl))

         do d = xdim, zdim
            if (allocated(cgl%cg%q_i_mbc(d, dom%nb)%mbc)) nrq = nrq + 2 * count(cgl%cg%q_i_mbc(d, dom%nb)%mbc(:) /= INVALID)
         enddo

         cgl => cgl%nxt
      enddo
      call inflate_req(nrq)

   end subroutine update_req

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
!!$      class(cg_list), intent(inout) :: this !< object invoking type-bound procedure
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

end module gc_list
