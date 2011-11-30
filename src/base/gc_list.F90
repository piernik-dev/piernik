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

   use constants, only: ndims
   use grid_cont, only: grid_container

   implicit none

   private
   public :: cg_list_global, cg_list_level, cg_list_patch, cg_list, cg_list_element

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

!      procedure :: init_el
      procedure :: init_new                          !< A constructor for an empty list
      generic, public :: init => init_new!, init_el  !< All constructors of a new list

      procedure :: add_new                           !< Add new element to the list
      generic, public :: add => add_new              !< All methods of adding a new element

      procedure :: del_lnk                           !< Destroy the element
      procedure :: del_lst                           !< Destroy the list
      generic, public :: delete => del_lnk, del_lst  !< All methods of destroying

      procedure :: un_link                           !< Un-link the element
      procedure :: get_extremum                      !< Find munimum or maximum value over a s list
      procedure :: print_list                        !< Print the list and associated cg ID
!> \todo merge lists

   end type cg_list

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
   type, extends(cg_list) :: cg_list_global
   !> \todo store the information about registered entries also here (name, restart_mode, index, dim4)
    contains
      procedure :: reg_var        !< Add a variable (cg%q or cg%w) to all grid containers
   end type cg_list_global

   !>
   !! \brief A list of all cg of the same resolution.
   !!
   !! \details For positive refinement levels the list may be composed of several disconnected subsets of cg ("islands: made of one or more cg: cg_list_patch).
   !<
   type, extends(cg_list) :: cg_list_level
      integer :: lev                           !< level number (relative to base level). For printing, debug, and I/O use only. No arithmetic should depend on it.
      type(cg_list_level), pointer :: coarser  !< coarser level cg set or null()
      type(cg_list_level), pointer :: finer    !< finer level cg set or null()
    contains
!      procedure :: prolong                    !< interpolate the grid data to this%finer level
!      procedure :: restrict                   !< interpolate the grid data from this%coarser level
! fine-coarse boundary exchanges may also belong to this type
   end type cg_list_level

   !>
   !! \brief A list of grid containers that cover single box (or rectangle) on a certain resolution level
   !!
   !! \details This set would be a result of base domain or patch decomposition
   !<
   type, extends(cg_list) :: cg_list_patch
      integer(kind=4), dimension(ndims) :: n_d                !< number of grid cells
      integer(kind=8), dimension(ndims) :: off                !< offset (with respect to the base level, counted on own level)
      type(cg_list_patch), pointer :: parent                  !< Parent patch (or null()). \todo Consider relaxing this restriction and allow multi-parent patches
      type(cg_list_patch), dimension(:), pointer :: children  !< refined patches
      !> \todo consider creating neigbour list (or (ndims, LO:HI) lists)
   end type cg_list_patch

contains

!> \brief a constructor for an empty list
   subroutine init_new(this)

      implicit none

      class(cg_list), intent(out) :: this !< object invoking type-bound procedure

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

      !> \todo update the lists of registered variables

   end subroutine add_new

!> \brief destroy the element
   subroutine del_lnk(this, cgle)

      use constants,  only: INVALID
      use dataio_pub, only: warn

      implicit none

      class(cg_list), intent(inout) :: this !< object invoking type-bound procedure
      type(cg_list_element), pointer, intent(inout) :: cgle !< the element to be eradicated

      if (.not. associated(cgle)) then
         call warn("[gc_list:delete] tried to remove null() element")
         return
      endif

      call this%un_link(cgle)
      if (associated(cgle%cg)) then
         if (cgle%cg%grid_n > INVALID) then ! dirty trick
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

!> \brief remove the element from the list, keep ith content
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

      use dataio_pub, only: warn, msg

      implicit none

      class(cg_list), intent(inout) :: this !< object invoking type-bound procedure

      type(cg_list_element), pointer :: cur
      integer :: cnt

      if (.not. associated(this%first)) then
         write(msg,'(a)')"[gc_list:print_list] Empty list"
         call warn(msg)
         if (this%cnt /= 0) then
            write(msg,'(a,i6)')"[gc_list:print_list] Empty list length /= 0 : ",this%cnt
            call warn(msg)
         endif
         if (associated(this%last)) then
            call warn("[gc_list:print_list] Tail without head")
            if (associated(this%last%cg)) then
               write(msg,'(a,i7)')"Last element #",this%last%cg%grid_n
               call warn(msg)
            endif
         endif
         return
      endif

      if (.not. associated(this%last)) call warn("[gc_list:print_list] Head without tail")

      if (associated(this%first%cg)) then
         write(msg,'(a,i7)')"First element #",this%first%cg%grid_n
         call warn(msg)
      endif
      if (associated(this%last%cg)) then
         write(msg,'(a,i7)')"Last element #",this%last%cg%grid_n
         call warn(msg)
      endif

      cnt = 0
      cur => this%first
      do while (associated(cur))
         cnt = cnt + 1
         if (associated(cur%cg)) then
            write(msg,'(i5,a,i7)')cnt,"-th element #",cur%cg%grid_n
            call warn(msg)
         endif
         cur => cur%nxt
      enddo

      if (cnt /= this%cnt) then
         write(msg, '(2(a,i5))')"[gc_list:print_list] this%cnt = ",this%cnt," /= ",cnt
         call warn(msg)
      endif

      cur => this%last
      do while (associated(cur))
         cnt = cnt - 1
         if (associated(cur%cg)) then
            write(msg,'(i5,a,i7)')cnt,"-th element #",cur%cg%grid_n
            call warn(msg)
         endif
         cur => cur%prv
      enddo
    end subroutine print_list

!>
!! \brief Use this routine to add a variable (cg%q or cg%w) to all grid containers.
!!
!! \details Register a rank-3 array of given name in each grid container (in cg%q) and decide what to do with it on restart.
!! When dim4 is present then create a rank-4 array instead.(in cg%w)
!<

   subroutine reg_var(this, name, restart_mode, dim4)

      use dataio_pub, only: die, msg

      implicit none

      class(cg_list_global),     intent(in) :: this !< object invoking type-bound procedure
      character(len=*),          intent(in) :: name          !< Name of the variable to be registered
      integer(kind=4),           intent(in) :: restart_mode  !< Write to the restar if not AT_IGNORE. Several write modes can be supported.
      integer(kind=4), optional, intent(in) :: dim4          !< If present then register the variable in the cg%w array.

      type(cg_list_element), pointer :: cgl

      cgl => this%first
      do while (associated(cgl))
         if (present(dim4)) then
            if (dim4<=0) then
               write(msg,'(3a)')"[gc_list:reg_var] dim4<=0 for variable'",name,"'"
               call die(msg)
            endif
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

      class(cg_list),  intent(in)  :: this !< object invoking type-bound procedure
      integer(kind=4), intent(in)  :: ind     !< Index in cg%q(:)
      integer(kind=4), intent(in)  :: minmax  !< minimum or maximum ?
      type(value),     intent(out) :: prop    !< precise location of the extremum to be found


      type(grid_container), pointer :: cg
      type(cg_list_element), pointer :: cgl
      integer, parameter :: tag1 = 11, tag2 = tag1 + 1
      integer, dimension(MINL:MAXL), parameter :: op = [ MPI_MINLOC, MPI_MAXLOC ]
      real, dimension(:,:,:), pointer :: tab
      enum, bind(C)
         enumerator :: I_V, I_P !< value and proc
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
