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

   use constants, only: dsetnamelen
   use grid_cont, only: grid_container

   implicit none

   private
   public :: cg_list, cg_list_element, ind_val, all_cg

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
!> \todo merge lists

   end type cg_list

   !>
   !! \brief Common properties of 3D and 4D named arrays
   !!
   !! \todo consider upgrading vital from logical to integer describing prolongation/restriction order or INVALID..
   !<
   type :: na_var
      character(len=dsetnamelen)                 :: name          !< a user-provided id for the array
      logical                                    :: vital         !< fields that are subject of automatic prolongation and restriction (e.g. state variables)
      integer(kind=4)                            :: restart_mode  !< AT_IGNORE: do not write to restart, AT_NO_B write without ext. boundaries, AT_OUT_B write with ext. boundaries
      integer(kind=4)                            :: ord_prolong   !< Prolongation order for the variable
      integer(kind=4), allocatable, dimension(:) :: position      !< VAR_CENTER by default, also possible VAR_CORNER and VAR_[XYZ]FACE
      integer(kind=4)                            :: dim4          !< <=0 for 3D arrays, >0 for 4D arrays
      logical                                    :: multigrid     !< .true. for variables that may exist below base level (e.g. work fields for multigrid solver)
   end type na_var

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
   !!
   !! Unfortunately we cannot move this type to a separate file because of dependencies in add_new
   !<
   type, extends(cg_list) :: cg_list_global
      type(na_var), dimension(:), allocatable :: q_lst !< information about registered 3D named arrays
      type(na_var), dimension(:), allocatable :: w_lst !< information about registered 4D named arrays
      !< indices of the most commonly used arrays stored in cg%w and cg%q
      integer :: fi                                    !< fluid           : cg%w(all_cg%fi)
      integer :: bi                                    !< magnetic field  : cg%w(all_cg%bi)
      integer :: wai                                   !< auxiliary array : cg%q(all_cg%wai)
    contains
      procedure :: reg_var        !< Add a variable (cg%q or cg%w) to all grid containers
      procedure :: check_na       !< Check if all named arrays are consistently registered
      procedure :: check_for_dirt !< Check all named arrays for constants:big_float
      procedure :: print_vars     !< write a summary on registered fields. Can be useful for debugging
      procedure :: ind
      procedure :: ind_4d
      procedure :: exists
      procedure :: exists_4d
   end type cg_list_global

   !> \brief Index - value pairs
   type ind_val
      integer :: ind  !< index in cg%q
      real    :: val  !< value for multiplication
   end type ind_val

   type(cg_list_global) :: all_cg                                     !< all grid containers; \todo restore protected

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
      integer :: i

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

      ! adding a cg to all_cg list implies registering all known named arrays for this cg
      select type(this)
         type is (cg_list_global)
            if (allocated(this%q_lst)) then
               do i = lbound(this%q_lst, dim=1), ubound(this%q_lst, dim=1)
                  call add_na(this%last%cg, this%q_lst(i)%multigrid)
               enddo
            endif
            if (allocated(this%w_lst)) then
               do i = lbound(this%w_lst, dim=1), ubound(this%w_lst, dim=1)
                  call add_na_4d(this%last%cg, this%w_lst(i)%dim4 )
               enddo
            endif
      end select

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
!! \brief Use this routine to add a variable (cg%q or cg%w) to all grid containers.
!!
!! \details Register a rank-3 array of given name in each grid container (in cg%q) and decide what to do with it on restart.
!! When dim4 is present then create a rank-4 array instead.(in cg%w)
!<

   subroutine reg_var(this, name, vital, restart_mode, ord_prolong, dim4, position, multigrid)

      use constants,   only: INVALID, VAR_CENTER, AT_NO_B, AT_IGNORE, O_INJ
      use dataio_pub,  only: die, warn, msg

      implicit none

      class(cg_list_global),                   intent(inout) :: this          !< object invoking type-bound procedure
      character(len=*),                        intent(in)    :: name          !< Name of the variable to be registered
      logical,                       optional, intent(in)    :: vital         !< .false. for arrays that don't need to be prolonged or restricted automatically
      integer(kind=4),               optional, intent(in)    :: restart_mode  !< Write to the restart if not AT_IGNORE. Several write modes can be supported.
      integer(kind=4),               optional, intent(in)    :: ord_prolong   !< Prolongation order for the variable
      integer(kind=4),               optional, intent(in)    :: dim4          !< If present then register the variable in the cg%w array.
      integer(kind=4), dimension(:), optional, intent(in)    :: position      !< If present then use this value instead of VAR_CENTER
      logical,                       optional, intent(in)    :: multigrid     !< If present and .true. then allocate cg%q(:)%arr and cg%w(:)%arr also below base level

      type(cg_list_element), pointer :: cgl
      logical :: mg, vit
      integer :: nvar
      integer(kind=4) :: op, d4, rm
      integer(kind=4), allocatable, dimension(:) :: pos

      vit = .false.
      if (present(vital)) vit = vital

      rm = AT_IGNORE
      if (present(restart_mode)) rm = restart_mode

      op = O_INJ
      if (present(ord_prolong)) op = ord_prolong

      mg = .false.
      if (present(multigrid)) mg = multigrid

      if (present(dim4)) then
         if (this%exists_4d(name)) then
            write(msg, '(3a)')"[gc_list:reg_var] A rank-4 array '",trim(name),"' was already registered."
            call die(msg)
         endif
         if (mg) call die("[gc_list:reg_var] there are no rank-4 multigrid arrays yet")
         d4 = dim4
         nvar = dim4
      else
         if (this%exists(name)) then
            write(msg, '(3a)')"[gc_list:reg_var] A rank-3 array '",trim(name),"' was already registered."
            call die(msg)
         endif
         d4 = int(INVALID, kind=4)
         nvar = 1
      endif

      if (allocated(pos)) call die("[gc_list:reg_var] pos(:) already allocated")
      allocate(pos(nvar))
      pos(:) = VAR_CENTER
      if (present(position)) then
         if (any(size(position) == [1, nvar])) then
            pos = position  !> \deprecated BEWARE: lhs reallocation
         else
            write(msg,'(2(a,i3))')"[gc_list:reg_var] position should be an array of 1 or ",nvar," values. Got ",size(position)
            call die(msg)
         endif
      endif
      if (any(pos(:) /= VAR_CENTER) .and. rm == AT_NO_B) then
         write(msg,'(3a)')"[gc_list:reg_var] no boundaries for restart with non cel-centered variable '",name,"' may result in loss of information in the restart files."
         call warn(msg)
      endif

      if (present(dim4)) then
         call add2lst(this%w_lst, name, vit, rm, op, pos, d4, mg)
      else
         call add2lst(this%q_lst, name, vit, rm, op, pos, d4, mg)
      endif

      cgl => this%first
      do while (associated(cgl))
         if (present(dim4)) then
            if (dim4<=0) then
               write(msg,'(3a)')"[gc_list:reg_var] dim4<=0 for variable'",name,"'"
               call die(msg)
            endif
            call add_na_4d(cgl%cg, dim4)
         else
            call add_na(cgl%cg, mg)
         endif
         cgl => cgl%nxt
      enddo

      deallocate(pos)

   end subroutine reg_var

!> \brief Add a named array properties to the list

   subroutine add2lst(lst, name, vital, restart_mode, ord_prolong, position, dim4, multigrid)

      implicit none

      type(na_var), dimension(:), allocatable, intent(inout) :: lst           !< the list to which we want to appent an entry
      character(len=*),                        intent(in)    :: name          !< Name of the variable to be registered
      logical,                                 intent(in)    :: vital         !< .false. for arrays that don't need to be prolonged or restricted automatically
      integer(kind=4),                         intent(in)    :: restart_mode  !< Write to the restart if not AT_IGNORE. Several write modes can be supported.
      integer(kind=4),                         intent(in)    :: ord_prolong   !< Prolongation order for the variable
      integer(kind=4), dimension(:),           intent(in)    :: position      !< VAR_CENTER in most cases, VAR_[XYZ]FACE for magnetic field
      integer(kind=4),                         intent(in)    :: dim4          !< If not INVALID then the variable is in the cg%w array.
      logical,                                 intent(in)    :: multigrid     !< If .true. then cg%q(:)%arr and cg%w(:)%arr are allocated also below base level

      type(na_var), dimension(:), allocatable :: tmp

      if (.not. allocated(lst)) then
         allocate(lst(1))
      else
         allocate(tmp(lbound(lst(:),dim=1):ubound(lst(:), dim=1) + 1))
         tmp(:ubound(lst(:), dim=1)) = lst(:)
         call move_alloc(from=tmp, to=lst)
      endif
      lst(ubound(lst(:), dim=1)) = na_var(name, vital, restart_mode, ord_prolong, position, dim4, multigrid)

   end subroutine add2lst

!> \brief Check if all named arrays are consistently registered

   subroutine check_na(this)

      use constants,  only: INVALID, base_level_id
      use dataio_pub, only: die

      implicit none

      class(cg_list_global), intent(in) :: this          !< object invoking type-bound procedure

      integer :: i
      type(cg_list_element), pointer :: cgl
      logical :: bad

      cgl => this%first
      do while (associated(cgl))
         if (allocated(this%q_lst) .neqv. allocated(cgl%cg%q)) then
            call die("[gc_list:check_na] allocated(all_cg%q_lst) .neqv. allocated(cgl%cg%q)")
         else if (allocated(this%q_lst)) then
            if (size(this%q_lst) /= size(cgl%cg%q)) then
               call die("[gc_list:check_na] size(all_cg%q_lst) /= size(cgl%cg%q)")
            else
               do i = lbound(this%q_lst, dim=1), ubound(this%q_lst, dim=1)
                  if (this%q_lst(i)%dim4 /= INVALID) call die("[gc_list:check_na] all_cg%q_lst(i) /= cgl%cg%q(i)")
                  if (associated(cgl%cg%q(i)%arr) .and. cgl%cg%level_id < base_level_id .and. .not. this%q_lst(i)%multigrid) &
                       call die("[gc_list:check_na] non-multigrid cgl%cg%q(i) allocated on coarse level")
               enddo
            endif
         endif
         if (allocated(this%w_lst) .neqv. allocated(cgl%cg%w)) then
            call die("[gc_list:check_na] allocated(all_cg%w_lst) .neqv. allocated(cgl%cg%w)")
         else if (allocated(this%w_lst)) then
            if (size(this%w_lst) /= size(cgl%cg%w)) then
               call die("[gc_list:check_na] size(all_cg%w_lst) /= size(cgl%cg%w)")
            else
               do i = lbound(this%w_lst, dim=1), ubound(this%w_lst, dim=1)
                  bad = .false.
                  if (associated(cgl%cg%w(i)%arr)) bad = this%w_lst(i)%dim4 /= size(cgl%cg%w(i)%arr, dim=1) .and. cgl%cg%level_id >= base_level_id
                  if (this%w_lst(i)%dim4 <= 0 .or. bad) call die("[gc_list:check_na] all_cg%w_lst(i) /= cgl%cg%w(i)")
                  if (associated(cgl%cg%w(i)%arr) .and. cgl%cg%level_id < base_level_id) call die("[gc_list:check_na] cgl%cg%w(i) allocated on coarse level")
               enddo
            endif
         endif
         cgl => cgl%nxt
      enddo

   end subroutine check_na

!> \brief Check values of all named arrays for big_float

   subroutine check_for_dirt(this)

      use constants,  only: big_float
      use dataio_pub, only: warn, msg

      implicit none

      class(cg_list_global), intent(in) :: this          !< object invoking type-bound procedure

      integer :: i
      type(cg_list_element), pointer :: cgl

      cgl => this%first
      do while (associated(cgl))
         do i = lbound(this%q_lst, dim=1), ubound(this%q_lst, dim=1)
            if (cgl%cg%q(i)%check()) then
               write(msg,'(3a,I12,a)') "[gc_list:check_for_dirt] Array ", trim(this%q_lst(i)%name), " has ", &
                  & count(cgl%cg%q(i)%arr >= big_float), " wrong values."
               call warn(msg)
            endif
         enddo
         do i = lbound(this%w_lst, dim=1), ubound(this%w_lst, dim=1)
            if (cgl%cg%w(i)%check()) then
               write(msg,'(3a,I12,a)') "[gc_list:check_for_dirt] Array ", trim(this%w_lst(i)%name), " has ", &
                  & count(cgl%cg%w(i)%arr >= big_float), " wrong values."
               call warn(msg)
            endif
         enddo
         cgl => cgl%nxt
      enddo

   end subroutine check_for_dirt
!>
!! \brief Register a new 3D entry in current cg with given name. Called from cg_list_global::reg_var
!!
!! \warning This routine should not be called directly from user routines
!<
   subroutine add_na(this, multigrid)

      use constants,   only: base_level_id
      use grid_cont,   only: grid_container
      use named_array, only: named_array3d

      implicit none

      type(grid_container), intent(inout) :: this
      logical,              intent(in)    :: multigrid     !< If .true. then cg%q(:)%arr and cg%w(:)%arr are allocated also below base level

      type(named_array3d), allocatable, dimension(:) :: tmp

      if (.not. allocated(this%q)) then
         allocate(this%q(1))
      else
         allocate(tmp(lbound(this%q(:),dim=1):ubound(this%q(:), dim=1) + 1))
         tmp(:ubound(this%q(:), dim=1)) = this%q(:)
         call move_alloc(from=tmp, to=this%q)
      endif

      if (multigrid .or. this%level_id >= base_level_id) call this%q(ubound(this%q(:), dim=1))%init(this%n_(:))

   end subroutine add_na

!>
!! \brief Register a new 4D entry in current cg with given name. Called from cg_list_global::reg_var
!!
!! \warning This routine should not be called directly from user routines
!! \deprecated Almost duplicated code with add_na
!<
   subroutine add_na_4d(this, n)

      use constants,   only: xdim, zdim, ndims, base_level_id
      use domain,      only: dom
      use grid_cont,   only: grid_container, set_mpi_types
      use named_array, only: named_array4d

      implicit none

      type(grid_container), intent(inout) :: this
      integer(kind=4), intent(in) :: n

      type(named_array4d), allocatable, dimension(:) :: tmp
      integer :: iw, d, ib, g

      if (.not. allocated(this%w)) then
         allocate(this%w(1))
      else
         allocate(tmp(lbound(this%w(:),dim=1):ubound(this%w(:), dim=1) + 1))
         tmp(:ubound(this%w(:), dim=1)) = this%w(:)
         do iw = lbound(this%w(:), dim=1), ubound(this%w(:), dim=1) ! prevent memory leak
            if (allocated(this%w(iw)%w_i_mbc)) deallocate(this%w(iw)%w_i_mbc)
            if (allocated(this%w(iw)%w_o_mbc)) deallocate(this%w(iw)%w_o_mbc)
         enddo
         call move_alloc(from=tmp, to=this%w)
      endif

      if (this%level_id >= base_level_id) then
         iw = ubound(this%w(:), dim=1)
         allocate(this%w(iw)%w_i_mbc(ndims, dom%nb), this%w(iw)%w_o_mbc(ndims, dom%nb))
         do d = xdim, zdim
            do ib = 1, dom%nb
               if (allocated(this%i_bnd(d, ib)%seg)) then
                  allocate(this%w(iw)%w_i_mbc(d, ib)%mbc(lbound(this%i_bnd(d, ib)%seg, dim=1):ubound(this%i_bnd(d, ib)%seg, dim=1)), &
                       &   this%w(iw)%w_o_mbc(d, ib)%mbc(lbound(this%i_bnd(d, ib)%seg, dim=1):ubound(this%i_bnd(d, ib)%seg, dim=1)))

                  do g = lbound(this%i_bnd(d, ib)%seg, dim=1), ubound(this%i_bnd(d, ib)%seg, dim=1)
                     call set_mpi_types([n, this%n_(:)], this%i_bnd(d, ib)%seg(g)%se(:,:), this%w(iw)%w_i_mbc(d, ib)%mbc(g))
                     call set_mpi_types([n, this%n_(:)], this%o_bnd(d, ib)%seg(g)%se(:,:), this%w(iw)%w_o_mbc(d, ib)%mbc(g))
                  enddo
               endif
            enddo
         enddo
         call this%w(iw)%init( [n, this%n_(:)] )
      endif

   end subroutine add_na_4d

!>
!! \brief Get the index of a named 3d array of given name.
!!
!! \warning OPT The indices aren't updated so cache them, whenever possible
!<
   function ind(this, name) result(rind)

      use dataio_pub, only: die, msg, warn

      implicit none

      class(cg_list_global), intent(inout) :: this
      character(len=*), intent(in) :: name

      integer :: rind, i

      rind = 0

      do i = lbound(this%q_lst, dim=1, kind=4), ubound(this%q_lst, dim=1, kind=4)
         if (trim(name) ==  this%q_lst(i)%name) then
            if (rind /= 0) then
               write(msg, '(2a)') "[gc_list:ind] multiple entries with the same name: ", trim(name)
               call die(msg)
            endif
            rind = i
         endif
      enddo

      if (rind == 0) then
         write(msg, '(2a)') "[gc_list:ind] requested entry not found: ", trim(name)
         call warn(msg)
      endif

   end function ind

!>
!! \brief Get the index of a named 4d array of given name.
!!
!! \details This method is provided for convenience only. Use ptr whenever possible.
!<
   function ind_4d(this, name) result(rind)

      use dataio_pub, only: die, msg, warn

      implicit none

      class(cg_list_global), intent(inout) :: this
      character(len=*), intent(in) :: name

      integer :: rind, i

      rind = 0

      do i = lbound(this%w_lst, dim=1, kind=4), ubound(this%w_lst, dim=1, kind=4)
         if (trim(name) ==  this%w_lst(i)%name) then
            if (rind /= 0) then
               write(msg, '(2a)') "[gc_list:ind_4d] multiple entries with the same name: ", trim(name)
               call die(msg)
            endif
            rind = i
         endif
      enddo

      if (rind == 0) then
         write(msg, '(2a)') "[gc_list:ind_4d] requested entry not found: ", trim(name)
         call warn(msg)
      endif

   end function ind_4d

!> \brief Check if a 3D array of given name is already registered

   function exists(this, name)

      use dataio_pub, only: die, msg

      implicit none

      class(cg_list_global), intent(inout) :: this
      character(len=*), intent(in) :: name

      logical :: exists
      integer :: i

      exists = .false.

      if (allocated(this%q_lst)) then
         do i = lbound(this%q_lst, dim=1), ubound(this%q_lst, dim=1)
            if (trim(name) ==  this%q_lst(i)%name) then
               if (exists) then
                  write(msg, '(2a)') "[gc_list:exists] multiple entries with the same name: ", trim(name)
                  call die(msg)
               endif
               exists = .true.
            endif
         enddo
      endif

   end function exists

!> \brief Check if a 4D array of given name is already registered

   function exists_4d(this, name) result(exists)

      use dataio_pub, only: die, msg

      implicit none

      class(cg_list_global), intent(inout) :: this
      character(len=*), intent(in) :: name

      logical :: exists
      integer :: i

      exists = .false.

      if (allocated(this%w_lst)) then
         do i = lbound(this%w_lst, dim=1), ubound(this%w_lst, dim=1)
            if (trim(name) ==  this%w_lst(i)%name) then
               if (exists) then
                  write(msg, '(2a)') "[gc_list:exists] multiple entries with the same name: ", trim(name)
                  call die(msg)
               endif
               exists = .true.
            endif
         enddo
      endif

   end function exists_4d

   subroutine print_vars(this)

      use dataio_pub, only: printinfo, msg
      use mpisetup,   only: slave

      implicit none

      class(cg_list_global), intent(inout) :: this

      integer :: i

      if (slave) return

      write(msg,'(a,i2,a)')"[gc_list:print_vars] Found ",size(this%q_lst)," rank-3 arrays:"
      call printinfo(msg)
      do i = lbound(this%q_lst(:), dim=1), ubound(this%q_lst(:), dim=1)
         write(msg,'(3a,l2,a,i2,a,l2,2(a,i2))')"'", this%q_lst(i)%name, "', vital=", this%q_lst(i)%vital, ", restart_mode=", this%q_lst(i)%restart_mode, &
              &                                ", multigrid=", this%q_lst(i)%multigrid, ", ord_prolong=", this%q_lst(i)%ord_prolong, ", position=", this%q_lst(i)%position(:)
         call printinfo(msg)
      enddo

      write(msg,'(a,i2,a)')"[gc_list:print_vars] Found ",size(this%w_lst)," rank-4 arrays:"
      call printinfo(msg)
      do i = lbound(this%w_lst(:), dim=1), ubound(this%w_lst(:), dim=1)
         write(msg,'(3a,l2,a,i2,a,l2,2(a,i2),a,100i2)')"'", this%w_lst(i)%name, "', vital=", this%w_lst(i)%vital, ", restart_mode=", this%w_lst(i)%restart_mode, &
              &                                        ", multigrid=", this%w_lst(i)%multigrid, ", ord_prolong=", this%w_lst(i)%ord_prolong, &
              &                                        ", components=", this%w_lst(i)%dim4, ", position=", this%w_lst(i)%position(:)
         call printinfo(msg)
      enddo

   end subroutine print_vars

!>
!! \brief Find munimum or maximum value over a specified list of grid containers
!!
!! \details It should be possible to find an extremum over a given level or leaf blocks or something
!<
   subroutine get_extremum(this, ind, minmax, prop, dir)

      use constants,  only: MINL, MAXL, I_ONE, ndims, xdim, ydim, zdim
      use dataio_pub, only: msg, warn, die
      use domain,     only: dom
      use grid_cont,  only: grid_container
      use mpi,        only: MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_STATUS_IGNORE, MPI_2DOUBLE_PRECISION, MPI_MINLOC, MPI_MAXLOC, MPI_IN_PLACE
      use mpisetup,   only: comm, ierr, master, proc, FIRST
      use types,      only: value

      implicit none

      class(cg_list),            intent(in)  :: this    !< object invoking type-bound procedure
      integer,                   intent(in)  :: ind     !< Index in cg%q(:)
      integer(kind=4),           intent(in)  :: minmax  !< minimum or maximum ?
      type(value),               intent(out) :: prop    !< precise location of the extremum to be found
      integer(kind=4), optional, intent(in)  :: dir     !< order the cell size in dir direction


      type(grid_container),   pointer :: cg
      type(cg_list_element),  pointer :: cgl
      real, dimension(:,:,:), pointer :: tab
      integer,                       parameter :: tag1 = 11, tag2 = tag1 + 1, tag3 = tag2 + 1
      integer, dimension(MINL:MAXL), parameter :: op = [ MPI_MINLOC, MPI_MAXLOC ]
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

         tab => cg%q(ind)%span(cg%ijkse)
         select case (minmax)
            case (MINL)
               if (minval(tab) < prop%val) then
                  prop%val = minval(tab)
                  prop%loc = minloc(tab) + dom%nb
               endif
            case (MAXL)
               if (maxval(tab) > prop%val) then
                  prop%val = maxval(tab)
                  prop%loc = maxloc(tab) + dom%nb
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
         if (present(dir))      prop%assoc        = cg%dl(dir)
      endif

      if (prop%proc /= 0) then
         if (proc == prop%proc) then ! slave
            call MPI_Send (prop%loc,    ndims, MPI_INTEGER,          FIRST, tag1, comm, ierr)
            call MPI_Send (prop%coords, ndims, MPI_DOUBLE_PRECISION, FIRST, tag2, comm, ierr)
            if (present(dir)) call MPI_Send (prop%assoc, I_ONE, MPI_DOUBLE_PRECISION, FIRST, tag3, comm, ierr)
         endif
         if (master) then
            call MPI_Recv (prop%loc,    ndims, MPI_INTEGER,          prop%proc, tag1, comm, MPI_STATUS_IGNORE, ierr)
            call MPI_Recv (prop%coords, ndims, MPI_DOUBLE_PRECISION, prop%proc, tag2, comm, MPI_STATUS_IGNORE, ierr)
            if (present(dir)) call MPI_Recv (prop%assoc, I_ONE, MPI_DOUBLE_PRECISION, prop%proc, tag3, comm, MPI_STATUS_IGNORE, ierr)
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
      use mpisetup,   only: comm, ierr

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
      call MPI_Allreduce(MPI_IN_PLACE, avg, I_ONE, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, vol, I_ONE, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr) !! \todo calculate this in some init routine
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
      use mpisetup,   only: comm, ierr

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
      call MPI_Allreduce(MPI_IN_PLACE, norm, I_ONE, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
      norm = sqrt(norm)

   end function norm_sq

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
