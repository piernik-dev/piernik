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

!> \brief This module contains list of all grid containers and related methods

module cg_list_global

   use constants, only: dsetnamelen
   use gc_list,   only: cg_list

   implicit none

   private
   public :: all_cg

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

   !>
   !! \brief A list of grid containers that are supposed to have the same variables registered
   !!
   !! \details The main purpose of this type is to provide a type for the set of all grid containers with methods and properties
   !! that should not be available for any arbitrarily composed subset of grid containers. Typically there will be only one variable
   !! of this type available in the code: all_cg.
   !!
   !! It should be possible to use this type for more fancy, multi-domain  grid configurations such as:
   !! - Yin-Yang grid (covering sphere with two domains shaped as parts of the sphere to avoid polar singularities),
   !! - Cylindrical grid with cartesian core covering singularity at the axis.
   !! - Simulations with mixed dimensionality (e.g. 2d grid for dust particles and 3d grid for gas) should probably also use separate cg_list
   !! for their data (and additional routine for coupling the two grid sets).
   !<
   type, extends(cg_list) :: cg_list_glob

      type(na_var), dimension(:), allocatable :: q_lst !< information about registered 3D named arrays
      type(na_var), dimension(:), allocatable :: w_lst !< information about registered 4D named arrays
      integer(kind=4) :: ord_prolong_nb                !< Maximum number of boundary cells required for prolongation

      !< indices of the most commonly used arrays stored in cg%w and cg%q
      integer :: fi                                    !< fluid           : cg%w(all_cg%fi)
      integer :: bi                                    !< magnetic field  : cg%w(all_cg%bi)
      integer :: wai                                   !< auxiliary array : cg%q(all_cg%wai)

    contains
      procedure :: reg_var        !< Add a variable (cg%q or cg%w) to all grid containers
      procedure :: check_na       !< Check if all named arrays are consistently registered
      procedure :: check_for_dirt !< Check all named arrays for constants:big_float
      procedure :: print_vars     !< Write a summary on registered fields. Can be useful for debugging
      procedure :: ind            !< Get the index of a named 3d array of given name.
      procedure :: ind_4d         !< Get the index of a named 4d array of given name.
      procedure :: exists         !< Check if a 3D array of given name is already registered
      procedure :: exists_4d      !< Check if a 4D array of given name is already registered
      procedure :: update_req     !< Update mpisetup::req(:)
   end type cg_list_glob

   type(cg_list_glob) :: all_cg   !< all grid containers; \todo restore protected

contains

!>
!! \brief Use this routine to add a variable (cg%q or cg%w) to all grid containers.
!!
!! \details Register a rank-3 array of given name in each grid container (in cg%q) and decide what to do with it on restart.
!! When dim4 is present then create a rank-4 array instead.(in cg%w)
!<

   subroutine reg_var(this, name, vital, restart_mode, ord_prolong, dim4, position, multigrid)

      use constants,  only: INVALID, VAR_CENTER, AT_NO_B, AT_IGNORE, O_INJ
      use dataio_pub, only: die, warn, msg
      use gc_list,    only: cg_list_element

      implicit none

      class(cg_list_glob),                     intent(inout) :: this          !< object invoking type-bound procedure
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
            write(msg, '(3a)')"[cg_list_global:reg_var] A rank-4 array '",trim(name),"' was already registered."
            call die(msg)
         endif
         if (mg) call die("[cg_list_global:reg_var] there are no rank-4 multigrid arrays yet")
         d4 = dim4
         nvar = dim4
      else
         if (this%exists(name)) then
            write(msg, '(3a)')"[cg_list_global:reg_var] A rank-3 array '",trim(name),"' was already registered."
            call die(msg)
         endif
         d4 = int(INVALID, kind=4)
         nvar = 1
      endif

      if (allocated(pos)) call die("[cg_list_global:reg_var] pos(:) already allocated")
      allocate(pos(nvar))
      pos(:) = VAR_CENTER
      if (present(position)) then
         if (any(size(position) == [1, nvar])) then
            pos = position  !> \deprecated BEWARE: lhs reallocation
         else
            write(msg,'(2(a,i3))')"[cg_list_global:reg_var] position should be an array of 1 or ",nvar," values. Got ",size(position)
            call die(msg)
         endif
      endif
      if (any(pos(:) /= VAR_CENTER) .and. rm == AT_NO_B) then
         write(msg,'(3a)')"[cg_list_global:reg_var] no boundaries for restart with non cel-centered variable '",name,"' may result in loss of information in the restart files."
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
               write(msg,'(3a)')"[cg_list_global:reg_var] dim4<=0 for variable'",name,"'"
               call die(msg)
            endif
            call cgl%cg%add_na_4d(dim4)
         else
            call cgl%cg%add_na(mg)
         endif
         cgl => cgl%nxt
      enddo

      deallocate(pos)

   end subroutine reg_var

!> \brief Add a named array properties to the list

   subroutine add2lst(lst, name, vital, restart_mode, ord_prolong, position, dim4, multigrid)

      use constants,  only: I_ZERO, I_ONE, I_TWO, O_INJ, O_LIN, O_I2, O_D2, O_I3, O_I4, O_D3, O_D4
      use dataio_pub, only: die
      use domain,     only: dom

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
         !! \warning slight memory leak here, e.g. in use at exit: 572 bytes in 115 blocks, perhaps associated with na_var%position
      endif
      lst(ubound(lst(:), dim=1)) = na_var(name, vital, restart_mode, ord_prolong, position, dim4, multigrid)

      select case (ord_prolong)
         case (O_INJ)
            all_cg%ord_prolong_nb = max(all_cg%ord_prolong_nb, I_ZERO)
         case (O_LIN, O_I2, O_D2)
            all_cg%ord_prolong_nb = max(all_cg%ord_prolong_nb, I_ONE)
         case (O_I3, O_I4, O_D3, O_D4)
            all_cg%ord_prolong_nb = max(all_cg%ord_prolong_nb, I_TWO)
         case default
            call die("[cg_list_global:add2lst] Unknown prolongation order")
      end select

      if (all_cg%ord_prolong_nb > dom%nb) call die("[cg_list_global:add2lst] Insufficient number of guardcells for requested prolongation stencil")

   end subroutine add2lst

!> \brief Check if all named arrays are consistently registered

   subroutine check_na(this)

      use constants,  only: INVALID, base_level_id
      use dataio_pub, only: msg, die
      use gc_list,    only: cg_list_element

      implicit none

      class(cg_list_glob), intent(in) :: this          !< object invoking type-bound procedure

      integer :: i
      type(cg_list_element), pointer :: cgl
      logical :: bad

      cgl => this%first
      do while (associated(cgl))
         if (allocated(this%q_lst) .neqv. allocated(cgl%cg%q)) then
            write(msg,'(2(a,l2))')"[cg_list_global:check_na] allocated(all_cg%q_lst) .neqv. allocated(cgl%cg%q):",allocated(this%q_lst)," .neqv. ",allocated(cgl%cg%q)
            call die(msg)
         else if (allocated(this%q_lst)) then
            if (size(this%q_lst) /= size(cgl%cg%q)) then
               write(msg,'(2(a,i5))')"[cg_list_global:check_na] size(all_cg%q_lst) /= size(cgl%cg%q)",size(this%q_lst)," /= ",size(cgl%cg%q)
               call die(msg)
            else
               do i = lbound(this%q_lst, dim=1), ubound(this%q_lst, dim=1)
                  if (this%q_lst(i)%dim4 /= INVALID) then
                     write(msg,'(3a,i10)')"[cg_list_global:check_na] all_cg%q_lst(",i,"), named '",this%q_lst(i)%name,"' has dim4 set to ",this%q_lst(i)%dim4
                     call die(msg)
                  endif
                  if (associated(cgl%cg%q(i)%arr) .and. cgl%cg%level_id < base_level_id .and. .not. this%q_lst(i)%multigrid) then
                     write(msg,'(a,i3,3a)')"[cg_list_global:check_na] non-multigrid cgl%cg%q(",i,"), named '",this%q_lst(i)%name,"' allocated on coarse level"
                     call die(msg)
                  endif
               enddo
            endif
         endif
         if (allocated(this%w_lst) .neqv. allocated(cgl%cg%w)) then
            write(msg,'(2(a,l2))')"[cg_list_global:check_na] allocated(all_cg%w_lst) .neqv. allocated(cgl%cg%w)",allocated(this%w_lst)," .neqv. ",allocated(cgl%cg%w)
            call die(msg)
         else if (allocated(this%w_lst)) then
            if (size(this%w_lst) /= size(cgl%cg%w)) then
               write(msg,'(2(a,i5))')"[cg_list_global:check_na] size(all_cg%w_lst) /= size(cgl%cg%w)",size(this%w_lst)," /= ",size(cgl%cg%w)
               call die(msg)
            else
               do i = lbound(this%w_lst, dim=1), ubound(this%w_lst, dim=1)
                  bad = .false.
                  if (associated(cgl%cg%w(i)%arr)) bad = this%w_lst(i)%dim4 /= size(cgl%cg%w(i)%arr, dim=1) .and. cgl%cg%level_id >= base_level_id
                  if (this%w_lst(i)%dim4 <= 0 .or. bad) then
                     write(msg,'(a,i3,2a,2(a,i7))')"[cg_list_global:check_na] all_cg%w_lst(",i,") named '",this%w_lst(i)%name,"' has inconsistent dim4: ",&
                          &         this%w_lst(i)%dim4," /= ",size(cgl%cg%w(i)%arr, dim=1)
                     call die(msg)
                  endif
                  if (associated(cgl%cg%w(i)%arr) .and. cgl%cg%level_id < base_level_id) then
                     write(msg,'(a,i3,3a)')"[cg_list_global:check_na] cgl%cg%w(",i,"), named '",this%w_lst(i)%name,"' allocated on coarse level"
                     call die(msg)
                  endif
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
      use gc_list,    only: cg_list_element

      implicit none

      class(cg_list_glob), intent(in) :: this          !< object invoking type-bound procedure

      integer :: i
      type(cg_list_element), pointer :: cgl

      cgl => this%first
      do while (associated(cgl))
         do i = lbound(this%q_lst, dim=1), ubound(this%q_lst, dim=1)
            if (cgl%cg%q(i)%check()) then
               write(msg,'(3a,I12,a)') "[cg_list_global:check_for_dirt] Array ", trim(this%q_lst(i)%name), " has ", &
                  & count(cgl%cg%q(i)%arr >= big_float), " wrong values."
               call warn(msg)
            endif
         enddo
         do i = lbound(this%w_lst, dim=1), ubound(this%w_lst, dim=1)
            if (cgl%cg%w(i)%check()) then
               write(msg,'(3a,I12,a)') "[cg_list_global:check_for_dirt] Array ", trim(this%w_lst(i)%name), " has ", &
                  & count(cgl%cg%w(i)%arr >= big_float), " wrong values."
               call warn(msg)
            endif
         enddo
         cgl => cgl%nxt
      enddo

   end subroutine check_for_dirt

!>
!! \brief Get the index of a named 3d array of given name.
!!
!! \warning OPT The indices aren't updated so cache them, whenever possible
!<
   function ind(this, name) result(rind)

      use dataio_pub, only: die, msg, warn

      implicit none

      class(cg_list_glob), intent(inout) :: this
      character(len=*),    intent(in)    :: name

      integer :: rind, i

      rind = 0

      do i = lbound(this%q_lst, dim=1, kind=4), ubound(this%q_lst, dim=1, kind=4)
         if (trim(name) ==  this%q_lst(i)%name) then
            if (rind /= 0) then
               write(msg, '(2a)') "[cg_list_global:ind] multiple entries with the same name: ", trim(name)
               call die(msg)
            endif
            rind = i
         endif
      enddo

      if (rind == 0) then
         write(msg, '(2a)') "[cg_list_global:ind] requested entry not found: ", trim(name)
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

      class(cg_list_glob), intent(inout) :: this
      character(len=*),    intent(in)    :: name

      integer :: rind, i

      rind = 0

      do i = lbound(this%w_lst, dim=1, kind=4), ubound(this%w_lst, dim=1, kind=4)
         if (trim(name) ==  this%w_lst(i)%name) then
            if (rind /= 0) then
               write(msg, '(2a)') "[cg_list_global:ind_4d] multiple entries with the same name: ", trim(name)
               call die(msg)
            endif
            rind = i
         endif
      enddo

      if (rind == 0) then
         write(msg, '(2a)') "[cg_list_global:ind_4d] requested entry not found: ", trim(name)
         call warn(msg)
      endif

   end function ind_4d

!> \brief Check if a 3D array of given name is already registered

   function exists(this, name)

      use dataio_pub, only: die, msg

      implicit none

      class(cg_list_glob), intent(inout) :: this
      character(len=*),    intent(in)    :: name

      logical :: exists
      integer :: i

      exists = .false.

      if (allocated(this%q_lst)) then
         do i = lbound(this%q_lst, dim=1), ubound(this%q_lst, dim=1)
            if (trim(name) ==  this%q_lst(i)%name) then
               if (exists) then
                  write(msg, '(2a)') "[cg_list_global:exists] multiple entries with the same name: ", trim(name)
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

      class(cg_list_glob), intent(inout) :: this
      character(len=*),    intent(in)    :: name

      logical :: exists
      integer :: i

      exists = .false.

      if (allocated(this%w_lst)) then
         do i = lbound(this%w_lst, dim=1), ubound(this%w_lst, dim=1)
            if (trim(name) ==  this%w_lst(i)%name) then
               if (exists) then
                  write(msg, '(2a)') "[cg_list_global:exists] multiple entries with the same name: ", trim(name)
                  call die(msg)
               endif
               exists = .true.
            endif
         enddo
      endif

   end function exists_4d

!> \brief Update mpisetup::req(:)

   subroutine update_req(this)

      use constants, only: INVALID, xdim, zdim
      use domain,    only: dom
      use gc_list,   only: cg_list_element
      use mpisetup,  only: inflate_req

      implicit none

      class(cg_list_glob), intent(inout) :: this

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

!> \brief Summarize all registered fields and their properties

   subroutine print_vars(this)

      use dataio_pub, only: printinfo, msg
      use mpisetup,   only: slave

      implicit none

      class(cg_list_glob), intent(inout) :: this

      integer :: i

      if (slave) return

      write(msg,'(a,i2,a)')"[cg_list_global:print_vars] Found ",size(this%q_lst)," rank-3 arrays:"
      call printinfo(msg)
      do i = lbound(this%q_lst(:), dim=1), ubound(this%q_lst(:), dim=1)
         write(msg,'(3a,l2,a,i2,a,l2,2(a,i2))')"'", this%q_lst(i)%name, "', vital=", this%q_lst(i)%vital, ", restart_mode=", this%q_lst(i)%restart_mode, &
              &                                ", multigrid=", this%q_lst(i)%multigrid, ", ord_prolong=", this%q_lst(i)%ord_prolong, ", position=", this%q_lst(i)%position(:)
         call printinfo(msg)
      enddo

      write(msg,'(a,i2,a)')"[cg_list_global:print_vars] Found ",size(this%w_lst)," rank-4 arrays:"
      call printinfo(msg)
      do i = lbound(this%w_lst(:), dim=1), ubound(this%w_lst(:), dim=1)
         write(msg,'(3a,l2,a,i2,a,l2,2(a,i2),a,100i2)')"'", this%w_lst(i)%name, "', vital=", this%w_lst(i)%vital, ", restart_mode=", this%w_lst(i)%restart_mode, &
              &                                        ", multigrid=", this%w_lst(i)%multigrid, ", ord_prolong=", this%w_lst(i)%ord_prolong, &
              &                                        ", components=", this%w_lst(i)%dim4, ", position=", this%w_lst(i)%position(:)
         call printinfo(msg)
      enddo

   end subroutine print_vars

end module cg_list_global
