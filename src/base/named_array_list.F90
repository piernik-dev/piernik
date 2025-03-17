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
!! \brief Lists of 3D and 4D named arrays and associated routines
!!
!! \details The qna(:) and wna(:) lists contain only descriptions of registers named arrays.
!! Named arrays themselves are defined in "named_array" module.
!!
!! \todo Add array of names attributed to 4D array components and make it available for data dumps
!! \todo Split this type into generic part, 3D part and 4D part.
!!       Name of the 3D named arrays should be scalar, for 4D arrays should be 1D vectors.
!<

module named_array_list

   use constants, only: dsetnamelen, INVALID

   implicit none

   private
   public :: na_var_4d, na_var, na_var_list, na_var_list_q, na_var_list_w, qna, wna

   !> \brief Common properties of 3D and 4D named arrays
   type, abstract :: na_var_base
      character(len=dsetnamelen)                 :: name          !< a user-provided id for the array
      logical                                    :: vital         !< fields that are subject of automatic prolongation and restriction (e.g. state variables)
      integer(kind=4)                            :: restart_mode  !< AT_BACKUP, AT_IGNORE: do not write to restart
                                                                  !< AT_NO_B write without ext. boundaries
                                                                  !< AT_OUT_B write with ext. boundaries
      integer(kind=4)                            :: ord_prolong   !< Prolongation order for the variable
      logical                                    :: multigrid     !< .true. for variables that may exist below base level (e.g. work fields for multigrid solver)
   end type na_var_base

   type, extends(na_var_base) :: na_var
   private
      type(na_var), allocatable, dimension(:) :: tmp  !< temporary array for list extension
   contains
      procedure :: copy_var
      generic :: assignment(=) => copy_var
   end type na_var

   type, extends(na_var_base) :: na_var_4d
      integer(kind=4)                                     :: dim4  !< should be >0 for 4D arrays
      type(na_var_4d), allocatable, dimension(:), private :: tmp   !< temporary array for list extension
   contains
      procedure :: copy_var_4d
      generic :: assignment(=) => copy_var_4d
   end type na_var_4d

   !> \brief the generic list of named arrays with supporting routines
   type, abstract :: na_var_list
      class(na_var_base), allocatable, dimension(:) :: lst
   contains
      procedure, private ::  find_ind                            !< Get the index of a named array of given name. Don't die when can't find requested field.
      procedure :: ind                                           !< Get the index of a named array of given name.
      procedure :: exists                                        !< Check if a named array of given name is already registered
      procedure :: print_vars                                    !< Write a summary on registered fields. Can be useful for debugging
      procedure :: get_reslst                                    !< List of fields that are stored in the restart file
   end type na_var_list

   !> \brief the most commonly used 3D named array is wa, thus we add a shortcut here
   type, extends(na_var_list) :: na_var_list_q
      integer(kind=4) :: wai = INVALID                           !< auxiliary array : cg%q(qna%wai)
   contains
      procedure :: add2lst => add2lst_q                          !< Add a 3D array to the list
   end type na_var_list_q

   !> \brief the most commonly used 4D named arrays are u and b, thus we add shortcuts here
   type, extends(na_var_list) :: na_var_list_w
      integer(kind=4) :: fi = INVALID                            !< fluid           : cg%w(wna%fi)
      integer(kind=4) :: bi = INVALID                            !< magnetic field  : cg%w(wna%bi)
   contains
      procedure :: add2lst => add2lst_w                          !< Add a 4D array to the list
      procedure :: get_dim4                                      !< Get dim4 value for given array index
   end type na_var_list_w

   type(na_var_list_q), target :: qna !< list of registered 3D named arrays
   type(na_var_list_w), target :: wna !< list of registered 4D named arrays

contains

!>
!! \brief Get the index of a named array of given name.
!!
!! \warning OPT The indices aren't updated so cache them, whenever possible
!<
   function find_ind(this, name) result(rind)

      use constants,  only: INVALID
      use dataio_pub, only: die, msg

      implicit none

      class(na_var_list), intent(in) :: this
      character(len=*),   intent(in) :: name

      integer(kind=4) :: rind, i

      rind = INVALID
      if (allocated(this%lst)) then
         do i = lbound(this%lst, dim=1, kind=4), ubound(this%lst, dim=1, kind=4)
            if (trim(name) == this%lst(i)%name) then
               if (rind /= INVALID) then
                  write(msg, '(3a)') "[named_array_list:find_ind] multiple entries with the same name '", trim(name),"'"
                  call die(msg)
               endif
               rind = i
            endif
         enddo
      endif

   end function find_ind

!>
!! \brief Get the index of a named array of given name.
!!
!! \warning OPT The indices aren't updated so cache them, whenever possible
!<
   function ind(this, name) result(rind)

      use constants,  only: INVALID
      use dataio_pub, only: die, msg

      implicit none

      class(na_var_list), intent(in) :: this
      character(len=*),   intent(in) :: name

      integer(kind=4) :: rind

      rind = this%find_ind(name)

      if (rind == INVALID) then
         write(msg, '(3a)') "[named_array_list:ind] requested entry not found: '", trim(name), "'"
         call die(msg)
      endif

   end function ind

!> \brief Check if a named array of given name is already registered

   logical function exists(this, name)

      use constants, only: INVALID

      implicit none

      class(na_var_list), intent(in) :: this
      character(len=*),   intent(in) :: name

      exists = (this%find_ind(name) /= INVALID)

   end function exists

!> \brief Add a named array properties to the list

   subroutine add2lst_q(this, element)
      use constants,  only: wa_n
      use dataio_pub, only: die, msg

      implicit none

      class(na_var_list_q), intent(inout) :: this
      type(na_var),         intent(in)    :: element

      if (this%exists(element%name)) then
         write(msg, '(3a)')"[named_array_list:add2lst_q] An array '",trim(element%name),"' was already registered in this list."
         call die(msg)
      endif

      if (len_trim(element%name) < 1) call die("[named_array_list:add2lst_q] empty names not allowed")

      if (.not. allocated(this%lst)) then
         allocate(na_var :: this%lst(1))
         select type (lst => this%lst)
            type is (na_var)
               lst = element
         end select
      else
         block
            type(na_var), allocatable :: tmp(:)
            allocate(tmp(ubound(this%lst(:), dim=1) + 1))
            select type (lst => this%lst)
               type is (na_var)
                  tmp(1:ubound(this%lst(:), dim=1)) = lst
                  tmp(ubound(this%lst(:), dim=1) + 1) = element
            end select
            call move_alloc(from=tmp, to=this%lst)
         end block
      endif

      if (element%name == wa_n) this%wai = ubound(this%lst(:), dim=1, kind=4)

   end subroutine add2lst_q

   subroutine add2lst_w(this, element)
      use constants,  only: fluid_n, mag_n
      use dataio_pub, only: die, msg

      implicit none

      class(na_var_list_w), intent(inout) :: this
      type(na_var_4d),      intent(in)    :: element  !< entry to be added

      if (this%exists(element%name)) then
         write(msg, '(3a)')"[named_array_list:add2lst_w] An array '",trim(element%name),"' was already registered in this list."
         call die(msg)
      endif

      if (len_trim(element%name) < 1) call die("[named_array_list:add2lst_w] empty names not allowed")

      if (.not. allocated(this%lst)) then
         allocate(na_var_4d :: this%lst(1))
         select type (lst => this%lst)
            type is (na_var_4d)
               lst = element
         end select
      else
         block
            type(na_var_4d), allocatable :: tmp(:)
            allocate(tmp(ubound(this%lst(:), dim=1) + 1))
            select type (lst => this%lst)
               type is (na_var_4d)
                  tmp(1:ubound(this%lst(:), dim=1))   = lst
                  tmp(ubound(this%lst(:), dim=1) + 1) = element
            end select
            call move_alloc(from=tmp, to=this%lst)
         end block
      endif

      if (element%name == fluid_n) this%fi = ubound(this%lst(:), dim=1, kind=4)
      if (element%name == mag_n)   this%bi = ubound(this%lst(:), dim=1, kind=4)

   end subroutine add2lst_w

!> \brief Find out which fields (cg%q and cg%w arrays) are stored in the restart file

   subroutine get_reslst(this, lst)

      use constants, only: AT_IGNORE
      use func,      only: append_int4_to_array

      implicit none

      class(na_var_list),                         intent(in)  :: this
      integer(kind=4), dimension(:), allocatable, intent(out) :: lst

      integer(kind=4) :: i

      allocate(lst(0))  ! we rely on its existence
      do i = lbound(this%lst(:), dim=1, kind=4), ubound(this%lst(:), dim=1, kind=4)
         if (this%lst(i)%restart_mode > AT_IGNORE) call append_int4_to_array(lst, i)
      enddo

   end subroutine get_reslst

!> \brief Summarize all registered fields and their properties

   subroutine print_vars(this, verbosity)

      use constants,  only: V_LOG
      use dataio_pub, only: printinfo, msg, die
      use mpisetup,   only: slave

      implicit none

      class(na_var_list),        intent(in) :: this
      integer(kind=4), optional, intent(in) :: verbosity  !< verbosity level

      integer :: i
      integer(kind=4) :: v

      v = V_LOG
      if (present(verbosity)) v = verbosity

      if (slave) return

      select type (this)
         type is (na_var_list_q)
            write(msg,'(a,i2,a)')"[named_array_list:print_vars] Found ",size(this%lst(:))," rank-3 arrays:"
         type is (na_var_list_w)
            write(msg,'(a,i2,a)')"[named_array_list:print_vars] Found ",size(this%lst(:))," rank-4 arrays:"
         class default
            call die("[named_array_list:print_vars] Unknown type of named array list")
      end select
      call printinfo(msg, v)

      do i = lbound(this%lst(:), dim=1), ubound(this%lst(:), dim=1)
         select type (lst => this%lst)
            type is (na_var)
               write(msg,'(3a,l2,a,i2,a,l2,a,i2)')"'", this%lst(i)%name, "', vital=", this%lst(i)%vital, ", restart_mode=", this%lst(i)%restart_mode, &
                    &                                ", multigrid=", this%lst(i)%multigrid, ", ord_prolong=", this%lst(i)%ord_prolong
            type is (na_var_4d)
               write(msg,'(3a,l2,a,i2,a,l2,a,i2,a,i5)')"'", this%lst(i)%name, "', vital=", this%lst(i)%vital, ", restart_mode=", this%lst(i)%restart_mode, &
                    &                                    ", multigrid=", this%lst(i)%multigrid, ", ord_prolong=", this%lst(i)%ord_prolong, &
                    &                                    ", components=", lst(i)%dim4
         end select
         call printinfo(msg, v)
      enddo

   end subroutine print_vars

   subroutine copy_var(this, other)

      implicit none

      class(na_var), intent(out) :: this   !< object to be copied to
      type(na_var),  intent(in)  :: other  !< object to be copied from

      this%name = other%name
      this%vital = other%vital
      this%restart_mode = other%restart_mode
      this%ord_prolong = other%ord_prolong
      this%multigrid = other%multigrid

   end subroutine copy_var

   subroutine copy_var_4d(this, other)

      implicit none

      class(na_var_4d), intent(out) :: this   !< object to be copied to
      type(na_var_4d),  intent(in)  :: other  !< object to be copied from

      this%name = other%name
      this%vital = other%vital
      this%restart_mode = other%restart_mode
      this%ord_prolong = other%ord_prolong
      this%dim4 = other%dim4
      this%multigrid = other%multigrid

   end subroutine copy_var_4d

   pure function get_dim4(this, iv) result(d4)

      use constants,  only: INVALID

      implicit none

      class(na_var_list_w), intent(in) :: this
      integer(kind=4),      intent(in) :: iv

      integer(kind=4) :: d4

      select type (lst => this%lst)
         type is (na_var_4d)
            d4 = lst(iv)%dim4
         class default
            d4 = INVALID
      end select

   end function get_dim4

end module named_array_list
