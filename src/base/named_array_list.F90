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
!!       Name and position of the 3D named arrays should be scalar, for 4D arrays should be 1D vectors,
!!       dim4 should occur only in 4D type
!<

module named_array_list

   use constants, only: dsetnamelen, INVALID

   implicit none

   private
   public :: na_var, na_var_list, na_var_list_q, na_var_list_w, qna, wna

   !> \brief Common properties of 3D and 4D named arrays
   type :: na_var
      character(len=dsetnamelen)                 :: name          !< a user-provided id for the array
      logical                                    :: vital         !< fields that are subject of automatic prolongation and restriction (e.g. state variables)
      integer(kind=4)                            :: restart_mode  !< AT_BACKUP, AT_IGNORE: do not write to restart
                                                                  !< AT_NO_B write without ext. boundaries
                                                                  !< AT_OUT_B write with ext. boundaries
                                                                  !< \todo position /= VAR_CENTER should automatically force AT_OUT_B if AT_IGNORE was not chosen
      integer(kind=4)                            :: ord_prolong   !< Prolongation order for the variable
      integer(kind=4), allocatable, dimension(:) :: position      !< VAR_CENTER by default, also possible VAR_CORNER and VAR_[XYZ]FACE
      integer(kind=4)                            :: dim4          !< <=0 for 3D arrays, >0 for 4D arrays
      logical                                    :: multigrid     !< .true. for variables that may exist below base level (e.g. work fields for multigrid solver)
   end type na_var

   !> \brief the generic list of named arrays with supporting routines
   type, abstract :: na_var_list
      type(na_var), dimension(:), allocatable :: lst
   contains
      procedure, private ::  find_ind                            !< Get the index of a named array of given name. Don't die when can't find requested field.
      procedure :: ind                                           !< Get the index of a named array of given name.
      procedure :: exists                                        !< Check if a named array of given name is already registered
      procedure :: print_vars                                    !< Write a summary on registered fields. Can be useful for debugging
      procedure :: add2lst                                       !< Add a named array properties to the list
      procedure :: get_reslst                                    !< List of fields that are stored in the restart file
   end type na_var_list

   ! types with indices of the most commonly used arrays stored in cg%w and cg%q

   !> \brief the most commonly used 3D named array is wa, thus we add a shortcut here
   type, extends(na_var_list) :: na_var_list_q
      integer(kind=4) :: wai = INVALID                           !< auxiliary array : cg%q(qna%wai)
   end type na_var_list_q

   !> \brief the most commonly used 4D named arrays are u and b, thus we add shortcuts here
   type, extends(na_var_list) :: na_var_list_w
      integer(kind=4) :: fi = INVALID                            !< fluid           : cg%w(wna%fi)
      integer(kind=4) :: bi = INVALID                            !< magnetic field  : cg%w(wna%bi)
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

   subroutine add2lst(this, element)

      use constants,  only: fluid_n, mag_n, wa_n
      use dataio_pub, only: die, msg

      implicit none

      class(na_var_list), intent(inout) :: this
      type(na_var),       intent(in)    :: element

      type(na_var), dimension(:), allocatable :: tmp
      integer :: i

      if (this%exists(element%name)) then
         write(msg, '(3a)')"[named_array_list:add2lst] An array '",trim(element%name),"' was already registered in this list."
         call die(msg)
      endif

      if (len_trim(element%name) < 1) call die("[named_array_list:add2lst] empty names not allowed")

      if (.not. allocated(this%lst)) then
         allocate(this%lst(1))
      else
         allocate(tmp(lbound(this%lst(:),dim=1):ubound(this%lst(:), dim=1) + 1))
         tmp(:ubound(this%lst(:), dim=1)) = this%lst(:)
         ! manually deallocate arrays inside user-types, as it seems that move_alloc is unable to do that
         do i = lbound(this%lst(:), dim=1), ubound(this%lst(:), dim=1)
            if (allocated(this%lst(i)%position)) deallocate(this%lst(i)%position)
         enddo
         call move_alloc(from=tmp, to=this%lst)
      endif
      this%lst(ubound(this%lst(:), dim=1)) = element

      select type(this)
         type is (na_var_list_w)
            if (element%name == fluid_n) this%fi  = ubound(this%lst(:), dim=1, kind=4)
            if (element%name == mag_n)   this%bi  = ubound(this%lst(:), dim=1, kind=4)
         type is (na_var_list_q)
            if (element%name == wa_n)    this%wai = ubound(this%lst(:), dim=1, kind=4)
      end select

   end subroutine add2lst

!> \brief Find out which fields (cg%q and cg%w arrays) are stored in the restart file

   subroutine get_reslst(this, lst)

      use constants, only: AT_IGNORE
      use func,      only: append_int_to_array

      implicit none

      class(na_var_list),                 intent(in)  :: this
      integer, dimension(:), allocatable, intent(out) :: lst

      integer :: i

      allocate(lst(0))  ! we rely on its existence
      do i = lbound(this%lst(:), dim=1), ubound(this%lst(:), dim=1)
         if (this%lst(i)%restart_mode > AT_IGNORE) call append_int_to_array(lst, i)
      enddo

   end subroutine get_reslst

!> \brief Summarize all registered fields and their properties

   subroutine print_vars(this, verbosity)

      use constants,  only: INVALID, V_LOG
      use dataio_pub, only: printinfo, warn, msg
      use mpisetup,   only: slave

      implicit none

      class(na_var_list),        intent(in) :: this
      integer(kind=4), optional, intent(in) :: verbosity  !< verbosity level

      integer :: i, d3
      integer(kind=4) :: v

      v = V_LOG
      if (present(verbosity)) v = verbosity

      if (slave) return

      d3 = count(this%lst(:)%dim4 == INVALID)

      if (d3 /= 0) then
         write(msg,'(a,i2,a)')"[named_array_list:print_vars] Found ",size(this%lst(:))," rank-3 arrays:"
         call printinfo(msg, v)
      endif
      if (count(this%lst(:)%dim4 /= INVALID) /= 0) then
         write(msg,'(a,i2,a)')"[named_array_list:print_vars] Found ",size(this%lst(:))," rank-4 arrays:"
         call printinfo(msg, v)
         if (d3 /=0) call warn("[named_array_list:print_vars] Both rank-3 and rank-4 named arrays are present in the same list!")
      endif

      do i = lbound(this%lst(:), dim=1), ubound(this%lst(:), dim=1)
         if (this%lst(i)%dim4 == INVALID) then
            write(msg,'(3a,l2,a,i2,a,l2,2(a,i2))')"'", this%lst(i)%name, "', vital=", this%lst(i)%vital, ", restart_mode=", this%lst(i)%restart_mode, &
                 &                                ", multigrid=", this%lst(i)%multigrid, ", ord_prolong=", this%lst(i)%ord_prolong, ", position=", this%lst(i)%position(:)
         else
            write(msg,'(3a,l2,a,i2,a,l2,a,i2,a,i3,a)')"'", this%lst(i)%name, "', vital=", this%lst(i)%vital, ", restart_mode=", this%lst(i)%restart_mode, &
                 &                                    ", multigrid=", this%lst(i)%multigrid, ", ord_prolong=", this%lst(i)%ord_prolong, &
                 &                                    ", components=", this%lst(i)%dim4, ", position="
            if (any(this%lst(i)%position /= 0)) then
               ! theoretically we can print even about (len(msg) - len_trim(msg))/2 ~= 452 entries for position
               if (size(this%lst(i)%position) <= 400) then  ! hardcoded integer here and in the formats below
                  write(msg(len_trim(msg)+1:), '(400i2)') this%lst(i)%position(:)
               else
                  write(msg(len_trim(msg)+1:), '(400i2,a)') this%lst(i)%position(:400), " ... ??? ..."
               endif
            else
               write(msg(len_trim(msg)+1:), '(a,i4,a)') "[ ", size(this%lst(i)%position), " * 0 ]"
            endif
         endif
         call printinfo(msg, v)
      enddo

   end subroutine print_vars

end module named_array_list
