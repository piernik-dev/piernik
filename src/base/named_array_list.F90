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
!! \brief Lists of 3D and 4D named arrays and associated routines
!!
!! \details The qna(:) and wna(:) lists contain only descriptions of registeres named arrays.
!! Named arrays themselves are defined in "named_array" module.
!<

module named_array_list

   use constants, only: dsetnamelen

   implicit none

   private
   public :: na_var, qna, wna

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

   !> \brief the generic list of named arrays with supporting routines
   type :: na_var_list
      type(na_var), dimension(:), allocatable :: lst
    contains
      procedure :: ind                                           !< Get the index of a named array of given name.
      procedure :: exists                                        !< Check if a named array of given name is already registered
      procedure :: print_vars                                    !< Write a summary on registered fields. Can be useful for debugging
      procedure :: add2lst                                       !< Add a named array properties to the list
   end type na_var_list

   ! types with indices of the most commonly used arrays stored in cg%w and cg%q

   !> \brief the most commonly used 3D named array is wa, thus we add a shortcut here
   type, extends(na_var_list) :: na_var_list_q
      integer :: wai                                   !< auxiliary array : cg%q(qna%wai)
   end type na_var_list_q

   !> \brief the most commonly used 4D named arraya are u and b, thus we add shortcuts here
   type, extends(na_var_list) :: na_var_list_w
      integer :: fi                                    !< fluid           : cg%w(wna%fi)
      integer :: bi                                    !< magnetic field  : cg%w(wna%bi)
   end type na_var_list_w

   type(na_var_list_q) :: qna !< list of registered 3D named arrays
   type(na_var_list_w) :: wna !< list of registered 4D named arrays

contains

!>
!! \brief Get the index of a named array of given name.
!!
!! \warning OPT The indices aren't updated so cache them, whenever possible
!<
   function ind(this, name) result(rind)

      use dataio_pub,  only: die, msg, warn

      implicit none

      class(na_var_list), intent(inout) :: this
      character(len=*),   intent(in)    :: name

      integer :: rind, i

      rind = 0

      if (allocated(this%lst)) then
         do i = lbound(this%lst, dim=1), ubound(this%lst, dim=1)
            if (trim(name) == this%lst(i)%name) then
               if (rind /= 0) then
                  write(msg, '(2a)') "[named_array_list:ind] multiple entries with the same name: ", trim(name)
                  call die(msg)
               endif
               rind = i
            endif
         enddo
      endif

      if (rind == 0) then
         write(msg, '(2a)') "[named_array_list:ind] requested entry not found: ", trim(name)
         call warn(msg)
      endif

   end function ind

!> \brief Check if a named array of given name is already registered

   function exists(this, name)

      use dataio_pub,  only: die, msg

      implicit none

      class(na_var_list), intent(inout) :: this
      character(len=*),   intent(in)    :: name

      logical :: exists
      integer :: i

      exists = .false.

      if (allocated(this%lst)) then
         do i = lbound(this%lst(:), dim=1), ubound(this%lst(:), dim=1)
            if (trim(name) ==  this%lst(i)%name) then
               if (exists) then
                  write(msg, '(2a)') "[named_array_list:exists] multiple entries with the same name: ", trim(name)
                  call die(msg)
               endif
               exists = .true.
            endif
         enddo
      endif

   end function exists

!> \brief Add a named array properties to the list

   subroutine add2lst(this, element)

      use constants,  only: fluid_n, mag_n, wa_n
      use dataio_pub, only: die, msg

      implicit none

      class(na_var_list), intent(inout) :: this
      type(na_var),       intent(in)    :: element

      type(na_var), dimension(:), allocatable :: tmp

      if (this%exists(element%name)) then
         write(msg, '(3a)')"[named_array_list:add2lst] An array '",trim(element%name),"' was already registered in this list."
         call die(msg)
      endif

      if (.not. allocated(this%lst)) then
         allocate(this%lst(1))
      else
         allocate(tmp(lbound(this%lst(:),dim=1):ubound(this%lst(:), dim=1) + 1))
         tmp(:ubound(this%lst(:), dim=1)) = this%lst(:)
         call move_alloc(from=tmp, to=this%lst)
         !! \warning slight memory leak here, e.g. in use at exit: 572 bytes in 115 blocks, perhaps associated with na_var%position
      endif
      this%lst(ubound(this%lst(:), dim=1)) = element

      select type(this)
         type is (na_var_list_w)
            if (element%name == fluid_n) this%fi  = ubound(this%lst(:), dim=1)
            if (element%name == mag_n)   this%bi  = ubound(this%lst(:), dim=1)
         type is (na_var_list_q)
            if (element%name == wa_n)    this%wai = ubound(this%lst(:), dim=1)
      end select

   end subroutine add2lst

!> \brief Summarize all registered fields and their properties

   subroutine print_vars(this)

      use constants,  only: INVALID
      use dataio_pub, only: printinfo, warn, msg
      use mpisetup,   only: slave

      implicit none

      class(na_var_list), intent(inout) :: this

      integer :: i, d3

      if (slave) return

      d3 = count(this%lst(:)%dim4 == INVALID)

      if (d3 /= 0) then
         write(msg,'(a,i2,a)')"[named_array_list:print_vars] Found ",size(this%lst(:))," rank-3 arrays:"
         call printinfo(msg)
      endif
      if (count(this%lst(:)%dim4 /= INVALID) /= 0) then
         write(msg,'(a,i2,a)')"[named_array_list:print_vars] Found ",size(this%lst(:))," rank-4 arrays:"
         call printinfo(msg)
         if (d3 /=0) call warn("[named_array_list:print_vars] Both rank-3 and rank-4 named arrays are present in the same list!")
      endif

      do i = lbound(this%lst(:), dim=1), ubound(this%lst(:), dim=1)
         if (this%lst(i)%dim4 == INVALID) then
            write(msg,'(3a,l2,a,i2,a,l2,2(a,i2))')"'", this%lst(i)%name, "', vital=", this%lst(i)%vital, ", restart_mode=", this%lst(i)%restart_mode, &
                 &                                ", multigrid=", this%lst(i)%multigrid, ", ord_prolong=", this%lst(i)%ord_prolong, ", position=", this%lst(i)%position(:)
         else
            write(msg,'(3a,l2,a,i2,a,l2,2(a,i2),a,100i2)')"'", this%lst(i)%name, "', vital=", this%lst(i)%vital, ", restart_mode=", this%lst(i)%restart_mode, &
                 &                                        ", multigrid=", this%lst(i)%multigrid, ", ord_prolong=", this%lst(i)%ord_prolong, &
                 &                                        ", components=", this%lst(i)%dim4, ", position=", this%lst(i)%position(:)
         endif
         call printinfo(msg)
      enddo

   end subroutine print_vars

end module named_array_list
