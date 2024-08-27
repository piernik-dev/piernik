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

!> \brief Unified refinement criteria for user-provided routine

module unified_ref_crit_user

   use unified_ref_crit, only: urc

   implicit none

   private
   public :: urc_user, mark_urc_user

!>
!! \brief Things that should be common for all user criteria.
!!
!! All necessary parameters are supposed to be visible in the initproblem,
!! where the user routine is supposed to be defined.
!!
!! \warning (spaghetti) urc_filter also defines plotfield in the same role
!<

   type, extends(urc) :: urc_user
      logical :: plotfield  !< create a 3D array to keep the value of refinement criterion
      procedure(mark_urc_user), pointer :: mark_u
   contains
      procedure :: mark => mark_user  !< a routine that takes a cg and leaves suggestions on refining
   end type urc_user

   interface

!> \brief Mark refinements on given grid container

      subroutine mark_urc_user(this, cg)

         use grid_cont, only: grid_container
         import urc_user

         implicit none

         class(urc_user),               intent(inout) :: this  !< an object invoking the type-bound procedure
         type(grid_container), pointer, intent(inout) :: cg    !< current grid piece

      end subroutine mark_urc_user

   end interface

   interface urc_user
      procedure :: init
   end interface urc_user

contains

!> \brief Mark refinements on given grid container

   subroutine mark_user(this, cg)

      use dataio_pub, only: die
      use grid_cont,  only: grid_container

      implicit none

      class(urc_user),               intent(inout) :: this  !< an object invoking the type-bound procedure
      type(grid_container), pointer, intent(inout) :: cg    !< current grid piece

      if (associated(this%mark_u)) then
         call this%mark_u(cg)
      else
         call die("[URC_user:mark_user] .not. associated(this%mark_u)")
      endif

   end subroutine mark_user

!> \brief A simple constructor for user-provided refinement criteria

   function init(user_mark, plotfield) result(this)

      use constants,  only: V_INFO
      use dataio_pub, only: printinfo, msg
      use mpisetup,   only: master

      implicit none

      logical, intent(in)      :: plotfield  !< create an array to keep the value of refinement criterion
      procedure(mark_urc_user) :: user_mark  !< user-provided routine

      type(urc_user) :: this  !< an object to be constructed

      this%plotfield = plotfield
      this%mark_u => user_mark

      if (master) then
         write(msg, '(5a,2g13.5,a)')"[URC user]  Initializing user-provided criteria"
         if (plotfield)  write(msg(len_trim(msg)+1:), '(a)') ", with plotfield"
         call printinfo(msg, V_INFO)
      endif

   end function init

end module unified_ref_crit_user
