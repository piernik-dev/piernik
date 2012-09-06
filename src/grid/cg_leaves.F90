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
!! \brief This module contains list of "leaves"
!!
!! \details EXPLAIN MORE
!<

module cg_leaves

   use cg_list_bnd, only: cg_list_bnd_T

   implicit none

   private
   public :: leaves

   type, extends(cg_list_bnd_T) :: cg_leaves_T
    contains
      procedure :: update !> Select grids that should be included on leaves list
   end type cg_leaves_T

   !>
   !! \deprecated A it is much easied to complete boundary exchanges on whole levels, the leaves list contains all grids from the base level upwards.
   !! Thus the variable name could be a bit misleading.
   !!
   !! \todo exclude base level and some higher levels if these are fully covered by finer grids (does it have side effects on visualization?)
   !! Perhaps it will be useful to keep few similar lists with slightly different inclusion criteria
   !<
   type(cg_leaves_T) :: leaves   !< grid containers not fully covered by finer grid containers

contains

!> \brief Select grids that should be included on leaves list

   subroutine update(this)

      use cg_list,          only: cg_list_element
      use cg_list_level,    only: base_lev, cg_list_level_T
      use constants,        only: I_ONE
      use dataio_pub,       only: msg, printinfo
      use list_of_cg_lists, only: all_lists
      use mpi,              only: MPI_INTEGER, MPI_SUM
      use mpisetup,         only: master, comm, mpi_err

      implicit none

      class(cg_leaves_T), intent(inout) :: this          !< object invoking type-bound procedure

      type(cg_list_level_T), pointer :: curl
      type(cg_list_element), pointer :: cgl

      integer :: g_cnt

      call leaves%delete

      call all_lists%register(this, "leaves")

      msg = "[cg_leaves:update] Leaves on levels: "
      curl => base_lev
      do while (associated(curl))
         cgl => curl%first
         do while (associated(cgl))
            call this%add(cgl%cg)
            cgl => cgl%nxt
         enddo
         call MPI_Allreduce(curl%cnt, g_cnt, I_ONE, MPI_INTEGER, MPI_SUM, comm, mpi_err)
         write(msg(len_trim(msg)+1:),'(i6)') g_cnt
         curl => curl%finer
      enddo
      call MPI_Allreduce(leaves%cnt, g_cnt, I_ONE, MPI_INTEGER, MPI_SUM, comm, mpi_err)
      write(msg(len_trim(msg)+1:), '(a,i7,a)')",      Total: ",g_cnt, " leaves"
      if (master) call printinfo(msg)

   end subroutine update

end module cg_leaves
