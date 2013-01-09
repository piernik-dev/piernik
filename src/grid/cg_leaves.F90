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

#if defined(__INTEL_COMPILER)
   !! \deprecated remove this clause as soon as Intel Compiler gets required
   !! features and/or bug fixes, it's needed for 12.1, fixed in 13.0 but the
   !! latter is broken and we cannot use it yet
   use cg_list,            only: cg_list_T   ! QA_WARN intel
#endif /* __INTEL_COMPILER */
   use cg_level_connected, only: cg_level_connected_T
   use cg_list_bnd,        only: cg_list_bnd_T

   implicit none

   private
   public :: leaves, cg_leaves_T

   !>
   !! \brief Special list of grid containers that does not include fully-covered multigrid levels
   !!
   !! \todo Exclude also non-multigrid levels when fully covered
   !<

   type, extends(cg_list_bnd_T) :: cg_leaves_T
      type(cg_level_connected_T), pointer :: coarsest_leaves
   contains
      procedure :: update           !< Select grids that should be included on leaves list
      procedure :: arr3d_boundaries !< Wrapper routine to set up all guardcells (internal, external and fine-coarse) for given rank-3 arrays
                                    !< \warning the name 'arr3d_boundaries' intentionally collides with cg_level_connested::arr3d_boundaries
                                    !< \todo unite this routine with cg_level_connected::arr3d_boundaries somehow
      !< \todo fix bnd_routines contents for this type as soon as possible
      procedure :: arr4d_boundaries !< Wrapper routine to set up all guardcells (internal, external and fine-coarse) for given rank-4 arrays
      procedure :: internal_bnd_3d  !< Wrapper routine to set up internal boundaries for for given rank-3 arrays
      procedure :: internal_bnd_4d  !< Wrapper routine to set up internal boundaries for for given rank-4 arrays
      procedure :: external_bnd_3d  !< Wrapper routine to set up external boundaries for for given rank-3 arrays
      procedure :: external_bnd_4d  !< Wrapper routine to set up external boundaries for for given rank-4 arrays
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

      use cg_level_connected, only: base_lev, cg_level_connected_T
      use cg_list,            only: cg_list_element
      use constants,          only: pSUM
      use dataio_pub,         only: msg, printinfo
      use list_of_cg_lists,   only: all_lists
      use mpisetup,           only: master, piernik_MPI_Allreduce

      implicit none

      class(cg_leaves_T), intent(inout)   :: this          !< object invoking type-bound procedure

      type(cg_level_connected_T), pointer :: curl
      type(cg_list_element),      pointer :: cgl

      integer :: g_cnt

      call leaves%delete

      call all_lists%register(this, "leaves")

      msg = "[cg_leaves:update] Leaves on levels: "
      curl => base_lev
      this%coarsest_leaves => curl !> \todo Start from first not fully covered level
      do while (associated(curl))
         cgl => curl%first
         do while (associated(cgl))
            call this%add(cgl%cg)
            cgl => cgl%nxt
         enddo
         g_cnt = curl%cnt
         call piernik_MPI_Allreduce(g_cnt, pSUM)
         write(msg(len_trim(msg)+1:),'(i6)') g_cnt
         curl => curl%finer
      enddo
      g_cnt = leaves%cnt
      call piernik_MPI_Allreduce(g_cnt, pSUM)
      write(msg(len_trim(msg)+1:), '(a,i7,a)')",      Total: ",g_cnt, " leaves"
      if (master) call printinfo(msg)

   end subroutine update

!> \brief This routine sets up all guardcells (internal, external and fine-coarse) for given rank-3 arrays

   subroutine arr3d_boundaries(this, ind, area_type, bnd_type)

      use cg_level_connected, only: cg_level_connected_T

      implicit none

      class(cg_leaves_T),        intent(in) :: this       !< the list on which to perform the boundary exchange
      integer,                   intent(in) :: ind        !< Negative value: index of cg%q(:) 3d array
      integer(kind=4), optional, intent(in) :: area_type  !< defines how do we treat boundaries
      integer(kind=4), optional, intent(in) :: bnd_type   !< Override default boundary type on external boundaries (useful in multigrid solver).
                                                          !< Note that BND_PER, BND_MPI, BND_SHE and BND_COR aren't external and cannot be overridden

      type(cg_level_connected_T), pointer   :: curl

      curl => this%coarsest_leaves
      do while (associated(curl))
         call curl%arr3d_boundaries(ind, area_type, bnd_type)
         curl => curl%finer
      enddo

   end subroutine arr3d_boundaries

!>
!! \brief This routine sets up all guardcells (internal, external and fine-coarse) for given rank-4 arrays
!! \todo if there is any routine for base_lev then place it here as a wrapper
!<
   subroutine arr4d_boundaries(this) !, ind, nb, area_type, bnd_type, corners)

      use cg_level_connected, only: cg_level_connected_T
      use dataio_pub,         only: die

      implicit none

      class(cg_leaves_T),        intent(in) :: this       !< the list on which to perform the boundary exchange
!      integer,                   intent(in) :: ind        !< Negative value: index of cg%q(:) 3d array
!      integer(kind=4), optional, intent(in) :: nb         !< number of grid cells to exchange
!      integer(kind=4), optional, intent(in) :: area_type  !< defines how do we treat boundaries
!      integer(kind=4), optional, intent(in) :: bnd_type   !< Override default boundary type on external boundaries (useful in multigrid solver).
!                                                          !< Note that BND_PER, BND_MPI, BND_SHE and BND_COR aren't external and cannot be overridden
!      logical,         optional, intent(in) :: corners    !< When present and .true. then call internal_boundaries_3d for each direction separately

      type(cg_level_connected_T), pointer   :: curl

      call die("[cg_leaves::arr4d_boundaries] This routine has not been implemented yet.")

      curl => this%coarsest_leaves
!      do while (associated(curl))
!         call curl%arr4d_boundaries(ind, nb, area_type, bnd_type, corners)
!         curl => curl%finer
!      enddo

   end subroutine arr4d_boundaries

!>
!! \brief Wrapper routine to set up internal boundaries for for given rank-3 arrays
!! \todo make it completed
!<
   subroutine internal_bnd_3d(this, ind, dim)

      use cg_level_connected, only: cg_level_connected_T
      use dataio_pub,         only: die

      implicit none

      class(cg_leaves_T),        intent(in) :: this       !< the list on which to perform the boundary exchange
      integer,                   intent(in) :: ind        !< Negative value: index of cg%q(:) 3d array
      integer(kind=4), optional, intent(in) :: dim        !< do the internal boundaries only in the specified dimension

      type(cg_level_connected_T), pointer   :: curl

      curl => this%coarsest_leaves
      do while (associated(curl))
         call curl%internal_boundaries_3d(ind, dim)
         curl => curl%finer
         if (associated(curl)) call die("[cg_leaves::internal_bnd_3d] This routine does not work with finer levels yet")
      enddo

   end subroutine internal_bnd_3d

!>
!! \brief Wrapper routine to set up internal boundaries for for given rank-4 arrays
!! \todo make it completed
!<
   subroutine internal_bnd_4d(this, ind, dim)

      use cg_level_connected, only: cg_level_connected_T
      use dataio_pub,         only: die

      implicit none

      class(cg_leaves_T),        intent(in) :: this       !< the list on which to perform the boundary exchange
      integer,                   intent(in) :: ind        !< Negative value: index of cg%q(:) 3d array
      integer(kind=4), optional, intent(in) :: dim        !< do the internal boundaries only in the specified dimension

      type(cg_level_connected_T), pointer   :: curl

      curl => this%coarsest_leaves
      do while (associated(curl))
         call curl%internal_boundaries_4d(ind, dim)
         curl => curl%finer
         if (associated(curl)) call die("[cg_leaves::internal_bnd_4d] This routine does not work with finer levels yet")
      enddo

   end subroutine internal_bnd_4d

!>
!! \brief Wrapper routine to set up external boundaries for for given rank-3 arrays
!! \todo make it completed
!<
   subroutine external_bnd_3d(this, ind, area_type, bnd_type)

      use cg_level_connected, only: cg_level_connected_T
      use dataio_pub,         only: die

      implicit none

      class(cg_leaves_T),        intent(in) :: this       !< the list on which to perform the boundary exchange
      integer,                   intent(in) :: ind        !< Negative value: index of cg%q(:) 3d array
      integer(kind=4), optional, intent(in) :: area_type  !< defines how do we treat boundaries
      integer(kind=4), optional, intent(in) :: bnd_type   !< Override default boundary type on external boundaries (useful in multigrid solver).
                                                          !< Note that BND_PER, BND_MPI, BND_SHE and BND_COR aren't external and cannot be overridden
      type(cg_level_connected_T), pointer   :: curl

      curl => this%coarsest_leaves
      do while (associated(curl))
         call curl%external_boundaries(ind, area_type, bnd_type)
         curl => curl%finer
         if (associated(curl)) call die("[cg_leaves::external_bnd_3d] This routine does not work with finer levels yet")
      enddo

   end subroutine external_bnd_3d

!>
!! \brief Wrapper routine to set up external boundaries for for given rank-4 arrays
!! \todo make it completed
!<
   subroutine external_bnd_4d(this) !, ind, area_type, bnd_type)

      use cg_level_connected, only: cg_level_connected_T
      use dataio_pub,         only: die

      implicit none

      class(cg_leaves_T),        intent(in) :: this       !< the list on which to perform the boundary exchange
!      integer,                   intent(in) :: ind        !< Negative value: index of cg%q(:) 3d array
!      integer(kind=4), optional, intent(in) :: area_type  !< defines how do we treat boundaries
!      integer(kind=4), optional, intent(in) :: bnd_type   !< Override default boundary type on external boundaries (useful in multigrid solver).
                                                          !< Note that BND_PER, BND_MPI, BND_SHE and BND_COR aren't external and cannot be overridden
      type(cg_level_connected_T), pointer   :: curl

      call die("[cg_leaves::external_bnd_4d] This routine has not been implemented yet.")
      curl => this%coarsest_leaves
!      do while (associated(curl))
!         call curl%external_boundaries(ind, area_type, bnd_type)
!         curl => curl%finer
!      enddo

   end subroutine external_bnd_4d

end module cg_leaves
