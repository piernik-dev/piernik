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
   use cg_list,            only: cg_list_t   ! QA_WARN intel
#endif /* __INTEL_COMPILER */
   use cg_level_connected, only: cg_level_connected_t
   use cg_list_bnd,        only: cg_list_bnd_t

   implicit none

   private
   public :: leaves, cg_leaves_t

   !>
   !! \brief Special list of grid containers that does not include fully-covered multigrid levels
   !!
   !! \todo Exclude also non-multigrid levels when fully covered
   !<

   type, extends(cg_list_bnd_t) :: cg_leaves_t ! cg_list_bnd_t is required for calling bnd_u and bnd_b
      type(cg_level_connected_t), private, pointer :: coarsest_leaves
   contains
      procedure :: update                  !< Select grids that should be included on leaves list
      procedure :: balance_and_update      !< Rebalance if required and update
      procedure :: leaf_arr3d_boundaries   !< Wrapper routine to set up all guardcells (internal, external and fine-coarse) for given rank-3 arrays
      procedure :: leaf_arr4d_boundaries   !< Wrapper routine to set up all guardcells (internal, external and fine-coarse) for given rank-4 arrays
   end type cg_leaves_t

   !>
   !! \deprecated A it is much easier to complete boundary exchanges on whole levels, the leaves list contains all grids from the base level upwards.
   !! Thus the variable name could be a bit misleading.
   !!
   !! \todo exclude base level and some higher levels if these are fully covered by finer grids (does it have side effects on visualization?)
   !! Perhaps it will be useful to keep few similar lists with slightly different inclusion criteria
   !<
   type(cg_leaves_t) :: leaves   !< grid containers not fully covered by finer grid containers

contains

!> \brief Select grids that should be included on leaves list

   subroutine update(this, str)

!      use cg_level_base,      only: base  ! can't use it because of cyclic dependency
      use cg_level_finest,    only: finest
      use cg_level_connected, only: cg_level_connected_t
      use cg_list,            only: cg_list_element
      use constants,          only: pSUM, pMAX, base_level_id, refinement_factor, base_level_id, INVALID, tmr_amr, PPP_AMR
      use dataio_pub,         only: msg, printinfo
      use domain,             only: dom
      use list_of_cg_lists,   only: all_lists
      use mpisetup,           only: master, piernik_MPI_Allreduce, nproc
      use ppp,                only: ppp_main
      use timer,              only: set_timer

      implicit none

      class(cg_leaves_t),         intent(inout) :: this  !< object invoking type-bound procedure
      character(len=*), optional, intent(in)    :: str   !< optional string identifier to show the progress of updating refinement

      type(cg_level_connected_t), pointer :: curl
      type(cg_list_element),      pointer :: cgl

      integer :: g_cnt, g_max, sum_max, ih, is, b_cnt
      integer, save :: prev_is = 0
      character(len=len(msg)), save :: prev_msg
      real :: lf
      character(len=*), parameter :: leaves_label = "leaves_update"

      call ppp_main%start(leaves_label, PPP_AMR)

      call leaves%delete
      call all_lists%register(this, "leaves")

      msg = "[cg_leaves:update] Grids on levels: "
      if (present(str)) msg(len_trim(msg)+1:) = str
      ih = len_trim(msg) + 1

      sum_max = 0
      curl => finest%level
      do while (associated(curl))
         if (curl%l%id == base_level_id) this%coarsest_leaves => curl !> \todo Find first not fully covered level
         curl => curl%coarser
      enddo
      b_cnt = INVALID
      curl => this%coarsest_leaves  ! base%level
      do while (associated(curl))
         cgl => curl%first
         do while (associated(cgl))
            call this%add(cgl%cg)
            cgl => cgl%nxt
         enddo
         g_max = curl%cnt
         call piernik_MPI_Allreduce(g_max, pMAX)
         sum_max = sum_max + g_max * nproc
         g_cnt = curl%cnt
         call piernik_MPI_Allreduce(g_cnt, pSUM)
         write(msg(len_trim(msg)+1:),'(i6)') g_cnt
         if (associated(curl, this%coarsest_leaves)) b_cnt = g_cnt
         call curl%vertical_prep
         curl => curl%finer
      enddo
      g_cnt = leaves%cnt
      call piernik_MPI_Allreduce(g_cnt, pSUM)
      write(msg(len_trim(msg)+1:), '(a,i7,a,f6.3)')", Sum: ",g_cnt, ", load balance: ",g_cnt/real(sum_max)
      is = len_trim(msg)
      write(msg(len_trim(msg)+1:), '(a,f7.3)') ", dt_wall= ", set_timer(tmr_amr)
      if (finest%level%l%id > base_level_id) then
         write(msg(len_trim(msg)+1:), '(a)')", leaves/finest: "
         lf = g_cnt/(b_cnt * real(refinement_factor)**(dom%eff_dim * finest%level%l%id))
         if (lf >= 0.0001) then
            write(msg(len_trim(msg)+1:), '(" ",f8.6)') lf
         else
            write(msg(len_trim(msg)+1:), '(" ",e9.2)') lf
         endif
      endif
      if (master .and. (msg(ih:is) /= prev_msg(ih:prev_is))) call printinfo(msg)
      prev_msg = msg
      prev_is = is

      call ppp_main%stop(leaves_label, PPP_AMR)

   end subroutine update

!> \brief Rebalance if required and update.

   subroutine balance_and_update(this, str)

      use cg_level_finest,    only: finest
      use cg_level_connected, only: cg_level_connected_t
      use constants,          only: PPP_AMR
      use ppp,                only: ppp_main

      implicit none

      class(cg_leaves_t),         intent(inout) :: this  !< object invoking type-bound procedure
      character(len=*), optional, intent(in)    :: str   !< optional string identifier to show the progress of updating refinement

      type(cg_level_connected_t), pointer :: curl
      character(len=*), parameter :: bu_label = "leaves_balance_and_update"

      call ppp_main%start(bu_label, PPP_AMR)
      curl => finest%level
      do while (associated(curl)) ! perhaps it is worth to limit to the base level
         call curl%balance_old
         call curl%sync_ru
         curl => curl%coarser
      enddo

      call this%update(str)
      call ppp_main%stop(bu_label, PPP_AMR)

   end subroutine balance_and_update

!> \brief This routine sets up all guardcells (internal, external and fine-coarse) for given rank-3 arrays

   subroutine leaf_arr3d_boundaries(this, ind, area_type, bnd_type, dir, nocorners)

      use cg_level_connected, only: cg_level_connected_t

      implicit none

      class(cg_leaves_t),        intent(in) :: this       !< the list on which to perform the boundary exchange
      integer(kind=4),           intent(in) :: ind        !< index of cg%q(:) 3d array
      integer(kind=4), optional, intent(in) :: area_type  !< defines how do we treat boundaries
      integer(kind=4), optional, intent(in) :: bnd_type   !< Override default boundary type on external boundaries (useful in multigrid solver).
                                                          !< Note that BND_PER, BND_MPI, BND_SHE and BND_COR aren't external and cannot be overridden
      integer(kind=4), optional, intent(in) :: dir        !< select only this direction
      logical,         optional, intent(in) :: nocorners  !< .when .true. then don't care about proper edge and corner update

      type(cg_level_connected_t), pointer   :: curl

      curl => this%coarsest_leaves
      do while (associated(curl))
         ! OPT this results in duplicated calls to level_3d_boundaries for levels from this%coarsest_leaves to finest%level%coarser
         !> \todo implement it with lower level routines to remove this duplication
         call curl%arr3d_boundaries(ind, area_type=area_type, bnd_type=bnd_type, dir=dir, nocorners=nocorners)
         curl => curl%finer
      enddo

   end subroutine leaf_arr3d_boundaries

!> \brief This routine sets up all guardcells (internal, external and fine-coarse) for given rank-4 arrays

   subroutine leaf_arr4d_boundaries(this, ind, area_type, dir, nocorners)

      use cg_level_connected, only: cg_level_connected_t

      implicit none

      class(cg_leaves_t),        intent(in) :: this       !< the list on which to perform the boundary exchange
      integer(kind=4),           intent(in) :: ind        !< index of cg%w(:) 4d array
      integer(kind=4), optional, intent(in) :: area_type  !< defines how do we treat boundaries
      integer(kind=4), optional, intent(in) :: dir        !< select only this direction
      logical,         optional, intent(in) :: nocorners  !< .when .true. then don't care about proper edge and corner update

      type(cg_level_connected_t), pointer   :: curl

      curl => this%coarsest_leaves
      do while (associated(curl))
         ! OPT this results in duplicated calls to level_4d_boundaries for levels from this%coarsest_leaves to finest%level%coarser
         !> \todo implement it with lower level routines to remove this duplication
         call curl%arr4d_boundaries(ind, area_type=area_type, dir=dir, nocorners=nocorners)
         curl => curl%finer
      enddo

   end subroutine leaf_arr4d_boundaries

end module cg_leaves
