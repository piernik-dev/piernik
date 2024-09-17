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
   use cg_list,            only: cg_list_element_ptr
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
      type(cg_level_connected_t), private, pointer         :: coarsest_leaves  !< For optimization
      type(cg_list_element_ptr), dimension(:), allocatable :: up_to_level      !< An array of pointers to get cg leaves no finer than a given level
   contains
      procedure :: update                  !< Select grids that should be included on leaves list
      procedure :: balance_and_update      !< Rebalance if required and update
      procedure :: leaf_arr3d_boundaries   !< Wrapper routine to set up all guardcells (internal, external and fine-coarse) for given rank-3 arrays on leaves
      procedure :: leaf_arr4d_boundaries   !< Wrapper routine to set up all guardcells (internal, external and fine-coarse) for given rank-4 arrays on leaves
      procedure :: prioritized_cg          !< Return a leaves list with different ordering of cg to optimize fine->coarse flux transfer
      procedure :: leaf_only_cg            !< Return a leaves list without fully covered cg
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
      use constants,          only: pSUM, pMAX, base_level_id, refinement_factor, INVALID, tmr_amr, PPP_AMR, V_VERBOSE
      use dataio_pub,         only: msg, printinfo
      use domain,             only: dom
      use list_of_cg_lists,   only: all_lists
      use mpisetup,           only: master, piernik_MPI_Allreduce, nproc
      use ppp,                only: ppp_main
      use timer,              only: set_timer
#ifdef NBODY
      use  func,              only: operator(.notequals.)
      use particle_utils,     only: global_balance_particles
#endif /* NBODY */

      implicit none

      class(cg_leaves_t), intent(inout) :: this  !< object invoking type-bound procedure
      character(len=*),   intent(in)    :: str   !< optional string identifier to show the progress of updating refinement

      type(cg_level_connected_t), pointer :: curl
      type(cg_list_element),      pointer :: cgl

      integer :: g_cnt, g_max, sum_max, ih, is, b_cnt, i
      integer, save :: prev_is = 0
      character(len=len(msg)), save :: prev_msg
      real :: lf
      character(len=*), parameter :: leaves_label = "leaves_update"
#ifdef NBODY
      real :: lb_part
      real, save :: prev_lb_part = -1.
#endif /* NBODY */

      call ppp_main%start(leaves_label, PPP_AMR)

      call leaves%delete
      call all_lists%register(this, "leaves")

      msg = "[cg_leaves:update] Grids on levels:" // str
      ih = len_trim(msg) + 1

      curl => finest%level
      do while (associated(curl))
         if (curl%l%id == base_level_id) this%coarsest_leaves => curl
         !> \todo Find first not fully covered level, but current IO implementation depends on leaves
         !!       as a complete set of cg from base level to finest level, so be careful.
         curl => curl%coarser
      enddo

      ! Create leaves sorted from finest to coarser levels to easily obtain pointers to
      ! cg subsets that doesn't contain cgs finer than given level
      if (allocated(this%up_to_level)) deallocate(this%up_to_level)
      allocate(this%up_to_level(base_level_id:finest%level%l%id))
      curl => finest%level
      do while (associated(curl))
         nullify(this%up_to_level(curl%l%id)%p)
         cgl => curl%first
         do while (associated(cgl))
            call this%add(cgl%cg)
            if (associated(curl%first%cg, this%last%cg)) this%up_to_level(curl%l%id)%p => this%last
            ! On top levels it will be unassociated if proc doesn't have that fine cg, fixed below
            cgl => cgl%nxt
         enddo
         curl => curl%coarser
         if (associated(curl)) then
            if (curl%l%id < this%coarsest_leaves%l%id) nullify(curl)
         endif
      enddo
      do i = lbound(this%up_to_level, dim=1) + 1, ubound(this%up_to_level, dim=1)
         ! Set up pointers for top levels in case ihese weren't associated due to lack of that fine cgs
         if (.not. associated(this%up_to_level(i)%p)) this%up_to_level(i)%p => this%up_to_level(i-1)%p
      enddo

      sum_max = 0
      b_cnt = INVALID
      curl => this%coarsest_leaves  ! base%level
      do while (associated(curl))
         g_max = curl%cnt
         call piernik_MPI_Allreduce(g_max, pMAX)
         sum_max = sum_max + g_max * nproc
         g_cnt = curl%cnt
         call piernik_MPI_Allreduce(g_cnt, pSUM)
         write(msg(len_trim(msg)+1:),'(i6)') g_cnt
         if (associated(curl, this%coarsest_leaves)) b_cnt = g_cnt
         curl => curl%finer
      enddo

      g_cnt = leaves%cnt
      call piernik_MPI_Allreduce(g_cnt, pSUM)
      write(msg(len_trim(msg)+1:), '(a,i7,a,f6.3)')", Sum: ",g_cnt, ", cg load balance: ",g_cnt/real(sum_max)
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
      ! Consider a higher verbosity when str == ' ( derefine )'
      if (master .and. (msg(ih:is) /= prev_msg(ih:prev_is))) call printinfo(msg, V_VERBOSE)
      prev_msg = msg
      prev_is = is

#ifdef NBODY
      lb_part = global_balance_particles()
      if (master .and. (lb_part .notequals. prev_lb_part)) then
         write(msg, '(a,f7.4)')"[cg_leaves:update] particles load balance: ", lb_part
         call printinfo(msg, V_VERBOSE)
      endif
      prev_lb_part = lb_part
#endif /* NBODY */

      call ppp_main%stop(leaves_label, PPP_AMR)

   end subroutine update

!> \brief Rebalance if required and update grid structure.

   subroutine balance_and_update(this, str)

      use cg_level_finest,    only: finest
      use cg_level_connected, only: cg_level_connected_t
      use rebalance,          only: rebalance_all
      use constants,          only: PPP_AMR
      use ppp,                only: ppp_main

      implicit none

      class(cg_leaves_t), intent(inout) :: this  !< object invoking type-bound procedure
      character(len=*),   intent(in)    :: str   !< optional string identifier to show the progress of updating refinement

      type(cg_level_connected_t), pointer :: curl
      character(len=*), parameter :: bu_label = "leaves_balance_and_update"

      call ppp_main%start(bu_label, PPP_AMR)

      call rebalance_all

      curl => finest%level
      do while (associated(curl))
         call curl%check_update_all  ! here the buffers for horizontal communication are updated
         call curl%sync_ru  ! no need to update this%recently_changed here (was done in check_update_all)
         curl => curl%coarser
      enddo

      call this%update(str)
      call this%coarsest_leaves%update_verticals  ! update prolongation/restriction buffers

      call ppp_main%stop(bu_label, PPP_AMR)

   end subroutine balance_and_update

!> \brief This routine sets up all guardcells (internal, external and fine-coarse) for given rank-3 arrays on leaves.

   subroutine leaf_arr3d_boundaries(this, ind, area_type, bnd_type, dir, nocorners)

      use cg_level_connected, only: cg_level_connected_t
      use constants,          only: PPP_AMR, O_INJ
      use global,             only: dirty_debug
      use named_array_list,   only: qna
      use ppp,                only: ppp_main

      implicit none

      class(cg_leaves_t),        intent(in) :: this       !< the list on which to perform the boundary exchange
      integer(kind=4),           intent(in) :: ind        !< index of cg%q(:) 3d array
      integer(kind=4), optional, intent(in) :: area_type  !< defines how do we treat boundaries
      integer(kind=4), optional, intent(in) :: bnd_type   !< Override default boundary type on external boundaries (useful in multigrid solver).
                                                          !< Note that BND_PER, BND_MPI, BND_SHE and BND_COR aren't external and cannot be overridden
      integer(kind=4), optional, intent(in) :: dir        !< select only this direction
      logical,         optional, intent(in) :: nocorners  !< .when .true. then don't care about proper edge and corner update

      type(cg_level_connected_t), pointer   :: curl
      character(len=*), parameter :: l3b_label = "leaf:arr3d_boundaries", l3bp_label = "leaf:arr3d_boundaries:prolong"
      logical :: nc

      call ppp_main%start(l3b_label)

      nc = .false.
      if (present(nocorners)) nc=nocorners
      if (qna%lst(ind)%ord_prolong /= O_INJ) nc = .false.

      curl => this%coarsest_leaves
      do while (associated(curl))
         if (dirty_debug) call curl%dirty_boundaries(ind)
         call curl%level_3d_boundaries(ind, area_type=area_type, bnd_type=bnd_type, dir=dir, nocorners=nc)
         ! corners are required on all levels except for finest anyway if prolongation is greater than injection
         curl => curl%finer
      enddo

      call ppp_main%start(l3bp_label, PPP_AMR)
      curl => this%coarsest_leaves
      do while (associated(curl))
         ! here we can use any high order prolongation without destroying conservation
         call curl%prolong_bnd_from_coarser(ind, dir=dir, nocorners=nocorners)

         if (present(bnd_type)) call curl%external_boundaries(ind, area_type=area_type, bnd_type=bnd_type)
         ! this is needed in multigrid residual to enforce strict BND_NEGREF on fine corners at domain boundaries

         curl => curl%finer
      enddo
      call ppp_main%stop(l3bp_label, PPP_AMR)

      call ppp_main%stop(l3b_label)

   end subroutine leaf_arr3d_boundaries

!> \brief This routine sets up all guardcells (internal, external and fine-coarse) for given rank-4 arrays on leaves.

   subroutine leaf_arr4d_boundaries(this, ind, area_type, dir, nocorners, no_fc)

      use cg_level_connected, only: cg_level_connected_t
      use cg_level_finest,    only: finest
      use constants,          only: O_INJ
      use named_array_list,   only: wna
      use ppp,                only: ppp_main

      implicit none

      class(cg_leaves_t),        intent(in) :: this       !< the list on which to perform the boundary exchange
      integer(kind=4),           intent(in) :: ind        !< index of cg%w(:) 4d array
      integer(kind=4), optional, intent(in) :: area_type  !< defines how do we treat boundaries
      integer(kind=4), optional, intent(in) :: dir        !< select only this direction
      logical,         optional, intent(in) :: nocorners  !< when .true. then don't care about proper edge and corner update
      logical,         optional, intent(in) :: no_fc      !< when .true. then skip prolong_bnd_from_coarser and do only intra-level internal boundaries

      type(cg_level_connected_t), pointer   :: curl
      character(len=*), parameter :: l4b_label = "leaf:arr4d_boundaries"
      logical :: nc, do_fc

      call ppp_main%start(l4b_label)

      nc = .false.
      if (present(nocorners)) nc=nocorners
      if (wna%lst(ind)%ord_prolong /= O_INJ) nc = .false.

      do_fc = .true.
      if (present(no_fc)) do_fc = .not. no_fc

      curl => finest%level
      do while (associated(curl))
         if (this%coarsest_leaves%l%id <= curl%l%id) &
              call curl%level_4d_boundaries(ind, area_type=area_type, dir=dir, nocorners=nocorners)
         curl => curl%coarser
      enddo

      if (do_fc) then
         curl => this%coarsest_leaves
         do while (associated(curl))
            call curl%prolong_bnd_from_coarser(ind, arr4d=.true., dir=dir, nocorners=nc)
            ! corners are required on all levels except for finest anyway if prolongation order is greater than injection
            curl => curl%finer
         enddo
      endif

      call ppp_main%stop(l4b_label)

   end subroutine leaf_arr4d_boundaries

!>
!! \brief Change the order of cg to optimize fine->coarse flux transfer.
!!
!! * Put these cg which have something to send, but aren't waiting for receives.
!! * Then put remaining cg which have something to send.
!! * Then put remaining cg which have any leaf (active) cells.
!! * Fully covered cg should optionally be appended at the end just in case.
!<

   function prioritized_cg(this, dir, covered_too) result(sorted_leaves)

      use cg_list,        only: cg_list_element
      use cg_list_dataop, only: cg_list_dataop_t

      implicit none

      class(cg_leaves_t), intent(in) :: this  !< object invoking type-bound procedure
      integer(kind=4),    intent(in) :: dir
      logical, optional,  intent(in) :: covered_too

      type(cg_list_dataop_t), pointer :: sorted_leaves
      type(cg_list_element),  pointer :: cgl

      allocate(sorted_leaves)
      call sorted_leaves%init_new("sorted_leaves")

      cgl => this%first
      do while (associated(cgl))
         cgl%cg%processed = .false.
         if (cgl%cg%is_sending_fc_flux(dir) .and. .not. cgl%cg%is_receiving_fc_flux(dir)) then
            call sorted_leaves%add(cgl%cg)
            cgl%cg%processed = .true.
         endif
         cgl => cgl%nxt
      enddo

      cgl => this%first
      do while (associated(cgl))
         if (.not. cgl%cg%processed) then
            if (cgl%cg%is_sending_fc_flux(dir)) then
               call sorted_leaves%add(cgl%cg)
               cgl%cg%processed = .true.
            endif
         endif
         cgl => cgl%nxt
      enddo

      cgl => this%first
      do while (associated(cgl))
         if (.not. cgl%cg%processed) then
            if (cgl%cg%has_leaves()) then
               call sorted_leaves%add(cgl%cg)
               cgl%cg%processed = .true.
            endif
         endif
         cgl => cgl%nxt
      enddo

      if (covered_too) then
         cgl => this%first
         do while (associated(cgl))
            if (.not. cgl%cg%processed) then
               call sorted_leaves%add(cgl%cg)
               cgl%cg%processed = .true.
            endif
            cgl => cgl%nxt
         enddo
      endif

   end function prioritized_cg

!> \brief Return a leaves list without fully covered cg

   function leaf_only_cg(this) result(selected_leaves)

      use cg_list,        only: cg_list_element
      use cg_list_dataop, only: cg_list_dataop_t

      implicit none

      class(cg_leaves_t), intent(in) :: this  !< object invoking type-bound procedure

      type(cg_list_dataop_t), pointer :: selected_leaves
      type(cg_list_element), pointer :: cgl

      allocate(selected_leaves)
      call selected_leaves%init_new("non-covered_leaves")

      cgl => this%first
      do while (associated(cgl))
         if (cgl%cg%has_leaves()) call selected_leaves%add(cgl%cg)
         cgl => cgl%nxt
      enddo

   end function leaf_only_cg

end module cg_leaves
