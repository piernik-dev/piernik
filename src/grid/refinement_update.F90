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

!> \brief This module contains routines for managing refinement and derefinement on the domain

module refinement_update

   implicit none

   private
   public :: update_refinement

contains

!> \brief Apply all (de)refinement criteria: user, automatic and by primitive geometric shapes

   subroutine scan_for_refinements

      use all_boundaries,        only: all_bnd, all_bnd_vital_q
      use cg_leaves,             only: leaves
      use cg_list,               only: cg_list_element
      use constants,             only: I_ONE, pSUM
      use dataio_pub,            only: msg, printinfo
      use mpisetup,              only: master, piernik_MPI_Allreduce
      use unified_ref_crit_list, only: urc_list
      use user_hooks,            only: problem_refine_derefine

      implicit none

      type(cg_list_element), pointer :: cgl
      logical, parameter :: verbose = .true.  ! be verbose for now, later we may want to be able to make it quiet
      enum, bind(C)
         enumerator :: PROBLEM
         enumerator :: URC
      end enum
      integer, dimension(PROBLEM:URC) :: cnt

      cnt = 0

      call prepare_ref

      ! We have to guarantee up-to-date guardcells on all vital fields
      call all_bnd ! \todo find a way to minimize calling this - perhaps manage a flag that says whether the boundaries are up to date or not
      call all_bnd_vital_q

      if (associated(problem_refine_derefine)) &
           call problem_refine_derefine ! call user routine first, so it cannot alter flags set by automatic routines

      if (verbose) then
         cnt(PROBLEM) = count_ref_flags()
         call piernik_MPI_Allreduce(cnt(PROBLEM), pSUM)
      endif

      call urc_list%all_mark(leaves%first)

      if (verbose) then
         cnt(URC) = count_ref_flags()
         call piernik_MPI_Allreduce(cnt(URC), pSUM)
      endif

      call parents_prevent_derefinement
      call sanitize_all_ref_flags  ! last call has to be done regardless of verbosity
      ! \todo count the change in derefinement flags here?

      if (verbose) then
         if (cnt(ubound(cnt, dim=1)) > 0) then
            write(msg,'(a,2i6,a)')"[refinement_update:scan_for_refinements] User routine and URC marked ", &
                 &                cnt(PROBLEM), cnt(PROBLEM+I_ONE:URC)-cnt(PROBLEM:URC-I_ONE), " block(s) for refinement, respectively."
            if (master) call printinfo(msg)
         endif
      endif

      ! \todo convert this to URC routine
      cgl => leaves%first
      do while (associated(cgl))
         call cgl%cg%refinemap2SFC_list  ! modifies cg%flag%map so no parent correction is possible afterwards
         cgl => cgl%nxt
      enddo

   contains

      !> \brief Initialize refinement flags and other data structures

      subroutine prepare_ref

         use cg_level_connected, only: cg_level_connected_t
         use cg_level_finest,    only: finest
         use cg_list,            only: cg_list_element

         implicit none

         type(cg_level_connected_t), pointer :: curl
         type(cg_list_element), pointer :: cgl

         curl => finest%level
         do while (associated(curl))
            call curl%deallocate_patches

            ! formerly cg_list::clear_ref_flags
            cgl => curl%first
            do while (associated(cgl))
               call cgl%cg%flag%init

               ! Mark everything for derefinement by default.
               ! It requires correct propagation of refinement requests from parent blocks as
               ! derefinement inhibitions on their appropriate children.
               cgl%cg%flag%derefine = .true.

               cgl => cgl%nxt
            enddo

            curl => curl%coarser
         enddo

      end subroutine prepare_ref

      !>
      !! \brief Sanitize refinement requests
      !!
      !! \todo convert this to URC routine
      !<

      subroutine sanitize_all_ref_flags

         use cg_leaves, only: leaves
         use cg_list,   only: cg_list_element

         implicit none

         type(cg_list_element), pointer :: cgl

         !> \todo communicate refines from coarse to fine blocks to prevent oscillations that might occur when there is derefinement request on fine, when coarse requests refinement

         cgl => leaves%first
         do while (associated(cgl))
            associate (cg => cgl%cg)
               cg%flag%refine = cg%flag%refine .or. cg%flag%get(int(cg%ijkse, kind=8), cg%leafmap)
               call cg%flag%sanitize(cg%l%id)
            end associate
            cgl => cgl%nxt
         enddo

      end subroutine sanitize_all_ref_flags

      !> \brief Count refinement flags everywhere

      function count_ref_flags() result(cnt)

         use cg_leaves, only: leaves
         use cg_list,   only: cg_list_element

         implicit none

         integer :: cnt !< returned counter

         type(cg_list_element), pointer :: cgl

         call sanitize_all_ref_flags

         cnt = 0
         cgl => leaves%first
         do while (associated(cgl))
            if (cgl%cg%flag%refine) cnt = cnt + 1
            cgl => cgl%nxt
         enddo

      end function count_ref_flags

   end subroutine scan_for_refinements

!>
!! \brief Clear derefine flag on child whenever there is a refine flag below, on its parent.
!!
!! Here we do loop over levels because we use mpi tags for prolongation to identify target grids.
!! These tags are guaranteed to be unique for a given pair of consecutive levels but not for the whole domain.
!! Once we implement a better structure to identify parent-child relation, we should switch to leaves-based list here (as it was implemented in 1a048b08).
!! Also consolidation of MPI communication is advised.
!<

   subroutine parents_prevent_derefinement

      use cg_level_base,      only: base
      use cg_level_connected, only: cg_level_connected_t

      implicit none

      type(cg_level_connected_t), pointer :: curl

      curl => base%level
      do while (associated(curl))
         if (associated(curl%finer)) call parents_prevent_derefinement_lev(curl)
         curl => curl%finer
      enddo

   end subroutine parents_prevent_derefinement

!> \brief Clear derefine flag on child whenever there is a refine flag below, on its parent (single-level)

   subroutine parents_prevent_derefinement_lev(lev)

      use cg_level_connected, only: cg_level_connected_t
      use cg_list,            only: cg_list_element
      use constants,          only: xdim, zdim, LO, HI, I_ONE, I_ZERO
      use dataio_pub        , only: die, warn
      use domain,             only: dom
      use mpi,                only: MPI_INTEGER, MPI_STATUS_IGNORE
      use mpisetup,           only: FIRST, LAST, comm, mpi_err, proc, req, status, inflate_req

      implicit none

      type(cg_level_connected_t), pointer, intent(in) :: lev

      type(cg_list_element), pointer :: cgl
      integer :: i, g, nr
      integer(kind=8), dimension(xdim:zdim, LO:HI) :: se
      integer, parameter :: perimeter = 2
      integer(kind=4), dimension(FIRST:LAST) :: gscnt, grcnt
      type :: pt  ! (proc, tag) pair
         integer(kind=4) :: proc
         integer(kind=4) :: tag  ! or grid_id, if tag fails
      end type pt
      type(pt), dimension(:), allocatable :: pt_list
      integer :: pt_cnt
      integer(kind=4) :: rtag
      integer(kind=4), dimension(:,:), pointer :: mpistatus

      if (perimeter > dom%nb) call die("[refinement_update:parents_prevent_derefinement_lev] perimeter > dom%nb")
      if (.not. associated(lev%finer)) then
         call warn("[refinement_update:parents_prevent_derefinement_lev] .not. associated(lev%finer)")
         return
      endif

      gscnt = 0
      allocate(pt_list(lev%cnt * 2**dom%eff_dim))  ! ToDo make it more flexible
      pt_cnt = 0
      cgl => lev%first
      do while (associated(cgl))
         associate (cg => cgl%cg)
            ! look for refineflag .and. .not. leafmap + preimeter
            if (cg%flag%get(int(cg%ijkse, kind=8), .not. cg%leafmap)) then
               if (allocated(cg%ri_tgt%seg)) then
                  do g = lbound(cg%ri_tgt%seg(:), dim=1), ubound(cg%ri_tgt%seg(:), dim=1)
                     ! clear own derefine flags (single-thread test)
                     se(:, LO) = cg%ri_tgt%seg(g)%se(:, LO) - perimeter * dom%D_
                     se(:, HI) = cg%ri_tgt%seg(g)%se(:, HI) + perimeter * dom%D_
                     if (cg%flag%get(se)) then

!!$                     if (associated(cg%ri_tgt%seg(g)%local)) then
!!$                        cg%ri_tgt%seg(g)%local%flag%derefine = .false.
!!$                     else
                        associate (fproc => int(cg%ri_tgt%seg(g)%proc, kind=4), ftag => int(cg%ri_tgt%seg(g)%tag, kind=4))
                           if (fproc == proc) then
                              call disable_derefine_by_tag(lev%finer, ftag)  ! beware: O(leaves%cnt^2)
                           else
                              ! create a list of foreign blocks that need not be derefined (proc, grid_id or SFC_id)
                              ! use prolongation structures
                              gscnt(fproc) = gscnt(fproc) + I_ONE
                              pt_cnt = pt_cnt + 1
                              if (pt_cnt > size(pt_list)) call die("[refinement_update:parents_prevent_derefinement_lev] pt_cnt > size(pt_list)")
                              pt_list(pt_cnt) = pt(fproc, ftag)
                           endif
                        end associate
                     endif
                  enddo
               endif
            endif
         end associate
         cgl => cgl%nxt

      enddo

      ! communicate how many ids to which threads
      call MPI_Alltoall(gscnt, I_ONE, MPI_INTEGER, grcnt, I_ONE, MPI_INTEGER, comm, mpi_err)

      ! Apparently gscnt/grcnt represent quite sparse matrix, so we better do nonblocking point-to-point than MPI_AlltoAllv
      nr = 0
      if (pt_cnt > 0) then
         do g = lbound(pt_list, dim=1), pt_cnt
            nr = nr + 1
            if (nr > size(req, dim=1)) call inflate_req
            call MPI_Isend(pt_list(g)%tag, I_ONE, MPI_INTEGER, pt_list(g)%proc, I_ZERO, comm, req(nr), mpi_err)
            ! OPT: Perhaps it will be more efficient to allocate arrays according to gscnt and send tags in bunches
         enddo
      endif

      do g = lbound(grcnt, dim=1), ubound(grcnt, dim=1)
         if (grcnt(g) /= 0) then
            if (g == proc) call die("[refinement_update:parents_prevent_derefinement] MPI_Recv from self")  ! this is not an error but it should've been handled as locat thing
            do i = 1, grcnt(g)
               call MPI_Recv(rtag, I_ONE, MPI_INTEGER, g, I_ZERO, comm, MPI_STATUS_IGNORE, mpi_err)
               call disable_derefine_by_tag(lev%finer, rtag)  ! beware: O(leaves%cnt^2)
            enddo
         endif
      enddo

      if (nr > 0) then
         mpistatus => status(:, :nr)
         call MPI_Waitall(nr, req(:nr), mpistatus, mpi_err)
      endif

      deallocate(pt_list)

   contains

      !>
      !! \brief Disable derefinement on child if we know only the tag for outgoing restriction
      !!
      !! Beware: slow: O(leaves%cnt^2)
      !<

      subroutine disable_derefine_by_tag(flev, tag)

         use cg_level_connected, only: cg_level_connected_t

         implicit none

         type(cg_level_connected_t), pointer, intent(in) :: flev

         integer(kind=4), intent(in) :: tag

         integer :: g
         logical :: done
         type(cg_list_element), pointer :: cgl

         done = .false.

         cgl => flev%first
         do while (associated(cgl) .and. .not. done)
            if (allocated(cgl%cg%ro_tgt%seg)) then
               do g = lbound(cgl%cg%ro_tgt%seg(:), dim=1), ubound(cgl%cg%ro_tgt%seg(:), dim=1)
                  if (cgl%cg%ro_tgt%seg(g)%tag == tag) then
                     cgl%cg%flag%derefine = .false.
                     done = .true.
                     exit
                  endif
               enddo
            endif
            cgl => cgl%nxt
         enddo

      end subroutine disable_derefine_by_tag

   end subroutine parents_prevent_derefinement_lev

!> \brief Update the refinement topology

!#define DEBUG_DUMPS

   subroutine update_refinement(act_count, refinement_fixup_only)

      use all_boundaries,        only: all_bnd, all_bnd_vital_q
      use cg_leaves,             only: leaves
      use cg_list,               only: cg_list_element
      use cg_level_base,         only: base
      use cg_level_connected,    only: cg_level_connected_t
      use cg_level_finest,       only: finest
      use cg_list_global,        only: all_cg
      use constants,             only: pLOR, pLAND, pSUM, tmr_amr
      use dataio_pub,            only: warn, die
      use global,                only: nstep
      use grid_cont,             only: grid_container
      use list_of_cg_lists,      only: all_lists
      use mpisetup,              only: piernik_MPI_Allreduce!, proc
      use refinement,            only: n_updAMR, emergency_fix
      use timer,                 only: set_timer
      use unified_ref_crit_list, only: urc_list
#ifdef GRAV
      use gravity,               only: update_gp
#endif /* GRAV */
#ifdef DEBUG_DUMPS
      use data_hdf5,             only: write_hdf5
#endif /* DEBUG_DUMPS */

      implicit none

      integer, optional, intent(out) :: act_count             !< counts number of blocks refined or deleted
      logical, optional, intent(in)  :: refinement_fixup_only !< When present and .true. then do not check refinement criteria, do only correction, if necessary.

      integer :: nciter
      integer, parameter :: nciter_max = 100 ! should be more than refinement levels
      logical :: some_refined, derefined
      type(cg_list_element), pointer :: cgl, aux
      type(cg_level_connected_t), pointer :: curl
      type(grid_container),  pointer :: cg
      logical :: correct, full_update
      real :: ts  !< time for runtime profiling

      ts =  set_timer(tmr_amr, .true.)

      if (present(act_count)) act_count = 0

      full_update = .true.
      if (present(refinement_fixup_only)) full_update = .not. refinement_fixup_only

      call finest%level%restrict_to_base ! implies update of leafmap

      if (n_updAMR <= 0) then
         ! call print_time("[refinement_update] No refinement update")  ! ToDo: enable only when detailed profiling is required
         return
      endif
      if (mod(nstep, n_updAMR) /= 0 .and. .not. emergency_fix) then
         ! call print_time("[refinement_update] No refinement update")  ! ToDo: enable only when detailed profiling is required
         return
      endif

      call all_cg%prevent_prolong
      call all_cg%set_is_old
      curl => finest%level
      do while (associated(curl))
         curl%recently_changed = .false.
         curl => curl%coarser
      enddo

      if (full_update) then
         call scan_for_refinements

         ! do the refinements first
         curl => base%level
         do while (associated(curl))
            cgl => curl%first
            do while (associated(cgl))
               if (any(cgl%cg%leafmap)) then
                  if (cgl%cg%flag%refine) then
                     call refine_one_grid(curl, cgl)
                     if (present(act_count)) act_count = act_count + 1
                     cgl%cg%flag%refine = .false.
                  endif
               endif
               cgl => cgl%nxt
            enddo

            call finest%equalize

            !> \todo merge small blocks into larger ones
            if (associated(curl%finer)) then
               call curl%finer%init_all_new_cg
               call curl%finer%deallocate_patches
               call curl%finer%sync_ru
               call curl%prolong
            endif

            curl => curl%finer
         enddo

         call all_cg%mark_orphans

         call leaves%update(" (  refine  ) ")
      endif

      ! fix the structures, mark grids for refinement (and unmark derefinements) due to refinement restrictions.
      ! With clever use of SFC properties this part can be done together with refine stage above - everything can be determined in sanitize_all_ref_flags for regular refinement update.
      ! For refinement update due to domain expansion, everything can be fixed in the expansion routine.
      correct = fix_refinement()
      call piernik_MPI_Allreduce(correct, pLAND)
      nciter = 0
      do while (.not. correct)

         call all_cg%prevent_prolong

         curl => finest%level
         do while (associated(curl))
            call curl%deallocate_patches
            curl => curl%coarser
         enddo

         curl => finest%level%coarser
         do while (associated(curl) .and. curl%l%id >= base%level%l%id)
            some_refined = .false.
            cgl => curl%first
            do while (associated(cgl))
               if (cgl%cg%flag%refine) then
                  if (finest%level%l%id <= cgl%cg%l%id) call warn("[refinement_update:update_refinement] growing too fast!")
                  if (associated(curl%finer)) then
                     call refine_one_grid(curl, cgl)
                     if (present(act_count)) act_count = act_count + 1
                     some_refined = .true.
                  else
                     call warn("[refinement_update:update_refinement] nowhere to add!")
                  endif
                  cgl%cg%flag%refine = .false.
               endif
               cgl => cgl%nxt
            enddo

            call piernik_MPI_Allreduce(some_refined, pLOR)
            if (some_refined) then
               call curl%finer%init_all_new_cg
               call curl%finer%sync_ru
               call curl%deallocate_patches
               !call finest%equalize
               call curl%prolong
               call all_cg%mark_orphans
            endif

            curl => curl%coarser
         enddo

         ! sync structure before trying to fix it
         call leaves%update(" (correcting) ")
         ! \todo implement 1-pass correcting
         correct = fix_refinement()
         call piernik_MPI_Allreduce(correct, pLAND)
         if (.not. correct) nciter = nciter + 1
         if (nciter > nciter_max) then
            call die("[refinement_update:update_refinement] Cannot converge to correct structure. Too many iterations.")
            correct = .true.
         endif
      enddo

      ! Now try to derefine any excess of refinement.
      ! Derefinement saves memory and CPU usage, but it is of lowest priority here.
      ! Just go through derefinement stage once and don't try to do too much here at once.
      ! Any excess of refinement will be handled in next call to this routine anyway.
      if (full_update) then

         curl => finest%level
         do while (associated(curl) .and. .not. associated(curl, base%level))
            cgl => curl%first
            derefined = .false.
            do while (associated(cgl))
               aux => cgl ! Auxiliary pointer makes it easier to loop over the list of grids when some of the elements are disappearing
               cgl => cgl%nxt
               if (aux%cg%flag%derefine) then
                  if (aux%cg%is_old) then ! forgetting fresh blocks is not good because their data hasn't been properly initialized yet
                     if (all(aux%cg%leafmap)) then
                        cg => aux%cg
                        call all_lists%forget(cg)
                        if (present(act_count)) act_count = act_count + 1
                        curl%recently_changed = .true.
                        derefined = .true.
                     endif
                  endif
               endif
            enddo
            call piernik_MPI_Allreduce(derefined, pLOR)
            if (derefined) call curl%init_all_new_cg ! no new cg to really initialize, but some other routines need to be called to refresh datastructures
            !> \todo replace call curl%init_all_new_cg by something cheaper
            call curl%sync_ru
            curl => curl%coarser
         enddo
         ! sync structure
         call leaves%balance_and_update(" ( derefine ) ")
         !call fix_refinement

      endif

      ! Check refinement topology and crash if anything got broken
      if (.not. fix_refinement()) call die("[refinement_update:update_refinement] Refinement defects still present")

      call all_bnd

      call urc_list%plot_mark(leaves%first)

      !> \todo call the update of cs_i2 and other vital variables if and only if something has changed
      !> \todo add another flag to named_array_list::na_var so the user can also specify fields that need boundary updates on fine/coarse boundaries
      call all_bnd_vital_q
#ifdef GRAV
      call update_gp
#endif /* GRAV */

      call all_cg%enable_prolong
      if (present(act_count)) call piernik_MPI_Allreduce(act_count, pSUM)

#ifdef DEBUG_DUMPS
      call write_hdf5
#endif /* DEBUG_DUMPS */

      emergency_fix = .false.

      if (full_update) then
         call print_time("[refinement_update] Finishing (full update)")
      else
         call print_time("[refinement_update] Finishing (fixup only)")
      endif

   contains

      subroutine print_time(str)

         use constants,  only: tmr_amr
         use dataio_pub, only: msg, printinfo
         use mpisetup,   only: master
         use timer,      only: set_timer

         implicit none

         character(len=*), intent(in) :: str

         write(msg, '(2a,f7.3)') str, ", dt_wall= ", set_timer(tmr_amr)
         if (master) call printinfo(msg)

      end subroutine print_time

   end subroutine update_refinement

!> \brief Refine a single grid piece. Pay attention whether it is already refined

   subroutine refine_one_grid(curl, cgl)

      use cg_level_connected, only: cg_level_connected_t
      use cg_level_finest,    only: finest
      use cg_list,            only: cg_list_element
      use refinement,         only: bsize
      use dataio_pub,         only: warn, die

      implicit none

      type(cg_level_connected_t), pointer, intent(inout) :: curl
      type(cg_list_element),      pointer, intent(in)    :: cgl

      integer :: b

      if (.not. any(cgl%cg%leafmap)) then
         call warn("[refinement_update:refine_one_grid] Attempting to refine a grid that is completely refined")
         return
      endif

      if (.not. associated(curl%finer)) call finest%add_finer
      if (size(cgl%cg%flag%SFC_refine_list) > 0) then ! we require detailed map!
         do b = lbound(cgl%cg%flag%SFC_refine_list, dim=1), ubound(cgl%cg%flag%SFC_refine_list, dim=1)
            if (cgl%cg%flag%SFC_refine_list(b)%level == curl%finer%l%id) then
               call curl%finer%add_patch(int(bsize, kind=8), cgl%cg%flag%SFC_refine_list(b)%off)
            else
               call die("[refinement_update:refine_one_grid] wrong level!")
            endif
         enddo
         call cgl%cg%flag%init ! it is safer to forget it now
      else
         call die("[refinement_update:refine_one_grid] populating flag%SFC_refine_list before refining a grid is strongly encouragrequired")
      endif

   end subroutine refine_one_grid

!>
!! \brief Apply some rules to fix refinement defects
!!
!! \details Important! Requires updated cg_leafmap!
!<

   logical function fix_refinement() result(correct)

      use cg_level_finest,  only: finest
      use cg_list,          only: cg_list_element
!      use cg_list_global,   only: all_cg
      use cg_leaves,        only: leaves
      use constants,        only: xdim, ydim, zdim, I_ONE!, I_TWO
      use dataio_pub,       only: die
      use domain,           only: dom
      use named_array_list, only: qna
!      use cg_level_base, only: base

      implicit none

      type(cg_list_element), pointer :: cgl
      integer(kind=4) ::  i, j, k
      integer :: lleaf, lnear, range
      enum, bind(C)
         enumerator :: INSIDE   = -1
         enumerator :: OUTSIDE  =  0
         enumerator :: BOUNDARY =  1
      end enum

      correct = .true.

      !> \todo check for excess of refinement levels

      ! Put a level number to the working array, restrict it and exchange internal boundaries
      cgl => leaves%first
      do while (associated(cgl))
         cgl%cg%wa = cgl%cg%l%id
         cgl => cgl%nxt
      enddo
      call finest%level%restrict_to_base_q_1var(qna%wai)

!!$      curl => base%level
!!$      do while (associated(curl))
!!$         call curl%clear_boundaries(ind, value=10.)
!!$         call curl%prolong_bnd_from_coarser(ind)
!!$         call curl%level_3d_boundaries(ind, nb, area_type, bnd_type, corners)
!!$!         call curl%arr3d_boundaries(ind, nb, area_type, bnd_type, corners)
!!$         curl => curl%finer
!!$      enddo
      call leaves%leaf_arr3d_boundaries(qna%wai)

      range = 1
      ! range = min((dom%nb+I_ONE)/I_TWO + all_cg%ord_prolong_nb, dom%nb)
      ! (dom%nb+1)/2 + all_cg%ord_prolong_nb is a range of influence of coarse to fine levels - it can be suitable if we were looking for too low levels in the neighbourhood
      ! ATM we're looking for high levels, so range = 1 seems to be appropriate

      ! detect high refinements near leafmap and alter refinement flags appropriately
      cgl => leaves%first
      do while (associated(cgl))

         ! this is perhaps a bit suboptimal: todo let cgl%cg%flag%SFC_refine_list(:) grow and implement INVALID entries
         if (allocated(cgl%cg%flag%SFC_refine_list)) deallocate(cgl%cg%flag%SFC_refine_list)
         allocate(cgl%cg%flag%SFC_refine_list(0))
         call cgl%cg%flag%clear

         if (any(cgl%cg%leafmap)) then
            cgl%cg%prolong_xyz = OUTSIDE
            lleaf = -huge(1)

            ! find the level of refinement on the leaf part and mark it with a negative value
            do k = lbound(cgl%cg%leafmap, dim=3, kind=4), ubound(cgl%cg%leafmap, dim=3, kind=4)
               do j = lbound(cgl%cg%leafmap, dim=2, kind=4), ubound(cgl%cg%leafmap, dim=2, kind=4)
                  do i = lbound(cgl%cg%leafmap, dim=1, kind=4), ubound(cgl%cg%leafmap, dim=1, kind=4)
                     if (cgl%cg%leafmap(i, j, k)) then
                        if (lleaf == -huge(1)) then
                           lleaf = int(cgl%cg%wa(i, j, k))
                        else
                           if (lleaf /= int(cgl%cg%wa(i, j, k)) .or. lleaf /= cgl%cg%l%id) call die("[refinement_update:fix_refinement] Inconsistent level map")
                        endif
                        cgl%cg%prolong_xyz(i,j,k) = INSIDE
                     endif
                  enddo
               enddo
            enddo

            ! find the border of the leaf map and mark it with positive value
            do k = lbound(cgl%cg%leafmap, dim=3, kind=4)-dom%D_z, ubound(cgl%cg%leafmap, dim=3, kind=4)+dom%D_z
               do j = lbound(cgl%cg%leafmap, dim=2, kind=4)-dom%D_y, ubound(cgl%cg%leafmap, dim=2, kind=4)+dom%D_y
                  do i = lbound(cgl%cg%leafmap, dim=1, kind=4)-dom%D_x, ubound(cgl%cg%leafmap, dim=1, kind=4)+dom%D_x
                     if (int(cgl%cg%prolong_xyz(i, j, k)) == OUTSIDE) then
                        if (dom%has_dir(xdim)) then
                           if (int(cgl%cg%prolong_xyz(i+1, j, k)) /= 0) cgl%cg%prolong_xyz(i+1, j, k) = BOUNDARY
                           if (int(cgl%cg%prolong_xyz(i-1, j, k)) /= 0) cgl%cg%prolong_xyz(i-1, j, k) = BOUNDARY
                        endif
                        if (dom%has_dir(ydim)) then
                           if (int(cgl%cg%prolong_xyz(i, j+1, k)) /= 0) cgl%cg%prolong_xyz(i, j+1, k) = BOUNDARY
                           if (int(cgl%cg%prolong_xyz(i, j-1, k)) /= 0) cgl%cg%prolong_xyz(i, j-1, k) = BOUNDARY
                        endif
                        if (dom%has_dir(zdim)) then
                           if (int(cgl%cg%prolong_xyz(i, j, k+1)) /= 0) cgl%cg%prolong_xyz(i, j, k+1) = BOUNDARY
                           if (int(cgl%cg%prolong_xyz(i, j, k-1)) /= 0) cgl%cg%prolong_xyz(i, j, k-1) = BOUNDARY
                        endif
                     endif
                  enddo
               enddo
            enddo

            !find highest refinement level witin a "range" cells from the border of the leaf map
! dirty check
!!$            if (maxval(cgl%cg%wa) > 10) then
!!$               do k = lbound(cgl%cg%wa, dim=3), ubound(cgl%cg%wa, dim=3)
!!$                  do j = lbound(cgl%cg%wa, dim=2), ubound(cgl%cg%wa, dim=2)
!!$                     do i = lbound(cgl%cg%wa, dim=1), ubound(cgl%cg%wa, dim=1)
!!$                        if (cgl%cg%wa(i,j,k) > 10) write(msg,*)"ur:fr WA(",i,j,k,") = ",cgl%cg%wa(i,j,k)
!!$                     enddo
!!$                  enddo
!!$               enddo
!!$            endif

            do k = lbound(cgl%cg%leafmap, dim=3, kind=4), ubound(cgl%cg%leafmap, dim=3, kind=4)
               do j = lbound(cgl%cg%leafmap, dim=2, kind=4), ubound(cgl%cg%leafmap, dim=2, kind=4)
                  do i = lbound(cgl%cg%leafmap, dim=1, kind=4), ubound(cgl%cg%leafmap, dim=1, kind=4)
                     if (int(cgl%cg%prolong_xyz(i, j, k)) == BOUNDARY) then
                        lnear = int(min(huge(I_ONE)/10.,maxval(cgl%cg%wa(i-range*dom%D_x:i+range*dom%D_x, j-range*dom%D_y:j+range*dom%D_y, k-range*dom%D_z:k+range*dom%D_z))))
                        if (lnear > lleaf) then
                           cgl%cg%flag%derefine = .false.
                           if (lnear > lleaf+1 .and. lnear <= finest%level%l%id) then
                              call cgl%cg%flag%set(i, j, k)
                              correct = .false.
                           endif
                        endif

                     endif
                  enddo
               enddo
            enddo
            call cgl%cg%refinemap2SFC_list
         endif

         call cgl%cg%flag%sanitize(cgl%cg%l%id)

         cgl => cgl%nxt
      enddo

   end function fix_refinement

end module refinement_update
