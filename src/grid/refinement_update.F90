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
      use cg_cost_data,          only: I_REFINE
      use cg_leaves,             only: leaves
      use cg_list,               only: cg_list_element
      use constants,             only: I_ONE, pSUM, PPP_AMR, PPP_PROB
      use dataio_pub,            only: msg, printinfo
      use mpisetup,              only: master, piernik_MPI_Allreduce
      use ppp,                   only: ppp_main
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
      character(len=*), parameter :: scan_ref_label = "scan_for_refinement", u_ref_label = "user_scan_for_refinement", urc_label = "URC_scan_for_refinement"

      call ppp_main%start(scan_ref_label, PPP_AMR)
      cnt = 0
      call prepare_ref

      ! We have to guarantee up-to-date guardcells on all vital fields
      call all_bnd ! \todo find a way to minimize calling this - perhaps manage a flag that says whether the boundaries are up to date or not
      call all_bnd_vital_q

      ! Call user routine first, so it cannot alter flags set by automatic routines
      if (associated(problem_refine_derefine)) then
         call ppp_main%start(u_ref_label, PPP_AMR + PPP_PROB)
         call problem_refine_derefine
         call ppp_main%stop(u_ref_label, PPP_AMR + PPP_PROB)
      endif

      if (verbose) then
         cnt(PROBLEM) = count_ref_flags()
         call piernik_MPI_Allreduce(cnt(PROBLEM), pSUM)
      endif

      ! Apply all Unified Refinement Criteria (everything that works locally, cg-wise)
      call ppp_main%start(urc_label, PPP_AMR)
      call urc_list%all_mark(leaves%first)
      call ppp_main%stop(urc_label, PPP_AMR)

      if (verbose) then
         cnt(URC) = count_ref_flags()
         call piernik_MPI_Allreduce(cnt(URC), pSUM)
      endif

      ! ToDo: exploit dot for making corrections
      ! if dot is complete then correction can be done in one step

      ! Transform refinement requests on non-leaf parts of the grid into derefinement inhibitions on child grids to avoid refinement flickering
      call parents_prevent_derefinement
      call sanitize_ref_flags  ! it can be converted to URC routine but it must be guaranteed that it goes after all refinement mark routines
      ! \todo count the change in derefinement flags here too?

      if (verbose) then
         if (cnt(ubound(cnt, dim=1)) > 0) then
            write(msg,'(a,2i6,a)')"[refinement_update:scan_for_refinements] User routine and URC marked ", &
                 &                cnt(PROBLEM), cnt(PROBLEM+I_ONE:URC)-cnt(PROBLEM:URC-I_ONE), " block(s) for refinement, respectively."
            if (master) call printinfo(msg)
         endif
      endif

      ! it can be converted to URC routine but it must be guaranteed that it goes after parents_prevent_derefinement because it modifies cg%flag%map
      cgl => leaves%first
      do while (associated(cgl))
         call cgl%cg%costs%start

         call cgl%cg%refinemap2SFC_list

         call cgl%cg%costs%stop(I_REFINE)
         cgl => cgl%nxt
      enddo
      call ppp_main%stop(scan_ref_label, PPP_AMR)

   contains

      !> \brief Initialize refinement flags and other data structures

      subroutine prepare_ref

         use cg_cost_data,       only: I_REFINE
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
               call cgl%cg%costs%start

               call cgl%cg%flag%init

               ! Mark everything for derefinement by default.
               ! It requires correct propagation of refinement requests from parent blocks as derefinement inhibitions on their appropriate children later
               cgl%cg%flag%derefine = .true.

               call cgl%cg%costs%stop(I_REFINE)
               cgl => cgl%nxt
            enddo

            curl => curl%coarser
         enddo

      end subroutine prepare_ref

      !> \brief Sanitize refinement requests

      subroutine sanitize_ref_flags

         use cg_cost_data, only: I_REFINE
         use cg_leaves,    only: leaves
         use cg_list,      only: cg_list_element
         use refinement,   only: level_min, level_max

         implicit none

         type(cg_list_element), pointer :: cgl

         cgl => leaves%first
         do while (associated(cgl))
            call cgl%cg%costs%start

            associate (cg => cgl%cg)
               if (cg%flag%get(int(cg%ijkse, kind=8), cg%leafmap)) cg%flag%derefine = .false.
               if (cg%l%id >= level_max) call cg%flag%clear
               if (cg%l%id <  level_min) call cg%flag%set
               if (cg%l%id >  level_max) cg%flag%derefine = .true.
               if (cg%l%id <= level_min) cg%flag%derefine = .false.
            end associate

            call cgl%cg%costs%stop(I_REFINE)
            cgl => cgl%nxt
         enddo

      end subroutine sanitize_ref_flags

      !> \brief Count refinement flags everywhere

      integer function count_ref_flags() result(cnt)

         use cg_cost_data, only: I_REFINE
         use cg_leaves,    only: leaves
         use cg_list,      only: cg_list_element

         implicit none

         type(cg_list_element), pointer :: cgl

         call sanitize_ref_flags

         cnt = 0
         cgl => leaves%first
         do while (associated(cgl))
            call cgl%cg%costs%start

            if (cgl%cg%flag%get(int(cgl%cg%ijkse, kind=8), cgl%cg%leafmap)) cnt = cnt + 1

            call cgl%cg%costs%stop(I_REFINE)
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
!! Also consolidation of MPI communication is advised (OPT).
!<

   subroutine parents_prevent_derefinement

      use cg_level_base,      only: base
      use cg_level_connected, only: cg_level_connected_t
      use constants,          only: PPP_AMR
      use ppp,                only: ppp_main

      implicit none

      type(cg_level_connected_t), pointer :: curl
      character(len=*), parameter :: ppd_label = "refinement_parent_protection"

      curl => base%level
      do while (associated(curl))
         if (associated(curl%finer)) then
            call ppp_main%start(ppd_label, PPP_AMR)
            call parents_prevent_derefinement_lev(curl)
            call ppp_main%stop(ppd_label, PPP_AMR)
         endif
         curl => curl%finer
      enddo

   end subroutine parents_prevent_derefinement

!> \brief Clear derefine flag on child whenever there is a refine flag below, on its parent (single-level)

   subroutine parents_prevent_derefinement_lev(lev)

      use cg_cost_data,       only: I_REFINE
      use cg_level_connected, only: cg_level_connected_t
      use cg_list,            only: cg_list_element
      use constants,          only: xdim, zdim, LO, HI, I_ONE, I_ZERO
      use dataio_pub,         only: die, warn
      use domain,             only: dom
      use MPIF,               only: MPI_INTEGER, MPI_STATUS_IGNORE, MPI_COMM_WORLD
      use MPIFUN,             only: MPI_Alltoall, MPI_Isend, MPI_Recv, MPI_Comm_dup, MPI_Comm_free
      use mpisetup,           only: FIRST, LAST, err_mpi, proc, req, inflate_req
      use ppp_mpi,            only: piernik_Waitall
#ifdef MPIF08
      use MPIF,       only: MPI_Comm
#endif /* MPIF08 */

      implicit none

      type(cg_level_connected_t), pointer, intent(in) :: lev

      type(cg_list_element), pointer :: cgl
      integer :: i
      integer(kind=4) :: nr, g
      integer(kind=8), dimension(xdim:zdim, LO:HI) :: se
      integer, parameter :: perimeter = 2
      integer(kind=4), dimension(FIRST:LAST) :: gscnt, grcnt
      type :: pt
         integer(kind=4) :: proc
         integer(kind=4) :: tag
      end type pt
      type(pt), dimension(:), allocatable :: pt_list
      integer(kind=4) :: pt_cnt
      integer(kind=4) :: rtag
#ifdef MPIF08
      type(MPI_Comm)  :: ref_comm
#else /* !MPIF08 */
      integer(kind=4) :: ref_comm
#endif /* !MPIF08 */

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
            call cg%costs%start

            if (cg%flag%get()) then
               if (allocated(cg%ri_tgt%seg)) then
                  do g = lbound(cg%ri_tgt%seg(:), dim=1, kind=4), ubound(cg%ri_tgt%seg(:), dim=1, kind=4)
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
                              ! here we use prolongation structures
                              gscnt(fproc) = gscnt(fproc) + I_ONE
                              pt_cnt = pt_cnt + I_ONE
                              if (pt_cnt > size(pt_list)) call die("[refinement_update:parents_prevent_derefinement_lev] pt_cnt > size(pt_list)")
                              pt_list(pt_cnt) = pt(fproc, ftag)
                           endif
                        end associate
                     endif
                  enddo
               endif
            endif

            call cg%costs%stop(I_REFINE)
         end associate
         cgl => cgl%nxt

      enddo

      ! communicate how many ids to which threads
      call MPI_Alltoall(gscnt, I_ONE, MPI_INTEGER, grcnt, I_ONE, MPI_INTEGER, MPI_COMM_WORLD, err_mpi)

      ! Apparently gscnt/grcnt represent quite sparse matrix, so we better do nonblocking point-to-point than MPI_AlltoAllv
      call MPI_Comm_dup(MPI_COMM_WORLD, ref_comm, err_mpi)
      nr = 0
      if (pt_cnt > 0) then
         do g = lbound(pt_list, dim=1, kind=4), pt_cnt
            nr = nr + I_ONE
            if (nr > size(req, dim=1)) call inflate_req
            call MPI_Isend(pt_list(g)%tag, I_ONE, MPI_INTEGER, pt_list(g)%proc, I_ZERO, ref_comm, req(nr), err_mpi)
            ! OPT: Perhaps it will be more efficient to allocate arrays according to gscnt and send tags in bunches
         enddo
      endif

      do g = lbound(grcnt, dim=1, kind=4), ubound(grcnt, dim=1, kind=4)
         if (grcnt(g) /= 0) then
            if (g == proc) call die("[refinement_update:parents_prevent_derefinement] MPI_Recv from self")  ! this is not an error but it should've been handled as local thing
            do i = 1, grcnt(g)
               call MPI_Recv(rtag, I_ONE, MPI_INTEGER, g, I_ZERO, ref_comm, MPI_STATUS_IGNORE, err_mpi)
               call disable_derefine_by_tag(lev%finer, rtag)  ! beware: O(leaves%cnt^2)
            enddo
         endif
      enddo

      call piernik_Waitall(nr, "prevent_derefinement")
      call MPI_Comm_free(ref_comm, err_mpi)

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
         integer(kind=4),                     intent(in) :: tag

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

!> \brief PPP-aware wrapper for update the refinement topology

   subroutine update_refinement(act_count, refinement_fixup_only)

      use constants, only: PPP_AMR
      use ppp,       only: ppp_main

      implicit none

      integer, optional, intent(out) :: act_count             !< counts number of blocks refined or deleted
      logical, optional, intent(in)  :: refinement_fixup_only !< When present and .true. then do not check refinement criteria, do only correction, if necessary.

      character(len=*), parameter :: ref_label = "refinement"

      call ppp_main%start(ref_label, PPP_AMR)
      call update_refinement_wrapped(act_count, refinement_fixup_only)
      call ppp_main%stop(ref_label, PPP_AMR)

   end subroutine update_refinement

!> \brief Update the refinement topology

   subroutine update_refinement_wrapped(act_count, refinement_fixup_only)

      use all_boundaries,        only: all_bnd, all_bnd_vital_q
      use cg_leaves,             only: leaves
      use cg_list,               only: cg_list_element
      use cg_level_base,         only: base
      use cg_level_connected,    only: cg_level_connected_t
      use cg_level_finest,       only: finest
      use cg_list_global,        only: all_cg
      use constants,             only: pLOR, pLAND, pSUM, pMAX, tmr_amr, PPP_AMR
      use dataio_pub,            only: die
      use global,                only: nstep
      use grid_cont,             only: grid_container
      use list_of_cg_lists,      only: all_lists
      use mpisetup,              only: piernik_MPI_Allreduce
      use ppp,                   only: ppp_main
      use refinement,            only: n_updAMR, emergency_fix
      use timer,                 only: set_timer
      use unified_ref_crit_list, only: urc_list
#ifdef GRAV
      use gravity,               only: update_gp
#endif /* GRAV */

      implicit none

      integer, optional, intent(out) :: act_count             !< counts number of blocks refined or deleted
      logical, optional, intent(in)  :: refinement_fixup_only !< When present and .true. then do not check refinement criteria, do only correction, if necessary.

      integer :: nciter, top_level
      integer, parameter :: nciter_max = 100 ! should be more than refinement levels
      logical :: some_refined, derefined, ctop_exists, fu
      type(cg_list_element), pointer :: cgl, aux
      type(cg_level_connected_t), pointer :: curl
      type(grid_container),  pointer :: cg
      logical :: correct, full_update
      real :: ts  !< time for runtime profiling
      character(len=*), parameter :: newref_label = "refine", deref_label = "derefinement", prol_label = "prolong_new"

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

      ! Maybe a bit overkill
      fu = full_update
      call piernik_MPI_Allreduce(fu, pLAND)
      if (fu .neqv. full_update) call die("[refinement_update:update_refinement] inconsistent full_update")

      if (full_update) then
         call scan_for_refinements

         ! do the refinements first
         curl => base%level
         do while (associated(curl))

            cgl => curl%first
            do while (associated(cgl))
               if (cgl%cg%has_leaves()) then
                  if (cgl%cg%flag%pending_blocks()) then
                     call refine_one_grid(curl, cgl)
                     if (present(act_count)) act_count = act_count + 1
                  endif
               endif
               cgl => cgl%nxt
            enddo

            call finest%equalize

            top_level = finest%level%l%id
            call piernik_MPI_Allreduce(top_level, pMAX)
            if (top_level /= finest%level%l%id) call die("[refinement_update:update_refinement] inconsistent top level (r)")

            call ppp_main%start(prol_label, PPP_AMR)
            !> \todo merge small blocks into larger ones
            if (associated(curl%finer)) then
               call curl%finer%init_all_new_cg
               call curl%finer%deallocate_patches
               call curl%finer%sync_ru
               call curl%prolong
            endif
            call ppp_main%stop(prol_label, PPP_AMR)

            curl => curl%finer
         enddo

         call all_cg%mark_orphans

         call leaves%update(" (  refine  ) ")
      endif

      ! fix the structures, mark grids for refinement (and unmark derefinements) due to refinement restrictions.
      ! With clever use of SFC properties this part can be done together with refine stage above - everything can be determined in sanitize_ref_flags for regular refinement update.
      ! For refinement update due to domain expansion, everything can be fixed in the expansion routine.
      correct = fix_refinement()
      call piernik_MPI_Allreduce(correct, pLAND)
      nciter = 0
      do while (.not. correct)
         ! exploit dot and iterate only when dot is not complete
         ! perform quick dot-based checks instead of wa workaround
         ! warn if dot%is_complete and any corrections were required

         call all_cg%prevent_prolong

         curl => finest%level
         do while (associated(curl))
            call curl%deallocate_patches
            curl => curl%coarser
         enddo

         call ppp_main%start(newref_label, PPP_AMR)

         ! Maybe a bit overkill - call finest%equalize should've take care of that
         top_level = finest%level%l%id
         call piernik_MPI_Allreduce(top_level, pMAX)
         if (top_level /= finest%level%l%id) call die("[refinement_update:update_refinement] inconsistent top level (c)")

         curl => finest%level%coarser
         ctop_exists = associated(curl)
         call piernik_MPI_Allreduce(ctop_exists, pLAND)
         if (ctop_exists .neqv. associated(curl)) call die("[refinement_update:update_refinement] inconsistent coarser than top level")

         do while (associated(curl))
            if (curl%l%id >= base%level%l%id) then

               some_refined = .false.
               cgl => curl%first
               do while (associated(cgl))
                  if (cgl%cg%flag%pending_blocks()) then
                     if (finest%level%l%id <= cgl%cg%l%id) call die("[refinement_update:update_refinement] growing too fast!")
                     if (associated(curl%finer)) then
                        call refine_one_grid(curl, cgl)
                        if (present(act_count)) act_count = act_count + 1
                        some_refined = .true.
                     else
                        call die("[refinement_update:update_refinement] nowhere to add!")
                     endif
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
            endif

            curl => curl%coarser
         enddo
         call ppp_main%stop(newref_label, PPP_AMR)

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
      ! Derefinement saves memory and CPU usage, but it is of lowest priority.
      ! Just go through derefinement stage once and don't try to do too much here at once.
      ! Any excess of refinement will be handled in the next call to this routine anyway.
      if (full_update) then

         call ppp_main%start(deref_label, PPP_AMR)
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
         call ppp_main%stop(deref_label, PPP_AMR)

         ! sync structure
         call leaves%balance_and_update(" ( derefine ) ")

      endif

      ! Check refinement topology and crash if anything got broken
      if (.not. fix_refinement()) call die("[refinement_update:update_refinement] Refinement defects still present")

      call all_bnd

      call urc_list%plot_mark(leaves%first)

      !> \todo call the update of cs_i2 and other vital variables if and only if something has changed
      !> \todo add another flag to named_array_list::na_var so the user can also specify fields that need boundary updates on fine/coarse boundaries
      call all_bnd_vital_q

      call all_cg%enable_prolong

#ifdef GRAV
      call update_gp
#endif /* GRAV */

      if (present(act_count)) call piernik_MPI_Allreduce(act_count, pSUM)

      emergency_fix = .false.

      if (full_update) then
         call print_time("[refinement_update] Finishing (full update)")
      else
         call print_time("[refinement_update] Finishing (fixup only)")
      endif

      call log_req

   contains

      !>
      !! \brief Print how much time it took to execute a given stage
      !! Call only for stages that don't call cg_leaves:update.
      !<

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

      !> \brief Notify when the size of the req arrays has increased

      subroutine log_req

         use constants,  only: pMAX
         use dataio_pub, only: printinfo, msg
         use mpisetup,   only: master, req, req2, piernik_MPI_Allreduce

         implicit none

         integer, save :: oldsize = 0
         integer :: cursize

         cursize = 0
         if (allocated(req)) cursize = cursize + size(req)
         if (allocated(req2)) cursize = cursize + size(req2)

         call piernik_MPI_Allreduce(cursize, pMAX)
         if (cursize > oldsize) then
            if (master) then
               write(msg, *) cursize
               call printinfo("[refinement_update] Total number of simultaneous nonblocking MPI requests increased to " // trim(adjustl(msg)))
            endif
            oldsize = cursize
         endif

      end subroutine log_req

   end subroutine update_refinement_wrapped

!> \brief Refine a single grid piece. Pay attention whether it is already refined

   subroutine refine_one_grid(curl, cgl)

      use cg_cost_data,       only: I_REFINE
      use cg_level_connected, only: cg_level_connected_t
      use cg_level_finest,    only: finest
      use cg_list,            only: cg_list_element
      use refinement,         only: bsize
      use dataio_pub,         only: warn, die

      implicit none

      type(cg_level_connected_t), pointer, intent(inout) :: curl
      type(cg_list_element),      pointer, intent(in)    :: cgl

      integer :: b

      if (.not. cgl%cg%has_leaves()) then
         call warn("[refinement_update:refine_one_grid] Attempting to refine a grid that is completely refined")
         return
      endif

      if (.not. associated(curl%finer)) call finest%add_finer

      call cgl%cg%costs%start
      associate (flag => cgl%cg%flag)
         if (flag%pending_blocks()) then
            do b = lbound(flag%SFC_refine_list, dim=1), ubound(flag%SFC_refine_list, dim=1)
               if (flag%SFC_refine_list(b)%level == curl%finer%l%id) then
                  call curl%finer%add_patch(int(bsize, kind=8), flag%SFC_refine_list(b)%off)
               else
                  call die("[refinement_update:refine_one_grid] wrong level!")
               endif
            enddo
            call flag%init ! it is safer to forget it now
            ! flag%init clears flag%derefine as a side-effect but it should be already cleared by sanitization anyway
         else
            call die("[refinement_update:refine_one_grid] populating flag%SFC_refine_list before refining a grid is required")
         endif
      end associate
      call cgl%cg%costs%stop(I_REFINE)

   end subroutine refine_one_grid

!>
!! \brief Apply some rules to fix refinement defects such as more than one level change across a face, edge or corner.
!!
!! \details Important: requires updated cg_leafmap!
!!
!! OPT: use %dot for direct identification if possible without invoking prolongation and restriction
!<

   logical function fix_refinement() result(correct)

      use cg_cost_data,     only: I_REFINE
      use cg_level_finest,  only: finest
      use cg_list,          only: cg_list_element
!      use cg_list_global,   only: all_cg
      use cg_leaves,        only: leaves
      use constants,        only: xdim, ydim, zdim, PPP_AMR, I_ONE!, I_TWO
      use dataio_pub,       only: die
      use domain,           only: dom
      use ppp,              only: ppp_main
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
      character(len=*), parameter :: fix_ref_label = "fix_refinement", frp_label = "fix_refinement_wa"

      call ppp_main%start(fix_ref_label, PPP_AMR)
      correct = .true.
      if (dom%nb < 2) call die("[refinement_update:fix_refinement] dom%nb >= 2 required")

      !> \todo check for excess of refinement levels

      ! Put a level number to the working array, restrict it and exchange internal boundaries
      cgl => leaves%first
      do while (associated(cgl))
         call cgl%cg%costs%start

         cgl%cg%wa = cgl%cg%l%id

         call cgl%cg%costs%stop(I_REFINE)
         cgl => cgl%nxt
      enddo

      call ppp_main%start(frp_label, PPP_AMR)
      call finest%level%restrict_to_base_q_1var(qna%wai)
      call leaves%leaf_arr3d_boundaries(qna%wai)
      call ppp_main%stop(frp_label, PPP_AMR)

      range = 1
      ! range = min((dom%nb+I_ONE)/I_TWO + all_cg%ord_prolong_nb, dom%nb)
      ! (dom%nb+1)/2 + all_cg%ord_prolong_nb is a range of influence of coarse to fine levels - it can be suitable if we were looking for too low levels in the neighbourhood
      ! ATM we're looking for high levels, so range = 1 seems to be appropriate

      ! detect high refinements near leafmap and alter refinement flags appropriately
      cgl => leaves%first
      do while (associated(cgl))

         associate (cg => cgl%cg)
            call cg%costs%start

            call cg%flag%reset_blocks
            call cg%flag%clear

            if (cg%has_leaves()) then
               cg%prolong_xyz = OUTSIDE
               lleaf = -huge(1)

               ! find the level of refinement on the leaf part and mark it with a negative value
               do k = lbound(cg%leafmap, dim=3, kind=4), ubound(cg%leafmap, dim=3, kind=4)
                  do j = lbound(cg%leafmap, dim=2, kind=4), ubound(cg%leafmap, dim=2, kind=4)
                     do i = lbound(cg%leafmap, dim=1, kind=4), ubound(cg%leafmap, dim=1, kind=4)
                        if (cg%leafmap(i, j, k)) then
                           if (lleaf == -huge(1)) then
                              lleaf = int(cg%wa(i, j, k))
                           else
                              if (lleaf /= int(cg%wa(i, j, k)) .or. lleaf /= cg%l%id) call die("[refinement_update:fix_refinement] Inconsistent level map")
                           endif
                           cg%prolong_xyz(i,j,k) = INSIDE
                        endif
                     enddo
                  enddo
               enddo

               ! find the border of the leaf map and mark it with positive value (requires dom%nb >= 2)
               do k = lbound(cg%leafmap, dim=3, kind=4)-dom%D_z, ubound(cg%leafmap, dim=3, kind=4)+dom%D_z
                  do j = lbound(cg%leafmap, dim=2, kind=4)-dom%D_y, ubound(cg%leafmap, dim=2, kind=4)+dom%D_y
                     do i = lbound(cg%leafmap, dim=1, kind=4)-dom%D_x, ubound(cg%leafmap, dim=1, kind=4)+dom%D_x
                        if (int(cg%prolong_xyz(i, j, k)) == OUTSIDE) then
                           if (dom%has_dir(xdim)) then
                              if (int(cg%prolong_xyz(i+1, j, k)) /= 0) cg%prolong_xyz(i+1, j, k) = BOUNDARY
                              if (int(cg%prolong_xyz(i-1, j, k)) /= 0) cg%prolong_xyz(i-1, j, k) = BOUNDARY
                           endif
                           if (dom%has_dir(ydim)) then
                              if (int(cg%prolong_xyz(i, j+1, k)) /= 0) cg%prolong_xyz(i, j+1, k) = BOUNDARY
                              if (int(cg%prolong_xyz(i, j-1, k)) /= 0) cg%prolong_xyz(i, j-1, k) = BOUNDARY
                           endif
                           if (dom%has_dir(zdim)) then
                              if (int(cg%prolong_xyz(i, j, k+1)) /= 0) cg%prolong_xyz(i, j, k+1) = BOUNDARY
                              if (int(cg%prolong_xyz(i, j, k-1)) /= 0) cg%prolong_xyz(i, j, k-1) = BOUNDARY
                           endif
                        endif
                     enddo
                  enddo
               enddo

               !find highest refinement level within a "range" cells from the border of the leaf map
! dirty check
!!$            if (maxval(cg%wa) > 10) then
!!$               do k = lbound(cg%wa, dim=3), ubound(cg%wa, dim=3)
!!$                  do j = lbound(cg%wa, dim=2), ubound(cg%wa, dim=2)
!!$                     do i = lbound(cg%wa, dim=1), ubound(cg%wa, dim=1)
!!$                        if (cg%wa(i,j,k) > 10) write(msg,*)"ur:fr WA(",i,j,k,") = ",cg%wa(i,j,k)
!!$                     enddo
!!$                  enddo
!!$               enddo
!!$            endif

               do k = lbound(cg%leafmap, dim=3, kind=4), ubound(cg%leafmap, dim=3, kind=4)
                  do j = lbound(cg%leafmap, dim=2, kind=4), ubound(cg%leafmap, dim=2, kind=4)
                     do i = lbound(cg%leafmap, dim=1, kind=4), ubound(cg%leafmap, dim=1, kind=4)
                        if (int(cg%prolong_xyz(i, j, k)) == BOUNDARY) then
                           lnear = int(min(huge(I_ONE)/10.,maxval(cg%wa(i-range*dom%D_x:i+range*dom%D_x, j-range*dom%D_y:j+range*dom%D_y, k-range*dom%D_z:k+range*dom%D_z))))
                           if (lnear > lleaf) then
                              cg%flag%derefine = .false.
                              if (lnear > lleaf+1 .and. lnear <= finest%level%l%id) then
                                 call cg%flag%set(i, j, k)
                                 correct = .false.
                              endif
                           endif

                        endif
                     enddo
                  enddo
               enddo
               call cg%refinemap2SFC_list
            endif

            call cg%costs%stop(I_REFINE)
         end associate

         cgl => cgl%nxt
      enddo
      call ppp_main%stop(fix_ref_label, PPP_AMR)

   end function fix_refinement

end module refinement_update
