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
#include "macros.h"

!> \brief This module contains routines for managing refinement and derefinement on the domain

module refinement_update

   implicit none

   private
   public :: update_refinement

contains

   subroutine scan_for_refinements

      use cg_level_connected,    only: cg_level_connected_T
      use cg_level_finest,       only: finest
      use cg_list_global,        only: all_cg
      use refinement_primitives, only: mark_all_primitives
      use user_hooks,            only: problem_refine_derefine
#ifdef VERBOSE
      use constants,             only: pSUM
      use dataio_pub,            only: msg, printinfo
      use mpisetup,              only: master, piernik_MPI_Allreduce
#endif

      implicit none

      type(cg_level_connected_T), pointer :: curl
      enum, bind(C)
         enumerator :: PROBLEM
         !enumerator :: SHOCKS
         enumerator :: PRIMITIVES
      end enum
      integer, dimension(PROBLEM:PRIMITIVES) :: cnt

      curl => finest%level
      do while (associated(curl))
         if (allocated(curl%patches)) deallocate(curl%patches)
         curl => curl%coarser
      enddo
      call all_cg%clear_ref_flags
      cnt = 0

      if (associated(problem_refine_derefine)) call problem_refine_derefine ! call user routine first, so it cannot alter flags set by automatic routines

#ifdef VERBOSE
      cnt(PROBLEM) = all_cg%count_ref_flags()
      call piernik_MPI_Allreduce(cnt(PROBLEM), pSUM)
#endif

      ! call mark_shocks !> \todo implement automatic refinement criteria
#ifdef VERBOSE
!      cnt(SHOCKS) = all_cg%count_ref_flags()
!      call piernik_MPI_Allreduce(cnt(SHOCKS), pSUM)
#endif

      call mark_all_primitives
#ifdef VERBOSE
      cnt(PRIMITIVES) = all_cg%count_ref_flags()
      call piernik_MPI_Allreduce(cnt(PRIMITIVES), pSUM)
      if (cnt(ubound(cnt, dim=1)) > 0) then
         write(msg,'(2(a,i6),a)')"[refinement_update:scan_for_refinements] User-defined routine marked ", cnt(PROBLEM), " block(s) for refinement, primitives marked ",cnt(PRIMITIVES)," block(s)"
      else
         write(msg,'(a)')"[refinement_update:scan_for_refinements] No blocks marked for refinement"
      endif
      if (master) call printinfo(msg)
#endif

   end subroutine scan_for_refinements

!>
!! \brief Update the refinement topology
!!
!! \deprecated Now it works only on whole domains
!!
!! \deprecated lot of duplicated code with grid::dom2cg and grid::init_grid
!<

!#define DEBUG_DUMPS

   subroutine update_refinement(act_count, refinement_fixup_only)

      use all_boundaries,     only: all_bnd
      use cg_leaves,          only: leaves
      use cg_list,            only: cg_list_element
      use cg_level_base,      only: base
      use cg_level_connected, only: cg_level_connected_T
      use cg_level_finest,    only: finest
      use cg_list_global,     only: all_cg
      use constants,          only: pLOR, pLAND, pSUM, cs_i2_n
      use dataio_pub,         only: warn, die
      use global,             only: nstep
      use grid_cont,          only: grid_container
      use list_of_cg_lists,   only: all_lists
      use mpisetup,           only: piernik_MPI_Allreduce!, proc
      use named_array_list,   only: qna
      use refinement,         only: n_updAMR, emergency_fix
#ifdef GRAV
      use gravity,            only: update_gp
#endif /* GRAV */
#ifdef DEBUG_DUMPS
      use data_hdf5,          only: write_hdf5
#endif /* DEBUG_DUMPS */

      implicit none

      integer, optional, intent(out) :: act_count             !< counts number of blocks refined or deleted
      logical, optional, intent(in)  :: refinement_fixup_only !< When present and .true. then do not check refinement criteria, do only correction, if necessary.

      integer :: nciter
      integer, parameter :: nciter_max = 100 ! should be more than refinement levels
      logical :: some_refined, derefined
      type(cg_list_element), pointer :: cgl, aux
      type(cg_level_connected_T), pointer :: curl
      type(grid_container),  pointer :: cg
      logical :: correct, full_update

      if (present(act_count)) act_count = 0

      full_update = .true.
      if (present(refinement_fixup_only)) full_update = .not. refinement_fixup_only

      call finest%level%restrict_to_base ! implies update of leafmap

      if (n_updAMR <= 0) return
      if (mod(nstep, n_updAMR) /= 0 .and. .not. emergency_fix) return

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
                  call cgl%cg%refine_flags%sanitize(cgl%cg%level_id)
                  if (cgl%cg%refine_flags%refine) then
                     call refine_one_grid(curl, cgl)
                     if (present(act_count)) act_count = act_count + 1
                     cgl%cg%refine_flags%refine = .false.
                  endif
               endif
               cgl => cgl%nxt
            enddo

            call finest%equalize !> \todo Try to split this loop here and place this call between the loops

            !> \todo merge small blocks into larger ones
            if (associated(curl%finer)) then
               call curl%finer%init_all_new_cg
               if (allocated(curl%finer%patches)) deallocate(curl%finer%patches)
               call curl%finer%sync_ru
               call curl%prolong
            endif

            curl => curl%finer
         enddo

         call all_cg%mark_orphans

         call leaves%update(" (  refine  ) ")
      endif

      ! fix the structures, mark grids for refinement (and unmark derefinements) due to refinement restrictions
      call all_cg%clear_ref_flags
      call fix_refinement(correct)
      call piernik_MPI_Allreduce(correct, pLAND)
      nciter = 0
      do while (.not. correct)

         call all_cg%prevent_prolong

         curl => finest%level
         do while (associated(curl))
            if (allocated(curl%patches)) deallocate(curl%patches)
            curl => curl%coarser
         enddo

         curl => finest%level%coarser
         do while (associated(curl) .and. curl%level_id >= base%level%level_id)
            some_refined = .false.
            cgl => curl%first
            do while (associated(cgl))
               call cgl%cg%refine_flags%sanitize(cgl%cg%level_id)
               if (cgl%cg%refine_flags%refine) then
                  if (finest%level%level_id <= cgl%cg%level_id) call warn("[refinement_update:update_refinement] growing too fast!")
!                  write(msg,*)"addp ^",curl%level_id," ^^",curl%level_id+1," @[]",cgl%cg%my_se(:, LO)*refinement_factor, " []",cgl%cg%n_b(:)*refinement_factor
                  if (associated(curl%finer)) then
                     call refine_one_grid(curl, cgl)
                     if (present(act_count)) act_count = act_count + 1
                     some_refined = .true.
                  else
                     call warn("[refinement_update:update_refinement] nowhere to add!")
                  endif
                  cgl%cg%refine_flags%refine = .false.
               endif
               cgl => cgl%nxt
            enddo

            call piernik_MPI_Allreduce(some_refined, pLOR)
            if (some_refined) then
               call curl%finer%init_all_new_cg
               call curl%finer%sync_ru
               if (allocated(curl%patches)) deallocate(curl%patches)
               !call finest%equalize
               call curl%prolong
               call all_cg%mark_orphans
            endif

            curl => curl%coarser
         enddo

         ! sync structure before trying to fix it
         call leaves%update(" (correcting) ")
         call all_cg%clear_ref_flags
         call fix_refinement(correct)
         call piernik_MPI_Allreduce(correct, pLAND)
         if (.not. correct) nciter = nciter + 1
         if (nciter > nciter_max) then
            call die("[refinement_update:update_refinement] Cannot converge to correct structure. Too many iterations.")
            correct = .true.
         endif
      enddo

      ! Now try to derefine any excess of refinement
      if (full_update) call scan_for_refinements
      call fix_refinement(correct)
      if (.not. correct) call die("[refinement_update:update_refinement] Refinement defects still present")

      ! Derefinement saves memory and CPU usage, but is not of highest priority.
      ! Just do it once and hope that any massive excess of refinement will be handled in next call to this routine
      if (full_update) then

         curl => finest%level
         do while (associated(curl) .and. .not. associated(curl, base%level))
            cgl => curl%first
            derefined = .false.
            do while (associated(cgl))
               call cgl%cg%refine_flags%sanitize(cgl%cg%level_id)
               aux => cgl ! Auxiliary pointer makes it easier to loop over the list of grids when some of the elements are disappearing
               cgl => cgl%nxt
               if (aux%cg%refine_flags%derefine) then
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
         call fix_refinement

      endif

      call all_bnd
      !> \todo call the update of cs_i2 if and only if something has changed
      !> \todo add another flag to named_array_list::na_var so the user can also specify fields that need boundary updates on fine/coarse boundaries
      if (qna%exists(cs_i2_n)) call leaves%leaf_arr3d_boundaries(qna%ind(cs_i2_n))
#ifdef GRAV
      call update_gp
#endif /* GRAV */

      call all_cg%enable_prolong
      if (present(act_count)) call piernik_MPI_Allreduce(act_count, pSUM)

#ifdef DEBUG_DUMPS
      call write_hdf5
#endif /* DEBUG_DUMPS */

      emergency_fix = .false.

   end subroutine update_refinement

!> \brief Refine a single grid piece. Pay attention wheter it is already refined

   subroutine refine_one_grid(curl, cgl)

      use cg_level_connected, only: cg_level_connected_T
      use cg_level_finest,    only: finest
      use cg_list,            only: cg_list_element
      use constants,          only: refinement_factor, LO, HI, ndims
      use dataio_pub,         only: warn
      use mergebox,           only: wmap

      implicit none

      type(cg_level_connected_T), pointer, intent(inout) :: curl
      type(cg_list_element),      pointer, intent(in)    :: cgl

      type(wmap) :: lmap
      integer(kind=8), dimension(ndims, LO:HI)  :: box_8   !< temporary storage
      integer :: b

      if (.not. any(cgl%cg%leafmap)) then
         call warn("[refinement_update:refine_one_grid] Attempting to refine a grid that is completely refined")
         return
      endif

      if (.not. associated(curl%finer)) call finest%add_finer

      if (.not. all(cgl%cg%leafmap)) then ! decompose the partially refined grid container into boxes that contain all leafcells
         box_8 = int(cgl%cg%ijkse, kind=8)
         call lmap%init(box_8)
         lmap%map(cgl%cg%is:cgl%cg%ie, cgl%cg%js:cgl%cg%je, cgl%cg%ks:cgl%cg%ke) = cgl%cg%leafmap(:,:,:)
         call lmap%find_boxes
         do b = lbound(lmap%blist%blist, dim=1), ubound(lmap%blist%blist, dim=1)
            call curl%finer%add_patch(int((lmap%blist%blist(b)%b(:,HI)-lmap%blist%blist(b)%b(:,LO)+1)*refinement_factor, kind=8), lmap%blist%blist(b)%b(:,LO)*refinement_factor)
         enddo
         call lmap%cleanup
      else
         call curl%finer%add_patch(int(cgl%cg%n_b(:)*refinement_factor, kind=8), cgl%cg%my_se(:, LO)*refinement_factor)
      endif

   end subroutine refine_one_grid

!>
!! \brief Apply some rules to fix refinement defects
!!
!! \details Important! Requires updated cg_leafmap!
!<

   subroutine fix_refinement(correct)

      use cg_level_finest,  only: finest
      use cg_list,          only: cg_list_element
!      use cg_list_global,   only: all_cg
      use cg_leaves,        only: leaves
      use constants,        only: xdim, ydim, zdim, I_ONE!, I_TWO, LO, HI
      use dataio_pub,       only: die, warn, msg!, printinfo
      use domain,           only: dom
      use refinement,       only: allow_face_rstep, allow_corner_rstep
      use named_array_list, only: qna
!      use cg_level_base, only: base
#ifdef DEBUG_DUMPS
      use constants,        only: pLOR
      use mpisetup,         only: piernik_MPI_Allreduce
      use data_hdf5,        only: write_hdf5
#endif /* DEBUG_DUMPS */

      implicit none

      logical, optional, intent(out) :: correct !< Flag to tell that corrections are required.

      type(cg_list_element), pointer :: cgl
      integer ::  i, j, k
      integer :: lleaf, lnear, range
      logical :: failed
      enum, bind(C)
         enumerator :: INSIDE   = -1
         enumerator :: OUTSIDE  =  0
         enumerator :: BOUNDARY =  1
      end enum

      if (present(correct)) correct = .true.
      failed = .false.

      if (allow_face_rstep .and. allow_corner_rstep) return
      !> \todo also check for excess or refinement levels

      ! Put a level number to the working array, restrict it and exchange internal boundaries
      cgl => leaves%first
      do while (associated(cgl))
         cgl%cg%wa = cgl%cg%level_id
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
         if (any(cgl%cg%leafmap)) then
            cgl%cg%prolong_xyz = OUTSIDE
            lleaf = -huge(1)

            ! find the level of refinement on the leaf part and mark it with a negative value
            do k = lbound(cgl%cg%leafmap, dim=3), ubound(cgl%cg%leafmap, dim=3)
               do j = lbound(cgl%cg%leafmap, dim=2), ubound(cgl%cg%leafmap, dim=2)
                  do i = lbound(cgl%cg%leafmap, dim=1), ubound(cgl%cg%leafmap, dim=1)
                     if (cgl%cg%leafmap(i, j, k)) then
                        if (lleaf == -huge(1)) then
                           lleaf = int(cgl%cg%wa(i, j, k))
                        else
                           if (lleaf /= int(cgl%cg%wa(i, j, k)) .or. lleaf /= cgl%cg%level_id) call die("[refinement_update:fix_refinement] Inconsistent level map")
                        endif
                        cgl%cg%prolong_xyz(i,j,k) = INSIDE
                     endif
                  enddo
               enddo
            enddo

            ! find the border of the leaf map and mark it with positive value
            do k = lbound(cgl%cg%leafmap, dim=3)-dom%D_z, ubound(cgl%cg%leafmap, dim=3)+dom%D_z
               do j = lbound(cgl%cg%leafmap, dim=2)-dom%D_y, ubound(cgl%cg%leafmap, dim=2)+dom%D_y
                  do i = lbound(cgl%cg%leafmap, dim=1)-dom%D_x, ubound(cgl%cg%leafmap, dim=1)+dom%D_x
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
            lnear = -huge(1)
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

            do k = lbound(cgl%cg%leafmap, dim=3), ubound(cgl%cg%leafmap, dim=3)
               do j = lbound(cgl%cg%leafmap, dim=2), ubound(cgl%cg%leafmap, dim=2)
                  do i = lbound(cgl%cg%leafmap, dim=1), ubound(cgl%cg%leafmap, dim=1)
                     if (int(cgl%cg%prolong_xyz(i, j, k)) == BOUNDARY) &
                          lnear = max(lnear, int(min(huge(I_ONE)/10.,maxval(cgl%cg%wa(i-range*dom%D_x:i+range*dom%D_x, j-range*dom%D_y:j+range*dom%D_y, k-range*dom%D_z:k+range*dom%D_z)))))
                  enddo
               enddo
            enddo

            if (lnear > lleaf) then
               cgl%cg%refine_flags%derefine = .false.
               if (lnear > lleaf+1 .and. lnear <= finest%level%level_id) then
                  cgl%cg%refine_flags%refine = .true.
                  if (present(correct)) then
                     correct = .false.
                  else
                     write(msg,'(a,i3,a,6i5,a,i3)')"[refinement_update:fix_refinement] neighbour level ^",lnear," at [",cgl%cg%my_se,"] ^",cgl%cg%level_id
                     call warn(msg)
                     failed = .true.
                  endif
               endif
            endif

         else
            cgl%cg%refine_flags%refine = .false.
            cgl%cg%refine_flags%derefine = .false.
         endif
         cgl => cgl%nxt
      enddo

#ifdef DEBUG_DUMPS
      call piernik_MPI_Allreduce(failed, pLOR)
#endif
      if (failed) then
#ifdef DEBUG_DUMPS
         call write_hdf5
#endif /* DEBUG_DUMPS */
         call die("[refinement_update:fix_refinement] Refinement defects found.")
      endif
!!$      call leaves%corners2wa(qna%wai)

   end subroutine fix_refinement

end module refinement_update

!!$!> \brief Add a new level - refine whole domain
!!$
!!$   subroutine refine_domain
!!$
!!$#if defined(__INTEL_COMPILER)
!!$   !! \deprecated remove this clause as soon as Intel Compiler gets required
!!$   !! features and/or bug fixes
!!$      use cg_level_connected, only: cg_level_connected_T  ! QA_WARN INTEL
!!$#endif /* __INTEL_COMPILER */
!!$      use cg_level_finest, only: finest
!!$      use dataio_pub,      only: msg, printinfo
!!$      use mpisetup,        only: master
!!$
!!$      implicit none
!!$
!!$      if (master) then
!!$         write(msg, '(a,i3)')"[refinement_update:refine_domain] refining level ",finest%level%level_id
!!$         call printinfo(msg)
!!$      endif
!!$
!!$      !> \todo Check if finest is a complete level first
!!$
!!$      call finest%add_finer
!!$      call finest%level%add_patch
!!$      call finest%level%init_all_new_cg
!!$      call finest%level%coarser%prolong
!!$
!!$   end subroutine refine_domain
!!$
!!$!> \brief Mark finest level for derefinement
!!$
!!$   subroutine derefine_domain
!!$
!!$      use cg_level_finest, only: finest
!!$      use cg_list,         only: cg_list_element
!!$      use dataio_pub,      only: msg, printinfo
!!$      use mpisetup,        only: master
!!$
!!$      implicit none
!!$
!!$      type(cg_list_element), pointer :: cgl
!!$
!!$      if (master) then
!!$         write(msg, '(a,i3)')"[refinement_update:derefine_domain] derefining level ",finest%level%level_id
!!$         call printinfo(msg)
!!$      endif
!!$      call finest%level%restrict
!!$
!!$      cgl => finest%level%first
!!$      do while (associated(cgl))
!!$         cgl%cg%refine_flags%refine   = .false.
!!$         cgl%cg%refine_flags%derefine = .true.
!!$         cgl => cgl%nxt
!!$      enddo
!!$
!!$   end subroutine derefine_domain
