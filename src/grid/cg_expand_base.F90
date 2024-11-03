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
!! \brief This module implements expansion of the base level
!!
!! \todo Domain shrinking, crawling, uplifting etc. should also be implemented here, when needed
!!      procedure :: shrink      !< delete one line of blocks in some directions
!!      procedure, private :: shrink_side !< delete one line of blocks in selected direction
!!      procedure :: lift        !< move one level higher
!!      procedure :: plunge      !< move one level lower
!!      procedure :: condense    !< shrink on all sides (halve size) and plunge
!!      procedure :: inflate     !< expand on all sides to double size and lift
!<

module cg_expand_base

   implicit none

   private
   public :: expand_base

contains

!> \brief Add one line of blocks in some directions (wrapper for expand_side)

   subroutine expand_base(sides)

      use all_boundaries,    only: all_fluid_boundaries
      use allreduce,         only: piernik_MPI_Allreduce
      use cg_level_coarsest, only: coarsest
      use cg_level_finest,   only: finest
      use constants,         only: xdim, zdim, LO, HI, pLOR, V_ESSENTIAL
      use dataio_pub,        only: msg, printinfo
      use domain,            only: dom
      use mpisetup,          only: master
#ifdef MULTIGRID
      use multigrid,         only: init_multigrid
#endif /* MULTIGRID */

      implicit none

      logical, dimension(xdim:zdim, LO:HI), intent(in)    :: sides  !< logical mask of sides to be extended

      integer :: d, lh
      logical :: changed, ch, side

      changed = .false.
      do d = xdim, zdim
         if (dom%has_dir(d)) then
            do lh = LO, HI
               side = sides(d, lh)
               call piernik_MPI_Allreduce(side, pLOR)
               if (side) then
                  ch = expand_side(d, lh)
                  changed = changed .or. ch
               endif
            enddo
         endif
      enddo

      if (changed) then
         if (master) then
            write(msg, '(a,3i8,a,i3)')"[cg_expand_base:expand_base] Effective resolution is [", finest%level%l%n_d(:), " ] at level ", finest%level%l%id
            call printinfo(msg, V_ESSENTIAL) ! As long as the restart file does not automagically recognize changed parameters, this message should be easily visible
         endif
         call coarsest%delete_coarser_than_base
#ifdef MULTIGRID
         call init_multigrid
#endif /* MULTIGRID */
         call all_fluid_boundaries
         ! the cg%gp and cg%cs_iso2 are updated in refinement_update::update_refinement which should be called right after domain expansion to fix refinement structure
      endif

   end subroutine expand_base

!> \brief Add one line of blocks in selected direction

   logical function expand_side(d, lh)

      use cg_leaves,          only: leaves
      use cg_level_base,      only: base
      use cg_level_connected, only: cg_level_connected_t
      use cg_list,            only: cg_list_element
      use cg_list_dataop,     only: expanded_domain
      use constants,          only: xdim, zdim, LO, HI, refinement_factor, &
           &                        BND_OUT, BND_OUTD, BND_OUTH, BND_OUTHD, BND_PER, BND_REF, BND_MPI, BND_FC
      use dataio_pub,         only: die, warn, msg
      use domain,             only: dom
      use mpisetup,           only: master
      use refinement,         only: emergency_fix, bsize
      use user_hooks,         only: late_initial_conditions

      implicit none

      integer,                intent(in)    :: d      !< direction to be expanded
      integer,                intent(in)    :: lh     !< side to be expanded

      integer(kind=8), dimension(xdim:zdim) :: e_size, e_off, new_n_d, new_off
      type(cg_list_element),  pointer :: cgl
      type(cg_level_connected_t), pointer :: curl

      expand_side = .false.

      if (.not. dom%has_dir(d)) call die("[cg_expand_base:expand_side] Non-existing direction")
      if (bsize(d) < dom%nb) call die("[cg_expand_base:expand_side] Invalid AMR::bsize")
      if (.not. associated(late_initial_conditions)) call die("[cg_expand_base:expand_side] You must provide a routine to initialize vital arrays for current time")

      if (any(dom%bnd(d, lh) == [ BND_PER, BND_REF ])) return
      if (all(dom%bnd(d, lh) /= [ BND_OUT, BND_OUTD, BND_OUTH, BND_OUTHD ])) then
         write(msg, '(a,i2,a,2i2,a)')"[cg_expand_base:expand_side] Not sure if boundary type ", dom%bnd(d, lh), " at direction (", d, lh, ") is expandable. You were warned."
         call warn(msg)
      endif

      expand_side = .true.

      call base%level%set_is_old
      cgl => leaves%first
      do while (associated(cgl))
         if (cgl%cg%ext_bnd(d, lh)) then
            cgl%cg%ext_bnd(d, lh) = .false.
            if (cgl%cg%l%id == base%level%l%id) then
               cgl%cg%bnd(d, lh) = BND_MPI
            else
               cgl%cg%bnd(d, lh) = BND_FC
            endif
         endif
         cgl => cgl%nxt
      enddo

      e_size = base%level%l%n_d
      e_size(d) = bsize(d)

      e_off = base%level%l%off
      select case (lh)
         case (LO)
            e_off(d) = base%level%l%off(d) - bsize(d)
         case (HI)
            e_off(d) = base%level%l%off(d) + base%level%l%n_d(d)
      end select

      curl => base%level
      do while (associated(curl))
         new_n_d = curl%l%n_d
         new_n_d(d) = curl%l%n_d(d) + bsize(d)*refinement_factor**(curl%l%id-base%level%l%id)
         new_off = curl%l%off
         new_off(d) = min(curl%l%off(d), e_off(d)*refinement_factor**(curl%l%id-base%level%l%id))
         call curl%l%update(curl%l%id, new_n_d, new_off)
         call curl%refresh_SFC_id
         curl => curl%finer
      enddo
      ! multigrid levels are destroyed and re-created in expand_base

      call dom%modify_side(d, lh, bsize(d))
      if (master) call base%level%add_patch(e_size, e_off)
      call base%level%init_all_new_cg

      cgl => base%level%first
      do while (associated(cgl))
         if (.not. cgl%cg%is_old) call expanded_domain%add(cgl%cg)
         cgl => cgl%nxt
      enddo

      call base%level%sync_ru
      call leaves%update(" (  expand  ) ") !cannot call balance here as it will mess up the expanded_domain list
      call base%level%deallocate_patches
      call late_initial_conditions
      emergency_fix = .true.

   end function expand_side

end module cg_expand_base
