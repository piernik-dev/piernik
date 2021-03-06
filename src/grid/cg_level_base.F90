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

!> \brief This module implements the base level and related methods

module cg_level_base

   use cg_level_connected, only: cg_level_connected_t

   implicit none

   private
   public :: base

   abstract interface
      subroutine no_args
         implicit none
      end subroutine no_args
   end interface

   !! \brief The pointer of the base level and a method to initialize it
   !> \todo Domainshrinking, expanding and crawling should also be implemented here
   type :: cg_level_base_t
      type(cg_level_connected_t), pointer :: level            !< The base level
      procedure(no_args), nopass, pointer :: init_multigrid   !< a pointer to multigrid:init_multigrid or null()
   contains
      procedure          :: set                               !< initialize the base level
      procedure          :: expand                            !< add one line of blocks in some directions
      procedure, private :: expand_side                       !< add one line of blocks in selected direction
!      procedure :: shrink      !< delete one line of blocks in some directions
!      procedure, private :: shrink_side !< delete one line of blocks in selected direction
!      procedure :: lift        !< move one level higher
!      procedure :: plunge      !< move one level lower
!      procedure :: condense    !< shrink on all sides (halve size) and plunge
!      procedure :: inflate     !< expand on all sides to double size and lift
   end type cg_level_base_t

   type(cg_level_base_t), pointer :: base                     !< base level grid containers

contains

!> \brief Initialize the base level

   subroutine set(this, n_d)

      use constants,          only: base_level_id, ndims
      use dataio_pub,         only: die
      use domain,             only: dom
      use list_of_cg_lists,   only: all_lists

      implicit none

      class(cg_level_base_t),            intent(inout) :: this   !< object invoking type bound procedure
      integer(kind=4), dimension(ndims), intent(in)    :: n_d    !< size of global base grid in cells

      ! Multigrid and refinement work properly with non-0, even offset.
      ! Offset value equal to k*2**n, where k is odd will allow at most n levels of coarsening.
      ! Odd offsets or domain sizes prevent creation of coarse levels.
      !> \todo Find the limit that comes from multigrid: maximum refinement should not depend on base level offset
      ! Offset of the base domain may change after the domain gets expanded, shrinked or resized.

      if (any(n_d(:) < 1)) call die("[cg_level_base:set] non-positive base grid sizes")
      if (any(dom%has_dir(:) .neqv. (n_d(:) > 1))) call die("[cg_level_base:set] base grid size incompatible with has_dir masks")

      this%init_multigrid => null()  ! it will be set by init_multigrid, if it is included and called
      allocate(this%level)
      call this%level%init_level

      call this%level%l%init(base_level_id, int(n_d, kind=8), dom%off)

      call all_lists%register(this%level, "Base level")

   end subroutine set

!> \brief Add one line of blocks in some directions (wrapper for expand_side)

   subroutine expand(this, sides)

      use all_boundaries,    only: all_fluid_boundaries
      use cg_level_coarsest, only: coarsest
      use cg_level_finest,   only: finest
      use constants,         only: xdim, zdim, LO, HI, pLOR
      use dataio_pub,        only: msg, warn
      use domain,            only: dom
      use mpisetup,          only: piernik_MPI_Allreduce, master

      implicit none

      class(cg_level_base_t),               intent(inout) :: this   !< object invoking type bound procedure
      logical, dimension(xdim:zdim, LO:HI), intent(in)    :: sides  !< logical mask of sides to be extended

      integer :: d, lh
      logical :: changed, side

      changed = .false.
      do d = xdim, zdim
         if (dom%has_dir(d)) then
            do lh = LO, HI
               side = sides(d, lh)
               call piernik_MPI_Allreduce(side, pLOR)
               if (side) then
                  call this%expand_side(d, lh)
                  changed = .true.
               endif
            enddo
         endif
      enddo

      if (changed) then
         if (master) then
            write(msg, '(a,3i8,a,i3)')"[cg_level_base:expand] Effective resolution is [", finest%level%l%n_d(:), " ] at level ", finest%level%l%id
            call warn(msg) ! As long as the restart file does not automagically recognize changed parameters, this message should be easily visible
         endif
         call coarsest%delete_coarser_than_base
         if (associated(this%init_multigrid)) call this%init_multigrid
         call all_fluid_boundaries
         ! the cg%gp and cg%cs_iso2 are updated in refinement_update::update_refinement which should be called right after domain expansion to fix refinement structure
      endif

   end subroutine expand

!> \brief Add one line of blocks in selected direction

   subroutine expand_side(this, d, lh)

      use cg_leaves,          only: leaves
      use cg_level_connected, only: cg_level_connected_t
      use cg_list,            only: cg_list_element
      use cg_list_dataop,     only: expanded_domain
      use constants,          only: xdim, zdim, LO, HI, BND_MPI, BND_FC, refinement_factor
      use dataio_pub,         only: die
      use domain,             only: dom
      use list_of_cg_lists,   only: all_lists
      use mpisetup,           only: master
      use refinement,         only: emergency_fix, bsize
      use user_hooks,         only: late_initial_conditions

      implicit none

      class(cg_level_base_t), intent(inout) :: this   !< object invoking type bound procedure
      integer,                intent(in)    :: d      !< direction to be expanded
      integer,                intent(in)    :: lh     !< side to be expanded

      integer(kind=8), dimension(xdim:zdim) :: e_size, e_off, new_n_d, new_off
      type(cg_list_element),  pointer :: cgl
      type(cg_level_connected_t), pointer :: curl

      if (.not. dom%has_dir(d)) call die("[cg_level_base:expand_side] Non-existing direction")
      if (bsize(d) < dom%nb) call die("[cg_level_base:expand_side] Invalid AMR::bsize")
      if (.not. associated(late_initial_conditions)) call die("[cg_level_base:expand_side] You must provide a routine to initialize vital arrays for current time")

      call this%level%set_is_old
      cgl => leaves%first
      do while (associated(cgl))
         if (cgl%cg%ext_bnd(d, lh)) then
            cgl%cg%ext_bnd(d, lh) = .false.
            if (cgl%cg%l%id == this%level%l%id) then
               cgl%cg%bnd(d, lh) = BND_MPI
            else
               cgl%cg%bnd(d, lh) = BND_FC
            endif
         endif
         cgl => cgl%nxt
      enddo

      e_size = this%level%l%n_d
      e_size(d) = bsize(d)

      e_off = this%level%l%off
      select case (lh)
         case (LO)
            e_off(d) = this%level%l%off(d) - bsize(d)
         case (HI)
            e_off(d) = this%level%l%off(d) + this%level%l%n_d(d)
      end select

      curl => this%level
      do while (associated(curl))
         new_n_d = curl%l%n_d
         new_n_d(d) = curl%l%n_d(d) + bsize(d)*refinement_factor**(curl%l%id-this%level%l%id)
         new_off = curl%l%off
         new_off(d) = min(curl%l%off(d), e_off(d)*refinement_factor**(curl%l%id-this%level%l%id))
         call curl%l%update(curl%l%id, new_n_d, new_off)
         call curl%refresh_SFC_id
         curl => curl%finer
      enddo
      ! multigrid levels are destroyed and re-created in this%expand

      call dom%modify_side(d, lh, bsize(d))
      if (master) call this%level%add_patch(e_size, e_off)
      call this%level%init_all_new_cg

      call all_lists%register(expanded_domain, "e-dom")
      cgl => this%level%first
      do while (associated(cgl))
         if (.not. cgl%cg%is_old) call expanded_domain%add(cgl%cg)
         cgl => cgl%nxt
      enddo

      call this%level%sync_ru
      call leaves%update(" (  expand  ) ") !cannot call balance here as it will mess up the expanded_domain list
      call this%level%deallocate_patches
      call late_initial_conditions
      emergency_fix = .true.

   end subroutine expand_side

end module cg_level_base
