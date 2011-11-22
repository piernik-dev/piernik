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

!>
!! \brief (DW) Module containing lists of grid containers for computational mesh and initialization and cleanup routines
!<

module grid

   use gc_list,   only: cg_list_global, cg_list_level, cg_list_patch, cg_list

   implicit none

   private
   public :: init_grid, cleanup_grid, all_cg, base_lev, leaves

   type(cg_list_global), protected :: all_cg                          !< all grid containers
   type(cg_list_level), protected  :: base_lev                        !< base level grid containers
   type(cg_list), protected  :: leaves                                !< grid containers not fully covered by finer grid containers
   integer, parameter :: NBD = 1                                      !< at the moment the base domain may be composed of only one patch
   type(cg_list_patch), dimension(NBD), target, protected :: base_dom !< base level patches; \todo relax the NBD=1 restriction if we want something like L-shaped or more complex domains

contains

!>
!! \brief Routine that allocates all grid containers and most important field arrays inside gc's
!<
   subroutine init_grid

      use constants,   only: PIERNIK_INIT_DOMAIN, AT_NO_B, AT_OUT_B, AT_IGNORE, INVALID, LONG, &
           &                 ndims, xdim, zdim, fluid_n, uh_n, mag_n, wa_n, u0_n, b0_n,cs_i2_n
      use dataio_pub,  only: printinfo, die, code_progress
      use domain,      only: pdom, is_multicg
      use fluidindex,  only: flind
      use gc_list,     only: cg_list_element
      use global,      only: repeat_step
      use grid_cont,   only: grid_container
      use mpisetup,    only: proc, inflate_req

      implicit none

      integer :: nrq, d
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
      type(cg_list_patch), pointer :: pbd

      if (code_progress < PIERNIK_INIT_DOMAIN) call die("[grid:init_grid] domain not initialized.")

#ifdef VERBOSE
      call printinfo("[grid:init_grid]: commencing...")
#endif /* VERBOSE */

      ! Create the empty main lists for base level only.
      ! Refinement lists will be added by iterating the initproblem::init_prob routine, in restart_hdf5::read_restart_hdf5 or in not_yet_implemented::refinement_update
      ! Underground levels will be added in multigrid::init_multigrid
      call all_cg%init
      call base_lev%init
      call leaves%init
      do d = lbound(base_dom, dim=1), ubound(base_dom, dim=1) ! currently we have only one base patch
         call base_dom(d)%init
      enddo

      pbd => base_dom(NBD)
      call dom2cg(pdom%n_d(:), [ 0_LONG, 0_LONG, 0_LONG ], 0, pdom%pse(proc), pbd)

#ifdef VERBOSE
      call printinfo("[grid:init_grid]: all_cg finished. \o/")
#endif /* VERBOSE */

      call all_cg%reg_var(wa_n, AT_IGNORE)                 ! BEWARE: magic string across multiple files
      call all_cg%reg_var(fluid_n, AT_NO_B, flind%all)     !< Main array of all fluids' components, "u"
      call all_cg%reg_var(uh_n, AT_IGNORE, flind%all)      !< Main array of all fluids' components (for t += dt/2)
      call all_cg%reg_var(mag_n, AT_OUT_B, ndims)          !< Main array of magnetic field's components, "b"
      if (repeat_step) then
         call all_cg%reg_var(u0_n, AT_IGNORE, flind%all)   !< Copy of main array of all fluids' components
         call all_cg%reg_var(b0_n, AT_IGNORE, ndims)       !< Copy of main array of magnetic field's components
      endif

      nrq = 0
      cgl => all_cg%first
      do while (associated(cgl))
         cg => cgl%cg

         cg%u  => cg%get_na_ptr_4d(fluid_n)
         cg%b  => cg%get_na_ptr_4d(mag_n)
         cg%wa => cg%get_na_ptr(wa_n)

         if (allocated(cg%w)) then
            do d = xdim, zdim
               if (allocated(cg%w(1)%w_i_mbc(d, cg%nb)%mbc)) nrq = nrq + 2 * count(cg%w(1)%w_i_mbc(d, cg%nb)%mbc(:) /= INVALID) ! w(1) is probably fluid, but it can be any registered field
            enddo
         endif

         cgl => cgl%nxt
      enddo
      call inflate_req(nrq)

#ifdef ISO
      if (is_multicg) call die("[grid:init_cs_iso2] multiple grid pieces per procesor not fully implemented yet") !nontrivial maxval

      call all_cg%reg_var(cs_i2_n, AT_NO_B) ! BEWARE: magic string across multiple files

      cgl => all_cg%first
      do while (associated(cgl))
         cgl%cg%cs_iso2 => cgl%cg%get_na_ptr(cs_i2_n)
         cgl%cg%cs_iso2(:,:,:) = maxval(flind%all_fluids(:)%cs2)   ! set cs2 with sane values
         cgl => cgl%nxt
      enddo
#endif /* ISO */

#ifdef VERBOSE
      call printinfo("[grid:init_grid]: cg finished. \o/")
#endif /* VERBOSE */

   end subroutine init_grid

!> \brief Create new cg according to domain decomposition data and add them to appropriate lists
   subroutine dom2cg(n_d, offset, level, local_decomposition, patch)

      use constants, only: ndims
      use domain,    only: cuboids, pdom

      implicit none

      integer(kind=4), dimension(ndims), intent(in) :: n_d
      integer(kind=8), dimension(ndims), intent(in) :: offset
      integer, intent(in) :: level
      type(cuboids), pointer, intent(in) :: local_decomposition
      type(cg_list_patch), pointer, intent(inout) :: patch

      integer :: g

      do g = lbound(local_decomposition%sel(:,:,:), dim=1), ubound(local_decomposition%sel(:,:,:), dim=1)

         ! create the new element in the patch and initialize it
         call patch%add
         call patch%last%cg%init(pdom, g)

         ! add to the other lists
         call all_cg%add(patch%last%cg)

         call base_lev%add(patch%last%cg)

         call leaves%add(patch%last%cg)

      enddo

      base_lev%lev = level
      patch%n_d = n_d
      patch%off = offset
      patch%parent => null()
      patch%children => null()

   end subroutine dom2cg

!> \brief Update the list of leaves
! subroutine update leaves
! end subroutine update leaves

!> \brief deallocate everything
   subroutine cleanup_grid

      use gc_list, only: cg_list_element

      implicit none

      type(cg_list_element), pointer :: cgl, erase

      cgl => all_cg%first
      do while (associated(cgl))

         call cgl%cg%cleanup
         erase => cgl
         cgl => cgl%nxt

         call all_cg%un_link(erase)
         deallocate(erase%cg)
         deallocate(erase)
      enddo

!!$      deallocate(levels)

   end subroutine cleanup_grid

end module grid
