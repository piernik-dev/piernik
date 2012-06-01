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

   use cg_list_lev, only: cg_list_level, cg_list_patch
   use gc_list,     only: cg_list

   implicit none

   private
   public :: init_grid, cleanup_grid, base_lev, leaves, top_lev

   type(cg_list_level), target  :: base_lev                           !< base level grid containers \todo restore "protected"
   type(cg_list), protected  :: leaves                                !< grid containers not fully covered by finer grid containers
   integer, parameter :: NBD = 1                                      !< at the moment the base domain may be composed of only one patch
   type(cg_list_patch), dimension(NBD), target, protected :: base_dom !< base level patches; \todo relax the NBD=1 restriction if we want something like L-shaped or more complex domains
   type(cg_list_level), pointer :: top_lev                            !< finest level of refinement

contains

!>
!! \brief Routine that allocates all grid containers and most important field arrays inside gcs
!<
   subroutine init_grid

      use constants,  only: PIERNIK_INIT_DOMAIN, AT_NO_B, AT_OUT_B, VAR_XFACE, VAR_YFACE, VAR_ZFACE, INVALID, &
           &                ndims, xdim, zdim, fluid_n, uh_n, mag_n, wa_n, u0_n, b0_n, base_level_id, base_level_offset
#ifdef ISO
      use constants,  only: cs_i2_n
#endif /* ISO */
      use dataio_pub, only: printinfo, die, code_progress
      use domain,     only: pdom
#ifdef ISO
      use domain,     only: is_multicg
#endif /* ISO */
      use fluidindex, only: flind
      use gc_list,    only: cg_list_element, all_cg
      use global,     only: repeat_step
      use grid_cont,  only: grid_container
      use mpisetup,   only: inflate_req

      implicit none

      integer :: nrq, d
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
      type(cg_list_patch),   pointer :: pbd
      integer(kind=4), parameter, dimension(ndims) :: xyz_face = [ VAR_XFACE, VAR_YFACE, VAR_ZFACE ]

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
      call dom2cg(pdom%n_d(:), base_level_offset, base_level_id, pbd)

!#ifdef VERBOSE
      call base_lev%print_segments
!#endif /* VERBOSE */

#ifdef VERBOSE
      call printinfo("[grid:init_grid]: all_cg finished. \o/")
#endif /* VERBOSE */

      call all_cg%reg_var(wa_n,                                                           multigrid=.true.)  !! Auxiliary array. Multigrid required only for CR diffusion
      call all_cg%reg_var(fluid_n, vital = .true., restart_mode = AT_NO_B,  dim4 = flind%all)                !! Main array of all fluids' components, "u"
      call all_cg%reg_var(uh_n,                                             dim4 = flind%all)                !! Main array of all fluids' components (for t += dt/2)
      call all_cg%reg_var(mag_n,   vital = .true., restart_mode = AT_OUT_B, dim4 = ndims, position=xyz_face) !! Main array of magnetic field's components, "b"
      if (repeat_step) then
         call all_cg%reg_var(u0_n,                                          dim4 = flind%all)                !! Copy of main array of all fluids' components
         call all_cg%reg_var(b0_n,                                          dim4 = ndims, position=xyz_face) !! Copy of main array of magnetic field's components
      endif

      all_cg%fi  = all_cg%ind_4d(fluid_n)
      all_cg%bi  = all_cg%ind_4d(mag_n)
      all_cg%wai = all_cg%ind(wa_n)

      nrq = 0
      cgl => all_cg%first
      do while (associated(cgl))
         cg => cgl%cg

         cg%u  => cg%w(all_cg%fi)%arr
         cg%b  => cg%w(all_cg%bi)%arr
         cg%wa => cg%q(all_cg%wai)%arr

         if (allocated(cg%w)) then
            do d = xdim, zdim
               if (allocated(cg%w(1)%w_i_mbc(d, pdom%nb)%mbc)) nrq = nrq + 2 * count(cg%w(1)%w_i_mbc(d, pdom%nb)%mbc(:) /= INVALID) ! w(1) is probably fluid, but it can be any registered field
            enddo
         endif

         cgl => cgl%nxt
      enddo
      call inflate_req(nrq)

#ifdef ISO
      if (is_multicg) call die("[grid:init_cs_iso2] multiple grid pieces per procesor not fully implemented yet") !nontrivial maxval

      call all_cg%reg_var(cs_i2_n, vital = .true., restart_mode = AT_NO_B)

      cgl => all_cg%first
      do while (associated(cgl))
         cgl%cg%cs_iso2 => cgl%cg%q(all_cg%ind(cs_i2_n))%arr
         cgl%cg%cs_iso2(:,:,:) = maxval(flind%all_fluids(:)%cs2)   ! set cs2 with sane values
         cgl => cgl%nxt
      enddo
#endif /* ISO */

#ifdef VERBOSE
      call printinfo("[grid:init_grid]: cg finished. \o/")
#endif /* VERBOSE */

   end subroutine init_grid

!> \brief Create new cg according to domain decomposition data and add them to appropriate lists
   subroutine dom2cg(n_d, offset, level, patch)

      use constants,     only: ndims, I_ONE, base_level_id
      use decomposition, only: divide_domain
      use dataio_pub,    only: warn, die
      use domain,        only: is_mpi_noncart, is_multicg, is_refined, is_uneven
      use gc_list,       only: all_cg
      use mpi,           only: MPI_IN_PLACE, MPI_COMM_NULL, MPI_LOGICAL, MPI_LOR
      use mpisetup,      only: proc, comm, mpi_err, master
      use types,         only: cdd

      implicit none

      integer(kind=4), dimension(ndims), intent(in)    :: n_d
      integer(kind=8), dimension(ndims), intent(in)    :: offset
      integer(kind=4),                   intent(in)    :: level
      type(cg_list_patch), pointer,      intent(inout) :: patch

      integer :: g

      patch%n_d = n_d
      patch%off = offset
      patch%parent => null()
      patch%children => null()

      if (level == base_level_id) then
         base_lev%lev = level
         base_lev%n_d = n_d
         patch%list_level => base_lev
      else
         patch%list_level => null()
         call die("[grid:dom2cg] Only base level is currently supported")
      endif

      if (.not. divide_domain(patch)) then
         if (master) call die("[grid:dom2cg] Domain decomposition failed")
      endif

      !\todo Analyze the decomposition and set up [ is_uneven, is_mpi_noncart, is_refined, ... ]
      is_multicg = (ubound(base_lev%pse(proc)%sel(:, :, :), dim=1) > 1)
      call MPI_Allreduce(MPI_IN_PLACE, is_multicg, I_ONE, MPI_LOGICAL, MPI_LOR, comm, mpi_err)
      if (is_multicg .and. cdd%comm3d /= MPI_COMM_NULL) call die("[grid:dom2cg] is_multicg cannot be used with comm3d")
      if (is_refined) then
         is_mpi_noncart = .true.
         is_multicg = .true.
      endif
      if (is_mpi_noncart) is_uneven = .true.

      ! bnd_[xyz][lr] now become altered according to local topology of processes
      if (is_multicg .and. master) call warn("[grid:dom2cg] Multiple blocks per process not fully implemented yet")
      if (is_refined) call die("[grid:dom2cg] Refinements are not implemented")

      do g = lbound(base_lev%pse(proc)%sel(:,:,:), dim=1), ubound(base_lev%pse(proc)%sel(:,:,:), dim=1)

         if (level == base_level_id) then
            call base_lev%init_new_cg(g)
         else
            call die("[grid:dom2cg] can only create base level yet")
         endif

         ! add to the other lists
         call all_cg%add(base_lev%last%cg)
         call patch%add(base_lev%last%cg)
         call leaves%add(base_lev%last%cg)

      enddo

      top_lev => base_lev

   end subroutine dom2cg

! \brief Update the list of leaves
! subroutine update leaves
! end subroutine update leaves

!> \brief deallocate everything
   subroutine cleanup_grid

      use gc_list, only: cg_list_element, all_cg

      implicit none

      integer :: d
      type(cg_list_element), pointer :: cgle

      call leaves%delete
      call base_lev%delete
      if (allocated(base_lev%pse)) deallocate(base_lev%pse)
      do d = lbound(base_dom, dim=1), ubound(base_dom, dim=1) ! currently we have only one base patch
         call base_dom(d)%delete
      enddo

      ! manually deallocate all grid containers first
      cgle => all_cg%first
      do while (associated(cgle))
         deallocate(cgle%cg)
         cgle => cgle%nxt
      enddo
      if (allocated(all_cg%q_lst)) deallocate(all_cg%q_lst)
      if (allocated(all_cg%w_lst)) deallocate(all_cg%w_lst)
      call all_cg%delete

   end subroutine cleanup_grid

end module grid
