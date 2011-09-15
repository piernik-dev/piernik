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
!! \brief (DW) Module containing routines to specify required computational mesh.
!! \date January/February 2006
!!
!!
!! In this module a following namelist of parameters is specified:
!! \copydetails grid::init_grid
!<

module grid

   use gc_list,   only: cg_list

   implicit none

   private
   public :: init_grid, grid_mpi_boundaries_prep, cleanup_grid, all_cg !, base, leafs, levels

   type(cg_list), protected :: all_cg    !< all grid containers
!!$   type(cg_list), protected  :: base   !< base level grid containers
!!$   type(cg_list), protected  :: leafs  !< grid containers not covered by other grid containers
!!$   type(cg_list), dimension(:), allocatable, protected  :: levels !< grid containers grouped by levels

contains

!>
!! \brief Routine that allocates all grid containers and most important field arrays inside gc's
!<
   subroutine init_grid

      use constants,   only: PIERNIK_INIT_DOMAIN, AT_NO_B, ndims, zdim
      use dataio_pub,  only: printinfo, die, code_progress
      use diagnostics, only: my_allocate
      use domain,      only: dom
      use fluidindex,  only: flind
      use gc_list,     only: cg_list_element
      use global,      only: repeat_step
      use grid_cont,   only: grid_container
      use mpisetup,    only: proc

      implicit none

      integer :: g
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
      integer(kind=4), dimension(:), allocatable :: ind_arr

      if (code_progress < PIERNIK_INIT_DOMAIN) call die("[grid:init_grid] domain not initialized.")

#ifdef VERBOSE
      call printinfo("[grid:init_grid]: commencing...")
#endif /* VERBOSE */

      call all_cg%init
!!$      call base%init
!!$      call leafs%init
!!$      allocate(levels(1))
!!$      call levels(1)%init

      do g = 1, ubound(dom%pse(proc)%sel(:,:,:), dim=1)
         call all_cg%add
         call all_cg%last%cg%init(dom, g)
      enddo

#ifdef VERBOSE
      call printinfo("[grid:init_grid]: all_cg finished. \o/")
#endif /* VERBOSE */

      cgl => all_cg%first
      do while (associated(cgl))
         cg => cgl%cg

         allocate(ind_arr(ndims+1))

         ind_arr = [ flind%all, cg%n_(:) ]
         call cg%u%init(ind_arr)
         call cg%uh%init(ind_arr)
         if (repeat_step) call cg%u0%init(ind_arr)

         ind_arr = [ ndims, cg%n_(:) ]
         call cg%b%init(ind_arr)
         if (repeat_step) call cg%b0%init(ind_arr)

         deallocate(ind_arr)

         call cgl%cg%add_na("wa") ! BEWARE: magic string across multiple files
         cgl%cg%wa => cgl%cg%get_na_ptr("wa")
#ifdef GRAV
         allocate(ind_arr(1))
         ind_arr = [cg%n_(zdim)]
         call my_allocate(cg%dprof, ind_arr, "dprof")
         deallocate(ind_arr)
#endif /* GRAV */
         cgl => cgl%nxt
      enddo

#ifdef ISO
      if (all_cg%cnt > 1) call die("[grid:init_cs_iso2] multiple grid pieces per procesor not fully implemented yet") !nontrivial maxval

      cgl => all_cg%first
      do while (associated(cgl))
         call cgl%cg%add_na("cs_iso2", AT_NO_B) ! BEWARE: magic string across multiple files
         cgl%cg%cs_iso2 => cgl%cg%get_na_ptr("cs_iso2")
         cgl%cg%cs_iso2(:,:,:) = maxval(flind%all_fluids(:)%cs2)   ! set cs2 with sane values
         cgl => cgl%nxt
      enddo
#endif /* ISO */

#ifdef VERBOSE
      call printinfo("[grid:init_grid]: cg finished. \o/")
#endif /* VERBOSE */

   end subroutine init_grid

!>
!! \brief deallocate everything
!<

   subroutine cleanup_grid

      use gc_list, only: cg_list_element

      implicit none

      type(cg_list_element), pointer :: cgl

      cgl => all_cg%first
      do while (associated(cgl))

         call cgl%cg%u%clean()
         call cgl%cg%u0%clean()
         call cgl%cg%uh%clean()

         call cgl%cg%b%clean()
         call cgl%cg%b0%clean()

#ifdef GRAV
         if (allocated(cgl%cg%dprof)) deallocate(cgl%cg%dprof)
#endif /* GRAV */

         call cgl%cg%cleanup
         deallocate(cgl%cg)
         cgl => cgl%nxt

         if (associated(cgl)) call all_cg%del(cgl%prv)
      enddo
      deallocate(all_cg%last)

!!$      deallocate(levels)

   end subroutine cleanup_grid

!>
!! \brief Set up subsets of u,b and sgp arrays for MPI communication
!!
!! \todo this should be called from cg%init only
!<

   subroutine grid_mpi_boundaries_prep

      use constants,  only: PIERNIK_INIT_GRID, FLUID, LO, HI
      use dataio_pub, only: die, code_progress
      use gc_list,    only: cg_list_element
      use mpisetup,   only: proc

      implicit none

      type(cg_list_element), pointer :: cgl

      if (code_progress < PIERNIK_INIT_GRID) call die("[grid:grid_mpi_boundaries_prep] grid or fluids not initialized.")
      if (all_cg%cnt > 1) call die("[grid:grid_mpi_boundaries_prep] Multiple blocks per process not implemented yet")

      cgl => all_cg%first
      do while (associated(cgl))
         ! find neighbours and set up the MPI containers

         call cgl%cg%mpi_bnd_types

         cgl => cgl%nxt
      enddo

   end subroutine grid_mpi_boundaries_prep

end module grid
