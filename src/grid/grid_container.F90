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

!> \brief Module providing the full, usable grid container type and its methods that don't fit to any abstract subtypes of grid container

module grid_cont

   use constants,         only: LO, HI
   use grid_cont_bnd,     only: segment
   use grid_cont_prolong, only: grid_container_prolong_t
#if defined(GRAV) && defined(NBODY)
   use particle_types,    only: particle_set
#endif /* GRAV && NBODY */

   implicit none

   private
   public :: grid_container

   !< \brief target list container for prolongations, restrictions and boundary exchanges
   type :: tgt_list
      type(segment), dimension(:), allocatable :: seg  !< a segment of data to be received or sent
   end type tgt_list

   !> \brief Everything required for autonomous computation of a single sweep on a portion of the domain on a single process
   type, extends(grid_container_prolong_t) :: grid_container

      ! Prolongation and restriction

      ! these lists are initialized and maintained in cg_level_connected
      ! perhaps parts of vertical_bf_prep can be moved here
      ! vertical_b_prep has to stay in cg_level_connected because we don't have dot here
      type(tgt_list) :: ri_tgt        !< description of incoming restriction data (this should be a linked list)
      type(tgt_list) :: ro_tgt        !< description of outgoing restriction data
      type(tgt_list) :: pi_tgt        !< description of incoming prolongation data
      type(tgt_list) :: po_tgt        !< description of outgoing prolongation data
      type(tgt_list) :: pib_tgt       !< description of incoming boundary prolongation data
      type(tgt_list) :: pob_tgt       !< description of outgoing boundary prolongation data
      type(tgt_list) :: rif_tgt       !< description of fluxes incoming from fine grid
      type(tgt_list) :: rof_tgt       !< description of fluxes outgoing to coarse grid

      ! Particles
#if defined(GRAV) && defined(NBODY)
      type(particle_set) :: pset                                  !< set of particles that belong to this grid part
#endif /* GRAV && NBODY */

      ! Misc
      integer(kind=8) :: SFC_id                                   !< position of the grid on space-filling curve
      integer :: membership                                       !< How many cg lists use this grid piece?
      logical :: ignore_prolongation                              !< When .true. do not upgrade interior with incoming prolonged values
      logical :: is_old                                           !< .true. if a given grid existed prior to  upgrade_refinement call
      logical :: processed                                        !< for use in sweeps.F90

   contains

      procedure          :: init_gc                               !< Initialization
      procedure          :: cleanup                               !< Deallocate all internals
      procedure          :: update_leafmap                        !< Check if the grid container has any parts covered by finer grids and update appropriate map

   end type grid_container

contains

!> \brief This method sets up remaining grid container variables and arrays.

   subroutine init_gc(this, my_se, grid_id, l)

      use constants,        only: PIERNIK_INIT_DOMAIN, LO
      use dataio_pub,       only: die, code_progress
      use level_essentials, only: level_t
      use ordering,         only: SFC_order

      implicit none

      class(grid_container), target,   intent(inout) :: this  ! intent(out) would silently clear everything, that was already set
                                                              ! (also the fields in types derived from grid_container)
      integer(kind=8), dimension(:,:), intent(in)    :: my_se    !< my segment
      integer,                         intent(in)    :: grid_id  !< ID which should be unique across level
      class(level_t), pointer,         intent(in)    :: l        !< level essential data

      if (code_progress < PIERNIK_INIT_DOMAIN) call die("[grid_container:init_gc] MPI not initialized.")

      call this%init_gc_base(my_se, grid_id, l)
      call this%init_gc_bnd
      call this%add_all_na
      call this%init_gc_prolong
#ifdef NBODY
      call this%pset%init()
#endif /* NBODY */

      this%membership = 1
      this%SFC_id     = SFC_order(this%my_se(:, LO) - l%off)

      call this%flag%init
      this%ignore_prolongation = .false.
      this%is_old = .false.
      this%has_previous_timestep = .false.

   end subroutine init_gc

!> \brief Routines that deallocates all internals of the grid container

   subroutine cleanup(this)

      implicit none

      class(grid_container), intent(inout) :: this !< object invoking type-bound procedure

      integer :: b, g
      integer, parameter :: nseg = 4*2
      type(tgt_list), dimension(nseg) :: rpio_tgt

      call this%cleanup_base
      call this%cleanup_na
      call this%cleanup_bnd
      call this%cleanup_prolong
#ifdef NBODY
      call this%pset%cleanup
#endif /* NBODY */

      rpio_tgt(1:nseg) = [ this%ri_tgt,  this%ro_tgt,  this%pi_tgt,  this%po_tgt, &
           &               this%pib_tgt, this%pob_tgt, this%rif_tgt, this%rof_tgt ]
      do b = 1, nseg
         if (allocated(rpio_tgt(b)%seg)) then
            do g = lbound(rpio_tgt(b)%seg, dim=1), ubound(rpio_tgt(b)%seg, dim=1)
               if (allocated(rpio_tgt(b)%seg(g)%buf)) deallocate(rpio_tgt(b)%seg(g)%buf)
            enddo
            deallocate(rpio_tgt(b)%seg)
         endif
      enddo

   end subroutine cleanup

!> \brief Check if the grid container has any parts covered by finer grids and update appropriate map

   subroutine update_leafmap(this)

      use constants, only: xdim, ydim, zdim, LO, HI

      implicit none

      class(grid_container), intent(inout) :: this !< object invoking type-bound procedure

      integer(kind=8), dimension(xdim:zdim, LO:HI) :: se
      integer :: g

      this%leafmap = .true.
      if (allocated(this%ri_tgt%seg)) then
         do g = lbound(this%ri_tgt%seg(:), dim=1), ubound(this%ri_tgt%seg(:), dim=1)
            se(:, :) = this%ri_tgt%seg(g)%se(:, :)
            this%leafmap(se(xdim, LO):se(xdim, HI), se(ydim, LO):se(ydim, HI), se(zdim, LO):se(zdim, HI)) = .false.
         enddo
      endif

   end subroutine update_leafmap

end module grid_cont
