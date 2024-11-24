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

   use cg_cost,           only: cg_cost_t
   use cg_cost_data,      only: cg_cost_data_t
   use constants,         only: LO, HI
   use grid_cont_bseg,    only: tgt_list
   use grid_cont_prolong, only: grid_container_prolong_t
#if defined(GRAV) && defined(NBODY)
   use constants,         only: ndims
   use particle_types,    only: particle_set
#endif /* GRAV && NBODY */

   implicit none

   private
   public :: grid_container

   !> \brief Everything required for autonomous computation of a single sweep on a portion of the domain on a single process
   type, extends(grid_container_prolong_t) :: grid_container

      ! Prolongation and restriction

      ! these lists are initialized and maintained in cg_level_connected
      ! perhaps parts of vertical_bf_prep can be moved here
      ! vertical_b_prep has to stay in cg_level_connected because we don't have dot here
      type(tgt_list) :: ri_tgt   !< description of incoming restriction data (this should be a linked list)
      type(tgt_list) :: ro_tgt   !< description of outgoing restriction data
      type(tgt_list) :: pi_tgt   !< description of incoming prolongation data
      type(tgt_list) :: po_tgt   !< description of outgoing prolongation data
      type(tgt_list) :: pib_tgt  !< description of incoming boundary prolongation data
      type(tgt_list) :: pob_tgt  !< description of outgoing boundary prolongation data
      type(tgt_list) :: rif_tgt  !< description of fluxes incoming from fine grid
      type(tgt_list) :: rof_tgt  !< description of fluxes outgoing to coarse grid

      ! Particles
#if defined(GRAV) && defined(NBODY)
      type(particle_set) :: pset                        !< set of particles that belong to this grid part
      real, dimension(:,:), allocatable :: psave        !< for timestep retry
      real, dimension(ndims, LO:HI) :: bnd_in, bnd_out  !< coordinates for qualifying particles in is_part_in_cg
      integer, dimension(LO:HI, LO:HI, LO:HI) :: chld_pcnt  !< particle counters in children
#endif /* GRAV && NBODY */

      ! Misc
      integer(kind=8) :: SFC_id          !< position of the grid on space-filling curve
      integer :: membership              !< How many cg lists use this grid piece?
      logical :: ignore_prolongation     !< When .true. do not upgrade interior with incoming prolonged values
      logical :: is_old                  !< .true. if a given grid existed prior to  upgrade_refinement call
      logical :: processed               !< for use in sweeps.F90
      type(cg_cost_t) :: costs           !< accumulate cg costs here for better work balance
      type(cg_cost_data_t) :: old_costs  !< accumulated cg costs from previous step for better work balance

   contains

      procedure :: init_gc               !< Initialization
      procedure :: cleanup               !< Deallocate all internals
      procedure :: update_leafmap        !< Check if the grid container has any parts covered by finer grids and update appropriate map
      procedure :: print_tgt             !< Print all tgt_lists (for debugging only)
      procedure :: is_sending_fc_flux    !< Returns .true. if this block has fine hydro flux to be sent to some coarse block in specified direction
      procedure :: is_receiving_fc_flux  !< Returns .true. if this block expects fine hydro flux to be received from some fine block in specified direction
      procedure :: count_particles       !< Returns the number of phy particles in this%pset
      procedure :: count_all_particles   !< Returns the number of all particles in this%pset

   end type grid_container

contains

!> \brief This method sets up remaining grid container variables and arrays.

   subroutine init_gc(this, my_se, grid_id, l)

      use constants,        only: PIERNIK_INIT_DOMAIN, LO, PPP_AMR, PPP_CG
      use dataio_pub,       only: die, code_progress
      use level_essentials, only: level_t
      use ordering,         only: SFC_order
      use ppp,              only: ppp_main
#if defined(GRAV) && defined(NBODY)
      use constants,        only: LEFT, RIGHT, xdim, ydim, zdim
      use particle_types,   only: npb
#endif /* GRAV && NBODY */
      implicit none

      class(grid_container), target,   intent(inout) :: this     ! intent(out) would silently clear everything, that was already set
                                                                 ! (also the fields in types derived from grid_container)
      integer(kind=8), dimension(:,:), intent(in)    :: my_se    !< my segment
      integer,                         intent(in)    :: grid_id  !< ID which should be unique across level
      class(level_t), pointer,         intent(in)    :: l        !< level essential data

      character(len=*), parameter :: na_label = "add_all_na"

      if (code_progress < PIERNIK_INIT_DOMAIN) call die("[grid_container:init_gc] MPI not initialized.")

      call this%init_gc_base(my_se, grid_id, l)

      call this%init_gc_ref
      call this%init_gc_fcflx

      call ppp_main%start(na_label, PPP_AMR + PPP_CG)
      call this%add_all_na
      call ppp_main%stop(na_label, PPP_AMR + PPP_CG)

      call this%init_gc_prolong
#ifdef NBODY
      call this%pset%init()

      this%bnd_out(:,LO) = [this%coord(LEFT, xdim)%r(this%ijkse(xdim,LO)-npb), this%coord(LEFT, ydim)%r(this%ijkse(ydim,LO)-npb), this%coord(LEFT, zdim)%r(this%ijkse(zdim,LO)-npb)]
      this%bnd_out(:,HI) = [this%coord(RIGHT,xdim)%r(this%ijkse(xdim,HI)+npb), this%coord(RIGHT,ydim)%r(this%ijkse(ydim,HI)+npb), this%coord(RIGHT,zdim)%r(this%ijkse(zdim,HI)+npb)]

      this%bnd_in(:,LO)  = [this%coord(LEFT, xdim)%r(this%ijkse(xdim,LO)+npb), this%coord(LEFT, ydim)%r(this%ijkse(ydim,LO)+npb), this%coord(LEFT, zdim)%r(this%ijkse(zdim,LO)+npb)]
      this%bnd_in(:,HI)  = [this%coord(RIGHT,xdim)%r(this%ijkse(xdim,HI)-npb), this%coord(RIGHT,ydim)%r(this%ijkse(ydim,HI)-npb), this%coord(RIGHT,zdim)%r(this%ijkse(zdim,HI)-npb)]
#endif /* NBODY */

      this%membership = 1
      this%SFC_id     = SFC_order(this%my_se(:, LO) - l%off)

      call this%flag%init
      this%ignore_prolongation = .false.
      this%is_old = .false.
      this%has_previous_timestep = .false.

      call this%costs%reset
      call this%old_costs%init

   end subroutine init_gc

!> \brief Routines that deallocates all internals of the grid container

   subroutine cleanup(this)

      implicit none

      class(grid_container), intent(inout) :: this  !< object invoking type-bound procedure

      call this%cleanup_base
      call this%cleanup_na
      call this%cleanup_bseg
      call this%cleanup_ref
      call this%cleanup_fcflx
      call this%cleanup_prolong
#ifdef NBODY
      call this%pset%cleanup
#endif /* NBODY */

      call this%ri_tgt%cleanup
      call this%ro_tgt%cleanup
      call this%pi_tgt%cleanup
      call this%po_tgt%cleanup
      call this%pib_tgt%cleanup
      call this%pob_tgt%cleanup
      call this%rif_tgt%cleanup
      call this%rof_tgt%cleanup

   end subroutine cleanup

!> \brief Check if the grid container has any parts covered by finer grids and update appropriate map

   subroutine update_leafmap(this)

      use constants, only: xdim, ydim, zdim, LO, HI

      implicit none

      class(grid_container), intent(inout) :: this  !< object invoking type-bound procedure

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

!> \brief Returns .true. if this block has fine hydro flux to be sent to some coarse block in specified direction

   pure logical function is_sending_fc_flux(this, cdim)

      use constants, only: LO, HI

      implicit none

      class(grid_container), intent(in) :: this
      integer(kind=4),       intent(in) :: cdim

      integer :: g

      is_sending_fc_flux = .false.

      if (allocated(this%rof_tgt%seg)) then
         associate ( seg => this%rof_tgt%seg )
            do g = lbound(seg, dim=1), ubound(seg, dim=1)
               if (seg(g)%se(cdim, LO) == seg(g)%se(cdim, HI)) then
                  is_sending_fc_flux = .true.
                  exit
               endif
            enddo
         end associate
      endif

   end function is_sending_fc_flux

!> \brief Returns .true. if this block expects fine hydro flux to be received from some fine block in specified direction

   pure logical function is_receiving_fc_flux(this, cdim)

      use constants, only: LO, HI

      implicit none

      class(grid_container), intent(in) :: this
      integer(kind=4),       intent(in) :: cdim

      integer :: g

      is_receiving_fc_flux = .false.

      if (allocated(this%rif_tgt%seg)) then
         associate ( seg => this%rif_tgt%seg )
            do g = lbound(seg, dim=1), ubound(seg, dim=1)
               if (seg(g)%se(cdim, LO) == seg(g)%se(cdim, HI)) then
                  is_receiving_fc_flux = .true.
                  exit
               endif
            enddo
         end associate
      endif

   end function is_receiving_fc_flux

!> \brief Print all tgt_lists (for debugging only)

   subroutine print_tgt(this)

      use constants, only: LO, HI, cwdlen
      use mpisetup,  only: proc

      implicit none

      class(grid_container), intent(in) :: this  !< object invoking type-bound procedure

      character(len=cwdlen) :: msg

      write(msg, '(a,i5,a,i2,a,i7,a,2(3i7,a))')"@", proc, " ^", this%l%id, " #", this%grid_id, " [", this%ijkse(:, LO), "]..[", this%ijkse(:, HI), "] "

      call pr_tgt(this%ri_tgt,  "ri " // trim(msg))
      call pr_tgt(this%ro_tgt,  "ro " // trim(msg))
      call pr_tgt(this%pi_tgt,  "pi " // trim(msg))
      call pr_tgt(this%po_tgt,  "po " // trim(msg))
      call pr_tgt(this%rif_tgt, "rif" // trim(msg))
      call pr_tgt(this%rof_tgt, "rof" // trim(msg))
      call pr_tgt(this%pib_tgt, "pib" // trim(msg))
      call pr_tgt(this%pob_tgt, "pob" // trim(msg))

   contains

      subroutine pr_tgt(tgt, label)

         use constants,  only: LO, HI
         use dataio_pub, only: printinfo, msg

         implicit none

         type(tgt_list),   intent(in) :: tgt   ! segment list to print
         character(len=*), intent(in) :: label ! identifier of segment type

         integer :: i
         logical, parameter :: stdout = .true.

         if (.not. allocated(tgt%seg)) then
            call printinfo(label // " not allocated", stdout)
         else if (size(tgt%seg) == 0) then
            call printinfo(label // " no entries", stdout)
         else
            do i = lbound(tgt%seg, 1), ubound(tgt%seg, 1)
               write(msg, '(2a,i5,2(a,3i7),a,i9)') label,  " -> @", tgt%seg(i)%proc, " [", tgt%seg(i)%se(:, LO), "]..[", tgt%seg(i)%se(:, HI), "] t", tgt%seg(i)%tag
               call printinfo(msg, stdout)
            enddo
         endif

      end subroutine pr_tgt

   end subroutine print_tgt

!> \brief Returns the number of phy particles in this%pset

   integer(kind=4) function count_particles(this) result(n_part)

#if !defined(GRAV) || !defined(NBODY)
      use constants, only: I_ZERO
#endif

      implicit none

      class(grid_container), intent(in) :: this  !< an object invoking the type-bound procedure

#if defined(GRAV) && defined(NBODY)
      n_part = this%pset%count()
#else
      n_part = I_ZERO * this%maxxyz  ! suppress compiler warnings
#endif

   end function count_particles

!> \brief Returns the number of phy particles in this%pset

   integer(kind=4) function count_all_particles(this) result(n_part)

#if !defined(GRAV) || !defined(NBODY)
      use constants, only: I_ZERO
#endif

      implicit none

      class(grid_container), intent(in) :: this  !< an object invoking the type-bound procedure

#if defined(GRAV) && defined(NBODY)
      n_part = this%pset%cnt
#else
      n_part = I_ZERO * this%maxxyz  ! suppress compiler warnings
#endif

   end function count_all_particles

end module grid_cont
