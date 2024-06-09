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
!! \brief Unified refinement criterion based on particle count
!!
!! Putting refinement on cg where a lot of particles reside can be justified in two ways:
!! * If there are a lot of particles, then maybe something worth resolving is going on.
!! * Refining crowded cg may let improve load balancing a bit because chances are that the resulting particle sets will be less numerous.
!<

module unified_ref_crit_nbody

   use unified_ref_crit_filter, only: urc_filter

   implicit none

   private
   public :: urc_nbody

!> \brief The type for refinement criterion based on particle count

   type, extends(urc_filter) :: urc_nbody
   contains
      procedure :: mark => mark_nbody
   end type urc_nbody

   interface urc_nbody
      procedure :: init
   end interface urc_nbody

contains

!> \brief A simple constructor for particle refinement

   function init(nbody_ref) result(this)

      implicit none

      integer(kind=4), intent(in) :: nbody_ref  !< maximum number of particles per cg

      type(urc_nbody) :: this  !< an object to be constructed

      this%ref_thr   = nbody_ref
      this%plotfield = .false.

   end function init

!>
!! \brief implementation of refinement criterion based on the number of particles
!!
!! Mark for refinement only selected octants (quadrants, halves in lower dimensions), where the particle count exceeds the threshold
!!
!! This routine has to conform to unified_ref_crit_user::mark_urc_user
!<

   subroutine mark_nbody(this, cg)

      use grid_cont,      only: grid_container
#if defined(GRAV) && defined(NBODY)
      use constants,      only: ndims, xdim, ydim, zdim, LO, HI, half, I_ZERO, I_ONE, I_TWO
      !use dataio_pub,     only: warn
      use particle_types, only: particle
#endif /* GRAV && NBODY */

      implicit none

      class(urc_nbody),              intent(inout) :: this  !< an object invoking the type-bound procedure
      type(grid_container), pointer, intent(inout) :: cg    !< current grid piece

#if defined(GRAV) && defined(NBODY)
      integer, dimension(LO:HI, LO:HI, LO:HI) :: oct_cnt
      integer, dimension(ndims) :: ijk
      integer :: i, j, k
      type(particle), pointer :: part

      if (cg%pset%cnt + sum(cg%chld_pcnt) >= this%ref_thr/I_TWO**ndims) then
         oct_cnt = I_ZERO
         part => cg%pset%first
         do while (associated(part))
            ijk = LO
            if (part%pdata%phy) then
               where (part%pdata%pos > half*(cg%fbnd(:, LO) + cg%fbnd(:, HI))) ijk = HI
               oct_cnt(ijk(xdim), ijk(ydim), ijk(zdim)) = oct_cnt(ijk(xdim), ijk(ydim), ijk(zdim)) + I_ONE
            endif
            part => part%nxt
         enddo

         ! mark the octants with big particle count
         ! To prevent refinement blinking we have to set flags under those children who have enough particles.
         ! The cg%chld_pcnt has to be updated at this point (non-local operation).
         do i = LO, HI
            do j = LO, HI
               do k = LO, HI
                  ! Looks like some particles remain nonprolonged during IC. Later it works correctly.
                  ! This bug is also visible on the nbdn field of the 0-th dump
                  if (oct_cnt(i, j, k) + cg%chld_pcnt(i, j, k) > this%ref_thr) then
                     call cg%flag%set(cg%ijkse(xdim, i), cg%ijkse(ydim, j), cg%ijkse(zdim, k))
                     oct_cnt(i, j, k) = I_ZERO
                  endif
               enddo
            enddo
         enddo

         ! if (sum(oct_cnt) > this%ref_thr) call warn("[URC_nbody:mark_nbody] Too many particles left in non-refined region")

      endif
#else /*  !(GRAV && NBODY) */
      if (.false.) this%ref_thr = this%ref_thr + cg%grid_id
#endif /* GRAV && NBODY */

   end subroutine mark_nbody

end module unified_ref_crit_nbody
