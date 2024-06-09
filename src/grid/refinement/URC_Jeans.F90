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
!! \brief Unified refinement criterion based on Jeans length
!!
!! The importance of proper resolution of Jeans length was presented in
!! Truelove et al., ApJ, 489, L179-L183, Bibcode: 1997ApJ...489L.179T
!!
!! Their suggestion of 4 cells per Jeans length was well justified on computers available in 1997.
!! Some decades later don't hesitate to ask for more, unless you're using your phone for computing.
!<

module unified_ref_crit_Jeans

   use unified_ref_crit_filter, only: urc_filter

   implicit none

   private
   public :: urc_jeans

!> \brief The type for Jeans length-based refinement criterion.

   type, extends(urc_filter) :: urc_jeans
   contains
      procedure :: mark => mark_Jeans
   end type urc_jeans

   interface urc_jeans
      procedure :: init
   end interface urc_jeans

contains

!> \brief A simple constructor for Jeans refinement

   function init(jeans_ref, jeans_plot) result(this)

      implicit none

      real,    intent(in) :: jeans_ref   !< minimum resolution in cells per Jeans wavelength
      logical, intent(in) :: jeans_plot  !< create an array to keep the value of Jeans resolution

      type(urc_jeans) :: this  !< an object to be constructed

      this%ref_thr   = jeans_ref
      this%plotfield = jeans_plot

   end function init

!>
!! \brief implementation of user-provided Jeans length refinement criterion
!!
!! This routine has to conform to unified_ref_crit_user::mark_urc_user
!!
!! \todo implement predictive marks based on velocity, something like:
!! where (cg%flag%map) cg%flag%map(x + v*2*n_updAMR) = .true.
!! count cases where internal refine flag goes beyond guardcells and suggest reducing n_updAMR
!<

   subroutine mark_Jeans(this, cg)

      use constants,  only: pi, GEO_XYZ, INVALID, dirtyH
      use dataio_pub, only: die
      use domain,     only: dom
      use fluidindex, only: iarr_all_sg, flind
      use func,       only: ekin !, emag
      use grid_cont,  only: grid_container
      use units,      only: newtong
#ifdef MAGNETIC
      use dataio_pub, only: warn
#endif /* MAGNETIC */

      implicit none

      class(urc_jeans),              intent(inout) :: this  !< an object invoking the type-bound procedure
      type(grid_container), pointer, intent(inout) :: cg    !< current grid piece

      real, dimension(:,:,:), pointer :: p3d
      integer :: f
      logical :: any_has_energy
#ifdef MAGNETIC
      logical, save :: warned = .false.

      if (.not. warned) call warn("[URC_Jeans:mark_Jeans] magnetic pressure is ignored (unimplemented yet)")
      warned = .true.
#endif /* MAGNETIC */

      if (dom%geometry_type /= GEO_XYZ) call die("[URC_Jeans:mark_Jeans] unsupported (non-cartesian) geometry")

      if (this%iplot /= INVALID) then
         p3d => cg%q(this%iplot)%arr
      else
         p3d => cg%wa
      endif

      ! assume that fluid and cs2 have updated boundaries
      if (associated(cg%cs_iso2)) then
         p3d(:,:,:) = sqrt(pi/newtong) / maxval(cg%dl) * &
              sqrt(cg%cs_iso2(:,:,:) / sum(cg%u(iarr_all_sg, :, :, :), dim=1))
      else
         any_has_energy = .false.
         do f = 1, flind%fluids
            if (flind%all_fluids(f)%fl%has_energy) any_has_energy = .true.
         enddo

         if (any_has_energy) then
            p3d(:,:,:) = 0.
         else
            p3d(:,:,:) = dirtyH
         endif

         ! l_J = sqrt( pi gam (gam - 1) (e_i - e_k + e_mag) / G ) / rho
         ! find the effective pressure of all components (non-selfgravitation as well)
         do f = 1, flind%fluids
            associate (fl => flind%all_fluids(f)%fl)
               if (fl%has_energy) p3d = p3d + fl%gam * fl%gam_1 * ( &
                    cg%u(fl%ien, :, :, :) - &
                    ekin(cg%u(fl%imx, :,:,:), cg%u(fl%imy, :,:,:), cg%u(fl%imz, :,:,:), cg%u(fl%idn, :,:,:)) &
                    )
            end associate
         enddo
         ! sum of ion/neu/dst
         p3d(:,:,:) = sqrt(pi/newtong) / maxval(cg%dl) * sqrt(p3d)/ sum(cg%u(iarr_all_sg, :, :, :), dim=1)
      endif

      call cg%flag%set(p3d < this%ref_thr)

   end subroutine mark_Jeans

end module unified_ref_crit_Jeans
