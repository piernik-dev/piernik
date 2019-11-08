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
!    along with PIERNIK.  If not, see http://www.gnu.org/licenses/.
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
!! \brief Unified refinement criteriun based on Jeans length
!!
!! The importance of proper resolution of Jeans length was presented in
!! Truelove et al., ApJ, 489, L179-L183, Bibcode: 1997ApJ...489L.179T
!!
!! Their suggestion of 4 cell per length was well justified on computers available in 1997.
!! Some decades later don't hesitate to ask for more, unless you're using your phone for computing.
!<

module unified_ref_crit_Jeans

   use unified_ref_crit_filter, only: urc_filter

   implicit none

   private
   public :: urc_Jeans

!> \brief Things that should be common for all refinement criteria based on filters that decide whether local conditions deserve refinement or derefinement.

   type, extends(urc_filter) :: urc_Jeans
   contains
      procedure :: mark => mark_Jeans
   end type urc_Jeans

   interface urc_Jeans
      procedure :: init
   end interface urc_Jeans

contains

!> \brief A simple constructor for Jeans refinement

   function init(jeans_ref, jeans_plot) result(this)

      use constants, only: refinement_factor

      implicit none

      real,    intent(in) :: jeans_ref   !< minimum resolution in cells per Jeans wavelength
      logical, intent(in) :: jeans_plot  !< create an array to keep the value of Jeans resolution

      type(urc_Jeans) :: this  !< an object to be constructed

      real, parameter :: safe_deref = refinement_factor * 1.25  !< if it proves to be not save then implement it as a problem.par parameter

      this%ref_thr   = jeans_ref
      this%deref_thr = safe_deref * jeans_ref
      this%plotfield = jeans_plot

   end function init

!>
!! \brief implementation of user-provided Jeans length refinement criterion
!!
!! This routine has to conform to unified_ref_crit_user::mark_urc_user
!<

   subroutine mark_Jeans(this, cg)

      use constants,  only: pi, GEO_XYZ, INVALID
      use dataio_pub, only: die
      use domain,     only: dom
      use fluidindex, only: iarr_all_sg
      use grid_cont,  only: grid_container
      use units,      only: newtong

      implicit none

      class(urc_Jeans),              intent(inout) :: this  !< an object invoking the type-bound procedure
      type(grid_container), pointer, intent(inout) :: cg    !< current grid piece

      real, dimension(:,:,:), pointer :: p3d

      if (dom%geometry_type /= GEO_XYZ) call die("[unified_ref_crit_Jeans:mark_Jeans] unsupported (non-cartesian) geometry")

      if (this%iplot /= INVALID) then
         p3d => cg%q(this%iplot)%arr
      else
         p3d => cg%wa
      endif

      p3d(:,:,:) = sqrt(pi/newtong) / maxval(cg%dl) * &  ! assumes that fluid and cs2 have updated boundaries
           sqrt(cg%cs_iso2(:,:,:) / sum(cg%u(iarr_all_sg, :, :, :), dim=1))

      where (p3d < this%ref_thr)
         cg%refinemap = .true.
      endwhere
      if (any(p3d < this%deref_thr)) cg%refine_flags%derefine = .false.

   end subroutine mark_Jeans

end module unified_ref_crit_Jeans
