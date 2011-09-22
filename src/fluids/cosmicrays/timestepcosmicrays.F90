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
!>
!! \brief (MH) Computation of %timestep for diffusive Cosmic Ray transport
!<

module timestepcosmicrays

! pulled by COSM_RAYS

   implicit none

   private
   public :: dt_crs, timestep_crs

   real, save :: dt_crs = huge(1.)

contains

!>
!! \details On a static grid with simple domain decompositions the fFollowing subroutine evaluates some constants, there is no need to run it
!! more than once, apart from wasting CPU cycles.
!<
   subroutine timestep_crs(cg)

      use constants,           only: one, half
      use grid,                only: all_cg
      use grid_cont,           only: grid_container
      use initcosmicrays,      only: cfl_cr, K_crs_paral, K_crs_perp, use_split
#ifdef MULTIGRID
      use multigrid_diffusion, only: diff_explicit, diff_tstep_fac, diff_dt_crs_orig
#endif /* MULTIGRID */

      implicit none

      type(grid_container), pointer, intent(in) :: cg

      logical, save :: frun = .true.
      real :: dt

      if (all_cg%cnt == 1 .and. .not. frun) return
      ! with multiple cg% there are few cg%dxmn to be checked
      ! with AMR minval(cg%dxmn) may change with time

      if (maxval(K_crs_paral+K_crs_perp) <= 0) then
         dt_crs = huge(one)
      else
         dt = cfl_cr * half/maxval(K_crs_paral+K_crs_perp)
         if (cg%dxmn < sqrt(huge(one))/dt) then
            dt = dt * cg%dxmn**2
#ifdef MULTIGRID
            diff_dt_crs_orig = min(dt_crs, dt)
            if (.not. (use_split .or. diff_explicit)) dt = dt * diff_tstep_fac ! enlarge timestep for non-explicit diffusion
#endif /* MULTIGRID */
            dt_crs = min(dt_crs, dt)
         endif
      endif

      frun = .false.

   end subroutine timestep_crs

end module timestepcosmicrays
