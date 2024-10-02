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
!! \brief Computation of %timestep for diffusive Cosmic Ray transport
!<

module timestepcosmicrays

! pulled by COSM_RAYS

   implicit none

   private
   public :: dt_crs, timestep_crs

   real :: dt_crs = huge(1.)

contains

!>
!! \details On a static grid with simple domain decompositions the following subroutine evaluates some constants, there is no need to run it
!! more than once, apart from wasting CPU cycles.
!<
   subroutine timestep_crs(dt)

      use allreduce,           only: piernik_MPI_Allreduce
      use cg_leaves,           only: leaves
      use cg_list,             only: cg_list_element
      use constants,           only: pMIN
      use domain,              only: is_multicg
      use initcosmicrays,      only: def_dtcrs, K_crs_valid, diff_max_lev
#ifdef MULTIGRID
      use multigrid_diffusion, only: diff_explicit, diff_tstep_fac, diff_dt_crs_orig
#endif /* MULTIGRID */
#ifdef CRESP
      use timestep_cresp,      only: dt_cre, cresp_timestep
#endif /* CRESP */

      implicit none

      real, intent(inout)            :: dt
      type(cg_list_element), pointer :: cgl

      logical, save                  :: frun = .true.

#ifdef CRESP
      call cresp_timestep
      dt = min(dt, dt_cre)
#endif /* CRESP */

      if (.not. K_crs_valid) return

      if ((is_multicg .or. frun)) then
      ! with multiple cg% there are few cg%dxmn to be checked
      ! with AMR minval(cg%dxmn) may change with time

         dt_crs = huge(1.)
         cgl => leaves%first
         do while (associated(cgl))
            if (cgl%cg%l%id <= diff_max_lev) &
                 dt_crs = min(dt_crs, def_dtcrs * cgl%cg%dxmn2)
            cgl => cgl%nxt
         enddo
         call piernik_MPI_Allreduce(dt_crs, pMIN)

#ifdef MULTIGRID
         diff_dt_crs_orig = dt_crs
#ifdef CRESP
         diff_dt_crs_orig = min(dt_crs, dt_cre)
#endif /* CRESP */
         if (.not. diff_explicit) dt_crs = diff_dt_crs_orig * diff_tstep_fac ! enlarge timestep for non-explicit diffusion
#endif /* MULTIGRID */
         frun = .false.
      endif

      dt = min(dt, dt_crs)

   end subroutine timestep_crs

end module timestepcosmicrays
