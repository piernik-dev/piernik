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
#include "piernik.def"
!>
!! \brief (MH) Computation of %timestep for diffusive Cosmic Ray transport
!!
!<
module timestepcosmicrays
   implicit none
   private
   public :: dt_crs, timestep_crs

   real   :: dt_crs

contains
   ! Following subroutine evaluates some constants, there is no need to run it
   ! more than once, apart from wasting CPU cycles...
   subroutine timestep_crs
      use grid,           only: dxmn
      use initcosmicrays, only: cfl_cr, K_crs_paral, K_crs_perp

      implicit none
      logical, save :: frun = .true.

      if (.not.frun) return

      if (maxval(K_crs_paral+K_crs_perp) <= 0) then
         dt_crs = huge(1.0)
      else
         dt_crs = cfl_cr * 0.5/maxval(K_crs_paral+K_crs_perp)
         if (dxmn < sqrt(huge(1.0))/dt_crs) then
            dt_crs = dt_crs * dxmn**2
         else
            dt_crs = huge(1.0)
         endif
      endif

      frun = .false.
   end subroutine timestep_crs

end module timestepcosmicrays
