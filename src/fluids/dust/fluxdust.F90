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
#define RNG 2:nm
!/*
!>
!! \brief (MH/JD) [R] Computation of %fluxes for the dust fluid
!!
!!The flux functions for dust are given by
!!
!!\f[
!!  \vec{F}{(\vec{u})} =
!!  \left(\begin{array}{c}
!!    \rho v_x \\
!!    \rho v_x^2 \\
!!    \rho v_x v_y \\
!!    \rho v_x v_z
!!  \end{array}\right),
!!  \qquad
!!  \vec{G}{(\vec{u})} =
!!  \left(\begin{array}{c}
!!    \rho v_y \\
!!    \rho v_y v_x \\
!!    \rho v_y^2 \\
!!    \rho v_y v_z
!!  \end{array}\right),
!!\qquad
!!  \vec{H}{(\vec{u})} =
!!  \left(\begin{array}{c}
!!    \rho v_z \\
!!    \rho v_z v_x\\
!!    \rho v_z v_y \\
!!    \rho v_z^2
!!  \end{array}\right),
!!\f]
!!
!<
!*/
module fluxdust
! pulled by DUST
   implicit none
   private
   public :: flux_dst

contains

   subroutine flux_dst(fluxd, cfrd, uud, n, vx, ps, bb, cs_iso2)

      use fluidindex, only: idn, imx, imy, imz
#ifdef GLOBAL_FR_SPEED
      use timestep,   only: c_all
#endif /* GLOBAL_FR_SPEED */

      implicit none
      integer(kind=4), intent(in)                  :: n         !< number of cells in the current sweep
      real, dimension(:,:), intent(inout), pointer :: fluxd     !< flux of dust
      real, dimension(:,:), intent(inout), pointer :: cfrd      !< freezing speed for dust
      real, dimension(:,:), intent(in),    pointer :: uud       !< part of u for dust
      real, dimension(:),   intent(inout), pointer :: vx        !< velocity of dust fluid for current sweep
      real, dimension(:),   intent(inout), pointer :: ps        !< pressure of dust fluid for current sweep
      real, dimension(:,:), intent(in),    pointer :: bb        !< magnetic field x,y,z-components table
      real, dimension(:),   intent(in),    pointer :: cs_iso2   !< local isothermal sound speed squared (optional)

      ! locals
!      real               :: minvx, maxvx, amp
      real, dimension(size(vx)) :: absvx
      integer                   :: nm

      nm = n-1

      ps(:)  = 0.0
      vx(RNG)=uud(imx,RNG)/uud(idn,RNG) ; vx(1) = vx(2); vx(n) = vx(nm)

      fluxd(idn,RNG)=uud(imx,RNG)
      fluxd(imx,RNG)=uud(imx,RNG)*vx(RNG)
      fluxd(imy,RNG)=uud(imy,RNG)*vx(RNG)
      fluxd(imz,RNG)=uud(imz,RNG)*vx(RNG)

      fluxd(:,1) = fluxd(:,2); fluxd(:,n) = fluxd(:,nm)

#ifdef LOCAL_FR_SPEED

      ! The freezing speed is now computed locally (in each cell)
      !  as in Trac & Pen (2003). This ensures much sharper shocks,
      !  but sometimes may lead to numerical instabilities
!     minvx = minval(vx(RNG))
!     maxvx = maxval(vx(RNG))
!     amp   = (maxvx-minvx)*0.5
!      cfrd(1,RNG) = max(sqrt(vx(RNG)**2+cfr_smooth*amp),small)
      absvx = abs(vx)
      cfrd(idn,RNG) = max( absvx(1:n-2), absvx(2:nm), absvx(3:n) )

      cfrd(idn,1) = cfrd(idn,2); cfrd(idn,n) = cfrd(idn,nm)

      cfrd(imx,:) = cfrd(idn,:)
      cfrd(imy,:) = cfrd(idn,:)
      cfrd(imz,:) = cfrd(idn,:)

#endif /* LOCAL_FR_SPEED */

#ifdef GLOBAL_FR_SPEED
      ! The freezing speed is now computed globally
      !  (c=const for the whole domain) in subroutine 'timestep'

      !  cfrd(:,:) = flind%dst%snap%c
      cfrd(:,:) = c_all
#endif /* GLOBAL_FR_SPEED */
      return
      if (.false.) write(0,*) bb, cs_iso2

   end subroutine flux_dst

end module fluxdust
