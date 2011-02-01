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
#define RNG 2:n-1
!/*
!>
!! \brief (MH/JD) [R] Computation of %fluxes for the ionized fluid
!!
!!The flux functions for ionized fluid are given by
!!\f[
!!  \vec{F}{(\vec{u})} =
!!  \left(\begin{array}{c}
!!    \rho v_x \\
!!    \rho v_x^2 + p_* - B_x^2 \\
!!    \rho v_x v_y - B_x B_y\\
!!    \rho v_x v_z - B_x B_z\\
!!    (e + p)v_x - \vec{B} \cdot \vec{v} \; B_x
!!  \end{array}\right),
!!  \qquad
!!  \vec{G}{(\vec{u})} =
!!  \left(\begin{array}{c}
!!    \rho v_y \\
!!    \rho v_y v_x  - B_y B_x\\
!!    \rho v_y^2 + p_* - B_y^2 \\
!!    \rho v_y v_z  - B_y B_z\\
!!    (e + p)v_y - \vec{B} \cdot \vec{v} \; B_y
!!  \end{array}\right),
!!\qquad
!!  \vec{H}{(\vec{u})} =
!!  \left(\begin{array}{c}
!!    \rho v_z \\
!!    \rho v_z v_x  - B_z B_x\\
!!    \rho v_z v_y  - B_z B_y\\
!!    \rho v_z^2 + p_* - B_z^2 \\
!!    (e + p)v_z - \vec{B} \cdot \vec{v} \; B_z
!!  \end{array}\right),
!!\f]
!!where \f$p_* = p + B^2/2\f$, \f$e= e_{th} + \frac{1}{2} \rho v^2 + B^2/2\f$,  are the total pressure and total energy density,
!!while \f$e_{th}\f$ is thermal energy density and  \f$e_{mag} = B^2/2\f$ is the magnetic energy density.
!<
!*/
module fluxionized
! pulled by IONIZED
   implicit none
   private
   public :: flux_ion

contains

   subroutine flux_ion(fluxi, cfri, uui, n, vx, ps, bb, cs_iso2)

      use constants,       only: small
      use dataio_pub,      only: die
      use fluidindex,      only: idn, imx, imy, imz, ien, flind, ibx, iby, ibz
      use mpisetup,        only: cfr_smooth

      implicit none
      integer, intent(in)                        :: n           !< number of cells in the current sweep
      real, dimension(:,:), intent(in),  pointer :: uui         !< part of u for ionized fluid
      real, dimension(:,:), intent(out), pointer :: fluxi       !< flux of ionized fluid
      real, dimension(:,:), intent(out), pointer :: cfri        !< freezing speed for ionized fluid
      real, dimension(:,:), intent(in),  pointer :: bb          !< magnetic field x,y,z-components table
      real, dimension(:),   intent(out), pointer :: vx          !< velocity of ionized fluid for current sweep
      real, dimension(:),   intent(out), pointer :: ps          !< gas pressure of ionized fluid for current sweep
      real, dimension(:),   intent(in),  pointer :: cs_iso2     !< local isothermal sound speed (optional)

      real, dimension(n)               :: p           !< thermal pressure of ionized fluid
      real, dimension(n)               :: pmag        !< pressure of magnetic field
#ifdef LOCAL_FR_SPEED
      real                             :: minvx       !<
      real                             :: maxvx       !<
      real                             :: amp         !<
#endif /* LOCAL_FR_SPEED */

      fluxi = 0.0; cfri = 0.0; vx = 0.0; ps = 0.0; p = 0.0; pmag = 0.0

#ifdef MAGNETIC
      pmag(RNG)=0.5*( bb(ibx,RNG)**2 + bb(iby,RNG)**2 +bb(ibz,RNG)**2 )
#else /* !MAGNETIC */
      pmag(:) = 0.0
#endif /* !MAGNETIC */
      vx(RNG)=uui(imx,RNG)/uui(idn,RNG)

#ifndef ISO_LOCAL
      if (associated(cs_iso2)) call die("[fluxionized:flux_ion] cs_iso2 should not be present")
#endif /* !ISO_LOCAL */

#ifdef ISO
#ifdef ISO_LOCAL
      p(RNG) = cs_iso2(RNG) * uui(idn,RNG)
#else /* !ISO_LOCAL */
      p(RNG) = flind%ion%cs2 * uui(idn,RNG)
#endif /* !ISO_LOCAL */
      ps(RNG)= p(RNG) + pmag(RNG)
#else /* !ISO */
      ps(RNG)=(uui(ien,RNG) - &
      0.5*( uui(imx,RNG)**2 + uui(imy,RNG)**2 + uui(imz,RNG)**2 ) &
          / uui(idn,RNG))*(flind%ion%gam_1) + (2.0-flind%ion%gam)*pmag(RNG)
      p(RNG) = ps(RNG)- pmag(RNG)
#endif /* !ISO */

      fluxi(idn,RNG)=uui(imx,RNG)
      fluxi(imx,RNG)=uui(imx,RNG)*vx(RNG)+ps(RNG) - bb(ibx,RNG)**2
      fluxi(imy,RNG)=uui(imy,RNG)*vx(RNG)-bb(iby,RNG)*bb(ibx,RNG)
      fluxi(imz,RNG)=uui(imz,RNG)*vx(RNG)-bb(ibz,RNG)*bb(ibx,RNG)
#ifndef ISO
      fluxi(ien,RNG)=(uui(ien,RNG)+ps(RNG))*vx(RNG)-bb(ibx,RNG)*(bb(ibx,RNG)*uui(imx,RNG) &
                +bb(iby,RNG)*uui(imy,RNG)+bb(ibz,RNG)*uui(imz,RNG))/uui(idn,RNG)
#endif /* !ISO */

#ifdef LOCAL_FR_SPEED

      ! The freezing speed is now computed locally (in each cell)
      !  as in Trac & Pen (2003). This ensures much sharper shocks,
      !  but sometimes may lead to numerical instabilities
      minvx = minval(vx(RNG))
      maxvx = maxval(vx(RNG))
      amp   = 0.5*(maxvx-minvx)
#ifdef ISO
      cfri(1,RNG) = sqrt(vx(RNG)**2+cfr_smooth*amp) + max(sqrt( abs(2.0*pmag(RNG) +              p(RNG))/uui(idn,RNG)),small)
#else /* !ISO */
      cfri(1,RNG) = sqrt(vx(RNG)**2+cfr_smooth*amp) + max(sqrt( abs(2.0*pmag(RNG) + flind%ion%gam*p(RNG))/uui(idn,RNG)),small)
#endif /* !ISO */
      !> \deprecated BEWARE: that is the cause of fast decreasing of timestep in galactic disk problem
      !>
      !! \todo  find why is it so
      !! if such a treatment is OK then should be applied also in both cases of neutral and ionized gas
      !!    do i = 2,n-1
      !!       cfri(1,i) = maxval( [c_fr(i-1), c_fr(i), c_fr(i+1)] )
      !!    enddo
      !<

      cfri(1,1) = cfri(1,2)
      cfri(1,n) = cfri(1,n-1)
      cfri = spread(cfri(1,:),1,flind%ion%all)
#endif /* LOCAL_FR_SPEED */

#ifdef GLOBAL_FR_SPEED
      ! The freezing speed is now computed globally
      !  (c=const for the whole domain) in subroutine 'timestep'
      cfri(:,:) = flind%ion%snap%c
#endif /* GLOBAL_FR_SPEED */

   end subroutine flux_ion

end module fluxionized
