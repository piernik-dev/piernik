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
!    Initial implemetation of PIERNIK code was based on TVD split MHD code by
!    Ue-Li Pen
!        see: Pen, Arras & Wong (2003) for algorithm and
!             http://www.cita.utoronto.ca/~pen/MHD
!             for original source code "mhd.f90"
!
!    For full list of developers see $PIERNIK_HOME/license/pdt.txt
!
#include "piernik.def"
#define RNG 2:n-1

!>
!! \brief (MH/JD) (doxy comments ready) Computation of %fluxes for the ionized fluid
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

module fluxionized
  implicit none

  contains
!==========================================================================================

  subroutine flux_ion(fluxi,cfri,vx,uui,bb,n)

    use constants,       only : small
    use fluidindex,      only : nmag
    use fluidindex,      only : ibx,iby,ibz
    use fluidindex,      only : idn,imx,imy,imz,ien
    use fluidindex,      only : nvar_ion

    use initionized,     only : gamma_ion, cs_iso_ion2
    use timestepionized, only : c_ion
    use constants,       only : small

    implicit none
    integer :: n                            !< number of cells in the current sweep

! locals
    real, dimension(nvar_ion,n):: fluxi     !< flux of ionized fluid
    real, dimension(nvar_ion,n):: uui       !< part of u for ionized fluid
    real, dimension(nvar_ion,n):: cfri      !< freezing speed for ionized fluid
    real, dimension(nmag,n):: bb            !< magnetic field
    real, dimension(n) :: vx                !< velocity for current sweep
    real, dimension(n) :: ps                !< total pressure of ionized fluid
    real, dimension(n) :: p                 !< thermal pressure of ionized fluid
    real, dimension(n) :: pmag              !< pressure of magnetic field

    fluxi   = 0.0
    cfri    = 0.0
    vx      = 0.0

    pmag(RNG)=0.5*( bb(ibx,RNG)**2 + bb(iby,RNG)**2 +bb(ibz,RNG)**2 )
    vx(RNG)=uui(imx,RNG)/uui(idn,RNG)

#ifdef ISO
    p(RNG) = cs_iso_ion2*uui(idn,RNG)
    ps(RNG)= p(RNG) + pmag(RNG)
#else /* ISO */
    ps(RNG)=(uui(ien,RNG) - &
      0.5*( uui(imx,RNG)**2 + uui(imy,RNG)**2 + uui(imz,RNG)**2 ) &
          / uui(idn,RNG))*(gamma_ion-1.0) + (2.0-gamma_ion)*pmag(RNG)
    p(RNG) = ps(RNG)- pmag(RNG)
#endif /* ISO */

    fluxi(idn,RNG)=uui(imx,RNG)
    fluxi(imx,RNG)=uui(imx,RNG)*vx(RNG)+ps(RNG) - bb(ibx,RNG)**2
    fluxi(imy,RNG)=uui(imy,RNG)*vx(RNG)-bb(iby,RNG)*bb(ibx,RNG)
    fluxi(imz,RNG)=uui(imz,RNG)*vx(RNG)-bb(ibz,RNG)*bb(ibx,RNG)
#ifndef ISO
    fluxi(ien,RNG)=(uui(ien,RNG)+ps(RNG))*vx(RNG)-bb(ibx,RNG)*(bb(ibx,RNG)*uui(imx,RNG) &
                +bb(iby,RNG)*uui(imy,RNG)+bb(ibz,RNG)*uui(imz,RNG))/uui(idn,RNG)
#endif /* ISO */

#ifdef LOCAL_FR_SPEED

!       The freezing speed is now computed locally (in each cell)
!       as in Trac & Pen (2003). This ensures much sharper shocks,
!       but sometimes may lead to numerical instabilities
#ifdef ISO
    cfri(1,RNG) = abs(vx(RNG)) &
               + max(sqrt( abs(2.0*pmag(RNG) + p(RNG))/uui(idn,RNG)),small)
#else /* ISO */
    cfri(1,RNG) = abs(vx(RNG)) &
               + max(sqrt( abs(2.0*pmag(RNG) + gamma_ion*p(RNG) &
                )/uui(idn,RNG)),small)
#endif /* ISO */
    cfri(1,1) = cfri(1,2)
    cfri(1,n) = cfri(1,n-1)
    cfri = spread(cfri(1,:),1,nvar_ion)
#endif /* LOCAL_FR_SPEED */

#ifdef GLOBAL_FR_SPEED
!       The freezing speed is now computed globally
!       (c=const for the whole domain) in sobroutine 'timestep'
    cfri(:,:) = c_ion
#endif /* GLOBAL_FR_SPEED */

  end subroutine flux_ion


end module fluxionized
