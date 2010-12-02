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
!! \brief (MH/JD) (doxy comments ready) Computation of %fluxes for the neutral fluid
!!
!!The flux functions for neutral fluid are given by
!!
!!\f[
!!  \vec{F}{(\vec{u})} =
!!  \left(\begin{array}{c}
!!    \rho v_x \\
!!    \rho v_x^2+p \\
!!    \rho v_x v_y\\
!!    \rho v_x v_z\\
!!    (e + p)v_x
!!  \end{array}\right),
!!  \qquad
!!  \vec{G}{(\vec{u})} =
!!  \left(\begin{array}{c}
!!    \rho v_y \\
!!    \rho v_y v_x\\
!!    \rho v_y^2+p\\
!!    \rho v_y v_z\\
!!    (e + p)v_y
!!  \end{array}\right),
!!\qquad
!!  \vec{H}{(\vec{u})} =
!!  \left(\begin{array}{c}
!!    \rho v_z \\
!!    \rho v_z v_x\\
!!    \rho v_z v_y \\
!!    \rho v_z^2+p \\
!!    (e + p)v_z
!!  \end{array}\right),
!!\f]
!<
!*/
module fluxneutral
! pulled by NEUTRAL
   implicit none
   private
   public :: flux_neu

contains
!==========================================================================================

   subroutine flux_neu(fluxn,cfrn,uun,n,vx,bb,cs_iso2)

      use constants,       only: small
      use fluidindex,      only: idn, imx, imy, imz, ien, nvar
      use mpisetup,        only: cfr_smooth, smallp
      use timestep,        only: c_all

      implicit none
      integer, intent(in)                        :: n         !< number of cells in the current sweep
      real, dimension(:,:), intent(out), pointer :: fluxn     !< flux of neutral fluid
      real, dimension(:,:), intent(in),  pointer :: uun       !< part of u for neutral fluid
      real, dimension(:,:), intent(out), pointer :: cfrn      !< freezing speed for neutral fluid
      real, dimension(:,:), intent(in),  pointer :: bb        !< \copydoc fluxes::interface::flux_interface::bb
      real, dimension(:),   intent(out), pointer :: vx        !< velocity of neutral fluid for current sweep
      real, dimension(:),   intent(in),  pointer :: cs_iso2   !< \copydoc fluxes::interface::flux_interface::cs_iso2

      ! locals
      real, dimension(n) :: p         !< pressure of neutral fluid
#ifdef LOCAL_FR_SPEED
      real               :: minvx     !<
      real               :: maxvx     !<
      real               :: amp       !<
#endif /* LOCAL_FR_SPEED */

      fluxn   = 0.0; cfrn    = 0.0; vx      = 0.0; p       = 0.0

      vx(RNG) = uun(imx,RNG)/uun(idn,RNG)
#ifdef ISO
      p(RNG)  = nvar%neu%cs2*uun(idn,RNG)
#else /* !ISO */
      p(RNG)  = (uun(ien,RNG)  &
      - 0.5*( uun(imx,RNG)**2 + uun(imy,RNG)**2 + uun(imz,RNG)**2 ) &
          / uun(idn,RNG))*(nvar%neu%gam_1)
      p(RNG)  = max(p(RNG),smallp)
#endif /* !ISO */

      fluxn(idn,RNG)=uun(imx,RNG)
      fluxn(imx,RNG)=uun(imx,RNG)*vx(RNG)+p(RNG)
      fluxn(imy,RNG)=uun(imy,RNG)*vx(RNG)
      fluxn(imz,RNG)=uun(imz,RNG)*vx(RNG)
#ifndef ISO
      fluxn(ien,RNG)=(uun(ien,RNG)+p(RNG))*vx(RNG)
#endif /* !ISO */

#ifdef LOCAL_FR_SPEED

      !       The freezing speed is now computed locally (in each cell)
      !       as in Trac & Pen (2003). This ensures much sharper shocks,
      !       but sometimes may lead to numerical instabilities
      minvx = minval(vx(RNG))
      maxvx = maxval(vx(RNG))
      amp   = 0.5*(maxvx-minvx)
      !    c_fr  = 0.0
#ifdef ISO
      cfrn(1,RNG) = sqrt(vx(RNG)**2+cfr_smooth*amp) + max(sqrt( abs(             p(RNG))/uun(idn,RNG)),small)
#else /* !ISO */
      cfrn(1,RNG) = sqrt(vx(RNG)**2+cfr_smooth*amp) + max(sqrt( abs(nvar%neu%gam*p(RNG))/uun(idn,RNG)),small)
#endif /* !ISO */
      !BEWARE: that is the cause of fast decreasing of timestep in galactic disk problem
      !TODO: find why is it so
      !if such a treatment is OK then should be applied also in both cases of neutral and ionized gas
      !    do i = 2,n-1
      !       cfrn(1,i) = maxval( [c_fr(i-1), c_fr(i), c_fr(i+1)] )
      !    enddo

      cfrn(1,1) = cfrn(1,2)
      cfrn(1,n) = cfrn(1,n-1)
      cfrn = spread(cfrn(1,:),1,nvar%neu%all)
#endif /* LOCAL_FR_SPEED */

#ifdef GLOBAL_FR_SPEED
      !       The freezing speed is now computed globally
      !       (c=const for the whole domain) in sobroutine 'timestep'

      !    cfrn(:,:) = nvar%neu%snap%c   ! check which c_xxx is better
      cfrn(:,:) = c_all
#endif /* GLOBAL_FR_SPEED */

   end subroutine flux_neu

end module fluxneutral
