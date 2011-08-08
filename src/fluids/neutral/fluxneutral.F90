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
!! \brief (MH/JD) [R] Computation of %fluxes for the neutral fluid
!!
!!The flux functions for neutral fluid are given by
!!
!!\f[
!!  \vec{F}{(\vec{u})} =
!!  \left(\begin{array}{c}
!!    \rho v_x \\
!!    \rho v_x^2 + p \\
!!    \rho v_x v_y\\
!!    \rho v_x v_z\\
!!    (e + p)v_x
!!  \end{array}\right),
!!  \qquad
!!  \vec{G}{(\vec{u})} =
!!  \left(\begin{array}{c}
!!    \rho v_y \\
!!    \rho v_y v_x\\
!!    \rho v_y^2 + p\\
!!    \rho v_y v_z\\
!!    (e + p)v_y
!!  \end{array}\right),
!!\qquad
!!  \vec{H}{(\vec{u})} =
!!  \left(\begin{array}{c}
!!    \rho v_z \\
!!    \rho v_z v_x\\
!!    \rho v_z v_y \\
!!    \rho v_z^2 + p \\
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
!
! OPT: This routine may cost as much as 30% of rtvd. It seems that all the data fit well a 512kB L2 cache, but Ir:Dr:Dw is like 8:2:1
! OPT: \todo Try an explicit loop over RNG to check if we're better than the compiler
! OPT: similar treatment may be helpful for fluxionized.F90, fluxdust.F90 and fluxcosmicrays.F90
!
   subroutine flux_neu(fluxn, cfrn, uun, n, vx, p, bb, cs_iso2)

      use constants,  only: small
      use func,       only: ekin
      use fluidindex, only: idn, imx, imy, imz, ien, flind
      use global,     only: cfr_smooth, smallp
#ifdef GLOBAL_FR_SPEED
      use timestep,   only: c_all
#endif /* GLOBAL_FR_SPEED */

      implicit none

      integer(kind=4), intent(in)                  :: n         !< number of cells in the current sweep
      real, dimension(:,:), intent(inout), pointer :: fluxn     !< flux of neutral fluid
      real, dimension(:,:), intent(inout), pointer :: cfrn      !< freezing speed for neutral fluid
      real, dimension(:,:), intent(in),    pointer :: uun       !< part of u for neutral fluid
      real, dimension(:),   intent(inout), pointer :: vx        !< velocity of neutral fluid for current sweep
      real, dimension(:),   intent(inout), pointer :: p         !< pressure of neutral fluid for current sweep
      real, dimension(:,:), intent(in),    pointer :: bb        !< magnetic field x,y,z-components table
      real, dimension(:),   intent(in),    pointer :: cs_iso2   !< isothermal sound speed squared

      ! locals
      integer            :: nm
#ifdef LOCAL_FR_SPEED
      integer            :: i
      real               :: minvx     !<
      real               :: maxvx     !<
      real               :: amp       !<
#endif /* LOCAL_FR_SPEED */

      nm = n-1
      vx(RNG) = uun(imx,RNG)/uun(idn,RNG) ; vx(1) = vx(2); vx(n) = vx(nm)
#ifdef ISO
      p(RNG)  = cs_iso2(RNG)*uun(idn,RNG) ; p(1) = p(2); p(n) = p(nm)
#else /* !ISO */
      p(RNG)  = (uun(ien,RNG) - ekin(uun(imx,RNG),uun(imy,RNG),uun(imz,RNG),uun(idn,RNG)) )*(flind%neu%gam_1)
      p(RNG)  = max(p(RNG),smallp)
#endif /* !ISO */

      fluxn(idn,RNG)=uun(imx,RNG)
      fluxn(imx,RNG)=uun(imx,RNG)*vx(RNG)+p(RNG)
      fluxn(imy,RNG)=uun(imy,RNG)*vx(RNG)
      fluxn(imz,RNG)=uun(imz,RNG)*vx(RNG)
#ifndef ISO
      fluxn(ien,RNG)=(uun(ien,RNG)+p(RNG))*vx(RNG)
#endif /* !ISO */
      fluxn(:,1) = fluxn(:,2) ; fluxn(:,n) = fluxn(:,nm)

#ifdef LOCAL_FR_SPEED

      ! The freezing speed is now computed locally (in each cell)
      !  as in Trac & Pen (2003). This ensures much sharper shocks,
      !  but sometimes may lead to numerical instabilities
      minvx = minval(vx(RNG))
      maxvx = maxval(vx(RNG))
      amp   = 0.5*(maxvx-minvx)
      !    c_fr  = 0.0
#ifdef ISO
      cfrn(1,RNG) = sqrt(vx(RNG)**2+cfr_smooth*amp) + max(sqrt( abs(              p(RNG))/uun(idn,RNG)),small)
#else /* !ISO */
      cfrn(1,RNG) = sqrt(vx(RNG)**2+cfr_smooth*amp) + max(sqrt( abs(flind%neu%gam*p(RNG))/uun(idn,RNG)),small)
#endif /* !ISO */
      !> \deprecated BEWARE: that is the cause of fast decreasing of timestep in galactic disk problem
      !>
      !! \todo find why is it so
      !! if such a treatment is OK then should be applied also in both cases of neutral and ionized gas
      !!    do i = 2,nm
      !!       cfrn(1,i) = maxval( [c_fr(i-1), c_fr(i), c_fr(i+1)] )
      !!    enddo
      !<

      cfrn(1,1) = cfrn(1,2);  cfrn(1,n) = cfrn(1,nm)
      do i = 2, flind%neu%all
         cfrn(i,:) = cfrn(1,:)
      enddo
#endif /* LOCAL_FR_SPEED */

#ifdef GLOBAL_FR_SPEED
      ! The freezing speed is now computed globally
      !  (c=const for the whole domain) in subroutine 'timestep'

      !    cfrn(:,:) = flind%neu%snap%c   ! check which c_xxx is better
      cfrn(:,:) = c_all
#endif /* GLOBAL_FR_SPEED */
      return
      if (.false.) write(0,*) bb, cs_iso2

   end subroutine flux_neu

end module fluxneutral
