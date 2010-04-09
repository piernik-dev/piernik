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
!! \brief (MH/JD) (doxy comments ready) Computation of %fluxes for the dust fluid
!!
!!The flux functions for dust are given by:
!!\f[
!!  \vec{F}{(\vec{u})} =
!!  \left(\begin{array}{c}
!!    \rho v_x \\
!!    (\rho v_x) v_x\\
!!    (\rho v_y) v_x\\
!!    (\rho v_z) v_x
!!  \end{array}\right),
!!  \qquad
!!  \vec{G}{(\vec{u})} =
!!  \left(\begin{array}{c}
!!    \rho v_y \\
!!    (\rho v_x) v_y \\
!!    (\rho v_y) v_y \\
!!    (\rho v_z) v_y
!!  \end{array}\right),
!!\qquad
!!  \vec{H}{(\vec{u})} =
!!  \left(\begin{array}{c}
!!    \rho v_z \\
!!    (\rho v_x) v_z\\
!!    (\rho v_y) v_z \\
!!    (\rho v_z) v_z
!!  \end{array}\right),
!!\f]
!!
!<


module fluxdust
  implicit none

  contains
!==========================================================================================

  subroutine flux_dst(fluxd,cfrd,uud,n)

    use mpisetup,        only : cfr_smooth
    use constants,       only : small
    use fluidindex,      only : idn,imx,imy,imz,ien
    use fluidindex,      only : nvar_dst

!    use timestepdust, only : c_dst
    use timestep, only : c_all


    implicit none
    integer n                                !< number of cells in the current sweep

! locals
    real :: minvx, maxvx, amp
    real, dimension(nvar_dst,n):: fluxd      !< flux for dust
    real, dimension(nvar_dst,n):: uud        !< part of u for dust
    real, dimension(nvar_dst,n):: cfrd       !< freezing speed for dust
    real, dimension(n) :: p                  !< pressure
    real, dimension(n) :: vx                 !< velocity for current sweep

    fluxd   = 0.0
    cfrd    = 0.0
    vx      = 0.0

    where(uud(idn,RNG) > 0.0)
       vx(RNG)=uud(imx,RNG)/uud(idn,RNG)
    elsewhere
       vx(RNG) = 0.0
    endwhere

    where(abs(vx) < 1.e-20 .and. vx /= 0.0) vx = sign(1.e-20,vx)

    fluxd(idn,RNG)=uud(imx,RNG)
    fluxd(imx,RNG)=uud(imx,RNG)*vx(RNG)
    fluxd(imy,RNG)=uud(imy,RNG)*vx(RNG)
    fluxd(imz,RNG)=uud(imz,RNG)*vx(RNG)

#ifdef LOCAL_FR_SPEED
!       The freezing speed is now computed locally (in each cell)
!       as in Trac & Pen (2003). This ensures much sharper shocks,
!       but sometimes may lead to numerical instabilities

    minvx = minval(vx(RNG))
    maxvx = maxval(vx(RNG))
    amp   = (maxvx-minvx)*0.5
    cfrd(1,RNG) = max(sqrt(vx(RNG)**2+cfr_smooth*amp),small)

    cfrd(1,1) = cfrd(1,2)
    cfrd(1,n) = cfrd(1,n-1)
    cfrd = spread(cfrd(1,:),1,nvar_dst)
#endif /* LOCAL_FR_SPEED */

#ifdef GLOBAL_FR_SPEED
!       The freezing speed is now computed globally
!       (c=const for the whole domain) in sobroutine 'timestep'
    cfrd(:,:) = c_dst
#endif /* GLOBAL_FR_SPEED */

  end subroutine flux_dst


end module fluxdust
