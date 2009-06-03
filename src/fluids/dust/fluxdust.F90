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

module fluxdust
  implicit none

  contains
!==========================================================================================

  subroutine flux_dst(fluxd,cfrd,uud,n)
   
    use mpisetup,        only : cfr_smooth 
    use constants,       only : small
    use fluidindex,      only : idn,imx,imy,imz,ien
    use fluidindex,      only : nvar_dst

    use timestepdust, only : c_dst

    implicit none
    integer n
    
! locals
    real :: minvx, maxvx, amp
    real, dimension(nvar_dst,n):: fluxd,uud,cfrd
    real, dimension(n) :: vx,p  
    
    fluxd   = 0.0
    cfrd    = 0.0
    vx      = 0.0

    where(uud(idn,RNG) > 0.0) 
       vx(RNG)=uud(imx,RNG)/uud(idn,RNG)
    elsewhere
       vx(RNG) = 0.0
    endwhere

    fluxd(idn,RNG)=uud(imx,RNG)
    fluxd(imx,RNG)=uud(imx,RNG)*vx(RNG)
    fluxd(imy,RNG)=uud(imy,RNG)*vx(RNG)
    fluxd(imz,RNG)=uud(imz,RNG)*vx(RNG)
    
#ifndef ISO
    fluxd(ien,RNG)=(uud(ien,RNG))*vx(RNG) 
#endif /* ISO */

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
