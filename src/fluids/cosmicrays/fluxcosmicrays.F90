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
!! \brief (MH) Computation of advection fluxes of Cosmic Rays
!!
!!
!<

module fluxcosmicrays
  implicit none

  contains 
!==========================================================================================

  subroutine flux_crs(fluxc,vion,uuc,n)
  
    use constants,       only : small
    use fluidindex,      only : nvar_crs,icr   
  
    implicit none
    integer n

! locals
    real, dimension(nvar_crs,n):: fluxc,uuc
    real, dimension(n) :: vion 
    
    fluxc   = 0.0

    fluxc(icr,RNG)= uuc(icr,RNG)*vion(RNG)

  end subroutine flux_crs


end module fluxcosmicrays
