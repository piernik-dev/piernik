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
!! \brief Computation of %fluxes for the tracer fluid
!!
!!The flux functions for dust are given by:
!!\f[
!!  \vec{F}{(\vec{u})} =
!!  \left(\begin{array}{c}
!!    \rho v_x \\
!!  \end{array}\right),
!!  \qquad
!!  \vec{G}{(\vec{u})} =
!!  \left(\begin{array}{c}
!!    \rho v_y \\
!!  \end{array}\right),
!!\qquad
!!  \vec{H}{(\vec{u})} =
!!  \left(\begin{array}{c}
!!    \rho v_z \\
!!  \end{array}\right),
!!\f]
!!
!<
!*/
module fluxtracer
! pulled by TRACER
   implicit none
   private
   public :: flux_tracer

contains
!==========================================================================================

   subroutine flux_tracer(fluxt,uut,vx)

      implicit none
      real, dimension(:), intent(out), pointer  :: fluxt   !< flux for tracer
      real, dimension(:), intent(in),  pointer  :: uut     !< part of u for tracer
      real, dimension(:), intent(in),  pointer  :: vx      !< velocity field of fluid for current sweep

      integer :: n, nm

      n  = size(fluxt,1); nm = n-1

      fluxt(RNG) = uut(RNG)*vx(RNG)
      fluxt(1)   = fluxt(2); fluxt(n) = fluxt(nm)
      return
   end subroutine flux_tracer

end module fluxtracer
