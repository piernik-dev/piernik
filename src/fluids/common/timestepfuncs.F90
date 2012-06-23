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

!>
!! \brief (KK)
!<
module timestepfuncs

! pulled by ANY

   implicit none

   private
   public :: timestep_fluid

contains

   subroutine timestep_fluid(cg, fl, dt, c_fl)

      use constants,     only: big, xdim, ydim, zdim, ndims, GEO_RPZ, LO
      use domain,        only: dom
      use fluidtypes,    only: component_fluid
      use global,        only: cfl
      use grid_cont,     only: grid_container

      implicit none

      type(grid_container), pointer, intent(in) :: cg !< current grid container
      type(component_fluid), pointer, intent(in) :: fl
      real, intent(out) :: dt !< resulting timestep
      real, intent(out) :: c_fl !< maximum speed at which information travels in fluid

      ! locals
      real, dimension(ndims) :: c !< maximum velocity for all directions
      real, dimension(ndims) :: dt_dim !< timestep for the current cg
      real, dimension(ndims) :: v !< velocity
      integer :: i, j, k, d

      c_fl = 0.0; c(:) = 0.0

      do k= cg%ks, cg%ke
         do j = cg%js, cg%je
            do i = cg%is, cg%ie

               if (cg%u(fl%idn,i,j,k) > 0.0) then
                  v = abs(cg%u(fl%imx:fl%imz, i, j, k) / cg%u(fl%idn, i, j, k))
               else
                  v = 0.0
               endif
               c = max(c, v + fl%get_cs(cg, fl, i, j, k))
               c_fl = max(c_fl, maxval(c))

            enddo
         enddo
      enddo

      do d = xdim, zdim
         if (dom%has_dir(d) .and. c(d) > 0.0) then
            dt_dim(d) = cg%dl(d)/c(d)
            if (dom%geometry_type == GEO_RPZ .and. d == ydim) dt_dim(d) = dt_dim(d) * dom%edge(xdim, LO)
         else
            dt_dim(d) = big
         endif
      enddo

      dt = cfl * minval(dt_dim)
      fl%c  = c_fl

   end subroutine timestep_fluid

end module timestepfuncs
