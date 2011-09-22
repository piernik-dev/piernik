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
   public :: compute_c_max, compute_dt

contains

   subroutine compute_c_max(fl, cs, i, j, k, c, c_max, cg)

      use constants,  only: ndims
      use fluidtypes, only: component_fluid
      use grid_cont,  only: grid_container

      implicit none

      type(component_fluid), pointer, intent(in) :: fl
      real, intent(in)                           :: cs
      integer, intent(in)                        :: i, j, k
      real, dimension(ndims), intent(inout)      :: c
      real, intent(inout)                        :: c_max
      type(grid_container), pointer, intent(in)  :: cg

      real, dimension(ndims)                     :: v

      if ( cg%u%arr(fl%idn,i,j,k) > 0.0) then
         v(:) = abs(cg%u%arr(fl%imx:fl%imz, i, j, k)/cg%u%arr(fl%idn, i, j, k))
      else
         v(:) = 0.0
      endif

      c(:) = max( c(:), v(:) + cs )
      c_max = max(c_max, maxval(c(:)))

   end subroutine compute_c_max

   subroutine compute_dt(c, dt_out, cg)

      use constants,  only: big, xdim, ydim, zdim, ndims, GEO_RPZ, LO
      use domain,     only: dom, geometry_type, has_dir
      use global,     only: cfl
      use grid_cont,  only: grid_container

      implicit none

      real, dimension(ndims), intent(in)  :: c
      real, intent(out) :: dt_out
      type(grid_container), pointer, intent(in) :: cg

      real, dimension(ndims) :: dt_proc !< timestep for the current cg
      integer :: d

      do d = xdim, zdim
         if (has_dir(d) .and. c(d) /= 0.) then
            dt_proc(d) = cg%dl(d)/c(d)
            if (geometry_type == GEO_RPZ .and. d == ydim) dt_proc(d) = dt_proc(d) * dom%edge(xdim, LO)
         else
            dt_proc(d) = big
         endif
      enddo

      dt_out = cfl * minval(dt_proc)

   end subroutine compute_dt

end module timestepfuncs
