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
!! \brief (MH) Timestep computation for the dust fluid
!<
module timestepdust
! pulled by ANY
   implicit none

   private
   public :: timestep_dst

contains

   subroutine timestep_dst(cg, dt, c_dst)

      use constants,     only: ndims
      use fluidtypes,    only: component_fluid
      use grid_cont,     only: grid_container
      use fluidindex,    only: flind
      use timestepfuncs, only: compute_c_max, compute_dt

      implicit none

      type(grid_container), pointer, intent(in) :: cg !< current grid container
      real, intent(out)                         :: dt !< resulting timestep
      real, intent(out)                         :: c_dst !< maximum speed at which information travels in dust

      real, dimension(ndims) :: c !< maximum velocity for all directions

! locals

      integer                        :: i, j, k
      class(component_fluid), pointer :: fl

      c_dst = 0.0; c(:) = 0.0
      fl => flind%dst

      do k= cg%ks, cg%ke
         do j = cg%js, cg%je
            do i = cg%is, cg%ie
               call compute_c_max(fl, 0.0, i, j, k, c(:), c_dst, cg)
            enddo
         enddo
      enddo
      call compute_dt(c(:), dt, cg)
      fl%c  = c_dst

   end subroutine timestep_dst

end module timestepdust
