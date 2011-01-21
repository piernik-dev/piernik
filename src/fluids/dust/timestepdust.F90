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
!!
!!
!<
module timestepdust
! pulled by DUST
   implicit none
   private
   public :: dt_dst, c_dst, timestep_dst
   real   :: dt_dst                !< final timestep for dust
   real   :: c_dst                 !< maximum speed at which information travels in dust

contains

   subroutine timestep_dst
      use types,         only: component_fluid
      use arrays,        only: u
      use grid,          only: cg
      use fluidindex,    only: flind
      use timestepfuncs, only: compute_c_max, compute_dt

      implicit none

      real :: cx                  !< maximum velocity for X direction
      real :: cy                  !< maximum velocity for Y direction
      real :: cz                  !< maximum velocity for Z direction

! locals

      real                           :: c_max
      integer                        :: i, j, k
      type(component_fluid), pointer :: fl

      c_max = 0.0; cx = 0.0; cy = 0.0; cz = 0.0
      fl => flind%dst

      do k= cg%ks, cg%ke
         do j = cg%js, cg%je
            do i = cg%is, cg%ie
               call compute_c_max(fl,0.0,i,j,k,cx,cy,cz,c_max)
            enddo
         enddo
      enddo
      call compute_dt(fl,cx,cy,cz,c_max,c_dst,dt_dst)

   end subroutine timestep_dst

end module timestepdust
