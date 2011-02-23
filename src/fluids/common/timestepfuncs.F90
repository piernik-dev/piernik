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

   subroutine compute_c_max(fl,cs,i,j,k,cx,cy,cz,c_max)

      use arrays,      only: u
      use types,       only: component_fluid

      implicit none

      type(component_fluid), pointer, intent(in) :: fl
      real, intent(in)                           :: cs
      integer, intent(in)                        :: i, j, k
      real, intent(inout)                        :: cx, cy, cz, c_max
      real                                       :: vx, vy, vz

      if ( u(fl%idn,i,j,k) > 0.0) then
         vx = abs(u(fl%imx,i,j,k)/u(fl%idn,i,j,k))
         vy = abs(u(fl%imy,i,j,k)/u(fl%idn,i,j,k))
         vz = abs(u(fl%imz,i,j,k)/u(fl%idn,i,j,k))
      else
         vx = 0.0; vy = 0.0; vz = 0.0
      endif

      cx    = max(cx,vx+cs)
      cy    = max(cy,vy+cs)
      cz    = max(cz,vz+cs)
      c_max = max(c_max,cx,cy,cz)

   end subroutine compute_c_max

   subroutine compute_dt(fl,cx,cy,cz,c_max,c_out,dt_out)

      use types,     only: component_fluid
      use grid,      only: cg
      use constants, only: big, xdim, ydim, zdim, GEO_RPZ
      use mpi,       only: MPI_DOUBLE_PRECISION, MPI_MIN, MPI_MAX
      use mpisetup,  only: comm, ierr, cfl, has_dir, dom, geometry_type

      implicit none

      type(component_fluid), pointer, intent(inout) :: fl
      real, intent(in)  :: cx, cy, cz, c_max
      real, intent(out) :: c_out, dt_out
      real :: dt_proc             !< minimum timestep for the current processor
      real :: dt_all              !< minimum timestep for all the processors
      real :: c_max_all           !< maximum speed for the fluid for all the processors
      real :: dt_proc_x           !< timestep computed for X direction for the current processor
      real :: dt_proc_y           !< timestep computed for Y direction for the current processor
      real :: dt_proc_z           !< timestep computed for Z direction for the current processor

      if (has_dir(xdim) .and. cx /= 0) then
         dt_proc_x = cg%dx/cx
      else
         dt_proc_x = big
      endif
      if (has_dir(ydim) .and. cy /= 0) then
         dt_proc_y = cg%dy/cy
         if (geometry_type == GEO_RPZ) dt_proc_y = dt_proc_y * dom%xmin
      else
         dt_proc_y = big
      endif
      if (has_dir(zdim) .and. cz /= 0) then
         dt_proc_z = cg%dz/cz
      else
         dt_proc_z = big
      endif

      dt_proc   = min(dt_proc_x, dt_proc_y, dt_proc_z)

      call MPI_Reduce(c_max, c_max_all, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierr)
      call MPI_Bcast(c_max_all, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      call MPI_Reduce(dt_proc, dt_all, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, comm, ierr)
      call MPI_Bcast(dt_all, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      c_out  = c_max_all
      dt_out = cfl*dt_all

      fl%snap%c  = c_max_all
      fl%snap%dt = dt_out

   end subroutine compute_dt

end module timestepfuncs
