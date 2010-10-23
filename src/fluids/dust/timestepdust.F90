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
#include "piernik.def"

!>
!! \brief (MH) Timestep computation for the dust fluid
!!
!!
!<
module timestepdust

   real :: dt_dst                !< final timestep for dust
   real :: c_dst                 !< maximum speed at which information travels in dust

contains

   subroutine timestep_dst
      use types,       only: component_fluid
      use arrays,      only: u, b
      use constants,   only: big, small
      use grid,        only: dx, dy, dz, nb, ks, ke, is, ie, js, je, nxd, nyd, nzd
      use mpisetup,    only: ierr, comm, cfl, MPI_MAX, MPI_MIN, MPI_DOUBLE_PRECISION
      use fluidindex,  only: nvar

      implicit none

      real :: dt_proc = 0.0       !< minimum timestep for the current processor
      real :: dt_all = 0.0        !< minimum timestep for all the processors
      real :: c_max_all = 0.0     !< maximum speed for the fluid for all the processors
      real :: dt_proc_x = 0.0     !< timestep computed for X direction for the current processor
      real :: dt_proc_y = 0.0     !< timestep computed for Y direction for the current processor
      real :: dt_proc_z = 0.0     !< timestep computed for Z direction for the current processor
      real :: cx = 0.0            !< maximum velocity for X direction
      real :: cy = 0.0            !< maximum velocity for Y direction
      real :: cz = 0.0            !< maximum velocity for Z direction
      real :: vx = 0.0            !< velocity in X direction computed for current cell
      real :: vy = 0.0            !< velocity in Y direction computed for current cell
      real :: vz = 0.0            !< velocity in Z direction computed for current cell
      real :: cs = 0.0            !< speed of sound

! locals

      real                           :: c_max = 0
      integer                        :: i, j, k
      type(component_fluid), pointer :: fl

      fl => nvar%dst

      do k=ks,ke
         do j=js,je
            do i=is,ie
               if (u(fl%idn,i,j,k) > 0.0) then
                  vx = abs(u(fl%imx,i,j,k)/u(fl%idn,i,j,k))
                  vy = abs(u(fl%imy,i,j,k)/u(fl%idn,i,j,k))
                  vz = abs(u(fl%imz,i,j,k)/u(fl%idn,i,j,k))
               else
                  vx = 0.0
                  vy = 0.0
                  vz = 0.0
               endif

               cx    = max(cx,vx+cs)
               cy    = max(cy,vy+cs)
               cz    = max(cz,vz+cs)
               c_max = max(c_max,cx,cy,cz)

            enddo
         enddo
      enddo

      if (nxd /= 1 .and. cx /= 0) then
         dt_proc_x = dx/cx
      else
         dt_proc_x = big
      endif
      if (nyd /= 1 .and. cy /= 0) then
         dt_proc_y = dy/cy
      else
         dt_proc_y = big
      endif
      if (nzd /= 1 .and. cz /= 0) then
         dt_proc_z = dz/cz
      else
         dt_proc_z = big
      endif

      dt_proc   = min(dt_proc_x, dt_proc_y, dt_proc_z)

      call MPI_Reduce(c_max, c_max_all, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierr)
      call MPI_Bcast(c_max_all, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      call MPI_Reduce(dt_proc, dt_all, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, comm, ierr)
      call MPI_Bcast(dt_all, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)
    
      c_dst  = c_max_all
      dt_dst = cfl*dt_all

   end subroutine timestep_dst

end module timestepdust
