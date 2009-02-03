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

module timestepdust

  real :: dt_dst,c_dst

contains

  subroutine timestep_dst
    use mpisetup
    use constants,   only : big    
    use grid, only      : dx,dy,dz,nb,ks,ke,is,ie,js,je,nxd,nyd,nzd
    use arrays, only    : u,b
    use initdust, only  : idnd,imxd,imyd,imzd
    use constants, only : small

    implicit none

    real dt_dst_proc, dt_dst_all, c_max_all
    real dt_dst_proc_x, dt_dst_proc_y, dt_dst_proc_z
    real cx, cy, cz, vx, vy, vz
    

! locals

    real p
    integer i,j,k


    cx    = 0.0
    cy    = 0.0
    cz    = 0.0
    c_dst     = 0.0

    do k=ks,ke
      do j=js,je
        do i=is,ie

          vx=abs(u(imxd,i,j,k)/u(idnd,i,j,k))
          vy=abs(u(imyd,i,j,k)/u(idnd,i,j,k))
          vz=abs(u(imzd,i,j,k)/u(idnd,i,j,k))

          cx=max(cx,vx)
          cy=max(cy,vy)
          cz=max(cz,vz)
          c_dst =max(c_dst,cx,cy,cz)

        end do
      end do
    end do

    
    if(nxd /= 1) then
       dt_dst_proc_x = dx/max(cx,small)  
    else 
       dt_dst_proc_x = big 
    endif
    if(nyd /= 1) then 
       dt_dst_proc_y = dy/max(cy,small) 
    else 
       dt_dst_proc_y = big 
    endif
    if(nzd /= 1) then 
       dt_dst_proc_z = dz/max(cz,small) 
    else 
       dt_dst_proc_z = big 
    endif
    
    
    dt_dst_proc   = min(dt_dst_proc_x, dt_dst_proc_y, dt_dst_proc_z)

    call MPI_REDUCE(c_dst, c_max_all, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierr)
    call MPI_BCAST(c_max_all, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)

    c_dst = c_max_all

    call MPI_REDUCE(dt_dst_proc, dt_dst_all, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, comm, ierr)
    call MPI_BCAST(dt_dst_all, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)
    dt_dst = cfl*dt_dst_all
    
!    write(*,*) 'timestep_dst:', dt_dst

  end subroutine timestep_dst

!-------------------------------------------------------------------------------
end module timestepdust

