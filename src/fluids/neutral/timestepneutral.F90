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

module timestepneutral

  real :: dt_neu,c_neu

contains

  subroutine timestep_neu
    use mpisetup
    use grid, only     : dx,dy,dz,nb,ks,ke,is,ie,js,je
    use arrays, only   : u,b
    use start, only  : cfl
    use initneutral, only : gamma_neu, cs_iso_neu,cs_iso_neu2   
    use initneutral, only : idnn,imxn,imyn,imzn
#ifndef ISO
    use initneutral, only : ienn
#endif /* ISO */


    implicit none

    real dt_neu_proc, dt_neu_all, c_max_all
    real dt_neu_proc_x, dt_neu_proc_y, dt_neu_proc_z
    real cx, cy, cz, vx, vy, vz, cs
    

! locals

    real p
    integer i,j,k


    cx    = 0.0
    cy    = 0.0
    cz    = 0.0
    c_neu = 0.0

    do k=ks,ke
      do j=js,je
        do i=is,ie

          vx=abs(u(imxn,i,j,k)/u(idnn,i,j,k))
          vy=abs(u(imyn,i,j,k)/u(idnn,i,j,k))
          vz=abs(u(imzn,i,j,k)/u(idnn,i,j,k))

#ifdef ISO
            p = cs_iso_neu2*u(idnn,i,j,k)
            cs = sqrt(abs( p/u(idnn,i,j,k)) )
#else /* ISO */
            p=(u(ienn,i,j,k)-sum(u(imxn:imzn,i,j,k)**2,1) &
              /u(idnn,i,j,k)/2.)*(gamma_neu-1.)

            cs = sqrt(abs(  (gamma_neu*p)/u(idnn,i,j,k)) )
#endif /* ISO */

          cx=max(cx,vx+cs)
          cy=max(cy,vy+cs)
          cz=max(cz,vz+cs)
          c_neu =max(c_neu,cx,cy,cz)

        end do
      end do
    end do

    dt_neu_proc_x = dx/cx
    dt_neu_proc_y = dy/cy
    dt_neu_proc_z = dz/cz
    dt_neu_proc   = min(dt_neu_proc_x, dt_neu_proc_y, dt_neu_proc_z)

    call MPI_REDUCE(c_neu, c_max_all, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierr)
    call MPI_BCAST(c_max_all, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)

    c_neu = c_max_all

    call MPI_REDUCE(dt_neu_proc, dt_neu_all, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, comm, ierr)
    call MPI_BCAST(dt_neu_all, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)
    dt_neu = cfl*dt_neu_all
    
!    write(*,*) 'timestep_neu:', dt_neu

  end subroutine timestep_neu

!-------------------------------------------------------------------------------
end module timestepneutral

