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

module timestepionized

  real :: dt_ion,c_ion

 contains 

  subroutine timestep_ion
    use mpisetup
    use grid, only   : dx,dy,dz,nb,ks,ke,is,ie,js,je
    use arrays, only : u,b
    use start, only  : cfl
    use initionized, only : gamma_ion, cs_iso_ion2   
    use initionized, only : idni,imxi,imyi,imzi
#ifndef ISO
    use initionized, only : ieni
#endif /* ISO */

    implicit none

    real dt_ion_proc, dt_ion_all, c_max_all
    real dt_ion_proc_x, dt_ion_proc_y, dt_ion_proc_z
    real cx, cy, cz, vx, vy, vz, cf
    

! locals
    real pmag
    real ps,p
    integer i,j,k


    cx    = 0.0
    cy    = 0.0
    cz    = 0.0
    c_ion     = 0.0

    do k=ks,ke
      do j=js,je
        do i=is,ie

          vx=abs(u(imxi,i,j,k)/u(idni,i,j,k))
          vy=abs(u(imyi,i,j,k)/u(idni,i,j,k))
          vz=abs(u(imzi,i,j,k)/u(idni,i,j,k))

          pmag = sum(b(:,i,j,k)**2,1)/2.

#ifdef ISO
            p = cs_iso_ion2*u(idni,i,j,k)
            ps =p+pmag
            cf = sqrt(abs(  (2.*pmag+p)/u(idni,i,j,k)) )
#else /* ISO */
            ps=(u(ieni,i,j,k)-sum(u(imxi:imzi,i,j,k)**2,1) &
             /u(idni,i,j,k)/2.)*(gamma_ion-1.)+(2.-gamma_ion)*pmag
            p=ps-pmag
            cf = sqrt(abs(  (2.*pmag+gamma_ion*p)/u(idni,i,j,k)) )
#endif /* ISO */

          cx=max(cx,vx+cf)
          cy=max(cy,vy+cf)
          cz=max(cz,vz+cf)
          c_ion =max(c_ion,cx,cy,cz)

        end do
      end do
    end do

    dt_ion_proc_x = dx/cx
    dt_ion_proc_y = dy/cy
    dt_ion_proc_z = dz/cz
    dt_ion_proc   = min(dt_ion_proc_x, dt_ion_proc_y, dt_ion_proc_z)

    call MPI_REDUCE(c_ion, c_max_all, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierr)
    call MPI_BCAST(c_max_all, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)

    c_ion = c_max_all

    call MPI_REDUCE(dt_ion_proc, dt_ion_all, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, comm, ierr)
    call MPI_BCAST(dt_ion_all, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)
    dt_ion = cfl*dt_ion_all
    
!    write(*,*) 'timestep_ion:', dt_ion

  end subroutine timestep_ion

!-------------------------------------------------------------------------------
end module timestepionized

