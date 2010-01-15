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

!> 
!! \brief (MH/JD) Timestep computation for the neutral fluid (doxy comments ready)
!!
!! %Timestep for the neutral fluid is set as the minimum %timestep for all of the MPI blocks times the Courant number.
!! To compute the %timestep in each MPI block, the fastest speed at which information travels in each direction is computed as
!! \f{equation}
!! c_x=\max\limits_{i,j,k}{\left(v_x^{i,j,k}+c_s^{i,j,k}\right)},
!! \f}
!! where \f$v_x^{i,j,k}\f$ is the maximum speed in \f$x\f$ direction for the cell \f$(i,j,k)\f$ and \f$c_s^{i,j,k}\f$ is the speed of sound for
!! neutral fluid computed as \f$c_s^{i,j,k}=\sqrt{\left|\frac{\gamma p}{\rho^{i,j,k}}\right|}\f$, where \f$p\f$ stands for pressure,
!! \f$\gamma\f$ is adiabatic index for neutral fluid and \f$\rho^{i,j,k}\f$ is the neutral fluid density in the cell
!! \f$(i,j,k)\f$. For directions \f$y, z\f$ the computations are made in similar way.
!!
!! %Timestep for each MPI block is then computed as
!! \f{equation}
!! dt=\min{\left(\left|\frac{dx}{c_x}\right|,\left|\frac{dy}{c_y}\right|,\left|\frac{dz}{c_z}\right|\right)},
!! \f}
!! where \f$dx\f$, \f$dy\f$ and \f$dz\f$ are the cell lengths in each direction.
!!
!! Information about the computed %timesteps is exchanged between MPI blocks in order to choose the minimum %timestep for the fluid.
!! The final %timestep is multiplied by the Courant number specified in parameters of each task.
!<
module timestepneutral

  real :: dt_neu                                                        !< final timestep for neutral fluid
  real :: c_neu                                                         !< maximum speed at which information travels in the neutral fluid

contains

  subroutine timestep_neu
    use mpisetup
    use constants,   only : big    
    use grid, only        : dx,dy,dz,nb,ks,ke,is,ie,js,je,nxd,nyd,nzd
    use arrays, only      : u,b
    use initneutral, only : gamma_neu, cs_iso_neu,cs_iso_neu2   
    use initneutral, only : idnn,imxn,imyn,imzn
#ifndef ISO
    use initneutral, only : ienn
#endif /* ISO */


    implicit none

    real dt_neu_proc                   !< minimum timestep for the neutral fluid for the current processor
    real dt_neu_all                    !< minimum timestep for the neutral fluid for all the processors
    real c_max_all                     !< maximum speed for the neutral fluid for all the processors
    real dt_neu_proc_x                 !< timestep computed for X direction for the current processor
    real dt_neu_proc_y                 !< timestep computed for Y direction for the current processor
    real dt_neu_proc_z                 !< timestep computed for Z direction for the current processor
    real cx                            !< maximum velocity for X direction
    real cy                            !< maximum velocity for Y direction
    real cz                            !< maximum velocity for Z direction
    real vx                            !< velocity in X direction computed for current cell
    real vy                            !< velocity in Y direction computed for current cell
    real vz                            !< velocity in Z direction computed for current cell
    real cs                            !< speed of sound for the neutral fluid


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
            cs = sqrt(cs_iso_neu2)           
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

    if(nxd /= 1 .and. cx /= 0) then
       dt_neu_proc_x = dx/cx
    else 
       dt_neu_proc_x = big 
    endif
    if(nyd /= 1 .and. cy /= 0) then 
       dt_neu_proc_y = dy/cy  
    else 
       dt_neu_proc_y = big 
    endif
    if(nzd /= 1 .and. cz /= 0) then 
       dt_neu_proc_z = dz/cz  
    else 
       dt_neu_proc_z = big 
    endif
    
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

