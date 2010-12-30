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
!! \brief (MH/JD) (doxy comments ready) %Timestep computation for the ionized fluid
!!
!! %Timestep for the ionized fluid is set as the minimum %timestep for all of the MPI blocks times the Courant number.
!! To compute the %timestep in each MPI block, the fastest speed at which information travels in each direction is computed as
!! \f{equation}
!! c_x=\max\limits_{i,j,k}{\left(v_x^{i,j,k}+c_f^{i,j,k}\right)},
!! \f}
!! where \f$v_x^{i,j,k}\f$ is the maximum speed in \f$x\f$ direction for the cell \f$(i,j,k)\f$ and \f$c_f^{i,j,k}\f$ is the speed of sound for
!! ionized fluid computed as \f$c_f^{i,j,k}=\sqrt{\left|\frac{2p_{mag}+\gamma p}{\rho^{i,j,k}}\right|}\f$, where \f$p\f$ stands for pressure,
!! \f$p_{mag}\f$ is pressure of magnetic field, \f$\gamma\f$ is the adiabatic index of the ionized fluid and \f$\rho^{i,j,k}\f$ is fluid density in the cell
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
module timestepionized
! pulled by IONIZED
   implicit none

   private
   public :: dt_ion, c_ion, timestep_ion

   real   :: dt_ion             !< final timestep for ionized fluids
   real   :: c_ion              !< maximum speed at which information travels in the ionized fluid

contains

   subroutine timestep_ion

      use types,         only: component_fluid
      use arrays,        only: u, b
      use grid,          only: D_x, D_y, D_z, cg
      use fluidindex,    only: nvar, ibx, iby, ibz
      use timestepfuncs, only: compute_c_max, compute_dt

      implicit none

      real :: cx                  !< maximum velocity for X direction
      real :: cy                  !< maximum velocity for Y direction
      real :: cz                  !< maximum velocity for Z direction
      real :: cs                  !< speed of sound

! locals
      real                           :: pmag, bx, by, bz, ps, p, c_max
      integer                        :: i, j, k
      type(component_fluid), pointer :: fl

      cx = 0.0; cy = 0.0; cz = 0.0; cs = 0.0;
      pmag = 0.0; bx = 0.0; by = 0.0; bz = 0.0; ps = 0.0; p = 0.0; c_max = 0.0
      fl => nvar%ion

      do k = cg%ks, cg%ke
         do j = cg%js, cg%je
            do i = cg%is, cg%ie

#ifdef MAGNETIC
               bx = (b(ibx,i,j,k) + b(ibx, i+D_x, j,     k    ))/(1.+D_x)
               by = (b(iby,i,j,k) + b(iby, i,     j+D_y, k    ))/(1.+D_y)
               bz = (b(ibz,i,j,k) + b(ibz, i,     j,     k+D_z))/(1.+D_z)

               pmag = 0.5*(bx*bx + by*by + bz*bz)
#else /* !MAGNETIC */
               ! all_mag_boundaries has not been called so we cannot trust b(ibx, cg%ie+D_x:), b(iby,:cg%je+D_y and b(ibz,:,:, cg%ke+D_z
               pmag = 0.
#endif /* !MAGNETIC */

#ifdef ISO
               p  = fl%cs2*u(fl%idn,i,j,k)
               ps = p + pmag
               cs = sqrt(abs(  (2.*pmag+p)/u(fl%idn,i,j,k)) )
#else /* !ISO */
               ps = (u(fl%ien,i,j,k)-sum(u(fl%imx:fl%imz,i,j,k)**2,1) &
                     /u(fl%idn,i,j,k)*0.5)*(fl%gam_1)+(2.-fl%gam)*pmag
               p  = ps - pmag
               cs = sqrt(abs(  (2.*pmag+fl%gam*p)/u(fl%idn,i,j,k)) )
#endif /* !ISO */
               call compute_c_max(fl,cs,i,j,k,cx,cy,cz,c_max)
            enddo
         enddo
      enddo
      call compute_dt(fl,cx,cy,cz,c_max,c_ion,dt_ion)

   end subroutine timestep_ion

end module timestepionized
