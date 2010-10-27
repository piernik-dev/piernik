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
!! \brief (MH/JD) (doxy comments ready) %Timestep computation for the ionized fluid
!!
!! %Timestep for the ionized fluid is set as the minimum %timestep for all of the MPI blocks times the Courant number.
!! To compute the %timestep in each MPI block, the fastest speed at which information travels in each direction is computed as
!! \f{equation}
!! c_x=\max\limits_{i,j,k}{\left(v_x^{i,j,k}+c_f^{i,j,k}\right)},
!! \f}
!! where \f$v_x^{i,j,k}\f$ is the maximum speed in \f$x\f$ direction for the cell \f$(i,j,k)\f$ and \f$c_f^{i,j,k}\f$ is the speed of sound for
!! ionized fluid computed as \f$c_f^{i,j,k}=\sqrt{\left|\frac{2p_{mag}+\gamma p}{\rho^{i,j,k}}\right|}\f$, where \f$p\f$ stands for pressure,
!! \f$p_{mag}\f$ is pressure of magnetic field, \f$\gamma\f$ is adiabatic index for ionized fluid and \f$\rho^{i,j,k}\f$ is fluid density in the cell
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

   real :: dt_ion             !< final timestep for ionized fluids
   real :: c_ion              !< maximum speed at which information travels in the ionized fluid

contains

   subroutine timestep_ion
      use types,         only: component_fluid
      use arrays,        only: u, b
      use grid,          only: ks, ke, is, ie, js, je
      use fluidindex,    only: nvar, ibx, iby, ibz
      use timestepfuncs, only: compute_c_max, compute_dt
#ifdef ISO_LOCAL
      use arrays,      only: cs_iso2_arr
#endif /* ISO_LOCAL */

      implicit none

      real :: cx = 0.0            !< maximum velocity for X direction
      real :: cy = 0.0            !< maximum velocity for Y direction
      real :: cz = 0.0            !< maximum velocity for Z direction
      real :: cs = 0.0            !< speed of sound

! locals
      real                           :: pmag = 0.0, bx = 0.0, by = 0.0, bz = 0.0, ps = 0.0, p = 0.0, c_max = 0.0
      integer                        :: i, j, k, ip, jp, kp
      type(component_fluid), pointer :: fl

      fl => nvar%ion

      do k = ks, ke
         kp = mod(k,ke)+1
         do j = js, je
            jp = mod(j,je)+1
            do i = is, ie
               ip = mod(i,ie)+1

               bx = 0.5*(b(ibx,i,j,k) + b(ibx,ip,j,k))
               by = 0.5*(b(iby,i,j,k) + b(iby,i,jp,k))
               bz = 0.5*(b(ibz,i,j,k) + b(ibz,i,j,kp))

               pmag = 0.5*(bx*bx + by*by + bz*bz)

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
