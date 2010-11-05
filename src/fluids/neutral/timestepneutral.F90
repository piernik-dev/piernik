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
!! \brief (MH/JD) (doxy comments ready) Timestep computation for the neutral fluid
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
! pulled by NEUTRAL
   implicit none
   private
   public :: dt_neu, c_neu, timestep_neu
   real   :: dt_neu               !< final timestep for neutral fluid
   real   :: c_neu                !< maximum speed at which information travels in the neutral fluid

contains

   subroutine timestep_neu
      use types,         only: component_fluid
      use arrays,        only: u
      use grid,          only: ks, ke, is, ie, js, je
      use fluidindex,    only: nvar
      use timestepfuncs, only: compute_c_max, compute_dt

      implicit none

      real :: cx                  !< maximum velocity for X direction
      real :: cy                  !< maximum velocity for Y direction
      real :: cz                  !< maximum velocity for Z direction
      real :: cs                  !< speed of sound

! locals

      real                           :: p, c_max
      integer                        :: i, j, k
      type(component_fluid), pointer :: fl

      cx = 0.0; cy = 0.0; cz = 0.0; cs = 0.0; p = 0.0; c_max = 0.0

      fl => nvar%neu

      do k = ks, ke
         do j = js, je
            do i = is, ie
#ifdef ISO
               p  = fl%cs2*u(fl%idn,i,j,k)
               cs = sqrt(fl%cs2)
#else /* !ISO */
               p  = (u(fl%ien,i,j,k)-sum(u(fl%imx:fl%imz,i,j,k)**2,1) &
                     /u(fl%idn,i,j,k)/2.)*(fl%gam_1)

               cs = sqrt(abs(  (fl%gam*p)/u(fl%idn,i,j,k)) )
#endif /* !ISO */
               call compute_c_max(fl,cs,i,j,k,cx,cy,cz,c_max)
            enddo
         enddo
      enddo
      call compute_dt(fl,cx,cy,cz,c_max,c_neu,dt_neu)

   end subroutine timestep_neu

end module timestepneutral
