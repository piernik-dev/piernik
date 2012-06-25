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
!! \brief Timestep computation for the neutral fluid
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
! pulled by ANY
   implicit none

   private
   public :: timestep_neu

contains

   subroutine timestep_neu(cg, dt, c_neu)

      use constants,     only: ndims
      use fluidindex,    only: flind
      use fluidtypes,    only: component_fluid
      use grid_cont,     only: grid_container
      use timestepfuncs, only: compute_c_max, compute_dt
#ifndef ISO
      use constants,     only: half
#endif /* !ISO */

      implicit none

      type(grid_container), pointer, intent(in) :: cg !< current grid container
      real, intent(out)                         :: dt !< resulting timestep
      real, intent(out)                         :: c_neu !< maximum speed at which information travels in the neutral fluid

      real, dimension(ndims) :: c !< maximum velocity for all directions
      real :: cs                  !< speed of sound

! locals

      real                           :: p
      integer                        :: i, j, k
      class(component_fluid), pointer :: fl

      c(:) = 0.0; cs = 0.0; p = 0.0; c_neu = 0.0

      fl => flind%neu

      do k = cg%ks, cg%ke
         do j = cg%js, cg%je
            do i = cg%is, cg%ie

#ifdef ISO
               p  = cg%cs_iso2(i,j,k)*cg%u(fl%idn,i,j,k)
               cs = sqrt(cg%cs_iso2(i,j,k))
#else /* !ISO */
               p  = (cg%u(fl%ien,i,j,k)-half*sum(cg%u(fl%imx:fl%imz,i,j,k)**2,1)/cg%u(fl%idn,i,j,k))*(fl%gam_1)

               cs = sqrt(abs(  (fl%gam*p)/cg%u(fl%idn,i,j,k)) )
#endif /* !ISO */
               call compute_c_max(fl, cs, i, j, k, c(:), c_neu, cg)
            enddo
         enddo
      enddo
      call compute_dt(c(:), dt, cg)
      fl%c  = c_neu

   end subroutine timestep_neu

end module timestepneutral
