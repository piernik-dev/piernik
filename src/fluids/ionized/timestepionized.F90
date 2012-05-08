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
!! \brief (MH/JD) [R] %Timestep computation for the ionized fluid
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
! pulled by ANY
   implicit none

   private
   public :: timestep_ion

contains

   subroutine timestep_ion(cg, dt, c_ion)

      use constants,     only: zero, two, half, ndims
#ifdef MAGNETIC
      use constants,     only: xdim, ydim, zdim
      use domain,        only: dom
#endif /* MAGNETIC */
      use fluidindex,    only: flind
      use fluidtypes,    only: component_fluid
      use grid_cont,     only: grid_container
      use timestepfuncs, only: compute_c_max, compute_dt

      implicit none

      type(grid_container), pointer, intent(in) :: cg !< current grid container
      real, intent(out)                         :: dt !< resulting timestep
      real, intent(out)                         :: c_ion !< maximum speed at which information travels in the ionized fluid

      real, dimension(ndims) :: c !< maximum velocity for all directions
      real :: cs                  !< speed of sound

! locals
      real                           :: pmag, bx, by, bz, ps, p
      integer                        :: i, j, k
      type(component_fluid), pointer :: fl

      c(:) = 0.0; cs = 0.0;
      pmag = 0.0; bx = 0.0; by = 0.0; bz = 0.0; ps = 0.0; p = 0.0; c_ion = 0.0
      fl => flind%ion

      do k = cg%ks, cg%ke
         do j = cg%js, cg%je
            do i = cg%is, cg%ie

#ifdef MAGNETIC
               bx = (cg%b(xdim,i,j,k) + cg%b(xdim, i+dom%D_x, j,         k        ))/(1.+dom%D_x)
               by = (cg%b(ydim,i,j,k) + cg%b(ydim, i,         j+dom%D_y, k        ))/(1.+dom%D_y)
               bz = (cg%b(zdim,i,j,k) + cg%b(zdim, i,         j,         k+dom%D_z))/(1.+dom%D_z)

               pmag = half*(bx*bx + by*by + bz*bz)
#else /* !MAGNETIC */
               ! all_mag_boundaries has not been called so we cannot trust cg%b(xdim, cg%ie+dom%D_x:), cg%b(ydim,:cg%je+dom%D_y and cg%b(zdim,:,:, cg%ke+dom%D_z
               pmag = zero
#endif /* !MAGNETIC */

#ifdef ISO
               p  = cg%cs_iso2(i,j,k)*cg%u(fl%idn,i,j,k)
               ps = p + pmag
               cs = sqrt(abs(  (two*pmag+p)/cg%u(fl%idn,i,j,k)) )
#else /* !ISO */
               ps = (cg%u(fl%ien,i,j,k)-sum(cg%u(fl%imx:fl%imz,i,j,k)**2,1) &
                     /cg%u(fl%idn,i,j,k)*half)*(fl%gam_1)+(two-fl%gam)*pmag
               p  = ps - pmag
               cs = sqrt(abs(  (two*pmag+fl%gam*p)/cg%u(fl%idn,i,j,k)) )
#endif /* !ISO */
               call compute_c_max(fl, cs, i, j, k, c(:), c_ion, cg)
            enddo
         enddo
      enddo
      call compute_dt(c(:), dt, cg)
      fl%c  = c_ion

   end subroutine timestep_ion

end module timestepionized
