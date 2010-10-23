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
      use types,       only: component_fluid
      use arrays,      only: u, b
      use constants,   only: big
      use grid,        only: dx, dy, dz, nb, ks, ke, is, ie, js, je, nxd, nyd, nzd
      use mpisetup,    only: MPI_DOUBLE_PRECISION, MPI_MIN, MPI_MAX, comm, ierr, cfl
      use fluidindex,  only: nvar, ibx, iby, ibz
#ifdef ISO_LOCAL
      use arrays,      only: cs_iso2_arr
#endif /* ISO_LOCAL */

      implicit none

      real :: dt_proc             !< minimum timestep for the current processor
      real :: dt_all              !< minimum timestep for all the processors
      real :: c_max_all           !< maximum speed for the fluid for all the processors
      real :: dt_proc_x           !< timestep computed for X direction for the current processor
      real :: dt_proc_y           !< timestep computed for Y direction for the current processor
      real :: dt_proc_z           !< timestep computed for Z direction for the current processor
      real :: cx                  !< maximum velocity for X direction
      real :: cy                  !< maximum velocity for Y direction
      real :: cz                  !< maximum velocity for Z direction
      real :: vx                  !< velocity in X direction computed for current cell
      real :: vy                  !< velocity in Y direction computed for current cell
      real :: vz                  !< velocity in Z direction computed for current cell
      real :: cs                  !< speed of sound

! locals
      real                           :: pmag, bx, by, bz, ps, p, c_max
      integer                        :: i, j, k, ip, jp, kp
      type(component_fluid), pointer :: fl

      fl => nvar%ion

      cx    = 0.0
      cy    = 0.0
      cz    = 0.0
      c_max = 0.0

      do k = ks, ke
         kp = mod(k,ke)+1
         do j = js, je
            jp = mod(j,je)+1
            do i = is, ie
               ip = mod(i,ie)+1

               vx = abs(u(fl%imx,i,j,k)/u(fl%idn,i,j,k))
               vy = abs(u(fl%imy,i,j,k)/u(fl%idn,i,j,k))
               vz = abs(u(fl%imz,i,j,k)/u(fl%idn,i,j,k))

               bx = 0.5*(b(ibx,i,j,k) + b(ibx,ip,j,k))
               by = 0.5*(b(iby,i,j,k) + b(iby,i,jp,k))
               bz = 0.5*(b(ibz,i,j,k) + b(ibz,i,j,kp))

               pmag = 0.5*(bx*bx + by*by + bz*bz)

#ifdef ISO
               p  = fl%cs2*u(fl%idn,i,j,k)
               ps = p + pmag
               cs = sqrt(abs(  (2.*pmag+p)/u(fl%idn,i,j,k)) )
#else /* ISO */
               ps = (u(fl%ien,i,j,k)-sum(u(fl%imx:fl%imz,i,j,k)**2,1) &
                     /u(fl%idn,i,j,k)*0.5)*(fl%gam-1.)+(2.-fl%gam)*pmag
               p  = ps - pmag
               cs = sqrt(abs(  (2.*pmag+fl%gam*p)/u(fl%idn,i,j,k)) )
#endif /* ISO */

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

      c_ion  = c_max_all
      dt_ion = cfl*dt_all

   end subroutine timestep_ion

end module timestepionized
