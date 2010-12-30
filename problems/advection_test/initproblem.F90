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
#include "macros.h"

module initproblem

   implicit none

   private
   public :: read_problem_par, init_prob

   real :: pulse_size, pulse_vel_x, pulse_vel_y, pulse_vel_z, pulse_amp, &
        &  pulse_x_min, pulse_x_max, pulse_y_min, pulse_y_max, pulse_z_min, pulse_z_max

   namelist /PROBLEM_CONTROL/  pulse_size, pulse_vel_x, pulse_vel_y, pulse_vel_z, pulse_amp

contains

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use dataio_pub,    only: ierrh, par_file, namelist_errh, compare_namelist      ! QA_WARN required for diff_nml
      use dataio_pub,    only: die
      use grid,          only: cg
      use mpisetup,      only: ierr, rbuff, master, slave, buffer_dim, comm
      use mpi,           only: MPI_DOUBLE_PRECISION

      implicit none

      ! namelist default parameter values
      pulse_size   = 0.5                   !< "fill factor" in each direction
      pulse_vel_x  = 0.0                   !< pulse velocity in x-direction
      pulse_vel_y  = 0.0                   !< pulse velocity in y-direction
      pulse_vel_z  = 0.0                   !< pulse velocity in z-direction
      pulse_amp    = 2.0                   !< pulse relative amplitude

      if (master) then

         diff_nml(PROBLEM_CONTROL)

         rbuff(1) = pulse_size
         rbuff(2) = pulse_vel_x
         rbuff(3) = pulse_vel_y
         rbuff(4) = pulse_vel_z
         rbuff(5) = pulse_amp

      endif

      call MPI_Bcast(rbuff, buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      if (slave) then

         pulse_size  = rbuff(1)
         pulse_vel_x = rbuff(2)
         pulse_vel_y = rbuff(3)
         pulse_vel_z = rbuff(4)
         pulse_amp   = rbuff(5)

      endif

      if (pulse_size <= 0. .or. pulse_size >= 1.) call die("[initproblem:read_problem_par] Pulse width should be between 0. and 1.")

      pulse_x_min = (cg%xmax+cg%xmin)/2. - (cg%xmax-cg%xmin)*pulse_size/2.
      pulse_x_max = (cg%xmax+cg%xmin)/2. + (cg%xmax-cg%xmin)*pulse_size/2.
      pulse_y_min = (cg%ymax+cg%ymin)/2. - (cg%ymax-cg%ymin)*pulse_size/2.
      pulse_y_max = (cg%ymax+cg%ymin)/2. + (cg%ymax-cg%ymin)*pulse_size/2.
      pulse_z_min = (cg%zmax+cg%zmin)/2. - (cg%zmax-cg%zmin)*pulse_size/2.
      pulse_z_max = (cg%zmax+cg%zmin)/2. + (cg%zmax-cg%zmin)*pulse_size/2.

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine init_prob

      use arrays,        only: u, b
      use grid,          only: cg
      use fluidindex,    only: nvar
      use mpisetup,      only: smalld, smallei

      implicit none

      integer :: i, j, k
      real    :: pulse_low_density, pulse_pressure

      !BEWARE: hardcoded magic numbers
      pulse_low_density = smalld * 1e5
      pulse_pressure = smallei * nvar%neu%gam_1 * 1e2

      b(:, :, :, :) = 0.
      u(nvar%neu%idn, :, :, :) = pulse_low_density

      ! Initialize density with uniform sphere
      do i = cg%is, cg%ie
         if (cg%x(i) > pulse_x_min .and. cg%x(i) < pulse_x_max) then
            do j = cg%js, cg%je
               if (cg%y(j) > pulse_y_min .and. cg%y(j) < pulse_y_max) then
                  do k = cg%ks, cg%ke
                     if (cg%z(k) > pulse_z_min .and. cg%z(k) < pulse_z_max) then
                        u(nvar%neu%idn, i, j, k) = pulse_low_density * pulse_amp
                     endif
                  enddo
               endif
            enddo
         endif
      enddo

      where (u(nvar%neu%idn, :, :, :) < smalld) u(nvar%neu%idn, :, :, :) = smalld

      u(nvar%neu%imx, :, :, :) = pulse_vel_x * u(nvar%neu%idn, :, :, :)
      u(nvar%neu%imy, :, :, :) = pulse_vel_y * u(nvar%neu%idn, :, :, :)
      u(nvar%neu%imz, :, :, :) = pulse_vel_z * u(nvar%neu%idn, :, :, :)
      do k = cg%ks, cg%ke
         do j = cg%js, cg%je
            do i = cg%is, cg%ie
               u(nvar%neu%ien,i,j,k) = max(smallei,                                             &
                    &              pulse_pressure / nvar%neu%gam_1        + &
                    &              0.5 * sum(u(nvar%neu%imx:nvar%neu%imz,i,j,k)**2,1) / u(nvar%neu%idn,i,j,k))
            enddo
         enddo
      enddo

   end subroutine init_prob

end module initproblem
