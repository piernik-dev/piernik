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
#include "macros.h"

module initproblem

  use problem_pub, only: problem_name, run_id

  real :: pulse_size, pulse_vel_x, pulse_vel_y, pulse_vel_z, pulse_amp, pulse_low_density, pulse_pressure, &
       &  pulse_x_min, pulse_x_max, pulse_y_min, pulse_y_max, pulse_z_min, pulse_z_max

  namelist /PROBLEM_CONTROL/  problem_name, run_id, pulse_size, pulse_vel_x, pulse_vel_y, pulse_vel_z, pulse_amp

contains

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use grid,        only : xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, dz, nxd, nyd, nzd
      use errh,        only : namelist_errh, die
      use mpisetup,    only : ierr, rbuff, cbuff, proc, buffer_dim, comm, smalld, smallei, MPI_CHARACTER, MPI_DOUBLE_PRECISION
      use initneutral, only : gamma_neu
      use dataio_public, only : ierrh, msg, par_file
      use func,        only : compare_namelist

      implicit none


      ! namelist default parameter values
      problem_name = 'selfgrav_clump'      !< The default problem name
      run_id       = '_'                   !< Auxiliary run identifier
      pulse_size   = 0.5                   !< "fill factor" in each direction
      pulse_vel_x  = 0.0                   !< pulse velocity in x-direction
      pulse_vel_y  = 0.0                   !< pulse velocity in y-direction
      pulse_vel_z  = 0.0                   !< pulse velocity in z-direction
      pulse_amp    = 2.0                   !< pulse relative amplitude

      if (proc == 0) then

         diff_nml(PROBLEM_CONTROL)

         cbuff(1) = problem_name
         cbuff(2) = run_id

         rbuff(1) = pulse_size
         rbuff(2) = pulse_vel_x
         rbuff(3) = pulse_vel_y
         rbuff(4) = pulse_vel_z
         rbuff(5) = pulse_amp

      end if

      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      if (proc /= 0) then

         problem_name = cbuff(1)
         run_id       = cbuff(2)(1:3)

         pulse_size  = rbuff(1)
         pulse_vel_x = rbuff(2)
         pulse_vel_y = rbuff(3)
         pulse_vel_z = rbuff(4)
         pulse_amp   = rbuff(5)

      endif

      if (pulse_size <= 0. .or. pulse_size >= 1.) call die("[initproblem:read_problem_par] Pulse width should be between 0. and 1.")

      pulse_x_min = (xmax+xmin)/2. - (xmax-xmin)*pulse_size/2.
      pulse_x_max = (xmax+xmin)/2. + (xmax-xmin)*pulse_size/2.
      pulse_y_min = (ymax+ymin)/2. - (ymax-ymin)*pulse_size/2.
      pulse_y_max = (ymax+ymin)/2. + (ymax-ymin)*pulse_size/2.
      pulse_z_min = (zmax+zmin)/2. - (zmax-zmin)*pulse_size/2.
      pulse_z_max = (zmax+zmin)/2. + (zmax-zmin)*pulse_size/2.

      !BEWARE: hardcoded values
      pulse_low_density = smalld * 1e5
      pulse_pressure = smallei * (gamma_neu-1.0) * 1e5

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine init_prob

      use arrays,        only : u, b
      use grid,          only : x, y, z, is, ie, js, je, ks, ke
      use initneutral,   only : gamma_neu, idnn, imxn, imyn, imzn, ienn
      use mpisetup,      only : smalld, smallei

      implicit none

      integer                   :: i, j, k

      b(:,    :, :, :) = 0.
      u(idnn, :, :, :) = pulse_low_density

      ! Initialize density with uniform sphere
      do i = is, ie
         if (x(i) > pulse_x_min .and. x(i) < pulse_x_max) then
            do j = js, je
               if (y(j) > pulse_y_min .and. y(j) < pulse_y_max) then
                  do k = ks, ke
                     if (z(k) > pulse_z_min .and. z(k) < pulse_z_max) then
                        u(idnn, i, j, k) = pulse_low_density * pulse_amp
                     end if
                  end do
               end if
            end do
         end if
      end do

      where (u(idnn, :, :, :) < smalld) u(idnn, :, :, :) = smalld

      u(imxn, :, :, :) = pulse_vel_x * u(idnn,:,:,:)
      u(imyn, :, :, :) = pulse_vel_y * u(idnn,:,:,:)
      u(imzn, :, :, :) = pulse_vel_z * u(idnn,:,:,:)
      do k = ks, ke
         do j = js, je
            do i = is, ie
               u(ienn,i,j,k) = max(smallei,                                             &
                    &              pulse_pressure / (gamma_neu-1.0)        + &
                    &              0.5 * sum(u(imxn:imzn,i,j,k)**2,1) / u(idnn,i,j,k) + &
                    &              0.5 * sum(b(:,i,j,k)**2,1))
            enddo
         enddo
      enddo

   end subroutine init_prob

end module initproblem
