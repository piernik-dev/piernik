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

   real, dimension(:,:,:), allocatable :: inid

contains

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use dataio_pub,    only: ierrh, par_file, namelist_errh, compare_namelist, cmdl_nml      ! QA_WARN required for diff_nml
      use dataio_pub,    only: die, user_vars_hdf5
      use mpisetup,      only: ierr, rbuff, master, slave, buffer_dim, comm, dom
      use mpi,           only: MPI_DOUBLE_PRECISION
      use types,         only: finalize_problem, cleanup_problem

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

      pulse_x_min = (dom%xmax+dom%xmin)/2. - (dom%xmax-dom%xmin)*pulse_size/2.
      pulse_x_max = (dom%xmax+dom%xmin)/2. + (dom%xmax-dom%xmin)*pulse_size/2.
      pulse_y_min = (dom%ymax+dom%ymin)/2. - (dom%ymax-dom%ymin)*pulse_size/2.
      pulse_y_max = (dom%ymax+dom%ymin)/2. + (dom%ymax-dom%ymin)*pulse_size/2.
      pulse_z_min = (dom%zmax+dom%zmin)/2. - (dom%zmax-dom%zmin)*pulse_size/2.
      pulse_z_max = (dom%zmax+dom%zmin)/2. + (dom%zmax-dom%zmin)*pulse_size/2.

      finalize_problem => finalize_problem_adv
      cleanup_problem  => cleanup_adv
      user_vars_hdf5   => inid_var_hdf5

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine init_prob

      use arrays,        only: u, b
      use grid,          only: cg
      use fluidindex,    only: flind
      use mpisetup,      only: smalld, smallei
      use diagnostics,   only: my_allocate

      implicit none

      integer :: i, j, k
      real    :: pulse_low_density, pulse_pressure

      !BEWARE: hardcoded magic numbers
      pulse_low_density = smalld * 1e5
      pulse_pressure = smallei * flind%neu%gam_1 * 1e2

      b(:, :, :, :) = 0.
      u(flind%neu%idn, :, :, :) = pulse_low_density

      ! Initialize density with uniform sphere
      do i = cg%is, cg%ie
         if (cg%x(i) > pulse_x_min .and. cg%x(i) < pulse_x_max) then
            do j = cg%js, cg%je
               if (cg%y(j) > pulse_y_min .and. cg%y(j) < pulse_y_max) then
                  do k = cg%ks, cg%ke
                     if (cg%z(k) > pulse_z_min .and. cg%z(k) < pulse_z_max) then
                        u(flind%neu%idn, i, j, k) = pulse_low_density * pulse_amp
                     endif
                  enddo
               endif
            enddo
         endif
      enddo

      where (u(flind%neu%idn, :, :, :) < smalld) u(flind%neu%idn, :, :, :) = smalld

      call my_allocate(inid, [cg%nxb, cg%nyb, cg%nzb], "inid")
      inid(1:cg%nxb, 1:cg%nyb, 1:cg%nzb) = u(flind%neu%idn, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke)

      u(flind%neu%imx, :, :, :) = pulse_vel_x * u(flind%neu%idn, :, :, :)
      u(flind%neu%imy, :, :, :) = pulse_vel_y * u(flind%neu%idn, :, :, :)
      u(flind%neu%imz, :, :, :) = pulse_vel_z * u(flind%neu%idn, :, :, :)
      do k = cg%ks, cg%ke
         do j = cg%js, cg%je
            do i = cg%is, cg%ie
               u(flind%neu%ien,i,j,k) = max(smallei,                                             &
                    &              pulse_pressure / flind%neu%gam_1        + &
                    &              0.5 * sum(u(flind%neu%imx:flind%neu%imz,i,j,k)**2,1) / u(flind%neu%idn,i,j,k))
            enddo
         enddo
      enddo

   end subroutine init_prob

!-----------------------------------------------------------------------------

   subroutine inid_var_hdf5(var, tab, ierrh)

      implicit none

      character(len=*), intent(in)                    :: var
      real(kind=4), dimension(:,:,:), intent(inout)   :: tab
      integer, intent(inout)                          :: ierrh

      ierrh = 0
      select case (trim(var))
         case ("inid")
            tab(:,:,:) = inid(:,:,:)
            case default
            ierrh = -1
      end select

   end subroutine inid_var_hdf5

!-----------------------------------------------------------------------------

   subroutine finalize_problem_adv

      use arrays,        only: u
      use dataio_pub,    only: msg, printinfo
      use grid,          only: cg
      use mpisetup,      only: master, comm3d, ierr
      use mpi,           only: MPI_DOUBLE_PRECISION, MPI_SUM, MPI_MIN, MPI_MAX, MPI_IN_PLACE
      use fluidindex,    only: flind

      implicit none

      integer            :: i, j, k
      real, dimension(2) :: norm, dev
      real               :: dini

      norm(:) = 0.
      dev(1) = huge(1.0)
      dev(2) = -dev(1)

      do k = cg%ks, cg%ke
         do j = cg%js, cg%je
            do i = cg%is, cg%ie
               dini =  inid(i-cg%is+1, j-cg%js+1, k-cg%ks+1)
               norm(1) = norm(1) + (dini - u(flind%neu%idn, i, j, k))**2
               norm(2) = norm(2) + dini**2
               dev(1) = min(dev(1), (dini - u(flind%neu%idn, i, j, k))/dini)
               dev(2) = max(dev(2), (dini - u(flind%neu%idn, i, j, k))/dini)
            enddo
         enddo
      enddo
      call MPI_Allreduce(MPI_IN_PLACE, norm,   2, MPI_DOUBLE_PRECISION, MPI_SUM, comm3d, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, dev(1), 1, MPI_DOUBLE_PRECISION, MPI_MIN, comm3d, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, dev(2), 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm3d, ierr)

      if (master) then
         write(msg,'(a,f12.6,a,2f12.6)')"[initproblem:finalize_problem_adv] L2 error norm = ", sqrt(norm(1)/norm(2)), ", min and max error = ", dev(1:2)
         call printinfo(msg)
      endif

   end subroutine finalize_problem_adv

!-----------------------------------------------------------------------------

   subroutine cleanup_adv

      use diagnostics, only: my_deallocate

      implicit none

      call my_deallocate(inid)

   end subroutine cleanup_adv

end module initproblem
