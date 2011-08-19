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

   use constants, only: varlen, ndims, LO, HI

   implicit none

   private
   public :: read_problem_par, init_prob, problem_pointers

   real :: pulse_size, pulse_vel_x, pulse_vel_y, pulse_vel_z, pulse_amp
   real, dimension(ndims, LO:HI) :: pulse_edge

   namelist /PROBLEM_CONTROL/  pulse_size, pulse_vel_x, pulse_vel_y, pulse_vel_z, pulse_amp

   real, dimension(:,:,:), allocatable, target :: inid
   character(len=varlen), parameter :: inid_n = "inid"

contains

!-----------------------------------------------------------------------------

   subroutine problem_pointers

      use dataio_user, only: user_vars_hdf5, problem_write_restart, problem_read_restart
      use types,       only: finalize_problem, cleanup_problem

      implicit none

      finalize_problem      => finalize_problem_adv
      cleanup_problem       => cleanup_adv
      user_vars_hdf5        => inid_var_hdf5
      problem_write_restart => write_IC_to_restart
      problem_read_restart  => read_IC_from_restart

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use dataio_pub, only: ierrh, par_file, namelist_errh, compare_namelist, cmdl_nml      ! QA_WARN required for diff_nml
      use dataio_pub, only: warn
      use domain,     only: dom
      use mpisetup,   only: ierr, rbuff, master, slave, buffer_dim, comm, proc, have_mpi, LAST
      use mpi,        only: MPI_DOUBLE_PRECISION

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

      if (pulse_size <= 0. .or. pulse_size >= 1.) then
         pulse_size = 0.5
         if (master) call warn("[initproblem:read_problem_par] Pulse width was invalid. Adjusted to 0.5.")
      endif
      if (pulse_amp <= 0.) then
         if (have_mpi) then
            pulse_amp = 1. + proc/real(LAST)
            pulse_size = 1.
         else
            pulse_amp = 2.
         endif
      endif

      pulse_edge(:, LO) = dom%C_(:) - dom%L_(:) * pulse_size/2.
      pulse_edge(:, HI) = dom%C_(:) + dom%L_(:) * pulse_size/2.

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine init_prob

      use constants,   only: xdim, ydim, zdim
      use diagnostics, only: my_allocate
      use fluidindex,  only: flind
      use global,      only: smalld, smallei
      use grid,        only: cga
      use grid_cont,   only: cg_list_element, grid_container

      implicit none

      integer :: i, j, k
      real    :: pulse_low_density, pulse_pressure
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      !BEWARE: hardcoded magic numbers
      pulse_low_density = smalld * 1e5
      pulse_pressure = smallei * flind%neu%gam_1 * 1e2

      call cga%get_root(cgl)
      do while (associated(cgl))
         cg => cgl%cg

      cg%b%arr(:, :, :, :) = 0.
      cg%u%arr(flind%neu%idn, :, :, :) = pulse_low_density

      ! Initialize density with uniform sphere
      do k = cg%ks, cg%ke
         if (cg%z(k) > pulse_edge(zdim, LO) .and. cg%z(k) < pulse_edge(zdim, HI)) then
            do j = cg%js, cg%je
               if (cg%y(j) > pulse_edge(ydim, LO) .and. cg%y(j) < pulse_edge(ydim, HI)) then
                  do i = cg%is, cg%ie
                     if (cg%x(i) > pulse_edge(xdim, LO) .and. cg%x(i) < pulse_edge(xdim, HI)) then
                        cg%u%arr(flind%neu%idn, i, j, k) = pulse_low_density * pulse_amp
                     endif
                  enddo
               endif
            enddo
         endif
      enddo

      where (cg%u%arr(flind%neu%idn, :, :, :) < smalld) cg%u%arr(flind%neu%idn, :, :, :) = smalld

      call my_allocate(inid, cg%n_(:), inid_n)
      inid(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) = cg%u%arr(flind%neu%idn, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke)

      cg%u%arr(flind%neu%imx, :, :, :) = pulse_vel_x * cg%u%arr(flind%neu%idn, :, :, :)
      cg%u%arr(flind%neu%imy, :, :, :) = pulse_vel_y * cg%u%arr(flind%neu%idn, :, :, :)
      cg%u%arr(flind%neu%imz, :, :, :) = pulse_vel_z * cg%u%arr(flind%neu%idn, :, :, :)
      do k = cg%ks, cg%ke
         do j = cg%js, cg%je
            do i = cg%is, cg%ie
               cg%u%arr(flind%neu%ien,i,j,k) = max(smallei,                                             &
                    &              pulse_pressure / flind%neu%gam_1        + &
                    &              0.5 * sum(cg%u%arr(flind%neu%imx:flind%neu%imz,i,j,k)**2,1) / cg%u%arr(flind%neu%idn,i,j,k))
            enddo
         enddo
      enddo
         cgl => cgl%nxt
      enddo

   end subroutine init_prob

!-----------------------------------------------------------------------------

   subroutine inid_var_hdf5(var, tab, ierrh, cg)

      use grid_cont, only: grid_container

      implicit none

      character(len=*), intent(in)                    :: var
      real(kind=4), dimension(:,:,:), intent(inout)   :: tab
      integer, intent(inout)                          :: ierrh
      type(grid_container), pointer, intent(in) :: cg

      ierrh = 0
      select case (trim(var))
         case (inid_n)
            tab(:,:,:) = real(inid(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), 4)
            case default
            ierrh = -1
      end select

   end subroutine inid_var_hdf5

!-----------------------------------------------------------------------------

   subroutine write_IC_to_restart(file_id, cg)

      use constants,   only: AT_NO_B
      use dataio_pub,  only: warn, die
      use dataio_hdf5, only: write_arr_to_restart
      use grid,        only: cga
      use grid_cont,   only: grid_container
      use hdf5,        only: HID_T

      implicit none

      integer(HID_T), intent(in) :: file_id
      type(grid_container), pointer, intent(in) :: cg

      real, dimension(:,:,:), pointer :: p3d

      if (ubound(cga%cg_all(:), dim=1) > 1) call die("[initproblem:write_IC_to_restart] multiple grid pieces per procesor not implemented yet") !nontrivial inid

      if (allocated(inid)) then
         if (associated(p3d)) nullify(p3d)
         p3d => inid
         call write_arr_to_restart(file_id, p3d, AT_NO_B, inid_n, cg)
      else
         call warn("[initproblem:write_IC_to_restart] Cannot store inid(:,:,:) in the restart file because it mysteriously deallocated.")
      endif

   end subroutine write_IC_to_restart

!-----------------------------------------------------------------------------

   subroutine read_IC_from_restart(file_id, cg)

      use constants,   only: AT_NO_B
      use dataio_pub,  only: warn, die
      use dataio_hdf5, only: read_arr_from_restart
      use diagnostics, only: my_allocate
      use grid,        only: cga
      use grid_cont,   only: grid_container
      use hdf5,        only: HID_T

      implicit none

      integer(HID_T), intent(in) :: file_id
      type(grid_container), pointer, intent(in) :: cg

      real, dimension(:,:,:), pointer :: p3d

      if (ubound(cga%cg_all(:), dim=1) > 1) call die("[initproblem:read_IC_from_restart] multiple grid pieces per procesor not implemented yet") !nontrivial inid

      if (associated(p3d)) nullify(p3d)
      if (.not.allocated(inid)) call my_allocate(inid, cg%n_(:), inid_n)
      p3d => inid(:,:,:)

      if (allocated(inid) .and. associated(p3d)) then
         call read_arr_from_restart(file_id, p3d, AT_NO_B, inid_n, cg)
         nullify(p3d)
      else
         call warn("[initproblem:read_IC_from_restart] Cannot read inid(:,:,:). It's really weird...")
      endif

   end subroutine read_IC_from_restart

!-----------------------------------------------------------------------------

   subroutine finalize_problem_adv

      use dataio_pub, only: msg, printinfo, warn, die
      use fluidindex, only: flind
      use grid,       only: cga
      use grid_cont,  only: cg_list_element, grid_container
      use mpi,        only: MPI_DOUBLE_PRECISION, MPI_SUM, MPI_MIN, MPI_MAX, MPI_IN_PLACE
      use mpisetup,   only: master, comm, ierr

      implicit none

      integer            :: i, j, k
      real, dimension(2) :: norm, dev
      real               :: dini
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      if (ubound(cga%cg_all(:), dim=1) > 1) call die("[initproblem:finalize_problem_adv] multiple grid pieces per procesor not implemented yet") !nontrivial inid

      if (.not. allocated(inid)) then
         if (master) call warn("[initproblem:finalize_problem_adv] Cannot compare results with the initial conditions. Perhaps there is no 'inid' array in the restart file?")
         return
      endif

      norm(:) = 0.
      dev(1) = huge(1.0)
      dev(2) = -dev(1)

      call cga%get_root(cgl)
      do while (associated(cgl))
         cg => cgl%cg

      do k = cg%ks, cg%ke
         do j = cg%js, cg%je
            do i = cg%is, cg%ie
               dini =  inid(i, j, k)
               norm(1) = norm(1) + (dini - cg%u%arr(flind%neu%idn, i, j, k))**2
               norm(2) = norm(2) + dini**2
               dev(1) = min(dev(1), (dini - cg%u%arr(flind%neu%idn, i, j, k))/dini)
               dev(2) = max(dev(2), (dini - cg%u%arr(flind%neu%idn, i, j, k))/dini)
            enddo
         enddo
      enddo
         cgl => cgl%nxt
      enddo

      call MPI_Allreduce(MPI_IN_PLACE, norm,   2, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, dev(1), 1, MPI_DOUBLE_PRECISION, MPI_MIN, comm, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, dev(2), 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm, ierr)

      if (master) then
         write(msg,'(a,f12.6,a,2f12.6)')"[initproblem:finalize_problem_adv] L2 error norm = ", sqrt(norm(1)/norm(2)), ", min and max error = ", dev(1:2)
         call printinfo(msg)
      endif

   end subroutine finalize_problem_adv

!-----------------------------------------------------------------------------

   subroutine cleanup_adv

      use diagnostics, only: my_deallocate

      implicit none

      if (allocated(inid)) call my_deallocate(inid)

   end subroutine cleanup_adv

end module initproblem
