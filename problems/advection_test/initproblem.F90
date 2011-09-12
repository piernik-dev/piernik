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

   use constants, only: dsetnamelen, ndims, LO, HI

   implicit none

   private
   public :: read_problem_par, init_prob, problem_pointers

   ! namelist parameters
   real :: pulse_size, pulse_vel_x, pulse_vel_y, pulse_vel_z, pulse_amp
   namelist /PROBLEM_CONTROL/  pulse_size, pulse_vel_x, pulse_vel_y, pulse_vel_z, pulse_amp

   ! other private data
   real, dimension(ndims, LO:HI) :: pulse_edge
   character(len=dsetnamelen), parameter :: inid_n = "init_dens"

contains

!-----------------------------------------------------------------------------

   subroutine problem_pointers

      use dataio_user, only: user_vars_hdf5, problem_read_restart
      use user_hooks,  only: finalize_problem

      implicit none

      finalize_problem      => finalize_problem_adv
      user_vars_hdf5        => inid_var_hdf5
      problem_read_restart  => register_user_var

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use dataio_pub, only: ierrh, par_file, namelist_errh, compare_namelist, cmdl_nml, lun, getlun      ! QA_WARN required for diff_nml
      use dataio_pub, only: warn
      use domain,     only: dom
      use mpisetup,   only: ierr, rbuff, master, slave, buffer_dim, comm, proc, have_mpi, LAST, FIRST
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

      call MPI_Bcast(rbuff, buffer_dim, MPI_DOUBLE_PRECISION, FIRST, comm, ierr)

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

      use constants,   only: xdim, ydim, zdim, INT4
      use fluidindex,  only: flind
      use global,      only: smalld, smallei
      use grid,        only: all_cg
      use gc_list,     only: cg_list_element
      use grid_cont,   only: grid_container

      implicit none

      integer :: i, j, k
      real    :: pulse_low_density, pulse_pressure
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      !BEWARE: hardcoded magic numbers
      pulse_low_density = smalld * 1e5
      pulse_pressure = smallei * flind%neu%gam_1 * 1e2

      cgl => all_cg%first
      do while (associated(cgl))
         cg => cgl%cg

         cg%b%arr(:, :, :, :) = 0.
         cg%u%arr(flind%neu%idn, :, :, :) = pulse_low_density

         ! Initialize density with uniform box
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

         ! Make uniform, completely boring flow
         cg%u%arr(flind%neu%imx, :, :, :) = pulse_vel_x * cg%u%arr(flind%neu%idn, :, :, :)
         cg%u%arr(flind%neu%imy, :, :, :) = pulse_vel_y * cg%u%arr(flind%neu%idn, :, :, :)
         cg%u%arr(flind%neu%imz, :, :, :) = pulse_vel_z * cg%u%arr(flind%neu%idn, :, :, :)

         ! Set up the internal energy
         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               do i = cg%is, cg%ie
                  cg%u%arr(flind%neu%ien,i,j,k) = max(smallei,                                             &
                       &              pulse_pressure / flind%neu%gam_1        + &
                       &              0.5 * sum(cg%u%arr(flind%neu%imx:flind%neu%imz,i,j,k)**2,1) / cg%u%arr(flind%neu%idn,i,j,k))
               enddo
            enddo
         enddo

         ! Save the initial density
         call register_user_var(0_INT4, cg)
         i = cg%get_na_ind(inid_n)
         cg%q(i)%arr(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) = cg%u%arr(flind%neu%idn, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke)

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
      type(grid_container), pointer, intent(in)       :: cg

      real, dimension(:,:,:), pointer :: inid

      ierrh = 0
      inid => cg%get_na_ptr(var)
      if (associated(inid)) then
         tab(:,:,:) = real(inid(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), 4)
      else
         ierrh = -1
      endif

   end subroutine inid_var_hdf5

!-----------------------------------------------------------------------------

   subroutine register_user_var(file_id, cg)

      use constants,   only: AT_NO_B
      use grid_cont,   only: grid_container
      use hdf5,        only: HID_T

      implicit none

      integer(HID_T), intent(in) :: file_id
      type(grid_container), pointer, intent(in) :: cg

      call cg%add_na(inid_n, AT_NO_B)

      if (.false.) write(*,*) file_id ! QA_WARN suppress compiler warnings on unused files

   end subroutine register_user_var

!-----------------------------------------------------------------------------

   subroutine finalize_problem_adv

      use constants,  only: I_ONE, I_TWO
      use dataio_pub, only: msg, printinfo, warn
      use fluidindex, only: flind
      use grid,       only: all_cg
      use gc_list,    only: cg_list_element
      use grid_cont,  only: grid_container
      use mpi,        only: MPI_DOUBLE_PRECISION, MPI_SUM, MPI_MIN, MPI_MAX, MPI_IN_PLACE
      use mpisetup,   only: master, comm, ierr

      implicit none

      integer            :: i, j, k
      real, dimension(2) :: norm, dev
      real               :: dini
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg
      real, dimension(:,:,:), pointer :: inid

      norm(:) = 0.
      dev(1) = huge(1.0)
      dev(2) = -dev(1)

      cgl => all_cg%first
      do while (associated(cgl))
         cg => cgl%cg

         inid => cg%get_na_ptr(inid_n)
         if (.not. associated(inid))then
            if (master) call warn("[initproblem:finalize_problem_adv] Cannot compare results with the initial conditions.")
            return
         endif

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

      call MPI_Allreduce(MPI_IN_PLACE, norm,   I_TWO, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, dev(1), I_ONE, MPI_DOUBLE_PRECISION, MPI_MIN, comm, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, dev(2), I_ONE, MPI_DOUBLE_PRECISION, MPI_MAX, comm, ierr)

      if (master) then
         write(msg,'(a,f12.6,a,2f12.6)')"[initproblem:finalize_problem_adv] L2 error norm = ", sqrt(norm(1)/norm(2)), ", min and max error = ", dev(1:2)
         call printinfo(msg)
      endif

   end subroutine finalize_problem_adv

end module initproblem
