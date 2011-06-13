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
   public :: read_problem_par, init_prob, problem_pointers

   integer            :: norm_step
   real               :: d0, p0, bx0, by0, bz0, x0, y0, z0, r0, beta_cr, amp_cr
   real, dimension(:,:,:), allocatable :: aecr1

   namelist /PROBLEM_CONTROL/ d0, p0, bx0, by0, bz0, x0, y0, z0, r0, &
        &                     beta_cr, amp_cr, norm_step

contains

!-----------------------------------------------------------------------------

   subroutine problem_pointers

      use dataio_user, only: user_vars_hdf5
      use types,       only: problem_customize_solution, finalize_problem, cleanup_problem

      implicit none

      problem_customize_solution => check_norm
      finalize_problem           => check_norm
      cleanup_problem            => cleanup_crtest
      user_vars_hdf5             => crtest_analytic_ecr1

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use dataio_pub,  only: ierrh, par_file, namelist_errh, compare_namelist, cmdl_nml      ! QA_WARN required for diff_nml
      use dataio_pub,  only: die
      use diagnostics, only: my_allocate
      use domain,      only: dom
      use grid,        only: cga
      use grid_cont,   only: cg_list_element, grid_container
      use mpi,         only: MPI_INTEGER, MPI_DOUBLE_PRECISION
      use mpisetup,    only: ibuff, rbuff, buffer_dim, comm, ierr, master, slave

      implicit none

      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      d0           = 1.0e5     !< density
      p0           = 1.0       !< pressure
      bx0          =   0.      !< Magnetic field component x
      by0          =   0.      !< Magnetic field component y
      bz0          =   0.      !< Magnetic field component z
      x0           = 0.0       !< x-position of the blob
      y0           = 0.0       !< y-position of the blob
      z0           = 0.0       !< z-position of the blob
      r0           = 5.* minval([dom%Lx, dom%lY, dom%Lz]/dom%n_d) !< radius of the blob

      beta_cr      = 0.0       !< ambient level
      amp_cr       = 1.0       !< amplitude of the blob

      norm_step    = 10        !< how often to compute the norm (in steps)

      if (master) then

         diff_nml(PROBLEM_CONTROL)

         rbuff(1) = d0
         rbuff(2) = p0
         rbuff(3) = bx0
         rbuff(4) = by0
         rbuff(5) = bz0
         rbuff(6) = x0
         rbuff(7) = y0
         rbuff(8) = z0
         rbuff(9) = r0
         rbuff(10)= beta_cr
         rbuff(11)= amp_cr

         ibuff(1) = norm_step

      endif

      call MPI_Bcast(ibuff, buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_Bcast(rbuff, buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      if (slave) then

         d0           = rbuff(1)
         p0           = rbuff(2)
         bx0          = rbuff(3)
         by0          = rbuff(4)
         bz0          = rbuff(5)
         x0           = rbuff(6)
         y0           = rbuff(7)
         z0           = rbuff(8)
         r0           = rbuff(9)
         beta_cr      = rbuff(10)
         amp_cr       = rbuff(11)

         norm_step    = ibuff(1)

      endif

      if (r0 == 0.) call die("[initproblem:read_problem_par] r0 == 0")

      if (ubound(cga%cg_all(:), dim=1) > 1) call die("[initproblem:read_problem_par] multiple grid pieces per procesor not implemented yet") !nontrivial aecr1
      !> \todo provide mechanism for rank-3 user arrays in grid_container

      cgl => cga%cg_leafs%cg_l(1)
      do while (associated(cgl))
         cg => cgl%cg

         call my_allocate(aecr1, [cg%nxb, cg%nyb, cg%nzb], "aecr1")
         aecr1(:,:,:) = 0.
         cgl => cgl%nxt
      enddo

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine cleanup_crtest

      use diagnostics, only: my_deallocate

      implicit none

      call my_deallocate(aecr1)

   end subroutine cleanup_crtest

!-----------------------------------------------------------------------------

   subroutine init_prob

      use constants,      only: xdim, ydim, zdim
      use dataio_pub,     only: die, warn
      use domain,         only: has_dir, dom
      use fluidindex,     only: ibx, iby, ibz
      use grid,           only: cga
      use grid_cont,      only: cg_list_element, grid_container
      use initcosmicrays, only: gamma_crs, iarr_crs, ncrn, ncre, K_crn_paral, K_crn_perp
      use initionized,    only: idni, imxi, imzi, ieni, gamma_ion

      implicit none

      integer :: i, j, k
      integer, parameter :: icr = 1 !< Only first CR component
      integer :: iecr
      real    :: cs_iso, r2
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      iecr = -1
      if (ncrn+ncre >= icr) then
         iecr = iarr_crs(icr)
      else
         call die("[initproblem:init_prob] No CR components defined.")
      endif

! Uniform equilibrium state

      cs_iso = sqrt(p0/d0)

      if (.not.has_dir(xdim)) bx0 = 0. ! ignore B field in nonexistent direction to match the analytical solution
      if (.not.has_dir(ydim)) by0 = 0.
      if (.not.has_dir(zdim)) bz0 = 0.

      if ((bx0**2 + by0**2 + bz0**2 == 0.) .and. (K_crn_paral(icr) /= 0. .or. K_crn_perp(icr) /= 0.)) then
         call warn("[initproblem:init_prob] No magnetic field is set, K_crn_* also have to be 0.")
         K_crn_paral(icr) = 0.
         K_crn_perp(icr)  = 0.
      endif

      cgl => cga%cg_leafs%cg_l(1)
      do while (associated(cgl))
         cg => cgl%cg

         cg%b%arr(ibx, :, :, :) = bx0
         cg%b%arr(iby, :, :, :) = by0
         cg%b%arr(ibz, :, :, :) = bz0
         cg%u%arr(idni, :, :, :) = d0
         cg%u%arr(imxi:imzi, :, :, :) = 0.0

#ifndef ISO
         do k = 1, cg%nz
            do j = 1, cg%ny
               do i = 1, cg%nx
                  cg%u%arr(ieni,i,j,k) = p0/(gamma_ion-1.0) + &
                       &                 0.5*sum(cg%u%arr(imxi:imzi,i,j,k)**2,1)/cg%u%arr(idni,i,j,k) + &
                       &                 0.5*sum(cg%b%arr(:,i,j,k)**2,1)
               enddo
            enddo
         enddo
#endif /* !ISO */

#ifdef COSM_RAYS
         cg%u%arr(iecr, :, :, :)      =  beta_cr*cs_iso**2 * cg%u%arr(idni, :, :, :)/(gamma_crs(icr)-1.0)

! Explosions
         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               do i = cg%is, cg%ie
                  r2 = (cg%x(i)-x0)**2+(cg%y(j)-y0)**2+(cg%z(k)-z0)**2
                  if (cg%x(i)> 2*x0-dom%xmax .and. cg%y(j) > 2*y0-dom%ymax) &
                       cg%u%arr(iecr, i, j, k)= cg%u%arr(iecr, i, j, k) + amp_cr*exp(-r2/r0**2)
               enddo
            enddo
         enddo
#endif /* COSM_RAYS */

         cgl => cgl%nxt
      enddo

      call check_norm

   end subroutine init_prob

!-----------------------------------------------------------------------------

   subroutine compute_analytic_ecr1

      use dataio_pub,     only: die
      use grid,           only: cga
      use grid_cont,      only: cg_list_element, grid_container
      use initcosmicrays, only: iarr_crs, ncrn, ncre, K_crn_paral, K_crn_perp
      use global,         only: t

      implicit none

      integer            :: i, j, k
      real               :: r_par2, r_perp2, delx, dely, delz, magb, ampt, r0_par2, r0_perp2, bxn, byn, bzn
      integer            :: iecr
      integer, parameter :: icr = 1 !< Only first CR component
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      if (ubound(cga%cg_all(:), dim=1) > 1) call die("[initproblem:compute_analytic_ecr1] multiple grid pieces per procesor not implemented yet") !nontrivial aecr1

      iecr = -1

      if (ncrn+ncre >= icr) then
         iecr = iarr_crs(icr)
      else
         call die("[initproblem:compute_analytic_ecr1] No CR components defined.")
      endif

      magb = sqrt(bx0**2 + by0**2 + bz0**2)
      if (magb > 0.) then
         bxn = bx0 / magb
         byn = by0 / magb
         bzn = bz0 / magb
      else
         bxn = 0.0
         byn = 0.0
         bzn = 0.0
      endif

      r0_par2  = r0**2 + 4 * (K_crn_paral(icr) + K_crn_perp(icr)) * t
      r0_perp2 = r0**2 + 4 * K_crn_perp(icr) * t

      if (r0_par2 == 0. .or. r0_perp2 == 0.) call die("[initproblem:compute_analytic_ecr1] r0_par2 == 0. .or. r0_perp2 == 0.")

      ampt     = amp_cr * r0**2 / sqrt(r0_par2 * r0_perp2)

      cgl => cga%cg_leafs%cg_l(1)
      do while (associated(cgl))
         cg => cgl%cg

         do k = cg%ks, cg%ke
            delz = cg%z(k) - z0
            do j = cg%js, cg%je
               dely = cg%y(j) - y0
               do i = cg%is, cg%ie
                  delx = cg%x(i) - x0

                  r_par2 = (bxn*delx + byn*dely + bzn*delz)**2 ! square of the distance form the center of the bump in direction parallel to the magnetic field
                  r_perp2 = delx**2 + dely**2 + delz**2 - r_par2
                  aecr1(i-cg%is+1, j-cg%js+1, k-cg%ks+1) = ampt * exp( - r_par2/r0_par2 - r_perp2/r0_perp2)

               enddo
            enddo
         enddo

         cgl => cgl%nxt
      enddo

   end subroutine compute_analytic_ecr1

!-----------------------------------------------------------------------------

   subroutine check_norm

      use dataio_pub,     only: code_progress, halfstep, msg, die, printinfo
      use constants,      only: PIERNIK_FINISHED
      use global,         only: t, nstep
      use grid,           only: cga
      use grid_cont,      only: cg_list_element, grid_container
      use initcosmicrays, only: iarr_crs, ncrn, ncre
      use mpisetup,       only: master, comm, ierr
      use mpi,            only: MPI_DOUBLE_PRECISION, MPI_SUM, MPI_MIN, MPI_MAX, MPI_IN_PLACE

      implicit none

      integer            :: i, j, k
      real, dimension(2) :: norm, dev
      real               :: crt
      integer            :: iecr
      integer, parameter :: icr = 1 !< Only first CR component
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      if (ubound(cga%cg_all(:), dim=1) > 1) call die("[] multiple grid pieces per procesor not implemented yet") !nontrivial aecr1

      iecr = -1

      if (ncrn+ncre >= icr) then
         iecr = iarr_crs(icr)
      else
         call die("[initproblem:check_norm] No CR components defined.")
      endif

      if (code_progress < PIERNIK_FINISHED .and. (mod(nstep, norm_step) /=0 .or. halfstep)) return

      call compute_analytic_ecr1

      norm(:) = 0.
      dev(1) = huge(1.0)
      dev(2) = -dev(1)

      cgl => cga%cg_leafs%cg_l(1)
      do while (associated(cgl))
         cg => cgl%cg
         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               do i = cg%is, cg%ie
                  crt = aecr1(i-cg%is+1, j-cg%js+1, k-cg%ks+1)
                  norm(1) = norm(1) + (crt - cg%u%arr(iecr, i, j, k))**2
                  norm(2) = norm(2) + crt**2
                  dev(1) = min(dev(1), (crt - cg%u%arr(iecr, i, j, k)))
                  dev(2) = max(dev(2), (crt - cg%u%arr(iecr, i, j, k)))
               enddo
            enddo
         enddo
         cgl => cgl%nxt
      enddo

      call MPI_Allreduce(MPI_IN_PLACE, norm,   2, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, dev(1), 1, MPI_DOUBLE_PRECISION, MPI_MIN, comm, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, dev(2), 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm, ierr)

      if (master) then
         if (norm(2) /= 0) then
            write(msg,'(a,f12.5,a,2f12.5)')"[initproblem:check_norm] L2 error norm = ", sqrt(norm(1)/norm(2)), " min and max error = ", dev(1:2)
         else
            write(msg,'(a,2f12.5)')"[initproblem:check_norm] Cannot compute L2 error norm, min and max error = ", dev(1:2)
         endif
         call printinfo(msg)
      endif

   end subroutine check_norm

!-----------------------------------------------------------------------------

   subroutine crtest_analytic_ecr1(var, tab, ierrh, cg)

      use grid_cont,      only: grid_container
      use initcosmicrays, only: iarr_crs

      implicit none

      character(len=*), intent(in)                    :: var
      real(kind=4), dimension(:,:,:), intent(inout)   :: tab
      integer, intent(inout)                          :: ierrh
      type(grid_container), pointer, intent(in)       :: cg

      call compute_analytic_ecr1

      ierrh = 0
      select case (trim(var))
         case ("acr1")
            tab(:,:,:) = real(aecr1(:,:,:), 4)
         case ("err1")
            tab(:,:,:) = real(aecr1(:,:,:) - cg%u%arr(iarr_crs(1), cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), 4)
         case default
            ierrh = -1
      end select

   end subroutine crtest_analytic_ecr1

end module initproblem
