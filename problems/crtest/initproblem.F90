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

   use constants, only: dsetnamelen

   implicit none

   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   integer(kind=4)    :: norm_step
   real               :: d0, p0, bx0, by0, bz0, x0, y0, z0, r0, beta_cr, amp_cr
   character(len=dsetnamelen), parameter :: aecr1_n = "aecr"

   namelist /PROBLEM_CONTROL/ d0, p0, bx0, by0, bz0, x0, y0, z0, r0, beta_cr, amp_cr, norm_step

contains

!-----------------------------------------------------------------------------

   subroutine problem_pointers

      use dataio_user, only: user_vars_hdf5
      use user_hooks,  only: problem_customize_solution, finalize_problem

      implicit none

      problem_customize_solution => check_norm_wrapper
      finalize_problem           => check_norm
      user_vars_hdf5             => crtest_analytic_ecr1

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use cg_list_global, only: all_cg
      use constants,      only: I_ONE, I_TEN, AT_NO_B
      use dataio_pub,     only: nh      ! QA_WARN required for diff_nml
      use dataio_pub,     only: die
      use domain,         only: dom
      use mpisetup,       only: ibuff, rbuff, master, slave, piernik_MPI_Bcast

      implicit none

      d0           = 1.0e5     !< density
      p0           = 1.0       !< pressure
      bx0          =   0.      !< Magnetic field component x
      by0          =   0.      !< Magnetic field component y
      bz0          =   0.      !< Magnetic field component z
      x0           = 0.0       !< x-position of the blob
      y0           = 0.0       !< y-position of the blob
      z0           = 0.0       !< z-position of the blob
      r0           = 5.* minval(dom%L_(:)/dom%n_d) !< radius of the blob

      beta_cr      = 0.0       !< ambient level
      amp_cr       = 1.0       !< amplitude of the blob

      norm_step    = I_TEN     !< how often to compute the norm (in steps)

      if (master) then

         diff_nml(PROBLEM_CONTROL)

         rbuff(1)  = d0
         rbuff(2)  = p0
         rbuff(3)  = bx0
         rbuff(4)  = by0
         rbuff(5)  = bz0
         rbuff(6)  = x0
         rbuff(7)  = y0
         rbuff(8)  = z0
         rbuff(9)  = r0
         rbuff(10) = beta_cr
         rbuff(11) = amp_cr

         ibuff(1)  = norm_step

      endif

      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         d0        = rbuff(1)
         p0        = rbuff(2)
         bx0       = rbuff(3)
         by0       = rbuff(4)
         bz0       = rbuff(5)
         x0        = rbuff(6)
         y0        = rbuff(7)
         z0        = rbuff(8)
         r0        = rbuff(9)
         beta_cr   = rbuff(10)
         amp_cr    = rbuff(11)

         norm_step = int(ibuff(1), kind=4)

      endif

      if (r0 == 0.) call die("[initproblem:read_problem_par] r0 == 0")

      call all_cg%reg_var(aecr1_n, restart_mode = AT_NO_B)

      if (norm_step <= 0) norm_step = huge(I_ONE)

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine problem_initial_conditions

      use cg_leaves,      only: leaves
      use cg_list,        only: cg_list_element
      use constants,      only: xdim, ydim, zdim, HI
      use dataio_pub,     only: die
      use domain,         only: dom
      use fluidindex,     only: flind
      use fluidtypes,     only: component_fluid
      use func,           only: ekin, emag
      use grid_cont,      only: grid_container
      use initcosmicrays, only: gamma_crs, iarr_crs, ncrn, ncre

      implicit none

      class(component_fluid), pointer :: fl
      integer                         :: i, j, k, iecr
      integer, parameter              :: icr = 1 !< Only first CR component
      real                            :: cs_iso, r2
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg

      fl => flind%ion

      iecr = -1
      if (ncrn+ncre >= icr) then
         iecr = iarr_crs(icr)
      else
         call die("[initproblem:problem_initial_conditions] No CR components defined.")
      endif

! Uniform equilibrium state

      cs_iso = sqrt(p0/d0)

      if (.not.dom%has_dir(xdim)) bx0 = 0. ! ignore B field in nonexistent direction to match the analytical solution
      if (.not.dom%has_dir(ydim)) by0 = 0.
      if (.not.dom%has_dir(zdim)) bz0 = 0.

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         cg%b(xdim, :, :, :) = bx0
         cg%b(ydim, :, :, :) = by0
         cg%b(zdim, :, :, :) = bz0
         cg%u(fl%idn, :, :, :) = d0
         cg%u(fl%imx:fl%imz, :, :, :) = 0.0

#ifndef ISO
         cg%u(fl%ien,:,:,:) = p0/fl%gam_1 + emag(cg%b(xdim,:,:,:), cg%b(ydim,:,:,:), cg%b(zdim,:,:,:)) + &
              &               ekin(cg%u(fl%imx,:,:,:), cg%u(fl%imy,:,:,:), cg%u(fl%imz,:,:,:), cg%u(fl%idn,:,:,:))
#endif /* !ISO */

#ifdef COSM_RAYS
         cg%u(iecr,:,:,:) = beta_cr*fl%cs2 * cg%u(fl%idn,:,:,:)/(gamma_crs(icr)-1.0)

! Explosions
         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               do i = cg%is, cg%ie
                  r2 = (cg%x(i)-x0)**2+(cg%y(j)-y0)**2+(cg%z(k)-z0)**2
                  if (cg%x(i)> 2*x0-dom%edge(xdim, HI) .and. cg%y(j) > 2*y0-dom%edge(ydim, HI)) &
                       cg%u(iecr, i, j, k)= cg%u(iecr, i, j, k) + amp_cr*exp(-r2/r0**2)
               enddo
            enddo
         enddo
#endif /* COSM_RAYS */

         cgl => cgl%nxt
      enddo

      call check_norm

   end subroutine problem_initial_conditions

!-----------------------------------------------------------------------------

   subroutine compute_analytic_ecr1

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use dataio_pub,       only: die
      use global,           only: t
      use grid_cont,        only: grid_container
      use initcosmicrays,   only: iarr_crs, ncrn, ncre, K_crn_paral, K_crn_perp
      use named_array_list, only: qna

      implicit none

      integer                        :: i, j, k, iecr
      real                           :: r_par2, r_perp2, delx, dely, delz, magb, ampt, r0_par2, r0_perp2, bxn, byn, bzn
      integer, parameter             :: icr = 1 !< Only first CR component
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg

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

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         if (.not. qna%exists(aecr1_n)) call die("[initproblem:compute_analytic_ecr1] cannot find aecr1")

         do k = cg%ks, cg%ke
            delz = cg%z(k) - z0
            do j = cg%js, cg%je
               dely = cg%y(j) - y0
               do i = cg%is, cg%ie
                  delx = cg%x(i) - x0

                  r_par2 = (bxn*delx + byn*dely + bzn*delz)**2 ! square of the distance form the center of the bump in direction parallel to the magnetic field
                  r_perp2 = delx**2 + dely**2 + delz**2 - r_par2
                  cg%q(qna%ind(aecr1_n))%arr(i, j, k) = ampt * exp( - r_par2/r0_par2 - r_perp2/r0_perp2)

               enddo
            enddo
         enddo

         cgl => cgl%nxt
      enddo

   end subroutine compute_analytic_ecr1

!-----------------------------------------------------------------------------

   subroutine check_norm_wrapper(forward)
      implicit none
      logical, intent(in) :: forward
      call check_norm
      return
      if (.false. .and. forward) d0 = 0.0 ! suppress compiler warnings on unused arguments
   end subroutine check_norm_wrapper

!-----------------------------------------------------------------------------
   subroutine check_norm

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: PIERNIK_FINISHED, pSUM, pMIN, pMAX
      use dataio_pub,       only: code_progress, halfstep, msg, die, printinfo
      use global,           only: nstep
      use grid_cont,        only: grid_container
      use initcosmicrays,   only: iarr_crs, ncrn, ncre
      use mpisetup,         only: master, piernik_MPI_Allreduce
      use named_array_list, only: qna

      implicit none

      integer                        :: i, j, k, iecr
      real, dimension(2)             :: norm, dev
      real                           :: crt
      integer, parameter             :: icr = 1 !< Only first CR component
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg

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

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         if (.not. qna%exists(aecr1_n)) call die("[initproblem:check_norm] cannot find aecr1")
         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               do i = cg%is, cg%ie
                  crt = cg%q(qna%ind(aecr1_n))%arr(i, j, k)
                  norm(1) = norm(1) + (crt - cg%u(iecr, i, j, k))**2
                  norm(2) = norm(2) + crt**2
                  dev(1) = min(dev(1), (crt - cg%u(iecr, i, j, k)))
                  dev(2) = max(dev(2), (crt - cg%u(iecr, i, j, k)))
               enddo
            enddo
         enddo
         cgl => cgl%nxt
      enddo

      call piernik_MPI_Allreduce(norm,   pSUM)
      call piernik_MPI_Allreduce(dev(1), pMIN)
      call piernik_MPI_Allreduce(dev(2), pMAX)

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

      use dataio_pub,       only: die
      use grid_cont,        only: grid_container
      use initcosmicrays,   only: iarr_crs
      use named_array_list, only: qna

      implicit none

      character(len=*),               intent(in)    :: var
      real(kind=4), dimension(:,:,:), intent(inout) :: tab
      integer,                        intent(inout) :: ierrh
      type(grid_container), pointer,  intent(in)    :: cg

      call compute_analytic_ecr1

      if (.not. qna%exists(aecr1_n)) call die("[initproblem:crtest_analytic_ecr1] cannot find aecr1")

      ierrh = 0
      select case (trim(var))
         case ("acr1")
            tab(:,:,:) = real(cg%q(qna%ind(aecr1_n))%span(cg%ijkse), 4)
         case ("err1")
            tab(:,:,:) = real(cg%q(qna%ind(aecr1_n))%span(cg%ijkse) - cg%u(iarr_crs(1), cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), 4)
         case default
            ierrh = -1
      end select

   end subroutine crtest_analytic_ecr1

end module initproblem
