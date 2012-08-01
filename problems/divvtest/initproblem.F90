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
   public :: read_problem_par, init_prob, problem_pointers

   integer(kind=4) :: norm_step
   real            :: t_sn
   integer         :: n_sn
   real            :: d0, p0, bx0, by0, bz0, x0, y0, z0, r0, beta_cr, amp_cr, &
                      c_exp_x, c_exp_y, c_exp_z, c_rot_x, c_rot_y, c_rot_z

   character(len=dsetnamelen), parameter :: dvel_n = "dvel"

   namelist /PROBLEM_CONTROL/  d0, p0, bx0, by0, bz0, &
                               x0, y0, z0, r0, &
                               beta_cr, amp_cr, &
                               norm_step, &
                               c_exp_x, c_exp_y, c_exp_z, c_rot_x, c_rot_y, c_rot_z

contains
!-----------------------------------------------------------------------------
   subroutine problem_pointers

      use dataio_user, only: user_vars_hdf5, user_reg_var_restart

      implicit none

      user_vars_hdf5             => dvel_var_hdf5
      user_reg_var_restart       => register_user_var

   end subroutine problem_pointers
!-----------------------------------------------------------------------------
   subroutine register_user_var

      use cg_list_global, only: all_cg
      use constants,      only: AT_NO_B

      implicit none

      call all_cg%reg_var(dvel_n, restart_mode = AT_NO_B)

   end subroutine register_user_var
!-----------------------------------------------------------------------------
   subroutine read_problem_par

      use dataio_pub, only: ierrh, par_file, namelist_errh, compare_namelist, cmdl_nml, lun ! QA_WARN required for diff_nml
      use dataio_pub, only: die
      use domain,     only: dom
      use mpi,        only: MPI_INTEGER, MPI_DOUBLE_PRECISION
      use mpisetup,   only: ibuff, rbuff, buffer_dim, comm, mpi_err, master, slave, FIRST

      implicit none

      t_sn = 0.0

      d0           = 1.0e5     !< density
      p0           = 1.0       !< pressure
      bx0          =   0.      !< Magnetic field component x
      by0          =   0.      !< Magnetic field component y
      bz0          =   0.      !< Magnetic field component z
      x0           = 0.0       !< x-position of the blob
      y0           = 0.0       !< y-position of the blob
      z0           = 0.0       !< z-position of the blob
      r0           = 5.* minval(dom%L_(:)/dom%n_d(:), mask=dom%has_dir(:))  !< radius of the blob

      beta_cr      = 0.0       !< ambient level
      amp_cr       = 1.0       !< amplitude of the blob

      norm_step    = 10        !< how often to compute the norm (in steps)

      c_exp_x = 0.0
      c_exp_y = 0.0
      c_exp_z = 0.0
      c_rot_x = 0.0
      c_rot_y = 0.0
      c_rot_z = 0.0

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
         rbuff(10) = beta_cr
         rbuff(11) = amp_cr
         rbuff(12) = c_exp_x
         rbuff(13) = c_exp_y
         rbuff(14) = c_exp_z
         rbuff(15) = c_rot_x
         rbuff(16) = c_rot_y
         rbuff(17) = c_rot_x

         ibuff(1) = norm_step
      endif

      call MPI_Bcast(ibuff, buffer_dim, MPI_INTEGER,          FIRST, comm, mpi_err)
      call MPI_Bcast(rbuff, buffer_dim, MPI_DOUBLE_PRECISION, FIRST, comm, mpi_err)

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
         c_exp_x      = rbuff(12)
         c_exp_y      = rbuff(13)
         c_exp_z      = rbuff(14)
         c_rot_x      = rbuff(15)
         c_rot_y      = rbuff(16)
         c_rot_z      = rbuff(17)

         norm_step = int(ibuff(1), kind=4)
      endif

      if (r0 == 0.0) call die("[initproblem:read_problem_par] r0 == 0")

   end subroutine read_problem_par
!-----------------------------------------------------------------------------
   subroutine init_prob

      use constants,      only: xdim, ydim, zdim, big, pi
      use cg_list,        only: cg_list_element
#ifdef COSM_RAYS_SOURCES
      use cr_data,        only: icr_H1, icr_C12
#endif /* COSM_RAYS_SOURCES */
      use crhelpers,      only: divv_n, div_v
      use dataio_pub,     only: warn, die
      use domain,         only: dom, is_multicg
      use fluidindex,     only: flind
      use fluidtypes,     only: component_fluid
      use func,           only: emag, ekin
      use grid,           only: leaves
      use grid_cont,      only: grid_container
      use initcosmicrays, only: iarr_crn, iarr_crs, gamma_crn, K_crn_paral, K_crn_perp
      use mpisetup,       only: master
      use named_array,    only: qna

      implicit none

      class(component_fluid), pointer :: fl
      integer :: i, j, k, icr
      real    :: cs_iso
      real    :: xsn, ysn, zsn, r
      real    :: r2, maxv
      integer :: ipm, jpm, kpm
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg
#ifndef COSM_RAYS_SOURCES
      integer, parameter :: icr_H1 = 1, icr_C12 = 2
#endif /* !COSM_RAYS_SOURCES */

      call register_user_var

      fl => flind%ion

      if (.not. qna%exists(dvel_n)) then
         if (master) call warn("[initproblem: ] Cannot access dvel.")
      endif


      ! BEWARE: temporary fix
      xsn = x0
      ysn = y0
      zsn = z0

! Uniform equilibrium state

      cs_iso = sqrt(p0/d0)

      if (.not.dom%has_dir(xdim)) bx0 = 0. ! ignore B field in nonexistent direction
      if (.not.dom%has_dir(ydim)) by0 = 0.
      if (.not.dom%has_dir(zdim)) bz0 = 0.

      if ((bx0**2 + by0**2 + bz0**2 == 0.) .and. (any(K_crn_paral(:) /= 0.) .or. any(K_crn_perp(:) /= 0.))) then
         call warn("[initproblem:init_prob] No magnetic field is set, K_crn_* also have to be 0.")
         K_crn_paral(:) = 0.
         K_crn_perp(:)  = 0.
      endif

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         cg%b(xdim, :, :, :) = bx0
         cg%b(ydim, :, :, :) = by0
         cg%b(zdim, :, :, :) = bz0
         cg%u(fl%idn, :, :, :) = d0
         cg%u(fl%imx:fl%imz, :, :, :) = 0.0
         cg%q(qna%ind(dvel_n))%arr(:, :, :) = big


         do k = 1, cg%n_(zdim)
            do j = 1, cg%n_(ydim)
               do i = 1, cg%n_(xdim)

                  cg%u(fl%imx,i,j,k) = c_exp_x * cg%u(fl%idn,i,j,k) * cg%x(i)
                  cg%u(fl%imy,i,j,k) = c_exp_y * cg%u(fl%idn,i,j,k) * cg%y(j)
                  cg%u(fl%imz,i,j,k) = c_exp_z * cg%u(fl%idn,i,j,k) * cg%z(k)

                  r = sqrt( (cg%x(i))**2 + (cg%y(j))**2 + (cg%z(k))**2 )

                  cg%u(fl%imx,i,j,k) = cg%u(fl%imx,i,j,k) + c_rot_z * cg%u(fl%idn,i,j,k) * (-cg%y(j)/r) * sin(pi*r)
                  cg%u(fl%imy,i,j,k) = cg%u(fl%imy,i,j,k) + c_rot_z * cg%u(fl%idn,i,j,k) * ( cg%x(i)/r) * sin(pi*r)
                  cg%u(fl%imz,i,j,k) = cg%u(fl%imz,i,j,k) + 0.0
               enddo
            enddo
         enddo

         call div_v(flind%ion%pos, cg)
         cg%q(qna%ind(dvel_n))%arr(:, :, :) = cg%q(qna%ind(divv_n))%arr(:, :, :)

#ifndef ISO
         do k = 1, cg%n_(zdim)
            do j = 1, cg%n_(ydim)
               do i = 1, cg%n_(xdim)
                  cg%u(fl%ien,i,j,k) = p0 / fl%gam_1 + &
                       &               ekin(cg%u(fl%imx,i,j,k), cg%u(fl%imy,i,j,k), cg%u(fl%imz,i,j,k), cg%u(fl%idn,i,j,k)) + &
                       &               emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))
               enddo
            enddo
         enddo
#endif /* !ISO */

#ifdef COSM_RAYS
         do icr = 1, flind%crs%all
            cg%u(iarr_crs(icr), :, :, :) =  beta_cr*fl%cs2 * cg%u(fl%idn, :, :, :)/(gamma_crn(icr)-1.0)
         enddo

! Explosions
         do icr = 1, flind%crn%all
            do k = cg%ks, cg%ke
               do j = cg%js, cg%je
                  do i = cg%is, cg%ie
                     do ipm = -1, 1
                        do jpm = -1, 1
                           do kpm = -1, 1

                              r2 = (cg%x(i)-xsn+real(ipm)*dom%L_(xdim))**2+(cg%y(j)-ysn+real(jpm)*dom%L_(ydim))**2+(cg%z(k)-zsn+real(kpm)*dom%L_(zdim))**2
                              if (icr == icr_H1) then
                                 cg%u(iarr_crn(icr), i, j, k) = cg%u(iarr_crn(icr), i, j, k) + amp_cr*exp(-r2/r0**2)
                              elseif (icr == icr_C12) then
                                 cg%u(iarr_crn(icr), i, j, k) = cg%u(iarr_crn(icr), i, j, k) + amp_cr*0.1*exp(-r2/r0**2) ! BEWARE: magic number
                              else
                                 cg%u(iarr_crn(icr), i, j, k) = 0.0
                              endif

                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo

         cgl => cgl%nxt
      enddo

      if (is_multicg) call die("[initproblem:init_prob] multiple grid pieces per procesor not implemented yet") !nontrivial maxv
#endif /* COSM_RAYS */

   end subroutine init_prob
!-----------------------------------------------------------------------------
!
! This routine provides the "dvel"  variable values to be dumped to the .h5 file
! "dvel" is the velocity divergence
!
   subroutine dvel_var_hdf5(var, tab, ierrh, cg)

      use dataio_pub,  only: die
      use grid_cont,   only: grid_container
      use named_array, only: qna

      implicit none

      character(len=*), intent(in)                    :: var
      real(kind=4), dimension(:,:,:), intent(inout)   :: tab
      integer, intent(inout)                          :: ierrh
      type(grid_container), pointer, intent(in)       :: cg

      if (.not. qna%exists(dvel_n)) call die("[initproblem:dvel_vars] Cannot find dvel_n")

      ierrh = 0
      select case (trim(var))
         case (dvel_n)
            tab(:,:,:) = real(cg%q(qna%ind(dvel_n))%span(cg%ijkse), 4)
         case default
            ierrh = -1
      end select

   end subroutine dvel_var_hdf5
!-----------------------------------------------------------------------------
end module initproblem
