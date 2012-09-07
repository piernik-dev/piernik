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

   use constants, only: ndims

   implicit none

   private
   public :: read_problem_par, init_prob, problem_pointers

   integer(kind=4) :: norm_step
   real            :: t_sn
   integer         :: n_sn
   real            :: d0, p0, r0, beta_cr, amp_cr
   real, dimension(ndims) :: c_exp, c_rot, b0, sn_pos


   namelist /PROBLEM_CONTROL/  d0, p0, b0, sn_pos, r0, beta_cr, amp_cr, norm_step, c_exp, c_rot

contains
!-----------------------------------------------------------------------------
   subroutine problem_pointers

      use dataio_user, only: user_vars_hdf5

      implicit none

      user_vars_hdf5             => dvel_var_hdf5

   end subroutine problem_pointers
!-----------------------------------------------------------------------------
   subroutine read_problem_par

      use dataio_pub, only: ierrh, par_file, namelist_errh, compare_namelist, cmdl_nml, lun ! QA_WARN required for diff_nml
      use dataio_pub, only: die
      use domain,     only: dom
      use mpisetup,   only: ibuff, rbuff, master, slave, piernik_MPI_Bcast

      implicit none

      t_sn = 0.0

      d0           = 1.0e5     !< density
      p0           = 1.0       !< pressure
      b0           = [ 0., 0., 0. ] !< Magnetic field
      sn_pos       = [ 0., 0., 0. ] !< position of blob
      r0           = 5.* minval(dom%L_(:)/dom%n_d(:), mask=dom%has_dir(:))  !< radius of the blob

      beta_cr      = 0.0       !< ambient level
      amp_cr       = 1.0       !< amplitude of the blob

      norm_step    = 10        !< how often to compute the norm (in steps)

      c_exp = [ 0.0, 0.0, 0.0 ]
      c_rot = [ 0.0, 0.0, 0.0 ]

      if (master) then

         diff_nml(PROBLEM_CONTROL)

         rbuff(1)     = d0
         rbuff(2)     = p0
         rbuff(3:5)   = b0
         rbuff(6:8)   = sn_pos
         rbuff(9)     = r0
         rbuff(10)    = beta_cr
         rbuff(11)    = amp_cr
         rbuff(12:14) = c_exp
         rbuff(15:17) = c_rot

         ibuff(1) = norm_step
      endif

      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         d0           = rbuff(1)
         p0           = rbuff(2)
         b0           = rbuff(3:5)
         sn_pos       = rbuff(6:8)
         r0           = rbuff(9)
         beta_cr      = rbuff(10)
         amp_cr       = rbuff(11)
         c_exp        = rbuff(12:14)
         c_rot        = rbuff(15:17)

         norm_step = int(ibuff(1), kind=4)
      endif

      if (r0 == 0.0) call die("[initproblem:read_problem_par] r0 == 0")

   end subroutine read_problem_par
!-----------------------------------------------------------------------------
   subroutine init_prob

      use cg_leaves,      only: leaves
      use cg_list,        only: cg_list_element
      use constants,      only: xdim, ydim, zdim, pi, I_ONE
      use crhelpers,      only: div_v
      use dataio_pub,     only: warn
      use domain,         only: dom
      use fluidindex,     only: flind
      use fluidtypes,     only: component_fluid
      use func,           only: emag, ekin
      use grid_cont,      only: grid_container
      use initcosmicrays, only: iarr_crn, iarr_crs, gamma_crn, K_crn_paral, K_crn_perp
#ifdef COSM_RAYS_SOURCES
      use cr_data,        only: icr_H1, icr_C12, cr_table
#endif /* COSM_RAYS_SOURCES */

      implicit none

      class(component_fluid), pointer  :: fl
      integer                          :: i, j, k, icr, ipm, jpm, kpm
      real                             :: cs_iso, r, r2
      type(cg_list_element),  pointer  :: cgl
      type(grid_container),   pointer  :: cg
#ifndef COSM_RAYS_SOURCES
      integer, parameter               :: icr_H1 = 1, icr_C12 = 2
      integer, parameter, dimension(2) :: cr_table = [1,2]
#endif /* !COSM_RAYS_SOURCES */

      fl => flind%ion

! Uniform equilibrium state

      cs_iso = sqrt(p0/d0)

      where (.not. dom%has_dir)
         b0 = 0.0  ! ignore B field in nonexistent direction
      endwhere

      if (sum(b0**2) == 0.0 .and. (any(K_crn_paral(:) /= 0.) .or. any(K_crn_perp(:) /= 0.))) then
         call warn("[initproblem:init_prob] No magnetic field is set, K_crn_* also have to be 0.")
         K_crn_paral(:) = 0.
         K_crn_perp(:)  = 0.
      endif

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         cg%b(xdim, :, :, :) = b0(xdim)
         cg%b(ydim, :, :, :) = b0(ydim)
         cg%b(zdim, :, :, :) = b0(zdim)
         cg%u(fl%idn, :, :, :) = d0
         cg%u(fl%imx:fl%imz, :, :, :) = 0.0

         do k = lbound(cg%u, zdim+I_ONE), ubound(cg%u, zdim+I_ONE)
            do j = lbound(cg%u, ydim+I_ONE), ubound(cg%u, ydim+I_ONE)
               do i = lbound(cg%u, zdim+I_ONE), ubound(cg%u, xdim+I_ONE)

                  cg%u(fl%imx,i,j,k) = c_exp(xdim) * cg%u(fl%idn,i,j,k) * cg%x(i)
                  cg%u(fl%imy,i,j,k) = c_exp(ydim) * cg%u(fl%idn,i,j,k) * cg%y(j)
                  cg%u(fl%imz,i,j,k) = c_exp(zdim) * cg%u(fl%idn,i,j,k) * cg%z(k)

                  r = sqrt( (cg%x(i))**2 + (cg%y(j))**2 + (cg%z(k))**2 )

                  cg%u(fl%imx,i,j,k) = cg%u(fl%imx,i,j,k) + c_rot(zdim) * cg%u(fl%idn,i,j,k) * (-cg%y(j)/r) * sin(pi*r)
                  cg%u(fl%imy,i,j,k) = cg%u(fl%imy,i,j,k) + c_rot(zdim) * cg%u(fl%idn,i,j,k) * ( cg%x(i)/r) * sin(pi*r)
                  cg%u(fl%imz,i,j,k) = cg%u(fl%imz,i,j,k) + 0.0
#ifndef ISO
                  cg%u(fl%ien,i,j,k) = p0 / fl%gam_1 + &
                       &               ekin(cg%u(fl%imx,i,j,k), cg%u(fl%imy,i,j,k), cg%u(fl%imz,i,j,k), cg%u(fl%idn,i,j,k)) + &
                       &               emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))
#endif /* !ISO */
               enddo
            enddo
         enddo

         call div_v(flind%ion%pos, cg)

#ifdef COSM_RAYS
         do icr = lbound(iarr_crs, 1), ubound(iarr_crs, 1)
            cg%u(iarr_crs(icr), :, :, :) =  beta_cr*fl%cs2 * cg%u(fl%idn, :, :, :) / (gamma_crn(icr)-1.0)
         enddo

! Explosions
         do icr = 1, flind%crn%all
            do k = cg%ks, cg%ke
               do j = cg%js, cg%je
                  do i = cg%is, cg%ie
                     do ipm = -1, 1
                        do jpm = -1, 1
                           do kpm = -1, 1

                              r2 = (cg%x(i) - sn_pos(xdim) + real(ipm) * dom%L_(xdim))**2 + &
                                 & (cg%y(j) - sn_pos(ydim) + real(jpm) * dom%L_(ydim))**2 + &
                                 & (cg%z(k) - sn_pos(zdim) + real(kpm) * dom%L_(zdim))**2
                              if (icr == cr_table(icr_H1)) then
                                 cg%u(iarr_crn(icr), i, j, k) = cg%u(iarr_crn(icr), i, j, k) + amp_cr*exp(-r2/r0**2)
                              elseif (icr == cr_table(icr_C12)) then
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

#endif /* COSM_RAYS */

   end subroutine init_prob
!-----------------------------------------------------------------------------
!
! This routine provides the "dvel"  variable values to be dumped to the .h5 file
! "dvel" is the velocity divergence
!
   subroutine dvel_var_hdf5(var, tab, ierrh, cg)

      use constants,   only: dsetnamelen
      use crhelpers,   only: divv_n
      use grid_cont,   only: grid_container
      use named_array_list, only: qna

      implicit none

      character(len=*), intent(in)                    :: var
      real(kind=4), dimension(:,:,:), intent(inout)   :: tab
      integer, intent(inout)                          :: ierrh
      type(grid_container), pointer, intent(in)       :: cg
      character(len=dsetnamelen), parameter :: dvel_n = "dvel"

      ierrh = 0
      select case (trim(var))
         case (dvel_n, divv_n)
            tab(:,:,:) = real(cg%q(qna%ind(divv_n))%span(cg%ijkse), 4)
         case default
            ierrh = -1
      end select

   end subroutine dvel_var_hdf5
!-----------------------------------------------------------------------------
end module initproblem
