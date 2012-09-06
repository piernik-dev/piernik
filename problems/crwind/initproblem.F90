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

! Initial condition for the cosmic ray driven dynamo
! Based on Parker instability setup
! Written by: M. Hanasz, February 2006
! Modified by M.Hanasz for CR-driven dynamo

   implicit none

   private
   public :: read_problem_par, init_prob, problem_pointers

   real :: d0, alpha, bxn,byn,bzn, amp_cr, beta_cr                           !< galactic disk specific parameters
   real :: x0, y0, z0                                                        !< parameters for a single supernova exploding at t=0

   namelist /PROBLEM_CONTROL/  d0, bxn, byn, bzn, x0, y0, z0, alpha, amp_cr, beta_cr

contains

!-----------------------------------------------------------------------------

   subroutine problem_pointers

      use user_hooks, only: problem_customize_solution

      implicit none

      problem_customize_solution => supernovae_wrapper

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use dataio_pub, only: ierrh, par_file, namelist_errh, compare_namelist, cmdl_nml, lun      ! QA_WARN required for diff_nml
      use mpisetup,   only: rbuff, master, slave, piernik_MPI_Bcast

      implicit none

      d0     = 1.0
      bxn    = 0.0
      byn    = 1.0
      bzn    = 0.0
      x0     = 0.0
      y0     = 0.0
      z0     = 0.0

      if (master) then

         diff_nml(PROBLEM_CONTROL)

         rbuff(1)  = d0
         rbuff(2)  = bxn
         rbuff(3)  = byn
         rbuff(4)  = bzn
         rbuff(5)  = x0
         rbuff(6)  = y0
         rbuff(7)  = z0
         rbuff(8)  = amp_cr
         rbuff(9)  = beta_cr
         rbuff(10) = alpha

      endif

      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         d0        = rbuff(1)
         bxn       = rbuff(2)
         byn       = rbuff(3)
         bzn       = rbuff(4)
         x0        = rbuff(5)
         y0        = rbuff(6)
         z0        = rbuff(7)
         amp_cr    = rbuff(8)
         beta_cr   = rbuff(9)
         alpha     = rbuff(10)

      endif

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine init_prob

      use cg_list,        only: cg_list_element
      use cg_leaves,      only: leaves
      use constants,      only: xdim, ydim, zdim
      use domain,         only: dom
      use fluidindex,     only: flind
      use fluidtypes,     only: component_fluid
      use func,           only: ekin, emag
      use global,         only: smalld
      use grid_cont,      only: grid_container
      use hydrostatic,    only: hydrostatic_zeq_densmid, set_default_hsparams, dprof
      use snsources,      only: r_sn
#ifdef SHEAR
      use shear,          only: qshear, omega
#endif /* SHEAR */
#ifdef COSM_RAYS
      use initcosmicrays, only: gamma_crs, iarr_crs
#endif /* COSM_RAYS */

      implicit none

      class(component_fluid), pointer :: fl
      integer                         :: i, j, k
      real                            :: b0, csim2
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg

!   Secondary parameters
      fl => flind%ion

      b0 = sqrt(2.*alpha*d0*fl%cs2)

      csim2 = fl%cs2*(1.0+alpha)

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         if (associated(cg%cs_iso2)) cg%cs_iso2(:,:,:) = fl%cs2

         call set_default_hsparams(cg)
         call hydrostatic_zeq_densmid(1, 1, d0, csim2)

         do k = 1, cg%n_(zdim)
            do j = 1, cg%n_(ydim)
               do i = 1, cg%n_(xdim)
                  cg%u(fl%idn,i,j,k)   = max(smalld, dprof(k))

                  cg%u(fl%imx,i,j,k) = 0.0
                  cg%u(fl%imy,i,j,k) = 0.0
                  cg%u(fl%imz,i,j,k) = 0.0
#ifdef SHEAR
                  cg%u(fl%imy,i,j,k) = -qshear*omega*cg%x(i)*cg%u(fl%idn,i,j,k)
#endif /* SHEAR */

#ifndef ISO
                  cg%u(fl%ien,i,j,k) = fl%cs2/(fl%gam_1) * cg%u(fl%idn,i,j,k) + ekin(cg%u(fl%imx,i,j,k), cg%u(fl%imy,i,j,k), cg%u(fl%imz,i,j,k), cg%u(fl%idn,i,j,k))
#endif /* !ISO */
#ifdef COSM_RAYS
                  cg%u(iarr_crs,i,j,k) =  beta_cr*fl%cs2 * cg%u(fl%idn,i,j,k)/( gamma_crs - 1.0 )
#ifdef GALAXY
! Single SN explosion in x0,y0,z0 at t = 0 if amp_cr /= 0
                  cg%u(iarr_crs,i,j,k)= cg%u(iarr_crs,i,j,k) &
                       + amp_cr*exp(-((cg%x(i)- x0              )**2 + (cg%y(j)- y0              )**2 + (cg%z(k)-z0)**2)/r_sn**2) &
                       + amp_cr*exp(-((cg%x(i)-(x0+dom%L_(xdim)))**2 + (cg%y(j)- y0              )**2 + (cg%z(k)-z0)**2)/r_sn**2) &
                       + amp_cr*exp(-((cg%x(i)- x0              )**2 + (cg%y(j)-(y0+dom%L_(ydim)))**2 + (cg%z(k)-z0)**2)/r_sn**2) &
                       + amp_cr*exp(-((cg%x(i)-(x0+dom%L_(xdim)))**2 + (cg%y(j)-(y0+dom%L_(ydim)))**2 + (cg%z(k)-z0)**2)/r_sn**2)
#endif /* GALAXY */
#endif /* COSM_RAYS */
               enddo
            enddo
         enddo

         do k = 1, cg%n_(zdim)
            do j = 1, cg%n_(ydim)
               do i = 1, cg%n_(xdim)
                  cg%b(xdim,i,j,k)   = b0*sqrt(cg%u(fl%idn,i,j,k)/d0)* bxn/sqrt(bxn**2+byn**2+bzn**2)
                  cg%b(ydim,i,j,k)   = b0*sqrt(cg%u(fl%idn,i,j,k)/d0)* byn/sqrt(bxn**2+byn**2+bzn**2)
                  cg%b(zdim,i,j,k)   = b0*sqrt(cg%u(fl%idn,i,j,k)/d0)* bzn/sqrt(bxn**2+byn**2+bzn**2)
#ifndef ISO
                  cg%u(fl%ien,i,j,k)   = cg%u(fl%ien,i,j,k) + emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))
#endif /* !ISO */
               enddo
            enddo
         enddo

         cgl => cgl%nxt
      enddo

   end subroutine init_prob

   subroutine supernovae_wrapper(forward)

      use snsources, only: random_sn

      implicit none

      logical, intent(in) :: forward

      if (forward) call random_sn

   end subroutine supernovae_wrapper

end module initproblem
