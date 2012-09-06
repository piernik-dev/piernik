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

      use gravity,    only: grav_accel
      use user_hooks, only: problem_customize_solution

      implicit none

      grav_accel                 => galactic_grav_accel
      problem_customize_solution => supernovae_wrapper

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use dataio_pub, only: ierrh, par_file, namelist_errh, compare_namelist, cmdl_nml, lun      ! QA_WARN required for diff_nml
      use mpisetup,   only: rbuff, master, slave, piernik_MPI_Bcast
#ifdef GRAV
      use gravity,    only: grav_pot_3d, user_grav
#endif /* GRAV */

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

      if (user_grav) grav_pot_3d => my_grav_pot_3d

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine init_prob

      use cg_list,        only: cg_list_element
      use cg_list_bnd,    only: leaves
      use constants,      only: xdim, ydim, zdim
      use fluidindex,     only: flind
      use fluidtypes,     only: component_fluid
      use func,           only: ekin, emag
      use global,         only: smalld
      use grid_cont,      only: grid_container
      use hydrostatic,    only: hydrostatic_zeq_densmid, set_default_hsparams, dprof
#ifdef SHEAR
      use shear,          only: qshear, omega
#endif /* SHEAR */
#ifdef COSM_RAYS
      use initcosmicrays, only: gamma_crn, iarr_crn
#endif /* COSM_RAYS */
#ifdef GRAV
      use gravity,        only: grav_pot_3d
#endif /* GRAV */
      implicit none

      class(component_fluid), pointer :: fl
      integer                         :: i, j, k
      real                            :: b0, csim2
      real, dimension(3)              :: sn_pos
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg

      call grav_pot_3d

      sn_pos = [x0,y0,z0]

!   Secondary parameters
      fl => flind%ion

      b0 = sqrt(2.*alpha*d0*fl%cs2)

      csim2 = fl%cs2*(1.0+alpha)

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

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
                  cg%u(iarr_crn,i,j,k)  = 0.0
                  cg%u(iarr_crn(1),i,j,k) = beta_cr*fl%cs2 * cg%u(fl%idn,i,j,k)/( gamma_crn(1) - 1.0 )
!#ifdef GALAXY
!! Single SN explosion in x0,y0,z0 at t = 0 if amp_cr /= 0
!
!
!               cg%u(iarr_crn(icr_H1),i,j,k)= cg%u(iarr_crn(icr_H1),i,j,k) &
!                     + amp_cr*exp(-((x(i)- x0    )**2 + (y(j)- y0    )**2 + (z(k)-z0)**2)/r_sn**2) &
!                     + amp_cr*exp(-((x(i)-(x0+Lx))**2 + (y(j)- y0    )**2 + (z(k)-z0)**2)/r_sn**2) &
!                     + amp_cr*exp(-((x(i)- x0    )**2 + (y(j)-(y0+Ly))**2 + (z(k)-z0)**2)/r_sn**2) &
!                     + amp_cr*exp(-((x(i)-(x0+Lx))**2 + (y(j)-(y0+Ly))**2 + (z(k)-z0)**2)/r_sn**2)
!
!               cg%u(iarr_crn(icr_C12),i,j,k)= cg%u(iarr_crn(icr_C12),i,j,k) &
!                     + 0.1*amp_cr*exp(-((x(i)- x0    )**2 + (y(j)- y0    )**2 + (z(k)-z0)**2)/r_sn**2) &
!                     + 0.1*amp_cr*exp(-((x(i)-(x0+Lx))**2 + (y(j)- y0    )**2 + (z(k)-z0)**2)/r_sn**2) &
!                     + 0.1*amp_cr*exp(-((x(i)- x0    )**2 + (y(j)-(y0+Ly))**2 + (z(k)-z0)**2)/r_sn**2) &
!                     + 0.1*amp_cr*exp(-((x(i)-(x0+Lx))**2 + (y(j)-(y0+Ly))**2 + (z(k)-z0)**2)/r_sn**2)
!!
!
!#endif /* GALAXY */
#endif /* COSM_RAYS */
               enddo
            enddo
         enddo

         cgl => cgl%nxt
      enddo

#ifdef COSM_RAYS

      call cr_sn_beware(sn_pos)

#endif /* COSM_RAYS */


      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

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

   subroutine my_grav_pot_3d

      use dataio_pub, only: die, warn
      use gravity,    only: grav_accel, grav_accel2pot
      use mpisetup,   only: master

      implicit none

      logical, save  :: frun = .true.

      if (.not.frun) return

      if (associated(grav_accel)) then
         if (master) call warn("[initproblem:my_grav_pot_3d]: using 'grav_accel' defined by user")
         call grav_accel2pot
      else
         call die("[initproblem:my_grav_pot_3d]: GRAV is defined, but 'gp' is not initialized")
      endif
      frun = .false.

   end subroutine my_grav_pot_3d

!--------------------------------------------------------------------------
!>
!! \brief Routine that compute values of gravitational acceleration
!! \param sweep string of characters that points out the current sweep direction
!! \param i1 integer, number of column in the first direction after one pointed out by sweep
!! \param i2 integer, number of column in the second direction after one pointed out by sweep
!! \param xsw 1D position array in the direction pointed out by sweep
!! \param n number of elements of xsw array
!! \param grav 1D array of gravitational acceleration values computed for positions from xsw and returned by the routine
!! \n\n
!! one type of %gravity is implemented here: \n\n
!! local Galactic %gravity only in z-direction (see <a href="http://cdsads.u-strasbg.fr/abs/1998ApJ...497..759F">Ferriere K., 1998, Astrophys. Journal, 497, 759</a>)\n
!! \f[
!! F_z = 3.23 \cdot 10^8 \cdot \left[\left(-4.4 \cdot 10^{-9} \cdot exp\left(-\frac{(r_{gc}-r_{gc_{}Sun})}{(4.9kpc)}\right) \cdot \frac{z}{\sqrt{(z^2+(0.2kpc)^2)}}\right)
!! -\left( 1.7 \cdot 10^{-9} \cdot \frac{(r_{gc_{}Sun}^2 + (2.2kpc)^2)}{(r_{gc}^2 + (2.2kpc)^2)} \cdot \frac{z}{1kpc}\right) \right]
!! \f]
!! where \f$r_{gc}\f$ is galactocentric radius and \f$r_{gcSun}\f$ is the galactocentric radius of Sun.
!<

   subroutine galactic_grav_accel(sweep, i1,i2, xsw, n, grav)

      use constants, only: zdim
      use gravity,   only: r_gc
      use units,     only: r_gc_sun, kpc

      implicit none

      integer(kind=4),   intent(in)  :: sweep
      integer,           intent(in)  :: i1, i2
      integer(kind=4),   intent(in)  :: n
      real, dimension(n),intent(in)  :: xsw
      real, dimension(n),intent(out) :: grav

      if (.false.) grav(1) = i1+i2 ! suppress compiler warning on unused argument

      if (sweep == zdim) then
         grav = 3.23e8 * (  &
           (-4.4e-9 * exp(-(r_gc-r_gc_sun)/(4.9*kpc)) * xsw/sqrt(xsw**2+(0.2*kpc)**2)) &
           -( 1.7e-9 * (r_gc_sun**2 + (2.2*kpc)**2)/(r_gc**2 + (2.2*kpc)**2)*xsw/kpc) )
!          -Om*(Om+G) * Z * (kpc ?) ! in the transition region between rigid
!                                   ! and flat rotation F'98: eq.(36)
      else
         grav=0.0
      endif

   end subroutine galactic_grav_accel

!------------------------------------------------------------------------------
!BEWARE!
!>
!! \brief Routine that inserts an amount of cosmic ray energy around the position of supernova
!! \param pos real, dimension(3), array of supernova position components
!! \author M. Hanasz
!!
!! BEWARE: the code is very similar to snsources:cr_sn . Merge?
!!
!<
   subroutine cr_sn_beware(pos)

      use cg_list,        only: cg_list_element
      use cg_list_bnd,    only: leaves
      use constants,      only: xdim, ydim, zdim, ndims
      use cr_data,        only: icr_H1, icr_C12, icr_N14, icr_O16, primary_C12, primary_N14, primary_O16
      use domain,         only: dom
      use fluidindex,     only: flind
      use grid_cont,      only: grid_container
      use initcosmicrays, only: iarr_crn
      use snsources,      only: r_sn

      implicit none

      real, dimension(ndims), intent(in) :: pos

      integer                            :: i, j, k, ipm, jpm, icr
      real                               :: decr, xsn, ysn, zsn, ysna !, ysni, ysno
      type(cg_list_element), pointer     :: cgl
      type(grid_container),  pointer     :: cg

      xsn = pos(xdim)
      ysn = pos(ydim)
      zsn = pos(zdim)

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         do k=1, cg%n_(zdim)
            do j=1, cg%n_(ydim)
               do i=1, cg%n_(xdim)

                  do ipm=-1,1

                     ysna = ysn
! ToDo: when implementing SHEAR, use select case construct or an assignment like this: ysna = ysn_array(ipm)
!                  if (ipm == -1) ysna = ysn   ! if (SHEAR) => ysna = ysno
!                  if (ipm == 0) ysna = ysn
!                  if (ipm == 1) ysna = ysn   ! if (SHEAR) => ysna = ysni

                     do jpm=-1,1

!                     decr = amp_ecr_sn * ethu  &
                        decr = amp_cr  &
                             * exp(-((cg%x(i)-xsn+real(ipm)*dom%L_(xdim))**2  &
                             + (cg%y(j)-ysna+real(jpm)*dom%L_(ydim))**2  &
                             + (cg%z(k)-zsn)**2)/r_sn**2)
!                     cg%u(iarr_crn,i,j,k) = cg%u(iarr_crn,i,j,k) + max(decr,1e-10) * [1., primary_C12*12., primary_N14*14., primary_O16*16.]
                        do icr = 1, flind%crn%all
                           select case (icr)
                              case (icr_H1)
                                 cg%u(iarr_crn(icr),i,j,k) = cg%u(iarr_crn(icr),i,j,k) + decr
                              case (icr_C12)
                                 cg%u(iarr_crn(icr),i,j,k) = cg%u(iarr_crn(icr),i,j,k) + primary_C12*12*decr
                              case (icr_N14)
                                 cg%u(iarr_crn(icr),i,j,k) = cg%u(iarr_crn(icr),i,j,k) + primary_N14*14*decr
                              case (icr_O16)
                                 cg%u(iarr_crn(icr),i,j,k) = cg%u(iarr_crn(icr),i,j,k) + primary_O16*16*decr
                           end select
                        enddo

                     enddo ! jpm
                  enddo ! ipm

               enddo ! i
            enddo ! j
         enddo ! k
         cgl => cgl%nxt
      enddo

   end subroutine cr_sn_beware

end module initproblem
