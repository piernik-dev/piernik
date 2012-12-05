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
#if defined GALAXY && defined SN_SRC
#define SN_GALAXY
#endif /* GALAXY && SN_SRC */
#if defined COSM_RAYS_SOURCES || defined SN_GALAXY
#define CRS_GALAXY
#endif /* COSM_RAYS_SOURCES || SN_GALAXY */
#if defined COSM_RAYS && defined SN_SRC
#define CR_SN
#endif /* COSM_RAYS && SN_SRC */

module initproblem

! Initial condition for the cosmic ray driven dynamo
! Based on Parker instability setup
! Written by: M. Hanasz, February 2006
! Modified by M.Hanasz for CR-driven dynamo

   use constants, only: ndims

   implicit none

   private
   public :: read_problem_par, init_prob, problem_pointers

   real :: d0, alpha, bxn, byn, bzn, amp_cr, beta_cr                         !< galactic disk specific parameters
   real :: x0, y0, z0                                                        !< parameters for a single supernova exploding at t=0
   real, dimension(ndims) :: b_n, sn_pos

   namelist /PROBLEM_CONTROL/  d0, bxn, byn, bzn, x0, y0, z0, alpha, amp_cr, beta_cr

contains

!-----------------------------------------------------------------------------

   subroutine problem_pointers

!      use gravity,    only: grav_accel
      use user_hooks, only: problem_customize_solution

      implicit none

!      grav_accel                 => galactic_grav_accel
      problem_customize_solution => supernovae_wrapper

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use dataio_pub, only: nh      ! QA_WARN required for diff_nml
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

      sn_pos = [x0,  y0,  z0 ]
      b_n    = [bxn, byn, bzn]

#ifdef GRAV
!      if (user_grav) grav_pot_3d => my_grav_pot_3d
      if (user_grav) grav_pot_3d => galactic_grav_pot_3d
#endif /* GRAV */

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine init_prob

      use cg_list,        only: cg_list_element
      use cg_leaves,      only: leaves
      use constants,      only: xdim, ydim, zdim, LO, HI
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
#ifdef CRS_GALAXY
      use cr_data,        only: eCRSP, cr_table
#endif /* CRS_GALAXY */
#ifdef SN_GALAXY
      use cr_data,        only: icr_H1, icr_C12
      use domain,         only: dom
      use snsources,      only: r_sn
#endif /* SN_GALAXY */
#endif /* COSM_RAYS */
#ifdef GRAV
      use gravity,        only: grav_pot_3d
#endif /* GRAV */
#ifdef COSM_RAYS_SOURCES
      use cr_data,        only: cr_sigma, icr_N14, icr_O16
#endif /* COSM_RAYS_SOURCES */
      implicit none

      class(component_fluid), pointer :: fl
      integer                         :: i, j, k
      real                            :: b0, csim2
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg
#ifdef SN_GALAXY
      real                            :: decr, x1, x2, y1, y2, z1
#endif /* SN_GALAXY */

#ifdef COSM_RAYS_SOURCES
! really workaround for the gold
      if (eCRSP(icr_N14)) cr_sigma(cr_table(icr_N14),:) = 0.0
      if (eCRSP(icr_O16)) cr_sigma(cr_table(icr_O16),:) = 0.0
#endif /* COSM_RAYS_SOURCES */

      call grav_pot_3d

!   Secondary parameters
      fl => flind%ion

      b0 = sqrt(2.*alpha*d0*fl%cs2)

      csim2 = fl%cs2*(1.0+alpha)

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         call set_default_hsparams(cg)
         i = cg%lhn(xdim,LO)
         j = cg%lhn(ydim,LO)
         call hydrostatic_zeq_densmid(i, j, d0, csim2)

         do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
            do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
               do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
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
#ifdef SN_GALAXY
! Single SN explosion in x0,y0,z0 at t = 0 if amp_cr /= 0
                  if (any([eCRSP(icr_H1), eCRSP(icr_C12)])) then
                     x1 = (cg%x(i)- x0)**2 ; x2 = (cg%x(i)-(x0+dom%L_(xdim)))**2
                     y1 = (cg%y(j)- y0)**2 ; y2 = (cg%y(j)-(y0+dom%L_(ydim)))**2
                     z1 = (cg%z(k)- z0)**2
                     decr = amp_cr * (exp(-(x1 + y1 + z1)/r_sn**2) + exp(-(x2 + y1 + z1)/r_sn**2) + exp(-(x1 + y2 + z1)/r_sn**2) + exp(-(x2 + y2 + z1)/r_sn**2))
                  endif
                  if (eCRSP(icr_H1 )) cg%u(iarr_crn(cr_table(icr_H1 )),i,j,k)= cg%u(iarr_crn(cr_table(icr_H1 )),i,j,k) +     decr
                  if (eCRSP(icr_C12)) cg%u(iarr_crn(cr_table(icr_C12)),i,j,k)= cg%u(iarr_crn(cr_table(icr_C12)),i,j,k) + 0.1*decr
#endif /* SN_GALAXY */
#endif /* COSM_RAYS */
               enddo
            enddo
         enddo

         cgl => cgl%nxt
      enddo

#ifdef CR_SN
      call cr_sn_beware(sn_pos)
#endif /* CR_SN */


      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
            do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
               do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
                  cg%b(:,i,j,k) = b0*sqrt(cg%u(fl%idn,i,j,k)/d0) * b_n/sqrt(sum(b_n**2))
#ifndef ISO
                  cg%u(fl%ien,i,j,k) = cg%u(fl%ien,i,j,k) + emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))
#endif /* !ISO */
               enddo
            enddo
         enddo

         cgl => cgl%nxt
      enddo

   end subroutine init_prob

   subroutine supernovae_wrapper(forward)
#ifdef SN_SRC
      use snsources, only: random_sn
#endif /* SN_SRC */
      implicit none

      logical, intent(in) :: forward
#ifdef SN_SRC
      if (forward) call random_sn
#endif /* SN_SRC */
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

   subroutine galactic_grav_pot_3d

      use axes_M,    only: axes
      use cg_leaves, only: leaves
      use cg_list,   only: cg_list_element
      use gravity,   only: grav_type
      use grid_cont, only: grid_container

      implicit none

      type(axes)                     :: ax
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         call ax%allocate_axes(cg%lhn)
         ax%x(:) = cg%x(:)
         ax%y(:) = cg%y(:)
         ax%z(:) = cg%z(:)

         call galactic_grav_pot(cg%gp, ax, cg%lhn)
         grav_type => galactic_grav_pot

         call ax%deallocate_axes

         cgl => cgl%nxt
      enddo

   end subroutine galactic_grav_pot_3d

   subroutine galactic_grav_pot(gp, ax, lhn, flatten)

      use axes_M,    only: axes
      use constants, only: ndims, LO, HI, zdim, half
      use gravity,   only: r_gc
      use units,     only: r_gc_sun, kpc

      implicit none

      real, dimension(:,:,:), pointer                     :: gp
      type(axes),                              intent(in) :: ax
      integer(kind=4), dimension(ndims,LO:HI), intent(in) :: lhn
      logical,                       optional, intent(in) :: flatten
      integer                                             :: k

      real, parameter :: f1 = 3.23e8, f2 = -4.4e-9, f3 = 1.7e-9
      real            :: r1, r22, r32, s4, s5

      r1  = 4.9*kpc
      r22 = (0.2*kpc)**2
      r32 = (2.2*kpc)**2
      s4  = f2 * exp(-(r_gc-r_gc_sun)/(r1))
      s5  = f3 * (r_gc_sun**2 + r32)/(r_gc**2 + r32)

!      grav = f1 * ((s4 * xsw/sqrt(xsw**2+r22)) - (s5 * xsw/kpc) )
!!          -Om*(Om+G) * Z * (kpc ?) ! in the transition region between rigid and flat rotation F'98: eq.(36)

      do k = lhn(zdim,LO), lhn(zdim,HI)
         gp(:,:,k) = -f1 * (s4 * sqrt(ax%z(k)**2+r22) - s5 * half * ax%z(k)**2 / kpc)
      enddo
      return

      if (.false. .and. present(flatten)) k = 0 ! suppress compiler warnings

   end subroutine galactic_grav_pot

!------------------------------------------------------------------------------
#ifdef CR_SN
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
      use cg_leaves,      only: leaves
      use constants,      only: xdim, ydim, zdim, ndims, LO, HI
      use domain,         only: dom
      use grid_cont,      only: grid_container
#ifdef COSM_RAYS_SOURCES
      use cr_data,        only: cr_table, cr_primary, eCRSP, icr_H1, icr_C12 !, icr_N14, icr_O16
      use initcosmicrays, only: iarr_crn
#endif /* COSM_RAYS_SOURCES */
      use snsources,      only: r_sn
#ifdef SHEAR
      use snsources,      only: sn_shear
#endif /* SHEAR */

      implicit none

      real, dimension(ndims), intent(in) :: pos
      integer                            :: i, j, k, ipm, jpm
      real                               :: decr, xsn, ysn, zsn, ysna, zr
      type(cg_list_element), pointer     :: cgl
      type(grid_container),  pointer     :: cg
#ifdef SHEAR
      real, dimension(3)                 :: ysnoi
#endif /* SHEAR */

      xsn = pos(xdim)
      ysn = pos(ydim)
      zsn = pos(zdim)

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

#ifdef SHEAR
         ysnoi(2) = ysn
         call sn_shear(cg, ysnoi)
#endif /* !SHEAR */

         do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
            zr = (cg%z(k)-zsn)**2
            do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
               do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)

                  decr = 0.0
                  do ipm=-1,1
#ifdef SHEAR
                     ysna = ysnoi(ipm+2)
#else /* !SHEAR */
                     ysna = ysn
#endif /* !SHEAR */
                     do jpm=-1,1

!                     decr = amp_ecr_sn * ethu  &
                        decr = decr + exp(-((cg%x(i)-xsn +real(ipm)*dom%L_(xdim))**2  &
                             +              (cg%y(j)-ysna+real(jpm)*dom%L_(ydim))**2 + zr)/r_sn**2)

                     enddo ! jpm
                  enddo ! ipm
                  decr = decr * amp_cr
#ifdef COSM_RAYS_SOURCES
!                     cg%u(iarr_crn,i,j,k) = cg%u(iarr_crn,i,j,k) + max(decr,1e-10) * [1., primary_C12*12., primary_N14*14., primary_O16*16.]
                  if (eCRSP(icr_H1 )) cg%u(iarr_crn(cr_table(icr_H1 )),i,j,k) = cg%u(iarr_crn(cr_table(icr_H1 )),i,j,k) + decr
                  if (eCRSP(icr_C12)) cg%u(iarr_crn(cr_table(icr_C12)),i,j,k) = cg%u(iarr_crn(cr_table(icr_C12)),i,j,k) + cr_primary(cr_table(icr_C12))*12*decr
                  !> \deprecated BEWARE: following lines are inconsistent with the gold for some reason
!                  if (eCRSP(icr_N14)) cg%u(iarr_crn(cr_table(icr_N14)),i,j,k) = cg%u(iarr_crn(cr_table(icr_N14)),i,j,k) + cr_primary(cr_table(icr_N14))*14*decr
!                  if (eCRSP(icr_O16)) cg%u(iarr_crn(cr_table(icr_O16)),i,j,k) = cg%u(iarr_crn(cr_table(icr_O16)),i,j,k) + cr_primary(cr_table(icr_O16))*16*decr
#endif /* COSM_RAYS_SOURCES */
               enddo ! i
            enddo ! j
         enddo ! k
         cgl => cgl%nxt
      enddo

   end subroutine cr_sn_beware
#endif /* CR_SN */
end module initproblem
