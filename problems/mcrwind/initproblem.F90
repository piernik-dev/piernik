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

module initproblem

! Initial condition for the cosmic ray driven dynamo
! Based on Parker instability setup
! Written by: M. Hanasz, February 2006
! Modified by M.Hanasz for CR-driven dynamo

   use constants, only: ndims

   implicit none

   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   real :: d0, alpha, bxn, byn, bzn, amp_cr, beta_cr                         !< galactic disk specific parameters
   real :: x0, y0, z0                                                        !< parameters for a single supernova exploding at t=0
   real, dimension(ndims) :: b_n, sn_pos


   namelist /PROBLEM_CONTROL/  d0, bxn, byn, bzn, x0, y0, z0, alpha, amp_cr, beta_cr

contains

!-----------------------------------------------------------------------------

   subroutine problem_pointers

#ifdef GRAV
      use gravity,    only: grav_pot_3d
#endif /* GRAV */
      use user_hooks, only: problem_customize_solution

      implicit none

#ifdef GRAV
      grav_pot_3d => galactic_grav_pot_3d
#endif /* GRAV */
      problem_customize_solution => supernovae_wrapper

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use dataio_pub, only: nh
      use mpisetup,   only: rbuff, master, slave, piernik_MPI_Bcast
#if defined(COSM_RAYS) && defined(SN_SRC)
      use snsources,  only: amp_ecr_sn
#endif /* COSM_RAYS && SN_SRC */

      implicit none

      d0      = 1.0
      bxn     = 0.0
      byn     = 1.0
      bzn     = 0.0
      x0      = 0.0
      y0      = 0.0
      z0      = 0.0
      alpha   = 0.0
      amp_cr  = 0.0
      beta_cr = 0.0

      if (master) then

         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=PROBLEM_CONTROL)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=""
         read(unit=nh%lun, nml=PROBLEM_CONTROL, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "PROBLEM_CONTROL")
         read(nh%cmdl_nml,nml=PROBLEM_CONTROL, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "PROBLEM_CONTROL", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=PROBLEM_CONTROL)
         close(nh%lun)
         call nh%compare_namelist()

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

#if defined(COSM_RAYS) && defined(SN_SRC)
      if (amp_cr < 0.) amp_cr = amp_ecr_sn
#endif /* COSM_RAYS && SN_SRC */

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine problem_initial_conditions

      use cg_leaves,      only: leaves
      use cg_list,        only: cg_list_element
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
      use initcosmicrays, only: gamma_cr_1, iarr_crn, iarr_crs, nspc
#ifdef SN_SRC
      use snsources,      only: cr_sn
#endif /* SN_SRC */
#endif /* COSM_RAYS */

      implicit none


      !print *, icr_prim
      !print *, icr_sec
      !stop
        !print *, 'nspc : ', nspc
        !print *, 'icr_H1 : ', icr_H1, ' ',eH1(PRIM)
        !print *, 'icr_C12 : ', icr_C12, ' ',eC12(PRIM)
        !print *, 'icr_N14 : ', icr_N14, ' ',eN14(PRIM)
        !print *, 'icr_O16 : ', icr_O16, ' ',eO16(PRIM)

      class(component_fluid), pointer :: fl
      integer                         :: i, j, k
      real                            :: b0, csim2
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg

!   Secondary parameters
      fl => flind%ion

      b0 = sqrt(2. * alpha * d0 * fl%cs2)
      csim2 = fl%cs2 * (1.0 + alpha)

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         call set_default_hsparams(cg)
         i = cg%lhn(xdim,LO)
         j = cg%lhn(ydim,LO)
         call hydrostatic_zeq_densmid(i, j, d0, csim2)

         cg%u(fl%imx, RNG) = 0.0
         cg%u(fl%imy, RNG) = 0.0
         cg%u(fl%imz, RNG) = 0.0
#ifdef COSM_RAYS
         cg%u(iarr_crs, RNG)  = 0.0
#endif /* COSM_RAYS */

         do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
            cg%u(fl%idn,:,:,k) = max(smalld, dprof(k))
            do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
               do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)

#ifdef SHEAR
                  cg%u(fl%imy,i,j,k) = -qshear * omega * cg%x(i) * cg%u(fl%idn,i,j,k)
#endif /* SHEAR */
                  cg%b(:,i,j,k) = b0 * sqrt(cg%u(fl%idn,i,j,k) / d0) * b_n / sqrt(sum(b_n**2))
#ifndef ISO
                  cg%u(fl%ien,i,j,k) = fl%cs2 / fl%gam_1 * cg%u(fl%idn,i,j,k) + ekin(cg%u(fl%imx,i,j,k), cg%u(fl%imy,i,j,k), cg%u(fl%imz,i,j,k), cg%u(fl%idn,i,j,k)) + &
                                     & emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))
#endif /* !ISO */
#ifdef COSM_RAYS
#ifndef CRESP
                  cg%u(iarr_crn(1),i,j,k) = beta_cr * fl%cs2 * cg%u(fl%idn,i,j,k) / gamma_cr_1
#endif /* CRESP */
#endif /* COSM_RAYS */
               enddo
            enddo
         enddo

         cgl => cgl%nxt
      enddo

#if defined(COSM_RAYS) && defined(SN_SRC)
      call cr_sn(sn_pos, amp_cr)
#endif /* COSM_RAYS && SN_SRC */

   end subroutine problem_initial_conditions

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

#ifdef GRAV
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
!            -Om*(Om+G) * Z * (kpc ?) ! in the transition region between rigid
!                                    ! and flat rotation F'98: eq.(36)
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

      grav_type => galactic_grav_pot

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         if (.not. cg%is_old) then
            call ax%allocate_axes(cg%lhn)
            ax%x(:) = cg%x(:)
            ax%y(:) = cg%y(:)
            ax%z(:) = cg%z(:)

            call galactic_grav_pot(cg%gp, ax, cg%lhn)

            call ax%deallocate_axes
         endif

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
#endif /* GRAV */

end module initproblem
