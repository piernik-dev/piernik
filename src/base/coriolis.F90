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

!>
!! \brief Implementation of the Coriolis forces for rotating grid
!! \todo Drop the CORIOLIS preprocessor symbol in favour of something more elegant if it does not incur too much overhead on rtvd.
!<

module coriolis
! pulled by CORIOLIS

   implicit none

   private
   public  :: init_coriolis, coriolis_force, set_omega, coriolis_omega

   real, protected :: coriolis_omega              !< Rotational angular velocity
   logical, save   :: omega_uninitialized = .true.

contains

!>
!! \brief Reset coriolis_omega and perform basic checks
!<

   subroutine init_coriolis

      use constants,  only: PIERNIK_INIT_GRID, GEO_XYZ
      use dataio_pub, only: die, code_progress
      use domain,     only: dom

      implicit none

      if (code_progress < PIERNIK_INIT_GRID) call die("[coriolis:init_coriolis] grid not initialized.") ! this is a weak check, the real dependency is init_geometry at the moment

      if (omega_uninitialized) coriolis_omega = 0.

      if (dom%geometry_type /= GEO_XYZ) call die("[coriolis:init_coriolis] Only cartesian geometry is implemented")
#if !(defined GRAV || defined SHEAR )
      call die("[coriolis:init_coriolis] Check how and under what conditions the rtvd::relaxing_tvd handles additional source terms")
#endif /* !(GRAV || SHEAR ) */

   end subroutine init_coriolis

!>
!! \brief Compute the Coriolis acceleration for a given row of cells.
!!
!! The equations for the non-inertial acceleration in the x and y directions are given by:
!! \f{equation}
!! acc_x = 2 * \Omega * v_y
!! acc_y = -2 * \Omega * v_x
!! \f}
!!
!! \details This is a low-order estimate of the Coriolis accelerations, because this routine uses density and velocity fields
!! from the beginning of the time step. This is a simple approach, but ignores any changes due to other accelerations during the time step.
!!
!! \todo add cylindrical geometry support (check carefully what exactly is u(iarr_all_my(:), :)/u(iarr_all_dn(:), :) and figure out whether any geometrical factors are needed)
!<

   function coriolis_force(sweep, u) result(rotacc)

      use fluidindex, only: flind, iarr_all_dn, iarr_all_my
      use constants,  only: xdim, ydim !, zdim

      implicit none

      integer(kind=4), intent(in)              :: sweep  !< string of characters that points out the current sweep direction
      real, dimension(:,:), intent(in)         :: u      !< current fluid state vector
      real, dimension(size(u,1), flind%fluids) :: rotacc !< an array for Coriolis accelerations

      ! Coriolis force for corotating coords
      select case (sweep)
         case (xdim)
            rotacc(:,:) = +2.0 * coriolis_omega * u(:, iarr_all_my)/u(:, iarr_all_dn)
         case (ydim)
            rotacc(:,:) = -2.0 * coriolis_omega * u(:, iarr_all_my)/u(:, iarr_all_dn)
!         case (zdim) !no z-component of the Coriolis force
         case default
            rotacc(:,:) = 0.0
      end select

   end function coriolis_force

!> \brief Provides a way to set the protected variable coriolis_omega

   subroutine set_omega(omega_in)

      use dataio_pub, only: warn, printinfo, msg
      use mpisetup,   only: master

      implicit none

      real, intent(in) :: omega_in

      if (.not. omega_uninitialized) then
         write(msg,"(2(a,g12.5))")"[coriolis:set_omega] coriolis_omega was already set to ",coriolis_omega,". Resetting to ",omega_in
         call warn(msg)
      else
         write(msg,"(a,g12.5)")"[coriolis:set_omega] Setting coriolis_omega to ",omega_in
         if (master) call printinfo(msg)
      endif
      coriolis_omega = omega_in
      omega_uninitialized = .false.

   end subroutine set_omega

end module coriolis
