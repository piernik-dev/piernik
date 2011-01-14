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

!<
!! \brief Implementation of the Coriolis forces for rotating grid
!>

module coriolis
! pulled by CORIOLIS
! \todo Drop the CORIOLIS preprocessor symbol in favour of something more elegant if it does not incur too much overhead on rtvd.

   implicit none

   private
   public  :: init_coriolis, coriolis_force, set_omega, coriolis_omega

   real, protected :: coriolis_omega !< Rotational angular velocity

contains

!<
!! \brief Reset coriolis_omega and perform basic checks
!>

   subroutine init_coriolis

      use dataio_pub, only: die
      use grid,       only: geometry

      implicit none

      coriolis_omega = 0.

      if (geometry /= "cartesian") call die("[coriolis:init_coriolis] Only cartesian geometry is implemented")
#if !(defined GRAV || defined SHEAR || defined FLUID_INTERACTIONS)
      call die("coriolis:init_coriolis] Check how and under what conditions the rtvd::relaxing_tvd handles additional source terms")
#endif

   end subroutine init_coriolis

!<
!! \brief Compute the Coriolis acceleration for a given row of cells.
!!
!! \details This is a low-order estimate of the Coriolis accelerations, because this routine uses density and velocity fields
!! from the beginning of the time step. This is a simple approach, but ignores any changes due to other accelerations during the time step.
!>

   subroutine coriolis_force(sweep, u, rotacc)

      use fluidindex, only: nvar, iarr_all_dn, iarr_all_my

      implicit none

      character(len=*), intent(in)                   :: sweep  !< string of characters that points out the current sweep direction
      real, dimension(:,:), intent(in)               :: u
      real, dimension(nvar%fluids, size(u,2)), intent(inout) :: rotacc !< an array for Coriolis accelerations

      ! Coriolis force for corotating coords
      select case (sweep)
         case ('xsweep')
            rotacc(:,:) = rotacc(:,:) + 2.0 * coriolis_omega * u(iarr_all_my(:), :)/u(iarr_all_dn(:), :)
         case ('ysweep')
            rotacc(:,:) = rotacc(:,:) - 2.0 * coriolis_omega * u(iarr_all_my(:), :)/u(iarr_all_dn(:), :)
!         case ('zsweep') !no z-component of the Coriolis force
      end select

   end subroutine coriolis_force

!< \brief Provides a way to set the protected variable coriolis_omega

   subroutine set_omega(omega_in)

      implicit none

      real, intent(in) :: omega_in

      coriolis_omega = omega_in

   end subroutine set_omega

end module coriolis
