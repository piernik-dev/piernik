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

!>
!! \brief Initialization of the dust fluid
!!
!!
!!
!! In this module following namelist of parameters is specified:
!! \copydetails initdust::init_dust
!<

module initdust
! pulled by ANY
   use fluidtypes, only: component_fluid
   implicit none

   private
   public :: init_dust, dust_fluid, cleanup_dust

   logical :: selfgrav    !< true if dust is selfgravitating

   type, extends(component_fluid) :: dust_fluid
      contains
         procedure, nopass :: get_tag
         procedure, pass   :: get_cs => dust_cs
         procedure, pass   :: compute_flux => flux_dust
         procedure, pass   :: initialize_indices => initialize_dust_indices
   end type dust_fluid

contains

   subroutine initialize_dust_indices(this, flind)
      use constants,  only: DST
      use fluidtypes, only: var_numbers
      implicit none
      class(dust_fluid), intent(inout) :: this
      type(var_numbers), intent(inout) :: flind

      call this%set_fluid_index(flind, .false., selfgrav, .false., 0.0, -1.0, DST)

   end subroutine initialize_dust_indices

   real function dust_cs(this, i, j, k, u, b, cs_iso2)
      implicit none
      class(dust_fluid),                 intent(in) :: this
      integer,                           intent(in) :: i, j, k
      real, dimension(:,:,:,:), pointer, intent(in) :: u       !< pointer to array of fluid properties
      real, dimension(:,:,:,:), pointer, intent(in) :: b       !< pointer to array of magnetic fields (used for ionized fluid with MAGNETIC #defined)
      real, dimension(:,:,:),   pointer, intent(in) :: cs_iso2 !< pointer to array of isothermal sound speeds (used when ISO was #defined)
      dust_cs = 0.0
      if (.false.) print *, u(:, i, j, k), b(:, i, j, k), cs_iso2(i, j, k), this%cs
   end function dust_cs

   function get_tag() result(tag)
      use constants, only: idlen
      implicit none
      character(len=idlen) :: tag

      tag = "DST"
   end function get_tag

!>
!! \brief Routine to set parameters from namelist FLUID_DUST
!!
!! \n \n
!! @b FLUID_DUST
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>selfgrav  </td><td>.false.</td><td>logical   </td><td>\copydoc initdust::selfgrav  </td></tr>
!! </table>
!! The list is active while \b "DUST" is defined.
!! \n \n
!<
   subroutine init_dust

      use dataio_pub,     only: nh  ! QA_WARN required for diff_nml
      use mpisetup,       only: lbuff, master, slave, piernik_MPI_Bcast

      implicit none

      namelist /FLUID_DUST/ selfgrav

      selfgrav = .false.

      if (master) then

         diff_nml(FLUID_DUST)

         lbuff(1)   = selfgrav

      endif

      call piernik_MPI_Bcast(lbuff)

      if (slave) then

         selfgrav    = lbuff(1)

      endif

   end subroutine init_dust

   subroutine cleanup_dust

      implicit none

   end subroutine cleanup_dust

#define RNG 2:nm
!/*
!>
!! \brief Computation of %fluxes for the dust fluid
!!
!!The flux functions for dust are given by
!!
!!\f[
!!  \vec{F}{(\vec{u})} =
!!  \left(\begin{array}{c}
!!    \rho v_x \\
!!    \rho v_x^2 \\
!!    \rho v_x v_y \\
!!    \rho v_x v_z
!!  \end{array}\right),
!!  \qquad
!!  \vec{G}{(\vec{u})} =
!!  \left(\begin{array}{c}
!!    \rho v_y \\
!!    \rho v_y v_x \\
!!    \rho v_y^2 \\
!!    \rho v_y v_z
!!  \end{array}\right),
!!\qquad
!!  \vec{H}{(\vec{u})} =
!!  \left(\begin{array}{c}
!!    \rho v_z \\
!!    \rho v_z v_x\\
!!    \rho v_z v_y \\
!!    \rho v_z^2
!!  \end{array}\right),
!!\f]
!!
!<
!*/
   subroutine flux_dust(this, flux, cfr, uu, n, vx, ps, bb, cs_iso2)

      use constants, only: idn, imx, imy, imz
#ifdef GLOBAL_FR_SPEED
      use timestep_pub, only: c_all
#endif /* GLOBAL_FR_SPEED */

      implicit none
      class(dust_fluid), intent(in)                :: this
      integer(kind=4), intent(in)                  :: n         !< number of cells in the current sweep
      real, dimension(:,:), intent(inout), pointer :: flux     !< flux of dust
      real, dimension(:,:), intent(inout), pointer :: cfr      !< freezing speed for dust
      real, dimension(:,:), intent(in),    pointer :: uu       !< part of u for dust
      real, dimension(:),   intent(inout), pointer :: vx        !< velocity of dust fluid for current sweep
      real, dimension(:),   intent(inout), pointer :: ps        !< pressure of dust fluid for current sweep
      real, dimension(:,:), intent(in),    pointer :: bb        !< magnetic field x,y,z-components table
      real, dimension(:),   intent(in),    pointer :: cs_iso2   !< local isothermal sound speed squared (optional)

      ! locals
!      real               :: minvx, maxvx, amp
#ifdef LOCAL_FR_SPEED
      real, dimension(size(vx)) :: absvx
#endif /* LOCAL_FR_SPEED */
      integer                   :: nm

      nm = n-1

      ps(:)  = 0.0
      vx(RNG)=uu(RNG, imx)/uu(RNG, idn) ; vx(1) = vx(2); vx(n) = vx(nm)

      flux(RNG, idn)=uu(RNG, imx)
      flux(RNG, imx)=uu(RNG, imx)*vx(RNG)
      flux(RNG, imy)=uu(RNG, imy)*vx(RNG)
      flux(RNG, imz)=uu(RNG, imz)*vx(RNG)

      flux(1, :) = flux(2, :); flux(n, :) = flux(nm, :)

#ifdef LOCAL_FR_SPEED

      ! The freezing speed is now computed locally (in each cell)
      !  as in Trac & Pen (2003). This ensures much sharper shocks,
      !  but sometimes may lead to numerical instabilities
!     minvx = minval(vx(RNG))
!     maxvx = maxval(vx(RNG))
!     amp   = (maxvx-minvx)*half
!      cfr(RNG, 1) = max(sqrt(vx(RNG)**2+cfr_smooth*amp),small)
      absvx = abs(vx)
      cfr(RNG, idn) = max( absvx(1:n-2), absvx(2:nm), absvx(3:n) )

      cfr(1, idn) = cfr(2, idn); cfr(n, idn) = cfr(nm, idn)

      cfr(:, imx) = cfr(:, idn)
      cfr(:, imy) = cfr(:, idn)
      cfr(:, imz) = cfr(:, idn)

#endif /* LOCAL_FR_SPEED */

#ifdef GLOBAL_FR_SPEED
      ! The freezing speed is now computed globally
      !  (c=const for the whole domain) in subroutine 'timestep'

      !  cfr(:,:) = flind%dst%c
      cfr(:,:) = c_all
#endif /* GLOBAL_FR_SPEED */
      return
      if (.false.) write(0,*) bb, cs_iso2, this%all

   end subroutine flux_dust
end module initdust
