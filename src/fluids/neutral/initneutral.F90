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
!! \brief Initialization of the neutral fluid
!!
!!
!!
!! In this module following namelist of parameters is specified:
!! \copydetails initneutral::init_neutral
!! \deprecated This module should not export any variables
!<

module initneutral
! pulled by ANY
   use fluidtypes, only: component_fluid
   implicit none

   private
   public :: init_neutral, cleanup_neutral, neutral_fluid

   real    :: gamma             !< adiabatic index for the neutral gas component
   real    :: cs_iso            !< isothermal sound speed (p = cs_iso<sup>2</sup>\f$\rho\f$), active only if neutral gas is \ref isothermal
   logical :: selfgrav          !< true if neutral gas is selfgravitating

   type, extends(component_fluid) :: neutral_fluid
      contains
         procedure, nopass :: get_tag
         procedure, pass   :: get_cs => neu_cs
         procedure, pass   :: compute_flux => flux_neu
         procedure, pass   :: initialize_indices => initialize_neu_indices
   end type neutral_fluid

contains

   subroutine initialize_neu_indices(this, flind)
      use constants,  only: NEU
      use fluidtypes, only: var_numbers
      implicit none
      class(neutral_fluid), intent(inout) :: this
      type(var_numbers),    intent(inout) :: flind

      logical :: has_energy
#ifdef ISO
      has_energy = .false.
#else /* !ISO */
      has_energy = .true.
#endif /* !ISO */

      call this%set_fluid_index(flind, .false., selfgrav, has_energy, cs_iso, gamma, NEU)

   end subroutine initialize_neu_indices

   real function neu_cs(this, i, j, k, u, b, cs_iso2)
#ifndef ISO
      use func,      only: ekin
#endif /* !ISO */
      implicit none
      class(neutral_fluid),              intent(in) :: this
      integer,                           intent(in) :: i, j, k
      real, dimension(:,:,:,:), pointer, intent(in) :: u       !< pointer to array of fluid properties
      real, dimension(:,:,:,:), pointer, intent(in) :: b       !< pointer to array of magnetic fields (used for ionized fluid with MAGNETIC #defined)
      real, dimension(:,:,:),   pointer, intent(in) :: cs_iso2 !< pointer to array of isothermal sound speeds (used when ISO was #defined)

      real :: p
#ifdef ISO
      p  = cs_iso2(i, j, k) * u(this%idn, i, j, k)
      neu_cs = sqrt(cs_iso2(i, j, k))
#else /* !ISO */
      p  = (u(this%ien, i, j, k) - &
         &   ekin(u(this%imx, i, j, k), u(this%imy, i, j, k), u(this%imz, i, j, k), u(this%idn, i, j, k)) &
         & ) * this%gam_1
      neu_cs = sqrt(abs((this%gam * p) / u(this%idn, i, j, k)))
#endif /* !ISO */
      if (.false.) print *, u(:, i, j, k), b(:, i, j, k), cs_iso2(i, j, k), this%cs
   end function neu_cs

   function get_tag() result(tag)
      use constants, only: idlen
      implicit none
      character(len=idlen) :: tag

      tag = "NEU"
   end function get_tag

!>
!! \brief Routine to set parameters from namelist FLUID_NEUTRAL
!!
!! \n \n
!! @b FLUID_NEUTRAL
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>gamma     </td><td>5./3.   </td><td>real value </td><td>\copydoc initneutral::gamma    </td></tr>
!! <tr><td>cs_iso    </td><td>1.0     </td><td>real value </td><td>\copydoc initneutral::cs_iso   </td></tr>
!! <tr><td>selfgrav  </td><td>.false. </td><td>logical    </td><td>\copydoc initneutral::selfgrav </td></tr>
!! </table>
!! The list is active while \b "NEUTRAL" is defined.
!! \n \n
!<
   subroutine init_neutral

      use dataio_pub, only: nh ! QA_WARN required for diff_nml
      use mpisetup,   only: rbuff, lbuff, master, slave, piernik_MPI_Bcast

      implicit none

      namelist /FLUID_NEUTRAL/ gamma, cs_iso, selfgrav

      gamma    = 5./3.
      cs_iso   = 1.0
      selfgrav = .false.

      if (master) then

         diff_nml(FLUID_NEUTRAL)

         lbuff(1)  = selfgrav

         rbuff(1)  = gamma
         rbuff(2)  = cs_iso

      endif

      call piernik_MPI_Bcast(rbuff)
      call piernik_MPI_Bcast(lbuff)

      if (slave) then

         selfgrav = lbuff(1)

         gamma    = rbuff(1)
         cs_iso   = rbuff(2)

      endif

   end subroutine init_neutral

   subroutine cleanup_neutral

      implicit none

   end subroutine cleanup_neutral

!==========================================================================================
!
! OPT: This routine may cost as much as 30% of rtvd. It seems that all the data fit well a 512kB L2 cache, but Ir:Dr:Dw is like 8:2:1
! OPT: \todo Try an explicit loop over RNG to check if we're better than the compiler
! OPT: similar treatment may be helpful for fluxionized.F90, fluxdust.F90 and fluxcosmicrays.F90
!
!/*
!>
!! \brief Computation of %fluxes for the neutral fluid
!!
!!The flux functions for neutral fluid are given by
!!
!!\f[
!!  \vec{F}{(\vec{u})} =
!!  \left(\begin{array}{c}
!!    \rho v_x \\
!!    \rho v_x^2 + p \\
!!    \rho v_x v_y\\
!!    \rho v_x v_z\\
!!    (e + p)v_x
!!  \end{array}\right),
!!  \qquad
!!  \vec{G}{(\vec{u})} =
!!  \left(\begin{array}{c}
!!    \rho v_y \\
!!    \rho v_y v_x\\
!!    \rho v_y^2 + p\\
!!    \rho v_y v_z\\
!!    (e + p)v_y
!!  \end{array}\right),
!!\qquad
!!  \vec{H}{(\vec{u})} =
!!  \left(\begin{array}{c}
!!    \rho v_z \\
!!    \rho v_z v_x\\
!!    \rho v_z v_y \\
!!    \rho v_z^2 + p \\
!!    (e + p)v_z
!!  \end{array}\right),
!!\f]
!<
!*/
#define RNG 2:nm
   subroutine flux_neu(this, flux, cfr, uu, n, vx, ps, bb, cs_iso2)

      use constants,    only: idn, imx, imy, imz
      use func,         only: ekin
#ifdef LOCAL_FR_SPEED
      use constants,    only: small, half
      use global,       only: cfr_smooth
#endif /* LOCAL_FR_SPEED */
#ifdef GLOBAL_FR_SPEED
      use timestep_pub, only: c_all
#endif /* GLOBAL_FR_SPEED */
#ifndef ISO
      use constants,    only: ien
      use global,       only: smallp
#endif /* !ISO */

      implicit none
      class(neutral_fluid), intent(in)             :: this
      integer(kind=4),      intent(in)             :: n         !< number of cells in the current sweep
      real, dimension(:,:), intent(inout), pointer :: flux      !< flux of neutral fluid
      real, dimension(:,:), intent(inout), pointer :: cfr       !< freezing speed for neutral fluid
      real, dimension(:,:), intent(in),    pointer :: uu        !< part of u for neutral fluid
      real, dimension(:),   intent(inout), pointer :: vx        !< velocity of neutral fluid for current sweep
      real, dimension(:),   intent(inout), pointer :: ps        !< pressure of neutral fluid for current sweep
      real, dimension(:,:), intent(in),    pointer :: bb        !< magnetic field x,y,z-components table
      real, dimension(:),   intent(in),    pointer :: cs_iso2   !< isothermal sound speed squared

      ! locals
      integer            :: nm
#ifdef LOCAL_FR_SPEED
      integer            :: i
      real               :: minvx     !<
      real               :: maxvx     !<
      real               :: amp       !<
#endif /* LOCAL_FR_SPEED */

      nm = n-1
      vx(RNG) = uu(RNG, imx)/uu(RNG, idn) ; vx(1) = vx(2); vx(n) = vx(nm)
#ifdef ISO
      ps(RNG) = cs_iso2(RNG) * uu(RNG, idn) ; ps(1) = ps(2); ps(n) = ps(nm)
#else /* !ISO */
      ps(RNG) = (uu(RNG, ien) - ekin(uu(RNG, imx),uu(RNG, imy),uu(RNG, imz),uu(RNG, idn)) )*(this%gam_1)
      ps(RNG) = max(ps(RNG), smallp)
#endif /* !ISO */

      flux(RNG, idn)=uu(RNG, imx)
      flux(RNG, imx)=uu(RNG, imx)*vx(RNG)+ps(RNG)
      flux(RNG, imy)=uu(RNG, imy)*vx(RNG)
      flux(RNG, imz)=uu(RNG, imz)*vx(RNG)
#ifndef ISO
      flux(RNG, ien)=(uu(RNG, ien)+ps(RNG))*vx(RNG)
#endif /* !ISO */
      flux(1, :) = flux(2, :) ; flux(n, :) = flux(nm, :)

#ifdef LOCAL_FR_SPEED

      ! The freezing speed is now computed locally (in each cell)
      !  as in Trac & Pen (2003). This ensures much sharper shocks,
      !  but sometimes may lead to numerical instabilities
      minvx = minval(vx(RNG))
      maxvx = maxval(vx(RNG))
      amp   = half*(maxvx-minvx)
      !    c_fr  = 0.0
#ifdef ISO
      cfr(RNG, 1) = sqrt(vx(RNG)**2+cfr_smooth*amp) + max(sqrt( abs(         ps(RNG))/uu(RNG, idn)),small)
#else /* !ISO */
      cfr(RNG, 1) = sqrt(vx(RNG)**2+cfr_smooth*amp) + max(sqrt( abs(this%gam*ps(RNG))/uu(RNG, idn)),small)
#endif /* !ISO */
      !> \deprecated BEWARE: that is the cause of fast decreasing of timestep in galactic disk problem
      !>
      !! \todo find why is it so
      !! if such a treatment is OK then should be applied also in both cases of neutral and ionized gas
      !!    do i = 2,nm
      !!       cfr(1,i) = maxval( [c_fr(i-1), c_fr(i), c_fr(i+1)] )
      !!    enddo
      !<

      cfr(1,1) = cfr(2,1);  cfr(n, 1) = cfr(nm, 1)
      do i = 2, this%all
         cfr(:, i) = cfr(:, 1)
      enddo
#endif /* LOCAL_FR_SPEED */

#ifdef GLOBAL_FR_SPEED
      ! The freezing speed is now computed globally
      !  (c=const for the whole domain) in subroutine 'timestep'

      !    cfr(:,:) = this%c   ! check which c_xxx is better
      cfr(:,:) = c_all
#endif /* GLOBAL_FR_SPEED */
      return
      if (.false.) write(0,*) bb, cs_iso2
#if defined(LOCAL_FR_SPEED) || defined(ISO)
      if (.false.) print *, this%all
#endif /* defined(LOCAL_FR_SPEED) || defined(ISO) */

   end subroutine flux_neu

end module initneutral
