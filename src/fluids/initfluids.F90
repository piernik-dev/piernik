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
#include "piernik.def"

!>
!! \brief (MH) Module to initialize all fluids, fluid components, tracers
!! and other dependent variables which are relevant for the current problem (doxy
!! comments ready).
!!
!! The Relaxing TVD scheme by Pen, Arras \& Wong (2003) was extended in PIERNIK
!! code to treat multiple fluids, by concatenation of the array \a u of
!! conservative variables of individual fluids
!!
!! \f[
!! \vec{u} =\big( \underbrace{\rho^{i}, m_x^{i}, m_y^{i}, m_z^{i}, e^{i}}_{\textrm{\scriptsize ionized gas}},
!! \underbrace{\rho^{n}, m_x^{n}, m_y^{n}, m_z^{n}, e^{n}}_{\textrm{\scriptsize neutral gas}},
!! \underbrace{\rho^{d}, m_x^{d}, m_y^{d}, m_z^{d}}_{\textrm{\scriptsize dust}}  \big)^T,
!! \f]
!! representing ionized gas, neutral gas, and dust treated as a pressureless fluid.
!! The Cosmic Ray (CR) gas  component, described by the diffusion--advection equation,
!! is incorporated in the same way.
!!
!! The general purpose of the multi-fluid framework is to simplify all code
!! modifications when a new fluid, component or variable is added.
!! We want to avoid, for example, modifications of boundary conditions routines,
!! when a new fluid is added. We formulate boundary condition routines in
!! a general manner and apply them to all variables of particular type (such as
!! gas density, x-momentum, ...) in one instant, through the array indexes:
!! iarr_all_dn, iarr_all_mx, etc ... (see comments to fluidindex module).
!!
!! \par DEFINITIONS
!!
!! \n\b Fluid \b component: an ingredient characterized by mass density, momenta,
!!    and optionally energy density.
!!    Examples: ionized fluid, neutral fluid, dust fluid.
!!    Variable \b "nfluid" (in fluidindex) counts fluids.
!!
!! \n\b Non-fluid \b component: Some constituents of the system are described by
!!      a system of equations, which does not involve momentum equation.
!!   These constituents are named non-fluid components
!!   Examples: CR energy density, or a set variables describing
!!   several CR energy bins, CR nuclear species or tracer or helicity.
!!   There is no specific counter for non-fluid components.
!!
!! \n\b Components: all fluid and non-fluid components. Variable \b "ncomponents"
!!      (in fluidindex) counts components.
!!
!! \n\b Non-isothermal \b fluid: an isothermal fluid does not need energy equation,
!!    thus \b non-isothermal fluids are distinguished as those, which do not
!!   engage the energy equation. Variable \b "nadiab" (in fluidindex) counts
!!      independent energy equations used for fluids.
!!
!! \n\b Fluid \b variable: Single quantity, such as gas density, x,y,z-momentum
!!   density, and energy density.
!!
!! \n\b Non-fluid \b variable: Single quantity, such as CR energy density,
!!      tracer, etc ... Variable \b "nvar" (in fluidindex) counts all fluid and
!!    non-fluid variables. There is no specific counter for non-fluid variables.
!!
!! \n All these constituents are organized in the module fluidindex into
!! the array of conservative variables \a u(ivar,:,:,:), where \a ivar is the
!! index of particular variable.
!!
!! \par The module initfluids is organized as follows:
!! \n (1)  All fluids defined in "piernik.def" are initialized subsequently.
!! \n (2)  The routine fluidindex is invoked to construct %arrays of indexes for
!!         each fluid, e.g.. ionized_index, neutral_index, etc.
!!         See fluidindex for more details.
!! \n (3)  Physical parameters common for all fluids are computed if necessary.
!!
!! \todo Change variable name "nadiab" to "nenerg". Reason: we are not limited
!!       to isothermal and  adiabatic equos. Energy equation will be used
!!       for non-adiabatic fluids in presence of cooling and heating.
!!
!! \todo Subdivide different fluids into species
!!
!! \warning check if cs_iso and cs_neu are correctly defined (end of init_fluids
!!
!!  subroutine) for your purposes (if used).
!<

module initfluids

  implicit none

  real, allocatable :: gamma(:)          !< array containing adiabatic indices of all fluids, indexed by   ifluid = i_ion, i_neu, etc.
  real              :: cs_iso            !< isothermal sound speed for a mixture of fluids
  real              :: cs_iso2           !< square of isothermal sound speed for a mixture of fluids

  contains

   subroutine init_fluids

      use fluidindex,      only: fluid_index, nvar
#ifdef VERBOSE
      use dataio_public,   only: printinfo
#endif /* VERBOSE */
#if defined NEUTRAL && defined IONIZED
      use dataio_public,   only: warn
#endif /* defined NEUTRAL && defined IONIZED  */
#ifdef IONIZED
      use initionized,     only: init_ionized, gamma_ion, cs_iso_ion, cs_iso_ion2
#endif /* IONIZED */
#ifdef NEUTRAL
      use initneutral,     only: init_neutral, gamma_neu, cs_iso_neu, cs_iso_neu2
#endif /* NEUTRAL */
#ifdef DUST
      use initdust,        only: init_dust
#endif /* DUST */
#ifdef COSM_RAYS
      use initcosmicrays,  only: init_cosmicrays
#endif /* COSM_RAYS */

      implicit none
#ifdef VERBOSE
      call printinfo("[initfluids:init_fluids]: commencing...")
#endif /* VERBOSE */

#ifdef IONIZED
      call init_ionized
#endif /* IONIZED */
#ifdef NEUTRAL
      call init_neutral
#endif /* NEUTRAL */
#ifdef DUST
      call init_dust
#endif /* DUST */
#ifdef COSM_RAYS
      call init_cosmicrays
#endif /* COSM_RAYS */

      call fluid_index

      allocate(gamma(nvar%fluids))

#if defined NEUTRAL && defined IONIZED
      if (cs_iso_neu /= cs_iso_ion) &
         call warn("[initfluids:init_fluids]: 'cs_iso_neu' and 'cs_iso_ion' should be equal")
#endif /* defined NEUTRAL && defined IONIZED  */

#ifdef IONIZED
      gamma(nvar%ion%pos) = gamma_ion
      cs_iso   = cs_iso_ion
      cs_iso2  = cs_iso_ion2
#endif /* IONIZED */

#ifdef NEUTRAL
      gamma(nvar%neu%pos) = gamma_neu
      cs_iso  = cs_iso_neu
      cs_iso2 = cs_iso_neu2
#endif /* NEUTRAL  */
#ifdef VERBOSE
      call printinfo("[initfluids:init_fluids]: finished. \o/")
#endif /* VERBOSE */

   end subroutine init_fluids

   subroutine cleanup_fluids

      use fluidindex, only: cleanup_fluid_index

      implicit none

      if (allocated(gamma)) deallocate(gamma)
      call cleanup_fluid_index

   end subroutine cleanup_fluids

end module initfluids
