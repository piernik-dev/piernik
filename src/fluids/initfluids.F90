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

!>
!! \brief Module to initialize all fluids, fluid components, tracers
!! and other dependent variables which are relevant for the current problem.
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
!!    Variable \b "flind\%fluids" (in fluidindex) counts fluids.
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
!!   engage the energy equation. Variable \b "flind\%energ" (in fluidindex) counts
!!      independent energy equations used for fluids.
!!
!! \n\b Fluid \b variable: Single quantity, such as gas density, x,y,z-momentum
!!   density, and energy density.
!!
!! \n\b Non-fluid \b variable: Single quantity, such as CR energy density,
!!      tracer, etc ... Variable \b "flind" (in fluidindex) counts all fluid and
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
!! \todo Subdivide different fluids into species
!!
!<

module initfluids
! pulled by ANY
   implicit none
   private
   public :: init_fluids, cleanup_fluids, sanitize_smallx_checks

contains

   subroutine init_fluids

      use constants,      only: PIERNIK_INIT_GLOBAL
      use dataio_pub,     only: die, code_progress, warn
      use fluidindex,     only: fluid_index, flind
      use fluids_pub,     only: has_dst, has_ion, has_neu, cs2_max
      use fluxes,         only: set_limiter
      use global,         only: limiter
      use initdust,       only: init_dust
      use initionized,    only: init_ionized
      use initneutral,    only: init_neutral
      use mass_defect,    only: init_magic_mass
#ifdef COSM_RAYS
      use initcosmicrays, only: init_cosmicrays
#endif /* COSM_RAYS */
#ifdef TRACER
      use inittracer,     only: init_tracer
#endif /* TRACER */
#ifdef VERBOSE
      use dataio_pub,     only: printinfo
#endif /* VERBOSE */

      implicit none

      integer :: ifl

      if (code_progress < PIERNIK_INIT_GLOBAL) call die("[initfluids:init_fluids] MPI not initialized.") ! limiter, init_ionized, init_neutral, init_dust, init_cosmicrays

#ifdef VERBOSE
      call printinfo("[initfluids:init_fluids]: commencing...")
#endif /* VERBOSE */

      if (has_ion) call init_ionized
      if (has_neu) call init_neutral
      if (has_dst) call init_dust
#ifdef COSM_RAYS
      call init_cosmicrays
#endif /* COSM_RAYS */
#ifdef TRACER
      call init_tracer
#endif /* TRACER */

      call fluid_index    ! flind has valid values afterwards

      cs2_max = 0.0
      do ifl = lbound(flind%all_fluids, dim=1), ubound(flind%all_fluids, dim=1)
         cs2_max = max(cs2_max, flind%all_fluids(ifl)%fl%cs2)
      enddo
      ! All processes should have same fluids, so piernik_MPI_Allreduce shouldn't be required
      !call piernik_MPI_Allreduce(cs_max, pMAX)

      if (has_neu .and. has_ion) then
         if (flind%ion%cs2 /= flind%neu%cs2) &
            call warn("[initfluids:init_fluids]: flind%neu%cs2 and flind%ion%cs should be equal")
      endif

      !> \todo find a better place for the following (somewhere between calling fluid_index and reading restart)
      call init_magic_mass

      call set_limiter(limiter)
#ifdef VERBOSE
      call printinfo("[initfluids:init_fluids]: finished. \o/")
#endif /* VERBOSE */

   end subroutine init_fluids

   subroutine cleanup_fluids
      use fluids_pub,     only: has_ion
      use initionized,    only: cleanup_ionized
      use mass_defect,    only: cleanup_magic_mass
#ifdef COSM_RAYS
      use initcosmicrays, only: cleanup_cosmicrays
#endif /* COSM_RAYS */

      implicit none

      if (has_ion) call cleanup_ionized
#ifdef COSM_RAYS
      call cleanup_cosmicrays
#endif /* COSM_RAYS */

      call cleanup_magic_mass

   end subroutine cleanup_fluids

!>
!! \brief Find sane values for smalld and smallp
!! \warning Use span(cg%ijkse) or guarantee that all boundaries with corners have all guardcells with proper values
!<

   subroutine sanitize_smallx_checks

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: big_float, DST, xdim, ydim, zdim, cs_i2_n, pMAX, pMIN
      use dataio_pub,       only: warn, msg
      use fluidindex,       only: flind
      use fluidtypes,       only: component_fluid
      use func,             only: ekin, emag
      use global,           only: smalld, smallp
      use grid_cont,        only: grid_container
      use mpisetup,         only: master, piernik_MPI_Allreduce
      use named_array_list, only: qna, wna

      implicit none

      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg
      class(component_fluid), pointer :: fl
      integer                         :: i
      real, pointer, dimension(:,:,:) :: dn, mx, my, mz, en, bx, by, bz
      real, parameter                 :: safety_factor = 1.e-4
      real, parameter                 :: max_dens_span = 5.0
      real                            :: maxdens, span, mindens, minpres

      maxdens = 0.0
      mindens = smalld
      minpres = smallp
      ! collect the extrema
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         bx => cg%w(wna%bi)%span(xdim,cg%ijkse)
         by => cg%w(wna%bi)%span(ydim,cg%ijkse)
         bz => cg%w(wna%bi)%span(zdim,cg%ijkse)

         if (smalld >= big_float) then
            do i = lbound(flind%all_fluids,1), ubound(flind%all_fluids,1)
               dn => cg%w(wna%fi)%span(flind%all_fluids(i)%fl%idn,cg%ijkse)
               maxdens = max( maxval(dn), maxdens )
               mindens = min( minval(dn), mindens )
            enddo
         endif

         if (smallp >= big_float) then
            do i = lbound(flind%all_fluids,1), ubound(flind%all_fluids,1)
               fl => flind%all_fluids(i)%fl
               if (fl%tag == DST) cycle
               dn => cg%w(wna%fi)%span(fl%idn,cg%ijkse)
               mx => cg%w(wna%fi)%span(fl%imx,cg%ijkse)
               my => cg%w(wna%fi)%span(fl%imy,cg%ijkse)
               mz => cg%w(wna%fi)%span(fl%imz,cg%ijkse)
               if (fl%has_energy) then
                  en => cg%w(wna%fi)%span(fl%ien,cg%ijkse)
                  if (fl%is_magnetized) then
                     minpres = min( minval( en - ekin(mx,my,mz,dn) - emag(bx,by,bz))/fl%gam_1, minpres )
                  else
                     minpres = min( minval( en - ekin(mx,my,mz,dn))/fl%gam_1, minpres )
                  endif
               else
                  minpres = min( minval( cg%q(qna%ind(cs_i2_n))%span(cg%ijkse)*dn ), minpres )
               endif
            enddo
         endif
         cgl => cgl%nxt
      enddo

      ! reduce across processes
      if (smalld >= big_float) then
         smalld = mindens
         call piernik_MPI_Allreduce(smalld,  pMIN)
         call piernik_MPI_Allreduce(maxdens, pMAX)
         span = 0
         if (maxdens > mindens) span = log10(maxdens/mindens)
         mindens = mindens * safety_factor
         if (master) then
            write(msg,'(A,ES11.4)') "[initfluids:sanitize_smallx_checks] adjusted smalld to ", smalld
            call warn(msg)
            if (span > max_dens_span) then
               write(msg,'(A,I3,A)') "[initfluids:sanitize_smallx_checks] density spans over ", int(span), " orders of magnitude!"
               call warn(msg)
            endif
         endif
      endif

      if (smallp >= big_float) then
         smallp = minpres * safety_factor
         call piernik_MPI_Allreduce(smallp, pMIN)
         if (smallp < 0.) then
            write(msg,'(A,ES11.4,A)') "[initfluids:sanitize_smallx_checks] Negative smallp detected! smallp=",smallp," may indicate nonphysical initial conditions."
            if (master) call warn(msg)
            smallp = tiny(1.)
         endif
         if (master) then
            write(msg,'(A,ES11.4)') "[initfluids:sanitize_smallx_checks] adjusted smallp to ", smallp
            call warn(msg)
         endif
      endif

      if (associated(dn)) nullify(dn)
      if (associated(mx)) nullify(mx)
      if (associated(my)) nullify(my)
      if (associated(mz)) nullify(mz)
      if (associated(en)) nullify(en)
      if (associated(bx)) nullify(bx)
      if (associated(by)) nullify(by)
      if (associated(bz)) nullify(bz)
      if (associated(fl)) nullify(fl)

   end subroutine sanitize_smallx_checks

end module initfluids
