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
!! \brief This module collect all sources to be added to an integration scheme
!<
module sources

   implicit none

   private
   public :: all_sources, care_for_positives, init_sources, prepare_sources, timestep_sources

contains

!/*
!>
!! \brief Subroutine computes any scheme sources (yet, now it is based on rtvd scheme)
!!
!! \todo Do not pass i1 and i2, pass optional pointer to gravacc instead
!<
!*/
   subroutine init_sources

      use interactions, only: init_interactions
#ifdef CORIOLIS
      use coriolis,     only: init_coriolis
#endif /* CORIOLIS */
#ifdef NON_INERTIAL
      use non_inertial, only: init_non_inertial
#endif /* NON_INERTIAL */
#ifdef SHEAR
      use shear,        only: init_shear
#endif /* SHEAR */
#ifdef SN_SRC
      use snsources,    only: init_snsources
#endif /* SN_SRC */
#ifdef THERM
      use thermal,      only: init_thermal
#endif /* THERM */

      implicit none

      call init_interactions                 ! requires flind and units

#ifdef CORIOLIS
      call init_coriolis                     ! depends on geometry
#endif /* CORIOLIS */

#ifdef NON_INERTIAL
      call init_non_inertial                 ! depends on geometry
#endif /* NON_INERTIAL */

#ifdef SHEAR
      call init_shear                        ! depends on fluids
#endif /* SHEAR */

#ifdef SN_SRC
      call init_snsources                    ! depends on grid and fluids/cosmicrays
#endif /* SN_SRC */

#ifdef THERM
      call init_thermal
#endif /* THERM */

   end subroutine init_sources

!/*
!>
!! \brief Subroutine computes any scheme sources (yet, now it is based on rtvd scheme)
!!
!! \todo Do not pass i1 and i2, pass optional pointer to gravacc instead
!<
!*/
   subroutine prepare_sources(cg)

      use grid_cont,  only: grid_container
#if defined(COSM_RAYS) && defined(IONIZED)
      use crhelpers,  only: div_v
      use fluidindex, only: flind
#endif /* COSM_RAYS && IONIZED */

      implicit none

      type(grid_container), pointer, intent(in) :: cg                 !< current grid piece

#if defined(COSM_RAYS) && defined(IONIZED)
      call div_v(flind%ion%pos, cg)
#endif /* COSM_RAYS && IONIZED */
      if (.false. .and. cg%is_old) return ! to supress compiler warnings

   end subroutine prepare_sources

!/*
!>
!! \brief Subroutine computes any scheme sources (yet, now it is based on rtvd scheme)
!!
!! \todo Do not pass i1 and i2, pass optional pointer to gravacc instead
!<
!*/
   subroutine all_sources(n, u, u1, bb, cg, istep, sweep, i1, i2, coeffdt, vel_sweep)

      use fluidindex,       only: flind, nmag
      use grid_cont,        only: grid_container
      use gridgeometry,     only: geometry_source_terms_exec
#ifdef BALSARA
      use interactions,     only: balsara_implicit_interactions
#else /* !BALSARA */
      use interactions,     only: fluid_interactions_exec
#endif /* !BALSARA */
#ifdef GRAV
      use gravity,          only: grav_src_exec
#endif /* GRAV */
#ifdef COSM_RAYS
      use sourcecosmicrays, only: src_gpcr_exec
#ifdef COSM_RAYS_SOURCES
      use sourcecosmicrays, only: src_crn_exec
#endif /* COSM_RAYS_SOURCES */
#endif /* COSM_RAYS */
#ifdef CORIOLIS
      use coriolis,         only: coriolis_force
#endif /* CORIOLIS */
#ifdef NON_INERTIAL
      use non_inertial,     only: non_inertial_force
#endif /* NON_INERTIAL */
#ifdef SHEAR
      use shear,            only: shear_acc
#endif /* SHEAR */
#ifdef THERM
      use thermal,          only: src_thermal_exec
#endif /* THERM */

      implicit none

      integer(kind=4),                          intent(in)    :: n                  !< array size
      real, dimension(n, flind%all),            intent(in)    :: u                  !< vector of conservative variables
      real, dimension(n, flind%all),            intent(inout) :: u1                 !< updated vector of conservative variables (after one timestep in second order scheme)
      real, dimension(n, nmag),                 intent(in)    :: bb                 !< local copy of magnetic field
      type(grid_container), pointer,            intent(in)    :: cg                 !< current grid piece
      integer,                                  intent(in)    :: istep              !< stage in the time integration scheme
      integer(kind=4),                          intent(in)    :: sweep              !< direction (x, y or z) we are doing calculations for
      integer,                                  intent(in)    :: i1                 !< coordinate of sweep in the 1st remaining direction
      integer,                                  intent(in)    :: i2                 !< coordinate of sweep in the 2nd remaining direction
      real,                                     intent(in)    :: coeffdt            !< time step times scheme coefficient
      real, dimension(n, flind%fluids), target, intent(in)    :: vel_sweep          !< velocity in the direction of current sweep

!locals
      real, dimension(n, flind%all)                           :: usrc, newsrc       !< u array update from sources
      real, dimension(:,:), pointer                           :: vx

      vx   => vel_sweep

      usrc = 0.0

      call geometry_source_terms_exec(u, bb, sweep, i1, i2, cg, newsrc)  ! n safe
      usrc(:,:) = usrc(:,:) + newsrc(:,:)

#ifndef BALSARA
      call get_updates_from_acc(n, u, usrc, fluid_interactions_exec(n, u, vx))  ! n safe
#else /* !BALSARA */
      call balsara_implicit_interactions(u1, vx, istep, sweep, i1, i2, cg) ! n safe
#endif /* !BALSARA */
#ifdef SHEAR
      call get_updates_from_acc(n, u, usrc, shear_acc(sweep,u)) ! n safe
#endif /* SHEAR */
#ifdef CORIOLIS
      call get_updates_from_acc(n, u, usrc, coriolis_force(sweep,u)) ! n safe
#endif /* CORIOLIS */
#ifdef NON_INERTIAL
      call get_updates_from_acc(n, u, usrc, non_inertial_force(sweep, u, cg))
#endif /* NON_INERTIAL */

#ifdef GRAV
      call grav_src_exec(n, u, cg, sweep, i1, i2, istep, newsrc)
      usrc(:,:) = usrc(:,:) + newsrc(:,:)
#endif /* !GRAV */

#if defined COSM_RAYS && defined IONIZED
      call src_gpcr_exec(u, n, newsrc, sweep, i1, i2, cg, vx)
      usrc(:,:) = usrc(:,:) + newsrc(:,:)
#ifdef COSM_RAYS_SOURCES
      call src_crn_exec(u, n, newsrc, coeffdt) ! n safe
      usrc(:,:) = usrc(:,:) + newsrc(:,:)
#endif /* COSM_RAYS_SOURCES */
#endif /* COSM_RAYS && IONIZED */
#ifdef THERM
      call src_thermal_exec(u, n, bb, newsrc)
      usrc(:,:) = usrc(:,:) + newsrc(:,:)
#endif /* THERM */

! --------------------------------------------------

      u1(:,:) = u1(:,:) + usrc(:,:) * coeffdt

      return
      if (.false.) write(0,*) bb, istep

   end subroutine all_sources

!/*
!>
!! \brief Subroutine computes any scheme sources (yet, now it is based on rtvd scheme)
!!
!! \todo Do not pass i1 and i2, pass optional pointer to gravacc instead
!<
!*/
   subroutine get_updates_from_acc(n, u, usrc, acc)

      use fluidindex, only: iarr_all_dn, iarr_all_mx, flind
#ifndef ISO
      use fluidindex, only: iarr_all_en
#endif /* !ISO */

      implicit none

      integer(kind=4),               intent(in)    :: n                  !< array size
      real, dimension(n, flind%all), intent(in)    :: u                  !< vector of conservative variables
      real, dimension(n, flind%all), intent(inout) :: usrc               !< u array update from sources
      real, dimension(n, flind%fluids), intent(in) :: acc                !< acceleration

      usrc(:, iarr_all_mx) = usrc(:, iarr_all_mx) + acc(:,:) * u(:, iarr_all_dn)
#ifndef ISO
      usrc(:, iarr_all_en) = usrc(:, iarr_all_en) + acc(:,:) * u(:, iarr_all_mx)
#endif /* !ISO */

   end subroutine get_updates_from_acc

!==========================================================================================

!/*
!>
!! \brief Subroutine collects dt limits estimations from sources
!<
!*/
   subroutine timestep_sources(dt)

#ifndef BALSARA
      use timestepinteractions, only: timestep_interactions
#endif /* !BALSARA */
#ifdef THERM
      use timestepthermal,      only: timestep_thermal
#endif /* THERM */

      implicit none

      real, intent(inout) :: dt

#ifndef BALSARA
         dt = min(dt, timestep_interactions())
#endif /* !BALSARA */
#ifdef THERM
         dt = min(dt, timestep_thermal())
#endif /* THERM */

      return
      if (.false. .and. dt < 0) return

   end subroutine timestep_sources

!==========================================================================================
   subroutine care_for_positives(n, u1, bb, cg, sweep, i1, i2)

      use fluidindex, only: flind, nmag
      use grid_cont,  only: grid_container

      implicit none

      integer(kind=4),               intent(in)    :: n                  !< array size
      real, dimension(n, flind%all), intent(inout) :: u1                 !< updated vector of conservative variables (after one timestep in second order scheme)
      real, dimension(n, nmag),      intent(in)    :: bb                 !< local copy of magnetic field
      type(grid_container), pointer, intent(in)    :: cg                 !< current grid piece
      integer(kind=4),               intent(in)    :: sweep              !< direction (x, y or z) we are doing calculations for
      integer,                       intent(in)    :: i1                 !< coordinate of sweep in the 1st remaining direction
      integer,                       intent(in)    :: i2                 !< coordinate of sweep in the 2nd remaining direction
!locals
      logical                                      :: full_dim

      full_dim = n > 1

      call limit_minimal_density(n, u1, cg, sweep, i1, i2)
      call limit_minimal_intener(n, bb, u1)
#ifdef COSM_RAYS
      if (full_dim) call limit_minimal_ecr(n, u1)
#endif /* COSM_RAYS */

   end subroutine care_for_positives

!==========================================================================================
   subroutine limit_minimal_density(n, u1, cg, sweep, i1, i2)

      use constants,   only: GEO_XYZ, GEO_RPZ, xdim, ydim, zdim, zero
      use dataio_pub,  only: msg, die, warn
      use domain,      only: dom
      use fluidindex,  only: flind, iarr_all_dn
      use global,      only: smalld, use_smalld, dn_negative, disallow_negatives
      use grid_cont,   only: grid_container
      use mass_defect, only: local_magic_mass

      implicit none

      integer(kind=4),               intent(in)    :: n                  !< array size
      real, dimension(n, flind%all), intent(inout) :: u1                 !< updated vector of conservative variables (after one timestep in second order scheme)
      type(grid_container), pointer, intent(in)    :: cg                 !< current grid piece
      integer(kind=4),               intent(in)    :: sweep              !< direction (x, y or z) we are doing calculations for
      integer,                       intent(in)    :: i1                 !< coordinate of sweep in the 1st remaining direction
      integer,                       intent(in)    :: i2                 !< coordinate of sweep in the 2nd remaining direction
      logical                                      :: dnneg

!locals

      integer :: ifl

      dnneg = any(u1(:, iarr_all_dn) < zero)
      dn_negative = dn_negative .or. dnneg
      if (use_smalld) then
         ! This is needed e.g. for outflow boundaries in presence of perp. gravity
         select case (dom%geometry_type)
            case (GEO_XYZ)
               local_magic_mass(:) = local_magic_mass(:) - sum(u1(dom%nb+1:n-dom%nb, iarr_all_dn), dim=1) * cg%dvol
               u1(:, iarr_all_dn) = max(u1(:, iarr_all_dn),smalld)
               local_magic_mass(:) = local_magic_mass(:) + sum(u1(dom%nb+1:n-dom%nb, iarr_all_dn), dim=1) * cg%dvol
            case (GEO_RPZ)
               select case (sweep)
                  case (xdim)
                     do ifl = lbound(iarr_all_dn, dim=1), ubound(iarr_all_dn, dim=1)
                        local_magic_mass(ifl) = local_magic_mass(ifl) - sum(u1(dom%nb+1:n-dom%nb, iarr_all_dn(ifl)) * cg%x(cg%is:cg%ie)) * cg%dvol
                     enddo
                     u1(:, iarr_all_dn) = max(u1(:, iarr_all_dn),smalld)
                     do ifl = lbound(iarr_all_dn, dim=1), ubound(iarr_all_dn, dim=1)
                        local_magic_mass(ifl) = local_magic_mass(ifl) + sum(u1(dom%nb+1:n-dom%nb, iarr_all_dn(ifl)) * cg%x(cg%is:cg%ie)) * cg%dvol
                     enddo
                  case (ydim)
                     local_magic_mass(:) = local_magic_mass(:) - sum(u1(dom%nb+1:n-dom%nb, iarr_all_dn), dim=1) * cg%dvol * cg%x(i2)
                     u1(:, iarr_all_dn) = max(u1(:, iarr_all_dn),smalld)
                     local_magic_mass(:) = local_magic_mass(:) + sum(u1(dom%nb+1:n-dom%nb, iarr_all_dn), dim=1) * cg%dvol * cg%x(i2)
                  case (zdim)
                     local_magic_mass(:) = local_magic_mass(:) - sum(u1(dom%nb+1:n-dom%nb, iarr_all_dn), dim=1) * cg%dvol * cg%x(i1)
                     u1(:, iarr_all_dn) = max(u1(:, iarr_all_dn),smalld)
                     local_magic_mass(:) = local_magic_mass(:) + sum(u1(dom%nb+1:n-dom%nb, iarr_all_dn), dim=1) * cg%dvol * cg%x(i1)
               end select
            case default
               call die("[sources:limit_minimal_density] Unsupported geometry")
         end select
      else
         if (dnneg) then
            write(msg,'(3A,I4,1X,I4,A)') "[sources:limit_minimal_density] negative density in sweep ",sweep,"( ", i1, i2, " )"
            if (disallow_negatives) then
               call warn(msg)
            else
               call die(msg)
            endif
         endif
      endif

   end subroutine limit_minimal_density

!==========================================================================================
   subroutine limit_minimal_intener(n, bb, u1)

      use constants,  only: xdim, ydim, zdim, zero
      use fluidindex, only: flind, nmag
      use fluidtypes, only: component_fluid
      use func,       only: emag, ekin
      use global,     only: smallei, use_smallei, ei_negative, disallow_negatives

      implicit none

      integer(kind=4),               intent(in)    :: n                  !< array size
      real, dimension(n, nmag),      intent(in)    :: bb                 !< local copy of magnetic field
      real, dimension(n, flind%all), intent(inout) :: u1                 !< updated vector of conservative variables (after one timestep in second order scheme)

!locals

      real, dimension(n)              :: kin_ener, int_ener, mag_ener
      class(component_fluid), pointer :: pfl
      integer                         :: ifl

      do ifl = 1, flind%fluids
         pfl => flind%all_fluids(ifl)%fl
         if (pfl%has_energy) then
            kin_ener = ekin(u1(:, pfl%imx), u1(:, pfl%imy), u1(:, pfl%imz), u1(:, pfl%idn))
            if (pfl%is_magnetized) then
               mag_ener = emag(bb(:, xdim), bb(:, ydim), bb(:, zdim))
               int_ener = u1(:, pfl%ien) - kin_ener - mag_ener
            else
               int_ener = u1(:, pfl%ien) - kin_ener
            endif

            if (disallow_negatives) ei_negative = ei_negative .or. (any(int_ener < zero))
            if (use_smallei) int_ener = max(int_ener, smallei)

            u1(:, pfl%ien) = int_ener + kin_ener
            if (pfl%is_magnetized) u1(:, pfl%ien) = u1(:, pfl%ien) + mag_ener
         endif
      enddo

   end subroutine limit_minimal_intener

#ifdef COSM_RAYS
   subroutine limit_minimal_ecr(n, u1)

      use constants,      only: zero
      use fluidindex,     only: flind
      use global,         only: cr_negative, disallow_CRnegatives
      use initcosmicrays, only: iarr_crs, smallecr, use_smallecr

      implicit none

      integer(kind=4),               intent(in)    :: n                  !< array size
      real, dimension(n, flind%all), intent(inout) :: u1                 !< updated vector of conservative variables (after one timestep in second order scheme)

      if (disallow_CRnegatives) cr_negative = cr_negative .or. (any(u1(:, iarr_crs(:)) < zero))
      if (use_smallecr) u1(:, iarr_crs(:)) = max(smallecr, u1(:, iarr_crs(:)))

   end subroutine limit_minimal_ecr
#endif /* COSM_RAYS */

!==========================================================================================
end module sources
