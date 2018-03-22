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
!! \brief This module implements relaxing TVD scheme
!!
!! The implementation was based on TVD split MHD code by Pen et al. (2003).
!<
module sources

! pulled by RTVD

   implicit none

   private
   public  :: all_sources, prepare_sources

contains

!/*
!>
!! \brief Subroutine computes any scheme sources (yet, now it is based on rtvd scheme)
!!
!! \todo Do not pass i1 and i2, pass optional pointer to gravacc instead
!<
!*/
   subroutine prepare_sources(cg)

      use grid_cont,  only: grid_container
#ifdef COSM_RAYS
      use crhelpers,  only: div_v
      use fluidindex, only: flind
#endif /* COSM_RAYS */

      implicit none

      type(grid_container), pointer, intent(in)    :: cg                 !< current grid piece

#ifdef COSM_RAYS
                        call div_v(flind%ion%pos, cg)
#endif /* COSM_RAYS */
      if (.false. .and. cg%is_old) return ! to supress compiler warnings

   end subroutine prepare_sources

!/*
!>
!! \brief Subroutine computes any scheme sources (yet, now it is based on rtvd scheme)
!!
!! \todo Do not pass i1 and i2, pass optional pointer to gravacc instead
!<
!*/
   subroutine all_sources(n, u, u0, istep, sweep, i1, i2, dx, dt, cg, u1, pressure, vel_sweep)

      use constants,        only: one, zero, half
      use fluidindex,       only: iarr_all_dn, iarr_all_mx, flind
#ifndef ISO
      use fluidindex,       only: iarr_all_en
#endif /* !ISO */
      use global,           only: integration_order
      use grid_cont,        only: grid_container
      use gridgeometry,     only: geometry_source_terms
#ifdef BALSARA
      use interactions,     only: balsara_implicit_interactions
#else /* !BALSARA */
      use interactions,     only: fluid_interactions
#endif /* !BALSARA */
#ifdef GRAV
      use gravity,          only: grav_pot2accel
#endif /* GRAV */
#ifdef COSM_RAYS
      use initcosmicrays,   only: iarr_crs, smallecr
      use sourcecosmicrays, only: src_gpcr
#ifdef COSM_RAYS_SOURCES
      use initcosmicrays,   only: iarr_crn
      use sourcecosmicrays, only: src_crn
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

      implicit none

      integer(kind=4),               intent(in)    :: n                  !< array size
      real, dimension(n, flind%all), intent(inout) :: u                  !< vector of conservative variables
      real, dimension(n, flind%all), intent(in)    :: u0                 !< vector of conservative variables
      integer,                       intent(in)    :: istep              !< step number in the time integration scheme
      integer(kind=4),               intent(in)    :: sweep              !< direction (x, y or z) we are doing calculations for
      integer,                       intent(in)    :: i1                 !< coordinate of sweep in the 1st remaining direction
      integer,                       intent(in)    :: i2                 !< coordinate of sweep in the 2nd remaining direction
      real,                          intent(in)    :: dx                 !< cell length
      real,                          intent(in)    :: dt                 !< time step
      type(grid_container), pointer, intent(in)    :: cg                 !< current grid piece
      real, dimension(n, flind%all), intent(inout) :: u1                 !< updated vector of conservative variables (after one timestep in second order scheme)
      real, dimension(n, flind%fluids), intent(in) :: pressure           !< gas pressure
      real, dimension(n, flind%fluids), target, intent(in) :: vel_sweep          !< velocity in the direction of current sweep

#ifdef GRAV
      integer                                      :: ind                !< fluid index
      real, dimension(n)                           :: gravacc            !< acceleration caused by gravitation
#endif /* GRAV */

!locals
      real, dimension(n, flind%all)                 :: usrc               !< u array update from sources
      real, dimension(n, flind%fluids)              :: acc                !< acceleration
      real, dimension(n, flind%fluids)              :: geosrc             !< source terms caused by geometry of coordinate system
      real, dimension(n, flind%fluids), target      :: density            !< gas density
      real, dimension(:,:),            pointer      :: dens, vx

#ifdef COSM_RAYS
      real, dimension(n)                            :: grad_pcr
      real, dimension(n, flind%crs%all)             :: decr
#ifdef COSM_RAYS_SOURCES
      real, dimension(n, flind%crn%all)             :: srccrn
#endif /* COSM_RAYS_SOURCES */
#endif /* COSM_RAYS */
      logical                                       :: full_dim

      real, dimension(2,2), parameter               :: rk2coef = reshape( [ one, half, zero, one ], [ 2, 2 ] )

      full_dim = n > 1

      vx   => vel_sweep
      dens => density

      density(:,:) = u(:, iarr_all_dn)

      geosrc = geometry_source_terms(u, pressure, sweep, cg)  ! n safe

      u1(:, iarr_all_mx) = u1(:, iarr_all_mx) + rk2coef(integration_order,istep)*geosrc(:,:)*dt ! n safe

      usrc = 0.0
#ifndef BALSARA
      acc = fluid_interactions(dens, vx)  ! n safe
      usrc(:, iarr_all_mx) = usrc(:, iarr_all_mx) + acc(:,:) * u(:, iarr_all_dn)
#ifndef ISO
      usrc(:, iarr_all_en) = usrc(:, iarr_all_en) + acc(:,:) * u(:, iarr_all_mx)
#endif /* !ISO */
#else /* !BALSARA */
      call balsara_implicit_interactions(u1, u0, vx, istep, sweep, i1, i2, cg) ! n safe
#endif /* !BALSARA */
#ifdef SHEAR
      acc = shear_acc(sweep,u) ! n safe
      usrc(:, iarr_all_mx) = usrc(:, iarr_all_mx) + acc(:,:) * u(:, iarr_all_dn)
#ifndef ISO
      usrc(:, iarr_all_en) = usrc(:, iarr_all_en) + acc(:,:) * u(:, iarr_all_mx)
#endif /* !ISO */
#endif /* SHEAR */
#ifdef CORIOLIS
      acc = coriolis_force(sweep,u) ! n safe
      usrc(:, iarr_all_mx) = usrc(:, iarr_all_mx) + acc(:,:) * u(:, iarr_all_dn)
#ifndef ISO
      usrc(:, iarr_all_en) = usrc(:, iarr_all_en) + acc(:,:) * u(:, iarr_all_mx)
#endif /* !ISO */
#endif /* CORIOLIS */
#ifdef NON_INERTIAL
      acc = non_inertial_force(sweep, u, cg)
      usrc(:, iarr_all_mx) = usrc(:, iarr_all_mx) + acc(:,:) * u(:, iarr_all_dn)
#ifndef ISO
      usrc(:, iarr_all_en) = usrc(:, iarr_all_en) + acc(:,:) * u(:, iarr_all_mx)
#endif /* !ISO */
#endif /* NON_INERTIAL */

      if (full_dim) then
#ifdef GRAV
         call grav_pot2accel(sweep, i1, i2, n, gravacc, istep, cg)

         do ind = 1, flind%fluids
            usrc(:, iarr_all_mx(ind)) = usrc(:, iarr_all_mx(ind)) + gravacc(:) * u(:, iarr_all_dn(ind))
#ifndef ISO
            usrc(:, iarr_all_en(ind)) = usrc(:, iarr_all_en(ind)) + gravacc(:) * u(:, iarr_all_mx(ind))
#endif /* !ISO */
         enddo
#endif /* !GRAV */
      endif

! --------------------------------------------------

#if defined COSM_RAYS && defined IONIZED
      if (full_dim) then
         call src_gpcr(u, n, dx, decr, grad_pcr, sweep, i1, i2, cg)
         usrc(:,                iarr_crs(:)) = usrc(:,               iarr_crs(:)) + decr(:,:)
         usrc(:, iarr_all_mx(flind%ion%pos)) = usrc(:, iarr_all_mx(flind%ion%pos)) + grad_pcr
#ifndef ISO
         usrc(:, iarr_all_en(flind%ion%pos)) = usrc1(:, iarr_all_en(flind%ion%pos)) + vx(:, flind%ion%pos) * grad_pcr
#endif /* !ISO */
      endif
#ifdef COSM_RAYS_SOURCES
      call src_crn(u, n, srccrn, rk2coef(integration_order, istep) * dt) ! n safe
      usrc(:, iarr_crn) = usrc(:, iarr_crn) + srccrn(:,:)
#endif /* COSM_RAYS_SOURCES */
#endif /* COSM_RAYS && IONIZED */
      u1(:,:) = u1(:,:) + rk2coef(integration_order, istep)*usrc(:,:)*dt
#if defined COSM_RAYS && defined IONIZED
      if (full_dim) u1(:, iarr_crs(:)) = max(smallecr, u1(:, iarr_crs(:)))
#endif /* COSM_RAYS && IONIZED */

   end subroutine all_sources

!==========================================================================================
end module sources
