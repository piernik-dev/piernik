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
module rtvd_sources

! pulled by RTVD

   implicit none

   private
   public  :: rtvd_sources_proc

contains
!/*
!>
!! \brief Subroutine computes any scheme sources (yet, now it is based on rtvd scheme)
!!
!! \todo Do not pass i1 and i2, pass optional pointer to gravacc instead
!<
!*/
   subroutine rtvd_sources_proc(n, u, u0, cs_iso2, istep, sweep, i1, i2, dx, dt, cg, u1, full_dim, pressure, vel_sweep)

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
      real, dimension(:), pointer,   intent(in)    :: cs_iso2            !< square of local isothermal sound speed
      integer,                       intent(in)    :: istep              !< step number in the time integration scheme
      integer(kind=4),               intent(in)    :: sweep              !< direction (x, y or z) we are doing calculations for
      integer,                       intent(in)    :: i1                 !< coordinate of sweep in the 1st remaining direction
      integer,                       intent(in)    :: i2                 !< coordinate of sweep in the 2nd remaining direction
      real,                          intent(in)    :: dx                 !< cell length
      real,                          intent(in)    :: dt                 !< time step
      type(grid_container), pointer, intent(in)    :: cg                 !< current grid piece
      real, dimension(n, flind%all), intent(inout) :: u1                 !< updated vector of conservative variables (after one timestep in second order scheme)
      logical,                       intent(in)    :: full_dim
      real, dimension(n, flind%fluids), intent(in) :: pressure           !< gas pressure
      real, dimension(n, flind%fluids), target, intent(in) :: vel_sweep          !< velocity in the direction of current sweep

#ifdef GRAV
      integer                                      :: ind                !< fluid index
      real, dimension(n)                           :: gravacc            !< acceleration caused by gravitation
#endif /* GRAV */

!locals
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

      real, dimension(2,2), parameter              :: rk2coef = reshape( [ one, half, zero, one ], [ 2, 2 ] )

      vx   => vel_sweep
      dens => density

      density(:,:) = u(:, iarr_all_dn)

         geosrc = geometry_source_terms(u, pressure, sweep, cg)  ! n safe

         u1(:, iarr_all_mx) = u1(:, iarr_all_mx) + rk2coef(integration_order,istep)*geosrc(:,:)*dt ! n safe

         acc = 0.0
#ifndef BALSARA
         acc = acc + fluid_interactions(dens, vx)  ! n safe
#else /* !BALSARA */
         call balsara_implicit_interactions(u1, u0, vx, cs_iso2, dt, istep) ! n safe
#endif /* !BALSARA */
#ifdef SHEAR
         acc = acc + shear_acc(sweep,u) ! n safe
#endif /* SHEAR */
#ifdef CORIOLIS
         acc = acc + coriolis_force(sweep,u) ! n safe
#endif /* CORIOLIS */
#ifdef NON_INERTIAL
         acc = acc + non_inertial_force(sweep, u, cg)
#endif /* NON_INERTIAL */

         if (full_dim) then
#ifdef GRAV
            call grav_pot2accel(sweep, i1, i2, n, gravacc, istep, cg)

            do ind = 1, flind%fluids
               acc(:, ind) =  acc(:, ind) + gravacc(:)
            enddo
#endif /* !GRAV */

            acc(n, :) = acc(n-1, :)
            acc(1, :) = acc(2, :)
         endif

         u1(:, iarr_all_mx) = u1(:, iarr_all_mx) + rk2coef(integration_order,istep)*acc(:,:)*u(:, iarr_all_dn)*dt
#ifndef ISO
         u1(:, iarr_all_en) = u1(:, iarr_all_en) + rk2coef(integration_order,istep)*acc(:,:)*u(:, iarr_all_mx)*dt
#endif /* !ISO */

! --------------------------------------------------

#if defined COSM_RAYS && defined IONIZED
         if (full_dim) then
            call src_gpcr(u, n, dx, decr, grad_pcr, sweep, i1, i2, cg)
            u1(:,                iarr_crs(:)) = u1(:,               iarr_crs(:)) + rk2coef(integration_order,istep) * decr(:,:) * dt
            u1(:,                iarr_crs(:)) = max(smallecr, u1(:, iarr_crs(:)))
            u1(:, iarr_all_mx(flind%ion%pos)) = u1(:, iarr_all_mx(flind%ion%pos)) + rk2coef(integration_order,istep) * grad_pcr * dt
#ifndef ISO
            u1(:, iarr_all_en(flind%ion%pos)) = u1(:, iarr_all_en(flind%ion%pos)) + rk2coef(integration_order,istep) * vx(:, flind%ion%pos) * grad_pcr * dt
#endif /* !ISO */
         endif
#ifdef COSM_RAYS_SOURCES
         call src_crn(u, n, srccrn, rk2coef(integration_order, istep) * dt) ! n safe
         u1(:, iarr_crn) = u1(:, iarr_crn) +  rk2coef(integration_order, istep)*srccrn(:,:)*dt
#endif /* COSM_RAYS_SOURCES */
#endif /* COSM_RAYS && IONIZED */

   end subroutine rtvd_sources_proc

!==========================================================================================
end module rtvd_sources
