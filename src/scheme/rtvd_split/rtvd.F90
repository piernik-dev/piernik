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
!! \brief (JD) (doxy comments ready) This module implements relaxing TVD scheme
!!
!! The implementation was based on TVD split MHD code by Pen et al. (2003).
!<
module rtvd ! split orig
! pulled by ANY
   implicit none
   private
   public  :: tvdb, relaxing_tvd
   contains
!/*
!>
!! \brief Subroutine computes magnetic field evolution
!!
!! The original RTVD MHD scheme incorporates magnetic field evolution via the Constrained transport (CT)
!! algorithm by Evans & Hawley (1988). The idea behind the CT scheme is to integrate numerically the induction equation
!! \f{equation}
!! \frac{\partial\vec{B}}{\partial t} = \nabla\times (\vec{v}  \times \vec{B}),
!! \f}
!! in a manner ensuring that the condition
!! \f{equation}
!! \nabla \cdot \vec{B} = 0,
!! \f}
!! is fulfilled to the machine accuracy.
!!
!! The divergence-free evolution of magnetic-field on a discrete computational
!! grid can be realized if the last equation holds for the discrete representation of the initial
!! condition, and that subsequent updates of \f$B\f$ do not change the total magnetic
!! flux threading cell faces.
!!
!! Time variations of magnetic flux threading surface \f$S\f$ bounded by contour \f$C\f$, due to Stokes theorem, can be written as
!! \f{equation}
!! \frac{\partial\Phi_S}{\partial t} = \frac{\partial}{\partial t} \int_S \vec{B} \cdot \vec{d\sigma} = \int_S \nabla \times ( \vec{v} \times \vec{B})\cdot \vec{d \sigma} =  \oint_{C} (\vec{v} \times \vec{B}) \cdot \vec{dl},
!! \f}
!! where \f$\vec{E}= \vec{v} \times \vec{B}\f$ is electric field, named also electromotive force (EMF).
!!
!! In a discrete representation variations of magnetic %fluxes, in %timestep \f$\Delta t\f$, threading faces of cell \f$(i,j,k)\f$ are
!! given by
!! \f{eqnarray}
!! \frac{\Phi^{x,n+1}_{i+1/2,j,k} - \Phi^{x,n}_{i+1/2,j,k}}{\Delta t} &=&
!! {E}^y_{i+1/2,j,k-1/2} \Delta y + {E}^z_{i+1/2,j+1/2,k} \Delta z     \nonumber\\
!! &-&{E}^y_{i+1/2,j,k+1/2} \Delta y - {E}^z_{i+1/2,j-1/2,k} \Delta z,  \nonumber
!! \f}
!! \f{eqnarray}
!! \frac{\Phi^{y,n+1}_{i,j+1/2,k} - \Phi^{y,n}_{i,j+1/2,k}}{\Delta t} &=&
!! {E}^x_{i,j+1/2,k+1/2} \Delta x + {E}^z_{i-1/2,j+1/2,k} \Delta z      \nonumber\\
!! &-&{E}^x_{i,j+1/2,k-1/2} \Delta x - {E}^z_{i+1/2,j+1/2,k} \Delta z,  \nonumber
!! \f}
!! \f{eqnarray}
!! \frac{\Phi^{z,n+1}_{i,j,k+1/2} - \Phi^{z,n}_{i,j,k+1/2}}{\Delta t} &=&
!! {E}^x_{i,j-1/2,k+1/2} \Delta x + {E}^y_{i+1/2,j,k+1/2} \Delta y      \nonumber\\
!! &-&{E}^x_{i,j+1/2,k+1/2} \Delta x - {E}^y_{i-1/2,j,k+1/2} \Delta y,  \nonumber
!! \f}
!! Variations of magnetic flux threading remaining three faces of cell \f$(i,j,k)\f$ can be written in a similar manner. We note that each EMF
!! contribution appears twice with opposite sign. Thus, the total change of magnetic flux, piercing all cell-faces, vanishes to machine accuracy.
!!
!! To present the idea in a slightly different way, following Pen et al. (2003),  we consider the three components of the induction
!! equation written in the explicit form:
!! \f{eqnarray}
!! \partial_t B_x &=& \partial_y (\mathrm{v_x B_y}) + \partial_z (v_x B_z) - \partial_y (v_y B_x) - \partial_z (v_z B_x),\\
!! \partial_t B_y &=& \partial_x (v_y B_x) + \partial_z (v_y B_z) - \partial_x (\mathrm{v_x B_y}) - \partial_z (v_z B_y),\\
!! \partial_t B_z &=& \partial_x (v_z B_x) + \partial_y (v_z B_y) - \partial_x (v_x B_z) - \partial_y (v_y B_z).
!! \f}
!!
!!We note that each combination of \f$v_a B_b\f$ appears twice in these equations. Let us consider \f$v_x B_y\f$.
!!Once  \f$v_x B_y\f$ is computed for numerical integration of the second equation, it should be also used for
!!integration of the first equation, to ensure cancellation of electromotive forces contributing to the total change of magnetic flux threading
!!cell faces.
!!
!!The scheme proposed by Pen et al. (2003), consists of the following steps:
!!\n (1) Computation of the edge--centered EMF component \f$ v_x B_y\f$.
!!\n (2) Update of $B_y$, according to the equation \f$\partial_t B_y = \partial_x (v_x B_y)\f$.
!!\n (3) Update of $B_x$, according to the equation \f$\partial_t B_x = \partial_y (v_x B_y)\f$.
!!
!!Analogous procedure applies to remaining EMF components.
!<
!*/
   subroutine tvdb(vibj,b,vg,n,dt,idi)
      use constants, only: big
      implicit none
      integer, intent(in) :: n       !< array size
      real, intent(in)    :: dt      !< time step
      real, intent(in)    :: idi     !< cell length, depends on direction x, y or z
      real, dimension(n)  :: vibj    !< face-centered electromotive force components (b*vg)
      real, dimension(n)  :: b       !< magnetic field
      real, dimension(n)  :: vg      !< velocity in the center of cell boundary
! locals
      real, dimension(n)  :: b1      !< magnetic field
      real, dimension(n)  :: vibj1   !< face-centered electromotive force (EMF) components (b*vg)
      real, dimension(n)  :: vh      !< velocity interpolated to the cell edges
      real :: dti                    !< dt/di
      real :: v                      !< auxiliary variable to compute EMF
      real :: w                      !< EMF component
      real :: dw                     !< The second-order correction to EMF component
      real :: dwm                    !< face centered EMF interpolated to left cell-edge
      real :: dwp                    !< face centered EMF interpolated to right cell-edge
      integer :: i                   !< auxiliary array indicator
      integer :: ip                  !< i+1
      integer :: ipp                 !< i+2
      integer :: im                  !< i-1

  ! unlike the B field, the vibj lives on the right cell boundary
      vh = 0.0

! velocity interpolation to the cell boundaries

      vh(1:n-1) =(vg(1:n-1)+ vg(2:n))*0.5;     vh(n) = vh(n-1)

      dti = dt*idi

! face-centered EMF components computation, depending on the sign of vh, the components are upwinded  to cell edges, leading to 1st order EMF

      where (vh > 0.)
         vibj1=b*vg
      elsewhere
         vibj1=eoshift(b*vg,1,boundary=big)
      endwhere

! values of magnetic field computation in Runge-Kutta half step

      b1(2:n) = b(2:n) -(vibj1(2:n)-vibj1(1:n-1))*dti*0.5;    b1(1) = b(2)

      do i = 3, n-3
         ip  = i  + 1
         ipp = ip + 1
         im  = i  - 1
         v   = vh(i)

! recomputation of EMF components (w) with b1 and face centered EMF interpolation to cell-edges (dwp, dwm), depending on the sign of v.

         if (v > 0.0) then
            w=vg(i)*b1(i)
            dwp=(vg(ip)*b1(ip)-w)*0.5
            dwm=(w-vg(im)*b1(im))*0.5
         else
            w=vg(ip)*b1(ip)
            dwp=(w-vg(ipp)*b1(ipp))*0.5
            dwm=(vg(i)*b1(i)-w)*0.5
         endif

! the second-order corrections to the EMF components computation with the aid of the van Leer monotonic interpolation and 2nd order EMF computation

         dw=0.0
         if (dwm*dwp > 0.0) dw=2.0*dwm*dwp/(dwm+dwp)
         vibj(i)=(w+dw)*dt
      enddo

   end subroutine tvdb
!/*
!>
!! \brief Subroutine implements the Relaxing TVD scheme for conserved physical quantities
!!
!! The main point of the Relaxing TVD scheme is to decompose vectors of conservative variables \f$u\f$ and %fluxes \f$F\f$
!! into left-moving and right-moving waves:
!! \f{equation}
!! u=u^L+u^P,
!! \f}
!! where
!! \f{equation}
!! u^L=\frac{1}{2}\left(U-\frac{F}{c}\right), u^P=\frac{1}{2}\left(U+\frac{F}{c}\right).
!! \f}
!! The symbol \f$c\f$ stands for freezing speed - a function that satisfies \f$c\ge \max\left(|v\pm c_f|\right)\f$,
!! \f$v\f$ is the fluid velocity and \f$c_f\f$ is the fast magnetosonic speed.
!! The %fluxes then can be written as
!! \f{equation}
!! F^L=-cu^L, F^P=cu^P,
!! \f}
!! and their sum is the flux of \f$u\f$
!! \f{equation}
!! F=F^L+F^P.
!! \f}
!!
!! To make time integration PIERNIK uses Runge-Kutta scheme. We can choose between the first and the second order accuracy.
!! The second order scheme consists of the following steps:
!! \n (1) First order %fluxes \f$\vec{F}_{i\pm 1/2}^{(1)L,R}\f$   are
!!    calculated together with source terms \f$\vec{S}_i(\vec{u})\f$ at \f$t^n\f$.
!! \n (2) Fluxes and source terms derived in the first step are used to calculate
!!    \f$\vec{u}^{n+1/2}\f$ at \f$t^{n+1/2}\f$.
!! \n (3) Second order %fluxes \f$\vec{F}_{i+1/2}^{(2)}\f$, obtained via monotonic, upwind
!!    interpolation to cell boundaries,  and source terms \f$\vec{S}\f$ are evaluated for
!!    \f$\vec{u}^{n+1/2}\f$ at \f$t^{n+1/2}\f$.
!! \n (4) Update of \f$\Delta\vec{u}\f$, corresponding to the full
!!    %timestep \f$\Delta t\f$ is done, using %fluxes and source terms calculated
!!    at the intermediate time step \f$t^{n+1/2}\f$.
!!
!! To achieve second order spatial accuracy, a monotone upwind interpolation of %fluxes onto cell boundaries is made,
!! with the aid of a flux limiter, to obtain the Monotone Upwind Scheme for Conservation Laws (MUSCL) (step (3) in Runge-Kutta scheme).
!! The monotone interpolation is used to avoid spurious oscillations in discrete solutions of higher order.
!! The Total Variation Diminishing property of the numerical scheme is related to the measure of the overall amount
!! of oscillations, called Total Variation.
!<
!*/
   subroutine relaxing_tvd(n, u, bb, sweep, i1, i2, dx, dt)

      use dataio_pub,       only: msg, die
      use fluidindex,       only: iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz, ibx, iby, ibz, nvar, nmag
      use fluxes,           only: flimiter, all_fluxes
      use mpisetup,         only: smalld, integration_order, use_smalld, local_magic_mass
      use gridgeometry,     only: gc, geometry_source_terms
#ifndef ISO
      use fluidindex,       only: iarr_all_en
      use mpisetup,         only: smallei
#endif /* !ISO */
#ifdef GRAV
      use gravity,          only: grav_pot2accel
#endif /* GRAV */
#ifdef SHEAR
      use grid,             only: x
      use shear,            only: qshear, omega, global_gradP
#endif /* SHEAR */
#ifdef COSM_RAYS
      use arrays,           only: divvel
      use initcosmicrays,   only: iarr_crs, iarr_crn, gamma_crs, cr_active, smallecr
#ifdef COSM_RAYS_SOURCES
      use sourcecosmicrays, only: src_crn
#endif /* COSM_RAYS_SOURCES */
#endif /* COSM_RAYS */
#ifdef FLUID_INTERACTIONS
      use initdust,         only: dragc_gas_dust
      use interactions,     only: fluid_interactions
#endif /* FLUID_INTERACTIONS */
#ifdef ISO_LOCAL
      use arrays,           only: cs_iso2_arr
#endif /* ISO_LOCAL */

      implicit none

      integer,                     intent(in)  :: n                  !< array size
      real, dimension(nvar%all,n), intent(out) :: u                  !< vector of conservative variables
      real, dimension(nmag,n),     intent(in)  :: bb                 !< local copy of magnetic field
      character(len=*),            intent(in)  :: sweep              !< direction (x, y or z) we are doing calculations for
      integer,                     intent(in)  :: i1                 !< coordinate of sweep in the 1st remaining direction
      integer,                     intent(in)  :: i2                 !< coordinate of sweep in the 2nd remaining direction
      real,                        intent(in)  :: dx                 !< cell length
      real,                        intent(in)  :: dt                 !< time step

      integer                        :: istep              !< step number in the time integration scheme
#if defined GRAV || defined SHEAR || defined FLUID_INTERACTIONS
      integer                        :: ind                !< fluid index
#endif /* defined GRAV || defined SHEAR || defined FLUID_INTERACTIONS */

      real                           :: dtx                !< dt/dx
      real, dimension(nvar%all,n)    :: cfr                !< freezing speed
#if defined GRAV || defined SHEAR || defined FLUID_INTERACTIONS
      real, dimension(nvar%fluids,n) :: acc                !< acceleration
#endif /* defined GRAV || defined SHEAR || defined FLUID_INTERACTIONS */
!locals
      real, dimension(nvar%all,n)    :: u0                 !< initial state of conservative variables
      real, dimension(nvar%all,n)    :: w                  !< auxiliary vector to calculate fluxes
      real, dimension(nvar%all,n)    :: fr                 !< flux of the right-moving waves
      real, dimension(nvar%all,n)    :: fl                 !< flux of the left-moving waves
      real, dimension(nvar%all,n)    :: fu                 !< sum of fluxes of right- and left-moving waves
      real, dimension(nvar%all,n)    :: dfrp               !< second order correction of right-moving waves flux on the right cell boundary
      real, dimension(nvar%all,n)    :: dfrm               !< second order correction of right-moving waves flux on the left cell boundary
      real, dimension(nvar%all,n)    :: dflm               !< second order correction of left-moving waves flux on the left cell boundary
      real, dimension(nvar%all,n)    :: dflp               !< second order correction of left-moving waves flux on the right cell boundary
      real, dimension(nvar%all,n)    :: u1                 !< updated vector of conservative variables (after one timestep in second order scheme)
      real, dimension(nvar%fluids,n) :: rotacc             !< acceleration caused by rotation
      real, dimension(nvar%fluids,n) :: fricacc            !< acceleration caused by friction
      real, dimension(nvar%fluids,n) :: geosrc             !< source terms caused by geometry of coordinate system
      real, dimension(nvar%fluids,n) :: pressure           !< gas pressure
      real, dimension(2)             :: df                 !< marker
      real, dimension(n)             :: gravacc            !< acceleration caused by gravitation
#ifdef ISO_LOCAL
      real, dimension(n)             :: cs_iso2            !< square of local isothermal sound speed (optional for ISO_LOCAL)
#endif /* ISO_LOCAL */

#ifdef SHEAR
      real, dimension(nvar%fluids,n) :: vy0
#endif /* SHEAR */

#ifdef COSM_RAYS
       integer                       :: icr
#ifdef COSM_RAYS_SOURCES
       real                          :: srccrn(nvar%crn%all,n)
#endif /* COSM_RAYS_SOURCES */
#endif /* COSM_RAYS */

#ifndef ISO
      real, dimension(nvar%fluids,n) :: ekin,eint
#if defined IONIZED && defined MAGNETIC
      real, dimension(n)             :: emag
#endif /* IONIZED && MAGNETIC */
#endif /* !ISO */
#ifdef COSM_RAYS
      real, dimension(n)             :: divv,grad_pcr,ecr
      real, dimension(n)             :: decr
#endif /* COSM_RAYS */
#ifdef FLUID_INTERACTIONS
      real, dimension(nvar%fluids,n) :: epsa, vx0
#endif /* FLUID_INTERACTIONS */
#ifdef FLUID_INTERACTIONS_DW
      real, dimension(nvar%all,n)    :: dintr
#endif /* FLUID_INTERACTIONS_DW */

      real, dimension(2,2), parameter:: rk2coef = RESHAPE( (/1.0,0.5,0.0,1.0/),(/2,2/))

#if !defined(ISO_LOCAL) && !defined(GRAV) && ! defined(FLUID_INTERACTIONS_DW) && !(defined COSM_RAYS && defined IONIZED)
      integer                        :: dummy
      if (.false.) dummy = i1 + i2 ! suppress compiler warnings on unused arguments
#endif /* !ISO_LOCAL && !GRAV && !FLUID_INTERACTIONS_DW && !(COSM_RAYS && IONIZED) */
#if !defined(ISO_LOCAL) && !defined(SHEAR) && !defined(GRAV) && ! defined(FLUID_INTERACTIONS_DW) && !(defined COSM_RAYS && defined IONIZED)
      if (.false.) dummy = len(sweep)
#endif /* !ISO_LOCAL && !SHEAR && !GRAV && !FLUID_INTERACTIONS_DW && !(COSM_RAYS && IONIZED) */

      w         = 0.0
      cfr       = 0.0
      dtx       = dt / dx

      u1 = u
      u0 = u

#ifdef ISO_LOCAL
      if (sweep .eq. 'xsweep') then
         cs_iso2(:) =  cs_iso2_arr(:,i1,i2)
      else if (sweep .eq. 'ysweep')  then
         cs_iso2(:) =  cs_iso2_arr(i2,:,i1)
      else
         cs_iso2(:) =  cs_iso2_arr(i1,i2,:)
      endif
#endif /* ISO_LOCAL */

      do istep=1,integration_order

! Fluxes calculation for cells centers
#ifdef ISO_LOCAL
         call all_fluxes(n, w, cfr, u1, bb, cs_iso2)
#else /* !ISO_LOCAL */
         call all_fluxes(n, w, cfr, u1, bb, pressure)
#endif /* !ISO_LOCAL */
! Right and left fluxes decoupling

         fl = (u1*cfr-w)*0.5
!        fr = (u1*cfr+w)*0.5
         fr = fl + w

         fl(:,1:n-1) = fl(:,2:n)                         ; fl(:,n)   = fl(:,n-1)

         if (istep == 2) then

! Second order flux corrections

            dfrp(:,1:n-1) = 0.5*(fr(:,2:n) - fr(:,1:n-1)); dfrp(:,n) = dfrp(:,n-1)
            dfrm(:,2:n)   = dfrp(:,1:n-1);                 dfrm(:,1) = dfrm(:,2)

! Flux limiter application

            call flimiter(fr,dfrm,dfrp)

            dflp(:,1:n-1) = 0.5*(fl(:,1:n-1) - fl(:,2:n)); dflp(:,n) = dflp(:,n-1)
            dflm(:,2:n)   = dflp(:,1:n-1);                 dflm(:,1) = dflm(:,2)
            call flimiter(fl,dflm,dflp)
         endif

! u update

         fu = fr - fl
         u1(:,2:n) = u0(:,2:n) - rk2coef(integration_order,istep) * gc(1,:,2:n) * dtx * ( gc(2,:,2:n)*fu(:,2:n) - gc(3,:,2:n)*fu(:,1:n-1) )
         u1(:,1)   = u1(:,2)

         if (use_smalld) then
            ! This is needed e.g. for outflow boundaries in presence of perp. gravity
            u1(iarr_all_dn,:) = max(u1(iarr_all_dn,:),smalld)
            local_magic_mass = local_magic_mass + sum( abs(u1(iarr_all_dn,2:n-1) - u0(iarr_all_dn,2:n-1)) )
         else
            if (any(u1(iarr_all_dn,:) < 0.0)) then
               write(msg,'(3A,I4,1X,I4,A)') "[rtvd:relaxing_tvd] negative density in sweep ",sweep,"( ", i1, i2, " )"
               call die(msg)
            endif
         endif

! Source terms -------------------------------------
         geosrc = geometry_source_terms(u,pressure,sweep)
#ifndef GRAV
         u1(iarr_all_mx,:) = u1(iarr_all_mx,:) + rk2coef(integration_order,istep)*geosrc(:,:)*dt
#endif /* !GRAV */

#ifdef FLUID_INTERACTIONS
#ifdef SHEAR
         df = global_gradP
#else /* !SHEAR */
         df = 0.0
#endif /* !SHEAR */
         epsa(1,:) = dragc_gas_dust * u(iarr_all_dn(2),:)  / u(iarr_all_dn(1),:)
         epsa(2,:) = dragc_gas_dust
         where (u(iarr_all_dn,:) > 0.0)
            vx0(:,:)  = u(iarr_all_mx,:)/u(iarr_all_dn,:)
         elsewhere
            vx0(:,:)  = 0.0
         endwhere

         do ind = 1, nvar%fluids
            if (ind == 1) then
               fricacc(ind,:) = - epsa(ind,:) * (vx0(1,:) - vx0(2,:))
            else
               fricacc(ind,:) = - epsa(ind,:) * (vx0(2,:) - vx0(1,:))
            endif
         enddo

#else /* !FLUID_INTERACTIONS */
         fricacc(:,:) = 0.0
         df = 0.0
#endif /* !FLUID_INTERACTIONS */

#ifdef SHEAR
         where (u(iarr_all_dn,:) > 0.0)
            vy0(:,:)  = u(iarr_all_my,:)/u(iarr_all_dn,:)
         elsewhere
            vy0(:,:)  = 0.0
         endwhere
         do ind = 1, nvar%fluids
!            if (sweep .eq. 'xsweep') then
!               rotacc(ind,:) =  2.0*omega*(vy0(ind,:) + qshear*omega*x(:))
!            else if (sweep .eq. 'ysweep')  then
!               rotacc(ind,:) = - 2.0*omega*vy0(ind,:)          ! with global shear
!            else
!               rotacc(ind,:) = 0.0
!            endif
            if (sweep .eq. 'xsweep') then
               rotacc(ind,:) =  2.0*omega*vy0(ind,:) + df(ind)  ! global_gradient
            else if (sweep .eq. 'ysweep')  then
               rotacc(ind,:) = (qshear - 2.0)*omega*vy0(ind,:)  ! with respect to global shear (2.5D)
            else
               rotacc(ind,:) = 0.0
            endif
         enddo
#else /* !SHEAR */
         rotacc(:,:) = 0.0
#endif /* !SHEAR */

#ifdef GRAV
         call grav_pot2accel(sweep,i1,i2, n, gravacc, istep)
#else /* !GRAV */
         gravacc = 0.0
#endif /* !GRAV */

#if defined GRAV || defined SHEAR || defined FLUID_INTERACTIONS
         acc     =  rotacc + fricacc
         do ind = 1, nvar%fluids
            acc(ind,:) =  acc(ind,:) + gravacc(:)
         enddo

         acc(:,n)   = acc(:,n-1); acc(:,1) = acc(:,2)

         u1(iarr_all_mx,:) = u1(iarr_all_mx,:) + rk2coef(integration_order,istep)*(acc(:,:)*u(iarr_all_dn,:)+geosrc(:,:))*dt
#ifndef ISO
         u1(iarr_all_en,:) = u1(iarr_all_en,:) + rk2coef(integration_order,istep)*acc(:,:)*u(iarr_all_mx,:)*dt
#endif /* !ISO */

#endif /* defined GRAV || defined SHEAR || defined FLUID_INTERACTIONS */

#ifdef FLUID_INTERACTIONS_DW
         !!!! old, more general way of adding interactions
         call fluid_interactions(sweep,i1,i2, n, dintr, u)
         u1 = u1 + rk2coef(integration_order,istep)*dintr*dt
#endif /* FLUID_INTERACTIONS_DW */

#if defined COSM_RAYS && defined IONIZED
         select case (sweep)
            case ('xsweep')
               divv = divvel(:,i1,i2)
            case ('ysweep')
               divv = divvel(i2,:,i1)
            case ('zsweep')
               divv = divvel(i1,i2,:)
         end select

         grad_pcr(:) = 0
         if (nvar%crn%all > 0) then ! BEWARE: quick hack
            do icr = 1, 1 !nvar_crs  !<BEWARE TEMPORARY!
               decr(:)                = -(gamma_crs(icr)-1.)*u1(iarr_crs(icr),:)*divv(:)*dt
               u1  (iarr_crs(icr),:)  = u1(iarr_crs(icr),:) + rk2coef(integration_order,istep)*decr(:)
               u1  (iarr_crs(icr),:)  = max(smallecr,u1(iarr_crs(icr),:))

               ecr                    = u1(iarr_crs(icr),:)
               grad_pcr(2:n-1) = grad_pcr(2:n-1) + cr_active*(gamma_crs(icr) -1.)*(ecr(3:n)-ecr(1:n-2))/(2.*dx)

            enddo
         endif
         grad_pcr(1:2)   = 0.0 ; grad_pcr(n-1:n) = 0.0

#ifndef ISO
         u1(iarr_all_en(nvar%ion%pos),:) = u1(iarr_all_en(nvar%ion%pos),:) &
                              - rk2coef(integration_order,istep)*u1(iarr_all_mx(nvar%ion%pos),:)/u1(iarr_all_dn(nvar%ion%pos),:)*grad_pcr*dt
#endif /* !ISO */
         u1(iarr_all_mx(nvar%ion%pos),:) = u1(iarr_all_mx(nvar%ion%pos),:) - rk2coef(integration_order,istep)*grad_pcr*dt

#ifdef COSM_RAYS_SOURCES
         call src_crn(u1,n, srccrn)
         u1  (iarr_crn,:)  = u1(iarr_crn,:) + rk2coef(integration_order,istep)*srccrn(:,:)*dt

#endif /* COSM_RAYS_SOURCES */

#endif /* COSM_RAYS && IONIZED */

#if defined IONIZED || defined NEUTRAL
#ifndef ISO
         ekin = 0.5*( u1(iarr_all_mx,:)**2 + u1(iarr_all_my,:)**2 &
                +u1(iarr_all_mz,:)**2) /u1(iarr_all_dn,:)
         eint = u1(iarr_all_en,:)-ekin
#if defined IONIZED && defined MAGNETIC
         emag = 0.5*(bb(ibx,:)*bb(ibx,:) + bb(iby,:)*bb(iby,:) + bb(ibz,:)*bb(ibz,:))
         eint(nvar%ion%pos,:) = eint(nvar%ion%pos,:) - emag
#endif /* IONIZED && MAGNETIC */

         eint = max(eint,smallei)

         u1(iarr_all_en,:) = eint+ekin
#if defined IONIZED && defined MAGNETIC
         u1(iarr_all_en(nvar%ion%pos),:) = u1(iarr_all_en(nvar%ion%pos),:)+emag
#endif /* IONIZED && MAGNETIC */
#endif /* !ISO */
#endif /*  defined IONIZED || defined NEUTRAL  */

         u(:,:) = u1(:,:)
      enddo

      return

   end subroutine relaxing_tvd

!==========================================================================================
end module rtvd
