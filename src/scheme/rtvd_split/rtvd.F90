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
!    Initial implemetation of PIERNIK code was based on TVD split MHD code by
!    Ue-Li Pen
!        see: Pen, Arras & Wong (2003) for algorithm and
!             http://www.cita.utoronto.ca/~pen/MHD
!             for original source code "mhd.f90"
!
!    For full list of developers see $PIERNIK_HOME/license/pdt.txt
!
#include "piernik.def"
!>
!! \brief (JD) This module implements relaxing TVD scheme (doxy comments ready)
!! 
!! The implementation was based on TVD split MHD code by Pen et al. (2003).
!<
module rtvd ! split orig

   contains
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
!! In a discrete representation variations of magnetic fluxes, in timestep \f$\Delta t\f$, threading faces of cell \f$(i,j,k)\f$ are
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
!! A particular implementation of the CT scheme depends on actual centering of fluid variables. In the RTVD scheme all fluid variables  (\f$\rho\f$, 
!! \f$m_x\f$, \f$m_y\f$, \f$m_z\f$ and \f$e\f$) are  cell-centered, and  magnetic field components (\f$B_x\f$, \f$B_y\f$, \f$B_z\f$) are 
!! face-centered. For that reason interpolation is necessary. For the detailed scheme used in our code see sections prepared for Lecture Notes in 
!! Physics, Springer, that are placed in doc/evora subdirectory. 
!<
   subroutine tvdb(vibj,b,vg,n,dt,di)
      use constants, only : big
      implicit none
      integer, intent(in) :: n       !< array size
      real, intent(in)    :: dt      !< time step
      real, intent(in)    :: di      !< cell length, depends on direction x, y or z
      real, dimension(n)  :: vibj    !< face-centered electromotive force components (b*vg)
      real, dimension(n)  :: b       !< magnetic field
      real, dimension(n)  :: vg      !< velocity in the center of cell boundary
! locals
      real, dimension(n)  :: b1      !< magnetic field
      real, dimension(n)  :: vibj1   !< face-centered electromotive force (EMF) components (b*vg)
      real, dimension(n)  :: vh      !< velocity interpolated to the cell edges 
      real :: dti		     !< dt/di
      real :: v			     !< auxiliary variable to compute EMF
      real :: w			     !< EMF component
      real :: dw                     !< The second-order correction to EMF component 
      real :: dwm                    !< face centered EMF interpolated to left cell-edge
      real :: dwp                    !< face centered EMF interpolated to right cell-edge
      integer :: i                   !< auxiliary array indicator
      integer :: ip                  !< i+1
      integer :: ipp                 !< i+2
      integer :: im                  !< i-1

  ! unlike the B field, the vibj lives on the right cell boundary
      vh = 0.0
      vh(1:n-1) =(vg(1:n-1)+ vg(2:n))*0.5;     vh(n) = vh(n-1)

      dti = dt/di
      where(vh > 0.)
         vibj1=b*vg
      elsewhere
         vibj1=eoshift(b*vg,1,boundary=big)
      end where

      b1(2:n) = b(2:n) -(vibj1(2:n)-vibj1(1:n-1))*dti*0.5;    b1(1) = b(2)

      do i = 3, n-3
         ip  = i  + 1
         ipp = ip + 1
         im  = i  - 1
         v   = vh(i)
         if (v > 0.0) then
            w=vg(i)*b1(i)
            dwp=(vg(ip)*b1(ip)-w)*0.5
            dwm=(w-vg(im)*b1(im))*0.5
         else
            w=vg(ip)*b1(ip)
            dwp=(w-vg(ipp)*b1(ipp))*0.5
            dwm=(vg(i)*b1(i)-w)*0.5
         end if
         dw=0.0
         if(dwm*dwp > 0.0) dw=2.0*dwm*dwp/(dwm+dwp)
         vibj(i)=(w+dw)*dt
      enddo

   end subroutine tvdb

!>
!! \brief Subroutine implements the Relaxing TVD scheme for conserved physical quantities
!!
!! The main point of the Relaxing TVD scheme is to decompose vectors of conservative variables \f$u\f$ and fluxes \f$F\f$ 
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
!! The fluxes then can be written as
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
!! \n (1) First order fluxes \f$\vec{F}_{i\pm 1/2}^{(1)L,R}\f$   are 
!!    calculated together with source terms \f$\vec{S}_i(\vec{u})\f$ at \f$t^n\f$.
!! \n (2) Fluxes and source terms derived in the first step are used to calculate 
!!    \f$\vec{u}^{n+1/2}\f$ at \f$t^{n+1/2}\f$.
!! \n (3) Second order fluxes \f$\vec{F}_{i+1/2}^{(2)}\f$, obtained via monotonic, upwind
!!    interpolation to cell boundaries,  and source terms \f$\vec{S}\f$ are evaluated for 
!!    \f$\vec{u}^{n+1/2}\f$ at \f$t^{n+1/2}\f$.
!! \n (4) Update of \f$\Delta\vec{u}\f$, corresponding to the full 
!!    timestep \f$\Delta t\f$ is done, using fluxes and source terms calculated 
!!    at the intermediate time step \f$t^{n+1/2}\f$. 
!!
!! To achieve second order spatial accuracy, a monotone upwind interpolation of fluxes onto cell boundaries is made, 
!! with the aid of a flux limiter, to obtain the Monotone Upwind Scheme for Conservation Laws (MUSCL) (step (3) in Runge-Kutta scheme).
!! The monotone interpolation is used to avoid spurious oscillations in discrete solutions of higher order.
!! The Total Variation Diminishing property of the numerical scheme is related to the measure of the overall amount
!! of oscillations, called Total Variation.
!<
   subroutine relaxing_tvd(u,bb,sweep,i1,i2,dx,n,dt)

#ifdef IONIZED
      use fluidindex,      only : i_ion
#endif /* IONIZED */

      use mpisetup,        only : smalld, integration_order
      use fluxes,          only : flimiter,all_fluxes
      use fluidindex,      only : nvar,nmag,nfluid
      use fluidindex,      only : ibx,iby,ibz
      use fluidindex,      only : iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz
#ifndef ISO
      use fluidindex,      only : iarr_all_en
      use mpisetup,        only : smallei
#endif /* ISO */

#ifdef GRAV
      use gravity,         only : grav_pot2accel
#endif /* GRAV */
#ifdef SHEAR
      use grid,            only : x
      use shear,           only : qshear, omega
#ifdef NEUTRAL
      use initneutral,     only : global_gradP_neu
#endif /* NEUTRAL */
#endif /* SHEAR */
#ifdef COSM_RAYS
      use initcosmicrays,  only : gamma_cr, cr_active, smallecr
      use initcosmicrays,  only : iecr
      use arrays,          only : divvel
#endif /* COSM_RAYS */
#ifdef FLUID_INTERACTIONS
      use initdust,        only : dragc_gas_dust 
      use interactions
#endif /* FLUID_INTERACTIONS */

      implicit none

      integer                   :: i1                 !< coordinate of sweep in the 1st remaining direction
      integer                   :: i2                 !< coordinate of sweep in the 2nd remaining direction
      integer                   :: n                  !< array size
      integer                   :: istep              !< step number in the time integration scheme
      integer                   :: ind                !< fluid index
      character(len=6)          :: sweep              !< direction (x, y or z) we are doing calculations for

      real                      :: dt                 !< time step
      real                      :: dx                 !< cell length
      real                      :: dtx                !< dt/dx
      real, dimension(nvar,n)   :: u                  !< vector of conservative vatiables
      real, dimension(nvar,n)   :: cfr                !< freezing speed
      real, dimension(nmag,n)   :: bb                 !< local copy of magnetic field
      real, dimension(nfluid,n) :: acc                !< acceleration
!locals
      real, dimension(nvar,n)   :: w                  !< auxiliary vector to calculate fluxes
      real, dimension(nvar,n)   :: fr                 !< flux of the right-moving waves
      real, dimension(nvar,n)   :: fl                 !< flux of the left-moving waves
      real, dimension(nvar,n)   :: dfrp               !< second order correction of right-moving waves flux on the right cell boundary
      real, dimension(nvar,n)   :: dfrm               !< second order correction of right-moving waves flux on the left cell boundary
      real, dimension(nvar,n)   :: dflm               !< second order correction of left-moving waves flux on the left cell boundary
      real, dimension(nvar,n)   :: dflp               !< second order correction of left-moving waves flux on the right cell boundary
      real, dimension(nvar,n)   :: dulf               !< second order correction of the vector of conservative variables for the left-moving waves
      real, dimension(nvar,n)   :: durf               !< second order correction of the vector of conservative variables for the right-moving waves
      real, dimension(nvar,n)   :: ul0                !< left moving wave in first order scheme
      real, dimension(nvar,n)   :: ur0                !< right moving wave in first order scheme
      real, dimension(nvar,n)   :: u1                 !< uptaded vector of conservative variables (after one timestep in second order scheme)
      real, dimension(nvar,n)   :: ul1                !< left moving wave (after one timestep in second order scheme)
      real, dimension(nvar,n)   :: ur1                !< right moving wave (after one timestep in second order scheme)
      real, dimension(nfluid,n) :: rotacc             !< acceleration caused by rotation
      real, dimension(nfluid,n) :: fricacc            !< acceleration caused by friction
      real, dimension(2)        :: df                 !< marker                 
      real, dimension(n)        :: gravacc            !< acceleration caused by gravitation

#ifdef SHEAR
      real, dimension(nfluid,n) :: vy0
#endif /* SHEAR */

#ifdef COSM_RAYS
      real, dimension(n)       :: vx
#endif /* COSM_RAYS */

#ifndef ISO
      real, dimension(nfluid,n):: ekin,eint
      real, dimension(n)       :: emag
#endif /* ISO */
#ifdef COSM_RAYS
      real, dimension(n)       :: divv,decr,grad_pcr,ecr
#endif /* COSM_RAYS */
#ifdef FLUID_INTERACTIONS
      real, dimension(nvar,n)  :: dintr
      real, dimension(nfluid,n):: epsa, vx0
#endif /* FLUID_INTERACTIONS */

      real, dimension(2,2), parameter  :: rk2coef = RESHAPE( (/1.0,0.5,0.0,1.0/),(/2,2/)) 

      w         = 0.0
      cfr       = 0.0
      dtx       = dt / dx

      u1 = u

      do istep=1,integration_order

! Fluxes calculation for cells centers

         call all_fluxes(w,cfr,u1,bb,n)

! Right and left fluxes decoupling

         fr = (u1*cfr+w)*0.5
         fl = (u1*cfr-w)*0.5

         if(istep == 1) then

! Right moving waves construction for first-order scheme

            ur0 = fr/cfr

! Left moving waves construction for first-order scheme

            ul0 = fl/cfr
         endif

         fl(:,1:n-1) = fl(:,2:n)                         ; fl(:,n)   = fl(:,n-1)

         if(istep == 2) then

! Second order flux corrections

            dfrp(:,1:n-1) = 0.5*(fr(:,2:n) - fr(:,1:n-1)); dfrp(:,n) = dfrp(:,n-1)
            dfrm(:,2:n)   = dfrp(:,1:n-1);                 dfrm(:,1) = dfrm(:,2)

! Flux limiter application

            call flimiter(fr,dfrm,dfrp,nvar,n)

            dflp(:,1:n-1) = 0.5*(fl(:,1:n-1) - fl(:,2:n)); dflp(:,n) = dflp(:,n-1)
            dflm(:,2:n)   = dflp(:,1:n-1);                 dflm(:,1) = dflm(:,2)
            call flimiter(fl,dflm,dflp,nvar,n)
         endif

! u corrections

         durf(:,2:n) = dtx*(fr(:,2:n) - fr(:,1:n-1));     durf(:,1) = durf(:,2)
         dulf(:,2:n) = dtx*(fl(:,2:n) - fl(:,1:n-1));     dulf(:,1) = dulf(:,2)

! u update

         ur1(:,:) = ur0 - rk2coef(integration_order,istep)*durf
         ul1(:,:) = ul0 + rk2coef(integration_order,istep)*dulf

         u1 = ul1 + ur1
         u1(iarr_all_dn(1),:) = max(u1(iarr_all_dn(1),:), smalld)

! Source terms -------------------------------------
#ifdef FLUID_INTERACTIONS
#ifdef SHEAR
         df = (/global_gradP_neu,0.0/)        ! znacznik1
#else 
         df = 0.0
#endif
         epsa(1,:) = dragc_gas_dust * u(iarr_all_dn(2),:)  / u(iarr_all_dn(1),:)
         epsa(2,:) = dragc_gas_dust
         where(u(iarr_all_dn,:) > 0.0) 
            vx0(:,:)  = u(iarr_all_mx,:)/u(iarr_all_dn,:)
         elsewhere
            vx0(:,:)  = 0.0
         endwhere


         do ind = 1, nfluid
            if(ind == 1) then
               fricacc(ind,:) = - epsa(ind,:) * (vx0(1,:) - vx0(2,:))
            else
               fricacc(ind,:) = - epsa(ind,:) * (vx0(2,:) - vx0(1,:))
            endif
         enddo

#else /* FLUID_INTERACTIONS */
         fricacc(:,:) = 0.0
         df = 0.0
#endif /* FLUID_INTERACTIONS */

#ifdef SHEAR
         where(u(iarr_all_dn,:) > 0.0) 
            vy0(:,:)  = u(iarr_all_my,:)/u(iarr_all_dn,:)
         elsewhere
            vy0(:,:)  = 0.0
         endwhere
         do ind = 1, nfluid
!            if(sweep .eq. 'xsweep') then
!               rotacc(ind,:) =  2.0*omega*(vy0(ind,:) + qshear*omega*x(:))
!            else if(sweep .eq. 'ysweep')  then
!               rotacc(ind,:) = - 2.0*omega*vy0(ind,:)          ! with global shear
!            else
!               rotacc(ind,:) = 0.0
!            endif
            if(sweep .eq. 'xsweep') then
               rotacc(ind,:) =  2.0*omega*vy0(ind,:) + df(ind)  ! global_gradient
            else if(sweep .eq. 'ysweep')  then
               rotacc(ind,:) = (qshear - 2.0)*omega*vy0(ind,:)  ! with respect to global shear (2.5D)
            else
               rotacc(ind,:) = 0.0
            endif
         enddo
#else  /* SHEAR */
         rotacc(:,:) = 0.0
#endif /* SHEAR */

#ifdef GRAV
         call grav_pot2accel(sweep,i1,i2, n, gravacc)
#else /* GRAV */
         gravacc = 0.0
#endif /* GRAV */


#if defined GRAV || defined SHEAR || defined FLUID_INTERACTIONS
         acc     =  rotacc + fricacc
         do ind = 1, nfluid
            acc(ind,:) =  acc(ind,:) + gravacc(:)
         enddo

         acc(:,n)   = acc(:,n-1); acc(:,1) = acc(:,2)

         !!!! BEWARE: May not be necessary anymore
         where(u1(iarr_all_dn,:) < 0.0)
            acc(:,:) = 0.0
         endwhere

         u1(iarr_all_mx,:) = u1(iarr_all_mx,:) + rk2coef(integration_order,istep)*acc(:,:)*u(iarr_all_dn,:)*dt
#ifndef ISO
         u1(iarr_all_en,:) = u1(iarr_all_en,:) + rk2coef(integration_order,istep)*acc(:,:)*u(iarr_all_mx,:)*dt
#endif /* ISO */

#endif /* defined GRAV || defined SHEAR || defined FLUID_INTERACTIONS */

#ifdef FLUID_INTERACTIONS_DW
         !!!! old, more general way of adding interactions
         call fluid_interactions(sweep,i1,i2, n, dintr, u)
         u1 = u1 + rk2coef(integration_order,istep)*dintr*dt
#endif /* FLUID_INTERACTIONS_DW */

#if defined COSM_RAYS && defined IONIZED
         select case (sweep)
            case('xsweep')
               divv = divvel(:,i1,i2)
            case('ysweep')
               divv = divvel(i2,:,i1)
            case('zsweep')
               divv = divvel(i1,i2,:)
         end select

         decr(:)    = -(gamma_cr-1.)*u1(iecr,:)*divv(:)*dt
         u1(iecr,:) = u1(iecr,:) + rk2coef(integration_order,istep)*decr(:)
         u1(iecr,:) = max(smallecr,u1(iecr,:))

         vx  = u1(iarr_all_mx(i_ion),:)/u1(iarr_all_dn(i_ion),:)
         ecr = u1(iecr,:)

         grad_pcr(2:n-1) = cr_active*(gamma_cr -1.)*(ecr(3:n)-ecr(1:n-2))/(2.*dx)
         grad_pcr(1:2)=0.0 ; grad_pcr(n-1:n) = 0.0

#ifndef ISO
         u1(iarr_all_en(i_ion),:) = u1(iarr_all_en(i_ion),:) &
                              - rk2coef(integration_order,istep)*u1(iarr_all_mx(i_ion),:)/u1(iarr_all_dn(i_ion),:)*grad_pcr*dt
#endif /* ISO */
         u1(iarr_all_mx(i_ion),:) = u1(iarr_all_mx(i_ion),:) - rk2coef(integration_order,istep)*grad_pcr*dt

#endif /* COSM_RAYS && IONIZED */

#ifndef ISO
         ekin = 0.5*( u1(iarr_all_mx,:)**2 + u1(iarr_all_my,:)**2 &
                +u1(iarr_all_mz,:)**2) /u1(iarr_all_dn,:)
         eint = u1(iarr_all_en,:)-ekin
#if defined IONIZED && defined MAGNETIC
         emag = 0.5*(bb(ibx,:)*bb(ibx,:) + bb(iby,:)*bb(iby,:) + bb(ibz,:)*bb(ibz,:))
         eint(i_ion,:) = eint(i_ion,:) - emag
#endif /* IONIZED && MAGNETIC */

         eint = max(eint,smallei)

         u1(iarr_all_en,:) = eint+ekin
#if defined IONIZED && defined MAGNETIC
         u1(iarr_all_en(i_ion),:) = u1(iarr_all_en(i_ion),:)+emag
#endif /* IONIZED && MAGNETIC */
#endif /* ISO */

         u(:,:) = u1(:,:)
      enddo

      return

   end subroutine relaxing_tvd

!==========================================================================================
end module rtvd
