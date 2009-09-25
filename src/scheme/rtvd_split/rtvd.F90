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
module rtvd ! split orig

   contains

   subroutine tvdb(vibj,b,vg,n,dt,di)
      use constants, only : big
      implicit none
      integer, intent(in) :: n
      real, intent(in)    :: dt,di
      real, dimension(n)  :: vibj,b,vg
! locals
      real, dimension(n)  :: b1,vibj1,vh
      real :: dti, v, w, dw, dwm, dwp
      integer :: i, ip, ipp, im

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

      integer                   :: i1,i2, n, istep, ind
      character(len=6)          :: sweep

      real                      :: dt,dx,dtx
      real, dimension(nvar,n)   :: u,cfr
      real, dimension(nmag,n)   :: bb
      real, dimension(nfluid,n) :: accr
!locals
      real, dimension(nvar,n)   :: w,fr,fl,dfrp,dfrm,dflm,dflp,dulf,durf
      real, dimension(nvar,n)   :: ul0,ur0,u1,ul1,ur1

      real, dimension(nfluid,n) :: rotaccr, fricaccr
      real, dimension(2)        :: df
      real, dimension(n)        :: gravaccr

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

         call all_fluxes(w,cfr,u1,bb,n)

         fr = (u1*cfr+w)*0.5
         fl = (u1*cfr-w)*0.5
         if(istep == 1) then
            ur0 = fr/cfr
            ul0 = fl/cfr
         endif

         fl(:,1:n-1) = fl(:,2:n)                         ; fl(:,n)   = fl(:,n-1)

         if(istep == 2) then
            dfrp(:,1:n-1) = 0.5*(fr(:,2:n) - fr(:,1:n-1)); dfrp(:,n) = dfrp(:,n-1)
            dfrm(:,2:n)   = dfrp(:,1:n-1);                 dfrm(:,1) = dfrm(:,2)
            call flimiter(fr,dfrm,dfrp,nvar,n)

            dflp(:,1:n-1) = 0.5*(fl(:,1:n-1) - fl(:,2:n)); dflp(:,n) = dflp(:,n-1)
            dflm(:,2:n)   = dflp(:,1:n-1);                 dflm(:,1) = dflm(:,2)
            call flimiter(fl,dflm,dflp,nvar,n)
         endif

         durf(:,2:n) = dtx*(fr(:,2:n) - fr(:,1:n-1));     durf(:,1) = durf(:,2)
         dulf(:,2:n) = dtx*(fl(:,2:n) - fl(:,1:n-1));     dulf(:,1) = dulf(:,2)

         ur1(:,:) = ur0 - rk2coef(integration_order,istep)*durf
         ul1(:,:) = ul0 + rk2coef(integration_order,istep)*dulf

         u1 = ul1 + ur1
         u1(iarr_all_dn(1),:) = max(u1(iarr_all_dn(1),:), smalld)

! Source terms -------------------------------------
#ifdef FLUID_INTERACTIONS
#ifdef SHEAR
         df = (/global_gradP_neu,0.0/)
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
               fricaccr(ind,:) = - epsa(ind,:) * (vx0(1,:) - vx0(2,:))
            else
               fricaccr(ind,:) = - epsa(ind,:) * (vx0(2,:) - vx0(1,:))
            endif
         enddo

#else /* FLUID_INTERACTIONS */
         fricaccr(:,:) = 0.0
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
!               rotaccr(ind,:) =  2.0*omega*(vy0(ind,:) + qshear*omega*x(:))
!            else if(sweep .eq. 'ysweep')  then
!               rotaccr(ind,:) = - 2.0*omega*vy0(ind,:)          ! with global shear
!            else
!               rotaccr(ind,:) = 0.0
!            endif
            if(sweep .eq. 'xsweep') then
               rotaccr(ind,:) =  2.0*omega*vy0(ind,:) + df(ind)  ! global_gradient
            else if(sweep .eq. 'ysweep')  then
               rotaccr(ind,:) = (qshear - 2.0)*omega*vy0(ind,:)  ! with respect to global shear (2.5D)
            else
               rotaccr(ind,:) = 0.0
            endif
         enddo
#else  /* SHEAR */
         rotaccr(:,:) = 0.0
#endif /* SHEAR */

#ifdef GRAV
         call grav_pot2accel(sweep,i1,i2, n, gravaccr)
#else /* GRAV */
         gravaccr = 0.0
#endif /* GRAV */


#if defined GRAV || defined SHEAR || defined FLUID_INTERACTIONS
         accr     =  rotaccr + fricaccr
         do ind = 1, nfluid
            accr(ind,:) =  accr(ind,:) + gravaccr(:)
         enddo

         accr(:,n)   = accr(:,n-1); accr(:,1) = accr(:,2)

         !!!! BEWARE: May not be necessary anymore
         where(u1(iarr_all_dn,:) < 0.0)
            accr(:,:) = 0.0
         endwhere

         u1(iarr_all_mx,:) = u1(iarr_all_mx,:) + rk2coef(integration_order,istep)*accr(:,:)*u(iarr_all_dn,:)*dt
#ifndef ISO
         u1(iarr_all_en,:) = u1(iarr_all_en,:) + rk2coef(integration_order,istep)*accr(:,:)*u(iarr_all_mx,:)*dt
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
