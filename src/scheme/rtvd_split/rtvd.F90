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
module rtvd ! split orig

! pulled by ANY

   implicit none

   private
   public  :: relaxing_tvd

contains
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
!!
!! \todo Do not pass i1 and i2, pass optional pointer to gravacc instead
!<
!*/

! OPT: 15% of CPU time spent in the relaxing_tvd routine is attributed to the entry point. The rest is more or less evenly distributed across all array operations
! OPT: \todo try to pass pointers instead of arrays, or assemble the arrays here
! OPT: n is usually short enough for all the data to fit L2 cache (checked on 512kB)
! OPT: we may also try to work on bigger parts of the u(:,:,:,:) at a time , but the exact amount may depend on size of the L2 cache
! OPT: try an explicit loop over n to see if better pipelining can be achieved

   subroutine relaxing_tvd(n, u0, u1, vel_sweep, bb, cs_iso2, istep, dtx, eflx)

      use constants,    only: half, GEO_XYZ, RK2_1, RK2_2
      use dataio_pub,   only: die
      use domain,       only: dom
      use fluidindex,   only: flind, nmag
      use fluxes,       only: flimiter, all_fluxes
      use fluxtypes,    only: ext_fluxes
      use gridgeometry, only: gc, GC1, GC2, GC3

      implicit none

      integer(kind=4),                  intent(in)    :: n                  !< array size
      real, dimension(n, flind%all),    intent(in)    :: u0                 !< vector of conservative variables
      real, dimension(n, flind%all),    intent(inout) :: u1                 !< vector of conservative variables
      real, dimension(n, flind%fluids), intent(in)    :: vel_sweep          !< velocity in the direction of current sweep
      real, dimension(n, nmag),         intent(in)    :: bb                 !< local copy of magnetic field
      real, dimension(:), pointer,      intent(in)    :: cs_iso2            !< square of local isothermal sound speed
      integer,                          intent(in)    :: istep              !< stage in the time integration scheme
      real,                             intent(in)    :: dtx                !< RK_coeff * time step / dx
      type(ext_fluxes),                 intent(inout) :: eflx               !< external fluxes

      real, dimension(n, flind%all)                   :: cfr                !< freezing speed
!locals
      real, dimension(n, flind%all)                   :: w                  !< auxiliary vector to calculate fluxes
      real, dimension(n, flind%all)                   :: fr                 !< flux of the right-moving waves
      real, dimension(n, flind%all)                   :: fl                 !< flux of the left-moving waves
      real, dimension(n, flind%all)                   :: fu                 !< sum of fluxes of right- and left-moving waves
      real, dimension(n, flind%all)                   :: dfp                !< second order correction of left/right-moving waves flux on the right cell boundary
      real, dimension(n, flind%all)                   :: dfm                !< second order correction of left/right-moving waves flux on the left cell boundary
      logical                                         :: full_dim

      !OPT: try to avoid these explicit initializations of u0(:,:)
      full_dim = n > 1

      if (full_dim) then
         ! Fluxes calculation for cells centers
         call all_fluxes(n, w, cfr, u1, bb, vel_sweep, cs_iso2)
         ! Right and left fluxes decoupling

         ! original code
         ! fl(1:n-1, :) = (cfr(2:n, :)*u1(2:n, :) - wl(2:n, :)) * 0.5
         ! fr          = (cfr*u1 + w) * 0.5
         ! following is equivalent but faster

         fl(1:n-1, :) = (u1(2:n, :) * cfr(2:n, :) - w(2:n, :)) * half
         fl(n, :) = fl(n-1, :)
         fr(2:n, :) = fl(1:n-1, :) + w(2:n, :)
         fr(1, :) = fr(2, :)

         select case (istep)
         case (RK2_2)
            ! Second order flux corrections

            dfp(1:n-1, :) = half * (fr(2:n, :) - fr(1:n-1, :)); dfp(n, :) = dfp(n-1, :)
            dfm(2:n, :)   = dfp(1:n-1, :);                      dfm(1, :) = dfm(2, :)
            ! Flux limiter application
            call flimiter(fr, dfm, dfp)

            dfp(1:n-1, :) = half * (fl(1:n-1, :) - fl(2:n, :)); dfp(n, :) = dfp(n-1, :)
            dfm(2:n, :)   = dfp(1:n-1, :);                      dfm(1, :) = dfm(2, :)
            call flimiter(fl, dfm, dfp)
            !OPT 60% of D1mr and 40% D1mw occurred in few above lines (D1mr = 0.1% Dr, D1mw = 0.5% Dw)
            ! That ^^ should be fixed now, please confirm
         case (RK2_1)
         case default
            call die("[rtvd:relaxing_tvd] Unsupported substep")
         end select

         ! u update
         fu = fr - fl

         ! Setting uflx to 0 can be used to enforce reflecting boundary
         ! To enforce diode boundaries one can either:
         ! * Set fl or fr to 0 (depending on side)
         ! or
         ! * Set f0 to 0 only when it would produce incoming flux.
         ! I don't remember which approach was already (unsuccesfully) tested
         ! \todo remove transpositions by changing index order in eflx
         if (associated(eflx%li)) fu(eflx%li%index, :) = eflx%li%uflx
         if (associated(eflx%ri)) fu(eflx%ri%index, :) = eflx%ri%uflx
         if (associated(eflx%lo)) eflx%lo%uflx = fu(eflx%lo%index, :)
         if (associated(eflx%ro)) eflx%ro%uflx = fu(eflx%ro%index, :)

         if (dom%geometry_type == GEO_XYZ) then
            u1(2:n, :) = u0(2:n, :) -                  dtx * (                fu(2:n, :) -                fu(1:n-1, :) )
         else
            u1(2:n, :) = u0(2:n, :) - gc(GC1,2:n, :) * dtx * ( gc(GC2,2:n, :)*fu(2:n, :) - gc(GC3,2:n, :)*fu(1:n-1, :) )
         endif
         u1(1, :)   = u1(2, :)
      endif ! (n > 1)

   end subroutine relaxing_tvd

end module rtvd
