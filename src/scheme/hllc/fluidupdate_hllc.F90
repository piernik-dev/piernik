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
!! \brief module for HLLC scheme (SPLIT MUSCL HANCOCK)
!!
!! \deprecated BEWARE: this module only care about neutral fluid
!<
module fluidupdate_hllc   ! SPLIT MUSCL HANCOCK
! pulled by ANY

   implicit none
   private
   public :: fluid_update_simple

contains

   subroutine fluid_update_simple

      use cg_list,     only: cg_list_element
      use cg_leaves,   only: leaves
      use constants,   only: xdim, zdim
      use dataio_pub,  only: halfstep, die
      use domain,      only: dom, is_multicg
      use global,      only: dt, dtm, t
      use user_hooks,  only: problem_customize_solution

      implicit none

      logical, save                  :: first_run = .true.
      type(cg_list_element), pointer :: cgl
      integer(kind=4)                :: ddim

      !> \todo figure out what the problem is and enable multicg and AMR as well
      if (is_multicg) call die("[fluidupdate_hllc:fluid_update_simple] something here is not compatible with multiple blocks per process yet")
#ifdef MAGNETIC
      call die("[fluidupdate_hllc:fluid_update_simple] Magnetic field is not compatible with HLLC")
#endif /* MAGNETIC */

      halfstep = .false.
      if (first_run) then
         dtm = 0.0
      else
         dtm = dt
      endif

      t=t+dt

      cgl => leaves%first
      do while (associated(cgl))
         do ddim = xdim, zdim, 1
            if (dom%has_dir(ddim)) call sweep(cgl%cg,dt,ddim)
         enddo
         if (associated(problem_customize_solution)) call problem_customize_solution(.true.)
         cgl => cgl%nxt
      enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      t=t+dt
      dtm = dt
      halfstep = .true.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      cgl => leaves%first
      do while (associated(cgl))
         do ddim = zdim, xdim, -1
            if (dom%has_dir(ddim)) call sweep(cgl%cg,dt,ddim)
         enddo
         if (associated(problem_customize_solution)) call problem_customize_solution(.false.)
         cgl => cgl%nxt
      enddo

      if (first_run) first_run = .false.

   end subroutine fluid_update_simple
!---------------------------------------------------------------------------
   subroutine sweep(cg,dt,ddim)

      use constants,        only: pdims, xdim, zdim, cs_i2_n, ORTHO1, ORTHO2, LO, HI
      use all_boundaries,   only: all_fluid_boundaries
      use dataio_pub,       only: warn
      use fluidindex,       only: iarr_all_swp
      use grid_cont,        only: grid_container
      use named_array_list, only: qna, wna

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      real,                          intent(in) :: dt
      integer(kind=4),               intent(in) :: ddim

      integer                                   :: i1, i2, i_cs_iso2
      real, dimension(size(cg%u,1),cg%n_(ddim)) :: u1d
      real, dimension(xdim:zdim,   cg%n_(ddim)) :: b1d
      real, dimension(:,:), pointer             :: pu
      real, dimension(:),   pointer             :: cs2
      logical, save                             :: firstcall = .true.

      cs2 => null()
      if (qna%exists(cs_i2_n)) then
         i_cs_iso2 = qna%ind(cs_i2_n) !> \deprecated BEWARE: magic strings across multiple files
      else
         i_cs_iso2 = -1
      endif

      b1d=0.
      if (firstcall) call warn("[fluidupdate_hllc:sweep] magnetic field unimplemented yet. Forcing to be 0")
      firstcall = .false.

      do i2 = cg%lhn(pdims(ddim, ORTHO2), LO), cg%lhn(pdims(ddim, ORTHO2), HI)
         do i1 = cg%lhn(pdims(ddim, ORTHO1), LO), cg%lhn(pdims(ddim, ORTHO1), HI)
            pu => cg%w(wna%fi)%get_sweep(ddim, i1, i2)
            if (i_cs_iso2 > 0) cs2 => cg%q(i_cs_iso2)%get_sweep(ddim, i1, i2)

            u1d(iarr_all_swp(ddim,:),:) = pu(:,:)

            call sweep1d_mh(u1d,b1d,cs2,dt/cg%dl(ddim))

            pu(:,:) = u1d(iarr_all_swp(ddim,:),:)
         enddo
      enddo
      call all_fluid_boundaries
   end subroutine sweep
!---------------------------------------------------------------------------
   function calculate_slope_vanleer(u) result(dq)

      implicit none

      real, dimension(:,:), intent(in)     :: u

      real, dimension(size(u,1),size(u,2)) :: dlft, drgt, dcen, dq
      integer :: n

      n = size(u,2)

      dlft(:,2:n)   = (u(:,2:n) - u(:,1:n-1)) ; dlft(:,1) = dlft(:,2)    ! (14.38)
      drgt(:,1:n-1) = dlft(:,2:n) ;             drgt(:,n) = drgt(:,n-1)

      dcen = dlft*drgt

      where (dcen>0.0)
         dq = 2.0*dcen / (dlft+drgt)       ! (14.54) ?
      elsewhere
         dq = 0.0
      endwhere

   end function calculate_slope_vanleer
!---------------------------------------------------------------------------
#if 0

! An unused function.
! It is an alternative to calculate_slope_vanleer.
! It may be worth making an option if we ever bring back the HLLC solver to regular use.

   function calculate_slope_moncen(u) result(dq)

      use constants, only: half, one

      implicit none

      real, dimension(:,:), intent(in)     :: u

      real, dimension(size(u,1),size(u,2)) :: dlft,drgt,dcen,dlim, dq
      integer :: n
      real :: sl

      sl = one

      n = size(u,2)

      dlft(:,2:n)   = sl*(u(:,2:n)   - u(:,1:n-1)) ; dlft(:,1) = dlft(:,2)
      drgt(:,1:n-1) = dlft(:,2:n) ;                   drgt(:,n) = drgt(:,n-1)

      dcen = half*(dlft+drgt)/sl

      where (dlft*drgt<=0.0)
         dlim = 0.0
      elsewhere
         dlim = min(abs(dlft),abs(drgt))
      endwhere
      dq = sign(1.0, dcen) * min(dlim,abs(dcen))

   end function calculate_slope_moncen
#endif /* 0 */
!---------------------------------------------------------------------------
   subroutine sweep1d_mh(u,b,cs2,dtodx)

      use constants,  only: half
      use fluidindex, only: flind
      use fluidtypes, only: component_fluid

      implicit none

      real,                        intent(in)    :: dtodx
      real, dimension(:), pointer, intent(in)    :: cs2
      real, dimension(:,:),        intent(in)    :: b
      real, dimension(:,:),        intent(inout) :: u

      class(component_fluid), pointer              :: fl
      real, dimension(size(u,1),size(u,2)), target :: flux, ql, qr, qgdn
      real, dimension(size(u,1),size(u,2)), target :: du, ul, ur, u_l, u_r
      real, dimension(:,:), pointer                :: p_ql, p_qr, p_q, p_flux
      integer :: nx, p

      nx = size(u,2)

      du = calculate_slope_vanleer(u)
      ul = u - half*du   ! (14.33)
      ur = u + half*du

      flux = compute_flux(ul,b,cs2) - compute_flux(ur,b,cs2)    !> \todo interpolate b?

      u_l = ur + half*dtodx*flux   ! (14.34) + (14.35)
      u_r(:,1:nx-1) = ul(:,2:nx) + half*dtodx*flux(:,2:nx); u_r(:,nx) = u_r(:,nx-1)

      ql = utoq(u_l,b)
      qr = utoq(u_r,b)

      do p = 1, flind%fluids
         fl     => flind%all_fluids(p)%fl
         p_ql   => ql(fl%beg:fl%end,:)
         p_qr   => qr(fl%beg:fl%end,:)
         p_q    => qgdn(fl%beg:fl%end,:)
         p_flux => flux(fl%beg:fl%end,:)
         call riemann_hllc(p_ql, p_qr, p_q, p_flux, nx, fl%gam, cs2)
      enddo

      u(:,2:nx) = u(:,2:nx) + dtodx*(flux(:,1:nx-1) - flux(:,2:nx))
      u(:,1)  = u(:,2); u(:,nx) = u(:,nx-1)

   end subroutine sweep1d_mh
!---------------------------------------------------------------------------
   function utoq(u,b) result(q)

      use constants,  only: half, two
      use fluidindex, only: flind
      use fluidtypes, only: component_fluid
      use func,       only: ekin

      implicit none

      real, dimension(:,:), intent(in)           :: u, b

      real, dimension(size(u,1),size(u,2))       :: q
      integer :: p
      class(component_fluid), pointer :: fl

      do p = 1, flind%fluids
         fl => flind%all_fluids(p)%fl

         q(fl%idn,:) = u(fl%idn,:)
         q(fl%imx,:) = u(fl%imx,:)/u(fl%idn,:)
         q(fl%imy,:) = u(fl%imy,:)/u(fl%idn,:)
         q(fl%imz,:) = u(fl%imz,:)/u(fl%idn,:)

         if (fl%has_energy) then
            q(fl%ien,:) = (u(fl%ien,:) - ekin(u(fl%imx,:), u(fl%imy,:), u(fl%imz,:), u(fl%idn,:))) * fl%gam_1
            if (fl%is_magnetized) q(fl%ien,:) = q(fl%ien,:) + (two-fl%gam)*half*sum(b(:,:)**2,dim=1)
         endif
      enddo
   end function utoq
!---------------------------------------------------------------------------
   function compute_flux(u,b,cs2) result(f)

      use constants,  only: half, two, xdim, ydim, zdim
      use fluidindex, only: flind
      use fluidtypes, only: component_fluid
      use func,       only: ekin

      implicit none

      real, dimension(:,:),        intent(in) :: u, b
      real, dimension(:), pointer, intent(in) :: cs2

      real, dimension(size(u,1), size(u,2))  :: f
      real, dimension(size(u,2))             :: vx, p
      integer :: ip
      class(component_fluid), pointer :: fl

      do ip = 1, flind%fluids
         fl => flind%all_fluids(ip)%fl

         vx = u(fl%imx,:) / u(fl%idn,:)
         if (fl%has_energy) then
            p = (u(fl%ien,:) - ekin(u(fl%imx,:), u(fl%imy,:), u(fl%imz,:), u(fl%idn,:))) * fl%gam_1
            if (fl%is_magnetized) p = p + (two-fl%gam)*half*sum(b**2,dim=1)
         else
            if (associated(cs2)) then
               p = cs2*u(fl%idn,:)
            else
               p = 0.
            endif

         endif

         f(fl%idn,:) = u(fl%imx,:)
         f(fl%imx,:) = u(fl%imx,:)*vx + p
         f(fl%imy,:) = u(fl%imy,:)*vx
         f(fl%imz,:) = u(fl%imz,:)*vx
         if (fl%is_magnetized) then
            f(fl%imx,:) = f(fl%imx,:) - b(xdim,:)*b(xdim,:)
            f(fl%imy,:) = f(fl%imy,:) - b(ydim,:)*b(xdim,:)
            f(fl%imz,:) = f(fl%imz,:) - b(zdim,:)*b(xdim,:)
         endif

         if (fl%has_energy) then
            f(fl%ien,:) = (u(fl%ien,:) + p(:))*vx(:)
            if (fl%is_magnetized) f(fl%ien,:) = f(fl%ien,:) - b(xdim,:)*(b(xdim,:)*u(fl%imx,:)+b(ydim,:)*u(fl%imy,:)+b(zdim,:)*u(fl%imz,:))/u(fl%idn,:)
         endif
      enddo

      return

   end function compute_flux
!---------------------------------------------------------------------------
!>
!! \brief HLLC Riemann solver (Toro)
!<
   subroutine riemann_hllc(qleft,qright,qgdnv,fgdnv, n, gamma, cs2)

      use constants,  only: zero, half, idn, imx, imy, imz, ien
      use global,     only: smalld
#ifndef ISO
      use constants,  only: one
#endif /* !ISO */

      implicit none

      integer,                       intent(in)    :: n
      real,                          intent(in)    :: gamma
      real, dimension(:,:), pointer, intent(in)    :: qleft,qright
      real, dimension(:,:), pointer, intent(inout) :: qgdnv,fgdnv
      real, dimension(:),            intent(in)    :: cs2

      real, parameter    :: smallc = 1.e-8
      real, dimension(n) :: SL,SR
      real, dimension(n) :: rl,Pl,ul,etotl,ptotl
      real, dimension(n) :: rr,Pr,ur,etotr,ptotr
      real, dimension(n) :: cfastl,rcl,rstarl
      real, dimension(n) :: cfastr,rcr,rstarr
      real, dimension(n) :: etotstarl,etotstarr
      real, dimension(n) :: ustar,Ptotstar
      real, dimension(n) :: ro,uo,Ptoto,etoto
      real    :: smallp, entho
      integer :: ivar
      logical :: has_e ! has_energy, .false. for dust
#ifndef ISO
      real, dimension(n) :: ekinl, ekinr
#endif /* !ISO */
#ifdef __INTEL_COMPILER
      integer :: i
#endif /* __INTEL_COMPILER */

      ! constants
      smallp = 1.e-7   !> \deprecated BEWARE

      has_e = (size(qleft, dim=1) >= ien)
#ifndef ISO
      entho = one/(gamma-one)
#else /* ISO */
      entho = 1.e25
#endif /* ISO */
      ! Left variables
      rl = max(qleft (idn,:), smalld)
      ul =     qleft (imx,:)
#ifndef ISO
      if (has_e) then
         Pl = max(qleft (ien,:),rl(:)*smallp)
      else
         Pl = rl(:)*smallp
      endif

      ekinl = half * rl * ( ul*ul + qleft(imy,:)**2 + qleft(imz,:)**2 )
      etotl = Pl*entho + ekinl
#else /* ISO */
      Pl = cs2 * rl
      etotl = zero
#endif /* ISO */
      Ptotl = Pl

      ! Right variables
      rr = max(qright(idn,:), smalld)
      ur =     qright(imx,:)
#ifndef ISO
      if (has_e) then
         Pr = max(qright(ien,:), rr*smallp)
      else
         Pr = rr*smallp
      endif

      ekinr = half * rr * ( ur*ur + qright(imy,:)**2 + qright(imz,:)**2 )
      etotr = Pr*entho + ekinr
#else /* ISO */
      PR = cs2 * rr
      etotr = zero
#endif /* ISO */
      Ptotr = Pr

      ! Compute average velocity
      qgdnv(imx,:) = half*( qleft(imx,:) + qright(imx,:) )

      ! Find the largest eigenvalues in the normal direction to the interface

#if 0
      ! strange FPExceptions here with gfortran 4.8.3 20140911 (Red Hat 4.8.3-7) when gamma*Pl/rl = -1.e-9 and there is -O3 (switching to -O2 fixes it)
      cfastl=sqrt(max(gamma*Pl/rl,smallc**2))
      cfastr=sqrt(max(gamma*Pr/rr,smallc**2))
#else /* !0 */
      cfastl=max(gamma*Pl/rl,smallc**2)
      cfastl=sqrt(cfastl)
      cfastr=max(gamma*Pr/rr,smallc**2)
      cfastr=sqrt(cfastr)
#endif /* !0 */

      ! Compute HLL wave speed
      SL=min(ul,ur)-max(cfastl,cfastr)
      SR=max(ul,ur)+max(cfastl,cfastr)

      ! Compute lagrangian sound speed
      rcl=rl*(ul-SL)
      rcr=rr*(SR-ur)

      ! Compute acoustic star state
      ustar   =(rcr*ur   +rcl*ul   +  (Ptotl-Ptotr))/(rcr+rcl)
      Ptotstar=(rcr*Ptotl+rcl*Ptotr+rcl*rcr*(ul-ur))/(rcr+rcl)

      ! Left star region variables
      rstarl=rl*(SL-ul)/(SL-ustar)
#ifndef ISO
      etotstarl=((SL-ul)*etotl-Ptotl*ul+Ptotstar*ustar)/(SL-ustar)
#else /* ISO */
      etotstarl=zero
#endif /* ISO */
      ! Right star region variables
      rstarr=rr*(SR-ur)/(SR-ustar)
#ifndef ISO
      etotstarr=((SR-ur)*etotr-Ptotr*ur+Ptotstar*ustar)/(SR-ustar)
#else /* ISO */
      etotstarr=zero
#endif /* ISO */

      ! Sample the solution at x/t=0
      where (SL>zero)
         ro=rl
         uo=ul
         Ptoto=Ptotl
         etoto=etotl
      elsewhere (ustar>zero)
         ro=rstarl
         uo=ustar
         Ptoto=Ptotstar
         etoto=etotstarl
      elsewhere (SR>zero)
         ro=rstarr
         uo=ustar
         Ptoto=Ptotstar
         etoto=etotstarr
      elsewhere
         ro=rr
         uo=ur
         Ptoto=Ptotr
         etoto=etotr
      endwhere

      ! Compute the Godunov flux
      fgdnv(idn,:) = ro*uo
      fgdnv(imx,:) = ro*uo*uo+Ptoto
#ifndef ISO
      if (has_e) fgdnv(ien,:) = (etoto+Ptoto)*uo
#endif /* !ISO */

      do ivar = imy,imz

      ! BEWARE the version with WHERE had huge, unexplained memory leaks when compiled with some Intel compilers
#ifndef __INTEL_COMPILER
         where (fgdnv(idn,:)>zero)
            fgdnv(ivar,:) = fgdnv(idn,:)*qleft (ivar,:)
         elsewhere
            fgdnv(ivar,:) = fgdnv(idn,:)*qright(ivar,:)
         endwhere
#else /* !__INTEL_COMPILER */
         do i = lbound(fgdnv(:,:),2), ubound(fgdnv(:,:),2)
            if (fgdnv(idn,i)>zero) then
               fgdnv(ivar,i) = fgdnv(idn,i)*qleft (ivar,i)
            else
               fgdnv(ivar,i) = fgdnv(idn,i)*qright(ivar,i)
            endif
         enddo
#endif /* !__INTEL_COMPILER */
      enddo
      return
#ifndef ISO
      if (.false.) write(0,*) cs2
#endif /* !ISO */
   end subroutine riemann_hllc
!---------------------------------------------------------------------------
end module fluidupdate_hllc
