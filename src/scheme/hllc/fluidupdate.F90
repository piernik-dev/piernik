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
module fluidupdate   ! SPLIT MUSCL HANCOCK
! pulled by HLLC
!BEWARE: this module only care about neutral fluid

  implicit none
  private
  public :: fluid_update

contains
   subroutine fluid_update
      use dataio_pub,     only: halfstep
      use global,         only: dt, dtm, t
      use user_hooks,     only: problem_customize_solution
      use grid,           only: cga
      use grid_cont,      only: cg_list_element
      use constants,      only: xdim, zdim
      use domain,         only: has_dir

      implicit none

      logical, save                  :: first_run = .true.
      type(cg_list_element), pointer :: cgl
      integer                        :: ddim

      halfstep = .false.
      if (first_run) then
         dtm = 0.0
      else
         dtm = dt
      endif

      t=t+dt

      call cga%get_root(cgl)
      do while (associated(cgl))
         do ddim = xdim, zdim, 1
            if (has_dir(ddim)) call sweep(cgl%cg,dt,ddim)
         enddo
         if (associated(problem_customize_solution)) call problem_customize_solution
         cgl => cgl%nxt
      enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      t=t+dt
      dtm = dt
      halfstep = .true.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call cga%get_root(cgl)
      do while (associated(cgl))
         do ddim = zdim, xdim, -1
            if (has_dir(ddim)) call sweep(cgl%cg,dt,ddim)
         enddo
         if (associated(problem_customize_solution)) call problem_customize_solution
         cgl => cgl%nxt
      enddo

      if (first_run) first_run = .false.

   end subroutine fluid_update
!---------------------------------------------------------------------------
   subroutine sweep(cg,dt,ddim)
      use fluidboundaries, only: all_fluid_boundaries
      use fluidindex,      only: iarr_all_swp, ibx, ibz
      use grid_cont,       only: grid_container
      use constants,       only: ndims, LO, HI, pdims

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      real,    intent(in)                       :: dt
      integer, intent(in)                       :: ddim

      integer :: i1, i2
      integer, dimension(ndims,LO:HI)       :: i

      real, dimension(size(cg%u%arr,1),cg%n_(ddim)) :: u1d
      real, dimension(:,:), pointer                 :: pu
      real, dimension(ibx:ibz,         cg%n_(ddim)) :: b1d

      do i2 = 1, cg%n_(pdims(ddim,3))
         do i1 = 1, cg%n_(pdims(ddim,2))
            pu => cg%u%get_sweep(ddim,i1,i2)
            u1d(iarr_all_swp(ddim,:),:) = pu(:,:)
            call sweep1d_mh(u1d,b1d,cg%cs_iso2%get_sweep(ddim,i1,i2),dt/cg%dl(ddim))
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
   function calculate_slope_moncen(u) result(dq)
      use constants,     only: half, one
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
!---------------------------------------------------------------------------
   subroutine sweep1d_mh(u,b,cs2,dtodx)
      use constants,    only: half
      use fluidtypes,   only: component_fluid
      use fluidindex,   only: flind
      implicit none
      real, intent(in)                           :: dtodx
      real, dimension(:), intent(in), pointer    :: cs2
      real, dimension(:,:), intent(in)           :: b
      real, dimension(:,:), intent(inout)        :: u

      type(component_fluid), pointer             :: fl

      real, dimension(size(u,1),size(u,2)), target :: flux, ql, qr, qgdn
      real, dimension(size(u,1),size(u,2)), target :: du, ul, ur, u_l, u_r
      real, dimension(:,:), pointer :: p_ql, p_qr, p_q, p_flux
      integer :: nx, p

      nx = size(u,2)

      du = calculate_slope_vanleer(u)
      ul = u - half*du   ! (14.33)
      ur = u + half*du

      flux = compute_flux(ul,b,cs2) - compute_flux(ur,b,cs2)    ! interpolate b?

      u_l = ur + half*dtodx*flux   ! (14.34) + (14.35)
      u_r(:,1:nx-1) = ul(:,2:nx) + half*dtodx*flux(:,2:nx); u_r(:,nx) = u_r(:,nx-1)

      ql = utoq(u_l,b)
      qr = utoq(u_r,b)

      do p = 1, flind%fluids
         fl     => flind%all_fluids(p)
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
      use constants,    only: half, two
      use fluidtypes,   only: component_fluid
      use fluidindex,   only: flind

      implicit none

      real, dimension(:,:), intent(in)           :: u, b
      real, dimension(size(u,1),size(u,2))       :: q

      integer :: p
      type(component_fluid), pointer :: fl

      do p = 1, flind%fluids
         fl => flind%all_fluids(p)

         q(fl%idn,:) = u(fl%idn,:)
         q(fl%imx,:) = u(fl%imx,:)/u(fl%idn,:)
         q(fl%imy,:) = u(fl%imy,:)/u(fl%idn,:)
         q(fl%imz,:) = u(fl%imz,:)/u(fl%idn,:)

         if (fl%has_energy) then
            q(fl%ien,:) = (u(fl%ien,:) - half*( u(fl%imx,:)**2 + u(fl%imy,:)**2 + u(fl%imz,:)**2) / u(fl%idn,:)) * fl%gam_1
            if (fl%is_magnetized) q(fl%ien,:) = q(fl%ien,:) + (two-fl%gam)*half*sum(b(:,:)**2,dim=1)
         endif
      enddo
   end function utoq
!---------------------------------------------------------------------------
   function compute_flux(u,b,cs2) result(f)
      use constants,    only: half, two
      use fluidindex,   only: flind, ibx, iby, ibz
      use fluidtypes,   only: component_fluid
      implicit none
      real, dimension(:,:), intent(in)       :: u, b
      real, dimension(:),   intent(in)       :: cs2
      real, dimension(size(u,1), size(u,2))  :: f
      real, dimension(size(u,2))             :: vx, p

      integer :: ip
      type(component_fluid), pointer :: fl

      do ip = 1, flind%fluids
         fl => flind%all_fluids(ip)

         vx = u(fl%imx,:) / u(fl%idn,:)
         if (fl%has_energy) then
            p = (u(fl%ien,:) - half*( u(fl%imx,:)**2 + u(fl%imy,:)**2 + u(fl%imz,:)**2) / u(fl%idn,:)) * flind%neu%gam_1
            if (fl%is_magnetized) p = p + (two-fl%gam)*half*sum(b**2,dim=1)
         else
            p = cs2*u(fl%idn,:)
         endif

         f(fl%idn,:) = u(fl%imx,:)
         f(fl%imx,:) = u(fl%imx,:)*vx + p
         f(fl%imy,:) = u(fl%imy,:)*vx
         f(fl%imz,:) = u(fl%imz,:)*vx
         if (fl%is_magnetized) then
            f(fl%imx,:) = f(fl%imx,:) - b(ibx,:)*b(ibx,:)
            f(fl%imy,:) = f(fl%imy,:) - b(iby,:)*b(ibx,:)
            f(fl%imz,:) = f(fl%imz,:) - b(ibz,:)*b(ibx,:)
         endif

         if (fl%has_energy) then
            f(fl%ien,:) = (u(fl%ien,:) + p(:))*vx(:)
            if (fl%is_magnetized) f(fl%ien,:) = f(fl%ien,:) - b(ibx,:)*(b(ibx,:)*u(fl%imx,:)+b(iby,:)*u(fl%imy,:)+b(ibz,:)*u(fl%imz,:))/u(fl%idn,:)
         endif
      enddo

      return

   end function compute_flux
!---------------------------------------------------------------------------
   subroutine riemann_hllc(qleft,qright,qgdnv,fgdnv, n, gamma, cs2)

      use constants,  only: zero, one, half
      use global,     only: smalld
      use fluidindex, only: idn, imx, imy, imz, ien

      implicit none

      real, parameter    :: smallc = 1.e-8

      ! HLLC Riemann solver (Toro)
      integer, intent(in) :: n
      real, intent(in)    :: gamma
      real, dimension(:,:), intent(in), pointer     :: qleft,qright
      real, dimension(:,:), intent(inout), pointer  :: qgdnv,fgdnv
      real, dimension(:), intent(in)                :: cs2

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

#ifndef ISO
      real, dimension(n) :: ekinl, ekinr
#endif
#ifdef __IFORT__
      integer :: i
#endif

      ! constants
      smallp = 1.e-7   ! BEWARE

#ifndef ISO
      entho = one/(gamma-one)
#else
      entho = 1.e25
#endif
      ! Left variables
      rl = max(qleft (idn,:), smalld)
      ul =     qleft (imx,:)
#ifndef ISO
      Pl = max(qleft (ien,:),rl(:)*smallp)

      ekinl = half * rl * ( ul*ul + qleft(imy,:)**2 + qleft(imz,:)**2 )
      etotl = Pl*entho + ekinl
#else
      Pl = cs2 * rl
      etotl = zero
#endif
      Ptotl = Pl

      ! Right variables
      rr = max(qright(idn,:), smalld)
      ur =     qright(imx,:)
#ifndef ISO
      Pr = max(qright(ien,:), rr*smallp)

      ekinr = half * rr * ( ur*ur + qright(imy,:)**2 + qright(imz,:)**2 )
      etotr = Pr*entho + ekinr
#else
      PR = cs2 * rr
      etotr = zero
#endif
      Ptotr = Pr

      ! Compute average velocity
      qgdnv(imx,:) = half*( qleft(imx,:) + qright(imx,:) )

      ! Find the largest eigenvalues in the normal direction to the interface
      cfastl=sqrt(max(gamma*Pl/rl,smallc**2))
      cfastr=sqrt(max(gamma*Pr/rr,smallc**2))

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
#else
      etotstarl=zero
#endif
      ! Right star region variables
      rstarr=rr*(SR-ur)/(SR-ustar)
#ifndef ISO
      etotstarr=((SR-ur)*etotr-Ptotr*ur+Ptotstar*ustar)/(SR-ustar)
#else
      etotstarr=zero
#endif

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
      fgdnv(ien,:) = (etoto+Ptoto)*uo
#endif /* ISO */

      do ivar = imy,imz

      ! BEWARE the version with WHERE had huge, unexplained memory leaks when compiled with some Intel compilers
      ! The __IFORT__ macro has to be defined manually, e.g. in appropriate compiler.in file
#ifndef __IFORT__
          where (fgdnv(idn,:)>zero)
             fgdnv(ivar,:) = fgdnv(idn,:)*qleft (ivar,:)
          elsewhere
             fgdnv(ivar,:) = fgdnv(idn,:)*qright(ivar,:)
          endwhere
#else
          do i = lbound(fgdnv(:,:),2), ubound(fgdnv(:,:),2)
             if (fgdnv(idn,i)>zero) then
                fgdnv(ivar,i) = fgdnv(idn,i)*qleft (ivar,i)
             else
                fgdnv(ivar,i) = fgdnv(idn,i)*qright(ivar,i)
             endif
          enddo
#endif
      enddo
      return
#ifndef ISO
      if (.false.) write(0,*) cs2
#endif /* ISO */
   end subroutine riemann_hllc
!---------------------------------------------------------------------------
end module fluidupdate
