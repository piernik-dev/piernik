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
!! \brief COMMENT ME
!!
!!    HLLD Riemann solver for ideal magnetohydrodynamics
!!    Varadarajan Parthasarathy, CAMK, Warszawa. 2015.
!!    Dr. Artur Gawryszczak, CAMK, Warszawa.
!!
!!    Energy fix up routines for CT and its related comments are not used in the current version.
!!    The algorithm is simply present for experimental purposes.
!<

module solvecg

! pulled by RIEMANN

   implicit none

   private
   public  :: solve_cg

contains

   subroutine solve_cg(cg, ddim, istep, fargo_vel)

      use constants,        only: pdims, xdim, zdim, ORTHO1, ORTHO2, LO, HI, psi_n, INVALID, GEO_XYZ, I_ZERO, I_ONE, first_stage
      use dataio_pub,       only: die
      use domain,           only: dom, is_refined
      use fluidindex,       only: flind, iarr_all_dn, iarr_all_mx, iarr_all_swp, iarr_mag_swp
      use global,           only: dt, force_cc_mag, use_fargo, integration_order
      use grid_cont,        only: grid_container
      use named_array_list, only: wna, qna
      use sources,          only: prepare_sources, all_sources
#ifdef COSM_RAYS
      use sources,          only: limit_minimal_ecr
#endif /* COSM_RAYS */

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer(kind=4),               intent(in) :: ddim
      integer,                       intent(in) :: istep     ! stage in the time integration scheme
      integer(kind=4), optional,     intent(in) :: fargo_vel

      real, dimension(cg%n_(ddim), size(cg%u,1)) :: u1d
      real, dimension(cg%n_(ddim), xdim:zdim)    :: b_cc1d
      real, dimension(cg%n_(ddim), 1)            :: psi_d ! artificial rank-2 to conform to flux limiter interface
      real, dimension(:,:), pointer              :: pu, pb
      integer                                    :: i1, i2
      integer                                    :: bi
      real, dimension(:), pointer                :: ppsi
      integer                                    :: psii
      real, dimension(size(u1d,1),size(u1d,2))           :: u0
      real, dimension(size(u1d,1), flind%fluids), target :: vx
      integer(kind=4)                            :: nmag, i

      ! is_multicg should be safe
      if (.false.) write(0,*) present(fargo_vel) ! suppress compiler warning on unused argument
      if (istep /= first_stage(integration_order)) call die("[solve_cg:solve_cg] istep >1 not implemented yet")
      if (use_fargo) call die("[solve_cg:solve_cg] Fargo is not yet enabled for Riemann")
      if (is_refined) call die("[solve_cg:solve_cg] This Rieman solver is not compatible with mesh refinements yet!")
      if (dom%geometry_type /= GEO_XYZ) call die("[solve_cg:solve_cg] Non-cartesian geometry is not implemented yet in this Riemann solver.")
      nmag = I_ZERO
      do i = 1, flind%fluids
         if (flind%all_fluids(i)%fl%is_magnetized) nmag = nmag + I_ONE
      enddo
      if (nmag > 1) call die("[solve_cg:solve_cg] At most one magnetized fluid is implemented")

      if (force_cc_mag) then
         bi = wna%bi
      else
         bi = wna%bcci
      endif

      psii = INVALID
      if (qna%exists(psi_n)) psii = qna%ind(psi_n)
      psi_d = 0.
      nullify(ppsi)

      call prepare_sources(cg)

      do i2 = cg%ijkse(pdims(ddim, ORTHO2), LO), cg%ijkse(pdims(ddim, ORTHO2), HI)
         do i1 = cg%ijkse(pdims(ddim, ORTHO1), LO), cg%ijkse(pdims(ddim, ORTHO1), HI)
            pu => cg%w(wna%fi)%get_sweep(ddim,i1,i2)
            u1d(:, iarr_all_swp(ddim,:)) = transpose(pu(:,:))
            pb => cg%w(bi)%get_sweep(ddim,i1,i2)
            b_cc1d(:, iarr_mag_swp(ddim,:)) = transpose(pb(:,:))

            if (psii /= INVALID) then
               ppsi => cg%q(psii)%get_sweep(ddim,i1,i2)
               psi_d(:, 1) = ppsi(:)
               call solve(u1d, b_cc1d, dt/cg%dl(ddim), psi_d)
               ppsi(:) = psi_d(:,1)
            else
               call solve(u1d, b_cc1d, dt/cg%dl(ddim), psi_d)
            endif

            ! This is lowest order implementation of CR
            ! It agrees with the implementation in RTVD in the limit of small CR energy amounts
            ! ToDo: integrate it into h_solver schemes

            ! transposition for compatibility with RTVD-based routines
            u0 = transpose(pu)
            vx = u1d(:, iarr_all_mx) / u1d(:, iarr_all_dn) ! this may also be useful for gravitational acceleration
            call all_sources(size(u1d, 1, kind=4), u0, u1d, b_cc1d, cg, istep, ddim, i1, i2, dt, vx)

            ! Beware: this is bypassing integration scheme, so the source terms are applied in lowest order fashion.
            ! See the results of Jeans test with RTVD and RIEMANN for comparision.

#if defined COSM_RAYS && defined IONIZED
            if (size(u1d, 1) > 1) call limit_minimal_ecr(size(u1d, 1), u1d)
#endif /* COSM_RAYS && IONIZED */

            pu(:,:) = transpose(u1d(:,iarr_all_swp(ddim,:)))
            pb(:,:) = transpose(b_cc1d(:, iarr_mag_swp(ddim,:))) ! ToDo figure out how to manage CT energy fixup without extra storage
         enddo
      enddo

   end subroutine solve_cg

!! k-th interface is between k-th cell and (k+1)-th cell
!! We don't calculate n-th interface because it is as incomplete as 0-th interface

#undef RK_HIGH
!! Disabled rk3 as high orders require more careful approach to guardcells and valid ranges

   subroutine solve(u, b_cc, dtodx, psi)

      use constants,  only: half
#ifdef RK_HIGH
      use constants,  only: oneq, twot
#endif /* RK_HIGH */
      use dataio_pub, only: die
      use global,     only: h_solver
      use hlld,           only: riemann_wrap
      use interpolations, only: interpol

      implicit none

      real, dimension(:,:), intent(inout) :: u
      real, dimension(:,:), intent(inout) :: b_cc
      real,                 intent(in)    :: dtodx
      real, dimension(:,:), intent(inout) :: psi

      ! left and right states at interfaces 1 .. n-1
      real, dimension(size(u,   1)-1, size(u,   2)), target :: ql, qr
      real, dimension(size(b_cc,1)-1, size(b_cc,2)), target :: b_cc_l, b_cc_r
      real, dimension(size(psi, 1)-1, size(psi, 2)), target :: psi_l, psi_r

      ! fluxes through interfaces 1 .. n-1
      real, dimension(size(u,   1)-1, size(u,   2)), target :: flx
      real, dimension(size(b_cc,1)-1, size(b_cc,2)), target :: mag_cc
      real, dimension(size(psi, 1)-1, size(psi, 2)), target :: psi_cc

      ! left and right MUSCL fluxes at interfaces 1 .. n-1
      real, dimension(size(u,   1)-1,size(u,   2))          :: flx_l, flx_r
      real, dimension(size(b_cc,1)-1,size(b_cc,2))          :: bclflx, bcrflx
      real, dimension(size(psi, 1)-1,size(psi, 2))          :: psilflx, psirflx

      ! left and right states for cells 2 .. n-1
      real, dimension(2:size(u,   1)-1, size(u,   2))       :: u1    !, du2, du3
      real, dimension(2:size(b_cc,1)-1, size(b_cc,2))       :: b1    !, db2, db3
      real, dimension(2:size(psi, 1)-1, size(psi, 2))       :: psi1  !, dpsi2, dpsi3

#ifdef RK_HIGH
      ! left and right states for cells (RK3)
      ! " to be checked "
      real, dimension(2:size(u,   1)-2, size(u,  2 ))        :: u2
      real, dimension(2:size(b_cc,1)-2, size(b_cc,2))        :: b2
      real, dimension(2:size(psi, 1)-2, size(psi, 2))        :: psi2
#endif /* RK_HIGH */

      ! updates required for higher order of integration will likely have shorter length

!!$     ! fluxes through interfaces 2 .. n-2
!!$     real, dimension(size(u,   1), 2:size(u,   2)-2), target :: flx1
!!$     real, dimension(size(b_cc,1), 2:size(b_cc,2)-2), target :: mag_cc1
!!$     real, dimension(size(psi, 1), 2:size(psi, 2)-2), target :: psi_cc1

      integer, parameter                                    :: in = 1  ! index for cells
      integer                                               :: nx

      nx  = size(u, in)
      if (size(b_cc, in) /= nx) call die("[solve_cg:solve] size b_cc and u mismatch")
      mag_cc = huge(1.)

      ! RK(N) with N .GE. 3 could be helpful for WENO3 ( this statement to be tested )
      !>
      !! Reference:Relativistic Hydrodynamics, L. Rezzolla, O. Zanotti
      !! ---------------------------------------------------------------------------
      !! L (or dtodx)--> discretization of spatial differential operator (Eq. 9.135)
      !! ---------------------------------------------------------------------------
      !! RK2 (Eq. 9.140)
      !! u^(1)   = u^(n) + \Delta t L(u^(n))
      !! u^(n+1) = 1/2 ( u^(n) + u^(1) + \Delta t L(u^(1)  )
      !! ---------------------------------------------------------------------------
      !! RK3 (Eq. 9.141)
      !! u^(1)   = u(n) + \Delta t L(u^(n))
      !! u^(2)   = 1/4 ( 3 u^(n) + u^(1) + \Delta t L(u^(1) ) )
      !! u^(n+1) = 1/3 u^(n) + 2/3 u^(2) + 2/3 \Delta t (u^(2))
      !! ---------------------------------------------------------------------------
      !<
      select case (h_solver)
         case ("rk2")
            call interpol(u, b_cc, psi, ql, qr, b_cc_l, b_cc_r, psi_l, psi_r)
            call riemann_wrap(ql, qr, b_cc_l, b_cc_r, psi_l, psi_r, flx, mag_cc, psi_cc) ! Now we advance the left and right states by a timestep.
            call du_db(u1, b1, psi1, half * dtodx)
            call interpol(u1, b1, psi1, ql(2:nx-2,:), qr(2:nx-2,:), b_cc_l(2:nx-2,:), b_cc_r(2:nx-2,:), psi_l(2:nx-2,:), psi_r(2:nx-2,:))
            call riemann_wrap(ql(2:nx-2,:), qr(2:nx-2,:), b_cc_l(2:nx-2,:), b_cc_r(2:nx-2,:), psi_l(2:nx-2,:), psi_r(2:nx-2,:), flx(2:nx-2,:), mag_cc(2:nx-2,:), psi_cc(2:nx-2,:)) ! second call for Riemann problem uses states evolved to half timestep
            call update  ! ToDo tell explicitly what range to update
#if RK_HIGH
         case ("rk3") ! " to be checked "
            call interpol(u, b_cc, psi, ql, qr, b_cc_l, b_cc_r, psi_l, psi_r)
            call riemann_wrap(ql, qr, b_cc_l, b_cc_r, psi_l, psi_r, flx, mag_cc, psi_cc) ! Now we advance the left and right states by a timestep.
            call du_db(u1, b1, psi1, oneq * dtodx)
            call interpol(u1, b1, psi1, ql(2:nx-2,:), qr(2:nx-2,:), b_cc_l(2:nx-2,:), b_cc_r(2:nx-2,:), psi_l(2:nx-2,:), psi_r(2:nx-2,:))
            call riemann_wrap(ql(2:nx-2,:), qr(2:nx-2,:), b_cc_l(2:nx-2,:), b_cc_r(2:nx-2,:), psi_l(2:nx-2,:), psi_r(2:nx-2),:, flx(2:nx-2,:), mag_cc(2:nx-2,:), psi_cc(2:nx-2,:)) ! second call for Riemann problem uses states evolved to half timestep
            call du_db(u2,b2,psi2, twot * dtodx)
            call interpol(u1, b1, psi1, ql(2:nx-2,:), qr(2:nx-2,:), b_cc_l(2:nx-2,:), b_cc_r(2:nx-2,:), psi_l(2:nx-2,:), psi_r(2:nx-2,:))
            call riemann_wrap(ql(2:nx-2,:), qr(2:nx-2,:), b_cc_l(2:nx-2,:), b_cc_r(2:nx-2,:), psi_l(2:nx-2,:), psi_r(2:nx-2,:), flx(2:nx-2,:), mag_cc(2:nx-2,:),:, psi_cc(2:nx-2,:))
            call update

            !! rk3 and most other possible schemes will need more sophisticated update routine.

            !! Beware: weno3, ppm (and other high spatial order of reconstruction) in conjunction with
            !! high order temporal integration schemes will need reliable, automatic estimate of required
            !! minimal number of guardcells.

#endif /* RK_HIGH */
         case ("muscl")
            call interpol(u,b_cc,psi,ql,qr,b_cc_l,b_cc_r,psi_l,psi_r)
            call musclflx(ql, b_cc_l, psi_l, flx_l, bclflx, psilflx)
            call musclflx(qr, b_cc_r, psi_r, flx_r, bcrflx, psirflx)
            call ulr_fluxes_qlr
            call riemann_wrap(ql(2:nx-2,:), qr(2:nx-2,:), b_cc_l(2:nx-2,:), b_cc_r(2:nx-2,:), psi_l(2:nx-2,:), psi_r(2:nx-2,:), flx(2:nx-2,:), mag_cc(2:nx-2,:), psi_cc(2:nx-2,:)) ! 2:nx-1 should be possible here
            call update
         case default
            call die("[sweeps:solve_cg:solve] No recognized solver")
      end select

   contains

      ! some shortcuts

!>
!! \brief calculate the change of cell state due to face fluxes (1-D)
!!
!! For cells 1 .. n we have 1 .. n-1 fluxes  (at all internal faces),
!! so we can reliably calculate the update for cells 2 .. (n-1).
!<

      subroutine du_db(u_new, b_new, psi_new, fac)

         use constants,  only: DIVB_HDC
         use global,     only: divB_0_method

         implicit none

         real, dimension(2:size(u,   1)-1, size(u,   2)), intent(out) :: u_new
         real, dimension(2:size(b_cc,1)-1, size(b_cc,2)), intent(out) :: b_new
         real, dimension(2:size(psi, 1)-1, size(psi, 2)), intent(out) :: psi_new
         real, intent(in) :: fac

         ! shape(flx) = shape(u) - [ 0, 1 ] = [ n_variables, nx-1 ]
         u_new = u(2:nx-1, :) + fac * (flx(:nx-2, :) - flx(2:, :))
         b_new = b_cc(2:nx-1, :) + fac * (mag_cc(:nx-2, :) - mag_cc(2:, :))

         if (divB_0_method == DIVB_HDC) then
            psi_new = psi(2:nx-1, :) + fac * (psi_cc(:nx-2, :) - psi_cc(2:, :))
         else
            psi_new = 0.
         endif

      end subroutine du_db

      subroutine ulr_fluxes_qlr

         use constants, only: DIVB_HDC
         use global,    only: divB_0_method

         implicit none

         ql(2:nx-2, :) = ql(2:nx-2, :) + half * dtodx * (flx_r(1:nx-3, :) - flx_l(2:nx-2, :))
         qr(2:nx-2, :) = qr(2:nx-2, :) + half * dtodx * (flx_r(2:nx-2, :) - flx_l(3:nx-1, :))

         b_cc_l(2:nx-2, :) = b_cc_l(2:nx-2, :) + half * dtodx * (bcrflx(1:nx-3, :) - bclflx(2:nx-2, :))
         b_cc_r(2:nx-2, :) = b_cc_r(2:nx-2, :) + half * dtodx * (bcrflx(2:nx-2, :) - bclflx(3:nx-1, :))

         if (divB_0_method == DIVB_HDC) then
            psi_l(2:nx-2, :) = psi_l(2:nx-2, :) + half * dtodx * (psirflx(1:nx-3, :) - psilflx(2:nx-2, :))
            psi_r(2:nx-2, :) = psi_r(2:nx-2, :) + half * dtodx * (psirflx(2:nx-2, :) - psilflx(3:nx-1, :))
         endif

      end subroutine ulr_fluxes_qlr

      subroutine musclflx(q, b_cc, psi, qf, b_ccf, psif)

         use constants,  only: half, xdim, ydim, zdim, DIVB_HDC, zero
         use fluidindex, only: flind
         use fluidtypes, only: component_fluid
         use global,     only: divB_0_method
         use hdc,        only: chspeed

         implicit none

         real, dimension(:,:), intent(in)  :: q
         real, dimension(:,:), intent(in)  :: b_cc
         real, dimension(:,:), intent(in)  :: psi
         real, dimension(:,:), intent(out) :: qf
         real, dimension(:,:), intent(out) :: b_ccf, psif

         class(component_fluid), pointer   :: fl
         real                              :: en
         integer                           :: ip, i

#ifdef ISO
#  error Isothermal EOS is not implemented yet in musclflx
#endif /* ISO */

         qf    = zero
         b_ccf = zero
         psif  = zero

         do ip = 1, flind%fluids

            fl => flind%all_fluids(ip)%fl

            do i = 1, size(q, 1)

               qf(i,fl%idn) = q(i,fl%idn)*q(i,fl%imx)
               if (fl%has_energy) then
                  if (fl%is_magnetized) then
                     qf(i,fl%imx) = q(i,fl%idn)*q(i,fl%imx)*q(i,fl%imx) + (q(i,fl%ien) + half*sum(b_cc(i,xdim:zdim)**2)) - b_cc(i,xdim)**2
                  else
                     qf(i,fl%imx) = q(i,fl%idn)*q(i,fl%imx)*q(i,fl%imx) + q(i,fl%ien)
                  endif
               else
                  qf(i,fl%imx) = q(i,fl%idn)*q(i,fl%imx)*q(i,fl%imx)
               endif
               if (fl%is_magnetized) then
                  qf(i,fl%imy) = q(i,fl%idn)*q(i,fl%imx)*q(i,fl%imy) - b_cc(i,xdim)*b_cc(i,ydim)
                  qf(i,fl%imz) = q(i,fl%idn)*q(i,fl%imx)*q(i,fl%imz) - b_cc(i,xdim)*b_cc(i,zdim)
                  b_ccf(i,ydim)  = b_cc(i,ydim)*q(i,fl%imx) - b_cc(i,xdim)*q(i,fl%imy)
                  b_ccf(i,zdim)  = b_cc(i,zdim)*q(i,fl%imx) - b_cc(i,xdim)*q(i,fl%imz)
                  if (divB_0_method .eq. DIVB_HDC) then
                     b_ccf(i,xdim) = psi(i,1)
                     psif(i,1)   = (chspeed**2)*b_cc(i,xdim)
                  endif
               else
                  qf(i,fl%imy) = q(i,fl%idn)*q(i,fl%imx)*q(i,fl%imy)
                  qf(i,fl%imz) = q(i,fl%idn)*q(i,fl%imx)*q(i,fl%imz)
               endif
               if (fl%has_energy) then
                  if (fl%is_magnetized) then
                     en = (q(i,fl%ien)/(fl%gam_1)) + half*q(i,fl%idn)*sum(q(i,fl%imx:fl%imz)**2) + half*sum(b_cc(i,xdim:zdim)**2)
                     qf(i,fl%ien) = (en + (q(i,fl%ien) + half*sum(b_cc(i,xdim:zdim)**2)))*q(i,fl%imx) - b_cc(i,xdim)*dot_product(q(i,fl%imx:fl%imz),b_cc(i,xdim:zdim))
                  else
                     en = (q(i,fl%ien)/(fl%gam_1)) + half*q(i,fl%idn)*sum(q(i,fl%imx:fl%imz)**2)
                     qf(i,fl%ien) = (en + (q(i,fl%ien)))*q(i,fl%imx)
                  endif
               endif

            enddo

         enddo

      end subroutine musclflx

      subroutine update !(weights)

         use constants,  only: xdim, ydim, zdim, DIVB_HDC
         use fluidindex, only: flind
         use global,     only: divB_0_method

         implicit none

!       real, optional, dimension(:), intent(in) :: weights

!       real, dimension(:), allocatable :: w
         integer :: iend  !< last component of any fluid (i.e. exclude CR or tracers here)


         iend = flind%all_fluids(flind%fluids)%fl%end

!!$       if (present(weights)) then
!!$          allocate(w(size(weights)))
!!$          w = weights/sum(weights)
!!$       else
!!$          allocate(w(1))
!!$          w(1) = 1.
!!$       endif

         u(3:nx-2, :iend) = u(3:nx-2, :iend) + dtodx * (flx(2:nx-3, :iend) - flx(3:nx-2, :iend)) ! * w(1)
!       if (size(w)>=2) u(2:nx, :iend) = u(2:nx, :iend) + w(2) * du1(2:nx, :iend)
!       if (size(w)>=3) u(2:nx, :iend) = u(2:nx, :iend) + w(3) * du2(2:nx, :iend)
!       if (size(w)>=4) u(2:nx, :iend) = u(2:nx, :iend) + w(4) * du3(2:nx, :iend)
!       u(1, :iend) = u(2, :iend)
!       u(nx, :iend) = u(nx-1, :iend)

         b_cc(3:nx-2, ydim:zdim) = b_cc(3:nx-2, ydim:zdim) + dtodx * (mag_cc(2:nx-3, ydim:zdim) - mag_cc(3:nx-2, ydim:zdim)) ! * w(1)
!       if (size(w)>=2)  b_cc(2:nx, ydim:zdim) = b_cc(2:nx, ydim:zdim) + w(2) * db1(2:nx, ydim:zdim)
!       if (size(w)>=3)  b_cc(2:nx, ydim:zdim) = b_cc(2:nx, ydim:zdim) + w(3) * db2(2:nx, ydim:zdim)
!       if (size(w)>=4)  b_cc(2:nx, ydim:zdim) = b_cc(2:nx, ydim:zdim) + w(4) * db3(2:nx, ydim:zdim)

         if (divB_0_method == DIVB_HDC) then

            b_cc(3:nx-2, xdim) = b_cc(3:nx-2, xdim) + dtodx * (mag_cc(2:nx-3, xdim) - mag_cc(3:nx-2, xdim)) ! * w(1)
!          if (size(w)>=2)  b_cc(2:nx, xdim) = b_cc(2:nx, xdim) + w(2) * db1(2:nx, xdim)
!          if (size(w)>=3)  b_cc(2:nx, xdim) = b_cc(2:nx, xdim) + w(3) * db2(2:nx, xdim)
!          if (size(w)>=4)  b_cc(2:nx, xdim) = b_cc(2:nx, xdim) + w(4) * db3(2:nx, xdim)

            psi(3:nx-2, 1) = psi(3:nx-2, 1) + dtodx * (psi_cc(2:nx-3, 1) - psi_cc(3:nx-2, 1)) ! * w(1)
!          if (size(w)>=2)  psi_cc(2:nx, 1) = psi_cc(2:nx, 1) + w(2) * dpsi1(2:nx, 1)
!          if (size(w)>=3)  psi_cc(2:nx, 1) = psi_cc(2:nx, 1) + w(3) * dpsi2(2:nx, 1)
!          if (size(w)>=4)  psi_cc(2:nx, 1) = psi_cc(2:nx, 1) + w(4) * dpsi3(2:nx, 1)
!          psi(1, :) = psi(2, :)
!          psi(nx, :) = psi(nx-1, :)

         endif

!       b_cc(1, :)  = b_cc(2, :)
!       b_cc(nx, :) = b_cc(nx-1, :)

!       deallocate(w)

      end subroutine update

   end subroutine solve

end module solvecg
