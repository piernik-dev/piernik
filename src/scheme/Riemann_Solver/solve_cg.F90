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

      use bfc_bcc,          only: interpolate_mag_field
      use constants,        only: pdims, xdim, zdim, ORTHO1, ORTHO2, LO, HI, psi_n, uh_n, magh_n, psih_n, INVALID, GEO_XYZ, I_ZERO, I_ONE, rk_coef
      use dataio_pub,       only: die
      use domain,           only: dom, is_refined
      use fluidindex,       only: flind, iarr_all_dn, iarr_all_mx, iarr_all_swp, iarr_mag_swp
      use global,           only: dt, force_cc_mag, use_fargo
      use grid_cont,        only: grid_container
      use hlld,             only: psidim
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

      real, dimension(cg%n_(ddim), size(cg%u,1)) :: u
      real, dimension(cg%n_(ddim), xdim:zdim)    :: b
      real, dimension(cg%n_(ddim)), target       :: psi
      real, dimension(:,:), pointer              :: pu, pu0, pb, pb0
      integer                                    :: i1, i2
      real, dimension(:), pointer                :: ppsi, ppsi0
      integer                                    :: psii, uhi, bhi, psihi
      real, dimension(size(u,1),size(u,2))       :: u0, u1
      real, dimension(size(b,1),size(b,2)+1)     :: b0, b1  ! Bx, By, Bz, psi
      real, dimension(size(u,1), flind%fluids), target :: vx
      integer(kind=4)                            :: nmag, i

      ! is_multicg should be safe
      if (.false.) write(0,*) present(fargo_vel) ! suppress compiler warning on unused argument
      if (use_fargo) call die("[solve_cg:solve_cg] Fargo is not yet enabled for Riemann")
      if (is_refined) call die("[solve_cg:solve_cg] This Rieman solver is not compatible with mesh refinements yet!")
      if (dom%geometry_type /= GEO_XYZ) call die("[solve_cg:solve_cg] Non-cartesian geometry is not implemented yet in this Riemann solver.")

      uhi = wna%ind(uh_n)
      bhi = wna%ind(magh_n)
      nmag = I_ZERO
      do i = 1, flind%fluids
         if (flind%all_fluids(i)%fl%is_magnetized) nmag = nmag + I_ONE
      enddo
      if (nmag > 1) call die("[solve_cg:solve_cg] At most one magnetized fluid is implemented")

      psii  = INVALID
      psihi = INVALID
      if (qna%exists(psi_n)) then
         psii = qna%ind(psi_n)
         psihi = qna%ind(psih_n)
      endif
      psi = 0.

      ppsi => psi ! suppress compiler complains on possibly uninitialized pointer

      call prepare_sources(cg)

      do i2 = cg%ijkse(pdims(ddim, ORTHO2), LO), cg%ijkse(pdims(ddim, ORTHO2), HI)
         do i1 = cg%ijkse(pdims(ddim, ORTHO1), LO), cg%ijkse(pdims(ddim, ORTHO1), HI)

            ! transposition for compatibility with RTVD-based routines
            pu0 => cg%w(uhi)%get_sweep(ddim,i1,i2)
            u0(:, iarr_all_swp(ddim,:)) = transpose(pu0(:,:))
            pu => cg%w(wna%fi)%get_sweep(ddim,i1,i2)
            u(:, iarr_all_swp(ddim,:)) = transpose(pu(:,:))

            pb => cg%w(wna%bi)%get_sweep(ddim,i1,i2)
            if (force_cc_mag) then
               pb0 => cg%w(bhi)%get_sweep(ddim,i1,i2)
               b0(:, iarr_mag_swp(ddim,:)) = transpose(pb0(:,:))
               b(:, iarr_mag_swp(ddim,:)) = transpose(pb(:,:))
            else
               ! For CT we have small inconsequence here: we don't call magfield
               ! and we discard transverse magnetic fluxes after first stage.
               ! Same applies to RTVD + CT.
               ! Beware: staggered grid will perhaps require magnetic boundary
               ! exchange with corners every time.
               b0(:, xdim:zdim) = interpolate_mag_field(ddim, cg, i1, i2, bhi)
               b(:, :) = interpolate_mag_field(ddim, cg, i1, i2, wna%bi)
            endif

            u1 = u
            b1(:, xdim:zdim) = b
            vx = u(:, iarr_all_mx) / u(:, iarr_all_dn) ! this may also be useful for gravitational acceleration
            if (psii /= INVALID) then
               ppsi0 => cg%q(psihi)%get_sweep(ddim,i1,i2)
               b0(:, psidim) = ppsi0(:)
               ppsi => cg%q(psii)%get_sweep(ddim,i1,i2)
               psi = ppsi
               b1(:, psidim) = psi

               call solve(u0, b0, u1, b1, rk_coef(istep) * dt/cg%dl(ddim))

            else
               call solve(u0, b0(:, xdim:zdim), u1, b1(:, xdim:zdim), rk_coef(istep) * dt/cg%dl(ddim))
            endif

            call all_sources(size(u, 1, kind=4), u, u1, b, cg, istep, ddim, i1, i2, rk_coef(istep) * dt, vx)
            ! See the results of Jeans test with RTVD and RIEMANN for estimate of accuracy.

#if defined COSM_RAYS && defined IONIZED
            if (size(u, 1) > 1) call limit_minimal_ecr(size(u, 1), u)
#endif /* COSM_RAYS && IONIZED */
            u(:,:) = u1(:,:)
            b(:,:) = b1(:, xdim:zdim)
            if (psii /= INVALID) psi = b1(:, psidim)

            pu(:,:) = transpose(u(:, iarr_all_swp(ddim,:)))
            if (force_cc_mag) pb(:,:) = transpose(b(:, iarr_mag_swp(ddim,:))) ! ToDo figure out how to manage CT energy fixup without extra storage
            if (psii /= INVALID) ppsi = psi
         enddo
      enddo

   end subroutine solve_cg

!! k-th interface is between k-th cell and (k+1)-th cell
!! We don't calculate n-th interface because it is as incomplete as 0-th interface

   subroutine solve(u, b_cc, u1, b1, dtodx)

      use constants,  only: DIVB_HDC, xdim, ydim, zdim
      use dataio_pub, only: die
      use global,     only: divB_0_method
      use hlld,           only: riemann_wrap, psidim
      use interpolations, only: interpol

      implicit none

      real, dimension(:,:), intent(in) :: u
      real, dimension(:,:), intent(in) :: b_cc
      real, dimension(:,:), intent(inout) :: u1
      real, dimension(:,:), intent(inout) :: b1
      real,                 intent(in)    :: dtodx

      ! left and right states at interfaces 1 .. n-1
      real, dimension(size(u,   1)-1, size(u,   2)), target :: ql, qr
      real, dimension(size(b_cc,1)-1, size(b_cc,2)), target :: b_cc_l, b_cc_r

      ! fluxes through interfaces 1 .. n-1
      real, dimension(size(u,   1)-1, size(u,   2)), target :: flx
      real, dimension(size(b_cc,1)-1, size(b_cc,2)), target :: mag_cc

      ! updates required for higher order of integration will likely have shorter length

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

      call interpol(u1, b1, ql, qr, b_cc_l, b_cc_r)
      call riemann_wrap(ql, qr, b_cc_l, b_cc_r, flx, mag_cc) ! Now we advance the left and right states by a timestep.

      u1(2:nx-1, :) = u(2:nx-1, :) + dtodx * (flx(:nx-2, :) - flx(2:, :))
      b1(2:nx-1, ydim:zdim) = b_cc(2:nx-1, ydim:zdim) + dtodx * (mag_cc(:nx-2, ydim:zdim) - mag_cc(2:, ydim:zdim))
      if (divB_0_method == DIVB_HDC) then
         b1(2:nx-1, xdim) = b_cc(2:nx-1, xdim) + dtodx * (mag_cc(:nx-2, xdim) - mag_cc(2:, xdim))
         b1(2:nx-1, psidim) = b_cc(2:nx-1, psidim) + dtodx * (mag_cc(:nx-2, psidim) - mag_cc(2:, psidim))
      else
         b1(2:nx-1, xdim) = b_cc(2:nx-1, xdim)
      endif

   end subroutine solve

end module solvecg
