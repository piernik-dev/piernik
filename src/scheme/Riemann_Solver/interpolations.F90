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
!--------------------------------------------------------------------------------------------------------------

#include "piernik.h"

!>
!!  \brief This module implements interpolation of left and right face states for cell-centered vectors of conservative variables.
!<

module interpolations

! pulled by ANY

   implicit none

   private
   public :: set_interpolations, interpol

   interface

      subroutine interpolation(prim_var, prim_var_l, prim_var_r, f_limiter)

         use fluxlimiters, only: limiter

         implicit none

         real, dimension(:,:),        intent(in)  :: prim_var
         real, dimension(:,:),        intent(out) :: prim_var_l
         real, dimension(:,:),        intent(out) :: prim_var_r
         procedure(limiter), pointer, intent(in)  :: f_limiter

      end subroutine interpolation

   end interface

   procedure(interpolation), pointer :: interp => null()

contains

!>
!! \brief Convert a vector of conservative variables to primitive ones.
!<

   function utoq(u, b_cc) result(q)

      use constants,  only: half, xdim, zdim, I_ONE
      use fluidindex, only: flind
      use fluidtypes, only: component_fluid
      use func,       only: ekin

      implicit none

      real, dimension(:,:),           intent(in) :: u
      real, dimension(:,:), optional, intent(in) :: b_cc

      real, dimension(size(u, 1), size(u, 2)) :: q
      integer                                 :: p
      class(component_fluid), pointer         :: fl

      do p = 1, flind%fluids
         fl => flind%all_fluids(p)%fl

         q(:, fl%idn) =  u(:, fl%idn)
         q(:, fl%imx) =  u(:, fl%imx)/u(:, fl%idn)
         q(:, fl%imy) =  u(:, fl%imy)/u(:, fl%idn)
         q(:, fl%imz) =  u(:, fl%imz)/u(:, fl%idn)
         ! J.CoPhy 208 (2005), Pg 317, Eq. 2. Gas pressure: p = (gamma-1)*(e-half*rho*v^2-half*B^2) and Total pressure: p_T = p + half*B^2. (1) and (2) are markers for HD and MHD.
         if (fl%has_energy) then
            q(:, fl%ien) =  fl%gam_1*(u(:, fl%ien) - ekin(u(:, fl%imx), u(:, fl%imy), u(:, fl%imz), u(:, fl%idn))) ! Primitive variable for gas pressure (p) without magnetic fields. (1)
            if (fl%is_magnetized) then
               q(:, fl%ien) =  q(:, fl%ien) - half*fl%gam_1*sum(b_cc(:, xdim:zdim)**2, dim=2) ! Primitive variable for gas pressure (p) with magnetic fields. The requirement of total pressure is dealt in the fluxes and hlld routines. (2)
            endif
         endif

      enddo

      associate (iend => flind%all_fluids(flind%fluids)%fl%end)
         if (iend < flind%all) q(:, iend + I_ONE:) = u(:, iend + I_ONE:)
      end associate

   end function utoq

!<
!! \brief Apply chosen interpolation scheme to obtain estimates of left and right state for the Riemann solver.
!>

   subroutine interpol(u, ql, qr, bcc, bccl, bccr)

      use fluxlimiters, only: flimiter, blimiter

      implicit none

      real, dimension(:,:), intent(in)     :: u
      real, dimension(:,:), intent(out)    :: ql
      real, dimension(:,:), intent(out)    :: qr

      real, dimension(:,:), intent(in), optional     :: bcc
      real, dimension(:,:), intent(out), optional    :: bccl
      real, dimension(:,:), intent(out), optional    :: bccr

      real, dimension(size(u, 1), size(u, 2)) :: q

      q = utoq(u, bcc)
      call interp(q,   ql,   qr,   flimiter)
      if (present(bcc)) call interp(bcc, bccl, bccr, blimiter)

   end subroutine interpol

!>
!! \brief Interpret and set desired interpolation scheme.
!<

   subroutine set_interpolations

      use dataio_pub,   only: msg, die
      use global,       only: interpol_str, limiter, limiter_b
      use fluxlimiters, only: set_limiters

      implicit none

      if (associated(interp)) call die("[interpolations:set_interpolations] interp already associated")

      select case (interpol_str)
         case ('linear', 'LINEAR', 'lin', '1')
            interp => linear
            call set_limiters(limiter, limiter_b)
         case ('weno3', 'WENO3', 'wo3', '2')
            interp => weno3
         case default
            write(msg, '(3a)') "[interpolations:set_interpolations] unknown interpolation '", interpol_str, "'"
            call die(msg)
            interp => null()
      end select

   end subroutine set_interpolations

!>
!! \brief Linear interpolation scheme for estimating left and right states on cell interfaces.
!!
!! Note that the arrays with left and right states are one cell shorter because we don't have enough knowledge to compute complete interface states at both vector ends.
!!
!! OPT: it looks like vectors longer than 1000 elements tend to be slower than shorter ones
!<

   subroutine linear(q, ql, qr, f_limiter)

      use fluxlimiters, only: limiter
      use domain,       only: dom
      use constants,    only: half, GEO_XYZ
      use dataio_pub,   only: die, msg

      implicit none

      real, dimension(:,:),        intent(in)  :: q
      real, dimension(:,:),        intent(out) :: ql
      real, dimension(:,:),        intent(out) :: qr
      procedure(limiter), pointer, intent(in)  :: f_limiter

      real, dimension(size(q, 1), size(q, 2)) :: dq_interp
      integer                                 :: n
      integer, parameter                      :: in = 1  ! index for cells

      if (dom%geometry_type /= GEO_XYZ) call die("[interpolations:linear] non-cartesian geometry not implemented yet.")
      if (size(q, in) - size(ql, in) /= 1) then
         write(msg, '(2(a,2i7),a)')"[interpolations:linear] face vector of wrong length: ", size(q, in), size(ql, in), " (expecting: ", size(q, in), size(q, in)-1, ")"
         call die(msg)
      endif
      if (any(shape(ql) /= shape(qr))) call die("[interpolations:linear] face vectors of different lengths")

      n = size(q, in)

      dq_interp = half * f_limiter(q) ! interpolate by half of cell

      ! interpolate from i-th cell center to its right face to get left state at interface i
      ql = q(:n-1, :) + dq_interp(:n-1, :)
      ! the right state at interface i comes from cell (i+1) and we have to obtain interpolation towards its left face
      qr = q(2:, :) - dq_interp(2:, :)

   end subroutine linear

!>
!! \brief Weighted Essentially Non-oscillatory 3rd order (WENO3) interpolation.
!! Based on Yamaleev N. K., Carpenter M. H., Journal of Computational Physics 228 (2009) 3025-3047.
!<

   subroutine weno3(q, ql, qr, f_limiter)

      use constants,    only: GEO_XYZ, onet, twot, one, two
      use dataio_pub,   only: die, msg
      use domain,       only: dom, is_refined
      use fluxlimiters, only: limiter
      use global,       only: w_epsilon

      implicit none

      real, dimension(:,:),        intent(in)  :: q
      real, dimension(:,:),        intent(out) :: ql
      real, dimension(:,:),        intent(out) :: qr
      procedure(limiter), pointer, intent(in)  :: f_limiter

      ! WENO3 definitions

      real                                     :: w0, w1
      real                                     :: alpha0, alpha1
      real                                     :: beta0, beta1
      real                                     :: flux0, flux1
      real                                     :: tau
      real                                     :: d0,d1
      !real                                     :: Delta_xi, Delta_xi2
      integer                                  :: n
      integer, parameter                       :: in = 1  ! index for cells
      integer                                  :: i, v

      if (is_refined .and. dom%nb <= 3) call die("[interpolations:weno3] AMR and WENO3 lead to inaccurate results when dom%nb <= 3")
      if (dom%geometry_type /= GEO_XYZ) call die("[interpolations:weno3] non-cartesian geometry not implemented yet.")
      if (size(q, in) - size(ql, in) /= 1) then
         write(msg, '(2(a,2i7),a)')"[interpolations:weno3] face vector of wrong length: ", size(q, in), size(ql, in), " (expecting: ", size(q, in), size(q, in)-1, ")"
         call die(msg)
      endif
      if (any(shape(ql) /= shape(qr))) call die("[interpolations:weno3] face vectors of different lengths")

      n = size(q, in)

      ! Eq. 19
      d0 = twot
      d1 = onet

      do v = lbound(q,2), ubound(q,2)
         do i = 2, n-1

            ! Eq. 20
            beta0 = (q(i+1, v) - q(i, v))**2  ! beta_0(j) = (u(j+1) - u(j))**2
            beta1 = (q(i, v) - q(i-1, v))**2  ! beta_1(j) = (u(j) - u(j-1))**2

            ! Eq. 22
            tau = (q(i+1, v) - two*q(i, v) + q(i-1, v))**2  ! tau(j) = (u(j+1) - 2 * u(j) + u(j-1))**2

            !>
            !! The WENO scheme is self-similar. The same applies to ESWENO.
            !! The grid spacing \Delta x is replaced with the grid spacing
            !! in the computational domain \Delta xi = 1/j, where j is the
            !! total number of gird cells.
            !<

            ! \Delta xi is mentioned before Eq. 63
            ! Delta_xi = one/n
            ! Delta_xi2 = Delta_xi*Delta_xi

            !>
            !! Sec 4.3
            !! The value of the tuning parameter "epsilon" is based on
            !! truncation error analysis. In Eq. 21 (or Eq. 19) this
            !! term epsilon is added with beta_r. Hence, it should be
            !! scaled consistently. The terms beta_r are ~ u_\xi^2 \Delta xi^2
            !! near smooth regions and ~ u^2 near unresolved regions.
            !! Therefore, epsilon can be chosen as:
            !! epsilon = max( L1norm(u_0^2), L1norm(u_0_\xi^2))*\Delta xi^2,
            !! where \xi != \xi_d, and L1_norm = sum ( abs(x(:)) ). The terms
            !! u_0 is the "initial condition", and u_0_\xi^2 is the "initial condition"
            !! discarding the set of points \xi_d where it is discontinuous.
            !! The term epsilon is calculated only once, and same value is used
            !! during the entire calculation to reduce CPU time. But, this could be
            !! relaxed, and possibly epsilon may be calculated at every time step.
            !<

            ! Eq. 65
            !w_epsilon   = max( sum( abs( q**2 ) ), one  )*Delta_xi2 ! Ad-hoc for making it work, but not correct!
            !w_epsilon   = max( sum( abs( u0^2  )  )  , sum( abs( u0_xi^2 ) ) )*Delta_xi2 ! Actual formulation, \xi != \xi_d --> points of discontinuity.

            ! For simplicity here we use just constant small epsilon from NUMERICAL_SETUP.

            ! Eq. 21 improved version compared to Eq. 19
            alpha0 = d0*(one + tau/(w_epsilon + beta0))  ! alpha_r(j) = d_r(j) * ( 1. + tau(j) / (epsilon + beta_r(j)) ) , r = 0, 1
            alpha1 = d1*(one + tau/(w_epsilon + beta1))

            ! Eq. 18
            w0 = alpha0/(alpha0 + alpha1)  ! w_r(j) = alpha_r(j) / (alpha_0(j) + alpha_1(j)) , r = 0, 1
            w1 = alpha1/(alpha0 + alpha1)  ! w_r(j) is positioned at right face of u(j); w^r_{j+1/2}  in the paper

            ! Left state interpolation, Eq. 15

            flux0 = 0.5*(q(i, v) + q(i+1, v))     ! f_0(j) = 1/2 * (  q(j) + q(j+1)) ; face at position j+1/2
            flux1 = 0.5*(-q(i-1, v) + 3.0*q(i, v))  ! f_1(j) = 1/2 * (3*q(j) - q(j-1))

            ! WENO3 flux, Eq. 14
            ql(i, v) = w0*flux0 + w1*flux1  ! f_W(j) = w_0(j) * f_0(j) + w_1(j) * f_1(j) ; face at position j+1/2

            ! Right state interpolation, Eq. 15
            ! we have to construct the right state symmetrically

            alpha0 = d0*(one + tau/(w_epsilon + beta1))  ! alpha_r(j) = d_r(j) * ( 1. + tau(j) / (epsilon + beta_r(1-j)) ) , r = 0, 1
            alpha1 = d1*(one + tau/(w_epsilon + beta0))  ! we use opposite beta to interpolate from the opposite side

            w0 = alpha0/(alpha0 + alpha1)  ! w_r(j) = alpha_r(j) / (alpha_0(j) + alpha_1(j)) , r = 0, 1
            w1 = alpha1/(alpha0 + alpha1)

            flux0 = 0.5*(q(i, v) + q(i-1, v))          ! f_0(j) = 1/2 * (  q(j) + q(j-1)) ; face at position j+1/2 from the other side
            flux1 = 0.5*(-q(i+1, v) + 3.0*q(i, v))  ! f_1(j) = 1/2 * (3*q(j) - q(j+1))

            ! WENO3 flux, Eq. 14, already shifted
            qr(i-1, v) = w0*flux0 + w1*flux1  ! similar as for left state, but weights are constructed from symmetric stencil

            if (.false.) qr = f_limiter(q)  ! suppress compiler warning on argument needed for other interpolation scheme

         enddo
      enddo
      ! Q&D: fix for FPE
      ! ToDo: handle it properly
      ql(1, :) = q(1, :)
      qr(n-1, :) = q(n, :)

   end subroutine weno3

end module interpolations
