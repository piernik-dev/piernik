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

#include "piernik.def"

!>
!!  \brief This module implements interpolation of left and right face states for cell-centered vectors of conservative variables.
!<

module interpolations
! pulled by RIEMANN

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
!! \brief Convert a vector of conservarive variables to primitive ones.
!<

   function utoq(u, b_cc) result(q)

     use constants,  only: half, xdim, zdim
     use fluidindex, only: flind
     use fluidtypes, only: component_fluid
     use func,       only: ekin

     implicit none

     real, dimension(:,:),   intent(in)    :: u , b_cc

     real, dimension(size(u, 1), size(u, 2))  :: q
     integer  :: p
     class(component_fluid), pointer       :: fl

     do p = 1, flind%fluids
        fl => flind%all_fluids(p)%fl

        q(fl%idn, :) =  u(fl%idn, :)
        q(fl%imx, :) =  u(fl%imx, :)/u(fl%idn, :)
        q(fl%imy, :) =  u(fl%imy, :)/u(fl%idn, :)
        q(fl%imz, :) =  u(fl%imz, :)/u(fl%idn, :)
        ! J.CoPhy 208 (2005), Pg 317, Eq. 2. Gas pressure: p = (gamma-1)*(e-half*rho*v^2-half*B^2) and Total pressure: p_T = p + half*B^2. (1) and (2) are markers for HD and MHD.
        if (fl%has_energy) then
            q(fl%ien, :) =  fl%gam_1*(u(fl%ien, :) - ekin(u(fl%imx, :), u(fl%imy, :), u(fl%imz, :), u(fl%idn, :))) ! Primitive variable for gas pressure (p) without magnetic fields. (1)
            if (fl%is_magnetized) then
               q(fl%ien, :) =  q(fl%ien, :) - half*fl%gam_1*sum(b_cc(xdim:zdim, :)**2, dim=1) ! Primitive variable for gas pressure (p) with magnetic fields. The requirement of total pressure is dealt in the fluxes and hlld routines. (2)
            endif
        endif

     enddo

   end function utoq

!<
!! \brief Apply chosen interpolation scheme to obtain estimates of left and right state for the Riemann solver.
!>

  subroutine interpol(u, bcc, psi, ql, qr, bccl, bccr, psil, psir)

    use fluxlimiters, only: flimiter, blimiter

    implicit none

    real, dimension(:,:), intent(in)     :: u
    real, dimension(:,:), intent(out)    :: ql
    real, dimension(:,:), intent(out)    :: qr

    real, dimension(:,:), intent(in)     :: bcc
    real, dimension(:,:), intent(out)    :: bccl
    real, dimension(:,:), intent(out)    :: bccr

    real, dimension(:,:), intent(in)     :: psi
    real, dimension(:,:), intent(out)    :: psil
    real, dimension(:,:), intent(out)    :: psir

    real, dimension(size(u, 1), size(u, 2)) :: q

    q = utoq(u, bcc)
    call interp(q,   ql,   qr,   flimiter)
    call interp(bcc, bccl, bccr, blimiter)
    call interp(psi, psil, psir, blimiter)

  end subroutine interpol

!>
!! \brief Interpret and set desired interpolation scheme.
!<

  subroutine set_interpolations

    use dataio_pub, only: msg, die
    use global,     only: interpol_str

    implicit none

    if (associated(interp)) call die("[interpolations:set_interpolations] interp already associated")

    select case (interpol_str)
    case ('linear', 'LINEAR', 'lin', '1')
       interp => linear
    case('weno3', 'WENO3', 'wo3', '2')
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
    integer, parameter                      :: in = 2  ! index for cells

    if (dom%geometry_type /= GEO_XYZ) call die("[interpolations:linear] non-cartesian geometry not implemented yet.")
    if (size(q, in) - size(ql, in) /= 1) then
       write(msg, '(2(a,2i7),a)')"[interpolations:linear] face vector of wrong length: ", size(q, in), size(ql, in), " (expecting: ", size(q, in), size(q, in)-1, ")"
       call die(msg)
    endif
    if (any(shape(ql) /= shape(qr))) call die("[interpolations:linear] face vectors of different lengths")

    n = size(q, in)

    dq_interp = half * f_limiter(q) ! interpolate by half of cell

    ! interpolate from i-th cell center to its right face to get left state at interface i
    ql = q(:, :n-1) + dq_interp(:, :n-1)
    ! the right state at interface i comes from cell (i+1) and we have to obtain interpolation towards its left face
    qr = q(:, 2:) - dq_interp(:, 2:)

  end subroutine linear

!>
!! \brief Weighted Essentially Non-oscillatory 3rd order (WENO3) interpolation.
!! Based on Yamaleev N. K., Carpenter M. H., Journal of Computational Physics 228 (2009) 3025-3047.
!<
#undef EPS

  subroutine weno3(q, ql, qr, f_limiter)

#ifdef EPS
    use cg_leaves,    only: leaves
    use cg_list,      only: cg_list_element
#endif
    use domain,       only: dom
    use fluxlimiters, only: limiter
    use domain,       only: dom
    use constants,    only: GEO_XYZ
    use dataio_pub,   only: die, msg
    use constants,    only: one, two
    
    implicit none

    real, dimension(:,:),        intent(in)  :: q
    real, dimension(:,:),        intent(out) :: ql
    real, dimension(:,:),        intent(out) :: qr
    procedure(limiter), pointer, intent(in)  :: f_limiter

    ! WENO3 definitions
    ! Artur check this routine very carefully!

    real, dimension(size(q,1),size(q,2))     :: w0, w1
    real, dimension(size(q,1),size(q,2))     :: alpha0, alpha1
    real, dimension(size(q,1),size(q,2))     :: beta0, beta1
    real, dimension(size(q,1),size(q,2))     :: flux0l, flux1l, flux0r, flux1r
    real, dimension(size(q,1),size(q,2))     :: tau
    real                                     :: d0,d1
    real                                     :: epsilon

#ifdef EPS
    type(cg_list_element), pointer           :: cgl
#endif
    integer                                  :: n
    integer, parameter                       :: in = 2  ! index for cells

    
    if (dom%geometry_type /= GEO_XYZ) call die("[interpolations:linear] non-cartesian geometry not implemented yet.")
    if (size(q, in) - size(ql, in) /= 1) then
       write(msg, '(2(a,2i7),a)')"[interpolations:linear] face vector of wrong length: ", size(q, in), size(ql, in), " (expecting: ", size(q, in), size(q, in)-1, ")"
       call die(msg)
    endif
    if (any(shape(ql) /= shape(qr))) call die("[interpolations:linear] face vectors of different lengths")

    n = size(q, in)

    ! Eq. 19
    d0 = two/3.0
    d1 = one/3.0

    ! remove this, or make it dirty
    beta0  = 0.
    beta1  = 0.
    tau    = 0.
    flux0l = 0.
    flux0r = 0.
    flux1l = 0.
    flux1r = 0.

    ! Eq. 20
    beta0(:, :n-1) = (q(:, 2:) - q(:, :n-1))**2  ! beta_0(j) = (u(j+1) - u(j))**2
    beta1(:, 2:)   = (q(:, 2:) - q(:, :n-1))**2  ! beta_1(j) = (u(j) - u(j-1))**2

    ! Eq. 22
    tau(:, 2:n-1) = (q(:, 3:) - two*q(:, 2:n-1) + q(:, :n-2))**2  ! tau(j) = (u(j+1) - 2 * u(j) + u(j-1))**2

    epsilon = 1e-6 ! if not for this declaration, epsilon goes unints in alpha0 or alpha1

#ifdef EPS
    ! let's debug other things first
    ! epsilon depends on dx
    cgl => leaves%first
    do while(associated(cgl))
       ! AJG: I think this should be local parameter, depending only on current block and perhaps also on current direction
       epsilon = minval(cgl%cg%dl,mask=dom%has_dir) ! check for correctness
       epsilon = epsilon*epsilon
       cgl => cgl%nxt
    end do
#endif

    ! Eq. 21 improved version compared to Eq. 19
    alpha0 = d0*(one + tau/(epsilon + beta0))  ! alpha_r(j) = d_r(j) * ( 1. + tau(j) / (epsilon + beta_r(j)) ) , r = 0, 1
    alpha1 = d1*(one + tau/(epsilon + beta1))

    ! Eq. 18
    w0 = alpha0/(alpha0 + alpha1)  ! w_r(j) = alpha_r(j) / (alpha_0(j) + alpha_1(j)) , r = 0, 1
    w1 = alpha1/(alpha0 + alpha1)  ! w_r(j) is positioned at right face of u(j); w^r_{j+1/2}  in the paper

    ! Left state interpolation, Eq. 15
    
    flux0l(:, :n-1) = 0.5*(q(:, :n-1) + q(:, 2:))     ! f_0(j) = 1/2 * (  q(j) + q(j+1)) ; face at position j+1/2
    flux1l(:, 2:) = 0.5*(-q(:, :n-1) + 3.0*q(:, 2:))  ! f_1(j) = 1/2 * (3*q(j) - q(j-1))

    ! WENO3 flux, Eq. 14
    ql = w0(:, :n-1)*flux0l(:, :n-1) + w1(:, :n-1)*flux1l(:, :n-1)  ! f_W(j) = w_0(j) * f_0(j) + w_1(j) * f_1(j) ; face at position j+1/2

    ! Right state interpolation, Eq. 15


    ! it is messed up!
    ! don't know yet whether we can recycle weights for the left state
    ! most likely we have to construct them symmetrically

    flux0r(:, 2:) = 0.5*(q(:, 2:) + q(:,:n-1))          ! f_0(j) = 1/2 * (  q(j) + q(j-1)) ; face at position j+1/2
    flux1r(:, :n-1) = 0.5*(-q(:, 2:) + 3.0*q(:, :n-1))  ! f_1(j) = 1/2 * (3*q(j) - q(j+1))

    ! WENO3 flux, Eq. 14
    qr = w0(:, :n-1)*flux0r (:, :n-1)+ w1(:, :n-1)*flux1r(:, :n-1)

    ! Shift right state

    qr(:,:n-2) = qr(:, 2:n-1)

    ! Update interpolation for first and last points, may be redundant
    ql(:,1) = q(:,1)
    qr(:,n-1) = ql(:,n-1)
    
    if (.false.) qr = f_limiter(q)  ! suppress compiler worning on argument needed for other interpolation scheme

  end subroutine weno3

end module interpolations
