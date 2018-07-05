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

  subroutine weno3(q, ql, qr, f_limiter)

    use cg_leaves,    only: leaves
    use cg_list,      only: cg_list_element
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
   
    real, dimension(size(q,1),size(q,2))     :: w0, w1
    real, dimension(size(q,1),size(q,2))     :: alpha0, alpha1
    real, dimension(size(q,1),size(q,2))     :: beta0, beta1
    real, dimension(size(q,1),size(q,2))     :: flux0l, flux1l, flux0r, flux1r
    real, dimension(size(q,1),size(q,2))     :: tau
    real                                     :: d0,d1
    real                                     :: epsilon

    type(cg_list_element), pointer           :: cgl
    integer                                  :: n, ip1, im1
    integer, parameter                       :: in = 2  ! index for cells

    
    if (dom%geometry_type /= GEO_XYZ) call die("[interpolations:linear] non-cartesian geometry not implemented yet.")
    if (size(q, in) - size(ql, in) /= 1) then
       write(msg, '(2(a,2i7),a)')"[interpolations:linear] face vector of wrong length: ", size(q, in), size(ql, in), " (expecting: ", size(q, in), size(q, in)-1, ")"
       call die(msg)
    endif
    if (any(shape(ql) /= shape(qr))) call die("[interpolations:linear] face vectors of different lengths")

    n = size(q, in)

    ip1 = min(n,n+1)
    im1 = max(1,n-1)
    
    ! Eq. 19
    d0 = two/3.0
    d1 = one/3.0

    ! Eq. 20
    
    beta0 = (q(:,1:ip1) - q(:,1:n))
    beta0 = beta0*beta0
    beta1 = (q(:,1:n) - q(:,1:im1))
    beta1 = beta1*beta1
    
    ! Eq. 22
    tau = (q(:,1:ip1) - two*q(:,1:n) + q(:,1:im1))
    tau = tau*tau
   
    epsilon = 0.0 ! if not for this declaration, epsilon goes unints in alpha0 or alpha1
    ! epsilon depends on dx
    cgl => leaves%first
    do while(associated(cgl))
       epsilon = minval(cgl%cg%dl,mask=dom%has_dir) ! check for correctness
       epsilon = epsilon*epsilon
       cgl => cgl%nxt
    end do ! please check if this is correct. Constant value for epislon like 1.d-6 can be used for test.
       

    ! Eq. 21 improved version compared to Eq. 19
    
    alpha0 = d0*(one + tau/(epsilon + beta0))
    alpha1 = d1*(one + tau/(epsilon + beta1))

    ! Eq. 18
    w0 = alpha0/(alpha0 + alpha1)
    w1 = alpha1/(alpha0 + alpha1)
    
    ! Left state interpolation, Eq. 15
    
    flux0l = 0.5*(q(:,1:n) + q(:,1:ip1))
    flux1l = 0.5*(-q(:,1:im1) + 3.0*q(:,1:n))

    ! WENO3 flux, Eq. 14
    ql = w0*flux0l + w1*flux1l

    ! Right state interpolation, Eq. 15

    flux0r = 0.5*(q(:,1:n) + q(:,1:im1))
    flux1r = 0.5*(-q(:,1:ip1) + 3.0*q(:,1:n))

    ! WENO3 flux, Eq. 14
    qr = w0*flux0r + w1*flux1r
       
    ! Shift right state
    
    !qr(:,1:im1) = qr(:,2:n)

    ! Update interpolation for first and last points, may be redundant
    !ql(:,1) = q(:,1)
    !qr(:,n) = q(:,n)
       
    if (.false.) qr = f_limiter(q)  ! suppress compiler worning on argument needed for other interpolation scheme

       
  end subroutine weno3

end module interpolations
