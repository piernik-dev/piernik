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

module unsplit_state_update

!>
!! \brief Applies multidimensional unsplit flux updates to conserved fields.
!!
!! Updates fluid and magnetic states from their directional flux arrays, and
!! advances the divergence-cleaning psi field from its face fluxes.
!<

   implicit none

   private
   public :: apply_flux, update_psi

contains

!>
!! \brief Apply the unsplit flux divergence to the fluid or magnetic state.
!!
!! Uses the fluid flux arrays when \p mag is false and the magnetic flux arrays
!! when it is true. Intermediate Runge-Kutta stages are written to the
!! corresponding half-step state; the final stage updates the primary state.
!! \param cg Grid block whose state and flux arrays are updated.
!! \param istep Current Runge-Kutta stage.
!! \param mag Select magnetic fields instead of fluid fields.
!<
   subroutine apply_flux(cg, istep, mag)

      use constants,        only: xdim, ydim, zdim, last_stage, rk_coef, uh_n, I_ONE, ndims, magh_n
      use domain,           only: dom
      use global,           only: integration_order, dt
      use grid_cont,        only: grid_container
      use named_array_list, only: wna

      implicit none

      type :: fxptr
         real, pointer :: flx(:,:,:,:)
      end type fxptr

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: istep
      logical,                       intent(in) :: mag

      logical                     :: active(ndims)
      integer                     :: L0(ndims), U0(ndims), L(ndims), U(ndims), shift(ndims)
      integer                     :: afdim, uhi, bhi
      real, pointer               :: T(:,:,:,:)
      type(fxptr)                 :: F(ndims)

      T => null()
      active = [ dom%has_dir(xdim), dom%has_dir(ydim), dom%has_dir(zdim) ]

      if (mag) then
         F(xdim)%flx => cg%bfx   ;  F(ydim)%flx => cg%bgy   ;  F(zdim)%flx => cg%bhz

         L0 = [ lbound(cg%w(wna%bi)%arr, 2), lbound(cg%w(wna%bi)%arr, 3), lbound(cg%w(wna%bi)%arr, 4) ]
         U0 = [ ubound(cg%w(wna%bi)%arr, 2), ubound(cg%w(wna%bi)%arr, 3), ubound(cg%w(wna%bi)%arr, 4) ]

         bhi = wna%ind(magh_n)
         if (istep == last_stage(integration_order) .or. integration_order == I_ONE) then
            T => cg%w(wna%bi)%arr
         else
            cg%w(bhi)%arr(:,:,:,:) = cg%w(wna%bi)%arr(:,:,:,:)
            T => cg%w(bhi)%arr
         endif
      else
         F(xdim)%flx => cg%fx   ;  F(ydim)%flx => cg%gy   ;  F(zdim)%flx => cg%hz

         L0 = [ lbound(cg%w(wna%fi)%arr, 2), lbound(cg%w(wna%fi)%arr, 3), lbound(cg%w(wna%fi)%arr, 4) ]
         U0 = [ ubound(cg%w(wna%fi)%arr, 2), ubound(cg%w(wna%fi)%arr, 3), ubound(cg%w(wna%fi)%arr, 4) ]

         uhi = wna%ind(uh_n)
         if (istep == last_stage(integration_order) .or. integration_order == I_ONE) then
            T => cg%w(wna%fi)%arr
         else
            cg%w(uhi)%arr(:,:,:,:) = cg%w(wna%fi)%arr(:,:,:,:)
            T => cg%w(uhi)%arr
         endif
      endif

      do afdim = xdim, zdim
         if (.not. active(afdim)) cycle

         call bounds_for_flux(L0, U0, active, afdim, L, U)

         shift = 0
         shift(afdim) = I_ONE
         T(:, L(xdim):U(xdim), L(ydim):U(ydim), L(zdim):U(zdim)) = T(:, L(xdim):U(xdim), L(ydim):U(ydim), L(zdim):U(zdim)) &
              + dt / cg%dl(afdim) * rk_coef(istep) * ( &
              F(afdim)%flx(:, L(xdim):U(xdim), L(ydim):U(ydim), L(zdim):U(zdim)) - &
              F(afdim)%flx(:, L(xdim)+shift(xdim):U(xdim)+shift(xdim), &
              &               L(ydim)+shift(ydim):U(ydim)+shift(ydim), &
              &               L(zdim)+shift(zdim):U(zdim)+shift(zdim)) )
      enddo

   end subroutine apply_flux

!>
!! \brief Apply the unsplit flux divergence to the divergence-cleaning psi field.
!!
!! Uses the psi interface fluxes and the same Runge-Kutta stage convention as
!! \ref apply_flux.
!! \param cg Grid block whose psi field and fluxes are updated.
!! \param istep Current Runge-Kutta stage.
!<
   subroutine update_psi(cg, istep)

      use constants,        only: xdim, ydim, zdim, last_stage, rk_coef, I_ONE, ndims, psi_n, psih_n
      use domain,           only: dom
      use global,           only: integration_order, dt
      use grid_cont,        only: grid_container
      use named_array_list, only: qna

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: istep

      logical         :: active(ndims)
      integer         :: L0(ndims), U0(ndims), L(ndims), U(ndims), shift(ndims)
      integer         :: afdim, psihi, psii
      real, pointer   :: TP(:,:,:)

      TP => null()
      active = [ dom%has_dir(xdim), dom%has_dir(ydim), dom%has_dir(zdim) ]

      psii = qna%ind(psi_n)
      psihi = qna%ind(psih_n)
      L0 = [ lbound(cg%q(psii)%arr, 1), lbound(cg%q(psii)%arr, 2), lbound(cg%q(psii)%arr, 3) ]
      U0 = [ ubound(cg%q(psii)%arr, 1), ubound(cg%q(psii)%arr, 2), ubound(cg%q(psii)%arr, 3) ]

      if (istep == last_stage(integration_order) .or. integration_order == I_ONE) then
         TP => cg%q(psii)%arr
      else
         cg%q(psihi)%arr(:,:,:) = cg%q(psii)%arr(:,:,:)
         TP => cg%q(psihi)%arr
      endif

      do afdim = xdim, zdim
         if (.not. active(afdim)) cycle
         call bounds_for_flux(L0, U0, active, afdim, L, U)
         shift = 0
         shift(afdim) = I_ONE
         TP(L(xdim):U(xdim), L(ydim):U(ydim), L(zdim):U(zdim)) = TP(L(xdim):U(xdim), L(ydim):U(ydim), L(zdim):U(zdim)) &
              + dt / cg%dl(afdim) * rk_coef(istep) * ( &
              cg%psiflx(afdim, L(xdim):U(xdim), L(ydim):U(ydim), L(zdim):U(zdim)) - &
              cg%psiflx(afdim, L(xdim)+shift(xdim):U(xdim)+shift(xdim), &
              &                 L(ydim)+shift(ydim):U(ydim)+shift(ydim), &
                  &                 L(zdim)+shift(zdim):U(zdim)+shift(zdim)) )
      enddo

   end subroutine update_psi

!>
!! \brief Derive the writable cell bounds for one directional flux update.
!!
!! Excludes outer guard cells in active directions and trims transverse bounds
!! so that the update uses only cells covered by the directional fluxes.
!<
   subroutine bounds_for_flux(L0, U0, active, afdim, L, U)

      use constants, only: xdim, zdim, I_ONE, ndims
      use domain,    only: dom

      implicit none

      integer, intent(in)  :: L0(ndims), U0(ndims)
      logical, intent(in)  :: active(ndims)
      integer, intent(in)  :: afdim
      integer, intent(out) :: L(ndims), U(ndims)

      integer :: d, nb_1

      L = L0
      U = U0
      nb_1 = dom%nb - I_ONE

      do d = xdim, zdim
         if (active(d)) then
            L(d) = L(d) + I_ONE
            U(d) = U(d) - I_ONE
            if (d /= afdim) then
               L(d) = L(d) + nb_1
               U(d) = U(d) - nb_1
            endif
         endif
      enddo

   end subroutine bounds_for_flux

end module unsplit_state_update
