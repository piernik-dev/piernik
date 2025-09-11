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

module unsplit_mag_modules

! pulled by ANY

   implicit none

   private
   public  :: solve_cg_ub

contains

   subroutine solve_cg_ub(cg,istep)
      use grid_cont,        only: grid_container
      use named_array_list, only: wna, qna
      use constants,        only: pdims, ORTHO1, ORTHO2, I_ONE, LO, HI, magh_n, uh_n, &
                                  psi_n, psih_n, psidim, cs_i2_n, first_stage, xdim, ydim, zdim, I_ONE
      use global,           only: integration_order
      use domain,           only: dom
      use fluidindex,       only: iarr_all_swp, iarr_mag_swp
      use fluxtypes,        only: ext_fluxes
      use unsplit_source,   only: apply_source
      use diagnostics,      only: my_allocate, my_deallocate

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: istep

      integer                                    :: i1, i2
      integer(kind=4)                            :: uhi, bhi, psii, psihi, ddim
      real, dimension(:,:),allocatable           :: u
      real, dimension(:,:),allocatable           :: b
      real, dimension(:,:),allocatable           :: b_psi                ! This will carry both b and psi so it will have one extra size in dim=2
      real, dimension(:,:), pointer              :: pu, pb
      real, dimension(:), pointer                :: ppsi
      real, dimension(:,:), pointer              :: pflux, pbflux,apsiflux
      real, dimension(:), pointer                :: ppsiflux
      real, dimension(:),   pointer              :: cs2
      real, dimension(:,:),allocatable           :: flux
      real, dimension(:,:),allocatable           :: bflux
      real, dimension(:,:),allocatable           :: tflux                 ! to temporarily store transpose of flux
      real, dimension(:,:),allocatable           :: tbflux                ! to temporarily store transpose of bflux
      type(ext_fluxes)                           :: eflx
      integer                                    :: i_cs_iso2

      uhi = wna%ind(uh_n)
      bhi = wna%ind(magh_n)

      psii = qna%ind(psi_n)
      psihi = qna%ind(psih_n)

      if (qna%exists(cs_i2_n)) then
         i_cs_iso2 = qna%ind(cs_i2_n)
      else
         i_cs_iso2 = -1
      endif
      cs2 => null()

      do ddim = xdim, zdim
         if (.not. dom%has_dir(ddim)) cycle

         call my_allocate(u, [cg%n_(ddim), size(cg%u,1, kind=4)])
         call my_allocate(b, [cg%n_(ddim), size(cg%b,1, kind=4)])
         call my_allocate(b_psi,  [size(b, 1, kind=4),         size(b, 2, kind=4) + I_ONE ])
         call my_allocate(flux,   [size(u, 1, kind=4) - I_ONE, size(u, 2, kind=4)])
         call my_allocate(tflux,  [size(u, 2, kind=4),         size(u, 1, kind=4)])
         call my_allocate(bflux,  [size(b, 1, kind=4) - I_ONE, size(b_psi, 2, kind=4)])
         call my_allocate(tbflux, [size(b_psi, 2, kind=4),     size(b, 1, kind=4)])

         do i2 = cg%ijkse(pdims(ddim, ORTHO2), LO), cg%ijkse(pdims(ddim, ORTHO2), HI)
            do i1 = cg%ijkse(pdims(ddim, ORTHO1), LO), cg%ijkse(pdims(ddim, ORTHO1), HI)

               if (ddim==xdim) then
                  pflux => cg%w(wna%xflx)%get_sweep(xdim, i1, i2)
                  pbflux => cg%w(wna%xbflx)%get_sweep(xdim, i1, i2)
                  apsiflux => cg%w(wna%psiflx)%get_sweep(xdim, i1, i2)
                  ppsiflux => apsiflux(xdim,:)
               else if (ddim==ydim) then
                  pflux => cg%w(wna%yflx)%get_sweep(ydim, i1, i2)
                  pbflux => cg%w(wna%ybflx)%get_sweep(ydim, i1, i2)
                  apsiflux => cg%w(wna%psiflx)%get_sweep(ydim, i1, i2)
                  ppsiflux => apsiflux(ydim,:)
               else if (ddim==zdim) then
                  pflux => cg%w(wna%zflx)%get_sweep(zdim, i1, i2)
                  pbflux => cg%w(wna%zbflx)%get_sweep(zdim, i1, i2)
                  apsiflux => cg%w(wna%psiflx)%get_sweep(zdim, i1, i2)
                  ppsiflux => apsiflux(zdim,:)
               endif
               pu   => cg%w(uhi)%get_sweep(ddim, i1, i2)
               pb   => cg%w(bhi)%get_sweep(ddim, i1, i2)
               ppsi => cg%q(psihi)%get_sweep(ddim, i1, i2)
               if (istep == first_stage(integration_order) .or. integration_order < 2 ) then
                  pu   => cg%w(wna%fi)%get_sweep(ddim, i1, i2)
                  pb   => cg%w(wna%bi)%get_sweep(ddim, i1, i2)
                  ppsi => cg%q(psii)%get_sweep(ddim, i1, i2)
               endif

               u(:, iarr_all_swp(ddim,:)) = transpose(pu(:,:))
               b(:, iarr_mag_swp(ddim,:)) = transpose(pb(:,:))

               b_psi(:, xdim:zdim) = b(:,:) ; b_psi(:,psidim) = ppsi(:)

               if (i_cs_iso2 > 0) cs2 => cg%q(i_cs_iso2)%get_sweep(ddim, i1, i2)

               call cg%set_fluxpointers(ddim, i1, i2, eflx)

               call solve(u, b_psi ,cs2, eflx, flux, bflux)

               call cg%save_outfluxes(ddim, i1, i2, eflx)

               tflux(:,2:) = transpose(flux(:, iarr_all_swp(ddim,:)))
               tflux(:,1) = 0.0
               pflux(:,:) = tflux

               tbflux(:,2:) = transpose(bflux(:, iarr_mag_swp(ddim,:)))
               tbflux(:,1) = 0
               tbflux(psidim,2:) = bflux(:,psidim)
               pbflux(:,:) = tbflux(xdim:zdim,:)
               ppsiflux(:) =  tbflux(psidim,:)

            enddo
         enddo

         call my_deallocate(u); call my_deallocate(flux); call my_deallocate(tflux)
         call my_deallocate(b); call my_deallocate(b_psi); call my_deallocate(tbflux)
         call my_deallocate(bflux)

      enddo

      call apply_flux(cg,istep,.true.)
      call apply_flux(cg,istep,.false.)
      call update_psi(cg,istep)
      call apply_source(cg,istep)
      nullify(cs2)

   end subroutine solve_cg_ub

   subroutine solve(ui, bi, cs2, eflx, flx, bflx)

      use constants,      only: DIVB_HDC
      use fluxtypes,      only: ext_fluxes
      use global,         only: divB_0_method
      use hlld,           only: riemann_wrap
      use interpolations, only: interpol
      use dataio_pub,     only: die

      implicit none

      real, dimension(:,:),        intent(in)    :: ui      !< cell-centered initial fluid states
      real, dimension(:,:),        intent(in)    :: bi      !< cell-centered initial magnetic field states (including psi field when necessary)
      real, dimension(:,:),        intent(inout) :: flx     !< cell-centered intermediate fluid states
      real, dimension(:,:),        intent(inout) :: bflx    !< cell-centered intermediate magnetic field states (including psi field when necessary)
      real, dimension(:), pointer, intent(in)    :: cs2     !< square of local isothermal sound speed
      type(ext_fluxes),            intent(inout) :: eflx    !< external fluxes

      ! left and right states at interfaces 1 .. n-1
      real, dimension(size(ui, 1)-1, size(ui, 2)), target :: ql, qr
      real, dimension(size(bi, 1)-1, size(bi, 2)), target :: bl, br

      ! updates required for higher order of integration will likely have shorter length

      bflx = huge(1.)

      call interpol(ui, ql, qr, bi, bl, br)
      call riemann_wrap(ql, qr, bl, br, cs2, flx, bflx) ! Now we advance the left and right states by a timestep.

      if (associated(eflx%li)) flx(eflx%li%index, :) = eflx%li%uflx
      if (associated(eflx%ri)) flx(eflx%ri%index, :) = eflx%ri%uflx
      if (associated(eflx%lo)) eflx%lo%uflx = flx(eflx%lo%index, :)
      if (associated(eflx%ro)) eflx%ro%uflx = flx(eflx%ro%index, :)

      if (divB_0_method == DIVB_HDC) then
         if (associated(eflx%li)) bflx(eflx%li%index, :) = eflx%li%bflx
         if (associated(eflx%ri)) bflx(eflx%ri%index, :) = eflx%ri%bflx
         if (associated(eflx%lo)) eflx%lo%bflx = bflx(eflx%lo%index, :)
         if (associated(eflx%ro)) eflx%ro%bflx = bflx(eflx%ro%index, :)
      else
         call die("[unsplit_mag_modules:solve] Unplit method is only implemented with Hyperbolic Divergence Cleaning")
      endif

   end subroutine solve

   subroutine apply_flux(cg, istep, mag)

      use domain,             only: dom
      use grid_cont,          only: grid_container
      use global,             only: integration_order, dt
      use named_array_list,   only: wna
      use constants,          only: xdim, ydim, zdim, last_stage, rk_coef, uh_n, I_ONE, ndims, magh_n

      implicit none

      type :: fxptr
         real, pointer :: flx(:,:,:,:)
      end type fxptr

      type(grid_container), pointer, intent(in)   :: cg
      integer,                       intent(in)   :: istep
      logical,                       intent(in)   :: mag

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
         if (istep==last_stage(integration_order) .or. integration_order==I_ONE) then
            T => cg%w(wna%bi)%arr
         else
            cg%w(bhi)%arr(:,:,:,:) = cg%w(wna%bi)%arr(:,:,:,:)
            T => cg%w(bhi)%arr
         endif
      else
         F(xdim)%flx => cg%fx   ;  F(ydim)%flx => cg%gy   ;  F(zdim)%flx => cg%hz

         L0 = [ lbound(cg%w(wna%fi)%arr,2), lbound(cg%w(wna%fi)%arr,3), lbound(cg%w(wna%fi)%arr,4) ]
         U0 = [ ubound(cg%w(wna%fi)%arr,2), ubound(cg%w(wna%fi)%arr,3), ubound(cg%w(wna%fi)%arr,4) ]

         uhi = wna%ind(uh_n)
         if (istep==last_stage(integration_order) .or. integration_order==I_ONE) then
            T => cg%w(wna%fi)%arr
         else
            cg%w(uhi)%arr(:,:,:,:) = cg%w(wna%fi)%arr(:,:,:,:)
            T => cg%w(uhi)%arr
         endif
      endif
      do afdim = xdim, zdim
         if (.not. active(afdim)) cycle

         call bounds_for_flux(L0,U0,active,afdim,L,U)

         shift = 0 ;  shift(afdim) = I_ONE
         T(:, L(xdim):U(xdim), L(ydim):U(ydim), L(zdim):U(zdim)) = T(:, L(xdim):U(xdim), L(ydim):U(ydim), L(zdim):U(zdim)) &
              + dt / cg%dl(afdim) * rk_coef(istep) * ( &
              F(afdim)%flx(:, L(xdim):U(xdim), L(ydim):U(ydim), L(zdim):U(zdim)) - &
              F(afdim)%flx(:, L(xdim)+shift(xdim):U(xdim)+shift(xdim), &
              &               L(ydim)+shift(ydim):U(ydim)+shift(ydim), &
              &               L(zdim)+shift(zdim):U(zdim)+shift(zdim)) )
      enddo

   end subroutine apply_flux

   subroutine update_psi(cg, istep)

      use domain,             only: dom
      use grid_cont,          only: grid_container
      use global,             only: integration_order, dt
      use named_array_list,   only: qna
      use constants,          only: xdim, ydim, zdim, last_stage, rk_coef, I_ONE, ndims, psi_n, psih_n

      implicit none

      type(grid_container), pointer, intent(in)   :: cg
      integer,                       intent(in)   :: istep

      logical                     :: active(ndims)
      integer                     :: L0(ndims), U0(ndims), L(ndims), U(ndims), shift(ndims)
      integer                     :: afdim, psihi, psii
      real, pointer               :: TP(:,:,:)

      TP => null()

      active = [ dom%has_dir(xdim), dom%has_dir(ydim), dom%has_dir(zdim) ]

      psii = qna%ind(psi_n)
      psihi = qna%ind(psih_n)
      L0 = [ lbound(cg%q(psii)%arr, 1), lbound(cg%q(psii)%arr, 2), lbound(cg%q(psii)%arr, 3) ]
      U0 = [ ubound(cg%q(psii)%arr, 1), ubound(cg%q(psii)%arr, 2), ubound(cg%q(psii)%arr, 3) ]

      if (istep==last_stage(integration_order) .or. integration_order==I_ONE) then
         TP => cg%q(psii)%arr
      else
         cg%q(psihi)%arr(:,:,:) = cg%q(psii)%arr(:,:,:)
         TP => cg%q(psihi)%arr
      endif

      do afdim = xdim, zdim
         if (.not. active(afdim)) cycle
         call bounds_for_flux(L0,U0,active,afdim,L,U)
         shift = 0 ;  shift(afdim) = I_ONE
         TP(L(xdim):U(xdim), L(ydim):U(ydim), L(zdim):U(zdim)) = TP(L(xdim):U(xdim), L(ydim):U(ydim), L(zdim):U(zdim)) &
              + dt / cg%dl(afdim) * rk_coef(istep) * ( &
              cg%psiflx(afdim,L(xdim):U(xdim), L(ydim):U(ydim), L(zdim):U(zdim)) - &
              cg%psiflx(afdim,L(xdim)+shift(xdim):U(xdim)+shift(xdim), &
              &               L(ydim)+shift(ydim):U(ydim)+shift(ydim), &
              &               L(zdim)+shift(zdim):U(zdim)+shift(zdim)) )
      enddo

   end subroutine update_psi

   subroutine bounds_for_flux(L0, U0, active, afdim, L, U)

      use constants, only: xdim, zdim, I_ONE, ndims
      use domain,    only: dom

      implicit none

      integer, intent(in)  :: L0(ndims), U0(ndims)   ! original bounds
      logical, intent(in)  :: active(ndims)          ! dom%has_dir flags
      integer, intent(in)  :: afdim                  ! direction we are updating (1,2,3)
      integer, intent(out) :: L(ndims), U(ndims)     ! returned bounds

      integer :: d,nb_1

      L = L0 ;  U = U0                     ! start from raw array bounds
      nb_1 = dom%nb - I_ONE

      do d = xdim, zdim
         if (active(d)) then               ! remove outer 1-cell ghosts
            L(d) = L(d) + I_ONE
            U(d) = U(d) - I_ONE
            if (d /= afdim) then           ! shrink transverse dirs by 3 extra
               L(d) = L(d) + nb_1
               U(d) = U(d) - nb_1
            endif
         endif
      enddo

   end subroutine bounds_for_flux

end module unsplit_mag_modules
