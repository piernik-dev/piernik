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

module solvecg_unsplit

! pulled by ANY

   implicit none

   private
   public  :: solve_cg_unsplit

contains

! This routine has to conform to the interface defined in sweeps::sweep

   subroutine solve_cg_unsplit(cg, istep)

      use constants,             only: mag_n, GEO_XYZ
      use dataio_pub,            only: die
      use domain,                only: dom
      use fluidindex,            only: flind
      use grid_cont,             only: grid_container
      use named_array_list,      only: wna
      use sources,               only: prepare_sources
      use unsplit_mag_modules,   only: solve_cg_ub

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: istep     ! stage in the time integration scheme
      integer :: nmag, i

      if (dom%geometry_type /= GEO_XYZ) call die("[solve_cg_unsplit:solve_cg_unsplit] Non-cartesian geometry is not implemented yet in this Unsplit solver.")

      call prepare_sources(cg)

      if (wna%exists(mag_n)) then
         !call die("[solve_cg_unsplit:solve_cg_unsplit] Magnetic field is still unsafe for the flux named arrays")
         nmag = 0
         do i = 1, flind%fluids
            if (flind%all_fluids(i)%fl%is_magnetized) nmag = nmag + 1
         enddo
         if (nmag > 1) call die("[solve_cg_unsplit:solve_cg_unsplit] At most one magnetized fluid is implemented")
            call solve_cg_ub(cg, istep)
      else
         call solve_cg_u(cg, istep)
      endif

      cg%processed = .true.

   end subroutine solve_cg_unsplit

   subroutine solve_cg_u(cg, istep)

      use grid_cont,        only: grid_container
      use named_array_list, only: wna, qna
      use constants,        only: pdims, ORTHO1, ORTHO2, I_ONE, LO, HI, uh_n, cs_i2_n, first_stage, xdim, ydim, zdim
      use global,           only: integration_order
      use domain,           only: dom
      use fluidindex,       only: iarr_all_swp
      use fluxtypes,        only: ext_fluxes
      use unsplit_source,   only: apply_source
      use diagnostics,      only: my_allocate, my_deallocate

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: istep

      integer                                    :: i1, i2
      integer(kind=4)                            :: uhi, ddim
      real, dimension(:,:),allocatable           :: u
      real, dimension(:,:), pointer              :: pu
      real, dimension(:,:), pointer              :: pflux
      real, dimension(:),   pointer              :: cs2
      real, dimension(:,:),allocatable           :: flux
      real, dimension(:,:),allocatable           :: tflux
      type(ext_fluxes)                           :: eflx
      integer                                    :: i_cs_iso2

      uhi = wna%ind(uh_n)

      if (qna%exists(cs_i2_n)) then
         i_cs_iso2 = qna%ind(cs_i2_n)
      else
         i_cs_iso2 = -1
      endif
      cs2 => null()

      do ddim = xdim, zdim
         if (.not. dom%has_dir(ddim)) cycle
         call my_allocate(u, [cg%n_(ddim), size(cg%u, 1, kind=4)])
         call my_allocate(flux, [size(u, 1, kind=4)-I_ONE, size(u, 2, kind=4)])
         call my_allocate(tflux, [size(u, 2, kind=4), size(u, 1, kind=4)])
         do i2 = cg%ijkse(pdims(ddim, ORTHO2), LO), cg%ijkse(pdims(ddim, ORTHO2), HI)
            do i1 = cg%ijkse(pdims(ddim, ORTHO1), LO), cg%ijkse(pdims(ddim, ORTHO1), HI)
               if (ddim == xdim) then
                  pflux => cg%w(wna%xflx)%get_sweep(xdim, i1, i2)
               else if (ddim == ydim) then
                  pflux => cg%w(wna%yflx)%get_sweep(ydim, i1, i2)
               else if (ddim == zdim) then
                  pflux => cg%w(wna%zflx)%get_sweep(zdim, i1, i2)
               endif

               pu => cg%w(uhi)%get_sweep(ddim, i1, i2)
               if (istep == first_stage(integration_order) .or. integration_order < 2 ) pu => cg%w(wna%fi)%get_sweep(ddim, i1, i2)

               u(:, iarr_all_swp(ddim, :)) = transpose(pu(:,:))

               if (i_cs_iso2 > 0) cs2 => cg%q(i_cs_iso2)%get_sweep(ddim, i1, i2)

               call cg%set_fluxpointers(ddim, i1, i2, eflx)
               call solve_u(u,cs2, eflx, flux)

               call cg%save_outfluxes(ddim, i1, i2, eflx)
               tflux(:, 2:) = transpose(flux(:, iarr_all_swp(ddim,:)))
               tflux(:, 1) = 0.0
               pflux(:, :) = tflux
            enddo
         enddo

         call my_deallocate(u); call my_deallocate(flux); call my_deallocate(tflux)

      enddo
      call apply_flux(cg, istep)
      call apply_source(cg, istep)
      nullify(cs2)

   end subroutine solve_cg_u

   subroutine solve_u(ui, cs2, eflx, flx)

      use fluxtypes,      only: ext_fluxes
      use hlld,           only: riemann_wrap_u
      use interpolations, only: interpol
      use dataio_pub,     only: die

      implicit none

      real, dimension(:,:),        intent(in)    :: ui       !< cell-centered intermediate fluid states
      real, dimension(:), pointer, intent(in)    :: cs2     !< square of local isothermal sound speed
      type(ext_fluxes),            intent(inout) :: eflx    !< external fluxes
      real, dimension(:,:),        intent(inout) :: flx     !< Output flux of a 1D chain of a domain at a fixed ortho location of that dimension

      ! left and right states at interfaces 1 .. n-1
      real, dimension(size(ui, 1)-1, size(ui, 2)), target :: ql, qr

      ! updates required for higher order of integration will likely have shorter length
      if (size(flx,dim = 1) /= size(ui, 1) - 1 .or. size(flx, dim=2) /= size(ui, 2)  ) then
         call die("[solve_cg_unsplit:solve_u] flux array dimension does not match the expected dimensions")
      endif

      call interpol(ui, ql, qr)
      call riemann_wrap_u(ql, qr, cs2, flx) ! Now we advance the left and right states by a timestep.

      if (associated(eflx%li)) flx(eflx%li%index, :) = eflx%li%uflx
      if (associated(eflx%ri)) flx(eflx%ri%index, :) = eflx%ri%uflx
      if (associated(eflx%lo)) eflx%lo%uflx = flx(eflx%lo%index, :)
      if (associated(eflx%ro)) eflx%ro%uflx = flx(eflx%ro%index, :)

   end subroutine solve_u

   subroutine apply_flux(cg, istep)
      use domain,             only: dom
      use grid_cont,          only: grid_container
      use global,             only: integration_order, dt
      use named_array_list,   only: wna
      use constants,          only: xdim, ydim, zdim, last_stage, rk_coef, uh_n, I_ONE, ndims

      implicit none

      type :: fxptr
         real, pointer :: flx(:,:,:,:)
      end type fxptr

      type(grid_container), pointer, intent(in)   :: cg
      integer,                       intent(in)   :: istep

      logical                     :: active(ndims)
      integer                     :: L0(ndims), U0(ndims), L(ndims), U(ndims), shift(ndims)
      integer                     :: afdim, uhi
      real, pointer               :: T(:,:,:,:)
      type(fxptr)                 :: F(ndims)

      T=> null()
      active = [ dom%has_dir(xdim), dom%has_dir(ydim), dom%has_dir(zdim) ]

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

      do afdim = xdim, zdim
         if (.not. active(afdim)) cycle

         call bounds_for_flux(L0, U0, active, afdim, L, U)

         shift = 0 ;  shift(afdim) = I_ONE

         T(:, L(xdim):U(xdim), L(ydim):U(ydim), L(zdim):U(zdim)) = T(:, L(xdim):U(xdim), L(ydim):U(ydim), L(zdim):U(zdim)) &
              + dt / cg%dl(afdim) * rk_coef(istep) * ( &
              F(afdim)%flx(:, L(xdim):U(xdim), L(ydim):U(ydim), L(zdim):U(zdim)) - &
              F(afdim)%flx(:, L(xdim)+shift(xdim):U(xdim)+shift(xdim), &
              &               L(ydim)+shift(ydim):U(ydim)+shift(ydim), &
              &               L(zdim)+shift(zdim):U(zdim)+shift(zdim)) )
      enddo

   end subroutine apply_flux

   subroutine bounds_for_flux(L0, U0, active, afdim, L, U)

      use constants, only: xdim, zdim, I_ONE, ndims
      use domain,    only: dom

      implicit none

      integer, intent(in)  :: L0(ndims), U0(ndims)   ! original bounds
      logical, intent(in)  :: active(ndims)          ! dom%has_dir flags
      integer, intent(in)  :: afdim                  ! direction we are updating (1,2,3)
      integer, intent(out) :: L(ndims), U(ndims)     ! returned bounds

      integer :: d, nb_1

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

end module solvecg_unsplit
