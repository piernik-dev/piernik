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
!! \brief HLLD Riemann solver for ideal magnetohydrodynamics
!!
!! Varadarajan Parthasarathy, CAMK, Warszawa. 2015.
!! Dr. Artur Gawryszczak, CAMK, Warszawa.
!!
!! RK(N) with N .GE. 3 could be helpful for WENO3 ( this statement to be tested )
!!
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
!!
!! Energy fix up routines for CT and its related comments are not used in the current version.
!! The algorithm is simply present for experimental purposes.
!<

module solvecg_riemann

! pulled by ANY

   implicit none

   private
   public  :: solve_cg_riemann

contains

! This routine has to conform to the interface defined in sweeps::sweep

   subroutine solve_cg_riemann(cg, ddim, istep, fargo_vel)

      use constants,        only: mag_n, GEO_XYZ
      use dataio_pub,       only: die
      use domain,           only: dom
      use fluidindex,       only: flind
      use grid_cont,        only: grid_container
      use global,           only: use_fargo
      use named_array_list, only: wna
      use sources,          only: prepare_sources

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer(kind=4),               intent(in) :: ddim
      integer,                       intent(in) :: istep     ! stage in the time integration scheme
      integer(kind=4), optional,     intent(in) :: fargo_vel

      integer :: nmag, i

      if (.false.) write(0,*) present(fargo_vel) ! suppress compiler warning on unused argument
      if (use_fargo) call die("[solve_cg_riemann:solve_cg_riemann] Fargo is not yet enabled for Riemann")
      if (dom%geometry_type /= GEO_XYZ) call die("[solve_cg_riemann:solve_cg_riemann] Non-cartesian geometry is not implemented yet in this Riemann solver.")

      call prepare_sources(cg)

      if (wna%exists(mag_n)) then
         nmag = 0
         do i = 1, flind%fluids
            if (flind%all_fluids(i)%fl%is_magnetized) nmag = nmag + 1
         enddo
         if (nmag > 1) call die("[solve_cg_riemann:solve_cg_riemann] At most one magnetized fluid is implemented")
         call solve_cg_ub(cg, ddim, istep)
      else
         call solve_cg_u(cg, ddim, istep)
      endif

      cg%processed = .true.

   end subroutine solve_cg_riemann

!>
!! \brief Apply MHD update + source terms to a single grid container, rely on properly updated guardcells.
!!
!! \warning:
!! * exceptionally slow when blocks are 120^3 (sometimes an order of magnitude)
!! * significantly slower on 2520^2 when compared to 1260^2
!! * optimal vector length could be around 180 cells (from 90 to 360 are also good)
!<

   subroutine solve_cg_ub(cg, ddim, istep)

      use bfc_bcc,          only: interpolate_mag_field
      use constants,        only: pdims, xdim, zdim, ORTHO1, ORTHO2, LO, HI, psi_n, uh_n, magh_n, psih_n, INVALID, &
           &                      rk_coef, psidim, cs_i2_n, first_stage
      use fluidindex,       only: flind, iarr_all_dn, iarr_all_mx, iarr_all_swp, iarr_mag_swp
      use fluxtypes,        only: ext_fluxes
      use global,           only: dt, cc_mag, integration_order
      use grid_cont,        only: grid_container
      use named_array_list, only: wna, qna
      use sources,          only: internal_sources, care_for_positives

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer(kind=4),               intent(in) :: ddim
      integer,                       intent(in) :: istep     ! stage in the time integration scheme

      real, dimension(cg%n_(ddim), size(cg%u,1)) :: u
      real, dimension(cg%n_(ddim), xdim:zdim)    :: b
      real, dimension(cg%n_(ddim)), target       :: psi
      real, dimension(:,:), pointer              :: pu, pu0, pb, pb0
      real, dimension(:),   pointer              :: cs2
      integer                                    :: i1, i2
      real, dimension(:), pointer                :: ppsi, ppsi0
      integer(kind=4)                            :: psii, uhi, bhi, psihi
      real, dimension(size(u,1),size(u,2))       :: u0, u1
      real, dimension(size(b,1),size(b,2)+1)     :: b0, b1  ! Bx, By, Bz, psi
      real, dimension(size(u,1), flind%fluids), target :: vx
      type(ext_fluxes)                           :: eflx
      integer                                    :: i_cs_iso2

      uhi = wna%ind(uh_n)
      bhi = wna%ind(magh_n)

      psii  = INVALID
      psihi = INVALID
      if (qna%exists(psi_n)) then
         psii = qna%ind(psi_n)
         psihi = qna%ind(psih_n)
      endif
      psi = 0.

      ppsi => psi ! suppress compiler complains on possibly uninitialized pointer

      if (qna%exists(cs_i2_n)) then
         i_cs_iso2 = qna%ind(cs_i2_n)
      else
         i_cs_iso2 = -1
      endif
      cs2 => null()

      do i2 = cg%ijkse(pdims(ddim, ORTHO2), LO), cg%ijkse(pdims(ddim, ORTHO2), HI)
         do i1 = cg%ijkse(pdims(ddim, ORTHO1), LO), cg%ijkse(pdims(ddim, ORTHO1), HI)

            ! transposition for compatibility with RTVD-based routines
            pu0 => cg%w(uhi)%get_sweep(ddim,i1,i2)
            pu => cg%w(wna%fi)%get_sweep(ddim,i1,i2)
            if (istep == first_stage(integration_order)) pu0 = pu
            ! such copy is a bit faster than whole copy of u and we don't have to modify all the source routines

            u0(:, iarr_all_swp(ddim,:)) = transpose(pu0(:,:))
            u(:, iarr_all_swp(ddim,:)) = transpose(pu(:,:))

            pb => cg%w(wna%bi)%get_sweep(ddim,i1,i2)
            pb0 => cg%w(bhi)%get_sweep(ddim,i1,i2)
            if (istep == first_stage(integration_order)) pb0 = pb

            if (cc_mag) then
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

            if (i_cs_iso2 > 0) cs2 => cg%q(i_cs_iso2)%get_sweep(ddim,i1,i2)

            call cg%set_fluxpointers(ddim, i1, i2, eflx)
            u1 = u
            b1(:, xdim:zdim) = b
            vx = u(:, iarr_all_mx) / u(:, iarr_all_dn) ! this may also be useful for gravitational acceleration
            if (psii > INVALID) then
               ppsi0 => cg%q(psihi)%get_sweep(ddim,i1,i2)
               ppsi => cg%q(psii)%get_sweep(ddim,i1,i2)
               if (istep == first_stage(integration_order)) ppsi0 = ppsi

               b0(:, psidim) = ppsi0(:)
               b1(:, psidim) = ppsi(:)

               call solve(u0, b0, u1, b1, cs2, rk_coef(istep) * dt/cg%dl(ddim), eflx)

            else
               call solve(u0, b0(:, xdim:zdim), u1, b1(:, xdim:zdim), cs2, rk_coef(istep) * dt/cg%dl(ddim), eflx)
            endif

            call internal_sources(size(u, 1, kind=4), u, u1, b, cg, istep, ddim, i1, i2, rk_coef(istep) * dt, vx)
            ! See the results of Jeans test with RTVD and RIEMANN for estimate of accuracy.

            call care_for_positives(size(u, 1, kind=4), u1, b1, cg, ddim, i1, i2)

            call cg%save_outfluxes(ddim, i1, i2, eflx)
            pu(:,:) = transpose(u1(:, iarr_all_swp(ddim,:)))
            if (cc_mag) pb(:,:) = transpose(b1(:, iarr_mag_swp(ddim,:))) ! ToDo figure out how to manage CT energy fixup without extra storage
            if (psii /= INVALID) ppsi = b1(:, psidim)
         enddo
      enddo

      nullify(cs2)

   end subroutine solve_cg_ub

!! \warning:
!! * exceptionally slow when blocks are 120^3 (sometimes an order of magnitude)
!! * kinda slow when blocks are 20^2 and 504^2
!! * significantly slower on 2520^2 when compared to 1260^2
!! * optimal vector length could be around 180 cells (from 90 to 630 are also good)

   subroutine solve_cg_u(cg, ddim, istep)

      use constants,        only: pdims, ORTHO1, ORTHO2, LO, HI, uh_n, rk_coef, cs_i2_n, first_stage
      use fluidindex,       only: flind, iarr_all_dn, iarr_all_mx, iarr_all_swp
      use fluxtypes,        only: ext_fluxes
      use global,           only: dt, integration_order
      use grid_cont,        only: grid_container
      use named_array_list, only: wna, qna
      use sources,          only: internal_sources, care_for_positives

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer(kind=4),               intent(in) :: ddim
      integer,                       intent(in) :: istep     ! stage in the time integration scheme

      real, dimension(cg%n_(ddim), size(cg%u,1)) :: u
      real, dimension(:,:), pointer              :: pu, pu0
      real, dimension(:),   pointer              :: cs2
      integer                                    :: i1, i2
      integer(kind=4)                            :: uhi
      real, dimension(size(u,1),size(u,2))       :: u0, u1
      real, dimension(size(u,1), flind%fluids), target :: vx
      type(ext_fluxes)                           :: eflx
      real, dimension(1, 1) :: b ! ugly
      integer                                    :: i_cs_iso2

      b = 0.
      uhi = wna%ind(uh_n)
      if (qna%exists(cs_i2_n)) then
         i_cs_iso2 = qna%ind(cs_i2_n)
      else
         i_cs_iso2 = -1
      endif
      cs2 => null()

      do i2 = cg%ijkse(pdims(ddim, ORTHO2), LO), cg%ijkse(pdims(ddim, ORTHO2), HI)
         do i1 = cg%ijkse(pdims(ddim, ORTHO1), LO), cg%ijkse(pdims(ddim, ORTHO1), HI)

            ! transposition for compatibility with RTVD-based routines
            pu0 => cg%w(uhi)%get_sweep(ddim,i1,i2)
            pu => cg%w(wna%fi)%get_sweep(ddim,i1,i2)
            if (istep == first_stage(integration_order)) pu0 = pu

            u0(:, iarr_all_swp(ddim,:)) = transpose(pu0(:,:))
            u(:, iarr_all_swp(ddim,:)) = transpose(pu(:,:))
            if (i_cs_iso2 > 0) cs2 => cg%q(i_cs_iso2)%get_sweep(ddim,i1,i2)

            call cg%set_fluxpointers(ddim, i1, i2, eflx)
            u1 = u
            vx = u(:, iarr_all_mx) / u(:, iarr_all_dn) ! this may also be useful for gravitational acceleration

            call solve_u(u0, u1, cs2, rk_coef(istep) * dt/cg%dl(ddim), eflx)

            call internal_sources(size(u, 1, kind=4), u, u1, b, cg, istep, ddim, i1, i2, rk_coef(istep) * dt, vx)
            ! See the results of Jeans test with RTVD and RIEMANN for estimate of accuracy.

            call care_for_positives(size(u, 1, kind=4), u1, b, cg, ddim, i1, i2)

            call cg%save_outfluxes(ddim, i1, i2, eflx)
            pu(:,:) = transpose(u1(:, iarr_all_swp(ddim,:)))
         enddo
      enddo

      nullify(cs2)

   end subroutine solve_cg_u

!>
!! \brief Make an Euler step of length dtodx from state [u0, b0] to [u1, b1]
!!
!! This is a basic block that can be used for higher order Runge-Kutta schemes too.
!!
!! k-th interface is between k-th cell and (k+1)-th cell
!! We don't calculate n-th interface because it is as incomplete as 0-th interface
!<

   subroutine solve(u0, b0, u1, b1, cs2, dtodx, eflx)

      use constants,      only: DIVB_HDC, xdim, ydim, zdim
      use fluxtypes,      only: ext_fluxes
      use global,         only: divB_0_method
      use hlld,           only: riemann_wrap
      use interpolations, only: interpol

      implicit none

      real, dimension(:,:),        intent(in)    :: u0     !< cell-centered initial fluid states
      real, dimension(:,:),        intent(in)    :: b0     !< cell-centered initial magnetic field states (including psi field when necessary)
      real, dimension(:,:),        intent(inout) :: u1     !< cell-centered intermediate fluid states
      real, dimension(:,:),        intent(inout) :: b1     !< cell-centered intermediate magnetic field states (including psi field when necessary)
      real, dimension(:), pointer, intent(in)    :: cs2    !< square of local isothermal sound speed
      real,                        intent(in)    :: dtodx  !< timestep advance: RK-factor * timestep / cell length
      type(ext_fluxes),            intent(inout) :: eflx   !< external fluxes

      ! left and right states at interfaces 1 .. n-1
      real, dimension(size(u0, 1)-1, size(u0, 2)), target :: ql, qr
      real, dimension(size(b0, 1)-1, size(b0, 2)), target :: bl, br

      ! fluxes through interfaces 1 .. n-1
      real, dimension(size(u0, 1)-1, size(u0, 2)), target :: flx
      real, dimension(size(b0, 1)-1, size(b0, 2)), target :: mag_flx

      ! updates required for higher order of integration will likely have shorter length

      integer, parameter :: in = 1  ! index for cells

      mag_flx = huge(1.)

      call interpol(u1, ql, qr, b1, bl, br)
      call riemann_wrap(ql, qr, bl, br, cs2, flx, mag_flx) ! Now we advance the left and right states by a timestep.

      if (associated(eflx%li)) flx(eflx%li%index, :) = eflx%li%uflx
      if (associated(eflx%ri)) flx(eflx%ri%index, :) = eflx%ri%uflx
      if (associated(eflx%lo)) eflx%lo%uflx = flx(eflx%lo%index, :)
      if (associated(eflx%ro)) eflx%ro%uflx = flx(eflx%ro%index, :)

      if (divB_0_method == DIVB_HDC) then
         if (associated(eflx%li)) mag_flx(eflx%li%index, :) = eflx%li%bflx
         if (associated(eflx%ri)) mag_flx(eflx%ri%index, :) = eflx%ri%bflx
         if (associated(eflx%lo)) eflx%lo%bflx = mag_flx(eflx%lo%index, :)
         if (associated(eflx%ro)) eflx%ro%bflx = mag_flx(eflx%ro%index, :)
      endif

      associate (nx => size(u0, in))
         u1(2:nx-1, :) = u0(2:nx-1, :) + dtodx * (flx(:nx-2, :) - flx(2:, :))
         if (divB_0_method == DIVB_HDC) then
            b1(2:nx-1, :) = b0(2:nx-1, :) + dtodx * (mag_flx(:nx-2, :) - mag_flx(2:, :))
         else
            b1(2:nx-1, xdim) = b0(2:nx-1, xdim)
            b1(2:nx-1, ydim:zdim) = b0(2:nx-1, ydim:zdim) + dtodx * (mag_flx(:nx-2, ydim:zdim) - mag_flx(2:, ydim:zdim))
            ! no psidim for CT
         endif
      end associate

   end subroutine solve

   subroutine solve_u(u0, u1, cs2, dtodx, eflx)

      use fluxtypes,      only: ext_fluxes
      use hlld,           only: riemann_wrap_u
      use interpolations, only: interpol

      implicit none

      real, dimension(:,:),        intent(in)    :: u0     !< cell-centered initial fluid states
      real, dimension(:,:),        intent(inout) :: u1     !< cell-centered intermediate fluid states
      real, dimension(:), pointer, intent(in)    :: cs2    !< square of local isothermal sound speed
      real,                        intent(in)    :: dtodx  !< timestep advance: RK-factor * timestep / cell length
      type(ext_fluxes),            intent(inout) :: eflx   !< external fluxes

      ! left and right states at interfaces 1 .. n-1
      real, dimension(size(u0, 1)-1, size(u0, 2)), target :: ql, qr

      ! fluxes through interfaces 1 .. n-1
      real, dimension(size(u0, 1)-1, size(u0, 2)), target :: flx

      ! updates required for higher order of integration will likely have shorter length

      integer, parameter :: in = 1  ! index for cells

      call interpol(u1, ql, qr)
      call riemann_wrap_u(ql, qr, cs2, flx) ! Now we advance the left and right states by a timestep.

      if (associated(eflx%li)) flx(eflx%li%index, :) = eflx%li%uflx
      if (associated(eflx%ri)) flx(eflx%ri%index, :) = eflx%ri%uflx
      if (associated(eflx%lo)) eflx%lo%uflx = flx(eflx%lo%index, :)
      if (associated(eflx%ro)) eflx%ro%uflx = flx(eflx%ro%index, :)

      associate (nx => size(u0, in))
         u1(2:nx-1, :) = u0(2:nx-1, :) + dtodx * (flx(:nx-2, :) - flx(2:, :))
      end associate

   end subroutine solve_u

end module solvecg_riemann
