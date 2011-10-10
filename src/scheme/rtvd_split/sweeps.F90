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
#include "piernik.h"

module sweeps     ! split sweeps

! pulled by RTVD

   implicit none

   private
   public  :: sweep
#if defined SHEAR && defined FLUID_INTERACTIONS
   public  :: source_terms_y
#endif /* SHEAR && FLUID_INTERACTIONS */

contains

#if defined SHEAR && defined FLUID_INTERACTIONS
   subroutine source_terms_y

      use constants,       only: xdim, zdim, half, one
      use dataio_pub,      only: die
      use domain,          only: is_multicg
      use fluidboundaries, only: all_fluid_boundaries
      use fluidindex,      only: iarr_all_dn, iarr_all_mx, iarr_all_my, flind
      use global,          only: dt
      use grid,            only: all_cg
      use grid_cont,       only: grid_container
      use interactions,    only: dragc_gas_dust
      use shear,           only: omega, qshear

      implicit none

      real, dimension(:,:,:), allocatable :: vxr, v_r, rotaccr
      real, dimension(:,:), allocatable   :: epsa
      real, dimension(:,:,:), allocatable :: u1
      integer :: ind,i
      real, parameter, dimension(2) :: fac = [half, one]
      type(grid_container), pointer :: cg

      cg => all_cg%first%cg
      if (is_multicg) call die("[sweeps:source_terms_y] multiple grid pieces per procesor not implemented yet") !nontrivial u1

      allocate(vxr(size(iarr_all_my), cg%n_(xdim), cg%n_(zdim)), v_r(size(iarr_all_my), cg%n_(xdim), cg%n_(zdim)), rotaccr(size(iarr_all_my), cg%n_(xdim), cg%n_(zdim)))
      allocate(epsa(cg%n_(xdim), cg%n_(zdim)))
      allocate(u1(flind%all, cg%n_(xdim), cg%n_(zdim)))

      u1(:,:,:) = cg%u%arr(:,:,1,:)

      do i = 1,2
         where (u1(iarr_all_dn,:,:) > 0.0)
            vxr(:,:,:) = u1(iarr_all_mx,:,:) / u1(iarr_all_dn,:,:)
            v_r(:,:,:) = u1(iarr_all_my,:,:) / u1(iarr_all_dn,:,:)
         elsewhere
            vxr(:,:,:) = 0.0
            v_r(:,:,:) = 0.0
         endwhere
         epsa(:,:) =  u1(iarr_all_dn(2),:,:) / u1(iarr_all_dn(1),:,:)

         do ind=1,size(iarr_all_my)
            if (ind == 1) then
               rotaccr(ind,:,:) = - dragc_gas_dust * epsa(:,:) * (v_r(1,:,:) - v_r(2,:,:))
            else
               rotaccr(ind,:,:) = - dragc_gas_dust * 1.0    * (v_r(2,:,:) - v_r(1,:,:))
            endif
!           rotaccr(ind,:,:) = rotaccr(ind,:,:) - 2.0*omega*vxr(ind,:,:)
            rotaccr(ind,:,:) = rotaccr(ind,:,:) + (qshear-2.0)*omega*vxr(ind,:,:)
         enddo

         where (u1(iarr_all_dn,:,:) > 0.0)
            u1(iarr_all_my,:,:) = u1(iarr_all_my,:,:) + fac(i)*dt*rotaccr(:,:,:)*cg%u%arr(iarr_all_dn,:,1,:)
         endwhere
         cg%u%arr(iarr_all_my,:,1,:) = u1(iarr_all_my,:,:)
      enddo

      deallocate(vxr, v_r, rotaccr)
      deallocate(epsa)
      deallocate(u1)

      call all_fluid_boundaries

   end subroutine source_terms_y
#endif /* SHEAR && FLUID_INTERACTIONS */
!------------------------------------------------------------------------------------------
   function interpolate_mag_field(cdim, cg, i1, i2) result (b)

      use constants,       only: pdims, xdim, ydim, zdim, half
      use fluidindex,      only: iarr_mag_swp, nmag
      use domain,          only: D_
      use grid_cont,       only: grid_container

      implicit none

      integer(kind=4), intent(in)                  :: cdim
      type(grid_container), pointer, intent(inout) :: cg
      integer, intent(in)                          :: i1, i2
      real, dimension(nmag, cg%n_(cdim))           :: b

      real, dimension(:), pointer :: pb, pb1
      integer(kind=4)                         :: ibx, iby, ibz
      integer                                 :: i1p, i2p

      !> OPTIMIZE ME

      ibx = iarr_mag_swp(cdim,xdim)
      iby = iarr_mag_swp(cdim,ydim)
      ibz = iarr_mag_swp(cdim,zdim)

      i1p = i1+D_(pdims(cdim,ydim))
      i2p = i2+D_(pdims(cdim,zdim))

      pb => cg%b%get_sweep(cdim,ibx,i1,i2)
      b(ibx,1:cg%n_(cdim)-1) = half*( pb(1:cg%n_(cdim)-1)+pb(2:cg%n_(cdim)) )
      b(ibx,  cg%n_(cdim)  ) = b(ibx,  cg%n_(cdim)-1)

      pb  => cg%b%get_sweep(cdim,iby,i1,i2)
      if (cdim == xdim) then
         pb1 => cg%b%get_sweep(cdim,iby,i1p,i2)
      else
         pb1 => cg%b%get_sweep(cdim,iby,i1,i2p)
      endif
      b(iby,:) = half*(pb + pb1)

      pb  => cg%b%get_sweep(cdim,ibz,i1,i2)
      if (cdim == xdim) then
         pb1 => cg%b%get_sweep(cdim,ibz,i1,i2p)
      else
         pb1 => cg%b%get_sweep(cdim,ibz,i1p,i2)
      endif
      b(ibz,:) = half*(pb + pb1)

      b( iarr_mag_swp(cdim,:),:) = b(:,:)
      nullify(pb,pb1)

   end function interpolate_mag_field
!------------------------------------------------------------------------------------------
   subroutine sweep(cdim)

      use constants,       only: pdims, LO, HI, ydim, zdim
      use fluidboundaries, only: all_fluid_boundaries
      use fluidindex,      only: flind, iarr_all_swp, nmag
      use gc_list,         only: cg_list_element
      use global,          only: dt, integration_order
      use grid,            only: all_cg
      use grid_cont,       only: grid_container
      use gridgeometry,    only: set_geo_coeffs
      use rtvd,            only: relaxing_tvd
#ifdef COSM_RAYS
      use crhelpers,       only: div_v, set_div_v1d
#endif /* COSM_RAYS */

      implicit none

      integer(kind=4), intent(in)                  :: cdim

      real, dimension(:,:), allocatable :: b
      real, dimension(:,:), allocatable :: u, u0
      real, dimension(:,:), pointer     :: pu, pu0
      real, dimension(:), pointer       :: div_v1d => null(), cs2
      integer                           :: i1, i2, uhi
      integer                           :: istep
      integer                           :: i_cs_iso2
      type(cg_list_element), pointer    :: cgl
      type(grid_container), pointer     :: cg

      do istep = 1, integration_order
         cgl => all_cg%first
         do while (associated(cgl))
            cg => cgl%cg

            uhi = cg%get_na_ind_4d("uh")

            if (allocated(b)) deallocate(b)
            if (allocated(u)) deallocate(u)
            if (allocated(u0)) deallocate(u0)

            allocate(b(nmag, cg%n_(cdim)), u(flind%all, cg%n_(cdim)), u0(flind%all, cg%n_(cdim)))

            b(:,:) = 0.0
            u(:,:) = 0.0

            if (istep == 1) then
#ifdef COSM_RAYS
               call div_v(flind%ion%pos, cg)
#endif /* COSM_RAYS */
               cg%w(uhi)%arr = cg%u%arr
            endif

            cs2 => null()
            if (cg%exists("cs_iso2")) then
               i_cs_iso2 = cg%get_na_ind("cs_iso2") ! BEWARE: magic strings across multiple files
            else
               i_cs_iso2 = -1
            endif

            do i2 = cg%ijkse(pdims(cdim,zdim),LO), cg%ijkse(pdims(cdim,zdim),HI)
               do i1 = cg%ijkse(pdims(cdim,ydim),LO), cg%ijkse(pdims(cdim,ydim),HI)

#ifdef MAGNETIC
                  b = interpolate_mag_field(cdim, cg, i1, i2)
#endif /* MAGNETIC */

                  call set_geo_coeffs(cdim, flind, i1, i2, cg)
#ifdef COSM_RAYS
                  call set_div_v1d(div_v1d, cdim, i1, i2, cg)
#endif /* COSM_RAYS */

                  pu => cg%u%get_sweep(cdim,i1,i2)
                  pu0 => cg%w(uhi)%get_sweep(cdim,i1,i2)
                  if (i_cs_iso2 > 0) cs2 => cg%q(i_cs_iso2)%get_sweep(cdim,i1,i2)

                  u (iarr_all_swp(cdim,:),:) = pu(:,:)
                  u0(iarr_all_swp(cdim,:),:) = pu0(:,:)

                  call relaxing_tvd(cg%n_(cdim), u, u0, b, div_v1d, cs2, istep, cdim, i1, i2, cg%dl(cdim), dt, cg)
                  pu(:,:) = u(iarr_all_swp(cdim,:),:)
                  nullify(pu,pu0,cs2)
               enddo
            enddo

            cgl => cgl%nxt
         enddo

         call all_fluid_boundaries    ! \todo : call only x for istep=1, call all for istep=2
      enddo

      if (allocated(b)) deallocate(b)
      if (allocated(u)) deallocate(u)
      if (allocated(u0)) deallocate(u0)

   end subroutine sweep
!------------------------------------------------------------------------------------------
end module sweeps
