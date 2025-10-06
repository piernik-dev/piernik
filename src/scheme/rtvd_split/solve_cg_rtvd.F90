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
!! \brief Module that implements 1D hydro calls for a single cg for RTVD
!<

module solvecg_rtvd

! pulled by ANY

   implicit none

   private
   public  :: solve_cg_rtvd

contains

!>
!! \brief Apply MHD update + source terms to a single grid container, rely on properly updated guardcells, handle local fine-coarse fluxes.
!!
!! This routine has to conform to the interface defined in sweeps::sweep
!<

   subroutine solve_cg_rtvd(cg, cdim, istep, fargo_vel)

      use bfc_bcc,            only: interpolate_mag_field
      use cg_level_connected, only: cg_level_connected_t
      use constants,          only: pdims, LO, HI, uh_n, cs_i2_n, ORTHO1, ORTHO2, VEL_CR, VEL_RES, ydim, rk_coef, first_stage
      use dataio_pub,         only: die
      use domain,             only: dom
      use find_lev,           only: find_level
      use fluidindex,         only: flind, iarr_all_swp, nmag, iarr_all_dn, iarr_all_mx
      use fluxtypes,          only: ext_fluxes
      use global,             only: dt, use_fargo, integration_order
      use grid_cont,          only: grid_container
      use gridgeometry,       only: set_geo_coeffs
      use named_array_list,   only: qna, wna
      use rtvd,               only: relaxing_tvd
      use sources,            only: internal_sources, care_for_positives
#ifdef MAGNETIC
      use fluidindex,         only: iarr_mag_swp
#endif /* MAGNETIC */

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer(kind=4),               intent(in) :: cdim
      integer,                       intent(in) :: istep     ! stage in the time integration scheme
      integer(kind=4), optional,     intent(in) :: fargo_vel

      real, dimension(:,:), allocatable :: b, u, u0, u1, vx
      integer                           :: i1, i2, uhi, ifl
      logical                           :: full_dim
      type(ext_fluxes)                  :: eflx
      real, dimension(:,:),  pointer    :: pu, pu0
      integer                           :: i_cs_iso2
#ifdef MAGNETIC
      real, dimension(:,:),  pointer    :: pb
#endif /* MAGNETIC */
      real, dimension(:),    pointer    :: cs2
      logical :: apply_sources
      type(cg_level_connected_t), pointer :: curl

      uhi = wna%ind(uh_n)
      full_dim = dom%has_dir(cdim)
      if (qna%exists(cs_i2_n)) then
         i_cs_iso2 = qna%ind(cs_i2_n)
      else
         i_cs_iso2 = -1
      endif
      allocate( b(cg%n_(cdim), nmag), u(cg%n_(cdim), flind%all), u0(cg%n_(cdim), flind%all), u1(cg%n_(cdim), flind%all), vx(cg%n_(cdim), flind%fluids))
      !OPT for AMR it may be worthwhile to move it to global scope

      !> \todo OPT: use cg%leafmap to skip lines fully covered by finer grids
      ! it should be also possible to compute only parts of lines that aren't covered by finer grids
      curl => find_level(cg%l%id)

      cs2 => null()
      do i2 = cg%ijkse(pdims(cdim, ORTHO2), LO), cg%ijkse(pdims(cdim, ORTHO2), HI)
         do i1 = cg%ijkse(pdims(cdim, ORTHO1), LO), cg%ijkse(pdims(cdim, ORTHO1), HI)

#ifdef MAGNETIC
            if (full_dim) then
               b(:,:) = interpolate_mag_field(cdim, cg, i1, i2, wna%bi)
            else
               pb => cg%w(wna%bi)%get_sweep(cdim, i1, i2)   ! BEWARE: is it correct for 2.5D ?
               b(:, iarr_mag_swp(cdim,:))  = transpose(pb(:,:))
            endif
#endif /* MAGNETIC */

            call set_geo_coeffs(cdim, flind, i1, i2, cg)

            pu                     => cg%w(wna%fi   )%get_sweep(cdim, i1, i2)
            pu0                    => cg%w(uhi      )%get_sweep(cdim, i1, i2)
            if (i_cs_iso2 > 0) cs2 => cg%q(i_cs_iso2)%get_sweep(cdim, i1, i2)

            u (:, iarr_all_swp(cdim,:)) = transpose(pu (:,:))

            if (istep == first_stage(integration_order)) pu0 = pu
            ! such copy is a bit faster than whole copy of u and we don't have to modify all the source routines

            u0(:, iarr_all_swp(cdim,:)) = transpose(pu0(:,:))
            if (use_fargo .and. cdim == ydim) then
               if (fargo_vel == VEL_RES) then
                  do ifl = 1, flind%fluids
                     vx(:, ifl) = u(:, iarr_all_mx(ifl)) / u(:, iarr_all_dn(ifl)) - curl%omega_mean(i2, ifl) * cg%x(i2)
                  enddo
                  apply_sources = .true.
               elseif (fargo_vel == VEL_CR) then
                  do ifl = 1, flind%fluids
                     vx(:, ifl) = curl%omega_cr(i2, ifl) * cg%x(i2)
                  enddo
                  apply_sources = .false.
               else
                  call die("[solve_cg_rtvd:solve_cg_rtvd] Unknown FARGO_VEL")
                  apply_sources = .false.
               endif
            else
               apply_sources = .true.
               vx(:,:) = u(:,iarr_all_mx(:)) / u(:,iarr_all_dn(:))
               if (full_dim) then
                  vx(1,:) = vx(2,:)
                  vx(cg%n_(cdim),:) = vx(cg%n_(cdim)-1,:)
               endif
            endif

            call cg%set_fluxpointers(cdim, i1, i2, eflx)
            !OPT: try to avoid these explicit initializations of u1(:,:)
            u1 = u

            call relaxing_tvd(cg%n_(cdim), u0, u1, vx, b, cs2, istep, rk_coef(istep) * dt / cg%dl(cdim), eflx)
            ! RTVD needs istep only to do something in 2nd stage of RK2
! Source terms -------------------------------------
            if (apply_sources) call internal_sources(cg%n_(cdim), u, u1, b, cg, istep, cdim, i1, i2, rk_coef(istep) * dt, vx)
            ! istep is important only for balsara and selfgravity

            call care_for_positives(cg%n_(cdim), u1, b, cg, cdim, i1, i2)
            call cg%save_outfluxes(cdim, i1, i2, eflx)

            pu(:,:) = transpose(u1(:, iarr_all_swp(cdim,:)))
            nullify(pu,pu0,cs2)
         enddo
      enddo

      deallocate(b, u, u0, u1, vx)

      cg%processed = .true.

   end subroutine solve_cg_rtvd

end module solvecg_rtvd
