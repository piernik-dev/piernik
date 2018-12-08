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
!! \brief Module that implements a single sweep
!!
!! \details This one is the most difficult to aggregate messages carrying flux data on fine/coarse interfaces as
!! there are complicated dependencies between grids. It is possible to calculate which fine grids should be
!! computed first in ourder to make critical fluxes available as early as possible.
!!
!! OPT: some fluxes can be copied locally without invoking MPI
!<

module sweeps

! pulled by RTVD

   implicit none

   private
   public  :: sweep

contains
!------------------------------------------------------------------------------------------
   function interpolate_mag_field(cdim, cg, i1, i2) result (b)

      use constants,        only: pdims, xdim, ydim, zdim, half, ORTHO1, ORTHO2
      use domain,           only: dom
      use fluidindex,       only: iarr_mag_swp, nmag
      use grid_cont,        only: grid_container
      use named_array_list, only: wna

      implicit none

      integer(kind=4),               intent(in)    :: cdim
      type(grid_container), pointer, intent(inout) :: cg
      integer,                       intent(in)    :: i1, i2
      real, dimension(cg%n_(cdim), nmag)           :: b

      real, dimension(:), pointer                  :: pb, pb1
      integer(kind=4)                              :: ibx, iby, ibz
      integer                                      :: i1p, i2p

      !> OPTIMIZE ME

      ibx = iarr_mag_swp(cdim,xdim)
      iby = iarr_mag_swp(cdim,ydim)
      ibz = iarr_mag_swp(cdim,zdim)

      i1p = i1+dom%D_(pdims(cdim, ORTHO1))
      i2p = i2+dom%D_(pdims(cdim, ORTHO2))

      pb => cg%w(wna%bi)%get_sweep(cdim,ibx,i1,i2)
      b(1:cg%n_(cdim)-1, ibx) = half*( pb(1:cg%n_(cdim)-1)+pb(2:cg%n_(cdim)) )
      b(cg%n_(cdim),     ibx) = b(cg%n_(cdim)-1, ibx)

      pb  => cg%w(wna%bi)%get_sweep(cdim,iby,i1,i2)
      if (cdim == xdim) then
         pb1 => cg%w(wna%bi)%get_sweep(cdim,iby,i1p,i2)
      else
         pb1 => cg%w(wna%bi)%get_sweep(cdim,iby,i1,i2p)
      endif
      b(:, iby) = half*(pb + pb1)

      pb  => cg%w(wna%bi)%get_sweep(cdim,ibz,i1,i2)
      if (cdim == xdim) then
         pb1 => cg%w(wna%bi)%get_sweep(cdim,ibz,i1,i2p)
      else
         pb1 => cg%w(wna%bi)%get_sweep(cdim,ibz,i1p,i2)
      endif
      b(:, ibz) = half*(pb + pb1)

      b(:, iarr_mag_swp(cdim,:)) = b(:,:)
      nullify(pb,pb1)

   end function interpolate_mag_field
!------------------------------------------------------------------------------------------
   !>
   !! TODO: comment me and change name if necessary
   !<
   integer function compute_nr_recv(cdim) result(nr)

      use constants,        only: LO, HI, I_ONE
      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use mpi,              only: MPI_DOUBLE_PRECISION
      use mpisetup,         only: comm, mpi_err, req, inflate_req

      implicit none
      integer(kind=4), intent(in)       :: cdim

      type(cg_list_element), pointer    :: cgl
      integer(kind=8), dimension(LO:HI) :: jc
      integer :: g

      nr = 0
      cgl => leaves%first
      do while (associated(cgl))
         cgl%cg%processed = .false.
         cgl%cg%finebnd(cdim, LO)%uflx(:, :, :) = 0. !> \warning overkill
         cgl%cg%finebnd(cdim, HI)%uflx(:, :, :) = 0.
         if (allocated(cgl%cg%rif_tgt%seg)) then
            associate ( seg => cgl%cg%rif_tgt%seg )
               do g = lbound(seg, dim=1), ubound(seg, dim=1)
                  jc = seg(g)%se(cdim, :)
                  if (jc(LO) == jc(HI)) then
                     nr = nr + I_ONE
                     if (nr > size(req, dim=1)) call inflate_req
                     call MPI_Irecv(seg(g)%buf, size(seg(g)%buf(:, :, :)), MPI_DOUBLE_PRECISION, seg(g)%proc, seg(g)%tag, comm, req(nr), mpi_err)
                     seg(g)%req => req(nr)
                  endif
               enddo
            end associate
         endif
         cgl => cgl%nxt
      enddo
   end function compute_nr_recv
!------------------------------------------------------------------------------------------
   !>
   !! TODO: comment me and change name if necessary
   !<
   subroutine recv_cg_finebnd(cdim, cg, all_received)
      use constants,        only: LO, HI, INVALID, ORTHO1, ORTHO2, pdims
      use dataio_pub,       only: die
      use grid_cont,        only: grid_container
      use mpi,              only: MPI_STATUS_IGNORE
      use mpisetup,         only: mpi_err
      implicit none
      integer(kind=4), intent(in)                  :: cdim
      type(grid_container), pointer, intent(inout) :: cg
      logical, intent(out)                         :: all_received

      integer :: g, lh
      logical :: received
      integer(kind=8), dimension(LO:HI) :: j1, j2, jc

      all_received = .true.
      if (allocated(cg%rif_tgt%seg)) then
         associate ( seg => cg%rif_tgt%seg )
         do g = lbound(seg, dim=1), ubound(seg, dim=1)
            jc = seg(g)%se(cdim, :)
            if (jc(LO) == jc(HI)) then
               call MPI_Test(seg(g)%req, received, MPI_STATUS_IGNORE, mpi_err)
               if (received) then
                  jc = seg(g)%se(cdim, :) !> \warning: partially duplicated code (see below)
                  j1 = seg(g)%se(pdims(cdim, ORTHO1), :)
                  j2 = seg(g)%se(pdims(cdim, ORTHO2), :)
                  if (jc(LO) /= jc(HI)) call die("[sweeps:sweep] layer too thick (Recv)")
                  if (all(cg%finebnd(cdim, LO)%index(j1(LO):j1(HI), j2(LO):j2(HI)) == jc(LO))) then
                     lh = LO
                  else if (all(cg%finebnd(cdim, HI)%index(j1(LO):j1(HI), j2(LO):j2(HI)) == jc(LO))) then
                     lh = HI
                  else
                     call die("[sweeps:sweep] Cannot determine side (Recv)")
                     lh = INVALID
                  endif
                  ! cg%finebnd(cdim, lh)%uflx(:, j1(LO):j1(HI), j2(LO):j2(HI)) = &
                  !     cg%finebnd(cdim, lh)%uflx(:, j1(LO):j1(HI), j2(LO):j2(HI)) + seg(g)%buf(:, :, :)
                  ! for more general decompositions with odd-offset patches it might be necessary to do sum, but it need to be debugged first
                  cg%finebnd(cdim, lh)%uflx(:, j1(LO):j1(HI), j2(LO):j2(HI)) = seg(g)%buf(:, :, :)
               else
                  all_received = .false.
               endif
            endif
         enddo
         end associate
      endif
      return
   end subroutine recv_cg_finebnd
!------------------------------------------------------------------------------------------
   !>
   !! TODO: comment me and change name if necessary
   !<
   subroutine send_cg_coarsebnd(cdim, cg, nr)
      use constants,        only: pdims, LO, HI, ORTHO1, ORTHO2, I_ONE, INVALID
      use dataio_pub,       only: die
      use domain,           only: dom
      use grid_cont,        only: grid_container
      use grid_helpers,     only: f2c_o
      use mpi,              only: MPI_DOUBLE_PRECISION
      use mpisetup,         only: comm, mpi_err, req, inflate_req

      implicit none
      integer(kind=4), intent(in)                  :: cdim
      type(grid_container), pointer, intent(inout) :: cg
      integer, intent(inout)                       :: nr

      integer :: g, lh
      integer(kind=8), dimension(LO:HI) :: j1, j2, jc
      integer(kind=8) :: j, k


      if (allocated(cg%rof_tgt%seg)) then
         associate ( seg => cg%rof_tgt%seg )
         do g = lbound(seg, dim=1), ubound(seg, dim=1)
            jc = seg(g)%se(cdim, :) !> \warning: partially duplicated code (see above)
            if (jc(LO) == jc(HI)) then
               j1 = seg(g)%se(pdims(cdim, ORTHO1), :)
               j2 = seg(g)%se(pdims(cdim, ORTHO2), :)
               if (all(cg%coarsebnd(cdim, LO)%index(j1(LO):j1(HI), j2(LO):j2(HI)) == jc(LO))) then
                  lh = LO
               else if (all(cg%coarsebnd(cdim, HI)%index(j1(LO):j1(HI), j2(LO):j2(HI)) == jc(LO))) then
                  lh = HI
               else
                  call die("[sweeps:sweep] Cannot determine side (Send)")
                  lh = INVALID
               endif

               seg(g)%buf(:, :, :) = 0.
               do j = j1(LO), j1(HI)
                  do k = j2(LO), j2(HI)
                     seg(g)%buf(:, f2c_o(j), f2c_o(k)) = seg(g)%buf(:, f2c_o(j), f2c_o(k)) + cg%coarsebnd(cdim, lh)%uflx(:, j, k)
                  enddo
               enddo
               seg(g)%buf = 1/2.**(dom%eff_dim-1) * seg(g)%buf

               nr = nr + I_ONE
               if (nr > size(req, dim=1)) call inflate_req
               call MPI_Isend(seg(g)%buf, size(seg(g)%buf(:, :, :)), MPI_DOUBLE_PRECISION, seg(g)%proc, seg(g)%tag, comm, req(nr), mpi_err)
               seg(g)%req => req(nr)
            endif
         enddo
         end associate
      endif
   end subroutine send_cg_coarsebnd
!------------------------------------------------------------------------------------------
   subroutine solve_cg(cg, cdim, istep, fargo_vel)

      use cg_level_connected, only: cg_level_connected_T, find_level
      use constants,          only: pdims, LO, HI, uh_n, cs_i2_n, ORTHO1, ORTHO2, VEL_CR, VEL_RES, ydim, rk_coef, first_stage
      use dataio_pub,         only: die
      use domain,             only: dom
      use fluidindex,         only: flind, iarr_all_swp, nmag, iarr_all_dn, iarr_all_mx
      use fluxtypes,          only: ext_fluxes
      use global,             only: dt, integration_order, use_fargo
      use grid_cont,          only: grid_container
      use gridgeometry,       only: set_geo_coeffs
      use named_array_list,   only: qna, wna
      use rtvd,               only: relaxing_tvd
      use sources,            only: prepare_sources, all_sources, care_for_positives
#ifdef MAGNETIC
      use fluidindex,         only: iarr_mag_swp
#endif /* MAGNETIC */

      implicit none

      type(grid_container), pointer, intent(inout) :: cg
      integer(kind=4),               intent(in)    :: cdim
      integer,                       intent(in)    :: istep     ! stage in the time integration scheme
      integer(kind=4), optional,     intent(in)    :: fargo_vel

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
      type(cg_level_connected_T), pointer :: curl

      uhi = wna%ind(uh_n)
      full_dim = dom%has_dir(cdim)
      if (qna%exists(cs_i2_n)) then
         i_cs_iso2 = qna%ind(cs_i2_n)
      else
         i_cs_iso2 = -1
      endif
      call eflx%init
      allocate( b(cg%n_(cdim), nmag), u(cg%n_(cdim), flind%all), u0(cg%n_(cdim), flind%all), u1(cg%n_(cdim), flind%all), vx(cg%n_(cdim), flind%fluids))
      !OPT for AMR it may be worthwhile to move it to global scope

      b(:,:) = 0.0
      u(:,:) = 0.0

      if (istep == first_stage(integration_order)) then
         call prepare_sources(cg)
         cg%w(uhi)%arr = cg%u
      endif

      !> \todo OPT: use cg%leafmap to skip lines fully covered by finer grids
      ! it should be also possible to compute only parts of lines that aren't covered by finer grids
      curl => find_level(cg%l%id)

      cs2 => null()
      do i2 = cg%ijkse(pdims(cdim, ORTHO2), LO), cg%ijkse(pdims(cdim, ORTHO2), HI)
         do i1 = cg%ijkse(pdims(cdim, ORTHO1), LO), cg%ijkse(pdims(cdim, ORTHO1), HI)

#ifdef MAGNETIC
            if (full_dim) then
               b(:,:) = interpolate_mag_field(cdim, cg, i1, i2)
            else
               pb => cg%w(wna%bi)%get_sweep(cdim, i1, i2)   ! BEWARE: is it correct for 2.5D ?
               b(:, iarr_mag_swp(cdim,:))  = transpose(pb(:,:))
            endif
#endif /* MAGNETIC */

            call set_geo_coeffs(cdim, flind, i1, i2, cg)

            pu                     => cg%w(wna%fi   )%get_sweep(cdim,i1,i2)
            pu0                    => cg%w(uhi      )%get_sweep(cdim,i1,i2)
            if (i_cs_iso2 > 0) cs2 => cg%q(i_cs_iso2)%get_sweep(cdim,i1,i2)

            u (:, iarr_all_swp(cdim,:)) = transpose(pu (:,:))
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
                  call die("[sweeps:sweep] Unknown FARGO_VEL")
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
            if (apply_sources) call all_sources(cg%n_(cdim), u, u1, b, cg, istep, cdim, i1, i2, rk_coef(istep) * dt, vx)
            ! istep is important only for balsara and selfgravity

            call care_for_positives(cg%n_(cdim), u1, b, cg, cdim, i1, i2)
            u(:,:) = u1(:,:)
            call cg%save_outfluxes(cdim, i1, i2, eflx)

            pu(:,:) = transpose(u(:, iarr_all_swp(cdim,:)))
            nullify(pu,pu0,cs2)
         enddo
      enddo

      deallocate(b, u, u0, u1, vx)

      cg%processed = .true.

   end subroutine solve_cg

!>
!! \brief Call all boundaries, try to avoid unnecessary parts.
!!
!! For some reasons dir=cdim affect mcrwind tests if sweeps_mgu
!! \todo Find out why. Is it related to position of magnetic field components?
!!
!! \todo Once it gets simplified enough merge it back to sweep.
!<

   subroutine update_boundaries(cdim, istep)

      use all_boundaries, only: all_fluid_boundaries
      use constants,      only: first_stage
      use domain,         only: dom
      use global,         only: sweeps_mgu, integration_order

      implicit none

      integer(kind=4), intent(in) :: cdim
      integer,         intent(in) :: istep

      if (dom%has_dir(cdim)) then
         if (sweeps_mgu) then
            if (istep == first_stage(integration_order)) then
               call all_fluid_boundaries(nocorners = .true., dir = cdim)
            else
               call all_fluid_boundaries(nocorners = .true.)
            endif
         else
            if (istep == first_stage(integration_order)) then
               call all_fluid_boundaries(nocorners = .true.)
            else
               call all_fluid_boundaries
            endif
         endif
      endif

   end subroutine update_boundaries
!------------------------------------------------------------------------------------------
   subroutine sweep(cdim, fargo_vel)

      use cg_leaves,          only: leaves
      use cg_list,            only: cg_list_element
      use constants,          only: ydim, first_stage, last_stage
      use dataio_pub,         only: die
      use global,             only: integration_order, use_fargo
      use grid_cont,          only: grid_container
      use mpisetup,           only: mpi_err, req, status

      implicit none

      integer(kind=4),           intent(in) :: cdim
      integer(kind=4), optional, intent(in) :: fargo_vel

      integer                        :: istep
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
      logical                        :: all_processed, all_received
      integer                        :: blocks_done
      integer                        :: g, nr, nr_recv
      integer(kind=4), dimension(:,:), pointer :: mpistatus

      if (use_fargo .and. cdim == ydim .and. .not. present(fargo_vel)) &
           call die("[sweeps:sweep] FARGO velocity keyword not present in y sweep")

      if (integration_order < lbound(first_stage, 1) .or. integration_order > ubound(first_stage, 1)) &
           call die("[sweeps:sweep] unknown integration_order")

      do istep = first_stage(integration_order), last_stage(integration_order)
         nr_recv = compute_nr_recv(cdim)
         nr = nr_recv
         all_processed = .false.

         do while (.not. all_processed)
            all_processed = .true.
            blocks_done = 0
            cgl => leaves%first
            do while (associated(cgl))
               cg => cgl%cg

               if (.not. cg%processed) then
                  call recv_cg_finebnd(cdim, cg, all_received)

                  if (all_received) then
                     call solve_cg(cg, cdim, istep, fargo_vel)
                     call send_cg_coarsebnd(cdim, cg, nr)
                     blocks_done = blocks_done + 1
                  else
                     all_processed = .false.
                  endif
               endif
               cgl => cgl%nxt
            enddo

            if (.not. all_processed .and. blocks_done == 0) then
               if (nr_recv > 0) then
                  mpistatus => status(:, :nr_recv)
                  call MPI_Waitany(nr_recv, req(:nr_recv), g, mpistatus, mpi_err)
                  ! g is the number of completed operations
               endif
            endif
         enddo

         if (nr > 0) then
            mpistatus => status(:, :nr)
            call MPI_Waitall(nr, req(:nr), mpistatus, mpi_err)
         endif

         call update_boundaries(cdim, istep)
      enddo

   end subroutine sweep

end module sweeps
