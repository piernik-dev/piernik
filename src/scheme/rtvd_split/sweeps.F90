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
      real, dimension(nmag, cg%n_(cdim))           :: b

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
      b(ibx,1:cg%n_(cdim)-1) = half*( pb(1:cg%n_(cdim)-1)+pb(2:cg%n_(cdim)) )
      b(ibx,  cg%n_(cdim)  ) = b(ibx,  cg%n_(cdim)-1)

      pb  => cg%w(wna%bi)%get_sweep(cdim,iby,i1,i2)
      if (cdim == xdim) then
         pb1 => cg%w(wna%bi)%get_sweep(cdim,iby,i1p,i2)
      else
         pb1 => cg%w(wna%bi)%get_sweep(cdim,iby,i1,i2p)
      endif
      b(iby,:) = half*(pb + pb1)

      pb  => cg%w(wna%bi)%get_sweep(cdim,ibz,i1,i2)
      if (cdim == xdim) then
         pb1 => cg%w(wna%bi)%get_sweep(cdim,ibz,i1,i2p)
      else
         pb1 => cg%w(wna%bi)%get_sweep(cdim,ibz,i1p,i2)
      endif
      b(ibz,:) = half*(pb + pb1)

      b( iarr_mag_swp(cdim,:),:) = b(:,:)
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
      use func,             only: f2c_o
      use grid_cont,        only: grid_container
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
   subroutine sweep(cdim)

      use all_boundaries,   only: all_fluid_boundaries
      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: pdims, LO, HI, uh_n, cs_i2_n, ORTHO1, ORTHO2
      use domain,           only: dom
      use fluidindex,       only: flind, iarr_all_swp, nmag
      use fluxtypes,        only: ext_fluxes
      use global,           only: dt, integration_order
      use grid_cont,        only: grid_container
      use gridgeometry,     only: set_geo_coeffs
      use mpisetup,         only: mpi_err, req, status
      use named_array_list, only: qna, wna
      use rtvd,             only: relaxing_tvd
#ifdef COSM_RAYS
      use crhelpers,        only: div_v, set_div_v1d
#endif /* COSM_RAYS */
#ifdef MAGNETIC
      use fluidindex,       only: iarr_mag_swp
#endif /* MAGNETIC */

      implicit none

      integer(kind=4), intent(in)       :: cdim

      integer                           :: i1, i2, uhi
      integer                           :: istep
      integer                           :: i_cs_iso2
      logical                           :: full_dim
      real, dimension(:,:), allocatable :: b
      real, dimension(:,:), allocatable :: u, u0
      real, dimension(:,:),  pointer    :: pu, pu0
#ifdef MAGNETIC
      real, dimension(:,:),  pointer    :: pb
#endif /* MAGNETIC */
      real, dimension(:),    pointer    :: div_v1d => null(), cs2
      type(cg_list_element), pointer    :: cgl
      type(grid_container),  pointer    :: cg
      type(ext_fluxes)                  :: eflx
      logical                           :: all_processed, all_received
      integer                           :: blocks_done
      integer                           :: g, nr, nr_recv
      integer(kind=4), dimension(:, :), pointer :: mpistatus
      integer :: cn_

      cn_ = 0
      full_dim = dom%has_dir(cdim)
      uhi = wna%ind(uh_n)
      if (qna%exists(cs_i2_n)) then
         i_cs_iso2 = qna%ind(cs_i2_n)
      else
         i_cs_iso2 = -1
      endif

      call eflx%init

      nr = 0
      do istep = 1, integration_order
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

                     if (cn_ /= cg%n_(cdim)) then
                        if (allocated(b))  deallocate(b)
                        if (allocated(u))  deallocate(u)
                        if (allocated(u0)) deallocate(u0)
                     endif
                     if (.not. allocated(u)) allocate(b(nmag, cg%n_(cdim)), u(flind%all, cg%n_(cdim)), u0(flind%all, cg%n_(cdim)))
                     cn_ = cg%n_(cdim)

                     b(:,:) = 0.0
                     u(:,:) = 0.0

                     if (istep == 1) then
#ifdef COSM_RAYS
                        call div_v(flind%ion%pos, cg)
#endif /* COSM_RAYS */
                        cg%w(uhi)%arr = cg%u
                     endif

                     !> \todo OPT: use cg%leafmap to skip lines fully covered by finer grids
                     ! it should be also possible to compute only parts of lines that aren't covered by finer grids
                     cs2 => null()
                     do i2 = cg%ijkse(pdims(cdim, ORTHO2), LO), cg%ijkse(pdims(cdim, ORTHO2), HI)
                        do i1 = cg%ijkse(pdims(cdim, ORTHO1), LO), cg%ijkse(pdims(cdim, ORTHO1), HI)

#ifdef MAGNETIC
                           if (full_dim) then
                              b(:,:) = interpolate_mag_field(cdim, cg, i1, i2)
                           else
                              pb => cg%w(wna%bi)%get_sweep(cdim, i1, i2)   ! BEWARE: is it correct for 2.5D ?
                              b(iarr_mag_swp(cdim,:),:)  = pb(:,:)
                           endif
#endif /* MAGNETIC */

                           call set_geo_coeffs(cdim, flind, i1, i2, cg)
#ifdef COSM_RAYS
                           call set_div_v1d(div_v1d, cdim, i1, i2, cg)
#endif /* COSM_RAYS */

                           pu                     => cg%w(wna%fi   )%get_sweep(cdim,i1,i2)
                           pu0                    => cg%w(uhi      )%get_sweep(cdim,i1,i2)
                           if (i_cs_iso2 > 0) cs2 => cg%q(i_cs_iso2)%get_sweep(cdim,i1,i2)

                           u (iarr_all_swp(cdim,:),:) = pu (:,:)
                           u0(iarr_all_swp(cdim,:),:) = pu0(:,:)

                           call cg%set_fluxpointers(cdim, i1, i2, eflx)
                           call relaxing_tvd(cg%n_(cdim), u, u0, b, div_v1d, cs2, istep, cdim, i1, i2, cg%dl(cdim), dt, cg, eflx)
                           call cg%save_outfluxes(cdim, i1, i2, eflx)

                           pu(:,:) = u(iarr_all_swp(cdim,:),:)
                           nullify(pu,pu0,cs2)
                        enddo
                     enddo

                     call send_cg_coarsebnd(cdim, cg, nr)

                     deallocate(b, u, u0)

                     cg%processed = .true.
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
               endif
            endif
         enddo

         if (nr > 0) then
            mpistatus => status(:, :nr)
            call MPI_Waitall(nr, req(:nr), mpistatus, mpi_err)
         endif

         if (full_dim) then
            if (istep == 1) then
               call all_fluid_boundaries(nocorners = .true.) !dir = cdim)
               ! For some weird reasons dir=cdim here affect mcrwind tests. \todo Find out why.
            else
               call all_fluid_boundaries!(nocorners = .true.)
               ! For some weird reasons nocorners here affect mcrwind tests. \todo Find out why.
            end if
         end if
      enddo

      if (allocated(b))  deallocate(b)
      if (allocated(u))  deallocate(u)
      if (allocated(u0)) deallocate(u0)

   end subroutine sweep

end module sweeps
