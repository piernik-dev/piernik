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

module internal_bnd
! pulled by ANY

   implicit none

   private
   public :: arr3d_boundaries, internal_boundaries

contains

!-----------------------------------------------------------------------------
!
! This routine exchanges guardcells for BND_MPI and BND_PER boundaries.in u(:,:,:,:), b(:,:,:,:) and rank-3 arrays passed through argument list
! (preferably only pointers to the actual arrays are passed).
! The corners should be properly updated if this%[io]_bnd(:, ind) was set up appropriately (MPI_Waitall is called separately for each dimension).
!

   subroutine internal_boundaries(ind, nb, pa3d, pa4d)

      use constants,  only: FLUID, MAG, CR, ARR, LO, HI, xdim, ydim, zdim, I_ONE
      use dataio_pub, only: die, warn
      use domain,     only: has_dir, cdd
      use grid,       only: all_cg
      use grid_cont,  only: grid_container
      use gc_list,    only: cg_list_element
      use mpi,        only: MPI_COMM_NULL
      use mpisetup,   only: comm, ierr, proc, req, status

      implicit none

      integer(kind=4), intent(in) :: ind   !< second index in [io]_bnd arrays
      integer, optional, intent(in) :: nb !< number of grid cells to exchange (not implemented for comm3d)
      real, optional, pointer, dimension(:,:,:)   :: pa3d
      real, optional, pointer, dimension(:,:,:,:) :: pa4d

      integer :: g, d, n
      integer(kind=4) :: nr
      integer(kind=8), dimension(xdim:zdim, LO:HI) :: ise, ose
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

!BEWARE: MPI_Waitall should be called after all grid containers post Isends and Irecvs
! This routine thus cannot be a metod of the grid_container type

      if (ind < minval([FLUID, MAG, CR, ARR]) .or. ind > maxval([FLUID, MAG, CR, ARR])) call die("[grid_container:internal_boundaries] wrong index")

      if (cdd%comm3d /= MPI_COMM_NULL) then
         call warn("[internal_bnd:internal_boundaries] comm3d is implemented somewhere else.")
         return
         ! ToDo: move comm3d variants here
      endif

      select case (ind)
         case (FLUID, MAG, CR)
            if (.not. present(pa4d)) call die("[grid_container:internal_boundaries] pa4d not provided")
            if (.not. associated(pa4d)) call die("[grid_container:internal_boundaries] pa4d == null()")
         case (ARR)
            if (.not. present(pa3d)) call die("[grid_container:internal_boundaries] pa3d not provided")
            if (.not. associated(pa3d)) call die("[grid_container:internal_boundaries] pa3d == null()")
         case default
            call die("[internal_bnd:internal_boundaries] wrong array")
            return
      end select
      if (present(pa3d) .and. present(pa4d)) call die("[grid_container:internal_boundaries] Both pa3d and pa4d are present")

      cgl => all_cg%first
      do while (associated(cgl))
         cg => cgl%cg

         n = cg%nb
         if (present(nb)) then
            n = nb
            if (n<=0 .or. n>cg%nb) call die("[internal_bnd:internal_boundaries] wrong number of guardcell layers")
         endif

         do d = xdim, zdim
            nr = 0
            if (has_dir(d)) then

               if (allocated(cg%i_bnd(d, ind, n)%seg)) then
                  if (.not. allocated(cg%o_bnd(d, ind, n)%seg)) call die("[internal_bnd:internal_boundaries] cg%i_bnd without cg%o_bnd")
                  if (ubound(cg%i_bnd(d, ind, n)%seg(:), dim=1) /= ubound(cg%o_bnd(d, ind, n)%seg(:), dim=1)) call die("[internal_bnd:internal_boundaries] cg%i_bnd differs in number of entries from cg%o_bnd")
                  do g = 1, ubound(cg%i_bnd(d, ind, n)%seg(:), dim=1)
                     if (proc == cg%i_bnd(d, ind, n)%seg(g)%proc) then
                        ise = cg%i_bnd(d, ind, n)%seg(g)%se
                        ose(:,:) = ise(:,:)
                        if (ise(d, LO) < cg%n_b(d)) then
                           ose(d, :) = ise(d, :) + cg%n_b(d)
                        else
                           ose(d, :) = ise(d, :) - cg%n_b(d)
                        endif
                        ! boundaries are always paired
                        if (ind == ARR) then
                           pa3d     (ise(xdim, LO):ise(xdim,HI), ise(ydim, LO):ise(ydim, HI), ise(zdim, LO):ise(zdim, HI)) = &
                                pa3d(ose(xdim, LO):ose(xdim,HI), ose(ydim, LO):ose(ydim, HI), ose(zdim, LO):ose(zdim, HI))
                        else
                           pa4d     (:, ise(xdim, LO):ise(xdim,HI), ise(ydim, LO):ise(ydim, HI), ise(zdim, LO):ise(zdim, HI)) = &
                                pa4d(:, ose(xdim, LO):ose(xdim,HI), ose(ydim, LO):ose(ydim, HI), ose(zdim, LO):ose(zdim, HI))
                        endif
                     else
                        ! BEWARE: Here we assume, that we have at most one chunk to communicate with a given process on a single side od the domain.
                        ! This will not be true when we allow many blocks per process and tag will need to be modified to include g or seg(g)%lh should become seg(g)%tag
                        nr = nr + I_ONE
                        if (ind == ARR) then
                           call MPI_Irecv(pa3d(1, 1, 1), I_ONE, cg%i_bnd(d, ind, n)%seg(g)%mbc, cg%i_bnd(d, ind, n)%seg(g)%proc, cg%i_bnd(d, ind, n)%seg(g)%tag, comm, req(nr), ierr)
                        else
                           call MPI_Irecv(pa4d(1, 1, 1, 1), I_ONE, cg%i_bnd(d, ind, n)%seg(g)%mbc, cg%i_bnd(d, ind, n)%seg(g)%proc, cg%i_bnd(d, ind, n)%seg(g)%tag, comm, req(nr), ierr)
                        endif
                     endif
                  enddo
               else
                  if (allocated(cg%o_bnd(d, ind, n)%seg)) call die("[grid_container:internal_boundaries] cg%o_bnd without cg%i_bnd")
               endif
               if (allocated(cg%o_bnd(d, ind, n)%seg)) then
                  do g = 1, ubound(cg%o_bnd(d, ind, n)%seg(:), dim=1)
                     if (proc /= cg%o_bnd(d, ind, n)%seg(g)%proc) then
                        nr = nr + I_ONE
                        ! for noncartesian division some y-boundary corner cells are independent from x-boundary face cells, (similarly for z-direction).
                        if (ind == ARR) then
                           call MPI_Isend(pa3d(1, 1, 1), I_ONE, cg%o_bnd(d, ind, n)%seg(g)%mbc, cg%o_bnd(d, ind, n)%seg(g)%proc, cg%o_bnd(d, ind, n)%seg(g)%tag, comm, req(nr), ierr)
                        else
                           call MPI_Isend(pa4d(1, 1, 1, 1), I_ONE, cg%o_bnd(d, ind, n)%seg(g)%mbc, cg%o_bnd(d, ind, n)%seg(g)%proc, cg%o_bnd(d, ind, n)%seg(g)%tag, comm, req(nr), ierr)
                        endif
                     endif
                  enddo
               endif
               if (ubound(cg%i_bnd(d, ind, n)%seg(:), dim=1) /= ubound(cg%o_bnd(d, ind, n)%seg(:), dim=1)) call die("g:ib u/=u")
               if (nr>0) call MPI_Waitall(nr, req(:nr), status(:,:nr), ierr)
            endif
         enddo

         cgl => cgl%nxt
      enddo

   end subroutine internal_boundaries

!-----------------------------------------------------------------------------
!
! This routine sets up all guardcells (internal and external) for given rank-3 arrays.
!

   subroutine arr3d_boundaries(pa3d, nb, area_type, dname)

      use constants,  only: ARR, xdim, ydim, zdim, LO, HI, BND, BLK, BND_PER, BND_MPI, BND_SHE, BND_COR, AT_NO_B, I_ONE
      use dataio_pub, only: die, msg
      use domain,     only: has_dir, cdd, is_multicg
      use gc_list,    only: cg_list_element
      use grid,       only: all_cg
      use grid_cont,  only: grid_container
      use mpi,        only: MPI_REQUEST_NULL, MPI_IN_PLACE, MPI_LOGICAL, MPI_LOR, MPI_COMM_NULL
      use mpisetup,   only: ierr, comm, proc, req, status

      implicit none

      real, dimension(:,:,:), pointer, intent(inout) :: pa3d
      integer, optional, intent(in) :: nb !< number of grid cells to exchange (not implemented for comm3d)
      integer(kind=4), intent(in), optional          :: area_type
      character(len=*), intent(in), optional         :: dname

      integer :: i, d, n
      integer(kind=4) :: lh
      logical :: dodie, do_permpi
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      dodie = .false.

      !> \todo fill corners with big_float ?

      if (is_multicg) call die("[grid:arr3d_boundaries] multiple grid pieces per procesor not implemented yet") !nontrivial MPI_Waitall should be outside do while (associated(cgl)) loop

      if (cdd%comm3d == MPI_COMM_NULL) then

         do_permpi = .true.
         if (present(area_type)) then
            if (area_type /= AT_NO_B) do_permpi = .false.
         endif

         if (do_permpi) call internal_boundaries(ARR, nb=nb, pa3d=pa3d)

      endif

      cgl => all_cg%first
      do while (associated(cgl))
         cg => cgl%cg

         n = cg%nb
         if (present(nb)) then
            n = nb
            if (n<=0 .or. n>cg%nb) call die("[internal_bnd:arr3d_boundaries] wrong number of guardcell layers")
         endif

         req(:) = MPI_REQUEST_NULL

         do d = xdim, zdim
            if (has_dir(d)) then
               do lh = LO, HI

                  select case (cg%bnd(d, lh))
                     case (BND_PER)
                        if (cdd%comm3d /= MPI_COMM_NULL) then
                           if (present(area_type)) then
                              if (area_type /= AT_NO_B) cycle
                           endif
                           do i = 1, ceiling(n/real(cg%n_b(d))) ! Repeating is important for domains that are narrower than their guardcells (e.g. cg%n_b(d) = 2)
                              select case (2*d+lh)
                                 case (2*xdim+LO)
                                    pa3d(1:cg%nb, :, :) = pa3d(cg%ieb:cg%ie, :, :) ! local copy is cheap (and don't occur so often in large runs) so don't boyher with the value of n
                                 case (2*ydim+LO)
                                    pa3d(:, 1:cg%nb, :) = pa3d(:, cg%jeb:cg%je, :)
                                 case (2*zdim+LO)
                                    pa3d(:, :, 1:cg%nb) = pa3d(:, :, cg%keb:cg%ke)
                                 case (2*xdim+HI)
                                    pa3d(cg%ie+1:cg%n_(xdim), :, :) = pa3d(cg%is:cg%isb, :, :)
                                 case (2*ydim+HI)
                                    pa3d(:, cg%je+1:cg%n_(ydim), :) = pa3d(:, cg%js:cg%jsb, :)
                                 case (2*zdim+HI)
                                    pa3d(:, :, cg%ke+1:cg%n_(zdim)) = pa3d(:, :, cg%ks:cg%ksb)
                              end select
                           enddo
                        endif
                     case (BND_MPI)
                        if (cdd%comm3d /= MPI_COMM_NULL) then
                           if (cdd%psize(d) > 1) then
                              call MPI_Isend(pa3d(1, 1, 1), I_ONE, cg%mbc(ARR, d, lh, BLK, n), cdd%procn(d, lh), int(2*d+(LO+HI-lh), kind=4), cdd%comm3d, req(4*(d-xdim)+1+2*(lh-LO)), ierr)
                              call MPI_Irecv(pa3d(1, 1, 1), I_ONE, cg%mbc(ARR, d, lh, BND, n), cdd%procn(d, lh), int(2*d+       lh,  kind=4), cdd%comm3d, req(4*(d-xdim)+2+2*(lh-LO)), ierr)
                           else
                              call die("[grid:arr3d_boundaries] bnd_[xyz][lr] == 'mpi' && cdd%psize([xyz]dim) <= 1")
                           endif
                        endif
                     case (BND_SHE) !> \todo move appropriate code from poissonsolver::poisson_solve or do nothing. or die until someone really needs SHEAR.
                        write(msg,*) "[grid:arr3d_boundaries] 'she' not implemented for ",dname
                        dodie = .true.
                     case (BND_COR)
                        if (present(area_type)) then
                           if (area_type /= AT_NO_B) cycle
                        endif
                        write(msg,*) "[grid:arr3d_boundaries] 'cor' not implemented for ", dname
                        dodie = .true.
                     case default ! Set gradient == 0 on the external boundaries
                        if (present(area_type)) then
                           if (area_type /= AT_NO_B) cycle
                        endif
                        do i = 1, cg%nb
                           select case (2*d+lh)
                              case (2*xdim+LO)
                                 pa3d(i, :, :) = pa3d(cg%is, :, :)
                              case (2*ydim+LO)
                                 pa3d(:, i, :) = pa3d(:, cg%js, :)
                              case (2*zdim+LO)
                                 pa3d(:, :, i) = pa3d(:, :, cg%ks)
                              case (2*xdim+HI)
                                 pa3d(cg%ie+i, :, :) = pa3d(cg%ie, :, :)
                              case (2*ydim+HI)
                                 pa3d(:, cg%je+i, :) = pa3d(:, cg%je, :)
                              case (2*zdim+HI)
                                 pa3d(:, :, cg%ke+i) = pa3d(:, :, cg%ke)
                           end select
                        enddo
                  end select

               enddo
            endif
            !> \warning outside xdim-zdim loop MPI_Waitall may change the operations order and as a result may leave mpi-corners uninitiallized
            if (cdd%comm3d /= MPI_COMM_NULL) call MPI_Waitall(size(req(:)), req(:), status(:,:), ierr)
         enddo

         cgl => cgl%nxt
      enddo

      call MPI_Allreduce(MPI_IN_PLACE, dodie, I_ONE, MPI_LOGICAL, MPI_LOR, comm, ierr)
      if (dodie) call die(msg)

   end subroutine arr3d_boundaries

end module internal_bnd
