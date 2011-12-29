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

!>
!! \brief Update external boundaries
!!
!! \details This module contains subroutines that are responsible for preparing external boundary cells for all simple boundary types for cg%q(:) and cg%w(:) arrays.
!! (i.e. no communication, dedicated arrays etc). Fancy, specialized boundary conditions should be defined somewhere else, in appropriate modules.
!!
!! Note that this routine may not properly update some layers of guardcells when number of guardcell layers exceedes number of active cells.
!! Appropriate checks should be made in divide_domain routine.
!!
!! \todo integrate here as much stuff from fluidboundaries, magboundaries, etc. as possible.
!<

module external_bnd
! pulled by ANY

   implicit none

   private
   public :: arr3d_boundaries

contains

!-----------------------------------------------------------------------------
!
! This routine sets up all guardcells (internal and external) for given rank-3 arrays.
!

   subroutine arr3d_boundaries(ind, nb, area_type)

      use constants,    only: ARR, xdim, ydim, zdim, LO, HI, BND, BLK, BND_PER, BND_MPI, BND_SHE, BND_COR, AT_NO_B, I_ONE
      use dataio_pub,   only: die, msg
      use domain,       only: dom
      use gc_list,      only: cg_list_element
      use grid,         only: leaves, all_cg
      use grid_cont,    only: grid_container
      use internal_bnd, only: internal_boundaries_3d
      use mpi,          only: MPI_REQUEST_NULL, MPI_IN_PLACE, MPI_LOGICAL, MPI_LOR, MPI_COMM_NULL
      use mpisetup,     only: ierr, comm, proc, req, status
      use types,        only: cdd

      implicit none

      integer(kind=4), intent(in) :: ind  !> Negative value: index of cg%q(:) 3d array
      integer, optional, intent(in) :: nb !> number of grid cells to exchange (not implemented for comm3d)
      integer(kind=4), intent(in), optional          :: area_type

      integer :: i, d, n
      integer(kind=4) :: lh
      logical :: dodie, do_permpi
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg
      real, dimension(:,:,:), pointer :: pa3d

      dodie = .false.

      !> \todo fill corners with big_float ?

      if (cdd%comm3d == MPI_COMM_NULL) then

         do_permpi = .true.
         if (present(area_type)) then
            if (area_type /= AT_NO_B) do_permpi = .false.
         endif

         if (do_permpi) call internal_boundaries_3d(all_cg, ind, nb=nb) ! should be more selective (modified leaves?)

      endif

      if (ind > ubound(all_cg%q_lst(:), dim=1) .or. ind < lbound(all_cg%q_lst(:), dim=1)) call die("[internal_bnd:arr3d_boundaries] wrong 3d index")

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         n = dom%nb
         if (present(nb)) then
            n = nb
            if (n<=0 .or. n>dom%nb) call die("[internal_bnd:arr3d_boundaries] wrong number of guardcell layers")
         endif

         req(:) = MPI_REQUEST_NULL

         pa3d =>cg%q(ind)%arr

         do d = xdim, zdim
            if (dom%has_dir(d)) then
               do lh = LO, HI

                  select case (cg%bnd(d, lh))
                     case (BND_PER)
                        if (cdd%comm3d /= MPI_COMM_NULL) then
                           if (present(area_type)) then
                              if (area_type /= AT_NO_B) cycle
                           endif
                           select case (2*d+lh)
                              case (2*xdim+LO)
                                 pa3d(1:dom%nb, :, :) = pa3d(cg%ieb:cg%ie, :, :) ! local copy is cheap (and don't occur so often in large runs) so don't boyher with the value of n
                              case (2*ydim+LO)
                                 pa3d(:, 1:dom%nb, :) = pa3d(:, cg%jeb:cg%je, :)
                              case (2*zdim+LO)
                                 pa3d(:, :, 1:dom%nb) = pa3d(:, :, cg%keb:cg%ke)
                              case (2*xdim+HI)
                                 pa3d(cg%ie+1:cg%n_(xdim), :, :) = pa3d(cg%is:cg%isb, :, :)
                              case (2*ydim+HI)
                                 pa3d(:, cg%je+1:cg%n_(ydim), :) = pa3d(:, cg%js:cg%jsb, :)
                              case (2*zdim+HI)
                                 pa3d(:, :, cg%ke+1:cg%n_(zdim)) = pa3d(:, :, cg%ks:cg%ksb)
                           end select
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
                        write(msg,*) "[grid:arr3d_boundaries] 'she' not implemented for ",all_cg%q_lst(ind)%name
                        dodie = .true.
                     case (BND_COR)
                        if (present(area_type)) then
                           if (area_type /= AT_NO_B) cycle
                        endif
                        write(msg,*) "[grid:arr3d_boundaries] 'cor' not implemented for ", all_cg%q_lst(ind)%name
                        dodie = .true.
                     case default ! Set gradient == 0 on the external boundaries
                        if (present(area_type)) then
                           if (area_type /= AT_NO_B) cycle
                        endif
                        do i = 1, dom%nb
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

end module external_bnd
