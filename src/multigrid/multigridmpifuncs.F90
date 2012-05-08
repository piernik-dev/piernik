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
!! \brief This module is responsible for setting boundaries, either by MPI communication or by estimating external boundaries.
!<

module multigridmpifuncs
! pulled by MULTIGRID

   implicit none

   private
   public :: mpi_multigrid_bnd

contains


!!$ ============================================================================
!!
!! Boundary conditions
!!
!> \brief Routine for inter-process and periodic boundary conditions.
!> \details mpi_multigrid_bnd provides communication between local domains to couple solution on the global computational domain
!!

   subroutine mpi_multigrid_bnd(curl, iv, ng, mode, corners)

      use constants,     only: xdim, ydim, zdim, LO, HI, BND, BLK, ARR, INT4, I_ONE, I_FOUR
      use dataio_pub,    only: die
      use domain,        only: is_mpi_noncart, dom
      use cg_list_lev,   only: cg_list_level
      use internal_bnd,  only: internal_boundaries_3d
      use mpi,           only: MPI_REQUEST_NULL, MPI_COMM_NULL
      use mpisetup,      only: ierr, have_mpi, req, status
      use types,         only: cdd

      implicit none

      type(cg_list_level), pointer, intent(in) :: curl    !< level which we are doing communication at
      integer,                      intent(in) :: iv      !< variable which we want to communicate
      integer(kind=4),              intent(in) :: ng      !< number of guardcells to exchange
      integer(kind=4),              intent(in) :: mode    !< what to do with external boundaries
      logical, optional,            intent(in) :: corners !< if .true. then don't forget aboutpay close attention to corners

      integer(kind=4), parameter :: dreq = I_FOUR
      logical :: cor
      integer(kind=4) :: d, doff

      if (.not. associated(curl)) call die("[multigridmpifuncs:mpi_multigrid_bnd] Invalid level")
      if (ng > dom%nb .or. ng <= 0) call die("[multigridmpifuncs:mpi_multigrid_bnd] Too many or <0 guardcells requested.")

      if (present(corners)) then
         cor = corners
      else
         cor = .false.
      endif

      ! Set the external boundary, where appropriate
      call multigrid_ext_bnd(curl, iv, ng, mode, cor)

      if (cdd%comm3d == MPI_COMM_NULL) then

         do d = xdim, zdim
            call internal_boundaries_3d(curl, iv, dim=d)
         enddo

      else
         if (have_mpi .and. is_mpi_noncart) call die("[multigridmpifuncs:mpi_multigrid_bnd] is_mpi_noncart is not implemented") !procxl, procxr, procyl, procyr, proczl, proczr, psize,

         !! \deprecated cannot call arr3d_boundaries, because it would destroy external boundary

         req(:) = MPI_REQUEST_NULL

         do d = xdim, zdim
            if (dom%has_dir(d)) then
               doff = dreq*(d-xdim)
               if (cdd%psize(d) > 1) then ! \todo remove psize(:), try to rely on offsets or boundary types
                  if (.not. curl%first%cg%ext_bnd(d, LO)) call MPI_Isend(curl%first%cg%q(iv)%arr(1, 1, 1), I_ONE, curl%first%cg%mbc(ARR, d, LO, BLK, ng), cdd%procn(d, LO), 17_INT4+doff, cdd%comm3d, req(1+doff), ierr)
                  if (.not. curl%first%cg%ext_bnd(d, HI)) call MPI_Isend(curl%first%cg%q(iv)%arr(1, 1, 1), I_ONE, curl%first%cg%mbc(ARR, d, HI, BLK, ng), cdd%procn(d, HI), 19_INT4+doff, cdd%comm3d, req(2+doff), ierr)
                  if (.not. curl%first%cg%ext_bnd(d, LO)) call MPI_Irecv(curl%first%cg%q(iv)%arr(1, 1, 1), I_ONE, curl%first%cg%mbc(ARR, d, LO, BND, ng), cdd%procn(d, LO), 19_INT4+doff, cdd%comm3d, req(3+doff), ierr)
                  if (.not. curl%first%cg%ext_bnd(d, HI)) call MPI_Irecv(curl%first%cg%q(iv)%arr(1, 1, 1), I_ONE, curl%first%cg%mbc(ARR, d, HI, BND, ng), cdd%procn(d, HI), 17_INT4+doff, cdd%comm3d, req(4+doff), ierr)
               else
                  if (curl%first%cg%ext_bnd(d, LO) .neqv. curl%first%cg%ext_bnd(d, HI)) call die("[multigridmpifuncs:mpi_multigrid_bnd] inconsiztency in ext_bnd(:)")
                  if (.not. curl%first%cg%ext_bnd(d, LO)) then
                     select case (d)
                        case (xdim)
                           curl%first%cg%q(iv)%arr(curl%first%cg%is-ng:curl%first%cg%is-1,  :, :) = curl%first%cg%q(iv)%arr(curl%first%cg%ie-ng+1:curl%first%cg%ie,      :, :)
                           curl%first%cg%q(iv)%arr(curl%first%cg%ie+1 :curl%first%cg%ie+ng, :, :) = curl%first%cg%q(iv)%arr(curl%first%cg%is     :curl%first%cg%is+ng-1, :, :)
                        case (ydim)
                           curl%first%cg%q(iv)%arr(:, curl%first%cg%js-ng:curl%first%cg%js-1,  :) = curl%first%cg%q(iv)%arr(:, curl%first%cg%je-ng+1:curl%first%cg%je,      :)
                           curl%first%cg%q(iv)%arr(:, curl%first%cg%je+1 :curl%first%cg%je+ng, :) = curl%first%cg%q(iv)%arr(:, curl%first%cg%js     :curl%first%cg%js+ng-1, :)
                        case (zdim)
                           curl%first%cg%q(iv)%arr(:, :, curl%first%cg%ks-ng:curl%first%cg%ks-1) = curl%first%cg%q(iv)%arr(:, :, curl%first%cg%ke-ng+1:curl%first%cg%ke)
                           curl%first%cg%q(iv)%arr(:, :, curl%first%cg%ke+1 :curl%first%cg%ke+ng) = curl%first%cg%q(iv)%arr(:, :, curl%first%cg%ks     :curl%first%cg%ks+ng-1)
                     end select
                  endif
               endif
               if (cor) call MPI_Waitall(dreq, req(1+doff:4+doff), status(:,1+doff:4+doff), ierr)
            endif
         enddo

!>
!! \todo Make a benchmark of a massively parallel run to determine difference in execution between calling MPI_Waitall for each direction and calling it once.
!! If the difference is small then set cor permanently to .true.
!<
         if (.not. cor) call MPI_Waitall(size(req(:)), req(:), status(:,:), ierr)

      endif

   end subroutine mpi_multigrid_bnd

!!$ ============================================================================
!!
!> \brief Set external boundary (not required for periodic box) on domain faces.
!! \details In multigrid typically mirror boundaries are in use. Extrapolate isolated boundaries at exit.
!!
!! dom%has_dir() is not checked here because ext_bnd() should be set to .false. on non-existing directions in 1D and 2D setups
!<

   subroutine multigrid_ext_bnd(curl, iv, ng, mode, cor)

      use constants,     only: LO, HI, xdim, ydim, zdim
      use dataio_pub,    only: die, msg, warn
      use gc_list,       only: cg_list_element
      use cg_list_lev,   only: cg_list_level
      use grid_cont,     only: grid_container
      use multigridvars, only: extbnd_donothing, extbnd_zero, extbnd_extrapolate, extbnd_mirror, extbnd_antimirror

      implicit none

      type(cg_list_level), pointer, intent(in) :: curl  !< level which we are preparing the guardcells at
      integer,                      intent(in) :: iv    !< variable which we want to set
      integer(kind=4),              intent(in) :: ng    !< number of guardcells to set
      integer(kind=4),              intent(in) :: mode  !< what to do with external boundaries
      logical,                      intent(in) :: cor   !< if .true. then don't forget about corners \deprecated BEWARE: not implemented properly

      integer :: i
      logical, save :: warned = .false.
      type(grid_container),  pointer :: cg
      type(cg_list_element), pointer :: cgl

      cgl => curl%first
      do while (associated(cgl))
         cg => cgl%cg

         if (cor) then
            if (.not. warned) then
               call warn("[multigridmpifuncs:multigrid_ext_bnd] some corners may not be properly filled. FIX ME!")
               warned = .true.
            endif
         endif

         !> \todo implement mixed BC
         select case (mode)         !> \deprecated BEWARE: some cylindrical factors may be helpful
            case (extbnd_donothing) ! remember to initialize everything first!
               return
            case (extbnd_extrapolate) !> \deprecated mixed-tybe BC: free flux; BEWARE: it is not protected from inflow
               do i = 1, ng
                  if (cg%ext_bnd(xdim, LO)) cg%q(iv)%arr(cg%is-i, :, :) = (1+i) * cg%q(iv)%arr(cg%is, :, :) - i * cg%q(iv)%arr(cg%is+1, :, :)
                  if (cg%ext_bnd(xdim, HI)) cg%q(iv)%arr(cg%ie+i, :, :) = (1+i) * cg%q(iv)%arr(cg%ie, :, :) - i * cg%q(iv)%arr(cg%ie-1, :, :)
                  if (cg%ext_bnd(ydim, LO)) cg%q(iv)%arr(:, cg%js-i, :) = (1+i) * cg%q(iv)%arr(:, cg%js, :) - i * cg%q(iv)%arr(:, cg%js+1, :)
                  if (cg%ext_bnd(ydim, HI)) cg%q(iv)%arr(:, cg%je+i, :) = (1+i) * cg%q(iv)%arr(:, cg%je, :) - i * cg%q(iv)%arr(:, cg%je-1, :)
                  if (cg%ext_bnd(zdim, LO)) cg%q(iv)%arr(:, :, cg%ks-i) = (1+i) * cg%q(iv)%arr(:, :, cg%ks) - i * cg%q(iv)%arr(:, :, cg%ks+1)
                  if (cg%ext_bnd(zdim, HI)) cg%q(iv)%arr(:, :, cg%ke+i) = (1+i) * cg%q(iv)%arr(:, :, cg%ke) - i * cg%q(iv)%arr(:, :, cg%ke-1)
               enddo
            case (extbnd_zero) ! homogenous Dirichlet BC with 0 at first guardcell row
               if (cg%ext_bnd(xdim, LO)) cg%q(iv)%arr(:cg%is, :, :) = 0.
               if (cg%ext_bnd(xdim, HI)) cg%q(iv)%arr(cg%ie:, :, :) = 0.
               if (cg%ext_bnd(ydim, LO)) cg%q(iv)%arr(:, :cg%js, :) = 0.
               if (cg%ext_bnd(ydim, HI)) cg%q(iv)%arr(:, cg%je:, :) = 0.
               if (cg%ext_bnd(zdim, LO)) cg%q(iv)%arr(:, :, :cg%ks) = 0.
               if (cg%ext_bnd(zdim, HI)) cg%q(iv)%arr(:, :, cg%ke:) = 0.
            case (extbnd_mirror) ! reflecting BC (homogenous Neumamnn)
               do i = 1, ng
                  if (cg%ext_bnd(xdim, LO)) cg%q(iv)%arr(cg%is-i, :, :) = cg%q(iv)%arr(cg%is+i-1, :, :)
                  if (cg%ext_bnd(xdim, HI)) cg%q(iv)%arr(cg%ie+i, :, :) = cg%q(iv)%arr(cg%ie-i+1, :, :)
                  if (cg%ext_bnd(ydim, LO)) cg%q(iv)%arr(:, cg%js-i, :) = cg%q(iv)%arr(:, cg%js+i-1, :)
                  if (cg%ext_bnd(ydim, HI)) cg%q(iv)%arr(:, cg%je+i, :) = cg%q(iv)%arr(:, cg%je-i+1, :)
                  if (cg%ext_bnd(zdim, LO)) cg%q(iv)%arr(:, :, cg%ks-i) = cg%q(iv)%arr(:, :, cg%ks+i-1)
                  if (cg%ext_bnd(zdim, HI)) cg%q(iv)%arr(:, :, cg%ke+i) = cg%q(iv)%arr(:, :, cg%ke-i+1)
               enddo
            case (extbnd_antimirror) ! homogenous Dirichlet BC with 0 at external faces
               do i = 1, ng
                  if (cg%ext_bnd(xdim, LO)) cg%q(iv)%arr(cg%is-i, :, :) = - cg%q(iv)%arr(cg%is+i-1, :, :)
                  if (cg%ext_bnd(xdim, HI)) cg%q(iv)%arr(cg%ie+i, :, :) = - cg%q(iv)%arr(cg%ie-i+1, :, :)
                  if (cg%ext_bnd(ydim, LO)) cg%q(iv)%arr(:, cg%js-i, :) = - cg%q(iv)%arr(:, cg%js+i-1, :)
                  if (cg%ext_bnd(ydim, HI)) cg%q(iv)%arr(:, cg%je+i, :) = - cg%q(iv)%arr(:, cg%je-i+1, :)
                  if (cg%ext_bnd(zdim, LO)) cg%q(iv)%arr(:, :, cg%ks-i) = - cg%q(iv)%arr(:, :, cg%ks+i-1)
                  if (cg%ext_bnd(zdim, HI)) cg%q(iv)%arr(:, :, cg%ke+i) = - cg%q(iv)%arr(:, :, cg%ke-i+1)
               enddo
            case default
               write(msg, '(a,i3,a)')"[multigridmpifuncs:multigrid_ext_bnd] boundary type ",mode," not implemented"
               call die(msg)
         end select

         cgl => cgl%nxt
      enddo

   end subroutine multigrid_ext_bnd

end module multigridmpifuncs
