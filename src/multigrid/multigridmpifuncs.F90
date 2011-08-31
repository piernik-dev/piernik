! $Id$
!
! PIERNIK Code Copyright (C) 2006 Michal Hanasz
!
!    This file is part of PIERNIK code.
!
!    PIERNIK is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published
!    by
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
!    Initial implementation of PIERNIK code was based on TVD split MHD
!    code by
!    Ue-Li Pen
!        see: Pen, Arras & Wong (2003) for algorithm and
!             http://www.cita.utoronto.ca/~pen/MHD
!             for original source code "mhd.f90"
!
!    For full list of developers see $PIERNIK_HOME/license/pdt.txt
!

#include "piernik.h"

!!$ ============================================================================
!>
!! \brief This module is responsible for setting boundaries, either by MPI communication or by estimating external boundaries.
!<

module multigridmpifuncs
! pulled by MULTIGRID
   implicit none

   private

   public :: mpi_multigrid_prep, mpi_multigrid_bnd

contains

!!$ ============================================================================
!> \brief Initialize MPI shortcuts for communication. Called from init_multigrid.

   subroutine mpi_multigrid_prep

      use constants,     only: xdim, ydim, zdim, LO, HI, BND, BLK, ndims, INVALID, I_ZERO
      use dataio_pub,    only: warn, die
      use domain,        only: is_overlap, has_dir, cdd
      use mpi,           only: MPI_DOUBLE_PRECISION, MPI_ORDER_FORTRAN, MPI_COMM_NULL
      use mpisetup,      only: ierr, proc, FIRST, LAST, procmask
      use multigridvars, only: plvl, base, pr_segment

      implicit none

      integer(kind=4):: d, ib, lh, hl
      integer :: g, j
      integer(kind=4), dimension(ndims) :: sizes, subsizes, starts
      logical :: sharing
      integer(kind=8), dimension(xdim:zdim) :: ijks, per
      integer(kind=8), dimension(xdim:zdim, LO:HI) :: coarsened, b_layer, bp_layer, poff
      type(pr_segment), pointer :: seg
      type(plvl), pointer :: curl

      curl => base
      do while (associated(curl))
         if (ubound(curl%dom%pse(proc)%sel(:,:,:), dim=1) > 1) call die("[multigrid:mpi_multigrid_prep] Multiple blocks per process not implemented yet")

         ijks(:) = curl%ijkse(:, LO) - curl%off(:)  ! add this to convert absolute cell coordinates to local indices. (+mg_nb - off(:))

         ! find fine target for receiving restricted data or sending data to be prolonged
         if (associated(curl%finer)) then
            procmask(:) = 0
            do j = FIRST, LAST
               coarsened(:,:) = curl%finer%dom%pse(j)%sel(1, :, :)/2
               call is_overlap(curl%my_se(:, :), coarsened(:,:), sharing)
               if (sharing) procmask(j) = 1 ! we can store also neigh(:,:), face and corner as a bitmask, if necessary
            enddo
            allocate(curl%f_tgt%seg(count(procmask(:) /= 0)))

            g = 0
            do j = FIRST, LAST
               if (procmask(j) /= 0) then
                  g = g + 1
                  if (.not. allocated(curl%f_tgt%seg) .or. g>ubound(curl%f_tgt%seg, dim=1)) call die("m:im f_tgt g>")
                  seg => curl%f_tgt%seg(g)
                  if (allocated(seg%buf)) then
                     call warn("m:mmp i seg%buf a a")
                     deallocate(seg%buf)
                  endif
                  seg%proc = j
                  ! find cross-section of own segment with coarsened fine segment
                  seg%se(:, LO) = max(curl%my_se(:, LO), curl%finer%dom%pse(j)%sel(1, :, LO)/2) + ijks(:)
                  seg%se(:, HI) = min(curl%my_se(:, HI), curl%finer%dom%pse(j)%sel(1, :, HI)/2) + ijks(:)
                  if (j /= proc) allocate(seg%buf(seg%se(xdim, HI)-seg%se(xdim, LO) + 1, seg%se(ydim, HI)-seg%se(ydim, LO) + 1, seg%se(zdim, HI)-seg%se(zdim, LO) + 1))
                  ! not counted in mb_alloc
               endif
            enddo
         endif

         ! find coarse target for sending restricted data or receiving data to be prolonged
         !> \deprecated almost duplicated code
         if (associated(curl%coarser)) then
            procmask(:) = 0
            coarsened(:,:) = curl%my_se(:, :)/2
            do j = FIRST, LAST
               call is_overlap(coarsened(:,:), curl%coarser%dom%pse(j)%sel(1, :, :), sharing)
               if (sharing) procmask(j) = 1
            enddo
            allocate(curl%c_tgt%seg(count(procmask(:) /= 0)))

            g = 0
            do j = FIRST, LAST
               if (procmask(j) /= 0) then
                  g = g + 1
                  if (.not. allocated(curl%c_tgt%seg) .or. g>ubound(curl%c_tgt%seg, dim=1)) call die("m:im c_tgt g>")
                  seg => curl%c_tgt%seg(g)
                  if (allocated(seg%buf)) then
                     call warn("m:mmp o seg%buf a a")
                     deallocate(seg%buf)
                  endif
                  seg%proc = j
                  ! find cross-section of own segment with refined coarse segment
                  seg%se(:, LO) = max(curl%my_se(:, LO), curl%coarser%dom%pse(j)%sel(1, :, LO)*2  )
                  seg%se(:, HI) = min(curl%my_se(:, HI), curl%coarser%dom%pse(j)%sel(1, :, HI)*2+1)
                  if (j /= proc) allocate(seg%buf(seg%se(xdim, HI)/2-seg%se(xdim, LO)/2 + 1, &
                       &                          seg%se(ydim, HI)/2-seg%se(ydim, LO)/2 + 1, &
                       &                          seg%se(zdim, HI)/2-seg%se(zdim, LO)/2 + 1))
                  ! not counted in mb_alloc
                  seg%se(:, LO) = seg%se(:, LO) + ijks(:)
                  seg%se(:, HI) = seg%se(:, HI) + ijks(:)
               endif
            enddo
         endif

         !BEWARE: almost replicated code (grid:grid_mpi_boundaries_prep)
         ! find neighbours and set up the MPI containers
         if (cdd%comm3d == MPI_COMM_NULL) then

            ! assume that cuboids fill the domain and don't collide

            curl%mmbc(:, :, :, :) = INVALID

            per(:) = 0
            where (curl%dom%periodic(:)) per(:) = curl%dom%n_d(:)

            if (allocated(curl%i_bnd) .or. allocated(curl%o_bnd)) call die("[multigrid:mpi_multigrid_prep] curl%i_bnd or curl%o_bnd already allocated")
            allocate(curl%i_bnd(xdim:zdim, curl%nb), curl%o_bnd(xdim:zdim, curl%nb))

            do d = xdim, zdim
               if (has_dir(d) .and. .not. curl%empty) then

                  ! identify processes with interesting neighbour data
                  procmask(:) = 0
                  do lh = LO, HI
                     hl = LO+HI-lh ! HI for LO, LO for HI
                     b_layer(:,:) = curl%my_se(:, :)
                     b_layer(d, lh) = b_layer(d, lh) + lh-hl ! -1 for LO, +1 for HI
                     b_layer(d, hl) = b_layer(d, lh) ! boundary layer without corners
                     do j = FIRST, LAST
                        call is_overlap(b_layer(:,:), curl%dom%pse(j)%sel(1, :, :), sharing, per(:))
                        if (sharing) procmask(j) = procmask(j) + 1
                     enddo
                  enddo
                  do j = 1, curl%nb
                     allocate(curl%i_bnd(d, j)%seg(sum(procmask(:))))
                     allocate(curl%o_bnd(d, j)%seg(sum(procmask(:))))
                  enddo

                  ! set up segments to be sent or received
                  g = 0
                  do j = FIRST, LAST
                     if (procmask(j) /= 0) then
                        do lh = LO, HI
                           hl = LO+HI-lh
                           b_layer(:,:) = curl%my_se(:, :)
                           b_layer(d, lh) = b_layer(d, lh) + lh-hl
                           b_layer(d, hl) = b_layer(d, lh)
                           bp_layer(:, :) = b_layer(:, :)
                           where (per(:) > 0)
                              bp_layer(:, LO) = mod(b_layer(:, LO) + per(:), per(:))
                              bp_layer(:, HI) = mod(b_layer(:, HI) + per(:), per(:))
                           endwhere
                           call is_overlap(bp_layer(:,:), curl%dom%pse(j)%sel(1, :, :), sharing)

                           if (sharing) then
                              poff(:,:) = bp_layer(:,:) - b_layer(:,:) ! displacement due to periodicity
                              bp_layer(:, LO) = max(bp_layer(:, LO), curl%dom%pse(j)%sel(1, :, LO))
                              bp_layer(:, HI) = min(bp_layer(:, HI), curl%dom%pse(j)%sel(1, :, HI))
                              b_layer(:,:) = bp_layer(:,:) - poff(:,:)
                              g = g + 1
                              do ib = 1, curl%nb
                                 curl%i_bnd(d, ib)%seg(g)%mbc = INVALID
                                 curl%i_bnd(d, ib)%seg(g)%proc = j
                                 curl%i_bnd(d, ib)%seg(g)%se(:,LO) = b_layer(:, LO) + ijks(:)
                                 curl%i_bnd(d, ib)%seg(g)%se(:,HI) = b_layer(:, HI) + ijks(:)
                                 if (any(curl%i_bnd(d, ib)%seg(g)%se(d, :) < 0)) &
                                      curl%i_bnd(d, ib)%seg(g)%se(d, :) = curl%i_bnd(d, ib)%seg(g)%se(d, :) + curl%dom%n_d(d)
                                 if (any(curl%i_bnd(d, ib)%seg(g)%se(d, :) > curl%n_b(d) + 2*curl%nb)) &
                                      curl%i_bnd(d, ib)%seg(g)%se(d, :) = curl%i_bnd(d, ib)%seg(g)%se(d, :) - curl%dom%n_d(d)
                                 ! expand to cover corners (requires separate MPI_Waitall for each direction)
                                 ! \todo create separate %mbc for corner-less exchange with one MPI_Waitall (can scale better)
                                 where (has_dir(:d-1))
                                    curl%i_bnd(d, ib)%seg(g)%se(:d-1, LO) = curl%i_bnd(d, ib)%seg(g)%se(:d-1, LO) - ib
                                    curl%i_bnd(d, ib)%seg(g)%se(:d-1, HI) = curl%i_bnd(d, ib)%seg(g)%se(:d-1, HI) + ib
                                 endwhere
                                 curl%o_bnd(d, ib)%seg(g) = curl%i_bnd(d, ib)%seg(g)
                                 curl%i_bnd(d, ib)%seg(g)%tag = HI*d+lh-LO !> \todo add an id and number of local gc
                                 curl%o_bnd(d, ib)%seg(g)%tag = HI*d+hl-LO
                                 select case (lh)
                                    case (LO)
                                       curl%i_bnd(d, ib)%seg(g)%se(d, LO) = curl%i_bnd(d, ib)%seg(g)%se(d, HI) - (ib - 1)
                                       curl%o_bnd(d, ib)%seg(g)%se(d, LO) = curl%i_bnd(d, ib)%seg(g)%se(d, HI) + 1
                                       curl%o_bnd(d, ib)%seg(g)%se(d, HI) = curl%o_bnd(d, ib)%seg(g)%se(d, LO) + (ib - 1)
                                    case (HI)
                                       curl%i_bnd(d, ib)%seg(g)%se(d, HI) = curl%i_bnd(d, ib)%seg(g)%se(d, LO) + (ib - 1)
                                       curl%o_bnd(d, ib)%seg(g)%se(d, HI) = curl%i_bnd(d, ib)%seg(g)%se(d, LO) - 1
                                       curl%o_bnd(d, ib)%seg(g)%se(d, LO) = curl%o_bnd(d, ib)%seg(g)%se(d, HI) - (ib - 1)
                                 end select
                                 ! set MPI type only for non-local transfers
                                 call MPI_Type_create_subarray(ndims, curl%n_(:), &
                                      &                        int(curl%i_bnd(d, ib)%seg(g)%se(:, HI) - curl%i_bnd(d, ib)%seg(g)%se(:, LO) + 1, kind=4), &
                                      &                        int(curl%i_bnd(d, ib)%seg(g)%se(:, LO) - 1, kind=4), MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, &
                                      &                        curl%i_bnd(d, ib)%seg(g)%mbc, ierr)
                                 call MPI_Type_commit(curl%i_bnd(d, ib)%seg(g)%mbc, ierr)
                                 call MPI_Type_create_subarray(ndims, curl%n_(:), &
                                      &                        int(curl%o_bnd(d, ib)%seg(g)%se(:, HI) - curl%o_bnd(d, ib)%seg(g)%se(:, LO) + 1, kind=4), &
                                      &                        int(curl%o_bnd(d, ib)%seg(g)%se(:, LO) - 1, kind=4), MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, &
                                      &                        curl%o_bnd(d, ib)%seg(g)%mbc, ierr)
                                 call MPI_Type_commit(curl%o_bnd(d, ib)%seg(g)%mbc, ierr)
                              enddo
                           endif
                        enddo
                     endif
                  enddo
               endif
            enddo

         else

            ! cartesian decomposition
            do ib = 1, curl%nb
               sizes(:) = curl%n_(:)
               do d = xdim, zdim
                  if (has_dir(d) .and. .not. curl%empty) then
                     subsizes(:) = sizes(:)
                     subsizes(d) = ib
                     starts(:) = I_ZERO

                     starts(d) = curl%nb-ib
                     call MPI_Type_create_subarray(ndims, sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, curl%mmbc(d, LO, BND, ib), ierr)
                     call MPI_Type_commit(curl%mmbc(d, LO, BND, ib), ierr)

                     starts(d) = curl%nb
                     call MPI_Type_create_subarray(ndims, sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, curl%mmbc(d, LO, BLK, ib), ierr)
                     call MPI_Type_commit(curl%mmbc(d, LO, BLK, ib), ierr)

                     starts(d) = curl%n_b(d) + curl%nb - ib
                     call MPI_Type_create_subarray(ndims, sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, curl%mmbc(d, HI, BLK, ib), ierr)
                     call MPI_Type_commit(curl%mmbc(d, HI, BLK, ib), ierr)

                     starts(d) = curl%n_b(d) + curl%nb
                     call MPI_Type_create_subarray(ndims, sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, curl%mmbc(d, HI, BND, ib), ierr)
                     call MPI_Type_commit(curl%mmbc(d, HI, BND, ib), ierr)
                  endif
               enddo
            enddo

         endif

         curl => curl%finer
      enddo

   end subroutine mpi_multigrid_prep

!!$ ============================================================================
!!
!! Boundary conditions
!!
!> \brief Routine for inter-process and periodic boundary conditions.
!> \details mpi_multigrid_bnd provides communication between local domains to couple solution on the global computational domain
!!

   subroutine mpi_multigrid_bnd(curl, iv, ng, mode, corners)

      use constants,     only: xdim, ydim, zdim, LO, HI, BND, BLK, INT4, I_ONE, I_FOUR
      use dataio_pub,    only: die
      use domain,        only: is_mpi_noncart, cdd, has_dir
      use mpi,           only: MPI_REQUEST_NULL, MPI_COMM_NULL
      use mpisetup,      only: proc, comm, ierr, have_mpi, req, status
      use multigridvars, only: plvl, is_external, ngridvars

      implicit none

      type(plvl), pointer, intent(in) :: curl  !< level which we are doing communication at
      integer(kind=4), intent(in) :: iv        !< variable which we want to communicate
      integer, intent(in) :: ng                !< number of guardcells to exchange
      integer(kind=4), intent(in) :: mode      !< what to do with external boundaries
      logical, intent(in), optional :: corners !< if .true. then don't forget aboutpay close attention to corners

      integer(kind=4), parameter :: dreq = I_FOUR
      logical :: cor
      integer :: g
      integer(kind=4) :: d, doff
      integer(kind=8), dimension(:,:), pointer :: ise
      integer(kind=8), dimension(xdim:zdim, LO:HI) :: ose
      integer(kind=4) :: nr

      if (.not. associated(curl)) call die("[multigridmpifuncs:mpi_multigrid_bnd] Invalid level")
      if (iv < 1 .or. iv > ngridvars) call die("[multigridmpifuncs:mpi_multigrid_bnd] Invalid variable index.")
      if (ng > curl%nb .or. ng <= 0) call die("[multigridmpifuncs:mpi_multigrid_bnd] Too many or <0 guardcells requested.")

      if (present(corners)) then
         cor = corners
      else
         cor = .false.
      endif

      ! Set the external boundary, where appropriate
      if (any(is_external(:, :))) call multigrid_ext_bnd(curl, iv, ng, mode, cor)

      if (cdd%comm3d == MPI_COMM_NULL) then

         ! \todo call curl%internal_boundaries(ib, curl%mgvar)
         ! or    call curl%internal_boundaries(ib, curl%mgvar(:, :, :, iv)) (pointers in both cases)

         do d = xdim, zdim
            nr = 0
            if (has_dir(d)) then
               if (allocated(curl%i_bnd(d, ng)%seg)) then
                  do g = 1, ubound(curl%i_bnd(d, ng)%seg(:), dim=1)
                     if (proc == curl%i_bnd(d, ng)%seg(g)%proc) then
                        ise => curl%i_bnd(d, ng)%seg(g)%se
                        ose(:,:) = ise(:,:)
                        if (ise(d, LO) < curl%n_b(d)) then
                           ose(d, :) = ise(d, :) + curl%n_b(d)
                        else
                           ose(d, :) = ise(d, :) - curl%n_b(d)
                        endif
                        ! boundaries are always paired
                        curl%mgvar(ise(xdim, LO):ise(xdim,HI), ise(ydim, LO):ise(ydim, HI), ise(zdim, LO):ise(zdim, HI), iv) = &
                             curl%mgvar(ose(xdim, LO):ose(xdim,HI), ose(ydim, LO):ose(ydim, HI), ose(zdim, LO):ose(zdim, HI), iv)
                     else
                        ! BEWARE: Here we assume, that we have at most one chunk to communicate with a given process on a single side od the domain.
                        ! This will not be true when we allow many blocks per process and tag will need to be modified to include g or seg(g)%lh should become seg(g)%tag
                        nr = nr + I_ONE
                        call MPI_Irecv(curl%mgvar(1, 1, 1, iv), I_ONE, curl%i_bnd(d, ng)%seg(g)%mbc, curl%i_bnd(d, ng)%seg(g)%proc, curl%i_bnd(d, ng)%seg(g)%tag, comm, req(nr), ierr)
                     endif
                  enddo
               endif
               if (allocated(curl%o_bnd(d, ng)%seg)) then
                  do g = 1, ubound(curl%o_bnd(d, ng)%seg(:), dim=1)
                     if (proc /= curl%o_bnd(d, ng)%seg(g)%proc) then
                        nr = nr + I_ONE
                        ! if (cor) there should be MPI_Waitall for each d
                        ! for noncartesian division some y-boundary corner cells are independent from x-boundary face cells, (similarly for z-direction).
                        call MPI_Isend(curl%mgvar(1, 1, 1, iv), I_ONE, curl%o_bnd(d, ng)%seg(g)%mbc, curl%o_bnd(d, ng)%seg(g)%proc, curl%o_bnd(d, ng)%seg(g)%tag, comm, req(nr), ierr)
                     endif
                  enddo
               endif
               if (ubound(curl%i_bnd(d, ng)%seg(:), dim=1) /= ubound(curl%o_bnd(d, ng)%seg(:), dim=1)) call die("mmf u/=u")
               if (nr>0) call MPI_Waitall(nr, req(:nr), status(:,:nr), ierr)
            endif
         enddo

      else
         if (have_mpi .and. is_mpi_noncart) call die("[multigridmpifuncs:mpi_multigrid_bnd] is_mpi_noncart is not implemented") !procxl, procxr, procyl, procyr, proczl, proczr, psize,

         req(:) = MPI_REQUEST_NULL

         do d = xdim, zdim
            if (has_dir(d)) then
               doff = dreq*(d-xdim)
               if (cdd%psize(d) > 1) then ! \todo remove psize(:), try to rely on offsets or boundary types
                  if (.not. is_external(d, LO)) call MPI_Isend(curl%mgvar(1, 1, 1, iv), I_ONE, curl%mmbc(d, LO, BLK, ng), cdd%procn(d, LO), 17_INT4+doff, cdd%comm3d, req(1+doff), ierr)
                  if (.not. is_external(d, HI)) call MPI_Isend(curl%mgvar(1, 1, 1, iv), I_ONE, curl%mmbc(d, HI, BLK, ng), cdd%procn(d, HI), 19_INT4+doff, cdd%comm3d, req(2+doff), ierr)
                  if (.not. is_external(d, LO)) call MPI_Irecv(curl%mgvar(1, 1, 1, iv), I_ONE, curl%mmbc(d, LO, BND, ng), cdd%procn(d, LO), 19_INT4+doff, cdd%comm3d, req(3+doff), ierr)
                  if (.not. is_external(d, HI)) call MPI_Irecv(curl%mgvar(1, 1, 1, iv), I_ONE, curl%mmbc(d, HI, BND, ng), cdd%procn(d, HI), 17_INT4+doff, cdd%comm3d, req(4+doff), ierr)
               else
                  if (is_external(d, LO) .neqv. is_external(d, HI)) call die("[multigridmpifuncs:mpi_multigrid_bnd] inconsiztency in is_external(:)")
                  if (.not. is_external(d, LO)) then
                     select case (d)
                        case (xdim)
                           curl%mgvar(curl%is-ng:curl%is-1,  :, :, iv) = curl%mgvar(curl%ie-ng+1:curl%ie,      :, :, iv)
                           curl%mgvar(curl%ie+1 :curl%ie+ng, :, :, iv) = curl%mgvar(curl%is     :curl%is+ng-1, :, :, iv)
                        case (ydim)
                           curl%mgvar(:, curl%js-ng:curl%js-1,  :, iv) = curl%mgvar(:, curl%je-ng+1:curl%je,      :, iv)
                           curl%mgvar(:, curl%je+1 :curl%je+ng, :, iv) = curl%mgvar(:, curl%js     :curl%js+ng-1, :, iv)
                        case (zdim)
                           curl%mgvar(:, :, curl%ks-ng:curl%ks-1,  iv) = curl%mgvar(:, :, curl%ke-ng+1:curl%ke,      iv)
                           curl%mgvar(:, :, curl%ke+1 :curl%ke+ng, iv) = curl%mgvar(:, :, curl%ks     :curl%ks+ng-1, iv)
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
!! has_dir() is not checked here because is_external() should be set to .false. on non-existing directions in 1D and 2D setups
!<

   subroutine multigrid_ext_bnd(curl, iv, ng, mode, cor)

      use constants,       only: LO, HI, xdim, ydim, zdim
      use dataio_pub,      only: die, msg, warn
      use multigridvars,   only: extbnd_donothing, extbnd_zero, extbnd_extrapolate, extbnd_mirror, extbnd_antimirror, plvl, is_external

      implicit none

      type(plvl), pointer, intent(in) :: curl !< level which we are preparing the guardcells at
      integer(kind=4), intent(in) :: iv      !< variable which we want to set
      integer, intent(in) :: ng              !< number of guardcells to set
      integer(kind=4), intent(in) :: mode    !< what to do with external boundaries
      logical, intent(in) :: cor             !< if .true. then don't forget about corners \deprecated BEWARE: not implemented properly

      integer :: i
      logical, save :: warned = .false.

      if (curl%empty) return

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
               if (is_external(xdim, LO)) curl%mgvar(curl%is-i, :, :, iv) = (1+i) * curl%mgvar(curl%is, :, :, iv) - i * curl%mgvar(curl%is+1, :, :, iv)
               if (is_external(xdim, HI)) curl%mgvar(curl%ie+i, :, :, iv) = (1+i) * curl%mgvar(curl%ie, :, :, iv) - i * curl%mgvar(curl%ie-1, :, :, iv)
               if (is_external(ydim, LO)) curl%mgvar(:, curl%js-i, :, iv) = (1+i) * curl%mgvar(:, curl%js, :, iv) - i * curl%mgvar(:, curl%js+1, :, iv)
               if (is_external(ydim, HI)) curl%mgvar(:, curl%je+i, :, iv) = (1+i) * curl%mgvar(:, curl%je, :, iv) - i * curl%mgvar(:, curl%je-1, :, iv)
               if (is_external(zdim, LO)) curl%mgvar(:, :, curl%ks-i, iv) = (1+i) * curl%mgvar(:, :, curl%ks, iv) - i * curl%mgvar(:, :, curl%ks+1, iv)
               if (is_external(zdim, HI)) curl%mgvar(:, :, curl%ke+i, iv) = (1+i) * curl%mgvar(:, :, curl%ke, iv) - i * curl%mgvar(:, :, curl%ke-1, iv)
            enddo
         case (extbnd_zero) ! homogenous Dirichlet BC with 0 at first guardcell row
            if (is_external(xdim, LO)) curl%mgvar(:curl%is, :, :, iv) = 0.
            if (is_external(xdim, HI)) curl%mgvar(curl%ie:, :, :, iv) = 0.
            if (is_external(ydim, LO)) curl%mgvar(:, :curl%js, :, iv) = 0.
            if (is_external(ydim, HI)) curl%mgvar(:, curl%je:, :, iv) = 0.
            if (is_external(zdim, LO)) curl%mgvar(:, :, :curl%ks, iv) = 0.
            if (is_external(zdim, HI)) curl%mgvar(:, :, curl%ke:, iv) = 0.
         case (extbnd_mirror) ! reflecting BC (homogenous Neumamnn)
            do i = 1, ng
               if (is_external(xdim, LO)) curl%mgvar(curl%is-i, :, :, iv) = curl%mgvar(curl%is+i-1, :, :, iv)
               if (is_external(xdim, HI)) curl%mgvar(curl%ie+i, :, :, iv) = curl%mgvar(curl%ie-i+1, :, :, iv)
               if (is_external(ydim, LO)) curl%mgvar(:, curl%js-i, :, iv) = curl%mgvar(:, curl%js+i-1, :, iv)
               if (is_external(ydim, HI)) curl%mgvar(:, curl%je+i, :, iv) = curl%mgvar(:, curl%je-i+1, :, iv)
               if (is_external(zdim, LO)) curl%mgvar(:, :, curl%ks-i, iv) = curl%mgvar(:, :, curl%ks+i-1, iv)
               if (is_external(zdim, HI)) curl%mgvar(:, :, curl%ke+i, iv) = curl%mgvar(:, :, curl%ke-i+1, iv)
            enddo
         case (extbnd_antimirror) ! homogenous Dirichlet BC with 0 at external faces
            do i = 1, ng
               if (is_external(xdim, LO)) curl%mgvar(curl%is-i, :, :, iv) = - curl%mgvar(curl%is+i-1, :, :, iv)
               if (is_external(xdim, HI)) curl%mgvar(curl%ie+i, :, :, iv) = - curl%mgvar(curl%ie-i+1, :, :, iv)
               if (is_external(ydim, LO)) curl%mgvar(:, curl%js-i, :, iv) = - curl%mgvar(:, curl%js+i-1, :, iv)
               if (is_external(ydim, HI)) curl%mgvar(:, curl%je+i, :, iv) = - curl%mgvar(:, curl%je-i+1, :, iv)
               if (is_external(zdim, LO)) curl%mgvar(:, :, curl%ks-i, iv) = - curl%mgvar(:, :, curl%ks+i-1, iv)
               if (is_external(zdim, HI)) curl%mgvar(:, :, curl%ke+i, iv) = - curl%mgvar(:, :, curl%ke-i+1, iv)
            enddo
         case default
            write(msg, '(a,i3,a)')"[multigridmpifuncs:multigrid_ext_bnd] boundary type ",mode," not implemented"
            call die(msg)
      end select

   end subroutine multigrid_ext_bnd

end module multigridmpifuncs
