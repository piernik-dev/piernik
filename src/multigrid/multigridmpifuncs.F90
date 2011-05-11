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

      use constants,     only: xdim, ydim, zdim, LO, HI, BND, BLK, ndims, INVALID
      use dataio_pub,    only: die
      use mpi,           only: MPI_DOUBLE_PRECISION, MPI_ORDER_FORTRAN, MPI_COMM_NULL
      use mpisetup,      only: ierr, has_dir, comm3d, proc
      use multigridvars, only: level_min, level_max, lvl

      implicit none

      integer :: ib, l
      integer, dimension(ndims) :: sizes, subsizes, starts

      ! set up connections between levels
      do l = level_min, level_max

         lvl(l)%i_rst%proc = proc
         lvl(l)%i_rst%se(:,:) = reshape( [ lvl(l)%is, lvl(l)%js, lvl(l)%ks, lvl(l)%ie, lvl(l)%je, lvl(l)%ke ], shape(lvl(l)%i_rst%se(:,:)) )
         lvl(l)%i_rst%nextgrid => null()
         lvl(l)%o_rst = lvl(l)%i_rst

         if (l<level_max) lvl(l)%i_rst%nextgrid => lvl(l+1)
         if (l>level_min) lvl(l)%o_rst%nextgrid => lvl(l-1)

      enddo

      do l = level_min, level_max

         ! find neighbours and set up the MPI containers
         if (comm3d == MPI_COMM_NULL) then

            ! assume that cuboids fill the domain and don't collide

            lvl(l)%mmbc(:, :, :, :) = INVALID

            call die("mmf:imf comm3d == MPI_COMM_NULL not implemented yet")

         else

            ! assume cartesian decomposition

            if (.not. lvl(l)%empty) then

               do ib = 1, lvl(l)%nb

                  !> \todo provide a loop-friendly array lvl(l)%n_b(xdim:zdim) and add a loop "do d = xdim, zdim"
                  if (has_dir(xdim)) then          !! X direction
                     sizes    = [  lvl(l)%nx,   lvl(l)%ny, lvl(l)%nz ]
                     subsizes = [      ib,      lvl(l)%ny, lvl(l)%nz ]
                     starts   = [ lvl(l)%nb-ib,     0,         0     ]

                     call MPI_Type_create_subarray(ndims, sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, lvl(l)%mmbc(xdim, LO, BND, ib), ierr)
                     call MPI_Type_commit(lvl(l)%mmbc(xdim, LO, BND, ib), ierr)

                     starts(xdim) = lvl(l)%nb
                     call MPI_Type_create_subarray(ndims, sizes, subsizes,  starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, lvl(l)%mmbc(xdim, LO, BLK, ib), ierr)
                     call MPI_Type_commit(lvl(l)%mmbc(xdim, LO, BLK, ib), ierr)

                     starts(xdim) = lvl(l)%nxb + lvl(l)%nb - ib
                     call MPI_Type_create_subarray(ndims, sizes, subsizes,  starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, lvl(l)%mmbc(xdim, HI, BLK, ib), ierr)
                     call MPI_Type_commit(lvl(l)%mmbc(xdim, HI, BLK, ib), ierr)

                     starts(xdim) = lvl(l)%nxb + lvl(l)%nb
                     call MPI_Type_create_subarray(ndims, sizes, subsizes, starts,  MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, lvl(l)%mmbc(xdim, HI, BND, ib), ierr)
                     call MPI_Type_commit(lvl(l)%mmbc(xdim, HI, BND, ib), ierr)
                  endif

                  if (has_dir(ydim)) then         !! Y Direction
                     sizes    = [ lvl(l)%nx,  lvl(l)%ny,  lvl(l)%nz ]
                     subsizes = [ lvl(l)%nx,     ib,      lvl(l)%nz ]
                     starts   = [     0,     lvl(l)%nb-ib,    0     ]

                     call MPI_Type_create_subarray(ndims, sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, lvl(l)%mmbc(ydim, LO, BND, ib), ierr)
                     call MPI_Type_commit(lvl(l)%mmbc(ydim, LO, BND, ib), ierr)

                     starts(ydim) = lvl(l)%nb
                     call MPI_Type_create_subarray(ndims, sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, lvl(l)%mmbc(ydim, LO, BLK, ib), ierr)
                     call MPI_Type_commit(lvl(l)%mmbc(ydim, LO, BLK, ib), ierr)

                     starts(ydim) = lvl(l)%nyb + lvl(l)%nb - ib
                     call MPI_Type_create_subarray(ndims, sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, lvl(l)%mmbc(ydim, HI, BLK, ib), ierr)
                     call MPI_Type_commit(lvl(l)%mmbc(ydim, HI, BLK, ib), ierr)

                     starts(ydim) = lvl(l)%nyb + lvl(l)%nb
                     call MPI_Type_create_subarray(ndims, sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, lvl(l)%mmbc(ydim, HI, BND, ib), ierr)
                     call MPI_Type_commit(lvl(l)%mmbc(ydim, HI, BND, ib), ierr)
                  endif

                  if (has_dir(zdim)) then         !! Z Direction
                     sizes    = [ lvl(l)%nx, lvl(l)%ny,  lvl(l)%nz   ]
                     subsizes = [ lvl(l)%nx, lvl(l)%ny,      ib      ]
                     starts   = [     0,         0,     lvl(l)%nb-ib ]

                     call MPI_Type_create_subarray(ndims, sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, lvl(l)%mmbc(zdim, LO, BND, ib), ierr)
                     call MPI_Type_commit(lvl(l)%mmbc(zdim, LO, BND, ib), ierr)

                     starts(zdim) = lvl(l)%nb
                     call MPI_Type_create_subarray(ndims, sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, lvl(l)%mmbc(zdim, LO, BLK, ib), ierr)
                     call MPI_Type_commit(lvl(l)%mmbc(zdim, LO, BLK, ib), ierr)

                     starts(zdim) = lvl(l)%nzb + lvl(l)%nb - ib
                     call MPI_Type_create_subarray(ndims, sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, lvl(l)%mmbc(zdim, HI, BLK, ib), ierr)
                     call MPI_Type_commit(lvl(l)%mmbc(zdim, HI, BLK, ib), ierr)

                     starts(zdim) = lvl(l)%nzb + lvl(l)%nb
                     call MPI_Type_create_subarray(ndims, sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, lvl(l)%mmbc(zdim, HI, BND, ib), ierr)
                     call MPI_Type_commit(lvl(l)%mmbc(zdim, HI, BND, ib), ierr)
                  endif

               enddo
            endif

         endif
      enddo

   end subroutine mpi_multigrid_prep

!!$ ============================================================================
!!
!! Boundary conditions
!!
!> \brief Routine for inter-process and periodic boundary conditions.
!> \details mpi_multigrid_bnd provides communication between local domains to couple solution on the global computational domain
!!

   subroutine mpi_multigrid_bnd(lev, iv, ng, mode, corners)

      use dataio_pub,    only: die
      use mpisetup,      only: comm, comm3d, ierr, procn, proc, psize, has_dir
      use constants,     only: ndims, xdim, ydim, zdim, LO, HI, BND, BLK
      use mpi,           only: MPI_STATUS_SIZE, MPI_REQUEST_NULL, MPI_COMM_NULL
      use multigridvars, only: lvl, is_external, ngridvars, level_min, level_max

      implicit none

      integer, intent(in) :: lev               !< level which we are doing communication at
      integer, intent(in) :: iv                !< variable which we want to communicate
      integer, intent(in) :: ng                !< number of guardcells to exchange
      integer, intent(in) :: mode              !< what to do with external boundaries
      logical, intent(in), optional :: corners !< if .true. then don't forget aboutpay close attention to corners

      integer, parameter                        :: nreq = ndims*4
      integer, dimension(nreq)                  :: req3d
      integer, dimension(MPI_STATUS_SIZE, nreq) :: status3d
      logical                                   :: cor

      if (iv < 1 .or. iv > ngridvars) call die("[multigridmpifuncs:mpi_multigrid_bnd] Invalid variable index.")
      if (lev < level_min .or. lev > level_max) call die("[multigridmpifuncs:mpi_multigrid_bnd] Invalid level number.")

      if (ng > lvl(lev)%nb .or. ng <= 0) call die("[multigridmpifuncs:mpi_multigrid_bnd] Too many or <0 guardcells requested.")

      if (present(corners)) then
         cor = corners
      else
         cor = .false.
      endif

      ! Set the external boundary, where appropriate
      if (any(is_external(:, :))) call multigrid_ext_bnd(lev, iv, ng, mode, cor)

      if (comm3d == MPI_COMM_NULL) then

         call die("mmf:mmb comm3d == MPI_COMM_NULL not implemented yet")

      else

         req3d(:) = MPI_REQUEST_NULL

         if (has_dir(xdim)) then
            if (psize(xdim) > 1) then
               if (.not. is_external(xdim, LO)) call MPI_Isend (lvl(lev)%mgvar(1, 1, 1, iv), 1, lvl(lev)%mmbc(xdim, LO, BLK, ng), procn(xdim, LO), 15, comm3d, req3d(1), ierr)
               if (.not. is_external(xdim, HI)) call MPI_Isend (lvl(lev)%mgvar(1, 1, 1, iv), 1, lvl(lev)%mmbc(xdim, HI, BLK, ng), procn(xdim, HI), 25, comm3d, req3d(3), ierr)
               if (.not. is_external(xdim, LO)) call MPI_Irecv (lvl(lev)%mgvar(1, 1, 1, iv), 1, lvl(lev)%mmbc(xdim, LO, BND, ng), procn(xdim, LO), 25, comm3d, req3d(2), ierr)
               if (.not. is_external(xdim, HI)) call MPI_Irecv (lvl(lev)%mgvar(1, 1, 1, iv), 1, lvl(lev)%mmbc(xdim, HI, BND, ng), procn(xdim, HI), 15, comm3d, req3d(4), ierr)
            else
               if (.not. is_external(xdim, LO)) lvl(lev)%mgvar(lvl(lev)%is-ng:lvl(lev)%is-1,  :, :, iv) = lvl(lev)%mgvar(lvl(lev)%ie-ng+1:lvl(lev)%ie,      :, :, iv)
               if (.not. is_external(xdim, HI)) lvl(lev)%mgvar(lvl(lev)%ie+1 :lvl(lev)%ie+ng, :, :, iv) = lvl(lev)%mgvar(lvl(lev)%is     :lvl(lev)%is+ng-1, :, :, iv)
            endif
            if (cor) call MPI_Waitall(4, req3d(1:4), status3d(:,1:4), ierr)
         endif

         if (has_dir(ydim)) then
            if (psize(ydim) > 1) then
               if (.not. is_external(ydim, LO)) call MPI_Isend (lvl(lev)%mgvar(1, 1, 1, iv), 1, lvl(lev)%mmbc(ydim, LO, BLK, ng), procn(ydim, LO), 35, comm3d, req3d(5), ierr)
               if (.not. is_external(ydim, HI)) call MPI_Isend (lvl(lev)%mgvar(1, 1, 1, iv), 1, lvl(lev)%mmbc(ydim, HI, BLK, ng), procn(ydim, HI), 45, comm3d, req3d(6), ierr)
               if (.not. is_external(ydim, LO)) call MPI_Irecv (lvl(lev)%mgvar(1, 1, 1, iv), 1, lvl(lev)%mmbc(ydim, LO, BND, ng), procn(ydim, LO), 45, comm3d, req3d(7), ierr)
               if (.not. is_external(ydim, HI)) call MPI_Irecv (lvl(lev)%mgvar(1, 1, 1, iv), 1, lvl(lev)%mmbc(ydim, HI, BND, ng), procn(ydim, HI), 35, comm3d, req3d(8), ierr)
            else
               if (.not. is_external(ydim, LO)) lvl(lev)%mgvar(:, lvl(lev)%js-ng:lvl(lev)%js-1,  :, iv) = lvl(lev)%mgvar(:, lvl(lev)%je-ng+1:lvl(lev)%je,      :, iv)
               if (.not. is_external(ydim, HI)) lvl(lev)%mgvar(:, lvl(lev)%je+1 :lvl(lev)%je+ng, :, iv) = lvl(lev)%mgvar(:, lvl(lev)%js     :lvl(lev)%js+ng-1, :, iv)
            endif
            if (cor) call MPI_Waitall(4, req3d(5:8), status3d(:,5:8), ierr)
         endif

         if (has_dir(zdim)) then
            if (psize(zdim) > 1) then
               if (.not. is_external(zdim, LO)) call MPI_Isend (lvl(lev)%mgvar(1, 1, 1, iv), 1, lvl(lev)%mmbc(zdim, LO, BLK, ng), procn(zdim, LO), 55, comm3d, req3d(9), ierr)
               if (.not. is_external(zdim, HI)) call MPI_Isend (lvl(lev)%mgvar(1, 1, 1, iv), 1, lvl(lev)%mmbc(zdim, HI, BLK, ng), procn(zdim, HI), 65, comm3d, req3d(10), ierr)
               if (.not. is_external(zdim, LO)) call MPI_Irecv (lvl(lev)%mgvar(1, 1, 1, iv), 1, lvl(lev)%mmbc(zdim, LO, BND, ng), procn(zdim, LO), 65, comm3d, req3d(11), ierr)
               if (.not. is_external(zdim, HI)) call MPI_Irecv (lvl(lev)%mgvar(1, 1, 1, iv), 1, lvl(lev)%mmbc(zdim, HI, BND, ng), procn(zdim, HI), 55, comm3d, req3d(12), ierr)
            else
               if (.not. is_external(zdim, LO)) lvl(lev)%mgvar(:, :, lvl(lev)%ks-ng:lvl(lev)%ks-1,  iv) = lvl(lev)%mgvar(:, :, lvl(lev)%ke-ng+1:lvl(lev)%ke,      iv)
               if (.not. is_external(zdim, HI)) lvl(lev)%mgvar(:, :, lvl(lev)%ke+1 :lvl(lev)%ke+ng, iv) = lvl(lev)%mgvar(:, :, lvl(lev)%ks     :lvl(lev)%ks+ng-1, iv)
            endif
            if (cor) call MPI_Waitall(4, req3d(9:12), status3d(:,9:12), ierr)
         endif

!>
!! \todo Make a benchmark of a massively parallel run to determine difference in execution between calling MPI_Waitall for each direction and calling it once.
!! If the difference is small then set cor permanently to .true.
!<
         if (.not. cor) call MPI_Waitall(nreq, req3d(:), status3d(:,:), ierr)

      endif

   end subroutine mpi_multigrid_bnd

!!$ ============================================================================
!!
!> \brief Set external boundary (not required for periodic box) on domain faces.
!! \details In multigrid typically mirror boundaries are in use. Extrapolate isolated boundaries at exit.
!!
!! has_dir() is not checked here because is_external() should be set to .false. on non-existing directions in 1D and 2D setups
!<

   subroutine multigrid_ext_bnd(lev, iv, ng, mode, cor)

      use constants,       only: LO, HI, xdim, ydim, zdim
      use dataio_pub,      only: die, msg, warn
      use multigridvars,   only: extbnd_donothing, extbnd_zero, extbnd_extrapolate, extbnd_mirror, extbnd_antimirror, lvl, is_external

      implicit none

      integer, intent(in) :: lev             !< level which we are preparing the guardcells at
      integer, intent(in) :: iv              !< variable which we want to set
      integer, intent(in) :: ng              !< number of guardcells to set
      integer, intent(in) :: mode            !< what to do with external boundaries
      logical, intent(in) :: cor             !< if .true. then don't forget about corners \deprecated BEWARE: not implemented properly

      integer :: i
      logical, save :: warned = .false.

      if (lvl(lev)%empty) return

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
               if (is_external(xdim, LO)) lvl(lev)%mgvar(lvl(lev)%is-i, :, :, iv) = (1+i) * lvl(lev)%mgvar(lvl(lev)%is, :, :, iv) - i * lvl(lev)%mgvar(lvl(lev)%is+1, :, :, iv)
               if (is_external(xdim, HI)) lvl(lev)%mgvar(lvl(lev)%ie+i, :, :, iv) = (1+i) * lvl(lev)%mgvar(lvl(lev)%ie, :, :, iv) - i * lvl(lev)%mgvar(lvl(lev)%ie-1, :, :, iv)
               if (is_external(ydim, LO)) lvl(lev)%mgvar(:, lvl(lev)%js-i, :, iv) = (1+i) * lvl(lev)%mgvar(:, lvl(lev)%js, :, iv) - i * lvl(lev)%mgvar(:, lvl(lev)%js+1, :, iv)
               if (is_external(ydim, HI)) lvl(lev)%mgvar(:, lvl(lev)%je+i, :, iv) = (1+i) * lvl(lev)%mgvar(:, lvl(lev)%je, :, iv) - i * lvl(lev)%mgvar(:, lvl(lev)%je-1, :, iv)
               if (is_external(zdim, LO)) lvl(lev)%mgvar(:, :, lvl(lev)%ks-i, iv) = (1+i) * lvl(lev)%mgvar(:, :, lvl(lev)%ks, iv) - i * lvl(lev)%mgvar(:, :, lvl(lev)%ks+1, iv)
               if (is_external(zdim, HI)) lvl(lev)%mgvar(:, :, lvl(lev)%ke+i, iv) = (1+i) * lvl(lev)%mgvar(:, :, lvl(lev)%ke, iv) - i * lvl(lev)%mgvar(:, :, lvl(lev)%ke-1, iv)
            enddo
         case (extbnd_zero) ! homogenous Dirichlet BC with 0 at first guardcell row
            if (is_external(xdim, LO)) lvl(lev)%mgvar(:lvl(lev)%is, :, :, iv) = 0.
            if (is_external(xdim, HI)) lvl(lev)%mgvar(lvl(lev)%ie:, :, :, iv) = 0.
            if (is_external(ydim, LO)) lvl(lev)%mgvar(:, :lvl(lev)%js, :, iv) = 0.
            if (is_external(ydim, HI)) lvl(lev)%mgvar(:, lvl(lev)%je:, :, iv) = 0.
            if (is_external(zdim, LO)) lvl(lev)%mgvar(:, :, :lvl(lev)%ks, iv) = 0.
            if (is_external(zdim, HI)) lvl(lev)%mgvar(:, :, lvl(lev)%ke:, iv) = 0.
         case (extbnd_mirror) ! reflecting BC (homogenous Neumamnn)
            do i = 1, ng
               if (is_external(xdim, LO)) lvl(lev)%mgvar(lvl(lev)%is-i, :, :, iv) = lvl(lev)%mgvar(lvl(lev)%is+i-1, :, :, iv)
               if (is_external(xdim, HI)) lvl(lev)%mgvar(lvl(lev)%ie+i, :, :, iv) = lvl(lev)%mgvar(lvl(lev)%ie-i+1, :, :, iv)
               if (is_external(ydim, LO)) lvl(lev)%mgvar(:, lvl(lev)%js-i, :, iv) = lvl(lev)%mgvar(:, lvl(lev)%js+i-1, :, iv)
               if (is_external(ydim, HI)) lvl(lev)%mgvar(:, lvl(lev)%je+i, :, iv) = lvl(lev)%mgvar(:, lvl(lev)%je-i+1, :, iv)
               if (is_external(zdim, LO)) lvl(lev)%mgvar(:, :, lvl(lev)%ks-i, iv) = lvl(lev)%mgvar(:, :, lvl(lev)%ks+i-1, iv)
               if (is_external(zdim, HI)) lvl(lev)%mgvar(:, :, lvl(lev)%ke+i, iv) = lvl(lev)%mgvar(:, :, lvl(lev)%ke-i+1, iv)
            enddo
         case (extbnd_antimirror) ! homogenous Dirichlet BC with 0 at external faces
            do i = 1, ng
               if (is_external(xdim, LO)) lvl(lev)%mgvar(lvl(lev)%is-i, :, :, iv) = - lvl(lev)%mgvar(lvl(lev)%is+i-1, :, :, iv)
               if (is_external(xdim, HI)) lvl(lev)%mgvar(lvl(lev)%ie+i, :, :, iv) = - lvl(lev)%mgvar(lvl(lev)%ie-i+1, :, :, iv)
               if (is_external(ydim, LO)) lvl(lev)%mgvar(:, lvl(lev)%js-i, :, iv) = - lvl(lev)%mgvar(:, lvl(lev)%js+i-1, :, iv)
               if (is_external(ydim, HI)) lvl(lev)%mgvar(:, lvl(lev)%je+i, :, iv) = - lvl(lev)%mgvar(:, lvl(lev)%je-i+1, :, iv)
               if (is_external(zdim, LO)) lvl(lev)%mgvar(:, :, lvl(lev)%ks-i, iv) = - lvl(lev)%mgvar(:, :, lvl(lev)%ks+i-1, iv)
               if (is_external(zdim, HI)) lvl(lev)%mgvar(:, :, lvl(lev)%ke+i, iv) = - lvl(lev)%mgvar(:, :, lvl(lev)%ke-i+1, iv)
            enddo
         case default
            write(msg, '(a,i3,a)')"[multigridmpifuncs:multigrid_ext_bnd] boundary type ",mode," not implemented"
            call die(msg)
      end select

   end subroutine multigrid_ext_bnd

end module multigridmpifuncs
