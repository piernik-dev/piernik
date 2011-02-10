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
!> \brief Initialize MPI shortcuts for local domain communication

   subroutine mpi_multigrid_prep

      use multigridvars, only: level_min, level_max, lvl, NDIM
      use mpisetup,      only: ierr, xdim, ydim, zdim, has_dir
      use mpi,           only: MPI_DOUBLE_PRECISION, MPI_ORDER_FORTRAN

      implicit none

      integer                  :: i,ib,lnb
      integer (kind=4)         :: ord, old
      integer, dimension(:), allocatable :: sizes, subsizes, starts

      ord = MPI_ORDER_FORTRAN
      old = MPI_DOUBLE_PRECISION

      allocate( sizes(NDIM), subsizes(NDIM), starts(NDIM) )

      do i = level_min, level_max
         lnb = lvl(i)%nb

         do ib = 1, lnb

            if (has_dir(xdim)) then          !! X direction
               sizes    = [ lvl(i)%nx, lvl(i)%ny, lvl(i)%nz ]
               subsizes = [    ib,     lvl(i)%ny, lvl(i)%nz ]
               starts   = [  lnb-ib,       0,         0     ]

               call MPI_Type_create_subarray(NDIM, sizes, subsizes, starts, ord, old, lvl(i)%MPI_YZ_LEFT_BND(ib), ierr)
               call MPI_Type_commit(lvl(i)%MPI_YZ_LEFT_BND(ib), ierr)

               starts(xdim) = lnb
               call MPI_Type_create_subarray(NDIM, sizes, subsizes,  starts, ord, old, lvl(i)%MPI_YZ_LEFT_DOM(ib), ierr)
               call MPI_Type_commit(lvl(i)%MPI_YZ_LEFT_DOM(ib), ierr)

               starts(xdim) = lvl(i)%nxb + lnb - ib
               call MPI_Type_create_subarray(NDIM, sizes, subsizes,  starts, ord, old, lvl(i)%MPI_YZ_RIGHT_DOM(ib), ierr)
               call MPI_Type_commit(lvl(i)%MPI_YZ_RIGHT_DOM(ib), ierr)

               starts(xdim) = lvl(i)%nxb + lnb
               call MPI_Type_create_subarray(NDIM, sizes, subsizes, starts,  ord, old, lvl(i)%MPI_YZ_RIGHT_BND(ib), ierr)
               call MPI_Type_commit(lvl(i)%MPI_YZ_RIGHT_BND(ib), ierr)
            endif

            if (has_dir(ydim)) then         !! Y Direction
               sizes    = [ lvl(i)%nx, lvl(i)%ny, lvl(i)%nz ]
               subsizes = [ lvl(i)%nx,     ib,    lvl(i)%nz ]
               starts   = [     0,      lnb-ib,       0     ]

               call MPI_Type_create_subarray(NDIM, sizes, subsizes, starts, ord, old, lvl(i)%MPI_XZ_LEFT_BND(ib), ierr)
               call MPI_Type_commit(lvl(i)%MPI_XZ_LEFT_BND(ib), ierr)

               starts(ydim) = lnb
               call MPI_Type_create_subarray(NDIM, sizes, subsizes, starts, ord, old, lvl(i)%MPI_XZ_LEFT_DOM(ib), ierr)
               call MPI_Type_commit(lvl(i)%MPI_XZ_LEFT_DOM(ib), ierr)

               starts(ydim) = lvl(i)%nyb + lnb - ib
               call MPI_Type_create_subarray(NDIM, sizes, subsizes, starts, ord, old, lvl(i)%MPI_XZ_RIGHT_DOM(ib), ierr)
               call MPI_Type_commit(lvl(i)%MPI_XZ_RIGHT_DOM(ib), ierr)

               starts(ydim) = lvl(i)%nyb + lnb
               call MPI_Type_create_subarray(NDIM, sizes, subsizes, starts, ord, old, lvl(i)%MPI_XZ_RIGHT_BND(ib), ierr)
               call MPI_Type_commit(lvl(i)%MPI_XZ_RIGHT_BND(ib), ierr)
            endif

            if (has_dir(zdim)) then         !! Z Direction
               sizes    = [ lvl(i)%nx, lvl(i)%ny, lvl(i)%nz ]
               subsizes = [ lvl(i)%nx, lvl(i)%ny,     ib    ]
               starts   = [     0,         0,       lnb-ib  ]

               call MPI_Type_create_subarray(NDIM, sizes, subsizes, starts, ord, old, lvl(i)%MPI_XY_LEFT_BND(ib), ierr)
               call MPI_Type_commit(lvl(i)%MPI_XY_LEFT_BND(ib), ierr)

               starts(zdim) = lnb
               call MPI_Type_create_subarray(NDIM, sizes, subsizes, starts, ord, old, lvl(i)%MPI_XY_LEFT_DOM(ib), ierr)
               call MPI_Type_commit(lvl(i)%MPI_XY_LEFT_DOM(ib), ierr)

               starts(zdim) = lvl(i)%nzb + lnb - ib
               call MPI_Type_create_subarray(NDIM, sizes, subsizes, starts, ord, old, lvl(i)%MPI_XY_RIGHT_DOM(ib), ierr)
               call MPI_Type_commit(lvl(i)%MPI_XY_RIGHT_DOM(ib), ierr)

               starts(zdim) = lvl(i)%nzb + lnb
               call MPI_Type_create_subarray(NDIM, sizes, subsizes, starts, ord, old, lvl(i)%MPI_XY_RIGHT_BND(ib), ierr)
               call MPI_Type_commit(lvl(i)%MPI_XY_RIGHT_BND(ib), ierr)
            endif

         enddo
      enddo

      deallocate( sizes, subsizes, starts )

   end subroutine mpi_multigrid_prep

!!$ ============================================================================
!!
!! Boundary conditions
!!
!> \brief Routine for inter-process and periodic boundary conditions.
!> \details mpi_multigrid_bnd provides communication between local domains to couple solution on the global computational domain
!!

   subroutine mpi_multigrid_bnd(lev, iv, ng, mode, corners)

      use dataio_pub,         only: die
      use mpisetup,           only: comm3d, ierr, procxl, procxr, procyl, procyr, proczl, proczr, proc, psize, xdim, ydim, zdim, has_dir
      use mpi,                only: MPI_STATUS_SIZE, MPI_REQUEST_NULL
      use multigridvars,      only: NDIM, lvl, XLO, XHI, YLO, YHI, ZLO, ZHI, is_external, ngridvars, level_min, level_max

      implicit none

      integer, intent(in) :: lev               !< level which we are doing communication at
      integer, intent(in) :: iv                !< variable which we want to communicate
      integer, intent(in) :: ng                !< number of guardcells to exchange
      integer, intent(in) :: mode              !< what to do with external boundaries
      logical, intent(in), optional :: corners !< if .true. then don't forget aboutpay close attention to corners

      integer, parameter                        :: nreq = NDIM*4
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
      if (any(is_external(:))) call multigrid_ext_bnd(lev, iv, ng, mode, cor)

      req3d(:) = MPI_REQUEST_NULL

      if (has_dir(xdim)) then
         if (psize(xdim) > 1) then
            if (.not. is_external(XLO)) call MPI_Isend (lvl(lev)%mgvar(1, 1, 1, iv), 1, lvl(lev)%MPI_YZ_LEFT_DOM(ng),  procxl, 15, comm3d, req3d(1),  ierr)
            if (.not. is_external(XHI)) call MPI_Isend (lvl(lev)%mgvar(1, 1, 1, iv), 1, lvl(lev)%MPI_YZ_RIGHT_DOM(ng), procxr, 25, comm3d, req3d(3),  ierr)
            if (.not. is_external(XLO)) call MPI_Irecv (lvl(lev)%mgvar(1, 1, 1, iv), 1, lvl(lev)%MPI_YZ_LEFT_BND(ng),  procxl, 25, comm3d, req3d(2),  ierr)
            if (.not. is_external(XHI)) call MPI_Irecv (lvl(lev)%mgvar(1, 1, 1, iv), 1, lvl(lev)%MPI_YZ_RIGHT_BND(ng), procxr, 15, comm3d, req3d(4),  ierr)
         else
            if (.not. is_external(XLO)) lvl(lev)%mgvar(lvl(lev)%is-ng:lvl(lev)%is-1,  :, :, iv) = lvl(lev)%mgvar(lvl(lev)%ie-ng+1:lvl(lev)%ie,      :, :, iv)
            if (.not. is_external(XHI)) lvl(lev)%mgvar(lvl(lev)%ie+1 :lvl(lev)%ie+ng, :, :, iv) = lvl(lev)%mgvar(lvl(lev)%is     :lvl(lev)%is+ng-1, :, :, iv)
         endif
         if (cor) call MPI_Waitall(4, req3d(1:4), status3d(:,1:4), ierr)
      endif

      if (has_dir(ydim)) then
         if (psize(ydim) > 1) then
            if (.not. is_external(YLO)) call MPI_Isend (lvl(lev)%mgvar(1, 1, 1, iv), 1, lvl(lev)%MPI_XZ_LEFT_DOM(ng),  procyl, 35, comm3d, req3d(5),  ierr)
            if (.not. is_external(YHI)) call MPI_Isend (lvl(lev)%mgvar(1, 1, 1, iv), 1, lvl(lev)%MPI_XZ_RIGHT_DOM(ng), procyr, 45, comm3d, req3d(6),  ierr)
            if (.not. is_external(YLO)) call MPI_Irecv (lvl(lev)%mgvar(1, 1, 1, iv), 1, lvl(lev)%MPI_XZ_LEFT_BND(ng),  procyl, 45, comm3d, req3d(7),  ierr)
            if (.not. is_external(YHI)) call MPI_Irecv (lvl(lev)%mgvar(1, 1, 1, iv), 1, lvl(lev)%MPI_XZ_RIGHT_BND(ng), procyr, 35, comm3d, req3d(8),  ierr)
         else
            if (.not. is_external(YLO)) lvl(lev)%mgvar(:, lvl(lev)%js-ng:lvl(lev)%js-1,  :, iv) = lvl(lev)%mgvar(:, lvl(lev)%je-ng+1:lvl(lev)%je,      :, iv)
            if (.not. is_external(YHI)) lvl(lev)%mgvar(:, lvl(lev)%je+1 :lvl(lev)%je+ng, :, iv) = lvl(lev)%mgvar(:, lvl(lev)%js     :lvl(lev)%js+ng-1, :, iv)
         endif
         if (cor) call MPI_Waitall(4, req3d(5:8), status3d(:,5:8), ierr)
      endif

      if (has_dir(zdim)) then
         if (psize(zdim) > 1) then
            if (.not. is_external(ZLO)) call MPI_Isend (lvl(lev)%mgvar(1, 1, 1, iv), 1, lvl(lev)%MPI_XY_LEFT_DOM(ng),  proczl, 55, comm3d, req3d(9),  ierr)
            if (.not. is_external(ZHI)) call MPI_Isend (lvl(lev)%mgvar(1, 1, 1, iv), 1, lvl(lev)%MPI_XY_RIGHT_DOM(ng), proczr, 65, comm3d, req3d(10), ierr)
            if (.not. is_external(ZLO)) call MPI_Irecv (lvl(lev)%mgvar(1, 1, 1, iv), 1, lvl(lev)%MPI_XY_LEFT_BND(ng),  proczl, 65, comm3d, req3d(11), ierr)
            if (.not. is_external(ZHI)) call MPI_Irecv (lvl(lev)%mgvar(1, 1, 1, iv), 1, lvl(lev)%MPI_XY_RIGHT_BND(ng), proczr, 55, comm3d, req3d(12), ierr)
         else
            if (.not. is_external(ZLO)) lvl(lev)%mgvar(:, :, lvl(lev)%ks-ng:lvl(lev)%ks-1,  iv) = lvl(lev)%mgvar(:, :, lvl(lev)%ke-ng+1:lvl(lev)%ke,      iv)
            if (.not. is_external(ZHI)) lvl(lev)%mgvar(:, :, lvl(lev)%ke+1 :lvl(lev)%ke+ng, iv) = lvl(lev)%mgvar(:, :, lvl(lev)%ks     :lvl(lev)%ks+ng-1, iv)
         endif
         if (cor) call MPI_Waitall(4, req3d(9:12), status3d(:,9:12), ierr)
      endif

!>
!! \todo Make a benchmark of a massively parallel run to determine difference in execution between calling MPI_Waitall for each direction and calling it once.
!! If the difference is small then set cor permanently to .true.
!<
      if (.not. cor) call MPI_Waitall(nreq, req3d(:), status3d(:,:), ierr)

   end subroutine mpi_multigrid_bnd

!!$ ============================================================================
!!
!> \brief Set external boundary (not required for periodic box) on domain faces.
!! \details In multigrid typically mirror boundaries are in use. Extrapolate isolated boundaries at exit.
!!
!! has_dir() is not checked here because is_external() should be set to .false. on non-existing directions in 1D and 2D setups
!<

   subroutine multigrid_ext_bnd(lev, iv, ng, mode, cor)

      use dataio_pub,      only: die, msg, warn
      use multigridvars,   only: extbnd_donothing, extbnd_zero, extbnd_extrapolate, extbnd_mirror, extbnd_antimirror, &
           &                     XLO, XHI, YLO, YHI, ZLO, ZHI, lvl, is_external

      implicit none

      integer, intent(in) :: lev             !< level which we are preparing the guardcells at
      integer, intent(in) :: iv              !< variable which we want to set
      integer, intent(in) :: ng              !< number of guardcells to set
      integer, intent(in) :: mode            !< what to do with external boundaries
      logical, intent(in) :: cor             !< if .true. then don't forget about corners \deprecated BEWARE: not implemented properly

      integer :: i
      logical, save :: warned = .false.

      if (cor) then
         if (.not. warned) then
            call warn("[multigridmpifuncs:multigrid_ext_bnd] some corners may not be properly filled. FIX ME!")
            warned = .true.
         endif
      endif

      select case (mode)         !> \deprecated BEWARE: some cylindrical factors may be helpful
         case (extbnd_donothing) ! remember to initialize everything first!
            return
         case (extbnd_extrapolate) !> \deprecated mixed-tybe BC: free flux; BEWARE: it is not protected from inflow
            do i = 1, ng
               if (is_external(XLO)) lvl(lev)%mgvar(lvl(lev)%is-i, :, :, iv) = (1+i) * lvl(lev)%mgvar(lvl(lev)%is, :, :, iv) - i * lvl(lev)%mgvar(lvl(lev)%is+1, :, :, iv)
               if (is_external(XHI)) lvl(lev)%mgvar(lvl(lev)%ie+i, :, :, iv) = (1+i) * lvl(lev)%mgvar(lvl(lev)%ie, :, :, iv) - i * lvl(lev)%mgvar(lvl(lev)%ie-1, :, :, iv)
               if (is_external(YLO)) lvl(lev)%mgvar(:, lvl(lev)%js-i, :, iv) = (1+i) * lvl(lev)%mgvar(:, lvl(lev)%js, :, iv) - i * lvl(lev)%mgvar(:, lvl(lev)%js+1, :, iv)
               if (is_external(YHI)) lvl(lev)%mgvar(:, lvl(lev)%je+i, :, iv) = (1+i) * lvl(lev)%mgvar(:, lvl(lev)%je, :, iv) - i * lvl(lev)%mgvar(:, lvl(lev)%je-1, :, iv)
               if (is_external(ZLO)) lvl(lev)%mgvar(:, :, lvl(lev)%ks-i, iv) = (1+i) * lvl(lev)%mgvar(:, :, lvl(lev)%ks, iv) - i * lvl(lev)%mgvar(:, :, lvl(lev)%ks+1, iv)
               if (is_external(ZHI)) lvl(lev)%mgvar(:, :, lvl(lev)%ke+i, iv) = (1+i) * lvl(lev)%mgvar(:, :, lvl(lev)%ke, iv) - i * lvl(lev)%mgvar(:, :, lvl(lev)%ke-1, iv)
            enddo
         case (extbnd_zero) ! homogenous Dirichlet BC with 0 at first guardcell row
            if (is_external(XLO)) lvl(lev)%mgvar(:lvl(lev)%is, :, :, iv) = 0.
            if (is_external(XHI)) lvl(lev)%mgvar(lvl(lev)%ie:, :, :, iv) = 0.
            if (is_external(YLO)) lvl(lev)%mgvar(:, :lvl(lev)%js, :, iv) = 0.
            if (is_external(YHI)) lvl(lev)%mgvar(:, lvl(lev)%je:, :, iv) = 0.
            if (is_external(ZLO)) lvl(lev)%mgvar(:, :, :lvl(lev)%ks, iv) = 0.
            if (is_external(ZHI)) lvl(lev)%mgvar(:, :, lvl(lev)%ke:, iv) = 0.
         case (extbnd_mirror) ! reflecting BC (homogenous Neumamnn)
            do i = 1, ng
               if (is_external(XLO)) lvl(lev)%mgvar(lvl(lev)%is-i, :, :, iv) = lvl(lev)%mgvar(lvl(lev)%is+i-1, :, :, iv)
               if (is_external(XHI)) lvl(lev)%mgvar(lvl(lev)%ie+i, :, :, iv) = lvl(lev)%mgvar(lvl(lev)%ie-i+1, :, :, iv)
               if (is_external(YLO)) lvl(lev)%mgvar(:, lvl(lev)%js-i, :, iv) = lvl(lev)%mgvar(:, lvl(lev)%js+i-1, :, iv)
               if (is_external(YHI)) lvl(lev)%mgvar(:, lvl(lev)%je+i, :, iv) = lvl(lev)%mgvar(:, lvl(lev)%je-i+1, :, iv)
               if (is_external(ZLO)) lvl(lev)%mgvar(:, :, lvl(lev)%ks-i, iv) = lvl(lev)%mgvar(:, :, lvl(lev)%ks+i-1, iv)
               if (is_external(ZHI)) lvl(lev)%mgvar(:, :, lvl(lev)%ke+i, iv) = lvl(lev)%mgvar(:, :, lvl(lev)%ke-i+1, iv)
            enddo
          case (extbnd_antimirror) ! homogenous Dirichlet BC with 0 at external faces
            do i = 1, ng
               if (is_external(XLO)) lvl(lev)%mgvar(lvl(lev)%is-i, :, :, iv) = - lvl(lev)%mgvar(lvl(lev)%is+i-1, :, :, iv)
               if (is_external(XHI)) lvl(lev)%mgvar(lvl(lev)%ie+i, :, :, iv) = - lvl(lev)%mgvar(lvl(lev)%ie-i+1, :, :, iv)
               if (is_external(YLO)) lvl(lev)%mgvar(:, lvl(lev)%js-i, :, iv) = - lvl(lev)%mgvar(:, lvl(lev)%js+i-1, :, iv)
               if (is_external(YHI)) lvl(lev)%mgvar(:, lvl(lev)%je+i, :, iv) = - lvl(lev)%mgvar(:, lvl(lev)%je-i+1, :, iv)
               if (is_external(ZLO)) lvl(lev)%mgvar(:, :, lvl(lev)%ks-i, iv) = - lvl(lev)%mgvar(:, :, lvl(lev)%ks+i-1, iv)
               if (is_external(ZHI)) lvl(lev)%mgvar(:, :, lvl(lev)%ke+i, iv) = - lvl(lev)%mgvar(:, :, lvl(lev)%ke-i+1, iv)
            enddo
         case default
            write(msg, '(a,i3,a)')"[multigridmpifuncs:multigrid_ext_bnd] boundary type ",mode," not implemented"
            call die(msg)
      end select

   end subroutine multigrid_ext_bnd

end module multigridmpifuncs
