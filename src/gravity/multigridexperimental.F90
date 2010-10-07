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
!    Initial implemetation of PIERNIK code was based on TVD split MHD code by
!    Ue-Li Pen
!        see: Pen, Arras & Wong (2003) for algorithm and
!             http://www.cita.utoronto.ca/~pen/MHD
!             for original source code "mhd.f90"
!
!    For full list of developers see $PIERNIK_HOME/license/pdt.txt
!

!!$ ============================================================================
!
! Cell-face prolongation stencils for fast convergence on uniform grid:
!   -1./12., 7./12., 7./12., -1./12. (integral qubic)
! slightly slower, less wide stencil:
!            1./2.,  1./2.           (average; integral&direct linear)
!// Prolongation of cell faces from cell centers will be required for FFT local solver, red-black Gauss-Seidel relaxation don't use it.

!!$ ============================================================================
!
! Cell-centered prolongation stencils, for odd fine cells, for even fine cells reverse the order.
!  35./2048., -252./2048., 1890./2048., 420./2048., -45./2048. ; direct quartic
!               -7./128.,   105./128.,   35./128.,   -5./128.  ; direct cubic
!               -3./32.,     30./32.,     5./32.               ; direct quadratic
!                             1.                               ; injection (0-th order), same for direct and integral approach
!                             3./4.       1./4.                ; linear, same for direct and integral approach
!               -1./8.,       1.,         1./8.,               ; integral quadratic
!               -5./64.,     55./64.,    17./64.,    -3./64.   ; integral cubic
!   3./128.,   -11./64.,      1.,        11./64.,    -3./128.  ; integral quartic
!
! Genaral rule is that  so the second big positive coefficient will be assigned to closer neighbor of the coarse parent cell.
! Thus a single coarse contributes to fine cells in the following way:
!
! |          |        |          |        |       |      |         |         |         |        | fine level
! | -3./128., 3./128.,| -11./64., 11./64.,| 1.,    1.,   | 11./64., -11./64.,| 3./128., -3./128.| integral quartic coefficients
! |                   |                   |              |                   |                  | coarse level
!
! The term "n-th order integral interpolation" here means that the prolonged values satisfy the following condition:
! Integral over a cell of a n-th order polynomial fit to the nearest 5 points in each dimension on coarse level
! is equal to the sum of similar integrals over fine cells covering the coarse cell.
!
! The term "n-th order direct interpolation" here means that the prolonged values are n-th order polynomial fit
! to the nearest 5 points in each dimension on coarse level evaluated for fine cell centers.
!
! It seems that for 3D cartesian grid with isolated boundaries and relaxation implemented in approximate_solution
! pure injection gives highest norm reduction factors per V-cycle. For other boundary types, FFT implementation of
! approximate_solution, specific source distribution or some other multigrid scheme may give faster convergence than injection.
!
! Estimated prolongation costs for integral quartic stencil:
!  "gather" approach: loop over fine cells, each one collects weighted values from 5**3 coarse cells (125*n_fine multiplications
!  "scatter" approach: loop over coarse cells, each one contributes weighted values to 10**3 fine cells (1000*n_coarse multiplications, roughly equal to cost of gather)
!  "directionally split" approach: do the prolongation (either gather or scatter type) first in x direction (10*n_coarse multiplications -> 2*n_coarse intermediate cells
!                                  result), then in y direction (10*2*n_coarse multiplications -> 4*n_coarse intermediate cells result), then in z direction
!                                  (10*4*n_coarse multiplications -> 8*n_coarse = n_fine cells result). Looks like 70*n_coarse multiplications.
!                                  Will require two additional arrays for intermediate results.
!  "FFT" approach: do the convolution in Fourier space. Unfortunately it is not periodic box, so we cannot use power of 2 FFT sizes. No idea how fast or slow this can be.
!
! For AMR or nested grid low-order prolongation schemes (injection and linear interpolation at least) are known to produce artifacts
! on fine-coarse boundaries. For uniform grid the simplest operators are probably the fastest and give best V-cycle convergence rates.
!

#include "piernik.def"

!!$ ============================================================================
!!
!! This module contains experimental routines, not recommended for production runs.
!!

!module multigridexperimental

!contains

   subroutine prolong_level_hord(lev, iv)

      use mpisetup,           only: proc
      use errh,               only: die, warn
      use multigridvars,      only: ord_prolong
      use multigridmpifuncs,  only: mpi_multigrid_bnd
      use dataio_public,      only: msg

      implicit none

      integer, intent(in)      :: lev   !< level to prolong from
      integer, intent(in)      :: iv    !< variable to be prolonged

      logical, save :: firstcall = .true.

      if (firstcall) then
         if (proc == 0) then
            write(msg,'(a,i3,a)')"[multigridexperimental:prolong_level_hord] prolongation order ",ord_prolong," is experimental"
            call warn(msg)
         end if
         firstcall = .false.
      end if

      call mpi_multigrid_bnd(lev, iv, abs(ord_prolong/2), .false.) ! exchange guardcells with corners

      select case (ord_prolong)
      case (-4)
         call prolong_level4D(lev, iv)
      case (-2)
         call prolong_level2D(lev, iv)
      case (2)
         call prolong_level2I(lev, iv)
      case (4)
         call prolong_level4I(lev, iv)
      case default
         call die("[multigridexperimental:prolong_level_hord] Unsupported 'ord_prolong' value")
      end select

   end subroutine prolong_level_hord

!!$ ============================================================================
!!
!! 2nd order interpolation, integral version
!!

   subroutine prolong_level2I(lev, iv)

      use multigridvars, only: plvl, lvl, eff_dim, NDIM
      use errh,          only: die

      implicit none

      integer, intent(in)      :: lev   !< level to prolong from
      integer, intent(in)      :: iv    !< variable to be prolonged

      type(plvl), pointer :: coarse, fine

      real, parameter :: P0 = 1., P1 = 1./8.

      if (eff_dim<NDIM) call die("[multigrid:prolong_level2I] 1D and 2D not finished")

      coarse => lvl(lev)
      fine   => lvl(lev + 1)

      ! convolve with the prolongation operator
      fine%prolong_x(          fine%is    :fine%ie-1:2, coarse%js-1:coarse%je+1, coarse%ks-1:coarse%ke+1) = &  ! x-odd cells
           + P1 * coarse%mgvar(coarse%is-1:coarse%ie-1, coarse%js-1:coarse%je+1, coarse%ks-1:coarse%ke+1, iv) &
           + P0 * coarse%mgvar(coarse%is  :coarse%ie  , coarse%js-1:coarse%je+1, coarse%ks-1:coarse%ke+1, iv) &
           - P1 * coarse%mgvar(coarse%is+1:coarse%ie+1, coarse%js-1:coarse%je+1, coarse%ks-1:coarse%ke+1, iv)
      fine%prolong_x(          fine%is+1  :fine%ie:2,   coarse%js-1:coarse%je+1, coarse%ks-1:coarse%ke+1) = &  ! x-even cells
           - P1 * coarse%mgvar(coarse%is-1:coarse%ie-1, coarse%js-1:coarse%je+1, coarse%ks-1:coarse%ke+1, iv) &
           + P0 * coarse%mgvar(coarse%is  :coarse%ie+0, coarse%js-1:coarse%je+1, coarse%ks-1:coarse%ke+1, iv) &
           + P1 * coarse%mgvar(coarse%is+1:coarse%ie+1, coarse%js-1:coarse%je+1, coarse%ks-1:coarse%ke+1, iv)

      fine%prolong_xy(           fine%is:fine%ie, fine%js    :fine%je-1:2, coarse%ks-1:coarse%ke+1) = &    ! y-odd cells
           + P1 * fine%prolong_x(fine%is:fine%ie, coarse%js-1:coarse%je-1, coarse%ks-1:coarse%ke+1) &
           + P0 * fine%prolong_x(fine%is:fine%ie, coarse%js  :coarse%je  , coarse%ks-1:coarse%ke+1) &
           - P1 * fine%prolong_x(fine%is:fine%ie, coarse%js+1:coarse%je+1, coarse%ks-1:coarse%ke+1)
      fine%prolong_xy(           fine%is:fine%ie, fine%js+1  :fine%je:2,   coarse%ks-1:coarse%ke+1) = &    ! y-even cells
           - P1 * fine%prolong_x(fine%is:fine%ie, coarse%js-1:coarse%je-1, coarse%ks-1:coarse%ke+1) &
           + P0 * fine%prolong_x(fine%is:fine%ie, coarse%js  :coarse%je  , coarse%ks-1:coarse%ke+1) &
           + P1 * fine%prolong_x(fine%is:fine%ie, coarse%js+1:coarse%je+1, coarse%ks-1:coarse%ke+1)

      fine%mgvar(                 fine%is:fine%ie, fine%js:fine%je, fine%ks    :fine%ke-1:2, iv) = &   ! z-odd cells
           + P1 * fine%prolong_xy(fine%is:fine%ie, fine%js:fine%je, coarse%ks-1:coarse%ke-1) &
           + P0 * fine%prolong_xy(fine%is:fine%ie, fine%js:fine%je, coarse%ks  :coarse%ke  ) &
           - P1 * fine%prolong_xy(fine%is:fine%ie, fine%js:fine%je, coarse%ks+1:coarse%ke+1)
      fine%mgvar(                 fine%is:fine%ie, fine%js:fine%je, fine%ks+1  :fine%ke:2,   iv) = &   ! z-even cells
           - P1 * fine%prolong_xy(fine%is:fine%ie, fine%js:fine%je, coarse%ks-1:coarse%ke-1) &
           + P0 * fine%prolong_xy(fine%is:fine%ie, fine%js:fine%je, coarse%ks  :coarse%ke  ) &
           + P1 * fine%prolong_xy(fine%is:fine%ie, fine%js:fine%je, coarse%ks+1:coarse%ke+1)

   end subroutine prolong_level2I

!!$ ============================================================================
!!
!! 2nd order interpolation, direct version
!!

   subroutine prolong_level2D(lev, iv)

      use multigridvars, only: plvl, lvl, eff_dim, NDIM
      use errh,          only: die

      implicit none

      integer, intent(in)      :: lev   !< level to prolong from
      integer, intent(in)      :: iv    !< variable to be prolonged

      type(plvl), pointer :: coarse, fine
      real, parameter :: P_1 = -3./32., P0 = 30./32., P1 = 5./32.

      if (eff_dim<NDIM) call die("[multigrid:prolong_level2D] 1D and 2D not finished")

      coarse => lvl(lev)
      fine   => lvl(lev + 1)

      ! convolve with the prolongation operator
      fine%prolong_x(          fine%is    :fine%ie-1:2, coarse%js-1:coarse%je+1, coarse%ks-1:coarse%ke+1) = &  ! x-odd cells
           + P1 * coarse%mgvar(coarse%is-1:coarse%ie-1, coarse%js-1:coarse%je+1, coarse%ks-1:coarse%ke+1, iv) &
           + P0 * coarse%mgvar(coarse%is  :coarse%ie  , coarse%js-1:coarse%je+1, coarse%ks-1:coarse%ke+1, iv) &
           + P_1* coarse%mgvar(coarse%is+1:coarse%ie+1, coarse%js-1:coarse%je+1, coarse%ks-1:coarse%ke+1, iv)
      fine%prolong_x(          fine%is+1  :fine%ie:2,   coarse%js-1:coarse%je+1, coarse%ks-1:coarse%ke+1) = &  ! x-even cells
           + P_1* coarse%mgvar(coarse%is-1:coarse%ie-1, coarse%js-1:coarse%je+1, coarse%ks-1:coarse%ke+1, iv) &
           + P0 * coarse%mgvar(coarse%is  :coarse%ie+0, coarse%js-1:coarse%je+1, coarse%ks-1:coarse%ke+1, iv) &
           + P1 * coarse%mgvar(coarse%is+1:coarse%ie+1, coarse%js-1:coarse%je+1, coarse%ks-1:coarse%ke+1, iv)

      fine%prolong_xy(           fine%is:fine%ie, fine%js    :fine%je-1:2, coarse%ks-1:coarse%ke+1) = &    ! y-odd cells
           + P1 * fine%prolong_x(fine%is:fine%ie, coarse%js-1:coarse%je-1, coarse%ks-1:coarse%ke+1) &
           + P0 * fine%prolong_x(fine%is:fine%ie, coarse%js  :coarse%je  , coarse%ks-1:coarse%ke+1) &
           + P_1* fine%prolong_x(fine%is:fine%ie, coarse%js+1:coarse%je+1, coarse%ks-1:coarse%ke+1)
      fine%prolong_xy(           fine%is:fine%ie, fine%js+1  :fine%je:2,   coarse%ks-1:coarse%ke+1) = &    ! y-even cells
           + P_1* fine%prolong_x(fine%is:fine%ie, coarse%js-1:coarse%je-1, coarse%ks-1:coarse%ke+1) &
           + P0 * fine%prolong_x(fine%is:fine%ie, coarse%js  :coarse%je  , coarse%ks-1:coarse%ke+1) &
           + P1 * fine%prolong_x(fine%is:fine%ie, coarse%js+1:coarse%je+1, coarse%ks-1:coarse%ke+1)

      fine%mgvar(                 fine%is:fine%ie, fine%js:fine%je, fine%ks    :fine%ke-1:2, iv) = &   ! z-odd cells
           + P1 * fine%prolong_xy(fine%is:fine%ie, fine%js:fine%je, coarse%ks-1:coarse%ke-1) &
           + P0 * fine%prolong_xy(fine%is:fine%ie, fine%js:fine%je, coarse%ks  :coarse%ke  ) &
           + P_1* fine%prolong_xy(fine%is:fine%ie, fine%js:fine%je, coarse%ks+1:coarse%ke+1)
      fine%mgvar(                 fine%is:fine%ie, fine%js:fine%je, fine%ks+1  :fine%ke:2,   iv) = &   ! z-even cells
           + P_1* fine%prolong_xy(fine%is:fine%ie, fine%js:fine%je, coarse%ks-1:coarse%ke-1) &
           + P0 * fine%prolong_xy(fine%is:fine%ie, fine%js:fine%je, coarse%ks  :coarse%ke  ) &
           + P1 * fine%prolong_xy(fine%is:fine%ie, fine%js:fine%je, coarse%ks+1:coarse%ke+1)

   end subroutine prolong_level2D

!!$ ============================================================================
!!
!! 4th order interpolation, integral version
!!

   subroutine prolong_level4I(lev, iv)

      use multigridvars, only: plvl, lvl, eff_dim, NDIM
      use errh,          only: die

      implicit none

      integer, intent(in)      :: lev   !< level to prolong from
      integer, intent(in)      :: iv    !< variable to be prolonged

      type(plvl), pointer :: coarse, fine

      real, parameter :: P0 = 1., P1 = 11./64., P2 = 3./128.

      if (eff_dim<NDIM) call die("[multigrid:prolong_level4I] 1D and 2D not finished")

      coarse => lvl(lev)
      fine   => lvl(lev + 1)

      ! convolve with the prolongation operator
      fine%prolong_x(          fine%is    :fine%ie-1:2, coarse%js-2:coarse%je+2, coarse%ks-2:coarse%ke+2) = &  ! x-odd cells
           - P2 * coarse%mgvar(coarse%is-2:coarse%ie-2, coarse%js-2:coarse%je+2, coarse%ks-2:coarse%ke+2, iv) &
           + P1 * coarse%mgvar(coarse%is-1:coarse%ie-1, coarse%js-2:coarse%je+2, coarse%ks-2:coarse%ke+2, iv) &
           + P0 * coarse%mgvar(coarse%is  :coarse%ie  , coarse%js-2:coarse%je+2, coarse%ks-2:coarse%ke+2, iv) &
           - P1 * coarse%mgvar(coarse%is+1:coarse%ie+1, coarse%js-2:coarse%je+2, coarse%ks-2:coarse%ke+2, iv) &
           + P2 * coarse%mgvar(coarse%is+2:coarse%ie+2, coarse%js-2:coarse%je+2, coarse%ks-2:coarse%ke+2, iv)
      fine%prolong_x(          fine%is+1  :fine%ie:2,   coarse%js-2:coarse%je+2, coarse%ks-2:coarse%ke+2) = &  ! x-even cells
           + P2 * coarse%mgvar(coarse%is-2:coarse%ie-2, coarse%js-2:coarse%je+2, coarse%ks-2:coarse%ke+2, iv) &
           - P1 * coarse%mgvar(coarse%is-1:coarse%ie-1, coarse%js-2:coarse%je+2, coarse%ks-2:coarse%ke+2, iv) &
           + P0 * coarse%mgvar(coarse%is  :coarse%ie+0, coarse%js-2:coarse%je+2, coarse%ks-2:coarse%ke+2, iv) &
           + P1 * coarse%mgvar(coarse%is+1:coarse%ie+1, coarse%js-2:coarse%je+2, coarse%ks-2:coarse%ke+2, iv) &
           - P2 * coarse%mgvar(coarse%is+2:coarse%ie+2, coarse%js-2:coarse%je+2, coarse%ks-2:coarse%ke+2, iv)

      fine%prolong_xy(           fine%is:fine%ie, fine%js    :fine%je-1:2, coarse%ks-2:coarse%ke+2) = &    ! y-odd cells
           - P2 * fine%prolong_x(fine%is:fine%ie, coarse%js-2:coarse%je-2, coarse%ks-2:coarse%ke+2) &
           + P1 * fine%prolong_x(fine%is:fine%ie, coarse%js-1:coarse%je-1, coarse%ks-2:coarse%ke+2) &
           + P0 * fine%prolong_x(fine%is:fine%ie, coarse%js  :coarse%je  , coarse%ks-2:coarse%ke+2) &
           - P1 * fine%prolong_x(fine%is:fine%ie, coarse%js+1:coarse%je+1, coarse%ks-2:coarse%ke+2) &
           + P2 * fine%prolong_x(fine%is:fine%ie, coarse%js+2:coarse%je+2, coarse%ks-2:coarse%ke+2)
      fine%prolong_xy(           fine%is:fine%ie, fine%js+1  :fine%je:2,   coarse%ks-2:coarse%ke+2) = &    ! y-even cells
           + P2 * fine%prolong_x(fine%is:fine%ie, coarse%js-2:coarse%je-2, coarse%ks-2:coarse%ke+2) &
           - P1 * fine%prolong_x(fine%is:fine%ie, coarse%js-1:coarse%je-1, coarse%ks-2:coarse%ke+2) &
           + P0 * fine%prolong_x(fine%is:fine%ie, coarse%js  :coarse%je  , coarse%ks-2:coarse%ke+2) &
           + P1 * fine%prolong_x(fine%is:fine%ie, coarse%js+1:coarse%je+1, coarse%ks-2:coarse%ke+2) &
           - P2 * fine%prolong_x(fine%is:fine%ie, coarse%js+2:coarse%je+2, coarse%ks-2:coarse%ke+2)

      fine%mgvar(                 fine%is:fine%ie, fine%js:fine%je, fine%ks    :fine%ke-1:2, iv) = &   ! z-odd cells
           - P2 * fine%prolong_xy(fine%is:fine%ie, fine%js:fine%je, coarse%ks-2:coarse%ke-2) &
           + P1 * fine%prolong_xy(fine%is:fine%ie, fine%js:fine%je, coarse%ks-1:coarse%ke-1) &
           + P0 * fine%prolong_xy(fine%is:fine%ie, fine%js:fine%je, coarse%ks  :coarse%ke  ) &
           - P1 * fine%prolong_xy(fine%is:fine%ie, fine%js:fine%je, coarse%ks+1:coarse%ke+1) &
           + P2 * fine%prolong_xy(fine%is:fine%ie, fine%js:fine%je, coarse%ks+2:coarse%ke+2)
      fine%mgvar(                 fine%is:fine%ie, fine%js:fine%je, fine%ks+1  :fine%ke:2,   iv) = &   ! z-even cells
           + P2 * fine%prolong_xy(fine%is:fine%ie, fine%js:fine%je, coarse%ks-2:coarse%ke-2) &
           - P1 * fine%prolong_xy(fine%is:fine%ie, fine%js:fine%je, coarse%ks-1:coarse%ke-1) &
           + P0 * fine%prolong_xy(fine%is:fine%ie, fine%js:fine%je, coarse%ks  :coarse%ke  ) &
           + P1 * fine%prolong_xy(fine%is:fine%ie, fine%js:fine%je, coarse%ks+1:coarse%ke+1) &
           - P2 * fine%prolong_xy(fine%is:fine%ie, fine%js:fine%je, coarse%ks+2:coarse%ke+2)

   end subroutine prolong_level4I

!!$ ============================================================================
!!
!! 4th order interpolation, direct version
!!

   subroutine prolong_level4D(lev, iv)

      use multigridvars, only: plvl, lvl, eff_dim, NDIM
      use errh,          only: die

      implicit none

      integer, intent(in)      :: lev   !< level to prolong from
      integer, intent(in)      :: iv    !< variable to be prolonged

      type(plvl), pointer :: coarse, fine

      real, parameter :: P_2 = 35./2048., P_1 = -252./2048., P0 = 1890./2048., P1 = 420./2048., P2 = -45./2048.

      if (eff_dim<NDIM) call die("[multigrid:prolong_level4D] 1D and 2D not finished")

      coarse => lvl(lev)
      fine   => lvl(lev + 1)

      ! convolve with the prolongation operator
      fine%prolong_x(          fine%is    :fine%ie-1:2, coarse%js-2:coarse%je+2, coarse%ks-2:coarse%ke+2) = &  ! x-odd cells
           + P2 * coarse%mgvar(coarse%is-2:coarse%ie-2, coarse%js-2:coarse%je+2, coarse%ks-2:coarse%ke+2, iv) &
           + P1 * coarse%mgvar(coarse%is-1:coarse%ie-1, coarse%js-2:coarse%je+2, coarse%ks-2:coarse%ke+2, iv) &
           + P0 * coarse%mgvar(coarse%is  :coarse%ie  , coarse%js-2:coarse%je+2, coarse%ks-2:coarse%ke+2, iv) &
           + P_1* coarse%mgvar(coarse%is+1:coarse%ie+1, coarse%js-2:coarse%je+2, coarse%ks-2:coarse%ke+2, iv) &
           + P_2* coarse%mgvar(coarse%is+2:coarse%ie+2, coarse%js-2:coarse%je+2, coarse%ks-2:coarse%ke+2, iv)
      fine%prolong_x(          fine%is+1  :fine%ie:2,   coarse%js-2:coarse%je+2, coarse%ks-2:coarse%ke+2) = &  ! x-even cells
           + P_2* coarse%mgvar(coarse%is-2:coarse%ie-2, coarse%js-2:coarse%je+2, coarse%ks-2:coarse%ke+2, iv) &
           + P_1* coarse%mgvar(coarse%is-1:coarse%ie-1, coarse%js-2:coarse%je+2, coarse%ks-2:coarse%ke+2, iv) &
           + P0 * coarse%mgvar(coarse%is  :coarse%ie+0, coarse%js-2:coarse%je+2, coarse%ks-2:coarse%ke+2, iv) &
           + P1 * coarse%mgvar(coarse%is+1:coarse%ie+1, coarse%js-2:coarse%je+2, coarse%ks-2:coarse%ke+2, iv) &
           + P2 * coarse%mgvar(coarse%is+2:coarse%ie+2, coarse%js-2:coarse%je+2, coarse%ks-2:coarse%ke+2, iv)

      fine%prolong_xy(           fine%is:fine%ie, fine%js    :fine%je-1:2, coarse%ks-2:coarse%ke+2) = &    ! y-odd cells
           + P2 * fine%prolong_x(fine%is:fine%ie, coarse%js-2:coarse%je-2, coarse%ks-2:coarse%ke+2) &
           + P1 * fine%prolong_x(fine%is:fine%ie, coarse%js-1:coarse%je-1, coarse%ks-2:coarse%ke+2) &
           + P0 * fine%prolong_x(fine%is:fine%ie, coarse%js  :coarse%je  , coarse%ks-2:coarse%ke+2) &
           + P_1* fine%prolong_x(fine%is:fine%ie, coarse%js+1:coarse%je+1, coarse%ks-2:coarse%ke+2) &
           + P_2* fine%prolong_x(fine%is:fine%ie, coarse%js+2:coarse%je+2, coarse%ks-2:coarse%ke+2)
      fine%prolong_xy(           fine%is:fine%ie, fine%js+1  :fine%je:2,   coarse%ks-2:coarse%ke+2) = &    ! y-even cells
           + P_2* fine%prolong_x(fine%is:fine%ie, coarse%js-2:coarse%je-2, coarse%ks-2:coarse%ke+2) &
           + P_1* fine%prolong_x(fine%is:fine%ie, coarse%js-1:coarse%je-1, coarse%ks-2:coarse%ke+2) &
           + P0 * fine%prolong_x(fine%is:fine%ie, coarse%js  :coarse%je  , coarse%ks-2:coarse%ke+2) &
           + P1 * fine%prolong_x(fine%is:fine%ie, coarse%js+1:coarse%je+1, coarse%ks-2:coarse%ke+2) &
           + P2 * fine%prolong_x(fine%is:fine%ie, coarse%js+2:coarse%je+2, coarse%ks-2:coarse%ke+2)

      fine%mgvar(                 fine%is:fine%ie, fine%js:fine%je, fine%ks    :fine%ke-1:2, iv) = &   ! z-odd cells
           + P2 * fine%prolong_xy(fine%is:fine%ie, fine%js:fine%je, coarse%ks-2:coarse%ke-2) &
           + P1 * fine%prolong_xy(fine%is:fine%ie, fine%js:fine%je, coarse%ks-1:coarse%ke-1) &
           + P0 * fine%prolong_xy(fine%is:fine%ie, fine%js:fine%je, coarse%ks  :coarse%ke  ) &
           + P_1* fine%prolong_xy(fine%is:fine%ie, fine%js:fine%je, coarse%ks+1:coarse%ke+1) &
           + P_2* fine%prolong_xy(fine%is:fine%ie, fine%js:fine%je, coarse%ks+2:coarse%ke+2)
      fine%mgvar(                 fine%is:fine%ie, fine%js:fine%je, fine%ks+1  :fine%ke:2,   iv) = &   ! z-even cells
           + P_2* fine%prolong_xy(fine%is:fine%ie, fine%js:fine%je, coarse%ks-2:coarse%ke-2) &
           + P_1* fine%prolong_xy(fine%is:fine%ie, fine%js:fine%je, coarse%ks-1:coarse%ke-1) &
           + P0 * fine%prolong_xy(fine%is:fine%ie, fine%js:fine%je, coarse%ks  :coarse%ke  ) &
           + P1 * fine%prolong_xy(fine%is:fine%ie, fine%js:fine%je, coarse%ks+1:coarse%ke+1) &
           + P2 * fine%prolong_xy(fine%is:fine%ie, fine%js:fine%je, coarse%ks+2:coarse%ke+2)

   end subroutine prolong_level4D

!!$ ============================================================================
!!
!! 4th order Laplacian
!!
!! Significantly slows down convergence, does not seem to improve quality of solution in simple tests.
!!
!! L4 = [0, 1, -2, 1, 0] + L4_strength * 1./12. * [ -1, 4, -6, 4, -1 ]
!! For integrated face fluxes in the 4th order Laplacian estimate set L4_strength = 0.5
!! For simple 5-point L4 set L4_strength = 1.0
!!
!! There also exists more compact Mehrstellen scheme.
!!

   subroutine residual4(lev, src, soln, def)

      use mpisetup,           only: proc
      use multigridvars,      only: lvl, eff_dim, NDIM, L4_strength, grav_bnd, bnd_givenval
      use multigridmpifuncs,  only: mpi_multigrid_bnd
      use errh,               only: die, warn

      implicit none

      integer, intent(in) :: lev  !< level for which approximate the solution
      integer, intent(in) :: src  !< index of source in lvl()%mgvar
      integer, intent(in) :: soln !< index of solution in lvl()%mgvar
      integer, intent(in) :: def  !< index of defect in lvl()%mgvar

      real, parameter     :: L4_scaling = 1./12. ! with L4_strength = 1. this gives an L4 approximation for finite differences approach
      integer, parameter  :: L2w = 2             ! #layers of boundary cells for L2 operator

      real                :: c21, c41, c42 !, c20, c40
      real                :: L0, Lx1, Lx2, Ly1, Ly2, Lz1, Lz2, Lx, Ly, Lz

      logical, save       :: firstcall = .true.
      integer             :: i, j, k

      if (eff_dim<NDIM) call die("[multigrid:residual4] Only 3D is implemented")

      if (firstcall) then
         if (proc == 0) call warn("[multigridexperimental:residual4] residual order 4 is experimental.")
         firstcall = .false.
      end if

      call mpi_multigrid_bnd(lev, soln, 2, .false.) ! no corners required

      c21 = 1.
      c42 = - L4_scaling * L4_strength
      c41 = c21 + 4. * L4_scaling * L4_strength
      !c20 = -2.
      !c40 = c20 - 6. * L4_strength

      Lx1 = c41 * lvl(lev)%idx2
      Ly1 = c41 * lvl(lev)%idy2
      Lz1 = c41 * lvl(lev)%idz2
      Lx2 = c42 * lvl(lev)%idx2
      Ly2 = c42 * lvl(lev)%idy2
      Lz2 = c42 * lvl(lev)%idz2
!      L0  = c40 * (lvl(lev)%idx2 + lvl(lev)%idy2 + lvl(lev)%idz2 )
      L0 = -2. * (Lx1 + Lx2 + Ly1 + Ly2 + Lz1 + Lz2)

      lvl(     lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   lvl(lev)%ks  :lvl(lev)%ke,   def)        = &
           lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   lvl(lev)%ks  :lvl(lev)%ke,   src)        - &
           lvl(lev)%mgvar(lvl(lev)%is-2:lvl(lev)%ie-2, lvl(lev)%js  :lvl(lev)%je,   lvl(lev)%ks  :lvl(lev)%ke,   soln) * Lx2 - &
           lvl(lev)%mgvar(lvl(lev)%is+2:lvl(lev)%ie+2, lvl(lev)%js  :lvl(lev)%je,   lvl(lev)%ks  :lvl(lev)%ke,   soln) * Lx2 - &
           lvl(lev)%mgvar(lvl(lev)%is-1:lvl(lev)%ie-1, lvl(lev)%js  :lvl(lev)%je,   lvl(lev)%ks  :lvl(lev)%ke,   soln) * Lx1 - &
           lvl(lev)%mgvar(lvl(lev)%is+1:lvl(lev)%ie+1, lvl(lev)%js  :lvl(lev)%je,   lvl(lev)%ks  :lvl(lev)%ke,   soln) * Lx1 - &
           lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js-2:lvl(lev)%je-2, lvl(lev)%ks  :lvl(lev)%ke,   soln) * Ly2 - &
           lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js+2:lvl(lev)%je+2, lvl(lev)%ks  :lvl(lev)%ke,   soln) * Ly2 - &
           lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js-1:lvl(lev)%je-1, lvl(lev)%ks  :lvl(lev)%ke,   soln) * Ly1 - &
           lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js+1:lvl(lev)%je+1, lvl(lev)%ks  :lvl(lev)%ke,   soln) * Ly1 - &
           lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   lvl(lev)%ks-2:lvl(lev)%ke-2, soln) * Lz2 - &
           lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   lvl(lev)%ks+2:lvl(lev)%ke+2, soln) * Lz2 - &
           lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   lvl(lev)%ks-1:lvl(lev)%ke-1, soln) * Lz1 - &
           lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   lvl(lev)%ks+1:lvl(lev)%ke+1, soln) * Lz1 - &
           lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   lvl(lev)%ks  :lvl(lev)%ke,   soln) * L0

      ! WARNING: not optimized
      if (grav_bnd == bnd_givenval) then ! probably also in some other cases
         ! Use L2 Laplacian in two layers of cells next to the boundary because L4 seems to be incompatible with present image mass construction
         Lx = c21 * lvl(lev)%idx2
         Ly = c21 * lvl(lev)%idy2
         Lz = c21 * lvl(lev)%idz2
         L0 = -2. * (Lx + Ly + Lz)

         do k = lvl(lev)%ks, lvl(lev)%ke
            do j = lvl(lev)%js, lvl(lev)%je
               do i = lvl(lev)%is, lvl(lev)%ie
                  if ( i<lvl(lev)%is+L2w .or. i>lvl(lev)%ie-L2w .or. j<lvl(lev)%js+L2w .or. j>lvl(lev)%je-L2w .or. k<lvl(lev)%ks+L2w .or. k>lvl(lev)%ke-L2w) then
                     lvl(       lev)%mgvar(i,   j,   k,   def)   = lvl(lev)%mgvar(i,   j,   k,   src)        - &
                          ( lvl(lev)%mgvar(i-1, j,   k,   soln)  + lvl(lev)%mgvar(i+1, j,   k,   soln)) * Lx - &
                          ( lvl(lev)%mgvar(i,   j-1, k,   soln)  + lvl(lev)%mgvar(i,   j+1, k,   soln)) * Ly - &
                          ( lvl(lev)%mgvar(i,   j,   k-1, soln)  + lvl(lev)%mgvar(i,   j,   k+1, soln)) * Lz - &
                          & lvl(lev)%mgvar(i,   j,   k,   soln)  * L0
                  end if
               end do
            end do
         end do
      end if

   end subroutine residual4

!end module multigridexperimental
