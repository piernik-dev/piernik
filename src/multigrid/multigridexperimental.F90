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

!!$ ============================================================================
!>
!! \brief This module contains experimental routines, not recommended for production runs.
!! \details
!! <table border="1" cellpadding="4" cellspacing="0">
!!   <tr><td> Cell-face prolongation stencils for fast convergence on uniform grid </td>
!!       <td> -1./12. </td><td> 7./12. </td><td> 7./12. </td><td> -1./12. </td><td> integral cubic </td></tr>
!!   <tr><td> Slightly slower convergence, less wide stencil  </td>
!!       <td>         </td><td> 1./2.  </td><td> 1./2.  </td><td>         </td><td> average; integral and direct linear </td></tr>
!! </table>
!!\n Prolongation of cell faces from cell centers are required for FFT local solver, red-black Gauss-Seidel relaxation don't use it.
!!
!!\n Cell-centered prolongation stencils, for odd fine cells, for even fine cells reverse the order.
!! <table border="1" cellpadding="4" cellspacing="0">
!!   <tr><td> 35./2048. </td><td> -252./2048. </td><td> 1890./2048. </td><td> 420./2048. </td><td> -45./2048. </td><td> direct quartic </td></tr>
!!   <tr><td>           </td><td>   -7./128.  </td><td>  105./128.  </td><td>  35./128.  </td><td>  -5./128.  </td><td> direct cubic </td></tr>
!!   <tr><td>           </td><td>   -3./32.   </td><td>   30./32.   </td><td>   5./32.   </td><td>            </td><td> direct quadratic </td></tr>
!!   <tr><td>           </td><td>             </td><td>    1.       </td><td>            </td><td>            </td><td> injection (0-th order), direct and integral approach </td></tr>
!!   <tr><td>           </td><td>             </td><td>    3./4.    </td><td>   1./4.    </td><td>            </td><td> linear, direct and integral approach </td></tr>
!!   <tr><td>           </td><td>   -1./8.    </td><td>    1.       </td><td>   1./8.    </td><td>            </td><td> integral quadratic </td></tr>
!!   <tr><td>           </td><td>   -5./64.   </td><td>   55./64.   </td><td>  17./64.   </td><td>  -3./64.   </td><td> integral cubic </td></tr>
!!   <tr><td>   3./128. </td><td>  -11./64.   </td><td>    1.       </td><td>  11./64.   </td><td>  -3./128.  </td><td> integral quartic </td></tr>
!! </table>
!!
!!\n General rule is that the second big positive coefficient should be assigned to closer neighbor of the coarse parent cell.
!!\n Thus a single coarse contributes to fine cells in the following way:
!! <table border="1" cellpadding="4" cellspacing="0">
!!   <tr><td> fine level   </td>
!!       <td> -3./128. </td><td> 3./128. </td><td> -11./64. </td><td>  11./64. </td><td> 1. </td><td> 1. </td>
!!       <td> 11./64. </td><td> -11./64. </td><td> 3./128.  </td><td> -3./128. </td><td> integral quartic coefficients </td></tr>
!!   <tr><td> coarse level </td>
!!       <td colspan="2">                </td><td colspan="2">                 </td><td colspan="2"> 1.  </td>
!!       <td colspan="2">                </td><td colspan="2">                 </td><td>                               </td></tr>
!! </table>
!!\n
!!\n The term "n-th order integral interpolation" here means that the prolonged values satisfy the following condition:
!!\n Integral over a cell of a n-th order polynomial fit to the nearest 5 points in each dimension on coarse level
!! is equal to the sum of similar integrals over fine cells covering the coarse cell.
!!\n
!!\n The term "n-th order direct interpolation" here means that the prolonged values are n-th order polynomial fit
!! to the nearest 5 points in each dimension on coarse level evaluated for fine cell centers.
!!\n
!!\n It seems that for 3D Cartesian grid with isolated boundaries and relaxation implemented in approximate_solution
!! pure injection gives highest norm reduction factors per V-cycle. For other boundary types, FFT implementation of
!! approximate_solution, specific source distribution or some other multigrid scheme may give faster convergence than injection.
!!\n
!!\n Estimated prolongation costs for integral quartic stencil:
!!\n  "gather" approach: loop over fine cells, each one collects weighted values from 5**3 coarse cells (125*n_fine multiplications
!!\n  "scatter" approach: loop over coarse cells, each one contributes weighted values to 10**3 fine cells (1000*n_coarse multiplications, roughly equal to cost of gather)
!!\n  "directionally split" approach: do the prolongation (either gather or scatter type) first in x direction (10*n_coarse multiplications -> 2*n_coarse intermediate cells
!!                                  result), then in y direction (10*2*n_coarse multiplications -> 4*n_coarse intermediate cells result), then in z direction
!!                                  (10*4*n_coarse multiplications -> 8*n_coarse = n_fine cells result). Looks like 70*n_coarse multiplications.
!!                                  Will require two additional arrays for intermediate results.
!!\n  "FFT" approach: do the convolution in Fourier space. Unfortunately it is not periodic box, so we cannot use power of 2 FFT sizes. No idea how fast or slow this can be.
!!\n
!!\n For AMR or nested grid low-order prolongation schemes (injection and linear interpolation at least) are known to produce artifacts
!! on fine-coarse boundaries. For uniform grid the simplest operators are probably the fastest and give best V-cycle convergence rates.
!<

module multigridexperimental
! pulled by MULTIGRID

   implicit none

   private
   public :: prolong_level_hord

contains

!!$ ============================================================================
!!
!> \brief high-order order prolongation interpolation
!!

   subroutine prolong_level_hord(coarse, iv)

      use constants,         only: ndims
      use dataio_pub,        only: die, warn, msg
      use domain,            only: eff_dim
      use mpisetup,          only: master
      use multigridmpifuncs, only: mpi_multigrid_bnd
      use multigridvars,     only: ord_prolong, extbnd_antimirror, plvl

      implicit none

      type(plvl), pointer, intent(in) :: coarse   !< level to prolong from
      integer(kind=4), intent(in) :: iv    !< variable to be prolonged

      logical, save :: firstcall = .true.

      if (firstcall) then
         if (master) then
            write(msg,'(a,i3,a)')"[multigridexperimental:prolong_level_hord] prolongation order ",ord_prolong," is experimental"
            call warn(msg)
         endif
         firstcall = .false.
      endif

      call mpi_multigrid_bnd(coarse, iv, abs(ord_prolong/2), extbnd_antimirror) ! exchange guardcells with corners

      if (eff_dim<ndims) call die("[multigridexperimental:prolong_level_hord] 1D and 2D not finished")

      select case (ord_prolong)
      case (-4)
         call prolong_level4D(coarse, iv)
      case (-2)
         call prolong_level2D(coarse, iv)
      case (2)
         call prolong_level2I(coarse, iv)
      case (4)
         call prolong_level4I(coarse, iv)
      case default
         call die("[multigridexperimental:prolong_level_hord] Unsupported 'ord_prolong' value")
      end select

   end subroutine prolong_level_hord

!!$ ============================================================================
!>
!! \brief 2nd order interpolation, integral version
!<

   subroutine prolong_level2I(coarse, iv)

      use dataio_pub,    only: die
      use multigridvars, only: plvl, lvl

      implicit none

      type(plvl), pointer, intent(in) :: coarse  !< level to prolong from
      integer(kind=4), intent(in) :: iv    !< variable to be prolonged

      type(plvl), pointer :: fine

      real, parameter :: P0 = 1., P1 = 1./8.

      if (.not. associated(coarse)) call die("[multigridexperimental:prolong_level2I] coarse == null()")
      fine   => coarse%finer
      if (.not. associated(fine)) call die("[multigridexperimental:prolong_level2I] fine == null()")

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
!>
!! \brief 2nd order interpolation, direct version
!<

   subroutine prolong_level2D(coarse, iv)

      use dataio_pub,    only: die
      use multigridvars, only: plvl, lvl

      implicit none

      type(plvl), pointer, intent(in) :: coarse  !< level to prolong from
      integer(kind=4), intent(in) :: iv    !< variable to be prolonged

      type(plvl), pointer :: fine
      real, parameter :: P_1 = -3./32., P0 = 30./32., P1 = 5./32.

      if (.not. associated(coarse)) call die("[multigridexperimental:prolong_level2D] coarse == null()")
      fine   => coarse%finer
      if (.not. associated(fine)) call die("[multigridexperimental:prolong_level2D] fine == null()")

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
!>
!! \brief 4th order interpolation, integral version
!<

   subroutine prolong_level4I(coarse, iv)

      use dataio_pub,    only: die
      use multigridvars, only: plvl, lvl

      implicit none

      type(plvl), pointer, intent(in) :: coarse  !< level to prolong from
      integer(kind=4), intent(in) :: iv    !< variable to be prolonged

      type(plvl), pointer :: fine

      real, parameter :: P0 = 1., P1 = 11./64., P2 = 3./128.

      if (.not. associated(coarse)) call die("[multigridexperimental:prolong_level4I] coarse == null()")
      fine   => coarse%finer
      if (.not. associated(fine)) call die("[multigridexperimental:prolong_level4I] fine == null()")

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
!>
!! \brief 4th order interpolation, direct version
!<

   subroutine prolong_level4D(coarse, iv)

      use dataio_pub,    only: die
      use multigridvars, only: plvl, lvl

      implicit none

      type(plvl), pointer, intent(in) :: coarse  !< level to prolong from
      integer(kind=4), intent(in) :: iv    !< variable to be prolonged

      type(plvl), pointer :: fine

      real, parameter :: P_2 = 35./2048., P_1 = -252./2048., P0 = 1890./2048., P1 = 420./2048., P2 = -45./2048.

      if (.not. associated(coarse)) call die("[multigridexperimental:prolong_level4D] coarse == null()")
      fine   => coarse%finer
      if (.not. associated(fine)) call die("[multigridexperimental:prolong_level4D] fine == null()")

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

end module multigridexperimental
