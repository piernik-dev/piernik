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
!! \brief This module contains various routines (interpolation, boundaries and some global reduction)
!! that are useful for all flavours of multigrid solvers.
!<

module multigridbasefuncs
! pulled by MULTIGRID
   implicit none

   private

   public :: prolong_level

contains

!!$ ============================================================================
!>
!! \brief Multigrid elementary operators: prolongation, restriction, norm etc.
!<

   subroutine prolong_level(coarse, iv)

      use cg_list_lev,           only: cg_list_level
      use constants,             only: O_INJ
      use dataio_pub,            only: die
      use multigridhelpers,      only: dirty_debug, check_dirty, dirtyH
      use multigridvars,         only: roof, ord_prolong, is_mg_uneven

      implicit none

      type(cg_list_level), pointer, intent(in) :: coarse !< level to prolong from
      integer,                      intent(in) :: iv     !< variable to be prolonged

      type(cg_list_level), pointer :: fine

      if (associated(coarse, roof)) return ! can't prolong finest level

      if (.not. associated(coarse)) call die("[multigridbasefuncs:prolong_level] coarse == null()")
      fine   => coarse%finer
      if (.not. associated(fine)) call die("[multigridbasefuncs:prolong_level] fine == null()")

      !! \warning: need set_dirty_level routine
      if (dirty_debug) fine%first%cg%q(iv)%arr(:, :, :) = dirtyH

      call check_dirty(coarse, iv, "prolong-")

      if (ord_prolong == O_INJ .or. is_mg_uneven) then
         call coarse%prolong0_q_1var(iv)
      else
         call prolong_level_hord(coarse, iv) ! experimental part
      endif

      call check_dirty(fine, iv, "prolong+")

   end subroutine prolong_level

!>
!! \brief high-order order prolongation interpolation
!!
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
!! direct quadratic and cubic interpolations give best norm reduction factors per V-cycle (maclaurin problem).
!!  For other boundary types, FFT implementation of approximate_solution, specific source distribution or
!!  some other multigrid scheme may give faster convergence rate.
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
!!
!! \deprecated These routines assume simplest domain decomposition where fine grids cover exactly the same area as coarse grids
!!
!! \todo This will finally belong to cg_list_level type
!<

   subroutine prolong_level_hord(coarse, iv)

      use cg_list_lev,       only: cg_list_level
      use constants,         only: ndims, O_INJ, O_LIN, O_D2, O_D3, O_D4, O_I2, O_I3, O_I4, BND_NEGREF
      use dataio_pub,        only: die, warn, msg
      use domain,            only: dom, is_multicg
      use grid_cont,         only: grid_container
      use mpisetup,          only: master
      use multigridmpifuncs, only: mpi_multigrid_bnd
      use multigridvars,     only: ord_prolong

      implicit none

      type(cg_list_level), pointer, intent(in) :: coarse !< level to prolong from
      integer,                      intent(in) :: iv     !< variable to be prolonged

      logical, save :: firstcall = .true.
      real :: P_2, P_1, P0, P1, P2
      type(grid_container), pointer :: cg_f, cg_c

      if (firstcall) then
         if (master) then
            write(msg,'(a,i3,a)')"[multigridbasefuncs:prolong_level_hord] prolongation order ",ord_prolong," is experimental"
            call warn(msg)
         endif
         firstcall = .false.
      endif
      if (abs(ord_prolong) > 2*dom%nb) call die("[multigridbasefuncs:prolong_level_hord] not enough guardcells for given prolongation operator order")

      call mpi_multigrid_bnd(coarse, iv, int((abs(ord_prolong)+1)/2, kind=4), BND_NEGREF) ! exchange guardcells with corners

      if (dom%eff_dim<ndims) call die("[multigridbasefuncs:prolong_level_hord] 1D and 2D not finished")
      if (is_multicg) call die("[multigridbasefuncs:prolong_level_hord] multicg not implemented")

      select case (ord_prolong)
         case (O_D4)
            P_2 = 35./2048.; P_1 = -252./2048.; P0 = 1890./2048.; P1 = 420./2048.; P2 = -45./2048.
         case (O_D3)
            P_2 = 0.;        P_1 = -7./128.;    P0 = 105./128.;   P1 = 35./128.;   P2 = -5./128.
         case (O_D2)
            P_2 = 0.;        P_1 = -3./32.;     P0 = 30./32.;     P1 = 5./32.;     P2 = 0.
         case (O_LIN)
            P_2 = 0.;        P_1 = 0.;          P0 = 3./4.;       P1 = 1./4.;      P2 = 0.
         case (O_INJ)
            P_2 = 0.;        P_1 = 0.;          P0 = 1.;          P1 = 0.;         P2 = 0.
         case (O_I2)
            P_2 = 0.;        P_1 = -1./8.;      P0 = 1.;          P1 = 1./8.;      P2 = 0.
         case (O_I3)
            P_2 = 0.;        P_1 = -5./64.;     P0 = 55./64;      P1 = 17./64.;    P2 = -3./64.
         case (O_I4)
            P_2 = 3./128.;   P_1 = -11./64.;    P0 = 1.;          P1 = 11./64.;    P2 = -3./128.
            ! case 0 is handled through cg_list_level%prolong0_q_1var
         case default
            call die("[multigridbasefuncs:prolong_level_hord] Unsupported order")
            return
      end select

      ! this design doesn't allow multiple blocks or inter-process fine-coarse communication
      if (.not. associated(coarse)) call die("[multigridbasefuncs:prolong_level_hord] coarse == null()")
      if (.not. associated(coarse%finer)) call die("[multigridbasefuncs:prolong_level_hord] fine == null()")
      cg_c => coarse%first%cg
      cg_f => coarse%finer%first%cg

      ! convolve with the prolongation operator
      ! the two cases here are for optimization (see also call mpi_multigrid_bnd above)
      if (P_2 == 0. .and. P2 == 0.) then
         cg_f%prolong_x(         cg_f%is  :cg_f%ie-1:2, cg_c%js-1:cg_c%je+1, cg_c%ks-1:cg_c%ke+1) = &  ! x-odd cells
              + P1 * cg_c%q(iv)%arr(cg_c%is-1:cg_c%ie-1,   cg_c%js-1:cg_c%je+1, cg_c%ks-1:cg_c%ke+1)   &
              + P0 * cg_c%q(iv)%arr(cg_c%is  :cg_c%ie,     cg_c%js-1:cg_c%je+1, cg_c%ks-1:cg_c%ke+1)   &
              + P_1* cg_c%q(iv)%arr(cg_c%is+1:cg_c%ie+1,   cg_c%js-1:cg_c%je+1, cg_c%ks-1:cg_c%ke+1)
         cg_f%prolong_x(         cg_f%is+1:cg_f%ie:2,   cg_c%js-1:cg_c%je+1, cg_c%ks-1:cg_c%ke+1) = &  ! x-even cells
              + P_1* cg_c%q(iv)%arr(cg_c%is-1:cg_c%ie-1,   cg_c%js-1:cg_c%je+1, cg_c%ks-1:cg_c%ke+1)   &
              + P0 * cg_c%q(iv)%arr(cg_c%is  :cg_c%ie,     cg_c%js-1:cg_c%je+1, cg_c%ks-1:cg_c%ke+1)   &
              + P1 * cg_c%q(iv)%arr(cg_c%is+1:cg_c%ie+1,   cg_c%js-1:cg_c%je+1, cg_c%ks-1:cg_c%ke+1)

         cg_f%prolong_xy(           cg_f%is:cg_f%ie, cg_f%js  :cg_f%je-1:2, cg_c%ks-1:cg_c%ke+1) = &    ! y-odd cells
              + P1 * cg_f%prolong_x(cg_f%is:cg_f%ie, cg_c%js-1:cg_c%je-1,   cg_c%ks-1:cg_c%ke+1)   &
              + P0 * cg_f%prolong_x(cg_f%is:cg_f%ie, cg_c%js  :cg_c%je,     cg_c%ks-1:cg_c%ke+1)   &
              + P_1* cg_f%prolong_x(cg_f%is:cg_f%ie, cg_c%js+1:cg_c%je+1,   cg_c%ks-1:cg_c%ke+1)
         cg_f%prolong_xy(           cg_f%is:cg_f%ie, cg_f%js+1:cg_f%je:2,   cg_c%ks-1:cg_c%ke+1) = &    ! y-even cells
              + P_1* cg_f%prolong_x(cg_f%is:cg_f%ie, cg_c%js-1:cg_c%je-1,   cg_c%ks-1:cg_c%ke+1)   &
              + P0 * cg_f%prolong_x(cg_f%is:cg_f%ie, cg_c%js  :cg_c%je,     cg_c%ks-1:cg_c%ke+1)   &
              + P1 * cg_f%prolong_x(cg_f%is:cg_f%ie, cg_c%js+1:cg_c%je+1,   cg_c%ks-1:cg_c%ke+1)

         cg_f%q(iv)%arr(                cg_f%is:cg_f%ie, cg_f%js:cg_f%je, cg_f%ks  :cg_f%ke-1:2) = &   ! z-odd cells
              + P1 * cg_f%prolong_xy(cg_f%is:cg_f%ie, cg_f%js:cg_f%je, cg_c%ks-1:cg_c%ke-1  )   &
              + P0 * cg_f%prolong_xy(cg_f%is:cg_f%ie, cg_f%js:cg_f%je, cg_c%ks  :cg_c%ke    )   &
              + P_1* cg_f%prolong_xy(cg_f%is:cg_f%ie, cg_f%js:cg_f%je, cg_c%ks+1:cg_c%ke+1  )
         cg_f%q(iv)%arr(                cg_f%is:cg_f%ie, cg_f%js:cg_f%je, cg_f%ks+1:cg_f%ke:2  ) = &   ! z-even cells
              + P_1* cg_f%prolong_xy(cg_f%is:cg_f%ie, cg_f%js:cg_f%je, cg_c%ks-1:cg_c%ke-1  )   &
              + P0 * cg_f%prolong_xy(cg_f%is:cg_f%ie, cg_f%js:cg_f%je, cg_c%ks  :cg_c%ke    )   &
              + P1 * cg_f%prolong_xy(cg_f%is:cg_f%ie, cg_f%js:cg_f%je, cg_c%ks+1:cg_c%ke+1  )
      else
         cg_f%prolong_x(         cg_f%is  :cg_f%ie-1:2, cg_c%js-2:cg_c%je+2, cg_c%ks-2:cg_c%ke+2) = &  ! x-odd cells
              + P2 * cg_c%q(iv)%arr(cg_c%is-2:cg_c%ie-2,   cg_c%js-2:cg_c%je+2, cg_c%ks-2:cg_c%ke+2)   &
              + P1 * cg_c%q(iv)%arr(cg_c%is-1:cg_c%ie-1,   cg_c%js-2:cg_c%je+2, cg_c%ks-2:cg_c%ke+2)   &
              + P0 * cg_c%q(iv)%arr(cg_c%is  :cg_c%ie,     cg_c%js-2:cg_c%je+2, cg_c%ks-2:cg_c%ke+2)   &
              + P_1* cg_c%q(iv)%arr(cg_c%is+1:cg_c%ie+1,   cg_c%js-2:cg_c%je+2, cg_c%ks-2:cg_c%ke+2)   &
              + P_2* cg_c%q(iv)%arr(cg_c%is+2:cg_c%ie+2,   cg_c%js-2:cg_c%je+2, cg_c%ks-2:cg_c%ke+2)
         cg_f%prolong_x(         cg_f%is+1:cg_f%ie:2,   cg_c%js-2:cg_c%je+2, cg_c%ks-2:cg_c%ke+2) = &  ! x-even cells
              + P_2* cg_c%q(iv)%arr(cg_c%is-2:cg_c%ie-2,   cg_c%js-2:cg_c%je+2, cg_c%ks-2:cg_c%ke+2)   &
              + P_1* cg_c%q(iv)%arr(cg_c%is-1:cg_c%ie-1,   cg_c%js-2:cg_c%je+2, cg_c%ks-2:cg_c%ke+2)   &
              + P0 * cg_c%q(iv)%arr(cg_c%is  :cg_c%ie,     cg_c%js-2:cg_c%je+2, cg_c%ks-2:cg_c%ke+2)   &
              + P1 * cg_c%q(iv)%arr(cg_c%is+1:cg_c%ie+1,   cg_c%js-2:cg_c%je+2, cg_c%ks-2:cg_c%ke+2)   &
              + P2 * cg_c%q(iv)%arr(cg_c%is+2:cg_c%ie+2,   cg_c%js-2:cg_c%je+2, cg_c%ks-2:cg_c%ke+2)

         cg_f%prolong_xy(           cg_f%is:cg_f%ie, cg_f%js  :cg_f%je-1:2, cg_c%ks-2:cg_c%ke+2) = &    ! y-odd cells
              + P2 * cg_f%prolong_x(cg_f%is:cg_f%ie, cg_c%js-2:cg_c%je-2,   cg_c%ks-2:cg_c%ke+2)   &
              + P1 * cg_f%prolong_x(cg_f%is:cg_f%ie, cg_c%js-1:cg_c%je-1,   cg_c%ks-2:cg_c%ke+2)   &
              + P0 * cg_f%prolong_x(cg_f%is:cg_f%ie, cg_c%js  :cg_c%je,     cg_c%ks-2:cg_c%ke+2)   &
              + P_1* cg_f%prolong_x(cg_f%is:cg_f%ie, cg_c%js+1:cg_c%je+1,   cg_c%ks-2:cg_c%ke+2)   &
              + P_2* cg_f%prolong_x(cg_f%is:cg_f%ie, cg_c%js+2:cg_c%je+2,   cg_c%ks-2:cg_c%ke+2)
         cg_f%prolong_xy(           cg_f%is:cg_f%ie, cg_f%js+1:cg_f%je:2,   cg_c%ks-2:cg_c%ke+2) = &    ! y-even cells
              + P_2* cg_f%prolong_x(cg_f%is:cg_f%ie, cg_c%js-2:cg_c%je-2,   cg_c%ks-2:cg_c%ke+2)   &
              + P_1* cg_f%prolong_x(cg_f%is:cg_f%ie, cg_c%js-1:cg_c%je-1,   cg_c%ks-2:cg_c%ke+2)   &
              + P0 * cg_f%prolong_x(cg_f%is:cg_f%ie, cg_c%js  :cg_c%je,     cg_c%ks-2:cg_c%ke+2)   &
              + P1 * cg_f%prolong_x(cg_f%is:cg_f%ie, cg_c%js+1:cg_c%je+1,   cg_c%ks-2:cg_c%ke+2)   &
              + P2 * cg_f%prolong_x(cg_f%is:cg_f%ie, cg_c%js+2:cg_c%je+2,   cg_c%ks-2:cg_c%ke+2)

         cg_f%q(iv)%arr(                cg_f%is:cg_f%ie, cg_f%js:cg_f%je, cg_f%ks  :cg_f%ke-1:2) = &   ! z-odd cells
              + P2 * cg_f%prolong_xy(cg_f%is:cg_f%ie, cg_f%js:cg_f%je, cg_c%ks-2:cg_c%ke-2  )   &
              + P1 * cg_f%prolong_xy(cg_f%is:cg_f%ie, cg_f%js:cg_f%je, cg_c%ks-1:cg_c%ke-1  )   &
              + P0 * cg_f%prolong_xy(cg_f%is:cg_f%ie, cg_f%js:cg_f%je, cg_c%ks  :cg_c%ke    )   &
              + P_1* cg_f%prolong_xy(cg_f%is:cg_f%ie, cg_f%js:cg_f%je, cg_c%ks+1:cg_c%ke+1  )   &
              + P_2* cg_f%prolong_xy(cg_f%is:cg_f%ie, cg_f%js:cg_f%je, cg_c%ks+2:cg_c%ke+2  )
         cg_f%q(iv)%arr(                cg_f%is:cg_f%ie, cg_f%js:cg_f%je, cg_f%ks+1:cg_f%ke:2  ) = &   ! z-even cells
              + P_2* cg_f%prolong_xy(cg_f%is:cg_f%ie, cg_f%js:cg_f%je, cg_c%ks-2:cg_c%ke-2  )   &
              + P_1* cg_f%prolong_xy(cg_f%is:cg_f%ie, cg_f%js:cg_f%je, cg_c%ks-1:cg_c%ke-1  )   &
              + P0 * cg_f%prolong_xy(cg_f%is:cg_f%ie, cg_f%js:cg_f%je, cg_c%ks  :cg_c%ke    )   &
              + P1 * cg_f%prolong_xy(cg_f%is:cg_f%ie, cg_f%js:cg_f%je, cg_c%ks+1:cg_c%ke+1  )   &
              + P2 * cg_f%prolong_xy(cg_f%is:cg_f%ie, cg_f%js:cg_f%je, cg_c%ks+2:cg_c%ke+2  )
      endif

   end subroutine prolong_level_hord

end module multigridbasefuncs
