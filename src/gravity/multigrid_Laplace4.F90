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
#include "macros.h"

!> \brief Implementation of 4th order, non-compact (5 cell wide) Laplace operator

module multigrid_Laplace4
! pulled by MULTIGRID && GRAV

   implicit none

   private
   public :: residual4, approximate_solution_rbgs4

contains

!>
!! \brief 4th order Laplacian
!!
!! \details Significantly slows down convergence, does not seem to improve quality of solution in simple tests.
!!
!! L4 = [0, 1, -2, 1, 0] + L4_strength * 1./12. * [ -1, 4, -6, 4, -1 ] = 1./12. * [ -1, 16, -30, 16, -1 ]
!! For integrated face fluxes in the 4th order Laplacian estimate set L4_strength = 0.5
!! For simple 5-point L4 set L4_strength = 1.0
!!
!! There also exists more compact Mehrstellen 4th order scheme which outperforms this one in terms of accuracy
!!
!! As this operator is outperformed by the 4th order Mehrstellen operator, we don't plan any improvements here.
!! It can be either deleted or left for curious people.
!!
!! For curious people:
!! L6 = L4(L4_strength == 1.0) + L6_strength * 1./120. * [ 1, -6, 15, -20, 15, -6, 1 ] = 1./120. * [ 1, -16, 175, -320, 175, -16, 1 ]
!! For simple 7-point L6 assume L6_strength = 1.0
!! Don't expect much improvement by implementing any 6th order Laplacian unless other things in the code are similarly high order.
!<

   subroutine residual4(cg_llst, src, soln, def)

      use cg_list,       only: cg_list_element
      use cg_leaves,     only: cg_leaves_T
      use constants,     only: ndims, idm2, xdim, ydim, zdim, BND_NEGREF
      use dataio_pub,    only: die, warn
      use domain,        only: dom
      use grid_cont,     only: grid_container
      use mpisetup,      only: master
      use multigridvars, only: grav_bnd, bnd_givenval, L4_strength
      use named_array,   only: p3

      implicit none

      class(cg_leaves_T), intent(in) :: cg_llst !< pointer to a level for which we approximate the solution
      integer(kind=4),    intent(in) :: src     !< index of source in cg%q(:)
      integer(kind=4),    intent(in) :: soln    !< index of solution in cg%q(:)
      integer(kind=4),    intent(in) :: def     !< index of defect in cg%q(:)

      real, parameter     :: L4_scaling = 1./12. ! with L4_strength = 1. this gives an L4 approximation for finite differences approach
      integer, parameter  :: L2w = 2             ! #layers of boundary cells for L2 operator

      real                :: c21, c41, c42 !, c20, c40
      real                :: L0, Lx1, Lx2, Ly1, Ly2, Lz1, Lz2, Lx, Ly, Lz

      logical, save       :: firstcall = .true.
      integer             :: i, j, k
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg

      if (dom%eff_dim<ndims) call die("[multigrid_Laplace4:residual4] Only 3D is implemented")

      if (firstcall) then
         if (master) call warn("[multigrid_Laplace4:residual4] residual order 4 is experimental.")
         firstcall = .false.
      endif

      call cg_llst%arr3d_boundaries(soln, bnd_type = BND_NEGREF)

      c21 = 1.
      c42 = - L4_scaling * L4_strength
      c41 = c21 + 4. * L4_scaling * L4_strength
      !c20 = -2.
      !c40 = c20 - 6. * L4_strength

      cgl => cg_llst%first
      do while (associated(cgl))
         cg => cgl%cg

         Lx1 = c41 * cg%idx2
         Ly1 = c41 * cg%idy2
         Lz1 = c41 * cg%idz2
         Lx2 = c42 * cg%idx2
         Ly2 = c42 * cg%idy2
         Lz2 = c42 * cg%idz2
         ! L0  = c40 * (cg%idx2 + cg%idy2 + cg%idz2 )
         L0 = -2. * (Lx1 + Lx2 + Ly1 + Ly2 + Lz1 + Lz2)

         !> \deprecated BEWARE: cylindrical factors go here
         p3 => cg%q(def)%span(cg%ijkse)
         p3 = cg%q(src)%span(cg%ijkse) - &
              cg%q(soln)%span(cg%ijkse-2*idm2(xdim,:,:)) * Lx2 - cg%q(soln)%span(cg%ijkse+2*idm2(xdim,:,:)) * Lx2 - &
              cg%q(soln)%span(cg%ijkse-  idm2(xdim,:,:)) * Lx1 - cg%q(soln)%span(cg%ijkse+  idm2(xdim,:,:)) * Lx1 - &
              cg%q(soln)%span(cg%ijkse-2*idm2(ydim,:,:)) * Ly2 - cg%q(soln)%span(cg%ijkse+2*idm2(ydim,:,:)) * Ly2 - &
              cg%q(soln)%span(cg%ijkse-  idm2(ydim,:,:)) * Ly1 - cg%q(soln)%span(cg%ijkse+  idm2(ydim,:,:)) * Ly1 - &
              cg%q(soln)%span(cg%ijkse-2*idm2(zdim,:,:)) * Lz2 - cg%q(soln)%span(cg%ijkse+2*idm2(zdim,:,:)) * Lz2 - &
              cg%q(soln)%span(cg%ijkse-  idm2(zdim,:,:)) * Lz1 - cg%q(soln)%span(cg%ijkse+  idm2(zdim,:,:)) * Lz1 - cg%q(soln)%span(cg%ijkse)   * L0

         ! WARNING: not optimized
         if (grav_bnd == bnd_givenval) then ! probably also in some other cases
            ! Use L2 Laplacian in two layers of cells next to the boundary because L4 seems to be incompatible with present image mass construction
            Lx = c21 * cg%idx2
            Ly = c21 * cg%idy2
            Lz = c21 * cg%idz2
            L0 = -2. * (Lx + Ly + Lz)

            do k = cg%ks, cg%ke
               do j = cg%js, cg%je
                  do i = cg%is, cg%ie
                     if ( i<cg%is+L2w .or. i>cg%ie-L2w .or. j<cg%js+L2w .or. j>cg%je-L2w .or. k<cg%ks+L2w .or. k>cg%ke-L2w) then
                        cg%q(def)%arr        (i,   j,   k)   = cg%q(src)%arr (i,   j,   k)         - &
                             ( cg%q(soln)%arr(i-1, j,   k)   + cg%q(soln)%arr(i+1, j,   k))   * Lx - &
                             ( cg%q(soln)%arr(i,   j-1, k)   + cg%q(soln)%arr(i,   j+1, k))   * Ly - &
                             ( cg%q(soln)%arr(i,   j,   k-1) + cg%q(soln)%arr(i,   j,   k+1)) * Lz - &
                             & cg%q(soln)%arr(i,   j,   k)                                    * L0
                     endif
                  enddo
               enddo
            enddo
         endif
         cgl => cgl%nxt
      enddo

   end subroutine residual4

!>
!! \brief Relaxation.
!!
!! \warning For optimal convergence this Laplace operator requires specific relaxation scheme.
!!
!! \details Call relaxation for the residual2 operator while specific relaxation scheme is not implemented and expect severily limited convergence.
!! It is not planned to be implemented anytime soon because 4th order Mehrstellen operator is much better.
!<

   subroutine approximate_solution_rbgs4(curl, src, soln, nsmoo)

      use cg_level_connected, only: cg_level_connected_T
      use multigrid_Laplace2, only: approximate_solution_rbgs2

      implicit none

      type(cg_level_connected_T), pointer, intent(in) :: curl  !< pointer to a level for which we approximate the solution
      integer(kind=4),                     intent(in) :: src   !< index of source in cg%q(:)
      integer(kind=4),                     intent(in) :: soln  !< index of solution in cg%q(:)
      integer,                             intent(in) :: nsmoo !< number of smoothing repetitions

      call approximate_solution_rbgs2(curl, src, soln, nsmoo)

   end subroutine approximate_solution_rbgs4

end module multigrid_Laplace4
