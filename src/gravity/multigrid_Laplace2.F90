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

!> \brief Implementation of 2nd order, simplest Laplace operator

module multigrid_Laplace2
! pulled by MULTIGRID && GRAV

   implicit none

   private
   public :: residual2, approximate_solution_rbgs2

contains

!>
!! \brief 2nd order Laplacian
!!
!! \details
!! 1D: L2 = 1/dx**2 [ 1 -2 1 ]
!!
!!                  |        |           |    1   |
!! 2D: L2 = 1/dx**2 | 1 -2 1 | + 1/dy**2 |   -2   |
!!                  |        |           |    1   |
!!
!! 3D: similar
!<

   subroutine residual2(cg_llst, src, soln, def)

      use cg_list,       only: cg_list_element
      use cg_leaves,     only: cg_leaves_T
      use constants,     only: xdim, ydim, zdim, ndims, GEO_XYZ, GEO_RPZ, half, BND_NEGREF
      use dataio_pub,    only: die
      use domain,        only: dom
      use grid_cont,     only: grid_container
      use multigridvars, only: multidim_code_3D

      implicit none

      class(cg_leaves_T), intent(in) :: cg_llst !< pointer to a level for which we approximate the solution
      integer(kind=4),    intent(in) :: src     !< index of source in cg%q(:)
      integer(kind=4),    intent(in) :: soln    !< index of solution in cg%q(:)
      integer(kind=4),    intent(in) :: def     !< index of defect in cg%q(:)

      integer                         :: i, j, k
      real                            :: L0, Lx, Ly, Lz
      real, dimension(:), allocatable :: Lx1, Ly_a, L0_a
      type(cg_list_element), pointer  :: cgl
      type(grid_container),  pointer  :: cg

      call cg_llst%leaf_arr3d_boundaries(soln, bnd_type = BND_NEGREF)
      ! corners are required for non-cartesian decompositions because current implementation of arr3d_boundaries may use overlapping buffers at triple points

      ! Possible optimization candidate: reduce cache misses (secondary importance, cache-aware implementation required)
      ! Explicit loop over k gives here better performance than array operation due to less cache misses (at least on 32^3 and 64^3 arrays)
      cgl => cg_llst%first
      do while (associated(cgl))
         cg => cgl%cg

         ! Coefficients for a simplest 3-point Laplacian operator: [ 1, -2, 1 ]
         ! for 2D and 1D setups appropriate elements of [ Lx, Ly, Lz ] should be == 0.
         Lx = cg%idx2
         Ly = cg%idy2
         Lz = cg%idz2
         L0 = -2. * (Lx + Ly + Lz)

         select case (dom%geometry_type)
            case (GEO_XYZ)
               if (dom%eff_dim == ndims .and. .not. multidim_code_3D) then
                  do k = cg%ks, cg%ke
                     cg%q(def)%arr        (cg%is  :cg%ie,   cg%js  :cg%je,   k)         = &
                          & cg%q(src)%arr (cg%is  :cg%ie,   cg%js  :cg%je,   k)         - &
                          ( cg%q(soln)%arr(cg%is-1:cg%ie-1, cg%js  :cg%je,   k)         + &
                          & cg%q(soln)%arr(cg%is+1:cg%ie+1, cg%js  :cg%je,   k))   * Lx - &
                          ( cg%q(soln)%arr(cg%is  :cg%ie,   cg%js-1:cg%je-1, k)         + &
                          & cg%q(soln)%arr(cg%is  :cg%ie,   cg%js+1:cg%je+1, k))   * Ly - &
                          ( cg%q(soln)%arr(cg%is  :cg%ie,   cg%js  :cg%je,   k-1)       + &
                          & cg%q(soln)%arr(cg%is  :cg%ie,   cg%js  :cg%je,   k+1)) * Lz - &
                          & cg%q(soln)%arr(cg%is  :cg%ie,   cg%js  :cg%je,   k)    * L0
                  enddo
               else
                  ! In 3D this implementation can give a bit more cache misses, few times more writes and significantly more instructions executed than monolithic 3D above
                  do k = cg%ks, cg%ke
                     cg%q(def)%arr        (cg%is  :cg%ie,   cg%js  :cg%je,   k)    = &
                          & cg%q(src)%arr (cg%is  :cg%ie,   cg%js  :cg%je,   k)    - &
                          & cg%q(soln)%arr(cg%is  :cg%ie,   cg%js  :cg%je,   k)    * L0
                     if (dom%has_dir(xdim)) &
                          & cg%q(def)%arr (cg%is  :cg%ie,   cg%js  :cg%je,   k)    = &
                          & cg%q(def)%arr (cg%is  :cg%ie,   cg%js  :cg%je,   k)    - &
                          ( cg%q(soln)%arr(cg%is-1:cg%ie-1, cg%js  :cg%je,   k)    + &
                          & cg%q(soln)%arr(cg%is+1:cg%ie+1, cg%js  :cg%je,   k))   * Lx
                     if (dom%has_dir(ydim)) &
                          & cg%q(def)%arr (cg%is  :cg%ie,   cg%js  :cg%je,   k)    = &
                          & cg%q(def)%arr (cg%is  :cg%ie,   cg%js  :cg%je,   k)    - &
                          ( cg%q(soln)%arr(cg%is  :cg%ie,   cg%js-1:cg%je-1, k)    + &
                          & cg%q(soln)%arr(cg%is  :cg%ie,   cg%js+1:cg%je+1, k))   * Ly
                     if (dom%has_dir(zdim)) &
                          & cg%q(def)%arr (cg%is  :cg%ie,   cg%js  :cg%je,   k)    = &
                          & cg%q(def)%arr (cg%is  :cg%ie,   cg%js  :cg%je,   k)    - &
                          ( cg%q(soln)%arr(cg%is  :cg%ie,   cg%js  :cg%je,   k-1)  + &
                          & cg%q(soln)%arr(cg%is  :cg%ie,   cg%js  :cg%je,   k+1)) * Lz
                  enddo
               endif
            case (GEO_RPZ)
!               Lx = cg%idx2 ! already set
!               Lz = cg%idz2 ! already set
               !> \todo convert Lx1, Ly_a and L0_a into precomputed arrays
               allocate(Lx1(cg%is:cg%ie), Ly_a(cg%is:cg%ie), L0_a(cg%is:cg%ie))
               Ly_a(cg%is:cg%ie) = cg%idy2 * cg%inv_x(cg%is:cg%ie)**2 ! cylindrical factor
               Lx1 (cg%is:cg%ie) = half * (cg%idx * cg%inv_x(cg%is:cg%ie))
               L0_a(cg%is:cg%ie) = -2. * (Lx + Ly_a(cg%is:cg%ie) + Lz)
               do k = cg%ks, cg%ke
                  do j = cg%js, cg%je
                     do i = cg%is, cg%ie
                        cg%q(def)%arr(i, j, k) = cg%q(src)%arr (i,   j,   k)   - cg%q(soln)%arr(i,   j,   k)    * L0_a(i)
                        if (dom%has_dir(xdim))   cg%q(def)%arr (i,   j,   k)   = cg%q(def) %arr(i,   j,   k)         - &
                             &                  (cg%q(soln)%arr(i+1, j,   k)   + cg%q(soln)%arr(i-1, j,   k))   * Lx - &
                             &                  (cg%q(soln)%arr(i+1, j,   k)   - cg%q(soln)%arr(i-1, j,   k))   * Lx1(i)    ! cylindrical term
                        if (dom%has_dir(ydim))   cg%q(def)%arr (i,   j,   k)   = cg%q(def) %arr(i,   j,   k)         - &
                             &                  (cg%q(soln)%arr(i,   j+1, k)   + cg%q(soln)%arr(i,   j-1, k))   * Ly_a(i)
                        if (dom%has_dir(zdim))   cg%q(def)%arr (i,   j,   k)   = cg%q(def) %arr(i,   j,   k)         - &
                             &                  (cg%q(soln)%arr(i,   j,   k+1) + cg%q(soln)%arr(i,   j,   k-1)) * Lz
                     enddo
                  enddo
               enddo
               deallocate(Lx1, Ly_a, L0_a)
            case default
               call die("[multigrid_Laplace2:residual2] Unsupported geometry.")
         end select
         cgl => cgl%nxt
      enddo

    end subroutine residual2

!>
!! \brief Red-Black Gauss-Seidel relaxation for Laplace operator implemented in residual2
!!
!! \details  This relaxation can also be used for some other implementations of the Laplace operators during development stage. In such cases you may expect poor convergence.
!! 4th order operator residual4 uses this relaxation because it is not planned to implement specialized operator anytime soon.
!!
!! \todo Implement efficient use of guardcells in order to save some communication like in [7846]
!<

   subroutine approximate_solution_rbgs2(curl, src, soln, nsmoo)

      use cg_level_connected, only: cg_level_connected_T
      use cg_list,            only: cg_list_element
      use cg_list_dataop,     only: dirty_label
      use constants,          only: xdim, ydim, zdim, ndims, GEO_XYZ, GEO_RPZ, I_ONE, BND_NEGREF, LO
      use dataio_pub,         only: die
      use domain,             only: dom
      use global,             only: dirty_debug
      use grid_cont,          only: grid_container
      use multigridvars,      only: multidim_code_3D, Jacobi_damp, overrelax, overrelax_xyz

      implicit none

      type(cg_level_connected_T), pointer, intent(in) :: curl  !< pointer to a level for which we approximate the solution
      integer(kind=4),                     intent(in) :: src   !< index of source in cg%q(:)
      integer(kind=4),                     intent(in) :: soln  !< index of solution in cg%q(:)
      integer,                             intent(in) :: nsmoo !< number of smoothing repetitions

      integer, parameter :: RED_BLACK = 2 !< the checkerboard requires two sweeps

      integer :: n, i, j, k, i1, j1, k1, id, jd, kd
      integer(kind=8) :: ijko
      real    :: crx, crx1, cry, crz, cr
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg

      ! call curl%arr3d_boundaries(src) required when we want to eliminate some communication of soln at a cost of expanding relaxated area into guardcells

      do n = 1, RED_BLACK*nsmoo
         call curl%arr3d_boundaries(soln, bnd_type = BND_NEGREF)

         if (dirty_debug) then
            write(dirty_label, '(a,i5)')"relax2 soln- smoo=", n
            call curl%check_dirty(soln, dirty_label, expand=I_ONE)
         endif
         cgl => curl%first
         do while (associated(cgl))
            cg => cgl%cg

            ! Possible optimization: this is the most costly part of the RBGS relaxation (instruction count, read and write data, L1 and L2 read cache miss)
            ! do n = 1, nsmoo
            !    call curl%arr3d_boundaries(soln, bnd_type = BND_NEGREF)
            !    relax single layer of red cells at all faces
            !    call curl%arr3d_boundaries(soln, bnd_type = BND_NEGREF)
            !    relax interior cells (except for single layer of cells at all faces), first red, then 1-cell behind black one.
            !    relax single layer of black cells at all faces
            ! enddo
            ! OPT: try to relax without Red-Black and use cg%wa for temporary storage

            ! with explicit outer loops it is easier to describe a 3-D checkerboard :-)

            ijko = 0
            if (any(cg%my_se(:, LO) < 0)) ijko = -ndims*RED_BLACK*minval(cg%my_se(:, LO))
            if (dom%eff_dim==ndims .and. .not. multidim_code_3D) then
               do k = cg%ks, cg%ke
                  do j = cg%js, cg%je
                     i1 = cg%is + int(mod(ijko+n+cg%is+j+k, int(RED_BLACK, kind=8)), kind=4)
                     if (dom%geometry_type == GEO_RPZ) then
!!$                  cg%q(soln)%arr(i1  :cg%ie  :2, j,   k) = &
!!$                       cg%mg%rx * (cg%q(soln)%arr(i1-1:cg%ie-1:2, j,   k  ) + cg%q(soln)%arr(i1+1:cg%ie+1:2, j,   k))   + &
!!$                       cg%mg%ry * (cg%q(soln)%arr(i1  :cg%ie  :2, j-1, k  ) + cg%q(soln)%arr(i1  :cg%ie  :2, j+1, k))   + &
!!$                       cg%mg%rz * (cg%q(soln)%arr(i1  :cg%ie  :2, j,   k-1) + cg%q(soln)%arr(i1  :cg%ie  :2, j,   k+1)) - &
!!$                       cg%mg%r  *  cg%q(src)%arr( i1  :cg%ie  :2, j,   k  )  + &
!!$                       cg%mg%rx * (cg%q(soln)%arr(i1+1:cg%ie+1:2, j,   k  ) - cg%q(soln)%arr(i1-1:cg%ie-1:2, j,   k)) * fac(i1:cg%ie:2)
                        call die("[multigrid_Laplace2:approximate_solution_rbgs2] This variant of relaxation loop is not implemented for cylindrical coordinates.")
                     else
                        cg%q(soln)%arr(i1  :cg%ie  :2, j,   k) = &
                             cg%mg%rx * (cg%q(soln)%arr(i1-1:cg%ie-1:2, j,   k)   + cg%q(soln)%arr(i1+1:cg%ie+1:2, j,   k))   + &
                             cg%mg%ry * (cg%q(soln)%arr(i1  :cg%ie  :2, j-1, k)   + cg%q(soln)%arr(i1  :cg%ie  :2, j+1, k))   + &
                             cg%mg%rz * (cg%q(soln)%arr(i1  :cg%ie  :2, j,   k-1) + cg%q(soln)%arr(i1  :cg%ie  :2, j,   k+1)) - &
                             cg%mg%r  *  cg%q(src)%arr (i1  :cg%ie  :2, j,   k)
                     endif
                  enddo
               enddo
            else
               ! In 3D this variant significantly increases instruction count and also some data read
               i1 = cg%is; id = 1 ! mv to multigridvars, init_multigrid
               j1 = cg%js; jd = 1
               k1 = cg%ks; kd = 1
               if (dom%has_dir(xdim)) then
                  id = RED_BLACK
               else if (dom%has_dir(ydim)) then
                  jd = RED_BLACK
               else if (dom%has_dir(zdim)) then
                  kd = RED_BLACK
               endif

               if (kd == RED_BLACK) k1 = cg%ks + int(mod(ijko+n+cg%ks, int(RED_BLACK, kind=8)), kind=4)
               select case (dom%geometry_type)
                  case (GEO_XYZ)
                     do k = k1, cg%ke, kd
                        if (jd == RED_BLACK) j1 = cg%js + int(mod(ijko+n+cg%js+k, int(RED_BLACK, kind=8)), kind=4)
                        do j = j1, cg%je, jd
                           if (id == RED_BLACK) i1 = cg%is + int(mod(ijko+n+cg%is+j+k, int(RED_BLACK, kind=8)), kind=4)
                           cg%q(soln)%arr                           (i1  :cg%ie  :id, j,   k)   = &
                                & (1. - Jacobi_damp)* cg%q(soln)%arr(i1  :cg%ie  :id, j,   k)   - &
                                &       Jacobi_damp * cg%q(src)%arr (i1  :cg%ie  :id, j,   k)   * cg%mg%r
                           if (dom%has_dir(xdim))     cg%q(soln)%arr(i1  :cg%ie  :id, j,   k)   = cg%q(soln)%arr(i1  :cg%ie  :id, j,   k)    + &
                                &       Jacobi_damp *(cg%q(soln)%arr(i1-1:cg%ie-1:id, j,   k)   + cg%q(soln)%arr(i1+1:cg%ie+1:id, j,   k))   * cg%mg%rx
                           if (dom%has_dir(ydim))     cg%q(soln)%arr(i1  :cg%ie  :id, j,   k)   = cg%q(soln)%arr(i1  :cg%ie  :id, j,   k)    + &
                                &       Jacobi_damp *(cg%q(soln)%arr(i1  :cg%ie  :id, j-1, k)   + cg%q(soln)%arr(i1  :cg%ie  :id, j+1, k))   * cg%mg%ry
                           if (dom%has_dir(zdim))     cg%q(soln)%arr(i1  :cg%ie  :id, j,   k)   = cg%q(soln)%arr(i1  :cg%ie  :id, j,   k)    + &
                                &       Jacobi_damp *(cg%q(soln)%arr(i1  :cg%ie  :id, j,   k-1) + cg%q(soln)%arr(i1  :cg%ie  :id, j,   k+1)) * cg%mg%rz
                        enddo
                     enddo
                  case (GEO_RPZ)
                     do k = k1, cg%ke, kd
                        if (jd == RED_BLACK) j1 = cg%js + int(mod(ijko+n+cg%js+k, int(RED_BLACK, kind=8)), kind=4)
                        do j = j1, cg%je, jd
                           if (id == RED_BLACK) i1 = cg%is + int(mod(ijko+n+cg%is+j+k, int(RED_BLACK, kind=8)), kind=4)
                           do i = i1, cg%ie, id
                              cr  = overrelax / 2.
                              crx = cg%dvol**2 * cg%idx2 * cg%x(i)**2
                              cry = cg%dvol**2 * cg%idy2
                              crz = cg%dvol**2 * cg%idz2 * cg%x(i)**2
                              cr  = cr / (crx + cry + crz)
                              crx = overrelax_xyz(xdim)* crx * cr
                              cry = overrelax_xyz(ydim)* cry * cr
                              crz = overrelax_xyz(zdim)* crz * cr
                              cr  = cr * cg%dvol**2 * cg%x(i)**2

                              crx1 = 2. * cg%x(i) * cg%idx
                              if (crx1 /= 0.) crx1 = 1./crx1
                              cg%q(soln)%arr                           (i,   j,   k)   = &
                                   & (1. - Jacobi_damp)* cg%q(soln)%arr(i,   j,   k)   - &
                                   &       Jacobi_damp * cg%q(src)%arr (i,   j,   k)   * cr
                              if (dom%has_dir(xdim))     cg%q(soln)%arr(i,   j,   k)   = cg%q(soln)%arr(i,   j,   k)    + &
                                   &       Jacobi_damp *(cg%q(soln)%arr(i-1, j,   k)   + cg%q(soln)%arr(i+1, j,   k))   * crx + &
                                   &       Jacobi_damp *(cg%q(soln)%arr(i+1, j,   k)   - cg%q(soln)%arr(i-1, j,   k))   * crx * crx1
                              if (dom%has_dir(ydim))     cg%q(soln)%arr(i,   j,   k)   = cg%q(soln)%arr(i,   j,   k)    + &
                                   &       Jacobi_damp *(cg%q(soln)%arr(i,   j-1, k)   + cg%q(soln)%arr(i,   j+1, k))   * cry
                              if (dom%has_dir(zdim))     cg%q(soln)%arr(i,   j,   k)   = cg%q(soln)%arr(i,   j,   k)    + &
                                   &       Jacobi_damp *(cg%q(soln)%arr(i,   j,   k-1) + cg%q(soln)%arr(i,   j,   k+1)) * crz
                           enddo
                        enddo
                     enddo
                  case default
                     call die("[multigrid_Laplace2:approximate_solution_rbgs] Unsupported geometry.")
               end select
            endif
            cgl => cgl%nxt
         enddo

         if (dirty_debug) then
            write(dirty_label, '(a,i5)')"relax2 soln+ smoo=", n
            call curl%check_dirty(soln, dirty_label)
         endif

      enddo

   end subroutine approximate_solution_rbgs2


end module multigrid_Laplace2
