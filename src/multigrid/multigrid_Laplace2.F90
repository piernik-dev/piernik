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
   public :: residual2, approximate_solution_rbgs2, vT_A_v_2

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
      use constants,          only: xdim, ydim, zdim, ndims, GEO_XYZ, GEO_RPZ, BND_NEGREF, LO
      use dataio_pub,         only: die
      use domain,             only: dom
      use global,             only: dirty_debug
      use grid_cont,          only: grid_container
      use multigridvars,      only: multidim_code_3D, overrelax, set_relax_boundaries

      implicit none

      type(cg_level_connected_T), pointer, intent(in) :: curl  !< pointer to a level for which we approximate the solution
      integer(kind=4),                     intent(in) :: src   !< index of source in cg%q(:)
      integer(kind=4),                     intent(in) :: soln  !< index of solution in cg%q(:)
      integer(kind=4),                     intent(in) :: nsmoo !< number of smoothing repetitions

      integer, parameter :: RED_BLACK = 2 !< the checkerboard requires two sweeps

      integer :: n, i, j, k, i1, j1, k1, id, jd, kd
      integer(kind=4) :: b
      integer(kind=8) :: ijko
      real, dimension(:), allocatable :: crx, crx1, cry, crz, cr
      real :: cr0
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
      integer :: is, ie, js, je, ks, ke
      logical :: need_all_bnd_upd

      ! call curl%arr3d_boundaries(src) required when we want to eliminate some communication of soln at a cost of expanding relaxated area into guardcells

      allocate(crx(0), crx1(0), cry(0), crz(0), cr(0)) ! suppress compiler warnings
      cr0 = 1. - overrelax

      if (dom%nb > 1) call curl%internal_boundaries_3d(src)
      do n = 1, RED_BLACK*nsmoo
         need_all_bnd_upd = (mod(n-1, int(dom%nb)) == 0)
         if (need_all_bnd_upd) call curl%arr3d_boundaries(soln, bnd_type = BND_NEGREF)
         b = int(dom%nb - 1 - mod(n-1, int(dom%nb)), kind=4)

         if (dirty_debug) then
            write(dirty_label, '(a,i5)')"relax2 soln- smoo=", n
            call curl%check_dirty(soln, dirty_label)
         endif
         cgl => curl%first
         do while (associated(cgl))
            cg => cgl%cg
            if (dom%geometry_type == GEO_RPZ) then
               deallocate(crx, crx1, cry, crz, cr)
               is = lbound(cg%x, dim=1)
               ie = ubound(cg%x, dim=1)
               allocate(crx(is:ie), crx1(is:ie), cry(is:ie), crz(is:ie), cr(is:ie))
               cr  = overrelax / 2.
               crx = cg%dvol**2 * cg%idx2 * cg%x**2
               cry = cg%dvol**2 * cg%idy2
               crz = cg%dvol**2 * cg%idz2 * cg%x**2
               cr  = cr / (crx + cry + crz)
               crx = cr * crx
               cry = cr * cry
               crz = cr * crz
               cr  = cr * cg%dvol**2 * cg%x**2

               crx1 = 2. * cg%x(is:ie) * cg%idx
               where (crx1 /= 0.) crx1 = 1./crx1
            endif
            call set_relax_boundaries(cg, soln, is, ie, js, je, ks, ke, b, .not. need_all_bnd_upd)

            ! OPT: The "forall" construct would give more clear code that is 4 times slower than implementation with explicit loops to describe a 3-D checkerboard.
            ! OPT: The 8-colored implementation with loops replaced by array operations is ~20% slower

            ! This routine is really sensitive to tiny details such as which cells we do first (red or black).
            ! Before you optimize anything, make sure it does not change the results
            ijko = 0
            if (any(cg%lhn(:, LO) < 0)) ijko = -ndims*RED_BLACK*minval(cg%lhn(:, LO))
            if (dom%eff_dim==ndims .and. .not. multidim_code_3D) then
               do k = ks, ke
                  do j = js, je
                     i1 = is + int(mod(ijko+n+is+j+k, int(RED_BLACK, kind=8)), kind=4)
                     select case (dom%geometry_type)
                        case (GEO_RPZ)
                           cg%q(soln)%arr                      (i1  :ie  :2, j,   k) = &
                                cr0           *  cg%q(soln)%arr(i1  :ie  :2, j,   k) + &
                                crx( i1:ie:2) * (cg%q(soln)%arr(i1-1:ie-1:2, j,   k  ) + cg%q(soln)%arr(i1+1:ie+1:2, j,   k))   + &
                                cry( i1:ie:2) * (cg%q(soln)%arr(i1  :ie  :2, j-1, k  ) + cg%q(soln)%arr(i1  :ie  :2, j+1, k))   + &
                                crz( i1:ie:2) * (cg%q(soln)%arr(i1  :ie  :2, j,   k-1) + cg%q(soln)%arr(i1  :ie  :2, j,   k+1)) - &
                                cr(  i1:ie:2) *  cg%q(src )%arr(i1  :ie  :2, j,   k  )  + &
                                crx(i1:ie:2) * crx1(i1:ie:2) * &
                                &               (cg%q(soln)%arr(i1+1:ie+1:2, j,   k  ) - cg%q(soln)%arr(i1-1:ie-1:2, j,   k))
                        case (GEO_XYZ)
                           cg%q(soln)%arr                 (i1  :ie  :2, j,   k) = &
                                cr0      *  cg%q(soln)%arr(i1  :ie  :2, j,   k) + &
                                cg%mg%rx * (cg%q(soln)%arr(i1-1:ie-1:2, j,   k)   + cg%q(soln)%arr(i1+1:ie+1:2, j,   k))   + &
                                cg%mg%ry * (cg%q(soln)%arr(i1  :ie  :2, j-1, k)   + cg%q(soln)%arr(i1  :ie  :2, j+1, k))   + &
                                cg%mg%rz * (cg%q(soln)%arr(i1  :ie  :2, j,   k-1) + cg%q(soln)%arr(i1  :ie  :2, j,   k+1)) - &
                                cg%mg%r  *  cg%q(src )%arr(i1  :ie  :2, j,   k)
                        case default
                           call die("[multigrid_Laplace2:approximate_solution_rbgs] Unsupported geometry (3D).")
                     end select
                  enddo
               enddo
            else
               ! In 3D this variant significantly increases instruction count and also some data read
               i1 = is; id = 1 ! mv to multigridvars, init_multigrid
               j1 = js; jd = 1
               k1 = ks; kd = 1
               if (dom%has_dir(xdim)) then
                  id = RED_BLACK
               else if (dom%has_dir(ydim)) then
                  jd = RED_BLACK
               else if (dom%has_dir(zdim)) then
                  kd = RED_BLACK
               endif

               if (kd == RED_BLACK) k1 = ks + int(mod(ijko+n+ks, int(RED_BLACK, kind=8)), kind=4)
               select case (dom%geometry_type)
                  case (GEO_XYZ)
                     do k = k1, ke, kd
                        if (jd == RED_BLACK) j1 = js + int(mod(ijko+n+js+k, int(RED_BLACK, kind=8)), kind=4)
                        do j = j1, je, jd
                           if (id == RED_BLACK) i1 = is + int(mod(ijko+n+is+j+k, int(RED_BLACK, kind=8)), kind=4)
                           cg%q(soln)%arr                       (i1  :ie  :id, j,   k)   = &
                                &                 cg%q(soln)%arr(i1  :ie  :id, j,   k)   * cr0 - &
                                &                 cg%q(src)%arr (i1  :ie  :id, j,   k)   * cg%mg%r
                           if (dom%has_dir(xdim)) cg%q(soln)%arr(i1  :ie  :id, j,   k)   = cg%q(soln)%arr(i1  :ie  :id, j,   k)    + &
                                &                (cg%q(soln)%arr(i1-1:ie-1:id, j,   k)   + cg%q(soln)%arr(i1+1:ie+1:id, j,   k))   * cg%mg%rx
                           if (dom%has_dir(ydim)) cg%q(soln)%arr(i1  :ie  :id, j,   k)   = cg%q(soln)%arr(i1  :ie  :id, j,   k)    + &
                                &                (cg%q(soln)%arr(i1  :ie  :id, j-1, k)   + cg%q(soln)%arr(i1  :ie  :id, j+1, k))   * cg%mg%ry
                           if (dom%has_dir(zdim)) cg%q(soln)%arr(i1  :ie  :id, j,   k)   = cg%q(soln)%arr(i1  :ie  :id, j,   k)    + &
                                &                (cg%q(soln)%arr(i1  :ie  :id, j,   k-1) + cg%q(soln)%arr(i1  :ie  :id, j,   k+1)) * cg%mg%rz
                        enddo
                     enddo
                  case (GEO_RPZ)
                     do k = k1, ke, kd
                        if (jd == RED_BLACK) j1 = js + int(mod(ijko+n+js+k, int(RED_BLACK, kind=8)), kind=4)
                        do j = j1, je, jd
                           if (id == RED_BLACK) i1 = is + int(mod(ijko+n+is+j+k, int(RED_BLACK, kind=8)), kind=4)
                           do i = i1, ie, id
                              cg%q(soln)%arr                       (i,   j,   k)   = &
                                   &                 cg%q(soln)%arr(i,   j,   k)   * cr0 - &
                                   &                 cg%q(src )%arr(i,   j,   k)   * cr(i)
                              if (dom%has_dir(xdim)) cg%q(soln)%arr(i,   j,   k)   = cg%q(soln)%arr(i,   j,   k)    + &
                                   &                (cg%q(soln)%arr(i-1, j,   k)   + cg%q(soln)%arr(i+1, j,   k))   * crx(i) + &
                                   &                (cg%q(soln)%arr(i+1, j,   k)   - cg%q(soln)%arr(i-1, j,   k))   * crx(i) * crx1(i)
                              if (dom%has_dir(ydim)) cg%q(soln)%arr(i,   j,   k)   = cg%q(soln)%arr(i,   j,   k)    + &
                                   &                (cg%q(soln)%arr(i,   j-1, k)   + cg%q(soln)%arr(i,   j+1, k))   * cry(i)
                              if (dom%has_dir(zdim)) cg%q(soln)%arr(i,   j,   k)   = cg%q(soln)%arr(i,   j,   k)    + &
                                   &                (cg%q(soln)%arr(i,   j,   k-1) + cg%q(soln)%arr(i,   j,   k+1)) * crz(i)
                           enddo
                        enddo
                     enddo
                  case default
                     call die("[multigrid_Laplace2:approximate_solution_rbgs] Unsupported geometry (1D or 2D).")
               end select
            endif

            cgl => cgl%nxt
         enddo

         if (dirty_debug) then
            write(dirty_label, '(a,i5)')"relax2 soln+ smoo=", n
            call curl%check_dirty(soln, dirty_label)
         endif

      enddo

      deallocate(crx, crx1, cry, crz, cr)

   end subroutine approximate_solution_rbgs2

!>
!! \brief Compute var*Laplacian(var) for CG algorithm
!!
!! \details This implementation uses 2nd order approximation of the Laplace operator
!!
!! This is the first implementation that also serves for 4th order laplacians as well.
!<

   real function vT_A_v_2(var)

      use cg_leaves,  only: leaves
      use cg_list,    only: cg_list_element
      use constants,  only: GEO_XYZ, BND_NEGREF, ndims, pSUM
      use dataio_pub, only: die
      use domain,     only: dom
      use grid_cont,  only: grid_container
      use mpisetup,   only: piernik_MPI_Allreduce

      implicit none

      integer(kind=4), intent(in) :: var !< variable, for which the operation is done

      type(cg_list_element), pointer  :: cgl
      type(grid_container),  pointer  :: cg
      real                            :: L0, Lx, Ly, Lz

      if (dom%geometry_type /= GEO_XYZ) call die("[multigrid_Laplace2:vT_A_v_2] Implemented only for cartesian coords as yet")
      if (dom%eff_dim /= ndims) call die("[multigrid_Laplace2:vT_A_v_2] Implemented only for 3D")

      call leaves%leaf_arr3d_boundaries(var, bnd_type = BND_NEGREF)

      vT_A_v_2 = 0.
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         Lx = cg%idx2
         Ly = cg%idy2
         Lz = cg%idz2
         L0 = -2. * (Lx + Ly + Lz)

         vT_A_v_2 = vT_A_v_2 + sum( &
               & cg%q(var)%arr(cg%is  :cg%ie,   cg%js  :cg%je,   cg%ks  :cg%ke) * ( &
               ( cg%q(var)%arr(cg%is-1:cg%ie-1, cg%js  :cg%je,   cg%ks  :cg%ke)         + &
               & cg%q(var)%arr(cg%is+1:cg%ie+1, cg%js  :cg%je,   cg%ks  :cg%ke))   * Lx + &
               ( cg%q(var)%arr(cg%is  :cg%ie,   cg%js-1:cg%je-1, cg%ks  :cg%ke)         + &
               & cg%q(var)%arr(cg%is  :cg%ie,   cg%js+1:cg%je+1, cg%ks  :cg%ke))   * Ly + &
               ( cg%q(var)%arr(cg%is  :cg%ie,   cg%js  :cg%je,   cg%ks-1:cg%ke-1)       + &
               & cg%q(var)%arr(cg%is  :cg%ie,   cg%js  :cg%je,   cg%ks+1:cg%ke+1)) * Lz + &
               & cg%q(var)%arr(cg%is  :cg%ie,   cg%js  :cg%je,   cg%ks  :cg%ke)    * L0), mask=cg%leafmap)
         cgl => cgl%nxt
      enddo

      call piernik_MPI_Allreduce(vT_A_v_2, pSUM)

   end function vT_A_v_2

end module multigrid_Laplace2
