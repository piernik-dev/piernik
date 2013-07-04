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

!>
!! \brief Implementation of 4th order, compact, Mehrstellen-type Laplace operator
!!
!! \warning It seems that for isolated boundaries some additional modifications need to be done in order to improve accuracy.
!! Current implementation produces significant error in the first layer of cells next to the outer boundary
!<

module multigrid_Laplace4M
! pulled by MULTIGRID && GRAV

   implicit none

   private
   public :: residual_Mehrstellen, approximate_solution_relax4M

contains

!>
!! \brief Compute the residual using compact, 4th order Mehrstellen formula
!!
!! \details 2nd order Laplace operator contains some 4th order terms on the potential size. To improve accuracy we may add a little laplacian of the source on the RHS
!! to match these 4th order terms. Careful calculation leads to such numerical operators:
!!
!!                     |   1   |                    | 1   4 1 |
!! 2D: residuum = 1/12 | 1 8 1 | source - 1/6 h**-2 | 4 -20 4 | solution
!!                     |   1   |                    | 1   4 1 |
!!
!!                     |       | |   1   | |       |                    |   1   | | 1   2 1 | |   1   |
!! 3D: residuum = 1/12 |   1   | | 1 6 1 | |   1   | source - 1/6 h**-2 | 1 2 1 | | 2 -24 2 | | 1 2 1 | solution
!!                     |       | |   1   | |       |                    |   1   | | 1   2 1 | |   1   |
!!
!! where h = cg%dx = cg%dy = cg%dz
!! On an unequally spaced grid (cg%dx /= cg%dy /= cg%dz), the stencil convolved with the solution is a bit more complicated:
!!
!!                 |  1  -2  1 |               |  1  10  1 |
!! 2D: 1/12 dx**-2 | 10 -20 10 | + 1/12 dy**-2 | -2 -20 -2 |
!!                 |  1  -2  1 |               |  1  10  1 |
!!
!!                 |        | | 1  -2 1 | |        |               |   1   | |  1   8  1 | |   1   |               |   1   | |     -2    | |   1   |
!! 3D: 1/12 dx**-2 | 1 -2 1 | | 8 -16 8 | | 1 -2 1 | + 1/12 dy**-2 |  -2   | | -2 -16 -2 | |  -2   | + 1/12 dz**-2 | 1 8 1 | | -2 -16 -2 | | 1 8 1 |
!!                 |        | | 1  -2 1 | |        |               |   1   | |  1   8  1 | |   1   |               |   1   | |     -2    | |   1   |
!!
!! For even higher order of accuracy look for HODIE schemes which are more general then the one described above.
!!
!! \warning This implementation does not support cylindrical coordinates yet
!<

   subroutine residual_Mehrstellen(cg_llst, src, soln, def)

      use cg_list,       only: cg_list_element
      use cg_leaves,     only: cg_leaves_T
      use constants,     only: ndims, xdim, ydim, zdim, BND_NEGREF, LO, HI, GEO_XYZ
      use dataio_pub,    only: die
      use domain,        only: dom
      use grid_cont,     only: grid_container
      use multigridvars, only: grav_bnd, bnd_givenval, multidim_code_3D
      use named_array,   only: p3

      implicit none

      class(cg_leaves_T), intent(in) :: cg_llst !< pointer to a level for which we approximate the solution
      integer(kind=4),    intent(in) :: src     !< index of source in cg%q(:)
      integer(kind=4),    intent(in) :: soln    !< index of solution in cg%q(:)
      integer(kind=4),    intent(in) :: def     !< index of defect in cg%q(:)

      real                :: L0, Lx, Ly, Lz, Lxy, Lxz, Lyz, fac

      integer             :: i
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
      integer(kind=4), dimension(ndims,ndims,LO:HI) :: idm
      real :: src_lapl

      if (dom%geometry_type /= GEO_XYZ) call die("[multigrid_Laplace4M:residual_Mehrstellen] Unsupported geometry")

      src_lapl = 1./12.
      if (grav_bnd == bnd_givenval) src_lapl = 0.
      ! the contribution of outer potential is simulated by a single layer of cells with image of density and we don't want to operate on this structure with the Laplacian.
      ! This image density is supposed to be infinitesimally thin, which we obviously can't reproduce, so we modify the operator instead
      call cg_llst%leaf_arr3d_boundaries(soln, bnd_type=BND_NEGREF)
      if (src_lapl /= 0.) call cg_llst%leaf_arr3d_boundaries(src, bnd_type=BND_NEGREF, nocorners=.true.)

      idm = 0
      do i = xdim, zdim
         if (dom%has_dir(i)) idm(i, i, :) = 1
      enddo

      cgl => cg_llst%first
      do while (associated(cgl))
         cg => cgl%cg

         Lxy = 0. ; if (dom%has_dir(xdim) .and. dom%has_dir(ydim)) Lxy = (cg%idx2 + cg%idy2) / 12.
         Lxz = 0. ; if (dom%has_dir(xdim) .and. dom%has_dir(zdim)) Lxz = (cg%idx2 + cg%idz2) / 12.
         Lyz = 0. ; if (dom%has_dir(ydim) .and. dom%has_dir(zdim)) Lyz = (cg%idy2 + cg%idz2) / 12.
         Lx  = 0. ; if (dom%has_dir(xdim)) Lx = (cg%idx2 - 2. * (Lxy + Lxz))
         Ly  = 0. ; if (dom%has_dir(ydim)) Ly = (cg%idy2 - 2. * (Lxy + Lyz))
         Lz  = 0. ; if (dom%has_dir(zdim)) Lz = (cg%idz2 - 2. * (Lxz + Lyz))
         L0  = -(-2. * (Lx + Ly + Lz) - 4. * (Lxy + Lxz + Lyz))
         Lx = -Lx ; Ly = -Ly ; Lz = -Lz
         Lxy = -Lxy; Lxz = -Lxz; Lyz = -Lyz

         if (dom%eff_dim == ndims .and. .not. multidim_code_3D) then
            ! There's small speed up vs multidim_code_3D
            fac = (1.0 - src_lapl * 2.0 * dom%eff_dim)
            cg%q(def)%arr(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) = &
               fac * cg%q(src)%arr(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) + L0 * cg%q(soln)%arr(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) + &
                    (cg%q(src)%arr(cg%is-1:cg%ie-1, cg%js:cg%je,     cg%ks:cg%ke    ) +  cg%q(src)%arr(cg%is+1:cg%ie+1, cg%js:cg%je,     cg%ks:cg%ke    )) * src_lapl  &
                 + (cg%q(soln)%arr(cg%is-1:cg%ie-1, cg%js:cg%je,     cg%ks:cg%ke    ) + cg%q(soln)%arr(cg%is+1:cg%ie+1, cg%js:cg%je,     cg%ks:cg%ke    )) * Lx        &
                 +  (cg%q(src)%arr(cg%is:cg%ie,     cg%js-1:cg%je-1, cg%ks:cg%ke    ) +  cg%q(src)%arr(cg%is:cg%ie,     cg%js+1:cg%je+1, cg%ks:cg%ke    )) * src_lapl  &
                 + (cg%q(soln)%arr(cg%is:cg%ie,     cg%js-1:cg%je-1, cg%ks:cg%ke    ) + cg%q(soln)%arr(cg%is:cg%ie,     cg%js+1:cg%je+1, cg%ks:cg%ke    )) * Ly        &
                 + (cg%q(soln)%arr(cg%is-1:cg%ie-1, cg%js-1:cg%je-1, cg%ks:cg%ke    ) + cg%q(soln)%arr(cg%is-1:cg%ie-1, cg%js+1:cg%je+1, cg%ks:cg%ke    ) + &
                 &  cg%q(soln)%arr(cg%is+1:cg%ie+1, cg%js-1:cg%je-1, cg%ks:cg%ke    ) + cg%q(soln)%arr(cg%is+1:cg%ie+1, cg%js+1:cg%je+1, cg%ks:cg%ke    )) * Lxy       &
                 +  (cg%q(src)%arr(cg%is:cg%ie,     cg%js:cg%je,     cg%ks-1:cg%ke-1) +  cg%q(src)%arr(cg%is:cg%ie,     cg%js:cg%je,     cg%ks+1:cg%ke+1)) * src_lapl  &
                 + (cg%q(soln)%arr(cg%is:cg%ie,     cg%js:cg%je,     cg%ks-1:cg%ke-1) + cg%q(soln)%arr(cg%is:cg%ie,     cg%js:cg%je,     cg%ks+1:cg%ke+1)) * Lz        &
                 + (cg%q(soln)%arr(cg%is-1:cg%ie-1, cg%js:cg%je,     cg%ks-1:cg%ke-1) + cg%q(soln)%arr(cg%is-1:cg%ie-1, cg%js:cg%je,     cg%ks+1:cg%ke+1) + &
                 &  cg%q(soln)%arr(cg%is+1:cg%ie+1, cg%js:cg%je,     cg%ks-1:cg%ke-1) + cg%q(soln)%arr(cg%is+1:cg%ie+1, cg%js:cg%je,     cg%ks+1:cg%ke+1)) * Lxz       &
                 + (cg%q(soln)%arr(cg%is:cg%ie,     cg%js-1:cg%je-1, cg%ks-1:cg%ke-1) + cg%q(soln)%arr(cg%is:cg%ie,     cg%js+1:cg%je+1, cg%ks-1:cg%ke-1) + &
                 &  cg%q(soln)%arr(cg%is:cg%ie,     cg%js-1:cg%je-1, cg%ks+1:cg%ke+1) + cg%q(soln)%arr(cg%is:cg%ie,     cg%js+1:cg%je+1, cg%ks+1:cg%ke+1)) * Lyz
         else
            ! 2D and 1D runs do not need too much optimization as they will be rarely used in production runs with selfgravity
            p3 => cg%q(def)%span(cg%ijkse)
            p3 = (1.-src_lapl*2*dom%eff_dim)*cg%q(src)%span(cg%ijkse) + L0*cg%q(soln)%span(cg%ijkse)
            if (dom%has_dir(xdim)) p3 = p3 + &
                 src_lapl*(cg%q(src)%span(cg%ijkse-idm(xdim,:,:)) + cg%q(src)%span(cg%ijkse+idm(xdim,:,:))) &
                 + Lx * (cg%q(soln)%span(cg%ijkse-idm(xdim,:,:)) + cg%q(soln)%span(cg%ijkse+idm(xdim,:,:)))
            if (dom%has_dir(ydim)) p3 = p3 + &
                 src_lapl*(cg%q(src)%span(cg%ijkse-idm(ydim,:,:)) + cg%q(src)%span(cg%ijkse+idm(ydim,:,:))) &
                 + Ly * (cg%q(soln)%span(cg%ijkse-idm(ydim,:,:)) + cg%q(soln)%span(cg%ijkse+idm(ydim,:,:))) &
                 + Lxy * (cg%q(soln)%span(cg%ijkse-idm(xdim,:,:)-idm(ydim,:,:)) + cg%q(soln)%span(cg%ijkse-idm(xdim,:,:)+idm(ydim,:,:)) + &
                 &        cg%q(soln)%span(cg%ijkse+idm(xdim,:,:)-idm(ydim,:,:)) + cg%q(soln)%span(cg%ijkse+idm(xdim,:,:)+idm(ydim,:,:)) )
            if (dom%has_dir(zdim)) p3 = p3 + &
                 src_lapl*(cg%q(src)%span(cg%ijkse-idm(zdim,:,:)) + cg%q(src)%span(cg%ijkse+idm(zdim,:,:))) &
                    + Lz * (cg%q(soln)%span(cg%ijkse-idm(zdim,:,:)) + cg%q(soln)%span(cg%ijkse+idm(zdim,:,:))) &
                 + Lxz * (cg%q(soln)%span(cg%ijkse-idm(xdim,:,:)-idm(zdim,:,:)) + cg%q(soln)%span(cg%ijkse-idm(xdim,:,:)+idm(zdim,:,:)) + &
                 &        cg%q(soln)%span(cg%ijkse+idm(xdim,:,:)-idm(zdim,:,:)) + cg%q(soln)%span(cg%ijkse+idm(xdim,:,:)+idm(zdim,:,:)) ) &
                 + Lyz * (cg%q(soln)%span(cg%ijkse-idm(zdim,:,:)-idm(ydim,:,:)) + cg%q(soln)%span(cg%ijkse-idm(zdim,:,:)+idm(ydim,:,:)) + &
                 &        cg%q(soln)%span(cg%ijkse+idm(zdim,:,:)-idm(ydim,:,:)) + cg%q(soln)%span(cg%ijkse+idm(zdim,:,:)+idm(ydim,:,:)) )
         endif
         cgl => cgl%nxt
      enddo

   end subroutine residual_Mehrstellen

!> \brief Relaxation for Laplace operator implemented in residual_Mehrstellen

   subroutine approximate_solution_relax4M(curl, src, soln, nsmoo)

      use cg_level_coarsest,  only: coarsest
      use cg_level_connected, only: cg_level_connected_T
      use cg_list,            only: cg_list_element
      use cg_list_dataop,     only: dirty_label
      use constants,          only: xdim, ydim, zdim, ndims, GEO_XYZ, BND_NEGREF, pMAX
      use dataio_pub,         only: die, warn
      use domain,             only: dom
      use global,             only: dirty_debug
      use grid_cont,          only: grid_container
      use mpisetup,           only: piernik_MPI_Allreduce, master
      use multigridvars,      only: multidim_code_3D, set_relax_boundaries, coarsest_tol, nc_growth, copy_and_max
      use named_array_list,   only: qna

      implicit none

      type(cg_level_connected_T), pointer, intent(in) :: curl  !< pointer to a level for which we approximate the solution
      integer(kind=4),                     intent(in) :: src   !< index of source in cg%q(:)
      integer(kind=4),                     intent(in) :: soln  !< index of solution in cg%q(:)
      integer(kind=4),                     intent(in) :: nsmoo !< number of smoothing repetitions

      integer :: n
      integer(kind=4) :: b
      real    :: L0, Lx, Ly, Lz, Lxy, Lxz, Lyz, iL0
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
      integer :: is, ie, js, je, ks, ke
      logical :: need_all_bnd_upd
      real :: max_in, max_out
      integer :: ncheck

      if (dom%geometry_type /= GEO_XYZ) call die("[multigrid_Laplace4M:approximate_solution_relax4M] Relaxation for Mehrstellen not implemented for noncartesian grid")
!     call curl%arr3d_boundaries(src, bnd_type = BND_NEGREF) ! required when we use 7-point source term, not just 1-point
      ! Also required when we want to eliminate some communication of soln at a cost of expanding relaxated area into guardcells

      ncheck = 2*dom%nb ! first check for sonvergence of relaxation on coarsest level will be done at this n
      max_in = 0.

      ! Cannot use Red-Black for 4th order Mehrstellen relaxation due to data dependencies even if in some cases Red-Black gives better convergence.
      !> \todo try 4- or 8-color scheme.
      if (dom%nb > 1) call curl%arr3d_boundaries(src, bnd_type = BND_NEGREF)
      do n = 1, nsmoo
         need_all_bnd_upd = (mod(n-1, int(dom%nb)) == 0)
         if (need_all_bnd_upd) call curl%arr3d_boundaries(soln, bnd_type = BND_NEGREF)
         b = int(dom%nb - 1 - mod(n-1, int(dom%nb)), kind=4)

         if (dirty_debug) then
            write(dirty_label, '(a,i5)')"relax4M soln- smoo=", n
            call curl%check_dirty(soln, dirty_label)
         endif

         if (associated(curl, coarsest%level) .and. n==ncheck) then
            max_in = copy_and_max(curl, soln)
            max_out = 0.
         end if

         cgl => curl%first
         do while (associated(cgl))
            cg => cgl%cg

            Lxy = 0. ; if (dom%has_dir(xdim) .and. dom%has_dir(ydim)) Lxy = (cg%idx2 + cg%idy2) / 12.
            Lxz = 0. ; if (dom%has_dir(xdim) .and. dom%has_dir(zdim)) Lxz = (cg%idx2 + cg%idz2) / 12.
            Lyz = 0. ; if (dom%has_dir(ydim) .and. dom%has_dir(zdim)) Lyz = (cg%idy2 + cg%idz2) / 12.
            Lx  = 0. ; if (dom%has_dir(xdim)) Lx = cg%idx2 - 2. * (Lxy + Lxz)
            Ly  = 0. ; if (dom%has_dir(ydim)) Ly = cg%idy2 - 2. * (Lxy + Lyz)
            Lz  = 0. ; if (dom%has_dir(zdim)) Lz = cg%idz2 - 2. * (Lxz + Lyz)
            L0  = 2. * (Lx + Ly + Lz) + 4. * (Lxy + Lxz + Lyz)
            if (L0 /= 0.0) then
               iL0 = 1. / L0
               Lx = Lx * iL0
               Ly = Ly * iL0
               Lz = Lz * iL0
               Lxy = Lxy * iL0
               Lxz = Lxz * iL0
               Lyz = Lyz * iL0
            else
               iL0 = 0. !should never happen but the compiler complains otherwise
            endif
            call set_relax_boundaries(cg, soln, is, ie, js, je, ks, ke, b, .not. need_all_bnd_upd)

            if (dom%eff_dim == ndims .and. .not. multidim_code_3D) then
               ! Set multidim_code_3D to .true. if you want to see performance difference between these two variants of relaxation.
               ! Expect approximately 10-20% difference of the computational cost in favour of the 3D implementation
               ! OPT: In AMR applications (multiple blocks on same process) the most costly part of this routine can be call curl%arr3d_boundaries(soln),
               ! not the relaxation itself
               cg%wa(is  :ie  , js  :je  ,   ks:ke) = &
                    &  (cg%q(soln)%arr(is-1:ie-1, js  :je  , ks  :ke  ) + cg%q(soln)%arr(is+1:ie+1, js  :je  , ks  :ke  ))*Lx      &
                    +  (cg%q(soln)%arr(is  :ie  , js-1:je-1, ks  :ke  ) + cg%q(soln)%arr(is  :ie  , js+1:je+1, ks  :ke  ))*Ly      &
                    +  (cg%q(soln)%arr(is  :ie  , js  :je  , ks-1:ke-1) + cg%q(soln)%arr(is  :ie  , js  :je  , ks+1:ke+1))*Lz      &
                    + ((cg%q(soln)%arr(is-1:ie-1, js-1:je-1, ks  :ke  ) + cg%q(soln)%arr(is+1:ie+1, js-1:je-1, ks  :ke  ))         &
                    +  (cg%q(soln)%arr(is-1:ie-1, js+1:je+1, ks  :ke  ) + cg%q(soln)%arr(is+1:ie+1, js+1:je+1, ks  :ke  )))*Lxy    &
                    + ((cg%q(soln)%arr(is  :ie  , js-1:je-1, ks-1:ke-1) + cg%q(soln)%arr(is  :ie  , js+1:je+1, ks-1:ke-1))         &
                    +  (cg%q(soln)%arr(is  :ie  , js-1:je-1, ks+1:ke+1) + cg%q(soln)%arr(is  :ie  , js+1:je+1, ks+1:ke+1)))*Lyz    &
                    + ((cg%q(soln)%arr(is-1:ie-1, js  :je  , ks-1:ke-1) + cg%q(soln)%arr(is+1:ie+1, js  :je  , ks-1:ke-1))         &
                    +  (cg%q(soln)%arr(is-1:ie-1, js  :je  , ks+1:ke+1) + cg%q(soln)%arr(is+1:ie+1, js  :je  , ks+1:ke+1)))*Lxz    &
                    -   cg%q(src )%arr(is  :ie  , js  :je  , ks  :ke  ) * iL0
! For some weird reasons the formula that comes directly from the Mehrstellen operator worsens convergence
! Note that it requires enabling call curl%arr3d_boundaries(src, ...) above
!!$                 - ((12-2*dom%eff_dim)*cg%q(src)%arr (is  :ie  , j,   k  )                                     &
!!$                 &  + cg%q(src)%arr (is-1:ie-1, j,   k  ) + cg%q(src)%arr (is+1:ie+1, j,   k  )          &
!!$                 &  + cg%q(src)%arr (is  :ie  , j-1, k  ) + cg%q(src)%arr (is  :ie  , j+1, k  )          &
!!$                 &  + cg%q(src)%arr (is  :ie  , j,   k-1) + cg%q(src)%arr (is  :ie  , j,   k+1))/(12. * L0)
            else
               cg%wa(is  :ie  , js:je,   ks:ke) = &
                    -     cg%q(src )%arr(is  :ie  , js  :je,   ks  :ke  ) * iL0
               if (dom%has_dir(xdim)) &
                    &     cg%wa         (is  :ie  , js  :je,   ks  :ke  ) = cg%wa         (is  :ie  , js  :je,   ks  :ke  )          &
                    & +  (cg%q(soln)%arr(is-1:ie-1, js  :je,   ks  :ke  ) + cg%q(soln)%arr(is+1:ie+1, js  :je,   ks  :ke  ))*Lx
               if (dom%has_dir(ydim)) &
                    &     cg%wa         (is  :ie  , js  :je,   ks  :ke  ) = cg%wa         (is  :ie  , js  :je,   ks  :ke  )          &
                    & +  (cg%q(soln)%arr(is  :ie  , js-1:je-1, ks  :ke  ) + cg%q(soln)%arr(is  :ie  , js+1:je+1, ks  :ke  ))*Ly
               if (dom%has_dir(zdim)) &
                    &     cg%wa         (is  :ie  , js  :je,   ks  :ke  ) = cg%wa         (is  :ie  , js  :je,   ks  :ke  )          &
                    & +  (cg%q(soln)%arr(is  :ie  , js  :je,   ks-1:ke-1) + cg%q(soln)%arr(is  :ie  , js  :je,   ks+1:ke+1))*Lz
               if (dom%has_dir(xdim) .and. dom%has_dir(ydim)) &
                    &     cg%wa         (is  :ie  , js  :je,   ks  :ke  ) = cg%wa         (is  :ie  , js  :je,   ks  :ke  )          &
                    & + ((cg%q(soln)%arr(is-1:ie-1, js-1:je-1, ks  :ke  ) + cg%q(soln)%arr(is+1:ie+1, js-1:je-1, ks  :ke  ))         &
                    & + ( cg%q(soln)%arr(is-1:ie-1, js+1:je+1, ks  :ke  ) + cg%q(soln)%arr(is+1:ie+1, js+1:je+1, ks  :ke  )))*Lxy
               if (dom%has_dir(ydim) .and. dom%has_dir(zdim)) &
                    &     cg%wa         (is  :ie  , js  :je,   ks  :ke  ) = cg%wa         (is  :ie  , js  :je,   ks  :ke  )          &
                    & + ((cg%q(soln)%arr(is  :ie  , js-1:je-1, ks-1:ke-1) + cg%q(soln)%arr(is  :ie  , js+1:je+1, ks-1:ke-1))         &
                    & +  (cg%q(soln)%arr(is  :ie  , js-1:je-1, ks+1:ke+1) + cg%q(soln)%arr(is  :ie  , js+1:je+1, ks+1:ke+1)))*Lyz
               if (dom%has_dir(xdim) .and. dom%has_dir(zdim)) &
                    &     cg%wa         (is  :ie  , js  :je,   ks  :ke  ) = cg%wa         (is  :ie  , js  :je,   ks  :ke  )          &
                    & + ((cg%q(soln)%arr(is-1:ie-1, js  :je,   ks-1:ke-1) + cg%q(soln)%arr(is+1:ie+1, js  :je,   ks-1:ke-1))         &
                    & +  (cg%q(soln)%arr(is-1:ie-1, js  :je,   ks+1:ke+1) + cg%q(soln)%arr(is+1:ie+1, js  :je,   ks+1:ke+1)))*Lxz
            endif

            if (associated(curl, coarsest%level) .and. n == ncheck) &
                 max_out = max(max_out, maxval(abs(cg%prolong_xyz( cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) - cg%wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke))))

            cgl => cgl%nxt
         enddo
         call curl%q_copy(qna%wai, soln)

         if (dirty_debug) then
            write(dirty_label, '(a,i5)')"relax4M soln+ smoo=", n
            call curl%check_dirty(soln, dirty_label)
         endif

         if (associated(curl, coarsest%level) .and. n == ncheck) then
            call piernik_MPI_Allreduce(max_out, pMAX)
            if (coarsest_tol*max_in-max_out > 0.) exit
            ncheck = int(ncheck * nc_growth)
         end if

      enddo

      if (associated(curl, coarsest%level) .and. n > nsmoo .and. master) &
           call warn("[multigrid_Laplace4M:approximate_solution_relax4M] relaxation on coarsest level did not converge, consider increasing nsmoob")

   end subroutine approximate_solution_relax4M

end module multigrid_Laplace4M
