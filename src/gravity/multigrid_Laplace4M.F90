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
   public :: residual_Mehrstellen, approximate_solution_rbgs4M

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
      use constants,     only: ndims, xdim, ydim, zdim, BND_NEGREF, LO, HI
      use dataio_pub,    only: warn
      use domain,        only: dom
      use grid_cont,     only: grid_container
      use mpisetup,      only: master
      use named_array,   only: p3

      implicit none

      class(cg_leaves_T), intent(in) :: cg_llst !< pointer to a level for which we approximate the solution
      integer(kind=4),    intent(in) :: src     !< index of source in cg%q(:)
      integer(kind=4),    intent(in) :: soln    !< index of solution in cg%q(:)
      integer(kind=4),    intent(in) :: def     !< index of defect in cg%q(:)

      real                :: L0, Lx, Ly, Lz, Lxy, Lxz, Lyz

      logical, save       :: firstcall = .true.
      integer             :: i
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
      integer(kind=4), dimension(ndims,ndims,LO:HI) :: idm

      if (firstcall) then
         if (master) call warn("[multigrid_Laplace4M:residual_Mehrstellen] residual order 4 is experimental.")
         firstcall = .false.
      endif

      call cg_llst%arr3d_boundaries(soln, bnd_type = BND_NEGREF)
      call cg_llst%arr3d_boundaries(src,  bnd_type = BND_NEGREF)

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
         Lx  = 0. ; if (dom%has_dir(xdim)) Lx = cg%idx2 - 2. * (Lxy + Lxz)
         Ly  = 0. ; if (dom%has_dir(ydim)) Ly = cg%idy2 - 2. * (Lxy + Lyz)
         Lz  = 0. ; if (dom%has_dir(zdim)) Lz = cg%idz2 - 2. * (Lxz + Lyz)
         L0  = -2. * (Lx + Ly + Lz) - 4. * (Lxy + Lxz + Lyz)

         p3 => cg%q(def)%span(cg%ijkse)

         p3 = (12.-2*dom%eff_dim)/12.*cg%q(src)%span(cg%ijkse) - L0*cg%q(soln)%span(cg%ijkse)
         if (dom%has_dir(xdim)) p3 = p3 + &
              (cg%q(src)%span(cg%ijkse-idm(xdim,:,:)) + cg%q(src)%span(cg%ijkse+idm(xdim,:,:)))/12. &
              - Lx * (cg%q(soln)%span(cg%ijkse-idm(xdim,:,:)) + cg%q(soln)%span(cg%ijkse+idm(xdim,:,:)))
         if (dom%has_dir(ydim)) p3 = p3 + &
              (cg%q(src)%span(cg%ijkse-idm(ydim,:,:)) + cg%q(src)%span(cg%ijkse+idm(ydim,:,:)))/12. &
              - Ly * (cg%q(soln)%span(cg%ijkse-idm(ydim,:,:)) + cg%q(soln)%span(cg%ijkse+idm(ydim,:,:))) &
              - Lxy * (cg%q(soln)%span(cg%ijkse-idm(xdim,:,:)-idm(ydim,:,:)) + cg%q(soln)%span(cg%ijkse-idm(xdim,:,:)+idm(ydim,:,:)) + &
              &        cg%q(soln)%span(cg%ijkse+idm(xdim,:,:)-idm(ydim,:,:)) + cg%q(soln)%span(cg%ijkse+idm(xdim,:,:)+idm(ydim,:,:)) )
         if (dom%has_dir(zdim)) p3 = p3 + &
              (cg%q(src)%span(cg%ijkse-idm(zdim,:,:)) + cg%q(src)%span(cg%ijkse+idm(zdim,:,:)))/12. &
              - Lz * (cg%q(soln)%span(cg%ijkse-idm(zdim,:,:)) + cg%q(soln)%span(cg%ijkse+idm(zdim,:,:))) &
              - Lxz * (cg%q(soln)%span(cg%ijkse-idm(xdim,:,:)-idm(zdim,:,:)) + cg%q(soln)%span(cg%ijkse-idm(xdim,:,:)+idm(zdim,:,:)) + &
              &        cg%q(soln)%span(cg%ijkse+idm(xdim,:,:)-idm(zdim,:,:)) + cg%q(soln)%span(cg%ijkse+idm(xdim,:,:)+idm(zdim,:,:)) ) &
              - Lyz * (cg%q(soln)%span(cg%ijkse-idm(zdim,:,:)-idm(ydim,:,:)) + cg%q(soln)%span(cg%ijkse-idm(zdim,:,:)+idm(ydim,:,:)) + &
              &        cg%q(soln)%span(cg%ijkse+idm(zdim,:,:)-idm(ydim,:,:)) + cg%q(soln)%span(cg%ijkse+idm(zdim,:,:)+idm(ydim,:,:)) )
         ! OPT: Perhaps in 3D this would be more efficient without ifs.
         ! OPT: Try direct operations on cg%q(soln)%arr without calling span
         ! 2D and 1D runs do not need too much optimization as they will be rarely used in production runs with selfgravity

         cgl => cgl%nxt
      enddo

   end subroutine residual_Mehrstellen

!> \brief Relaxation for Laplace operator implemented in residual_Mehrstellen

   subroutine approximate_solution_rbgs4M(curl, src, soln, nsmoo)

      use cg_level_connected, only: cg_level_connected_T
      use cg_list,            only: cg_list_element
      use cg_list_dataop,     only: dirty_label
      use constants,          only: xdim, ydim, zdim, ndims, GEO_RPZ, I_ONE, BND_NEGREF
      use dataio_pub,         only: die
      use domain,             only: dom
      use global,             only: dirty_debug
      use grid_cont,          only: grid_container
      use multigridvars,      only: multidim_code_3D
      use named_array_list,   only: qna

      implicit none

      type(cg_level_connected_T), pointer, intent(in) :: curl !< pointer to a level for which we approximate the solution
      integer(kind=4),                     intent(in) :: src  !< index of source in cg%q(:)
      integer(kind=4),                     intent(in) :: soln !< index of solution in cg%q(:)
      integer,                             intent(in) :: nsmoo !< number of smoothing repetitions

      integer :: n
      real    :: L0, Lx, Ly, Lz, Lxy, Lxz, Lyz
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg

!     call curl%arr3d_boundaries(src, bnd_type = BND_NEGREF) ! required when we use 7-point source term, not just 1-point
      ! Also required when we want to eliminate some communication of soln at a cost of expanding relaxated area into guardcells

      ! Cannot use Red-Black for 4th order Mehrstellen relaxation due to data dependencies even if in some cases Red-Black gives better convergence
      do n = 1, nsmoo
         call curl%arr3d_boundaries(soln, bnd_type = BND_NEGREF)

         if (dirty_debug) then
            write(dirty_label, '(a,i5)')"relax4M soln- smoo=", n
            call curl%check_dirty(soln, dirty_label, expand=I_ONE)
         endif
         cgl => curl%first
         do while (associated(cgl))
            cg => cgl%cg

            if (dom%geometry_type == GEO_RPZ) call die("[multigrid_Laplace4M:approximate_solution_rbgs4M] Relaxation for Mehrstellen not implemented for noncartesian grid")
            Lxy = 0. ; if (dom%has_dir(xdim) .and. dom%has_dir(ydim)) Lxy = (cg%idx2 + cg%idy2) / 12.
            Lxz = 0. ; if (dom%has_dir(xdim) .and. dom%has_dir(zdim)) Lxz = (cg%idx2 + cg%idz2) / 12.
            Lyz = 0. ; if (dom%has_dir(ydim) .and. dom%has_dir(zdim)) Lyz = (cg%idy2 + cg%idz2) / 12.
            Lx  = 0. ; if (dom%has_dir(xdim)) Lx = cg%idx2 - 2. * (Lxy + Lxz)
            Ly  = 0. ; if (dom%has_dir(ydim)) Ly = cg%idy2 - 2. * (Lxy + Lyz)
            Lz  = 0. ; if (dom%has_dir(zdim)) Lz = cg%idz2 - 2. * (Lxz + Lyz)
            L0  = 2. * (Lx + Ly + Lz) + 4. * (Lxy + Lxz + Lyz)
            if (dom%eff_dim == ndims .and. .not. multidim_code_3D) then
               ! Set multidim_code_3D to .true. if you want to see performance difference between these two variants of relaxation.
               ! Expect approximately 10-20% difference of the computational cost in favour of the 3D implementation
               ! OPT: In AMR applications (multiple blocks on same process) the most costly part of this routine can be call curl%arr3d_boundaries(soln),
               ! not the relaxation itself
               cg%wa(cg%is  :cg%ie  , cg%js  :cg%je  ,   cg%ks:cg%ke) = &
                    &  (cg%q(soln)%arr(cg%is-1:cg%ie-1, cg%js  :cg%je  , cg%ks  :cg%ke  ) + cg%q(soln)%arr(cg%is+1:cg%ie+1, cg%js  :cg%je  , cg%ks  :cg%ke  ))*Lx/L0   &
                    +  (cg%q(soln)%arr(cg%is  :cg%ie  , cg%js-1:cg%je-1, cg%ks  :cg%ke  ) + cg%q(soln)%arr(cg%is  :cg%ie  , cg%js+1:cg%je+1, cg%ks  :cg%ke  ))*Ly/L0   &
                    +  (cg%q(soln)%arr(cg%is  :cg%ie  , cg%js  :cg%je  , cg%ks-1:cg%ke-1) + cg%q(soln)%arr(cg%is  :cg%ie  , cg%js  :cg%je  , cg%ks+1:cg%ke+1))*Lz/L0   &
                    + ((cg%q(soln)%arr(cg%is-1:cg%ie-1, cg%js-1:cg%je-1, cg%ks  :cg%ke  ) + cg%q(soln)%arr(cg%is+1:cg%ie+1, cg%js-1:cg%je-1, cg%ks  :cg%ke  ))         &
                    +  (cg%q(soln)%arr(cg%is-1:cg%ie-1, cg%js+1:cg%je+1, cg%ks  :cg%ke  ) + cg%q(soln)%arr(cg%is+1:cg%ie+1, cg%js+1:cg%je+1, cg%ks  :cg%ke  )))*Lxy/L0 &
                    + ((cg%q(soln)%arr(cg%is  :cg%ie  , cg%js-1:cg%je-1, cg%ks-1:cg%ke-1) + cg%q(soln)%arr(cg%is  :cg%ie  , cg%js+1:cg%je+1, cg%ks-1:cg%ke-1))         &
                    +  (cg%q(soln)%arr(cg%is  :cg%ie  , cg%js-1:cg%je-1, cg%ks+1:cg%ke+1) + cg%q(soln)%arr(cg%is  :cg%ie  , cg%js+1:cg%je+1, cg%ks+1:cg%ke+1)))*Lyz/L0 &
                    + ((cg%q(soln)%arr(cg%is-1:cg%ie-1, cg%js  :cg%je  , cg%ks-1:cg%ke-1) + cg%q(soln)%arr(cg%is+1:cg%ie+1, cg%js  :cg%je  , cg%ks-1:cg%ke-1))         &
                    +  (cg%q(soln)%arr(cg%is-1:cg%ie-1, cg%js  :cg%je  , cg%ks+1:cg%ke+1) + cg%q(soln)%arr(cg%is+1:cg%ie+1, cg%js  :cg%je  , cg%ks+1:cg%ke+1)))*Lxz/L0 &
                    -   cg%q(src )%arr(cg%is  :cg%ie  , cg%js  :cg%je  , cg%ks  :cg%ke  ) / L0
! For some weird reasons the formula that comes directly from the Mehrstellen operator worsens convergence
! Note that it requires enabling call curl%arr3d_boundaries(src, ...) above
!!$                 - ((12-2*dom%eff_dim)*cg%q(src)%arr (cg%is  :cg%ie  , j,   k  )                                     &
!!$                 &  + cg%q(src)%arr (cg%is-1:cg%ie-1, j,   k  ) + cg%q(src)%arr (cg%is+1:cg%ie+1, j,   k  )          &
!!$                 &  + cg%q(src)%arr (cg%is  :cg%ie  , j-1, k  ) + cg%q(src)%arr (cg%is  :cg%ie  , j+1, k  )          &
!!$                 &  + cg%q(src)%arr (cg%is  :cg%ie  , j,   k-1) + cg%q(src)%arr (cg%is  :cg%ie  , j,   k+1))/(12. * L0)
            else
               cg%wa(cg%is  :cg%ie  , cg%js:cg%je,   cg%ks:cg%ke) = &
                    -     cg%q(src )%arr(cg%is  :cg%ie  , cg%js  :cg%je,   cg%ks  :cg%ke  ) / L0
               if (dom%has_dir(xdim)) &
                    &     cg%wa         (cg%is  :cg%ie  , cg%js  :cg%je,   cg%ks  :cg%ke  ) = cg%wa         (cg%is  :cg%ie  , cg%js  :cg%je,   cg%ks  :cg%ke  )          &
                    & +  (cg%q(soln)%arr(cg%is-1:cg%ie-1, cg%js  :cg%je,   cg%ks  :cg%ke  ) + cg%q(soln)%arr(cg%is+1:cg%ie+1, cg%js  :cg%je,   cg%ks  :cg%ke  ))*Lx/L0
               if (dom%has_dir(ydim)) &
                    &     cg%wa         (cg%is  :cg%ie  , cg%js  :cg%je,   cg%ks  :cg%ke  ) = cg%wa         (cg%is  :cg%ie  , cg%js  :cg%je,   cg%ks  :cg%ke  )          &
                    & +  (cg%q(soln)%arr(cg%is  :cg%ie  , cg%js-1:cg%je-1, cg%ks  :cg%ke  ) + cg%q(soln)%arr(cg%is  :cg%ie  , cg%js+1:cg%je+1, cg%ks  :cg%ke  ))*Ly/L0
               if (dom%has_dir(zdim)) &
                    &     cg%wa         (cg%is  :cg%ie  , cg%js  :cg%je,   cg%ks  :cg%ke  ) = cg%wa         (cg%is  :cg%ie  , cg%js  :cg%je,   cg%ks  :cg%ke  )          &
                    & +  (cg%q(soln)%arr(cg%is  :cg%ie  , cg%js  :cg%je,   cg%ks-1:cg%ke-1) + cg%q(soln)%arr(cg%is  :cg%ie  , cg%js  :cg%je,   cg%ks+1:cg%ke+1))*Lz/L0
               if (dom%has_dir(xdim) .and. dom%has_dir(ydim)) &
                    &     cg%wa         (cg%is  :cg%ie  , cg%js  :cg%je,   cg%ks  :cg%ke  ) = cg%wa         (cg%is  :cg%ie  , cg%js  :cg%je,   cg%ks  :cg%ke  )          &
                    & + ((cg%q(soln)%arr(cg%is-1:cg%ie-1, cg%js-1:cg%je-1, cg%ks  :cg%ke  ) + cg%q(soln)%arr(cg%is+1:cg%ie+1, cg%js-1:cg%je-1, cg%ks  :cg%ke  ))         &
                    & + ( cg%q(soln)%arr(cg%is-1:cg%ie-1, cg%js+1:cg%je+1, cg%ks  :cg%ke  ) + cg%q(soln)%arr(cg%is+1:cg%ie+1, cg%js+1:cg%je+1, cg%ks  :cg%ke  )))*Lxy/L0
               if (dom%has_dir(ydim) .and. dom%has_dir(zdim)) &
                    &     cg%wa         (cg%is  :cg%ie  , cg%js  :cg%je,   cg%ks  :cg%ke  ) = cg%wa         (cg%is  :cg%ie  , cg%js  :cg%je,   cg%ks  :cg%ke  )          &
                    & + ((cg%q(soln)%arr(cg%is  :cg%ie  , cg%js-1:cg%je-1, cg%ks-1:cg%ke-1) + cg%q(soln)%arr(cg%is  :cg%ie  , cg%js+1:cg%je+1, cg%ks-1:cg%ke-1))         &
                    & +  (cg%q(soln)%arr(cg%is  :cg%ie  , cg%js-1:cg%je-1, cg%ks+1:cg%ke+1) + cg%q(soln)%arr(cg%is  :cg%ie  , cg%js+1:cg%je+1, cg%ks+1:cg%ke+1)))*Lyz/L0
               if (dom%has_dir(xdim) .and. dom%has_dir(zdim)) &
                    &     cg%wa         (cg%is  :cg%ie  , cg%js  :cg%je,   cg%ks  :cg%ke  ) = cg%wa         (cg%is  :cg%ie  , cg%js  :cg%je,   cg%ks  :cg%ke  )          &
                    & + ((cg%q(soln)%arr(cg%is-1:cg%ie-1, cg%js  :cg%je,   cg%ks-1:cg%ke-1) + cg%q(soln)%arr(cg%is+1:cg%ie+1, cg%js  :cg%je,   cg%ks-1:cg%ke-1))         &
                    & +  (cg%q(soln)%arr(cg%is-1:cg%ie-1, cg%js  :cg%je,   cg%ks+1:cg%ke+1) + cg%q(soln)%arr(cg%is+1:cg%ie+1, cg%js  :cg%je,   cg%ks+1:cg%ke+1)))*Lxz/L0
            endif
            cgl => cgl%nxt
         enddo
         call curl%q_copy(qna%wai, soln)

         if (dirty_debug) then
            write(dirty_label, '(a,i5)')"relax4M soln+ smoo=", n
            call curl%check_dirty(soln, dirty_label)
         endif
      enddo

   end subroutine approximate_solution_rbgs4M

end module multigrid_Laplace4M
