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

!>
!! \brief Variables and data structures required by multigrid routines.
!!
!! \details These variables are not meant to be accessed (or worse: altered) from the outside of multigrid routines (Fortran has no friends :-( )
!!
!! The multigrid solver relies typically at most on 2 guardcells, but the impact on performance is not big so we've decided to have it uniform with the rest of Piernik
!<

module multigridvars
! pulled by MULTIGRID

   use constants,   only: dsetnamelen

   implicit none

   public ! QA_WARN no secrets are kept here
   private :: dsetnamelen ! QA_WARN prevent re-exporting

   ! these 4 variables are required for basic use of the multigrid solver
   character(len=dsetnamelen), parameter :: source_n     = "source"     !< density field
   character(len=dsetnamelen), parameter :: solution_n   = "solution"   !< iterated solution (potential) fields
   character(len=dsetnamelen), parameter :: defect_n     = "defect"     !< defect field (effectively the density not accounted in current solution)
   character(len=dsetnamelen), parameter :: correction_n = "correction" !< correction to the potential to be applied at the end of V-cycle
   integer(kind=4) :: source, solution, defect, correction !< indices to the fields described above

   ! namelist parameters
   integer(kind=4)    :: ord_prolong                                  !< Prolongation operator order; allowed values are -4 -2, 0 (default), 2 and 4; -2 is often fast
   integer(kind=4)    :: ord_prolong_face_norm                        !< Face prolongation operator order in the direction normal to the face; allowed values are 0, 1  and 2
   integer(kind=4)    :: ord_prolong_face_par                         !< Face prolongation operator order in the directions parallel to the face; allowed values are -2 .. 2
   logical            :: stdout                                       !< print verbose messages to stdout
   logical            :: verbose_vcycle                               !< Print one line of log per V-cycle, summary otherwise

   ! miscellaneous
   real                    :: ts                                      !< time for runtime profiling
   real                    :: tot_ts                                  !< total multigrid time
   logical                 :: single_base                             !< .true. when the whole base level is located on a single cpu
   logical                 :: fft_full_relax                          !< Perform full or boundary relaxation after local FFT solve
   integer(kind=4)         :: nsmool                                  !< smoothing cycles per call
   integer(kind=4)         :: nsmoof                                  !< FFT iterations per call
   logical                 :: multidim_code_3D                        !< prefer code written for any 1D and 2D configuration even in 3D for benchmarking and debugging
   real                    :: overrelax                               !< overrealaxation factor (if < 1. then works as underrelaxation), provided for tests

   ! boundaries
   enum, bind(C)                                                      !< constants for enumerating multigrid boundary types
      enumerator :: bnd_periodic                                      !< periodic
      enumerator :: bnd_dirichlet                                     !< 0-value boundary type (uniform Dirichlet)
      enumerator :: bnd_isolated                                      !< isolated boundary type
      enumerator :: bnd_neumann                                       !< 0-gradient boundary type (uniform Neumann)
      enumerator :: bnd_givenval                                      !< given value boundary type (general Dirichlet)
      enumerator :: bnd_invalid = bnd_periodic - 1                    !< invalid
   end enum
   integer            :: grav_bnd                                     !< boundary type for computational domain

 contains

!>
!! \brief Put insane FP values into all multigrid working arrays
!!
!! \details If there are any uninitialized values used in the solver under certain circumstances, the dirtyH will most likely propagate and be easily detectable.
!! \deprecated remove this clause as soon as Intel Compiler gets required features and/or bug fixes
!<

   subroutine all_dirty
#if defined(__INTEL_COMPILER)
      use cg_list_bnd,    only: cg_list_bnd_T  ! QA_WARN intel
#endif /* __INTEL_COMPILER */
      use cg_list_global, only: all_cg
      use constants,      only: dirtyH
      use global,         only: dirty_debug

      implicit none

      call all_cg%set_dirty(source)
      call all_cg%set_dirty(solution)
      call all_cg%set_dirty(defect)
      call all_cg%set_dirty(correction)

      if (dirty_debug) call all_cg%reset_boundaries(dirtyH)

   end subroutine all_dirty

!> \brief Take care of boundaries of relaxated grid

   subroutine set_relax_boundaries(cg, ind, is, ie, js, je, ks, ke, b, need_bnd_upd)

      use constants, only: xdim, ydim, zdim, LO, HI
      use domain,    only: dom
      use grid_cont, only: grid_container

      implicit none

      type(grid_container), pointer, intent(inout) :: cg                     !< current grid container
      integer(kind=4),               intent(in)    :: ind                    !< index in cg%q(:)
      integer,                       intent(out)   :: is, ie, js, je, ks, ke !< indices in cg
      integer(kind=4),               intent(in)    :: b                      !< how far we look into boundary layer
      logical,                       intent(in)    :: need_bnd_upd           !< if .true. then update 1 layer of external boundaries

      ! calling curl%external_boundaries(ind, bnd_type = BND_NEGREF) is a bit overkill
      if (cg%ext_bnd(xdim, LO)) then
         is = cg%is
         if (need_bnd_upd) cg%q(ind)%arr(is-1, :, :) = - cg%q(ind)%arr(is, :, :)
      else
         is = cg%is-b*dom%D_(xdim)
      endif
      if (cg%ext_bnd(xdim, HI)) then
         ie = cg%ie
         if (need_bnd_upd) cg%q(ind)%arr(ie+1, :, :) = - cg%q(ind)%arr(ie, :, :)
      else
         ie = cg%ie+b*dom%D_(xdim)
      endif
      if (cg%ext_bnd(ydim, LO)) then
         js = cg%js
         if (need_bnd_upd) cg%q(ind)%arr(:, js-1, :) = - cg%q(ind)%arr(:, js, :)
      else
         js = cg%js-b*dom%D_(ydim)
      endif
      if (cg%ext_bnd(ydim, HI)) then
         je = cg%je
         if (need_bnd_upd) cg%q(ind)%arr(:, je+1, :) = - cg%q(ind)%arr(:, je, :)
      else
         je = cg%je+b*dom%D_(ydim)
      endif
      if (cg%ext_bnd(zdim, LO)) then
         ks = cg%ks
         if (need_bnd_upd) cg%q(ind)%arr(:, :, ks-1) = - cg%q(ind)%arr(:, :, ks)
      else
         ks = cg%ks-b*dom%D_(zdim)
      endif
      if (cg%ext_bnd(zdim, HI)) then
         ke = cg%ke
         if (need_bnd_upd) cg%q(ind)%arr(:, :, ke+1) = - cg%q(ind)%arr(:, :, ke)
      else
         ke = cg%ke+b*dom%D_(zdim)
      endif

   end subroutine set_relax_boundaries

end module multigridvars
