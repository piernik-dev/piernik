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

!> \brief Helper multigrid routines that depend at most on multigridvars and can be used everywhere else in multigrid.


module multigrid_helpers
! pulled by MULTIGRID

   implicit none

   public :: all_dirty, set_relax_boundaries, copy_and_max
   private 

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
      use multigridvars,  only: source, solution, defect, correction

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

!> \brief Copy solution to a temporary place and compute maximum value

   function copy_and_max(curl, soln) result(max_in)

      use cg_level_connected, only: cg_level_connected_T
      use cg_list,            only: cg_list_element
      use constants,          only: pMAX
      use mpisetup,           only: piernik_MPI_Allreduce

      implicit none

      type(cg_level_connected_T), pointer, intent(in) :: curl  !< pointer to a level for which we approximate the solution
      integer(kind=4),                     intent(in) :: soln  !< index of solution in cg%q(:)

      type(cg_list_element), pointer :: cgl
      real :: max_in

      max_in = 0.
      cgl => curl%first
      do while (associated(cgl))
         associate (cg => cgl%cg)
         cg%prolong_xyz(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) = cg%q(soln)%arr(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke)
         max_in = max(max_in, maxval(abs(cg%prolong_xyz( cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke))))
         end associate
         cgl => cgl%nxt
      enddo
      call piernik_MPI_Allreduce(max_in, pMAX)

   end function copy_and_max

end module multigrid_helpers
