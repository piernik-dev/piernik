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
#define RNG is:ie, js:je, ks:ke
module initproblem

! Initial condition for Sedov-Taylor explosion
! Written by: M. Hanasz, March 2006

   implicit none

   private
   public  :: read_problem_par, init_prob, problem_pointers

   integer(kind=4) :: n_sn
   real            :: d0, p0, bx0, by0, bz0, Eexpl, x0, y0, z0, r0, dt_sn, r, t_sn

   namelist /PROBLEM_CONTROL/ d0, p0, bx0, by0, bz0, Eexpl, x0, y0, z0, r0, n_sn, dt_sn

contains

!-----------------------------------------------------------------------------

   subroutine problem_pointers

      implicit none

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      implicit none

   end subroutine read_problem_par
!-----------------------------------------------------------------------------
   subroutine init_prob

      use constants,   only: xdim, ydim, zdim
      use cg_list,     only: cg_list_element
      use cg_leaves,   only: leaves
      use fluidindex,  only: flind
      use fluidtypes,  only: component_fluid
      use grid_cont,   only: grid_container
      use particle_types, only: pset
      use particle_integrators, only: hermit_4ord

      implicit none

      integer                         :: i, j, k, p
      class(component_fluid), pointer :: fl
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg

      do p = lbound(flind%all_fluids, dim=1), ubound(flind%all_fluids, dim=1)
         fl => flind%all_fluids(p)%fl

! Uniform equilibrium state

         cgl => leaves%first
         do while (associated(cgl))
            cg => cgl%cg

            do k = 1, cg%n_(zdim)
               do j = 1, cg%n_(ydim)
                  do i = 1, cg%n_(xdim)
                     cg%u(fl%idn,i,j,k) = 1.0
                     cg%u(fl%imx,i,j,k) = 0.0
                     cg%u(fl%imy,i,j,k) = 0.0
                     cg%u(fl%imz,i,j,k) = 0.0
                  enddo
               enddo
            enddo
            cgl => cgl%nxt
         enddo
      enddo

      call pset%add(1.0, [ 0.9700436,  -0.24308753,  0.0], [ 0.466203685,  0.43236573, 0.0])
      call pset%add(1.0, [-0.9700436,   0.24308753,  0.0], [ 0.466203685,  0.43236573, 0.0])
      call pset%add(1.0, [ 0.0,         0.0,         0.0], [-0.932407370, -0.86473146, 0.0])
      call hermit_4ord(pset, 0.0, 10.0)

   end subroutine init_prob
!-----------------------------------------------------------------------------
end module initproblem
