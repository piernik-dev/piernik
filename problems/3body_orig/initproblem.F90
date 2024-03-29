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
module initproblem

! Initial condition for Sedov-Taylor explosion
! Written by: M. Hanasz, March 2006

   implicit none

   private
   public  :: read_problem_par, problem_initial_conditions, problem_pointers

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
   subroutine problem_initial_conditions

      use cg_leaves,      only: leaves
      use cg_list,        only: cg_list_element
      use constants,      only: xdim, ydim, zdim, LO, HI, I_ONE, I_TWO, I_THREE
      use dataio_pub,     only: printinfo
      use fluidindex,     only: flind
      use particle_utils, only: add_part_in_proper_cg

      implicit none

      integer                         :: i, j, k, p
      type(cg_list_element),  pointer :: cgl
      logical, save                   :: first_run = .true.

      do p = lbound(flind%all_fluids, dim=1), ubound(flind%all_fluids, dim=1)
         cgl => leaves%first
         do while (associated(cgl))
            associate(cg => cgl%cg)
               do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
                  do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
                     do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
                        associate( fl => flind%all_fluids(p)%fl )
                           cg%u(fl%idn,i,j,k) = 1.0
                           cg%u(fl%imx,i,j,k) = 0.0
                           cg%u(fl%imy,i,j,k) = 0.0
                           cg%u(fl%imz,i,j,k) = 0.0
                        end associate
                     enddo
                  enddo
               enddo
            end associate
            cgl => cgl%nxt
         enddo
      enddo

      if (first_run) then
         call add_part_in_proper_cg(I_ONE,   1.0, [ 0.9700436,  -0.24308753,  0.0], [ 0.466203685,  0.43236573, 0.0], [0.0, 0.0, 0.0], 0.0)
         call add_part_in_proper_cg(I_TWO,   1.0, [-0.9700436,   0.24308753,  0.0], [ 0.466203685,  0.43236573, 0.0], [0.0, 0.0, 0.0], 0.0)
         call add_part_in_proper_cg(I_THREE, 1.0, [ 0.0,         0.0,         0.0], [-0.932407370, -0.86473146, 0.0], [0.0, 0.0, 0.0], 0.0)
         call printinfo('To see results type: gnuplot -p -e ''plot "nbody_out.log" u 2:3'' ')
         first_run = .false.
      endif

   end subroutine problem_initial_conditions
!-----------------------------------------------------------------------------
end module initproblem
