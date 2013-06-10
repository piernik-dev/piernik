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

!> \brief Multigrid routines and parameters useful for various variants of the solver

module multigrid_gravity_helper
! pulled by MULTIGRID && GRAV

   implicit none

   private
   public :: approximate_solution, &
        &    nsmoob

   ! namelist parameters
   integer(kind=4)    :: nsmoob                                       !< smoothing cycles on coarsest level when cannot use FFT. (a convergence check would be much better)

contains

!>
!! \brief This routine has to find an approximate solution for given source field and implemented differential operator
!!
!! \details First, it initializes the solution. It means that the solution is:
!! * set to 0 for coarsest level,
!! * prolonged from coarser level otherwise.
!! Then a relaxation is called. The number of relaxation passes is determined according to the level.
!!
!! If for some reason it is undesirable to initialize the solution, call approximate_solution_order directly.
!<

   subroutine approximate_solution(curl, src, soln)

      use cg_level_coarsest,  only: coarsest
      use cg_level_connected, only: cg_level_connected_T
      use constants,          only: BND_NEGREF
      use multigridvars,      only: nsmool
      use multigrid_Laplace,  only: approximate_solution_order
!!$      use constants,           only: fft_none
!!$      use multigrid_fftapprox, only: approximate_solution_fft

      implicit none

      type(cg_level_connected_T), pointer, intent(in) :: curl !< pointer to a level for which we approximate the solution
      integer(kind=4),                     intent(in) :: src  !< index of source in cg%q(:)
      integer(kind=4),                     intent(in) :: soln !< index of solution in cg%q(:)

      integer(kind=4) :: nsmoo

      call curl%check_dirty(src, "approx_soln src-")
!!$      if (curl%fft_type /= fft_none) then
!!$         call approximate_solution_fft(curl, src, soln)
!!$      else
      if (associated(curl, coarsest%level)) then
         !> \todo Implement automatic convergence check on coarsest level (not very important when we have a FFT solver for coarsest level)
         !> \todo Implement alternative bottom-solvers
         nsmoo = nsmoob
         call coarsest%level%set_q_value(soln, 0.)
      else
         nsmoo = nsmool
         call curl%coarser%prolong_q_1var(soln, bnd_type = BND_NEGREF)
         !> \warning when this is be incompatible with V-cycle or other scheme, use direct call to approximate_solution_order
      endif

      call approximate_solution_order(curl, src, soln, nsmoo)
!!$      endif

      call curl%check_dirty(soln, "approx_soln soln+")

   end subroutine approximate_solution

end module multigrid_gravity_helper
