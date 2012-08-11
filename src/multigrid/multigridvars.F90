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
   integer :: source, solution, defect, correction !< indices to the fields described above

   ! these constants should be moved to constants module

   ! namelist parameters
   integer(kind=4)    :: ord_prolong                                  !< Prolongation operator order; allowed values are -4 -2, 0 (default), 2 and 4; -2 is often fast
   integer(kind=4)    :: ord_prolong_face_norm                        !< Face prolongation operator order in the direction normal to the face; allowed values are 0, 1  and 2
   integer(kind=4)    :: ord_prolong_face_par                         !< Face prolongation operator order in the directions parallel to the face; allowed values are -2 .. 2
   logical            :: stdout                                       !< print verbose messages to stdout
   logical            :: verbose_vcycle                               !< Print one line of log per V-cycle, summary otherwise

   ! miscellaneous
   real                    :: ts                                      !< time for runtime profiling
   real                    :: tot_ts                                  !< total multigrid time
   logical                 :: is_mg_uneven                            !< .true. when domain shapes differ across procesors, even on the coarsest grids
   logical                 :: single_base                             !< .true. when the whole base level is located on a single cpu
   logical                 :: need_general_pf                         !< .false. only for most regular domain decomposition

end module multigridvars
