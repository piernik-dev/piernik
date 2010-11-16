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

!<
!! This module contains problem variables common to all problems.
!! It removes dependence of dataio_hdf5 on initproblem, which allows to hook user-provided I/O subroutines.
!>

module problem_pub

   use mpisetup, only: cbuff_len
   use types,    only: idlen

   implicit none

   public  ! QA_WARN nothing to hide here

   character(len=cbuff_len) :: problem_name                   !< The default problem name
   character(len=idlen)     :: run_id                         !< Auxiliary run identifier

   ! hack for tests
#ifdef JEANS_PROBLEM
   real :: jeans_d0
   integer :: jeans_mode
#endif /* JEANS_PROBLEM */

end module problem_pub
