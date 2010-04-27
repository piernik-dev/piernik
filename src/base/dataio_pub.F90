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
!    Initial implemetation of PIERNIK code was based on TVD split MHD code by
!    Ue-Li Pen
!        see: Pen, Arras & Wong (2003) for algorithm and
!             http://www.cita.utoronto.ca/~pen/MHD
!             for original source code "mhd.f90"
!
!    For full list of developers see $PIERNIK_HOME/license/pdt.txt
!
#include "piernik.def"

module dataio_public
   implicit none

   real                  :: tend                   !< simulation time to end
   real                  :: wend                   !< wall clock time to end (in hours)

   integer               :: nend                   !< number of the step to end simulation
   integer               :: nstep_start            !< number of start timestep
   integer               :: nhdf                   !< current number of hdf file

   character(len=128)    :: log_file               !< path to the current log file
   integer               :: log_lun = 3            !< luncher for log file

   logical               :: halfstep = .false.     !< true when X-Y-Z sweeps are done and Z-Y-X are not

end module dataio_public
