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

module cresp_variables ! & constants
! pulled by CRESP

   use constants, only: fpi, three

   implicit none

   public                                                ! QA_WARN no secrets are kept here

   real, parameter :: div_v        = 0.0
   real, parameter :: omega_d      = 0.0
   real, parameter :: clight_cresp = 1.0

   real, parameter :: fpcc         = fpi * clight_cresp
   real, parameter :: fpcc2        = fpi * clight_cresp**2
   real, parameter :: fp3cc        = fpi / three * clight_cresp

! these will most probably be in types and will be modified by the driver (piernik)
!
!   integer              :: taylor_coeff_2nd, taylor_coeff_3rd

!   real, parameter      :: cnst_c  = 1.0e0 ! speed of light
!   real                 :: p_lo_d, p_up_d, u_d_d, u_b_d   ! additional variables that are passed down to timestep by grid, but since grid is absent in driver module, they must be compensated for this way.
!   real, parameter      :: cnst_me = 1.0d0 ! mass of electron


end module cresp_variables
