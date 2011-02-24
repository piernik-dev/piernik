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
!! \brief Module containing numerical %constants, enumerations, string lengths, etc.
!! \details Module constants contains numerical %constants, enumerations, string lengths, etc.
!! This module should not depend on anything else in the Piernik code and should not contain any functions or subroutines.
!<
module constants

   use iso_fortran_env, only: error_unit, output_unit

   implicit none

   public                                                ! QA_WARN no secrets are kept here

   ! integers and rationals
   real, parameter :: zero       = 0.0                   !< zero
   real, parameter :: one        = 1.0                   !< one
   real, parameter :: two        = 2.0                   !< two
   real, parameter :: half       = 0.5                   !< a half
   real, parameter :: onet       = 1./3.                 !< one third
   real, parameter :: twot       = 2./3.                 !< two thirds
   real, parameter :: oneq       = 1./4.                 !< one fourth
   real, parameter :: thrq       = 3./4.                 !< three fourths

   ! irrational number approximations
   real, parameter :: pi         = 3.141592653589793238  !< Pi (Archimedes' constant)
   real, parameter :: dpi        = 2.*pi                 !< doubled Pi
   real, parameter :: fpi        = 4.*pi                 !< four Pi
   real, parameter :: e          = 2.718281828459045235  !< Napier's constant (base of Natural logarithm)

   ! some numerical representation extrema
   real, parameter :: big        = huge(real(1.0,4))     !< a constant used as the upper limit number
   real, parameter :: big_float  = huge(real(1.0,4))     !< replicated temporarily 'big' for compatibility \todo choose one and convert occurrences of the other one
   real, parameter :: small      = tiny(real(1.0,4))     !< a constant used as the lower limit number
   integer, parameter :: big_int = huge(int(1,4))

   ! dimensions
   integer, parameter :: xdim    = 1                     !< parameter assigned to the x-direction
   integer, parameter :: ydim    = xdim + 1              !< parameter assigned to the y-direction
   integer, parameter :: zdim    = ydim + 1              !< parameter assigned to the z-direction
   integer, parameter :: ndims   = zdim - xdim + 1       !< We live in a 3-dimensional world
   integer, parameter :: LO = 1                          !< index for low (left) boundary
   integer, parameter :: HI = LO + 1                     !< index for high (right) boundary
   ! \todo remove analogous vars from multigridvars

   ! string lengths
   integer, parameter :: cwdlen = 512                    !< allow for quite long CWD
   integer, parameter :: cbuff_len = 32                  !< length for problem parameters
   integer, parameter :: fplen = 24                      !< length of buffer for printed FP or integer number
   integer, parameter :: domlen = 16                     !< should be <= cbuff_len
   integer, parameter :: varlen = 4                      !< length of state variable names in hdf files
   integer, parameter :: bndlen = 4                      !< length of boundary names
   integer, parameter :: idlen  = 3                      !< COMMENT ME

   ! simulation state
   integer, parameter :: PIERNIK_START       = 1                       ! before initialization
   integer, parameter :: PIERNIK_INIT_MPI    = PIERNIK_START       + 1 ! initialized MPI
   integer, parameter :: PIERNIK_INIT_BASE   = PIERNIK_INIT_MPI    + 1 ! initialized most fundamental modules that depend only on MPI: constants, grid, fluids, etc.
   integer, parameter :: PIERNIK_INIT_ARRAYS = PIERNIK_INIT_BASE   + 1 ! initialized arrays
   integer, parameter :: PIERNIK_INIT_IO_IC  = PIERNIK_INIT_ARRAYS + 1 ! initialized all physics
   integer, parameter :: PIERNIK_INITIALIZED = PIERNIK_INIT_IO_IC  + 1 ! initialized I/O and IC, running
   integer, parameter :: PIERNIK_FINISHED    = PIERNIK_INITIALIZED + 1 ! finished simulation
   integer, parameter :: PIERNIK_CLEANUP     = PIERNIK_FINISHED    + 1 ! finished post-simulation computations and I/O

   ! grid geometry type
   integer, parameter :: GEO_XYZ     = 0                 !< cartesian grid with uniform cell spacing
   integer, parameter :: GEO_RPZ     = GEO_XYZ + 1       !< cylindrical grid with uniform spacing
   integer, parameter :: GEO_INVALID = GEO_XYZ - 1       !< non-recognized grid geometry

   ! boundary conditions type
   integer, parameter :: BND_MPI     = 0                 !< internal, processor-processor boundary
   integer, parameter :: BND_PER     = BND_MPI  + 1      !< periodic boudary
   integer, parameter :: BND_REF     = BND_PER  + 1      !< reflecting boundary
   integer, parameter :: BND_OUT     = BND_REF  + 1      !< free boundary
   integer, parameter :: BND_OUTD    = BND_OUT  + 1      !< one-way outflow boundary
   integer, parameter :: BND_OUTH    = BND_OUTD + 1      !< hydrostatic boundary
   integer, parameter :: BND_COR     = BND_OUTH + 1      !< corner boundary
   integer, parameter :: BND_SHE     = BND_COR  + 1      !< shear boundary
   integer, parameter :: BND_INVALID = BND_MPI  - 1      !< non-recognized boundary

   ! misc
   integer, parameter :: stdout = output_unit
   integer, parameter :: stderr = error_unit

   ! \todo: add enums for boundary condition types

end module constants
