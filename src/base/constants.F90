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

   ! precision
   integer, parameter :: LONG = selected_int_kind(16)    ! it is much more clear to type 0_LONG than 0_8
   ! \todo define 1- and 2-byte integers and short (4-byte), long double (10-byte) or quad_precision (16-byte) reals if needed,

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
   enum, bind(C)
      enumerator :: xdim = 1, ydim, zdim                 !< parameters assigned to the x-, y- and z-direction
   end enum
   integer, parameter :: ndims   = zdim - xdim + 1       !< We live in a 3-dimensional world

   enum, bind(C)
      enumerator :: LO = 1, HI                           !< indices for low (left) and high (right) boundaries
   end enum

   ! string lengths
   integer, parameter :: cwdlen = 512                    !< allow for quite long CWD
   integer, parameter :: cbuff_len = 32                  !< length for problem parameters
   integer, parameter :: fplen = 24                      !< length of buffer for printed FP or integer number
   integer, parameter :: domlen = 16                     !< should be <= cbuff_len
   integer, parameter :: varlen = 4                      !< length of state variable names in hdf files
   integer, parameter :: idlen  = 3                      !< COMMENT ME

   ! simulation state
   enum, bind(C)
      enumerator :: PIERNIK_START                        ! before initialization
      enumerator :: PIERNIK_INIT_MPI                     ! initialized MPI
      enumerator :: PIERNIK_INIT_BASE                    ! initialized most fundamental modules that depend only on MPI: constants, grid, fluids, etc.
      enumerator :: PIERNIK_INIT_ARRAYS                  ! initialized arrays
      enumerator :: PIERNIK_INIT_IO_IC                   ! initialized all physics
      enumerator :: PIERNIK_INITIALIZED                  ! initialized I/O and IC, running
      enumerator :: PIERNIK_FINISHED                     ! finished simulation
      enumerator :: PIERNIK_CLEANUP                      ! finished post-simulation computations and I/O
   end enum

   ! grid geometry type
   enum, bind(C)
       enumerator :: GEO_XYZ, GEO_RPZ                    !< cartesian (0) or cylindrical (1) grid with uniform cell spacing
       enumerator :: GEO_INVALID = GEO_XYZ - 1           !< non-recognized grid geometry (-1)
   end enum

   ! boundary conditions type
   enum, bind(C)
      enumerator :: BND_MPI                              !< internal, processor-processor boundary (by default equal to 0)
      enumerator :: BND_PER, BND_REF                     !< periodic boundary, reflecting boudary
      enumerator :: BND_OUT, BND_OUTD, BND_OUTH          !< free boundary, one-way outflow boundary, hydrostatic boundary
      enumerator :: BND_COR, BND_SHE                     !< corner boundary, shear boundary
      enumerator :: BND_INF                              !< 'inf', COMMENT ME
      enumerator :: BND_USER                             !< user boundaries (provided in read_problem_par)
      enumerator :: BND_INVALID = BND_MPI  - 1           !< non-recognized boundary
   end enum

   ! first index of cg%mbc(:,:,:) array
   enum, bind(C)
      enumerator :: FLUID = 1, MAG, ARR                  !< MPI container type for exchanging u(:,:,:,:),  b(ndims,:,:,:) and a rank-3 arrays
   end enum
   ! last index of cg%mbc(:,:,:) array
   enum, bind(C)
      enumerator :: BND = 1, BLK                         !< receiving and sending area
   end enum

   ! Handling boundary cells in the output
   enum, bind(C)
      enumerator :: AT_NO_B, AT_OUT_B, AT_ALL_B, AT_USER !< No boundary cells, external boundary cells and all boundary cells in the output, and user defined area type
   end enum

   ! Fluid type index, used in flind%tag
   enum, bind(C)
      enumerator :: ION = 1, NEU, DST
   end enum

   ! Domain decompositions
   enum, bind(C)
      enumerator :: DD_CART, DD_UE                       !< cartesian and uneven domain decompositions
   end enum

   ! misc
   enum, bind(C)
      enumerator :: MINL, MAXL                           !< constants for mpifind
   end enum
   integer, parameter :: stdout = output_unit
   integer, parameter :: stderr = error_unit
   integer, parameter :: INVALID = -1

   ! \todo: add enums for boundary condition types

end module constants
