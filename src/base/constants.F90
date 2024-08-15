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
!! \brief Module containing numerical %constants, enumerations, string lengths, etc.
!! \details Module constants contains numerical %constants, enumerations, string lengths, etc.
!! This module should not depend on anything else in the Piernik code and should not contain any functions or subroutines.
!<
module constants

   use iso_fortran_env, only: error_unit, output_unit

   implicit none

   public                                                ! QA_WARN no secrets are kept here

   ! precision
   integer, parameter :: FP_REAL   = selected_real_kind(5)   ! this should be 32-bit single precision
   integer, parameter :: FP_DOUBLE = selected_real_kind(12)  ! this should be 64-bit double precision
   integer, parameter :: FP_EXT    = selected_real_kind(17)  ! this should be 80-bit extended precision
   integer, parameter :: FP_QUAD   = selected_real_kind(30)  ! this should be 128-bit quad precision (don't expect hardware support in CPU)

   integer, parameter :: INT128 = selected_int_kind(30)      ! this should be 128-bit integer (don't expect hardware support in CPU)
   integer, parameter :: LONG   = selected_int_kind(16)      ! We need at least 8-byte integers to count to 10**16, it is much more clear to type 0_LONG than 0_8
   integer, parameter :: INT4   = selected_int_kind(9)       ! Assume that all MPI and HDF5 calls expect 4-byte integers
   ! \todo define:
   !     1-byte integer (selected_int_kind(1))
   !     2-byte integer (selected_int_kind(3))
   ! if needed

   ! integers and rationals
   real, parameter :: zero       = 0.0                   !< zero
   real, parameter :: one        = 1.0                   !< one
   real, parameter :: two        = 2.0                   !< two
   real, parameter :: three      = 3.0                   !< three
   real, parameter :: four       = 4.0                   !< four
   real, parameter :: five       = 5.0
   real, parameter :: eight      = 8.0                   !< eight
   real, parameter :: ten        = 10.0
   real, parameter :: half       = 0.5                   !< a half
   real, parameter :: onet       = 1./3.                 !< one third
   real, parameter :: twot       = 2./3.                 !< two thirds
   real, parameter :: oneq       = 1./4.                 !< one fourth
   real, parameter :: thrq       = 3./4.                 !< three fourths
   real, parameter :: onesth     = 1./6.                 !< one sixth
   real, parameter :: oneeig     = 1./8.                 !< one eighth
   real, parameter :: logten     = log(ten)              !< natural logarithm of ten

   enum, bind(C)
      enumerator :: idn = 1, imx, imy, imz, ien
   end enum

   enum, bind(C)
      enumerator :: I_ZERO = 0, I_ONE, I_TWO, I_THREE, I_FOUR, I_FIVE, I_SIX, I_SEVEN, I_EIGHT, I_NINE, I_TEN
   end enum

   ! enumerator for length/mass/time/velocity/magnetic field units
   enum, bind(C)
      enumerator :: U_LEN = 1, U_MASS, U_TIME, U_VEL, U_MAG, U_TEMP, U_ENER
   end enum

   ! irrational number approximations
   real, parameter :: pi         = 3.141592653589793238  !< Pi (Archimedes' constant)
   real, parameter :: dpi        = 2.*pi                 !< doubled Pi
   real, parameter :: fpi        = 4.*pi                 !< four Pi
   real, parameter :: e          = 2.718281828459045235  !< Napier's constant (base of Natural logarithm)

   ! some numerical representation extrema
   !> Warning: the valuses of big and small should be used only when exceeding the range od single precision is not desired
   !! e.g. when when using the default, 32-bit data output in .h5 files.
   !! In most cases it should be much safer to use tiny(1.) and huge(1.) (either directly or slightly scaled).
   !! Careless relying on SP huge and tiny may result in incorrect calculations when the involved values exceed single precision range.
   !! Such incorrect calculation may happen e.g. for extremely big or small cell sizes, extreme timestep length etc.
   real, parameter :: big        = huge(real(1.0,4))     !< a constant used as the upper limit number.
   real, parameter :: big_float  = huge(real(1.0,4))     !< replicated temporarily 'big' for compatibility \todo choose one and convert occurrences of the other one
   real, parameter :: dirtyH     = big                   !< If dirty_debug then pollute arrays with this insane value
   real, parameter :: dirtyH1    = 10.**int(log10(big))  !< this "round" dirty value makes it easier to detect which call contaminated the data
   real, parameter :: dirtyH1c   = 0.1 * dirtyH1         !< lowest value of marked dirty pollutions
   real, parameter :: dirtyL     = sqrt(dirtyH)          !< If dirty_debug then assume that the array got contaminated by dirtyH by checking against this value
   real, parameter :: small      = tiny(real(1.0,4))     !< a constant used as the lower limit number
   integer, parameter :: big_int = huge(int(1,4))

   !! Expected values of tiny, huge and epsilon for various lengths of real type (obtained from gfortran)
   !! x=          0._4 (SP)        0._8 (DP, the default)    0._10 (ext precision)          0._16 (quad precision, software implementation)
   !! tiny(x)     1.17549435E-38   2.2250738585072014E-308   3.36210314311209350626E-4932   3.36210314311209350626267781732175260E-4932
   !! huge(x)     3.40282347E+38   1.7976931348623157E+308   1.18973149535723176502E+4932   1.18973149535723176508575932662800702E+4932
   !! epsilon(x)  1.19209290E-07   2.2204460492503131E-016   1.08420217248550443401E-0019   1.92592994438723585305597794258492732E-0034

   ! dimensions
   enum, bind(C)
      enumerator :: xdim = 1, ydim, zdim                 !! parameters assigned to the x-, y- and z-direction
   end enum
   integer(kind=4), parameter :: ndims   = zdim - xdim + 1 !< We live in a 3-dimensional world
   integer(kind=4), parameter :: cor_dim = zdim+1        !< corner "direction", useful for selecting guardcell exchanges
   !>
   !! array of all positive permutations of xyzdim
   !<
   enum, bind(C)
      enumerator :: NORMAL = 1    ! Normal direction
      enumerator :: ORTHO1        ! First orthogonal to the NORMAL (x -> y -> z -> x cycle)
      enumerator :: ORTHO2        ! Second orthogonal to the NORMAL
   end enum
   integer(kind=4), dimension(xdim:zdim, NORMAL:ORTHO2) :: pdims = &
      reshape([xdim, ydim, zdim, ydim, zdim, xdim, zdim, xdim, ydim], [ndims, ndims])

   enum, bind(C)
      enumerator :: LO = 1, HI                           !! indices for low (left) and high (right) boundaries
   end enum

   ! string lengths
   integer(kind=4), parameter :: cwdlen = 512               !< allow for quite long CWD
   integer(kind=4), parameter :: fmt_len = 128              !< length of format string
   integer(kind=4), parameter :: fnamelen = 128             !< length of output filename
   integer(kind=4), parameter :: cbuff_len = 32             !< length for problem parameters
   integer(kind=4), parameter :: msg_len = cbuff_len        !< length for message file paths; ToDo restore cwdlen length and handle it properly
   integer(kind=4), parameter :: units_len = 5 * cbuff_len  !< length for unit strings
   integer(kind=4), parameter :: fplen = 24                 !< length of buffer for printed FP or integer number
   integer(kind=4), parameter :: domlen = 16                !< should be <= cbuff_len
   integer(kind=4), parameter :: dsetnamelen = cbuff_len    !< length of dataset name and state variable names in hdf files
   integer(kind=4), parameter :: idlen = 3                  !< COMMENT ME
   integer(kind=4), parameter :: singlechar = 1             !< a single character

   integer(kind=4), parameter :: rep_len = 112              !< default length of "---...", "+++..." and "===..." strings in the log

   ! simulation state
   enum, bind(C)
      enumerator :: PIERNIK_START                        ! before initialization
      enumerator :: PIERNIK_INIT_MPI                     ! initialized MPI
      enumerator :: PIERNIK_INIT_DOMAIN                  ! initialized domain
      enumerator :: PIERNIK_INIT_GLOBAL                  ! initialized global parameters
      enumerator :: PIERNIK_INIT_FLUIDS                  ! initialized fluid properties
      enumerator :: PIERNIK_INIT_GRID                    ! initialized grids
      enumerator :: PIERNIK_INIT_IO_IC                   ! initialized all physics
      enumerator :: PIERNIK_POST_IC                      ! initial conditions or restart completed
      enumerator :: PIERNIK_INITIALIZED                  ! initialized I/O and IC, running
      enumerator :: PIERNIK_FINISHED                     ! finished simulation
      enumerator :: PIERNIK_CLEANUP                      ! finished post-simulation computations and I/O
   end enum

   ! verbosity levels for printinfo
   enum, bind(C)
      enumerator :: V_SILENT                 ! messages that are supposed to be only in the log, not on stdout
      enumerator :: V_DEBUG                  ! debugging information
      enumerator :: V_VERBOSE                ! extra details
      enumerator :: V_INFO                   ! printed by default
      enumerator :: V_ESSENTIAL              ! to make the output very brief
      enumerator :: V_WARN                   ! redirect to warn()
      ! aliases
      enumerator :: V_LOG     = V_SILENT
      enumerator :: V_NORMAL  = V_INFO
      enumerator :: V_LACONIC = V_ESSENTIAL
      ! special aliases
      enumerator :: V_LOWEST = V_SILENT
      enumerator :: V_HIGHEST = V_WARN
   end enum

   character(len=*), dimension(V_LOWEST:V_HIGHEST), parameter :: v_name = [ &
        "log      ", &
        "debug    ", &
        "verbose  ", &
        "standard ", &
        "essential", &
        "warn     " ]

   ! grid geometry type
   enum, bind(C)
      enumerator :: GEO_XYZ, GEO_RPZ                    !! Cartesian (0) or cylindrical (1) grid with uniform cell spacing
      enumerator :: GEO_INVALID = GEO_XYZ - 1           !! non-recognized grid geometry (-1)
   end enum

   ! boundary conditions type
   enum, bind(C)
      enumerator :: BND_MPI                   !! internal, processor-processor boundary on the same level (by default equal to 0)
      enumerator :: BND_FC                    !! internal, processor-processor boundary with lower refinement level
      enumerator :: BND_MPI_FC                !! internal, processor-processor boundary partially on the same level, partially lower
      enumerator :: BND_PER                   !! periodic boundary
      enumerator :: BND_REF                   !! reflecting boundary
      enumerator :: BND_NEGREF                !! antireflecting boundary
      enumerator :: BND_ZERO                  !! zero-valued boundary
      enumerator :: BND_NONE                  !! do not touch external boundary
      enumerator :: BND_XTRAP                 !! linearly extrapolated boundary
      enumerator :: BND_OUT                   !! free boundary
      enumerator :: BND_OUTD                  !! one-way outflow boundary
      enumerator :: BND_OUTH                  !! hydrostatic boundary
      enumerator :: BND_OUTHD                 !! hydrostatic diode boundary (one-way outflow)
      enumerator :: BND_COR                   !! corner boundary
      enumerator :: BND_SHE                   !! shear boundary
      enumerator :: BND_USER                  !! user boundaries (provided in read_problem_par)
      enumerator :: BND_INVALID = BND_MPI - 1 !! non-recognized boundary
   end enum

   ! solver type
   enum, bind(C)
      enumerator :: RTVD_SPLIT    !! MHD RTVD, as it was implemented from the beginning of Piernik
      enumerator :: HLLC_SPLIT    !! non-magnetic (pure HD) HLLC as first attempt of something more precise than RTVD, lacks many features an ma be removed at some point
      enumerator :: RIEMANN_SPLIT !! MHD Riemann, implementations by Varadarajan Parthasarathy; HD variant is slower than HLLC_SPLIT
   end enum
   ! Perhaps it may make sense to create compatibility matrix for solvers.
   ! AMR, magnetic, FARGO, resistivity, ...


   ! enumerate stages of Runge-Kutta method in an unique way, so istep will contain information both about stage and method
   enum, bind(C)
      enumerator :: EULER = 10000  !! integration_order == 1
      enumerator :: RK2_1          !! halfstep RK2
      enumerator :: RK2_2          !! fullstep RK2
   end enum
   integer, parameter :: n_schemes = 2
   integer, dimension(n_schemes), parameter :: first_stage = [ EULER, RK2_1 ], &
        &                                      last_stage  = [ EULER, RK2_2 ]
   real, dimension(EULER:RK2_2), parameter :: rk_coef = [ one, &      !! EULER
        &                                                 half, one ] !! RK2

   ! 3D and 4D array names
   ! fluids
   character(len=dsetnamelen), parameter :: fluid_n = "fluid"   !< main fluid array
   character(len=dsetnamelen), parameter :: uh_n    = "uh"      !< auxiliary array for half-step values
   ! magnetic field
   character(len=dsetnamelen), parameter :: mag_n   = "mag"     !< main magnetic field array
   character(len=dsetnamelen), parameter :: magh_n  = "magh"    !< auxiliary array for half-step values
   ! gravitational potential
   character(len=dsetnamelen), parameter :: gp_n    = "gp"      !< static, external field, must be explicitly set to 0. if no external fields are applied
   character(len=dsetnamelen), parameter :: sgp_n   = "sgp"     !< current field from self-gravity
   character(len=dsetnamelen), parameter :: sgpm_n  = "sgpm"    !< previous field from self-gravity
   character(len=dsetnamelen), parameter :: gpot_n  = "gpot"    !< current sum of fields
   character(len=dsetnamelen), parameter :: hgpot_n = "hgpot"   !< sum of fields for half-step values
#ifdef NBODY
#ifdef NBODY_GRIDDIRECT
   character(len=dsetnamelen), parameter :: nbgp_n  = "nbgp"    !< current gravity field from nbody
#endif /* NBODY_GRIDDIRECT */
   character(len=dsetnamelen), parameter :: gp1b_n  = "gp1b"    !< current gravity field from 1 body
   character(len=dsetnamelen), parameter :: nbdn_n  = "nbdn"    !< density from particles
   character(len=dsetnamelen), parameter :: prth_n  = "prth"    !< histogram of particles on the grid
#endif /* NBODY */
   ! misc
   character(len=dsetnamelen), parameter :: wcu_n   = "wcu"     !< (resistivity) COMMENT ME
   character(len=dsetnamelen), parameter :: cs_i2_n = "cs_iso2" !< map of imposed isothermal sound speed
   character(len=dsetnamelen), parameter :: wcr_n   = "wcr"     !< auxiliary array for CR diffusion
   character(len=dsetnamelen), parameter :: wa_n    = "wa"      !< general-purpose auxiliary 3D array
   character(len=dsetnamelen), parameter :: psi_n   = "psi"     !< auxiliary 3D array for divergence cleaning
   character(len=dsetnamelen), parameter :: psih_n  = "psih"    !< auxiliary 3D array for divergence cleaning for half-step values

   ! timer names
   character(len=*), parameter :: tmr_fu  = "fluid_update"       !< main timer used to measure fluid_update step
   character(len=*), parameter :: tmr_hdf = "hdf_dump"           !< timer for I/O operations
   character(len=*), parameter :: tmr_mg  = "multigrid"          !< timer for gravity multigrid solver
   character(len=*), parameter :: tmr_mgd = "multigrid_diff"     !< timer for CR diffusion multigrid solver
   character(len=*), parameter :: tmr_amr = "refinement_update"  !< timer for refinement updates

   ! Handling boundary cells in the output (AT stands for Area Type)
   enum, bind(C)
      enumerator :: AT_BACKUP       !! backup field: no output AND no backup
      enumerator :: AT_IGNORE       !! no output for anything less or equal AT_IGNORE
      enumerator :: AT_NO_B         !! no boundary cells
      enumerator :: AT_OUT_B        !! external boundary cells
      enumerator :: AT_USER         !! user defined area type
   end enum

   ! Position of variable within a cell
   enum, bind(C)
      enumerator :: VAR_CENTER      !! cell-centered
      enumerator :: VAR_CORNER      !! corner (staggered-grid variable)
      enumerator :: VAR_XFACE       !! X-face
      enumerator :: VAR_YFACE       !! Y-face
      enumerator :: VAR_ZFACE       !! Z-face
   end enum

   ! Interpolation order
   enum, bind(C)
      enumerator :: O_INJ = 0  !! injection
      enumerator :: O_LIN = 1  !! linear
      enumerator :: O_I2  = 2  !! integral quadratic
      enumerator :: O_I3  = 3  !! integral cubic
      enumerator :: O_I4  = 4  !! integral quartic
      enumerator :: O_I5  = 5  !! integral quintic
      enumerator :: O_I6  = 6  !! integral sextic
      enumerator :: O_D2  = -2 !! direct quadratic
      enumerator :: O_D3  = -3 !! direct cubic
      enumerator :: O_D4  = -4 !! direct quartic
      enumerator :: O_D5  = -5 !! direct quintic
      enumerator :: O_D6  = -6 !! direct sextic
   end enum

   ! Fluid type index, used in flind%tag
   enum, bind(C)
      enumerator :: ION = 1, NEU, DST
   end enum

   ! Output type
   enum, bind(C)
      enumerator :: RES = 1, HDF, LOGF, TSL, INCEPTIVE, FINAL_DUMP, CHK
   end enum

   ! base level
   integer(kind=4), parameter :: base_level_id = 0 !< Base domain level id. Refinements are positively numbered, coarsened levels for use in multigrid solvers have negative numbers.
   integer(kind=4), parameter :: refinement_factor = 2 !< Resolution difference between consecutive levels. This is deeply hardwired into prolongation, restriction and such routines

   ! type of FFT transform used in multigrid
   enum, bind(C)
      enumerator :: fft_rcr  =  1 !< full
      enumerator :: fft_dst       !< discrete sine
      enumerator :: fft_none = -1 !< none
   end enum

   ! particle interpolation scheme
   enum, bind(C)
      enumerator :: I_NGP   ! Nearest grid point
      enumerator :: I_CIC   ! Cloud in cell
      enumerator :: I_TSC   ! Triangular shaped cloud
   end enum

   ! divB=0 constraining method
   enum, bind(C)
      enumerator :: DIVB_CT   ! Constrained Transport
      enumerator :: DIVB_HDC  ! Hyperbolic Divergence Cleaning (div(B) diffusion, GLM)
   end enum
   integer(kind=4), parameter :: psidim = zdim + 1

   ! -1, 0, 1
   enum, bind(C)
      enumerator :: IM = -1
      enumerator :: I0
      enumerator :: IP
   end enum

   ! coordinate positions;
   !> \warning gridgeometry depends on the sequence
   enum, bind(C)
      enumerator :: CENTER
      enumerator :: LEFT
      enumerator :: RIGHT
      enumerator :: INV_CENTER
   end enum

   ! MPI reduction type
   enum, bind(C)
      enumerator :: pSUM = 1, pMIN, pMAX, pLOR, pLAND
   end enum

   ! different velocities used by fargo algorithm
   enum, bind(C)
      enumerator :: VEL_RES, VEL_CR
   end enum

   ! misc
   logical, parameter :: has_B = &
#ifdef MAGNETIC
        .true.
#else /* !MAGNETIC */
        .false.
#endif /* !MAGNETIC */
   enum, bind(C)
      enumerator :: MINL, MAXL                           !< constants for func::get_extremum
      enumerator :: RD, WR                               !< constants for wd_{rd,wr} selection
   end enum
   integer, parameter :: stdout = output_unit
   integer, parameter :: stderr = error_unit
   integer(kind=4), parameter :: INVALID = -I_ONE

   integer(kind=4), dimension(ndims,ndims),       parameter :: idm  = reshape(int([ [1,0,0], [0,1,0], [0,0,1] ], kind=4),[ndims,ndims])   !< identity matrix 3x3
   integer(kind=4), dimension(ndims,ndims,LO:HI), parameter :: idm2 = reshape([idm,idm],[ndims,ndims,2_INT4])                             !< auxiliary matrix 3x3x2 based on identity matrix
   integer(kind=4), dimension(ndims),             parameter :: uv   = int([1,1,1], kind=4)                                                !< unity vector

   ! PPP timer categories
   enum, bind(C)                                      ! Event related to:
      !enumerator :: PPP_MAIN  = int(B"000000000000", kind=4)  ! general, unmasked, always enabled
      enumerator :: PPP_IO    = int(B"000000000001", kind=4)  ! I/O
      enumerator :: PPP_MG    = int(B"000000000010", kind=4)  ! multigrid
      enumerator :: PPP_GRAV  = int(B"000000000100", kind=4)  ! gravity
      enumerator :: PPP_CR    = int(B"000000001000", kind=4)  ! cosmic rays
      enumerator :: PPP_PART  = int(B"000000010000", kind=4)  ! particles
      enumerator :: PPP_MPI   = int(B"000000100000", kind=4)  ! MPI
      enumerator :: PPP_AMR   = int(B"000001000000", kind=4)  ! refinement
      enumerator :: PPP_CG    = int(B"000010000000", kind=4)  ! single cg
      enumerator :: PPP_MAG   = int(B"000100000000", kind=4)  ! magnetic field
      enumerator :: PPP_PROB  = int(B"001000000000", kind=4)  ! problem
      enumerator :: PPP_DEBUG = int(B"010000000000", kind=4)  ! debug
      enumerator :: PPP_AUX   = int(B"100000000000", kind=4)  ! auxiliary (unused by default)
   end enum

   ! OS type
   enum, bind(C)
      enumerator :: LINUX
      enumerator :: APPLE
      !enumerator :: ANDROID
      !enumerator :: WINDOWS
   end enum

end module constants
