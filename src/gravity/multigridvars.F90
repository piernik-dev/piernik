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

module multigridvars

   use mpisetup, only: cbuff_len

   implicit none

   ! multigrid constants
   integer, parameter :: source=1,          solution=source+1         !< Names of density and potential fields
   integer, parameter :: defect=solution+1, correction=defect+1       !< Names of auxiliary fields (related to source and solution, respectively)
   integer, parameter :: ngridvars = correction                       !< number of variables required for implementation of multigrid; here: 4
   integer, parameter :: level_min = 1, level_gb = level_min-1        !< Base (coarsest) level number and global-base level number
   integer, parameter :: mg_nb = 2                                    !< Number of guardcells in multigrid (simplest laplacian and relaxation require only 1)
   integer, parameter :: bnd_periodic=1, bnd_dirichlet=2              !< constants for enumerating multgrid boundary types: periodic, 0-value,
   integer, parameter :: bnd_isolated=3, bnd_neumann=4                !< isolated, 0-gradient
   integer, parameter :: bnd_givenval=5, bnd_invalid=-1               !< given value, invalid
   integer, parameter :: fft_none=-1, fft_rcr=1, fft_dst=fft_rcr+1    !< type of FFT transform: none, full, discrete sine

   ! these constants can be #defined and used in other source files as well
   integer, parameter :: LOW=1, HIGH=LOW+1                            !< indices for low and high boundary values (third index of plvl%bnd_[xyz] array)
   integer, parameter :: XDIR=1, YDIR=XDIR+1, ZDIR=YDIR+1, NDIM=ZDIR  !< directional indices and number of dimensions
   integer, parameter :: XLO=1,     XHI=XLO+1                         !< enumerated indices for low and high x boundary (used for is_external(:))
   integer, parameter :: YLO=XHI+1, YHI=YLO+1                         !< as above, low and high y
   integer, parameter :: ZLO=YHI+1, ZHI=ZLO+1                         !< as above, low and high y

   ! namelist parameters
   real               :: norm_tol                                     !< stop V-cycle iterations when the ratio of norms ||residual||/||source|| is below this value
   real               :: overrelax                                    !< overrealaxation factor (if < 1. then works as underrelaxation), provided for tests
   real               :: overrelax_x, overrelax_y, overrelax_z        !< overrealaxation factors for fine tuning convergence ratio when cell spacing is not equal in all 3 directions. Use with care, parience and lots of hope.
   real               :: Jacobi_damp                                  !< omega factor for damped Jacobi relaxation. Jacobi_damp == 1 gives undamped method. Try 0.5 in 1D.
   real               :: vcycle_abort                                 !< abort the V-cycle when lhs norm raises by this factor
   real               :: L4_strength                                  !< strength of the 4th order terms in the Laplace operator; 0.: 2nd, 1.: 4th direct, 0.5: 4th integral
   integer            :: level_max                                    !< Levels of multigrid refinement
   integer            :: max_cycles                                   !< Maximum allowed number ov V-cycles
   integer            :: nsmool                                       !< smoothing cycles per call
   integer            :: nsmoob                                       !< smoothing cycles on base level when gb_no_fft = .true. (a convergence check would be much better)
   integer            :: nsmoof                                       !< FFT iterations per call
   integer            :: ord_laplacian                                !< Laplace operator order; allowed values are 2 (default) and 4 (experimental, not fully implemented)
   integer            :: ord_prolong                                  !< Prolongation operator order; allowed values are -4 -2, 0 (default), 2 and 4; -2 is often fast
   integer            :: ord_prolong_face                             !< Face prolongation operator order; allowed values are -2 .. 2
   integer            :: ord_time_extrap                              !< Order of temporal extrapolation for solution recycling; -1 means 0-guess, 2 does parabolic interpolation
   logical            :: trust_fft_solution                           !< Bypas the V-cycle, when doing FFT on whole domain, make sure first that FFT is properly set up.
   logical            :: stdout                                       !< print verbose messages to stdout
   logical            :: verbose_vcycle                               !< Print one line of log per V-cycle, summary otherwise
   logical            :: gb_no_fft                                    !< Deny solving the base level with FFT. Can be very slow.
   logical            :: prefer_rbgs_relaxation                       !< Prefer relaxation over FFT local solver. Typically faster.
   ! \todo allow to perform one or more V-cycles with FFT method, the switch to the RBGS (may save one V-cycle in some cases)
   logical            :: fft_full_relax                               !< Perform full or boundary relaxation after local FFT solve
   logical            :: prefer_modified_norm                         !< Use norm of the unmodified source
   logical            :: gb_solve_gather                              !< Prefer MPI_Gather over Send/Recv when solving global base level (looks a bit faster on small domains)
   logical            :: fft_patient                                  !< Spend more time in init_multigrid to find faster fft plan
   logical            :: hdf5levels                                   !< Dump mgvar to the HDF5 file?
   character(len=cbuff_len) :: grav_bnd_str                                 !< Type of gravitational boundary conditions.

   ! single level container
   type :: plvl

      ! storage
      real, allocatable, dimension(:,:,:,:) :: mgvar                  !< main working array
      real, allocatable, dimension(:,:,:)   :: prolong_x, prolong_xy  !< auxiliary prolongation arrays
      real, allocatable, dimension(:)       :: x, y, z                !< coords of cell centers
      real, allocatable, dimension(:,:,:)   :: bnd_x, bnd_y, bnd_z    !< given boundary values for potential

      ! geometrical factors, cell counters, etc.
      integer :: nx, ny, nz, nb                                       !< x, y, z cell count (including boundary cells), boundary local cell count
      integer :: level                                                !< multigrid level, level_min == 1: coarsest, level_max: finest
      integer :: nxb, nyb, nzb                                        !< x, y and z local active cell count
      integer :: is, ie, js, je, ks, ke                               !< lowest and highest active cell indices
      real    :: dx, dy, dz, dvol, vol                                !< physical cell sizes, cell volume, processor domain volume
      real    :: dxy, dxz, dyz                                        !< cell surface area
      real    :: idx2, idy2, idz2                                     !< inverse of d{x,y,z} square
      real    :: dvol2                                                !< square of cell volume
      real    :: r, rx, ry, rz                                        !< geometric factors for relaxation (diffusion) used in approximate_solution_rbgs

      ! MPI datatype shortcuts
      integer, dimension(mg_nb) :: MPI_XZ_LEFT_BND, MPI_XZ_RIGHT_BND  !< MPI types for block boundary exchange
      integer, dimension(mg_nb) :: MPI_XZ_LEFT_DOM, MPI_XZ_RIGHT_DOM  !< To exchange ng guardcells use e.g. MPI_XZ_LEFT_BND(ng)
      integer, dimension(mg_nb) :: MPI_XY_LEFT_BND, MPI_XY_RIGHT_BND
      integer, dimension(mg_nb) :: MPI_XY_LEFT_DOM, MPI_XY_RIGHT_DOM
      integer, dimension(mg_nb) :: MPI_YZ_LEFT_BND, MPI_YZ_RIGHT_BND
      integer, dimension(mg_nb) :: MPI_YZ_LEFT_DOM, MPI_YZ_RIGHT_DOM

      ! data for FFT solver
      integer                                :: nxc                   !< first index (complex or real: fft(:,:,:) or fftr(:,:,:)) cell count
      integer                                :: fft_type              !< type of FFT to employ (depending on boundaries)
      complex, allocatable, dimension(:,:,:) :: fft                   !< a complex array for FFT operations (Fourier space)
      real,    allocatable, dimension(:,:,:) :: fftr                  !< a real array for FFT operations (Fourier space for sine transform)
      real,    allocatable, dimension(:,:,:) :: src                   !< an input array for FFT (real space data)
      real,    allocatable, dimension(:,:,:) :: Green3D               !< //Green's function (0.5 * gb_fft_norm / (kx + ky + kz))
      integer (kind = selected_int_kind(16)) :: planf, plani          !< FFT forward and inverse plans
      real                                   :: fft_norm              !< normalization factor

   end type plvl

   type(plvl), dimension(:), allocatable, target :: lvl               !< a stack of multigrid arrays
   type(plvl), pointer                           :: base, roof, gb    !< pointers to coarsest, finest and global-base levels, respectively

   ! dimensions
   integer                                 :: eff_dim=0               !< count number of dimensions (>1 cell in a direction)
   logical, dimension(NDIM)                :: has_dir                 !< Set to .true. when there is more than one cell on base level in particular direction
   integer                                 :: D_x, D_y, D_z              !< set to 1 when diven direction exists and use to construct dimensionally-safe indices for arrays

   ! global base-level FFT solver
   type cart   ! auxiliary type for rank-to-coordinates array
      integer, dimension(NDIM) :: lo, up, proc
   end type cart
   type(cart), dimension(:),       allocatable :: gb_cartmap          !< rank-to-coordinates array and ranges on gb_src for assembling base level;
   real,       dimension(:,:,:,:), allocatable :: gb_src_temp         !< Storage for collected base level if using gb_solve_gather

   ! constants from fftw3.f
   integer, parameter :: FFTW_MEASURE=0, FFTW_PATIENT=32, FFTW_ESTIMATE=64
   integer, parameter :: FFTW_RODFT01=8, FFTW_RODFT10=9

   integer            :: fftw_flags = FFTW_MEASURE                    !< or FFTW_PATIENT on request

   ! solution recycling
   integer, parameter :: nold_max=3                                   !< maximum implemented extrapolationorder
   integer :: nold                                                    !< number of old solutions kept for solution recycling
   type old_soln                                                      !< container for an old solution with its timestamp
      real, dimension(:,:,:), allocatable :: soln
      real :: time
   end type old_soln
   type soln_history                                                  !< container for a set of several old potential solutions
      type(old_soln), dimension(nold_max) :: old
      integer :: last                                                 !< index of the last stored potential
      logical :: valid                                                !< .true. when old(last) was properly initialized
   end type soln_history
   type(soln_history), target :: inner, outer                         !< storage for recycling the inner and outer potentials

   ! boundaries
   integer                     :: grav_bnd                            !< boundary type for computational domain
   logical, dimension(XLO:ZHI) :: is_external                         !< .true. for non-"mpi" local domain boundaries when gravity boundary is nonperiodic (even if the global domain is periodic)

   ! miscellaneous
   real                    :: ts, tot_ts                              !< time for runtime profiling, total multigrid time
   real                    :: norm_rhs_orig                           !< Norm of the unmodified source
   character(len=2)        :: cprefix                                 !< optional prefix for distinguishing inner and outer potential V-cycles in the log

   real, dimension(:,:), allocatable :: vcycle_factors                !< buffer for storing V-cycle convergence factors and execution times

end module multigridvars
