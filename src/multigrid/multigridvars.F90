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

module multigridvars
! pulled by MULTIGRID
   implicit none

   public ! QA_WARN no secrets are kept here

   ! multigrid constants
   integer, parameter :: source=1,          solution=source+1         !< Names of density and potential fields
   integer, parameter :: defect=solution+1, correction=defect+1       !< Names of auxiliary fields (related to source and solution, respectively)
   integer, parameter :: level_min = 1, level_gb = level_min-1        !< Base (coarsest) level number and global-base level number
   integer, parameter :: mg_nb = 2                                    !< Number of guardcells in multigrid (simplest laplacian and relaxation require only 1)

   ! these constants can be #defined and used in other source files as well
   integer, parameter :: LOW=1, HIGH=LOW+1                            !< indices for low and high boundary values (third index of plvl%bnd_[xyz] array)
   integer, parameter :: NDIM=3                                       !< number of dimensions
   integer, parameter :: XLO=1,     XHI=XLO+1                         !< enumerated indices for low and high x boundary (used for is_external(:))
   integer, parameter :: YLO=XHI+1, YHI=YLO+1                         !< as above, low and high y
   integer, parameter :: ZLO=YHI+1, ZHI=ZLO+1                         !< as above, low and high y

   ! namelist parameters
   integer            :: level_max                                    !< Levels of multigrid refinement
   integer            :: ord_prolong                                  !< Prolongation operator order; allowed values are -4 -2, 0 (default), 2 and 4; -2 is often fast
   integer            :: ord_prolong_face                             !< Face prolongation operator order; allowed values are -2 .. 2
   logical            :: stdout                                       !< print verbose messages to stdout
   logical            :: verbose_vcycle                               !< Print one line of log per V-cycle, summary otherwise

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
      real,    allocatable, dimension(:,:,:) :: Green3D               !< Green's function (0.5 * gb_fft_norm / (kx + ky + kz))
      integer (kind = selected_int_kind(16)) :: planf, plani          !< FFT forward and inverse plans
      real                                   :: fft_norm              !< normalization factor

   end type plvl

   type(plvl), dimension(:), allocatable, target :: lvl               !< a stack of multigrid arrays
   type(plvl), pointer                           :: base, roof, gb    !< pointers to coarsest, finest and global-base levels, respectively

   ! dimensions
   integer                                 :: eff_dim                 !< count number of dimensions (>1 cell in a direction)
   integer                                 :: ngridvars               !< number of variables required for implementation of multigrid

   ! boundaries
   integer, parameter :: bnd_periodic=1, bnd_dirichlet=2              !< constants for enumerating multigrid boundary types: periodic, 0-value,
   integer, parameter :: bnd_isolated=3, bnd_neumann=4                !< isolated, 0-gradient
   integer, parameter :: bnd_givenval=5, bnd_invalid=-1               !< given value, invalid
   logical, dimension(XLO:ZHI) :: is_external                         !< .true. for non-"mpi" local domain boundaries
   integer :: periodic_bnd_cnt, non_periodic_bnd_cnt                  !< counts periodic and non-periodic boundaries in existing directions
   integer, parameter :: extbnd_donothing = 0                         !< Do not touch external boundaries
   integer, parameter :: extbnd_zero = extbnd_donothing + 1           !< Fill external boundaries with zeroes
   integer, parameter :: extbnd_extrapolate = extbnd_zero + 1         !< Perform extrapolation in external boundaries
   integer, parameter :: extbnd_mirror = extbnd_extrapolate + 1       !< Zero-gradient, mirroring external boundaries
   integer, parameter :: extbnd_antimirror = - extbnd_mirror          !< mirroring external boundaries with opposite sign

   type cart   ! auxiliary type for rank-to-coordinates array
      integer, dimension(NDIM) :: lo, up, proc
   end type cart
   type(cart), dimension(:),       allocatable :: gb_cartmap          !< rank-to-coordinates array and ranges on gb_src for assembling base level;

   ! miscellaneous
   real                    :: ts, tot_ts                              !< time for runtime profiling, total multigrid time

   integer, parameter :: prefix_len = 3                               !< length of prefix for distinguishing V-cycles in the log
   type :: vcycle_stats
      real, allocatable, dimension(:) :: factor                       !< norm reduction factor
      real, allocatable, dimension(:) :: time                         !< time spent
      integer                         :: count                        !< number of executed V-cycles
      real                            :: norm_rhs                     !< norm of the source
      real                            :: norm_final                   !< norm of the defect relative to the source
      character(len=prefix_len)       :: cprefix                      !< prefix for distinguishing V-cycles in the log (e.g inner or outer potential, CR component)
   end type vcycle_stats

end module multigridvars
