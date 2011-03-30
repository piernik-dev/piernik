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

!!$ ============================================================================
!>
!! \brief Variables and data structures required by multigrid routines.
!!
!! \details These variables are not meant to be accessed (or worse: altered) from the outside of multigrid routines (Fortran has no friends :-( )
!<

module multigridvars
! pulled by MULTIGRID

   use constants, only: xdim, zdim, ndims, LO, HI, BND, BLK
   use grid, only: grid_container
   use types, only: domain_container, segment

   implicit none

   public ! QA_WARN no secrets are kept here
   private :: xdim, zdim, ndims, LO, HI, BND, BLK ! QA_WARN prevent re-exporting

   ! multigrid constants
   enum, bind(C)
      enumerator :: source = 1                                        !< Index of the density field
      enumerator :: solution                                          !< Index of the iterated solution (potential) fields
      enumerator :: defect                                            !< Index of the defect field (effectively the density not accounted in current solution)
      enumerator :: correction                                        !< Index of the correction to the potential to be applied at the end of V-cycle
   end enum

   integer, parameter :: level_min = 1                                !< Base (coarsest) level number
   integer, parameter :: level_gb = level_min-1                       !< Global-base level number (will be obsolete soon)
   integer, parameter :: mg_nb = 2                                    !< Number of guardcells in multigrid (simplest laplacian and relaxation require only 1)

   ! these constants should be moved to constants module

   ! namelist parameters
   integer            :: level_max                                    !< Levels of multigrid refinement
   integer            :: ord_prolong                                  !< Prolongation operator order; allowed values are -4 -2, 0 (default), 2 and 4; -2 is often fast
   integer            :: ord_prolong_face                             !< Face prolongation operator order; allowed values are -2 .. 2
   logical            :: stdout                                       !< print verbose messages to stdout
   logical            :: verbose_vcycle                               !< Print one line of log per V-cycle, summary otherwise

   ! dimensions
   integer                                 :: ngridvars               !< number of variables required for implementation of multigrid

   ! boundaries
   enum, bind(C)                                                      !< constants for enumerating multigrid boundary types
      enumerator :: bnd_periodic                                      !< periodic
      enumerator :: bnd_dirichlet                                     !< 0-value boundary type (uniform Dirichlet)
      enumerator :: bnd_isolated                                      !< isolated boundary type
      enumerator :: bnd_neumann                                       !< 0-gradient boundary type (uniform Neumann)
      enumerator :: bnd_givenval                                      !< given value boundary type (general Dirichlet)
      enumerator :: bnd_invalid = bnd_periodic - 1                    !< invalid
   end enum
   logical, dimension(xdim:zdim, LO:HI) :: is_external                !< .true. for non-"mpi" local domain boundaries
   enum, bind(C)
      enumerator :: extbnd_donothing                                  !< Do not touch external boundaries
      enumerator :: extbnd_zero                                       !< Fill external boundaries with zeroes
      enumerator :: extbnd_extrapolate                                !< Perform extrapolation in external boundaries
      enumerator :: extbnd_mirror                                     !< Zero-gradient, mirroring external boundaries
      enumerator :: extbnd_antimirror = - extbnd_mirror               !< mirroring external boundaries with opposite sign
   end enum

   type :: cart   ! auxiliary type for rank-to-coordinates array
      integer, dimension(ndims) :: lo, up, proc
   end type cart
   type(cart), dimension(:),       allocatable :: gb_cartmap          !< rank-to-coordinates array and ranges on gb_src for assembling base level;

   ! miscellaneous
   real                    :: ts                                      !< time for runtime profiling
   real                    :: tot_ts                                  !< total multigrid time

   integer, parameter :: prefix_len = 3                               !< length of prefix for distinguishing V-cycles in the log
   type :: vcycle_stats
      real, allocatable, dimension(:) :: factor                       !< norm reduction factor
      real, allocatable, dimension(:) :: time                         !< time spent
      integer                         :: count                        !< number of executed V-cycles
      real                            :: norm_rhs                     !< norm of the source
      real                            :: norm_final                   !< norm of the defect relative to the source
      character(len=prefix_len)       :: cprefix                      !< prefix for distinguishing V-cycles in the log (e.g inner or outer potential, CR component)
   end type vcycle_stats

   type, extends(segment) :: lvl_segment
     type(plvl), pointer :: nextgrid
   end type lvl_segment

   type, extends(grid_container) :: plvl                              !< single level container

      ! storage
      real, allocatable, dimension(:,:,:,:) :: mgvar                  !< main working array
      real, allocatable, dimension(:,:,:)   :: prolong_x, prolong_xy  !< auxiliary prolongation arrays
      real, allocatable, dimension(:,:,:)   :: bnd_x, bnd_y, bnd_z    !< given boundary values for potential; \todo consider converting it to bnd(:,:,:,xdim:zdim)

      ! geometrical factors, cell counters, etc.
      integer :: level                                                !< multigrid level, level_min == 1: coarsest, level_max: finest
      real    :: vol                                                  !< processor domain volume; BEWARE: for cylindrical geometry multiply by appropriate x(:) to get real volume
      real    :: dxy, dxz, dyz                                        !< cell surface area
      real    :: idx2, idy2, idz2                                     !< inverse of d{x,y,z} square
      real    :: dvol2                                                !< square of cell volume
      real    :: r, rx, ry, rz                                        !< geometric factors for relaxation (diffusion) used in approximate_solution_rbgs

      ! MPI datatype shortcut, similar to grid_container%mbc(ARR, :, :, :)
      integer, dimension(xdim:zdim, LO:HI, BND:BLK, mg_nb) :: mmbc    !< Multigrid MPI Boundary conditions Container for block boundary exchanges with 1 .. mg_nb layers

      type(lvl_segment) :: i_rst, o_rst                       !< description of incoming and outgoing restriction data (this will be a linked list)

      ! data for FFT solver
      integer                                :: nxc                   !< first index (complex or real: fft(:,:,:) or fftr(:,:,:)) cell count
      integer                                :: fft_type              !< type of FFT to employ (depending on boundaries)
      complex, allocatable, dimension(:,:,:) :: fft                   !< a complex array for FFT operations (Fourier space)
      real,    allocatable, dimension(:,:,:) :: fftr                  !< a real array for FFT operations (Fourier space for sine transform)
      real,    allocatable, dimension(:,:,:) :: src                   !< an input array for FFT (real space data)
      real,    allocatable, dimension(:,:,:) :: Green3D               !< Green's function (0.5 * gb_fft_norm / (kx + ky + kz))
      integer (kind = selected_int_kind(16)) :: planf, plani          !< FFT forward and inverse plans
      real                                   :: fft_norm              !< normalization factor

    contains

      procedure :: restrict_level

   end type plvl

   type(plvl), dimension(:), allocatable, target :: lvl               !< a stack of multigrid arrays
   type(plvl), pointer                           :: base              !< pointer to coarsest level
   type(plvl), pointer                           :: roof              !< pointer to finest level
   type(plvl), pointer                           :: gb                !< pointer to global-base level
   type(domain_container), dimension(:), allocatable, target :: dom_lvl       !< a stack of domains with various resolutions

contains

!!$ ============================================================================
!>
!! \brief Simplest restriction (averaging).
!! \todo implement high order restriction and test its influence on V-cycle convergence rate
!<

   subroutine restrict_level(this, iv)

      use constants,  only: xdim, ydim, zdim, LO, HI
      use dataio_pub, only: msg, warn, die
      use grid,       only: D_x, D_y, D_z
      use mpisetup,   only: proc

      implicit none

      class(plvl), intent(inout), target  :: this
      integer, intent(in)      :: iv

      class(plvl), pointer :: coarse
      integer(kind=8), dimension(:,:), pointer :: fse, cse ! shortcuts for fine segment and coarse segment

      if (iv < lbound(this%mgvar(:,:,:,:), dim=4) .or. iv > ubound(this%mgvar(:,:,:,:), dim=4)) call die("[multigridvars:restrict_level] Invalid variable index.")

      coarse => this%o_rst%nextgrid
      if (.not. associated(coarse)) then
         write(msg,'(a,i3)')"[multigridvars:restrict_level] no coarse level here: ", this%level
         call warn(msg) ! can't restrict base level
         return
      endif

!!$      call check_dirty(this%level, iv, "restrict_level-")

      fse => this%o_rst%se
      if (this%o_rst%proc == proc) then
      !> \deprecated BEWARE: unoptimized: some cells are used multiple times (1D and 2D speed-ups possible). Normalization factor will be: / ((1.+D_x)*(1.+D_y)*(1.+D_z))

         cse => coarse%i_rst%se
         coarse%mgvar(     cse(xdim, LO)    :cse(xdim, HI),             cse(ydim, LO)    :cse(ydim, HI),             cse(zdim, LO)    :cse(zdim, HI),             iv) = &
              ( this%mgvar(fse(xdim, LO)    :fse(xdim, HI)-D_x:(1+D_x), fse(ydim, LO)    :fse(ydim, HI)-D_y:(1+D_y), fse(zdim, LO)    :fse(zdim, HI)-D_z:(1+D_z), iv) + &
              & this%mgvar(fse(xdim, LO)+D_x:fse(xdim, HI)    :(1+D_x), fse(ydim, LO)    :fse(ydim, HI)-D_y:(1+D_y), fse(zdim, LO)    :fse(zdim, HI)-D_z:(1+D_z), iv) + &
              & this%mgvar(fse(xdim, LO)    :fse(xdim, HI)-D_x:(1+D_x), fse(ydim, LO)+D_y:fse(ydim, HI)    :(1+D_y), fse(zdim, LO)    :fse(zdim, HI)-D_z:(1+D_z), iv) + &
              & this%mgvar(fse(xdim, LO)+D_x:fse(xdim, HI)    :(1+D_x), fse(ydim, LO)+D_y:fse(ydim, HI)    :(1+D_y), fse(zdim, LO)    :fse(zdim, HI)-D_z:(1+D_z), iv) + &
              & this%mgvar(fse(xdim, LO)    :fse(xdim, HI)-D_x:(1+D_x), fse(ydim, LO)    :fse(ydim, HI)-D_y:(1+D_y), fse(zdim, LO)+D_z:fse(zdim, HI)    :(1+D_z), iv) + &
              & this%mgvar(fse(xdim, LO)+D_x:fse(xdim, HI)    :(1+D_x), fse(ydim, LO)    :fse(ydim, HI)-D_y:(1+D_y), fse(zdim, LO)+D_z:fse(zdim, HI)    :(1+D_z), iv) + &
              & this%mgvar(fse(xdim, LO)    :fse(xdim, HI)-D_x:(1+D_x), fse(ydim, LO)+D_y:fse(ydim, HI)    :(1+D_y), fse(zdim, LO)+D_z:fse(zdim, HI)    :(1+D_z), iv) + &
              & this%mgvar(fse(xdim, LO)+D_x:fse(xdim, HI)    :(1+D_x), fse(ydim, LO)+D_y:fse(ydim, HI)    :(1+D_y), fse(zdim, LO)+D_z:fse(zdim, HI)    :(1+D_z), iv) ) * 0.125
         !\todo add geometrical terms to improve convergence on cylindrical grids
      else
         call die("[multigridvars:restrict_level] cross-processor restriction not implemented yet")
      endif

!!$      call check_dirty(coarse%level, iv, "restrict_level+")

   end subroutine restrict_level

end module multigridvars
