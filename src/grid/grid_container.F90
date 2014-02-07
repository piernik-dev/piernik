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

!> \brief Module containing the grid container type and its associated methods

module grid_cont

   use constants,       only: xdim, zdim, ndims, LO, HI, CENTER, INV_CENTER
   use fluxtypes,       only: fluxarray, fluxpoint
   use named_array,     only: named_array4d, named_array3d
   use refinement_flag, only: ref_flag
   use real_vector,     only: real_vec_T

   implicit none

   private
   public :: grid_container, pr_segment, tgt_list, is_overlap, segment

   type(fluxpoint), target :: fpl, fpr, cpl, cpr

   !> \brief Specification of segment of data for boundary exchange
   type :: segment
      integer :: proc                                     !< target process
      integer(kind=8), dimension(xdim:zdim, LO:HI) :: se  !< range
      integer(kind=4) :: tag                              !< unique tag for data exchange
      real, allocatable, dimension(:,:,:)   :: buf        !< buffer for the 3D (scalar) data to be sent or received
      real, allocatable, dimension(:,:,:,:) :: buf4       !< buffer for the 4D (vector) data to be sent or received
      integer(kind=4), pointer :: req                     !< request ID, used for most asynchronous communication, such as fine-coarse flux exchanges
      integer(kind=8), dimension(xdim:zdim, LO:HI) :: se2 !< auxiliary range, used in cg_level_connected:vertical_bf_prep
      type(grid_container), pointer :: local              !< set this pointer to non-null when the exchange is local
   end type segment

   !> \brief coefficient-layer pair used for prolongation
   type :: c_layer
      integer :: layer                                                !< index of a layer with face-prolongation coefficient coeff
      real    :: coeff                                                !< coefficient for face prolongation
   end type c_layer

   !> \brief segment type for prolongation and restriction
   type, extends(segment) :: pr_segment
      type(c_layer), dimension(:), allocatable :: f_lay               !< face layers to contribute to the prolonged face value
   end type pr_segment

   !< \brief target list container for prolongations, restrictions and boundary exchanges
   type :: tgt_list
      type(segment), dimension(:), allocatable :: seg              !< a segment of data to be received or sent
   end type tgt_list

   !>
   !! \brief Multigrid-specific storage
   !!
   !! \details The FFT-related arrays here cannot be declared as named arrays in cg%q(:) because their size or type differ
   !!
   !! \todo Provide an initialization (procedure pointer) in multigrid gravity for these things when a cg is dynamically added
   !<
   type :: mg_arr
      ! storage
      real, allocatable, dimension(:,:,:)   :: bnd_x, bnd_y, bnd_z    !< given boundary values for potential; \todo make an array of pointers, indexed by (xdim:zdim)

      ! data for FFT multigrid self-gravity solver
      integer                                :: nxc                   !< first index (complex or real: fft(:,:,:) or fftr(:,:,:)) cell count
      complex, allocatable, dimension(:,:,:) :: fft                   !< a complex array for FFT operations (Fourier space)
      real,    allocatable, dimension(:,:,:) :: fftr                  !< a real array for FFT operations (Fourier space for sine transform)
      real,    allocatable, dimension(:,:,:) :: src                   !< an input array for FFT (real space data)
      real,    allocatable, dimension(:,:,:) :: Green3D               !< Green's function (0.5 * fft_norm / (kx + ky + kz))
      integer (kind = selected_int_kind(16)) :: planf, plani          !< FFT forward and inverse plans
      real                                   :: fft_norm              !< normalization factor

      ! geometrical factors
      !! \todo move to cg
      real    :: r, rx, ry, rz                                        !< geometric factors for relaxation (diffusion) used in approximate_solution_{rbgs,relax}*

   end type mg_arr

   !> \brief Array of boundary segments to exchange
   type :: bnd_list
      type(segment), dimension(:), allocatable :: seg !< segments
    contains
       procedure :: add_seg !< Add an new segment, reallocate if necessary
   end type bnd_list

   !> \brief Everything required for autonomous computation of a single sweep on a portion of the domain on a single process
   type :: grid_container

      ! Cell properties

      ! sizes
      real, dimension(ndims) :: dl                               !< array of %grid cell sizes in all directions, [ dx, dy, dz ]
      real :: dx                                                 !< length of the %grid cell in x-direction
      real :: dy                                                 !< length of the %grid cell in y-direction
      real :: dz                                                 !< length of the %grid cell in z-direction
      real :: dvol                                               !< volume of one %grid cell
      real :: dxy, dxz, dyz                                      !< cell surface area

      ! shortcuts
      real, dimension(ndims) :: idl                              !< array of inverted %grid cell sizes in all directions, 1./dl(:)
      real :: idx                                                !< inverted length of the %grid cell in x-direction
      real :: idy                                                !< inverted length of the %grid cell in y-direction
      real :: idz                                                !< inverted length of the %grid cell in z-direction
      real :: idx2, idy2, idz2                                   !< inverse of d{x,y,z} square
      real, dimension(ndims) :: idl2                             !< [ idx2, idy2, idz2 ]

      ! Grid properties

      ! cell count and position
      integer(kind=4), dimension(ndims) :: n_b                   !< [nxb, nyb, nzb]
      integer(kind=4) :: nxb                                     !< number of %grid cells in one block (without boundary cells) in x-direction
      integer(kind=4) :: nyb                                     !< number of %grid cells in one block (without boundary cells) in y-direction
      integer(kind=4) :: nzb                                     !< number of %grid cells in one block (without boundary cells) in z-direction
      integer(kind=4) :: level_id                                !< level id (number); do not use it without a good reason, use cg_level_T%lev where possible instead

      ! shortcuts
      !> \todo Change kind from 4 to 8 to allow really deep refinements (effective resolution > 2**31, perhaps the other requirement will be default integer  kind = 8)
      integer(kind=4) :: is                                      !< index of the first %grid cell of physical domain in x-direction
      integer(kind=4) :: ie                                      !< index of the last %grid cell of physical domain in x-direction
      integer(kind=4) :: js                                      !< index of the first %grid cell of physical domain in y-direction
      integer(kind=4) :: je                                      !< index of the last %grid cell of physical domain in y-direction
      integer(kind=4) :: ks                                      !< index of the first %grid cell of physical domain in z-direction
      integer(kind=4) :: ke                                      !< index of the last %grid cell of physical domain in z-direction
      integer(kind=4) :: isb, ieb, jsb, jeb, ksb, keb            !< auxiliary indices for exchanging boundary data, (e.g. is:isb -> ie+1:nx, ieb:ie -> 1:nb)
!      integer(kind=4) :: il1, ih1, jl1, jh1, kl1, kh1            !< index of the first guardcell, adjacent to the active interior of the grid, l for low/left, h for high/right side
!      integer(kind=4) :: iln, ihn, jln, jhn, kln, khn            !< index of the n-th guardcell, furthest from the active interior of the grid, l for low/left, h for high/right side
      integer(kind=4), dimension(ndims, LO:HI)  :: ijkse         !< [[is,  js,  ks ], [ie,  je,  ke ]]
      integer(kind=4), dimension(ndims, LO:HI)  :: ijkseb        !< [[isb, jsb, ksb], [ieb, jeb, keb]]
      integer(kind=4), dimension(ndims, LO:HI)  :: lh1           !< [[il1, jl1, kl1], [ih1, jh1, kh1]]
      integer(kind=4), dimension(ndims, LO:HI)  :: lhn           !< [[iln, jln, kln], [ihn, jhn, khn]]
      integer(kind=8), dimension(ndims)         :: h_cor1        !< offsets of the corner opposite to the one defined by off(:) + 1, a shortcut to be compared with dom%n_d(:) DEPRECATED will be equivalent to ijkse(:, HI)+1
      integer(kind=4), dimension(ndims)         :: n_            !< number of %grid cells in one block in x-, y- and z-directions (n_b(:) + 2 * nb)
      integer(kind=8), dimension(ndims, LO:HI) :: my_se          !< own segment. my_se(:,LO) = 0; my_se(:,HI) = dom%n_d(:) - 1 would cover entire domain on a base level
                                                                 !! my_se(:,LO) = 0; my_se(:,HI) = finest%level%n_d(:) -1 would cover entire domain on the most refined level
                                                                 !! DEPRECATED: will be equivalent to ijkse(:,:)
      integer(kind=8), dimension(ndims) :: level_off             !< offset of own level

      ! Physical size and coordinates

      real, dimension(ndims, LO:HI) :: fbnd                      !< current block boundary positions

      type(real_vec_T), dimension(CENTER:INV_CENTER, ndims) :: coord !< all coordinates (CENTER, LEFT, RIGHT, INV_CENTER)
      ! shortcuts
      real, pointer, dimension(:) :: x                             !< array of x-positions of %grid cells centers
      real, pointer, dimension(:) :: y                             !< array of x-positions of %grid cells centers
      real, pointer, dimension(:) :: z                             !< array of x-positions of %grid cells centers
      real, pointer, dimension(:) :: inv_x                         !< array of invert x-positions of %grid cells centers
      real, pointer, dimension(:) :: inv_y                         !< array of invert y-positions of %grid cells centers
      real, pointer, dimension(:) :: inv_z                         !< array of invert z-positions of %grid cells centers

      ! External boundary conditions and internal boundaries

      integer(kind=4), dimension(ndims, LO:HI)           :: bnd         !< type of boundary conditions coded in integers
      type(bnd_list),  dimension(:),         allocatable :: i_bnd       !< description of incoming boundary data, the shape is (xdim:zdim)
      type(bnd_list),  dimension(:),         allocatable :: o_bnd       !< description of outgoing boundary data, the shape is (xdim:zdim)
      logical,         dimension(xdim:zdim, LO:HI)       :: ext_bnd     !< .false. for BND_PER and BND_MPI
      type(fluxarray), dimension(xdim:zdim, LO:HI)       :: finebnd     !< indices and flux arrays for fine/coarse flux updates on coarse side
      type(fluxarray), dimension(xdim:zdim, LO:HI)       :: coarsebnd   !< indices and flux arrays for fine/coarse flux updates on fine side

      ! Prolongation and restriction

      type(tgt_list) :: ri_tgt                                    !< description of incoming restriction data (this should be a linked list)
      type(tgt_list) :: ro_tgt                                    !< description of outgoing restriction data
      type(tgt_list) :: pi_tgt                                    !< description of incoming prolongation data
      type(tgt_list) :: po_tgt                                    !< description of outgoing prolongation data
      type(tgt_list) :: pib_tgt                                   !< description of incoming boundary prolongation data
      type(tgt_list) :: pob_tgt                                   !< description of outgoing boundary prolongation data
      type(tgt_list) :: rif_tgt                                   !< description of fluxes incoming from fine grid
      type(tgt_list) :: rof_tgt                                   !< description of fluxes outgoing to coarse grid
      real, allocatable, dimension(:,:,:) :: prolong_, prolong_x, prolong_xy !< auxiliary prolongation arrays for intermediate results
      real, dimension(:,:,:), pointer ::  prolong_xyz             !< auxiliary prolongation array for final result. OPT: Valgrind indicates that operations on array allocated on pointer might be slower than on ordinary arrays due to poorer L2 cache utilization
      logical, allocatable, dimension(:,:,:) :: leafmap           !< .true. when a cell is not covered by finer cells, .false. otherwise
      logical, allocatable, dimension(:,:,:) :: refinemap         !< .true. when a cell triggers refinement criteria, .false. otherwise

      ! Non-cartesian geometrical factors

      real, allocatable, dimension(:,:,:) :: gc_xdim             !< array of geometrical coefficients in x-direction
      real, allocatable, dimension(:,:,:) :: gc_ydim             !< array of geometrical coefficients in y-direction
      real, allocatable, dimension(:,:,:) :: gc_zdim             !< array of geometrical coefficients in z-direction

      ! Registered variables

      type(named_array3d), allocatable, dimension(:) :: q        !< 3D arrays such as gravitational potential pr user-defined quantities or gravitational potential
      type(named_array4d), allocatable, dimension(:) :: w        !< 4D arrays such as u, vector fields (b) or other vector/multi-scalar user-defined quantities

      ! handy shortcuts to some entries in q(:)
      real, dimension(:,:,:), pointer :: gpot    => null()       !< Array for sum of gravitational potential at t += dt
      real, dimension(:,:,:), pointer :: hgpot   => null()       !< Array for sum of gravitational potential at t += 0.5*dt
      real, dimension(:,:,:), pointer :: gp      => null()       !< Array for gravitational potential from external fields
      real, dimension(:,:,:), pointer :: sgp     => null()       !< Array for gravitational potential from multigrid or FFT solver
      real, dimension(:,:,:), pointer :: sgpm    => null()       !< Array for gravitational potential from multigrid or FFT solver at previous timestep saved by source_terms_grav.
      real, dimension(:,:,:), pointer :: cs_iso2 => null()       !< COMMENT ME
      real, dimension(:,:,:), pointer :: wa      => null()       !< Temporary array used for different purposes, usually has dimension (grid::nx, grid::ny, grid::nz)

      ! handy shortcuts to some entries in w(:)
      real, dimension(:,:,:,:), pointer :: u     => null()       !< Main array of all fluids' components
      real, dimension(:,:,:,:), pointer :: b     => null()       !< Main array of magnetic field's components

      ! Misc
      type(mg_arr), pointer :: mg                                !< multigrid arrays
      real :: vol                                                !< volume of the grid; BEWARE: for cylindrical geometry it needs to be integrated over x(:) to get real volume
      real :: dxmn                                               !< the smallest length of the %grid cell (among dx, dy, and dz)
      integer(kind=4) :: maxxyz                                  !< maximum number of %grid cells in any direction
      integer :: grid_id                                         !< index of own segment in own level decomposition, e.g. my_se(:,:) = base%level%dot%gse(proc)%c(grid_id)%se(:,:)
      integer(kind=8) :: SFC_id                                  !< position of the grid on space-filling curve
      type(ref_flag) :: refine_flags                             !< refine or derefine this grid container?
      integer :: membership                                      !< How many cg lists use this grid piece?
      logical :: ignore_prolongation                             !< When .true. do not upgrade interior with incoming prolonged values
      logical :: is_old                                          !< .true. if a given grid existed prior to  upgrade_refinement call
      logical :: processed                                       !< for use in sweeps.F90
      logical :: has_previous_timestep                           !< used to prevent timestep retries on freshly created blocks

   contains

      procedure          :: init_gc                              !< Initialization
      procedure          :: cleanup                              !< Deallocate all internals
      procedure, private :: set_coords                           !< Calculate arrays of coordinates along a given direction
      procedure, private :: add_all_na                           !< Register all known named arrays for this cg, sey up shortcuts to the crucial fields
      procedure          :: add_na                               !< Register a new 3D entry in current cg with given name.
      procedure          :: add_na_4d                            !< Register a new 4D entry in current cg with given name.
      procedure          :: update_leafmap                       !< Check if the grid container has any parts covered by finer grids and update appropriate map
      procedure          :: set_fluxpointers
      procedure          :: save_outfluxes
      procedure          :: prolong                              !< perform prolongation of the data stored in this%prolong_
      procedure          :: refinemap2SFC_list                   !< create list of SFC indices to be created from refine flags

   end type grid_container

   interface is_overlap
      module procedure is_overlap_simple, is_overlap_per
   end interface

contains

!>
!! \brief Initialization of the grid container
!!
!! \details This method sets up the grid container variables, coordinates and allocates basic arrays.
!! Everything related to the interior of grid container should be set here.
!! Things that are related to communication with other grid containers or global properties are set up in cg_level::init_all_new_cg.
!<

   subroutine init_gc(this, n_d, off, my_se, grid_id, level_id)

      use constants,    only: PIERNIK_INIT_DOMAIN, xdim, ydim, zdim, ndims, big_float, LO, HI, I_ONE, I_TWO, BND_MPI, BND_COR, GEO_XYZ, GEO_RPZ, dpi
      use dataio_pub,   only: die, warn, code_progress
      use domain,       only: dom
      use grid_helpers, only: f2c
      use ordering,     only: SFC_order

      implicit none

      class(grid_container), target,   intent(inout) :: this ! intent(out) would silently clear everything, that was already set
                                                             ! (also the fields in types derived from grid_container)
      integer(kind=8), dimension(:),   intent(in) :: n_d     !< max resolution of my level
      integer(kind=8), dimension(:),   intent(in) :: off     !< offset of my level
      integer(kind=8), dimension(:,:), intent(in) :: my_se   !< my segment
      integer,                         intent(in) :: grid_id
      integer(kind=4),                 intent(in) :: level_id

      integer :: i
      integer(kind=8), dimension(ndims, LO:HI) :: rn

      if (code_progress < PIERNIK_INIT_DOMAIN) call die("[grid_container:init_gc] MPI not initialized.")

      this%level_off  = off
      this%membership = 1
      this%grid_id    = grid_id
      this%my_se(:,:) = my_se(:, :)
      this%h_cor1(:)  = this%my_se(:, HI) + I_ONE
      this%n_b(:)     = int(this%my_se(:, HI) - this%my_se(:, LO) + I_ONE, 4) ! Block 'physical' grid sizes
      this%level_id   = level_id
      this%SFC_id     = SFC_order(this%my_se(:, LO) - off)

      if (any(this%n_b(:) <= 0)) call die("[grid_container:init_gc] Mixed positive and non-positive grid sizes")

      ! Inherit the boundaries from the domain, then set MPI or SHEAR boundaries where applicable
      this%bnd(:,:) = dom%bnd(:,:)
      where (my_se(:, LO)   /= off(:)         ) this%bnd(:, LO) = BND_MPI
      where (this%h_cor1(:) /= off(:) + n_d(:)) this%bnd(:, HI) = BND_MPI
      ! For periodic boundaries do not set BND_MPI when local domain spans through the whole computational domain in given direction.
      where (dom%periodic(:) .and. this%h_cor1(:)    /= n_d(:)) this%bnd(:, LO) = BND_MPI
      where (dom%periodic(:) .and. this%my_se(:, LO) /= 0)      this%bnd(:, HI) = BND_MPI

      this%ext_bnd(:, :) = .false.
      do i = xdim, zdim
         if (dom%has_dir(i) .and. .not. dom%periodic(i)) then
            this%ext_bnd(i, LO) = (my_se(i, LO)   == off(i))
            this%ext_bnd(i, HI) = (this%h_cor1(i) == off(i) + n_d(i)) !! \warning not true on AMR
         endif
      enddo

      ! For shear boundaries and some domain decompositions it is possible that a boundary can be mixed 'per' with 'mpi'

!      call inflate_req
      ! write_plot_hdf5 requires nproc entries for the status array

      if (any(dom%bnd(xdim:ydim, :) == BND_COR)) call die("[grid_container:init_gc] BND_COR unimplemented")
      if (any(dom%bnd(zdim, :) == BND_COR)) call die("[grid_container:init_gc] Corner BC not allowed for z-direction")

#ifdef SHEAR_BND
      call die("[grid_container:init_gc] Shear-pediodic boundary conditions unimplemented")
      ! This is possible to be implemented
#endif /* SHEAR_BND */

      do i = xdim, zdim
         if (dom%has_dir(i)) then
            if (this%n_b(i) < 1) call die("[grid_init_gc] Too many CPUs for such a small grid.")
            if (this%n_b(i) < dom%nb) call warn("[grid_init_gc] domain size in some directions is < nb, which may result in incomplete boundary cell update")
         endif
      enddo

      where (dom%has_dir(:))
         this%n_(:)        = this%n_b(:) + I_TWO * dom%nb       ! Block total grid size with guardcells
         this%ijkse(:, LO) = int(this%my_se(:, LO), kind=4)
         this%ijkse(:, HI) = int(this%my_se(:, HI), kind=4)
         this%ijkseb(:,LO) = this%ijkse(:, LO) + dom%nb - I_ONE
         this%ijkseb(:,HI) = this%ijkse(:, HI) - dom%nb + I_ONE
         this%lh1(:,LO)    = this%ijkse(:, LO) - I_ONE
         this%lh1(:,HI)    = this%ijkse(:, HI) + I_ONE
         this%lhn(:,LO)    = this%ijkse(:, LO) - dom%nb
         this%lhn(:,HI)    = this%ijkse(:, HI) + dom%nb
         this%dl(:)        = dom%L_(:) / n_d(:)
         this%fbnd(:, LO)  = dom%edge(:, LO) + this%dl(:) * (this%my_se(:, LO) - off(:))
         this%fbnd(:, HI)  = dom%edge(:, LO) + this%dl(:) * (this%h_cor1(:) - off(:))
      elsewhere
         this%n_(:)        = 1
         this%ijkse(:, LO) = 0
         this%ijkse(:, HI) = 0
         this%ijkseb(:,LO) = this%ijkse(:, LO)
         this%ijkseb(:,HI) = this%ijkse(:, HI)
         this%lh1(:,LO)    = this%ijkse(:, LO)
         this%lh1(:,HI)    = this%ijkse(:, HI)
         this%lhn(:,LO)    = this%ijkse(:, LO)
         this%lhn(:,HI)    = this%ijkse(:, HI)
         this%dl(:)        = 1.0
         this%fbnd(:, LO)  = dom%edge(:, LO)
         this%fbnd(:, HI)  = dom%edge(:, HI)
      endwhere

      this%isb = this%ijkseb(xdim, LO)
      this%ieb = this%ijkseb(xdim, HI)
      this%jsb = this%ijkseb(ydim, LO)
      this%jeb = this%ijkseb(ydim, HI)
      this%ksb = this%ijkseb(zdim, LO)
      this%keb = this%ijkseb(zdim, HI)

      select case (dom%geometry_type)
         case (GEO_XYZ)
            this%vol = product(this%fbnd(:, HI)-this%fbnd(:, LO), mask=dom%has_dir(:))
            this%dvol = product(this%dl(:), mask=dom%has_dir(:))
         case (GEO_RPZ)
            if (.not. dom%has_dir(ydim)) then
               this%dl(ydim) = dpi
               this%fbnd(ydim, :) = [ 0., dpi ]
            endif
            this%vol = 1.
            if (dom%has_dir(xdim)) this%vol = this%vol * (this%fbnd(xdim, HI)**2 - this%fbnd(xdim, LO)**2)/2.
            this%vol = this%vol * (this%fbnd(ydim, HI) - this%fbnd(ydim, LO))
            if (dom%has_dir(zdim)) this%vol = this%vol * (this%fbnd(zdim, HI) - this%fbnd(zdim, LO))
            this%dvol = product(this%dl(:), mask=(dom%has_dir(:) .or. [.false., .true., .false.])) ! multiply by actual radius to get true cell volume
      end select

      this%maxxyz = maxval(this%n_(:), mask=dom%has_dir(:))

      call this%set_coords

      this%dxmn = minval(this%dl(:), mask=dom%has_dir(:))

      ! some shortcuts for convenience
      this%idl(:) = 1./this%dl(:)

      this%dx = this%dl(xdim)
      this%dy = this%dl(ydim)
      this%dz = this%dl(zdim)

      this%idx = 1./this%dx
      this%idy = 1./this%dy
      this%idz = 1./this%dz

      !> \deprecated this%n[xyz]b are almost unused. \todo remove it
      this%nxb = this%n_b(xdim)
      this%nyb = this%n_b(ydim)
      this%nzb = this%n_b(zdim)

      this%is = this%ijkse(xdim, LO)
      this%js = this%ijkse(ydim, LO)
      this%ks = this%ijkse(zdim, LO)
      this%ie = this%ijkse(xdim, HI)
      this%je = this%ijkse(ydim, HI)
      this%ke = this%ijkse(zdim, HI)

      ! copied from multigrid
      this%dxy = 1.
      this%dxz = 1.
      this%dyz = 1.

      if (dom%has_dir(xdim)) then
         this%idx2  = 1. / this%dx**2                      ! auxiliary invariants
         this%dxy   = this%dxy * this%dx
         this%dxz   = this%dxz * this%dx
      else
         this%idx2  = 0.
      endif

      if (dom%has_dir(ydim)) then
         this%idy2  = 1. / this%dy**2
         this%dxy   = this%dxy * this%dy
         this%dyz   = this%dyz * this%dy
      else
         this%idy2  = 0.
      endif

      if (dom%has_dir(zdim)) then
         this%idz2  = 1. / this%dz**2
         this%dxz   = this%dxz * this%dz
         this%dyz   = this%dyz * this%dz
      else
         this%idz2  = 0.
      endif

      this%idl2 = [ this%idx2, this%idy2, this%idz2 ]
      nullify(this%prolong_xyz)
      if (allocated(this%prolong_) .or. allocated(this%prolong_x) .or. allocated(this%prolong_xy)) &
           call die("[grid_container:init_gc] prolong_* arrays already allocated")
      ! size of coarsened grid with guardcells, additional cell is required only when even-sized grid has odd offset

      rn = int(this%ijkse, kind=8)
      rn = f2c(rn)
      where (dom%has_dir(:))
         rn(:, LO) = rn(:, LO) - dom%nb
         rn(:, HI) = rn(:, HI) + dom%nb
         ! +1 is because of some simplifications in cg_level::prolong_q_1var in treating grids with odd offsets
      endwhere
      allocate(this%prolong_   (      rn(xdim, LO):      rn(xdim, HI),       rn(ydim, LO):      rn(ydim, HI),       rn(zdim, LO):      rn(zdim, HI)), &
           &   this%prolong_x  (this%lhn(xdim, LO):this%lhn(xdim, HI),       rn(ydim, LO):      rn(ydim, HI),       rn(zdim, LO):      rn(zdim, HI)), &
           &   this%prolong_xy (this%lhn(xdim, LO):this%lhn(xdim, HI), this%lhn(ydim, LO):this%lhn(ydim, HI),       rn(zdim, LO):      rn(zdim, HI)), &
           &   this%prolong_xyz(this%lhn(xdim, LO):this%lhn(xdim, HI), this%lhn(ydim, LO):this%lhn(ydim, HI), this%lhn(zdim, LO):this%lhn(zdim, HI)))
      allocate(this%leafmap  (this%ijkse(xdim, LO):this%ijkse(xdim, HI), this%ijkse(ydim, LO):this%ijkse(ydim, HI), this%ijkse(zdim, LO):this%ijkse(zdim, HI)), &
           &   this%refinemap(this%ijkse(xdim, LO):this%ijkse(xdim, HI), this%ijkse(ydim, LO):this%ijkse(ydim, HI), this%ijkse(zdim, LO):this%ijkse(zdim, HI)))

      this%prolong_   (:, :, :) = big_float
      this%prolong_x  (:, :, :) = big_float
      this%prolong_xy (:, :, :) = big_float
      this%prolong_xyz(:, :, :) = big_float
      this%leafmap    (:, :, :) = .true.
      this%refinemap  (:, :, :) = .false.
      call this%refine_flags%init
      this%ignore_prolongation = .false.
      this%is_old = .false.
      this%has_previous_timestep = .false.

      do i = LO, HI
         call this%finebnd  (xdim, i)%fainit([ this%js, this%je ], [ this%ks, this%ke ])
         call this%finebnd  (ydim, i)%fainit([ this%ks, this%ke ], [ this%is, this%ie ])
         call this%finebnd  (zdim, i)%fainit([ this%is, this%ie ], [ this%js, this%je ])
         call this%coarsebnd(xdim, i)%fainit([ this%js, this%je ], [ this%ks, this%ke ])
         call this%coarsebnd(ydim, i)%fainit([ this%ks, this%ke ], [ this%is, this%ie ])
         call this%coarsebnd(zdim, i)%fainit([ this%is, this%ie ], [ this%js, this%je ])
      enddo
      do i = xdim, zdim
         this%finebnd  (i, LO)%index = this%lhn(i, LO) - 1
         this%finebnd  (i, HI)%index = this%lhn(i, HI) + 1
         this%coarsebnd(i, LO)%index = this%lhn(i, LO) - 1
         this%coarsebnd(i, HI)%index = this%lhn(i, HI) + 1
      enddo

      call this%add_all_na

   end subroutine init_gc

!> \brief Calculate arrays of coordinates along a given direction

   subroutine set_coords(this)

      use constants,  only: LO, HI, half, one, zero, xdim, ydim, zdim, CENTER, LEFT, RIGHT, INV_CENTER
      use dataio_pub, only: die
      use domain,     only: dom
      use func,       only: operator(.notequals.)

      implicit none

      class(grid_container), intent(inout) :: this !< grid container, where the arrays have to be set

      integer :: d, i

      do d = xdim, zdim
         do i = CENTER, INV_CENTER
            if (this%coord(i, d)%associated()) call die("[grid_container:set_coords] a coordinate already allocated")
            call this%coord(i, d)%allocate(this%lhn(d, LO), this%lhn(d, HI))
         enddo

         if (dom%has_dir(d)) then
            this%coord(CENTER, d)%r(:) = this%fbnd(d, LO) + this%dl(d) * ([(i, i=this%lhn(d, LO), this%lhn(d, HI))] - this%ijkse(d, LO) + half)
         else
            this%coord(CENTER, d)%r(:) = half*(this%fbnd(d, LO) + this%fbnd(d, HI))
         endif

         this%coord(LEFT,  d)%r(:) = this%coord(CENTER, d)%r(:) - half*this%dl(d)
         this%coord(RIGHT, d)%r(:) = this%coord(CENTER, d)%r(:) + half*this%dl(d)

         where ( this%coord(CENTER, d)%r(:).notequals.zero )
            this%coord(INV_CENTER, d)%r(:) = one/this%coord(CENTER, d)%r(:)
         elsewhere
            this%coord(INV_CENTER, d)%r(:) = zero
         endwhere

      enddo

      !--- Shortcuts --------------------
      ! left zone boundaries:  xl, yl, zl
      ! zone centers:          x,  y,  z
      ! right zone boundaries: xr, yr, zr

      this%x     => this%coord(CENTER,     xdim)%r
      this%y     => this%coord(CENTER,     ydim)%r
      this%z     => this%coord(CENTER,     zdim)%r

      this%inv_x => this%coord(INV_CENTER, xdim)%r
      this%inv_y => this%coord(INV_CENTER, ydim)%r
      this%inv_z => this%coord(INV_CENTER, zdim)%r

   end subroutine set_coords

!> \brief Routines that deallocates all internals of the grid container

   subroutine cleanup(this)

      use constants, only: xdim, zdim, CENTER, INV_CENTER, LO, HI

      implicit none

      class(grid_container), intent(inout) :: this
      integer :: d, g, b, cdim
      integer, parameter :: nseg = 4*2
      type(tgt_list), dimension(nseg) :: rpio_tgt

      if (associated(this%x))     nullify(this%x)
      if (associated(this%y))     nullify(this%y)
      if (associated(this%z))     nullify(this%z)
      if (associated(this%inv_x)) nullify(this%inv_x)
      if (associated(this%inv_y)) nullify(this%inv_y)
      if (associated(this%inv_z)) nullify(this%inv_z)
      do cdim = xdim, zdim
         do b = CENTER, INV_CENTER
            call this%coord(b, cdim)%deallocate()
         enddo
      enddo

      if (allocated(this%gc_xdim)) deallocate(this%gc_xdim)
      if (allocated(this%gc_ydim)) deallocate(this%gc_ydim)
      if (allocated(this%gc_zdim)) deallocate(this%gc_zdim)

      if (allocated(this%i_bnd)) then
         do d = lbound(this%i_bnd, dim=1), ubound(this%i_bnd, dim=1)
            if (allocated(this%i_bnd(d)%seg)) deallocate(this%i_bnd(d)%seg)
         enddo
         deallocate(this%i_bnd)
      endif
      if (allocated(this%o_bnd)) then
         do d = lbound(this%o_bnd, dim=1), ubound(this%o_bnd, dim=1)
            if (allocated(this%o_bnd(d)%seg)) deallocate(this%o_bnd(d)%seg)
         enddo
         deallocate(this%o_bnd)
      endif

      rpio_tgt(1:nseg) = [ this%ri_tgt,  this%ro_tgt,  this%pi_tgt,  this%po_tgt, &
           &               this%pib_tgt, this%pob_tgt, this%rif_tgt, this%rof_tgt ]
      do b = 1, nseg
         if (allocated(rpio_tgt(b)%seg)) then
            do g = lbound(rpio_tgt(b)%seg, dim=1), ubound(rpio_tgt(b)%seg, dim=1)
               if (allocated(rpio_tgt(b)%seg(g)%buf)) deallocate(rpio_tgt(b)%seg(g)%buf)
            enddo
            deallocate(rpio_tgt(b)%seg)
         endif
      enddo

      if (allocated(this%q)) then
         do g = lbound(this%q(:), dim=1), ubound(this%q(:), dim=1)
            call this%q(g)%clean
         enddo
         deallocate(this%q)
      endif

      if (allocated(this%w)) then
         do g = lbound(this%w(:), dim=1), ubound(this%w(:), dim=1)
            call this%w(g)%clean
         enddo
         deallocate(this%w)
      endif

      ! arrays not handled through named_array feature
      if (associated(this%prolong_xyz)) deallocate(this%prolong_xyz)
      if (allocated(this%prolong_xy))   deallocate(this%prolong_xy)
      if (allocated(this%prolong_x))    deallocate(this%prolong_x)
      if (allocated(this%prolong_))     deallocate(this%prolong_)
      if (allocated(this%leafmap))      deallocate(this%leafmap)
      if (allocated(this%refinemap))    deallocate(this%refinemap)
      do d = xdim, zdim
         do g = LO, HI
            call this%finebnd  (d, g)%facleanup
            call this%coarsebnd(d, g)%facleanup
         enddo
      enddo

      call fpl%fpcleanup
      call fpr%fpcleanup
      call cpl%fpcleanup
      call cpr%fpcleanup

   end subroutine cleanup

!>
!! \brief is_overlap_per checks if two given blocks placed within a periodic domain are overlapping.
!!
!! \details to handle shearing box which is divided in y-direction at the edges, one has to provide another subroutine (is_overlap_per_shear) and add it to interface is_overlap
!<

   logical function is_overlap_per(this, other, periods) result(share)

      use constants,  only: xdim, ydim, zdim, LO, HI
      use domain,     only: dom

      implicit none

      integer(kind=8), dimension(xdim:zdim, LO:HI), intent(in) :: this    !< this box
      integer(kind=8), dimension(xdim:zdim, LO:HI), intent(in) :: other   !< the other box
      integer(kind=8), dimension(xdim:zdim),        intent(in) :: periods !< where >0 then the direction is periodic with the given number of cells

      integer :: i, j, k
      integer(kind=8), dimension(xdim:zdim, LO:HI) :: oth

      share = .false.
      do i = -1, 1
         if ((dom%has_dir(xdim) .or. periods(xdim)>0) .or. i==0) then
            oth(xdim, :) = other(xdim, :) + i*periods(xdim)
            do j = -1, 1
               if ((dom%has_dir(ydim) .or. periods(ydim)>0) .or. j==0) then
                  oth(ydim, :) = other(ydim, :) + j*periods(ydim)
                  do k = -1, 1
                     if ((dom%has_dir(zdim) .or. periods(zdim)>0) .or. k==0) then
                        oth(zdim, :) = other(zdim, :) + k*periods(zdim)
                        share = is_overlap_simple(this, oth) .or. share
                     endif
                  enddo
               endif
            enddo
         endif
      enddo

   end function is_overlap_per

!>
!! \brief is_overlap_simple checks if two given blocks placed within a nonperiodic domain are overlapping.
!! This routine is not supposed to take care of periodic domain - use is_overlap_per when you check overlap for boxes that cross the periodic domain boundary
!<

   logical function is_overlap_simple(this, other) result(share)

      use constants,  only: xdim, zdim, LO, HI
      use domain,     only: dom

      implicit none

      integer(kind=8), dimension(:, :), intent(in) :: this  !< this box
      integer(kind=8), dimension(:, :), intent(in) :: other !< the other box

      integer :: d

      share = .true.
      do d = xdim, zdim
         if (dom%has_dir(d)) share = share .and. (other(d, LO) <= this(d, HI)) .and. (other(d, HI) >= this(d, LO))
      enddo

   end function is_overlap_simple

!> \brief Register all known named arrays for this cg, sey up shortcuts to the crucial fields

   subroutine add_all_na(this)

      use named_array_list, only: qna, wna
#ifdef ISO
      use constants,   only: cs_i2_n
      use fluids_pub,  only: cs2_max
#endif /* ISO */

      implicit none

      class(grid_container), intent(inout) :: this

      integer :: i

      if (allocated(qna%lst)) then
         do i = lbound(qna%lst(:), dim=1), ubound(qna%lst(:), dim=1)
            call this%add_na(qna%lst(i)%multigrid)
         enddo
      endif
      if (allocated(wna%lst)) then
         do i = lbound(wna%lst(:), dim=1), ubound(wna%lst(:), dim=1)
            call this%add_na_4d(wna%lst(i)%dim4)
         enddo
      endif

      ! shortcuts
      this%u  => this%w(wna%fi)%arr
      this%b  => this%w(wna%bi)%arr
      this%wa => this%q(qna%wai)%arr
#ifdef ISO
      this%cs_iso2 => this%q(qna%ind(cs_i2_n))%arr
      if (associated(this%cs_iso2)) this%cs_iso2(:,:,:) = cs2_max   ! set cs2 with sane values on non-multigrid grid pieces
#endif /* ISO */

   end subroutine add_all_na

!>
!! \brief Register a new 3D entry in current cg with given name. Called from cg_list_glob::reg_var
!!
!! \warning This routine should not be called directly from user routines
!<
   subroutine add_na(this, multigrid)

      use constants,   only: base_level_id
      use named_array, only: named_array3d

      implicit none

      class(grid_container), intent(inout) :: this
      logical,               intent(in)    :: multigrid     !< If .true. then cg%q(:)%arr and cg%w(:)%arr are allocated also below base level

      type(named_array3d), allocatable, dimension(:) :: tmp

      if (.not. allocated(this%q)) then
         allocate(this%q(1))
      else
         allocate(tmp(lbound(this%q(:),dim=1):ubound(this%q(:), dim=1) + 1))
         tmp(:ubound(this%q(:), dim=1)) = this%q(:)
         call move_alloc(from=tmp, to=this%q)
      endif

      if (multigrid .or. this%level_id >= base_level_id) call this%q(ubound(this%q(:), dim=1))%init(this%lhn(:, LO), this%lhn(:, HI))

   end subroutine add_na

!>
!! \brief Register a new 4D entry in current cg with given name. Called from cg_list_glob::reg_var
!!
!! \warning This routine should not be called directly from user routines
!! \deprecated Almost duplicated code with add_na
!<
   subroutine add_na_4d(this, n)

      use constants,   only: base_level_id, INT4
      use named_array, only: named_array4d

      implicit none

      class(grid_container), intent(inout) :: this
      integer(kind=4),       intent(in)    :: n            !< Length of the vector quantity to be stored (first dimension of the array)

      type(named_array4d), allocatable, dimension(:) :: tmp

      if (.not. allocated(this%w)) then
         allocate(this%w(1))
      else
         allocate(tmp(lbound(this%w(:),dim=1):ubound(this%w(:), dim=1) + 1))
         tmp(:ubound(this%w(:), dim=1)) = this%w(:)
         call move_alloc(from=tmp, to=this%w)
      endif

      if (this%level_id >= base_level_id) call this%w(ubound(this%w(:), dim=1))%init( [1_INT4, this%lhn(:, LO)], [n, this%lhn(:, HI)] ) !< \deprecated magic integer

   end subroutine add_na_4d

!> \brief Check if the grid container has any parts covered by finer grids and update appropriate map

   subroutine update_leafmap(this)

      use constants, only: xdim, ydim, zdim, LO, HI

      implicit none

      class(grid_container), intent(inout) :: this

      integer(kind=8), dimension(xdim:zdim, LO:HI) :: se
      integer :: g

      this%leafmap = .true.
      if (allocated(this%ri_tgt%seg)) then
         do g = lbound(this%ri_tgt%seg(:), dim=1), ubound(this%ri_tgt%seg(:), dim=1)
            se(:, :) = this%ri_tgt%seg(g)%se(:, :)
            this%leafmap(se(xdim, LO):se(xdim, HI), se(ydim, LO):se(ydim, HI), se(zdim, LO):se(zdim, HI)) = .false.
         enddo
      endif

   end subroutine update_leafmap

!> \brief Add a new segment, reallocate if necessary

   subroutine add_seg(this, proc, se, tag)

      use dataio_pub, only: die

      implicit none

      class(bnd_list),                              intent(inout) :: this
      integer,                                      intent(in)    :: proc
      integer(kind=8), dimension(xdim:zdim, LO:HI), intent(in)    :: se
      integer(kind=4),                              intent(in)    :: tag

      type(segment), dimension(:), allocatable :: tmp
      integer :: g

      if (tag <0) call die("[grid_container:add_seg] tag<0")

      if (allocated(this%seg)) then
         allocate(tmp(lbound(this%seg, dim=1):ubound(this%seg, dim=1)+1))
         tmp(:ubound(this%seg, dim=1)) = this%seg
         call move_alloc(from=tmp, to=this%seg)
      else
         allocate(this%seg(1))
      endif

      g = ubound(this%seg, dim=1)

      this%seg(g)%proc = proc
      this%seg(g)%se = se
      this%seg(g)%tag = tag
      nullify(this%seg(g)%local)

   end subroutine add_seg

   subroutine set_fluxpointers(this, cdim, i1, i2, eflx)

      use constants,  only: LO, HI
      use fluxtypes,  only: ext_fluxes

      implicit none

      class(grid_container), intent(in)    :: this
      integer(kind=4),       intent(in)    :: cdim
      integer,               intent(in)    :: i1, i2
      type(ext_fluxes),      intent(inout) :: eflx

      if (this%finebnd(cdim, LO)%index(i1, i2) >= this%ijkse(cdim, LO)) then
         fpl = this%finebnd(cdim, LO)%fa2fp(i1, i2)
         if (.not. allocated(fpl%uflx)) call fpl%fpinit
         eflx%li => fpl
      else
         nullify(eflx%li)
      endif
      if (this%finebnd(cdim, HI)%index(i1, i2) <= this%ijkse(cdim, HI)) then
         fpr = this%finebnd(cdim, HI)%fa2fp(i1, i2)
         if (.not. allocated(fpr%uflx)) call fpr%fpinit
         eflx%ri => fpr
      else
         nullify(eflx%ri)
      endif

      if (this%coarsebnd(cdim, LO)%index(i1, i2) >= this%ijkse(cdim, LO)) then
         cpl%index = this%coarsebnd(cdim, LO)%index(i1, i2)
         if (.not. allocated(cpl%uflx)) call cpl%fpinit
         eflx%lo => cpl
      else
         nullify(eflx%lo)
      endif
      if (this%coarsebnd(cdim, HI)%index(i1, i2) <= this%ijkse(cdim, HI)) then
         cpr%index = this%coarsebnd(cdim, HI)%index(i1, i2)
         if (.not. allocated(cpr%uflx)) call cpr%fpinit
         eflx%ro => cpr
      else
         nullify(eflx%ro)
      endif

   end subroutine set_fluxpointers

   subroutine save_outfluxes(this, cdim, i1, i2, eflx)

      use constants,  only: LO, HI
      use fluxtypes,  only: ext_fluxes

      implicit none

      class(grid_container), intent(inout) :: this
      integer(kind=4),       intent(in)    :: cdim
      integer,               intent(in)    :: i1, i2
      type(ext_fluxes),      intent(inout) :: eflx

      if (associated(eflx%lo)) call this%coarsebnd(cdim, LO)%fp2fa(eflx%lo, i1, i2)
      if (associated(eflx%ro)) call this%coarsebnd(cdim, HI)%fp2fa(eflx%ro, i1, i2)

   end subroutine save_outfluxes

!>
!! \brief perform high-order order prolongation interpolation of the data stored in this%prolong_
!!
!! \details
!! <table border="1" cellpadding="4" cellspacing="0">
!!   <tr><td> Cell-face prolongation stencils for fast convergence on uniform grid </td>
!!       <td> -1./12. </td><td> 7./12. </td><td> 7./12. </td><td> -1./12. </td><td> integral cubic </td></tr>
!!   <tr><td> Slightly slower convergence, less wide stencil  </td>
!!       <td>         </td><td> 1./2.  </td><td> 1./2.  </td><td>         </td><td> average; integral and direct linear </td></tr>
!! </table>
!!\n Prolongation of cell faces from cell centers are required for FFT local solver, red-black Gauss-Seidel relaxation don't use it.
!!
!!\n Cell-centered prolongation stencils, for odd fine cells, for even fine cells reverse the order.
!! <table border="1" cellpadding="4" cellspacing="0">
!!   <tr><td> 35./2048. </td><td> -252./2048. </td><td> 1890./2048. </td><td> 420./2048. </td><td> -45./2048. </td><td> direct quartic </td></tr>
!!   <tr><td>           </td><td>   -7./128.  </td><td>  105./128.  </td><td>  35./128.  </td><td>  -5./128.  </td><td> direct cubic </td></tr>
!!   <tr><td>           </td><td>   -3./32.   </td><td>   30./32.   </td><td>   5./32.   </td><td>            </td><td> direct quadratic </td></tr>
!!   <tr><td>           </td><td>             </td><td>    1.       </td><td>            </td><td>            </td><td> injection (0-th order), direct and integral approach </td></tr>
!!   <tr><td>           </td><td>             </td><td>    3./4.    </td><td>   1./4.    </td><td>            </td><td> linear, direct and integral approach </td></tr>
!!   <tr><td>           </td><td>   -1./8.    </td><td>    1.       </td><td>   1./8.    </td><td>            </td><td> integral quadratic </td></tr>
!!   <tr><td>           </td><td>   -5./64.   </td><td>   55./64.   </td><td>  17./64.   </td><td>  -3./64.   </td><td> integral cubic </td></tr>
!!   <tr><td>   3./128. </td><td>  -11./64.   </td><td>    1.       </td><td>  11./64.   </td><td>  -3./128.  </td><td> integral quartic </td></tr>
!! </table>
!!
!!\n General rule is that the second big positive coefficient should be assigned to closer neighbor of the coarse parent cell.
!!\n Thus a single coarse contributes to fine cells in the following way:
!! <table border="1" cellpadding="4" cellspacing="0">
!!   <tr><td> fine level   </td>
!!       <td> -3./128. </td><td> 3./128. </td><td> -11./64. </td><td>  11./64. </td><td> 1. </td><td> 1. </td>
!!       <td> 11./64. </td><td> -11./64. </td><td> 3./128.  </td><td> -3./128. </td><td> integral quartic coefficients </td></tr>
!!   <tr><td> coarse level </td>
!!       <td colspan="2">                </td><td colspan="2">                 </td><td colspan="2"> 1.  </td>
!!       <td colspan="2">                </td><td colspan="2">                 </td><td>                               </td></tr>
!! </table>
!!\n
!!\n The term "n-th order integral interpolation" here means that the prolonged values satisfy the following condition:
!!\n Integral over a cell of a n-th order polynomial fit to the nearest 5 points in each dimension on coarse level
!! is equal to the sum of similar integrals over fine cells covering the coarse cell.
!!\n
!!\n The term "n-th order direct interpolation" here means that the prolonged values are n-th order polynomial fit
!! to the nearest 5 points in each dimension on coarse level evaluated for fine cell centers.
!!\n
!!\n It seems that for 3D Cartesian grid with isolated boundaries and relaxation implemented in approximate_solution
!! direct quadratic and cubic interpolations give best norm reduction factors per V-cycle (maclaurin problem).
!!  For other boundary types, FFT implementation of approximate_solution, specific source distribution or
!!  some other multigrid scheme may give faster convergence rate.
!!\n
!!\n Estimated prolongation costs for integral quartic stencil:
!!\n  "gather" approach: loop over fine cells, each one collects weighted values from 5**3 coarse cells (125*n_fine multiplications
!!\n  "scatter" approach: loop over coarse cells, each one contributes weighted values to 10**3 fine cells (1000*n_coarse multiplications, roughly equal to cost of gather)
!!\n  "directionally split" approach: do the prolongation (either gather or scatter type) first in x direction (10*n_coarse multiplications -> 2*n_coarse intermediate cells
!!                                  result), then in y direction (10*2*n_coarse multiplications -> 4*n_coarse intermediate cells result), then in z direction
!!                                  (10*4*n_coarse multiplications -> 8*n_coarse = n_fine cells result). Looks like 70*n_coarse multiplications.
!!                                  Will require two additional arrays for intermediate results.
!!\n  "FFT" approach: do the convolution in Fourier space. Unfortunately it is not periodic box, so we cannot use power of 2 FFT sizes. No idea how fast or slow this can be.
!!\n
!!\n For AMR or nested grid low-order prolongation schemes (injection and linear interpolation at least) are known to produce artifacts
!! on fine-coarse boundaries. For uniform grid the simplest operators are probably the fastest and give best V-cycle convergence rates.
!<

   subroutine prolong(this, ind, cse, p_xyz)

      use constants,          only: xdim, ydim, zdim, zero, LO, HI, I_ZERO, I_ONE, I_TWO, I_THREE, O_INJ, O_LIN, O_D2, O_D3, O_D4, O_D5, O_D6, O_I2, O_I3, O_I4
      use dataio_pub,         only: die
      use domain,             only: dom
      use func,               only: operator(.notequals.)
      use grid_helpers,       only: c2f
      use named_array_list,   only: qna

      implicit none

      class(grid_container),                        intent(inout) :: this
      integer(kind=4),                              intent(in)    :: ind   !< index of cg%q(:) 3d array - variable to be prolonged
      integer(kind=8), dimension(xdim:zdim, LO:HI), intent(in)    :: cse   !< coarse segment
      logical,                                      intent(in)    :: p_xyz !< store the result in this%prolong_xyz when true, in this%q(ind)%arr otherwise

      integer :: stencil_range        !< how far to look for the data to be prolonged
      integer(kind=8), dimension(xdim:zdim) :: D
      integer(kind=8), dimension(xdim:zdim, LO:HI) :: fse ! fine segment
      real :: P_3, P_2, P_1, P0, P1, P2, P3 !< interpolation coefficients
      real, dimension(:,:,:), pointer :: pa3d

      if (p_xyz) then
         pa3d => this%prolong_xyz
      else
         pa3d => this%q(ind)%arr
      endif

      select case (qna%lst(ind)%ord_prolong)
         case (O_D6)
            P_3 = -231./65536. ; P_2 = 2002./65536.; P_1 = -9009./65536.; P0 = 60060./65536.; P1 = 15015./65536.; P2 = -2574./65536.; P3 = 273./65536.
         case (O_D5)
            P_3 = 0.;            P_2 = 77./8192.;    P_1 = -693./8192.;   P0 = 6930./8192.;   P1 = 2310./8192.;   P2 = -495./8192.;   P3 = 63./8192.
         case (O_D4)
            P_3 = 0.;            P_2 = 35./2048.;    P_1 = -252./2048.;   P0 = 1890./2048.;   P1 = 420./2048.;    P2 = -45./2048.;    P3 = 0.
         case (O_D3)
            P_3 = 0.;            P_2 = 0.;           P_1 = -7./128.;      P0 = 105./128.;     P1 = 35./128.;      P2 = -5./128.;      P3 = 0.
         case (O_D2)
            P_3 = 0.;            P_2 = 0.;           P_1 = -3./32.;       P0 = 30./32.;       P1 = 5./32.;        P2 = 0.;            P3 = 0.
         case (O_LIN)
            P_3 = 0.;            P_2 = 0.;           P_1 = 0.;            P0 = 3./4.;         P1 = 1./4.;         P2 = 0.;            P3 = 0.
         case (O_INJ)
            P_3 = 0.;            P_2 = 0.;           P_1 = 0.;            P0 = 1.;            P1 = 0.;            P2 = 0.;            P3 = 0.
         case (O_I2)
            P_3 = 0.;            P_2 = 0.;           P_1 = -1./8.;        P0 = 1.;            P1 = 1./8.;         P2 = 0.;            P3 = 0.
         case (O_I3)
            P_3 = 0.;            P_2 = 0.;           P_1 = -5./64.;       P0 = 55./64;        P1 = 17./64.;       P2 = -3./64.;       P3 = 0.
         case (O_I4)
            P_3 = 0.;            P_2 = 3./128.;      P_1 = -11./64.;      P0 = 1.;            P1 = 11./64.;       P2 = -3./128.;      P3 = 0.
         case default
            call die("[grid_container:prolong] Unsupported order")
            return
      end select

      ! this is just for optimization. Setting stencil_range = I_THREE should work correctly for all interpolations.
      stencil_range = I_ZERO
      if ((P_1.notequals.zero).or.(P1.notequals.zero)) stencil_range = I_ONE
      if ((P_2.notequals.zero).or.(P2.notequals.zero)) stencil_range = I_TWO
      if ((P_3.notequals.zero).or.(P3.notequals.zero)) stencil_range = I_THREE

      where (dom%has_dir(:))
         D(:) = 1
      elsewhere
         D(:) = 0
      endwhere

      !> \deprecated the comments below are quite old and may be outdated or inaccurate.
      ! When the grid offset is odd, the coarse data is shifted by half coarse cell (or one fine cell)
      ! odd(:) = int(mod(cg%off(:), int(refinement_factor, kind=8)), kind=4)
      ! When the grid offset is odd we need to apply mirrored prolongation stencil (swap even and odd stencils)
      ! when dom%nb is odd, one, most distant, layer of cells is not filled up

      fse = c2f(cse)

      ! Perform directional-split interpolation
      select case (stencil_range*dom%D_x) ! stencil_range or I_ZERO if .not. dom%has_dir(xdim)
         case (I_ZERO)
            this%prolong_x      (fse(xdim, LO):fse(xdim, HI):2, :, :) = &
                 this%prolong_  (cse(xdim, LO):cse(xdim, HI),   :, :)
            if (dom%has_dir(xdim)) &
                 this%prolong_x (fse(xdim, LO)+dom%D_x:fse(xdim, HI)+dom%D_x:2, :, :) = &
                 & this%prolong_(cse(xdim, LO):cse(xdim, HI),                   :, :)
         case (I_ONE)
            this%prolong_x          (fse(xdim, LO)        :fse(xdim, HI):2,         cse(ydim, LO)-dom%D_y:cse(ydim, HI)+dom%D_y, cse(zdim, LO)-dom%D_z:cse(zdim, HI)+dom%D_z) = &
                 +P1 * this%prolong_(cse(xdim, LO)-D(xdim):cse(xdim, HI)-D(xdim),   cse(ydim, LO)-dom%D_y:cse(ydim, HI)+dom%D_y, cse(zdim, LO)-dom%D_z:cse(zdim, HI)+dom%D_z) &
                 +P0 * this%prolong_(cse(xdim, LO)        :cse(xdim, HI),           cse(ydim, LO)-dom%D_y:cse(ydim, HI)+dom%D_y, cse(zdim, LO)-dom%D_z:cse(zdim, HI)+dom%D_z) &
                 +P_1* this%prolong_(cse(xdim, LO)+D(xdim):cse(xdim, HI)+D(xdim),   cse(ydim, LO)-dom%D_y:cse(ydim, HI)+dom%D_y, cse(zdim, LO)-dom%D_z:cse(zdim, HI)+dom%D_z)
            this%prolong_x          (fse(xdim, LO)+dom%D_x:fse(xdim, HI)+dom%D_x:2, cse(ydim, LO)-dom%D_y:cse(ydim, HI)+dom%D_y, cse(zdim, LO)-dom%D_z:cse(zdim, HI)+dom%D_z) = &
                 +P_1* this%prolong_(cse(xdim, LO)-D(xdim):cse(xdim, HI)-D(xdim),   cse(ydim, LO)-dom%D_y:cse(ydim, HI)+dom%D_y, cse(zdim, LO)-dom%D_z:cse(zdim, HI)+dom%D_z) &
                 +P0 * this%prolong_(cse(xdim, LO)        :cse(xdim, HI),           cse(ydim, LO)-dom%D_y:cse(ydim, HI)+dom%D_y, cse(zdim, LO)-dom%D_z:cse(zdim, HI)+dom%D_z) &
                 +P1 * this%prolong_(cse(xdim, LO)+D(xdim):cse(xdim, HI)+D(xdim),   cse(ydim, LO)-dom%D_y:cse(ydim, HI)+dom%D_y, cse(zdim, LO)-dom%D_z:cse(zdim, HI)+dom%D_z)
         case (I_TWO)
            this%prolong_x          (fse(xdim, LO)          :fse(xdim, HI):2,         cse(ydim, LO)-2*dom%D_y:cse(ydim, HI)+2*dom%D_y, cse(zdim, LO)-2*dom%D_z:cse(zdim, HI)+2*dom%D_z) = &
                 +P2 * this%prolong_(cse(xdim, LO)-2*D(xdim):cse(xdim, HI)-2*D(xdim), cse(ydim, LO)-2*dom%D_y:cse(ydim, HI)+2*dom%D_y, cse(zdim, LO)-2*dom%D_z:cse(zdim, HI)+2*dom%D_z) &
                 +P1 * this%prolong_(cse(xdim, LO)-  D(xdim):cse(xdim, HI)-  D(xdim), cse(ydim, LO)-2*dom%D_y:cse(ydim, HI)+2*dom%D_y, cse(zdim, LO)-2*dom%D_z:cse(zdim, HI)+2*dom%D_z) &
                 +P0 * this%prolong_(cse(xdim, LO)          :cse(xdim, HI),           cse(ydim, LO)-2*dom%D_y:cse(ydim, HI)+2*dom%D_y, cse(zdim, LO)-2*dom%D_z:cse(zdim, HI)+2*dom%D_z) &
                 +P_1* this%prolong_(cse(xdim, LO)+  D(xdim):cse(xdim, HI)+  D(xdim), cse(ydim, LO)-2*dom%D_y:cse(ydim, HI)+2*dom%D_y, cse(zdim, LO)-2*dom%D_z:cse(zdim, HI)+2*dom%D_z) &
                 +P_2* this%prolong_(cse(xdim, LO)+2*D(xdim):cse(xdim, HI)+2*D(xdim), cse(ydim, LO)-2*dom%D_y:cse(ydim, HI)+2*dom%D_y, cse(zdim, LO)-2*dom%D_z:cse(zdim, HI)+2*dom%D_z)
            this%prolong_x          (fse(xdim, LO)+dom%D_x  :fse(xdim, HI)+dom%D_x:2, cse(ydim, LO)-2*dom%D_y:cse(ydim, HI)+2*dom%D_y, cse(zdim, LO)-2*dom%D_z:cse(zdim, HI)+2*dom%D_z) = &
                 +P_2* this%prolong_(cse(xdim, LO)-2*D(xdim):cse(xdim, HI)-2*D(xdim), cse(ydim, LO)-2*dom%D_y:cse(ydim, HI)+2*dom%D_y, cse(zdim, LO)-2*dom%D_z:cse(zdim, HI)+2*dom%D_z) &
                 +P_1* this%prolong_(cse(xdim, LO)-  D(xdim):cse(xdim, HI)-  D(xdim), cse(ydim, LO)-2*dom%D_y:cse(ydim, HI)+2*dom%D_y, cse(zdim, LO)-2*dom%D_z:cse(zdim, HI)+2*dom%D_z) &
                 +P0 * this%prolong_(cse(xdim, LO)          :cse(xdim, HI),           cse(ydim, LO)-2*dom%D_y:cse(ydim, HI)+2*dom%D_y, cse(zdim, LO)-2*dom%D_z:cse(zdim, HI)+2*dom%D_z) &
                 +P1 * this%prolong_(cse(xdim, LO)+  D(xdim):cse(xdim, HI)+  D(xdim), cse(ydim, LO)-2*dom%D_y:cse(ydim, HI)+2*dom%D_y, cse(zdim, LO)-2*dom%D_z:cse(zdim, HI)+2*dom%D_z) &
                 +P2 * this%prolong_(cse(xdim, LO)+2*D(xdim):cse(xdim, HI)+2*D(xdim), cse(ydim, LO)-2*dom%D_y:cse(ydim, HI)+2*dom%D_y, cse(zdim, LO)-2*dom%D_z:cse(zdim, HI)+2*dom%D_z)
         case (I_THREE)
            this%prolong_x          (fse(xdim, LO)          :fse(xdim, HI):2,         cse(ydim, LO)-3*dom%D_y:cse(ydim, HI)+3*dom%D_y, cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) = &
                 +P3 * this%prolong_(cse(xdim, LO)-3*D(xdim):cse(xdim, HI)-3*D(xdim), cse(ydim, LO)-3*dom%D_y:cse(ydim, HI)+3*dom%D_y, cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) &
                 +P2 * this%prolong_(cse(xdim, LO)-2*D(xdim):cse(xdim, HI)-2*D(xdim), cse(ydim, LO)-3*dom%D_y:cse(ydim, HI)+3*dom%D_y, cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) &
                 +P1 * this%prolong_(cse(xdim, LO)-  D(xdim):cse(xdim, HI)-  D(xdim), cse(ydim, LO)-3*dom%D_y:cse(ydim, HI)+3*dom%D_y, cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) &
                 +P0 * this%prolong_(cse(xdim, LO)          :cse(xdim, HI),           cse(ydim, LO)-3*dom%D_y:cse(ydim, HI)+3*dom%D_y, cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) &
                 +P_1* this%prolong_(cse(xdim, LO)+  D(xdim):cse(xdim, HI)+  D(xdim), cse(ydim, LO)-3*dom%D_y:cse(ydim, HI)+3*dom%D_y, cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) &
                 +P_2* this%prolong_(cse(xdim, LO)+2*D(xdim):cse(xdim, HI)+2*D(xdim), cse(ydim, LO)-3*dom%D_y:cse(ydim, HI)+3*dom%D_y, cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) &
                 +P_3* this%prolong_(cse(xdim, LO)+3*D(xdim):cse(xdim, HI)+3*D(xdim), cse(ydim, LO)-3*dom%D_y:cse(ydim, HI)+3*dom%D_y, cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z)
            this%prolong_x          (fse(xdim, LO)+dom%D_x  :fse(xdim, HI)+dom%D_x:2, cse(ydim, LO)-3*dom%D_y:cse(ydim, HI)+3*dom%D_y, cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) = &
                 +P_3* this%prolong_(cse(xdim, LO)-3*D(xdim):cse(xdim, HI)-3*D(xdim), cse(ydim, LO)-3*dom%D_y:cse(ydim, HI)+3*dom%D_y, cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) &
                 +P_2* this%prolong_(cse(xdim, LO)-2*D(xdim):cse(xdim, HI)-2*D(xdim), cse(ydim, LO)-3*dom%D_y:cse(ydim, HI)+3*dom%D_y, cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) &
                 +P_1* this%prolong_(cse(xdim, LO)-  D(xdim):cse(xdim, HI)-  D(xdim), cse(ydim, LO)-3*dom%D_y:cse(ydim, HI)+3*dom%D_y, cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) &
                 +P0 * this%prolong_(cse(xdim, LO)          :cse(xdim, HI),           cse(ydim, LO)-3*dom%D_y:cse(ydim, HI)+3*dom%D_y, cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) &
                 +P1 * this%prolong_(cse(xdim, LO)+  D(xdim):cse(xdim, HI)+  D(xdim), cse(ydim, LO)-3*dom%D_y:cse(ydim, HI)+3*dom%D_y, cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) &
                 +P2 * this%prolong_(cse(xdim, LO)+2*D(xdim):cse(xdim, HI)+2*D(xdim), cse(ydim, LO)-3*dom%D_y:cse(ydim, HI)+3*dom%D_y, cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) &
                 +P3 * this%prolong_(cse(xdim, LO)+3*D(xdim):cse(xdim, HI)+3*D(xdim), cse(ydim, LO)-3*dom%D_y:cse(ydim, HI)+3*dom%D_y, cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z)
         case default
            call die("[grid_container:prolong] unsupported stencil size")
      end select

      select case (stencil_range*dom%D_y)
         case (I_ZERO)
            this%prolong_xy      (fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI):2, :) = &
                 this%prolong_x  (fse(xdim, LO):fse(xdim, HI), cse(ydim, LO):cse(ydim, HI),   :)
            if (dom%has_dir(ydim)) &
                 this%prolong_xy (fse(xdim, LO):fse(xdim, HI), fse(ydim, LO)+dom%D_y:fse(ydim, HI)+dom%D_y:2, :) = &
                 & this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO):cse(ydim, HI),                   :)
         case (I_ONE)
            this%prolong_xy           (fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI):2,               cse(zdim, LO)-dom%D_z:cse(zdim, HI)+dom%D_z) = &
                 + P1 * this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)-D(ydim):cse(ydim, HI)-D(ydim), cse(zdim, LO)-dom%D_z:cse(zdim, HI)+dom%D_z) &
                 + P0 * this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)        :cse(ydim, HI),         cse(zdim, LO)-dom%D_z:cse(zdim, HI)+dom%D_z) &
                 + P_1* this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)+D(ydim):cse(ydim, HI)+D(ydim), cse(zdim, LO)-dom%D_z:cse(zdim, HI)+dom%D_z)
            this%prolong_xy           (fse(xdim, LO):fse(xdim, HI), fse(ydim, LO)+dom%D_y:fse(ydim, HI)+dom%D_y:2, cse(zdim, LO)-dom%D_z:cse(zdim, HI)+dom%D_z) = &
                 + P_1* this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)-D(ydim):cse(ydim, HI)-D(ydim),   cse(zdim, LO)-dom%D_z:cse(zdim, HI)+dom%D_z) &
                 + P0 * this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)        :cse(ydim, HI),           cse(zdim, LO)-dom%D_z:cse(zdim, HI)+dom%D_z) &
                 + P1 * this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)+D(ydim):cse(ydim, HI)+D(ydim),   cse(zdim, LO)-dom%D_z:cse(zdim, HI)+dom%D_z)
         case (I_TWO)
            this%prolong_xy           (fse(xdim, LO):fse(xdim, HI), fse(ydim, LO)          :fse(ydim, HI):2,         cse(zdim, LO)-2*dom%D_z:cse(zdim, HI)+2*dom%D_z) = &
                 + P2 * this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)-2*D(ydim):cse(ydim, HI)-2*D(ydim), cse(zdim, LO)-2*dom%D_z:cse(zdim, HI)+2*dom%D_z) &
                 + P1 * this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)-  D(ydim):cse(ydim, HI)-  D(ydim), cse(zdim, LO)-2*dom%D_z:cse(zdim, HI)+2*dom%D_z) &
                 + P0 * this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)          :cse(ydim, HI),           cse(zdim, LO)-2*dom%D_z:cse(zdim, HI)+2*dom%D_z) &
                 + P_1* this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)+  D(ydim):cse(ydim, HI)+  D(ydim), cse(zdim, LO)-2*dom%D_z:cse(zdim, HI)+2*dom%D_z) &
                 + P_2* this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)+2*D(ydim):cse(ydim, HI)+2*D(ydim), cse(zdim, LO)-2*dom%D_z:cse(zdim, HI)+2*dom%D_z)
            this%prolong_xy           (fse(xdim, LO):fse(xdim, HI), fse(ydim, LO)+dom%D_y  :fse(ydim, HI)+dom%D_y:2, cse(zdim, LO)-2*dom%D_z:cse(zdim, HI)+2*dom%D_z) = &
                 + P_2* this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)-2*D(ydim):cse(ydim, HI)-2*D(ydim), cse(zdim, LO)-2*dom%D_z:cse(zdim, HI)+2*dom%D_z) &
                 + P_1* this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)-  D(ydim):cse(ydim, HI)-  D(ydim), cse(zdim, LO)-2*dom%D_z:cse(zdim, HI)+2*dom%D_z) &
                 + P0 * this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)          :cse(ydim, HI),           cse(zdim, LO)-2*dom%D_z:cse(zdim, HI)+2*dom%D_z) &
                 + P1 * this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)+  D(ydim):cse(ydim, HI)+  D(ydim), cse(zdim, LO)-2*dom%D_z:cse(zdim, HI)+2*dom%D_z) &
                 + P2 * this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)+2*D(ydim):cse(ydim, HI)+2*D(ydim), cse(zdim, LO)-2*dom%D_z:cse(zdim, HI)+2*dom%D_z)
         case (I_THREE)
            this%prolong_xy           (fse(xdim, LO):fse(xdim, HI), fse(ydim, LO)          :fse(ydim, HI):2,         cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) = &
                 + P3 * this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)-3*D(ydim):cse(ydim, HI)-3*D(ydim), cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) &
                 + P2 * this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)-2*D(ydim):cse(ydim, HI)-2*D(ydim), cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) &
                 + P1 * this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)-  D(ydim):cse(ydim, HI)-  D(ydim), cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) &
                 + P0 * this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)          :cse(ydim, HI),           cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) &
                 + P_1* this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)+  D(ydim):cse(ydim, HI)+  D(ydim), cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) &
                 + P_2* this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)+2*D(ydim):cse(ydim, HI)+2*D(ydim), cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) &
                 + P_3* this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)+3*D(ydim):cse(ydim, HI)+3*D(ydim), cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z)
            this%prolong_xy           (fse(xdim, LO):fse(xdim, HI), fse(ydim, LO)+dom%D_y  :fse(ydim, HI)+dom%D_y:2, cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) = &
                 + P_3* this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)-3*D(ydim):cse(ydim, HI)-3*D(ydim), cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) &
                 + P_2* this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)-2*D(ydim):cse(ydim, HI)-2*D(ydim), cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) &
                 + P_1* this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)-  D(ydim):cse(ydim, HI)-  D(ydim), cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) &
                 + P0 * this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)          :cse(ydim, HI),           cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) &
                 + P1 * this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)+  D(ydim):cse(ydim, HI)+  D(ydim), cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) &
                 + P2 * this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)+2*D(ydim):cse(ydim, HI)+2*D(ydim), cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z) &
                 + P3 * this%prolong_x(fse(xdim, LO):fse(xdim, HI), cse(ydim, LO)+3*D(ydim):cse(ydim, HI)+3*D(ydim), cse(zdim, LO)-3*dom%D_z:cse(zdim, HI)+3*dom%D_z)
         case default
            call die("[grid_container:prolong] unsupported stencil size")
      end select

      select case (stencil_range*dom%D_z)
         case (I_ZERO)
            pa3d                  (fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), fse(zdim, LO):fse(zdim, HI):2) = &
                 this%prolong_xy  (fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO):cse(zdim, HI))
            if (dom%has_dir(zdim)) &
                 pa3d             (fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), fse(zdim, LO)+dom%D_z:fse(zdim, HI)+dom%D_z:2) = &
                 & this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO):cse(zdim, HI))
         case (I_ONE)
            pa3d                       (fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), fse(zdim, LO)        :fse(zdim, HI):2) = &
                 + P1 * this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)-D(zdim):cse(zdim, HI)-D(zdim)) &
                 + P0 * this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)        :cse(zdim, HI)        ) &
                 + P_1* this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)+D(zdim):cse(zdim, HI)+D(zdim))
            pa3d                       (fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), fse(zdim, LO)+dom%D_z:fse(zdim, HI)+dom%D_z:2) = &
                 + P_1* this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)-D(zdim):cse(zdim, HI)-D(zdim)) &
                 + P0 * this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)        :cse(zdim, HI)        ) &
                 + P1 * this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)+D(zdim):cse(zdim, HI)+D(zdim))
         case (I_TWO)
            pa3d                       (fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), fse(zdim, LO)          :fse(zdim, HI):2) = &
                 + P2 * this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)-2*D(zdim):cse(zdim, HI)-2*D(zdim)) &
                 + P1 * this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)-  D(zdim):cse(zdim, HI)-  D(zdim)) &
                 + P0 * this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)          :cse(zdim, HI)          ) &
                 + P_1* this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)+  D(zdim):cse(zdim, HI)+  D(zdim)) &
                 + P_2* this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)+2*D(zdim):cse(zdim, HI)+2*D(zdim))
            pa3d                       (fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), fse(zdim, LO)+dom%D_z  :fse(zdim, HI)+dom%D_z:2) = &
                 + P_2* this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)-2*D(zdim):cse(zdim, HI)-2*D(zdim)) &
                 + P_1* this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)-  D(zdim):cse(zdim, HI)-  D(zdim)) &
                 + P0 * this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)          :cse(zdim, HI)          ) &
                 + P1 * this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)+  D(zdim):cse(zdim, HI)+  D(zdim)) &
                 + P2 * this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)+2*D(zdim):cse(zdim, HI)+2*D(zdim))
         case (I_THREE)
            pa3d                       (fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), fse(zdim, LO)          :fse(zdim, HI):2) = &
                 + P3 * this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)-3*D(zdim):cse(zdim, HI)-3*D(zdim)) &
                 + P2 * this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)-2*D(zdim):cse(zdim, HI)-2*D(zdim)) &
                 + P1 * this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)-  D(zdim):cse(zdim, HI)-  D(zdim)) &
                 + P0 * this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)          :cse(zdim, HI)          ) &
                 + P_1* this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)+  D(zdim):cse(zdim, HI)+  D(zdim)) &
                 + P_2* this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)+2*D(zdim):cse(zdim, HI)+2*D(zdim)) &
                 + P_3* this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)+3*D(zdim):cse(zdim, HI)+3*D(zdim))
            pa3d                       (fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), fse(zdim, LO)+dom%D_z  :fse(zdim, HI)+dom%D_z:2) = &
                 + P_3* this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)-3*D(zdim):cse(zdim, HI)-3*D(zdim)) &
                 + P_2* this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)-2*D(zdim):cse(zdim, HI)-2*D(zdim)) &
                 + P_1* this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)-  D(zdim):cse(zdim, HI)-  D(zdim)) &
                 + P0 * this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)          :cse(zdim, HI)          ) &
                 + P1 * this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)+  D(zdim):cse(zdim, HI)+  D(zdim)) &
                 + P2 * this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)+2*D(zdim):cse(zdim, HI)+2*D(zdim)) &
                 + P3 * this%prolong_xy(fse(xdim, LO):fse(xdim, HI), fse(ydim, LO):fse(ydim, HI), cse(zdim, LO)+3*D(zdim):cse(zdim, HI)+3*D(zdim))
         case default
            call die("[grid_container:prolong] unsupported stencil size")
      end select
      ! Alternatively, an FFT convolution may be employed after injection. No idea at what stencil size the FFT is faster. It is finite size for sure :-)

   end subroutine prolong

!< \brief Create list of SFC indices to be created from refine flags

   subroutine refinemap2SFC_list(this)

      use constants,    only: refinement_factor, xdim, ydim, zdim, I_ONE
      use dataio_pub,   only: die
      use domain,       only: AMR_bsize
      use grid_helpers, only: c2f_o

      implicit none

      class(grid_container), intent(inout) :: this

      integer :: i, j, k, ifs, ife, jfs, jfe, kfs, kfe
      enum, bind(C)
         enumerator :: NONE, REFINE, LEAF
      end enum
      integer :: type

      this%refinemap = this%refinemap .and. this%leafmap
      type = NONE
      if (any(this%refinemap)) then
         type = REFINE
      else if (this%refine_flags%refine) then
         type = LEAF
      end if

      if (type ==NONE) return

      do i = int(((this%is - this%level_off(xdim))*refinement_factor) / AMR_bsize(xdim)), int(((this%ie - this%level_off(xdim))*refinement_factor + I_ONE) / AMR_bsize(xdim))
         ifs = max(this%is, (i*AMR_bsize(xdim))/refinement_factor)
         ife = min(this%ie, ((i+I_ONE)*AMR_bsize(xdim)-I_ONE)/refinement_factor)

         do j = int(((this%js - this%level_off(ydim))*refinement_factor) / AMR_bsize(ydim)), int(((this%je - this%level_off(ydim))*refinement_factor + I_ONE) / AMR_bsize(ydim))
            jfs = max(this%js, (j*AMR_bsize(ydim))/refinement_factor)
            jfe = min(this%je, ((j+I_ONE)*AMR_bsize(ydim)-I_ONE)/refinement_factor)

            do k = int(((this%ks - this%level_off(zdim))*refinement_factor) / AMR_bsize(zdim)), int(((this%ke - this%level_off(zdim))*refinement_factor + I_ONE) / AMR_bsize(zdim))
               kfs = max(this%ks, (k*AMR_bsize(zdim))/refinement_factor)
               kfe = min(this%ke, ((k+I_ONE)*AMR_bsize(zdim)-I_ONE)/refinement_factor)

               select case (type)
                  case (REFINE)
                     if (any(this%refinemap(ifs:ife, jfs:jfe, kfs:kfe))) call this%refine_flags%add(this%level_id+1, c2f_o([i, j, k]*AMR_bsize-this%level_off))
                  case (LEAF)
                     if (all(this%leafmap(ifs:ife, jfs:jfe, kfs:kfe))) then
                        call this%refine_flags%add(this%level_id+1, c2f_o([i, j, k]*AMR_bsize-this%level_off))
                     else if (any(this%leafmap(ifs:ife, jfs:jfe, kfs:kfe))) then
                        call die("[grid_container:refinemap2SFC_list] cannot refine partially leaf parf of the grid")
                     end if
                  case default
                     call die("[grid_container:refinemap2SFC_list] invalid type")
               end select
            end do
         end do
      end do
      this%refinemap = .false.

   end subroutine refinemap2SFC_list

end module grid_cont
