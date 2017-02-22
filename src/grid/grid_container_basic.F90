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

!> \brief Module containing the grid container type and its most basic associated methods

module grid_cont_basic

   use constants,        only: xdim, zdim, ndims, LO, HI, CENTER, INV_CENTER
   use fluxtypes,        only: fluxarray
   use level_essentials, only: level_T
   use named_array,      only: named_array4d, named_array3d, named_array_fc, fc_3d_arr
   use real_vector,      only: real_vec_T

   implicit none

   private
   public :: grid_container_basic, tgt_list, segment, get_cs

   !> \brief Specification of segment of data for boundary exchange
   type :: segment
      integer :: proc                                     !< target process
      integer(kind=8), dimension(xdim:zdim, LO:HI) :: se  !< range
      integer(kind=4) :: tag                              !< unique tag for data exchange
      real, allocatable, dimension(:,:,:)   :: buf        !< buffer for the 3D (scalar) data to be sent or received
      real, allocatable, dimension(:,:,:,:) :: buf4       !< buffer for the 4D (vector) data to be sent or received
      integer(kind=4), pointer :: req                     !< request ID, used for most asynchronous communication, such as fine-coarse flux exchanges
      integer(kind=8), dimension(xdim:zdim, LO:HI) :: se2 !< auxiliary range, used in cg_level_connected:vertical_bf_prep
      class(grid_container_basic), pointer :: local       !< set this pointer to non-null when the exchange is local
   end type segment

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
   type :: grid_container_basic

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
      class(level_T), pointer :: l                               !< level essential data

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
      integer(kind=4), dimension(ndims, LO:HI)  :: lh_out        !< ijkse expanded at the external boundaries to include external guardcells for contexts where AT_OUT_B is used
      integer(kind=4), dimension(ndims)         :: n_            !< number of %grid cells in one block in x-, y- and z-directions (n_b(:) + 2 * nb)
      integer(kind=8), dimension(ndims, LO:HI)  :: my_se         !< own segment. my_se(:,LO) = 0; my_se(:,HI) = dom%n_d(:) - 1 would cover entire domain on a base level
                                                                 !! my_se(:,LO) = 0; my_se(:,HI) = finest%level%n_d(:) -1 would cover entire domain on the most refined level
                                                                 !! DEPRECATED: will be equivalent to ijkse(:,:)

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

      ! Non-cartesian geometrical factors

      real, allocatable, dimension(:,:,:) :: gc_xdim             !< array of geometrical coefficients in x-direction
      real, allocatable, dimension(:,:,:) :: gc_ydim             !< array of geometrical coefficients in y-direction
      real, allocatable, dimension(:,:,:) :: gc_zdim             !< array of geometrical coefficients in z-direction

      ! Registered variables

      type(named_array3d),  allocatable, dimension(:) :: q       !< 3D arrays such as gravitational potential pr user-defined quantities or gravitational potential
      type(named_array4d),  allocatable, dimension(:) :: w       !< 4D arrays such as u, vector fields (b) or other vector/multi-scalar user-defined quantities
      type(named_array_fc), allocatable, dimension(:) :: f       !< face-centered arrays such as magnetic field

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
      type(fc_3d_arr), dimension(:), pointer :: bf => null() !< Main array of face-centered magnetic field components

      ! Misc
      integer(kind=8) :: SFC_id                                  !< position of the grid on space-filling curve
      type(mg_arr), pointer :: mg                                !< multigrid arrays
      real :: vol                                                !< volume of the grid; BEWARE: for cylindrical geometry it needs to be integrated over x(:) to get real volume
      real :: dxmn                                               !< the smallest length of the %grid cell (among dx, dy, and dz)
      integer(kind=4) :: maxxyz                                  !< maximum number of %grid cells in any direction
      integer :: grid_id                                         !< index of own segment in own level decomposition, e.g. my_se(:,:) = base%level%dot%gse(proc)%c(grid_id)%se(:,:)
      integer :: membership                                      !< How many cg lists use this grid piece?

   contains

      procedure          :: init_gc_basic                        !< Initialization
      procedure          :: cleanup_basic                        !< Deallocate all internals
      procedure, private :: set_coords                           !< Calculate arrays of coordinates along a given direction
      procedure, private :: add_all_na                           !< Register all known named arrays for this cg, sey up shortcuts to the crucial fields
      procedure          :: add_na                               !< Register a new 3D entry in current cg with given name.
      procedure          :: add_na_4d                            !< Register a new 4D entry in current cg with given name.
      procedure          :: add_na_fc                            !< Register a new face_centered entry in current cg with given name.
      procedure          :: set_constant_b_field                 !< set constant magnetic field on whole block
      procedure          :: emag_point                           !< return energy asociated with magnetic field at specified point
      procedure          :: emag_range                           !< return energy asociated with magnetic field at specified range
      generic, public    :: emag => emag_point, emag_range

   end type grid_container_basic

contains

!>
!! \brief Initialization of the grid container basic
!!
!! \details This method sets up the grid container variables, coordinates and allocates basic arrays.
!! Everything related to the interior of grid container should be set here.
!! Things that are related to communication with other grid containers or global properties are set up in cg_level::init_all_new_cg.
!!
!<

   subroutine init_gc_basic(this, my_se, grid_id, l)

      use constants,        only: PIERNIK_INIT_DOMAIN, xdim, ydim, zdim, LO, HI, I_ONE, I_TWO, BND_MPI, BND_COR, GEO_XYZ, GEO_RPZ, dpi
      use dataio_pub,       only: die, warn, code_progress
      use domain,           only: dom
      use func,             only: operator(.equals.)
      use level_essentials, only: level_T
      use ordering,         only: SFC_order

      implicit none

      class(grid_container_basic), target, intent(inout) :: this  ! intent(out) would silently clear everything, that was already set
                                                                  ! (also the fields in types derived from grid_container)
      integer(kind=8), dimension(:,:),     intent(in) :: my_se    !< my segment
      integer,                             intent(in) :: grid_id  !< ID which should be unique across level
      class(level_T), pointer,             intent(in) :: l        !< level essential data

      integer :: i

      if (code_progress < PIERNIK_INIT_DOMAIN) call die("[grid_container_basic:init_gc] MPI not initialized.")

      this%l          => l
      this%membership = 1
      this%grid_id    = grid_id
      this%my_se(:,:) = my_se(:, :)
      this%n_b(:)     = int(this%my_se(:, HI) - this%my_se(:, LO) + I_ONE, 4) ! Block 'physical' grid sizes
      this%SFC_id     = SFC_order(this%my_se(:, LO) - l%off)

      if (any(this%n_b(:) <= 0)) call die("[grid_container_basic:init_gc] Mixed positive and non-positive grid sizes")

      ! Inherit the boundaries from the domain, then set MPI or SHEAR boundaries where applicable
      this%bnd(:,:) = dom%bnd(:,:)
      where (my_se(:, LO)         /= l%off(:)           ) this%bnd(:, LO) = BND_MPI
      where (my_se(:, HI) + I_ONE /= l%off(:) + l%n_d(:)) this%bnd(:, HI) = BND_MPI
      ! For periodic boundaries do not set BND_MPI when local domain spans through the whole computational domain in given direction.
      where (dom%periodic(:) .and. this%my_se(:, HI) + I_ONE /= l%n_d(:)) this%bnd(:, LO) = BND_MPI
      where (dom%periodic(:) .and. this%my_se(:, LO)         /= 0       ) this%bnd(:, HI) = BND_MPI

      this%ext_bnd(:, :) = .false.
      do i = xdim, zdim
         if (dom%has_dir(i) .and. .not. dom%periodic(i)) then
            this%ext_bnd(i, LO) = (my_se(i, LO)         == l%off(i))
            this%ext_bnd(i, HI) = (my_se(i, HI) + I_ONE == l%off(i) + l%n_d(i))
         endif
      enddo

      ! For shear boundaries and some domain decompositions it is possible that a boundary can be mixed 'per' with 'mpi'

!      call inflate_req
      ! write_plot_hdf5 requires nproc entries for the status array

      if (any(dom%bnd(xdim:ydim, :) == BND_COR)) call die("[grid_container_basic:init_gc] BND_COR unimplemented")
      if (any(dom%bnd(zdim, :) == BND_COR)) call die("[grid_container_basic:init_gc] Corner BC not allowed for z-direction")

#ifdef SHEAR_BND
      call die("[grid_container_basic:init_gc] Shear-pediodic boundary conditions unimplemented")
      ! This is possible to be implemented
#endif /* SHEAR_BND */

      do i = xdim, zdim
         if (dom%has_dir(i)) then
            if (this%n_b(i) < 1) call die("[grid_container_basic:init_gc] Too many CPUs for such a small grid.")
            if (this%n_b(i) < dom%nb) call warn("[grid_container_basic:init_gc] domain size in some directions is < nb, which may result in incomplete boundary cell update")
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
         this%dl(:)        = dom%L_(:) / l%n_d(:)
         this%fbnd(:, LO)  = dom%edge(:, LO) + this%dl(:) * (this%my_se(:, LO)         - l%off(:))
         this%fbnd(:, HI)  = dom%edge(:, LO) + this%dl(:) * (this%my_se(:, HI) + I_ONE - l%off(:))
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

      ! Compute indices that include external boundary cells
      ! Strangely, we ignore periodicity here, following what we had in restart_hdf5_v1::set_dims_for_restart

      this%lh_out = this%ijkse
      where (dom%has_dir(:) .and. (my_se(:, LO) == l%off(:)                   )) this%lh_out(:, LO) = this%lh_out(:, LO) - dom%nb
      where (dom%has_dir(:) .and. (my_se(:, HI) == l%off(:) + l%n_d(:) - I_ONE)) this%lh_out(:, HI) = this%lh_out(:, HI) + dom%nb
      !> \todo make sure the above works correctly with refinements

      if (any(this%dl .equals. 0.)) call die("[grid_container_basic:init_gc] found cell size equal to 0.")

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

   end subroutine init_gc_basic

!> \brief Calculate arrays of coordinates along a given direction

   subroutine set_coords(this)

      use constants,  only: LO, HI, half, one, zero, xdim, ydim, zdim, CENTER, LEFT, RIGHT, INV_CENTER
      use dataio_pub, only: die, warn
      use domain,     only: dom
      use func,       only: operator(.notequals.), operator(.equals.)

      implicit none

      class(grid_container_basic), intent(inout) :: this !< grid container, where the arrays have to be set

      integer :: d, i
      integer, parameter :: safety_warn_factor = 1000 ! warn if a cell size is smaller than this * epsilon(coordinates)

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

         ! Generally nobody should substract one cell coordinate from another in code solvers. One should use cell sizes instead.
         ! The problem may arise when initial conditions are comparing coordinates to set something on the left or right side of some line.
         ! When the cell size is too small compared to the coordinates, such line cannot be properly calculated
         ! Note that since we force real kind=8, we can use a named constant instead of epsilon
         if ( any(this%coord(CENTER, d)%r(:) .equals. this%coord(LEFT,  d)%r(:)) .or. &
              any(this%coord(CENTER, d)%r(:) .equals. this%coord(RIGHT, d)%r(:)) ) call die("[grid_container:set_coords] cannot distinguish between center and face coordinates of a cell")
         if ( any(abs(this%coord(CENTER, d)%r(:)-this%coord(LEFT,  d)%r(:)) < safety_warn_factor*epsilon(this%coord(CENTER, d)%r(:))*this%coord(CENTER, d)%r(:)) .or. &
              any(abs(this%coord(CENTER, d)%r(:)-this%coord(RIGHT, d)%r(:)) < safety_warn_factor*epsilon(this%coord(CENTER, d)%r(:))*this%coord(CENTER, d)%r(:))) &
              call warn("[grid_container:set_coords] cell sizes are much smaller than coordinates. Inaccuracies in setting the initial conditions may happen.")
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

!> \brief Routine that deallocates all internals of the grid container basic

   subroutine cleanup_basic(this)

      use constants, only: xdim, zdim, CENTER, INV_CENTER, LO, HI

      implicit none

      class(grid_container_basic), intent(inout) :: this !< object invoking type-bound procedure

      integer :: d, g, b, cdim

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

      if (allocated(this%f)) then
         do g = lbound(this%f(:), dim=1), ubound(this%f(:), dim=1)
            call this%f(g)%clean
         enddo
         deallocate(this%f)
      endif

      do d = xdim, zdim
         do g = LO, HI
            call this%finebnd  (d, g)%facleanup
            call this%coarsebnd(d, g)%facleanup
         enddo
      enddo

   end subroutine cleanup_basic

!> \brief Register all known named arrays for this cg, sey up shortcuts to the crucial fields

   subroutine add_all_na(this)

      use constants,        only: INVALID, base_level_id
      use dataio_pub,       only: die, warn
      use named_array_list, only: qna, wna, fna
#ifdef ISO
      use constants,   only: cs_i2_n
      use fluids_pub,  only: cs2_max
#endif /* ISO */

      implicit none

      class(grid_container_basic), intent(inout), target :: this  !< object invoking type-bound procedure

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

      if (allocated(fna%lst)) then
         do i = lbound(fna%lst(:), dim=1), ubound(fna%lst(:), dim=1)
            call this%add_na_fc()
         enddo
      endif

      ! shortcuts
      this%u  => this%w(wna%fi)%arr
      if (wna%bi /= INVALID) this%b  => this%w(wna%bi)%arr
      if (fna%bi /= INVALID) this%bf => this%f(fna%bi)%f_arr
      if (associated(this%b) .and. associated(this%bf)) call die("[grid_container:add_all_na] both positions of magnetic field can not be used at same time")
      ! if one ever needs other type of array for magnetic field storage (auxiliary one for interpolated values), it can be requested from user problem or solver but should not become a standard.
      if ((.not. associated(this%b)) .and. (.not. associated(this%bf)) .and. (this%l%id >= base_level_id)) call warn("[grid_container:add_all_na] no magnetic field was declared")
      !> \todo use this by default for nonmagnetic simulations  instead of having MAGNETIC preprocessor macro.
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

      class(grid_container_basic), intent(inout) :: this          !< object invoking type-bound procedure
      logical,                     intent(in)    :: multigrid     !< If .true. then cg%q(:)%arr and cg%w(:)%arr are allocated also below base level

      type(named_array3d), allocatable, dimension(:) :: tmp

      if (.not. allocated(this%q)) then
         allocate(this%q(1))
      else
         allocate(tmp(lbound(this%q(:),dim=1):ubound(this%q(:), dim=1) + 1))
         tmp(:ubound(this%q(:), dim=1)) = this%q(:)
         call move_alloc(from=tmp, to=this%q)
      endif

      if (multigrid .or. this%l%id >= base_level_id) call this%q(ubound(this%q(:), dim=1))%init(this%lhn(:, LO), this%lhn(:, HI))

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

      class(grid_container_basic), intent(inout) :: this         !< object invoking type-bound procedure
      integer(kind=4),             intent(in)    :: n            !< Length of the vector quantity to be stored (first dimension of the array)

      type(named_array4d), allocatable, dimension(:) :: tmp

      if (.not. allocated(this%w)) then
         allocate(this%w(1))
      else
         allocate(tmp(lbound(this%w(:),dim=1):ubound(this%w(:), dim=1) + 1))
         tmp(:ubound(this%w(:), dim=1)) = this%w(:)
         call move_alloc(from=tmp, to=this%w)
      endif

      if (this%l%id >= base_level_id) call this%w(ubound(this%w(:), dim=1))%init( [1_INT4, this%lhn(:, LO)], [n, this%lhn(:, HI)] ) !< \deprecated magic integer

   end subroutine add_na_4d

!>
!! \brief Register a new fc entry in current cg with given name. Called from cg_list_glob::reg_var
!!
!! \warning This routine should not be called directly from user routines
!<
   subroutine add_na_fc(this)

      use constants,   only: base_level_id
      use named_array, only: named_array_fc

      implicit none

      class(grid_container_basic), intent(inout) :: this

      type(named_array_fc), allocatable, dimension(:) :: tmp

      if (.not. allocated(this%f)) then
         allocate(this%f(1))
      else
         allocate(tmp(lbound(this%f(:),dim=1):ubound(this%f(:), dim=1) + 1))
         tmp(:ubound(this%f(:), dim=1)) = this%f(:)
         call move_alloc(from=tmp, to=this%f)
      endif

      if (this%l%id >= base_level_id) call this%f(ubound(this%f(:), dim=1))%init(this%lhn(:, LO), this%lhn(:, HI))

   end subroutine add_na_fc

!> \brief Add a new segment, reallocate if necessary

   subroutine add_seg(this, proc, se, tag)

      use dataio_pub, only: die

      implicit none

      class(bnd_list),                              intent(inout) :: this !< object invoking type-bound procedure
      integer,                                      intent(in)    :: proc !< process to be communicated
      integer(kind=8), dimension(xdim:zdim, LO:HI), intent(in)    :: se   !< segment definition
      integer(kind=4),                              intent(in)    :: tag  !< tag for MPI calls

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

!< \brief set constant magnetic field on whole block

   subroutine set_constant_b_field(this, b)

      use constants, only: xdim, ydim, zdim, I_ONE

      implicit none

      class(grid_container_basic), intent(inout) :: this !< object invoking type-bound procedure
      real, dimension(xdim:zdim),  intent(in)    :: b    !< the value of the magnetic field vector in whole block

      integer :: d

      if (associated(this%b)) then
         do d = xdim, zdim
            this%b(d, this%is:this%ie, this%js:this%je, this%ks:this%ke) = b(d)
         enddo
      endif

      if (associated(this%bf)) then
         this%bf(xdim)%arr(this%is:this%ie+I_ONE, this%js:this%je,       this%ks:this%ke      ) = b(xdim)
         this%bf(ydim)%arr(this%is:this%ie,       this%js:this%je+I_ONE, this%ks:this%ke      ) = b(ydim)
         this%bf(zdim)%arr(this%is:this%ie,       this%js:this%je,       this%ks:this%ke+I_ONE) = b(zdim)
      endif

   end subroutine set_constant_b_field

!< \brief return energy asociated with magnetic field at specified point

   function emag_point(this, ijk) result(e_mag)

      use constants,  only: xdim, ydim, zdim, I_ONE
      use dataio_pub, only: die
      use func,       only: emag

      implicit none

      class(grid_container_basic), intent(in)  :: this
      integer, dimension(:),       intent(in)  :: ijk

      real :: e_mag

      if (associated(this%b)) then
         e_mag = emag(this%b(xdim, ijk(xdim), ijk(ydim), ijk(zdim)), &
              &       this%b(ydim, ijk(xdim), ijk(ydim), ijk(zdim)), &
              &       this%b(zdim, ijk(xdim), ijk(ydim), ijk(zdim)))
      else if (associated(this%bf)) then
         e_mag = emag((this%bf(xdim)%arr(ijk(xdim), ijk(ydim), ijk(zdim)) + this%bf(xdim)%arr(ijk(xdim)+I_ONE, ijk(ydim), ijk(zdim)))/2., &
              &       (this%bf(ydim)%arr(ijk(xdim), ijk(ydim), ijk(zdim)) + this%bf(ydim)%arr(ijk(xdim), ijk(ydim)+I_ONE, ijk(zdim)))/2., &
              &       (this%bf(zdim)%arr(ijk(xdim), ijk(ydim), ijk(zdim)) + this%bf(zdim)%arr(ijk(xdim), ijk(ydim), ijk(zdim)+I_ONE))/2.)
      else
         call die("[grid_container:emag_point] no magnetic field declared here")
         e_mag = 0.
      endif

   end function emag_point

!< \brief return energy asociated with magnetic field at specified range

   function emag_range(this, ijk) result(e_mag)

      use constants,  only: xdim, ydim, zdim, LO, HI, I_ONE
      use dataio_pub, only: die
      use func,       only: emag

      implicit none

      class(grid_container_basic), intent(in)  :: this
      integer, dimension(:,:),     intent(in)  :: ijk

      real, dimension(ijk(xdim, LO):ijk(xdim, HI), ijk(ydim, LO):ijk(ydim, HI), ijk(zdim, LO):ijk(zdim, HI)) :: e_mag

      if (associated(this%b)) then
         e_mag = emag(this%b(xdim, ijk(xdim, LO):ijk(xdim, HI), ijk(ydim, LO):ijk(ydim, HI), ijk(zdim, LO):ijk(zdim, HI)), &
              &       this%b(ydim, ijk(xdim, LO):ijk(xdim, HI), ijk(ydim, LO):ijk(ydim, HI), ijk(zdim, LO):ijk(zdim, HI)), &
              &       this%b(zdim, ijk(xdim, LO):ijk(xdim, HI), ijk(ydim, LO):ijk(ydim, HI), ijk(zdim, LO):ijk(zdim, HI)))
      else if (associated(this%bf)) then
         e_mag = emag((this%bf(xdim)%arr(ijk(xdim, LO)      :ijk(xdim, HI),       ijk(ydim, LO):ijk(ydim, HI), ijk(zdim, LO):ijk(zdim, HI)) + &
              &        this%bf(xdim)%arr(ijk(xdim, LO)+I_ONE:ijk(xdim, HI)+I_ONE, ijk(ydim, LO):ijk(ydim, HI), ijk(zdim, LO):ijk(zdim, HI)))/2., &
              &       (this%bf(ydim)%arr(ijk(xdim, LO):ijk(xdim, HI), ijk(ydim, LO)      :ijk(ydim, HI),       ijk(zdim, LO):ijk(zdim, HI)) + &
              &        this%bf(ydim)%arr(ijk(xdim, LO):ijk(xdim, HI), ijk(ydim, LO)+I_ONE:ijk(ydim, HI)+I_ONE, ijk(zdim, LO):ijk(zdim, HI)))/2., &
              &       (this%bf(zdim)%arr(ijk(xdim, LO):ijk(xdim, HI), ijk(ydim, LO):ijk(ydim, HI), ijk(zdim, LO)      :ijk(zdim, HI)      ) + &
              &        this%bf(zdim)%arr(ijk(xdim, LO):ijk(xdim, HI), ijk(ydim, LO):ijk(ydim, HI), ijk(zdim, LO)+I_ONE:ijk(zdim, HI)+I_ONE))/2.)
      else
         call die("[grid_container:emag_range] no magnetic field declared here")
         e_mag = 0.
      endif

   end function emag_range


   real function get_cs(this, i, j, k, u, b, cs_iso2) ! ion_cs

      use constants,  only: two
      use fluidtypes, only: component_fluid
#ifndef ISO
      use func,      only: ekin
#endif /* !ISO */
#ifdef MAGNETIC
      use constants, only: xdim, ydim, zdim, half
      use domain,    only: dom
      use func,      only: emag
#else /* !MAGNETIC */
      use constants, only: zero
#endif /* !MAGNETIC */

      implicit none

      class(component_fluid),            intent(in) :: this
      integer,                           intent(in) :: i, j, k
      real, dimension(:,:,:,:), pointer, intent(in) :: u       !< pointer to array of fluid properties
      real, dimension(:,:,:,:), pointer, intent(in) :: b       !< pointer to array of magnetic fields (used for ionized fluid with MAGNETIC #defined)
      real, dimension(:,:,:),   pointer, intent(in) :: cs_iso2 !< pointer to array of isothermal sound speeds (used when ISO was #defined)

#ifdef MAGNETIC
      real :: bx, by, bz
#endif /* MAGNETIC */
      real :: pmag, p, ps

#ifdef MAGNETIC
      bx = half*(b(xdim,i,j,k) + b(xdim, i+dom%D_x, j,         k        ))
      by = half*(b(ydim,i,j,k) + b(ydim, i,         j+dom%D_y, k        ))
      bz = half*(b(zdim,i,j,k) + b(zdim, i,         j,         k+dom%D_z))

      pmag = emag(bx, by, bz)
#else /* !MAGNETIC */
      ! all_mag_boundaries has not been called so we cannot trust b(xdim, ie+dom%D_x:), b(ydim,:je+dom%D_y and b(zdim,:,:, ke+dom%D_z
      pmag = zero
#endif /* !MAGNETIC */

#ifdef ISO
      p  = cs_iso2(i, j, k) * u(this%idn, i, j, k)
      ps = p + pmag
      get_cs = sqrt(abs((two * pmag + p) / u(this%idn, i, j, k)))
#else /* !ISO */
      ps = (u(this%ien, i, j, k) - &
         &   ekin(u(this%imx, i, j, k), u(this%imy, i, j, k), u(this%imz, i, j, k), u(this%idn, i, j, k)) &
         & ) * (this%gam_1) + (two - this%gam) * pmag
      p  = ps - pmag
      get_cs = sqrt(abs((two * pmag + this%gam * p) / u(this%idn, i, j, k)))
#endif /* !ISO */
      if (.false.) print *, u(:, i, j, k), b(:, i, j, k), cs_iso2(i, j, k), this%cs
   end function get_cs

end module grid_cont_basic
