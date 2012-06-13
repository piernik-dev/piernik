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

!> \brief Module containing the grid container type and its associated methods

module grid_cont

   use constants,   only: xdim, zdim, ndims, LO, HI
   use types,       only: axes
   use named_array, only: named_array4d, named_array3d, mbc_list

   implicit none

   private
   public :: grid_container, segment, bnd_list, pr_segment, tgt_list, set_mpi_types, is_overlap

   !> \brief Specification of segment of data for boundary exchange, prolongation and restriction.
   type :: segment
      integer :: proc                                     !< target process
      integer(kind=8), dimension(xdim:zdim, LO:HI) :: se  !< range
   end type segment

   !> \brief coefficient-layer pair used for prolongation
   type :: c_layer
      integer :: layer                                                !< index of a layer with face-prolongation coefficient coeff
      real    :: coeff                                                !< coefficient for face prolongation
   end type c_layer

   !> \brief segment type for prolongation and restriction
   type, extends(segment) :: pr_segment
      real, allocatable, dimension(:,:,:) :: buf                      !< buffer for the coarse data (incoming prolongation and outgoing restriction) for each nonlocal operations
      type(c_layer), dimension(:), allocatable :: f_lay               !< face layers to contribute to the prolonged face value
      integer(kind=4) :: tag                                          !< unique tag for communication
   end type pr_segment                                                !< (not allocated for outgoing prolongation, incoming restriction and for local operations)

   !< \brief target list container for prolongations, restrictions and boundary exchanges
   type :: tgt_list
      type(pr_segment), dimension(:), allocatable :: seg              !< a segment of data to be received or sent
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
      integer                                :: fft_type              !< type of FFT to employ (depending on boundaries)
      complex, allocatable, dimension(:,:,:) :: fft                   !< a complex array for FFT operations (Fourier space)
      real,    allocatable, dimension(:,:,:) :: fftr                  !< a real array for FFT operations (Fourier space for sine transform)
      real,    allocatable, dimension(:,:,:) :: src                   !< an input array for FFT (real space data)
      real,    allocatable, dimension(:,:,:) :: Green3D               !< Green's function (0.5 * fft_norm / (kx + ky + kz))
      integer (kind = selected_int_kind(16)) :: planf, plani          !< FFT forward and inverse plans
      real                                   :: fft_norm              !< normalization factor

      ! geometrical factors
      !! \todo move to cg
      real    :: r, rx, ry, rz                                        !< geometric factors for relaxation (diffusion) used in approximate_solution_rbgs

      ! prolongation and restriction
      !! \todo move to cg, should be initialized by cg_list_level procedure
      type(tgt_list), dimension(xdim:zdim, LO:HI) :: pff_tgt, pfc_tgt !< description outgoing and incoming face prolongation data

   end type mg_arr

   !> \brief Segment type with additional parameters for boundary exchange
   type, extends(segment) :: bnd_segment
      integer(kind=4) :: tag                              !< unique tag for data exchange
   end type bnd_segment

   !> \brief Array of boundary segments to exchange
   type :: bnd_list
      type(bnd_segment), dimension(:), allocatable :: seg !< segments
   end type bnd_list

   !> \brief Everything required for autonomous computation of a single sweep on a portion of the domain on a single process
   type, extends(axes) :: grid_container

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
      real :: dvol2                                              !< square of one cell volume

      ! Grid properties

      ! cell count and position
      integer(kind=8), dimension(ndims) :: off                   !< offset of the local domain within computational domain
      integer(kind=4), dimension(ndims) :: n_b                   !< [nxb, nyb, nzb]
      integer(kind=4) :: nxb                                     !< number of %grid cells in one block (without boundary cells) in x-direction
      integer(kind=4) :: nyb                                     !< number of %grid cells in one block (without boundary cells) in y-direction
      integer(kind=4) :: nzb                                     !< number of %grid cells in one block (without boundary cells) in z-direction
      integer :: level_id                                        !< level id (number)

      ! shortcuts
      integer(kind=4) :: is                                      !< index of the first %grid cell of physical domain in x-direction
      integer(kind=4) :: ie                                      !< index of the last %grid cell of physical domain in x-direction
      integer(kind=4) :: js                                      !< index of the first %grid cell of physical domain in y-direction
      integer(kind=4) :: je                                      !< index of the last %grid cell of physical domain in y-direction
      integer(kind=4) :: ks                                      !< index of the first %grid cell of physical domain in z-direction
      integer(kind=4) :: ke                                      !< index of the last %grid cell of physical domain in z-direction
      integer(kind=4) :: isb, ieb, jsb, jeb, ksb, keb            !< auxiliary indices for exchanging boundary data, (e.g. is:isb -> ie+1:nx, ieb:ie -> 1:nb)
      integer(kind=4), dimension(ndims, LO:HI)  :: ijkse         !< [[is, js, ks], [ie, je, ke]]
      integer(kind=8), dimension(ndims) :: h_cor1                !< offsets of the corner opposite to the one defined by off(:) + 1, a shortcut to be compared with dom%n_d(:)
      integer(kind=4), dimension(ndims) :: n_                    !< number of %grid cells in one block in x-, y- and z-directions (n_b(:) + 2 * nb)
      integer(kind=8), dimension(ndims, LO:HI) :: my_se          !< own segment

      ! Physical size and coordinates

      real, dimension(ndims, LO:HI) :: fbnd                      !< current block boundary positions
      real, allocatable, dimension(:) :: inv_x                   !< array of invert x-positions of %grid cells centers
      real, allocatable, dimension(:) :: inv_y                   !< array of invert y-positions of %grid cells centers
      real, allocatable, dimension(:) :: inv_z                   !< array of invert z-positions of %grid cells centers
      real, allocatable, dimension(:) :: xl                      !< array of x-positions of %grid cells left borders
      real, allocatable, dimension(:) :: yl                      !< array of y-positions of %grid cells left borders
      real, allocatable, dimension(:) :: zl                      !< array of z-positions of %grid cells left borders
      real, allocatable, dimension(:) :: xr                      !< array of x-positions of %grid cells right borders
      real, allocatable, dimension(:) :: yr                      !< array of y-positions of %grid cells right borders
      real, allocatable, dimension(:) :: zr                      !< array of z-positions of %grid cells right borders

      ! External boundary conditions and internal boundaries

      integer(kind=4), dimension(ndims, LO:HI)  :: bnd           !< type of boundary conditions coded in integers
      integer(kind=4), allocatable, dimension(:,:,:,:,:) :: mbc  !< MPI Boundary conditions Container for comm3d-based communication
      type(bnd_list), dimension(:,:), allocatable :: i_bnd       !< description of incoming boundary data, the shape is (xdim:zdim, nb)
      type(bnd_list), dimension(:,:), allocatable :: o_bnd       !< description of outgoing boundary data, the shape is (xdim:zdim, nb)
      logical, dimension(xdim:zdim, LO:HI) :: ext_bnd            !< .false. for BND_PER and BND_MPI

      ! Prolongation and restriction

      type(tgt_list) :: f_tgt                                    !< description of incoming restriction and outgoing prolongation data (this should be a linked list)
      type(tgt_list) :: c_tgt                                    !< description of outgoing restriction and incoming prolongation data
      real, allocatable, dimension(:,:,:) :: prolong_, prolong_x, prolong_xy !< auxiliary prolongation arrays

      ! Non-cartesian geometrical factors

      real, allocatable, dimension(:,:,:) :: gc_xdim             !< array of geometrical coefficients in x-direction
      real, allocatable, dimension(:,:,:) :: gc_ydim             !< array of geometrical coefficients in y-direction
      real, allocatable, dimension(:,:,:) :: gc_zdim             !< array of geometrical coefficients in z-direction

      ! Registeed variables

      type(named_array3d), allocatable, dimension(:) :: q        !< 3D arrays such as gravitational potential pr user-defined quantities or gravitational potential
      type(named_array4d), allocatable, dimension(:) :: w        !< 4D arrays such as u, vector fields (b) or other vector/multi-scalar user-defined quantities

      type(mbc_list), dimension(:,:), allocatable :: q_i_mbc     !< MPI Boundary conditions Containers for incoming guardcell updates on the q arrays
      type(mbc_list), dimension(:,:), allocatable :: q_o_mbc     !< MPI Boundary conditions Containers for outgoing guardcell updates on the q arrays

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

      type(mg_arr) :: mg                                         !< multigrid arrays (without multigrid will remain unallocated)
      real :: vol                                                !< volume of the grid; BEWARE: for cylindrical geometry it need to be multiplied by appropriate x(:) to get real volume
      real :: dxmn                                               !< the smallest length of the %grid cell (among dx, dy, and dz)
      integer(kind=4) :: maxxyz                                  !< maximum number of %grid cells in any direction
      integer :: grid_id                                         !< index of own segment in own level decomposition, e.g. my_se(:,:) = base_lev%pse(proc)%sel(grid_id, :, :)

   contains

      procedure :: init                                          !< Initialization
      procedure :: cleanup                                       !< Deallocate all internals
      procedure :: set_axis                                      !< Calculate arrays of coordinates along a given direction
      procedure :: set_q_mbc                                     !< Initialize the communicators for q

   end type grid_container

   interface is_overlap
      module procedure is_overlap_simple, is_overlap_per
   end interface

contains

!>
!! \brief Initialization of the grid container
!!
!! \details This method sets up the grid container variables, coordinates and allocates basic arrays
!<

   subroutine init(this, n_d, my_se, grid_id, level_id)

      use constants,  only: PIERNIK_INIT_DOMAIN, xdim, ydim, zdim, ndims, FLUID, ARR, LO, HI, BND, BLK, INVALID, I_ONE, I_TWO, BND_MPI, BND_COR, big_float
#ifdef SHEAR_BND
#ifndef FFTW
      use constants,  only: BND_SHE
#endif /* !FTTW */
#endif /* SHEAR_BND */
      use dataio_pub, only: die, warn, printinfo, msg, code_progress
      use domain,     only: dom
      use mpi,        only: MPI_COMM_NULL
      use mpisetup,   only: nproc, inflate_req
      use types,      only: cdd
#ifdef DEBUG
      use mpisetup,   only: proc
#endif /* DEBUG */

      implicit none

      class(grid_container), target,   intent(inout) :: this ! intent(out) would silently clear everything, that was already set
                                                             ! (also the fields in types derived from grid_container)
      integer(kind=8), dimension(:),   intent(in) :: n_d     !< max resolution of my level
      integer(kind=8), dimension(:,:), intent(in) :: my_se !< my segment
      integer,                         intent(in) :: grid_id
      integer(kind=4),                 intent(in) :: level_id

      integer :: i

      if (code_progress < PIERNIK_INIT_DOMAIN) call die("[grid_container:init] MPI not initialized.")

      this%grid_id    = grid_id
      this%my_se(:,:) = my_se(:, :)
      this%off(:)     = this%my_se(:, LO)
      this%h_cor1(:)  = this%my_se(:, HI) + I_ONE
      this%n_b(:)     = int(this%h_cor1(:) - this%off(:), 4) ! Block 'physical' grid sizes
      this%level_id   = level_id

      if (any(this%n_b(:) <= 0)) call die("[grid_container:init] Mixed positive and non-positive grid sizes")

      ! Inherit the boundaries from the domain, then set MPI or SHEAR boundaries where applicable
      this%bnd(:,:) = dom%bnd(:,:)
      where (this%off(:)    /= 0)      this%bnd(:, LO) = BND_MPI
      where (this%h_cor1(:) /= n_d(:)) this%bnd(:, HI) = BND_MPI
      ! For periodic boundaries do not set BND_MPI when local domain spans through the whole computational domain in given direction.
      where (dom%periodic(:) .and. this%h_cor1(:) /= n_d(:)) this%bnd(:, LO) = BND_MPI
      where (dom%periodic(:) .and. this%off(:)    /= 0)      this%bnd(:, HI) = BND_MPI

      this%ext_bnd(:, :) = .false.
      do i = xdim, zdim
         if (dom%has_dir(i) .and. .not. dom%periodic(i)) then
            this%ext_bnd(i, LO) = (this%off(i)    == 0)
            this%ext_bnd(i, HI) = (this%h_cor1(i) == n_d(i)) !! \warning not true on AMR
         endif
      enddo

      ! For shear boundaries and some domain decompositions it is possible that a boundary can be mixed 'per' with 'mpi'

      if (cdd%comm3d == MPI_COMM_NULL) then
         call inflate_req(size([LO, HI]) * 2 * ndims * nproc) ! 2 = count([i_bnd, o_bnd])
         if (any(dom%bnd(:, :) == BND_COR)) call die("[grid_container:init] Corner BC not implemented without comm3d")
#ifdef SHEAR_BND
         call die("[grid_container:init] SHEAR_BND not implemented without comm3d")
#endif /* SHEAR_BND */
      else
         call inflate_req(max(size([LO, HI]) * size([BLK, BND]) * ndims, int(nproc))) ! just another way of defining '4 * 3' ;-)
         ! write_plot_hdf5 requires nproc entries for the status array

         if (any(dom%bnd(xdim:ydim, :) == BND_COR) .and. (cdd%psize(xdim) /= cdd%psize(ydim) .or. n_d(xdim) /= n_d(ydim))) then
            write(msg, '(a,4(i4,a))')"[grid_container:init] Corner BC require psize(xdim) equal to psize(ydim) and n_d(xdim) equal to n_d(ydim). Detected: [", &
                 &                   cdd%psize(xdim),",",cdd%psize(ydim), "] and [",n_d(xdim),",",n_d(ydim),"]"
            call die(msg)
         endif
         if (any(dom%bnd(zdim, :) == BND_COR)) call die("[grid_container:init] Corner BC not allowed for z-direction")

#ifdef SHEAR_BND
         if (cdd%psize(ydim) > 1) call die("[domain:initmpi] Shear-pediodic boundary conditions do not permit psize(ydim) > 1")
         ! This is possible to be implemented with mpi_noncart

#ifndef FFTW
         this%bnd(xdim, :) = BND_MPI
         if (cdd%pcoords(xdim) == 0)             this%bnd(xdim, LO) = BND_SHE
         if (cdd%pcoords(xdim) == psize(xdim)-1) this%bnd(xdim, HI) = BND_SHE
#endif /* !FFTW */
#endif /* SHEAR_BND */
#ifdef DEBUG
         do i = xdim, zdim
            write(msg,*) 'dir',i,': ',cdd%procn(i, LO), proc, cdd%procn(i, HI)
            call printinfo(msg)
         enddo
#endif /* DEBUG */
      endif

      !> \todo allocate this conditionally, only when comm3d is in use
      if (allocated(this%mbc)) call die("[grid_container:init] this%mbc already allocated")
      allocate(this%mbc(FLUID:ARR, xdim:zdim, LO:HI, BND:BLK, 1:dom%nb))
      this%mbc(:, :, :, :, :) = INVALID

      do i = xdim, zdim
         if (dom%has_dir(i)) then
            if (this%n_b(i) < 1) call die("[grid_init] Too many CPUs for such a small grid.")
            if (this%n_b(i) < dom%nb) call warn("[grid_init] domain size in some directions is < nb, which may result in incomplete boundary cell update")
         endif
      enddo

      where (dom%has_dir(:))
         this%n_(:) = this%n_b(:) + I_TWO * dom%nb       ! Block total grid size with guardcells
         this%ijkse(:, LO) = dom%nb + I_ONE
         this%ijkse(:, HI) = dom%nb + this%n_b(:)
         this%dl(:) = dom%L_(:) / n_d(:)
         this%fbnd(:, LO) = dom%edge(:, LO) + this%dl(:) * this%off(:)
         this%fbnd(:, HI) = dom%edge(:, LO) + this%dl(:) * this%h_cor1(:)
      elsewhere
         this%n_(:) = 1
         this%ijkse(:, LO) = 1
         this%ijkse(:, HI) = 1
         this%dl(:) = 1.0
         this%fbnd(:, LO) = dom%edge(:, LO)
         this%fbnd(:, HI) = dom%edge(:, HI)
      endwhere

      if (dom%has_dir(xdim)) then
         this%isb   = I_TWO*dom%nb
         this%ieb   = this%n_b(xdim)+I_ONE
      else
         this%isb   = 1
         this%ieb   = 1
      endif

      if (dom%has_dir(ydim)) then
         this%jsb   = I_TWO*dom%nb
         this%jeb   = this%n_b(ydim)+I_ONE
      else
         this%jsb   = 1
         this%jeb   = 1
      endif

      if (dom%has_dir(zdim)) then
         this%ksb   = I_TWO*dom%nb
         this%keb   = this%n_b(zdim)+I_ONE
      else
         this%ksb   = 1
         this%keb   = 1
      endif

      this%vol = product(this%fbnd(:, HI)-this%fbnd(:, LO), mask=dom%has_dir(:))
      this%dvol = product(this%dl(:), mask=dom%has_dir(:))

      this%maxxyz = maxval(this%n_(:), mask=dom%has_dir(:))

      do i = xdim, zdim
         call this%set_axis(i)
      enddo

      this%dxmn = minval(this%dl(:), mask=dom%has_dir(:))

      ! some shortcuts for convenience
      this%idl(:) = 1./this%dl(:)

      this%dx = this%dl(xdim)
      this%dy = this%dl(ydim)
      this%dz = this%dl(zdim)

      this%idx = 1./this%dx
      this%idy = 1./this%dy
      this%idz = 1./this%dz

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

      this%dvol2 = this%dvol**2

      if (allocated(this%prolong_) .or. allocated(this%prolong_x) .or. allocated(this%prolong_xy) ) call die("[grid_container:init] prolong_* arrays already allocated")
      allocate(this%prolong_  (this%nxb/2+2*dom%nb, this%nyb/2+2*dom%nb, this%nzb/2+2*dom%nb), &
           &   this%prolong_x (this%n_(xdim),       this%nyb/2+2*dom%nb, this%nzb/2+2*dom%nb), &
           &   this%prolong_xy(this%n_(xdim),       this%n_(ydim),       this%nzb/2+2*dom%nb))

      this%prolong_  (:, :, :) = big_float
      this%prolong_x (:, :, :) = big_float
      this%prolong_xy(:, :, :) = big_float

   end subroutine init

!> \brief Calculate arrays of coordinates along a given direction

   subroutine set_axis(this, d)

      use constants,  only: LO, HI, half, one, zero, xdim, ydim, zdim
      use dataio_pub, only: die
      use domain,     only: dom

      implicit none

      class(grid_container), intent(inout) :: this !< grid container, where the arrays have to be set
      integer,               intent(in)    :: d    !< direction

      real, dimension(:), allocatable :: a0, al, ar, ia
      integer :: i

      allocate(a0(this%n_(d)), al(this%n_(d)), ar(this%n_(d)), ia(this%n_(d)))

      if (dom%has_dir(d)) then
         a0(:) = dom%edge(d, LO) + this%dl(d) * ([(i, i=1, this%n_(d))] - half - dom%nb + this%off(d))
      else
         a0(:) = half*(this%fbnd(d, LO) + this%fbnd(d, HI))
      endif

      al(:) = a0(:) - half*this%dl(d)
      ar(:) = a0(:) + half*this%dl(d)

      where ( a0(:) /= zero )
         ia(:) = one/a0(:)
      elsewhere
         ia(:) = zero
      endwhere

!--- Assignments -----------------------------------------------------------
         ! left zone boundaries:  xl, yl, zl
         ! zone centers:          x,  y,  z
         ! right zone boundaries: xr, yr, zr

      select case (d)
         case (xdim)
            if (allocated(this%x) .or. allocated(this%xl) .or. allocated(this%xr) .or. allocated(this%inv_x)) call die("[grid_container:set_axis] x-coordinates already allocated")
            call move_alloc(a0, this%x)
            call move_alloc(al, this%xl)
            call move_alloc(ar, this%xr)
            call move_alloc(ia, this%inv_x)
         case (ydim)
            if (allocated(this%y) .or. allocated(this%yl) .or. allocated(this%yr) .or. allocated(this%inv_y)) call die("[grid_container:set_axis] y-coordinates already allocated")
            call move_alloc(a0, this%y)
            call move_alloc(al, this%yl)
            call move_alloc(ar, this%yr)
            call move_alloc(ia, this%inv_y)
         case (zdim)
            if (allocated(this%z) .or. allocated(this%zl) .or. allocated(this%zr) .or. allocated(this%inv_z)) call die("[grid_container:set_axis] z-coordinates already allocated")
            call move_alloc(a0, this%z)
            call move_alloc(al, this%zl)
            call move_alloc(ar, this%zr)
            call move_alloc(ia, this%inv_z)
         case default
            call die("[grid_container:set_axis] invalid direction")
      end select

   end subroutine set_axis

!> \brief Routines that deallocates all internals of the grid container

   subroutine cleanup(this)

      use constants, only: FLUID, ARR, xdim, zdim, LO, HI, BND, BLK, INVALID
      use domain,    only: dom
      use mpisetup,  only: mpi_err

      implicit none

      class(grid_container), intent(inout) :: this
      integer :: d, t, g, b
      integer, parameter :: nseg = 2
      type(tgt_list), dimension(nseg) :: fc_tgt

      if (this%grid_id <= INVALID) return ! very dirty workaround for unability to determine whether a given cg was already deallocated

      if (allocated(this%x))     deallocate(this%x)
      if (allocated(this%xl))    deallocate(this%xl)
      if (allocated(this%xr))    deallocate(this%xr)
      if (allocated(this%inv_x)) deallocate(this%inv_x)
      if (allocated(this%y))     deallocate(this%y)
      if (allocated(this%yl))    deallocate(this%yl)
      if (allocated(this%yr))    deallocate(this%yr)
      if (allocated(this%inv_y)) deallocate(this%inv_y)
      if (allocated(this%z))     deallocate(this%z)
      if (allocated(this%zl))    deallocate(this%zl)
      if (allocated(this%zr))    deallocate(this%zr)
      if (allocated(this%inv_z)) deallocate(this%inv_z)

      if (allocated(this%gc_xdim)) deallocate(this%gc_xdim)
      if (allocated(this%gc_ydim)) deallocate(this%gc_ydim)
      if (allocated(this%gc_zdim)) deallocate(this%gc_zdim)

      if (allocated(this%mbc)) then
         do d = xdim, zdim
            if (dom%has_dir(d)) then
               do t = FLUID, ARR
                  do b = 1, dom%nb
                     if (this%mbc(t, d, LO, BLK, b) /= INVALID) call MPI_Type_free(this%mbc(t, d, LO, BLK, b), mpi_err)
                     if (this%mbc(t, d, LO, BND, b) /= INVALID) call MPI_Type_free(this%mbc(t, d, LO, BND, b), mpi_err)
                     if (this%mbc(t, d, HI, BLK, b) /= INVALID) call MPI_Type_free(this%mbc(t, d, HI, BLK, b), mpi_err)
                     if (this%mbc(t, d, HI, BND, b) /= INVALID) call MPI_Type_free(this%mbc(t, d, HI, BND, b), mpi_err)
                  enddo
               enddo
            endif
         enddo
         deallocate(this%mbc)
      endif

      if (allocated(this%i_bnd)) then
         do d = xdim, zdim
            do b = lbound(this%i_bnd, dim=2), ubound(this%i_bnd, dim=2)
               if (allocated(this%i_bnd(d, b)%seg)) deallocate(this%i_bnd(d, b)%seg)
            enddo
         enddo
         deallocate(this%i_bnd)
      endif
      if (allocated(this%o_bnd)) then
         do d = xdim, zdim
            do b = lbound(this%o_bnd, dim=2), ubound(this%o_bnd, dim=2)
               if (allocated(this%o_bnd(d, b)%seg)) deallocate(this%o_bnd(d, b)%seg)
            enddo
         enddo
         deallocate(this%o_bnd)
      endif

      fc_tgt(1:nseg) = [ this%f_tgt, this%c_tgt ]
      do b = 1, nseg
         if (allocated(fc_tgt(b)%seg)) then
            do g = lbound(fc_tgt(b)%seg, dim=1), ubound(fc_tgt(b)%seg, dim=1)
               if (allocated(fc_tgt(b)%seg(g)%buf)) deallocate(fc_tgt(b)%seg(g)%buf)
            enddo
            deallocate(fc_tgt(b)%seg)
         endif
      enddo

      if (allocated(this%q)) then
         do g = lbound(this%q(:), dim=1), ubound(this%q(:), dim=1)
            call this%q(g)%clean
         enddo
         deallocate(this%q)
      endif

      do d = xdim, zdim
         do b = 1, dom%nb
            if (allocated(this%q_i_mbc)) then
               if (allocated(this%q_i_mbc(d, b)%mbc)) then
                  do g = lbound(this%q_i_mbc(d, b)%mbc, dim=1), ubound(this%q_i_mbc(d, b)%mbc, dim=1)
                     if (this%q_i_mbc(d, b)%mbc(g) /= INVALID) call MPI_Type_free(this%q_i_mbc(d, b)%mbc(g), mpi_err)
                  enddo
                  deallocate(this%q_i_mbc(d, b)%mbc)
               endif
            endif
            if (allocated(this%q_o_mbc)) then
               if (allocated(this%q_o_mbc(d, b)%mbc)) then
                  do g = lbound(this%q_o_mbc(d, b)%mbc, dim=1), ubound(this%q_o_mbc(d, b)%mbc, dim=1)
                     if (this%q_o_mbc(d, b)%mbc(g) /= INVALID) call MPI_Type_free(this%q_o_mbc(d, b)%mbc(g), mpi_err)
                  enddo
                  deallocate(this%q_o_mbc(d, b)%mbc)
               endif
            endif
         enddo
      enddo
      deallocate(this%q_i_mbc, this%q_o_mbc)


      if (allocated(this%w)) then
         do g = lbound(this%w(:), dim=1), ubound(this%w(:), dim=1)
            call this%w(g)%clean
         enddo
         deallocate(this%w)
      endif

      ! arrays not handled through named_array feature
      if (allocated(this%prolong_xy)) deallocate(this%prolong_xy)
      if (allocated(this%prolong_x))  deallocate(this%prolong_x)
      if (allocated(this%prolong_))   deallocate(this%prolong_)
      if (allocated(this%mg%bnd_x))   deallocate(this%mg%bnd_x)
      if (allocated(this%mg%bnd_y))   deallocate(this%mg%bnd_y)
      if (allocated(this%mg%bnd_z))   deallocate(this%mg%bnd_z)

      this%grid_id = INVALID

   end subroutine cleanup

!> \brief Create an MPI type for exchanging a segment of data

   subroutine set_mpi_types(sizes, se, mbc)

      use constants,  only: xdim, zdim, ndims
      use dataio_pub, only: die
      use mpi,        only: MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION
      use mpisetup,   only: mpi_err

      implicit none

      integer(kind=4), dimension(:),            intent(in)  :: sizes !< dimensions of the array
      integer(kind=8), dimension(ndims, LO:HI), intent(in)  :: se    !< segment to communicate
      integer(kind=4),                          intent(out) :: mbc   !< MPI Boundary conditions Container

      integer(kind=4), dimension(:), allocatable :: subsizes, starts

      allocate(subsizes(size(sizes)), starts(size(sizes)))

      if (size(sizes) == 1+ndims) then
         subsizes(1) = sizes(1)
         starts(1) = 0
      else if (size(sizes) /= ndims) then
         call die("[grid_container:set_mpi_types] Only 3D and 4D arrays are supported")
      endif

      subsizes(size(sizes)-zdim+xdim:size(sizes)) = int(se(:, HI) - se(:, LO) + 1, kind=4)
      starts  (size(sizes)-zdim+xdim:size(sizes)) = int(se(:, LO) - 1, kind=4)
      call MPI_Type_create_subarray(size(sizes), sizes, subsizes, starts,  MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, mbc, mpi_err)
      call MPI_Type_commit(mbc, mpi_err)

      deallocate(subsizes, starts)

   end subroutine set_mpi_types

!>
!! \brief is_overlap_per checks if two given blocks placed within a periodic domain are overlapping.
!!
!! \details to handle shearing box which is divided in y-direction at the edges, one has to provide another subroutine (is_overlap_per_shear) and add it to interface is_overlap
!<

   logical function is_overlap_per(this, other, periods) result(share)

      use constants,  only: xdim, ydim, zdim, ndims, LO, HI
      use domain,     only: dom

      implicit none

      integer(kind=8), dimension(xdim:zdim, LO:HI), intent(in) :: this !< object invoking type-bound procedure
      integer(kind=8), dimension(xdim:zdim, LO:HI), intent(in) :: other !< the other box
      integer(kind=8), dimension(xdim:zdim), intent(in)        :: periods !< where >0 then the direction is periodic with the given number of cells

      integer :: i, j, k
      integer(kind=8), dimension(xdim:zdim, LO:HI) :: oth

      share = .false.
      do i = -1, 1
         if ((dom%has_dir(xdim) .or. periods(xdim)>0) .or. i==0) then
            do j = -1, 1
               if ((dom%has_dir(ydim) .or. periods(ydim)>0) .or. j==0) then
                  do k = -1, 1
                     if ((dom%has_dir(zdim) .or. periods(zdim)>0) .or. k==0) then
                        oth(:,:) = other(:,:) + reshape([i*periods(xdim), j*periods(ydim), k*periods(zdim), i*periods(xdim), j*periods(ydim), k*periods(zdim)], [ndims, HI])
                        share = share .or. is_overlap_simple(this, oth)
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

      integer(kind=8), dimension(xdim:zdim, LO:HI), intent(in) :: this !< object invoking type-bound procedure
      integer(kind=8), dimension(xdim:zdim, LO:HI), intent(in) :: other !< the other box

      integer :: d

      share = .true.
      do d = xdim, zdim
         if (dom%has_dir(d)) share = share .and. (other(d, LO) <= this(d, HI)) .and. (other(d, HI) >= this(d, LO))
      enddo

   end function is_overlap_simple

!> \brief Initialize the communicators for q even if there are no q arrays at the moment

   subroutine set_q_mbc(this)

      use constants, only: ndims, xdim, zdim
      use domain,    only: dom

      implicit none

      class(grid_container), intent(inout) :: this
      integer :: d, ib, g

      allocate(this%q_i_mbc(ndims, dom%nb), this%q_o_mbc(ndims, dom%nb))
      do d = xdim, zdim
         do ib = 1, dom%nb
            if (allocated(this%i_bnd(d, ib)%seg)) then
               allocate(this%q_i_mbc(d, ib)%mbc(lbound(this%i_bnd(d, ib)%seg, dim=1):ubound(this%i_bnd(d, ib)%seg, dim=1)), &
                    &   this%q_o_mbc(d, ib)%mbc(lbound(this%i_bnd(d, ib)%seg, dim=1):ubound(this%i_bnd(d, ib)%seg, dim=1)))
               do g = lbound(this%i_bnd(d, ib)%seg, dim=1), ubound(this%i_bnd(d, ib)%seg, dim=1)
                  call set_mpi_types(this%n_(:), this%i_bnd(d, ib)%seg(g)%se(:,:), this%q_i_mbc(d, ib)%mbc(g))
                  call set_mpi_types(this%n_(:), this%o_bnd(d, ib)%seg(g)%se(:,:), this%q_o_mbc(d, ib)%mbc(g))
               enddo
            endif
         enddo
      enddo

   end subroutine set_q_mbc

end module grid_cont
