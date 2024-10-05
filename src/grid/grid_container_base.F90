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

!> \brief Module containing the basics of grid container as an abstract type and its associated methods

module grid_cont_base

   use constants,        only: xdim, zdim, ndims, LO, HI, CENTER, INV_CENTER
   use level_essentials, only: level_t
   use real_vector,      only: real_vec_t

   implicit none

   private
   public :: grid_container_base_t

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

   !> \brief Everything required for autonomous computation of a single sweep on a portion of the domain on a single process
   type, abstract :: grid_container_base_t

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
      real :: suminv                                             !< idx + idy + idz

      ! Grid properties

      ! cell count and position
      integer(kind=4), dimension(ndims) :: n_b                   !< [nxb, nyb, nzb]
      integer(kind=4) :: nxb                                     !< number of %grid cells in one block (without boundary cells) in x-direction
      integer(kind=4) :: nyb                                     !< number of %grid cells in one block (without boundary cells) in y-direction
      integer(kind=4) :: nzb                                     !< number of %grid cells in one block (without boundary cells) in z-direction
      class(level_t), pointer :: l                               !< level essential data

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

      ! External boundary conditions and internal boundaries

      integer(kind=4), dimension(ndims, LO:HI)           :: bnd         !< type of boundary conditions coded in integers
      logical,         dimension(xdim:zdim, LO:HI)       :: ext_bnd     !< .false. for BND_PER and BND_MPI

      ! Physical size and coordinates

      real, dimension(ndims, LO:HI) :: fbnd                      !< current block boundary positions

      type(real_vec_t), dimension(CENTER:INV_CENTER, ndims) :: coord !< all coordinates (CENTER, LEFT, RIGHT, INV_CENTER)
      ! shortcuts
      real, pointer, dimension(:) :: x                             !< array of x-positions of %grid cells centers
      real, pointer, dimension(:) :: y                             !< array of x-positions of %grid cells centers
      real, pointer, dimension(:) :: z                             !< array of x-positions of %grid cells centers
      real, pointer, dimension(:) :: inv_x                         !< array of invert x-positions of %grid cells centers
      real, pointer, dimension(:) :: inv_y                         !< array of invert y-positions of %grid cells centers
      real, pointer, dimension(:) :: inv_z                         !< array of invert z-positions of %grid cells centers

      ! Non-cartesian geometrical factors

      real, allocatable, dimension(:,:,:) :: gc_xdim             !< array of geometrical coefficients in x-direction
      real, allocatable, dimension(:,:,:) :: gc_ydim             !< array of geometrical coefficients in y-direction
      real, allocatable, dimension(:,:,:) :: gc_zdim             !< array of geometrical coefficients in z-direction

      ! Misc
      type(mg_arr), pointer :: mg                                !< multigrid arrays
      real :: vol                                                !< volume of the grid; BEWARE: for cylindrical geometry it needs to be integrated over x(:) to get real volume
      real :: dxmn                                               !< the smallest length of the %grid cell (among dx, dy, and dz)
      real :: dxmn2                                              !< squared dxmn
      integer(kind=4) :: maxxyz                                  !< maximum number of %grid cells in any direction
      integer :: grid_id                                         !< index of own segment in own level decomposition, e.g. my_se(:,:) = base%level%dot%gse(proc)%c(grid_id)%se(:,:)
      logical :: has_previous_timestep                           !< used to prevent timestep retries on freshly created blocks

   contains

      procedure          :: init_gc_base                         !< Initialization
      procedure          :: cleanup_base                         !< Deallocate all internals
      procedure, private :: set_coords                           !< Calculate arrays of coordinates along a given direction

   end type grid_container_base_t

contains

!> \brief This method sets up the grid container basic variables and coordinates.

   subroutine init_gc_base(this, my_se, grid_id, l)

      use constants,        only: xdim, ydim, zdim, LO, HI, I_ONE, I_TWO, BND_MPI, BND_COR, GEO_XYZ, GEO_RPZ, dpi
      use dataio_pub,       only: die, warn
      use domain,           only: dom
      use level_essentials, only: level_t

      implicit none

      class(grid_container_base_t), target, intent(inout) :: this     !< object invoking type-bound procedure (grid_container)
      integer(kind=8), dimension(:,:),      intent(in)    :: my_se    !< my segment
      integer,                              intent(in)    :: grid_id  !< ID which should be unique across level
      class(level_t), pointer,              intent(in)    :: l        !< level essential data

      integer :: i

      this%l          => l
      this%grid_id    = grid_id
      this%my_se(:,:) = my_se(:, :)
      this%n_b(:)     = int(this%my_se(:, HI) - this%my_se(:, LO) + I_ONE, 4) ! Block 'physical' grid sizes

      if (any(this%n_b(:) <= 0)) call die("[grid_container_base:init_gc_base] Mixed positive and non-positive grid sizes")

      ! Inherit the boundaries from the domain, then set MPI or SHEAR boundaries where applicable
      this%bnd(:,:) = dom%bnd(:,:)
      where (my_se(:, LO)         /= l%off(:)           ) this%bnd(:, LO) = BND_MPI
      where (my_se(:, HI) + I_ONE /= l%off(:) + l%n_d(:)) this%bnd(:, HI) = BND_MPI
      ! For periodic boundaries do not set BND_MPI when local domain spans through the whole computational domain in given direction.
      where (dom%periodic(:) .and. this%my_se(:, HI) + I_ONE /= l%n_d(:)) this%bnd(:, LO) = BND_MPI
      where (dom%periodic(:) .and. this%my_se(:, LO)         /= 0       ) this%bnd(:, HI) = BND_MPI

      this%ext_bnd(:, :) = .false.  ! slightly modified variant of l%has_ext_bnd would do the work
      do i = xdim, zdim
         if (dom%has_dir(i) .and. .not. dom%periodic(i)) then
            this%ext_bnd(i, LO) = (my_se(i, LO)         == l%off(i))
            this%ext_bnd(i, HI) = (my_se(i, HI) + I_ONE == l%off(i) + l%n_d(i))
         endif
      enddo

      ! For shear boundaries and some domain decompositions it is possible that a boundary can be mixed 'per' with 'mpi'

      if (any(dom%bnd(xdim:ydim, :) == BND_COR)) call die("[grid_container_base:init_gc_base] BND_COR unimplemented")
      if (any(dom%bnd(zdim, :) == BND_COR)) call die("[grid_container_base:init_gc_base] Corner BC not allowed for z-direction")

#ifdef SHEAR_BND
      call die("[grid_container_base:init_gc_base] Shear-pediodic boundary conditions unimplemented")
      ! This is possible to be implemented
#endif /* SHEAR_BND */

      do i = xdim, zdim
         if (dom%has_dir(i)) then
            if (this%n_b(i) < 1) call die("[grid_container_base:init_gc_base] Too many CPUs for such a small grid.")
            if (this%n_b(i) < dom%nb) call warn("[grid_container_base:init_gc_base] domain size in some directions is < nb, which may result in incomplete boundary cell update")
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

      if (any(this%dl < tiny(1.0))) call die("[grid_container_base:init_gc_base] found cell size < tiny()")

      this%isb = this%ijkseb(xdim, LO)
      this%ieb = this%ijkseb(xdim, HI)
      this%jsb = this%ijkseb(ydim, LO)
      this%jeb = this%ijkseb(ydim, HI)
      this%ksb = this%ijkseb(zdim, LO)
      this%keb = this%ijkseb(zdim, HI)

      ! Beware: 0-dimensional simulations should not rely on cg%vol or cg%dvol
      select case (dom%geometry_type)
         case (GEO_XYZ)
            this%vol = merge(0., product(this%fbnd(:, HI)-this%fbnd(:, LO), mask=dom%has_dir(:)), dom%eff_dim == 0)
            this%dvol = merge(0., product(this%dl(:), mask=dom%has_dir(:)), dom%eff_dim == 0)
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

      if (dom%eff_dim == 0) then
         this%maxxyz = 1
      else
         this%maxxyz = maxval(this%n_(:), mask=dom%has_dir(:))
      endif

      call this%set_coords

      if (dom%eff_dim == 0) then
         this%dxmn  = minval(dom%L_(:))
      else
         this%dxmn  = minval(this%dl(:), mask=dom%has_dir(:))
      endif
      this%dxmn2 = (this%dxmn)**2

      ! some shortcuts for convenience
      this%idl(:) = 1./this%dl(:)
      this%suminv = merge(sum(this%idl(:), mask=dom%has_dir(:)), 0., dom%eff_dim /= 0)

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

   end subroutine init_gc_base

!> \brief Calculate arrays of coordinates along a given direction

   subroutine set_coords(this)

      use constants,  only: LO, HI, half, one, zero, xdim, ydim, zdim, CENTER, LEFT, RIGHT, INV_CENTER
      use dataio_pub, only: die, warn
      use domain,     only: dom
      use func,       only: operator(.notequals.)

      implicit none

      class(grid_container_base_t), intent(inout) :: this  !< grid container, where the arrays have to be set

      integer :: d, i
      integer, parameter :: safety_warn_factor = 1000 ! warn if a cell size is smaller than this * epsilon(coordinates)

      do d = xdim, zdim
         do i = CENTER, INV_CENTER
            if (this%coord(i, d)%associated()) call die("[grid_container_base:set_coords] a coordinate already allocated")
            call this%coord(i, d)%allocate(this%lhn(d, LO), this%lhn(d, HI))
         enddo

         if (dom%has_dir(d)) then
            this%coord(CENTER, d)%r(:) = this%fbnd(d, LO) + this%dl(d) * ([(i, i=this%lhn(d, LO), this%lhn(d, HI))] - this%ijkse(d, LO) + half)
         else
            this%coord(CENTER, d)%r(:) = half*(this%fbnd(d, LO) + this%fbnd(d, HI))
         endif

         this%coord(LEFT,  d)%r(:) = this%coord(CENTER, d)%r(:) - half*this%dl(d)
         this%coord(RIGHT, d)%r(:) = this%coord(CENTER, d)%r(:) + half*this%dl(d)

         where (this%coord(CENTER, d)%r(:) .notequals. zero)
            this%coord(INV_CENTER, d)%r(:) = one/this%coord(CENTER, d)%r(:)
         elsewhere
            this%coord(INV_CENTER, d)%r(:) = zero
         endwhere

         ! Generally nobody should subtract one cell coordinate from another in code solvers. One should use cell sizes instead.
         ! The problem may arise when initial conditions are comparing coordinates to set something on the left or right side of some line.
         ! When the cell size is too small compared to the coordinates, such line cannot be properly calculated
         ! Note that since we force real kind=8, we can use a named constant instead of epsilon
         if (dom%has_dir(d)) then
            if ( any(abs(this%coord(CENTER, d)%r(:) - this%coord(LEFT,  d)%r(:)) < tiny(1.0)) .or. &
                 any(abs(this%coord(CENTER, d)%r(:) - this%coord(RIGHT, d)%r(:)) < tiny(1.0)) ) call die("[grid_container_base:set_coords] cannot distinguish between center and face coordinates of a cell")
            if ( any(abs(this%coord(CENTER, d)%r(:)-this%coord(LEFT,  d)%r(:)) < safety_warn_factor*epsilon(this%coord(CENTER, d)%r(:))*this%coord(CENTER, d)%r(:)) .or. &
                 any(abs(this%coord(CENTER, d)%r(:)-this%coord(RIGHT, d)%r(:)) < safety_warn_factor*epsilon(this%coord(CENTER, d)%r(:))*this%coord(CENTER, d)%r(:))) &
                 call warn("[grid_container_base:set_coords] cell sizes are much smaller than coordinates. Inaccuracies in setting the initial conditions may happen.")
         endif
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

   subroutine cleanup_base(this)

      use constants, only: xdim, zdim, CENTER, INV_CENTER

      implicit none

      class(grid_container_base_t), intent(inout) :: this  !< object invoking type-bound procedure

      integer :: b, cdim

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

   end subroutine cleanup_base

end module grid_cont_base
