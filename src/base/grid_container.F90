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
!! \brief Module containing the grid container type and its associated methods
!<
module grid_cont

   use constants, only: xdim, zdim, ndims, LO, HI, BND, BLK, FLUID, ARR
   use types,     only: axes, bnd_list, array4d, array3d

   implicit none

   private
   public :: grid_container, cg_list_element, cg_list, cg_set

   type, extends(axes) :: grid_container
      real    :: dx                             !< length of the %grid cell in x-direction
      real    :: dy                             !< length of the %grid cell in y-direction
      real    :: dz                             !< length of the %grid cell in z-direction
      real    :: idx                            !< inverted length of the %grid cell in x-direction
      real    :: idy                            !< inverted length of the %grid cell in y-direction
      real    :: idz                            !< inverted length of the %grid cell in z-direction
      real    :: dxmn                           !< the smallest length of the %grid cell (among dx, dy, and dz)
      real    :: dvol                           !< volume of one %grid cell
      real    :: vol                            !< volume of the grid; BEWARE: for cylindrical geometry it need to be multiplied by appropriate x(:) to get real volume
      real, dimension(ndims, LO:HI) :: fbnd     !< current block boundary positions

      real, dimension(ndims)          :: dl     !< array of %grid cell sizes in all directions
      real, dimension(ndims)          :: idl    !< array of inverted %grid cell sizes in all directions

      real, allocatable, dimension(:) :: inv_x  !< array of invert x-positions of %grid cells centers
      real, allocatable, dimension(:) :: inv_y  !< array of invert y-positions of %grid cells centers
      real, allocatable, dimension(:) :: inv_z  !< array of invert z-positions of %grid cells centers
      real, allocatable, dimension(:) :: xl     !< array of x-positions of %grid cells left borders
      real, allocatable, dimension(:) :: yl     !< array of y-positions of %grid cells left borders
      real, allocatable, dimension(:) :: zl     !< array of z-positions of %grid cells left borders
      real, allocatable, dimension(:) :: xr     !< array of x-positions of %grid cells right borders
      real, allocatable, dimension(:) :: yr     !< array of y-positions of %grid cells right borders
      real, allocatable, dimension(:) :: zr     !< array of z-positions of %grid cells right borders

      real, allocatable, dimension(:) :: dprof  !< Array used for storing density during calculation of hydrostatic equilibrium

      integer(kind=4) :: nx                             !< number of %grid cells in one block in x-direction
      integer(kind=4) :: ny                             !< number of %grid cells in one block in y-direction
      integer(kind=4) :: nz                             !< number of %grid cells in one block in z-direction
      integer(kind=4) :: nxb                            !< number of %grid cells in one block (without boundary cells) in x-direction
      integer(kind=4) :: nyb                            !< number of %grid cells in one block (without boundary cells) in y-direction
      integer(kind=4) :: nzb                            !< number of %grid cells in one block (without boundary cells) in z-direction
      integer(kind=4) :: is                             !< index of the first %grid cell of physical domain in x-direction
      integer(kind=4) :: ie                             !< index of the last %grid cell of physical domain in x-direction
      integer(kind=4) :: js                             !< index of the first %grid cell of physical domain in y-direction
      integer(kind=4) :: je                             !< index of the last %grid cell of physical domain in y-direction
      integer(kind=4) :: ks                             !< index of the first %grid cell of physical domain in z-direction
      integer(kind=4) :: ke                             !< index of the last %grid cell of physical domain in z-direction
      integer(kind=4) :: nb                             !< number of boundary cells surrounding the physical domain, same for all directions
      integer(kind=4) :: maxxyz                         !< maximum number of %grid cells in any direction
      integer(kind=4) :: isb, ieb, jsb, jeb, ksb, keb   !< auxiliary indices for exchanging boundary data, (e.g. is:isb -> ie+1:nx, ieb:ie -> 1:nb)

      logical :: empty                          !< .true. if there are no cells to process (e.g. some processes at base level in multigrid gravity)

      integer(kind=8), dimension(ndims) :: off  !< offset of the local domain within computational domain
      integer(kind=4), dimension(ndims) :: n_b  !< [nxb, nyb, nzb]
      integer(kind=4), dimension(ndims) :: n_   !< [nx,  ny,  nz ]

      integer(kind=4), dimension(ndims, LO:HI)  :: ijkse !< [[is, js, ks], [ie, je, ke]]
      integer, dimension(ndims, LO:HI)  :: bnd  !< type of boundary conditions coded in integers

      integer(kind=4), dimension(FLUID:ARR, xdim:zdim, LO:HI, BND:BLK) :: mbc !< MPI Boundary conditions Container

      !>
      !! description of incoming and outgoing boundary data,
      !! the shape is (xdim:zdim, FLUID:ARR) for cg (base grid container)) and (xdim:zdim, mg_nb) for plvl (multigrid level container)
      !<
      type(bnd_list), dimension(:, :), allocatable :: i_bnd, o_bnd

      real, allocatable, dimension(:,:,:) :: gc_xdim !< array of geometrical coefficients in x-direction
      real, allocatable, dimension(:,:,:) :: gc_ydim !< array of geometrical coefficients in y-direction
      real, allocatable, dimension(:,:,:) :: gc_zdim !< array of geometrical coefficients in z-direction

      type(array3d) :: cs_iso2
      type(array3d) :: wa                       !< Temporary array used for different purposes, usually has dimension (grid::nx, grid::ny, grid::nz)
      type(array3d) :: gpot                     !< Array for sum of gravitational potential at t += dt
      type(array3d) :: hgpot                    !< Array for sum of gravitational potential at t += 0.5*dt
      type(array3d) :: gp                       !< Array for gravitational potential from external fields
      type(array3d) :: sgp                      !< Array for gravitational potential from multigrid or FFT solver
      type(array3d) :: sgpm                     !< Array for gravitational potential from multigrid or FFT solver at previous timestep saved by source_terms_grav.

      type(array4d) :: u                        !< Main array of all fluids' components
      type(array4d) :: uh                       !< Main array of all fluids' components (for t += dt/2)
      type(array4d) :: u0                       !< Copy of main array of all fluids' components
      type(array4d) :: b0                       !< Copy of main array of magnetic field's components
      type(array4d) :: b                        !< Main array of magnetic field's components

   contains

      procedure :: init
      procedure :: cleanup
      procedure :: internal_boundaries

   end type grid_container

   ! the prv and nxt pointers are not elements of the grid_container type to allow membership in several lists simultaneously
   type cg_list_element
      type(grid_container), pointer :: cg        !< the current grid container
      type(cg_list_element), pointer :: prv, nxt !< pointers to previous and next grid container or null() at the end of the list
   end type cg_list_element

   type cg_list
      type(cg_list_element), dimension(:), pointer :: cg_l
   end type cg_list

   ! On an uniform, nowhere refined grid, cg_levels will have only one element, and all cg_all elements ill be members of cg_leafs, cg_base and cg_levels(1) lists
   type cg_set
      type(grid_container), dimension(:), allocatable :: cg_all    !< all grid containers
      type(cg_list), dimension(:), allocatable        :: cg_levels !< grid containers grouped by level
      type(cg_list)                                   :: cg_leafs  !< grid containers that are not fully covered by finer grids
      type(cg_list)                                   :: cg_base   !< grid containers on the base level

      contains
         procedure :: get_root
   end type cg_set

contains
!-----------------------------------------------------------------------------
!>
!! \brief sets pointer to grid root
!<
   subroutine get_root(this, cgp)
      implicit none
      class(cg_set), intent(in)                   :: this
      type(cg_list_element), pointer, intent(out) :: cgp

      cgp => this%cg_leafs%cg_l(1)
      return
   end subroutine get_root
!-----------------------------------------------------------------------------
   subroutine init(this, dom)

      use constants,  only: PIERNIK_INIT_DOMAIN, xdim, ydim, zdim, INVALID, INT4
      use dataio_pub, only: die, warn, code_progress
      use domain,     only: has_dir, translate_bnds_to_ints_dom, domain_container
      use mpisetup,   only: proc

      implicit none

      class(grid_container), intent(inout) :: this ! intent(out) would silently clear everything, that was already set (also the fields in types derived from grid_container)
      type(domain_container), intent(in) :: dom

      integer :: i

      if (code_progress < PIERNIK_INIT_DOMAIN) call die("[grid:init] MPI not initialized.")
      if (ubound(dom%pse(proc)%sel(:,:,:), dim=1) > 1) call die("[grid_container:init] Multiple blocks per process not implemented yet")

      this%nb = dom%nb
      this%dxmn = huge(1.0)

      this%off(:) = dom%pse(proc)%sel(1, :, LO)  ! Block offset on the dom% should be between 0 and nxd-nxb
      this%n_b(:) = int(dom%pse(proc)%sel(1, :, HI) - dom%pse(proc)%sel(1, :, LO) + 1, 4) ! Block 'physical' grid sizes

      if (all(this%n_b(:) == 0)) then
         this%empty = .true.
      else if (any(this%n_b(:) == 0)) then
         call die("[grid_init] Mixed positive and non-positive grid sizes")
      else
         this%empty = .false.
      endif

      this%nxb = this%n_b(xdim)
      this%nyb = this%n_b(ydim)
      this%nzb = this%n_b(zdim)

      this%bnd(:,:) = translate_bnds_to_ints_dom()
      this%mbc(:, :, :, :) = INVALID

      if (this%empty) then
         this%fbnd(:,:) = dom%edge(:,:)

         this%nx    = 0
         this%is    = this%nb + 1_INT4
         this%ie    = this%nb
         this%isb   = 0 ! ???
         this%ieb   = 0 ! ???
         this%dx    = 1.0

         this%ny    = 0
         this%js    = this%nb + 1_INT4
         this%je    = this%nb
         this%jsb   = 0
         this%jeb   = 0
         this%dy    = 1.0

         this%nz    = 0
         this%ks    = this%nb + 1_INT4
         this%ke    = this%nb
         this%ksb   = 0
         this%keb   = 0
         this%dz    = 1.0

         this%idx = 1./this%dx
         this%idy = 1./this%dy
         this%idz = 1./this%dz

         this%dl(xdim:zdim) = [ this%dx, this%dy, this%dz ]
         this%idl(:) = 1./this%dl(:)

         this%vol = 0.
         this%dvol = 0.
         this%maxxyz = 0

      else

         do i = xdim, zdim
            if (has_dir(i)) then
               if (this%n_b(i) < 1) call die("[grid_init] Too many CPUs for a small grid.")
               if (this%n_b(i) < this%nb) call warn("[grid_init] domain size in some directions is < nb, which may result in incomplete boundary cell update")
            endif
         enddo

         if (has_dir(xdim)) then
            this%nx    = this%nxb + 2_INT4 * this%nb       ! Block total grid sizes
            this%is    = this%nb + 1_INT4
            this%ie    = this%nb + this%nxb
            this%isb   = 2_INT4*this%nb
            this%ieb   = this%nxb+1_INT4
            this%dx    = dom%L_(xdim) / dom%n_d(xdim)
            this%dxmn  = min(this%dxmn, this%dx)
            this%fbnd(xdim, LO) = dom%edge(xdim, LO) + this%dx *  this%off(xdim)
            this%fbnd(xdim, HI) = dom%edge(xdim, LO) + this%dx * (this%off(xdim) + this%nxb)
         else
            this%nx    = 1
            this%is    = 1
            this%ie    = 1
            this%isb   = 1
            this%ieb   = 1
            this%dx    = 1.0
            this%fbnd(xdim, :) = dom%edge(xdim, :)
         endif

         if (has_dir(ydim)) then
            this%ny    = this%nyb + 2_INT4 * this%nb
            this%js    = this%nb + 1_INT4
            this%je    = this%nb + this%nyb
            this%jsb   = 2_INT4*this%nb
            this%jeb   = this%nyb+1_INT4
            this%dy    = dom%L_(ydim) / dom%n_d(ydim)
            this%dxmn  = min(this%dxmn, this%dy)
            this%fbnd(ydim, LO) = dom%edge(ydim, LO) + this%dy *  this%off(ydim)
            this%fbnd(ydim, HI) = dom%edge(ydim, LO) + this%dy * (this%off(ydim) + this%nyb)
         else
            this%ny    = 1
            this%js    = 1
            this%je    = 1
            this%jsb   = 1
            this%jeb   = 1
            this%dy    = 1.0
            this%fbnd(ydim, :) = dom%edge(ydim, :)
         endif

         if (has_dir(zdim)) then
            this%nz    = this%nzb + 2_INT4 * this%nb
            this%ks    = this%nb + 1_INT4
            this%ke    = this%nb + this%nzb
            this%ksb   = 2_INT4*this%nb
            this%keb   = this%nzb+1_INT4
            this%dz    = dom%L_(zdim) / dom%n_d(zdim)
            this%dxmn  = min(this%dxmn, this%dz)
            this%fbnd(zdim, LO) = dom%edge(zdim, LO) + this%dz *  this%off(zdim)
            this%fbnd(zdim, HI) = dom%edge(zdim, LO) + this%dz * (this%off(zdim) + this%nzb)
         else
            this%nz    = 1
            this%ks    = 1
            this%ke    = 1
            this%ksb   = 1
            this%keb   = 1
            this%dz    = 1.0
            this%fbnd(zdim, :) = dom%edge(zdim, :)
         endif

         this%n_(xdim:zdim) = [ this%nx, this%ny, this%nz ]

         this%vol = product(this%fbnd(:, HI)-this%fbnd(:, LO), mask=has_dir(:))

         this%idx = 1./this%dx
         this%idy = 1./this%dy
         this%idz = 1./this%dz

         this%dl(xdim:zdim) = [ this%dx, this%dy, this%dz ]
         this%idl(:) = 1./this%dl(:)

         this%dvol = product(this%dl(:))

         allocate(this%x(this%nx), this%xl(this%nx), this%xr(this%nx), this%inv_x(this%nx))
         allocate(this%y(this%ny), this%yl(this%ny), this%yr(this%ny), this%inv_y(this%ny))
         allocate(this%z(this%nz), this%zl(this%nz), this%zr(this%nz), this%inv_z(this%nz))
         this%maxxyz = int(maxval([size(this%x), size(this%y), size(this%z)]), kind=4)

!--- Assignments -----------------------------------------------------------
         ! left zone boundaries:  xl, yl, zl
         ! zone centers:          x,  y,  z
         ! right zone boundaries: xr, yr, zr

!--- x-grids --------------------------------------------------------------

         if (has_dir(xdim)) then
            this%x(:) = dom%edge(xdim, LO) + this%dx * ([(i, i=1, this%nx)] - 0.5 - this%nb + this%off(xdim))
         else
            this%x(:) = 0.5*(this%fbnd(xdim, LO) + this%fbnd(xdim, HI))
         endif
         this%xl(:) = this%x(:) - 0.5*this%dx
         this%xr(:) = this%x(:) + 0.5*this%dx

         where ( this%x /= 0.0 )
            this%inv_x(:) = 1./this%x(:)
         elsewhere
            this%inv_x(:) = 0.
         endwhere

!--- y-grids --------------------------------------------------------------

         if (has_dir(ydim)) then
            this%y(:) = dom%edge(ydim, LO) + this%dy * ([(i, i=1, this%ny)] - 0.5 - this%nb + this%off(ydim))
         else
            this%y(:) = 0.5*(this%fbnd(ydim, LO) + this%fbnd(ydim, HI))
         endif
         this%yl(:) = this%y(:) - 0.5*this%dy
         this%yr(:) = this%y(:) + 0.5*this%dy

         where ( this%y /= 0.0 )
            this%inv_y(:) = 1./this%y(:)
         elsewhere
            this%inv_y(:) = 0.
         endwhere

!--- z-grids --------------------------------------------------------------

         if (has_dir(zdim)) then
            this%z(:) = dom%edge(zdim, LO) + this%dz * ([(i, i=1, this%nz)] - 0.5 - this%nb + this%off(zdim))
         else
            this%z(:) = 0.5*(this%fbnd(zdim, LO) + this%fbnd(zdim, HI))
         endif
         this%zl(:) = this%z(:) - 0.5*this%dz
         this%zr(:) = this%z(:) + 0.5*this%dz

         where ( this%z /= 0.0 )
            this%inv_z(:) = 1./this%z(:)
         elsewhere
            this%inv_z(:) = 0.
         endwhere

      endif

      this%ijkse(:, LO) = [ this%is, this%js, this%ks ]
      this%ijkse(:, HI) = [ this%ie, this%je, this%ke ]

#ifdef ISO
      call this%cs_iso2%init(this%n_(:))
#endif /* ISO */
#ifdef GRAV
      call this%gpot%init(this%n_(:))
      call this%hgpot%init(this%n_(:))
      call this%gp%init(this%n_(:))
#ifdef SELF_GRAV
      call this%sgp%init(this%n_(:))
      call this%sgpm%init(this%n_(:))
#endif /* SELF_GRAV */
#endif /* GRAV */

   end subroutine init

!-----------------------------------------------------------------------------
!
! This routine exchanges guardcells for BND_MPI and BND_PER boundaries.in u(:,:,:,:), b(:,:,:,:) and rank-3 arrays passed through argument list
! (preferably only pointers to the actual arrays are passed).
! The corners should be properly updated if this%[io]_bnd(:, ind) was set up appropriately (MPI_Waitall is called separately for each dimension).
!

   subroutine internal_boundaries(this, ind, pa3d, pa4d)

      use constants,  only: FLUID, MAG, CR, ARR, LO, HI, xdim, ydim, zdim
      use dataio_pub, only: die
      use domain,     only: has_dir
      use mpisetup,   only: comm, ierr, proc, req, status

      implicit none

      class(grid_container) :: this
      integer(kind=4), intent(in) :: ind   !< second index in [io]_bnd arrays
      real, optional, pointer, dimension(:,:,:)   :: pa3d
      real, optional, pointer, dimension(:,:,:,:) :: pa4d

      integer :: g, tag, d, nr
      integer(kind=8), dimension(xdim:zdim, LO:HI) :: ise, ose

!BEWARE: MPI_Waitall should be called after all grid containers post Isends and Irecvs
! This routine thus cannot be a metod of the grid_container type

      if (ind < minval([FLUID, MAG, CR, ARR]) .or. ind > maxval([FLUID, MAG, CR, ARR])) call die("[grid:internal_boundaries] wrong index")

      select case (ind)
         case (FLUID, MAG, CR)
            if (.not. present(pa4d)) call die("[grid:internal_boundaries] pa4d not provided")
            if (.not. associated(pa4d)) call die("[grid:internal_boundaries] pa4d == null()")
         case (ARR)
            if (.not. present(pa3d)) call die("[grid:internal_boundaries] pa3d not provided")
            if (.not. associated(pa3d)) call die("[grid:internal_boundaries] pa3d == null()")
         case default
            call die("[grid:internal_boundaries] not implemented yet")
            return
      end select
      if (present(pa3d) .and. present(pa4d)) call die("[grid:internal_boundaries] Both pa3d and pa4d are present")

      do d = xdim, zdim
         nr = 0
         if (has_dir(d)) then

            if (allocated(this%i_bnd(d, ind)%seg)) then
               do g = 1, ubound(this%i_bnd(d, ind)%seg(:), dim=1)
                  if (proc == this%i_bnd(d, ind)%seg(g)%proc) then
                     ise = this%i_bnd(d, ind)%seg(g)%se
                     ose(:,:) = ise(:,:)
                     if (ise(d, LO) < this%n_b(d)) then
                        ose(d, :) = ise(d, :) + this%n_b(d)
                     else
                        ose(d, :) = ise(d, :) - this%n_b(d)
                     endif
                     ! boundaries are always paired
                     if (ind == ARR) then
                        pa3d     (ise(xdim, LO):ise(xdim,HI), ise(ydim, LO):ise(ydim, HI), ise(zdim, LO):ise(zdim, HI)) = &
                             pa3d(ose(xdim, LO):ose(xdim,HI), ose(ydim, LO):ose(ydim, HI), ose(zdim, LO):ose(zdim, HI))
                     else
                        pa4d     (:, ise(xdim, LO):ise(xdim,HI), ise(ydim, LO):ise(ydim, HI), ise(zdim, LO):ise(zdim, HI)) = &
                             pa4d(:, ose(xdim, LO):ose(xdim,HI), ose(ydim, LO):ose(ydim, HI), ose(zdim, LO):ose(zdim, HI))
                     endif
                  else
                     ! BEWARE: Here we assume, that we have at most one chunk to communicate with a given process on a single side od the domain.
                     ! This will not be true when we allow many blocks per process and tag will need to be modified to include g or seg(g)%lh should become seg(g)%tag
                     tag = this%i_bnd(d, ind)%seg(g)%lh + HI*d
                     nr = nr + 1
                     if (ind == ARR) then
                        call MPI_Irecv(pa3d(1, 1, 1), 1, this%i_bnd(d, ind)%seg(g)%mbc, this%i_bnd(d, ind)%seg(g)%proc, tag, comm, req(nr), ierr)
                     else
                        call MPI_Irecv(pa4d(1, 1, 1, 1), 1, this%i_bnd(d, ind)%seg(g)%mbc, this%i_bnd(d, ind)%seg(g)%proc, tag, comm, req(nr), ierr)
                     endif
                  endif
               enddo
            endif
            if (allocated(this%o_bnd(d, ind)%seg)) then
               do g = 1, ubound(this%o_bnd(d, ind)%seg(:), dim=1)
                  if (proc /= this%o_bnd(d, ind)%seg(g)%proc) then
                     tag = this%o_bnd(d, ind)%seg(g)%lh + HI*d
                     nr = nr + 1
                     ! for noncartesian division some y-boundary corner cells are independent from x-boundary face cells, (similarly for z-direction).
                     if (ind == ARR) then
                        call MPI_Isend(pa3d(1, 1, 1), 1, this%o_bnd(d, ind)%seg(g)%mbc, this%o_bnd(d, ind)%seg(g)%proc, tag, comm, req(nr), ierr)
                     else
                        call MPI_Isend(pa4d(1, 1, 1, 1), 1, this%o_bnd(d, ind)%seg(g)%mbc, this%o_bnd(d, ind)%seg(g)%proc, tag, comm, req(nr), ierr)
                     endif
                  endif
               enddo
            endif
            if (ubound(this%i_bnd(d, ind)%seg(:), dim=1) /= ubound(this%o_bnd(d, ind)%seg(:), dim=1)) call die("g:ib u/=u")
            if (nr>0) call MPI_Waitall(nr, req(:nr), status(:,:nr), ierr)
         endif
      enddo

   end subroutine internal_boundaries

!>
!! \brief Routines that deallocates directional meshes.
!<

   subroutine cleanup(this)

      use domain,    only: has_dir
      use mpisetup,  only: ierr
      use constants, only: FLUID, ARR, xdim, zdim, LO, HI, BND, BLK, INVALID

      implicit none

      class(grid_container) :: this
      integer :: d, t, g

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

      do d = xdim, zdim
         if (has_dir(d)) then
            do t = FLUID, ARR
               if (this%mbc(t, d, LO, BLK) /= INVALID) call MPI_Type_free(this%mbc(t, d, LO, BLK), ierr)
               if (this%mbc(t, d, LO, BND) /= INVALID) call MPI_Type_free(this%mbc(t, d, LO, BND), ierr)
               if (this%mbc(t, d, HI, BLK) /= INVALID) call MPI_Type_free(this%mbc(t, d, HI, BLK), ierr)
               if (this%mbc(t, d, HI, BND) /= INVALID) call MPI_Type_free(this%mbc(t, d, HI, BND), ierr)
            enddo
         endif
      enddo

      if (allocated(this%i_bnd)) then
         do d = xdim, zdim
            do t = 1, ubound(this%i_bnd, dim=2)
               if (allocated(this%i_bnd(d, t)%seg)) then
                  do g = 1, ubound(this%i_bnd(d, t)%seg(:), dim=1)
                     if (this%i_bnd(d, t)%seg(g)%mbc /= INVALID) call MPI_Type_free(this%i_bnd(d, t)%seg(g)%mbc, ierr)
                  enddo
                  deallocate(this%i_bnd(d, t)%seg)
               endif
            enddo
         enddo
         deallocate(this%i_bnd)
      endif
      if (allocated(this%o_bnd)) then
         do d = xdim, zdim
            do t = 1, ubound(this%o_bnd, dim=2)
               if (allocated(this%o_bnd(d, t)%seg)) then
                  do g = 1, ubound(this%o_bnd(d, t)%seg(:), dim=1)
                     if (this%o_bnd(d, t)%seg(g)%mbc /= INVALID) call MPI_Type_free(this%o_bnd(d, t)%seg(g)%mbc, ierr)
                  enddo
                  deallocate(this%o_bnd(d, t)%seg)
               endif
            enddo
         enddo
         deallocate(this%o_bnd)
      endif

      call this%cs_iso2%clean
      call this%gpot%clean
      call this%hgpot%clean
      call this%gp%clean
      call this%sgp%clean
      call this%sgpm%clean

   end subroutine cleanup

end module grid_cont
