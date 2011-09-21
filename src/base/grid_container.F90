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

   use constants, only: dsetnamelen, xdim, zdim, ndims, LO, HI
   use domain,    only: domain_container
   use types,     only: axes, array4d, array3d

   implicit none

   private
   public :: grid_container, segment, bnd_list

   ! specify segment of data for boundary exchange, prolongation and restriction.
   type :: segment
      integer :: proc                                     !< target process
      integer(kind=8), dimension(xdim:zdim, LO:HI) :: se  !< range
   end type segment

   !< segment type for boundary exchange
   type, extends(segment) :: bnd_segment
      integer(kind=4) :: mbc                              !< MPI Boundary conditions Container
      integer(kind=4) :: tag                              !< unique tag for data exchange
   end type bnd_segment

   type :: bnd_list
      type(bnd_segment), dimension(:), allocatable :: seg !< list of boundary segments to exchange
   end type bnd_list

   !< A named array for user-defined variables, scalar fields and similar
   type, extends(array3d):: named_array3d
      character(len=dsetnamelen) :: name !< a user-provided id for the array
      integer(kind=4) :: restart_mode !< \todo If not .true. then write names to the restart file
   end type named_array3d

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

      real, dimension(ndims)          :: dl     !< array of %grid cell sizes in all directions, [ dx, dy, dz ]
      real, dimension(ndims)          :: idl    !< array of inverted %grid cell sizes in all directions, 1./dl(:)

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

      integer(kind=8), dimension(ndims) :: off    !< offset of the local domain within computational domain
      integer(kind=4), dimension(ndims) :: n_b    !< [nxb, nyb, nzb]
      integer(kind=8), dimension(ndims) :: h_cor1 !< offsets of the corner opposite to the one defined by off(:) + 1, a shortcut to be compared with dom%n_d(:)
      integer(kind=4), dimension(ndims) :: n_     !< number of %grid cells in one block in x-, y- and z-directions (n_b(:) + 2 * nb)
      integer(kind=8), dimension(ndims, LO:HI) :: my_se !< own segment
      integer :: grid_n                           !< number of own segment: my_se(:,:) = dom%pse(proc)%sel(grid_n, :, :)

      integer(kind=4), dimension(ndims, LO:HI)  :: ijkse !< [[is, js, ks], [ie, je, ke]]
      integer, dimension(ndims, LO:HI)  :: bnd  !< type of boundary conditions coded in integers

      integer(kind=4), allocatable, dimension(:,:,:,:,:) :: mbc !< MPI Boundary conditions Container

      type(domain_container) :: dom                  !< contains domain decomposition of a domain this cg belongs to (BEWARE: antiparallel)

      !>
      !! description of incoming and outgoing boundary data,
      !! the shape is (xdim:zdim, FLUID:ARR) for cg (base grid container)) and (xdim:zdim, nb) for plvl (multigrid level container)
      !<
      type(bnd_list), dimension(:, :, :), allocatable :: i_bnd, o_bnd

      real, allocatable, dimension(:,:,:) :: gc_xdim !< array of geometrical coefficients in x-direction
      real, allocatable, dimension(:,:,:) :: gc_ydim !< array of geometrical coefficients in y-direction
      real, allocatable, dimension(:,:,:) :: gc_zdim !< array of geometrical coefficients in z-direction

      type(array4d) :: u                        !< Main array of all fluids' components
      type(array4d) :: uh                       !< Main array of all fluids' components (for t += dt/2)
      type(array4d) :: u0                       !< Copy of main array of all fluids' components
      type(array4d) :: b0                       !< Copy of main array of magnetic field's components
      type(array4d) :: b                        !< Main array of magnetic field's components

      !< Other 3D arrays (such as user-defined quantities or gravitational potential). The fourth index selects variable so it cannot be merged with u or b.
      type(named_array3d), allocatable, dimension(:) :: q
      ! handy shortcuts to some entries in q(:)
      real, dimension(:,:,:), pointer :: gpot    => null() !< Array for sum of gravitational potential at t += dt
      real, dimension(:,:,:), pointer :: hgpot   => null() !< Array for sum of gravitational potential at t += 0.5*dt
      real, dimension(:,:,:), pointer :: gp      => null() !< Array for gravitational potential from external fields
      real, dimension(:,:,:), pointer :: sgp     => null() !< Array for gravitational potential from multigrid or FFT solver
      real, dimension(:,:,:), pointer :: sgpm    => null() !< Array for gravitational potential from multigrid or FFT solver at previous timestep saved by source_terms_grav.
      real, dimension(:,:,:), pointer :: cs_iso2 => null()
      real, dimension(:,:,:), pointer :: wa      => null() !< Temporary array used for different purposes, usually has dimension (grid::nx, grid::ny, grid::nz)

   contains

      procedure :: init
      procedure :: cleanup
      procedure :: mpi_bnd_types
      procedure :: add_na
      procedure :: get_na_ptr
      procedure :: get_na_ind
      procedure :: exists

   end type grid_container

contains

!-----------------------------------------------------------------------------
!>
!! \brief Initialization of the grid container
!!
!! \details This method sets up the grid container variables, coordinates and allocates basic arrays
!<

   subroutine init(this, dom, grid_n)

      use constants,  only: PIERNIK_INIT_DOMAIN, xdim, ydim, zdim, ndims, FLUID, ARR, LO, HI, BND, BLK, INVALID, I_ONE, I_TWO, BND_MPI, BND_SHE, BND_COR
      use dataio_pub, only: die, warn, printinfo, msg, code_progress
      use domain,     only: has_dir, domain_container, cdd
      use mpi,        only: MPI_COMM_NULL
      use mpisetup,   only: proc, nproc, inflate_req

      implicit none

      class(grid_container), intent(inout), target :: this ! intent(out) would silently clear everything, that was already set (also the fields in types derived from grid_container)
      type(domain_container), intent(in) :: dom
      integer, intent(in) :: grid_n

      integer :: i
      real, dimension(:), pointer :: a0, al, ar, ia

      if (code_progress < PIERNIK_INIT_DOMAIN) call die("[grid_container:init] MPI not initialized.")

      this%dom = dom
      this%nb = dom%nb !> \todo make this one global constant back

      this%grid_n     = grid_n
      this%my_se(:,:) = dom%pse(proc)%sel(grid_n, :, :)
      this%off(:)     = this%my_se(:, LO)
      this%h_cor1(:)  = this%my_se(:, HI) + I_ONE
      this%n_b(:)     = int(this%h_cor1(:) - this%off(:), 4) ! Block 'physical' grid sizes

      if (all(this%n_b(:) == 0)) then
         this%empty = .true.
      else if (any(this%n_b(:) <= 0)) then
         call die("[grid_container:init] Mixed positive and non-positive grid sizes")
      else
         this%empty = .false.
      endif

      ! Inherit the boundaries from the domain, then set MPI or SHEAR boundaries where applicable
      this%bnd(:,:) = dom%bnd(:,:)
      where (this%off(:)    /= 0)          this%bnd(:, LO) = BND_MPI
      where (this%h_cor1(:) /= dom%n_d(:)) this%bnd(:, HI) = BND_MPI
      ! For periodic boundaries do not set BND_MPI when local domain spans through the whole computational domain in given direction.
      where (dom%periodic(:) .and. this%h_cor1(:) /= dom%n_d(:)) this%bnd(:, LO) = BND_MPI
      where (dom%periodic(:) .and. this%off(:)    /= 0)          this%bnd(:, HI) = BND_MPI

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

         if (any(dom%bnd(xdim:ydim, :) == BND_COR) .and. (cdd%psize(xdim) /= cdd%psize(ydim) .or. dom%n_d(xdim) /= dom%n_d(ydim))) then
            write(msg, '(a,4(i4,a))')"[grid_container:init] Corner BC require psize(xdim) equal to psize(ydim) and nxd equal to nyd. Detected: [", &
                 &                   cdd%psize(xdim),",",cdd%psize(ydim), "] and [",dom%n_d(xdim),",",dom%n_d(ydim),"]"
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

      !> \todo allocate this contitionally, only when comm3d is in use
      if (allocated(this%mbc)) call die("[grid_container:init] this%mbc already allocated")
      allocate(this%mbc(FLUID:ARR, xdim:zdim, LO:HI, BND:BLK, 1:this%nb))
      this%mbc(:, :, :, :, :) = INVALID

      if (this%empty) then

         this%fbnd(:,:) = dom%edge(:,:)
         this%n_(:) = 0
         this%ijkse(:, LO) = this%nb + I_ONE
         this%ijkse(:, HI) = this%nb
         this%dl(:) = 1.0

         this%isb   = 0 ! ???
         this%ieb   = 0
         this%jsb   = 0
         this%jeb   = 0
         this%ksb   = 0
         this%keb   = 0

         this%vol = 0.
         this%dvol = 0.
         this%maxxyz = 0

      else

         do i = xdim, zdim
            if (has_dir(i)) then
               if (this%n_b(i) < 1) call die("[grid_init] Too many CPUs for such a small grid.")
               if (this%n_b(i) < this%nb) call warn("[grid_init] domain size in some directions is < nb, which may result in incomplete boundary cell update")
            endif
         enddo

         where (has_dir(:))
            this%n_(:) = this%n_b(:) + I_TWO * this%nb       ! Block total grid size with guardcells
            this%ijkse(:, LO) = this%nb + I_ONE
            this%ijkse(:, HI) = this%nb + this%n_b(:)
            this%dl(:) = dom%L_(:) / dom%n_d(:)
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

         if (has_dir(xdim)) then
            this%isb   = I_TWO*this%nb
            this%ieb   = this%n_b(xdim)+I_ONE
         else
            this%isb   = 1
            this%ieb   = 1
         endif

         if (has_dir(ydim)) then
            this%jsb   = I_TWO*this%nb
            this%jeb   = this%n_b(ydim)+I_ONE
         else
            this%jsb   = 1
            this%jeb   = 1
         endif

         if (has_dir(zdim)) then
            this%ksb   = I_TWO*this%nb
            this%keb   = this%n_b(zdim)+I_ONE
         else
            this%ksb   = 1
            this%keb   = 1
         endif

         this%vol = product(this%fbnd(:, HI)-this%fbnd(:, LO), mask=has_dir(:))
         this%dvol = product(this%dl(:), mask=has_dir(:))

         this%maxxyz = maxval(this%n_(:), mask=has_dir(:))

!--- Assignments -----------------------------------------------------------
         ! left zone boundaries:  xl, yl, zl
         ! zone centers:          x,  y,  z
         ! right zone boundaries: xr, yr, zr

         allocate(this%x(this%n_(xdim)), this%xl(this%n_(xdim)), this%xr(this%n_(xdim)), this%inv_x(this%n_(xdim)))
         a0 => this%x; al => this%xl; ar => this%xr; ia => this%inv_x
         call set_axis(xdim, a0, al, ar, ia, this, dom)

         allocate(this%y(this%n_(ydim)), this%yl(this%n_(ydim)), this%yr(this%n_(ydim)), this%inv_y(this%n_(ydim)))
         a0 => this%y; al => this%yl; ar => this%yr; ia => this%inv_y
         call set_axis(ydim, a0, al, ar, ia, this, dom)

         allocate(this%z(this%n_(zdim)), this%zl(this%n_(zdim)), this%zr(this%n_(zdim)), this%inv_z(this%n_(zdim)))
         a0 => this%z; al => this%zl; ar => this%zr; ia => this%inv_z
         call set_axis(zdim, a0, al, ar, ia, this, dom)

      endif

      this%dxmn = minval(this%dl(:), mask=has_dir(:))

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

      call this%mpi_bnd_types

   end subroutine init

   subroutine set_axis(d, a0, al, ar, ia, cg, dom)

      use constants, only: LO, HI, half, one, zero
      use domain,    only: has_dir, domain_container

      implicit none

      integer(kind=4), intent(in) :: d                      !< direction
      real, dimension(:), pointer, intent(inout) :: a0, al, ar, ia !< arrays with coordinates
      class(grid_container), intent(in) :: cg
      type(domain_container), intent(in) :: dom

      integer :: i

      if (has_dir(d)) then
         a0(:) = dom%edge(d, LO) + cg%dl(d) * ([(i, i=1, cg%n_(d))] - half - cg%nb + cg%off(d))
      else
         a0(:) = half*(cg%fbnd(d, LO) + cg%fbnd(d, HI))
      endif

      al(:) = a0(:) - half*cg%dl(d)
      ar(:) = a0(:) + half*cg%dl(d)

      where ( a0(:) /= zero )
         ia(:) = one/a0(:)
      elsewhere
         ia(:) = zero
      endwhere

   end subroutine set_axis

!>
!! \brief Routines that deallocates all internals of the grid container
!<

   subroutine cleanup(this)

      use domain,    only: has_dir
      use mpisetup,  only: ierr
      use constants, only: FLUID, ARR, xdim, zdim, LO, HI, BND, BLK, INVALID

      implicit none

      class(grid_container) :: this
      integer :: d, t, g, b

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
            if (has_dir(d)) then
               do t = FLUID, ARR
                  do b = 1, this%nb
                     if (this%mbc(t, d, LO, BLK, b) /= INVALID) call MPI_Type_free(this%mbc(t, d, LO, BLK, b), ierr)
                     if (this%mbc(t, d, LO, BND, b) /= INVALID) call MPI_Type_free(this%mbc(t, d, LO, BND, b), ierr)
                     if (this%mbc(t, d, HI, BLK, b) /= INVALID) call MPI_Type_free(this%mbc(t, d, HI, BLK, b), ierr)
                     if (this%mbc(t, d, HI, BND, b) /= INVALID) call MPI_Type_free(this%mbc(t, d, HI, BND, b), ierr)
                  enddo
               enddo
            endif
         enddo
         deallocate(this%mbc)
      endif

      if (allocated(this%i_bnd)) then
         do d = xdim, zdim
            do t = lbound(this%i_bnd, dim=2), ubound(this%i_bnd, dim=2)
               do b = 1, ubound(this%i_bnd, dim=3)
                  if (allocated(this%i_bnd(d, t, b)%seg)) then
                     do g = 1, ubound(this%i_bnd(d, t, b)%seg(:), dim=1)
                        if (this%i_bnd(d, t, b)%seg(g)%mbc /= INVALID) call MPI_Type_free(this%i_bnd(d, t, b)%seg(g)%mbc, ierr)
                     enddo
                     deallocate(this%i_bnd(d, t, b)%seg)
                  endif
               enddo
            enddo
         enddo
         deallocate(this%i_bnd)
      endif
      if (allocated(this%o_bnd)) then
         do d = xdim, zdim
            do t = lbound(this%i_bnd, dim=2), ubound(this%o_bnd, dim=2)
                do b = 1, ubound(this%o_bnd, dim=3)
                   if (allocated(this%o_bnd(d, t, b)%seg)) then
                      do g = 1, ubound(this%o_bnd(d, t, b)%seg(:), dim=1)
                         if (this%o_bnd(d, t, b)%seg(g)%mbc /= INVALID) call MPI_Type_free(this%o_bnd(d, t, b)%seg(g)%mbc, ierr)
                      enddo
                      deallocate(this%o_bnd(d, t, b)%seg)
                   endif
                enddo
            enddo
         enddo
         deallocate(this%o_bnd)
      endif

      if (allocated(this%q)) then
         do g = 1, ubound(this%q(:), dim=1)
            call this%q(g)%clean
         enddo
         deallocate(this%q)
      endif

      if (allocated(this%dom%pse)) deallocate(this%dom%pse) ! this is a side effect of keepina a copy of dom instead of pointing it

   end subroutine cleanup

   subroutine mpi_bnd_types(this)

      use constants,  only: FLUID, ARR, xdim, zdim, ndims, LO, HI, BND, BLK, INVALID, I_ONE
      use dataio_pub, only: die
      use domain,     only: has_dir, dom, is_overlap, cdd
      use fluidindex, only: flind
      use mpi,        only: MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, MPI_COMM_NULL
      use mpisetup,   only: ierr, proc, FIRST, LAST, procmask

      implicit none

      class(grid_container), intent(inout) :: this

      integer(kind=4), dimension(:), allocatable :: sizes, subsizes, starts
      integer :: t, g, j, b
      integer(kind=4) :: d, hl, lh, ib
      integer(kind=4), parameter, dimension(FLUID:ARR) :: dims = [ I_ONE+ndims, I_ONE+ndims, I_ONE+ndims, ndims ] !< dimensionality of arrays
      integer(kind=4), dimension(FLUID:ARR) :: nc
      integer(kind=8), dimension(xdim:zdim) :: ijks, per
      integer(kind=8), dimension(xdim:zdim, LO:HI) :: b_layer, bp_layer, poff
      logical :: sharing

      nc = [ flind%all, ndims, max(flind%crs%all,I_ONE), I_ONE ]      !< number of fluids, magnetic field components, CRs, and 1 for a rank-3 array

      if (allocated(this%i_bnd) .or. allocated(this%o_bnd)) call die("[grid:grid_mpi_boundaries_prep] this%i_bnd or this%o_bnd already allocated")
      allocate(this%i_bnd(xdim:zdim, FLUID:ARR, this%nb), this%o_bnd(xdim:zdim, FLUID:ARR, this%nb))

      ! assume that cuboids fill the domain and don't collide

      ijks(:) = this%ijkse(:, LO) - this%off(:)
      per(:) = 0
      where (this%dom%periodic(:)) per(:) = this%dom%n_d(:)

      if (cdd%comm3d == MPI_COMM_NULL) then

         do d = xdim, zdim
            if (has_dir(d) .and. .not. this%empty) then

               ! identify processes with interesting neighbour data
               procmask(:) = 0
               do lh = LO, HI
                  hl = LO+HI-lh ! HI for LO, LO for HI
                  b_layer(:,:) = this%my_se(:, :)
                  b_layer(d, lh) = b_layer(d, lh) + lh-hl ! -1 for LO, +1 for HI
                  b_layer(d, hl) = b_layer(d, lh) ! boundary layer without corners
                  do j = FIRST, LAST
                     do b = 1, ubound(this%dom%pse(j)%sel(:, :, :), dim=1)
                        call is_overlap(b_layer(:,:), this%dom%pse(j)%sel(b, :, :), sharing, per(:))
                        if (sharing) procmask(j) = procmask(j) + 1
                     enddo
                  enddo
               enddo
               do j = FLUID, ARR
                  do ib = 1, this%nb
                     allocate(this%i_bnd(d, j, ib)%seg(sum(procmask(:))))
                     allocate(this%o_bnd(d, j, ib)%seg(sum(procmask(:))))
                  enddo
               enddo

               ! set up segments to be sent or received
               g = 0
               do j = FIRST, LAST
                  if (procmask(j) /= 0) then
                     do lh = LO, HI
                        hl = LO+HI-lh
                        do b = 1, ubound(dom%pse(j)%sel(:, :, :), dim=1)
                           b_layer(:,:) = this%my_se(:, :)
                           b_layer(d, lh) = b_layer(d, lh) + lh-hl
                           b_layer(d, hl) = b_layer(d, lh)
                           bp_layer(:, :) = b_layer(:, :)
                           where (per(:) > 0)
                              bp_layer(:, LO) = mod(b_layer(:, LO) + per(:), per(:))
                              bp_layer(:, HI) = mod(b_layer(:, HI) + per(:), per(:))
                           endwhere
                           !> \todo save b_layer(:,:) and bp_layer(:,:) and move above calculations outside the b loop

                           call is_overlap(bp_layer(:,:), this%dom%pse(j)%sel(b, :, :), sharing)

                           if (sharing) then
                              poff(:,:) = bp_layer(:,:) - b_layer(:,:) ! displacement due to periodicity
                              bp_layer(:, LO) = max(bp_layer(:, LO), this%dom%pse(j)%sel(b, :, LO))
                              bp_layer(:, HI) = min(bp_layer(:, HI), this%dom%pse(j)%sel(b, :, HI))
                              b_layer(:,:) = bp_layer(:,:) - poff(:,:)
                              g = g + 1
                              do t = FLUID, ARR
                                 do ib = 1, this%nb
                                    this%i_bnd(d, t, ib)%seg(g)%mbc = INVALID
                                    this%i_bnd(d, t, ib)%seg(g)%proc = j
                                    this%i_bnd(d, t, ib)%seg(g)%se(:,LO) = b_layer(:, LO) + ijks(:)
                                    this%i_bnd(d, t, ib)%seg(g)%se(:,HI) = b_layer(:, HI) + ijks(:)
                                    if (any(this%i_bnd(d, t, ib)%seg(g)%se(d, :) < 0)) &
                                         this%i_bnd(d, t, ib)%seg(g)%se(d, :) = this%i_bnd(d, t, ib)%seg(g)%se(d, :) + this%dom%n_d(d)
                                    if (any(this%i_bnd(d, t, ib)%seg(g)%se(d, :) > this%n_b(d) + 2*this%nb)) &
                                         this%i_bnd(d, t, ib)%seg(g)%se(d, :) = this%i_bnd(d, t, ib)%seg(g)%se(d, :) - this%dom%n_d(d)

                                    ! expand to cover corners (requires separate MPI_Waitall for each direction)
                                    ! \todo create separate %mbc for corner-less exchange with one MPI_Waitall (can scale better)
                                    where (has_dir(:d-1))
                                       this%i_bnd(d, t, ib)%seg(g)%se(:d-1, LO) = this%i_bnd(d, t, ib)%seg(g)%se(:d-1, LO) - ib
                                       this%i_bnd(d, t, ib)%seg(g)%se(:d-1, HI) = this%i_bnd(d, t, ib)%seg(g)%se(:d-1, HI) + ib
                                    endwhere
                                    this%o_bnd(d, t, ib)%seg(g) = this%i_bnd(d, t, ib)%seg(g)
                                    this%i_bnd(d, t, ib)%seg(g)%tag = int(b           + ubound(dom%pse(   j)%sel(:, :, :), dim=1) * (HI*d+lh-LO), kind=4)
                                    this%o_bnd(d, t, ib)%seg(g)%tag = int(this%grid_n + ubound(dom%pse(proc)%sel(:, :, :), dim=1) * (HI*d+hl-LO), kind=4)
                                    select case (lh)
                                       case (LO)
                                          this%i_bnd(d, t, ib)%seg(g)%se(d, LO) = this%i_bnd(d, t, ib)%seg(g)%se(d, HI) - (ib - 1)
                                          this%o_bnd(d, t, ib)%seg(g)%se(d, LO) = this%i_bnd(d, t, ib)%seg(g)%se(d, HI) + 1
                                          this%o_bnd(d, t, ib)%seg(g)%se(d, HI) = this%o_bnd(d, t, ib)%seg(g)%se(d, LO) + (ib - 1)
                                       case (HI)
                                          this%i_bnd(d, t, ib)%seg(g)%se(d, HI) = this%i_bnd(d, t, ib)%seg(g)%se(d, LO) + (ib - 1)
                                          this%o_bnd(d, t, ib)%seg(g)%se(d, HI) = this%i_bnd(d, t, ib)%seg(g)%se(d, LO) - 1
                                          this%o_bnd(d, t, ib)%seg(g)%se(d, LO) = this%o_bnd(d, t, ib)%seg(g)%se(d, HI) - (ib - 1)
                                    end select

                                    allocate(sizes(dims(t)), subsizes(dims(t)), starts(dims(t)))

                                    starts(:) = 0
                                    if (dims(t) == 1+ndims) then
                                       sizes(1) = nc(t)
                                       subsizes(1) = sizes(1)
                                    endif
                                    sizes   (dims(t)-zdim+xdim:dims(t)) = this%n_(:)
                                    subsizes(dims(t)-zdim+xdim:dims(t)) = int(this%i_bnd(d, t, ib)%seg(g)%se(:, HI) - this%i_bnd(d, t, ib)%seg(g)%se(:, LO) + 1, kind=4)
                                    starts  (dims(t)-zdim+xdim:dims(t)) = int(this%i_bnd(d, t, ib)%seg(g)%se(:, LO) - 1, kind=4)
                                    call MPI_Type_create_subarray(dims(t), sizes, subsizes, starts,  MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, this%i_bnd(d, t, ib)%seg(g)%mbc, ierr)
                                    call MPI_Type_commit(this%i_bnd(d, t, ib)%seg(g)%mbc, ierr)

                                    subsizes(dims(t)-zdim+xdim:dims(t)) = int(this%o_bnd(d, t, ib)%seg(g)%se(:, HI) - this%o_bnd(d, t, ib)%seg(g)%se(:, LO) + 1, kind=4)
                                    starts  (dims(t)-zdim+xdim:dims(t)) = int(this%o_bnd(d, t, ib)%seg(g)%se(:, LO) - 1, kind=4)
                                    call MPI_Type_create_subarray(dims(t), sizes, subsizes, starts,  MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, this%o_bnd(d, t, ib)%seg(g)%mbc, ierr)
                                    call MPI_Type_commit(this%o_bnd(d, t, ib)%seg(g)%mbc, ierr)

                                    deallocate(sizes, subsizes, starts)
                                 enddo
                              enddo
                           endif
                        enddo
                     enddo
                  endif
               enddo
            endif
         enddo

      else

         do d = xdim, zdim
            if (has_dir(d)) then
               do t = FLUID, ARR  ! fluid, Bfield, wcr, grav

                  allocate(sizes(dims(t)), subsizes(dims(t)), starts(dims(t)))
                  if (dims(t) == 1+ndims) sizes(1) = nc(t)
                  sizes(dims(t)-zdim+xdim:dims(t)) = this%n_(:)

                  do ib = 1, this%nb

                     subsizes(:) = sizes(:)
                     subsizes(dims(t)-zdim+d) = ib
                     starts(:) = 0

                     starts(dims(t)-zdim+d) = this%nb-ib
                     call MPI_Type_create_subarray(dims(t), sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, this%mbc(t, d, LO, BND, ib), ierr)
                     call MPI_Type_commit(this%mbc(t, d, LO, BND, ib), ierr)

                     starts(dims(t)-zdim+d) = this%nb
                     call MPI_Type_create_subarray(dims(t), sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, this%mbc(t, d, LO, BLK, ib), ierr)
                     call MPI_Type_commit(this%mbc(t, d, LO, BLK, ib), ierr)

                     starts(dims(t)-zdim+d) = this%n_b(d) + this%nb - ib
                     call MPI_Type_create_subarray(dims(t), sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, this%mbc(t, d, HI, BLK, ib), ierr)
                     call MPI_Type_commit(this%mbc(t, d, HI, BLK, ib), ierr)

                     starts(dims(t)-zdim+d) = this%ijkse(d, HI)
                     call MPI_Type_create_subarray(dims(t), sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, this%mbc(t, d, HI, BND, ib), ierr)
                     call MPI_Type_commit(this%mbc(t, d, HI, BND, ib), ierr)

                  enddo

                  deallocate(sizes, subsizes, starts)

               enddo
            endif
         enddo

      endif

   end subroutine mpi_bnd_types

!>
!! \brief Register a new entry in current cg with given name.
!!
!! \warning You should call it for every grid_container so each of them has the same list of named arrays.
!<

   subroutine add_na(this, name, res_m)

      use constants,  only: AT_IGNORE
      use dataio_pub, only: msg, die

      implicit none

      class(grid_container), intent(inout) :: this
      character(len=*), intent(in) :: name
      integer(kind=4), optional, intent(in) :: res_m

      type(named_array3d), allocatable, dimension(:) :: tmp
      integer :: i

      if (.not. allocated(this%q)) then
         allocate(this%q(1))
      else
         if (this%exists(name)) then
            write(msg, '(3a)')"[grid_container:add_na] Array '",trim(name),"' was already registered."
            call die(msg)
         endif
         allocate(tmp(ubound(this%q(:), dim=1) + 1))
         tmp(:ubound(this%q(:), dim=1)) = this%q(:)
         call move_alloc(from=tmp, to=this%q)
      endif

      i = ubound(this%q(:), dim=1)
      this%q(i)%name = name
      call this%q(i)%init(this%n_(:))

      this%q(i)%restart_mode = AT_IGNORE
      if (present(res_m)) this%q(i)%restart_mode = res_m

   end subroutine add_na

!>
!! \brief Get the pointer to a named array of given name.
!!
!! \details This is the preferred method to access the registered array
!<

   function get_na_ptr(this, name) result(ptr)

      use dataio_pub, only: die, warn

      implicit none

      class(grid_container), intent(inout) :: this
      character(len=*), intent(in) :: name

      real, dimension(:,:,:), pointer :: ptr
      integer :: i

      ptr => null()

      do i = 1, ubound(this%q, dim=1)
         if (trim(name) ==  this%q(i)%name) then
            if (associated(ptr)) call die("[grid_container:get_na_ptr] multiple entries with the same name")
            ptr => this%q(i)%arr
         endif
      enddo

      if (.not. associated(ptr)) call warn("[grid_container:get_na_ptr] requested entry not found")

   end function get_na_ptr

!>
!! \brief Get the index of a named array of given name.
!!
!! \details This method is provided for convenience only. Use get_na_ptr whenever possible.
!<

   function get_na_ind(this, name) result(ind)

      use dataio_pub, only: die, warn

      implicit none

      class(grid_container), intent(inout) :: this
      character(len=*), intent(in) :: name

      integer :: ind, i

      ind = 0

      do i = 1, ubound(this%q, dim=1)
         if (trim(name) ==  this%q(i)%name) then
            if (ind /= 0) call die("[grid_container:get_na_ind] multiple entries with the same name")
            ind = i
         endif
      enddo

      if (ind == 0) call warn("[grid_container:get_na_ind] requested entry not found")

   end function get_na_ind

!>
!! \brief Check if an array of given name is already registered
!<

   function exists(this, name)

      use dataio_pub, only: die

      implicit none

      class(grid_container), intent(inout) :: this
      character(len=*), intent(in) :: name

      logical :: exists
      integer :: i

      exists = .false.

      do i = 1, ubound(this%q, dim=1)
         if (trim(name) ==  this%q(i)%name) then
            if (exists) call die("[grid_container:exists] multiple entries with the same name")
            exists = .true.
         endif
      enddo

   end function exists

end module grid_cont
