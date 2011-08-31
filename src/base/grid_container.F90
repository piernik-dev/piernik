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
   use types,     only: axes, array4d, array3d

   implicit none

   private
   public :: grid_container, segment, bnd_list

   ! specify segment of data for boundary exchange, prolongation and restriction.
   type :: segment
      integer :: proc                                     !< target process
      integer :: ngc                                      !< number of grid container on the target process
      integer(kind=8), dimension(xdim:zdim, LO:HI) :: se  !< range
   end type segment

   !< segment type for boundary exchange
   type, extends(segment) :: bnd_segment
      integer(kind=4) :: mbc                              !< Multigrid MPI Boundary conditions Container
      integer(kind=4) :: lh                               !< low or high boundary; \todo store full tag here
   end type bnd_segment

   type :: bnd_list
      type(bnd_segment), dimension(:), allocatable :: seg !< list of boundary segments to exchange
   end type bnd_list

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

   end type grid_container

contains

!-----------------------------------------------------------------------------
   subroutine init(this, dom, grid_n)

      use constants,  only: PIERNIK_INIT_DOMAIN, xdim, ydim, zdim, INVALID, I_ONE, I_TWO, BND_MPI, BND_SHE, BND_COR
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
      if (ubound(dom%pse(proc)%sel(:,:,:), dim=1) > 1) call die("[grid_container:init] Multiple blocks per process not implemented yet")

      this%nb = dom%nb

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
         call inflate_req(size([LO, HI]) * 2 * nproc) ! 2 = count([i_bnd, o_bnd])
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

      this%mbc(:, :, :, :) = INVALID

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
