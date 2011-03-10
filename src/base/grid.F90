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
!! \brief (DW) Module containing routines to specify required computational mesh.
!! \date January/February 2006
!!
!!
!! In this module a following namelist of parameters is specified:
!! \copydetails grid::init_grid
!<
module grid

   use types,    only: grid_container

   implicit none

   private
   public :: cleanup_grid, init_grid, grid_mpi_boundaries_prep, arr3d_boundaries
   public :: total_ncells, cg, D_x, D_y, D_z

   integer, protected :: total_ncells                   !< total number of %grid cells
   integer, protected :: D_x                            !< set to 1 when x-direction exists, 0 otherwise. Use to construct dimensionally-safe indices for arrays
   integer, protected :: D_y                            !< set to 1 when y-direction exists, 0 otherwise.
   integer, protected :: D_z                            !< set to 1 when z-direction exists, 0 otherwise.
   type(grid_container), protected :: cg                !< A container for the grid. For AMR this will be a dynamically resized array

contains

!>
!! \brief Routine which sets numbers of cells for the domain, MPI blocks and initializes direction meshes (x,y,z).
!! Also compute domain maximum and minimum of coordinates, lengths of cells and coordinates of zone centers and left/right zone boundaries.
!!
!<
   subroutine init_grid

      use constants,  only: PIERNIK_INIT_MPI, xdim, ydim, zdim
      use dataio_pub, only: printinfo, die, code_progress
      use mpisetup,   only: dom, has_dir

      implicit none

      if (code_progress < PIERNIK_INIT_MPI) call die("[grid:init_grid] MPI not initialized.")

#ifdef VERBOSE
      call printinfo("[grid:init_grid]: commencing...")
#endif /* VERBOSE */

      call set_cg(cg, dom)

      D_x = 0; D_y = 0; D_z = 0
      if (has_dir(xdim)) D_x = 1
      if (has_dir(ydim)) D_y = 1
      if (has_dir(zdim)) D_z = 1
      total_ncells = product(dom%n_d(:))

#ifdef VERBOSE
      call printinfo("[grid:init_grid]: finished. \o/")
#endif /* VERBOSE */

   end subroutine init_grid

   subroutine set_cg(this, dom)

      use constants,  only: PIERNIK_INIT_MPI, xdim, ydim, zdim
      use dataio_pub, only: die, warn, code_progress
      use mpisetup,   only: psize, pcoords, has_dir, translate_bnds_to_ints_dom
      use types,      only: domain_container

      implicit none

      type(grid_container), intent(out) :: this
      type(domain_container), intent(in) :: dom

      integer :: i, j, k

      if (code_progress < PIERNIK_INIT_MPI) call die("[grid:set_cg] MPI not initialized.")

      this%nb = dom%nb
      this%dxmn = huge(1.0)

      this%off(:) = 0
      this%n_b(:) = 1
      where (has_dir(:))
         this%off(:) = (dom%n_d(:) * pcoords(:) ) / psize(:) ! Block offset on the dom% should be between 0 and nxd-nxb
         this%n_b(:) = (dom%n_d(:) * (pcoords(:)+1))/psize(:) - this%off(:)  ! Block 'physical' grid sizes
      endwhere

      do i = xdim, zdim
         if (has_dir(i)) then
            if (this%n_b(i) < 1) call die("[grid_set_cg] Too many CPUs for a small grid.")
            if (this%n_b(i) < this%nb) call warn("[grid_set_cg] domain size in some directions is < nb, which may result in incomplete boundary cell update")
         endif
      enddo

      this%nxb = this%n_b(xdim)
      this%nyb = this%n_b(ydim)
      this%nzb = this%n_b(zdim)

      if (has_dir(xdim)) then
         this%nx    = this%nxb + 2 * this%nb       ! Block total grid sizes
         this%is    = this%nb + 1
         this%ie    = this%nb + this%nxb
         this%isb   = 2*this%nb
         this%ieb   = this%nxb+1
         this%dx    = dom%Lx / dom%n_d(xdim)
         this%dxmn  = min(this%dxmn, this%dx)
         this%xminb = dom%xmin + this%dx *  this%off(xdim)
         this%xmaxb = dom%xmin + this%dx * (this%off(xdim) + this%nxb)
      else
         this%nx    = 1
         this%is    = 1
         this%ie    = 1
         this%isb   = 1
         this%ieb   = 1
         this%dx    = 1.0
         this%xminb = dom%xmin
         this%xmaxb = dom%xmax
      endif

      if (has_dir(ydim)) then
         this%ny    = this%nyb + 2 * this%nb
         this%js    = this%nb + 1
         this%je    = this%nb + this%nyb
         this%jsb   = 2*this%nb
         this%jeb   = this%nyb+1
         this%dy    = dom%Ly / dom%n_d(ydim)
         this%dxmn  = min(this%dxmn, this%dy)
         this%yminb = dom%ymin + this%dy *  this%off(ydim)
         this%ymaxb = dom%ymin + this%dy * (this%off(ydim) + this%nyb)
      else
         this%ny    = 1
         this%js    = 1
         this%je    = 1
         this%jsb   = 1
         this%jeb   = 1
         this%dy    = 1.0
         this%yminb = dom%ymin
         this%ymaxb = dom%ymax
      endif

      if (has_dir(zdim)) then
         this%nz    = this%nzb + 2 * this%nb
         this%ks    = this%nb + 1
         this%ke    = this%nb + this%nzb
         this%ksb   = 2*this%nb
         this%keb   = this%nzb+1
         this%dz    = dom%Lz / dom%n_d(zdim)
         this%dxmn  = min(this%dxmn, this%dz)
         this%zminb = dom%zmin + this%dz *  this%off(zdim)
         this%zmaxb = dom%zmin + this%dz * (this%off(zdim) + this%nzb)
      else
         this%nz    = 1
         this%ks    = 1
         this%ke    = 1
         this%ksb   = 1
         this%keb   = 1
         this%dz    = 1.0
         this%zminb = dom%zmin
         this%zmaxb = dom%zmax
      endif

      this%idx = 1./this%dx
      this%idy = 1./this%dy
      this%idz = 1./this%dz

      this%dl(xdim:zdim) = [ this%dx, this%dy, this%dz ]
      this%idl(:) = 1./this%dl(:)

      this%dvol = product(this%dl(:))

      allocate(this%x(this%nx), this%xl(this%nx), this%xr(this%nx), this%inv_x(this%nx))
      allocate(this%y(this%ny), this%yl(this%ny), this%yr(this%ny), this%inv_y(this%ny))
      allocate(this%z(this%nz), this%zl(this%nz), this%zr(this%nz), this%inv_z(this%nz))
      this%maxxyz = maxval([size(this%x), size(this%y), size(this%z)])

      this%bnd(:,:) = translate_bnds_to_ints_dom()

!--- Assignments -----------------------------------------------------------
      ! left zone boundaries:  xl, yl, zl
      ! zone centers:          x,  y,  z
      ! right zone boundaries: xr, yr, zr

!--- x-grids --------------------------------------------------------------

      if (has_dir(xdim)) then
         do i= 1, this%nx
            this%x(i)  = this%xminb + 0.5*this%dx + (i-this%nb-1)*this%dx
            this%xl(i) = this%x(i)  - 0.5*this%dx
            this%xr(i) = this%x(i)  + 0.5*this%dx
         enddo
      else
         this%x  =  0.5*(this%xminb + this%xmaxb)
         this%xl = -0.5*this%dx
         this%xr =  0.5*this%dx
      endif
      where ( this%x /= 0.0 )
         this%inv_x = 1./this%x
      elsewhere
         this%inv_x = 0.
      endwhere

!--- y-grids --------------------------------------------------------------

      if (has_dir(ydim)) then
         do j= 1, this%ny
            this%y(j)  = this%yminb + 0.5*this%dy + (j-this%nb-1)*this%dy
            this%yl(j) = this%y(j)  - 0.5*this%dy
            this%yr(j) = this%y(j)  + 0.5*this%dy
         enddo
      else
         this%y  =  0.5*(this%yminb + this%ymaxb)
         this%yl = -0.5*this%dy
         this%yr =  0.5*this%dy
      endif
      where ( this%y /= 0.0 )
         this%inv_y = 1./this%y
      elsewhere
         this%inv_y = 0.
      endwhere

!--- z-grids --------------------------------------------------------------

      if (has_dir(zdim)) then
         do k= 1, this%nz
            this%z(k)  = this%zminb + 0.5*this%dz + (k-this%nb-1) * this%dz
            this%zl(k) = this%z(k)  - 0.5*this%dz
            this%zr(k) = this%z(k)  + 0.5*this%dz
         enddo
      else
         this%z  =  0.5*(this%zminb + this%zmaxb)
         this%zl = -0.5*this%dz
         this%zr =  0.5*this%dz
      endif
      where ( this%z /= 0.0 )
         this%inv_z = 1./this%z
      elsewhere
         this%inv_z = 0.
      endwhere

   end subroutine set_cg

!>
!! \brief Set up subsets of u,b and sgp arrays for MPI communication
!<

   subroutine grid_mpi_boundaries_prep(numfluids)

      use dataio_pub, only: die, code_progress
      use constants,  only: PIERNIK_INIT_BASE, FLUID, ARR, xdim, zdim, ndims, LO, HI, BND, DOM
      use mpi,        only: MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION
      use mpisetup,   only: ierr, has_dir

      implicit none

      integer, intent(in) :: numfluids !< expect flind%all, here, cannot grab it directly because of cyclic deps in CR-based setups

      integer, dimension(:), allocatable :: sizes, subsizes, starts
      integer, parameter :: NOT_EXIST = -1
      integer :: d, t
      integer, dimension(FLUID:ARR) :: nc
      integer, parameter, dimension(FLUID:ARR) :: dim = [ 1+ndims, 1+ndims, ndims ] !< dimensionality of arrays
      integer, dimension(xdim:zdim) :: HBstart

      if (code_progress < PIERNIK_INIT_BASE) call die("[grid:grid_mpi_boundaries_prep] grid or fluids not initialized.")

      cg%mbc(:, :, :, :) = NOT_EXIST

      nc = [ numfluids, ndims, 1 ]      !< number of fluids, magnetic field components and 1 for rank-3 array
      HBstart = [ cg%ie, cg%je, cg%ke ]
      do d = xdim, zdim
         if (has_dir(d)) then
            do t = FLUID, ARR  ! fluid, Bfield, grav

               allocate(sizes(dim(t)), subsizes(dim(t)), starts(dim(t)))

               if (dim(t) == 1+ndims) sizes(1) = nc(t)
               sizes(dim(t)-zdim+xdim:dim(t)) = [ cg%nx, cg%ny, cg%nz ]
               subsizes(:) = sizes(:)
               subsizes(dim(t)-zdim+d) = cg%nb

               starts(:) = 0
               call MPI_Type_create_subarray(dim(t), sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, cg%mbc(t, d, LO, BND), ierr)
               call MPI_Type_commit(cg%mbc(t, d, LO, BND), ierr)

               starts(dim(t)-zdim+d) = cg%nb
               call MPI_Type_create_subarray(dim(t), sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, cg%mbc(t, d, LO, DOM), ierr)
               call MPI_Type_commit(cg%mbc(t, d, LO, DOM), ierr)

               starts(dim(t)-zdim+d) = cg%n_b(d)
               call MPI_Type_create_subarray(dim(t), sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, cg%mbc(t, d, HI, DOM), ierr)
               call MPI_Type_commit(cg%mbc(t, d, HI, DOM), ierr)

               starts(dim(t)-zdim+d) = HBstart(d)
               call MPI_Type_create_subarray(dim(t), sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, cg%mbc(t, d, HI, BND), ierr)
               call MPI_Type_commit(cg%mbc(t, d, HI, BND), ierr)

               deallocate(sizes, subsizes, starts)

            enddo
         endif
      enddo

   end subroutine grid_mpi_boundaries_prep


!-----------------------------------------------------------------------------

   subroutine arr3d_boundaries(pa3d)

      use dataio_pub,    only: die
      use mpi,           only: MPI_STATUS_SIZE, MPI_REQUEST_NULL
      use constants,     only: ARR, xdim, ydim, zdim, LO, HI, BND, DOM, BND_PER, BND_MPI, BND_SHE, BND_COR
      use mpisetup,      only: comm3d, ierr, has_dir, psize, procxl, procyl, proczl, procxr, procyr, proczr

      implicit none

      real, dimension(:,:,:), pointer, intent(inout) :: pa3d

      integer, parameter                          :: nreq = 3 * 4
      integer, dimension(nreq)                    :: req3d
      integer, dimension(MPI_STATUS_SIZE, nreq)   :: status3d
      integer                                     :: i, d, lh, p

      !> \todo fill corners with big_float ?

      req3d(:) = MPI_REQUEST_NULL

      do d = xdim, zdim
         if (has_dir(d)) then
            do lh = LO, HI

               select case (cg%bnd(d, lh))
                  case (BND_PER)
                     do i = 1, ceiling(cg%nb/real(cg%n_b(d))) ! Repeating is important for domains that are narrower than their guardcells (e.g. cg%n_b(d) = 2)
                        select case (2*d+lh)
                           case (2*xdim+LO)
                              pa3d(1:cg%nb, :, :) = pa3d(cg%ieb:cg%ie, :, :)
                           case (2*ydim+LO)
                              pa3d(:, 1:cg%nb, :) = pa3d(:, cg%jeb:cg%je, :)
                           case (2*zdim+LO)
                              pa3d(:, :, 1:cg%nb) = pa3d(:, :, cg%keb:cg%ke)
                           case (2*xdim+HI)
                              pa3d(cg%ie+1:cg%nx, :, :) = pa3d(cg%is:cg%isb, :, :)
                           case (2*ydim+HI)
                              pa3d(:, cg%je+1:cg%ny, :) = pa3d(:, cg%js:cg%jsb, :)
                           case (2*zdim+HI)
                              pa3d(:, :, cg%ke+1:cg%nz) = pa3d(:, :, cg%ks:cg%ksb)
                        end select
                     enddo
                  case (BND_MPI)
                     if (psize(d) > 1) then
                        select case (2*d+lh)
                           case (2*xdim+LO)
                              p = procxl
                           case (2*ydim+LO)
                              p = procyl
                           case (2*zdim+LO)
                              p = proczl
                           case (2*xdim+HI)
                              p = procxr
                           case (2*ydim+HI)
                              p = procyr
                           case (2*zdim+HI)
                              p = proczr
                        end select
                        call MPI_Isend(pa3d(1, 1, 1), 1, cg%mbc(ARR, d, lh, DOM),  p, 2*d+(LO+HI-lh), comm3d, req3d(4*(d-xdim)+1+2*(lh-LO)), ierr)
                        call MPI_Irecv(pa3d(1, 1, 1), 1, cg%mbc(ARR, d, lh, BND),  p, 2*d+       lh,  comm3d, req3d(4*(d-xdim)+2+2*(lh-LO)), ierr)
                     else
                        call die("[mpisetup:arr3d_boundaries] bnd_[xyz][lr] == 'mpi' && psize([xyz]dim) <= 1")
                     endif
                  case (BND_SHE, BND_COR) !> \todo move appropriate code from poissonsolver::poisson_solve or do nothing. or die until someone really needs SHEAR.
                     call die("[mpisetup:arr3d_boundaries] 'she' not implemented")
                  case default ! Set gradient == 0 on the external boundaries
                     do i = 1, cg%nb
                        select case (2*d+lh)
                           case (2*xdim+LO)
                              pa3d(i, :, :) = pa3d(cg%is, :, :)
                           case (2*ydim+LO)
                              pa3d(:, i, :) = pa3d(:, cg%js, :)
                           case (2*zdim+LO)
                              pa3d(:, :, i) = pa3d(:, :, cg%ks)
                           case (2*xdim+HI)
                              pa3d(cg%ie+i, :, :) = pa3d(cg%ie, :, :)
                           case (2*ydim+HI)
                              pa3d(:, cg%je+i, :) = pa3d(:, cg%je, :)
                           case (2*zdim+HI)
                              pa3d(:, :, cg%ke+i) = pa3d(:, :, cg%ke)
                        end select
                     enddo
               end select

            enddo
         endif
      enddo

      call MPI_Waitall(nreq, req3d(:), status3d(:,:), ierr)

   end subroutine arr3d_boundaries

!>
!! \brief Routines that deallocates directional meshes.
!<
   subroutine cleanup_grid

      use mpisetup,  only: has_dir, ierr
      use constants, only: FLUID, ARR, xdim, zdim, LO, HI, BND, DOM

      implicit none

      integer :: d, t

      if (allocated(cg%x))     deallocate(cg%x)
      if (allocated(cg%xl))    deallocate(cg%xl)
      if (allocated(cg%xr))    deallocate(cg%xr)
      if (allocated(cg%inv_x)) deallocate(cg%inv_x)
      if (allocated(cg%y))     deallocate(cg%y)
      if (allocated(cg%yl))    deallocate(cg%yl)
      if (allocated(cg%yr))    deallocate(cg%yr)
      if (allocated(cg%inv_y)) deallocate(cg%inv_y)
      if (allocated(cg%z))     deallocate(cg%z)
      if (allocated(cg%zl))    deallocate(cg%zl)
      if (allocated(cg%zr))    deallocate(cg%zr)
      if (allocated(cg%inv_z)) deallocate(cg%inv_z)

      do d = xdim, zdim
         if (has_dir(d)) then
            do t = FLUID, ARR
               call MPI_Type_free(cg%mbc(t, d, LO, BND),  ierr)
               call MPI_Type_free(cg%mbc(t, d, LO, DOM),  ierr)
               call MPI_Type_free(cg%mbc(t, d, HI, DOM), ierr)
               call MPI_Type_free(cg%mbc(t, d, HI, BND), ierr)
            enddo
         endif
      enddo

   end subroutine cleanup_grid

end module grid
