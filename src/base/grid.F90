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

   use mpisetup, only: geometry
   use types,    only: grid_container

   implicit none

   private
   public :: cleanup_grid, init_grid, grid_mpi_boundaries_prep
   public :: total_ncells, geometry, cg, D_x, D_y, D_z

   integer, protected :: total_ncells                   !< total number of %grid cells
   integer, protected :: D_x                            !< set to 1 when x-direction exists, 0 otherwise. Use to construct dimensionally-safe indices for arrays
   integer, protected :: D_y                            !< set to 1 when y-direction exists, 0 otherwise.
   integer, protected :: D_z                            !< set to 1 when z-direction exists, 0 otherwise.
   type(grid_container),     protected :: cg            !< A container for the grid. For AMR this will be a dynamically resized array

contains

!>
!! \brief Routine which sets numbers of cells for the domain, MPI blocks and initializes direction meshes (x,y,z).
!! Also compute domain maximum and minimum of coordinates, lengths of cells and coordinates of zone centers and left/right zone boundaries.
!!
!<
   subroutine init_grid

      use dataio_pub, only: par_file, ierrh, namelist_errh, compare_namelist, cmdl_nml  ! QA_WARN required for diff_nml
      use dataio_pub, only: printinfo, die, code_progress, PIERNIK_INIT_MPI
      use mpisetup,   only: psize, pcoords, comm, has_dir, xdim, ydim, zdim, ndims, dom, nb

      implicit none

      integer :: i, j, k

      if (code_progress < PIERNIK_INIT_MPI) call die("[grid:init_grid] MPI not initialized.")

#ifdef VERBOSE
      call printinfo("[grid:init_grid]: commencing...")
#endif /* VERBOSE */

      ! \todo Check if that statement is still necessary
      if ( (mod(dom%nxd, psize(xdim)) .ne. 0) .or. &
           (mod(dom%nyd, psize(ydim)) .ne. 0) .or. &
           (mod(dom%nzd, psize(zdim)) .ne. 0) ) then
         call die("One of: (mod(n_d,p_size) /= 0")
      endif

      if (has_dir(xdim)) then
         cg%nxb = dom%nxd / psize(xdim)     ! Block 'physical' grid sizes
         cg%nx  = cg%nxb + 2 * nb     ! Block total grid sizes
         cg%is  = nb + 1
         cg%ie  = nb + cg%nxb
         cg%isb = 2*nb
         cg%ieb = cg%nxb+1
         D_x = 1
      else
         cg%nxb    = 1
         cg%nx     = 1
         cg%is     = 1
         cg%ie     = 1
         cg%isb    = 1
         cg%ieb    = 1
         D_x    = 0
      endif

      if (has_dir(ydim)) then
         cg%nyb = dom%nyd / psize(ydim)
         cg%ny  = cg%nyb + 2 * nb
         cg%js  = nb + 1
         cg%je  = nb + cg%nyb
         cg%jsb = 2*nb
         cg%jeb = cg%nyb+1
         D_y = 1
      else
         cg%ny     = 1
         cg%nyb    = 1
         cg%js     = 1
         cg%je     = 1
         cg%jsb    = 1
         cg%jeb    = 1
         D_y    = 0
      endif

      if (has_dir(zdim)) then
         cg%nzb = dom%nzd / psize(zdim)
         cg%nz  = cg%nzb + 2 * nb
         cg%ks  = nb + 1
         cg%ke  = nb + cg%nzb
         cg%ksb = 2*nb
         cg%keb = cg%nzb+1
         D_z = 1
      else
         cg%nzb    = 1
         cg%nz     = 1
         cg%ks     = 1
         cg%ke     = 1
         cg%ksb    = 1
         cg%keb    = 1
         D_z    = 0
      endif
      cg%nb = nb

      allocate(cg%dl(ndims))
      allocate(cg%idl(ndims))
      allocate(cg%x(cg%nx), cg%xl(cg%nx), cg%xr(cg%nx), cg%inv_x(cg%nx))
      allocate(cg%y(cg%ny), cg%yl(cg%ny), cg%yr(cg%ny), cg%inv_y(cg%ny))
      allocate(cg%z(cg%nz), cg%zl(cg%nz), cg%zr(cg%nz), cg%inv_z(cg%nz))

      total_ncells = dom%nxd * dom%nyd * dom%nzd

      cg%maxxyz = maxval([size(cg%x), size(cg%y), size(cg%z)])

      cg%xminb = dom%xmin + real(pcoords(xdim)  )*dom%Lx/real(psize(xdim))
      cg%xmaxb = dom%xmin + real(pcoords(xdim)+1)*dom%Lx/real(psize(xdim))
      cg%yminb = dom%ymin + real(pcoords(ydim)  )*dom%Ly/real(psize(ydim))
      cg%ymaxb = dom%ymin + real(pcoords(ydim)+1)*dom%Ly/real(psize(ydim))
      cg%zminb = dom%zmin + real(pcoords(zdim)  )*dom%Lz/real(psize(zdim))
      cg%zmaxb = dom%zmin + real(pcoords(zdim)+1)*dom%Lz/real(psize(zdim))

      cg%dxmn = huge(1.0)
      if (has_dir(xdim)) then
         cg%dx = (cg%xmaxb-cg%xminb)/cg%nxb
         cg%dxmn = min(cg%dxmn, cg%dx)
      else
         cg%dx = 1.0
      endif
      cg%idx = 1./cg%dx
      if (has_dir(ydim)) then
         cg%dy = (cg%ymaxb-cg%yminb)/cg%nyb
         cg%dxmn = min(cg%dxmn, cg%dy)
      else
         cg%dy = 1.0
      endif
      cg%idy = 1./cg%dy
      if (has_dir(zdim)) then
         cg%dz = (cg%zmaxb-cg%zminb)/cg%nzb
         cg%dxmn = min(cg%dxmn, cg%dz)
      else
         cg%dz = 1.0
      endif
      cg%idz = 1./cg%dz

      cg%dl(xdim) = cg%dx
      cg%dl(ydim) = cg%dy
      cg%dl(zdim) = cg%dz
      cg%idl = 1./cg%dl

      cg%dvol = cg%dx*cg%dy*cg%dz

!--- Assignments -----------------------------------------------------------
    ! left zone boundaries:  xl, yl, zl
    ! zone centers:          x,  y,  z
    ! right zone boundaries: xr, yr, zr

!--- x-grids --------------------------------------------------------------

      if (has_dir(xdim)) then
         do i= 1, cg%nx
            cg%x(i)  = cg%xminb + 0.5*cg%dx + (i-nb-1)*cg%dx
            cg%xl(i) = cg%x(i)  - 0.5*cg%dx
            cg%xr(i) = cg%x(i)  + 0.5*cg%dx
         enddo
      else
         cg%x  =  0.5*(cg%xminb + cg%xmaxb)
         cg%xl = -0.5*cg%dx
         cg%xr =  0.5*cg%dx
      endif
      where ( cg%x /= 0.0 )
         cg%inv_x = 1./cg%x
      elsewhere
         cg%inv_x = 0.
      endwhere

!--- y-grids --------------------------------------------------------------

      if (has_dir(ydim)) then
         do j= 1, cg%ny
            cg%y(j)  = cg%yminb + 0.5*cg%dy + (j-nb-1)*cg%dy
            cg%yl(j) = cg%y(j)  - 0.5*cg%dy
            cg%yr(j) = cg%y(j)  + 0.5*cg%dy
         enddo
      else
         cg%y  =  0.5*(cg%yminb + cg%ymaxb)
         cg%yl = -0.5*cg%dy
         cg%yr =  0.5*cg%dy
      endif
      where ( cg%y /= 0.0 )
         cg%inv_y = 1./cg%y
      elsewhere
         cg%inv_y = 0.
      endwhere

!--- z-grids --------------------------------------------------------------

      if (has_dir(zdim)) then
         do k= 1, cg%nz
            cg%z(k)  = cg%zminb + 0.5*cg%dz + (k-nb-1) * cg%dz
            cg%zl(k) = cg%z(k)  - 0.5*cg%dz
            cg%zr(k) = cg%z(k)  + 0.5*cg%dz
         enddo
      else
         cg%z  =  0.5*(cg%zminb + cg%zmaxb)
         cg%zl = -0.5*cg%dz
         cg%zr =  0.5*cg%dz
      endif
      where ( cg%z /= 0.0 )
         cg%inv_z = 1./cg%z
      elsewhere
         cg%inv_z = 0.
      endwhere

#ifdef VERBOSE
      call printinfo("[grid:init_grid]: finished. \o/")
#endif /* VERBOSE */
   end subroutine init_grid

!>
!! \brief Set up subsets of u,b and sgp arrays for MPI communication
!<

   subroutine grid_mpi_boundaries_prep(numfluids)

      use dataio_pub, only: die, code_progress, PIERNIK_INIT_BASE
      use mpi,        only: MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION
      use mpisetup,   only: ierr, xdim, ydim, zdim, has_dir, ndims

      implicit none

      integer, intent(in) :: numfluids !< expect flind%all, here, cannot grab it directly because of cyclic deps in CR-based setups

      integer, dimension(:), allocatable :: sizes, subsizes, starts
      integer(kind=4) :: ord
      integer(kind=4) :: old
      integer, parameter :: dim4 = ndims +1 !< dimensionality of compund arrays: (q, x, y, z)
      integer, parameter :: dim3 = ndims    !< dimensionality of simple arrays: (x, y, z)
      integer, parameter :: NOT_EXIST = -1

      if (code_progress < PIERNIK_INIT_BASE) call die("[grid:grid_mpi_boundaries_prep] grid or fluids not initialized.")

      ord = MPI_ORDER_FORTRAN
      old = MPI_DOUBLE_PRECISION

!------------------------!
!   X dimension - fluid  !
!------------------------!
      if (has_dir(xdim)) then

         allocate(sizes(dim4), subsizes(dim4), starts(dim4))

         sizes    = [ numfluids, cg%nx, cg%ny, cg%nz ]
         subsizes = [ numfluids, cg%nb, cg%ny, cg%nz ]
         starts   = [ 0, 0, 0, 0 ]

         call MPI_Type_create_subarray(dim4, sizes, subsizes, starts, ord, old, cg%MPI_YZ_LEFT_BND,  ierr)
         call MPI_Type_commit(cg%MPI_YZ_LEFT_BND, ierr)

         starts(xdim + 1) = cg%nb
         call MPI_Type_create_subarray(dim4, sizes, subsizes, starts, ord, old, cg%MPI_YZ_LEFT_DOM,  ierr)
         call MPI_Type_commit(cg%MPI_YZ_LEFT_DOM, ierr)

         starts(xdim + 1) = cg%nxb
         call MPI_Type_create_subarray(dim4, sizes, subsizes, starts, ord, old, cg%MPI_YZ_RIGHT_DOM, ierr)
         call MPI_Type_commit(cg%MPI_YZ_RIGHT_DOM, ierr)

         starts(xdim + 1) = cg%ie
         call MPI_Type_create_subarray(dim4, sizes, subsizes, starts, ord, old, cg%MPI_YZ_RIGHT_BND, ierr)
         call MPI_Type_commit(cg%MPI_YZ_RIGHT_BND, ierr)

!------------------------!
!   X dimension - Bfield !
!------------------------!
         sizes    = [ ndims, cg%nx, cg%ny, cg%nz ]
         subsizes = [ ndims, cg%nb, cg%ny, cg%nz ]
         starts   = [ 0, 0, 0, 0 ]

         call MPI_Type_create_subarray(dim4, sizes, subsizes, starts, ord, old, cg%MAG_YZ_LEFT_BND,  ierr)
         call MPI_Type_commit(cg%MAG_YZ_LEFT_BND, ierr)

         starts(xdim + 1) = cg%nb
         call MPI_Type_create_subarray(dim4, sizes, subsizes, starts, ord, old, cg%MAG_YZ_LEFT_DOM,  ierr)
         call MPI_Type_commit(cg%MAG_YZ_LEFT_DOM, ierr)

         starts(xdim + 1) = cg%nxb
         call MPI_Type_create_subarray(dim4, sizes, subsizes, starts, ord, old, cg%MAG_YZ_RIGHT_DOM, ierr)
         call MPI_Type_commit(cg%MAG_YZ_RIGHT_DOM, ierr)

         starts(xdim + 1) = cg%ie
         call MPI_Type_create_subarray(dim4, sizes, subsizes, starts, ord, old, cg%MAG_YZ_RIGHT_BND, ierr)
         call MPI_Type_commit(cg%MAG_YZ_RIGHT_BND, ierr)

         deallocate(sizes, subsizes, starts)

!---------------------------------------!
!   X dimension - cg%nx*cg%ny*cg%nz array (grav) !
!---------------------------------------!
         allocate(sizes(dim3), subsizes(dim3), starts(dim3))

         sizes    = [ cg%nx, cg%ny, cg%nz ]
         subsizes = [ cg%nb, cg%ny, cg%nz ]
         starts   = [ 0, 0, 0 ]

         call MPI_Type_create_subarray(dim3, sizes, subsizes, starts, ord, old, cg%ARR_YZ_LEFT_BND,  ierr)
         call MPI_Type_commit(cg%ARR_YZ_LEFT_BND,  ierr)

         starts(xdim) = cg%nb
         call MPI_Type_create_subarray(dim3, sizes, subsizes, starts, ord, old, cg%ARR_YZ_LEFT_DOM,  ierr)
         call MPI_Type_commit(cg%ARR_YZ_LEFT_DOM,  ierr)

         starts(xdim) = cg%nxb
         call MPI_Type_create_subarray(dim3, sizes, subsizes, starts, ord, old, cg%ARR_YZ_RIGHT_DOM, ierr)
         call MPI_Type_commit(cg%ARR_YZ_RIGHT_DOM, ierr)

         starts(xdim) = cg%ie
         call MPI_Type_create_subarray(dim3, sizes, subsizes, starts, ord, old, cg%ARR_YZ_RIGHT_BND, ierr)
         call MPI_Type_commit(cg%ARR_YZ_RIGHT_BND, ierr)

         deallocate(sizes, subsizes, starts)

      else

         cg%MPI_YZ_LEFT_BND  = NOT_EXIST
         cg%MPI_YZ_RIGHT_BND = NOT_EXIST
         cg%MPI_YZ_LEFT_DOM  = NOT_EXIST
         cg%MPI_YZ_RIGHT_DOM = NOT_EXIST

         cg%MAG_YZ_LEFT_BND  = NOT_EXIST
         cg%MAG_YZ_RIGHT_BND = NOT_EXIST
         cg%MAG_YZ_LEFT_DOM  = NOT_EXIST
         cg%MAG_YZ_RIGHT_DOM = NOT_EXIST

         cg%ARR_YZ_LEFT_BND  = NOT_EXIST
         cg%ARR_YZ_RIGHT_BND = NOT_EXIST
         cg%ARR_YZ_LEFT_DOM  = NOT_EXIST
         cg%ARR_YZ_RIGHT_DOM = NOT_EXIST

      endif

!------------------------!
!   Y dimension - fluid  !
!------------------------!
      if (has_dir(ydim)) then

         allocate(sizes(dim4), subsizes(dim4), starts(dim4))

         sizes    = [ numfluids, cg%nx, cg%ny, cg%nz ]
         subsizes = [ numfluids, cg%nx, cg%nb, cg%nz ]
         starts   = [ 0, 0, 0, 0 ]

         call MPI_Type_create_subarray(dim4, sizes, subsizes, starts, ord, old, cg%MPI_XZ_LEFT_BND,  ierr)
         call MPI_Type_commit(cg%MPI_XZ_LEFT_BND, ierr)

         starts(ydim + 1) = cg%nb
         call MPI_Type_create_subarray(dim4, sizes, subsizes, starts, ord, old, cg%MPI_XZ_LEFT_DOM,  ierr)
         call MPI_Type_commit(cg%MPI_XZ_LEFT_DOM, ierr)

         starts(ydim + 1) = cg%nyb
         call MPI_Type_create_subarray(dim4, sizes, subsizes, starts, ord, old, cg%MPI_XZ_RIGHT_DOM, ierr)
         call MPI_Type_commit(cg%MPI_XZ_RIGHT_DOM, ierr)

         starts(ydim + 1) = cg%je
         call MPI_Type_create_subarray(dim4, sizes, subsizes, starts, ord, old, cg%MPI_XZ_RIGHT_BND, ierr)
         call MPI_Type_commit(cg%MPI_XZ_RIGHT_BND, ierr)

!------------------------!
!   Y dimension - Bfield !
!------------------------!
         sizes    = [ ndims, cg%nx, cg%ny, cg%nz ]
         subsizes = [ ndims, cg%nx, cg%nb, cg%nz ]
         starts   = [ 0, 0, 0, 0 ]

         call MPI_Type_create_subarray(dim4, sizes, subsizes, starts, ord, old, cg%MAG_XZ_LEFT_BND,  ierr)
         call MPI_Type_commit(cg%MAG_XZ_LEFT_BND, ierr)

         starts(ydim + 1) = cg%nb
         call MPI_Type_create_subarray(dim4, sizes, subsizes, starts, ord, old, cg%MAG_XZ_LEFT_DOM,  ierr)
         call MPI_Type_commit(cg%MAG_XZ_LEFT_DOM, ierr)

         starts(ydim + 1) = cg%nyb
         call MPI_Type_create_subarray(dim4, sizes, subsizes, starts, ord, old, cg%MAG_XZ_RIGHT_DOM, ierr)
         call MPI_Type_commit(cg%MAG_XZ_RIGHT_DOM, ierr)

         starts(ydim + 1) = cg%je
         call MPI_Type_create_subarray(dim4, sizes, subsizes, starts, ord, old, cg%MAG_XZ_RIGHT_BND, ierr)
         call MPI_Type_commit(cg%MAG_XZ_RIGHT_BND, ierr)

         deallocate(sizes, subsizes, starts)

!---------------------------------------!
!   Y dimension - cg%nx*cg%ny*cg%nz array (grav) !
!---------------------------------------!
         allocate(sizes(dim3), subsizes(dim3), starts(dim3))

         sizes    = [ cg%nx, cg%ny, cg%nz ]
         subsizes = [ cg%nx, cg%nb, cg%nz ]
         starts   = [ 0, 0, 0 ]

         call MPI_Type_create_subarray(dim3, sizes, subsizes, starts, ord, old, cg%ARR_XZ_LEFT_BND,  ierr)
         call MPI_Type_commit(cg%ARR_XZ_LEFT_BND,  ierr)

         starts(ydim) = cg%nb
         call MPI_Type_create_subarray(dim3, sizes, subsizes, starts, ord, old, cg%ARR_XZ_LEFT_DOM,  ierr)
         call MPI_Type_commit(cg%ARR_XZ_LEFT_DOM,  ierr)

         starts(ydim) = cg%nyb
         call MPI_Type_create_subarray(dim3, sizes, subsizes, starts, ord, old, cg%ARR_XZ_RIGHT_DOM, ierr)
         call MPI_Type_commit(cg%ARR_XZ_RIGHT_DOM, ierr)

         starts(ydim) = cg%je
         call MPI_Type_create_subarray(dim3, sizes, subsizes, starts, ord, old, cg%ARR_XZ_RIGHT_BND, ierr)
         call MPI_Type_commit(cg%ARR_XZ_RIGHT_BND, ierr)

         deallocate(sizes, subsizes, starts)

      else

         cg%MPI_XZ_LEFT_BND  = NOT_EXIST
         cg%MPI_XZ_RIGHT_BND = NOT_EXIST
         cg%MPI_XZ_LEFT_DOM  = NOT_EXIST
         cg%MPI_XZ_RIGHT_DOM = NOT_EXIST

         cg%MAG_XZ_LEFT_BND  = NOT_EXIST
         cg%MAG_XZ_RIGHT_BND = NOT_EXIST
         cg%MAG_XZ_LEFT_DOM  = NOT_EXIST
         cg%MAG_XZ_RIGHT_DOM = NOT_EXIST

         cg%ARR_XZ_LEFT_BND  = NOT_EXIST
         cg%ARR_XZ_RIGHT_BND = NOT_EXIST
         cg%ARR_XZ_LEFT_DOM  = NOT_EXIST
         cg%ARR_XZ_RIGHT_DOM = NOT_EXIST

      endif

!------------------------!
!   Z dimension - fluid  !
!------------------------!
      if (has_dir(zdim)) then

         allocate(sizes(dim4), subsizes(dim4), starts(dim4))

         sizes    = [ numfluids, cg%nx, cg%ny, cg%nz ]
         subsizes = [ numfluids, cg%nx, cg%ny, cg%nb ]
         starts   = [ 0, 0, 0, 0 ]

         call MPI_Type_create_subarray(dim4, sizes, subsizes, starts, ord, old, cg%MPI_XY_LEFT_BND,  ierr)
         call MPI_Type_commit(cg%MPI_XY_LEFT_BND, ierr)

         starts(zdim + 1) = cg%nb
         call MPI_Type_create_subarray(dim4, sizes, subsizes, starts, ord, old, cg%MPI_XY_LEFT_DOM,  ierr)
         call MPI_Type_commit(cg%MPI_XY_LEFT_DOM, ierr)

         starts(zdim + 1) = cg%nzb
         call MPI_Type_create_subarray(dim4, sizes, subsizes, starts, ord, old, cg%MPI_XY_RIGHT_DOM, ierr)
         call MPI_Type_commit(cg%MPI_XY_RIGHT_DOM, ierr)

         starts(zdim + 1) = cg%ke
         call MPI_Type_create_subarray(dim4, sizes, subsizes, starts, ord, old, cg%MPI_XY_RIGHT_BND, ierr)
         call MPI_Type_commit(cg%MPI_XY_RIGHT_BND, ierr)

!------------------------!
!   Z dimension - Bfield !
!------------------------!
         sizes    = [ ndims, cg%nx, cg%ny, cg%nz ]
         subsizes = [ ndims, cg%nx, cg%ny, cg%nb ]
         starts   = [ 0, 0, 0, 0 ]

         call MPI_Type_create_subarray(dim4, sizes, subsizes, starts, ord, old, cg%MAG_XY_LEFT_BND,  ierr)
         call MPI_Type_commit(cg%MAG_XY_LEFT_BND, ierr)

         starts(zdim + 1) = cg%nb
         call MPI_Type_create_subarray(dim4, sizes, subsizes, starts, ord, old, cg%MAG_XY_LEFT_DOM,  ierr)
         call MPI_Type_commit(cg%MAG_XY_LEFT_DOM, ierr)

         starts(zdim + 1) = cg%nzb
         call MPI_Type_create_subarray(dim4, sizes, subsizes, starts, ord, old, cg%MAG_XY_RIGHT_DOM, ierr)
         call MPI_Type_commit(cg%MAG_XY_RIGHT_DOM, ierr)

         starts(zdim + 1) = cg%ke
         call MPI_Type_create_subarray(dim4, sizes, subsizes, starts, ord, old, cg%MAG_XY_RIGHT_BND, ierr)
         call MPI_Type_commit(cg%MAG_XY_RIGHT_BND, ierr)

         deallocate(sizes, subsizes, starts)

!---------------------------------------!
!   Z dimension - cg%nx*cg%ny*cg%nz array (grav) !
!---------------------------------------!
         allocate(sizes(dim3), subsizes(dim3), starts(dim3))

         sizes    = [ cg%nx, cg%ny, cg%nz ]
         subsizes = [ cg%nx, cg%ny, cg%nb ]
         starts   = [ 0, 0, 0 ]

         call MPI_Type_create_subarray(dim3, sizes, subsizes, starts, ord, old, cg%ARR_XY_LEFT_BND,  ierr)
         call MPI_Type_commit(cg%ARR_XY_LEFT_BND,  ierr)

         starts(zdim) = cg%nb
         call MPI_Type_create_subarray(dim3, sizes, subsizes, starts, ord, old, cg%ARR_XY_LEFT_DOM,  ierr)
         call MPI_Type_commit(cg%ARR_XY_LEFT_DOM,  ierr)

         starts(zdim) = cg%nzb
         call MPI_Type_create_subarray(dim3, sizes, subsizes, starts, ord, old, cg%ARR_XY_RIGHT_DOM, ierr)
         call MPI_Type_commit(cg%ARR_XY_RIGHT_DOM, ierr)

         starts(zdim) = cg%ke
         call MPI_Type_create_subarray(dim3, sizes, subsizes, starts, ord, old, cg%ARR_XY_RIGHT_BND, ierr)
         call MPI_Type_commit(cg%ARR_XY_RIGHT_BND, ierr)

         deallocate(sizes, subsizes, starts)

      else

         cg%MPI_XY_LEFT_BND  = NOT_EXIST
         cg%MPI_XY_RIGHT_BND = NOT_EXIST
         cg%MPI_XY_LEFT_DOM  = NOT_EXIST
         cg%MPI_XY_RIGHT_DOM = NOT_EXIST

         cg%MAG_XY_LEFT_BND  = NOT_EXIST
         cg%MAG_XY_RIGHT_BND = NOT_EXIST
         cg%MAG_XY_LEFT_DOM  = NOT_EXIST
         cg%MAG_XY_RIGHT_DOM = NOT_EXIST

         cg%ARR_XY_LEFT_BND  = NOT_EXIST
         cg%ARR_XY_RIGHT_BND = NOT_EXIST
         cg%ARR_XY_LEFT_DOM  = NOT_EXIST
         cg%ARR_XY_RIGHT_DOM = NOT_EXIST

      endif

   end subroutine grid_mpi_boundaries_prep
!>
!! \brief Routines that deallocates directional meshes.
!<
   subroutine cleanup_grid

      use mpisetup, only: has_dir, xdim, ydim, zdim, ierr

      implicit none

      if (allocated(cg%dl))    deallocate(cg%dl)
      if (allocated(cg%idl))   deallocate(cg%idl)
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

      if (has_dir(xdim)) then
         call MPI_Type_free(cg%MPI_YZ_LEFT_BND,  ierr)
         call MPI_Type_free(cg%MPI_YZ_LEFT_DOM,  ierr)
         call MPI_Type_free(cg%MPI_YZ_RIGHT_DOM, ierr)
         call MPI_Type_free(cg%MPI_YZ_RIGHT_BND, ierr)
         call MPI_Type_free(cg%MAG_YZ_LEFT_BND,  ierr)
         call MPI_Type_free(cg%MAG_YZ_LEFT_DOM,  ierr)
         call MPI_Type_free(cg%MAG_YZ_RIGHT_DOM, ierr)
         call MPI_Type_free(cg%MAG_YZ_RIGHT_BND, ierr)
         call MPI_Type_free(cg%ARR_YZ_LEFT_BND,  ierr)
         call MPI_Type_free(cg%ARR_YZ_LEFT_DOM,  ierr)
         call MPI_Type_free(cg%ARR_YZ_RIGHT_DOM, ierr)
         call MPI_Type_free(cg%ARR_YZ_RIGHT_BND, ierr)
      endif

      if (has_dir(ydim)) then
         call MPI_Type_free(cg%MPI_XZ_LEFT_BND,  ierr)
         call MPI_Type_free(cg%MPI_XZ_LEFT_DOM,  ierr)
         call MPI_Type_free(cg%MPI_XZ_RIGHT_DOM, ierr)
         call MPI_Type_free(cg%MPI_XZ_RIGHT_BND, ierr)
         call MPI_Type_free(cg%MAG_XZ_LEFT_BND,  ierr)
         call MPI_Type_free(cg%MAG_XZ_LEFT_DOM,  ierr)
         call MPI_Type_free(cg%MAG_XZ_RIGHT_DOM, ierr)
         call MPI_Type_free(cg%MAG_XZ_RIGHT_BND, ierr)
         call MPI_Type_free(cg%ARR_XZ_LEFT_BND,  ierr)
         call MPI_Type_free(cg%ARR_XZ_LEFT_DOM,  ierr)
         call MPI_Type_free(cg%ARR_XZ_RIGHT_DOM, ierr)
         call MPI_Type_free(cg%ARR_XZ_RIGHT_BND, ierr)
      endif

      if (has_dir(zdim)) then
         call MPI_Type_free(cg%MPI_XY_LEFT_BND,  ierr)
         call MPI_Type_free(cg%MPI_XY_LEFT_DOM,  ierr)
         call MPI_Type_free(cg%MPI_XY_RIGHT_DOM, ierr)
         call MPI_Type_free(cg%MPI_XY_RIGHT_BND, ierr)
         call MPI_Type_free(cg%MAG_XY_LEFT_BND,  ierr)
         call MPI_Type_free(cg%MAG_XY_LEFT_DOM,  ierr)
         call MPI_Type_free(cg%MAG_XY_RIGHT_DOM, ierr)
         call MPI_Type_free(cg%MAG_XY_RIGHT_BND, ierr)
         call MPI_Type_free(cg%ARR_XY_LEFT_BND,  ierr)
         call MPI_Type_free(cg%ARR_XY_LEFT_DOM,  ierr)
         call MPI_Type_free(cg%ARR_XY_RIGHT_DOM, ierr)
         call MPI_Type_free(cg%ARR_XY_RIGHT_BND, ierr)
      endif

   end subroutine cleanup_grid

end module grid
