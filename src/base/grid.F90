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
!! \brief [DW] Module containing routines to specify required computational mesh.
!! \date January/February 2006
!!
!!
!! In this module two namelists of parameters are specified:
!! \copydetails grid::init_grid
!<
module grid

   use mpisetup, only: cbuff_len
   use types,    only: grid_container

   implicit none

   private
   public :: cleanup_grid, init_grid, total_ncells, geometry, cg, D_x, D_y, D_z

   integer, protected :: total_ncells                   !< total number of %grid cells
   integer, protected :: D_x, D_y, D_z                  !< set to 1 when given direction exists, 0 otherwise. Use to construct dimensionally-safe indices for arrays
   character(len=cbuff_len), protected :: geometry      !< define system of coordinates
   type(grid_container), protected :: cg             ! AMR: this will be a dynamically resized array

   contains

!>
!! \brief Routine which sets numbers of cells for the domain, MPI blocks and initializes direction meshes (x,y,z).
!! Also compute domain maximum and minimum of coordinates, lengths of cells and coordinates of zone centers and left/right zone boundaries.
!!
!! \n \n
!! @b DOMAIN_LIMITS
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>xmin</td><td></td><td>real</td><td>\copydoc grid::xmin</td></tr>
!! <tr><td>xmax</td><td></td><td>real</td><td>\copydoc grid::xmax</td></tr>
!! <tr><td>ymin</td><td></td><td>real</td><td>\copydoc grid::ymin</td></tr>
!! <tr><td>ymax</td><td></td><td>real</td><td>\copydoc grid::ymax</td></tr>
!! <tr><td>zmin</td><td></td><td>real</td><td>\copydoc grid::zmin</td></tr>
!! <tr><td>zmax</td><td></td><td>real</td><td>\copydoc grid::zmax</td></tr>
!!</table>
!! \n \n
!<
   subroutine init_grid

      use dataio_pub, only: par_file, ierrh, namelist_errh, compare_namelist  ! QA_WARN required for diff_nml
      use dataio_pub, only: printinfo, die
      use mpisetup,   only: ierr, rbuff, cbuff, master, slave, buffer_dim, psize, pxsize, pysize, pzsize, pcoords, comm, &
           &                has_dir, xdim, ydim, zdim, ndims, nxd, nyd, nzd, nb
      use mpi,        only: MPI_DOUBLE_PRECISION, MPI_CHARACTER

      implicit none

      integer :: i, j, k
      real    :: xmin, xmax, ymin, ymax, zmin, zmax

      namelist /DOMAIN_LIMITS/ xmin, xmax, ymin, ymax, zmin, zmax, geometry

#ifdef VERBOSE
      call printinfo("[grid:init_grid]: commencing...")
#endif /* VERBOSE */

      xmin = 0.; xmax = 1.
      ymin = 0.; ymax = 1.
      zmin = 0.; zmax = 1.
      geometry = "cartesian"

      if (master) then

         diff_nml(DOMAIN_LIMITS)

         rbuff(1)   = xmin
         rbuff(2)   = xmax
         rbuff(3)   = ymin
         rbuff(4)   = ymax
         rbuff(5)   = zmin
         rbuff(6)   = zmax

         cbuff(1)   = geometry

      endif

      call MPI_Bcast(cbuff, cbuff_len*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_Bcast(rbuff,           buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      if (slave) then

         xmin = rbuff(1)
         xmax = rbuff(2)
         ymin = rbuff(3)
         ymax = rbuff(4)
         zmin = rbuff(5)
         zmax = rbuff(6)

         geometry = cbuff(1)

      endif

      cg%xmin = xmin; cg%ymin = ymin; cg%zmin = zmin
      cg%xmax = xmax; cg%ymax = ymax; cg%zmax = zmax

      if ( (mod(nxd, pxsize) .ne. 0) .or. &
           (mod(nyd, pysize) .ne. 0) .or. &
           (mod(nzd, pzsize) .ne. 0) ) then
         call die("One of: (mod(n_d,p_size) /= 0")
      endif

      if (has_dir(xdim)) then
         cg%nxb = nxd / pxsize     ! Block 'physical' grid sizes
         cg%nx  = cg%nxb + 2 * nb     ! Block total grid sizes
         cg%nxt = nxd + 2 * nb     ! Domain total grid sizes
         cg%is  = nb + 1
         cg%ie  = nb + cg%nxb
         cg%isb = 2*nb
         cg%ieb = cg%nxb+1
         D_x = 1
      else
         cg%nxb    = 1
         cg%nx     = 1
         cg%nxt    = 1
         pxsize = 1
         cg%is     = 1
         cg%ie     = 1
         cg%isb    = 1
         cg%ieb    = 1
         D_x    = 0
      endif

      if (has_dir(ydim)) then
         cg%nyb = nyd / pysize
         cg%ny  = cg%nyb + 2 * nb
         cg%nyt = nyd + 2 * nb
         cg%js  = nb + 1
         cg%je  = nb + cg%nyb
         cg%jsb = 2*nb
         cg%jeb = cg%nyb+1
         D_y = 1
      else
         cg%ny     = 1
         cg%nyb    = 1
         cg%nyt    = 1
         pysize = 1
         cg%js     = 1
         cg%je     = 1
         cg%jsb    = 1
         cg%jeb    = 1
         D_y    = 0
      endif

      if (has_dir(zdim)) then
         cg%nzb = nzd / pzsize
         cg%nz  = cg%nzb + 2 * nb
         cg%nzt = nzd + 2 * nb
         cg%ks  = nb + 1
         cg%ke  = nb + cg%nzb
         cg%ksb = 2*nb
         cg%keb = cg%nzb+1
         D_z = 1
      else
         cg%nzb    = 1
         cg%nz     = 1
         cg%nzt    = 1
         pzsize = 1
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

      total_ncells = nxd * nyd * nzd

      cg%maxxyz = maxval([size(cg%x), size(cg%y), size(cg%z)])

      cg%xminb = cg%xmin + real(pcoords(1)  )*(cg%xmax-cg%xmin)/real(psize(1))
      cg%xmaxb = cg%xmin + real(pcoords(1)+1)*(cg%xmax-cg%xmin)/real(psize(1))
      cg%yminb = cg%ymin + real(pcoords(2)  )*(cg%ymax-cg%ymin)/real(psize(2))
      cg%ymaxb = cg%ymin + real(pcoords(2)+1)*(cg%ymax-cg%ymin)/real(psize(2))
      cg%zminb = cg%zmin + real(pcoords(3)  )*(cg%zmax-cg%zmin)/real(psize(3))
      cg%zmaxb = cg%zmin + real(pcoords(3)+1)*(cg%zmax-cg%zmin)/real(psize(3))

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

!--------------------------------------------------------------------------

      cg%Lx = cg%xmax - cg%xmin
      cg%Ly = cg%ymax - cg%ymin
      cg%Lz = cg%zmax - cg%zmin

      cg%Vol = 1.
      if (has_dir(xdim)) cg%Vol = cg%Vol * cg%Lx
      if (has_dir(ydim)) cg%Vol = cg%Vol * cg%Ly
      if (has_dir(zdim)) cg%Vol = cg%Vol * cg%Lz
      ! BEWARE: not true for non-cartesian geometry

#ifdef VERBOSE
      call printinfo("[grid:init_grid]: finished. \o/")
#endif /* VERBOSE */
   end subroutine init_grid
!>
!! \brief Routines that deallocates directional meshes.
!<
   subroutine cleanup_grid

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

   end subroutine cleanup_grid

end module grid
