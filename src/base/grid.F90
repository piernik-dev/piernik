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
!    Initial implemetation of PIERNIK code was based on TVD split MHD code by
!    Ue-Li Pen
!        see: Pen, Arras & Wong (2003) for algorithm and
!             http://www.cita.utoronto.ca/~pen/MHD
!             for original source code "mhd.f90"
!
!    For full list of developers see $PIERNIK_HOME/license/pdt.txt
!
#include "piernik.def"
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

   implicit none

   real    :: dx                             !< length of the %grid cell in x-direction
   real    :: dy                             !< length of the %grid cell in y-direction
   real    :: dz                             !< length of the %grid cell in z-direction
   real    :: dxmn                           !< the smallest length of the %grid cell (among dx, dy, and dz)
   real    :: dvol                           !< volume of one %grid cell
   integer :: nxd                            !< number of %grid cells in physical domain (without boundary cells) in x-direction (if equal to 1 then x-dimension is reduced to a point and boundary cells layer is not added)
   integer :: nyd                            !< number of %grid cells in physical domain (without boundary cells) in y-direction (if equal to 1 then y-dimension is reduced to a point and boundary cells layer is not added)
   integer :: nzd                            !< number of %grid cells in physical domain (without boundary cells) in z-direction (if equal to 1 then z-dimension is reduced to a point and boundary cells layer is not added)
   integer :: nb                             !< number of boundary cells surrounding the physical domain, same for all directions
   integer :: nx                             !< number of %grid cells in one block in x-direction
   integer :: ny                             !< number of %grid cells in one block in y-direction
   integer :: nz                             !< number of %grid cells in one block in z-direction
   integer :: nxb                            !< number of physical domain %grid cells in one block (without boundary cells) in x-direction
   integer :: nyb                            !< number of physical domain %grid cells in one block (without boundary cells) in y-direction
   integer :: nzb                            !< number of physical domain %grid cells in one block (without boundary cells) in z-direction
   integer :: nxt                            !< total number of %grid cells in the whole domain in x-direction
   integer :: nyt                            !< total number of %grid cells in the whole domain in y-direction
   integer :: nzt                            !< total number of %grid cells in the whole domain in z-direction
   integer :: is                             !< index of the first %grid cell of physical domain in x-direction
   integer :: ie                             !< index of the last %grid cell of physical domain in x-direction
   integer :: js                             !< index of the first %grid cell of physical domain in y-direction
   integer :: je                             !< index of the last %grid cell of physical domain in y-direction
   integer :: ks                             !< index of the first %grid cell of physical domain in z-direction
   integer :: ke                             !< index of the last %grid cell of physical domain in z-direction
   integer :: maxxyz                         !< maximum number of %grid cells in any direction

   real    :: xmin                           !< physical domain left x-boundary position
   real    :: xmax                           !< physical domain right x-boundary position
   real    :: ymin                           !< physical domain left y-boundary position
   real    :: ymax                           !< physical domain right y-boundary position
   real    :: zmin                           !< physical domain left z-boundary position
   real    :: zmax                           !< physical domain right z-boundary position
   real    :: xminb                          !< current block left x-boundary position
   real    :: xmaxb                          !< current block right x-boundary position
   real    :: yminb                          !< current block left y-boundary position
   real    :: ymaxb                          !< current block right y-boundary position
   real    :: zminb                          !< current block left z-boundary position
   real    :: zmaxb                          !< current block right z-boundary position
   real    :: Lx                             !< span of the physical domain in x-direction (xmax-xmin)
   real    :: Ly                             !< span of the physical domain in y-direction (ymax-ymin)
   real    :: Lz                             !< span of the physical domain in z-direction (zmax-zmin)
   integer,parameter  :: xdim=1              !< parameter assigned to x-direction
   integer,parameter  :: ydim=2              !< parameter assigned to y-direction
   integer,parameter  :: zdim=3              !< parameter assigned to z-direction

   real, allocatable, target :: dl(:)               !< array of %grid cell sizes in all directions
   real, allocatable, dimension(:), target :: x     !< array of x-positions of %grid cells centers
   real, allocatable, dimension(:), target :: y     !< array of y-positions of %grid cells centers
   real, allocatable, dimension(:), target :: z     !< array of z-positions of %grid cells centers
   real, allocatable, dimension(:), target :: xl    !< array of x-positions of %grid cells left borders
   real, allocatable, dimension(:), target :: yl    !< array of y-positions of %grid cells left borders
   real, allocatable, dimension(:), target :: zl    !< array of z-positions of %grid cells left borders
   real, allocatable, dimension(:), target :: xr    !< array of x-positions of %grid cells right borders
   real, allocatable, dimension(:), target :: yr    !< array of y-positions of %grid cells right borders
   real, allocatable, dimension(:), target :: zr    !< array of z-positions of %grid cells right borders

   contains
!>
!! \brief Routine which sets numbers of cells for the domain, MPI blocks and initializes direction meshes (x,y,z).
!!
!! \n \n
!! @b DOMAIN_SIZES
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>nxd</td><td>1</td><td>positive integer    </td><td>\copydoc grid::nxd</td></tr>
!! <tr><td>nyd</td><td>1</td><td>positive integer    </td><td>\copydoc grid::nyd</td></tr>
!! <tr><td>nzd</td><td>1</td><td>positive integer    </td><td>\copydoc grid::nzd</td></tr>
!! <tr><td>nb </td><td>4</td><td>non-negative integer</td><td>\copydoc grid::nb </td></tr>
!! </table>
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
   subroutine set_container_grid(cgrid)

      use types,           only: grid_container
      use dataio_public,   only: die, namelist_errh

      implicit none

      type(grid_container), intent(out) :: cgrid

      cgrid%dx = dx; cgrid%dy = dy; cgrid%dz = dz
      cgrid%dxmn = minval([dx,dy,dz])
      cgrid%dvol = dx*dy*dz
      cgrid%nb = nb
      cgrid%nxd = nxd; cgrid%nyd = nyd; cgrid%nzd = nzd
      cgrid%nx = nx; cgrid%ny = ny; cgrid%nz = nz
      cgrid%nxb = nxb; cgrid%nyb = nyb; cgrid%nzb = nzb
      cgrid%nxt = nxt; cgrid%nyt = nyt; cgrid%nzt = nzt
      cgrid%is = is; cgrid%js = js; cgrid%ks = ks
      cgrid%ie = ie; cgrid%je = je; cgrid%ks = ke

      cgrid%xmin = xmin; cgrid%ymin = ymin; cgrid%zmin = zmin
      cgrid%xmax = xmax; cgrid%ymax = ymax; cgrid%zmax = zmax
      cgrid%xminb = xminb; cgrid%yminb = yminb; cgrid%zminb = zminb
      cgrid%xmaxb = xmaxb; cgrid%ymaxb = ymaxb; cgrid%zmaxb = zmaxb
      cgrid%Lx = Lx; cgrid%Ly = Ly; cgrid%Lz = Lz

      !Association check is not required for uninitialized variables.
      cgrid%dl=>dl
      cgrid%x=>x
      cgrid%y=>y
      cgrid%z=>z
      cgrid%xl=>xl
      cgrid%yl=>yl
      cgrid%zl=>zl
      cgrid%xr=>xr
      cgrid%yr=>yr
      cgrid%zr=>zr

   end subroutine set_container_grid

   subroutine init_grid(cgrid)

      use types,         only: grid_container
      use mpisetup,      only: ierr, ibuff, rbuff, MPI_INTEGER, MPI_DOUBLE_PRECISION, proc, &
           &                buffer_dim, pxsize, pysize, pzsize, comm
      use dataio_public, only: par_file, ierrh, die, printinfo, namelist_errh, compare_namelist

      implicit none

      type(grid_container), intent(out) :: cgrid

      namelist /DOMAIN_SIZES/ nxd, nyd, nzd, nb
      namelist /DOMAIN_LIMITS/ xmin, xmax, ymin, ymax, zmin, zmax

#ifdef VERBOSE
      call printinfo("[grid:init_grid]: commencing...")
#endif /* VERBOSE */

      nxd  = 1
      nyd  = 1
      nzd  = 1
      nb   = 4

      if (proc == 0) then
         diff_nml(DOMAIN_SIZES)
         diff_nml(DOMAIN_LIMITS)
      endif

      nxd = max(1, nxd)
      nyd = max(1, nyd)
      nzd = max(1, nzd)

      if (proc == 0) then

         ibuff(1)   = nxd
         ibuff(2)   = nyd
         ibuff(3)   = nzd
         ibuff(4)   = nb

         rbuff(1)   = xmin
         rbuff(2)   = xmax
         rbuff(3)   = ymin
         rbuff(4)   = ymax
         rbuff(5)   = zmin
         rbuff(6)   = zmax

      endif

      call MPI_Bcast(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_Bcast(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      if (proc /= 0) then

         nxd  = ibuff(1)
         nyd  = ibuff(2)
         nzd  = ibuff(3)
         nb   = ibuff(4)

         xmin = rbuff(1)
         xmax = rbuff(2)
         ymin = rbuff(3)
         ymax = rbuff(4)
         zmin = rbuff(5)
         zmax = rbuff(6)

      endif

      if ((mod(nxd, pxsize) .ne. 0) .or. &
         (mod(nyd, pysize) .ne. 0) .or. &
         (mod(nzd, pzsize) .ne. 0)) then
            call die("One of: (mod(n_d,p_size) /= 0")
      endif

      nxb = nxd/pxsize     !
      nyb = nyd/pysize     ! Block 'physical' grid sizes
      nzb = nzd/pzsize     !

      nx=nxb+2*nb          !
      ny=nyb+2*nb          ! Block total grid sizes
      nz=nzb+2*nb          !

      nxt=nxd+2*nb         !
      nyt=nyd+2*nb         ! Domain total grid sizes
      nzt=nzd+2*nb         !

      if (nxd == 1) then
         nx     = 1
         nxb    = 1
         nxt    = 1
         pxsize = 1
         is     = 1
         ie     = 1
      else
         is = nb+1
         ie = nb+nxb
      endif

      if (nyd == 1) then
         ny     = 1
         nyb    = 1
         nyt    = 1
         pysize = 1
         js     = 1
         je     = 1
      else
         js = nb+1
         je = nb+nyb
      endif

      if (nzd == 1) then
         nz     = 1
         nzb    = 1
         nzt    = 1
         pzsize = 1
         ks     = 1
         ke     = 1
      else
         ks = nb+1
         ke = nb+nzb
      endif
      allocate(dl(3))
      allocate(x(nx), xl(nx), xr(nx))
      allocate(y(ny), yl(ny), yr(ny))
      allocate(z(nz), zl(nz), zr(nz))

      call grid_xyz

      call set_container_grid(cgrid)
#ifdef VERBOSE
      call printinfo("[grid:init_grid]: finished. \o/")
#endif /* VERBOSE */
   end subroutine init_grid
!>
!! \brief Routine that computes domain maximum and minimum of coordinates, lengths of cells and coordinates of zone centers and left/right zone boundaries.
!<
   subroutine grid_xyz

      use mpisetup, only: psize, pcoords

      implicit none

      integer :: i,j,k

      maxxyz = max(size(x),size(y))
      maxxyz = max(size(z),maxxyz)

      xminb = xmin + real(pcoords(1)  )*(xmax-xmin)/real(psize(1))
      xmaxb = xmin + real(pcoords(1)+1)*(xmax-xmin)/real(psize(1))
      yminb = ymin + real(pcoords(2)  )*(ymax-ymin)/real(psize(2))
      ymaxb = ymin + real(pcoords(2)+1)*(ymax-ymin)/real(psize(2))
      zminb = zmin + real(pcoords(3)  )*(zmax-zmin)/real(psize(3))
      zmaxb = zmin + real(pcoords(3)+1)*(zmax-zmin)/real(psize(3))

      dxmn = huge(1.0)
      if (nxd /= 1) then
         dx = (xmaxb-xminb)/nxb
         dxmn = min(dxmn,dx)
      else
         dx = 1.0
      endif
      if (nyd /= 1) then
         dy = (ymaxb-yminb)/nyb
         dxmn = min(dxmn,dy)
      else
         dy = 1.0
      endif
      if (nzd /= 1) then
         dz = (zmaxb-zminb)/nzb
         dxmn = min(dxmn,dz)
      else
         dz = 1.0
      endif

      dl(xdim) = dx
      dl(ydim) = dy
      dl(zdim) = dz

      dvol = dx*dy*dz

!--- Asignments -----------------------------------------------------------
    ! left zone boundaries:  xl, yl, zl
    ! zone centers:          x,  y,  z
    ! right zone boundaries: xr, yr, zr

!--- x-grids --------------------------------------------------------------

      if (nxd /= 1) then
         do i= 1, nx
            x(i)  = xminb + 0.5*dx + (i-nb-1)*dx
            xl(i) = x(i)  - 0.5*dx
            xr(i) = x(i)  + 0.5*dx
         enddo
      else
         x  =  0.5*(xminb + xmaxb)
         xl = -0.5*dx
         xr =  0.5*dx
      endif

!--- y-grids --------------------------------------------------------------

      if (nyd /= 1) then
         do j= 1, ny
            y(j)  = yminb + 0.5*dy + (j-nb-1)*dy
            yl(j) = y(j)  - 0.5*dy
            yr(j) = y(j)  + 0.5*dy
         enddo
      else
         y  =  0.5*(yminb + ymaxb)
         yl = -0.5*dy
         yr =  0.5*dy
      endif

!--- z-grids --------------------------------------------------------------

      if (nzd /= 1) then
         do k= 1, nz
            z(k)  = zminb + 0.5*dz + (k-nb-1) * dz
            zl(k) = z(k)  - 0.5*dz
            zr(k) = z(k)  + 0.5*dz
         enddo
      else
         z  =  0.5*(zminb + zmaxb)
         zl = -0.5*dz
         zr =  0.5*dz
      endif
!--------------------------------------------------------------------------

      Lx = xmax - xmin
      Ly = ymax - ymin
      Lz = zmax - zmin

   end subroutine grid_xyz
!>
!! \brief Routines that deallocates directional meshes.
!<
   subroutine cleanup_grid

      implicit none

      if (allocated(dl)) deallocate(dl)
      if (allocated(x))  deallocate(x)
      if (allocated(xl)) deallocate(xl)
      if (allocated(xr)) deallocate(xr)
      if (allocated(y))  deallocate(y)
      if (allocated(yl)) deallocate(yl)
      if (allocated(yr)) deallocate(yr)
      if (allocated(z))  deallocate(z)
      if (allocated(zl)) deallocate(zl)
      if (allocated(zr)) deallocate(zr)

   end subroutine cleanup_grid

end module grid
