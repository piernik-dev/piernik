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
!>
!! \brief [DW] Module containing routines to specify required computational mesh.
!! \author M. Hanasz
!! \date January/February 2006
!!
!!
!! In this module two namelists of parameters are specified:
!!
!! @b DOMAIN_SIZES
!!
!! \f[
!! \begin{tabular}{ | p{3cm} | p{3cm} | p{4cm} | p{8cm} | } \hline &&&\\
!! {\bf parameter} & {\bf default value} & {\bf possible values} & {\bf description} \\ &&&\\ \hline \hline &&&\\
!! nxd & 1 & positive integer & number of grid cells in physical domain (without boundary cells) in x-direction (if equal to 1 then x-dimension is reduced to a point and boundary cells layer is not added) \\ &&&\\ \hline &&&\\
!! nyd & 1 & positive integer & number of grid cells in physical domain (without boundary cells) in y-direction (if equal to 1 then y-dimension is reduced to a point and boundary cells layer is not added) \\ &&&\\ \hline &&&\\
!! nzd & 1 & positive integer & number of grid cells in physical domain (without boundary cells) in z-direction (if equal to 1 then z-dimension is reduced to a point and boundary cells layer is not added) \\ &&&\\ \hline &&&\\
!! nb & 4 & non-negative integer & number of boundary cells surrounding the physical domain, same for all directions \\ &&&\\ \hline
!! \end{tabular} \f]
!!
!! @b DOMAIN_LIMITS
!!
!! \f[
!! \begin{tabular}{ | p{3cm} | p{3cm} | p{4cm} | p{8cm} | } \hline &&&\\
!! {\bf parameter} & {\bf default value} & {\bf possible values} & {\bf description} \\ &&&\\ \hline \hline &&&\\
!! xmin &  & real & physical domain left x-boundary position \\ &&&\\ \hline &&&\\
!! xmax &  & real & physical domain right x-boundary position \\ &&&\\ \hline &&&\\
!! ymin &  & real & physical domain left y-boundary position \\ &&&\\ \hline &&&\\
!! ymax &  & real & physical domain right y-boundary position \\ &&&\\ \hline &&&\\
!! zmin &  & real & physical domain left z-boundary position \\ &&&\\ \hline &&&\\
!! zmax &  & real & physical domain right z-boundary position \\ &&&\\ \hline
!! \end{tabular} \f]
!!
!<
module grid

! Written by: M. Hanasz, January/February 2006

   implicit none
   real    :: dx                             !< length of the grid cell in x-direction
   real    :: dy                             !< length of the grid cell in y-direction
   real    :: dz                             !< length of the grid cell in z-direction
   real    :: dxmn                           !< the smallest length of the grid cell (among dx, dy, and dz)
   real    :: dvol                           !< volume of one grid cell
   integer :: nxd                            !< number of grid cells in physical domain (without boundary cells) in x-direction
   integer :: nyd                            !< number of grid cells in physical domain (without boundary cells) in y-direction
   integer :: nzd                            !< number of grid cells in physical domain (without boundary cells) in z-direction
   integer :: nb                             !< number cells in a boundary layer
   integer :: nx                             !< number of grid cells in one block in x-direction
   integer :: ny                             !< number of grid cells in one block in y-direction
   integer :: nz                             !< number of grid cells in one block in z-direction
   integer :: nxb                            !< number of physical domain grid cells in one block (without boundary cells) in x-direction
   integer :: nyb                            !< number of physical domain grid cells in one block (without boundary cells) in y-direction
   integer :: nzb                            !< number of physical domain grid cells in one block (without boundary cells) in z-direction
   integer :: nxt                            !< total number of grid cells in the whole domain in x-direction
   integer :: nyt                            !< total number of grid cells in the whole domain in y-direction
   integer :: nzt                            !< total number of grid cells in the whole domain in z-direction
   integer :: is                             !< index of the first grid cell of physical domain in x-direction
   integer :: ie                             !< index of the last grid cell of physical domain in x-direction
   integer :: js                             !< index of the first grid cell of physical domain in y-direction
   integer :: je                             !< index of the last grid cell of physical domain in y-direction
   integer :: ks                             !< index of the first grid cell of physical domain in z-direction
   integer :: ke                             !< index of the last grid cell of physical domain in z-direction
   integer :: maxxyz                         !< maximum number of grid cells in any direction

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

   real, allocatable :: dl(:)                !< array of grid cell sizes in all directions
   real, allocatable, dimension(:)  :: x     !< array of x-positions of grid cells centers
   real, allocatable, dimension(:)  :: y     !< array of y-positions of grid cells centers
   real, allocatable, dimension(:)  :: z     !< array of z-positions of grid cells centers
   real, allocatable, dimension(:)  :: xl    !< array of x-positions of grid cells left borders
   real, allocatable, dimension(:)  :: yl    !< array of y-positions of grid cells left borders
   real, allocatable, dimension(:)  :: zl    !< array of z-positions of grid cells left borders
   real, allocatable, dimension(:)  :: xr    !< array of x-positions of grid cells right borders
   real, allocatable, dimension(:)  :: yr    !< array of y-positions of grid cells right borders
   real, allocatable, dimension(:)  :: zr    !< array of z-positions of grid cells right borders


   contains
!>
!! \brief Routine which sets numbers of cells for the domain, MPI blocks and initializes direction meshes (x,y,z).
!<
   subroutine init_grid
      use errh, only : namelist_errh
      use mpisetup
      implicit none
      integer :: ierrh
      character(LEN=100) :: par_file, tmp_log_file

      namelist /DOMAIN_SIZES/ nxd, nyd, nzd, nb
      namelist /DOMAIN_LIMITS/ xmin, xmax, ymin, ymax, zmin, zmax

      nxd  = 1
      nyd  = 1
      nzd  = 1
      nb   = 4

      if(proc == 0) then
         par_file = trim(cwd)//'/problem.par'
         tmp_log_file = trim(cwd)//'/tmp.log'
         open(1,file=par_file)
            read(unit=1,nml=DOMAIN_SIZES,iostat=ierrh)
            call namelist_errh(ierrh,'DOMAIN_SIZES')
         close(1)
         open(1,file=par_file)
            read(unit=1,nml=DOMAIN_LIMITS,iostat=ierrh)
            call namelist_errh(ierrh,'DOMAIN_LIMITS')
         close(1)
         open(3, file='tmp.log', position='append')
           write(3,nml=DOMAIN_SIZES)
           write(3,nml=DOMAIN_LIMITS)
           write(3,*)
         close(3)
      endif


      if(proc == 0) then

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

         call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
         call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      else

         call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
         call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

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

      if((mod(nxd, pxsize) .ne. 0) .or. &
         (mod(nyd, pysize) .ne. 0) .or. &
         (mod(nzd, pzsize) .ne. 0)) then
         call mpistop
         if (proc .eq. 0) then
            write(*,*) 'One of: (mod(nxd,pxsize) .or. mod(nyd,pysize) .or. mod(nzd,pzsize)) .ne. 0'
         endif
         stop
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

      if(nxd == 1) then
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

      if(nyd == 1) then
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

      if(nzd == 1) then
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

   end subroutine init_grid
!>
!! \brief Routine that computes domain maximum and minimum of coordinates, lengths of cells and coordinates of zone centers and left/right zone boundaries.
!<
   subroutine grid_xyz
      use mpisetup

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

      dxmn = 1.e20
      if(nxd /= 1) then
         dx = (xmaxb-xminb)/nxb
         dxmn = min(dxmn,dx)
      else
         dx = 1.0
      endif
      if(nyd /= 1) then
         dy = (ymaxb-yminb)/nyb
         dxmn = min(dxmn,dy)
      else
         dy = 1.0
      endif
      if(nzd /= 1) then
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

      if(nxd /= 1) then
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

      if(nyd /= 1) then
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

      if(nzd /= 1) then
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

      if(allocated(dl)) deallocate(dl)
      if(allocated(x))  deallocate(x)
      if(allocated(xl)) deallocate(xl)
      if(allocated(xr)) deallocate(xr)
      if(allocated(y))  deallocate(y)
      if(allocated(yl)) deallocate(yl)
      if(allocated(yr)) deallocate(yr)
      if(allocated(z))  deallocate(z)
      if(allocated(zl)) deallocate(zl)
      if(allocated(zr)) deallocate(zr)

   end subroutine cleanup_grid

end module grid
