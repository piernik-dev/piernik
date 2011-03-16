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

!!$ ============================================================================
!>
!! \brief This module contains various routines (interpolation, boundaries and some global reduction)
!! that are useful for all flavours of multigrid solvers.
!<

module multigridbasefuncs
! pulled by MULTIGRID
   implicit none

   private

   public :: prolong_level, restrict_all, norm_sq, substract_average, prolong_faces, zero_boundaries

contains

!!$ ============================================================================
!>
!! \brief Clear boundary values
!<

   subroutine zero_boundaries(lev)

      use multigridvars, only: lvl

      implicit none

      integer, intent(in) :: lev  !< level for which clear the boundary values

      lvl(lev)%bnd_x(:,:,:) = 0.
      lvl(lev)%bnd_y(:,:,:) = 0.
      lvl(lev)%bnd_z(:,:,:) = 0.

   end subroutine zero_boundaries

!!$ ============================================================================
!>
!! \brief Multigrid elementary operators: prolongation, restriction, norm etc.
!<

   subroutine prolong_level(lev, iv)

      use dataio_pub,            only: die
      use multigridhelpers,      only: dirty_debug, check_dirty, dirtyH
      use multigridvars,         only: plvl, lvl, level_min, level_max, ord_prolong, ngridvars
      use multigridexperimental, only: prolong_level_hord

      implicit none

      integer, intent(in)      :: lev   !< level to prolong from
      integer, intent(in)      :: iv    !< variable to be prolonged

      type(plvl), pointer :: coarse, fine

      if (lev >= level_max) return ! can't prolong finest level
      if (lev <  level_min) call die("[multigridbasefuncs:prolong_level] level <= 0.")
      if (iv < 1 .or. iv > ngridvars) call die("[multigridbasefuncs:prolong_level] Invalid variable index.")

      coarse => lvl(lev)
      fine   => lvl(lev + 1)

      if (dirty_debug) fine%mgvar(:, :, :, iv) = dirtyH

      call check_dirty(coarse%level, iv, "prolong-")

      if (ord_prolong == 0) then
         call prolong_level0(lev, iv)
      else
         call prolong_level_hord(lev, iv) ! experimental part
      endif

      call check_dirty(fine%level, iv, "prolong+")

   end subroutine prolong_level

!!$ ============================================================================
!>
!! \brief 0th order prolongation : injection
!<

   subroutine prolong_level0(lev, iv)

      use multigridvars, only: plvl, lvl
      use grid,          only: D_x, D_y, D_z

      implicit none

      integer, intent(in)      :: lev   !< level to prolong from
      integer, intent(in)      :: iv    !< variable to be prolonged

      type(plvl), pointer :: coarse, fine

      coarse => lvl(lev)
      fine   => lvl(lev + 1)

      ! No guardcells required here

      ! Possible optimization candidate: reduce L1 and L2 cache misses on both read and write (RBGS only, secondary importance)
      fine%mgvar       (fine%is    :fine%ie-D_x:(1+D_x), fine%js    :fine%je-D_y:(1+D_y), fine%ks    :fine%ke-D_z:(1+D_z), iv) = &
           coarse%mgvar(coarse%is  :coarse%ie,           coarse%js  :coarse%je,           coarse%ks  :coarse%ke,           iv)
      fine%mgvar       (fine%is+D_x:fine%ie    :(1+D_x), fine%js    :fine%je-D_y:(1+D_y), fine%ks    :fine%ke-D_z:(1+D_z), iv) = &
           fine%mgvar  (fine%is    :fine%ie-D_x:(1+D_x), fine%js    :fine%je-D_y:(1+D_y), fine%ks    :fine%ke-D_z:(1+D_z), iv)
      fine%mgvar       (fine%is    :fine%ie,             fine%js+D_y:fine%je    :(1+D_y), fine%ks    :fine%ke-D_z:(1+D_z), iv) = &
           fine%mgvar  (fine%is    :fine%ie,             fine%js    :fine%je-D_y:(1+D_y), fine%ks    :fine%ke-D_z:(1+D_z), iv)
      fine%mgvar       (fine%is    :fine%ie,             fine%js    :fine%je,             fine%ks+D_z:fine%ke    :(1+D_z), iv) = &
           fine%mgvar  (fine%is    :fine%ie,             fine%js    :fine%je,             fine%ks    :fine%ke-D_z:(1+D_z), iv)

   end subroutine prolong_level0

!!$ ============================================================================
!>
!! \brief Restriction operators
!<

   subroutine restrict_all(iv)

      use dataio_pub,         only: die
      use multigridhelpers,   only: check_dirty
      use multigridvars,      only: level_min, level_max, ngridvars, lvl

      implicit none

      integer, intent(in)      :: iv    !< variable to be restricted

      integer :: i

      if (iv < 1 .or. iv > ngridvars) call die("[multigridbasefuncs:restrict_all] Invalid variable index.")

      call check_dirty(level_max, iv, "restrict_all-")

      do i = level_max, level_min+1, -1
         call lvl(i)%restrict_level(iv)
      enddo

      call check_dirty(level_min, iv, "restrict_all+")

   end subroutine restrict_all

!!$ ============================================================================
!>
!! \brief Calculate L2 norm
!<

   subroutine norm_sq(iv, norm)

      use dataio_pub,    only: die
      use mpisetup,      only: comm3d, ierr, geometry_type
      use mpi,           only: MPI_DOUBLE_PRECISION, MPI_SUM
      use multigridvars, only: ngridvars, roof
      use constants,     only: GEO_XYZ, GEO_RPZ

      implicit none

      integer, intent(in)  :: iv   !< index of variable in lvl()%mgvar for which we want to find the norm
      real,    intent(out) :: norm !< the calculated norm

      real                 :: lsum
      integer              :: i

      if (iv <= 0 .or. iv > ngridvars) call die("[multigridbasefuncs:norm_sq] Invalid variable index")

      select case (geometry_type)
         case (GEO_XYZ)
            lsum = sum(roof%mgvar(roof%is:roof%ie, roof%js:roof%je, roof%ks:roof%ke, iv)**2) * roof%dvol
         case (GEO_RPZ)
            lsum = 0.
            do i = roof%is, roof%ie
               lsum = lsum + sum(roof%mgvar(i, roof%js:roof%je, roof%ks:roof%ke, iv)**2) * roof%dvol * roof%x(i)
            enddo
         case default
            call die("[multigridbasefuncs:norm_sq] Unsupported geometry.")
      end select
      call MPI_Allreduce(lsum, norm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm3d, ierr)
      norm = sqrt(norm)

   end subroutine norm_sq

!!$ ============================================================================
!>
!! \brief Compute the global average value and subtract it from the whole domain
!<

   subroutine substract_average(lev, iv)

      use dataio_pub,    only: die
      use mpisetup,      only: comm3d, ierr, geometry_type
      use mpi,           only: MPI_DOUBLE_PRECISION, MPI_SUM
      use multigridvars, only: lvl, level_min, level_max, ngridvars
      use constants,     only: GEO_XYZ, GEO_RPZ

      implicit none

      integer, intent(in) :: lev  !< level for which we want to subtract its average from
      integer, intent(in) :: iv   !< index of variable in lvl()%mgvar which we want to have zero average

      real                :: lsum, avg, vol
      integer             :: i

      if (lev < level_min .or. lev > level_max) call die("[multigridbasefuncs:substract_average] Invalid level number.")
      if (iv < 1 .or. iv > ngridvars) call die("[multigridbasefuncs:substract_average] Invalid variable index.")

      select case (geometry_type)
         case (GEO_XYZ)
            lsum = sum(lvl(lev)%mgvar(lvl(lev)%is:lvl(lev)%ie, lvl(lev)%js:lvl(lev)%je, lvl(lev)%ks:lvl(lev)%ke, iv)) * lvl(lev)%dvol
         case (GEO_RPZ)
            lsum = 0.
            do i = lvl(lev)%is, lvl(lev)%ie
               lsum = lsum + sum(lvl(lev)%mgvar(i, lvl(lev)%js:lvl(lev)%je, lvl(lev)%ks:lvl(lev)%ke, iv)) * lvl(lev)%dvol * lvl(lev)%x(i)
            enddo
         case default
            call die("[multigridbasefuncs:substract_average] Unsupported geometry.")
      end select
      call MPI_Allreduce(lsum, avg, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm3d, ierr)
      call MPI_Allreduce(lvl(lev)%vol, vol, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm3d, ierr) ! This probably can be calculated in init_multigrid
      avg = avg / vol

      lvl(lev)%mgvar(:, :, :, iv) = lvl(lev)%mgvar(:, :, :, iv) - avg

   end subroutine substract_average

!!$ ============================================================================
!>
!! \brief Prolong solution data at level (lev-1) to faces at level lev
!<

   subroutine prolong_faces(lev, soln)

      use grid,               only: D_x, D_y, D_z
      use mpisetup,           only: has_dir
      use constants,          only: xdim, ydim, zdim, LO, HI
      use dataio_pub,         only: die, warn
      use multigridhelpers,   only: check_dirty
      use multigridmpifuncs,  only: mpi_multigrid_bnd
      use multigridvars,      only: plvl, lvl, ord_prolong_face, level_min, level_max, extbnd_antimirror

      implicit none

      integer, intent(in) :: lev  !< level for which approximate the solution
      integer, intent(in) :: soln !< index of solution in lvl()%mgvar

      integer                       :: i, j, k
      type(plvl), pointer           :: coarse, fine
      real, parameter, dimension(3) :: p0  = [ 0.,       1.,     0.     ] ! injection
      real, parameter, dimension(3) :: p1  = [ 0.,       3./4.,  1./4.  ] ! 1D linear prolongation stencil
      real, parameter, dimension(3) :: p2i = [ -1./8.,   1.,     1./8.  ] ! 1D integral cubic prolongation stencil
      real, parameter, dimension(3) :: p2d = [ -3./32., 30./32., 5./32. ] ! 1D direct cubic prolongation stencil
      real, dimension(-1:1)         :: p
      real, dimension(-1:1,-1:1,2,2):: pp   ! 2D prolongation stencil
      real                          :: pp_norm

      if (lev < level_min .or. lev > level_max) call die("[multigridbasefuncs:prolong_faces] Invalid level")

      if (lev == level_min) then
         call warn("[multigridbasefuncs:prolong_faces] Cannot prolong anything to base level")
         return
      endif

      select case (ord_prolong_face)
         case (0)
            p(:) = p0(:)
         case (1,-1)
            p(:) = p1(:)
         case (2)
            p(:) = p2i(:)
         case (-2)
            p(:) = p2d(:)
         case default
            p(:) = p0(:)
      end select

      do i = -1, 1
         pp(i,:,1,1) = 0.5*p( i)*p(:)       ! 0.5 because of face averaging
         pp(i,:,1,2) = 0.5*p( i)*p(1:-1:-1) ! or use matmul()
         pp(i,:,2,1) = 0.5*p(-i)*p(:)
         pp(i,:,2,2) = 0.5*p(-i)*p(1:-1:-1)
      enddo

      call mpi_multigrid_bnd(lev-1, soln, 1, extbnd_antimirror) !> \deprecated BEWARE for higher prolongation order more guardcell are required
      call check_dirty(lev-1, soln, "prolong_faces", 1)

      coarse => lvl(lev - 1)
      fine   => lvl(lev)

      if (has_dir(xdim)) then
         pp_norm = 2.*sum(pp(-D_y:D_y, -D_z:D_z, 1, 1)) ! normalization is required for ord_prolong_face == 1 and -2
         do j = coarse%js, coarse%je
            do k = coarse%ks, coarse%ke
               fine%bnd_x(-fine%js+2*j,    -fine%ks+2*k,    LO) =sum(pp(-D_y:D_y, -D_z:D_z, 1, 1) * (coarse%mgvar(coarse%is,j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%is-1,j-D_y:j+D_y,k-D_z:k+D_z,soln))) / pp_norm
               fine%bnd_x(-fine%js+2*j+D_y,-fine%ks+2*k,    LO) =sum(pp(-D_y:D_y, -D_z:D_z, 2, 1) * (coarse%mgvar(coarse%is,j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%is-1,j-D_y:j+D_y,k-D_z:k+D_z,soln))) / pp_norm
               fine%bnd_x(-fine%js+2*j,    -fine%ks+2*k+D_z,LO) =sum(pp(-D_y:D_y, -D_z:D_z, 1, 2) * (coarse%mgvar(coarse%is,j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%is-1,j-D_y:j+D_y,k-D_z:k+D_z,soln))) / pp_norm
               fine%bnd_x(-fine%js+2*j+D_y,-fine%ks+2*k+D_z,LO) =sum(pp(-D_y:D_y, -D_z:D_z, 2, 2) * (coarse%mgvar(coarse%is,j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%is-1,j-D_y:j+D_y,k-D_z:k+D_z,soln))) / pp_norm
               fine%bnd_x(-fine%js+2*j,    -fine%ks+2*k,    HI)=sum(pp(-D_y:D_y, -D_z:D_z, 1, 1) * (coarse%mgvar(coarse%ie,j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%ie+1,j-D_y:j+D_y,k-D_z:k+D_z,soln))) / pp_norm
               fine%bnd_x(-fine%js+2*j+D_y,-fine%ks+2*k,    HI)=sum(pp(-D_y:D_y, -D_z:D_z, 2, 1) * (coarse%mgvar(coarse%ie,j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%ie+1,j-D_y:j+D_y,k-D_z:k+D_z,soln))) / pp_norm
               fine%bnd_x(-fine%js+2*j,    -fine%ks+2*k+D_z,HI)=sum(pp(-D_y:D_y, -D_z:D_z, 1, 2) * (coarse%mgvar(coarse%ie,j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%ie+1,j-D_y:j+D_y,k-D_z:k+D_z,soln))) / pp_norm
               fine%bnd_x(-fine%js+2*j+D_y,-fine%ks+2*k+D_z,HI)=sum(pp(-D_y:D_y, -D_z:D_z, 2, 2) * (coarse%mgvar(coarse%ie,j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%ie+1,j-D_y:j+D_y,k-D_z:k+D_z,soln))) / pp_norm
            enddo
         enddo
      endif

      if (has_dir(ydim)) then
         pp_norm = 2.*sum(pp(-D_x:D_x, -D_z:D_z, 1, 1))
         do i = coarse%is, coarse%ie
            do k = coarse%ks, coarse%ke
               fine%bnd_y(-fine%is+2*i,    -fine%ks+2*k,    LO) =sum(pp(-D_x:D_x, -D_z:D_z, 1, 1) * (coarse%mgvar(i-D_x:i+D_x,coarse%js,k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%js-1,k-D_z:k+D_z,soln))) / pp_norm
               fine%bnd_y(-fine%is+2*i+D_x,-fine%ks+2*k,    LO) =sum(pp(-D_x:D_x, -D_z:D_z, 2, 1) * (coarse%mgvar(i-D_x:i+D_x,coarse%js,k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%js-1,k-D_z:k+D_z,soln))) / pp_norm
               fine%bnd_y(-fine%is+2*i,    -fine%ks+2*k+D_z,LO) =sum(pp(-D_x:D_x, -D_z:D_z, 1, 2) * (coarse%mgvar(i-D_x:i+D_x,coarse%js,k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%js-1,k-D_z:k+D_z,soln))) / pp_norm
               fine%bnd_y(-fine%is+2*i+D_x,-fine%ks+2*k+D_z,LO) =sum(pp(-D_x:D_x, -D_z:D_z, 2, 2) * (coarse%mgvar(i-D_x:i+D_x,coarse%js,k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%js-1,k-D_z:k+D_z,soln))) / pp_norm
               fine%bnd_y(-fine%is+2*i,    -fine%ks+2*k,    HI)=sum(pp(-D_x:D_x, -D_z:D_z, 1, 1) * (coarse%mgvar(i-D_x:i+D_x,coarse%je,k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%je+1,k-D_z:k+D_z,soln))) / pp_norm
               fine%bnd_y(-fine%is+2*i+D_x,-fine%ks+2*k,    HI)=sum(pp(-D_x:D_x, -D_z:D_z, 2, 1) * (coarse%mgvar(i-D_x:i+D_x,coarse%je,k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%je+1,k-D_z:k+D_z,soln))) / pp_norm
               fine%bnd_y(-fine%is+2*i,    -fine%ks+2*k+D_z,HI)=sum(pp(-D_x:D_x, -D_z:D_z, 1, 2) * (coarse%mgvar(i-D_x:i+D_x,coarse%je,k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%je+1,k-D_z:k+D_z,soln))) / pp_norm
               fine%bnd_y(-fine%is+2*i+D_x,-fine%ks+2*k+D_z,HI)=sum(pp(-D_x:D_x, -D_z:D_z, 2, 2) * (coarse%mgvar(i-D_x:i+D_x,coarse%je,k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%je+1,k-D_z:k+D_z,soln))) / pp_norm
            enddo
         enddo
      endif

      if (has_dir(zdim)) then
         pp_norm = 2.*sum(pp(-D_x:D_x, -D_y:D_y, 1, 1))
         do i = coarse%is, coarse%ie
            do j = coarse%js, coarse%je
               fine%bnd_z(-fine%is+2*i,    -fine%js+2*j,    LO) =sum(pp(-D_x:D_x, -D_y:D_y, 1, 1) * (coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks,soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks-1,soln))) / pp_norm
               fine%bnd_z(-fine%is+2*i+D_x,-fine%js+2*j,    LO) =sum(pp(-D_x:D_x, -D_y:D_y, 2, 1) * (coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks,soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks-1,soln))) / pp_norm
               fine%bnd_z(-fine%is+2*i,    -fine%js+2*j+D_y,LO) =sum(pp(-D_x:D_x, -D_y:D_y, 1, 2) * (coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks,soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks-1,soln))) / pp_norm
               fine%bnd_z(-fine%is+2*i+D_x,-fine%js+2*j+D_y,LO) =sum(pp(-D_x:D_x, -D_y:D_y, 2, 2) * (coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks,soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks-1,soln))) / pp_norm
               fine%bnd_z(-fine%is+2*i,    -fine%js+2*j,    HI)=sum(pp(-D_x:D_x, -D_y:D_y, 1, 1) * (coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke,soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke+1,soln))) / pp_norm
               fine%bnd_z(-fine%is+2*i+D_x,-fine%js+2*j,    HI)=sum(pp(-D_x:D_x, -D_y:D_y, 2, 1) * (coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke,soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke+1,soln))) / pp_norm
               fine%bnd_z(-fine%is+2*i,    -fine%js+2*j+D_y,HI)=sum(pp(-D_x:D_x, -D_y:D_y, 1, 2) * (coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke,soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke+1,soln))) / pp_norm
               fine%bnd_z(-fine%is+2*i+D_x,-fine%js+2*j+D_y,HI)=sum(pp(-D_x:D_x, -D_y:D_y, 2, 2) * (coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke,soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke+1,soln))) / pp_norm
            enddo
         enddo
      endif

   end subroutine prolong_faces

end module multigridbasefuncs
