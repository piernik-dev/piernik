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

module multigridbasefuncs

   implicit none

contains

!!$ ============================================================================
!!
!! Clear boundary values
!!

   subroutine zero_boundaries(lev)

      use multigridvars, only : lvl

      implicit none

      integer, intent(in) :: lev  !< level for which clear the boundary values

      lvl(lev)%bnd_x(:,:,:) = 0.
      lvl(lev)%bnd_y(:,:,:) = 0.
      lvl(lev)%bnd_z(:,:,:) = 0.

   end subroutine zero_boundaries

!!$ ============================================================================
!!
!! Multigrid elementary operators: prolongation, restriction, residual, norm etc.
!!

   subroutine prolong_level(lev, iv)

      use errh,                  only: die
      use multigridvars,         only: plvl, lvl, level_min, level_max, ord_prolong, ngridvars
      use multigridhelpers,      only: dirty_debug, check_dirty, dirtyH
!      use multigridexperimental, only: prolong_level_hord

      implicit none

      integer, intent(in)      :: lev   !< level to prolong from
      integer, intent(in)      :: iv    !< variable to be prolonged

      type(plvl), pointer :: coarse, fine

      interface
         subroutine prolong_level_hord(lev, iv)
            implicit none
            integer, intent(in)      :: lev   !< level to prolong from
            integer, intent(in)      :: iv    !< variable to be prolonged
         end subroutine prolong_level_hord
      end interface

      if (lev >= level_max) return ! can't prolong finest level
      if (lev <  level_min) call die("[multigridbasefuncs:prolong_level] level <= 0.")
      if (iv < 1 .or. iv > ngridvars) call die("[multigridbasefuncs:prolong_level] Invalid variable index.")

      coarse => lvl(lev)
      fine   => lvl(lev + 1)

      if (dirty_debug) fine%mgvar(:, :, :, iv) = dirtyH

      if (ord_prolong == 0) then
         call prolong_level0(lev, iv)
      else
         call prolong_level_hord(lev, iv) ! experimental part
      end if

      call check_dirty(fine%level, iv, "prolong")

   end subroutine prolong_level

!!$ ============================================================================
!!
!! Oth order prolongation : injection
!!

   subroutine prolong_level0(lev, iv)

      use multigridvars, only: plvl, lvl, eff_dim, NDIM, XDIR, YDIR, ZDIR, has_dir, D_x, D_y, D_z
      use errh,          only: die

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
!!
!! Restriction operators
!!

   subroutine restrict_all(iv)

      use errh,               only: die
      use multigridhelpers,   only: check_dirty
      use multigridvars,      only: level_min, level_max, ngridvars

      implicit none

      integer, intent(in)      :: iv    !< variable to be restricted

      integer :: i

      if (iv < 1 .or. iv > ngridvars) call die("[multigridbasefuncs:restrict_all] Invalid variable index.")

      call check_dirty(level_max, iv, "restrict_all-")

      do i = level_max, level_min+1, -1
         call restrict_level(i, iv)
      end do

      call check_dirty(0, iv, "restrict_all+")

   end subroutine restrict_all

!!$ ============================================================================
!!
!! Simplest restriction (averaging).
!! \todo implement high order restriction and test its influence on V-cycle convergence rate
!!

   subroutine restrict_level(lev, iv)

      use errh,               only: die
      use multigridhelpers,   only: check_dirty
      use multigridvars,      only: plvl, lvl, level_min, level_max, ngridvars, eff_dim, NDIM, XDIR, YDIR, ZDIR, has_dir, D_x, D_y, D_z

      implicit none

      integer, intent(in)      :: lev   !< level to restrict from
      integer, intent(in)      :: iv    !< variable to be restricted

      type(plvl), pointer :: coarse, fine

      if (lev <= level_min) return ! can't restrict base level
      if (lev >  level_max) call die("[multigridbasefuncs:restrict_level] level>max.")
      if (iv < 1 .or. iv > ngridvars) call die("[multigridbasefuncs:restrict_level] Invalid variable index.")

      coarse => lvl(lev-1)
      fine   => lvl(lev)

      call check_dirty(fine%level, iv, "restrict_level-")

      ! BEWARE: unoptimized: some cells are used multiple times (1D and 2D speed-ups possible)
      coarse%mgvar(      coarse%is:coarse%ie,   coarse%js:coarse%je,   coarse%ks:coarse%ke,   iv) = &
           ( fine%mgvar( fine%is    :fine%ie-D_x:(1+D_x), fine%js    :fine%je-D_y:(1+D_y), fine%ks    :fine%ke-D_z:(1+D_z), iv) + &
           & fine%mgvar( fine%is+D_x:fine%ie    :(1+D_x), fine%js    :fine%je-D_y:(1+D_y), fine%ks    :fine%ke-D_z:(1+D_z), iv) + &
           & fine%mgvar( fine%is    :fine%ie-D_x:(1+D_x), fine%js+D_y:fine%je    :(1+D_y), fine%ks    :fine%ke-D_z:(1+D_z), iv) + &
           & fine%mgvar( fine%is+D_x:fine%ie    :(1+D_x), fine%js+D_y:fine%je    :(1+D_y), fine%ks    :fine%ke-D_z:(1+D_z), iv) + &
           & fine%mgvar( fine%is    :fine%ie-D_x:(1+D_x), fine%js    :fine%je-D_y:(1+D_y), fine%ks+D_z:fine%ke    :(1+D_z), iv) + &
           & fine%mgvar( fine%is+D_x:fine%ie    :(1+D_x), fine%js    :fine%je-D_y:(1+D_y), fine%ks+D_z:fine%ke    :(1+D_z), iv) + &
           & fine%mgvar( fine%is    :fine%ie-D_x:(1+D_x), fine%js+D_y:fine%je    :(1+D_y), fine%ks+D_z:fine%ke    :(1+D_z), iv) + &
           & fine%mgvar( fine%is+D_x:fine%ie    :(1+D_x), fine%js+D_y:fine%je    :(1+D_y), fine%ks+D_z:fine%ke    :(1+D_z), iv) ) *0.125  !/ ((1.+D_x)*(1.+D_y)*(1.+D_z))

      call check_dirty(coarse%level, iv, "restrict_level+")

   end subroutine restrict_level

!!$ ============================================================================
!!
!! Calculate the residuum for the Poisson equation.
!!

   subroutine residual(lev, src, soln, def)

      use errh,                  only: die
      use multigridvars,         only: ord_laplacian, level_min, level_max, ngridvars
!      use multigridexperimental, only: residual4

      implicit none

      integer, intent(in) :: lev  !< level for which approximate the solution
      integer, intent(in) :: src  !< index of source in lvl()%mgvar
      integer, intent(in) :: soln !< index of solution in lvl()%mgvar
      integer, intent(in) :: def  !< index of defect in lvl()%mgvar

      interface
         subroutine residual4(lev, src, soln, def)
            implicit none
            integer, intent(in) :: lev  !< level for which approximate the solution
            integer, intent(in) :: src  !< index of source in lvl()%mgvar
            integer, intent(in) :: soln !< index of solution in lvl()%mgvar
            integer, intent(in) :: def  !< index of defect in lvl()%mgvar
         end subroutine residual4
      end interface

      if (any( [ src, soln, def ] <= 0) .or. any( [ src, soln, def ] > ngridvars)) call die("[multigridbasefuncs:residual] Invalid variable index")
      if (lev < level_min .or. lev > level_max) call die("[multigridbasefuncs:residual] Invalid level number")

      select case (ord_laplacian)
      case (2)
         call residual2(lev, src, soln, def)
      case (4)
         call residual4(lev, src, soln, def)
      case default
         call die("[multigridbasefuncs:residual] The parameter 'ord_laplacian' must be 2 or 4")
      end select

   end subroutine residual

!!$ ============================================================================
!!
!! 2nd order Laplacian
!!

   subroutine residual2(lev, src, soln, def)

      use multigridvars,      only: lvl, eff_dim, NDIM, XDIR, YDIR, ZDIR, has_dir
      use multigridmpifuncs,  only: mpi_multigrid_bnd
      use multigridhelpers,   only: multidim_code_3D

      implicit none

      integer, intent(in) :: lev  !< level for which approximate the solution
      integer, intent(in) :: src  !< index of source in lvl()%mgvar
      integer, intent(in) :: soln !< index of solution in lvl()%mgvar
      integer, intent(in) :: def  !< index of defect in lvl()%mgvar

      real                :: L0, Lx, Ly, Lz
      integer :: k

      call mpi_multigrid_bnd(lev, soln, 1, .false.) ! no corners required

      ! Coefficients for a simplest 3-point Laplacian operator: [ 1, -2, 1 ]
      ! for 2D and 1D setups appropriate elements of [ Lx, Ly, Lz ] should be == 0.
      Lx = lvl(lev)%idx2
      Ly = lvl(lev)%idy2
      Lz = lvl(lev)%idz2
      L0 = -2. * (Lx + Ly + Lz)

      ! Possible optimization candidate: reduce cache misses (secondary importance, cache-aware implementation required)
      ! Explicit loop over k gives here better performance than array operation due to less cache misses (at least on 32^3 and 64^3 arrays)
      if (eff_dim == NDIM .and. .not. multidim_code_3D) then
         do k = lvl(lev)%ks, lvl(lev)%ke
            lvl(       lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   k,   def)        = &
                 & lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   k,   src)        - &
                 ( lvl(lev)%mgvar(lvl(lev)%is-1:lvl(lev)%ie-1, lvl(lev)%js  :lvl(lev)%je,   k,   soln)       + &
                 & lvl(lev)%mgvar(lvl(lev)%is+1:lvl(lev)%ie+1, lvl(lev)%js  :lvl(lev)%je,   k,   soln)) * Lx - &
                 ( lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js-1:lvl(lev)%je-1, k,   soln)       + &
                 & lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js+1:lvl(lev)%je+1, k,   soln)) * Ly - &
                 ( lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   k-1, soln)       + &
                 & lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   k+1, soln)) * Lz - &
                 & lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   k,   soln)  * L0
         end do
      else
         ! In 3D this implementation can give a bit more cache misses, few times more writes and significantly more instructions executed than monolithic 3D above
         do k = lvl(lev)%ks, lvl(lev)%ke
            lvl(       lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   k,   def)   = &
                 & lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   k,   src)   - &
                 & lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   k,   soln)  * L0
            if (has_dir(XDIR)) &
                 & lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   k,   def)   = &
                 & lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   k,   def)   - &
                 ( lvl(lev)%mgvar(lvl(lev)%is-1:lvl(lev)%ie-1, lvl(lev)%js  :lvl(lev)%je,   k,   soln)  + &
                 & lvl(lev)%mgvar(lvl(lev)%is+1:lvl(lev)%ie+1, lvl(lev)%js  :lvl(lev)%je,   k,   soln)) * Lx
            if (has_dir(YDIR)) &
                 & lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   k,   def)   = &
                 & lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   k,   def)   - &
                 ( lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js-1:lvl(lev)%je-1, k,   soln)  + &
                 & lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js+1:lvl(lev)%je+1, k,   soln)) * Ly
            if (has_dir(ZDIR)) &
                 & lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   k,   def)   = &
                 & lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   k,   def)   - &
                 ( lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   k-1, soln)  + &
                 & lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   k+1, soln)) * Lz
         end do
      end if

   end subroutine residual2

!!$ ============================================================================
!!
!! Calculate L2 norm
!!

   subroutine norm_sq(iv, norm)

      use errh,          only: die
      use mpisetup,      only: comm3d, ierr, MPI_DOUBLE_PRECISION, MPI_SUM
      use multigridvars, only: ngridvars, roof

      implicit none

      integer, intent(in)  :: iv   !< index of variable in lvl()%mgvar for which we want to find the norm
      real,    intent(out) :: norm !< the calculated norm

      real                 :: lsum

      if (iv <= 0 .or. iv > ngridvars) call die("[multigridbasefuncs:norm_sq] Invalid variable index")

      lsum = sum(roof%mgvar(roof%is:roof%ie, roof%js:roof%je, roof%ks:roof%ke, iv)**2) * roof%dvol
      call MPI_Allreduce(lsum, norm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm3d, ierr)
      norm = sqrt(norm)

   end subroutine norm_sq

!!$ ============================================================================
!!
!! Compute the global average value and substract it from the whole domain
!!

   subroutine substract_average(lev, iv)

      use errh,          only: die
      use mpisetup,      only: comm3d, ierr, MPI_DOUBLE_PRECISION, MPI_SUM
      use multigridvars, only: lvl, level_min, level_max, ngridvars

      implicit none

      integer, intent(in) :: lev  !< level for which we want to substract its average from
      integer, intent(in) :: iv   !< index of variable in lvl()%mgvar which we want to have zero average

      real                :: lsum, avg, vol

      if (lev < level_min .or. lev > level_max) call die("[multigridbasefuncs:substract_average] Invalid level number.")
      if (iv < 1 .or. iv > ngridvars) call die("[multigridbasefuncs:substract_average] Invalid variable index.")

      lsum = sum(lvl(lev)%mgvar(lvl(lev)%is:lvl(lev)%ie, lvl(lev)%js:lvl(lev)%je, lvl(lev)%ks:lvl(lev)%ke, iv)) * lvl(lev)%dvol
      call MPI_Allreduce(lsum, avg, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm3d, ierr)
      call MPI_Allreduce(lvl(lev)%vol, vol, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm3d, ierr) ! This probably can be calculated in init_multigrid
      avg = avg / vol

      lvl(lev)%mgvar(:, :, :, iv) = lvl(lev)%mgvar(:, :, :, iv) - avg

   end subroutine substract_average

end module multigridbasefuncs
