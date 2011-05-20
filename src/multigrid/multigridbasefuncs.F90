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
      use multigridvars,         only: plvl, lvl, base, roof, ord_prolong, ngridvars, is_mg_uneven
      use multigridexperimental, only: prolong_level_hord

      implicit none

      integer, intent(in)      :: lev   !< level to prolong from
      integer, intent(in)      :: iv    !< variable to be prolonged

      type(plvl), pointer :: coarse, fine

      if (lev >= roof%level) return ! can't prolong finest level
      if (lev <  base%level) call die("[multigridbasefuncs:prolong_level] level < base%level.")
      if (iv < 1 .or. iv > ngridvars) call die("[multigridbasefuncs:prolong_level] Invalid variable index.")

      coarse => lvl(lev)
      if (.not. associated(coarse)) call die("[multigridbasefuncs:prolong_level] coarse == null()")
      fine   => coarse%finer
      if (.not. associated(fine)) call die("[multigridbasefuncs:prolong_level] fine == null()")

      if (dirty_debug) fine%mgvar(:, :, :, iv) = dirtyH

      call check_dirty(coarse%level, iv, "prolong-")

      if (ord_prolong == 0 .or. is_mg_uneven) then
         call coarse%prolong_level0(iv)
      else
         call prolong_level_hord(lev, iv) ! experimental part
      endif

      call check_dirty(fine%level, iv, "prolong+")

   end subroutine prolong_level


!!$ ============================================================================
!>
!! \brief Restriction operators
!<

   subroutine restrict_all(iv)

      use dataio_pub,         only: die
      use multigridhelpers,   only: check_dirty
      use multigridvars,      only: roof, base, ngridvars, lvl, plvl

      implicit none

      integer, intent(in)      :: iv    !< variable to be restricted

      type(plvl), pointer :: curl

      if (iv < 1 .or. iv > ngridvars) call die("[multigridbasefuncs:restrict_all] Invalid variable index.")

      call check_dirty(roof%level, iv, "restrict_all-")

      curl => roof
      do while (associated(curl) .and. .not. associated(curl, base))
         call curl%restrict_level(iv)
         curl => curl%coarser
      enddo

      call check_dirty(base%level, iv, "restrict_all+")

   end subroutine restrict_all

!!$ ============================================================================
!>
!! \brief Calculate L2 norm
!<

   subroutine norm_sq(iv, norm)

      use dataio_pub,    only: die
      use mpisetup,      only: comm, ierr, geometry_type
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
      call MPI_Allreduce(lsum, norm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
      norm = sqrt(norm)

   end subroutine norm_sq

!!$ ============================================================================
!>
!! \brief Compute the global average value and subtract it from the whole domain
!<

   subroutine substract_average(lev, iv)

      use dataio_pub,    only: die
      use mpisetup,      only: comm, ierr, geometry_type
      use mpi,           only: MPI_DOUBLE_PRECISION, MPI_SUM
      use multigridvars, only: lvl, plvl, base, roof, ngridvars
      use constants,     only: GEO_XYZ, GEO_RPZ

      implicit none

      integer, intent(in) :: lev  !< level for which we want to subtract its average from
      integer, intent(in) :: iv   !< index of variable in lvl()%mgvar which we want to have zero average

      real                :: lsum, avg, vol
      integer             :: i
      type(plvl), pointer :: curl

      if (lev < base%level .or. lev > roof%level) call die("[multigridbasefuncs:substract_average] Invalid level number.")
      if (iv < 1 .or. iv > ngridvars) call die("[multigridbasefuncs:substract_average] Invalid variable index.")

      curl => lvl(lev)

      select case (geometry_type)
         case (GEO_XYZ)
            lsum = sum(curl%mgvar(curl%is:curl%ie, curl%js:curl%je, curl%ks:curl%ke, iv)) * curl%dvol
         case (GEO_RPZ)
            lsum = 0.
            do i = curl%is, curl%ie
               lsum = lsum + sum(curl%mgvar(i, curl%js:curl%je, curl%ks:curl%ke, iv)) * curl%dvol * curl%x(i)
            enddo
         case default
            call die("[multigridbasefuncs:substract_average] Unsupported geometry.")
      end select
      call MPI_Allreduce(lsum, avg, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
      call MPI_Allreduce(curl%vol, vol, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr) ! This probably can be calculated in init_multigrid
      avg = avg / vol

      curl%mgvar(:, :, :, iv) = curl%mgvar(:, :, :, iv) - avg

   end subroutine substract_average

!!$ ============================================================================
!>
!! \brief Prolong solution data at level (lev-1) to faces at level lev
!!
!! \details It looks that parallel prolongation order 0 is the best from range -2 .. 2. Quite good is also order 2.(tests were performed with the maclaurin problem at r4051).
!! There were no experiments with higher order prolongation, but it seems there is not much room for any improvements.
!! The tests were performed on uniform grid, where the interpolation affects only convergence factors.
!! On a refinement step the interpolation type affects also sulution by influencing the way how fine and coarse grids are coupled.
!!
!! \todo move to plvl or multigrid_gravity
!!
!! OPT: completely unoptimized
!<

   subroutine prolong_faces(lev, soln)

      use grid,               only: D_x, D_y, D_z
      use mpisetup,           only: has_dir, comm
      use constants,          only: xdim, ydim, zdim, LO, HI, LONG
      use dataio_pub,         only: die, warn
      use multigridhelpers,   only: check_dirty
      use multigridmpifuncs,  only: mpi_multigrid_bnd
      use multigridvars,      only: plvl, lvl, ord_prolong_face_norm, ord_prolong_face_par, base, roof, extbnd_antimirror, need_general_pf

      implicit none

      integer, intent(in) :: lev  !< level for which approximate the solution
      integer, intent(in) :: soln !< index of solution in lvl()%mgvar

      integer                       :: i, j, k
      type(plvl), pointer           :: coarse, fine
      integer, parameter            :: s_wdth  = 3           ! interpolation stencil width
      integer, parameter            :: s_rng = (s_wdth-1)/2  ! stencil range around 0
      real, parameter, dimension(s_wdth) :: p0  = [ 0.,       1.,     0.     ] ! injection
      real, parameter, dimension(s_wdth) :: p1  = [ 0.,       3./4.,  1./4.  ] ! 1D linear prolongation stencil
      real, parameter, dimension(s_wdth) :: p2i = [ -1./8.,   1.,     1./8.  ] ! 1D integral cubic prolongation stencil
      real, parameter, dimension(s_wdth) :: p2d = [ -3./32., 30./32., 5./32. ] ! 1D direct cubic prolongation stencil
      real, dimension(-s_rng:s_rng)                    :: p
      real, dimension(-s_rng:s_rng, -s_rng:s_rng, 2, 2):: pp   ! 2D prolongation stencil
      real                          :: pp_norm
      real :: opfn1, opfn3
      integer :: b_rng

      ! coefficients for intepolation perpendicular to the face can be expressed as
      ! 1/2*| 1 0 ... ] + a2*| 1 -1 0 ... ] + a4*| 1 -2 1 0 ... ] + a6*| 1 -3 3 -1 0 ... ] + a8*| 1 -4 6 -4 1 0 ... ] + ...
      ! 1./2.   * ( 0   0   1   1   0 0 ) ! average, direct linear (a2 =  0.,      ...)
      ! 1./16.  * ( 0  -1   9   9  -1 0 ) ! direct cubic           (a2 =  1./16.,  a4 = 0.,      ...)
      ! 1./256. * ( 3 -25 150 150 -25 3 ) ! direct quintic         (a2 = 19./256., a4 = 3./256., a6 = 0., ...)

      ! factor 2./8. (a2 = 1./8.) was found experimentally
      opfn1 = 1. + 2./8. ! an additional factor equal to 0.5 is applied later as normalization of  pp(:,:,:,:)
      opfn3 = -2./8.
      ! the optimal value of a2 seems to depend at least on ord_prolong_face_par
      ! for ord_prolong_face_par =  0 it is ~1.7 (best convergence, reduction form 8/9 cycles to 7 cycles in maclaurin test for norm_tol=1e-10 and nsmool=8)
      ! for ord_prolong_face_par =  2 it is ~2,2 (slightly worse convergence, 8 cycles instead of 9/10 cycles)
      ! for ord_prolong_face_par = -2 it is ~2.5 (convergence in 8 cycles instead of 10 cycles)
      ! for ord_prolong_face_par =  1 it is ~3.6 (worst convergence, 11 -> 9 cycles)

      if (lev < base%level .or. lev > roof%level) call die("[multigridbasefuncs:prolong_faces] Invalid level")

      if (lev == base%level) then
         call warn("[multigridbasefuncs:prolong_faces] Cannot prolong anything to base level")
         return
      endif

      fine => lvl(lev)
      if (.not. associated(fine)) call die("[multigridbasefuncs:prolong_faces] fine == null()")
      coarse => fine%coarser
      if (.not. associated(coarse)) call die("[multigridbasefuncs:prolong_faces] coarse == null()")

      if (need_general_pf) then
         call die("[multigridbasefuncs:prolong_faces] comm3d == MPI_COMM_NULL not implemented yet")
      else
         select case (ord_prolong_face_par)
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

         do i = -s_rng, s_rng
            pp(i,:,1,1) = 0.5*p( i)*p(:)       ! 0.5 because of face averaging
            pp(i,:,1,2) = 0.5*p( i)*p(1:-1:-1) ! or use matmul()
            pp(i,:,2,1) = 0.5*p(-i)*p(:)
            pp(i,:,2,2) = 0.5*p(-i)*p(1:-1:-1)
         enddo

         b_rng = s_rng
         if (ord_prolong_face_norm > 0) b_rng = max(b_rng, ord_prolong_face_norm+1)
         call mpi_multigrid_bnd(lev-1, soln, b_rng, extbnd_antimirror, corners=(ord_prolong_face_par/=0)) !> \deprecated BEWARE for higher prolongation order more guardcell are required
         call check_dirty(lev-1, soln, "prolong_faces", s_rng)

         if (ord_prolong_face_norm > 0) then
            if (has_dir(xdim)) then
               pp_norm = 2.*sum(pp(-D_y:D_y, -D_z:D_z, 1, 1)) ! normalization is required for ord_prolong_face_par == 1 and -2
               do j = coarse%js, coarse%je
                  do k = coarse%ks, coarse%ke
                     fine%bnd_x(-fine%js+2*j,    -fine%ks+2*k,    LO)=sum(pp(-D_y:D_y, -D_z:D_z, 1, 1) * ( &
                          opfn1*(coarse%mgvar(coarse%is,  j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%is-1,j-D_y:j+D_y,k-D_z:k+D_z,soln)) + &
                          opfn3*(coarse%mgvar(coarse%is+1,j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%is-2,j-D_y:j+D_y,k-D_z:k+D_z,soln)))) / pp_norm
                     fine%bnd_x(-fine%js+2*j+D_y,-fine%ks+2*k,    LO)=sum(pp(-D_y:D_y, -D_z:D_z, 2, 1) * ( &
                          opfn1*(coarse%mgvar(coarse%is,  j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%is-1,j-D_y:j+D_y,k-D_z:k+D_z,soln)) + &
                          opfn3*(coarse%mgvar(coarse%is+1,j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%is-2,j-D_y:j+D_y,k-D_z:k+D_z,soln)))) / pp_norm
                     fine%bnd_x(-fine%js+2*j,    -fine%ks+2*k+D_z,LO)=sum(pp(-D_y:D_y, -D_z:D_z, 1, 2) * ( &
                          opfn1*(coarse%mgvar(coarse%is  ,j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%is-1,j-D_y:j+D_y,k-D_z:k+D_z,soln)) + &
                          opfn3*(coarse%mgvar(coarse%is+1,j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%is-2,j-D_y:j+D_y,k-D_z:k+D_z,soln)))) / pp_norm
                     fine%bnd_x(-fine%js+2*j+D_y,-fine%ks+2*k+D_z,LO)=sum(pp(-D_y:D_y, -D_z:D_z, 2, 2) * ( &
                          opfn1*(coarse%mgvar(coarse%is,  j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%is-1,j-D_y:j+D_y,k-D_z:k+D_z,soln)) + &
                          opfn3*(coarse%mgvar(coarse%is+1,j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%is-2,j-D_y:j+D_y,k-D_z:k+D_z,soln)))) / pp_norm
                     fine%bnd_x(-fine%js+2*j,    -fine%ks+2*k,    HI)=sum(pp(-D_y:D_y, -D_z:D_z, 1, 1) * ( &
                          opfn1*(coarse%mgvar(coarse%ie,  j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%ie+1,j-D_y:j+D_y,k-D_z:k+D_z,soln)) + &
                          opfn3*(coarse%mgvar(coarse%ie-1,j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%ie+2,j-D_y:j+D_y,k-D_z:k+D_z,soln)))) / pp_norm
                     fine%bnd_x(-fine%js+2*j+D_y,-fine%ks+2*k,    HI)=sum(pp(-D_y:D_y, -D_z:D_z, 2, 1) * ( &
                          opfn1*(coarse%mgvar(coarse%ie,  j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%ie+1,j-D_y:j+D_y,k-D_z:k+D_z,soln)) + &
                          opfn3*(coarse%mgvar(coarse%ie-1,j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%ie+2,j-D_y:j+D_y,k-D_z:k+D_z,soln)))) / pp_norm
                     fine%bnd_x(-fine%js+2*j,    -fine%ks+2*k+D_z,HI)=sum(pp(-D_y:D_y, -D_z:D_z, 1, 2) * ( &
                          opfn1*(coarse%mgvar(coarse%ie,  j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%ie+1,j-D_y:j+D_y,k-D_z:k+D_z,soln)) + &
                          opfn3*(coarse%mgvar(coarse%ie-1,j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%ie+2,j-D_y:j+D_y,k-D_z:k+D_z,soln)))) / pp_norm
                     fine%bnd_x(-fine%js+2*j+D_y,-fine%ks+2*k+D_z,HI)=sum(pp(-D_y:D_y, -D_z:D_z, 2, 2) * ( &
                          opfn1*(coarse%mgvar(coarse%ie,  j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%ie+1,j-D_y:j+D_y,k-D_z:k+D_z,soln)) + &
                          opfn3*(coarse%mgvar(coarse%ie-1,j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%ie+2,j-D_y:j+D_y,k-D_z:k+D_z,soln)))) / pp_norm
                  enddo
               enddo
            endif

            if (has_dir(ydim)) then
               pp_norm = 2.*sum(pp(-D_x:D_x, -D_z:D_z, 1, 1))
               do i = coarse%is, coarse%ie
                  do k = coarse%ks, coarse%ke
                     fine%bnd_y(-fine%is+2*i,    -fine%ks+2*k,    LO)=sum(pp(-D_x:D_x, -D_z:D_z, 1, 1) * ( &
                          opfn1*(coarse%mgvar(i-D_x:i+D_x,coarse%js,  k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%js-1,k-D_z:k+D_z,soln)) + &
                          opfn3*(coarse%mgvar(i-D_x:i+D_x,coarse%js+1,k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%js-2,k-D_z:k+D_z,soln)))) / pp_norm
                     fine%bnd_y(-fine%is+2*i+D_x,-fine%ks+2*k,    LO)=sum(pp(-D_x:D_x, -D_z:D_z, 2, 1) * ( &
                          opfn1*(coarse%mgvar(i-D_x:i+D_x,coarse%js,  k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%js-1,k-D_z:k+D_z,soln)) + &
                          opfn3*(coarse%mgvar(i-D_x:i+D_x,coarse%js+1,k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%js-2,k-D_z:k+D_z,soln)))) / pp_norm
                     fine%bnd_y(-fine%is+2*i,    -fine%ks+2*k+D_z,LO)=sum(pp(-D_x:D_x, -D_z:D_z, 1, 2) * ( &
                          opfn1*(coarse%mgvar(i-D_x:i+D_x,coarse%js,  k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%js-1,k-D_z:k+D_z,soln)) + &
                          opfn3*(coarse%mgvar(i-D_x:i+D_x,coarse%js+1,k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%js-2,k-D_z:k+D_z,soln)))) / pp_norm
                     fine%bnd_y(-fine%is+2*i+D_x,-fine%ks+2*k+D_z,LO)=sum(pp(-D_x:D_x, -D_z:D_z, 2, 2) * ( &
                          opfn1*(coarse%mgvar(i-D_x:i+D_x,coarse%js,  k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%js-1,k-D_z:k+D_z,soln)) + &
                          opfn3*(coarse%mgvar(i-D_x:i+D_x,coarse%js+1,k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%js-2,k-D_z:k+D_z,soln)))) / pp_norm
                     fine%bnd_y(-fine%is+2*i,    -fine%ks+2*k,    HI)=sum(pp(-D_x:D_x, -D_z:D_z, 1, 1) * ( &
                          opfn1*(coarse%mgvar(i-D_x:i+D_x,coarse%je,  k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%je+1,k-D_z:k+D_z,soln)) + &
                          opfn3*(coarse%mgvar(i-D_x:i+D_x,coarse%je-1,k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%je+2,k-D_z:k+D_z,soln)))) / pp_norm
                     fine%bnd_y(-fine%is+2*i+D_x,-fine%ks+2*k,    HI)=sum(pp(-D_x:D_x, -D_z:D_z, 2, 1) * ( &
                          opfn1*(coarse%mgvar(i-D_x:i+D_x,coarse%je,  k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%je+1,k-D_z:k+D_z,soln)) + &
                          opfn3*(coarse%mgvar(i-D_x:i+D_x,coarse%je-1,k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%je+2,k-D_z:k+D_z,soln)))) / pp_norm
                     fine%bnd_y(-fine%is+2*i,    -fine%ks+2*k+D_z,HI)=sum(pp(-D_x:D_x, -D_z:D_z, 1, 2) * ( &
                          opfn1*(coarse%mgvar(i-D_x:i+D_x,coarse%je,  k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%je+1,k-D_z:k+D_z,soln)) + &
                          opfn3*(coarse%mgvar(i-D_x:i+D_x,coarse%je-1,k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%je+2,k-D_z:k+D_z,soln)))) / pp_norm
                     fine%bnd_y(-fine%is+2*i+D_x,-fine%ks+2*k+D_z,HI)=sum(pp(-D_x:D_x, -D_z:D_z, 2, 2) * ( &
                          opfn1*(coarse%mgvar(i-D_x:i+D_x,coarse%je,  k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%je+1,k-D_z:k+D_z,soln)) + &
                          opfn3*(coarse%mgvar(i-D_x:i+D_x,coarse%je-1,k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%je+2,k-D_z:k+D_z,soln)))) / pp_norm
                  enddo
               enddo
            endif

            if (has_dir(zdim)) then
               pp_norm = 2.*sum(pp(-D_x:D_x, -D_y:D_y, 1, 1))
               do i = coarse%is, coarse%ie
                  do j = coarse%js, coarse%je
                     fine%bnd_z(-fine%is+2*i,    -fine%js+2*j,    LO)=sum(pp(-D_x:D_x, -D_y:D_y, 1, 1) * ( &
                          opfn1*(coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks,  soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks-1,soln)) + &
                          opfn3*(coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks+1,soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks-2,soln)))) / pp_norm
                     fine%bnd_z(-fine%is+2*i+D_x,-fine%js+2*j,    LO)=sum(pp(-D_x:D_x, -D_y:D_y, 2, 1) * ( &
                          opfn1*(coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks,  soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks-1,soln)) + &
                          opfn3*(coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks+1,soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks-2,soln)))) / pp_norm
                     fine%bnd_z(-fine%is+2*i,    -fine%js+2*j+D_y,LO)=sum(pp(-D_x:D_x, -D_y:D_y, 1, 2) * ( &
                          opfn1*(coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks,  soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks-1,soln)) + &
                          opfn3*(coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks+1,soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks-2,soln)))) / pp_norm
                     fine%bnd_z(-fine%is+2*i+D_x,-fine%js+2*j+D_y,LO)=sum(pp(-D_x:D_x, -D_y:D_y, 2, 2) * ( &
                          opfn1*(coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks,  soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks-1,soln)) + &
                          opfn3*(coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks+1,soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks-2,soln)))) / pp_norm
                     fine%bnd_z(-fine%is+2*i,    -fine%js+2*j,    HI)=sum(pp(-D_x:D_x, -D_y:D_y, 1, 1) * ( &
                          opfn1*(coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke,  soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke+1,soln)) + &
                          opfn3*(coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke-1,soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke+2,soln)))) / pp_norm
                     fine%bnd_z(-fine%is+2*i+D_x,-fine%js+2*j,    HI)=sum(pp(-D_x:D_x, -D_y:D_y, 2, 1) * ( &
                          opfn1*(coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke,  soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke+1,soln)) + &
                          opfn3*(coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke-1,soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke+2,soln)))) / pp_norm
                     fine%bnd_z(-fine%is+2*i,    -fine%js+2*j+D_y,HI)=sum(pp(-D_x:D_x, -D_y:D_y, 1, 2) * ( &
                          opfn1*(coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke,  soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke+1,soln)) + &
                          opfn3*(coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke-1,soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke+2,soln)))) / pp_norm
                     fine%bnd_z(-fine%is+2*i+D_x,-fine%js+2*j+D_y,HI)=sum(pp(-D_x:D_x, -D_y:D_y, 2, 2) * ( &
                          opfn1*(coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke,  soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke+1,soln)) + &
                          opfn3*(coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke-1,soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke+2,soln)))) / pp_norm
                  enddo
               enddo
            endif

         else

            if (has_dir(xdim)) then
               pp_norm = 2.*sum(pp(-D_y:D_y, -D_z:D_z, 1, 1)) ! normalization is required for ord_prolong_face_par == 1 and -2
               do j = coarse%js, coarse%je
                  do k = coarse%ks, coarse%ke
                     fine%bnd_x(-fine%js+2*j,    -fine%ks+2*k,    LO)=sum(pp(-D_y:D_y, -D_z:D_z, 1, 1) * (coarse%mgvar(coarse%is,j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%is-1,j-D_y:j+D_y,k-D_z:k+D_z,soln))) / pp_norm
                     fine%bnd_x(-fine%js+2*j+D_y,-fine%ks+2*k,    LO)=sum(pp(-D_y:D_y, -D_z:D_z, 2, 1) * (coarse%mgvar(coarse%is,j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%is-1,j-D_y:j+D_y,k-D_z:k+D_z,soln))) / pp_norm
                     fine%bnd_x(-fine%js+2*j,    -fine%ks+2*k+D_z,LO)=sum(pp(-D_y:D_y, -D_z:D_z, 1, 2) * (coarse%mgvar(coarse%is,j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%is-1,j-D_y:j+D_y,k-D_z:k+D_z,soln))) / pp_norm
                     fine%bnd_x(-fine%js+2*j+D_y,-fine%ks+2*k+D_z,LO)=sum(pp(-D_y:D_y, -D_z:D_z, 2, 2) * (coarse%mgvar(coarse%is,j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%is-1,j-D_y:j+D_y,k-D_z:k+D_z,soln))) / pp_norm
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
                     fine%bnd_y(-fine%is+2*i,    -fine%ks+2*k,    LO)=sum(pp(-D_x:D_x, -D_z:D_z, 1, 1) * (coarse%mgvar(i-D_x:i+D_x,coarse%js,k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%js-1,k-D_z:k+D_z,soln))) / pp_norm
                     fine%bnd_y(-fine%is+2*i+D_x,-fine%ks+2*k,    LO)=sum(pp(-D_x:D_x, -D_z:D_z, 2, 1) * (coarse%mgvar(i-D_x:i+D_x,coarse%js,k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%js-1,k-D_z:k+D_z,soln))) / pp_norm
                     fine%bnd_y(-fine%is+2*i,    -fine%ks+2*k+D_z,LO)=sum(pp(-D_x:D_x, -D_z:D_z, 1, 2) * (coarse%mgvar(i-D_x:i+D_x,coarse%js,k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%js-1,k-D_z:k+D_z,soln))) / pp_norm
                     fine%bnd_y(-fine%is+2*i+D_x,-fine%ks+2*k+D_z,LO)=sum(pp(-D_x:D_x, -D_z:D_z, 2, 2) * (coarse%mgvar(i-D_x:i+D_x,coarse%js,k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%js-1,k-D_z:k+D_z,soln))) / pp_norm
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
                     fine%bnd_z(-fine%is+2*i,    -fine%js+2*j,    LO)=sum(pp(-D_x:D_x, -D_y:D_y, 1, 1) * (coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks,soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks-1,soln))) / pp_norm
                     fine%bnd_z(-fine%is+2*i+D_x,-fine%js+2*j,    LO)=sum(pp(-D_x:D_x, -D_y:D_y, 2, 1) * (coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks,soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks-1,soln))) / pp_norm
                     fine%bnd_z(-fine%is+2*i,    -fine%js+2*j+D_y,LO)=sum(pp(-D_x:D_x, -D_y:D_y, 1, 2) * (coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks,soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks-1,soln))) / pp_norm
                     fine%bnd_z(-fine%is+2*i+D_x,-fine%js+2*j+D_y,LO)=sum(pp(-D_x:D_x, -D_y:D_y, 2, 2) * (coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks,soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks-1,soln))) / pp_norm
                     fine%bnd_z(-fine%is+2*i,    -fine%js+2*j,    HI)=sum(pp(-D_x:D_x, -D_y:D_y, 1, 1) * (coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke,soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke+1,soln))) / pp_norm
                     fine%bnd_z(-fine%is+2*i+D_x,-fine%js+2*j,    HI)=sum(pp(-D_x:D_x, -D_y:D_y, 2, 1) * (coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke,soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke+1,soln))) / pp_norm
                     fine%bnd_z(-fine%is+2*i,    -fine%js+2*j+D_y,HI)=sum(pp(-D_x:D_x, -D_y:D_y, 1, 2) * (coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke,soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke+1,soln))) / pp_norm
                     fine%bnd_z(-fine%is+2*i+D_x,-fine%js+2*j+D_y,HI)=sum(pp(-D_x:D_x, -D_y:D_y, 2, 2) * (coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke,soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke+1,soln))) / pp_norm
                  enddo
               enddo
            endif
         endif
      endif

   end subroutine prolong_faces

end module multigridbasefuncs
