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

module multipole

   ! needed for global vars in this module
   use multigridvars, only : plvl, NDIM, LOW, HIGH

   implicit none

   integer, parameter        :: INSIDE = 1, OUTSIDE = INSIDE + 1 !< distinction between interior and exterior multipole expansion

   ! namelist parameters
   integer                   :: lmax                             !< Maximum l-order of multipole moments
   integer                   :: mmax                             !< Maximum m-order of multipole moments. Equal to lmax by default.
   integer                   :: ord_prolong_mpole                !< boundary prolongation operator order; allowed values are -2 .. 2
   integer                   :: coarsen_multipole                !< If > 0 then evaluate multipoles at level_max-coarsen_multipole level
   logical                   :: use_point_monopole               !< Don't evaluate multipole moments, use point-like mass approximation (crudest possible)

   ! radial discredization
   integer                   :: rqbin                            !< number of radial samples of multipoles
   real                      :: drq, rscale                      !< radial resolution of multipoles
   integer                   :: irmin, irmax                     !< minimum and maximum Q(:, :, r) indices in use

   type(plvl), pointer       :: lmpole                           !< pointer to the level where multipoles are evaluated
   real, dimension(0:NDIM)   :: CoM                              !< Total mass and center of mass coordinates

   ! boundaries
   real, dimension(LOW:HIGH) :: fbnd_x, fbnd_y, fbnd_z           !< coordinates at faces of local domain

   ! multipoles and auxiliary factors
   real, dimension(:,:,:),   allocatable :: Q                    !< The whole moment array with dependence on radius
   real, dimension(:,:,:),   allocatable :: k12                  !< array of Legendre recurrence factors
   real, dimension(:),       allocatable :: sfac, cfac           !< 0..m_max sized array of azimuthal sine and cosine factors
   real, dimension(:),       allocatable :: ofact                !< arrays of Legendre normalization factor (compressed)
   real, dimension(:),       allocatable :: rn, irn              !< 0..l_max sized array of positive and negative powers of r

contains

!!$ ============================================================================
!!
!! This function returns array index for "compressed" Q(:, inout, r) array replacing "plain" Q(l, m, inout, r) array
!! Originally there were separate indices for l and m multipole numbers in Q and ofact arrays. Since there are no Y_{l,m} harmonics for l<|m|, approximately half of the Q array
!! was left unused. For each m there is only lmax-|m|+1 valid Y_{l,m} harmonic contributions (for l = |m| .. lmax). The function lm(l, m) converts each valid (l,m) pair into an
!! unique index leaving no unused entries in the Q array.
!! Assume that real Y_{l,m} harmonic was stored in Q(l, 2*|m|, :, :) for even (cosine, positive-m) hamonic, or in Q(l, 2*|m|-1, :, :) for odd (sine, negative-m) harmonic.
!! Then consecutive entries in compressed Q(:, inout, r) array should be related to Y_{l,m} harmonics, starting from index 0 as follows
!! Y_{0,0}, Y_{1,0}, ..., Y_{lmax, 0}, Y_{1,-1}, Y_{2,-1}, ..., Y_{lmax,-1}, Y_{1,1}, Y_{2,1}, ..., Y_{lmax,1}, Y_{2,-2}, Y_{3,-2}, ..., Y_{lmax,-2},
!! Y_{2,2}, Y_{3,2}, ..., Y_{lmax,2}, ...,, Y(lmax-1,-(mmax-1)), Y(lmax,-(mmax-1)), Y(lmax-1,mmax-1), Y(lmax,mmax-1), Y(lmax,-mmax), Y(lmax,mmax)
!! Does it looks a bit cryptic? I agree.
!!

   elemental integer function lm(l, m)

      implicit none

      integer, intent(in) :: l, m

      lm = l + m*lmax - int((m-1)/2)*int(m/2)

   end function lm

!!$ ============================================================================
!!
!! Initialization routine, called from init_multigrid
!!

   subroutine init_multipole(mb_alloc, cgrid)

      use errh,          only: die
      use types,         only: grid_container
      use mpisetup,      only: proc
      use multigridvars, only: level_min, level_max, lvl, eff_dim

      implicit none

      real,                 intent(inout) :: mb_alloc               !< multigrid allocation counter
      type(grid_container), intent(in)    :: cgrid                  !< copy of grid variables

      integer, dimension(4) :: aerr = 0
      integer               :: l,m

      ! external face coordinates
      fbnd_x(LOW:HIGH) = [ cgrid%xminb, cgrid%xmaxb ]
      fbnd_y(LOW:HIGH) = [ cgrid%yminb, cgrid%ymaxb ]
      fbnd_z(LOW:HIGH) = [ cgrid%zminb, cgrid%zmaxb ]

      ! assume that Center of Mass is approximately in the center of computational domain by default
      CoM(0) = 1.
      CoM(XDIR) = (cgrid%xmax + cgrid%xmin)/2.
      CoM(YDIR) = (cgrid%ymax + cgrid%ymin)/2.
      CoM(ZDIR) = (cgrid%zmax + cgrid%zmin)/2.

      if (eff_dim /= NDIM) call die("[multipole:init_multipole] Only 3D is supported")

      !multipole moments
      if (mmax > lmax) then
         if (proc == 0) write(*,'(a)')"[multipole:init_multipole] Warning: mmax reduced to lmax"
         mmax = lmax
      end if
      if (mmax < 0) mmax = lmax

      if (coarsen_multipole < 0) coarsen_multipole = 0
      if (level_max - coarsen_multipole < level_min) then
         if (proc == 0) write(*,'(a)')"[multipole:init_multipole] Warning: too deep multipole coarsening, setting level_min."
         coarsen_multipole = level_max - level_min
      end if
      lmpole => lvl(level_max - coarsen_multipole)

      if (.not. use_point_monopole) then
         if (allocated(rn) .or. allocated(irn) .or. allocated(sfac) .or. allocated(cfac)) call die("[multipole:init_multipole] rn, irn, sfac or cfac already allocated")
         allocate(  rn(0:lmax), stat=aerr(1))
         allocate( irn(0:lmax), stat=aerr(2))
         allocate(sfac(0:mmax), stat=aerr(3))
         allocate(cfac(0:mmax), stat=aerr(4))
         if (any(aerr(1:4) /= 0)) call die("[multipole:init_multipole] Allocation error: rn, irn sfac or cfac")
         mb_alloc = mb_alloc + size(rn) + size(irn) + size(sfac) + size(cfac)

         drq = min(lmpole%dx, lmpole%dy, lmpole%dz) / 2.
         rqbin = int(sqrt((cgrid%xmax - cgrid%xmin)**2 + (cgrid%ymax - cgrid%ymin)**2 + (cgrid%zmax - cgrid%zmin)**2)/drq) + 1
         ! arithmetic average of the closest and farthest points of computational domain with respect to its center
         rscale = ( min((cgrid%xmax - cgrid%xmin),     (cgrid%ymax - cgrid%ymin),     (cgrid%zmax - cgrid%zmin)) + &
              &    sqrt((cgrid%xmax - cgrid%xmin)**2 + (cgrid%ymax - cgrid%ymin)**2 + (cgrid%zmax - cgrid%zmin)**2) )/4.
         if (allocated(k12) .or. allocated(ofact) .or. allocated(Q)) call die("[multipole:init_multipole] k12, ofact or Q already allocated")
         allocate(   k12(2, 1:lmax, 0:mmax), stat=aerr(1))
         allocate(ofact(0:lm(lmax, 2*mmax)), stat=aerr(2))
         allocate(    Q(0:lm(lmax, 2*mmax), INSIDE:OUTSIDE, 0:rqbin), stat=aerr(3))
         if (any(aerr(1:3) /= 0)) call die("[multipole:init_multipole] Allocation error: k12, ofact or Q")
         mb_alloc = mb_alloc + size(k12) + size(ofact) + size(Q)

         ofact(:) = 0. ! prevent FPE spurious exceptions in multipole:img_mass2moments
         do l = 1, lmax
            do m = 0, min(l, mmax)
               if (m == 0) then
                  ofact(l) = 2. ! lm(l,0)
               else                     ! ofact(l, 2*m) = ((-1)^m (2m-1)!!)**2 (l-m)! / (l+m)! ; Should work up to m=512 and even beyond
                  ofact(lm(l, 2*m))   = ofact(lm(l, 2*(m-1))) * dble(1 - 2*m)**2 / dble((l+m) * (l-m+1))
                  ofact(lm(l, 2*m-1)) = ofact(lm(l, 2*m))
               end if
               if (l>m) then
                  k12(1, l, m) = dble(2 * l - 1) / dble(l - m)
                  k12(2, l, m) = dble(l + m - 1) / dble(l - m)
               end if
            end do
         end do
         ofact(0:lmax) = 1. ! lm(0:lmax,0)

         cfac(0) = 1.e0
         sfac(0) = 0.e0

      end if

   end subroutine init_multipole

!!$ ============================================================================
!!
!! Multipole cleanup
!!

   subroutine cleanup_multipole
      implicit none

      if  (allocated(rn))    deallocate(rn)
      if  (allocated(irn))   deallocate(irn)
      if  (allocated(sfac))  deallocate(sfac)
      if  (allocated(cfac))  deallocate(cfac)
      if  (allocated(Q))     deallocate(Q)
      if  (allocated(k12))   deallocate(k12)
      if  (allocated(ofact)) deallocate(ofact)

   end subroutine cleanup_multipole

!!$ ============================================================================
!!
!! Multipole solver
!!
!! \todo impmprove multipole expansion on coarser grids (see. "A Scalable Parallel Poisson Solver
!! in Three Dimensions with Infinite-Domain Boundary Conditions" by McCorquodale, Colella, Balls and Baden).
!! Coarsening by one level would reduce the multipole costs by a factor of 4.
!!

   subroutine multipole_solver

      use multigridvars,      only: level_max, solution
      use multigridhelpers,   only: dirtyH, dirty_debug
      use multigridbasefuncs, only: zero_boundaries, restrict_level

      implicit none

      integer :: lev

      if (dirty_debug) then
         lmpole%bnd_x(:, :, :) = dirtyH
         lmpole%bnd_y(:, :, :) = dirtyH
         lmpole%bnd_z(:, :, :) = dirtyH
      else
         call zero_boundaries(lmpole%level)
      end if

      if (lmpole%level <  level_max) then
         do lev = level_max, lmpole%level + 1, -1
            call restrict_level(lev, solution) ! Overkill, only some layers next to external boundary are needed.
         end do                                ! An alternative: do potential2img_mass on the roof and restrict bnd_[xyz] data.
      end if
      call potential2img_mass

      if (use_point_monopole) then
         call find_img_CoM
         call isolated_monopole
      else
         ! results seems to be slightly better without find_img_CoM.
         ! With CoM or when it is known than CoM is close to the domain center one may try to save some CPU time by lowering mmax.
         call img_mass2moments
         call moments2bnd_potential
      end if

      if (lmpole%level <  level_max) then
         do lev = lmpole%level, level_max - 1
            call prolong_ext_bnd(lev)
         end do
      end if

   end subroutine multipole_solver

!!$ ============================================================================
!!
!! Set boundary potential from monopole source. Fill lmpole%bnd_[xyz] arrays with expected values of the gravitational potential at external face of computational domain.
!! This is a simplified approach that can be used for tests and as a fast replacement for the
!! multipole boundary solver for nearly spherically symmetric source distributions.
!! The isolated_monopole subroutine ignores the radial profile of the monopole
!!

   subroutine isolated_monopole

      use multigridvars,   only: XDIR, YDIR, ZDIR, LOW, HIGH, is_external, XLO, XHI, YLO, YHI, ZLO, ZHI
      use constants, only: newtong

      implicit none

      integer :: i, j, k
      real    :: r2

      if (is_external(XLO) .or. is_external(XHI)) then
         do j = lmpole%js, lmpole%je
            do k = lmpole%ks, lmpole%ke
               r2 = (lmpole%y(j)-CoM(YDIR))**2 + (lmpole%z(k) - CoM(ZDIR))**2
               if (is_external(XLO)) lmpole%bnd_x(j, k, LOW)  = - newtong * CoM(0) / sqrt(r2 + (fbnd_x(LOW) -CoM(XDIR))**2)
               if (is_external(XHI)) lmpole%bnd_x(j, k, HIGH) = - newtong * CoM(0) / sqrt(r2 + (fbnd_x(HIGH)-CoM(XDIR))**2)
            end do
         end do
      end if

      if (is_external(YLO) .or. is_external(YHI)) then
         do i = lmpole%is, lmpole%ie
            do k = lmpole%ks, lmpole%ke
               r2 = (lmpole%x(i)-CoM(XDIR))**2 + (lmpole%z(k) - CoM(ZDIR))**2
               if (is_external(YLO)) lmpole%bnd_y(i, k, LOW)  = - newtong * CoM(0) / sqrt(r2 + (fbnd_y(LOW) -CoM(YDIR))**2)
               if (is_external(YHI)) lmpole%bnd_y(i, k, HIGH) = - newtong * CoM(0) / sqrt(r2 + (fbnd_y(HIGH)-CoM(YDIR))**2)
            end do
         end do
      end if

      if (is_external(ZLO) .or. is_external(ZHI)) then
         do i = lmpole%is, lmpole%ie
            do j = lmpole%js, lmpole%je
               r2 = (lmpole%x(i)-CoM(XDIR))**2 + (lmpole%y(j) - CoM(YDIR))**2
               if (is_external(ZLO)) lmpole%bnd_z(i, j, LOW)  = - newtong * CoM(0) / sqrt(r2 + (fbnd_z(LOW) -CoM(ZDIR))**2)
               if (is_external(ZHI)) lmpole%bnd_z(i, j, HIGH) = - newtong * CoM(0) / sqrt(r2 + (fbnd_z(HIGH)-CoM(ZDIR))**2)
            end do
         end do
      end if

   end subroutine isolated_monopole

!!$ ============================================================================
!!
!! Find total mass and its center
!! This routine does the summation only on external boundaries
!!

   subroutine find_img_CoM

      use multigridvars,   only: is_external, XLO, XHI, YLO, YHI, ZLO, ZHI, XDIR, YDIR, ZDIR, LOW, HIGH
      use errh,            only: die
      use mpisetup,        only: comm3d, ierr, MPI_DOUBLE_PRECISION, MPI_SUM

      implicit none

      real, dimension(0:NDIM) :: lsum, dsum

      lsum(:) = 0.

      if (is_external(XLO)) then
         dsum(0)    =      sum( lmpole%bnd_x(lmpole%js:lmpole%je, lmpole%ks:lmpole%ke, LOW) )
         dsum(XDIR) = dsum(0) * fbnd_x(LOW)
         dsum(YDIR) = sum( sum( lmpole%bnd_x(lmpole%js:lmpole%je, lmpole%ks:lmpole%ke, LOW),  dim=2) * lmpole%y(lmpole%js:lmpole%je) )
         dsum(ZDIR) = sum( sum( lmpole%bnd_x(lmpole%js:lmpole%je, lmpole%ks:lmpole%ke, LOW),  dim=1) * lmpole%z(lmpole%ks:lmpole%ke) )
         lsum(:)    = lsum(:) + dsum(:) * lmpole%dyz
      end if

      if (is_external(XHI)) then
         dsum(0)    =      sum( lmpole%bnd_x(lmpole%js:lmpole%je, lmpole%ks:lmpole%ke, HIGH) )
         dsum(XDIR) = dsum(0) * fbnd_x(HIGH)
         dsum(YDIR) = sum( sum( lmpole%bnd_x(lmpole%js:lmpole%je, lmpole%ks:lmpole%ke, HIGH), dim=2) * lmpole%y(lmpole%js:lmpole%je) )
         dsum(ZDIR) = sum( sum( lmpole%bnd_x(lmpole%js:lmpole%je, lmpole%ks:lmpole%ke, HIGH), dim=1) * lmpole%z(lmpole%ks:lmpole%ke) )
         lsum(:)    = lsum(:) + dsum(:) * lmpole%dyz
      end if

      if (is_external(YLO)) then
         dsum(0)    =      sum( lmpole%bnd_y(lmpole%is:lmpole%ie, lmpole%ks:lmpole%ke, LOW) )
         dsum(XDIR) = sum( sum( lmpole%bnd_y(lmpole%is:lmpole%ie, lmpole%ks:lmpole%ke, LOW),  dim=2) * lmpole%x(lmpole%is:lmpole%ie) )
         dsum(YDIR) = dsum(0) * fbnd_y(LOW)
         dsum(ZDIR) = sum( sum( lmpole%bnd_y(lmpole%is:lmpole%ie, lmpole%ks:lmpole%ke, LOW),  dim=1) * lmpole%z(lmpole%ks:lmpole%ke) )
         lsum(:)    = lsum(:) + dsum(:) * lmpole%dxz
      end if

      if (is_external(YHI)) then
         dsum(0)    =      sum( lmpole%bnd_y(lmpole%is:lmpole%ie, lmpole%ks:lmpole%ke, HIGH) )
         dsum(XDIR) = sum( sum( lmpole%bnd_y(lmpole%is:lmpole%ie, lmpole%ks:lmpole%ke, HIGH), dim=2) * lmpole%x(lmpole%is:lmpole%ie) )
         dsum(YDIR) = dsum(0) * fbnd_y(HIGH)
         dsum(ZDIR) = sum( sum( lmpole%bnd_y(lmpole%is:lmpole%ie, lmpole%ks:lmpole%ke, HIGH), dim=1) * lmpole%z(lmpole%ks:lmpole%ke) )
         lsum(:)    = lsum(:) + dsum(:) * lmpole%dxz
      end if

      if (is_external(ZLO)) then
         dsum(0)    =      sum( lmpole%bnd_z(lmpole%is:lmpole%ie, lmpole%js:lmpole%je, LOW) )
         dsum(XDIR) = sum( sum( lmpole%bnd_z(lmpole%is:lmpole%ie, lmpole%js:lmpole%je, LOW),  dim=2) * lmpole%x(lmpole%is:lmpole%ie) )
         dsum(YDIR) = sum( sum( lmpole%bnd_z(lmpole%is:lmpole%ie, lmpole%js:lmpole%je, LOW),  dim=1) * lmpole%z(lmpole%ks:lmpole%ke) )
         dsum(ZDIR) = dsum(0) * fbnd_z(LOW)
         lsum(:)    = lsum(:) + dsum(:) * lmpole%dxy
      end if

      if (is_external(ZHI)) then
         dsum(0)    =      sum( lmpole%bnd_z(lmpole%is:lmpole%ie, lmpole%js:lmpole%je, HIGH) )
         dsum(XDIR) = sum( sum( lmpole%bnd_z(lmpole%is:lmpole%ie, lmpole%js:lmpole%je, HIGH), dim=2) * lmpole%x(lmpole%is:lmpole%ie) )
         dsum(YDIR) = sum( sum( lmpole%bnd_z(lmpole%is:lmpole%ie, lmpole%js:lmpole%je, HIGH), dim=1) * lmpole%z(lmpole%ks:lmpole%ke) )
         dsum(ZDIR) = dsum(0) * fbnd_z(HIGH)
         lsum(:)    = lsum(:) + dsum(:) * lmpole%dxy
      end if

      call MPI_Allreduce(lsum(0:NDIM), CoM(0:NDIM), 4, MPI_DOUBLE_PRECISION, MPI_SUM, comm3d, ierr)

      if (CoM(0) /= 0.) then
         CoM(XDIR:ZDIR) = CoM(XDIR:ZDIR) / CoM(0)
      else
         call die("[multipole:find_img_CoM] Total mass == 0")
      end if

   end subroutine find_img_CoM

!!$ ============================================================================
!!
!! Convert potential into image mass. This way we reduce a 3D problem to a 2D one.
!! There will be work imbalance here because different PEs may operate on different amount of external boundary data
!!

   subroutine potential2img_mass

      use multigridvars,   only: is_external, XLO, XHI, YLO, YHI, ZLO, ZHI, LOW, HIGH, solution

      implicit none

      real, parameter :: a1 = -2., a2 = (-2. - a1)/3. ! interpolation parameters;   <---- a1=-2 => a2=0
      ! a1==-2. is the simplest, low order choice, gives best agreement of total mass and CoM location when compared to 3-D integration

      if (is_external(XLO)) lmpole%bnd_x(             lmpole%js:lmpole%je, lmpole%ks:lmpole%ke, LOW) =    ( &
           &           a1 * lmpole%mgvar(lmpole%is,   lmpole%js:lmpole%je, lmpole%ks:lmpole%ke, solution) + &
           &           a2 * lmpole%mgvar(lmpole%is+1, lmpole%js:lmpole%je, lmpole%ks:lmpole%ke, solution) ) / lmpole%dx

      if (is_external(XHI)) lmpole%bnd_x(             lmpole%js:lmpole%je, lmpole%ks:lmpole%ke, HIGH) =   ( &
           &           a1 * lmpole%mgvar(lmpole%ie,   lmpole%js:lmpole%je, lmpole%ks:lmpole%ke, solution) + &
           &           a2 * lmpole%mgvar(lmpole%ie-1, lmpole%js:lmpole%je, lmpole%ks:lmpole%ke, solution) ) / lmpole%dx

      if (is_external(YLO)) lmpole%bnd_y(lmpole%is:lmpole%ie,              lmpole%ks:lmpole%ke, LOW) =    ( &
           &           a1 * lmpole%mgvar(lmpole%is:lmpole%ie, lmpole%js,   lmpole%ks:lmpole%ke, solution) + &
           &           a2 * lmpole%mgvar(lmpole%is:lmpole%ie, lmpole%js+1, lmpole%ks:lmpole%ke, solution) ) / lmpole%dy

      if (is_external(YHI)) lmpole%bnd_y(lmpole%is:lmpole%ie,              lmpole%ks:lmpole%ke, HIGH) =   ( &
           &           a1 * lmpole%mgvar(lmpole%is:lmpole%ie, lmpole%je,   lmpole%ks:lmpole%ke, solution) + &
           &           a2 * lmpole%mgvar(lmpole%is:lmpole%ie, lmpole%je-1, lmpole%ks:lmpole%ke, solution) ) / lmpole%dy

      if (is_external(ZLO)) lmpole%bnd_z(lmpole%is:lmpole%ie, lmpole%js:lmpole%je,              LOW) =    ( &
           &           a1 * lmpole%mgvar(lmpole%is:lmpole%ie, lmpole%js:lmpole%je, lmpole%ks,   solution) + &
           &           a2 * lmpole%mgvar(lmpole%is:lmpole%ie, lmpole%js:lmpole%je, lmpole%ks+1, solution) ) / lmpole%dz

      if (is_external(ZHI)) lmpole%bnd_z(lmpole%is:lmpole%ie, lmpole%js:lmpole%je,              HIGH) =   ( &
           &           a1 * lmpole%mgvar(lmpole%is:lmpole%ie, lmpole%js:lmpole%je, lmpole%ke,   solution) + &
           &           a2 * lmpole%mgvar(lmpole%is:lmpole%ie, lmpole%js:lmpole%je, lmpole%ke-1, solution) ) / lmpole%dz

   end subroutine potential2img_mass

!!$ ============================================================================
!!
!! Prolong boundaries wrapper
!!

   subroutine prolong_ext_bnd(lev)

      use multigridvars,   only: is_external
      use errh,            only: die

      implicit none

      integer, intent(in) :: lev !< level to prolong from

      if (abs(ord_prolong_mpole) > 2) call die("[multipole:prolong_ext_bnd] interpolation order too high")

      if (any(is_external(:))) then
         if (ord_prolong_mpole == 0) then
            call prolong_ext_bnd0(lev)
         else
            call prolong_ext_bnd2(lev)
         end if
      end if

   end subroutine prolong_ext_bnd

!!$ ============================================================================
!!
!! Prolong boundaries by injection.
!!

   subroutine prolong_ext_bnd0(lev)

      use multigridvars,   only: lvl, is_external, level_max, XLO, XHI, YLO, YHI, ZLO, ZHI, HIGH, LOW, eff_dim, NDIM
      use errh,            only: die

      implicit none

      integer, intent(in) :: lev !< level to prolong from

      type(plvl), pointer :: coarse, fine

      if (lev >= level_max) return

      if (eff_dim<NDIM) call die("[multigridmultipole:prolong_ext_bnd0] 1D and 2D not finished")

      coarse => lvl(lev)
      fine   => lvl(lev + 1)

      if (is_external(XLO)) then
         fine%bnd_x(fine%js  :fine%je-1:2, fine%ks  :fine%ke-1:2, LOW)  = coarse%bnd_x(coarse%js:coarse%je,   coarse%ks:coarse%ke,   LOW)
         fine%bnd_x(fine%js+1:fine%je  :2, fine%ks  :fine%ke-1:2, LOW)  = fine  %bnd_x(fine%js  :fine%je-1:2, fine%ks  :fine%ke-1:2, LOW)
         fine%bnd_x(fine%js  :fine%je,     fine%ks+1:fine%ke  :2, LOW)  = fine  %bnd_x(fine%js  :fine%je,     fine%ks  :fine%ke-1:2, LOW)
      end if

       if (is_external(XHI)) then
         fine%bnd_x(fine%js  :fine%je-1:2, fine%ks  :fine%ke-1:2, HIGH) = coarse%bnd_x(coarse%js:coarse%je,   coarse%ks:coarse%ke,   HIGH)
         fine%bnd_x(fine%js+1:fine%je  :2, fine%ks  :fine%ke-1:2, HIGH) = fine  %bnd_x(fine%js  :fine%je-1:2, fine%ks  :fine%ke-1:2, HIGH)
         fine%bnd_x(fine%js  :fine%je,     fine%ks+1:fine%ke  :2, HIGH) = fine  %bnd_x(fine%js  :fine%je,     fine%ks  :fine%ke-1:2, HIGH)
      end if

      if (is_external(YLO)) then
         fine%bnd_y(fine%is  :fine%ie-1:2, fine%ks  :fine%ke-1:2, LOW)  = coarse%bnd_y(coarse%is:coarse%ie,   coarse%ks:coarse%ke,   LOW)
         fine%bnd_y(fine%is+1:fine%ie  :2, fine%ks  :fine%ke-1:2, LOW)  = fine  %bnd_y(fine%is  :fine%ie-1:2, fine%ks  :fine%ke-1:2, LOW)
         fine%bnd_y(fine%is  :fine%ie,     fine%ks+1:fine%ke  :2, LOW)  = fine  %bnd_y(fine%is  :fine%ie,     fine%ks  :fine%ke-1:2, LOW)
      end if

       if (is_external(YHI)) then
         fine%bnd_y(fine%is  :fine%ie-1:2, fine%ks  :fine%ke-1:2, HIGH) = coarse%bnd_y(coarse%is:coarse%ie,   coarse%ks:coarse%ke,   HIGH)
         fine%bnd_y(fine%is+1:fine%ie  :2, fine%ks  :fine%ke-1:2, HIGH) = fine  %bnd_y(fine%is  :fine%ie-1:2, fine%ks  :fine%ke-1:2, HIGH)
         fine%bnd_y(fine%is  :fine%ie,     fine%ks+1:fine%ke  :2, HIGH) = fine  %bnd_y(fine%is  :fine%ie,     fine%ks  :fine%ke-1:2, HIGH)
      end if

      if (is_external(ZLO)) then
         fine%bnd_z(fine%is  :fine%ie-1:2, fine%js  :fine%je-1:2, LOW)  = coarse%bnd_z(coarse%is:coarse%ie,   coarse%js:coarse%je,   LOW)
         fine%bnd_z(fine%is+1:fine%ie  :2, fine%js  :fine%je-1:2, LOW)  = fine  %bnd_z(fine%is  :fine%ie-1:2, fine%js  :fine%je-1:2, LOW)
         fine%bnd_z(fine%is  :fine%ie,     fine%js+1:fine%je  :2, LOW)  = fine  %bnd_z(fine%is  :fine%ie,     fine%js  :fine%je-1:2, LOW)
      end if

       if (is_external(ZHI)) then
         fine%bnd_z(fine%is  :fine%ie-1:2, fine%js  :fine%je-1:2, HIGH) = coarse%bnd_z(coarse%is:coarse%ie,   coarse%js:coarse%je,   HIGH)
         fine%bnd_z(fine%is+1:fine%ie  :2, fine%js  :fine%je-1:2, HIGH) = fine  %bnd_z(fine%is  :fine%ie-1:2, fine%js  :fine%je-1:2, HIGH)
         fine%bnd_z(fine%is  :fine%ie,     fine%js+1:fine%je  :2, HIGH) = fine  %bnd_z(fine%is  :fine%ie,     fine%js  :fine%je-1:2, HIGH)
      end if

   end subroutine prolong_ext_bnd0

!!$ ============================================================================
!!
!! Prolong boundaries by linear or quadratic interpolation.
!!
!! some code is replicated from prolong_faces
!! \todo write something more general, a routine that takes arrays or pointers and does the 2D prolongation
!!

   subroutine prolong_ext_bnd2(lev)

      use multigridvars,   only: lvl, is_external, level_max, XLO, XHI, YLO, YHI, ZLO, ZHI, HIGH, LOW, eff_dim, NDIM
      use errh,            only: die

      implicit none

      integer, intent(in) :: lev !< level to prolong from

      type(plvl), pointer :: coarse, fine

      integer                       :: i, j, k
      real, parameter, dimension(3) :: p0  = [ 0.,       1.,     0.     ] ! injection
      real, parameter, dimension(3) :: p1  = [ 0.,       3./4.,  1./4.  ] ! 1D linear prolongation stencil
      real, parameter, dimension(3) :: p2i = [ -1./8.,   1.,     1./8.  ] ! 1D integral cubic prolongation stencil
      real, parameter, dimension(3) :: p2d = [ -3./32., 30./32., 5./32. ] ! 1D direct cubic prolongation stencil
      real, dimension(-1:1)         :: p
      real, dimension(-1:1,-1:1,2,2):: pp   ! 2D prolongation stencil

      if (lev >= level_max) return

      if (eff_dim<NDIM) call die("[multigridmultipole:prolong_ext_bnd2] 1D and 2D not finished")

      select case(ord_prolong_mpole)
         case(0)
            p(:) = p0(:)
         case(1,-1)
            p(:) = p1(:)
         case(2)
            p(:) = p2i(:)
         case(-2)
            p(:) = p2d(:)
         case default
            p(:) = p0(:)
      end select

      do i = -1, 1
         pp(i,:,1,1) = p( i)*p(:)
         pp(i,:,1,2) = p( i)*p(1:-1:-1)
         pp(i,:,2,1) = p(-i)*p(:)
         pp(i,:,2,2) = p(-i)*p(1:-1:-1)
      end do

      coarse => lvl(lev)
      fine   => lvl(lev + 1)

      !at edges and corners we can only inject
      call prolong_ext_bnd0(lev) ! overkill: replace this

      if (is_external(XLO)) then
         do j = coarse%js+1, coarse%je-1
            do k = coarse%ks+1, coarse%ke-1
               fine%bnd_x(-fine%js+2*j,  -fine%ks+2*k,  LOW) =sum(pp(:,:,1,1) * coarse%bnd_x(j-1:j+1,k-1:k+1,LOW))
               fine%bnd_x(-fine%js+2*j+1,-fine%ks+2*k,  LOW) =sum(pp(:,:,2,1) * coarse%bnd_x(j-1:j+1,k-1:k+1,LOW))
               fine%bnd_x(-fine%js+2*j,  -fine%ks+2*k+1,LOW) =sum(pp(:,:,1,2) * coarse%bnd_x(j-1:j+1,k-1:k+1,LOW))
               fine%bnd_x(-fine%js+2*j+1,-fine%ks+2*k+1,LOW) =sum(pp(:,:,2,2) * coarse%bnd_x(j-1:j+1,k-1:k+1,LOW))
            end do
         end do
      end if

      if (is_external(XHI)) then
         do j = coarse%js+1, coarse%je-1
            do k = coarse%ks+1, coarse%ke-1
               fine%bnd_x(-fine%js+2*j,  -fine%ks+2*k,  HIGH)=sum(pp(:,:,1,1) * coarse%bnd_x(j-1:j+1,k-1:k+1,HIGH))
               fine%bnd_x(-fine%js+2*j+1,-fine%ks+2*k,  HIGH)=sum(pp(:,:,2,1) * coarse%bnd_x(j-1:j+1,k-1:k+1,HIGH))
               fine%bnd_x(-fine%js+2*j,  -fine%ks+2*k+1,HIGH)=sum(pp(:,:,1,2) * coarse%bnd_x(j-1:j+1,k-1:k+1,HIGH))
               fine%bnd_x(-fine%js+2*j+1,-fine%ks+2*k+1,HIGH)=sum(pp(:,:,2,2) * coarse%bnd_x(j-1:j+1,k-1:k+1,HIGH))
            end do
         end do
      end if


      if (is_external(YLO)) then
         do i = coarse%is+1, coarse%ie-1
            do k = coarse%ks+1, coarse%ke-1
               fine%bnd_y(-fine%is+2*i,  -fine%ks+2*k,  LOW) =sum(pp(:,:,1,1) * coarse%bnd_y(i-1:i+1,k-1:k+1,LOW))
               fine%bnd_y(-fine%is+2*i+1,-fine%ks+2*k,  LOW) =sum(pp(:,:,2,1) * coarse%bnd_y(i-1:i+1,k-1:k+1,LOW))
               fine%bnd_y(-fine%is+2*i,  -fine%ks+2*k+1,LOW) =sum(pp(:,:,1,2) * coarse%bnd_y(i-1:i+1,k-1:k+1,LOW))
               fine%bnd_y(-fine%is+2*i+1,-fine%ks+2*k+1,LOW) =sum(pp(:,:,2,2) * coarse%bnd_y(i-1:i+1,k-1:k+1,LOW))
            end do
         end do
      end if

      if (is_external(YHI)) then
         do i = coarse%is+1, coarse%ie-1
            do k = coarse%ks+1, coarse%ke-1
               fine%bnd_y(-fine%is+2*i,  -fine%ks+2*k,  HIGH) =sum(pp(:,:,1,1) * coarse%bnd_y(i-1:i+1,k-1:k+1,HIGH))
               fine%bnd_y(-fine%is+2*i+1,-fine%ks+2*k,  HIGH) =sum(pp(:,:,2,1) * coarse%bnd_y(i-1:i+1,k-1:k+1,HIGH))
               fine%bnd_y(-fine%is+2*i,  -fine%ks+2*k+1,HIGH) =sum(pp(:,:,1,2) * coarse%bnd_y(i-1:i+1,k-1:k+1,HIGH))
               fine%bnd_y(-fine%is+2*i+1,-fine%ks+2*k+1,HIGH) =sum(pp(:,:,2,2) * coarse%bnd_y(i-1:i+1,k-1:k+1,HIGH))
            end do
         end do
      end if

      if (is_external(ZLO)) then
         do i = coarse%is+1, coarse%ie-1
            do j = coarse%js+1, coarse%je-1
               fine%bnd_z(-fine%is+2*i,  -fine%js+2*j,  LOW) =sum(pp(:,:,1,1) * coarse%bnd_z(i-1:i+1,j-1:j+1,LOW))
               fine%bnd_z(-fine%is+2*i+1,-fine%js+2*j,  LOW) =sum(pp(:,:,2,1) * coarse%bnd_z(i-1:i+1,j-1:j+1,LOW))
               fine%bnd_z(-fine%is+2*i,  -fine%js+2*j+1,LOW) =sum(pp(:,:,1,2) * coarse%bnd_z(i-1:i+1,j-1:j+1,LOW))
               fine%bnd_z(-fine%is+2*i+1,-fine%js+2*j+1,LOW) =sum(pp(:,:,2,2) * coarse%bnd_z(i-1:i+1,j-1:j+1,LOW))
            end do
         end do
      end if

      if (is_external(ZHI)) then
         do i = coarse%is+1, coarse%ie-1
            do j = coarse%js+1, coarse%je-1
               fine%bnd_z(-fine%is+2*i,  -fine%js+2*j,  HIGH) =sum(pp(:,:,1,1) * coarse%bnd_z(i-1:i+1,j-1:j+1,HIGH))
               fine%bnd_z(-fine%is+2*i+1,-fine%js+2*j,  HIGH) =sum(pp(:,:,2,1) * coarse%bnd_z(i-1:i+1,j-1:j+1,HIGH))
               fine%bnd_z(-fine%is+2*i,  -fine%js+2*j+1,HIGH) =sum(pp(:,:,1,2) * coarse%bnd_z(i-1:i+1,j-1:j+1,HIGH))
               fine%bnd_z(-fine%is+2*i+1,-fine%js+2*j+1,HIGH) =sum(pp(:,:,2,2) * coarse%bnd_z(i-1:i+1,j-1:j+1,HIGH))
            end do
         end do
      end if

   end subroutine prolong_ext_bnd2

!!$ ============================================================================
!!
!! Compute multipole moments for image mass
!! \todo distribute excess of work more evenly (important only for large number of PEs, ticket:43)
!!

   subroutine img_mass2moments

      use multigridvars,   only: is_external, XLO, XHI, YLO, YHI, ZLO, ZHI, LOW, HIGH, XDIR, YDIR, ZDIR
      use mpisetup, only: comm3d, ierr, MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_SUM, MPI_MIN, MPI_MAX, MPI_IN_PLACE

      implicit none

      integer :: i, j, k, r, rr

      ! reset the multipole data
      Q(:, :, :) = 0.

      irmax = 0
      irmin = rqbin

      ! scan
      if (is_external(XLO) .or. is_external(XHI)) then
         do j = lmpole%js, lmpole%je
            do k = lmpole%ks, lmpole%ke
               if (is_external(XLO)) call point2moments(lmpole%bnd_x(j, k, LOW) *lmpole%dyz, fbnd_x(LOW) -CoM(XDIR), lmpole%y(j)-CoM(YDIR),  lmpole%z(k)-CoM(ZDIR))
               if (is_external(XHI)) call point2moments(lmpole%bnd_x(j, k, HIGH)*lmpole%dyz, fbnd_x(HIGH)-CoM(XDIR), lmpole%y(j)-CoM(YDIR),  lmpole%z(k)-CoM(ZDIR))
            end do
         end do
      end if

      if (is_external(YLO) .or. is_external(YHI)) then
         do i = lmpole%is, lmpole%ie
            do k = lmpole%ks, lmpole%ke
               if (is_external(YLO)) call point2moments(lmpole%bnd_y(i, k, LOW) *lmpole%dxz, lmpole%x(i)-CoM(XDIR),  fbnd_y(LOW) -CoM(YDIR), lmpole%z(k)-CoM(ZDIR))
               if (is_external(YHI)) call point2moments(lmpole%bnd_y(i, k, HIGH)*lmpole%dxz, lmpole%x(i)-CoM(XDIR),  fbnd_y(HIGH)-CoM(YDIR), lmpole%z(k)-CoM(ZDIR))
            end do
         end do
      end if

      if (is_external(ZLO) .or. is_external(ZHI)) then
         do i = lmpole%is, lmpole%ie
            do j = lmpole%js, lmpole%je
               if (is_external(ZLO)) call point2moments(lmpole%bnd_z(i, j, LOW) *lmpole%dxy, lmpole%x(i)-CoM(XDIR),  lmpole%y(j)-CoM(YDIR),  fbnd_z(LOW) -CoM(ZDIR))
               if (is_external(ZHI)) call point2moments(lmpole%bnd_z(i, j, HIGH)*lmpole%dxy, lmpole%x(i)-CoM(XDIR),  lmpole%y(j)-CoM(YDIR),  fbnd_z(HIGH)-CoM(ZDIR))
            end do
         end do
      end if

      call MPI_Allreduce(MPI_IN_PLACE, irmin, 1, MPI_INTEGER, MPI_MIN, comm3d, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, irmax, 1, MPI_INTEGER, MPI_MAX, comm3d, ierr)

      ! integrate radially and apply normalization factor (the (4 \pi)/(2 l  + 1) terms cancel out)
      rr = max(1, irmin)
      Q(:, INSIDE, rr-1) = Q(:, INSIDE, rr-1) * ofact(:)
      do r = rr, irmax
         Q(:, INSIDE, r) = Q(:, INSIDE, r) * ofact(:) + Q(:, INSIDE, r-1)
      end do

      Q(:, OUTSIDE, irmax+1) = Q(:, OUTSIDE, irmax+1) * ofact(:)
      do r = irmax, rr, -1
         Q(:, OUTSIDE, r) = Q(:, OUTSIDE, r) * ofact(:) + Q(:, OUTSIDE, r+1)
      end do

      call MPI_Allreduce(MPI_IN_PLACE, Q(:, :, irmin:irmax), size(Q(:, :, irmin:irmax)), MPI_DOUBLE_PRECISION, MPI_SUM, comm3d, ierr)

   end subroutine img_mass2moments

!!$ ============================================================================
!!
!! Compute multipole moments for a single point
!!
!! \todo try to improve accuracy with linear interpolaion over radius
!!

   subroutine point2moments(mass, x, y, z)

      implicit none

      real, intent(in) :: mass    !< mass of the contributing point
      real, intent(in) :: x, y, z !< coordinates of the contributing point

      real    :: sin_th, cos_th
      real    :: Ql, Ql1, Ql2
      integer :: l, m, ir, m2s, m2c

      call geomfac4moments(mass, x, y, z, sin_th, cos_th, ir)

      ! monopole, the (0,0) moment; P_0 = 1.
      Q(0, INSIDE,  ir) = Q(0, INSIDE,  ir) +  rn(0)
      Q(0, OUTSIDE, ir) = Q(0, OUTSIDE, ir) + irn(0)

      ! axisymmetric (l,0) moments
      ! Legendre polynomial recurrence: l P_l = x (2l-1) P_{l-1} - (l-1) P_{l-2}, x \eqiv \cos(\theta)
      Ql2 = 0.
      Ql1 = 1.
      do l = 1, lmax
         Ql = cos_th * k12(1, l, 0) * Ql1 - k12(2, l, 0) * Ql2
         Q(l, INSIDE,  ir) = Q(l, INSIDE,  ir) +  rn(l) * Ql
         Q(l, OUTSIDE, ir) = Q(l, OUTSIDE, ir) + irn(l) * Ql
         Ql2 = Ql1
         Ql1 = Ql
      end do

      ! non-axisymmetric (l,m) moments for 1 <= m <= mmax, m <= l <= lmax.
      do m = 1, mmax

         m2s = lm(0, 2*m-1)
         m2c = lm(0, 2*m)
         ! The (m,m) moment
         ! Associated Legendre polynomial: P_m^m = (-1)^m (2m-1)!! (1-x^2)^{m/2}
         ! The (2m-1)!! factor is integrated in ofact(:) array, where it mostly cancels out, note that (2m-1)!! \simeq m! exp(m/sqrt(2)) so it grows pretty fast with m
         Ql1 = sin_th ** m
         Q(m2s+m, INSIDE,  ir) = Q(m2s+m, INSIDE,  ir) +  rn(m) * Ql1 * sfac(m)
         Q(m2c+m, INSIDE,  ir) = Q(m2c+m, INSIDE,  ir) +  rn(m) * Ql1 * cfac(m)
         Q(m2s+m, OUTSIDE, ir) = Q(m2s+m, OUTSIDE, ir) + irn(m) * Ql1 * sfac(m)
         Q(m2c+m, OUTSIDE, ir) = Q(m2c+m, OUTSIDE, ir) + irn(m) * Ql1 * cfac(m)

         ! BEWARE: most of computational cost of multipoles is here
         ! from (m+1,m) to (lmax,m)
         ! Associated Legendre polynomial: P_{m+1}^m = x (2m+1) P_m^m
         ! Associated Legendre polynomial recurrence: (l-m) P_l^m = x (2l-1) P_{l-1}^m - (l+m-1) P_{l-2}^m
         Ql2 = 0.
         do l = m + 1, lmax
            Ql = cos_th * k12(1, l, m) * Ql1 - k12(2, l, m) * Ql2
            Q(m2s+l, INSIDE,  ir) = Q(m2s+l, INSIDE,  ir) +  rn(l) * Ql * sfac(m)
            Q(m2c+l, INSIDE,  ir) = Q(m2c+l, INSIDE,  ir) +  rn(l) * Ql * cfac(m)
            Q(m2s+l, OUTSIDE, ir) = Q(m2s+l, OUTSIDE, ir) + irn(l) * Ql * sfac(m)
            Q(m2c+l, OUTSIDE, ir) = Q(m2c+l, OUTSIDE, ir) + irn(l) * Ql * cfac(m)
            Ql2 = Ql1
            Ql1 = Ql
         end do

      enddo

   end subroutine point2moments

!!$ ============================================================================
!!
!! Compute infinite-boundary potential from multipole moments
!!
!! \todo distribute excess of work more evenly (important only for large number of PEs, ticket:43)
!!

   subroutine moments2bnd_potential

      use multigridvars,   only: is_external, XLO, XHI, YLO, YHI, ZLO, ZHI, XDIR, YDIR, ZDIR, LOW, HIGH

      implicit none

      integer :: i, j, k

      if (is_external(XLO) .or. is_external(XHI)) then
         do j = lmpole%js, lmpole%je
            do k = lmpole%ks, lmpole%ke
               if (is_external(XLO)) call moments2pot(lmpole%bnd_x(j, k, LOW),  fbnd_x(LOW) -CoM(XDIR), lmpole%y(j)-CoM(YDIR),  lmpole%z(k)-CoM(ZDIR))
               if (is_external(XHI)) call moments2pot(lmpole%bnd_x(j, k, HIGH), fbnd_x(HIGH)-CoM(XDIR), lmpole%y(j)-CoM(YDIR),  lmpole%z(k)-CoM(ZDIR))
            end do
         end do
      end if

      if (is_external(YLO) .or. is_external(YHI)) then
         do i = lmpole%is, lmpole%ie
            do k = lmpole%ks, lmpole%ke
               if (is_external(YLO)) call moments2pot(lmpole%bnd_y(i, k, LOW),  lmpole%x(i)-CoM(XDIR),  fbnd_y(LOW) -CoM(YDIR), lmpole%z(k)-CoM(ZDIR))
               if (is_external(YHI)) call moments2pot(lmpole%bnd_y(i, k, HIGH), lmpole%x(i)-CoM(XDIR),  fbnd_y(HIGH)-CoM(YDIR), lmpole%z(k)-CoM(ZDIR))
            end do
         end do
      end if

      if (is_external(ZLO) .or. is_external(ZHI)) then
         do i = lmpole%is, lmpole%ie
            do j = lmpole%js, lmpole%je
               if (is_external(ZLO)) call moments2pot(lmpole%bnd_z(i, j, LOW),  lmpole%x(i)-CoM(XDIR),  lmpole%y(j)-CoM(YDIR),  fbnd_z(LOW) -CoM(ZDIR))
               if (is_external(ZHI)) call moments2pot(lmpole%bnd_z(i, j, HIGH), lmpole%x(i)-CoM(XDIR),  lmpole%y(j)-CoM(YDIR),  fbnd_z(HIGH)-CoM(ZDIR))
            end do
         end do
      end if

   end subroutine moments2bnd_potential

!!$ ============================================================================
!!
!! Compute potential from multipole moments at a single point
!!
!! \todo improve accuracy with linear interpolaion over radius
!!

   subroutine moments2pot(potential, x, y, z)

      use constants, only: newtong

      implicit none

      real, intent(in)  :: x, y, z   !< coordinates of the point
      real, intent(out) :: potential !< calculated potential at given point

      real :: sin_th, cos_th
      real :: Ql, Ql1, Ql2
      integer :: l, m, ir, m2s, m2c

      call geomfac4moments(-newtong, x, y, z, sin_th, cos_th, ir)

      ! monopole, the (0,0) moment; P_0 = 1.
      potential = &
           Q(0, INSIDE,  ir)   * irn(0) + &
           Q(0, OUTSIDE, ir+1) *  rn(0)
      ! ir+1 to prevent duplicate accounting contributions from ir bin; alternatively one can modify radial integration

      ! axisymmetric (l,0) moments
      Ql2 = 0.
      Ql1 = 1.
      do l = 1, lmax
         Ql = cos_th * k12(1, l, 0) * Ql1 - k12(2, l, 0) * Ql2
         potential = potential + Ql * ( &
              &      Q(l, INSIDE,  ir)   * irn(l) + &
              &      Q(l, OUTSIDE, ir+1) *  rn(l) )
         Ql2 = Ql1
         Ql1 = Ql
      end do

      ! non-axisymmetric (l,m) moments for 1 <= m <= mmax, m <= l <= lmax.
      do m = 1, mmax
         m2s = lm(0, 2*m-1)
         m2c = lm(0, 2*m)
         ! The (m,m) moment
         Ql1 = sin_th ** m
         potential = potential + Ql1 * ( &
              &      (Q(m2c+m,   INSIDE,  ir)   * irn(m)             + &
              &       Q(m2c+m,   OUTSIDE, ir+1) *  rn(m) ) * cfac(m) + &
              &      (Q(m2s+m, INSIDE,  ir)   * irn(m)             + &
              &       Q(m2s+m, OUTSIDE, ir+1) *  rn(m) ) * sfac(m) )

         ! BEWARE: lots of computational cost of multipoles is here
         ! from (m+1,m) to (lmax,m)
         Ql2 = 0.
         do l = m+1, lmax
            Ql = cos_th * k12(1, l, m) * Ql1 - k12(2, l, m) * Ql2
            potential = potential + Ql * ( &
                 &      (Q(m2c+l, INSIDE,  ir)   * irn(l)             + &
                 &       Q(m2c+l, OUTSIDE, ir+1) *  rn(l) ) * cfac(m) + &
                 &      (Q(m2s+l, INSIDE,  ir)   * irn(l)             + &
                 &       Q(m2s+l, OUTSIDE, ir+1) *  rn(l) ) * sfac(m) )
            Ql2 = Ql1
            Ql1 = Ql
         end do
      end do

   end subroutine moments2pot

!!$ ============================================================================
!!
!! This routine calculates various geometrical numbers required for multipole evaluation
!! It modifies rn(:), irn(:), cfac(:) and sfac(:) arrays. Scalars are passed through argument list.
!!
!! \todo return also fraction of the radial bin, (r/drq - ir) for radial interpolation
!!

   subroutine geomfac4moments(factor, x, y, z, sin_th, cos_th, ir)

      use errh, only: die

      implicit none

      real,    intent(in)  :: factor         !< scaling factor (e.g. mass) of the contributing point
      real,    intent(in)  :: x, y, z        !< coordinates of the contributing point
      real,    intent(out) :: sin_th, cos_th !< sine and cosine of the theta angle
      integer, intent(out) :: ir             !< radial index for the Q(:, :, r) array

      real    :: rxy, r, rinv
      real    :: sin_ph, cos_ph
      integer :: l, m

      ! radius and its projection onto XY plane
      rxy  = x**2 + y**2
      r    = sqrt(rxy + z**2)
      rxy  = sqrt(rxy)
      if (r /= 0.) then
         rinv = 1. / r
      else
         rinv = 0.
      end if

      !radial index for the Q(:, :, r) array
      ir = int(r / drq)
      if (ir > rqbin .or. ir < 0) call die("[multipole:geomfac4moments] radial index outside Q(:, :, r) range")
      irmax = max(irmax, ir)
      irmin = min(irmin, ir)

      ! azimuthal angle sine and cosine tables
      ! ph = atan2(y, x); cfac(m) = cos(m * ph); sfac(m) = sin(m * ph)
      ! cfac(0) and sfac(0) are set in init_multigrid
      if (rxy /= 0.) then
         cos_ph = x / rxy
         sin_ph = y / rxy
      else
         cos_ph = 1.
         sin_ph = 0.
      end if
      ! \todo Possible optimization: number of computed elements can be doubled on each loop iteration (should give better pipelining)
      do m = 1, mmax
         cfac(m) = cos_ph*cfac(m-1) - sin_ph*sfac(m-1)
         sfac(m) = cos_ph*sfac(m-1) + sin_ph*cfac(m-1)
      end do

      ! vertical angle
      cos_th = z   * rinv
      sin_th = rxy * rinv

      ! \todo check how much it would degrade solution and improve performance to move the multiplications by rn(:) and irn(:) to img_mass2moments (before and after integration)
      ! rn(l) = factor * r ** l; irn(l) = factor * r ** -(l+1)
      rn(0)  = factor
      irn(0) = factor * rinv
      ! scale r to reduce risk of occuring a Floating Point Overflow for high l_max and large r values. Leave the original value of r only in irn(0)
      r = r / rscale
      rinv = rinv * rscale
      ! \todo Possible optimization: number of computed elements can be doubled on each loop iteration (should give better pipelining)
      do l = 1, lmax
         rn(l)  =  rn(l-1) * r
         irn(l) = irn(l-1) * rinv
      end do

   end subroutine geomfac4moments

!!$ ============================================================================
!!
!! HEAVY_DEBUG marks routines that normally are never called, but at some point
!! were useful to test correctness or something.
!!

!#define HEAVY_DEBUG
#ifdef HEAVY_DEBUG

!!$ ============================================================================
!!
!! Quick test for correctness of the multipole solver.
!! cphi should agree well with phi with largest errors at r = sqrt(sum(p(1:3)**2))
!!

#error "The test_multipoles routine is outdated and was commented out"

   subroutine test_multipoles

      use errh,               only: die
      use constants,          only: newtong
      use multigridhelpers,   only: aux_par_R0, aux_par_R1, aux_par_R2, aux_par_I0, dirty_debug

      implicit none

      integer :: i, j, r, rr
      real :: phi, cphi
      real, dimension(0:3) :: p, x

!!$      ! reset the multipole data
!!$      Q(:, :, :, :) = 0.
!!$
!!$      irmax = 0
!!$      irmin = rqbin
!!$
!!$      p(0) = 1./newtong
!!$      p(1) = aux_par_R0
!!$      p(2) = aux_par_R1
!!$      p(3) = aux_par_R2
!!$
!!$      call point2moments(p(0), p(1), p(2), p(3))
!!$
!!$      ! integrate radially and apply normalization factor (the (4 \pi)/(2 l  + 1) terms cancel out)
!!$      irmin = 0
!!$      irmax = rqbin-1
!!$      rr = max(1, irmin)
!!$
!!$      Q(1:lmax, :, INSIDE, rr-1) = Q(1:lmax, :, INSIDE, rr-1) * ofact(1:lmax, :)
!!$      do r = rr, irmax
!!$         Q(0,      :, INSIDE, r) = Q(0,      :, INSIDE, r)                    + Q(0,      :, INSIDE, r-1)
!!$         Q(1:lmax, :, INSIDE, r) = Q(1:lmax, :, INSIDE, r) * ofact(1:lmax, :) + Q(1:lmax, :, INSIDE, r-1)
!!$      end do
!!$
!!$      Q(1:lmax, :, OUTSIDE, irmax+1) = Q(1:lmax, :, OUTSIDE, irmax+1) * ofact(1:lmax, :)
!!$      do r = irmax, rr, -1
!!$         Q(0,      :, OUTSIDE, r) = Q(0,      :, OUTSIDE, r)                    + Q(0,      :, OUTSIDE, r+1)
!!$         Q(1:lmax, :, OUTSIDE, r) = Q(1:lmax, :, OUTSIDE, r) * ofact(1:lmax, :) + Q(1:lmax, :, OUTSIDE, r+1)
!!$      end do
!!$
!!$      do i = -60, 60
!!$         do j = -60, 60
!!$            x(1) = 0.05 * i
!!$            x(2) = 0.05 * j
!!$            x(3) = 0.1 * aux_par_I0
!!$
!!$            cphi = - p(0)*newtong / sqrt(1e-290+sum((p(1:3)-x(1:3))**2))
!!$
!!$            call moments2pot(phi, x(1), x(2), x(3))
!!$            write(*,'(a,3f7.2,2(a,g15.5))')" xyz= ",x(1:3)," phi = ",phi," calc= ",cphi
!!$         end do
!!$      end do
!!$
!!$      call die("qniec")

   end subroutine test_multipoles

#endif /* HEAVY_DEBUG */

end module multipole
