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
!! \brief Multipole solver for isolated boundaries
!!
!! \details This solver estimates gravitational potential on external (domain) boundaries, which allows to mimic \f$\Phi(\infty) = 0\f$.
!! The calculated potential may then be used in second pass of the Poisson solver (this time run with an empty space) to calculate correction
!! to be added to the first-pass solution obtained with homogenous Dirichlet boundary conditions.
!<

module multipole
! pulled by MULTIGRID && GRAV
   ! needed for global vars in this module
   use multigridvars, only: plvl
   use constants,     only: ndims

   implicit none

   private
   public :: init_multipole, cleanup_multipole, multipole_solver
   public :: lmax, mmax, ord_prolong_mpole, coarsen_multipole, use_point_monopole, interp_pt2mom, interp_mom2pot ! initialized in multigrid_gravity

   integer, parameter        :: INSIDE = 1                       !< index for interior multipole expansion
   integer, parameter        :: OUTSIDE = INSIDE + 1             !< index for exterior multipole expansion

   ! namelist parameters for MULTIGRID_GRAVITY
   integer                   :: lmax                             !< Maximum l-order of multipole moments
   integer                   :: mmax                             !< Maximum m-order of multipole moments. Equal to lmax by default.
   integer                   :: ord_prolong_mpole                !< boundary prolongation operator order; allowed values are -2 .. 2
   integer                   :: coarsen_multipole                !< If > 0 then evaluate multipoles at roof%level-coarsen_multipole level
   logical                   :: use_point_monopole               !< Don't evaluate multipole moments, use point-like mass approximation (crudest possible)
   logical                   :: interp_pt2mom                    !< Distribute contribution from a cell between two adjacent radial bins (linear interpolation in radius)
   logical                   :: interp_mom2pot                   !< Compute the potential from moments from two adjacent radial bins (linear interpolation in radius)

   ! radial discretization
   integer                   :: rqbin                            !< number of radial samples of multipoles
   real                      :: drq                              !< radial resolution of multipoles
   real                      :: rscale                           !< scaling factor that limits risk of floating point overflow for high powers of radius (rn(:) and irn(:))
   integer                   :: irmin                            !< minimum Q(:, :, r) index in use
   integer                   :: irmax                            !< maximum Q(:, :, r) index in use

   type(plvl), pointer       :: lmpole                           !< pointer to the level where multipoles are evaluated
   real, dimension(0:ndims)  :: CoM                              !< Total mass and center of mass coordinates
   logical                   :: zaxis_inside                     !< true when z-axis belongs to the inner radial boundary in polar coordinates

   ! multipoles and auxiliary factors
   real, dimension(:,:,:), allocatable :: Q                      !< The whole moment array with dependence on radius
   real, dimension(:,:,:), allocatable :: k12                    !< array of Legendre recurrence factors
   real, dimension(:),     allocatable :: sfac                   !< 0..m_max sized array of azimuthal sine factors
   real, dimension(:),     allocatable :: cfac                   !< 0..m_max sized array of azimuthal cosine factors
   real, dimension(:),     allocatable :: ofact                  !< arrays of Legendre normalization factor (compressed)
   real, dimension(:),     allocatable :: rn                     !< 0..l_max sized array of positive powers of r
   real, dimension(:),     allocatable :: irn                    !< 0..l_max sized array of negative powers of r

contains

!!$ ============================================================================
!>
!! \brief This function returns array index for "compressed" Q(:, inout, r) array replacing "plain" Q(l, m, inout, r) array
!! \details Originally there were separate indices for l and m multipole numbers in Q and ofact arrays. Since there are no Y_{l,m} harmonics for l<|m|, approximately half of the Q array
!! was left unused. For each m there is only lmax-|m|+1 valid Y_{l,m} harmonic contributions (for l = |m| .. lmax). The function lm(l, m) converts each valid (l,m) pair into an
!! unique index leaving no unused entries in the Q array.
!! Assume that real Y_{l,m} harmonic was stored in Q(l, 2*|m|, :, :) for even (cosine, positive-m) harmonic, or in Q(l, 2*|m|-1, :, :) for odd (sine, negative-m) harmonic.
!! Then consecutive entries in compressed Q(:, inout, r) array should be related to Y_{l,m} harmonics, starting from index 0 as follows
!! Y_{0,0}, Y_{1,0}, ..., Y_{lmax, 0}, Y_{1,-1}, Y_{2,-1}, ..., Y_{lmax,-1}, Y_{1,1}, Y_{2,1}, ..., Y_{lmax,1}, Y_{2,-2}, Y_{3,-2}, ..., Y_{lmax,-2},
!! Y_{2,2}, Y_{3,2}, ..., Y_{lmax,2}, ...,, Y(lmax-1,-(mmax-1)), Y(lmax,-(mmax-1)), Y(lmax-1,mmax-1), Y(lmax,mmax-1), Y(lmax,-mmax), Y(lmax,mmax)
!! Does it looks a bit cryptic? I agree.
!<

   elemental integer function lm(l, m)

      implicit none

      integer, intent(in) :: l, m

      lm = l + m*lmax - int((m-1)/2)*int(m/2)

   end function lm

!!$ ============================================================================
!>
!! \brief Initialization routine, called from init_multigrid
!<

   subroutine init_multipole(mb_alloc)

      use constants,     only: small, pi, xdim, ydim, zdim, ndims, GEO_XYZ, GEO_RPZ, LO, HI
      use dataio_pub,    only: die, warn
      use domain,        only: dom, eff_dim, geometry_type
      use mpisetup,      only: master
      use multigridvars, only: roof

      implicit none

      real,                 intent(inout) :: mb_alloc               !< multigrid allocation counter

      integer, dimension(4) :: aerr
      integer               :: l,m

      ! assume that Center of Mass is approximately in the center of computational domain by default
      CoM(0) = 1.

      select case (geometry_type)
         case (GEO_XYZ)
            CoM(xdim:zdim) = dom%C_(xdim:zdim)
            zaxis_inside = .false.
         case (GEO_RPZ)
            if (dom%L_(ydim) >= (2.-small)*pi) then
               CoM(xdim) = 0.
               CoM(ydim) = 0.
            else
!!$               CoM(xdim) = 2./3. * (dom%edge(xdim, HI)**3-dom%edge(xdim, LO)**3)/(dom%edge(xdim, HI)**2-dom%edge(xdim, LO)**2)
!!$               if (dom%L_(ydim) /= 0.) CoM(xdim) = CoM(xdim) * sin(dom%L_(ydim)/2.)/(dom%L_(ydim)/2.)
!!$               CoM(ydim) = dom%C_(ydim)
               CoM(xdim) = 0.
               CoM(ydim) = 0.
            endif
            CoM(zdim) = dom%C_(zdim)
            zaxis_inside = dom%edge(xdim, LO) <= dom%L_(xdim)/dom%n_d(xdim)
            if (master) then
               if (zaxis_inside) call warn("[multipole:init_multipole] Setups with Z-axis at the edge of the domain may not work as expected yet.")
               if (use_point_monopole) call warn("[multipole:init_multipole] Point-like monopole is not implemented.")
            endif
            use_point_monopole = .false.
         case default
            call die("[multipole:init_multipole] Unsupported geometry.")
      end select

      if (eff_dim /= ndims) call die("[multipole:init_multipole] Only 3D is supported") !> \todo add support for 2D RZ

      !multipole moments
      if (mmax > lmax) then
         if (master) call warn("[multipole:init_multipole] mmax reduced to lmax")
         mmax = lmax
      endif
      if (mmax < 0) mmax = lmax

      lmpole => roof
      do l = 1, coarsen_multipole
         if (associated(lmpole%coarser)) then
            lmpole => lmpole%coarser
         else
            if (master) call warn("[multipole:init_multipole] too deep multipole coarsening.")
         endif
      enddo

      if (coarsen_multipole > 0) then
         if (interp_pt2mom) then
            call warn("[multipole:init_multipole] coarsen_multipole > 0 disables interp_pt2mom.")
            interp_pt2mom = .false.
         endif
         if (interp_mom2pot) then
            call warn("[multipole:init_multipole] coarsen_multipole > 0 disables interp_mom2pot.")
            interp_mom2pot = .false.
         endif
      endif

      if (.not. use_point_monopole) then
         if (allocated(rn) .or. allocated(irn) .or. allocated(sfac) .or. allocated(cfac)) call die("[multipole:init_multipole] rn, irn, sfac or cfac already allocated")
         allocate(  rn(0:lmax), stat=aerr(1))
         allocate( irn(0:lmax), stat=aerr(2))
         allocate(sfac(0:mmax), stat=aerr(3))
         allocate(cfac(0:mmax), stat=aerr(4))
         if (any(aerr(1:4) /= 0)) call die("[multipole:init_multipole] Allocation error: rn, irn sfac or cfac")
         mb_alloc = mb_alloc + size(rn) + size(irn) + size(sfac) + size(cfac)

         select case (geometry_type)
            case (GEO_XYZ)
               drq = min(lmpole%dx, lmpole%dy, lmpole%dz) / 2.
               rqbin = int(sqrt(sum(dom%L_(:)**2))/drq) + 1
               ! arithmetic average of the closest and farthest points of computational domain with respect to its center
               !>
               !!\todo check what happens if there are points that are really close to the domain center (maybe we should use a harmonic average?)
               !! Issue a warning or error if it is known that given lmax leads to FP overflows in rn(:) and irn(:)
               !<
               rscale = ( minval(dom%L_(:)) + sqrt(sum(dom%L_(:)**2)) )/4.
            case (GEO_RPZ)
               drq = min(lmpole%dx, dom%C_(xdim)*lmpole%dy, lmpole%dz) / 2.
               rqbin = int(sqrt((2.*dom%edge(xdim, HI))**2 + dom%L_(zdim)**2)/drq) + 1
               rscale = ( min(2.*dom%edge(xdim, HI), dom%L_(zdim)) + sqrt((2.*dom%edge(xdim, HI))**2 + dom%L_(zdim)**2) )/4.
            case default
               call die("[multipole:init_multipole] Unsupported geometry.")
         end select
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
                  ofact(lm(l, 2*m))   = ofact(lm(l, 2*(m-1))) * real(1 - 2*m)**2 / real((l+m) * (l-m+1))
                  ofact(lm(l, 2*m-1)) = ofact(lm(l, 2*m))
               endif
               if (l>m) then
                  k12(1, l, m) = real(2 * l - 1) / real(l - m)
                  k12(2, l, m) = real(l + m - 1) / real(l - m)
               endif
            enddo
         enddo
         ofact(0:lmax) = 1. ! lm(0:lmax,0)

         cfac(0) = 1.e0
         sfac(0) = 0.e0

      endif

   end subroutine init_multipole

!!$ ============================================================================
!>
!! \brief Multipole cleanup
!<

   subroutine cleanup_multipole

      implicit none

      if (allocated(rn))    deallocate(rn)
      if (allocated(irn))   deallocate(irn)
      if (allocated(sfac))  deallocate(sfac)
      if (allocated(cfac))  deallocate(cfac)
      if (allocated(Q))     deallocate(Q)
      if (allocated(k12))   deallocate(k12)
      if (allocated(ofact)) deallocate(ofact)

   end subroutine cleanup_multipole

!!$ ============================================================================
!>
!! \brief Multipole solver
!!
!! \todo improve multipole expansion on coarser grids
!! (see. "A Scalable Parallel Poisson Solver in Three Dimensions with Infinite-Domain Boundary Conditions" by McCorquodale, Colella, Balls and Baden).
!! Coarsening by one level would reduce the multipole costs by a factor of 4.
!<

   subroutine multipole_solver

      use multigridvars,      only: roof, solution, lvl, plvl
      use multigridhelpers,   only: dirtyH, dirty_debug
      use multigridbasefuncs, only: zero_boundaries

      implicit none

      type(plvl), pointer :: curl

      if (dirty_debug) then
         lmpole%bnd_x(:, :, :) = dirtyH
         lmpole%bnd_y(:, :, :) = dirtyH
         lmpole%bnd_z(:, :, :) = dirtyH
      else
         call zero_boundaries(lmpole)
      endif

      if (.not. associated(lmpole, roof)) then
         curl => roof
         do while (associated(curl) .and. .not. associated(curl, lmpole)) ! do lev = roof%level, lmpole%level + 1, -1
            call curl%restrict_level(solution)  ! Overkill, only some layers next to external boundary are needed.
            curl => curl%coarser
         enddo                                ! An alternative: do potential2img_mass on the roof and restrict bnd_[xyz] data.
      endif
      call potential2img_mass

      if (use_point_monopole) then
         call find_img_CoM
         call isolated_monopole
      else
         ! results seems to be slightly better without find_img_CoM.
         ! With CoM or when it is known than CoM is close to the domain center one may try to save some CPU time by lowering mmax.
         call img_mass2moments
         call moments2bnd_potential
      endif

      if (.not. associated(lmpole, roof)) then
         curl => lmpole
         do while (associated(curl) .and. .not. associated(curl, roof)) ! do lev = lmpole%level, roof%level - 1
            call prolong_ext_bnd(curl)
            curl => curl%finer
         enddo
      endif

   end subroutine multipole_solver

!!$ ============================================================================
!>
!! \brief Set boundary potential from monopole source. Fill lmpole%bnd_[xyz] arrays with expected values of the gravitational potential at external face of computational domain.
!! \details This is a simplified approach that can be used for tests and as a fast replacement for the
!! multipole boundary solver for nearly spherically symmetric source distributions.
!! The isolated_monopole subroutine ignores the radial profile of the monopole
!<

   subroutine isolated_monopole

      use dataio_pub,    only: die
      use constants,     only: xdim, ydim, zdim, LO, HI, GEO_XYZ !, GEO_RPZ
      use multigridvars, only: is_external
      use domain,        only: geometry_type
      use units,         only: newtong

      implicit none

      integer :: i, j, k, lh
      real    :: r2

      if (geometry_type /= GEO_XYZ) call die("[multigridmultipole:isolated_monopole] non-cartesian geometry not implemented yet")

      do lh = LO, HI
         if (is_external(xdim, lh)) then
            do j = lmpole%js, lmpole%je
               do k = lmpole%ks, lmpole%ke
                  r2 = (lmpole%y(j)-CoM(ydim))**2 + (lmpole%z(k) - CoM(zdim))**2
                  lmpole%bnd_x(j, k, lh) = - newtong * CoM(0) / sqrt(r2 + (lmpole%fbnd(xdim, lh)-CoM(xdim))**2)
               enddo
            enddo
         endif
         if (is_external(ydim, lh)) then
            do i = lmpole%is, lmpole%ie
               do k = lmpole%ks, lmpole%ke
                  r2 = (lmpole%x(i)-CoM(xdim))**2 + (lmpole%z(k) - CoM(zdim))**2
                  lmpole%bnd_y(i, k, lh) = - newtong * CoM(0) / sqrt(r2 + (lmpole%fbnd(ydim, lh)-CoM(ydim))**2)
               enddo
            enddo
         endif
         if (is_external(zdim, lh)) then
            do i = lmpole%is, lmpole%ie
               do j = lmpole%js, lmpole%je
                  r2 = (lmpole%x(i)-CoM(xdim))**2 + (lmpole%y(j) - CoM(ydim))**2
                  lmpole%bnd_z(i, j, lh) = - newtong * CoM(0) / sqrt(r2 + (lmpole%fbnd(zdim, lh)-CoM(zdim))**2)
               enddo
            enddo
         endif
      enddo

   end subroutine isolated_monopole

!!$ ============================================================================
!>
!! \brief Find total mass and its center
!! \details This routine does the summation only on external boundaries
!<

   subroutine find_img_CoM

      use constants,     only: ndims, xdim, ydim, zdim, LO, HI, GEO_XYZ !, GEO_RPZ
      use dataio_pub,    only: die
      use domain,        only: geometry_type
      use mpi,           only: MPI_DOUBLE_PRECISION, MPI_SUM
      use mpisetup,      only: comm, ierr
      use multigridvars, only: is_external

      implicit none

      real, dimension(0:ndims) :: lsum, dsum
      integer :: lh

      if (geometry_type /= GEO_XYZ) call die("[multigridmultipole:find_img_CoM] non-cartesian geometry not implemented yet")

      lsum(:) = 0.

      do lh = LO, HI
         if (is_external(xdim, lh)) then
            dsum(0)    =      sum( lmpole%bnd_x(lmpole%js:lmpole%je, lmpole%ks:lmpole%ke, lh) )
            dsum(xdim) = dsum(0) * lmpole%fbnd(xdim, lh)
            dsum(ydim) = sum( sum( lmpole%bnd_x(lmpole%js:lmpole%je, lmpole%ks:lmpole%ke, lh),  dim=2) * lmpole%y(lmpole%js:lmpole%je) )
            dsum(zdim) = sum( sum( lmpole%bnd_x(lmpole%js:lmpole%je, lmpole%ks:lmpole%ke, lh),  dim=1) * lmpole%z(lmpole%ks:lmpole%ke) )
            lsum(:)    = lsum(:) + dsum(:) * lmpole%dyz
         endif
         if (is_external(ydim, lh)) then
            dsum(0)    =      sum( lmpole%bnd_y(lmpole%is:lmpole%ie, lmpole%ks:lmpole%ke, lh) )
            dsum(xdim) = sum( sum( lmpole%bnd_y(lmpole%is:lmpole%ie, lmpole%ks:lmpole%ke, lh),  dim=2) * lmpole%x(lmpole%is:lmpole%ie) )
            dsum(ydim) = dsum(0) * lmpole%fbnd(ydim, lh)
            dsum(zdim) = sum( sum( lmpole%bnd_y(lmpole%is:lmpole%ie, lmpole%ks:lmpole%ke, lh),  dim=1) * lmpole%z(lmpole%ks:lmpole%ke) )
            lsum(:)    = lsum(:) + dsum(:) * lmpole%dxz
         endif
         if (is_external(zdim, lh)) then
            dsum(0)    =      sum( lmpole%bnd_z(lmpole%is:lmpole%ie, lmpole%js:lmpole%je, lh) )
            dsum(xdim) = sum( sum( lmpole%bnd_z(lmpole%is:lmpole%ie, lmpole%js:lmpole%je, lh),  dim=2) * lmpole%x(lmpole%is:lmpole%ie) )
            dsum(ydim) = sum( sum( lmpole%bnd_z(lmpole%is:lmpole%ie, lmpole%js:lmpole%je, lh),  dim=1) * lmpole%z(lmpole%ks:lmpole%ke) )
            dsum(zdim) = dsum(0) * lmpole%fbnd(zdim, lh)
            lsum(:)    = lsum(:) + dsum(:) * lmpole%dxy
         endif
      enddo

      call MPI_Allreduce(lsum(0:ndims), CoM(0:ndims), 4, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)

      if (CoM(0) /= 0.) then
         CoM(xdim:zdim) = CoM(xdim:zdim) / CoM(0)
      else
         call die("[multipole:find_img_CoM] Total mass == 0")
      endif

   end subroutine find_img_CoM

!!$ ============================================================================
!>
!! \brief Convert potential into image mass. This way we reduce a 3D problem to a 2D one.
!! \details There will be work imbalance here because different PEs may operate on different amount of external boundary data
!<

   subroutine potential2img_mass

      use constants,     only: GEO_RPZ, LO, HI, xdim, ydim, zdim
      use domain,        only: geometry_type
      use multigridvars, only: is_external, solution

      implicit none

      integer :: i
      real, parameter :: a1 = -2., a2 = (-2. - a1)/3. ! interpolation parameters;   <---- a1=-2 => a2=0
      ! a1 = -2. is the simplest, low order choice, gives best agreement of total mass and CoM location when compared to 3-D integration
      ! a1 = -1., a2 = -1./3. seems to do the best job,
      !> \todo: find out how and why

      !> \deprecated BEWARE: some cylindrical factors may be helpful
      if (is_external(xdim, LO)) then
         if (zaxis_inside .and. geometry_type == GEO_RPZ) then
            lmpole%bnd_x(                         lmpole%js:lmpole%je, lmpole%ks:lmpole%ke, LO) = 0. ! treat as internal
         else
            lmpole%bnd_x(                         lmpole%js:lmpole%je, lmpole%ks:lmpole%ke, LO) =   ( &
                 & a1 * lmpole%mgvar(lmpole%is,   lmpole%js:lmpole%je, lmpole%ks:lmpole%ke, solution) + &
                 & a2 * lmpole%mgvar(lmpole%is+1, lmpole%js:lmpole%je, lmpole%ks:lmpole%ke, solution) ) / lmpole%dx
         endif
      endif

      if (is_external(xdim, HI)) lmpole%bnd_x(lmpole%js:lmpole%je, lmpole%ks:lmpole%ke, HI) =   ( &
           &   a1 * lmpole%mgvar(lmpole%ie,   lmpole%js:lmpole%je, lmpole%ks:lmpole%ke, solution) + &
           &   a2 * lmpole%mgvar(lmpole%ie-1, lmpole%js:lmpole%je, lmpole%ks:lmpole%ke, solution) ) / lmpole%dx

      if (is_external(ydim, LO)) then
         lmpole%bnd_y(            lmpole%is:lmpole%ie,              lmpole%ks:lmpole%ke, LO) =    ( &
              & a1 * lmpole%mgvar(lmpole%is:lmpole%ie, lmpole%js,   lmpole%ks:lmpole%ke, solution) + &
              & a2 * lmpole%mgvar(lmpole%is:lmpole%ie, lmpole%js+1, lmpole%ks:lmpole%ke, solution) ) / lmpole%dy
         if (geometry_type == GEO_RPZ) then
            do i = lmpole%is, lmpole%ie
               if (lmpole%x(i) /= 0.) then ! cylindrical factor in the gradient operator
                  lmpole%bnd_y(i, lmpole%ks:lmpole%ke, LO) = lmpole%bnd_y(i, lmpole%ks:lmpole%ke, LO) / lmpole%x(i) !> \todo precompute sanitized 1/x(:) (convert plvl% into cg%)
               else
                  lmpole%bnd_y(i, lmpole%ks:lmpole%ke, LO) = 0.
               endif
            enddo
         endif
      endif

      if (is_external(ydim, HI)) then
         lmpole%bnd_y(            lmpole%is:lmpole%ie,              lmpole%ks:lmpole%ke, HI) =   ( &
              & a1 * lmpole%mgvar(lmpole%is:lmpole%ie, lmpole%je,   lmpole%ks:lmpole%ke, solution) + &
              & a2 * lmpole%mgvar(lmpole%is:lmpole%ie, lmpole%je-1, lmpole%ks:lmpole%ke, solution) ) / lmpole%dy
         if (geometry_type == GEO_RPZ) then
            do i = lmpole%is, lmpole%ie
               if (lmpole%x(i) /= 0.) then
                  lmpole%bnd_y(i, lmpole%ks:lmpole%ke, HI) = lmpole%bnd_y(i, lmpole%ks:lmpole%ke, HI) / lmpole%x(i)
               else
                  lmpole%bnd_y(i, lmpole%ks:lmpole%ke, HI) = 0.
               endif
            enddo
         endif
      endif

      if (is_external(zdim, LO)) lmpole%bnd_z(lmpole%is:lmpole%ie, lmpole%js:lmpole%je,              LO) =    ( &
           &                a1 * lmpole%mgvar(lmpole%is:lmpole%ie, lmpole%js:lmpole%je, lmpole%ks,   solution) + &
           &                a2 * lmpole%mgvar(lmpole%is:lmpole%ie, lmpole%js:lmpole%je, lmpole%ks+1, solution) ) / lmpole%dz

      if (is_external(zdim, HI)) lmpole%bnd_z(lmpole%is:lmpole%ie, lmpole%js:lmpole%je,              HI) =   ( &
           &                a1 * lmpole%mgvar(lmpole%is:lmpole%ie, lmpole%js:lmpole%je, lmpole%ke,   solution) + &
           &                a2 * lmpole%mgvar(lmpole%is:lmpole%ie, lmpole%js:lmpole%je, lmpole%ke-1, solution) ) / lmpole%dz

   end subroutine potential2img_mass

!!$ ============================================================================
!>
!! \brief Prolong boundaries wrapper
!<

   subroutine prolong_ext_bnd(coarse)

      use constants,     only: ndims
      use dataio_pub,    only: die
      use domain,        only: eff_dim
      use multigridvars, only: is_external, is_mg_uneven, plvl

      implicit none

      type(plvl), pointer, intent(in) :: coarse !< level to prolong from

      if (is_mg_uneven) call die("[multigridmultipole:prolong_ext_bnd0] uneven decomposition not implemented yet")
      if (eff_dim<ndims) call die("[multigridmultipole:prolong_ext_bnd0] 1D and 2D not finished")

      if (abs(ord_prolong_mpole) > 2) call die("[multipole:prolong_ext_bnd] interpolation order too high")

      !> \deprecated BEWARE: do we need cylindrical factors for prolongation?
      if (any(is_external(:, :))) then
         if (ord_prolong_mpole == 0) then
            call prolong_ext_bnd0(coarse)
         else
            call prolong_ext_bnd2(coarse)
         endif
      endif

   end subroutine prolong_ext_bnd

!!$ ============================================================================
!>
!! \brief Prolong boundaries by injection.
!<

   subroutine prolong_ext_bnd0(coarse)

      use multigridvars,   only: plvl, is_external
      use constants,       only: HI, LO, xdim, ydim, zdim
      use dataio_pub,      only: die

      implicit none

      type(plvl), pointer, intent(in) :: coarse !< level to prolong from

      type(plvl), pointer :: fine
      integer :: lh

      if (.not. associated(coarse)) call die("[multigridmultipole:prolong_ext_bnd0] coarse == null()")
      fine   => coarse%finer
      if (.not. associated(fine)) call die("[multigridmultipole:prolong_ext_bnd0] fine == null()")

      do lh = LO, HI
         if (is_external(xdim, lh)) then
            fine%bnd_x(fine%js  :fine%je-1:2, fine%ks  :fine%ke-1:2, lh) = coarse%bnd_x(coarse%js:coarse%je,   coarse%ks:coarse%ke,   lh)
            fine%bnd_x(fine%js+1:fine%je  :2, fine%ks  :fine%ke-1:2, lh) = fine  %bnd_x(fine%js  :fine%je-1:2, fine%ks  :fine%ke-1:2, lh)
            fine%bnd_x(fine%js  :fine%je,     fine%ks+1:fine%ke  :2, lh) = fine  %bnd_x(fine%js  :fine%je,     fine%ks  :fine%ke-1:2, lh)
         endif
         if (is_external(ydim, lh)) then
            fine%bnd_y(fine%is  :fine%ie-1:2, fine%ks  :fine%ke-1:2, lh) = coarse%bnd_y(coarse%is:coarse%ie,   coarse%ks:coarse%ke,   lh)
            fine%bnd_y(fine%is+1:fine%ie  :2, fine%ks  :fine%ke-1:2, lh) = fine  %bnd_y(fine%is  :fine%ie-1:2, fine%ks  :fine%ke-1:2, lh)
            fine%bnd_y(fine%is  :fine%ie,     fine%ks+1:fine%ke  :2, lh) = fine  %bnd_y(fine%is  :fine%ie,     fine%ks  :fine%ke-1:2, lh)
         endif
         if (is_external(zdim, lh)) then
            fine%bnd_z(fine%is  :fine%ie-1:2, fine%js  :fine%je-1:2, lh) = coarse%bnd_z(coarse%is:coarse%ie,   coarse%js:coarse%je,   lh)
            fine%bnd_z(fine%is+1:fine%ie  :2, fine%js  :fine%je-1:2, lh) = fine  %bnd_z(fine%is  :fine%ie-1:2, fine%js  :fine%je-1:2, lh)
            fine%bnd_z(fine%is  :fine%ie,     fine%js+1:fine%je  :2, lh) = fine  %bnd_z(fine%is  :fine%ie,     fine%js  :fine%je-1:2, lh)
         endif
      enddo

   end subroutine prolong_ext_bnd0

!!$ ============================================================================
!>
!! \brief Prolong boundaries by linear or quadratic interpolation.
!!
!! \details some code is replicated from prolong_faces
!! \todo write something more general, a routine that takes arrays or pointers and does the 2D prolongation
!<

   subroutine prolong_ext_bnd2(coarse)

      use multigridvars,   only: plvl, is_external
      use constants,       only: HI, LO, xdim, ydim, zdim
      use dataio_pub,      only: die

      implicit none

      type(plvl), pointer, intent(in) :: coarse !< level to prolong from

      type(plvl), pointer :: fine

      integer                       :: i, j, k, lh
      real, parameter, dimension(3) :: p0  = [ 0.,       1.,     0.     ] ! injection
      real, parameter, dimension(3) :: p1  = [ 0.,       3./4.,  1./4.  ] ! 1D linear prolongation stencil
      real, parameter, dimension(3) :: p2i = [ -1./8.,   1.,     1./8.  ] ! 1D integral cubic prolongation stencil
      real, parameter, dimension(3) :: p2d = [ -3./32., 30./32., 5./32. ] ! 1D direct cubic prolongation stencil
      real, dimension(-1:1)         :: p
      real, dimension(-1:1,-1:1,2,2):: pp   ! 2D prolongation stencil

      if (.not. associated(coarse)) call die("[multigridmultipole:prolong_ext_bnd0] coarse == null()")
      fine   => coarse%finer
      if (.not. associated(fine)) call die("[multigridmultipole:prolong_ext_bnd0] fine == null()")

      select case (ord_prolong_mpole)
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
         pp(i,:,1,1) = p( i)*p(:)
         pp(i,:,1,2) = p( i)*p(1:-1:-1)
         pp(i,:,2,1) = p(-i)*p(:)
         pp(i,:,2,2) = p(-i)*p(1:-1:-1)
      enddo

      !at edges and corners we can only inject
      call prolong_ext_bnd0(coarse) ! overkill: replace this

      do lh = LO, HI
         if (is_external(xdim, lh)) then
            do j = coarse%js+1, coarse%je-1
               do k = coarse%ks+1, coarse%ke-1
                  fine%bnd_x(-fine%js+2*j,  -fine%ks+2*k,  lh) =sum(pp(:,:,1,1) * coarse%bnd_x(j-1:j+1,k-1:k+1,lh))
                  fine%bnd_x(-fine%js+2*j+1,-fine%ks+2*k,  lh) =sum(pp(:,:,2,1) * coarse%bnd_x(j-1:j+1,k-1:k+1,lh))
                  fine%bnd_x(-fine%js+2*j,  -fine%ks+2*k+1,lh) =sum(pp(:,:,1,2) * coarse%bnd_x(j-1:j+1,k-1:k+1,lh))
                  fine%bnd_x(-fine%js+2*j+1,-fine%ks+2*k+1,lh) =sum(pp(:,:,2,2) * coarse%bnd_x(j-1:j+1,k-1:k+1,lh))
               enddo
            enddo
         endif

         if (is_external(ydim, lh)) then
            do i = coarse%is+1, coarse%ie-1
               do k = coarse%ks+1, coarse%ke-1
                  fine%bnd_y(-fine%is+2*i,  -fine%ks+2*k,  lh) =sum(pp(:,:,1,1) * coarse%bnd_y(i-1:i+1,k-1:k+1,lh))
                  fine%bnd_y(-fine%is+2*i+1,-fine%ks+2*k,  lh) =sum(pp(:,:,2,1) * coarse%bnd_y(i-1:i+1,k-1:k+1,lh))
                  fine%bnd_y(-fine%is+2*i,  -fine%ks+2*k+1,lh) =sum(pp(:,:,1,2) * coarse%bnd_y(i-1:i+1,k-1:k+1,lh))
                  fine%bnd_y(-fine%is+2*i+1,-fine%ks+2*k+1,lh) =sum(pp(:,:,2,2) * coarse%bnd_y(i-1:i+1,k-1:k+1,lh))
               enddo
            enddo
         endif

         if (is_external(zdim, lh)) then
            do i = coarse%is+1, coarse%ie-1
               do j = coarse%js+1, coarse%je-1
                  fine%bnd_z(-fine%is+2*i,  -fine%js+2*j,  lh) =sum(pp(:,:,1,1) * coarse%bnd_z(i-1:i+1,j-1:j+1,lh))
                  fine%bnd_z(-fine%is+2*i+1,-fine%js+2*j,  lh) =sum(pp(:,:,2,1) * coarse%bnd_z(i-1:i+1,j-1:j+1,lh))
                  fine%bnd_z(-fine%is+2*i,  -fine%js+2*j+1,lh) =sum(pp(:,:,1,2) * coarse%bnd_z(i-1:i+1,j-1:j+1,lh))
                  fine%bnd_z(-fine%is+2*i+1,-fine%js+2*j+1,lh) =sum(pp(:,:,2,2) * coarse%bnd_z(i-1:i+1,j-1:j+1,lh))
               enddo
            enddo
         endif
      enddo

   end subroutine prolong_ext_bnd2

!!$ ============================================================================
!>
!! \brief Compute multipole moments for image mass
!! \todo distribute excess of work more evenly (important only for large number of PEs, ticket:43)
!<

   subroutine img_mass2moments

      use constants,     only: xdim, ydim, zdim, GEO_XYZ, GEO_RPZ, LO, HI
      use dataio_pub,    only: die
      use domain,        only: geometry_type
      use multigridvars, only: is_external
      use mpi,           only: MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_SUM, MPI_MIN, MPI_MAX, MPI_IN_PLACE
      use mpisetup,      only: comm, ierr

      implicit none

      integer :: i, j, k, r, rr
      real, dimension(LO:HI) :: geofac

      if (geometry_type /= GEO_XYZ .and. any(CoM(xdim:zdim) /= 0.)) call die("[multipole:img_mass2moments] CoM not allowed for non-cartesian geometry")

      ! reset the multipole data
      Q(:, :, :) = 0.

      irmax = 0
      irmin = rqbin
      geofac(:) = 1.

      !OPT: try to exchange loops i < j < k -> k < j < i
      ! scan
      if (any(is_external(xdim, :))) then
         if (geometry_type == GEO_RPZ) geofac(:) = [ lmpole%fbnd(xdim, LO), lmpole%fbnd(xdim, HI) ]
         do j = lmpole%js, lmpole%je
            do k = lmpole%ks, lmpole%ke
               if (is_external(xdim, LO) .and. (geometry_type /= GEO_RPZ .or. .not. zaxis_inside)) &
                    &                     call point2moments(lmpole%bnd_x(j, k, LO)*lmpole%dyz*geofac(LO), lmpole%fbnd(xdim, LO)-CoM(xdim), lmpole%y(j)-CoM(ydim), lmpole%z(k)-CoM(zdim))
               if (is_external(xdim, HI)) call point2moments(lmpole%bnd_x(j, k, HI)*lmpole%dyz*geofac(HI), lmpole%fbnd(xdim, HI)-CoM(xdim), lmpole%y(j)-CoM(ydim), lmpole%z(k)-CoM(zdim))
            enddo
         enddo
      endif

      if (any(is_external(ydim, :))) then
         do i = lmpole%is, lmpole%ie
            do k = lmpole%ks, lmpole%ke
               if (is_external(ydim, LO)) call point2moments(lmpole%bnd_y(i, k, LO)*lmpole%dxz, lmpole%x(i)-CoM(xdim), lmpole%fbnd(ydim, LO)-CoM(ydim), lmpole%z(k)-CoM(zdim))
               if (is_external(ydim, HI)) call point2moments(lmpole%bnd_y(i, k, HI)*lmpole%dxz, lmpole%x(i)-CoM(xdim), lmpole%fbnd(ydim, HI)-CoM(ydim), lmpole%z(k)-CoM(zdim))
            enddo
         enddo
      endif

      if (any(is_external(zdim, :))) then
         do i = lmpole%is, lmpole%ie
            if (geometry_type == GEO_RPZ) geofac(LO) = lmpole%x(i)
            do j = lmpole%js, lmpole%je
               if (is_external(zdim, LO)) call point2moments(lmpole%bnd_z(i, j, LO)*lmpole%dxy*geofac(LO), lmpole%x(i)-CoM(xdim), lmpole%y(j)-CoM(ydim), lmpole%fbnd(zdim, LO)-CoM(zdim))
               if (is_external(zdim, HI)) call point2moments(lmpole%bnd_z(i, j, HI)*lmpole%dxy*geofac(LO), lmpole%x(i)-CoM(xdim), lmpole%y(j)-CoM(ydim), lmpole%fbnd(zdim, HI)-CoM(zdim))
            enddo
         enddo
      endif

      call MPI_Allreduce(MPI_IN_PLACE, irmin, 1, MPI_INTEGER, MPI_MIN, comm, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, irmax, 1, MPI_INTEGER, MPI_MAX, comm, ierr)

      ! integrate radially and apply normalization factor (the (4 \pi)/(2 l  + 1) terms cancel out)
      rr = max(1, irmin)
      Q(:, INSIDE, rr-1) = Q(:, INSIDE, rr-1) * ofact(:)
      do r = rr, irmax
         Q(:, INSIDE, r) = Q(:, INSIDE, r) * ofact(:) + Q(:, INSIDE, r-1)
      enddo

      Q(:, OUTSIDE, irmax+1) = Q(:, OUTSIDE, irmax+1) * ofact(:)
      do r = irmax, rr, -1
         Q(:, OUTSIDE, r) = Q(:, OUTSIDE, r) * ofact(:) + Q(:, OUTSIDE, r+1)
      enddo

      call MPI_Allreduce(MPI_IN_PLACE, Q(:, :, irmin:irmax), size(Q(:, :, irmin:irmax)), MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)

   end subroutine img_mass2moments

!!$ ============================================================================
!>
!! \brief Compute multipole moments for a single point
!!
!! \todo try to improve accuracy with linear interpolation over radius
!<

   subroutine point2moments(mass, x, y, z)

      implicit none

      real, intent(in) :: mass    !< mass of the contributing point
      real, intent(in) :: x       !< x coordinate of the contributing point
      real, intent(in) :: y       !< y coordinate of the contributing point
      real, intent(in) :: z       !< z coordinate of the contributing point

      real    :: sin_th, cos_th, del
      real    :: Ql, Ql1, Ql2
      integer :: l, m, ir, m2s, m2c

      call geomfac4moments(mass, x, y, z, sin_th, cos_th, ir, del)

      if (.not. interp_pt2mom) del = 0.

      ! monopole, the (0,0) moment; P_0 = 1.
      Q(0, INSIDE,  ir)   = Q(0, INSIDE,  ir)   +  rn(0) * (1.-del)
      Q(0, OUTSIDE, ir)   = Q(0, OUTSIDE, ir)   + irn(0) * (1.-del)
      if (del /= 0.) then
         Q(0, INSIDE,  ir+1) = Q(0, INSIDE,  ir+1) +  rn(0) * del
         Q(0, OUTSIDE, ir+1) = Q(0, OUTSIDE, ir+1) + irn(0) * del
      endif

      ! axisymmetric (l,0) moments
      ! Legendre polynomial recurrence: l P_l = x (2l-1) P_{l-1} - (l-1) P_{l-2}, x \eqiv \cos(\theta)
      Ql2 = 0.
      Ql1 = 1.
      do l = 1, lmax
         Ql = cos_th * k12(1, l, 0) * Ql1 - k12(2, l, 0) * Ql2
         Q(l, INSIDE,  ir)   = Q(l, INSIDE,  ir)   +  rn(l) * Ql * (1.-del)
         Q(l, OUTSIDE, ir)   = Q(l, OUTSIDE, ir)   + irn(l) * Ql * (1.-del)
         if (del /= 0.) then
            Q(l, INSIDE,  ir+1) = Q(l, INSIDE,  ir+1) +  rn(l) * Ql * del
            Q(l, OUTSIDE, ir+1) = Q(l, OUTSIDE, ir+1) + irn(l) * Ql * del
         endif
         Ql2 = Ql1
         Ql1 = Ql
      enddo

      ! non-axisymmetric (l,m) moments for 1 <= m <= mmax, m <= l <= lmax.
      do m = 1, mmax

         m2s = lm(0, 2*m-1)
         m2c = lm(0, 2*m)
         ! The (m,m) moment
         ! Associated Legendre polynomial: P_m^m = (-1)^m (2m-1)!! (1-x^2)^{m/2}
         ! The (2m-1)!! factor is integrated in ofact(:) array, where it mostly cancels out, note that (2m-1)!! \simeq m! exp(m/sqrt(2)) so it grows pretty fast with m
         Ql1 = sin_th ** m
         Q(m2s+m, INSIDE,  ir)   = Q(m2s+m, INSIDE,  ir)   +  rn(m) * Ql1 * sfac(m) * (1.-del)
         Q(m2c+m, INSIDE,  ir)   = Q(m2c+m, INSIDE,  ir)   +  rn(m) * Ql1 * cfac(m) * (1.-del)
         Q(m2s+m, OUTSIDE, ir)   = Q(m2s+m, OUTSIDE, ir)   + irn(m) * Ql1 * sfac(m) * (1.-del)
         Q(m2c+m, OUTSIDE, ir)   = Q(m2c+m, OUTSIDE, ir)   + irn(m) * Ql1 * cfac(m) * (1.-del)
         if (del /= 0.) then
            Q(m2s+m, INSIDE,  ir+1) = Q(m2s+m, INSIDE,  ir+1) +  rn(m) * Ql1 * sfac(m) * del
            Q(m2c+m, INSIDE,  ir+1) = Q(m2c+m, INSIDE,  ir+1) +  rn(m) * Ql1 * cfac(m) * del
            Q(m2s+m, OUTSIDE, ir+1) = Q(m2s+m, OUTSIDE, ir+1) + irn(m) * Ql1 * sfac(m) * del
            Q(m2c+m, OUTSIDE, ir+1) = Q(m2c+m, OUTSIDE, ir+1) + irn(m) * Ql1 * cfac(m) * del
         endif

         !>
         !! \deprecated BEWARE: most of computational cost of multipoles is here
         !! from (m+1,m) to (lmax,m)
         !! Associated Legendre polynomial: P_{m+1}^m = x (2m+1) P_m^m
         !! Associated Legendre polynomial recurrence: (l-m) P_l^m = x (2l-1) P_{l-1}^m - (l+m-1) P_{l-2}^m
         !<
         Ql2 = 0.
         do l = m + 1, lmax
            Ql = cos_th * k12(1, l, m) * Ql1 - k12(2, l, m) * Ql2
            Q(m2s+l, INSIDE,  ir)   = Q(m2s+l, INSIDE,  ir)   +  rn(l) * Ql * sfac(m) * (1.-del)
            Q(m2c+l, INSIDE,  ir)   = Q(m2c+l, INSIDE,  ir)   +  rn(l) * Ql * cfac(m) * (1.-del)
            Q(m2s+l, OUTSIDE, ir)   = Q(m2s+l, OUTSIDE, ir)   + irn(l) * Ql * sfac(m) * (1.-del)
            Q(m2c+l, OUTSIDE, ir)   = Q(m2c+l, OUTSIDE, ir)   + irn(l) * Ql * cfac(m) * (1.-del)
            if (del /= 0.) then
               Q(m2s+l, INSIDE,  ir+1) = Q(m2s+l, INSIDE,  ir+1) +  rn(l) * Ql * sfac(m) * del
               Q(m2c+l, INSIDE,  ir+1) = Q(m2c+l, INSIDE,  ir+1) +  rn(l) * Ql * cfac(m) * del
               Q(m2s+l, OUTSIDE, ir+1) = Q(m2s+l, OUTSIDE, ir+1) + irn(l) * Ql * sfac(m) * del
               Q(m2c+l, OUTSIDE, ir+1) = Q(m2c+l, OUTSIDE, ir+1) + irn(l) * Ql * cfac(m) * del
            endif
            Ql2 = Ql1
            Ql1 = Ql
         enddo

      enddo

   end subroutine point2moments

!!$ ============================================================================
!>
!! \brief Compute infinite-boundary potential from multipole moments
!!
!! \todo distribute excess of work more evenly (important only for large number of PEs, ticket:43)
!<

   subroutine moments2bnd_potential

      use dataio_pub,      only: die
      use constants,       only: xdim, ydim, zdim, GEO_XYZ, GEO_RPZ, LO, HI
      use multigridvars,   only: is_external
      use domain,          only: geometry_type

      implicit none

      integer :: i, j, k

      if (geometry_type /= GEO_XYZ .and. any(CoM(xdim:zdim) /= 0.)) call die("[multipole:img_mass2moments] CoM not allowed for non-cartesian geometry")

      if (any(is_external(xdim, :))) then
         do j = lmpole%js, lmpole%je
            do k = lmpole%ks, lmpole%ke
               if (is_external(xdim, LO) .and. (geometry_type /= GEO_RPZ .or. .not. zaxis_inside)) &
                    &                     call moments2pot(lmpole%bnd_x(j, k, LO), lmpole%fbnd(xdim, LO)-CoM(xdim), lmpole%y(j)-CoM(ydim), lmpole%z(k)-CoM(zdim))
               if (is_external(xdim, HI)) call moments2pot(lmpole%bnd_x(j, k, HI), lmpole%fbnd(xdim, HI)-CoM(xdim), lmpole%y(j)-CoM(ydim), lmpole%z(k)-CoM(zdim))
            enddo
         enddo
      endif

      if (any(is_external(ydim, :))) then
         do i = lmpole%is, lmpole%ie
            do k = lmpole%ks, lmpole%ke
               if (is_external(ydim, LO)) call moments2pot(lmpole%bnd_y(i, k, LO), lmpole%x(i)-CoM(xdim), lmpole%fbnd(ydim, LO)-CoM(ydim), lmpole%z(k)-CoM(zdim))
               if (is_external(ydim, HI)) call moments2pot(lmpole%bnd_y(i, k, HI), lmpole%x(i)-CoM(xdim), lmpole%fbnd(ydim, HI)-CoM(ydim), lmpole%z(k)-CoM(zdim))
            enddo
         enddo
      endif

      if (any(is_external(zdim, :))) then
         do i = lmpole%is, lmpole%ie
            do j = lmpole%js, lmpole%je
               if (is_external(zdim, LO)) call moments2pot(lmpole%bnd_z(i, j, LO), lmpole%x(i)-CoM(xdim), lmpole%y(j)-CoM(ydim), lmpole%fbnd(zdim, LO)-CoM(zdim))
               if (is_external(zdim, HI)) call moments2pot(lmpole%bnd_z(i, j, HI), lmpole%x(i)-CoM(xdim), lmpole%y(j)-CoM(ydim), lmpole%fbnd(zdim, HI)-CoM(zdim))
            enddo
         enddo
      endif

   end subroutine moments2bnd_potential

!!$ ============================================================================
!>
!! \brief Compute potential from multipole moments at a single point
!!
!! \todo improve accuracy with linear interpolation over radius
!<

   subroutine moments2pot(potential, x, y, z)

      use units, only: newtong

      implicit none

      real, intent(in)  :: x         !< x coordinate of the contributing point
      real, intent(in)  :: y         !< y coordinate of the contributing point
      real, intent(in)  :: z         !< z coordinate of the contributing point
      real, intent(out) :: potential !< calculated potential at given point

      real :: sin_th, cos_th, del
      real :: Ql, Ql1, Ql2
      integer :: l, m, ir, m2s, m2c

      call geomfac4moments(-newtong, x, y, z, sin_th, cos_th, ir, del)

      if (.not. interp_mom2pot) del = 0.

      ! monopole, the (0,0) moment; P_0 = 1.
      potential = (1.-del) * ( &
           Q(0, INSIDE,  ir)   * irn(0) + &
           Q(0, OUTSIDE, ir+1) *  rn(0) )
      if (del /= 0.) potential = potential + del * ( &
           Q(0, INSIDE,  ir-1) * irn(0) + &
           Q(0, OUTSIDE, ir)   *  rn(0) )
      ! ir+1 to prevent duplicate accounting contributions from ir bin; alternatively one can modify radial integration

      ! axisymmetric (l,0) moments
      Ql2 = 0.
      Ql1 = 1.
      do l = 1, lmax
         Ql = cos_th * k12(1, l, 0) * Ql1 - k12(2, l, 0) * Ql2
         potential = potential + Ql * (1.-del) * ( &
              &      Q(l, INSIDE,  ir)   * irn(l) + &
              &      Q(l, OUTSIDE, ir+1) *  rn(l) )
         if (del /= 0.) potential = potential + Ql  * del * ( &
              &      Q(l, INSIDE,  ir-1) * irn(l) + &
              &      Q(l, OUTSIDE, ir)   *  rn(l) )
         Ql2 = Ql1
         Ql1 = Ql
      enddo

      ! non-axisymmetric (l,m) moments for 1 <= m <= mmax, m <= l <= lmax.
      do m = 1, mmax
         m2s = lm(0, 2*m-1)
         m2c = lm(0, 2*m)
         ! The (m,m) moment
         Ql1 = sin_th ** m
         potential = potential + Ql1 * (1.-del) * ( &
              &      (Q(m2c+m, INSIDE,  ir)   * irn(m) + &
              &       Q(m2c+m, OUTSIDE, ir+1) *  rn(m) ) * cfac(m) + &
              &      (Q(m2s+m, INSIDE,  ir)   * irn(m) + &
              &       Q(m2s+m, OUTSIDE, ir+1) *  rn(m) ) * sfac(m) )
         if (del /= 0.) potential = potential + Ql1 * del * ( &
              &      (Q(m2c+m, INSIDE,  ir-1) * irn(m) + &
              &       Q(m2c+m, OUTSIDE, ir)   *  rn(m) ) * cfac(m) + &
              &      (Q(m2s+m, INSIDE,  ir-1) * irn(m) + &
              &       Q(m2s+m, OUTSIDE, ir)   *  rn(m) ) * sfac(m) )

         !> \deprecated BEWARE: lots of computational cost of multipoles is here
         ! from (m+1,m) to (lmax,m)
         Ql2 = 0.
         do l = m+1, lmax
            Ql = cos_th * k12(1, l, m) * Ql1 - k12(2, l, m) * Ql2
            potential = potential + Ql * (1.-del) * ( &
                 &      (Q(m2c+l, INSIDE,  ir)   * irn(l) + &
                 &       Q(m2c+l, OUTSIDE, ir+1) *  rn(l) ) * cfac(m) + &
                 &      (Q(m2s+l, INSIDE,  ir)   * irn(l) + &
                 &       Q(m2s+l, OUTSIDE, ir+1) *  rn(l) ) * sfac(m) )
            if (del /= 0.) potential = potential + Ql * del * ( &
                 &      (Q(m2c+l, INSIDE,  ir-1) * irn(l) + &
                 &       Q(m2c+l, OUTSIDE, ir)   *  rn(l) ) * cfac(m) + &
                 &      (Q(m2s+l, INSIDE,  ir-1) * irn(l) + &
                 &       Q(m2s+l, OUTSIDE, ir)   *  rn(l) ) * sfac(m) )
            Ql2 = Ql1
            Ql1 = Ql
         enddo
      enddo

   end subroutine moments2pot

!!$ ============================================================================
!>
!! \brief This routine calculates various geometrical numbers required for multipole evaluation
!! \details It modifies the rn(:), irn(:), cfac(:) and sfac(:) arrays. Scalars are passed through argument list.
!<

   subroutine geomfac4moments(factor, x, y, z, sin_th, cos_th, ir, delta)

      use constants,  only: GEO_XYZ, GEO_RPZ
      use dataio_pub, only: die, msg
      use domain,     only: geometry_type

      implicit none

      real,    intent(in)  :: factor         !< scaling factor (e.g. mass) of the contributing point
      real,    intent(in)  :: x              !< x coordinate of the contributing point
      real,    intent(in)  :: y              !< x coordinate of the contributing point
      real,    intent(in)  :: z              !< x coordinate of the contributing point
      real,    intent(out) :: sin_th         !< sine of the theta angle
      real,    intent(out) :: cos_th         !< cosine of the theta angle
      integer, intent(out) :: ir             !< radial index for the Q(:, :, r) array
      real,    intent(out) :: delta          !< fraction of the radial cell for interpolation between ir and ir+1

      real    :: rxy, r, rinv
      real    :: sin_ph, cos_ph
      integer :: l, m

      ! radius and its projection onto XY plane
      select case (geometry_type)
         case (GEO_XYZ)
            rxy = x**2 + y**2
         case (GEO_RPZ)
            rxy = x**2
         case default
            call die("[multigridmultipole:geomfac4moments] Unsupported geometry.")
            rxy = 0.
      end select
      r    = sqrt(rxy + z**2)
      rxy  = sqrt(rxy)
      if (r /= 0.) then
         rinv = 1. / r
      else
         rinv = 0.
      endif

      !radial index for the Q(:, :, r) array
      ir = int(r / drq)
      if (ir < rqbin) then
         delta = r/drq - ir
      else
         delta = 0
      endif
      if (ir > rqbin .or. ir < 0) then
         write(msg,'(2(a,i7),a)')"[multipole:geomfac4moments] radial index = ",ir," outside Q(:, :, ",rqbin,") range"
         call die(msg)
      endif
      irmax = max(irmax, ir)
      irmin = min(irmin, ir)

      ! azimuthal angle sine and cosine tables
      ! ph = atan2(y, x); cfac(m) = cos(m * ph); sfac(m) = sin(m * ph)
      ! cfac(0) and sfac(0) are set in init_multigrid
      if (rxy /= 0.) then
         select case (geometry_type)
            case (GEO_XYZ)
               cos_ph = x / rxy
               sin_ph = y / rxy
            case (GEO_RPZ)
               cos_ph = cos(y)
               sin_ph = sin(y)
            case default
               call die("[multigridmultipole:geomfac4moments] Unsupported geometry.")
               cos_ph = 0. ; sin_ph = 0.
         end select
      else
         cos_ph = 1.
         sin_ph = 0.
      endif
!> \todo Possible optimization: number of computed elements can be doubled on each loop iteration (should give better pipelining)
      do m = 1, mmax
         cfac(m) = cos_ph*cfac(m-1) - sin_ph*sfac(m-1)
         sfac(m) = cos_ph*sfac(m-1) + sin_ph*cfac(m-1)
      enddo

      ! vertical angle
      cos_th = z   * rinv
      sin_th = rxy * rinv

!> \todo check how much it would degrade solution and improve performance to move the multiplications by rn(:) and irn(:) to img_mass2moments (before and after integration)
      ! rn(l) = factor * r ** l; irn(l) = factor * r ** -(l+1)
      rn(0)  = factor
      irn(0) = factor * rinv
      ! scale r to reduce risk of occurring a Floating Point Overflow for high l_max and large r values. Leave the original value of r only in irn(0)
      r = r / rscale
      rinv = rinv * rscale
!> \todo Possible optimization: number of computed elements can be doubled on each loop iteration (should give better pipelining)
      do l = 1, lmax
         rn(l)  =  rn(l-1) * r
         irn(l) = irn(l-1) * rinv
      enddo

   end subroutine geomfac4moments

!!$ ============================================================================
!>
!! \brief HEAVY_DEBUG marks routines that normally are never called, but at some point
!! were useful to test correctness or something.
!<

!#define HEAVY_DEBUG
#ifdef HEAVY_DEBUG

!!$ ============================================================================
!>
!! \brief Quick test for correctness of the multipole solver.
!! \details cphi should agree well with phi with largest errors at r = sqrt(sum(p(1:3)**2))
!<

#error "The test_multipoles routine is outdated and was commented out"

!!$   subroutine test_multipoles
!!$
!!$      use dataio_pub,       only: die, printinfo, msg
!!$      use units,            only: newtong
!!$      use multigridhelpers, only: dirty_debug
!!$      use debug,            only: aux_R, aux_I
!!$
!!$      implicit none
!!$
!!$      integer :: i, j, r, rr
!!$      real :: phi, cphi
!!$      real, dimension(0:3) :: p, x
!!$
!!$      ! reset the multipole data
!!$      Q(:, :, :, :) = 0.
!!$
!!$      irmax = 0
!!$      irmin = rqbin
!!$
!!$      p(0) = 1./newtong
!!$      p(1:3) = aux_R(1:3)
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
!!$      enddo
!!$
!!$      Q(1:lmax, :, OUTSIDE, irmax+1) = Q(1:lmax, :, OUTSIDE, irmax+1) * ofact(1:lmax, :)
!!$      do r = irmax, rr, -1
!!$         Q(0,      :, OUTSIDE, r) = Q(0,      :, OUTSIDE, r)                    + Q(0,      :, OUTSIDE, r+1)
!!$         Q(1:lmax, :, OUTSIDE, r) = Q(1:lmax, :, OUTSIDE, r) * ofact(1:lmax, :) + Q(1:lmax, :, OUTSIDE, r+1)
!!$      enddo
!!$
!!$      do i = -60, 60
!!$         do j = -60, 60
!!$            x(1) = 0.05 * i
!!$            x(2) = 0.05 * j
!!$            x(3) = 0.1 * aux_I(1)
!!$
!!$            cphi = - p(0)*newtong / sqrt(1e-290+sum((p(1:3)-x(1:3))**2))
!!$
!!$            call moments2pot(phi, x(1), x(2), x(3))
!!$            write(msg,'(a,3f7.2,2(a,g15.5))')" xyz= ",x(1:3)," phi = ",phi," calc= ",cphi
!!$            call printinfo(msg, .false.)
!!$         enddo
!!$      enddo
!!$
!!$      call die("qniec")
!!$
!!$   end subroutine test_multipoles

#endif /* HEAVY_DEBUG */

end module multipole
