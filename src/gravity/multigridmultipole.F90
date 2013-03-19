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
   use constants,   only: ndims, xdim
#if defined(__INTEL_COMPILER)
      !! \deprecated remove this clause as soon as Intel Compiler gets required
      !! features and/or bug fixes
   use cg_list_bnd,   only: cg_list_bnd_T   ! QA_WARN intel
#endif /* __INTEL_COMPILER */
!   use cg_level_connected, only: cg_level_connected_T

   implicit none

   private
   public :: init_multipole, cleanup_multipole, multipole_solver
   public :: lmax, mmax, ord_prolong_mpole, coarsen_multipole, use_point_monopole, interp_pt2mom, interp_mom2pot ! initialized in multigrid_gravity

   integer, parameter        :: INSIDE = 1                       !< index for interior multipole expansion
   integer, parameter        :: OUTSIDE = INSIDE + 1             !< index for exterior multipole expansion

   ! namelist parameters for MULTIGRID_GRAVITY
   integer(kind=4)           :: lmax                             !< Maximum l-order of multipole moments
   integer(kind=4)           :: mmax                             !< Maximum m-order of multipole moments. Equal to lmax by default.
   integer(kind=4)           :: ord_prolong_mpole                !< boundary prolongation operator order; allowed values are -2 .. 2
   integer(kind=4)           :: coarsen_multipole                !< If > 0 then evaluate multipoles at finest%level%level_id-coarsen_multipole level
   logical                   :: use_point_monopole               !< Don't evaluate multipole moments, use point-like mass approximation (crudest possible)
   logical                   :: interp_pt2mom                    !< Distribute contribution from a cell between two adjacent radial bins (linear interpolation in radius)
   logical                   :: interp_mom2pot                   !< Compute the potential from moments from two adjacent radial bins (linear interpolation in radius)

   ! radial discretization
   integer                   :: rqbin                            !< number of radial samples of multipoles
   real                      :: drq                              !< radial resolution of multipoles
   real                      :: rscale                           !< scaling factor that limits risk of floating point overflow for high powers of radius (rn(:) and irn(:))
   integer                   :: irmin                            !< minimum Q(:, :, r) index in use
   integer                   :: irmax                            !< maximum Q(:, :, r) index in use

!> \todo OPT derive a special set of leaves and coarsened leaves that can be safely used here
!   type(cg_level_connected_T), pointer :: lmpole                 !< pointer to the level where multipoles are evaluated
   integer, parameter                  :: imass = xdim-1         !< index for mass in CoM(:)
   real, dimension(imass:ndims)        :: CoM                    !< Total mass and center of mass coordinates
   logical                             :: zaxis_inside           !< true when z-axis belongs to the inner radial boundary in polar coordinates

   ! multipoles and auxiliary factors
   real, dimension(:,:,:), allocatable :: Q                      !< The whole moment array with dependence on radius
   real, dimension(:,:,:), allocatable :: k12                    !< array of Legendre recurrence factors
   real, dimension(:),     allocatable :: sfac                   !< 0..m_max sized array of azimuthal sine factors
   real, dimension(:),     allocatable :: cfac                   !< 0..m_max sized array of azimuthal cosine factors
   real, dimension(:),     allocatable :: ofact                  !< arrays of Legendre normalization factor (compressed)
   real, dimension(:),     allocatable :: rn                     !< 0..l_max sized array of positive powers of r
   real, dimension(:),     allocatable :: irn                    !< 0..l_max sized array of negative powers of r

contains

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

!> \brief Initialization routine, called once, from init_multigrid

   subroutine init_multipole

      use constants,  only: ndims
      use dataio_pub, only: die, warn
      use domain,     only: dom
      use mpisetup,   only: master

      implicit none

      integer :: l, m

      if (dom%eff_dim /= ndims) call die("[multigridmultipole:init_multipole] Only 3D is supported") !> \todo add support for 2D RZ

      !fixup multipole moments
      if (mmax > lmax) then
         if (master) call warn("[multigridmultipole:init_multipole] mmax reduced to lmax")
         mmax = lmax
      endif
      if (mmax < 0) mmax = lmax

      if (.not. use_point_monopole) then

         if (allocated(k12) .or. allocated(ofact)) call die("[multigridmultipole:init_multipole] k12 or ofact already allocated")
         allocate(k12(2, 1:lmax, 0:mmax), ofact(0:lm(int(lmax), int(2*mmax))))

         if (allocated(sfac) .or. allocated(cfac)) call die("[multigridmultipole:init_multipole] sfac or cfac already allocated")
         allocate(sfac(0:mmax), cfac(0:mmax))

         if (allocated(rn) .or. allocated(irn)) call die("[multigridmultipole:init_multipole] rn or irn already allocated")
         allocate(rn(0:lmax), irn(0:lmax))

         ofact(:) = 0. ! prevent spurious FP exceptions in multipole:img_mass2moments
         do l = 1, lmax
            do m = 0, min(l, int(mmax))
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

!>
!! \brief Initialization routine, called once per entry to the multipole solver.
!!
!! \details This routine reinitializes everything that may need reinitialization due to AMR activity
!!
!! \todo OPT: Detect only changes in highest required level
!<

   subroutine refresh_multipole

      use cg_level_finest, only: finest
      use constants,       only: small, pi, xdim, ydim, zdim, GEO_XYZ, GEO_RPZ, LO, HI, pMIN
      use dataio_pub,      only: die, warn
      use domain,          only: dom
      use mpisetup,        only: master, piernik_MPI_Allreduce
      use particle_pub,    only: pset

      implicit none

      integer :: i !, l

      ! Multipole coarsening may significantly improve performance at a cost of accuracy
      ! It requires a lot of new code to work fully on irregularly refined grids made by AMR
!!$      lmpole => finest
!!$      !> \todo find finest grid at the external outer boundary (i.e. non-reflecting, not periodic, ...), which might be significantly coarser than the globally finest grid.
!!$      do l = 1, coarsen_multipole
!!$         call die("[multigridmultipole:refresh_multipole] multipole coarsening is temporrarily disabled")
!!$         if (associated(lmpole%coarser)) then
!!$            lmpole => lmpole%coarser
!!$         else
!!$            if (master) call warn("[multigridmultipole:refresh_multipole] too deep multipole coarsening.")
!!$         endif
!!$      enddo

      if (coarsen_multipole > 0) then
         if (interp_pt2mom) then
            call warn("[multigridmultipole:refresh_multipole] coarsen_multipole > 0 disables interp_pt2mom.")
            interp_pt2mom = .false.
         endif
         if (interp_mom2pot) then
            call warn("[multigridmultipole:refresh_multipole] coarsen_multipole > 0 disables interp_mom2pot.")
            interp_mom2pot = .false.
         endif
      endif

      ! assume that Center of Mass is approximately in the center of computational domain by default
      CoM(imass) = 1.

      select case (dom%geometry_type)
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
            zaxis_inside = dom%edge(xdim, LO) <= dom%L_(xdim)/finest%level%n_d(xdim) ! lmpole
            if (master) then
               if (zaxis_inside) call warn("[multigridmultipole:refresh_multipole] Setups with Z-axis at the edge of the domain may not work as expected yet.")
               if (use_point_monopole) call warn("[multigridmultipole:refresh_multipole] Point-like monopole is not implemented.")
            endif
            use_point_monopole = .false.
         case default
            call die("[multigridmultipole:refresh_multipole] Unsupported geometry.")
      end select

      if (.not. use_point_monopole) then

         if (associated(finest%level%first)) then
            select case (dom%geometry_type)
               !> \warning With refinement lmpole might no longer be a single level
               case (GEO_XYZ)
                  drq = minval(finest%level%first%cg%dl(:), mask=dom%has_dir(:)) / 2.
               case (GEO_RPZ)
                  drq = min(finest%level%first%cg%dx, dom%C_(xdim)*finest%level%first%cg%dy, finest%level%first%cg%dz) / 2.
               case default
                  call die("[multigridmultipole:refresh_multipole] Unsupported geometry.")
            end select
         else
            drq = maxval(dom%L_(:))
         endif
         call piernik_MPI_Allreduce(drq, pMIN)

         select case (dom%geometry_type)
            case (GEO_XYZ)
               rqbin = int(sqrt(sum(dom%L_(:)**2))/drq) + 1
               ! arithmetic average of the closest and farthest points of computational domain with respect to its center
               !>
               !!\todo check what happens if there are points that are really close to the domain center (maybe we should use a harmonic average?)
               !! Issue a warning or error if it is known that given lmax leads to FP overflows in rn(:) and irn(:)
               !<
               rscale = ( minval(dom%L_(:)) + sqrt(sum(dom%L_(:)**2)) )/4.
            case (GEO_RPZ)
               rqbin = int(sqrt((2.*dom%edge(xdim, HI))**2 + dom%L_(zdim)**2)/drq) + 1
               rscale = ( min(2.*dom%edge(xdim, HI), dom%L_(zdim)) + sqrt((2.*dom%edge(xdim, HI))**2 + dom%L_(zdim)**2) )/4.
            case default
               call die("[multigridmultipole:refresh_multipole] Unsupported geometry.")
         end select

         do i = lbound(pset%p, dim=1), ubound(pset%p, dim=1)
            !> \warning this is not optimal, when domain is placed far away form the origin
            !> \warning a factor of up to 2 may be required if we use CoM as the origin
            if (pset%p(i)%outside) rqbin = max(rqbin, int(sqrt(sum(pset%p(i)%pos**2))/drq) + 1)
         enddo

         if (allocated(Q)) deallocate(Q)
         allocate(Q(0:lm(int(lmax), int(2*mmax)), INSIDE:OUTSIDE, 0:rqbin))

      endif

   end subroutine refresh_multipole

!> \brief Multipole cleanup

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

!>
!! \brief Multipole solver
!!
!! \todo improve multipole expansion on coarser grids
!! (see. "A Scalable Parallel Poisson Solver in Three Dimensions with Infinite-Domain Boundary Conditions" by McCorquodale, Colella, Balls and Baden).
!! Coarsening by one level would reduce the multipole costs by a factor of 4.
!<

   subroutine multipole_solver

      use cg_leaves,          only: leaves
!!$      use cg_level_connected, only: finest !, cg_level_connected_T
      use constants,          only: dirtyH
!!$      use dataio_pub,         only: die
      use global,             only: dirty_debug
!!$      use multigridvars,      only: solution
#ifdef MACLAURIN_PROBLEM
      use problem_pub,        only: maclaurin2bnd_potential
#endif /* MACLAURIN_PROBLEM */

      implicit none

!!$      type(cg_level_connected_T), pointer :: curl

      call refresh_multipole

!      if (.not. associated(lmpole, finest)) call die("[multigridmultipole:multipole_solver] lmpole /= finest requires a lot of work")

      if (dirty_debug) then
         call leaves%reset_boundaries(dirtyH)
      else
         call leaves%reset_boundaries
      endif

!!$      if (.not. associated(lmpole, finest)) then
!!$         curl => finest
!!$         do while (associated(curl) .and. .not. associated(curl, lmpole)) ! do lev = finest%level%level_id, lmpole%first%cg%level_id + 1, -1
!!$            call curl%restrict_q_1var(solution)  ! Overkill, only some layers next to external boundary are needed.
!!$            curl => curl%coarser
!!$         enddo                                ! An alternative: do potential2img_mass on the finest and restrict bnd_[xyz] data.
!!$      endif
      call potential2img_mass

#ifndef MACLAURIN_PROBLEM
      if (use_point_monopole) then
         call find_img_CoM
         call isolated_monopole
      else
         ! results seems to be slightly better without find_img_CoM.
         ! With CoM or when it is known than CoM is close to the domain center one may try to save some CPU time by lowering mmax.
         call img_mass2moments
         call moments2bnd_potential
      endif
#else /* MACLAURIN_PROBLEM */
      call maclaurin2bnd_potential
#endif /* ! MACLAURIN_PROBLEM */

      !> \todo The approach with lmpole should be reworked completely. With AMR use only current level and reinitialize some things if refinements have changed
!!$      if (.not. associated(lmpole, finest)) then
!!$         curl => lmpole
!!$         do while (associated(curl) .and. .not. associated(curl, finest)) ! do lev = lmpole%first%cg%level_id, finest%level%level_id - 1
!!$            call prolong_ext_bnd(curl)
!!$            curl => curl%finer
!!$         enddo
!!$      endif

   end subroutine multipole_solver

!>
!! \brief Set boundary potential from monopole source. Fill lmpole%first%cg%mg%bnd_[xyz] arrays with expected values of the gravitational potential at external face of computational domain.
!! \details This is a simplified approach that can be used for tests and as a fast replacement for the
!! multipole boundary solver for nearly spherically symmetric source distributions.
!! The isolated_monopole subroutine ignores the radial profile of the monopole
!<

   subroutine isolated_monopole

      use cg_list,    only: cg_list_element
      use cg_leaves,  only: leaves
      use constants,  only: xdim, ydim, zdim, LO, HI, GEO_XYZ !, GEO_RPZ
      use dataio_pub, only: die
      use domain,     only: dom
      use grid_cont,  only: grid_container
      use units,      only: newtong

      implicit none

      integer :: i, j, k, lh
      real    :: r2
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      if (dom%geometry_type /= GEO_XYZ) call die("[multigridmultipole:isolated_monopole] non-cartesian geometry not implemented yet")

      cgl => leaves%first ! lmpole
      do while (associated(cgl))
         cg => cgl%cg
         do lh = LO, HI
            if (cg%ext_bnd(xdim, lh)) then
               do j = cg%js, cg%je
                  do k = cg%ks, cg%ke
                     r2 = (cg%y(j)-CoM(ydim))**2 + (cg%z(k) - CoM(zdim))**2
                     cg%mg%bnd_x(j, k, lh) = - newtong * CoM(imass) / sqrt(r2 + (cg%fbnd(xdim, lh)-CoM(xdim))**2)
                  enddo
               enddo
            endif
            if (cg%ext_bnd(ydim, lh)) then
               do i = cg%is, cg%ie
                  do k = cg%ks, cg%ke
                     r2 = (cg%x(i)-CoM(xdim))**2 + (cg%z(k) - CoM(zdim))**2
                     cg%mg%bnd_y(i, k, lh) = - newtong * CoM(imass) / sqrt(r2 + (cg%fbnd(ydim, lh)-CoM(ydim))**2)
                  enddo
               enddo
            endif
            if (cg%ext_bnd(zdim, lh)) then
               do i = cg%is, cg%ie
                  do j = cg%js, cg%je
                     r2 = (cg%x(i)-CoM(xdim))**2 + (cg%y(j) - CoM(ydim))**2
                     cg%mg%bnd_z(i, j, lh) = - newtong * CoM(imass) / sqrt(r2 + (cg%fbnd(zdim, lh)-CoM(zdim))**2)
                  enddo
               enddo
            endif
         enddo
         cgl => cgl%nxt
      enddo

   end subroutine isolated_monopole

!>
!! \brief Find total mass and its center
!!
!! \details This routine does the summation only on external boundaries
!<

   subroutine find_img_CoM

      use cg_leaves,    only: leaves
      use cg_list,      only: cg_list_element
      use constants,    only: ndims, xdim, ydim, zdim, LO, HI, GEO_XYZ, pSUM !, GEO_RPZ
      use dataio_pub,   only: die
      use domain,       only: dom
      use grid_cont,    only: grid_container
      use mpisetup,     only: piernik_MPI_Allreduce
      use particle_pub, only: pset
#ifdef DEBUG
      use dataio_pub,   only: msg, printinfo
      use mpisetup,     only: master
      use units,        only: fpiG
#endif /* DEBUG */

      implicit none

      real, dimension(imass:ndims) :: lsum, dsum
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg
      integer :: lh, i, d

      if (dom%geometry_type /= GEO_XYZ) call die("[multigridmultipole:find_img_CoM] non-cartesian geometry not implemented yet")

      lsum(:) = 0.

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         do lh = LO, HI
            if (cg%ext_bnd(xdim, lh)) then
               d = 1
               if (lh == HI) d = cg%n_b(xdim)
               dsum(imass)     =        sum( cg%mg%bnd_x(cg%js:cg%je, cg%ks:cg%ke, lh), mask=cg%leafmap(d, :, :) )
               dsum(xdim:zdim) = [ dsum(imass) * cg%fbnd(xdim, lh), &
                    &              sum( sum( cg%mg%bnd_x(cg%js:cg%je, cg%ks:cg%ke, lh), mask=cg%leafmap(d, :, :), dim=2) * cg%y(cg%js:cg%je) ), &
                    &              sum( sum( cg%mg%bnd_x(cg%js:cg%je, cg%ks:cg%ke, lh), mask=cg%leafmap(d, :, :), dim=1) * cg%z(cg%ks:cg%ke) ) ]
               lsum(:)         = lsum(:) + dsum(:) * cg%dyz
            endif
            if (cg%ext_bnd(ydim, lh)) then
               d = 1
               if (lh == HI) d = cg%n_b(ydim)
               dsum(imass)     =        sum( cg%mg%bnd_y(cg%is:cg%ie, cg%ks:cg%ke, lh), mask=cg%leafmap(:, d, :) )
               dsum(xdim:zdim) = [ sum( sum( cg%mg%bnd_y(cg%is:cg%ie, cg%ks:cg%ke, lh), mask=cg%leafmap(:, d, :), dim=2) * cg%x(cg%is:cg%ie) ), &
                    &              dsum(imass) * cg%fbnd(ydim, lh), &
                    &              sum( sum( cg%mg%bnd_y(cg%is:cg%ie, cg%ks:cg%ke, lh), mask=cg%leafmap(:, d, :), dim=1) * cg%z(cg%ks:cg%ke) ) ]
               lsum(:)         = lsum(:) + dsum(:) * cg%dxz
            endif
            if (cg%ext_bnd(zdim, lh)) then
               d = 1
               if (lh == HI) d = cg%n_b(zdim)
               dsum(imass)     =        sum( cg%mg%bnd_z(cg%is:cg%ie, cg%js:cg%je, lh), mask=cg%leafmap(:, :, d) )
               dsum(xdim:zdim) = [ sum( sum( cg%mg%bnd_z(cg%is:cg%ie, cg%js:cg%je, lh), mask=cg%leafmap(:, :, d), dim=2) * cg%x(cg%is:cg%ie) ), &
                    &              sum( sum( cg%mg%bnd_z(cg%is:cg%ie, cg%js:cg%je, lh), mask=cg%leafmap(:, :, d), dim=1) * cg%y(cg%js:cg%je) ), &
                    &              dsum(imass) * cg%fbnd(zdim, lh) ]
               lsum(:)         = lsum(:) + dsum(:) * cg%dxy
            endif
         enddo
         cgl => cgl%nxt
      enddo

      ! Add only those particles, which are placed outside the domain. Particles inside the domain were already mapped on the grid.
      !> \warning Do we need to use the fppiG factor here?
      do i = lbound(pset%p, dim=1), ubound(pset%p, dim=1)
         if (pset%p(i)%outside) lsum(:) = lsum(:) + [ pset%p(i)%mass, pset%p(i)%mass * pset%p(i)%pos(:) ]
      enddo

      CoM(imass:ndims) = lsum(imass:ndims)
      call piernik_MPI_Allreduce(CoM(imass:ndims), pSUM)

      if (CoM(imass) /= 0.) then
         CoM(xdim:zdim) = CoM(xdim:zdim) / CoM(imass)
      else
         call die("[multigridmultipole:find_img_CoM] Total mass == 0")
      endif
#ifdef DEBUG
      if (master) then
         write(msg, '(a,g14.6,a,3g14.6,a)')"[multigridmultipole:find_img_CoM] Total mass = ", CoM(imass)/fpiG," at (",CoM(xdim:zdim),")"
         call printinfo(msg)
      endif
#endif /* DEBUG */

   end subroutine find_img_CoM

!>
!! \brief Convert potential into image mass. This way we reduce a 3D problem to a 2D one.
!!
!! \details There will be work imbalance here because different PEs may operate on different amount of external boundary data.
!!
!! The value stored in bnd_[xyz] has the meaning of mass inside the cell next to the external boundary,
!! divided by the surface of the cell projected onto the external boundary surface (surface density).
!! The surface factors are taken into account in img_mass2moments routine.
!! The exact position of the mass depends on choice of a1 and a2 parameters.
!! The surface density is estimated using gradient of potential.
!! The Taylor expansion of the potential is done wrt. boundary position, normal direction points to the outside of the boundary
!! with the assumption that the potential at the external boundary equals 0.
!<

   subroutine potential2img_mass

      use cg_leaves,     only: leaves
      use cg_list,       only: cg_list_element
      use constants,     only: GEO_RPZ, LO, HI, xdim, ydim, zdim
      use domain,        only: dom
      use grid_cont,     only: grid_container
      use multigridvars, only: solution

      implicit none

      integer :: i
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg
      real, parameter :: a1 = -2., a2 = (-2. - a1)/3. ! interpolation parameters;   <---- a1=-2 => a2=0
      ! a1 = -2. is the simplest, 1st order choice, gives best agreement of total mass and CoM location when compared to 3-D integration
      ! a1 = -1., a2 = -1./3. seems to do the best job,
      ! a1 = -3./2., a2 = -1./6. seems to be 2nd order estimator

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         associate( &
            bnd_x => cg%mg%bnd_x, &
            bnd_y => cg%mg%bnd_y, &
            bnd_z => cg%mg%bnd_z, &
            soln => cg%q(solution)%arr &
         )
         !> \deprecated BEWARE: some cylindrical factors may be helpful
         if (cg%ext_bnd(xdim, LO)) then
            if (zaxis_inside .and. dom%geometry_type == GEO_RPZ) then
               bnd_x(                    cg%js:cg%je, cg%ks:cg%ke, LO) = 0. ! treat as internal
            else
               bnd_x(                    cg%js:cg%je, cg%ks:cg%ke, LO) =   ( &
                    & a1 * soln(cg%is,   cg%js:cg%je, cg%ks:cg%ke) + &
                    & a2 * soln(cg%is+1, cg%js:cg%je, cg%ks:cg%ke) ) * cg%idx
            endif
         endif

         if (cg%ext_bnd(xdim, HI)) bnd_x(   cg%js:cg%je, cg%ks:cg%ke, HI) =   ( &
              &   a1 * soln(cg%ie,   cg%js:cg%je, cg%ks:cg%ke) + &
              &   a2 * soln(cg%ie-1, cg%js:cg%je, cg%ks:cg%ke) ) * cg%idx

         if (cg%ext_bnd(ydim, LO)) then
            bnd_y           (cg%is:cg%ie,          cg%ks:cg%ke, LO) =    ( &
                 & a1 * soln(cg%is:cg%ie, cg%js,   cg%ks:cg%ke) + &
                 & a2 * soln(cg%is:cg%ie, cg%js+1, cg%ks:cg%ke) ) * cg%idy
            if (dom%geometry_type == GEO_RPZ) then
               do i = cg%is, cg%ie ! cg%inv_x(i) is sanitized for x(i) == 0.
                  bnd_y(i, cg%ks:cg%ke, LO) = bnd_y(i, cg%ks:cg%ke, LO) * cg%inv_x(i)
               enddo
            endif
         endif

         if (cg%ext_bnd(ydim, HI)) then
            bnd_y           (cg%is:cg%ie,          cg%ks:cg%ke, HI) =   ( &
                 & a1 * soln(cg%is:cg%ie, cg%je,   cg%ks:cg%ke) + &
                 & a2 * soln(cg%is:cg%ie, cg%je-1, cg%ks:cg%ke) ) * cg%idy
            if (dom%geometry_type == GEO_RPZ) then
               do i = cg%is, cg%ie
                  bnd_y(i, cg%ks:cg%ke, HI) = bnd_y(i, cg%ks:cg%ke, HI) * cg%inv_x(i)
               enddo
            endif
         endif

         if (cg%ext_bnd(zdim, LO)) bnd_z(cg%is:cg%ie, cg%js:cg%je, LO) =    ( &
              &         a1 * soln(cg%is:cg%ie, cg%js:cg%je, cg%ks) + &
              &         a2 * soln(cg%is:cg%ie, cg%js:cg%je, cg%ks+1) ) * cg%idz

         if (cg%ext_bnd(zdim, HI)) bnd_z(cg%is:cg%ie, cg%js:cg%je, HI) =   ( &
              &         a1 * soln(cg%is:cg%ie, cg%js:cg%je, cg%ke) + &
              &         a2 * soln(cg%is:cg%ie, cg%js:cg%je, cg%ke-1) ) * cg%idz
         cgl => cgl%nxt
         end associate
      enddo

   end subroutine potential2img_mass

!> \brief Prolong boundaries wrapper

   subroutine prolong_ext_bnd(coarse)

      use cg_level_connected, only: cg_level_connected_T
      use constants,     only: ndims, O_INJ, O_D2, O_I2
      use dataio_pub,    only: die
      use domain,        only: dom

      implicit none

      type(cg_level_connected_T), pointer, intent(in) :: coarse !< level to prolong from


      call die("[multigridmultipole:prolong_ext_bnd] prolong_ext_bnd is severely outdated and thus unusable at the moment")

      if (dom%eff_dim<ndims) call die("[multigridmultipole:prolong_ext_bnd] 1D and 2D not finished")
      if (abs(ord_prolong_mpole) > maxval(abs([O_D2, O_I2]))) call die("[multigridmultipole:prolong_ext_bnd] interpolation order too high")

      !> \deprecated BEWARE: do we need cylindrical factors for prolongation?
      if (ord_prolong_mpole == O_INJ) then
         call prolong_ext_bnd0(coarse)
      else
         call prolong_ext_bnd2(coarse)
      endif

   end subroutine prolong_ext_bnd

!> \brief Prolong boundaries by injection.

   subroutine prolong_ext_bnd0(coarse)

      use cg_level_connected, only: cg_level_connected_T
      use constants,     only: HI, LO, xdim, ydim, zdim
      use dataio_pub,    only: die
      use domain,        only: is_multicg

      implicit none

      type(cg_level_connected_T), pointer, intent(in) :: coarse !< level to prolong from

      type(cg_level_connected_T), pointer :: fine
      integer :: lh

      if (is_multicg) call die("[multigridmultipole:prolong_ext_bnd0] multicg not implemented yet") ! fine%first%cg%mg%bnd_[xyz]

      if (.not. associated(coarse)) call die("[multigridmultipole:prolong_ext_bnd0] coarse == null()")
      fine   => coarse%finer
      if (.not. associated(fine)) call die("[multigridmultipole:prolong_ext_bnd0] fine == null()")

      do lh = LO, HI
         if (fine%first%cg%ext_bnd(xdim, lh)) then
            fine%first%cg%mg%bnd_x(fine%first%cg%js  :fine%first%cg%je-1:2, fine%first%cg%ks  :fine%first%cg%ke-1:2, lh) = coarse%first%cg%mg%bnd_x(coarse%first%cg%js:coarse%first%cg%je,   coarse%first%cg%ks:coarse%first%cg%ke,   lh)
            fine%first%cg%mg%bnd_x(fine%first%cg%js+1:fine%first%cg%je  :2, fine%first%cg%ks  :fine%first%cg%ke-1:2, lh) = fine%first%cg%mg%bnd_x(  fine%first%cg%js  :fine%first%cg%je-1:2, fine%first%cg%ks  :fine%first%cg%ke-1:2, lh)
            fine%first%cg%mg%bnd_x(fine%first%cg%js  :fine%first%cg%je,     fine%first%cg%ks+1:fine%first%cg%ke  :2, lh) = fine%first%cg%mg%bnd_x(  fine%first%cg%js  :fine%first%cg%je,     fine%first%cg%ks  :fine%first%cg%ke-1:2, lh)
         endif
         if (fine%first%cg%ext_bnd(ydim, lh)) then
            fine%first%cg%mg%bnd_y(fine%first%cg%is  :fine%first%cg%ie-1:2, fine%first%cg%ks  :fine%first%cg%ke-1:2, lh) = coarse%first%cg%mg%bnd_y(coarse%first%cg%is:coarse%first%cg%ie,   coarse%first%cg%ks:coarse%first%cg%ke,   lh)
            fine%first%cg%mg%bnd_y(fine%first%cg%is+1:fine%first%cg%ie  :2, fine%first%cg%ks  :fine%first%cg%ke-1:2, lh) = fine%first%cg%mg%bnd_y(  fine%first%cg%is  :fine%first%cg%ie-1:2, fine%first%cg%ks  :fine%first%cg%ke-1:2, lh)
            fine%first%cg%mg%bnd_y(fine%first%cg%is  :fine%first%cg%ie,     fine%first%cg%ks+1:fine%first%cg%ke  :2, lh) = fine%first%cg%mg%bnd_y(  fine%first%cg%is  :fine%first%cg%ie,     fine%first%cg%ks  :fine%first%cg%ke-1:2, lh)
         endif
         if (fine%first%cg%ext_bnd(zdim, lh)) then
            fine%first%cg%mg%bnd_z(fine%first%cg%is  :fine%first%cg%ie-1:2, fine%first%cg%js  :fine%first%cg%je-1:2, lh) = coarse%first%cg%mg%bnd_z(coarse%first%cg%is:coarse%first%cg%ie,   coarse%first%cg%js:coarse%first%cg%je,   lh)
            fine%first%cg%mg%bnd_z(fine%first%cg%is+1:fine%first%cg%ie  :2, fine%first%cg%js  :fine%first%cg%je-1:2, lh) = fine%first%cg%mg%bnd_z(  fine%first%cg%is  :fine%first%cg%ie-1:2, fine%first%cg%js  :fine%first%cg%je-1:2, lh)
            fine%first%cg%mg%bnd_z(fine%first%cg%is  :fine%first%cg%ie,     fine%first%cg%js+1:fine%first%cg%je  :2, lh) = fine%first%cg%mg%bnd_z(  fine%first%cg%is  :fine%first%cg%ie,     fine%first%cg%js  :fine%first%cg%je-1:2, lh)
         endif
      enddo

   end subroutine prolong_ext_bnd0

!>
!! \brief Prolong boundaries by linear or quadratic interpolation.
!!
!! \details some code is replicated from prolong_faces
!! \todo write something more general, a routine that takes arrays or pointers and does the 2D prolongation
!<

   subroutine prolong_ext_bnd2(coarse)

      use cg_level_connected, only: cg_level_connected_T
      use constants,     only: HI, LO, xdim, ydim, zdim, O_INJ, O_LIN, O_D2, O_I2
      use dataio_pub,    only: die
      use domain,        only: is_multicg

      implicit none

      type(cg_level_connected_T), pointer, intent(in) :: coarse !< level to prolong from

      type(cg_level_connected_T), pointer :: fine

      integer                       :: i, j, k, lh
      real, parameter, dimension(3) :: p0  = [ 0.,       1.,     0.     ] ! injection
      real, parameter, dimension(3) :: p1  = [ 0.,       3./4.,  1./4.  ] ! 1D linear prolongation stencil
      real, parameter, dimension(3) :: p2i = [ -1./8.,   1.,     1./8.  ] ! 1D integral cubic prolongation stencil
      real, parameter, dimension(3) :: p2d = [ -3./32., 30./32., 5./32. ] ! 1D direct cubic prolongation stencil
      real, dimension(-1:1)         :: p
      real, dimension(-1:1,-1:1,2,2):: pp   ! 2D prolongation stencil

      if (is_multicg) call die("[multigridmultipole:prolong_ext_bnd2] multicg not implemented yet") ! fine%first%cg%
      if (.not. associated(coarse)) call die("[multigridmultipole:prolong_ext_bnd0] coarse == null()")
      fine   => coarse%finer
      if (.not. associated(fine)) call die("[multigridmultipole:prolong_ext_bnd0] fine == null()")

      select case (ord_prolong_mpole)
         case (O_INJ)
            p(:) = p0(:)
         case (O_LIN)
            p(:) = p1(:)
         case (O_I2)
            p(:) = p2i(:)
         case (O_D2)
            p(:) = p2d(:)
         case default
            call die("[multigridmultipole:prolong_ext_bnd2] invalid ord_prolong_mpole")
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
         if (fine%first%cg%ext_bnd(xdim, lh)) then
            do j = coarse%first%cg%js+1, coarse%first%cg%je-1
               do k = coarse%first%cg%ks+1, coarse%first%cg%ke-1
                  fine%first%cg%mg%bnd_x(-fine%first%cg%js+2*j,  -fine%first%cg%ks+2*k,  lh) =sum(pp(:,:,1,1) * coarse%first%cg%mg%bnd_x(j-1:j+1,k-1:k+1,lh))
                  fine%first%cg%mg%bnd_x(-fine%first%cg%js+2*j+1,-fine%first%cg%ks+2*k,  lh) =sum(pp(:,:,2,1) * coarse%first%cg%mg%bnd_x(j-1:j+1,k-1:k+1,lh))
                  fine%first%cg%mg%bnd_x(-fine%first%cg%js+2*j,  -fine%first%cg%ks+2*k+1,lh) =sum(pp(:,:,1,2) * coarse%first%cg%mg%bnd_x(j-1:j+1,k-1:k+1,lh))
                  fine%first%cg%mg%bnd_x(-fine%first%cg%js+2*j+1,-fine%first%cg%ks+2*k+1,lh) =sum(pp(:,:,2,2) * coarse%first%cg%mg%bnd_x(j-1:j+1,k-1:k+1,lh))
               enddo
            enddo
         endif

         if (fine%first%cg%ext_bnd(ydim, lh)) then
            do i = coarse%first%cg%is+1, coarse%first%cg%ie-1
               do k = coarse%first%cg%ks+1, coarse%first%cg%ke-1
                  fine%first%cg%mg%bnd_y(-fine%first%cg%is+2*i,  -fine%first%cg%ks+2*k,  lh) =sum(pp(:,:,1,1) * coarse%first%cg%mg%bnd_y(i-1:i+1,k-1:k+1,lh))
                  fine%first%cg%mg%bnd_y(-fine%first%cg%is+2*i+1,-fine%first%cg%ks+2*k,  lh) =sum(pp(:,:,2,1) * coarse%first%cg%mg%bnd_y(i-1:i+1,k-1:k+1,lh))
                  fine%first%cg%mg%bnd_y(-fine%first%cg%is+2*i,  -fine%first%cg%ks+2*k+1,lh) =sum(pp(:,:,1,2) * coarse%first%cg%mg%bnd_y(i-1:i+1,k-1:k+1,lh))
                  fine%first%cg%mg%bnd_y(-fine%first%cg%is+2*i+1,-fine%first%cg%ks+2*k+1,lh) =sum(pp(:,:,2,2) * coarse%first%cg%mg%bnd_y(i-1:i+1,k-1:k+1,lh))
               enddo
            enddo
         endif

         if (fine%first%cg%ext_bnd(zdim, lh)) then
            do i = coarse%first%cg%is+1, coarse%first%cg%ie-1
               do j = coarse%first%cg%js+1, coarse%first%cg%je-1
                  fine%first%cg%mg%bnd_z(-fine%first%cg%is+2*i,  -fine%first%cg%js+2*j,  lh) =sum(pp(:,:,1,1) * coarse%first%cg%mg%bnd_z(i-1:i+1,j-1:j+1,lh))
                  fine%first%cg%mg%bnd_z(-fine%first%cg%is+2*i+1,-fine%first%cg%js+2*j,  lh) =sum(pp(:,:,2,1) * coarse%first%cg%mg%bnd_z(i-1:i+1,j-1:j+1,lh))
                  fine%first%cg%mg%bnd_z(-fine%first%cg%is+2*i,  -fine%first%cg%js+2*j+1,lh) =sum(pp(:,:,1,2) * coarse%first%cg%mg%bnd_z(i-1:i+1,j-1:j+1,lh))
                  fine%first%cg%mg%bnd_z(-fine%first%cg%is+2*i+1,-fine%first%cg%js+2*j+1,lh) =sum(pp(:,:,2,2) * coarse%first%cg%mg%bnd_z(i-1:i+1,j-1:j+1,lh))
               enddo
            enddo
         endif
      enddo

   end subroutine prolong_ext_bnd2

!>
!! \brief Compute multipole moments for image mass
!!
!! \todo distribute excess of work more evenly (important only for large number of PEs, ticket:43)
!<

   subroutine img_mass2moments

      use cg_leaves,    only: leaves
      use cg_list,      only: cg_list_element
      use constants,    only: xdim, ydim, zdim, GEO_XYZ, GEO_RPZ, LO, HI, pSUM, pMIN, pMAX
      use dataio_pub,   only: die
      use domain,       only: dom
      use grid_cont,    only: grid_container
      use mpisetup,     only: piernik_MPI_Allreduce
      use particle_pub, only: pset
      use units,        only: fpiG

      implicit none

      integer :: i, j, k, r, rr
      real, dimension(LO:HI) :: geofac
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      if (dom%geometry_type /= GEO_XYZ .and. any(CoM(xdim:zdim) /= 0.)) call die("[multigridmultipole:img_mass2moments] CoM not allowed for non-cartesian geometry")

      ! reset the multipole data
      Q(:, :, :) = 0.

      irmax = 0
      irmin = rqbin
      geofac(:) = 1.

      !OPT: try to exchange loops i < j < k -> k < j < i
      ! scan
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         if (any(cg%ext_bnd(xdim, :))) then
            if (dom%geometry_type == GEO_RPZ) geofac(:) = [ cg%fbnd(xdim, LO), cg%fbnd(xdim, HI) ]
            do j = cg%js, cg%je
               do k = cg%ks, cg%ke
                  if (cg%leafmap(cg%is, j, k) .and. cg%ext_bnd(xdim, LO) .and. (dom%geometry_type /= GEO_RPZ .or. .not. zaxis_inside)) &
                       call point2moments(cg%mg%bnd_x(j, k, LO)*cg%dyz*geofac(LO), cg%fbnd(xdim, LO)-CoM(xdim), cg%y(j)-CoM(ydim), cg%z(k)-CoM(zdim))
                  if (cg%leafmap(cg%ie, j, k) .and. cg%ext_bnd(xdim, HI)) &
                       call point2moments(cg%mg%bnd_x(j, k, HI)*cg%dyz*geofac(HI), cg%fbnd(xdim, HI)-CoM(xdim), cg%y(j)-CoM(ydim), cg%z(k)-CoM(zdim))
               enddo
            enddo
         endif

         if (any(cg%ext_bnd(ydim, :))) then
            do i = cg%is, cg%ie
               do k = cg%ks, cg%ke
                  if (cg%leafmap(i, cg%js, k) .and. cg%ext_bnd(ydim, LO)) &
                       call point2moments(cg%mg%bnd_y(i, k, LO)*cg%dxz, cg%x(i)-CoM(xdim), cg%fbnd(ydim, LO)-CoM(ydim), cg%z(k)-CoM(zdim))
                  if (cg%leafmap(i, cg%je, k) .and. cg%ext_bnd(ydim, HI)) &
                       call point2moments(cg%mg%bnd_y(i, k, HI)*cg%dxz, cg%x(i)-CoM(xdim), cg%fbnd(ydim, HI)-CoM(ydim), cg%z(k)-CoM(zdim))
               enddo
            enddo
         endif

         if (any(cg%ext_bnd(zdim, :))) then
            do i = cg%is, cg%ie
               if (dom%geometry_type == GEO_RPZ) geofac(LO) = cg%x(i)
               do j = cg%js, cg%je
                  if (cg%leafmap(i, j, cg%ks) .and. cg%ext_bnd(zdim, LO)) &
                       call point2moments(cg%mg%bnd_z(i, j, LO)*cg%dxy*geofac(LO), cg%x(i)-CoM(xdim), cg%y(j)-CoM(ydim), cg%fbnd(zdim, LO)-CoM(zdim))
                  if (cg%leafmap(i, j, cg%ke) .and. cg%ext_bnd(zdim, HI)) &
                       call point2moments(cg%mg%bnd_z(i, j, HI)*cg%dxy*geofac(LO), cg%x(i)-CoM(xdim), cg%y(j)-CoM(ydim), cg%fbnd(zdim, HI)-CoM(zdim))
               enddo
            enddo
         endif
         cgl => cgl%nxt
      enddo

      ! Add only those particles, which are placed outside the domain. Particles inside the domain were already mapped on the grid.
      do i = lbound(pset%p, dim=1), ubound(pset%p, dim=1)
         if (pset%p(i)%outside) call point2moments(fpiG*pset%p(i)%mass, pset%p(i)%pos(xdim), pset%p(i)%pos(ydim), pset%p(i)%pos(zdim))
      enddo

      call piernik_MPI_Allreduce(irmin, pMIN)
      call piernik_MPI_Allreduce(irmax, pMAX)

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

      call piernik_MPI_Allreduce(Q(:, :, irmin:irmax), pSUM)

   end subroutine img_mass2moments

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

!>
!! \brief Compute infinite-boundary potential from multipole moments
!!
!! \todo distribute excess of work more evenly (important only for large number of PEs, ticket:43)
!<

   subroutine moments2bnd_potential

      use cg_leaves,  only: leaves
      use cg_list,    only: cg_list_element
      use constants,  only: xdim, ydim, zdim, GEO_XYZ, GEO_RPZ, LO, HI
      use dataio_pub, only: die
      use domain,     only: dom
      use grid_cont,  only: grid_container

      implicit none

      integer :: i, j, k
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      if (dom%geometry_type /= GEO_XYZ .and. any(CoM(xdim:zdim) /= 0.)) call die("[multigridmultipole:img_mass2moments] CoM not allowed for non-cartesian geometry")

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         if (any(cg%ext_bnd(xdim, :))) then
            do j = cg%js, cg%je
               do k = cg%ks, cg%ke
                  if (cg%ext_bnd(xdim, LO) .and. (dom%geometry_type /= GEO_RPZ .or. .not. zaxis_inside)) &
                       &                    cg%mg%bnd_x(j, k, LO) = moments2pot(cg%fbnd(xdim, LO)-CoM(xdim), cg%y(j)-CoM(ydim), cg%z(k)-CoM(zdim))
                  if (cg%ext_bnd(xdim, HI)) cg%mg%bnd_x(j, k, HI) = moments2pot(cg%fbnd(xdim, HI)-CoM(xdim), cg%y(j)-CoM(ydim), cg%z(k)-CoM(zdim))
               enddo
            enddo
         endif

         if (any(cg%ext_bnd(ydim, :))) then
            do i = cg%is, cg%ie
               do k = cg%ks, cg%ke
                  if (cg%ext_bnd(ydim, LO)) cg%mg%bnd_y(i, k, LO) = moments2pot(cg%x(i)-CoM(xdim), cg%fbnd(ydim, LO)-CoM(ydim), cg%z(k)-CoM(zdim))
                  if (cg%ext_bnd(ydim, HI)) cg%mg%bnd_y(i, k, HI) = moments2pot(cg%x(i)-CoM(xdim), cg%fbnd(ydim, HI)-CoM(ydim), cg%z(k)-CoM(zdim))
               enddo
            enddo
         endif

         if (any(cg%ext_bnd(zdim, :))) then
            do i = cg%is, cg%ie
               do j = cg%js, cg%je
                  if (cg%ext_bnd(zdim, LO)) cg%mg%bnd_z(i, j, LO) = moments2pot(cg%x(i)-CoM(xdim), cg%y(j)-CoM(ydim), cg%fbnd(zdim, LO)-CoM(zdim))
                  if (cg%ext_bnd(zdim, HI)) cg%mg%bnd_z(i, j, HI) = moments2pot(cg%x(i)-CoM(xdim), cg%y(j)-CoM(ydim), cg%fbnd(zdim, HI)-CoM(zdim))
               enddo
            enddo
         endif
         cgl => cgl%nxt
      enddo

   contains

!>
!! \brief Compute potential from multipole moments at a single point
!!
!! \todo improve accuracy with linear interpolation over radius
!<

      real function moments2pot(x, y, z) result(potential)

         use units, only: newtong

         implicit none

         real, intent(in)  :: x         !< x coordinate of the contributing point
         real, intent(in)  :: y         !< y coordinate of the contributing point
         real, intent(in)  :: z         !< z coordinate of the contributing point

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

      end function moments2pot

   end subroutine moments2bnd_potential

!>
!! \brief This routine calculates various geometrical numbers required for multipole evaluation
!!
!! \details It modifies the rn(:), irn(:), cfac(:) and sfac(:) arrays. Scalars are passed through argument list.
!<

   subroutine geomfac4moments(factor, x, y, z, sin_th, cos_th, ir, delta)

      use constants,  only: GEO_XYZ, GEO_RPZ
      use dataio_pub, only: die, msg
      use domain,     only: dom

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
      select case (dom%geometry_type)
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
         write(msg,'(2(a,i7),a)')"[multigridmultipole:geomfac4moments] radial index = ",ir," outside Q(:, :, ",rqbin,") range"
         call die(msg)
      endif
      irmax = max(irmax, ir)
      irmin = min(irmin, ir)

      ! azimuthal angle sine and cosine tables
      ! ph = atan2(y, x); cfac(m) = cos(m * ph); sfac(m) = sin(m * ph)
      ! cfac(0) and sfac(0) are set in init_multigrid
      if (rxy /= 0.) then
         select case (dom%geometry_type)
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

!>
!! \brief HEAVY_DEBUG marks routines that normally are never called, but at some point
!! were useful to test correctness or something.
!<

!#define HEAVY_DEBUG
#ifdef HEAVY_DEBUG

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
