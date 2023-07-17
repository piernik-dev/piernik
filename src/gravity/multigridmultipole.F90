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
!! \details By default this solver estimates gravitational potential on external (domain) boundaries, which allows to mimic \f$\Phi(\infty) = 0\f$.
!! The calculated potential may then be used in second pass of the Poisson solver (this time run with an empty space) to calculate correction
!! to be added to the first-pass solution obtained with homogeneous Dirichlet boundary conditions.
!!
!! One can also use purely monopole estimate of the potential at the boundaries, provide an user routine that will impose some values.
!!
!! There is also a 3D multipole solver that does not depend on differentiation of the gravitational potential atdomain boundaries.
!! Note that differentiation still takes place in init_source to convert given-value potential into density - this is the biggest source of solution inaccuracy when everything else is set up to high orders.
!!
!! The 3D multipole solver is typically slower from the default, image mass method but is way more accurate for density distributions that doesn't have much contribution near boundaries.
!! It is possible, however, to compure is purely on base level or even coarser levels without much loss of accuracy.
!!
!! Another advantage of the 3D solver is that it usually doesn't benefit much from high-order multipoles so setting l_max to 8 or even 6 usually should be sufficient.
!! In the same situation the image mass solver tends to produce long tails for high l on both small m and m close to m_max.
!!
!! ToDo
!! * Automatic limiting of l_max based on relative strength w.r.t monopole contribution (1e-6 seems to be reasonable threshold).
!! * (3D solver) Automatic choice of best level to integrate on, based on "distance" between Q calculated on different levels.
!! * "True" given-value boundaries in multigrid, without differentiation.
!<

module multipole
! pulled by MULTIGRID && SELF_GRAV

   ! needed for global vars in this module
   use constants,       only: cbuff_len
   use multipole_array, only: mpole_container

#if defined(__INTEL_COMPILER)
   !! \deprecated remove this clause as soon as Intel Compiler gets required
   !! features and/or bug fixes
   use cg_list_bnd,   only: cg_list_bnd_t   ! QA_WARN intel
#endif /* __INTEL_COMPILER */

   implicit none

   private
   public :: init_multipole, cleanup_multipole, multipole_solver, moments2pot, compute_mpole_potential
   public :: lmax, mmax, mpole_solver, singlepass  ! initialized in multigrid_gravity

   interface moments2pot
      module procedure moments2pot_xyz
      module procedure moments2pot_r
   end interface moments2pot

   type(mpole_container)    :: Q             !< The whole moment array with dependence on radius

   ! namelist parameters for MULTIGRID_GRAVITY
   integer(kind=4)          :: lmax          !< Maximum l-order of multipole moments
   integer(kind=4)          :: mmax          !< Maximum m-order of multipole moments. Equal to lmax by default.
   character(len=cbuff_len) :: mpole_solver  !< Pick one of: "monopole", "img_mass" (default), "3D"

   logical                  :: zaxis_inside  !< true when z-axis belongs to the inner radial boundary in polar coordinates
   logical                  :: singlepass    !< When .true. it allows for single-pass multigrid solve

   enum, bind(C)
      enumerator :: MONOPOLE, IMG_MASS, THREEDIM
   end enum
   integer                  :: solver        !< mpole_solver decoded into one of the above enums

!> \todo OPT derive a special set of leaves and coarsened leaves that can be safely used here

contains

!>
!! \brief Initialization routine, called once, from init_multigrid
!!
!! Single-pass multigrid allowed only for user and 3D solvers as these doesn't depend on Dirichlet solution.
!! BEWARE: if user solution decides to depend on Dirichlet solution, this has to be changed.
!<

   subroutine init_multipole

      use constants,     only: ndims, INVALID
      use dataio_pub,    only: die, warn
      use domain,        only: dom
      use mpisetup,      only: master
      use multigridvars, only: grav_bnd, bnd_isolated
      use user_hooks,    only: ext_bnd_potential

      implicit none

      if (dom%eff_dim /= ndims) call die("[multigridmultipole:init_multipole] Only 3D is supported") !> \todo add support for 2D RZ
      if (grav_bnd /= bnd_isolated) call die("[multigridmultipole:init_multipole] Only fully isolated boundaries are implemented")

      !fixup multipole moments
      if (mmax > lmax) then
         if (master) call warn("[multigridmultipole:init_multipole] mmax reduced to lmax")
         mmax = lmax
      endif
      if (mmax < 0) mmax = lmax

      singlepass = associated(ext_bnd_potential)
      solver = INVALID
      select case (mpole_solver)
         case ("monopole", "mono")
            solver = MONOPOLE
         case ("img_mass", "surface", "dOmega")
            solver = IMG_MASS
         case ("3D", "3d", "volume")
            solver = THREEDIM
            singlepass = .true.
         case default
            call die("[multigridmultipole:init_multipole] unknown solver '" // trim(mpole_solver) // "'")
      end select
      if (solver /= MONOPOLE) call Q%init_once(lmax, mmax)

   end subroutine init_multipole

!>
!! \brief Initialization routine, called once per entry to the multipole solver.
!!
!! \details This routine reinitializes everything that may need reinitialization due to AMR activity
!!
!! Since the worst quality of multipole expansion occurs at the radius of the contributing source,
!! typically it is safest choice to use center of domain as the origin of radial distribution of multipoles.
!!
!! At the sphere around the origin that contains the source, the multipole expansion reduces to
!! spherical harmonics and there are no factors containing powers of r to help with convergence.
!!
!! \todo OPT: Detect only changes in highest required level
!<

   subroutine refresh_multipole

      use cg_level_finest, only: finest
      use constants,       only: small, pi, xdim, ydim, zdim, GEO_XYZ, GEO_RPZ, LO
      use dataio_pub,      only: die, warn
      use domain,          only: dom
      use mpisetup,        only: master

      implicit none

      select case (dom%geometry_type)
         case (GEO_XYZ)
            Q%center(xdim:zdim) = dom%C_(xdim:zdim)
            zaxis_inside = .false.
         case (GEO_RPZ)
            if (dom%L_(ydim) >= (2.-small)*pi) then
               Q%center(xdim) = 0.
               Q%center(ydim) = 0.
            else
!!$               Q%center(xdim) = 2./3. * (dom%edge(xdim, HI)**3-dom%edge(xdim, LO)**3)/(dom%edge(xdim, HI)**2-dom%edge(xdim, LO)**2)
!!$               if (dom%L_(ydim).notequals.zero) Q%center(xdim) = Q%center(xdim) * sin(dom%L_(ydim)/2.)/(dom%L_(ydim)/2.)
!!$               Q%center(ydim) = dom%C_(ydim)
               Q%center(xdim) = 0.
               Q%center(ydim) = 0.
            endif
            Q%center(zdim) = dom%C_(zdim)
            zaxis_inside = dom%edge(xdim, LO) <= dom%L_(xdim)/finest%level%l%n_d(xdim)
            if (master) then
               if (zaxis_inside) call warn("[multigridmultipole:refresh_multipole] Setups with Z-axis at the edge of the domain may not work as expected yet.")
               if (solver == MONOPOLE) call warn("[multigridmultipole:refresh_multipole] Point-like monopole is not implemented.")
            endif
            solver = IMG_MASS
         case default
            call die("[multigridmultipole:refresh_multipole] Unsupported geometry.")
      end select

      if (solver /= MONOPOLE) call Q%refresh

   end subroutine refresh_multipole

!> \brief Multipole cleanup

   subroutine cleanup_multipole

      implicit none

      call Q%cleanup

   end subroutine cleanup_multipole

!>
!! \brief Multipole solver
!!
!! \todo improve multipole expansion on coarser grids
!! (see. "A Scalable Parallel Poisson Solver in Three Dimensions with Infinite-Domain Boundary Conditions" by McCorquodale, Colella, Balls and Baden).
!<

   subroutine multipole_solver

      use cg_leaves,     only: leaves
      use constants,     only: dirtyH1, PPP_GRAV
      use dataio_pub,    only: die
      use global,        only: dirty_debug
      use mg_monopole,   only: isolated_monopole, find_img_CoM
      ! use multigridvars, only: grav_bnd, bnd_isolated
      use ppp,           only: ppp_main
      use user_hooks,    only: ext_bnd_potential

      implicit none

      character(len=*), parameter :: mpole_label = "multipole_solver"

      call ppp_main%start(mpole_label, PPP_GRAV)

      ! Cannot use this check, because we do tricks in multigrid_solve_grav
      ! if (grav_bnd /= bnd_isolated) call die("[multigridmultipole:multipole_solver] Only fully isolated boundaries are implemented")

      if (associated(ext_bnd_potential)) then
         call ext_bnd_potential
         return
      endif

      call refresh_multipole

      if (dirty_debug) then
         call leaves%reset_boundaries(0.959*dirtyH1)
      else
         call leaves%reset_boundaries
      endif

      select case (solver)
         case (MONOPOLE)
            call potential2img_mass
            call find_img_CoM
            call isolated_monopole
         case (IMG_MASS)
            call potential2img_mass
            ! results seems to be slightly better without find_img_CoM.
            ! With CoM or when it is known than CoM is close to the domain center one may try to save some CPU time by lowering mmax.
            call Q%reset
            call img_mass2moments
#ifdef NBODY
            call particles2moments
#endif /* NBODY */
            call Q%red_int_norm
            ! OPT: automagically reduce lmax for a couple of steps if higher multipoles fall below some threshold
            call moments2bnd_potential
         case (THREEDIM)
            call Q%reset
            call domain2moments
#ifdef NBODY
            call particles2moments
#endif /* NBODY */
            call Q%red_int_norm
            ! OPT: automagically reduce lmax for a couple of steps if higher multipoles fall below some threshold
            call moments2bnd_potential
         case default
            call die("[multigridmultipole:multipole_solver] unimplemented solver")
      end select

      call ppp_main%stop(mpole_label, PPP_GRAV)

   end subroutine multipole_solver

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

      use cg_cost_data,  only: I_MULTIPOLE
      use cg_leaves,     only: leaves
      use cg_list,       only: cg_list_element
      use constants,     only: GEO_RPZ, LO, HI, xdim, ydim, zdim, PPP_GRAV
      use domain,        only: dom
      use grid_cont,     only: grid_container
      use multigridvars, only: solution
      use ppp,           only: ppp_main

      implicit none

      integer :: i
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg
      real, parameter :: a1 = -2., a2 = (-2. - a1)/3. ! interpolation parameters;   <---- a1=-2 => a2=0
      ! a1 = -2. is the simplest, 1st order choice, gives best agreement of total mass and CoM location when compared to 3-D integration
      ! a1 = -1., a2 = -1./3. seems to do the best job,
      ! a1 = -3./2., a2 = -1./6. seems to be 2nd order estimator
      character(len=*), parameter :: p2m_label = "multipole_pot2img"

      call ppp_main%start(p2m_label, PPP_GRAV)
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         call cg%costs%start

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

         call cg%costs%stop(I_MULTIPOLE)
         cgl => cgl%nxt
         end associate
      enddo
      call ppp_main%stop(p2m_label, PPP_GRAV)

   end subroutine potential2img_mass

!>
!! \brief Compute multipole moments for image mass
!!
!! \todo distribute excess of work more evenly (important only for large number of PEs, ticket:43)
!<

   subroutine img_mass2moments

      use cg_cost_data, only: I_MULTIPOLE
      use cg_leaves,    only: leaves
      use cg_list,      only: cg_list_element
      use constants,    only: xdim, ydim, zdim, GEO_XYZ, GEO_RPZ, LO, HI, PPP_GRAV
      use dataio_pub,   only: die
      use domain,       only: dom
      use grid_cont,    only: grid_container
      use ppp,          only: ppp_main

      implicit none

      integer :: i, j, k
      real, dimension(LO:HI) :: geofac
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg
      character(len=*), parameter :: m2m_label = "multipole_img2mom"

      call ppp_main%start(m2m_label, PPP_GRAV)
      if (dom%geometry_type /= GEO_XYZ .and. any(abs(Q%center(xdim:zdim)) > tiny(1.))) call die("[multigridmultipole:img_mass2moments] Q%center /= 0. not implemented for non-cartesian geometry")

      geofac(:) = 1.

      !OPT: try to exchange loops i < j < k -> k < j < i
      ! scan
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         call cg%costs%start

         if (any(cg%ext_bnd(xdim, :))) then
            if (dom%geometry_type == GEO_RPZ) geofac(:) = [ cg%fbnd(xdim, LO), cg%fbnd(xdim, HI) ]
            do j = cg%js, cg%je
               do k = cg%ks, cg%ke
                  if (cg%leafmap(cg%is, j, k) .and. cg%ext_bnd(xdim, LO) .and. (dom%geometry_type /= GEO_RPZ .or. .not. zaxis_inside)) &
                       call Q%point2moments(cg%mg%bnd_x(j, k, LO)*cg%dyz*geofac(LO), cg%fbnd(xdim, LO), cg%y(j), cg%z(k))
                  if (cg%leafmap(cg%ie, j, k) .and. cg%ext_bnd(xdim, HI)) &
                       call Q%point2moments(cg%mg%bnd_x(j, k, HI)*cg%dyz*geofac(HI), cg%fbnd(xdim, HI), cg%y(j), cg%z(k))
               enddo
            enddo
         endif

         if (any(cg%ext_bnd(ydim, :))) then
            do i = cg%is, cg%ie
               do k = cg%ks, cg%ke
                  if (cg%leafmap(i, cg%js, k) .and. cg%ext_bnd(ydim, LO)) &
                       call Q%point2moments(cg%mg%bnd_y(i, k, LO)*cg%dxz, cg%x(i), cg%fbnd(ydim, LO), cg%z(k))
                  if (cg%leafmap(i, cg%je, k) .and. cg%ext_bnd(ydim, HI)) &
                       call Q%point2moments(cg%mg%bnd_y(i, k, HI)*cg%dxz, cg%x(i), cg%fbnd(ydim, HI), cg%z(k))
               enddo
            enddo
         endif

         if (any(cg%ext_bnd(zdim, :))) then
            do i = cg%is, cg%ie
               if (dom%geometry_type == GEO_RPZ) geofac(LO) = cg%x(i)
               do j = cg%js, cg%je
                  if (cg%leafmap(i, j, cg%ks) .and. cg%ext_bnd(zdim, LO)) &
                       call Q%point2moments(cg%mg%bnd_z(i, j, LO)*cg%dxy*geofac(LO), cg%x(i), cg%y(j), cg%fbnd(zdim, LO))
                  if (cg%leafmap(i, j, cg%ke) .and. cg%ext_bnd(zdim, HI)) &
                       call Q%point2moments(cg%mg%bnd_z(i, j, HI)*cg%dxy*geofac(LO), cg%x(i), cg%y(j), cg%fbnd(zdim, HI))
               enddo
            enddo
         endif

         call cg%costs%stop(I_MULTIPOLE)
         cgl => cgl%nxt
      enddo
      call ppp_main%stop(m2m_label, PPP_GRAV)

   end subroutine img_mass2moments

!>
!! \brief Compute multipole moments for the whole domain
!!
!! \todo test with CoM (implement find_CoM)
!! \todo implement Richardson extrapolation between coarsened levels
!<

   subroutine domain2moments

      use cg_cost_data,       only: I_MULTIPOLE
      use cg_leaves,          only: leaves
      use cg_level_base,      only: base
      use cg_level_finest,    only: finest
      use cg_level_connected, only: cg_level_connected_t
      use cg_list,            only: cg_list_element
      use constants,          only: xdim, zdim, GEO_XYZ, base_level_id, PPP_GRAV
      use dataio_pub,         only: die, msg, warn
      use domain,             only: dom
      use grid_cont,          only: grid_container
      use multigridvars,      only: source
      use multipole_array,    only: mpole_level, mpole_level_auto
      use ppp,                only: ppp_main

      implicit none

      integer :: i, j, k
      type(cg_level_connected_t), pointer :: level
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg
      character(len=*), parameter :: d2m_label = "multipole_dom2mom"

      call ppp_main%start(d2m_label, PPP_GRAV)

      if (dom%geometry_type /= GEO_XYZ .and. any(abs(Q%center(xdim:zdim)) > tiny(1.))) call die("[multigridmultipole:domain2moments] Q%center /= 0. not implemented for non-cartesian geometry")
      if (dom%geometry_type /= GEO_XYZ) call die("[multigridmultipole:domain2moments] Noncartesian geometry haven't been tested. Verify it before use.")

      ! scan
      if (mpole_level <= mpole_level_auto) then
         level => finest%find_finest_bnd()
         cgl => level%first
      else if (mpole_level <= base_level_id) then
         level => base%level
         call finest%level%restrict_to_base_q_1var(source)
         do while (level%l%id > mpole_level)
            if (associated(level%coarser)) then
               call level%restrict_1var(source)
               level => level%coarser
            else
               write(msg, '(2(a,i3))')"[multigridmultipole:domain2moments] Coarsest level reached. Will use level ", level%l%id, " instead of ", mpole_level
               call warn(msg)
               exit
            endif
         enddo
         cgl => level%first
      else
         level => finest%level
         do while (level%l%id > mpole_level)
            call level%restrict_1var(source)
            level => level%coarser
            if (.not. associated(level)) call die("[multigridmultipole:domain2moments] Coarsest level reached for mpole_level > base_level_id")
         enddo
         cgl => leaves%first
      endif

      do while (associated(cgl))
         cg => cgl%cg
         call cg%costs%start

         if (cg%l%id <= level%l%id) then
            do k = cg%ks, cg%ke
               do j = cg%js, cg%je
                  do i = cg%is, cg%ie
                     ! if (dom%geometry_type == GEO_RPZ) geofac = cg%x(i)
                     if (cg%leafmap(i, j, k) .or. level%l%id == cg%l%id) &
                          call Q%point2moments(cg%dvol * cg%q(source)%arr(i, j, k), cg%x(i) , cg%y(j) , cg%z(k) )  ! * geofac for GEO_RPZ
                  enddo
               enddo
            enddo
         endif

         call cg%costs%stop(I_MULTIPOLE)
         cgl => cgl%nxt
      enddo
      call ppp_main%stop(d2m_label, PPP_GRAV)

   end subroutine domain2moments

!> \brief Compute multipole moments for the particles

#ifdef NBODY
   subroutine particles2moments

      use cg_cost_data,   only: I_PARTICLE
      use cg_leaves,      only: leaves
      use cg_list,        only: cg_list_element
      use constants,      only: xdim, ydim, zdim, PPP_GRAV, PPP_PART
      use grid_cont,      only: grid_container
      use ppp,            only: ppp_main
      use units,          only: fpiG
      use particle_types, only: particle

      implicit none

      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer  :: cg
      type(particle), pointer        :: pset
      character(len=*), parameter    :: p2m_label = "multipole_part2mom"

      call ppp_main%start(p2m_label, PPP_GRAV + PPP_PART)

      ! Add only those particles, which are placed outside the domain. Particles inside the domain were already mapped on the grid.
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         call cg%costs%start

         pset => cg%pset%first
         do while (associated(pset))
            if (pset%pdata%phy .and. pset%pdata%outside) &
                 call Q%point2moments(fpiG*pset%pdata%mass, pset%pdata%pos(xdim), pset%pdata%pos(ydim), pset%pdata%pos(zdim))

            ! When the particles that are outside the domain happen to be close to projected exension od the cg-cg boundary,
            ! they were counted twice. Since we use the phy flag for counting unique particles, it should work also here
            ! to select the only ghost copy of the outside particle for multipole contribution.

            ! WARNING: Particles that are too close to the outer boundary aren't fully mapped onto the grid.
            ! This may cause huge errors in the potential, even for "3D" solver because their mass is counted
            ! only partially in the mapping routine and the rest is ignored here.
            !
            ! Suggested fix: set up some buffer at the boundaries (at least one cell wide, few cells are advised)
            ! and make a "not_mapped" flag instead of "outside". The transition between mapped and not mapped
            ! particles may be smooth (with partial mappings allowed) but it must not allow for
            ! particle mappings extending beyond domain boundaries.
            !
            ! A "not_mapped" flag set in the mapping routine would fix this issue.

            pset => pset%nxt
         enddo

         call cg%costs%stop(I_PARTICLE)
         cgl => cgl%nxt
      enddo

      call ppp_main%stop(p2m_label, PPP_GRAV + PPP_PART)

   end subroutine particles2moments
#endif /* NBODY */

!>
!! \brief Compute infinite-boundary potential from multipole moments
!!
!! \todo distribute excess of work more evenly (important only for large number of PEs, ticket:43)
!<

   subroutine moments2bnd_potential

      use cg_cost_data, only: I_MULTIPOLE
      use cg_leaves,    only: leaves
      use cg_list,      only: cg_list_element
      use constants,    only: xdim, ydim, zdim, GEO_XYZ, GEO_RPZ, LO, HI, PPP_GRAV
      use dataio_pub,   only: die
      use domain,       only: dom
      use grid_cont,    only: grid_container
      use ppp,          only: ppp_main

      implicit none

      integer :: i, j, k
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg
      character(len=*), parameter :: m2p_label = "multipole_mom2pot"

      call ppp_main%start(m2p_label, PPP_GRAV)

      if (dom%geometry_type /= GEO_XYZ .and. any(abs(Q%center(xdim:zdim)) > tiny(1.))) call die("[multigridmultipole:img_mass2moments] Q%center /= 0. not implemented for non-cartesian geometry")

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         call cg%costs%start

         if (any(cg%ext_bnd(xdim, :))) then
            do j = cg%js, cg%je
               do k = cg%ks, cg%ke
                  if (cg%ext_bnd(xdim, LO) .and. (dom%geometry_type /= GEO_RPZ .or. .not. zaxis_inside)) &
                       &                    cg%mg%bnd_x(j, k, LO) = Q%moments2pot(cg%fbnd(xdim, LO), cg%y(j), cg%z(k))
                  if (cg%ext_bnd(xdim, HI)) cg%mg%bnd_x(j, k, HI) = Q%moments2pot(cg%fbnd(xdim, HI), cg%y(j), cg%z(k))
               enddo
            enddo
         endif

         if (any(cg%ext_bnd(ydim, :))) then
            do i = cg%is, cg%ie
               do k = cg%ks, cg%ke
                  if (cg%ext_bnd(ydim, LO)) cg%mg%bnd_y(i, k, LO) = Q%moments2pot(cg%x(i), cg%fbnd(ydim, LO), cg%z(k))
                  if (cg%ext_bnd(ydim, HI)) cg%mg%bnd_y(i, k, HI) = Q%moments2pot(cg%x(i), cg%fbnd(ydim, HI), cg%z(k))
               enddo
            enddo
         endif

         if (any(cg%ext_bnd(zdim, :))) then
            do i = cg%is, cg%ie
               do j = cg%js, cg%je
                  if (cg%ext_bnd(zdim, LO)) cg%mg%bnd_z(i, j, LO) = Q%moments2pot(cg%x(i), cg%y(j), cg%fbnd(zdim, LO))
                  if (cg%ext_bnd(zdim, HI)) cg%mg%bnd_z(i, j, HI) = Q%moments2pot(cg%x(i), cg%y(j), cg%fbnd(zdim, HI))
               enddo
            enddo
         endif

         call cg%costs%stop(I_MULTIPOLE)
         cgl => cgl%nxt
      enddo
      call ppp_main%stop(m2p_label, PPP_GRAV)

   end subroutine moments2bnd_potential

!> \brief Wrapper for Q%moments2pot with proper normalisation for (x, y, z) argument

   real function moments2pot_xyz(x, y, z) result(pot)

      use dataio_pub,    only: die
      use multigridvars, only: grav_bnd, bnd_isolated
      use units,         only: fpiG

      implicit none

      real, intent(in) :: x !< absolute x-coordinate of the point
      real, intent(in) :: y !< absolute y-coordinate of the point
      real, intent(in) :: z !< absolute z-coordinate of the point

      if (grav_bnd /= bnd_isolated) call die("[multigridmultipole:moments2pot_xyz] Only fully isolated boundaries are implemented")

      pot = Q%moments2pot(x, y, z)/fpiG

   end function moments2pot_xyz

!> \brief Wrapper for Q%moments2pot with proper normalisation for [x, y, z] argument

   real function moments2pot_r(r) result(pot)

      use constants,     only: xdim, ydim, zdim, ndims
      use dataio_pub,    only: die
      use multigridvars, only: grav_bnd, bnd_isolated
      use units,         only: fpiG

      implicit none

      real, dimension(ndims), intent(in) :: r !< [x, y, z] coordinate of the point

      if (grav_bnd /= bnd_isolated) call die("[multigridmultipole:moments2pot_r] Only fully isolated boundaries are implemented")

      pot = Q%moments2pot(r(xdim), r(ydim), r(zdim))/fpiG

   end function moments2pot_r

!>
!! \brief Compute multipole potential in whole computational domain
!!
!! \details This routine is not intended for regular use because of huge
!! computational cost per cell (O(lmax**2)). It is useful for understanding
!! and debugging the multipole solver and for showing properties of multipole
!! estimate of gravitational potential.
!<

   subroutine compute_mpole_potential(qvar)

      use cg_cost_data,  only: I_MULTIPOLE
      use cg_leaves,     only: leaves
      use cg_list,       only: cg_list_element
      use multigridvars, only: grav_bnd, bnd_isolated

      implicit none

      integer(kind=4), intent(in) :: qvar  !< index in cg%q to store the results

      integer :: i, j, k
      type(cg_list_element), pointer :: cgl

      call leaves%set_q_value(qvar, 0.)

      if (grav_bnd /= bnd_isolated) return

      cgl => leaves%first
      do while (associated(cgl))
         ! an ELEMENTAL implementation of moments2pot would allow simplifications of this routine.
         ! At least the loop over i may be worth stripping.
         associate (cg => cgl%cg)
            call cg%costs%start

            do k = cg%ks, cg%ke
               do j = cg%js, cg%je
                  do i = cg%is, cg%ie
                     cg%q(qvar)%arr(i, j, k) = moments2pot([cg%x(i), cg%y(j), cg%z(k)])
                  enddo
               enddo
            enddo

            call cg%costs%stop(I_MULTIPOLE)
         end associate
         cgl => cgl%nxt
      enddo

   end subroutine compute_mpole_potential

end module multipole
