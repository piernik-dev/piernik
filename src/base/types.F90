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
!! \brief Definitions of compound types and subroutine templates for user hooks
!<
module types

   use constants, only: ndims, LO, HI, xdim, zdim

   implicit none

   private
   public :: axes, domain_container, segment, tsl_container, value, problem_customize_solution, problem_grace_passed, finalize_problem, cleanup_problem, custom_emf_bnd

   type :: value
      real                      :: val
      integer, dimension(ndims) :: loc
      integer                   :: proc
   end type value

! AMR: There will be at least one domain container for the base grid. It will be possible to host one or more refined domains on the base container and on the refined containers.
!      The refined domains may cover whole parent domain or only a tiny part of it.
!      The multigrid solver will operate also on a stack of coarser domains - parents of the base domain. The coarser domains must be no smaller than the base domain.
   type :: domain_container
      ! primary parameters, read from /DOMAIN_SIZES/, /BOUNDARIES/ and /DOMAIN_LIMITS/ namelists
      real    :: xmin                           !< physical domain left x-boundary position
      real    :: xmax                           !< physical domain right x-boundary position
      real    :: ymin                           !< physical domain left y-boundary position
      real    :: ymax                           !< physical domain right y-boundary position
      real    :: zmin                           !< physical domain left z-boundary position
      real    :: zmax                           !< physical domain right z-boundary position
      integer, dimension(ndims) :: n_d          !< number of grid cells in physical domain in x-, y- and z-direction (where equal to 1, the dimension is reduced to a point with no boundary cells)
      integer                   :: nb           !< number of boundary cells surrounding the physical domain, same for all directions
      integer, dimension(ndims, LO:HI) :: bnd   !< type of boundary conditions coded in integers

      integer, dimension(:,:,:), allocatable :: se !< span of all domains (0:nproc-1, xdim:zdim, LO:HI); Use with care, because this is an antiparallel thing
      !! \todo add hooks to parent and a list of children domains

      ! derived parameters
      real    :: Lx                             !< span of the physical domain in x-direction (xmax-xmin)
      real    :: Ly                             !< span of the physical domain in y-direction (ymax-ymin)
      real    :: Lz                             !< span of the physical domain in z-direction (zmax-zmin)
      real    :: x0                             !< center of the physical domain in x-direction (xmax+xmin)/2.
      real    :: y0                             !< center of the physical domain in y-direction (ymax+ymin)/2.
      real    :: z0                             !< center of the physical domain in z-direction (zmax+zmin)/2.
      real    :: Vol                            !< total volume of the physical domain

      logical, dimension(ndims) :: periodic     !< .true. for periodic and shearing boundary pairs

      ! Do not use n[xyz]t components in the Piernik source tree without a good reason
      ! Avoid as a plague allocating buffers of that size because it negates benefits of parallelization
      ! \todo move them to another type, that extends domain_container?
      integer :: nxt                            !< total number of %grid cells in the whole domain in x-direction
      integer :: nyt                            !< total number of %grid cells in the whole domain in y-direction
      integer :: nzt                            !< total number of %grid cells in the whole domain in z-direction

    contains

      procedure :: set_derived

   end type domain_container

   ! specify segment of data for boundary exchange, prolongation and restriction.
   type :: segment
      integer :: proc                               !< target process
      integer, dimension(xdim:zdim, LO:HI) :: se    !< range
      integer, dimension(:), allocatable   :: mc    !< MPI type Containers
   end type segment

   type :: axes
      real, allocatable, dimension(:) :: x      !< array of x-positions of %grid cells centers
      real, allocatable, dimension(:) :: y      !< array of y-positions of %grid cells centers
      real, allocatable, dimension(:) :: z      !< array of z-positions of %grid cells centers
   end type axes

   type :: tsl_container
      logical :: dummy
#ifdef COSM_RAYS
      real :: encr_min, encr_max
#endif /* COSM_RAYS */

#ifdef RESISTIVE
      real :: etamax
#endif /* RESISTIVE */

#ifdef MAGNETIC
      real :: b_min, b_max, divb_max
#endif /* MAGNETIC */

#ifdef IONIZED
      real :: vai_max
#endif /* IONIZED */

#ifdef VARIABLE_GP
      real :: gpxmax, gpymax, gpzmax
#endif /* VARIABLE_GP */

   end type tsl_container

   ! User hooks

   interface
      subroutine no_args
         implicit none
      end subroutine no_args
      subroutine tab_args(tab)
         implicit none
         real, dimension(:,:,:), intent(inout) :: tab
      end subroutine tab_args
   end interface

   procedure(no_args),  pointer :: problem_customize_solution => NULL()
   procedure(no_args),  pointer :: problem_grace_passed       => NULL()
   procedure(no_args),  pointer :: finalize_problem           => NULL()
   procedure(no_args),  pointer :: cleanup_problem            => NULL()
   procedure(tab_args), pointer :: custom_emf_bnd             => NULL()

contains

   subroutine set_derived(this)

      use constants,  only: xdim, ydim, zdim, LO, HI, BND_PER, BND_SHE
      use dataio_pub, only: die

      implicit none

      class(domain_container), intent(inout) :: this
      logical, dimension(xdim:zdim) :: has_dir
      integer :: d

      has_dir(:) = this%n_d(:) > 1

      this%periodic(:) = .false.
      do d = xdim, zdim
         if ((any(this%bnd(d, LO:HI) == BND_PER) .or. (d==xdim .and. any(this%bnd(d, LO:HI) == BND_SHE))) .and. has_dir(d)) then
            this%periodic(d) = .true.
            if (this%bnd(d, LO) /= this%bnd(d, HI)) call die("[types:set_derived] Periodic BC do not match")
         endif
      enddo
      if (any(this%bnd(ydim:zdim, LO:HI) == BND_SHE)) call die("[types:set_derived] Shearing BC not allowed for y- and z-direction")

      !> \todo convert L[xyz], [xyz]0 and n[xyz]t to (xdim:zdim)-sized arrays

      ! auxiliary lengths
      this%Lx = this%xmax - this%xmin
      this%Ly = this%ymax - this%ymin
      this%Lz = this%zmax - this%zmin
      this%x0 = (this%xmax + this%xmin)/2.
      this%y0 = (this%ymax + this%ymin)/2.
      this%z0 = (this%zmax + this%zmin)/2.

      !volume and total grid sizes
      this%Vol = 1.
      if (has_dir(xdim)) then
         this%Vol = this%Vol * this%Lx
         this%nxt = this%n_d(xdim) + 2 * this%nb
      else
         this%nxt = 1
      endif

      if (has_dir(ydim)) then
         this%Vol = this%Vol * this%Ly
         this%nyt = this%n_d(ydim) + 2 * this%nb
      else
         this%nyt = 1
      endif

      if (has_dir(zdim)) then
         this%Vol = this%Vol * this%Lz
         this%nzt = this%n_d(zdim) + 2 * this%nb
      else
         this%nzt = 1
      endif
      !> \deprecated BEWARE: Vol computed above is not true for non-cartesian geometry

   end subroutine set_derived

end module types
