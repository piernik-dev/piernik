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
   public :: axes, domain_container, segment, bnd_list, tsl_container, value, array4d, array3d, &
        &    problem_customize_solution, problem_grace_passed, finalize_problem, cleanup_problem, custom_emf_bnd, at_user_settings

   type :: array4d
      real, dimension(:,:,:,:), pointer :: arr => null()
      contains
         procedure :: array4d_init           ! \todo check why private here does not  work as expected
         procedure :: array4d_associate
         procedure :: clean => array4d_clean
         procedure :: array4d_get_sweep
         procedure :: array4d_get_sweep_one_var
         procedure :: check => array4d_check_if_dirty
         generic, public :: init => array4d_init, array4d_associate
         generic, public :: get_sweep => array4d_get_sweep_one_var, array4d_get_sweep
   end type array4d

   type :: array3d
      real, dimension(:,:,:), pointer :: arr => null()
      contains
         procedure :: array3d_init
         procedure :: array3d_associate
         procedure :: clean => array3d_clean
         procedure :: check => array3d_check_if_dirty
         procedure :: get_sweep => array3d_get_sweep
         generic, public :: init => array3d_init, array3d_associate
   end type array3d

   type :: value
      real                      :: val
      real,    dimension(ndims) :: coords
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

      integer(kind=8), dimension(:,:,:), allocatable :: se !< span of all domains (0:nproc-1, xdim:zdim, LO:HI); Use with care, because this is an antiparallel thing
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

      integer :: pdiv_type                      !< 0 for cartesian decomposition of base level, /= 0 for more spthisticated patterns
      integer, dimension(xdim:zdim) :: pdiv     !< when pdiv_type == 0 it is like psize(:)

    contains

      procedure :: set_derived

   end type domain_container

   ! specify segment of data for boundary exchange, prolongation and restriction.
   type :: segment
      integer :: proc                                     !< target process
      integer(kind=8), dimension(xdim:zdim, LO:HI) :: se  !< range
   end type segment

   !< segment type for boundary exchange
   type, extends(segment) :: bnd_segment
      integer :: mbc                                      !< Multigrid MPI Boundary conditions Container
      integer :: lh                                       !< low or high boundary; \todo store full tag here
   end type bnd_segment

   type :: bnd_list
      type(bnd_segment), dimension(:), allocatable :: seg !< list of boundary segments to exchange
   end type bnd_list

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
      subroutine indx_args(ar,ll,lr,ch,lo)
         implicit none
         integer,         dimension(:), intent(out) :: ar, ll, lr, ch
         integer(kind=8), dimension(:), intent(out) :: lo
      end subroutine indx_args
   end interface

   procedure(no_args),  pointer :: problem_customize_solution => NULL()
   procedure(no_args),  pointer :: problem_grace_passed       => NULL()
   procedure(no_args),  pointer :: finalize_problem           => NULL()
   procedure(no_args),  pointer :: cleanup_problem            => NULL()
   procedure(tab_args), pointer :: custom_emf_bnd             => NULL()
   procedure(indx_args),pointer :: at_user_settings           => NULL()

contains

   subroutine array3d_init(this,nx,ny,nz)
      use constants, only: big_float
      implicit none
      class(array3d), intent(inout) :: this
      integer, intent(in)           :: nx,ny,nz

      if (.not.associated(this%arr)) allocate(this%arr(nx,ny,nz))
      this%arr = big_float
!     if (.not.associated(this%arr)) this%arr = reshape( [(big_float,i=1,nx*ny*nz)], [nx,ny,nz] )   ! lhs realloc
      return
   end subroutine array3d_init

   subroutine array3d_clean(this)
      implicit none
      class(array3d), intent(inout) :: this

      if (associated(this%arr)) deallocate(this%arr)
      return
   end subroutine array3d_clean

   logical function array3d_check_if_dirty(this)
      use constants, only: big_float
      implicit none
      class(array3d), intent(inout) :: this                  !! \todo i want to become polymorphic class(*) :/

      array3d_check_if_dirty = any( this%arr >= big_float )

   end function array3d_check_if_dirty

   subroutine array3d_associate(this,other)
      implicit none
      class(array3d), intent(inout) :: this
      real, allocatable, dimension(:,:,:), target :: other

      if (.not.associated(this%arr)) this%arr => other
      return
   end subroutine array3d_associate

   subroutine array4d_init(this,nn,nx,ny,nz)
      use constants, only: big_float
      implicit none
      class(array4d), intent(inout) :: this
      integer, intent(in)           :: nn,nx,ny,nz

      if (.not.associated(this%arr)) allocate(this%arr(nn,nx,ny,nz))
      this%arr = big_float
!     if (.not.associated(this%arr)) this%arr = reshape( [(big_float,i=1,nn*nx*ny*nz)], [nn,nx,ny,nz] )   ! lhs realloc
      return
   end subroutine array4d_init

   subroutine array4d_associate(this,other)
      implicit none
      class(array4d), intent(inout) :: this
      real, allocatable, dimension(:,:,:,:), target :: other

      if (.not.associated(this%arr)) this%arr => other
      return
   end subroutine array4d_associate

   subroutine array4d_clean(this)
      implicit none
      class(array4d), intent(inout) :: this                  !! Unlimited polymorphism at (1) not yet supported

      if (associated(this%arr)) deallocate(this%arr)
      return
   end subroutine array4d_clean

   logical function array4d_check_if_dirty(this)
      use constants, only: big_float
      implicit none
      class(array4d), intent(inout) :: this                  !! \todo i want to become polymorphic class(*) when I grow older

      array4d_check_if_dirty = any( this%arr >= big_float )

   end function array4d_check_if_dirty

   function array3d_get_sweep(this,ndim,i1,i2) result(p1d)
      use constants, only: xdim, ydim, zdim
      implicit none
      class(array3d), intent(inout) :: this
      real, dimension(:),  pointer  :: p1d
      integer, intent(in)           :: ndim, i1, i2

      if (.not.associated(this%arr)) then
         p1d => null()
      else
         select case (ndim)
            case (xdim)
               p1d => this%arr(:,i1,i2)
            case (ydim)
               p1d => this%arr(i2,:,i1)
            case (zdim)
               p1d => this%arr(i1,i2,:)
         end select
      endif
   end function array3d_get_sweep

   function array4d_get_sweep_one_var(this,ndim,nn,i1,i2) result(p1d)
      use constants, only: xdim, ydim, zdim
      implicit none
      class(array4d), intent(inout)      :: this
      integer, intent(in)                :: ndim, nn, i1, i2
      real, dimension(:),  pointer       :: p1d

      select case (ndim)
         case (xdim)
            p1d => this%arr(nn,:,i1,i2)
         case (ydim)
            p1d => this%arr(nn,i2,:,i1)
         case (zdim)
            p1d => this%arr(nn,i1,i2,:)
      end select
   end function array4d_get_sweep_one_var

   function array4d_get_sweep(this,ndim,i1,i2) result(p1d)
      use constants, only: xdim, ydim, zdim
      implicit none
      class(array4d), intent(inout)      :: this
      real, dimension(:,:),  pointer     :: p1d
      integer, intent(in)                :: ndim, i1, i2

      select case (ndim)
         case (xdim)
            p1d => this%arr(:,:,i1,i2)
         case (ydim)
            p1d => this%arr(:,i2,:,i1)
         case (zdim)
            p1d => this%arr(:,i1,i2,:)
      end select
   end function array4d_get_sweep

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
