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

!> \brief This module adds named arrays to the grid container

module grid_cont_na

   use grid_cont_base, only: grid_container_base_t
   use named_array,    only: named_array4d, named_array3d

   implicit none

   private
   public :: grid_container_na_t

   !> \brief This type adds named arrays and related routines to the grid container
   type, extends(grid_container_base_t), abstract :: grid_container_na_t

      ! Registered variables

      type(named_array3d), allocatable, dimension(:) :: q  !< 3D arrays such as gravitational potential pr user-defined quantities or gravitational potential
      type(named_array4d), allocatable, dimension(:) :: w  !< 4D arrays such as u, vector fields (b) or other vector/multi-scalar user-defined quantities

      ! handy shortcuts to some entries in q(:)
      real, dimension(:,:,:), pointer :: gpot    => null()  !< Array for sum of gravitational potential at t += dt
      real, dimension(:,:,:), pointer :: hgpot   => null()  !< Array for sum of gravitational potential at t += 0.5*dt
      real, dimension(:,:,:), pointer :: gp      => null()  !< Array for gravitational potential from external fields
      real, dimension(:,:,:), pointer :: sgp     => null()  !< Array for gravitational potential from multigrid or FFT solver
      real, dimension(:,:,:), pointer :: sgpm    => null()  !< Array for gravitational potential from multigrid or FFT solver at previous timestep saved by source_terms_grav.
      real, dimension(:,:,:), pointer :: cs_iso2 => null()  !< Array for sound speed for isothermal EOS (not associated for gamma EOS)
      real, dimension(:,:,:), pointer :: wa      => null()  !< Temporary array used for different purposes, usually has dimension (grid::nx, grid::ny, grid::nz)
#ifdef NBODY
      real, dimension(:,:,:), pointer :: prth    => null()  !< Array for histogram of particles
      real, dimension(:,:,:), pointer :: nbdn    => null()  !< Array of density from particles
      real, dimension(:,:,:), pointer :: gp1b    => null()  !< Array of gravitational potential from particles
#ifdef NBODY_GRIDDIRECT
      real, dimension(:,:,:), pointer :: nbgp    => null()  !< Array of gravitational potential from particles
#endif /* NBODY_GRIDDIRECT */
#endif /* NBODY */

      ! handy shortcuts to some entries in w(:)
      real, dimension(:,:,:,:), pointer :: u     => null()  !< Main array of all fluids' components
      real, dimension(:,:,:,:), pointer :: b     => null()  !< Main array of magnetic field's components

   contains

      procedure :: cleanup_na            !< Deallocate all internals
      procedure :: add_all_na            !< Register all known named arrays for this cg, set up shortcuts to the crucial fields
      procedure :: set_constant_b_field  !< set constant magnetic field on whole block

      ! These should be private procedures but we need them in cg_list_global:reg_var
      procedure :: add_na                !< Register a new 3D entry in current cg with given name.
      procedure :: add_na_4d             !< Register a new 4D entry in current cg with given name.

   end type grid_container_na_t

contains

!> \brief Routines that deallocates all internals of the grid container

   subroutine cleanup_na(this)

      implicit none

      class(grid_container_na_t), intent(inout) :: this !< object invoking type-bound procedure

      integer :: g

      if (allocated(this%q)) then
         do g = lbound(this%q(:), dim=1), ubound(this%q(:), dim=1)
            call this%q(g)%clean
         enddo
         deallocate(this%q)
      endif

      if (allocated(this%w)) then
         do g = lbound(this%w(:), dim=1), ubound(this%w(:), dim=1)
            call this%w(g)%clean
         enddo
         deallocate(this%w)
      endif

   end subroutine cleanup_na

!> \brief Register all known named arrays for this cg, set up shortcuts to the crucial fields

   subroutine add_all_na(this)

      use constants,        only: INVALID
      use memory_usage,     only: check_mem_usage
      use named_array_list, only: qna, wna
#ifdef ISO
      use constants,        only: cs_i2_n
      use fluids_pub,       only: cs2_max
#endif /* ISO */

      implicit none

      class(grid_container_na_t), intent(inout) :: this  !< object invoking type-bound procedure

      integer :: i

      if (allocated(qna%lst)) then
         do i = lbound(qna%lst(:), dim=1), ubound(qna%lst(:), dim=1)
            call this%add_na(qna%lst(i)%multigrid)
         enddo
      endif

      if (allocated(wna%lst)) then
         do i = lbound(wna%lst(:), dim=1), ubound(wna%lst(:), dim=1)
            call this%add_na_4d(wna%get_dim4(int(i, kind=4)))
         enddo
      endif

      call check_mem_usage

      ! shortcuts
      if (wna%fi > INVALID)  this%u  => this%w(wna%fi)%arr
      if (wna%bi > INVALID)  this%b  => this%w(wna%bi)%arr
      if (qna%wai > INVALID) this%wa => this%q(qna%wai)%arr

#ifdef ISO
      this%cs_iso2 => this%q(qna%ind(cs_i2_n))%arr
      if (associated(this%cs_iso2)) this%cs_iso2(:,:,:) = cs2_max   ! set cs2 with sane values on non-multigrid grid pieces
#endif /* ISO */

   end subroutine add_all_na

!>
!! \brief Register a new 3D entry in current cg with given name. Called from cg_list_glob::reg_var
!!
!! \warning This routine should not be called directly from user routines
!<
   subroutine add_na(this, multigrid)

      use constants,   only: base_level_id, LO, HI
      use named_array, only: named_array3d

      implicit none

      class(grid_container_na_t), intent(inout) :: this       !< object invoking type-bound procedure
      logical,                    intent(in)    :: multigrid  !< If .true. then cg%q(:)%arr and cg%w(:)%arr are allocated also below base level

      type(named_array3d), allocatable, dimension(:) :: tmp

      if (.not. allocated(this%q)) then
         allocate(this%q(1))
      else
         allocate(tmp(lbound(this%q(:),dim=1):ubound(this%q(:), dim=1) + 1))
         tmp(:ubound(this%q(:), dim=1)) = this%q(:)
         call move_alloc(from=tmp, to=this%q)
      endif

      if (multigrid .or. this%l%id >= base_level_id) call this%q(ubound(this%q(:), dim=1))%init(this%lhn(:, LO), this%lhn(:, HI))

   end subroutine add_na

!>
!! \brief Register a new 4D entry in current cg with given name. Called from cg_list_glob::reg_var
!!
!! \warning This routine should not be called directly from user routines
!<
   subroutine add_na_4d(this, n)

      use constants,   only: base_level_id, INT4, LO, HI
      use named_array, only: named_array4d

      implicit none

      class(grid_container_na_t), intent(inout) :: this  !< object invoking type-bound procedure
      integer(kind=4),            intent(in)    :: n     !< Length of the vector quantity to be stored (first dimension of the array)

      type(named_array4d), allocatable, dimension(:) :: tmp

      if (.not. allocated(this%w)) then
         allocate(this%w(1))
      else
         allocate(tmp(lbound(this%w(:),dim=1):ubound(this%w(:), dim=1) + 1))
         tmp(:ubound(this%w(:), dim=1)) = this%w(:)
         call move_alloc(from=tmp, to=this%w)
      endif

      if (this%l%id >= base_level_id) call this%w(ubound(this%w(:), dim=1))%init( [1_INT4, this%lhn(:, LO)], [n, this%lhn(:, HI)] ) !< \deprecated magic integer

   end subroutine add_na_4d

!< \brief set constant magnetic field on whole block

   subroutine set_constant_b_field(this, b)

      use constants, only: xdim, zdim

      implicit none

      class(grid_container_na_t), intent(inout) :: this !< object invoking type-bound procedure
      real, dimension(xdim:zdim), intent(in)    :: b    !< the value of the magnetic field vector in whole block

      integer :: d

      if (associated(this%b)) then
         do d = xdim, zdim
            this%b(d, this%is:this%ie, this%js:this%je, this%ks:this%ke) = b(d)
         enddo
      endif

   end subroutine set_constant_b_field

end module grid_cont_na
