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
#include "macros.h"

!>  \brief Types for selfgravitating particles

module particle_types
! pulled by GRAV
   use constants, only: ndims

   implicit none

   private
   public :: pset

   !>
   !! \brief simple particle: just mass and position
   !!
   !! \todo Extend it a bit
   !<
   type :: particle
      real                   :: mass       !< mass of the particle
      real, dimension(ndims) :: position   !< physical position
      real, dimension(ndims) :: vel        !< particle velocity
   end type particle

   !> \brief A list of particles and some associated methods

   type :: particle_set
      type(particle), allocatable, dimension(:) :: p !< the list of particles
   contains
      procedure :: init   !< initialize the list
      procedure :: remove !< remove a particle
      procedure :: merge  !< merge two particles
      procedure :: map    !< project particles onto grid
      procedure :: add_using_basic_types   !< add a particle
      procedure :: add_using_derived_type  !< add a particle
      generic, public :: add => add_using_basic_types, add_using_derived_type
   end type particle_set

   type(particle_set) :: pset !< default particle list

contains

!> \brief initialize the list with 0 elements

   subroutine init(this)

      use dataio_pub, only: die

      implicit none

      class(particle_set),    intent(inout) :: this     !< an object invoking the type-bound procedure

      if (allocated(this%p)) call die("[particle_types:init] already initialized")
      allocate(this%p(0))

   end subroutine init

!> \brief Add a particle to the list
!> \todo Consider using lhs-realloc

   subroutine add_using_derived_type(this, part)

      implicit none

      class(particle_set), intent(inout) :: this     !< an object invoking the type-bound procedure
      type(particle), intent(in) :: part             !< new particle

      type(particle), allocatable, dimension(:)  :: new_p
      allocate(new_p(size(this%p, dim=1) + 1))
      new_p(:size(this%p, dim=1)) = this%p(:)
      new_p(ubound(new_p, dim=1)) = part
      call move_alloc(from=new_p, to=this%p)

      ! this%p = [this%p, particle]  ! LHS-realloc

   end subroutine add_using_derived_type

!> \brief Add a particle to the list

   subroutine add_using_basic_types(this, mass, position, vel)

      use constants, only: ndims

      implicit none

      class(particle_set),    intent(inout) :: this     !< an object invoking the type-bound procedure
      real,                   intent(in)    :: mass     !< mass of the particle (negative values are allowed just in case someone wants to calculate electric potential)
      real, dimension(ndims), intent(in)    :: position !< physical position
      real, dimension(ndims), intent(in)    :: vel      !< particle velosity

      call this%add(particle(mass, position, vel))

   end subroutine add_using_basic_types

!> \brief Remove a partilce number id from the list
!> \todo Consider using lhs-realloc

   subroutine remove(this, id)

      use dataio_pub, only: msg, warn

      implicit none

      class(particle_set), intent(inout) :: this !< an object invoking the type-bound procedure
      integer,             intent(in)    :: id   !< position in the array of particles

      type(particle), allocatable, dimension(:)  :: new_p
      integer :: i

      if (id > size(this%p, dim=1)) then
         write(msg, '("[particle_set:remove] Id = ",I6," is greater than total number of particles: ",I6)') id, size(this%p, dim=1)
         call warn(msg)
         call warn("[particle_set:remove] I'll remove last particle instead")
      endif

      i = min(id, size(this%p, dim=1))  ! Sanitize id

      allocate(new_p(size(this%p, dim=1) - 1))
      new_p(:i-1) = this%p(:i-1)
      if (i < size(this%p, dim=1)) new_p(i:) = this%p(i+1:)
      call move_alloc(from=new_p, to=this%p)

      ! if (id >= size(this%p, dim=1) then
      !    ... copy warns ...
      !    this%p = this%p(:size(this%p, dim=1) - 1)  ! LHS-realloc
      ! else
      !    this%p = [this%p(:id-1), this%p(id+1:)]  ! LHS-realloc
      ! endif
   end subroutine remove

!>
!! \brief Merge two particles
!!
!! \todo consider implementation as overloading of the (+) operator
!<

   subroutine merge(this, id1, id2)

      implicit none

      class(particle_set), intent(inout) :: this !< an object invoking the type-bound procedure
      integer,             intent(in)    :: id1  !< position of the first particle in the array of particles (particle to be replaced by the merger)
      integer,             intent(in)    :: id2  !< position of the second partilce in the array of particles (particle to be removed)

      type(particle) :: merger

      merger%mass     = this%p(id1)%mass + this%p(id2)%mass
      merger%position = (this%p(id1)%mass*this%p(id1)%position + this%p(id2)%mass*this%p(id2)%position) / merger%mass ! CoM

      this%p(id1) = merger
      call this%remove(id2)

   end subroutine merge

   subroutine map(this, iv, factor)

      use cg_leaves, only: leaves
      use cg_list,   only: cg_list_element
      use constants, only: xdim, ydim, zdim, ndims, LO, HI

      implicit none

      class(particle_set), intent(inout) :: this   !< an object invoking the type-bound procedure
      integer,             intent(in)    :: iv     !< index in cg%q array, where we want the particles to be projected
      real,                intent(in)    :: factor !< typically fpiG

      type(cg_list_element), pointer :: cgl
      integer :: p
      integer(kind=8), dimension(ndims) :: ijkp

      cgl => leaves%first
      do while (associated(cgl))

         do p = lbound(this%p, dim=1), ubound(this%p, dim=1)
            ijkp(:) = floor((this%p(p)%position(:)-cgl%cg%fbnd(:, LO))/cgl%cg%dl(:)) + cgl%cg%ijkse(:, LO)
            if (all(ijkp >= cgl%cg%ijkse(:,LO)) .and. all(ijkp <= cgl%cg%ijkse(:,HI))) &
                 cgl%cg%q(iv)%arr(ijkp(xdim), ijkp(ydim), ijkp(zdim)) = cgl%cg%q(iv)%arr(ijkp(xdim), ijkp(ydim), ijkp(zdim)) + factor * this%p(p)%mass / cgl%cg%dvol
         enddo

         cgl => cgl%nxt
      enddo


   end subroutine map

end module particle_types
