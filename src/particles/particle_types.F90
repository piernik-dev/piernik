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

!>  \brief Types for selfgravitating particles

module particle_types
! pulled by GRAV
   use constants, only: ndims, dsetnamelen
  
   implicit none

   private
   public :: particle, particle_set, particle_data

   !>
   !! \brief simple particle: just mass and position
   !!
   !! \todo Extend it a bit
   !<

   type :: particle_data
     integer                :: pid        !< particle ID
     real                   :: mass       !< mass of the particle
     real, dimension(ndims) :: pos        !< physical position
     real, dimension(ndims) :: vel        !< particle velocity
#ifdef NBODY
     real, dimension(ndims) :: acc        !< acceleration of the particle
     real                   :: energy     !< total energy of particle
     logical                :: in, phy, out !< Flags to locate particle in the inner part of the domain or the outer part
#endif /* NBODY */
     logical                :: outside    !< this flag is true if the particle is outside the domain
   contains
     procedure :: is_outside              !< compute the outside flag
  end type particle_data

   type :: particle
      type(particle_data), pointer :: pdata !< list of particle data
      type(particle), pointer :: prv, nxt !< pointers to previous and next particle
   end type particle

   !> \brief A list of particles and some associated methods
   type :: particle_set
      type(particle), pointer :: first
      type(particle), pointer :: last
      integer(kind=4) :: cnt                  !< number of chain links
      character(len=dsetnamelen) :: label     !< name of the list for diagnostic and identification purposes
      !type(particle), allocatable, dimension(:) :: p !< the list of particles
      !type(particle), pointer :: p !< the list of particles
   contains
      procedure :: init        !< initialize the list
      !procedure :: print       !< print the list
      procedure :: cleanup     !< delete the list
      procedure :: remove      !< remove a particle
      !procedure :: merge_parts !< merge two particles
      procedure :: add_part_list   !< add a particle
      !procedure :: particle_with_id_exists  !< Check if particle no. "i" exists
      !generic, public :: exists => particle_with_id_exists
      generic, public :: add => add_part_list
   end type particle_set

contains

!> \brief compute the outside flag

   subroutine is_outside(this)

      use constants, only: LO, HI
      use domain,    only: dom

      implicit none

      class(particle_data), intent(inout) :: this     !< an object invoking the type-bound procedure

      this%outside = any(dom%has_dir(:) .and. (this%pos(:) < dom%edge(:, LO) .or. this%pos(:) >= dom%edge(:, HI)))
      ! Inequalities above must match the rounding function used in map routine (floor() includes bottom edge, but excludes top edge)

   end subroutine is_outside

!> \brief initialize the list with 0 elements

   subroutine init(this)

      use dataio_pub, only: die

      implicit none

      class(particle_set), intent(inout) :: this     !< an object invoking the type-bound procedure
      !character(len=*), intent(in)    :: label !< name of the list

      this%first => null()
      this%last  => null()
      this%cnt   =  0
      !this%label =  trim(label)
      !if (allocated(this%pdata)) call die("[particle_types:init] already initialized")
      !allocate(this%pdata(0))

   end subroutine init

!> \brief print the list

!   subroutine print(this)

!      use dataio_pub, only: msg, printinfo
!      use mpisetup,   only: slave

!      implicit none

!      class(particle_set), intent(inout) :: this     !< an object invoking the type-bound procedure

!      integer :: i

      !> \ todo communicate particles that aren't known to the master
!      if (slave) return
!
!      if (size(this%p) <= 0) return
!
!      call printinfo("[particle_types:print] Known particles:")
!      write(msg, '(a,a12,2(a,a36),2a)')" #number   : ","mass"," [ ","position"," ] [ ","velocity"," ]  is_outside"
!      call printinfo(msg)
!      do i = lbound(this%p, dim=1), ubound(this%p, dim=1)
!         write(msg, '(a,i7,a,g12.3,2(a,3g12.3),a,l2)')" # ",i," : ",this%p(i)%mass," [ ",this%p(i)%pos," ] [ ",this%p(i)%vel," ] ",this%p(i)%outside
!         call printinfo(msg)
!      enddo

!   end subroutine print

!> \brief delete the list

   subroutine cleanup(this)

      implicit none

      class(particle_set), intent(inout) :: this     !< an object invoking the type-bound procedure

      !Does nothing for now... is it necessary for a list?
      !if (allocated(this%p)) deallocate(this%p)

    end subroutine cleanup

!> \brief Add a particle to the list


    subroutine add_part_list(this, pid, mass, pos, vel, acc, energy, in, phy, out)

      use constants,  only: I_ONE
      use dataio_pub, only: die
      
      implicit none

      class(particle_set), intent(inout) :: this     !< an object invoking the type-bound procedure
      type(particle_data), pointer  :: part     !< new particle
      type(particle),      pointer       :: new
      integer,                intent(in) :: pid
      real, dimension(ndims), intent(in) :: pos, vel
#ifdef NBODY
      real, dimension(ndims), intent(in) :: acc
      real,                   intent(in) :: energy
#endif /* NBODY */
      real,                   intent(in) :: mass
      logical                            :: in,phy,out

      allocate(new)
      allocate(part)
      !if (present(part)) then
     ! if (associated(part)) then
            new%pdata => part
            new%pdata%pid = pid
            new%pdata%mass = mass
            new%pdata%pos = pos
            new%pdata%vel = vel
            new%pdata%acc = acc
            new%pdata%energy = energy
            new%pdata%in = in
            new%pdata%phy = phy
            new%pdata%out = out
            new%pdata%outside = .false.
      !   endif
      !else
      !   allocate(new%pdata)
      !endif
      new%nxt => null()
      
      if (.not. associated(this%first)) then ! the list was empty
         if (associated(this%last)) call die("[cg_list:add_new_part] last without first")
         this%first => new
         new%prv => null()
      else
         if (.not. associated(this%last)) call die("[cg_list:add_new_part] first without last")
         this%last%nxt => new
         new%prv => this%last
      endif
      
      this%last => new
      this%cnt = this%cnt + I_ONE
      
    end subroutine add_part_list
      

!> \brief Remove a particle number id from the list

   subroutine remove(this, pdata)

      use constants,  only: I_ONE
      use dataio_pub, only: die
      
      implicit none

      class(particle_set), intent(inout) :: this !< an object invoking the type-bound procedure
      type(particle), pointer, intent(out) :: pdata
      !integer,             intent(in)    :: id   !< position in the array of particles
      
      if (.not. associated(pdata)) call die("[particle removal] tried to remove null() element")
      if (.not. associated(this%first)) call die("[particle removal] this%cnt <=0 .and. associated(this%first)")

      if (associated(this%first, pdata)) this%first => this%first%nxt
      if (associated(this%last,  pdata)) this%last  => this%last%prv
      if (associated(pdata%prv)) pdata%prv%nxt => pdata%nxt
      if (associated(pdata%nxt)) pdata%nxt%prv => pdata%prv
      deallocate(pdata)
      this%cnt = this%cnt - I_ONE
      !if (.not. this%exists(id)) then
      !   write(msg, '("[particle_types:remove] Particle no Id = ",I6," does not exist")') id
      !   call die(msg)
      !endif

      !this%p = [this%p(:id-1), this%p(id+1:)]   ! LHS-realloc, please note  that if id == lbound(this%p, 1)
                                                ! or id == ubound(this%p, 1) it does not cause out of bound
                                                ! access in p
   end subroutine remove

!>
!! \brief Merge two particles
!!
!! \todo consider implementation as overloading of the (+) operator
!<

!   subroutine merge_parts(this, id1, id2)

!      implicit none

!      class(particle_set), intent(inout) :: this !< an object invoking the type-bound procedure
!      integer,             intent(in)    :: id1  !< position of the first particle in the array of particles (particle to be replaced by the merger)
!      integer,             intent(in)    :: id2  !< position of the second partilce in the array of particles (particle to be removed)
!
!      type(particle) :: merger
!
!      merger%mass = this%p(id1)%mass + this%p(id2)%mass
!      merger%pos  = (this%p(id1)%mass*this%p(id1)%pos + this%p(id2)%mass*this%p(id2)%pos) / merger%mass ! CoM
!
!      this%p(id1) = merger
!      call this%p(id1)%is_outside
!      call this%remove(id2)

!   end subroutine merge_parts

!   function particle_with_id_exists(this, id) result (tf)

!      implicit none

!      class(particle_set), intent(inout) :: this
!      integer,             intent(in)    :: id

!      logical :: tf

!      tf = allocated(this%p)
!      if (tf) tf = (id >= lbound(this%p, dim=1)) .and. (id <= ubound(this%p, dim=1))

!   end function particle_with_id_exists

end module particle_types
