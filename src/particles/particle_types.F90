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
   use constants, only: ndims

   implicit none

   public ! QA_WARN no secrets are kept here

   integer(kind=4), parameter :: npb = 2   !< number of cells between in and phy or between phy and out boundaries

   !> \brief enumerators for packung and unpacking the particles for teleportation to another cg
   enum, bind(C)
      enumerator :: P_ID=1, P_MASS, P_POS_X, P_POS_Y, P_POS_Z, P_VEL_X, P_VEL_Y, P_VEL_Z, P_ACC_X, P_ACC_Y, P_ACC_Z, P_ENER, P_TFORM, P_TDYN
   end enum
   integer(kind=4), parameter :: npf = P_TDYN  !< last enumerated element: the number of single particle fields

   !> \brief Particle type

   type :: particle_data
      integer(kind=4)        :: pid            !< particle ID
      real                   :: mass           !< mass of the particle
      real, dimension(ndims) :: pos            !< physical position
      real, dimension(ndims) :: vel            !< particle velocity
      real, dimension(ndims) :: acc            !< acceleration of the particle
      real                   :: energy         !< total energy of particle
      real                   :: tform          !< formation time of the particle
      real                   :: tdyn           !< dynamical time for SF
      logical                :: in, phy, out   !< Flags to locate particle in the inner part of the domain or the outer part
      logical                :: fin            !< this flag is true if the particle is located in a finest level cell
      logical                :: outside        !< this flag is true if the particle is outside the domain
   contains
      procedure :: is_outside                  !< compute the outside flag
   end type particle_data

   type :: particle
      type(particle_data), pointer :: pdata    !< list of particle data
      type(particle),      pointer :: prv, nxt !< pointers to previous and next particle
   end type particle

   !> \brief A list of particles and some associated methods
   type :: particle_set
      type(particle), pointer :: first
      type(particle), pointer :: last
      integer(kind=4)         :: cnt         !< number of chain links
   contains
      procedure :: init                      !< initialize the list
      procedure :: print                     !< print the list
      procedure :: count_phy                 !< count particles with phy flag set
      procedure :: cleanup                   !< delete the list
      procedure :: remove                    !< remove a particle
      !procedure :: merge_parts              !< merge two particles
      procedure :: add_part_list             !< add a particle (create from parameters)
      procedure :: add_part                  !< add a particle (copy from particle type)
      !procedure :: particle_with_id_exists  !< Check if particle no. "i" exists
      !generic, public :: exists => particle_with_id_exists
      generic, public :: add => add_part_list, add_part
      generic, public :: count => count_phy
   end type particle_set

contains

!> \brief compute the outside flag

   subroutine is_outside(this)  ! TO DO: CHECK if this is always the needed definition of outside (with boundary cells)

      use constants, only: LO, HI
      use domain,    only: dom

      implicit none

      class(particle_data), intent(inout) :: this     !< an object invoking the type-bound procedure

      this%outside = any(dom%has_dir(:) .and. (this%pos(:) < dom%edge(:, LO) .or. this%pos(:) >= dom%edge(:, HI)))
      ! Inequalities above must match the rounding function used in map routine (floor() includes bottom edge, but excludes top edge)

   end subroutine is_outside

!> \brief initialize the list with 0 elements

   subroutine init(this)

      use constants,  only: ndims
      use dataio_pub, only: die
      use domain,     only: dom

      implicit none

      class(particle_set), intent(inout) :: this     !< an object invoking the type-bound procedure
      !character(len=*), intent(in)    :: label !< name of the list

      this%first => null()
      this%last  => null()
      this%cnt   =  0

      if (dom%eff_dim /= ndims) call die("[particle_types:init] Only 3D is supported")
      ! Various functions related to particles aren't really ready for 2D grid. Perhaps it won't ever be needed.

   end subroutine init

!> \brief print the list

   subroutine print(this, remark)

      use dataio_pub, only: msg, printinfo, warn

      implicit none

      class(particle_set), intent(inout) :: this    !< an object invoking the type-bound procedure
      character(len=*),    intent(in)    :: remark  !< additional comment from the caller

      integer :: cnt
      type(particle), pointer :: pp

      if (this%cnt < 0) call warn("[particle_types:print] this%cnt < 0")
      if (associated(this%first) .neqv. associated(this%last)) call warn("[particle_types:print] associated(this%first) .neqv. associated(this%last)")
      if (this%cnt > 0 .and. .not. associated(this%first)) call warn("[particle_types:print] this%cnt > 0 .and. .not. associated(this%first)")
      if (this%cnt > 0 .and. .not. associated(this%last))  call warn("[particle_types:print] this%cnt > 0 .and. .not. associated(this%last)")
      if (this%cnt <= 0 .and. (associated(this%first) .or. associated(this%last))) call warn("[particle_types:print] <=0 .and. (associated(this%first) .or. associated(this%last))")

      cnt = this%cnt
      if (associated(this%first)) then
         pp => this%first
         do while (associated(pp))
            call prntline
            cnt = cnt - 1
            pp => pp%nxt
         enddo
      else if (associated(this%last)) then
         pp => this%last
         call warn(trim(remark) // "[particle_types:print] Going backward")
         do while (associated(pp))
            call prntline
            cnt = cnt - 1
            pp => pp%prv
         enddo
      else
         if (this%cnt > 0) then
            write(msg, '(a,i6,a)')trim(remark), this%cnt, " particles missing in the set"
            call printinfo(msg)
         endif
      endif

      if (cnt /= 0) call warn(trim(remark) // "[particle_types:print] particle counter messed up")

   contains

      subroutine prntline

         use mpisetup,  only: proc

         implicit none

         associate (p => pp%pdata)
            write(msg, '(a,i5,a,i8,a,g12.5,2(a,3g12.5),a,5l2)')trim(remark) // "@", proc, " #",p%pid," : ",p%mass," [ ",p%pos," ] [ ",p%vel," ] ",p%outside, p%in, p%out, p%phy, p%fin
         end associate
         call printinfo(msg)

      end subroutine prntline

   end subroutine print

!> \brief count particles with phy flag set

   integer(kind=4) function count_phy(this) result(n_part)

      use constants, only: I_ONE

      implicit none

      class(particle_set), intent(in) :: this  !< an object invoking the type-bound procedure

      type(particle), pointer :: pset

      n_part = 0
      pset => this%first
      do while (associated(pset))
         if (pset%pdata%phy) n_part = n_part + I_ONE
         pset => pset%nxt
      enddo

   end function count_phy

!> \brief delete the list

   subroutine cleanup(this)

      implicit none

      class(particle_set), intent(inout) :: this     !< an object invoking the type-bound procedure

      type(particle), pointer :: pp

      do while (associated(this%first))
         pp => this%last
         call this%remove(pp) ! cannot just pass this%last because it will change after un_link and wrong element will be deallocated
      enddo

   end subroutine cleanup

!> \brief Add a particle to the list (create from parameters)

   subroutine add_part_list(this, pid, mass, pos, vel, acc, energy, in, phy, out, fin, tform, tdyn)

      use constants,  only: I_ONE
      use dataio_pub, only: die

      implicit none

      class(particle_set), intent(inout) :: this     !< an object invoking the type-bound procedure
      type(particle_data), pointer       :: part     !< new particle
      type(particle),      pointer       :: new
      integer(kind=4),        intent(in) :: pid
      real, dimension(ndims), intent(in) :: pos, vel
      real, dimension(ndims), intent(in) :: acc
      real,                   intent(in) :: energy
      real,                   intent(in) :: mass, tform, tdyn
      logical,                intent(in) :: in, phy, out, fin

      allocate(new)
      allocate(part)
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
            new%pdata%fin = fin
            new%pdata%outside = .false.
            new%pdata%tform = tform
            new%pdata%tdyn = tdyn
            call new%pdata%is_outside()
      new%nxt => null()

      if (.not. associated(this%first)) then ! the list was empty
         if (associated(this%last)) call die("[particle_types:add_part_list] last without first")
         this%first => new
         new%prv => null()
      else
         if (.not. associated(this%last)) call die("[particle_types:add_part_list] first without last")
         this%last%nxt => new
         new%prv => this%last
      endif

      this%last => new
      this%cnt = this%cnt + I_ONE

   end subroutine add_part_list

!> \brief Add a particle (copy from particle type)

   subroutine add_part(this, pd)

      use constants,  only: I_ONE
      use dataio_pub, only: die

      implicit none

      class(particle_set),          intent(inout) :: this  !< an object invoking the type-bound procedure
      type(particle_data), pointer, intent(in)    :: pd

      type(particle), pointer :: new

      allocate(new)
      allocate(new%pdata)
      new%pdata = pd

      ! spaghetti warning: this is the same code as in add_part_list
      new%nxt => null()

      if (.not. associated(this%first)) then ! the list was empty
         if (associated(this%last)) call die("[particle_types:add_part] last without first")
         this%first => new
         new%prv => null()
      else
         if (.not. associated(this%last)) call die("[particle_types:add_part] first without last")
         this%last%nxt => new
         new%prv => this%last
      endif

      this%last => new
      this%cnt = this%cnt + I_ONE

   end subroutine add_part

!> \brief Remove a particle number id from the list

   subroutine remove(this, pset)

      use constants,  only: I_ONE
      use dataio_pub, only: die

      implicit none

      class(particle_set),     intent(inout) :: this !< an object invoking the type-bound procedure
      type(particle), pointer, intent(inout) :: pset

      if (.not. associated(pset)) call die("[particle_types:remove] tried to remove null() element")
      if (.not. associated(this%first)) call die("[particle_types:remove] .not. associated(this%first)")

      if (associated(this%first, pset)) this%first => this%first%nxt
      if (associated(this%last,  pset)) this%last  => this%last%prv
      if (associated(pset%prv)) pset%prv%nxt => pset%nxt
      if (associated(pset%nxt)) pset%nxt%prv => pset%prv
      deallocate(pset%pdata)
      deallocate(pset)
      this%cnt = this%cnt - I_ONE

      if (this%cnt < 0) call die("[particle_types:remove] this%cnt < 0")

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
!      integer,             intent(in)    :: id2  !< position of the second particle in the array of particles (particle to be removed)
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
