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
!! \brief Collect names of nodes. Provide structures that allow for quick translation between rank and node name.
!!
!! To be used in load balancing and possibly in MPI-3 (shared memory parallelism).
!!
!! To allow usage in IO routines do not use entities from dataio_pub.
!<

module procnames

   use MPIF, only: MPI_MAX_PROCESSOR_NAME

   implicit none

   private
   public :: pnames

   ! all MPI rank associated with particular node name
   type nodeproc_t
      character(len=MPI_MAX_PROCESSOR_NAME) :: nodename   !< local $HOSTNAME
      integer(kind=4), allocatable, dimension(:) :: proc  !< list of MPI ranks that belong to this%nodename
      real :: wtime                                       !< estimated average execution time per cg of group of local MPI processes
   end type nodeproc_t

   ! all connections between MPI ranks and nodes
   type procnamelist_t
      character(len=MPI_MAX_PROCESSOR_NAME), allocatable, dimension(:) :: procnames  !< node names associated with MPI ranks
      real, allocatable, dimension(:) :: wtime                                       !< estimated execution time per cg of MPI ranks
      logical, allocatable, dimension(:) :: exclude                                  !< When .true. then exclude given thread from computations
      type(nodeproc_t), allocatable, dimension(:) :: proc_on_node                    !< array of nodes and MPI ranks
      integer, allocatable, dimension(:) :: hostindex                                !< index in proc_on_node for each MPI rank
      integer(kind=4) :: maxnamelen                                                  !< length of longest hostname
      logical :: speed_avail                                                         !< .true. after host speeds were calculated at least once
   contains
      procedure :: init                !< Initialize the pnames structure
      procedure :: cleanup             !< Clean up the pnames structure
      procedure :: calc_hostspeed      !< Compute this%proc_on_node(:)%wtime from this%wtime
      procedure :: mark_for_exclusion  !< Mark underperforming processes for exclusion
      procedure :: enable_all          !< Unmark any exclusions
   end type procnamelist_t

   type(procnamelist_t) :: pnames

contains

!< \brief Initialize the pnames structure

   subroutine init(this)

      use constants, only: I_ZERO
      use mpisetup,  only: FIRST, LAST, err_mpi
      use MPIF,      only: MPI_COMM_WORLD, MPI_CHARACTER, MPI_Get_processor_name
      use MPIFUN,    only: MPI_Allgather

      implicit none

      class(procnamelist_t), intent(inout) :: this !< an object invoking the type-bound procedure

      character(len=MPI_MAX_PROCESSOR_NAME) :: myname
      integer(kind=4) :: mynamelen
      character(len=MPI_MAX_PROCESSOR_NAME), allocatable, dimension(:) :: nodenames  !< aux array for unique node names

      this%speed_avail = .false.

      allocate(this%procnames(FIRST:LAST), &
           &   this%wtime    (FIRST:LAST), &
           &   this%exclude  (FIRST:LAST), &
           &   this%hostindex(FIRST:LAST))

      call this%enable_all
      this%maxnamelen = I_ZERO

      call MPI_Get_processor_name(myname, mynamelen, err_mpi)
      call MPI_Allgather(myname,         MPI_MAX_PROCESSOR_NAME, MPI_CHARACTER, &
           &             this%procnames, MPI_MAX_PROCESSOR_NAME, MPI_CHARACTER, &
           &             MPI_COMM_WORLD, err_mpi)

      call find_unique
      call fill_proc_on_node
      deallocate(nodenames)

      this%wtime = 0.

   contains

      !> \brief Find unique node names

      subroutine find_unique

         use constants,  only: I_ONE

         implicit none

         integer :: i, j
         logical :: found

         allocate(nodenames(I_ONE))
         nodenames(I_ONE) = this%procnames(FIRST)
         do i = lbound(this%procnames, 1), ubound(this%procnames, 1)
            found = .false.
            do j = lbound(nodenames, 1), ubound(nodenames, 1)
               found = found .or. (this%procnames(i) == nodenames(j))  ! longer loop perhaps would benefit from use of exit statement
            enddo
            if (.not. found) nodenames = [ nodenames, this%procnames(i) ]  ! lhs reallocation
         enddo

         do j = lbound(nodenames, 1), ubound(nodenames, 1)
            this%maxnamelen = max(this%maxnamelen, len_trim(nodenames(j), kind=4))
         enddo

      end subroutine find_unique

      !> Connect unique node names with MPI ranks

      subroutine fill_proc_on_node

         use constants, only: INVALID

         implicit none

         integer(kind=4) :: i, j

         allocate(this%proc_on_node(size(nodenames)))
         do i = lbound(nodenames, 1, kind=4), ubound(nodenames, 1, kind=4)
            this%proc_on_node(i)%nodename = nodenames(i)
            allocate(this%proc_on_node(i)%proc(0))
         enddo

         do j = lbound(this%proc_on_node, 1, kind=4), ubound(this%proc_on_node, 1, kind=4)
            do i = lbound(this%procnames, 1, kind=4), ubound(this%procnames, 1, kind=4)
               if (this%procnames(i) == this%proc_on_node(j)%nodename) &
                    & this%proc_on_node(j)%proc = [ this%proc_on_node(j)%proc, i ]  ! lhs reallocation
            enddo
         enddo

         ! Set up this%hostindex to be able to quickly refer to node properties knowing own MPI rank (mpisetup::proc)
         this%hostindex = INVALID
         do j = lbound(this%proc_on_node, 1, kind=4), ubound(this%proc_on_node, 1, kind=4)
            do i = lbound(this%proc_on_node(j)%proc, 1, kind=4), ubound(this%proc_on_node(j)%proc, 1, kind=4)
               this%hostindex(this%proc_on_node(j)%proc(i)) = j
            enddo
         enddo

      end subroutine fill_proc_on_node

   end subroutine init

!< \brief Clean up the pnames structure

   subroutine cleanup(this)

      implicit none

      class(procnamelist_t), intent(inout) :: this !< an object invoking the type-bound procedure

      integer :: i

      if (.not. allocated(this%procnames)) return  !< in case it wasn't ever initialized

      do i = lbound(this%proc_on_node, 1), ubound(this%proc_on_node, 1)
         deallocate(this%proc_on_node(i)%proc)
      enddo
      deallocate(this%proc_on_node)
      deallocate(this%procnames)
      deallocate(this%wtime)
      deallocate(this%exclude)
      deallocate(this%hostindex)

   end subroutine cleanup

!< \brief Compute this%proc_on_node(:)%wtime from this%wtime

   subroutine calc_hostspeed(this)

      implicit none

      class(procnamelist_t), intent(inout) :: this    !< an object invoking the type-bound procedure

      integer :: host
      real :: avg

      do host = lbound(this%proc_on_node, 1), ubound(this%proc_on_node, 1)
         associate (h => this%proc_on_node(host))
            avg = 0.  ! Don't average on unoccupied/excluded threads
            if (count(this%wtime(h%proc(:)) > 0.) > 0) avg = sum(this%wtime(h%proc(:))) / count(this%wtime(h%proc(:)) > 0.)
            h%wtime = avg
         end associate
      enddo

      this%speed_avail = .true.

   end subroutine calc_hostspeed

!< \brief Mark underperforming processes for exclusion

   subroutine mark_for_exclusion(this, threshold)

      implicit none

      class(procnamelist_t), intent(inout) :: this       !< an object invoking the type-bound procedure
      real,                  intent(in)    :: threshold  !< mark for exclusion when a process is that much slower than average

      real, parameter :: fast_enough = 1.2  ! count slightly slower threads in the average but reject marauders
      real :: avg, fast_avg

      if (count(.not. this%exclude .and. this%wtime(:) > 0.) <= 0) return  ! this may occur right after restart

      ! average MHD cost per cg on active threads
      avg = sum(this%wtime(:), mask = .not. this%exclude .and. this%wtime(:) > 0.) / &
           &                    count(.not. this%exclude .and. this%wtime(:) > 0.)

      ! average MHD cost per cg on active threads that aren't lagging too much behind average
      fast_avg = sum(this%wtime(:), mask = (.not. this%exclude .and. this%wtime(:) > 0. .and. this%wtime(:) <= fast_enough * avg)) / &
           &                          count(.not. this%exclude .and. this%wtime(:) > 0. .and. this%wtime(:) <= fast_enough * avg)

      this%exclude = this%exclude .or. this%wtime(:) > fast_avg * threshold

   end subroutine mark_for_exclusion

!< \brief Unmark any exclusions

   subroutine enable_all(this)

      implicit none

      class(procnamelist_t), intent(inout) :: this !< an object invoking the type-bound procedure

      this%exclude = .false.

   end subroutine enable_all

end module procnames
