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

   use MPIF,  only: MPI_MAX_PROCESSOR_NAME
   use types, only: ema_t

   implicit none

   private
   public :: pnames

   ! all MPI rank associated with particular node name
   type nodeproc_t
      character(len=MPI_MAX_PROCESSOR_NAME) :: nodename   !< local $HOSTNAME
      integer(kind=4), allocatable, dimension(:) :: proc  !< list of MPI ranks that belong to this%nodename
      type(ema_t) :: speed                                !< estimated average speed of local MPI processes
   end type nodeproc_t

   ! all connections between MPI ranks and nodes
   type procnamelist_t
      character(len=MPI_MAX_PROCESSOR_NAME), allocatable, dimension(:) :: procnames  !< node names associated with MPI ranks
      type(ema_t), allocatable, dimension(:) :: speed                                !< estimated speed of MPI ranks
      type(nodeproc_t), allocatable, dimension(:) :: proc_on_node                    !< array of nodes and MPI ranks
      integer(kind=4) :: maxnamelen
   contains
      procedure :: init            !< Initialize the pnames structure
      procedure :: cleanup         !< Clean up the pnames structure
      procedure :: calc_hostspeed  !< Compute this%proc_on_node(:)%speed from this%speed
   end type procnamelist_t

   type(procnamelist_t) :: pnames

contains

!< \brief Initialize the pnames structure

   subroutine init(this)

      use constants, only: I_ZERO
      use mpisetup,  only: FIRST, LAST, err_mpi
      use MPIF,      only: MPI_COMM_WORLD, MPI_CHARACTER, MPI_Get_processor_name, MPI_Allgather

      implicit none

      class(procnamelist_t), intent(inout) :: this

      character(len=MPI_MAX_PROCESSOR_NAME) :: myname
      integer(kind=4) :: mynamelen
      character(len=MPI_MAX_PROCESSOR_NAME), allocatable, dimension(:) :: nodenames  !< aux array for unique node names

      allocate(this%procnames(FIRST:LAST))
      allocate(this%speed(FIRST:LAST))
      this%maxnamelen = I_ZERO

      call MPI_Get_processor_name(myname, mynamelen, err_mpi)
      call MPI_Allgather(myname,         MPI_MAX_PROCESSOR_NAME, MPI_CHARACTER, &
           &             this%procnames, MPI_MAX_PROCESSOR_NAME, MPI_CHARACTER, &
           &             MPI_COMM_WORLD, err_mpi)

      call find_unique
      call fill_proc_on_node
      deallocate(nodenames)

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

      end subroutine fill_proc_on_node

   end subroutine init

!< \brief Clean up the pnames structure

   subroutine cleanup(this)

      implicit none

      class(procnamelist_t), intent(inout) :: this

      integer :: i

      if (.not. allocated(this%procnames)) return  !< in case it wasn't ever initialized

      do i = lbound(this%proc_on_node, 1), ubound(this%proc_on_node, 1)
         deallocate(this%proc_on_node(i)%proc)
      enddo
      deallocate(this%proc_on_node)
      deallocate(this%procnames)
      deallocate(this%speed)

   end subroutine cleanup

!< \brief Compute this%proc_on_node(:)%speed from this%speed

   subroutine calc_hostspeed(this, factor)

      implicit none

      class(procnamelist_t), intent(inout) :: this
      real, optional,        intent(in)    :: factor

      integer :: host
      real :: avg

      do host = lbound(this%proc_on_node, 1), ubound(this%proc_on_node, 1)
         avg = sum(this%speed(this%proc_on_node(host)%proc(:))%avg) / size(this%proc_on_node(host)%proc(:))
         if (present(factor)) then
            call this%proc_on_node(host)%speed%add(avg, factor)
         else
            call this%proc_on_node(host)%speed%add(avg)
         endif
      enddo

   end subroutine calc_hostspeed

end module procnames
