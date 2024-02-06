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

!> \brief Module for providing event array for precise parallel profiling data

module ppp_events

   use constants, only: cbuff_len

   implicit none

   private
   public :: event, eventarray

   !> \brief a cheap single-event entry
   type :: event
      character(len=cbuff_len) :: label  !< label used to identify the event
      real(kind=8)             :: wtime  !< output from MPI_Wtime(); positive for start, negative for end
   end type event

   !> \brief a cheap array of events
   type :: eventarray
      type(event), dimension(:), allocatable :: ev_arr
   contains
      procedure :: init     !< allocate eventarray of given size
      procedure :: cleanup  !< deallocate eventarray
   end type eventarray

contains

!> \brief allocate eventarray of given size

   subroutine init(this, asize)

      use memory_usage, only: check_mem_usage
      use dataio_pub,   only: die

      implicit none

      class(eventarray), intent(inout) :: this   !< an object invoking the type-bound procedure
      integer,           intent(in)    :: asize  !< size of the event array

      if (allocated(this%ev_arr)) call die("[ppp_events:init] already allocated")
      allocate(this%ev_arr(asize))
      call check_mem_usage
!      this%ev_arr(:)%wtime = 0.

   end subroutine init

!> \brief deallocate eventarray

   subroutine cleanup(this)

      use dataio_pub, only: die

      implicit none

      class(eventarray), intent(inout) :: this   !< an object invoking the type-bound procedure

      if (.not. allocated(this%ev_arr)) call die("[ppp_events:cleanup] not allocated")
      deallocate(this%ev_arr)

   end subroutine cleanup

end module ppp_events
