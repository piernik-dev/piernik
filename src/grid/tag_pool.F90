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

!> \brief This module manages tag pool for cg_level_connected_T%vertical_b_prep

module tag_pool

   implicit none

   private
   public :: t_pool

   integer(kind=4), parameter :: POOL_SIZE = 2**16  !< number of entries in a single pool. This one gives 2**15 pools. Try to make it smaller if you ever run out of pools.
   integer(kind=4), parameter :: FREE = -2147483647 !< -2**31+1, special id, that means a given pool is free

   type :: tag_pool_T
      integer(kind=4), dimension(:), allocatable :: t_map
   contains
      procedure :: cleanup !< Deallocate everything
      procedure :: get     !< Give a pool of tags and mark it with given id
      procedure :: release !< Release all pools of tags with given id
   end type tag_pool_T

   type(tag_pool_T) :: t_pool

contains

!> \brief Deallocate the map of pools

   subroutine cleanup(this)

      implicit none

      class(tag_pool_T), intent(inout) :: this

      if (allocated(this%t_map)) deallocate(this%t_map)

   end subroutine cleanup

!> \brief Give a pool of tags and mark it with given id

   subroutine get(this, id, t_start, t_end)

      use constants,  only: INVALID, I_ONE
      use dataio_pub, only: die

      implicit none

      class(tag_pool_T), intent(inout) :: this
      integer(kind=4),   intent(in)    :: id
      integer(kind=4),   intent(out)   :: t_start
      integer(kind=4),   intent(out)   :: t_end

      integer(kind=4) :: i

      if (id == FREE) call die("[tag_pool:get] reserved id")

      if (.not. allocated(this%t_map)) then
         allocate(this%t_map(1))
         this%t_map = FREE
      endif

      t_start = INVALID
      t_end = INVALID
      if (.not. any(this%t_map == FREE)) this%t_map = [ this%t_map, FREE ] ! LHS realloc
      do i = ubound(this%t_map, dim=1, kind=4), lbound(this%t_map, dim=1, kind=4), -1
         if (this%t_map(i) == FREE) then
            t_start = (i-lbound(this%t_map, dim=1, kind=4)) * POOL_SIZE
            t_end = t_start + POOL_SIZE - I_ONE
            this%t_map(i) = id
            exit
         endif
      enddo
      if (t_start == INVALID .or. t_end == INVALID .or. t_start < 0 .or. t_end < 0) call die("[tag_pool:get] Cannot give valid pool of tags")

   end subroutine get

!> \brief Release all pools of tags with given id

   subroutine release(this, id)

      use dataio_pub, only: die

      implicit none

      class(tag_pool_T), intent(inout) :: this
      integer(kind=4),   intent(in)    :: id

      integer :: i

      if (id == FREE) call die("[tag_pool:release] reserved id")
      if (.not. allocated(this%t_map)) return

      do i = lbound(this%t_map, dim=1), ubound(this%t_map, dim=1)
         if (this%t_map(i) == id) this%t_map(i) = FREE
      enddo

   end subroutine release

end module tag_pool
