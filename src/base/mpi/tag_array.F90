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

!> \brief Type for array of tag arrays

module tag_array

   use sort_tags, only: sort_tags_t

   implicit none

   private
   public :: tag_arr

   type :: tag_arr  ! array of tag arrays
      type(sort_tags_t), allocatable, dimension(:) :: t  !< list of tags for each proc in current communication
   contains
      procedure          :: cleanup     !< find problems free the resources
      procedure          :: store_tag   !< store tags for inspection upon Waitall
      procedure, private :: are_unique  !< are the stored tags unique?
   end type tag_arr

contains

!> \brief clean up

   subroutine cleanup(this, label)

      implicit none

      class(tag_arr), target, intent(inout) :: this   !< an object invoking the type-bound procedure
      character(len=*),       intent(in)    :: label  !< identification

      integer :: i

      if (allocated(this%t)) then
         call this%are_unique(label)
         do i = lbound(this%t, 1), ubound(this%t, 1)
            if (allocated(this%t(i)%list)) deallocate(this%t(i)%list)
         enddo
         deallocate(this%t)
      endif

   end subroutine cleanup

!> brief Accumulate tags for inspection

   subroutine store_tag(this, tag, other_proc, n)

      use mpisetup,  only: FIRST, LAST
      use sort_tags, only: req_data

      implicit none

      class(tag_arr),  intent(inout) :: this        !< an object invoking the type-bound procedure
      integer(kind=4), intent(in)    :: tag         !< the tag
      integer(kind=4), intent(in)    :: other_proc  !< source or destination
      integer(kind=4), intent(in)    :: n           !< request number

      if (.not. allocated(this%t)) allocate(this%t(FIRST:LAST))
      call this%t(other_proc)%dump(req_data(tag, n))

   end subroutine store_tag

!> \brief Are the stored tags unique?

   subroutine are_unique(this, label)

      use constants,  only: V_DEBUG
      use dataio_pub, only: msg, printinfo, warn
      use mpisetup,   only: proc

      implicit none

      class(tag_arr),   intent(inout) :: this   !< an object invoking the type-bound procedure
      character(len=*), intent(in)    :: label  !< identification

      integer :: i, j, cnt

      cnt = 0
      if (allocated(this%t)) then
         do i = lbound(this%t, 1), ubound(this%t, 1)
            if (allocated(this%t(i)%list)) then
               if (this%t(i)%cnt > 1) then
                  call this%t(i)%sort
                  do j = lbound(this%t(i)%list, 1) + 1, this%t(i)%cnt
                     if (this%t(i)%list(j-1)%tag >= this%t(i)%list(j)%tag) then
                        write(msg, '(3(a,i0))')"[tag_array:are_unique] " // trim(label) // " @", proc, " : tag collision detected in communication with ", i, " tag=", this%t(i)%list(j-1)%tag
                        call printinfo(msg, V_DEBUG)
                        cnt = cnt + 1
                     endif
                  enddo
               endif
            endif
         enddo
      endif

      if (cnt /= 0) then
         write(msg, '(2(a,i0),a)')"[tag_array:are_unique] " // trim(label) // " @", proc, " : ", cnt, " tag collisions detected in communication"
         call warn(msg)
      endif

   end subroutine are_unique

end module tag_array
