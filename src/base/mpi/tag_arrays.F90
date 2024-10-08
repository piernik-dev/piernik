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

!> \brief Type for handling array of tag arrays separately for send and receive operations.

module tag_arrays

   use constants, only: cbuff_len
   use tag_array, only: tag_arr

   implicit none

   private
   public :: tag_arrs

   type :: tag_arrs
      type(tag_arr) :: ts  !< list of "send" tags for each proc in current communication
      type(tag_arr) :: tr  !< list of "recv" tags for each proc in current communication
      character(len = cbuff_len) :: label  !< identificator
   contains
      procedure :: cleanup    !< find problems and free the resources
      procedure :: store_tag  !< store tags for inspection upon Waitall
      procedure :: set_label  !< set the label
   end type tag_arrs

contains

!> \brief clean up

   subroutine cleanup(this)

      implicit none

      class(tag_arrs), intent(inout), target :: this  !< an object invoking the type-bound procedure

      call this%ts%cleanup(trim(this%label) // "_send")
      call this%tr%cleanup(trim(this%label) // "_recv")

   end subroutine cleanup

!> brief Accumulate tags for inspection

   subroutine store_tag(this, tag, other_proc, n, recv)

      implicit none

      class(tag_arrs), intent(inout) :: this        !< an object invoking the type-bound procedure
      integer(kind=4), intent(in)    :: tag         !< the tag
      integer(kind=4), intent(in)    :: other_proc  !< source or destination
      integer(kind=4), intent(in)    :: n           !< request number
      logical,         intent(in)    :: recv        !< send or receive?

      if (recv) then
         call this%tr%store_tag(tag, other_proc, n)
      else
         call this%ts%store_tag(tag, other_proc, n)
      endif

   end subroutine store_tag

!> \brief Set the label

   subroutine set_label(this, label)

      use constants,  only: cbuff_len
      use dataio_pub, only: die

      implicit none

      class(tag_arrs),  intent(inout) :: this   !< an object invoking the type-bound procedure
      character(len=*), intent(in)    :: label  !< identification

      if (len_trim(this%label) > 0) then
         call die("[tag_arrays:set_label] label already set: '" // label // "' vs '" // this%label // "'")
      else
         this%label = label(1:min(cbuff_len, len_trim(label, kind=4)))
      endif

   end subroutine set_label

end module tag_arrays
