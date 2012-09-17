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

!> \brief Provide support routines for initialization and cleanup of data stored in grid containers by optional modules (like multigrid)

module grid_container_ext

   use constants, only: dsetnamelen

   implicit none

   private
   public :: cg_ext, cg_extptrs

   interface

      subroutine cg_ext(cg)

         use grid_cont, only: grid_container

         implicit none

         type(grid_container), pointer,  intent(inout) :: cg

      end subroutine cg_ext

   end interface

   type :: ext_ptr
      procedure(cg_ext), nopass, pointer :: init    !< initialization of extended features
      procedure(cg_ext), nopass, pointer :: cleanup !< cleanup of extended features
      character(len=dsetnamelen) :: label           !< name of the extension
   end type ext_ptr

   type :: ext_ptr_array
      type(ext_ptr), dimension(:), allocatable :: ext
   contains
      procedure :: epa_init    !< Allocate zero elements so we can safely count from lbound to ubound
      procedure :: epa_cleanup !< Deallocate
      generic, public :: extend => extend_independent, extend_dependent
      procedure :: extend_independent !< Add a new entry that does not depend on another entries
      procedure :: extend_dependent   !< Add a new entry that must be called after an existing one
   end type ext_ptr_array

   type(ext_ptr_array) :: cg_extptrs

contains

!> \brief Allocate zero elements so we can safely count from lbound to ubound

   subroutine epa_init(this)

      use dataio_pub, only: die

      implicit none

      class(ext_ptr_array), intent(inout) :: this

      if (allocated(this%ext)) call die("[grid_container_ext:epa_init] already allocated")
      allocate(this%ext(0))

   end subroutine epa_init

!> \bried Deallocate

   subroutine epa_cleanup(this)

      implicit none

      class(ext_ptr_array), intent(inout) :: this

      deallocate(this%ext)

   end subroutine epa_cleanup

!> \brief Add a new entry that does not depend on another entries

   subroutine extend_independent(this, init, cleanup, label)

      implicit none

      class(ext_ptr_array),       intent(inout) :: this
      procedure(cg_ext), pointer, intent(in)    :: init
      procedure(cg_ext), pointer, intent(in)    :: cleanup
      character(len=*),           intent(in)    :: label

      call this%extend_dependent(init, cleanup, label, "")

   end subroutine extend_independent

!> Add a new entry that must be called after an existing one

   subroutine extend_dependent(this, init, cleanup, label, after_label)

      use dataio_pub, only: die

      implicit none

      class(ext_ptr_array),       intent(inout) :: this
      procedure(cg_ext), pointer, intent(in)    :: init
      procedure(cg_ext), pointer, intent(in)    :: cleanup
      character(len=*),           intent(in)    :: label
      character(len=*),           intent(in)    :: after_label

      type(ext_ptr), dimension(:), allocatable :: new_ext
      integer :: i

      if (label == after_label) call die("[grid_container_ext:extend_dependent] label == after_label")

      do i = lbound(this%ext, dim=1), ubound(this%ext, dim=1)
         if (this%ext(i)%label == label) call die("[grid_container_ext:extend_dependent] label already exists")
      enddo

      if (len(after_label) > 0) then
         if (size(this%ext)<1) call die("[grid_container_ext:extend_dependent] no elements in the list")
         do i = lbound(this%ext, dim=1), ubound(this%ext, dim=1)
            if (this%ext(i)%label == after_label) exit
         enddo
         if (this%ext(i)%label /= after_label) call die("[grid_container_ext:extend_dependent] after_label not found")
      else
         i = ubound(this%ext, dim=1)
      endif

      allocate(new_ext(size(this%ext)+1))
      new_ext(:i) = this%ext(:i)
      new_ext(i+1)%init    => init
      new_ext(i+1)%cleanup => cleanup
      new_ext(i+1)%label   =  label
      if (i<ubound(this%ext, dim=1)) new_ext(i+2:ubound(new_ext, dim=1)) = this%ext(i+1:)
      call move_alloc(from=new_ext, to=this%ext)

   end subroutine extend_dependent

end module grid_container_ext
