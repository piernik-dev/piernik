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

!> \brief A module with an array of pointers to local grid containers

module gcpa

   use grid_cont, only: grid_container

   implicit none

   private
   public :: gcpa_T

   !> a pointer to a grid container
   type :: gcp
      type(grid_container), pointer :: p
   end type gcp

   !> an array of pointers to local grid containers
   type :: gcpa_T
      type(gcp), dimension(:), allocatable :: l_pse ! auxiliary array used to convert entries in cg_list_balance_T%dot%gse into pointers to grid containers for local exchanges
   contains
      procedure :: init
      procedure :: cleanup
   end type gcpa_T

contains

!> \brief Set up an array to be able to convert from local grid_id to pointers to cg

   subroutine init(this, curl)

      use cg_list,         only: cg_list_element
      use cg_list_balance, only: cg_list_balance_T
      use dataio_pub,      only: die
      use mpisetup,        only: proc

      implicit none

      class(gcpa_T),            intent(inout) :: this !< object invoking type bound procedure
      class(cg_list_balance_T), intent(inout) :: curl !< current level

      type(cg_list_element), pointer :: cgl
      integer                        :: b

      allocate(this%l_pse(lbound(curl%dot%gse(proc)%c(:), dim=1):ubound(curl%dot%gse(proc)%c(:), dim=1)))
      ! OPT: the curl%dot%gse is sorted, so the setting of this%l_pse can be done in a bit faster, less safe way. Or do it fast first, then try the safe way to fill up, what is missing, if anything
      do b = lbound(this%l_pse, dim=1), ubound(this%l_pse, dim=1)
         this%l_pse(b)%p => null()
         cgl => curl%first
         do while (associated(cgl))
            if (all(cgl%cg%my_se == curl%dot%gse(proc)%c(b)%se)) then
               this%l_pse(b)%p => cgl%cg
               exit
            endif
            cgl => cgl%nxt
         enddo
         if (.not. associated(this%l_pse(b)%p)) call die("[gcpa:init] this%l_pse pointer not set")
      enddo

   end subroutine init

!> \brief deallocate

   subroutine cleanup(this)

      implicit none

      class(gcpa_T), intent(inout) :: this !< object invoking type bound procedure

      if (allocated(this%l_pse)) deallocate(this%l_pse)

   end subroutine cleanup

end module gcpa
