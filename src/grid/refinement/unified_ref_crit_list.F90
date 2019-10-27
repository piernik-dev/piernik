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
!    along with PIERNIK.  If not, see http://www.gnu.org/licenses/.
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

!> \brief Unified refinement criteria list

module unified_ref_crit_list

   use unified_ref_crit, only: urc

   implicit none

   private
   public :: urc_list

   type urc_list_t
      class(urc), pointer :: first => null()  !< here the list should start
      ! A pointer to last refinement criterion would be of some use only during initialisation.
   contains
      procedure :: add      !< add a single criterion to the list
      procedure :: init     !< initialize the list of refinement criteria with everything that is known at the beginning
      procedure :: cleanup  !< do a cleanup of all refinement criteria and deallocate them
      procedure :: mark     !< put all refinement marks or derefinement unmarks
   end type urc_list_t

   type(urc_list_t) :: urc_list  !< the list of refinement criteria to be applied cg-wise

contains

!> \brief Add a single criterion to the list.

   subroutine add(this, urc_p)

      implicit none

      class(urc_list_t),   intent(inout) :: this   !< an object invoking the type-bound procedure
      class(urc), pointer, intent(in)    :: urc_p  !< refinement criterion to be appended to the end of list

      class(urc), pointer :: p

      p => this%first
      if (associated(p)) then
         do while (associated(p))
            if (.not. associated(p%next)) then
               p%next => urc_p
               exit
            else
               p => p%next
            endif
         enddo
      else
         this%first => urc_p
      endif

   end subroutine add

!>
!! \brief Initialize the list of refinement criteria with everything that is known at the beginning.
!! The list can be expanded later, if necessary.
!<

   subroutine init(this)

      use constants,                          only: base_level_id
      use refinement,                         only: refine_points, refine_boxes
      use unified_ref_crit_geometrical_box,   only: urc_box
      use unified_ref_crit_geometrical_point, only: urc_point

      implicit none

      class(urc_list_t), intent(inout) :: this  !< an object invoking the type-bound procedure

      type(urc_box), pointer :: urcb
      type(urc_point), pointer :: urcp
      integer :: ip

      ! add user hook from initproblem

      ! add automatic criteria detecting shock waves

      ! add Jeans-length criterion

      ! add geometric primitives specified in problem.par
      do ip = lbound(refine_points, dim=1), ubound(refine_points, dim=1)
         if (refine_points(ip)%level > base_level_id) then
            allocate(urcp)
            urcp = urc_point(refine_points(ip))
            call this%add(urcp)
         endif
      enddo

      do ip = lbound(refine_boxes, dim=1), ubound(refine_boxes, dim=1)
         if (refine_boxes(ip)%level > base_level_id) then
            allocate(urcb)
            urcb = urc_box(refine_boxes(ip))
            call this%add(urcb)
         endif
      enddo

   end subroutine init

!> \brief Do a cleanup of all refinement criteria and deallocate them.

   subroutine cleanup(this)

      implicit none

      class(urc_list_t), intent(inout) :: this  !< an object invoking the type-bound procedure

      class(urc), pointer :: p, pn

      p => this%first
      do while (associated(p))
         pn => p%next
         deallocate(p)  ! this should also deallocate private data of p
         p => pn
         this%first => p
      enddo

   end subroutine cleanup

!> \brief Put all refinement marks or derefinement unmarks

   subroutine mark(this, cg)

      use grid_cont, only: grid_container

      implicit none

      class(urc_list_t),             intent(inout) :: this  !< an object invoking the type-bound procedure
      type(grid_container), pointer, intent(inout) :: cg    !< current grid piece

      class(urc), pointer :: p

      p => this%first
      do while (associated(p))
         call p%mark(cg)
         p => p%next
      enddo

   end subroutine mark

end module unified_ref_crit_list
