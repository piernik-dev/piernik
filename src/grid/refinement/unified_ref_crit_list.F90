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
      procedure :: add                !< append a single criterion to the list, or a chain of these
      procedure :: init               !< initialize the list of refinement criteria with everything that is known at the beginning
      procedure :: cnt                !< return the current number of defined refinement criteria
      procedure :: cleanup            !< do a cleanup of all refinement criteria and deallocate them
      procedure :: all_mark           !< check refinement criteria on a given list of cg
      procedure :: plot_mark          !< check refinement criteria on a given list of cg only for iplot set
      procedure :: create_plotfields  !< set up qna fields for refinement criteria
      procedure, private :: mark      !< put all refinement marks or derefinement unmarks
   end type urc_list_t

   type(urc_list_t) :: urc_list  !< the list of refinement criteria to be applied cg-wise

contains

!> \brief Append a single criterion to the list. A chain should get correctly appended as well.

   subroutine add(this, urc_p)

      use dataio_pub, only: warn
      use mpisetup,   only: master

      implicit none

      class(urc_list_t),   intent(inout) :: this   !< an object invoking the type-bound procedure
      class(urc), pointer, intent(in)    :: urc_p  !< refinement criterion to be appended to the end of list

      class(urc), pointer :: p

      if (.not. associated(urc_p)) then
         if (master) call warn("[unified_ref_crit_list:add] .not. associated(urc_p)")
         return
      endif

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
      use dataio_pub,                         only: msg, printinfo
      use mpisetup,                           only: master
      use refinement,                         only: refine_points, refine_boxes, refine_vars, inactive_name
      use unified_ref_crit_geometrical_box,   only: urc_box
      use unified_ref_crit_geometrical_point, only: urc_point
      use unified_ref_crit_var,               only: decode_urcv

      implicit none

      class(urc_list_t), intent(inout) :: this  !< an object invoking the type-bound procedure

      type(urc_box),   pointer :: urcb
      type(urc_point), pointer :: urcp
      integer :: ip

      ! add user hook from initproblem

      ! add automatic criteria detecting shock waves

      do ip = lbound(refine_vars, 1), ubound(refine_vars, 1)
         if (trim(refine_vars(ip)%rname) /= trim(inactive_name)) then
            call this%add(decode_urcv(refine_vars(ip)))
         endif
      enddo

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

      call this%create_plotfields

      write(msg, '(a,i3,a)') "[unified_ref_crit_list:init] ", this%cnt(), " criteria defined."
      if (master) call printinfo(msg)

   end subroutine init

!< \brief Set up qna fields for refinement criteria where needed

   subroutine create_plotfields(this)

      use cg_list_global,   only: all_cg
      use constants,        only: INVALID, dsetnamelen
      use dataio_pub,       only: printinfo, msg, warn
      use named_array_list, only: qna, wna
      use mpisetup,         only: master
      use unified_ref_crit_var, only: urc_var

      implicit none

      class(urc_list_t), intent(inout) :: this  !< an object invoking the type-bound procedure

      class(urc), pointer :: p
      integer :: max, i
      character(len=dsetnamelen) :: ref_n

      max = this%cnt()

      p => this%first
      do while (associated(p))
         select type (p)
            class is (urc_var)
               if (p%plotfield .and. p%iplot == INVALID) then
                  do i = 1, max  ! Beware: O(n^2)
                     write(ref_n, '(a,i2.2)') "ref_", i
                     if (.not. qna%exists(ref_n)) exit
                  enddo
                  call all_cg%reg_var(ref_n)
                  p%iplot = qna%ind(ref_n)
                  write(msg, '(3a)') "[unified_ref_crit_list:create_plotfields] refinement criterion of type '", trim(p%rname), "' for '"
                  if (p%ic /= INVALID) then
                     write(msg, '(3a,i3,a)') trim(msg), trim(wna%lst(p%iv)%name), "(", p%ic, ")"
                  else
                     write(msg, '(2a)') trim(msg), trim(qna%lst(p%iv)%name)
                  endif
                  write(msg, '(4a)') trim(msg), "' is stored in array '", trim(ref_n), "'"
                  if (master) call printinfo(msg)
               endif
            class default
               if (p%iplot /= INVALID) then
                  write(msg, '(3a)') "[unified_ref_crit_list:create_plotfields] some refinement criterion may be stored in array '", trim(qna%lst(p%iplot)%name), "'"
                  if (master) call warn(msg)
               endif
         end select
         p => p%next
      enddo

   end subroutine create_plotfields

!> \brief return the current number of defined refinement criteria

   function cnt(this)

      implicit none

      class(urc_list_t), intent(inout) :: this  !< an object invoking the type-bound procedure

      class(urc), pointer :: p
      integer :: cnt

      cnt = 0
      p => this%first
      do while (associated(p))
         cnt  = cnt + 1
         p => p%next
      enddo
   end function cnt

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

!> \brief Check refinement criteria on a given list of cg

   subroutine all_mark(this, first)

      use cg_list,   only: cg_list_element

      implicit none

      class(urc_list_t),              intent(inout) :: this        !< an object invoking the type-bound procedure
      type(cg_list_element), pointer, intent(in)    :: first         !< the list of cgs (usually leaves)

      type(cg_list_element), pointer :: cgl

      cgl => first
      do while (associated(cgl))
         call this%mark(cgl%cg)
         cgl => cgl%nxt
      enddo

   end subroutine all_mark

!>
!! \brief Check refinement criteria on a given list of cg only for iplot set
!!
!! Perhaps it would be more optimal to exchange the loops and check if there is anything to do at all
!<

   subroutine plot_mark(this, first)

      use cg_list,   only: cg_list_element
      use constants, only: INVALID

      implicit none

      class(urc_list_t),              intent(inout) :: this        !< an object invoking the type-bound procedure
      type(cg_list_element), pointer, intent(in)    :: first         !< the list of cgs (usually leaves)

      type(cg_list_element), pointer :: cgl
      class(urc), pointer :: p

      p => this%first
      do while (associated(p))
         if (p%iplot > INVALID) then
            cgl => first
            do while (associated(cgl))
               call p%mark(cgl%cg)
               cgl => cgl%nxt
            enddo
         endif
         p => p%next
      enddo

   end subroutine plot_mark

end module unified_ref_crit_list
