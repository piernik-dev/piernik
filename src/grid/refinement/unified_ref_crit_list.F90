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
      procedure :: init           !< initialize the list of refinement criteria with everything that is known at the beginning
      procedure :: cleanup        !< do a cleanup of all refinement criteria and deallocate them
      procedure :: all_mark       !< check refinement criteria on a given list of cg
      procedure :: plot_mark      !< check refinement criteria on a given list of cg only for iplot set
      procedure :: add_user_urcv  !< add field-based refinement criteria from initproblem
      procedure :: add_user_urc   !< add user-provided routine with refinement criteria from initproblem
      procedure, private :: add                !< append a single criterion to the list, or a chain of these
      procedure, private :: cnt                !< return the current number of defined refinement criteria
      procedure, private :: create_plotfields  !< set up qna fields for refinement criteria
      procedure, private :: summary            !< print summary of refinement criteria in use
      procedure, private :: mark               !< put all refinement marks
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
!!
!! User criteria in initproblem are supposed to be added from read_problem_par this way:
!!     call urc_list%add_user_urc(mark_user, do_plotfield)
!! Multiple routines can be added, if necessary.
!!
!! Field-based refinement criteria can be added from initproblem as well by:
!!     call urc_list%add_user_urcv(iv, ic, ref_thr, aux, rname, plotfield)
!<

   subroutine init(this)

      use constants,                          only: base_level_id, I_ONE
      use refinement,                         only: refine_points, refine_boxes, refine_zcyls, refine_vars, inactive_name, &
           &                                        jeans_ref, jeans_plot, nbody_ref
      use unified_ref_crit_geometrical_box,   only: urc_box
      use unified_ref_crit_geometrical_point, only: urc_point
      use unified_ref_crit_geometrical_zcyl,  only: urc_zcyl
      use unified_ref_crit_Jeans,             only: urc_jeans
      use unified_ref_crit_nbody,             only: urc_nbody
      use unified_ref_crit_var,               only: decode_urcv

      implicit none

      class(urc_list_t), intent(inout) :: this  !< an object invoking the type-bound procedure

      type(urc_box),   pointer :: urcb
      type(urc_point), pointer :: urcp
      type(urc_zcyl),  pointer :: urczc
      type(urc_jeans), pointer :: urcj
      type(urc_nbody), pointer :: urcn
      class(urc),      pointer :: p_urc

      integer :: ip

      ! add automatic criteria detecting shock waves
      do ip = lbound(refine_vars, 1), ubound(refine_vars, 1)
         if (trim(refine_vars(ip)%rname) /= trim(inactive_name)) then
            p_urc => decode_urcv(refine_vars(ip))
            call this%add(p_urc)
         endif
      enddo

      ! add Jeans-length criterion
      if (jeans_ref > 0.) then
         allocate(urcj)
         urcj = urc_jeans(jeans_ref, jeans_plot)
         p_urc => urcj
         call this%add(p_urc)
      endif

      ! add particle count criterion
      if (nbody_ref < huge(I_ONE) .and. nbody_ref > 0) then
         allocate(urcn)
         urcn = urc_nbody(nbody_ref)
         p_urc => urcn
         call this%add(p_urc)
      endif

      ! add geometric primitives specified in problem.par
      do ip = lbound(refine_points, dim=1), ubound(refine_points, dim=1)
         if (refine_points(ip)%level > base_level_id) then
            allocate(urcp)
            urcp = urc_point(refine_points(ip))
            p_urc => urcp
            call this%add(p_urc)
         endif
      enddo

      do ip = lbound(refine_boxes, dim=1), ubound(refine_boxes, dim=1)
         if (refine_boxes(ip)%level > base_level_id) then
            allocate(urcb)
            urcb = urc_box(refine_boxes(ip))
            p_urc => urcb
            call this%add(p_urc)
         endif
      enddo

      do ip = lbound(refine_zcyls, dim=1), ubound(refine_zcyls, dim=1)
         if (refine_zcyls(ip)%level > base_level_id) then
            allocate(urczc)
            urczc = urc_zcyl(refine_zcyls(ip))
            p_urc => urczc
            call this%add(p_urc)
         endif
      enddo

      call this%create_plotfields

   end subroutine init

!> \brief print summary of refinement criteria in use

   subroutine summary(this)

      use dataio_pub, only: msg, printinfo
      use mpisetup,   only: master

      implicit none

      class(urc_list_t), intent(inout) :: this   !< an object invoking the type-bound procedure

      logical, save :: first_run = .true.

      write(msg, '(a,i3,a)') "[unified_ref_crit_list:summary] ", this%cnt(), " criteria defined."
      if (master) then
         if (this%cnt() /= 0) then  ! unfortunately this%cnt() can not be pure due to pointer assignment
            if (first_run) then
               call printinfo(msg)
            else
               call printinfo(trim(msg) // " (update)")
            endif
         endif
      endif
      first_run = .false.

   end subroutine summary

!< \brief Add field-based refinement criteria from initproblem

   subroutine add_user_urcv(this, iv, ic, ref_thr, aux, rname, plotfield)

      use refinement,            only: ref_auto_param
      use unified_ref_crit_var,  only: urc_var

      implicit none

      class(urc_list_t), intent(inout) :: this      !< an object invoking the type-bound procedure
      integer(kind=4),   intent(in)    :: iv        !< field index in cg%q or cg%w array
      integer(kind=4),   intent(in)    :: ic        !< component index of 4D array or INVALID for 3D arrays
      real,              intent(in)    :: ref_thr   !< refinement threshold
      real,              intent(in)    :: aux       !< auxiliary parameter
      character(len=*),  intent(in)    :: rname     !< name of the refinement routine
      logical,           intent(in)    :: plotfield !< create an array to keep the value of refinement criterion

      type(urc_var), pointer :: urcv
      class(urc),    pointer :: p_urc

      allocate(urcv)
      urcv = urc_var(ref_auto_param("user", rname, ref_thr, aux, plotfield), iv, ic)
      p_urc => urcv
      call this%add(p_urc)

      call this%create_plotfields

   end subroutine add_user_urcv

!< \brief Add user-provided routine with refinement criteria from initproblem

   subroutine add_user_urc(this, user_mark, plotfield)

      use unified_ref_crit_user, only: urc_user, mark_urc_user

      implicit none

      class(urc_list_t), intent(inout) :: this       !< an object invoking the type-bound procedure
      logical,           intent(in)    :: plotfield  !< create an array to keep the value of refinement criterion
      procedure(mark_urc_user)         :: user_mark  !< user-provided routine

      type(urc_user), pointer :: urcu
      class(urc),     pointer :: p_urc

      allocate(urcu)
      urcu = urc_user(user_mark, plotfield)
      p_urc => urcu
      call this%add(p_urc)

      call this%create_plotfields

   end subroutine add_user_urc

!< \brief Set up qna fields for refinement criteria where needed

   subroutine create_plotfields(this)

      use constants,              only: INVALID, dsetnamelen
      use dataio_pub,             only: printinfo, msg, warn
      use named_array_list,       only: qna, wna
      use mpisetup,               only: master
      use unified_ref_crit_Jeans, only: urc_jeans
      use unified_ref_crit_nbody, only: urc_nbody
      use unified_ref_crit_user,  only: urc_user
      use unified_ref_crit_var,   only: urc_var

      implicit none

      class(urc_list_t), intent(inout) :: this  !< an object invoking the type-bound procedure

      class(urc), pointer :: p
      integer :: max
      character(len=dsetnamelen) :: ref_n

      max = this%cnt()

      p => this%first
      do while (associated(p))
         select type (p)
            class is (urc_var)
               if (p%plotfield .and. p%iplot == INVALID) then
                  p%iplot = new_ref_field()
                  write(msg, '(3a)') "[unified_ref_crit_list:create_plotfields] refinement criterion of type '", trim(p%rname), "' for '"
                  if (p%ic /= INVALID) then
                     write(msg(len_trim(msg)+1:), '(2a,i3,a)') trim(wna%lst(p%iv)%name), "(", p%ic, ")"
                  else
                     write(msg(len_trim(msg)+1:), '(a)') trim(qna%lst(p%iv)%name)
                  endif
                  write(msg(len_trim(msg)+1:), '(3a)') "' is stored in array '", trim(ref_n), "'"
                  if (master) call printinfo(msg)
               endif
            class is (urc_jeans)
               if (p%plotfield .and. p%iplot == INVALID) then
                  p%iplot = new_ref_field("nJ")
                  write(msg, '(3a)') "[unified_ref_crit_list:create_plotfields] Jeans refinement criterion is stored in array '", trim(ref_n), "'"
                  if (master) call printinfo(msg)
               endif
            class is (urc_nbody)
               ! unused at the moment of implementation
               if (p%plotfield .and. p%iplot == INVALID) then
                  p%iplot = new_ref_field("nparticles")
                  write(msg, '(3a)') "[unified_ref_crit_list:create_plotfields] particle count criterion is stored in array '", trim(ref_n), "'"
                  if (master) call printinfo(msg)
               endif
            class is (urc_user)
               if (p%plotfield .and. p%iplot == INVALID) then
                  p%iplot = new_ref_field()
                  write(msg, '(3a)') "[unified_ref_crit_list:create_plotfields] user-provided refinement criterion is stored in array '", trim(ref_n), "'"
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

      call this%summary

   contains

      integer function new_ref_field(suggested)

         use cg_list_global,   only: all_cg
         use dataio_pub,       only: warn
         use named_array_list, only: qna

         implicit none

         character(len=*), optional :: suggested

         integer :: i

         ref_n = ""

         if (present(suggested)) then
            if (.not. qna%exists(suggested)) then
               ref_n = suggested
               call all_cg%reg_var(ref_n)
            else
               call warn("[unified_ref_crit_list:create_plotfields:new_ref_field] suggested name '" // trim(suggested) // "' is already occupied")
            endif
         endif

         ! If max >= 100 this should fail in add2list with empty name
         if (len_trim(ref_n) <= 0) then
            do i = 1, max  ! Beware: O(n^2)
               write(ref_n, '(a,i2.2)') "ref_", i
               if (.not. qna%exists(ref_n)) exit
            enddo
            call all_cg%reg_var(ref_n)
         endif

         new_ref_field = qna%ind(ref_n)

      end function new_ref_field

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

!> \brief Put all refinement marks

   subroutine mark(this, cg)

      use grid_cont, only: grid_container

      implicit none

      class(urc_list_t),             intent(inout) :: this  !< an object invoking the type-bound procedure
      type(grid_container), pointer, intent(inout) :: cg    !< starting cg, can't use here leaves%first explicitly because of cyclic dependencies

      class(urc), pointer :: p

      p => this%first
      do while (associated(p))
         call p%mark(cg)
         p => p%next
      enddo

   end subroutine mark

!> \brief Check refinement criteria on a given list of cg

   subroutine all_mark(this, first)

      use cg_cost_data, only: I_REFINE
      use cg_list,      only: cg_list_element

      implicit none

      class(urc_list_t),              intent(inout) :: this   !< an object invoking the type-bound procedure
      type(cg_list_element), pointer, intent(in)    :: first  !< the list of cgs (usually leaves)

      type(cg_list_element), pointer :: cgl

      cgl => first
      do while (associated(cgl))
         call cgl%cg%costs%start

         call this%mark(cgl%cg)

         call cgl%cg%costs%stop(I_REFINE)
         cgl => cgl%nxt
      enddo

   end subroutine all_mark

!>
!! \brief Check refinement criteria on a given list of cg only for iplot set.
!!
!! ToDo: check whether loop order matters performance-wise.
!<

   subroutine plot_mark(this, first)

      use cg_cost_data, only: I_REFINE
      use cg_list,      only: cg_list_element
      use constants,    only: INVALID, PPP_AMR, PPP_IO
      use ppp,          only: ppp_main

      implicit none

      class(urc_list_t),              intent(inout) :: this   !< an object invoking the type-bound procedure
      type(cg_list_element), pointer, intent(in)    :: first  !< the list of cgs (usually leaves)

      type(cg_list_element), pointer :: cgl
      class(urc), pointer :: p
      character(len=*), parameter :: plot_label = "URC_map_plot"

      p => this%first
      do while (associated(p))
         if (p%iplot > INVALID) then
            call ppp_main%start(plot_label, PPP_AMR + PPP_IO)
            cgl => first
            do while (associated(cgl))
               call cgl%cg%costs%start

               call p%mark(cgl%cg)

               call cgl%cg%costs%stop(I_REFINE)
               cgl => cgl%nxt
            enddo
            call ppp_main%stop(plot_label, PPP_AMR + PPP_IO)
         endif
         p => p%next
      enddo

   end subroutine plot_mark

end module unified_ref_crit_list
