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

!> \brief This module contains variables related to refinement

module refinement_crit_list

   use refinement_filters, only: ref_crit

   implicit none

   private
   public :: refines2list, user_ref2list, cleanup_refinement, auto_refine_derefine

   type(ref_crit), dimension(:), allocatable :: ref_crit_list  !< definitions of user-supplied automatic refinement criteria, processed and checked

contains

!> \brief free the memory

   subroutine cleanup_refinement

     implicit none

     if (allocated(ref_crit_list)) deallocate(ref_crit_list)

   end subroutine cleanup_refinement

!> \brief convert refinement criteria parameters to a list

   subroutine refines2list

      use constants,  only: INVALID
      use dataio_pub, only: die, msg
      use refinement, only: refine_vars, inactive_name

      implicit none

      integer :: i, c
      integer(kind=4) :: iv
      integer(kind=4), dimension(:), allocatable :: ic

      do i = lbound(refine_vars, 1), ubound(refine_vars, 1)
         call identify_field(refine_vars(i)%rvar, iv, ic)
         if (iv /= INVALID .and. trim(refine_vars(i)%rname) /= trim(inactive_name)) then
            if (.not. allocated(ic)) call die("[refinement_crit_list:refines2list] .not. allocated(ic))")
            if (size(ic) <= 0) then
               write(msg,'(3a)')"[refinement_crit_list:refines2list] iv /= INVALID and  unrecognized field name '",trim(refine_vars(i)%rname),"'"
               call die(msg)
            else
               do c = lbound(ic, dim=1), ubound(ic, dim=1)
                  call user_ref2list(iv, ic(c), refine_vars(i)%ref_thr, refine_vars(i)%deref_thr, refine_vars(i)%aux, refine_vars(i)%rname, refine_vars(i)%plotfield)
               enddo
            endif
         endif
         if (allocated(ic)) deallocate(ic)
      enddo

   end subroutine refines2list

!> \brief Add a user-defined criteria to the list

   subroutine user_ref2list(iv, ic, ref_thr, deref_thr, aux, rname, plotfield)

      use cg_list_global,     only: all_cg
      use constants,          only: INVALID, dsetnamelen
      use dataio_pub,         only: warn, die, printinfo, msg
      use mpisetup,           only: master
      use named_array_list,   only: qna, wna
      use refinement_filters, only: refine_on_gradient, refine_on_relative_gradient, refine_on_second_derivative
      use refinement,         only: inactive_name

      implicit none

      integer(kind=4),  intent(in) :: iv        !< field index in cg%q or cg%w array
      integer(kind=4),  intent(in) :: ic        !< component index of 4D array or INVALID for 3D arrays
      real,             intent(in) :: ref_thr   !< refinement threshold
      real,             intent(in) :: deref_thr !< derefinement threshold
      real,             intent(in) :: aux       !< auxiliary parameter
      character(len=*), intent(in) :: rname     !< name of the refinement routine
      logical,          intent(in) :: plotfield !< create an array to keep the value of refinement criterion

      character(len=dsetnamelen) :: ref_n
      integer :: i

      if (iv == INVALID) then
         call warn("[refinement_crit_list:user_ref2list] invalid field. Ignored.")
         return
      endif
      if (.not. allocated(ref_crit_list)) allocate(ref_crit_list(0))
      ref_crit_list = [ ref_crit_list, ref_crit(iv, ic, INVALID, ref_thr, deref_thr, aux, null()) ]
      select case (trim(rname))
         case ("grad")
            ref_crit_list(ubound(ref_crit_list, dim=1))%refine => refine_on_gradient
         case ("relgrad")
            ref_crit_list(ubound(ref_crit_list, dim=1))%refine => refine_on_relative_gradient
         case ("Loechner", "second_order", "d2")
            ref_crit_list(ubound(ref_crit_list, dim=1))%refine => refine_on_second_derivative

!> \todo Implement Richardson extrapolation method, as described in M. Berger papers

         case (trim(inactive_name)) ! do nothing
         case default
            call die("[refinement_crit_list:user_ref2list] unknown refinement detection routine")
      end select

      if (plotfield) then
         do i = 1, ubound(ref_crit_list, dim=1)  ! Beware: O(n^2)
            write(ref_n, '(a,i2.2)') "ref_", i
            if (.not. qna%exists(ref_n)) exit
         enddo
         call all_cg%reg_var(ref_n)
         ref_crit_list(ubound(ref_crit_list, dim=1))%iplot = qna%ind(ref_n)
         write(msg, '(3a)') "[refinement_crit_list:user_ref2list] refinement criterion of type '", trim(rname), "' for '"
         if (ic /= INVALID) then
            write(msg, '(3a,i3,a)') trim(msg), trim(wna%lst(iv)%name), "(", ic, ")"
         else
            write(msg, '(2a)') trim(msg), trim(qna%lst(iv)%name)
         endif
         write(msg, '(4a)') trim(msg), "' is stored in array '", trim(ref_n), "'"
         if (master) call printinfo(msg)
      endif

      !> \todo try to detect doubled criteria

   end subroutine user_ref2list

!>
!! \brief Apply automatic refinement citeria
!!
!! \details The leaves argument should normally be cg_leaves::leaves. W cannot use it directly here because of circular dependencies
!! Note that if you pass only subset of leaves (i.e. single level or somehow filtered cg list), the (de)refinement will be marked only on that list.
!! The rest of the domain will stay unaffected or be corrected for refinement defects.
!<

   subroutine auto_refine_derefine(plots_only)

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: INVALID
      use dataio_pub,       only: die
      use named_array_list, only: qna, wna

      implicit none

      logical, optional, intent(in) :: plots_only

      integer :: i
      logical :: var3d, skip
      type(cg_list_element), pointer :: cgl
      real, dimension(:,:,:), pointer :: p3d

      if (.not. allocated(ref_crit_list)) return
      do i = lbound(ref_crit_list, dim=1), ubound(ref_crit_list, dim=1)

         skip = .false.
         if (present(plots_only)) then
            if (plots_only) skip = (ref_crit_list(i)%iplot == INVALID)
         endif
         if (skip) cycle

         var3d = (ref_crit_list(i)%ic == INVALID)

         if (var3d) then
            if (ref_crit_list(i)%iv<lbound(qna%lst, dim=1) .or. ref_crit_list(i)%iv>ubound(qna%lst, dim=1)) &
                 call die("[refinement_crit_list:auto_refine_derefine] 3D index out of range")
         else
            if (ref_crit_list(i)%iv<lbound(wna%lst, dim=1) .or. ref_crit_list(i)%iv>ubound(wna%lst, dim=1)) &
                 call die("[refinement_crit_list:auto_refine_derefine] 4D index out of range")
            if (ref_crit_list(i)%ic <= 0 .or. ref_crit_list(i)%ic > wna%lst(ref_crit_list(i)%iv)%dim4) &
                 call die("[refinement_crit_list:auto_refine_derefine] component out of range")
         endif

         cgl => leaves%first
         do while (associated(cgl))
            if (any(cgl%cg%leafmap)) then
               if (var3d) then
                  p3d => cgl%cg%q(ref_crit_list(i)%iv)%arr
               else
                  associate (a=>cgl%cg%w(ref_crit_list(i)%iv)%arr)
                     p3d(lbound(a, dim=2):, lbound(a, dim=3):, lbound(a, dim=4):) => cgl%cg%w(ref_crit_list(i)%iv)%arr(ref_crit_list(i)%ic, :, :, :)
                  end associate
               endif
               call ref_crit_list(i)%refine(cgl%cg, p3d)
            endif
            cgl => cgl%nxt
         enddo
      enddo

   end subroutine auto_refine_derefine

!> \brief Identify field name and return indices to cg%q or cg%w arrays

   subroutine identify_field(vname, iv, ic)

      use constants,        only: INVALID, cbuff_len
      use dataio_pub,       only: msg, warn
      use fluidindex,       only: iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz, iarr_all_en
      use named_array_list, only: qna, wna
      use refinement,       only: inactive_name

      implicit none

      character(len=cbuff_len),                   intent(in)  :: vname !< string specifying the field on
      integer(kind=4),                            intent(out) :: iv    !< field index in cg%q or cg%w array
      integer(kind=4), dimension(:), allocatable, intent(out) :: ic    !< component index array (cg%w(iv)%arr(ic,:,:,:)) or INVALID for 3D arrays

      iv = INVALID

      if (trim(vname) == trim(inactive_name)) return ! ignore this

      if (qna%exists(trim(vname))) then
         iv = qna%ind(trim(vname))
         if (iv /= INVALID) then
            allocate(ic(1))
            ic = INVALID
            return ! this is a 3d array name
         endif
      endif

      if (trim(vname) == "dens") then
         allocate(ic(lbound(iarr_all_dn, dim=1):ubound(iarr_all_dn, dim=1)))
         iv = wna%fi
         ic = iarr_all_dn
         return
      else if (trim(vname) == "velx") then
         allocate(ic(lbound(iarr_all_mx, dim=1):ubound(iarr_all_mx, dim=1)))
         iv = wna%fi
         ic = iarr_all_mx
         return
      else if (trim(vname) == "vely") then
         allocate(ic(lbound(iarr_all_my, dim=1):ubound(iarr_all_my, dim=1)))
         iv = wna%fi
         ic = iarr_all_my
         return
      else if (trim(vname) == "velz") then
         allocate(ic(lbound(iarr_all_mz, dim=1):ubound(iarr_all_mz, dim=1)))
         iv = wna%fi
         ic = iarr_all_mz
         return
      else if (trim(vname) == "ener") then
         allocate(ic(lbound(iarr_all_en, dim=1):ubound(iarr_all_en, dim=1)))
         iv = wna%fi
         ic = iarr_all_en
         return
      endif
      !> \todo identify here all {den,vl[xyz],ene}{d,n,i}
      !> \todo introduce possibility to operate on pressure or other indirect fields

      write(msg,'(3a)')"[refinement_crit_list:identify_field] Unidentified refinement variable: '",trim(vname),"'"
      call warn(msg)

      allocate(ic(0))

   end subroutine identify_field

end module refinement_crit_list
