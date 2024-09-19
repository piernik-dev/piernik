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

!>
!! \brief This module contains list methods that operate on various
!! grid container properties and aren't related to field operations.
!! Field operations should be implemented in cg_list_dataop_t.
!<

module cg_list_misc

   use cg_list, only: cg_list_t

   implicit none

   private
   public :: cg_list_misc_t

   !> \brief Expansion of list of grid containers, not for direct use.
   type, abstract, extends(cg_list_t) :: cg_list_misc_t
   contains
      procedure :: numbered_ascii_dump  !< Construct name of emergency ASCII dump
      procedure :: ascii_dump           !< Emergency routine for quick ASCII dumps
      procedure :: prevent_prolong      !< Mark grids as untouchable for prolongation
      procedure :: enable_prolong       !< Mark grids eligible for prolongation
      procedure :: set_is_old           !< Mark grids as existing in the previous timestep
      procedure :: reset_costs          !< Gather info about measured costs and print them to the log
   end type cg_list_misc_t

contains

!> \brief Construct name of emergency ASCII dump

   subroutine numbered_ascii_dump(this, qlst, basename, a)

      use dataio_pub, only: halfstep, msg
      use global,     only: nstep, do_ascii_dump
      use mpisetup,   only: proc

      implicit none

      class(cg_list_misc_t),         intent(inout) :: this     !< list for which do the dump (usually all_cg)
      integer(kind=4), dimension(:), intent(in)    :: qlst     !< list of scalar fields to be printed
      character(len=*),              intent(in)    :: basename !< first part of the filename
      integer, optional,             intent(in)    :: a        !< additional number

      integer                              :: l, n

      if (.not. do_ascii_dump) return

      n = 2 * nstep
      if (halfstep) n = n + 1

      if (present(a)) then
         write(msg, '(a,i4,i6,i3)') trim(basename), proc, n, a
      else
         write(msg, '(a,i4,i6)')    trim(basename), proc, n
      endif
      do l = 1, len_trim(msg)
         if (msg(l:l) == " ") msg(l:l) = "_"
      enddo
      call this%ascii_dump(trim(msg), qlst)

   end subroutine numbered_ascii_dump

!>
!! \brief Emergency routine for quick ASCII dumps
!!
!! \details Absolute integer coordinates also allow seamless concatenation of dumps made by all PEs.
!!
!! \warning This routine is intended only for debugging. It is strongly discouraged to use it for data dumps in production runs.
!<

   subroutine ascii_dump(this, filename, qlst)

      use cg_list,          only: cg_list_element
      use dataio_pub,       only: msg, printio
      use named_array_list, only: qna

      implicit none

      class(cg_list_misc_t),         intent(inout) :: this     !< list for which do the dump (usually all_cg)
      character(len=*),              intent(in)    :: filename !< name to write the emergency dump (should be different on each process)
      integer(kind=4), dimension(:), intent(in)    :: qlst     !< list of scalar fields to be printed

      integer, parameter                   :: fu=30
      integer                              :: i, j, k, q
      type(cg_list_element), pointer       :: cgl

      open(fu, file=filename, status="unknown")
      write(fu, '("#",a3,2a4,a6,3a20)', advance='no')"i", "j", "k", "level", "x(i)", "y(j)", "z(k)"
      do q = lbound(qlst(:), dim=1), ubound(qlst(:), dim=1)
         write(fu, '(a20)', advance='no') trim(qna%lst(qlst(q))%name)
      enddo
      write(fu, '(/)')

      cgl => this%first
      do while (associated(cgl))
         do i = cgl%cg%is, cgl%cg%ie
            do j = cgl%cg%js, cgl%cg%je
               do k = cgl%cg%ks, cgl%cg%ke
                  write(fu, '(3i4,i6,3es20.11e3)', advance='no') i, j, k, cgl%cg%l%id, cgl%cg%x(i), cgl%cg%y(j), cgl%cg%z(k)
                  do q = lbound(qlst(:), dim=1), ubound(qlst(:), dim=1)
                     write(fu, '(es20.11e3)', advance='no') cgl%cg%q(qlst(q))%arr(i, j, k)
                  enddo
                  write(fu, '()')
               enddo
               write(fu, '()')
            enddo
            write(fu, '()')
         enddo
         write(fu, '()')
         cgl => cgl%nxt
      enddo

      close(fu)

      write(msg,'(3a)') "[cg_list_misc:ascii_dump] Wrote dump '",filename,"'"
      call printio(msg)

   end subroutine ascii_dump

!> \brief Mark grids as untouchable for prolongation

   subroutine prevent_prolong(this)

      use cg_list, only: cg_list_element

      implicit none

      class(cg_list_misc_t), intent(in) :: this  !< object invoking type-bound procedure

      type(cg_list_element), pointer :: cgl

      cgl => this%first
      do while (associated(cgl))
         cgl%cg%ignore_prolongation = .true.
         cgl => cgl%nxt
      enddo

   end subroutine prevent_prolong

!> \brief Mark grids as eligible for prolongation

   subroutine enable_prolong(this)

      use cg_list, only: cg_list_element

      implicit none

      class(cg_list_misc_t), intent(in) :: this  !< object invoking type-bound procedure

      type(cg_list_element), pointer :: cgl

      cgl => this%first
      do while (associated(cgl))
         cgl%cg%ignore_prolongation = .false.
         cgl => cgl%nxt
      enddo

   end subroutine enable_prolong

!> \brief Mark grids as existing in the previous timestep

   subroutine set_is_old(this)

      use cg_list, only: cg_list_element

      implicit none

      class(cg_list_misc_t), intent(in) :: this  !< object invoking type-bound procedure

      type(cg_list_element), pointer :: cgl

      cgl => this%first
      do while (associated(cgl))
         cgl%cg%is_old = .true.
         cgl => cgl%nxt
      enddo

   end subroutine set_is_old

!< \brief Gather info about measured costs and print them to the log

   subroutine reset_costs(this)

      use cg_list, only: cg_list_element

      implicit none

      class(cg_list_misc_t), intent(in) :: this  !< object invoking type-bound procedure

      type(cg_list_element), pointer :: cgl

      ! clear the data before next stage
      cgl => this%first
      do while (associated(cgl))
         cgl%cg%old_costs%wtime = cgl%cg%costs%wtime
         call cgl%cg%costs%reset
         cgl => cgl%nxt
      enddo

   end subroutine reset_costs

end module cg_list_misc
