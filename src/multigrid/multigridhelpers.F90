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

!!$ ============================================================================
!>
!! \brief This module contains variables and routines useful mostly for developing and debugging
!<

module multigridhelpers
! pulled by MULTIGRID
   implicit none

   private

   public :: set_dirty, check_dirty, numbered_ascii_dump, zero_boundaries
   public :: do_ascii_dump, dirty_debug, dirty_label, dirtyH, dirtyL

   ! namelist parameters
   logical            :: do_ascii_dump                      !< to dump, or not to dump: that is a question (ascii)
   logical            :: dirty_debug                        !< Initialize everything with some insane values (dirtyH, defined below) and check if they can propagate
   integer, parameter    :: dl_len = 64                     !< length of label buffer
   character(len=dl_len) :: dirty_label                     !< buffer for label for check_dirty subroutine
   real, parameter    :: big_s =  huge(real(1.0,4))         !< largest single-precision number
   real, parameter    :: dirtyH = big_s                     !< If dirty_debug then initialize arrays with this insane value
   real, parameter    :: dirtyL = sqrt(big_s)               !< If dirty_debug then check if the solution got contaminated by dirtyH by checking if it contains anything above dirtyL

contains

!> \todo write also general set_dirty_array and check_dirty_array routines so there will be no need to use dirtyH and dirtyL outside multigridhelpers module.

!!$ ============================================================================
!>
!! \brief This routine pollutes selected multigrid variable with an insane value dirtyH.
!! \details If anything in the multigrid works by accident, through compiler-dependent initialization or unintentional relying on outdated values,
!! the insane value should pollute the solution in an easily visible way.
!<

   subroutine set_dirty(iv)

      use cg_list_global, only: all_cg
      use dataio_pub,     only: die

      implicit none

      integer, intent(in) :: iv   !< index of variable in cg%q(:) which we want to pollute

      if (.not. dirty_debug) return

      if (iv < lbound(all_cg%q_lst, dim=1) .or. iv > ubound(all_cg%q_lst, dim=1)) call die("[multigridhelpers:set_dirty] Invalid variable index.")

      call all_cg%set_q_value(iv, dirtyH)

   end subroutine set_dirty

!!$ ============================================================================
!>
!! \brief This routine checks for detectable traces of set_dirty calls.
!<

   subroutine check_dirty(curl, iv, label, expand)

      use cg_list_global, only: all_cg
      use constants,      only: ndims
      use dataio_pub,     only: die, warn, msg
      use domain,         only: dom
      use gc_list,        only: cg_list_element
      use cg_list_lev,    only: cg_list_level
      use mpisetup,       only: proc

      implicit none

      type(cg_list_level), pointer, intent(in) :: curl   !< level which we are checking
      integer,                      intent(in) :: iv     !< index of variable in cg%q(:) which we want to pollute
      character(len=*),             intent(in) :: label  !< label to indicate the origin of call
      integer(kind=4), optional,    intent(in) :: expand !< also check guardcells

      integer :: i, j, k, ng
      type(cg_list_element), pointer :: cgl

      if (.not. dirty_debug) return
      if (iv < lbound(all_cg%q_lst, dim=1) .or. iv > ubound(all_cg%q_lst, dim=1)) call die("[multigridhelpers:check_dirty] Invalid variable index.")

      if (present(expand) .and. dom%eff_dim==ndims) then ! for 1D and 2D one should define ng_x,ng_y and ng_z
         ng = min(dom%nb, expand)
      else
         ng = 0
      endif

      cgl => curl%first
      do while (associated(cgl))
         do k = cgl%cg%ks-ng*dom%D_z, cgl%cg%ke+ng*dom%D_z
            do j = cgl%cg%js-ng*dom%D_y, cgl%cg%je+ng*dom%D_y
               do i = cgl%cg%is-ng*dom%D_x, cgl%cg%ie+ng*dom%D_x
                  if (abs(cgl%cg%q(iv)%arr(i, j, k)) > dirtyL) then
                     ! if (count([i<cgl%cg%is .or. i>cgl%cg%ie, j<cgl%cg%js .or. j>cgl%cg%je, k<cgl%cg%ks .or. k>cgl%cg%ke]) <=1) then ! excludes corners
                     write(msg, '(3a,i4,a,i3,a,i5,3a,3(i3,a),g20.12)') "[multigridhelpers:check_dirty] ", trim(label), "@", proc, " lvl^", cgl%cg%level_id, " cg#", cgl%cg%grid_id, &
                          &                                            " '", trim(all_cg%q_lst(iv)%name), "'(", i, ",", j, ",", k, ") = ", cgl%cg%q(iv)%arr(i, j, k)
                     call warn(msg)
                     ! endif
                  endif
               enddo
            enddo
         enddo
         cgl => cgl%nxt
      enddo

   end subroutine check_dirty

!!$ ============================================================================
!>
!! \brief Emergency routine for quick ASCII dumps
!!
!! \details Absolute integer coordinates also allow seamless concatenation of dumps made by all PEs.
!<

   subroutine ascii_dump(filename)

      use cg_list_global, only: all_cg
      use constants,      only: xdim, ydim, zdim
      use dataio_pub,     only: msg, printio
      use gc_list,        only: cg_list_element
      use mpisetup,       only: master
      use multigridvars,  only: source, solution, defect, correction

      implicit none

      character(len=*), intent(in) :: filename !< name to write the emergency dump

      integer, parameter :: fu=30
      integer, dimension(4) :: qlst
      integer            :: i, j, k, q
      type(cg_list_element), pointer :: cgl

      if (.not. do_ascii_dump) return

      qlst = [ source, solution, defect, correction ]

      open(fu, file=filename, status="unknown")
      write(fu, '("#",a3,2a4,a6,3a20)', advance='no')"i", "j", "k", "level", "x(i)", "y(j)", "z(k)"
      do q = lbound(qlst, dim=1), ubound(qlst, dim=1)
         write(fu, '(a20)', advance='no') trim(all_cg%q_lst(qlst(q))%name)
      enddo
      write(fu, '(/)')

      cgl => all_cg%first
      do while (associated(cgl))
         do i = cgl%cg%is, cgl%cg%ie
            do j = cgl%cg%js, cgl%cg%je
               do k = cgl%cg%ks, cgl%cg%ke
                  write(fu, '(3i4,i6,3es20.11e3)', advance='no') i-cgl%cg%is+cgl%cg%off(xdim), j-cgl%cg%js+cgl%cg%off(ydim), k-cgl%cg%ks+cgl%cg%off(zdim), &
                       &                           cgl%cg%level_id, cgl%cg%x(i), cgl%cg%y(j), cgl%cg%z(k)
                  do q = lbound(qlst, dim=1), ubound(qlst, dim=1)
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

      if (master) then
         write(msg,'(3a)') "[multigridhelpers:ascii_dump] Wrote dump '",filename,"'"
         call printio(msg)
      endif

   end subroutine ascii_dump

!!$ ============================================================================
!>
!! \brief Construct name of emergency ASCII dump
!<

   subroutine numbered_ascii_dump(basename, a)

      use dataio_pub, only: halfstep, msg
      use global,     only: nstep
      use mpisetup,   only: proc

      implicit none

      character(len=*), intent(in)  :: basename !< first part of the filename
      integer, optional, intent(in) :: a        !< additional number

      integer             :: l, n

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
      call ascii_dump(trim(msg))

   end subroutine numbered_ascii_dump


!>
!! \brief Clear boundary values
!<

   subroutine zero_boundaries(curl)

      use gc_list,     only: cg_list_element
      use cg_list_lev, only: cg_list_level

      implicit none

      type(cg_list_level), pointer, intent(in) :: curl  !< level for which clear the boundary values

      type(cg_list_element), pointer :: cgl

      cgl => curl%first
      do while (associated(cgl))
         cgl%cg%mg%bnd_x(:,:,:) = 0.
         cgl%cg%mg%bnd_y(:,:,:) = 0.
         cgl%cg%mg%bnd_z(:,:,:) = 0.
         cgl => cgl%nxt
      enddo

   end subroutine zero_boundaries

end module multigridhelpers
