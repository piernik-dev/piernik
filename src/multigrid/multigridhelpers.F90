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

   public :: numbered_ascii_dump
   public :: do_ascii_dump, dirty_label

   ! namelist parameters
   logical            :: do_ascii_dump                      !< to dump, or not to dump: that is a question (ascii)
   integer, parameter    :: dl_len = 64                     !< length of label buffer
   character(len=dl_len) :: dirty_label                     !< buffer for label for check_dirty subroutine
   real, parameter    :: big_s =  huge(real(1.0,4))         !< largest single-precision number

contains

!> \todo write also general set_dirty_array and check_dirty_array routines so there will be no need to use dirtyH and dirtyL outside multigridhelpers module.

!!$ ============================================================================
!>
!! \brief Construct name of emergency ASCII dump
!<

   subroutine numbered_ascii_dump(basename, a)

      use cg_list_global, only: all_cg
      use dataio_pub,     only: halfstep, msg
      use global,         only: nstep
      use mpisetup,       only: proc
      use multigridvars,  only: source, solution, defect, correction

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
      call all_cg%ascii_dump(trim(msg), [ source, solution, defect, correction ])

   end subroutine numbered_ascii_dump

end module multigridhelpers
