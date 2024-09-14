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
!! \brief Module controlling the repeat_step flag .
!!
!! \details Mishandling this flag (i.e. having a non reduced state across processes) may lead to deadlocks.
!! These deadlocks would likely to occur in places unrelated to the origin of the problem.
!!
!! Perhaps this module can be expressed as a type, if more synchonized logical variables are needed in the code.
!<

module repeatstep

   implicit none

   private
   public :: repeat_step, set_repeat_step, clear_repeat_step, sync_repeat_step

   logical :: repeat = .false.     ! this flag marks whether a timestep has to be repeated or not
   logical :: need_sync = .false.  ! die if unsynchronized repeat is used anywhere

   logical, parameter :: ultra_defensive = &
#ifdef DEBUG
        .true.   ! enforce unconditional global reduction or call die()
#else /* !DEBUG */
        .false.  ! allow the programmer to skip global reduction when it is known to be safe from some other sources
#endif /* DEBUG */

contains

!> \brief Return repeat value, check for synchronisation first

   logical function repeat_step()

      use dataio_pub, only: die

      implicit none

      if (need_sync) call die("[repeatstep:repeat_step] The repeat flag was possibly unsynchronized.")
      repeat_step = repeat

   end function repeat_step

!> \brief Change the state of the repeat flag, if required

   subroutine set_repeat_step(repeat_flag)

      implicit none

      logical, intent(in) :: repeat_flag  ! if .true. then set the repeat flag in this module

      need_sync = ultra_defensive .or. need_sync .or. (repeat_flag .and. .not. repeat)
      repeat = repeat .or. repeat_flag

   end subroutine set_repeat_step

!> \brief Unset the repeat flag

   subroutine clear_repeat_step(sync_flag)

      implicit none

      logical, optional, intent(in) :: sync_flag  ! if .true. then also do the synchronisation

      need_sync = ultra_defensive .or. need_sync .or. repeat
      repeat = .false.

      if (present(sync_flag)) then
         if (sync_flag) call sync_repeat_step
      endif

   end subroutine clear_repeat_step

!> \brief Synchronize the repeat flag (MPI global reduction)

   subroutine sync_repeat_step

      use allreduce, only: piernik_MPI_Allreduce
      use constants, only: pLOR

      implicit none

      need_sync = .false.
      call piernik_MPI_Allreduce(repeat, pLOR)

   end subroutine sync_repeat_step

end module repeatstep
