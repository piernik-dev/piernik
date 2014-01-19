!
! PIERNIK Code Copyright (C) 2006-2013 Michal Hanasz
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
!! \brief Module for handling signals
!!
!! \details
!!    This module allow to catch signals sent to Piernik from external
!!    sources, i.e. PBS
!<
module signalhandler
! pulled by ANY

   implicit none

   private
   public :: SIGTERM, SIGINT, register_sighandler

   integer(kind=4), parameter :: SIGINT = 2
   integer(kind=4), parameter :: SIGTERM = 15

   interface
      integer function signal_handler(signum)
         implicit none
         integer, intent(in) :: signum
      end function
   end interface

contains

   subroutine register_sighandler(signum, func)
#ifdef __INTEL_COMPILER
      use ifport, only: signal
#endif /* __INTEL_COMPILER */
      implicit none
      integer(kind=4), intent(in) :: signum
      procedure(signal_handler) :: func
      integer :: err

#if defined(__INTEL_COMPILER)
      err = signal(signum, func, -1)
#elif defined(__GFORTRAN__)
      err = signal(signum, func)
#else /* !(__INTEL_COMPILER || __GFORTRAN__) */
      err = -1
#endif /* !(__INTEL_COMPILER || __GFORTRAN__) */

   end subroutine register_sighandler

end module signalhandler

#ifdef __TEST
program test
   use signalhandler, only: SIGTERM, register_handler
   implicit none

   call register_handler(SIGTERM, warning_sigterm)
   print *, "The current process ID is ", getpid()
   print *, "KILL ME!"
   call sleep(600)
   print *, 'You were either too slow or something went wrong'

contains

   integer function warning_sigterm(signum)
      implicit none
      integer, intent(in) :: signum
      print *, 'Process interrupted (SIGTERM), exiting...'
      print *, '=.='
      warning_sigterm = 0
   end function warning_sigterm

end program test
#endif /* __TEST */
