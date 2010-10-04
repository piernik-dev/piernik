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
!    Initial implemetation of PIERNIK code was based on TVD split MHD code by
!    Ue-Li Pen
!        see: Pen, Arras & Wong (2003) for algorithm and
!             http://www.cita.utoronto.ca/~pen/MHD
!             for original source code "mhd.f90"
!
!    For full list of developers see $PIERNIK_HOME/license/pdt.txt
!
#include "piernik.def"
!>
!! \brief (KK) Module responsible for error handling
!!
!<
module errh


   implicit none

   private

   public :: die, warn, printinfo, namelist_errh, msg

   character(len=512) :: msg
   integer, parameter :: T_ERR = 1, T_WARN = T_ERR + 1, T_INFO = T_WARN + 1

   include 'mpif.h'

   contains

!-----------------------------------------------------------------------------

   subroutine colormessage(nm, mode)

      use dataio_public, only : log_file, log_lun, dataio_initialized

      implicit none

      character(len=*), intent(in) :: nm
      integer, intent(in) :: mode

      character(len=4), parameter :: ansi_black  = char(27)//"[0m"
      character(len=7), parameter :: ansi_red    = char(27)//'[1;31m'
      character(len=7), parameter :: ansi_green  = char(27)//'[1;32m'
      character(len=7), parameter :: ansi_yellow = char(27)//'[1;33m'

      character(len=7)  :: ansicolor
      character(len=8) :: msg_type_str
      integer :: proc, ierr

      select case(mode)
         case (T_ERR)
            ansicolor = ansi_red
            msg_type_str = "Error"
         case (T_WARN)
            ansicolor = ansi_yellow
            msg_type_str = "Warning"
         case (T_INFO)
            ansicolor = ansi_green
            msg_type_str = "Info"
         case default
            ansicolor = ansi_black
            msg_type_str = ""
      end select

      call MPI_comm_rank(MPI_COMM_WORLD, proc, ierr)

      write(*,'(a,i5,3a)') trim(ansicolor)//trim(msg_type_str)//" @"//ansi_black, proc, ': "', trim(nm), '"'       ! QA_WARN

      if (dataio_initialized) then
         open(log_lun, file=log_file, position='append')
         if (proc == 0 .and. mode == T_ERR) write(log_lun,'(/,a,/)')"###############     Crashing     ###############"
         write(log_lun,'(2a,i5,3a)')trim(msg_type_str)," @", proc, ': "', trim(nm), '"'
         close(log_lun)
      else
         write(*,'(a,i5,a)') ansi_yellow//"[errh:colormessage] dataio_initialized == .false. @", proc, ansi_black  ! QA_WARN
      endif

   end subroutine colormessage

!-----------------------------------------------------------------------------

   subroutine printinfo(nm)

      implicit none

      character(len=*), intent(in) :: nm

      call colormessage(nm, T_INFO)

   end subroutine printinfo

!-----------------------------------------------------------------------------

   subroutine warn(nm)

      implicit none

      character(len=*), intent(in) :: nm

      call colormessage(nm, T_WARN)

   end subroutine warn

!-----------------------------------------------------------------------------
   !! BEWARE: routine is not finished, it should kill PIERNIK gracefully

   subroutine die(nm, allprocs)

      implicit none

      character(len=*), intent(in)  :: nm
      integer, optional, intent(in) :: allprocs

      integer :: ierr

      call colormessage(nm, T_ERR)

      if(present(allprocs)) then
         if(allprocs /= 0) then
            call MPI_Barrier(MPI_COMM_WORLD, ierr)
            call MPI_Finalize(ierr)
         endif
      endif

      stop

   end subroutine die

!-----------------------------------------------------------------------------

   subroutine namelist_errh(ierrh,nm)
      implicit none
      integer, intent(in) :: ierrh
      character(len=*), intent(in) :: nm

      select case (ierrh)
         case (19)
            write(*,*) "severe (19): Invalid reference to variable in ",trim(nm), " namelist"                        ! QA_WARN
            write(*,*) "One of the following conditions occurred: "                                                  ! QA_WARN
            write(*,*) "    * The variable was not a member of the namelist group."                                  ! QA_WARN
            write(*,*) "    * An attempt was made to subscript a scalar variable."                                   ! QA_WARN
            write(*,*) "    * A subscript of the array variable was out-of-bounds."                                  ! QA_WARN
            write(*,*) "    * An array variable was specified with too many or too few subscripts for the variable." ! QA_WARN
            write(*,*) "    * An attempt was made to specify a substring of a noncharacter variable or array name."  ! QA_WARN
            write(*,*) "    * A substring specifier of the character variable was out-of-bounds."                    ! QA_WARN
            write(*,*) "    * A subscript or substring specifier of the variable was not an integer constant."       ! QA_WARN
            write(*,*) "    * An attempt was made to specify a substring by using an unsubscripted array variable."  ! QA_WARN
            stop
         case (-1)
            write(*,*) "Namelist: ",trim(nm)," not found in problem.par"                                             ! QA_WARN
            stop
         case (5010)
            write(*,*) "One of the variables found in problem.par doesn't belong to ",trim(nm), " namelist"          ! QA_WARN
            stop
         case (239)
            write(*,*) "One of the variables found in problem.par doesn't belong to ",trim(nm), " namelist"          ! QA_WARN
            stop
         case (0)
         case default
            write(*,*) 'Unknown error (', ierrh,') in namelist ',trim(nm)                                            ! QA_WARN
            stop
      endselect

   end subroutine namelist_errh

!-----------------------------------------------------------------------------
end module errh
