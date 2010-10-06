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

module dataio_public

   use types, only : hdf

   implicit none

   public

   real               :: tend                   !< simulation time to end
   real               :: wend                   !< wall clock time to end (in hours)

   integer            :: nend                   !< number of the step to end simulation
   integer            :: nstep_start            !< number of start timestep
   integer            :: nhdf                   !< current number of hdf file
   integer            :: nres                   !< current number of restart file
   integer            :: nlog                   !< current number of log file
   integer            :: ntsl                   !< current number of timeslice file
   integer            :: nrestart               !< number of restart file to be read while restart is not set to ''
   integer            :: step_hdf               !< number of simulation timestep corresponding to values dumped in hdf file

   character(len=128) :: log_file               !< path to the current log file

   logical            :: halfstep = .false.     !< true when X-Y-Z sweeps are done and Z-Y-X are not

   !! ToDo:
   !!  Currently to use PGPLOT you need to:
   !!   1. set one of i{x,y,z} to positive value and zero to the others
   !!   2. use only one-element vars array in problem.par
   !!
   logical            :: vizit = .false.        !< performe "live" vizualization using pgplot (BEWARE: highly experimental)
   real               :: fmin                   !< minimum on pgplot scale
   real               :: fmax                   !< maximum on pgplot scale

   type(hdf)          :: chdf                   !< container for some vital simulation parameters
   character(len=16)  :: domain                 !< string to choose if boundaries have to be dumped in hdf files
   real               :: last_hdf_time          !< time in simulation of the last resent hdf file dump

   integer, parameter :: T_PLAIN  = 0, &        !< enum for message types
        &                T_ERR    = T_PLAIN + 1, &
        &                T_WARN   = T_ERR   + 1, &
        &                T_INFO   = T_WARN  + 1, &
        &                T_SILENT = T_INFO  + 1

   integer, parameter :: PIERNIK_START       = 1                       ! before initialization
   integer, parameter :: PIERNIK_INITIALIZED = PIERNIK_START       + 1 ! initialized, running
   integer, parameter :: PIERNIK_FINISHED    = PIERNIK_INITIALIZED + 1
   integer, parameter :: PIERNIK_CLEANUP     = PIERNIK_FINISHED    + 1
   integer            :: code_progress          !< rough estimate of code execution progress

   logical            :: skip_advection = .false. !< .true. will instruct fluidupdate:make_3sweeps to skip sweeps (used by maclaurin problem, replaces precompiler symbol __NO_FLUID_STEP)
   logical, save      :: dataio_initialized = .false.

   character(LEN=1024) :: msg                   !< buffer for messages

   include 'mpif.h'

contains

!-----------------------------------------------------------------------------

   subroutine colormessage(nm, mode)

      implicit none

      character(len=*), intent(in) :: nm
      integer, intent(in) :: mode

      integer, parameter  :: log_lun = 3            !< luncher for log file

      character(len=4), parameter :: ansi_black  = char(27)//"[0m"
      character(len=7), parameter :: ansi_red    = char(27)//'[1;31m'
      character(len=7), parameter :: ansi_green  = char(27)//'[1;32m'
      character(len=7), parameter :: ansi_yellow = char(27)//'[1;33m'

      character(len=7)  :: ansicolor
      character(len=7) :: msg_type_str
      integer :: proc, ierr

      select case(mode)
         case (T_ERR)
            ansicolor = ansi_red
            msg_type_str = "Error  "
         case (T_WARN)
            ansicolor = ansi_yellow
            msg_type_str = "Warning"
         case (T_INFO)
            ansicolor = ansi_green
            msg_type_str = "Info   "
         case (T_SILENT)
            ansicolor = ansi_black
            msg_type_str = ""
         case default ! T_PLAIN
            ansicolor = ansi_black
            msg_type_str = ""
      end select

      call MPI_comm_rank(MPI_COMM_WORLD, proc, ierr)

      if (mode /= T_SILENT) then
         if (mode == T_PLAIN) then
            write(*,'(a)') trim(nm)                                                                               ! QA_WARN
         else
            write(*,'(a,i5,2a)') trim(ansicolor)//msg_type_str//" @"//ansi_black, proc, ': ', trim(nm)            ! QA_WARN
         end if
      end if

      if (dataio_initialized) then
         open(log_lun, file=log_file, position='append')
         if (proc == 0 .and. mode == T_ERR) write(log_lun,'(/,a,/)')"###############     Crashing     ###############"
         write(log_lun,'(a,i5,2a)') msg_type_str//" @", proc, ': ', trim(nm)
         close(log_lun)
      else
         write(*,'(a,i5,a)') ansi_yellow//"[errh:colormessage] dataio_initialized == .false. @", proc, ansi_black ! QA_WARN
      endif

   end subroutine colormessage

   subroutine get_container(nstep)

      implicit none

      integer, intent(out) :: nstep

      nstep         = chdf%nstep
      nhdf          = chdf%nhdf
      ntsl          = chdf%ntsl
      nres          = chdf%nres
      nlog          = chdf%nlog
      step_hdf      = chdf%step_hdf
      last_hdf_time = chdf%last_hdf_time
      log_file      = chdf%log_file
      nrestart      = chdf%nrestart
      domain        = chdf%domain

   end subroutine get_container

   subroutine set_container_chdf(nstep)

      implicit none

      integer, intent(in) ::  nstep

      chdf%nstep          = nstep
      chdf%nhdf           = nhdf
      chdf%ntsl           = ntsl
      chdf%nres           = nres
      chdf%nlog           = nlog
      chdf%step_hdf       = step_hdf
      chdf%last_hdf_time  = last_hdf_time
      chdf%log_file       = log_file
      chdf%nrestart       = nrestart
      chdf%domain         = domain

   end subroutine set_container_chdf

end module dataio_public
