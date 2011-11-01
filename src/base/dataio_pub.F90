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
!>
!! \brief module providing common dataio parameters, variables, routines and interfaces
!<
module dataio_pub

   use constants, only: cbuff_len, domlen, idlen, cwdlen

   implicit none

   public  ! QA_WARN most variables are not secrets here
   private :: colormessage, T_PLAIN, T_ERR, T_WARN, T_INFO, T_IO, T_SILENT, ansi_red, ansi_green, ansi_yellow, ansi_blue, ansi_magenta, ansi_cyan ! QA_WARN no need to use these symbols outside dataio_pub
   private :: cbuff_len, domlen, idlen, cwdlen ! QA_WARN prevent re-exporting
   !mpisetup uses: ansi_white and ansi_black

   real, parameter :: piernik_hdf5_version = 1.17   !< output version
   real, parameter :: piernik_hdf5_version2 = 2.0   !< output version for multi-file, multi-domain I/O
   logical         :: use_v2_io                     !< prefer the new I/O format
   integer(kind=4) :: nproc_io                      !< how many processes do the I/O (v2 only)
   logical         :: can_i_write                   !< .true. for processes allowed to write

   ! Buffer lengths used only in I/O routines
   integer, parameter :: msglen = 1024          !< 1kB for a message ought to be enough for anybody ;-)
   integer, parameter :: ansilen = 7, ansirst=4 !< length of our ANSI color-strings

   ! Simulation control
   character(len=cbuff_len) :: problem_name     !< The problem name
   character(len=idlen)     :: run_id           !< Auxiliary run identifier
   character(len=idlen)     :: new_id           !< COMMENT ME
   real                  :: tend                   !< simulation time to end
   real                  :: wend                   !< wall clock time to end (in hours)

   integer               :: lun                    !< current free logical unit
   integer               :: nend                   !< number of the step to end simulation
   integer               :: nstep_start            !< number of start timestep
   integer(kind=4)       :: nhdf                   !< current number of hdf file
   integer(kind=4)       :: nres                   !< current number of restart file
   real                  :: next_t_log             !< when to print statistics to the log file
   real                  :: next_t_tsl             !< when to produce the timeslice file
   integer(kind=4)       :: nrestart               !< number of restart file to be read while restart is not set to ''
   integer (kind=4)      :: step_hdf               !< number of simulation timestep corresponding to values dumped in hdf file
   integer(kind=4)       :: step_res               !< COMMENT ME
   character(len=domlen) :: domain_dump         !< string to choose if boundaries have to be dumped in hdf files

   ! Buffers for global use
   character(len=cwdlen), save :: cwd = "."     !< path to the current working directory
   character(len=msglen), save :: cmdl_nml =" " !< buffer for namelist supplied via commandline
   character(len=cwdlen) :: log_file            !< path to the current log file
   character(len=cwdlen) :: tmp_log_file        !< path to the temporary log file
   character(len=cwdlen) :: par_file            !< path to the parameter file
   ! Handy variables
   integer(kind=4)       :: ierrh                  !< variable for iostat

   real               :: last_hdf_time                  !< time in simulation of the last resent hdf file dump
   integer            :: code_progress                  !< rough estimate of code execution progress

   ! storage for the problem.par
   integer, parameter :: maxparfilelen   = 128                         !< max length of line in problem.par file
   integer, parameter :: maxparfilelines = 256                         !< max number of lines in problem.par
   character(len=maxparfilelen), dimension(maxparfilelines) :: parfile !< contents of the parameter file
   integer, save :: parfilelines = 0                                   !< number of lines in the parameter file

   logical, save      :: halfstep = .false.             !< true when X-Y-Z sweeps are done and Z-Y-X are not
   logical, save      :: log_file_initialized = .false. !< logical to mark initialization of logfile
   integer(kind=4), save :: require_init_prob = 0       !< 1 will call initproblem::init_prob on restart
!>
!! \todo
!!  Currently to use PGPLOT you need to:
!!   1. set one of i{x,y,z} to positive value and zero to the others
!!   2. use only one-element vars array in problem.par
!<
   logical            :: vizit = .false.        !< perform "live" visualization using pgplot (\deprecated BEWARE: highly experimental)
   logical            :: multiple_h5files = .false. !< write one HDF5 file per proc
   real               :: fmin                   !< minimum on pgplot scale
   real               :: fmax                   !< maximum on pgplot scale

   ! stdout and logfile messages
   integer, parameter :: T_PLAIN  = 0, &        !< enum for message types
        &                T_ERR    = T_PLAIN  + 1, &
        &                T_WARN   = T_ERR    + 1, &
        &                T_INFO   = T_WARN   + 1, &
        &                T_IO     = T_INFO   + 1, &
        &                T_IO_NOA = T_IO     + 1, &
        &                T_SILENT = T_IO_NOA + 1
   character(len=msglen)  :: msg                !< buffer for messages
   character(len=ansirst) :: ansi_black
   character(len=ansilen) :: ansi_red, ansi_green, ansi_yellow, ansi_blue, ansi_magenta, ansi_cyan, ansi_white
   character(len=*),parameter :: tmr_hdf = "hdf_dump"
   real                       :: thdf !< hdf dump wallclock

contains
!-----------------------------------------------------------------------------
   subroutine colormessage(nm, mode)

      use constants, only: stdout, stderr, idlen
      use mpi,       only: MPI_COMM_WORLD

      implicit none

      character(len=*),  intent(in) :: nm
      integer,           intent(in) :: mode

      integer                       :: log_lun                !< luncher for log file
      character(len=ansilen)        :: ansicolor
      integer, parameter            :: msg_type_len = 7       !< length of the "Warning" word.
      character(len=msg_type_len)   :: msg_type_str
      integer(kind=4)               :: proc
      integer                       :: outunit
      logical, save                 :: frun = .true.
      character(len=idlen)          :: adv

      if (frun) then
         write(ansi_black,  '(A1,A3)') char(27),"[0m"
         write(ansi_red,    '(A1,A6)') char(27),"[1;31m"
         write(ansi_green,  '(A1,A6)') char(27),"[1;32m"
         write(ansi_yellow, '(A1,A6)') char(27),"[1;33m"
         write(ansi_blue,   '(A1,A6)') char(27),"[1;34m"
         write(ansi_magenta,'(A1,A6)') char(27),"[1;35m"
         write(ansi_cyan,   '(A1,A6)') char(27),"[1;36m"
         write(ansi_white,  '(A1,A6)') char(27),"[1;37m"
         frun = .false.
      endif

!      write(stdout,*) ansi_red, "Red ", ansi_green, "Green ", ansi_yellow, "Yellow ", ansi_blue, "Blue ", ansi_magenta, "Magenta ", ansi_cyan, "Cyan ", ansi_white, "White ", ansi_black
      adv = 'yes'
      if (mode == T_IO_NOA) adv = 'no'
      select case (mode)
         case (T_ERR)
            ansicolor = ansi_red
            outunit   = stderr
            msg_type_str = "Error  "
         case (T_WARN)
            ansicolor = ansi_yellow
            outunit   = stdout
            msg_type_str = "Warning"
         case (T_INFO)
            ansicolor = ansi_green
            outunit   = stdout
            msg_type_str = "Info   "
         case (T_IO, T_IO_NOA)
            ansicolor = ansi_blue
            outunit   = stdout
            msg_type_str = "I/O    "
         case (T_SILENT)
            ansicolor = ansi_black
            outunit   = stdout
            msg_type_str = ''
         case default ! T_PLAIN
            ansicolor = ansi_black
            outunit   = stdout
            msg_type_str = ''
      end select

      call MPI_Comm_rank(MPI_COMM_WORLD, proc, ierrh)

      if (mode /= T_SILENT) then
         if (mode == T_PLAIN) then
            write(outunit,'(a)') trim(nm)
         else
            write(outunit,'(a,a," @",a,i5,2a)', advance=adv) trim(ansicolor),msg_type_str,ansi_black, proc, ': ', trim(nm)
         endif
      endif

      if (log_file_initialized) then
         open(newunit=log_lun, file=log_file, position='append')
      else
         open(newunit=log_lun, file=tmp_log_file, status='unknown', position='append')   ! BEWARE: possible race condition
      endif
      if (proc == 0 .and. mode == T_ERR) write(log_lun,'(/,a,/)')"###############     Crashing     ###############"
      write(log_lun,'(2a,i5,2a)') msg_type_str," @", proc, ': ', trim(nm)
      close(log_lun)

   end subroutine colormessage
!-----------------------------------------------------------------------------
   subroutine printinfo(nm, to_stdout)

      implicit none

      character(len=*),  intent(in) :: nm
      logical, optional, intent(in) :: to_stdout

      if (present(to_stdout)) then
         if (to_stdout) then
            call colormessage(nm, T_PLAIN)
         else
            call colormessage(nm, T_SILENT)
         endif
      else
         call colormessage(nm, T_INFO)
      endif

   end subroutine printinfo
!-----------------------------------------------------------------------------
   subroutine printio(nm, noadvance)

      implicit none

      character(len=*), intent(in) :: nm
      logical, optional, intent(in) :: noadvance

      logical :: adv

      adv  = .true.
      if (present(noadvance)) then
         if (noadvance) adv = .false.
      endif

      if (adv) then
         call colormessage(nm, T_IO)
      else
         call colormessage(nm, T_IO_NOA)
      endif

   end subroutine printio
!-----------------------------------------------------------------------------
   subroutine warn(nm)

      implicit none

      character(len=*), intent(in) :: nm

      call colormessage(nm, T_WARN)

   end subroutine warn
!-----------------------------------------------------------------------------
   !> \deprecated BEWARE: routine is not finished, it should kill PIERNIK gracefully
   subroutine die(nm, allprocs)

      use mpi,    only: MPI_COMM_WORLD

      implicit none

      character(len=*), intent(in)  :: nm
      integer, optional, intent(in) :: allprocs

      integer :: ierr

      call colormessage(nm, T_ERR)

      if (present(allprocs)) then
         if (allprocs /= 0) then
            call MPI_Barrier(MPI_COMM_WORLD, ierr)
            call MPI_Finalize(ierr)
         endif
      endif

      call exit(-1)

   end subroutine die
!-----------------------------------------------------------------------------
   subroutine namelist_errh(ierrh,nm,skip_eof)

      implicit none

      integer(kind=4),   intent(in) :: ierrh
      character(len=*),  intent(in) :: nm
      logical, optional, intent(in) :: skip_eof

      select case (ierrh)
         case (19)
            call warn("One of the following conditions occurred: ")
            call warn("    * The variable was not a member of the namelist group.")
            call warn("    * An attempt was made to subscript a scalar variable.")
            call warn("    * A subscript of the array variable was out-of-bounds.")
            call warn("    * An array variable was specified with too many or too few subscripts for the variable.")
            call warn("    * An attempt was made to specify a substring of a noncharacter variable or array name.")
            call warn("    * A substring specifier of the character variable was out-of-bounds.")
            call warn("    * A subscript or substring specifier of the variable was not an integer constant.")
            call warn("    * An attempt was made to specify a substring by using an unsubscripted array variable.")
            write(msg,'(3a)') "severe (19): Invalid reference to variable in the ",trim(nm)," namelist"
            call die(msg)
         case (-1)
            if (present(skip_eof)) then
               if (skip_eof) return
            endif
            write(msg,'(3a)') "Namelist: ",trim(nm)," not found in problem.par"
            call die(msg)
         case (239, 5010)
            write(msg,'(3a)') "One of the variables found in problem.par doesn't belong to the ",trim(nm)," namelist"
            call die(msg)
         case (0)
         case default
            write(msg, *)'Unknown error (', ierrh,') in namelist ',trim(nm)
            if (present(skip_eof)) then
               if (skip_eof) return
            endif
            call die(msg)
      endselect

   end subroutine namelist_errh
!-----------------------------------------------------------------------------
   subroutine compare_namelist(nml_bef, nml_aft)

      use constants, only: PIERNIK_INIT_IO_IC
      use mpi,       only: MPI_COMM_WORLD

      implicit none

      character(len=*), intent(in)     :: nml_bef, nml_aft
      integer                          :: io
      character(len=maxparfilelen)     :: sa, sb
      integer                          :: lun_bef, lun_aft
      integer(kind=4)                  :: proc

      call MPI_Comm_rank(MPI_COMM_WORLD, proc, ierrh)
      if (proc > 0) call die("[dataio_pub:compare_namelist] This routine must not be called by many threads at once. Make sure that diff_nml macro is called only from rank 0.")

      if (code_progress > PIERNIK_INIT_IO_IC) call warn("[dataio_pub:compare_namelist] Late namelist")

      open(newunit=lun_bef, file=nml_bef, status='old')
      open(newunit=lun_aft, file=nml_aft, status='old')
      io = 0
      do
         read(lun_bef,'(a)', iostat=io) sa
         read(lun_aft,'(a)', iostat=io) sb
         if (io/=0) exit
         if ((sa/=sb)) then
            write(msg,'(a1,a)') '*',trim(sb)
         else
            write(msg,'(a1,a)') ' ',trim(sb)
         endif
         call printinfo(msg, .false.)
      enddo
      close(lun_aft, status="delete")
      close(lun_bef, status="delete")

   end subroutine compare_namelist
!-----------------------------------------------------------------------------
   integer function move_file(a,b) result (stat)

      use iso_fortran_env, only: iostat_end

      implicit none

      character(len=*), intent(in) :: a, b

      integer                      :: old_u, new_u
      integer                      :: io_stat

      open(newunit=old_u, file=a, status="old")
      open(newunit=new_u, file=b, status="unknown")
      do
         read(old_u, '(a)', iostat=io_stat) msg
         if (io_stat == iostat_end) exit
         write(new_u,'(a)') trim(msg)
      enddo
      close(new_u)
      close(old_u, status="delete")

      stat = 0
   end function move_file
!-----------------------------------------------------------------------------

!>
!! \brief Sanitize a file name
!!
!! \warning This routine is used only after reading the restart file.
!! This may result in inconsistent file naming, when problem_name or run_id contains characters such as '+'.
!<

   function fix_string(str) result (outstr)

      implicit none

      character(len=*), intent(in)  :: str
      character(len=len(str)) :: outstr

      integer            :: i

      outstr = repeat(" ", len(str))

      do i=1, len(str)
         outstr(i:i) = ''
         select case (str(i:i))
            case ('a':'z', 'A':'Z', '0':'9', '_', '-') ! restrict file names to a very safe character set
               outstr(i:i) = str(i:i)
            case default
               exit
         end select
      enddo

   end function fix_string

end module dataio_pub
