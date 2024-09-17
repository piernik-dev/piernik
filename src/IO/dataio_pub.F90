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

   use constants, only: cbuff_len, domlen, idlen, cwdlen, V_INFO

   implicit none

   public  ! QA_WARN most variables are not secrets here
   private :: err_mpi, colormessage, T_PLAIN, T_ERR, T_WARN, T_INFO, T_IO, T_SILENT, & ! QA_WARN no need to use these symbols outside dataio_pub
        &     ansi_red, ansi_green, ansi_yellow, ansi_blue, ansi_magenta, ansi_cyan, & ! QA_WARN
        &     namelist_handler_t, logbuffer                                            ! QA_WARN
   private :: cbuff_len, domlen, idlen, cwdlen, V_INFO ! QA_WARN prevent re-exporting
   !mpisetup uses: ansi_white and ansi_black

   real, parameter             :: piernik_hdf5_version = 1.19    !< output version

   ! v2 specific
   real, parameter             :: piernik_hdf5_version2 = 2.03   !< output version for multi-file, multi-domain I/O
   logical                     :: use_v2_io                      !< prefer the new I/O format
   logical                     :: gdf_strict                     !< adhere more strictly to GDF standard
   integer(kind=4)             :: nproc_io                       !< how many processes do the I/O (v2 only)
   logical                     :: can_i_write                    !< .true. for processes allowed to write
   logical                     :: enable_compression             !< set to .true. to enable automatic compression (test I/O performance before use, avoid on serial I/O)
   integer(kind=4)             :: gzip_level                     !< gzip compression strength: 1 - lowest and fast, 9 - best and slow
   logical                     :: h5_64bit                       !< single or double precision plotfiles

   ! Buffer lengths used only in I/O routines
   integer, parameter          :: msglen = 1024                  !< 1kB for a message ought to be enough for anybody ;-)
   integer, parameter          :: ansilen = 7, ansirst=4         !< length of our ANSI color-strings

   ! Simulation control
   character(len=cbuff_len)    :: problem_name                   !< The problem name
   character(len=idlen)        :: run_id                         !< Auxiliary run identifier
   character(len=idlen)        :: res_id                         !< Auxiliary run to restart identifier, yet different then current run_id of restarted simulation (e.g. to avoid overwriting of the output from the previous (pre-restart) simulation; if res_id = '' then run_id is used also for reading restart file)
   real                        :: tend                           !< simulation time to end
   real                        :: wend                           !< wall clock time to end (in hours)

   integer, target             :: lun                            !< current free logical unit
   integer                     :: tsl_lun                        !< logical unit number for timeslice file
   integer                     :: log_lun                        !< logical unit number for log file
   integer(kind=4)             :: nend                           !< number of the step to end simulation
   integer(kind=4), save       :: cbline = 0                     !< current buffer line
   integer                     :: nstep_start                    !< number of start timestep
   integer(kind=4)             :: nhdf                           !< current number of hdf file
   integer(kind=4)             :: nres                           !< current number of restart file
   integer(kind=4)             :: nrestart                       !< number of restart file to be read while restart is not set to ''
   character(len=domlen)       :: domain_dump                    !< string to choose if boundaries have to be dumped in hdf files

   ! Buffers for global use
   character(len=cwdlen), save :: wd_rd = "./"                   !< path to problem.par and/or restarts
   character(len=cwdlen), save :: wd_wr = "./"                   !< path where output is written
   character(len=cwdlen), save :: log_wr = "./"                  !< path where log is written
   character(len=msglen), save, target :: cmdl_nml =" "          !< buffer for namelist supplied via commandline
   character(len=cwdlen)       :: log_file                       !< path to the current log file
   character(len=cwdlen)       :: tsl_file                       !< path to the current tsl file
   character(len=cwdlen)       :: tmp_log_file                   !< path to the temporary log file
   character(len=cwdlen), target :: par_file                     !< path to the parameter file
   ! Handy variables
   integer(kind=4), target     :: ierrh                          !< variable for iostat error on reading namelists
   integer(kind=4)             :: err_mpi                        !< variable for error code in MPI calls (should we export it to mpisetup?)
   character(len=cwdlen), target :: errstr                       !< string for storing error messages

   real                        :: last_log_time                  !< time in simulation of the recent dump of statistics into a log file
   real                        :: last_tsl_time                  !< time in simulation of the recent timeslice dump
   real                        :: last_hdf_time                  !< time in simulation of the recent hdf file dump
   real                        :: last_res_time                  !< time in simulation of the recent res file dump
   integer                     :: code_progress                  !< rough estimate of code execution progress

   ! storage for the problem.par
   integer, parameter          :: maxparfilelen   = 500          !< max length of line in problem.par file
   integer, parameter          :: maxparfilelines = 256          !< max number of lines in problem.par
   integer(kind=4), parameter  :: bufferlines = 128              !< max number of lines in the log buffer
   character(len=maxparfilelen), allocatable, dimension(:) :: parfile !< contents of the parameter file
   character(len=msglen), allocatable, dimension(:) :: logbuffer !< buffer for log I/O
   integer, save               :: parfilelines = 0               !< number of lines in the parameter file
   integer(kind=4)             :: maxparlen                      !< max length of parfile lines
   integer(kind=4)             :: maxenvlen                      !< max length of env array

   logical, save               :: halfstep = .false.             !< true when X-Y-Z sweeps are done and Z-Y-X are not
   logical, save               :: log_file_initialized = .false. !< logical to mark initialization of logfile
   logical, save               :: log_file_opened = .false.      !< logical to mark opening of logfile
   logical, save               :: restarted_sim = .false.        !< logical to distinguish between new and restarted simulation
   integer(kind=4), save       :: require_problem_IC = 0         !< 1 will call initproblem::problem_initial_conditions on restart

   logical                     :: multiple_h5files = .false.     !< write one HDF5 file per proc

   !< enum for message types
   ! stdout and logfile messages
   enum, bind(C)
      enumerator ::  T_PLAIN = 0
      enumerator :: T_PLAIN_NOA, T_ERR, T_WARN, T_INFO, T_IO, T_IO_NOA, T_SILENT
   end enum
   character(len=msglen)       :: msg                            !< buffer for messages
   character(len=ansirst)      :: ansi_black
   character(len=ansilen)      :: ansi_red, ansi_green, ansi_yellow, ansi_blue, ansi_magenta, ansi_cyan, ansi_white
   real                        :: thdf                           !< hdf dump wallclock

   integer(kind=4)             :: piernik_verbosity = V_INFO

   ! Per suggestion of ZEUS sysops:
   ! http://www.fz-juelich.de/ias/jsc/EN/Expertise/Supercomputers/JUROPA/UserInfo/IO_Tuning.htm
#if defined(__INTEL_COMPILER)
   integer, parameter               :: io_par = 4
   character(len=io_par), parameter :: io_buffered = "yes"
   integer, parameter               :: io_blocksize = 1048576
   integer, parameter               :: io_buffno = 1
#endif /* __INTEL_COMPILER */

   interface
      subroutine namelist_errh_P(ierrh, nm, skip_eof)
         implicit none
         integer(kind=4),   intent(in) :: ierrh
         character(len=*),  intent(in) :: nm
         logical, optional, intent(in) :: skip_eof
      end subroutine namelist_errh_P

      subroutine compare_namelist_P(nml_bef, nml_aft)
         implicit none
         character(len=*), intent(in)     :: nml_bef, nml_aft
      end subroutine compare_namelist_P
   end interface

   type :: namelist_handler_t
      character(len=msglen), pointer :: cmdl_nml   !< buffer for namelist supplied via commandline
      character(len=cwdlen), pointer :: par_file   !< path to the parameter file
      character(len=cwdlen), pointer :: errstr     !< string for storing error messages
      character(len=cwdlen) :: tmp1
      character(len=cwdlen) :: tmp2
      integer(kind=4), pointer       :: ierrh      !< variable for iostat error on reading namelists
      integer, pointer               :: lun        !< current free logical unit
      procedure(namelist_errh_P), nopass, pointer    :: namelist_errh
      logical :: initialized = .false.
   contains
      procedure :: init => namelist_handler_t_init
      procedure :: compare_namelist
   end type namelist_handler_t

   type(namelist_handler_t) :: nh

   interface printinfo
      module procedure printinfo_legacy
      module procedure printinfo_v
   end interface printinfo

contains

   subroutine namelist_handler_t_init(this)

      use constants, only: cbuff_len
#if defined(__INTEL_COMPILER)
      use ifport,    only: getpid
#endif /* __INTEL_COMPILER */

      implicit none

      class(namelist_handler_t), intent(inout) :: this

      character(len=cwdlen) :: tmpdir
      character(len=cbuff_len) :: pid
      integer :: lchar_tmpdir

      call get_environment_variable("PIERNIK_TMPDIR", tmpdir)
      lchar_tmpdir = len_trim(tmpdir)
      if (lchar_tmpdir == 0) then
         tmpdir = "."
         lchar_tmpdir = 1
      else
         if (tmpdir(lchar_tmpdir:lchar_tmpdir) == '/') lchar_tmpdir = lchar_tmpdir - 1
      endif

      write(pid, '(i5)') getpid()

      ! The filenames here will be unique, unless multiple Piernik instances
      ! from multiple machines are started in the same shared directory.
      write(this%tmp1, '(a,"/temp1_' // trim(pid) // '.dat")') tmpdir(1:lchar_tmpdir)
      write(this%tmp2, '(a,"/temp2_' // trim(pid) // '.dat")') tmpdir(1:lchar_tmpdir)

      this%cmdl_nml => cmdl_nml
      this%par_file => par_file
      this%errstr => errstr
      this%ierrh => ierrh
      this%lun => lun

      this%namelist_errh => namelist_errh

      this%initialized = .true.
   end subroutine namelist_handler_t_init
!-----------------------------------------------------------------------------
   subroutine set_colors(enable)

      implicit none

      logical, intent(in) :: enable

      if (enable) then
         write(ansi_black,  '(A1,A3)') char(27),"[0m"
         write(ansi_red,    '(A1,A6)') char(27),"[1;31m"
         write(ansi_green,  '(A1,A6)') char(27),"[1;32m"
         write(ansi_yellow, '(A1,A6)') char(27),"[1;33m"
         write(ansi_blue,   '(A1,A6)') char(27),"[1;34m"
         write(ansi_magenta,'(A1,A6)') char(27),"[1;35m"
         write(ansi_cyan,   '(A1,A6)') char(27),"[1;36m"
         write(ansi_white,  '(A1,A6)') char(27),"[1;37m"
      else
         ansi_black = ''
         ansi_red = ''
         ansi_green = ''
         ansi_yellow = ''
         ansi_blue = ''
         ansi_magenta = ''
         ansi_cyan = ''
         ansi_white = ''
      endif

   end subroutine set_colors
!-----------------------------------------------------------------------------
   subroutine colormessage(nm, mode)

      use constants, only: stdout, stderr, idlen, I_ONE
      use MPIF,      only: MPI_COMM_WORLD, MPI_Comm_rank

      implicit none

      character(len=*),  intent(in) :: nm
      integer(kind=4),   intent(in) :: mode

      character(len=ansilen)        :: ansicolor
      integer, parameter            :: msg_type_len = len("Warning")
      character(len=msg_type_len)   :: msg_type_str
      integer(kind=4)               :: proc
      integer                       :: outunit
      character(len=idlen)          :: adv

!      write(stdout,*) ansi_red, "Red ", ansi_green, "Green ", ansi_yellow, "Yellow ", ansi_blue, "Blue ", ansi_magenta, "Magenta ", ansi_cyan, "Cyan ", ansi_white, "White ", ansi_black
      adv = 'yes'
      if ((mode == T_IO_NOA) .or. (mode == T_PLAIN_NOA)) adv = 'no'
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
         case default ! T_PLAIN, T_PLAIN_NOA
            ansicolor = ansi_black
            outunit   = stdout
            msg_type_str = ''
      end select

      call MPI_Comm_rank(MPI_COMM_WORLD, proc, err_mpi)

      if (mode /= T_SILENT) then
         if ((mode == T_PLAIN) .or. (mode == T_PLAIN_NOA)) then
            write(outunit,'(a)') trim(nm)
         else
            write(outunit,'(a,a," @",a,i5,2a)', advance=adv) trim(ansicolor), msg_type_str, ansi_black, proc, ': ', trim(nm)
         endif
      endif

      ! TODO: get msgs from other procs if necessary
      if (proc == 0) then
         if (log_file_initialized) then
            if (.not.log_file_opened) then
#if defined(__INTEL_COMPILER)
               open(newunit=log_lun, file=log_file, position='append', blocksize=io_blocksize, buffered=io_buffered, buffercount=io_buffno)
#else /* __INTEL_COMPILER */
               open(newunit=log_lun, file=log_file, position='append') !> \todo reconstruct asynchronous writing to log files
#endif /* !__INTEL_COMPILER */
               log_file_opened = .true.
            endif
         else
            ! BEWARE: possible race condition
#if defined(__INTEL_COMPILER)
            open(newunit=log_lun, file=tmp_log_file, status='unknown', position='append', blocksize=io_blocksize, buffered=io_buffered, buffercount=io_buffno)
#else /* __INTEL_COMPILER */
            open(newunit=log_lun, file=tmp_log_file, status='unknown', position='append')
#endif /* !__INTEL_COMPILER */
         endif
         if (proc == 0 .and. mode == T_ERR) write(log_lun,'(/,a,/)')"###############     Crashing     ###############"
         call allocate_text_buffers
         if (cbline < size(logbuffer)) then
            cbline = cbline + I_ONE
            write(logbuffer(cbline), '(2a,i5,2a)') msg_type_str," @", proc, ': ', trim(nm)
         endif
         if (cbline >= size(logbuffer)) call flush_to_log
         if (mode == T_ERR) call flush_to_log
         if (.not. log_file_initialized) close(log_lun)
      else
         if (mode == T_SILENT) &
            write(stderr,'(a,a," @",a,i5,2a)', advance=adv) trim(ansi_red), "not logged", ansi_black, proc, ': ', trim(nm)
      endif

   end subroutine colormessage
!-----------------------------------------------------------------------------
   subroutine flush_to_log
      implicit none
      integer :: line

      if (.not. allocated(logbuffer)) return

      do line = 1, min(cbline, size(logbuffer, kind=4))
         write(log_lun, '(a)') trim(logbuffer(line)) !> \todo reconstruct asynchronous writing to log files
      enddo
      cbline = 0

   end subroutine flush_to_log
!-----------------------------------------------------------------------------
   subroutine allocate_text_buffers

      implicit none

      if (.not. allocated(logbuffer)) allocate(logbuffer(bufferlines))
      if (.not. allocated(parfile)) allocate(parfile(maxparfilelines))

   end subroutine allocate_text_buffers
!-----------------------------------------------------------------------------
   subroutine cleanup_text_buffers

      implicit none

      if (cbline > 0) call flush_to_log

      if (allocated(logbuffer)) deallocate(logbuffer)
      if (allocated(parfile)) deallocate(parfile)

   end subroutine cleanup_text_buffers

   !> \brief Printinfo variant which prints a line of repeated characters

   subroutine print_char_line(c, verbosity)

      use constants, only: I_ONE, rep_len, V_LOG

      implicit none

      character(len=I_ONE),      intent(in) :: c          !< the character to repeat
      integer(kind=4), optional, intent(in) :: verbosity  !< verbosity level

      integer(kind=4) :: v

      v = V_LOG
      if (present(verbosity)) v = verbosity

      call printinfo(repeat(c, rep_len), v)

   end subroutine print_char_line

   !> \brief Printinfo variant which does not print colored tag (legacy)

   subroutine printinfo_legacy(nm, to_stdout)

      use constants, only: cbuff_len, V_LOG, V_DEBUG, V_INFO

      implicit none

      character(len=*), intent(in) :: nm         !< message
      logical,          intent(in) :: to_stdout  !< to print or only to log

      character(len=cbuff_len) :: l

      l = ""
      if (piernik_verbosity <= V_DEBUG) l = "<Legacy>"

!      call colormessage(nm, merge(T_PLAIN, T_SILENT, to_stdout))
      call printinfo_v(trim(l) // nm, merge(V_INFO, V_LOG, to_stdout), .true.)

   end subroutine printinfo_legacy

   !> \brief New printinfo implementation that does the filtering wrt. verbosity level

   subroutine printinfo_v(nm, verbosity, no_tag)

      use constants, only: V_DEBUG, V_INFO, V_WARN, V_LOWEST, V_HIGHEST

      implicit none

      character(len=*),          intent(in) :: nm         !< message
      integer(kind=4), optional, intent(in) :: verbosity  !< verbosity level
      logical,         optional, intent(in) :: no_tag     !< use plain formatting (no color tag, ignored for V_WARN)

      integer :: v, i
      logical :: plain
      character(len=V_HIGHEST-V_LOWEST+1) :: v_lev

      v = V_INFO
      if (present(verbosity)) v = verbosity

      if (piernik_verbosity <= V_DEBUG) then  ! add verbosity level info to the messages
         do i = V_LOWEST, V_HIGHEST-1
            v_lev(i-V_LOWEST+1:i-V_LOWEST+1) = merge(merge("+", ".", present(verbosity)), " ", i < v)
         enddo
         v_lev(V_HIGHEST-V_LOWEST+1:V_HIGHEST-V_LOWEST+1) = ":"
      else
         v_lev = ""
      endif

      if (v >= V_WARN) then  ! now it is also possible to implement several warning levels, if necessary
         call warn(nm)
      else                   ! emit standard green-tagged message if verbosity level allows, print everything to the logile
         plain = .false.
         if (present(no_tag)) plain = no_tag
         call colormessage(trim(v_lev) // nm, merge(merge(T_PLAIN, T_INFO, plain), T_SILENT, v >= piernik_verbosity))
      endif

   end subroutine printinfo_v
!-----------------------------------------------------------------------------
   subroutine printio(nm, noadvance)

      implicit none

      character(len=*),  intent(in) :: nm
      logical, optional, intent(in) :: noadvance

      logical :: adv

      adv  = .true.
      if (present(noadvance)) adv = .not. noadvance

      call colormessage(nm, merge(T_IO, T_IO_NOA, adv))

   end subroutine printio
!-----------------------------------------------------------------------------
   subroutine print_plain(nm, noadvance)

      implicit none

      character(len=*), intent(in) :: nm
      logical,          intent(in) :: noadvance

      call colormessage(nm, merge(T_PLAIN_NOA, T_PLAIN, noadvance))

   end subroutine print_plain
!-----------------------------------------------------------------------------
   subroutine warn(nm)

      implicit none

      character(len=*), intent(in) :: nm

      call colormessage(nm, T_WARN)

   end subroutine warn
!-----------------------------------------------------------------------------
   !> \deprecated BEWARE: routine is not finished, it should kill PIERNIK gracefully
   subroutine die(nm, allprocs)

      use MPIF,   only: MPI_COMM_WORLD, MPI_Barrier, MPI_Finalize
#if defined(__INTEL_COMPILER)
      use ifcore, only: tracebackqq
#endif /* __INTEL_COMPILER */

      implicit none

      character(len=*),  intent(in) :: nm
      integer, optional, intent(in) :: allprocs

      call colormessage(nm, T_ERR)

      if (present(allprocs)) then
         if (allprocs /= 0) then
            call MPI_Barrier(MPI_COMM_WORLD, err_mpi)
            call MPI_Finalize(err_mpi)
         endif
      endif
      call colormessage("Following backtrace is used for debugging, please attach it to your bug report", T_ERR)
#if defined(__INTEL_COMPILER)
      call colormessage("Be advised, that you need to compile the code with -g -traceback for the dump to be meaningful", T_ERR)
      call tracebackqq()
#else /* !__INTEL_COMPILER */
      call colormessage("Be advised, that you need to compile the code with -g for the dump to be meaningful", T_ERR)
      call abort()
#endif /* !__INTEL_COMPILER */
      call exit(-1)

   end subroutine die
!-----------------------------------------------------------------------------
   subroutine namelist_errh(ierrh,nm,skip_eof)

      use constants, only: V_VERBOSE

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
            call warn(msg)
            call die(errstr)
         case (-1)
            if (present(skip_eof)) then
               if (skip_eof) return
            endif
            write(msg,'(3a)') "Namelist: ",trim(nm)," not found in problem.par. Assuming defaults." ! Can happen also when there is no EOL after the namelist in problem.par
            call printinfo(msg, V_VERBOSE)
         case (239, 5010)
            write(msg,'(3a)') "A problem with one of the variables that belong to the ",trim(nm)," namelist was found in problem.par:"
            call warn(msg)
            call die(errstr)
         case (0)
         case default
            write(msg, *)'Unknown error (', ierrh,') in namelist ',trim(nm)
            call warn(msg)
            if (present(skip_eof)) then
               if (skip_eof) return
            endif
            call die(errstr)
      endselect

   end subroutine namelist_errh
!-----------------------------------------------------------------------------
   subroutine compare_namelist(this)

      use constants, only: PIERNIK_INIT_IO_IC, V_LOG
      use MPIF,      only: MPI_COMM_WORLD, MPI_Comm_rank

      implicit none
      class(namelist_handler_t), intent(inout) :: this
      integer                          :: io
      character(len=maxparfilelen)     :: sa, sb
      integer                          :: lun_bef, lun_aft
      integer(kind=4)                  :: proc

      call MPI_Comm_rank(MPI_COMM_WORLD, proc, err_mpi)
      if (proc > 0) call die("[dataio_pub:compare_namelist] This routine must not be called by many threads at once. Make sure that diff_nml macro is called only from rank 0.")

      if (code_progress > PIERNIK_INIT_IO_IC) call warn("[dataio_pub:compare_namelist] Late namelist")

      open(newunit=lun_bef, file=this%tmp1, status='old')
      open(newunit=lun_aft, file=this%tmp2, status='old')
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
         call printinfo(msg, V_LOG)
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
   subroutine close_txt_file(lfile, llun)

      implicit none

      integer,               intent(in) :: llun     !< logical unit number for txt file
      character(len=cwdlen), intent(in) :: lfile    !< path to txt file

      logical :: lopen

      inquire(file=lfile, opened=lopen)
      if (lopen) then
         flush(llun)
         close(llun)
      endif

   end subroutine close_txt_file

   subroutine close_logs

      use MPIF, only: MPI_COMM_WORLD, MPI_Comm_rank

      implicit none

      integer(kind=4) :: proc

      call MPI_Comm_rank(MPI_COMM_WORLD, proc, err_mpi)

      if (proc == 0) then
         call flush_to_log
         call close_txt_file(log_file, log_lun)
         call close_txt_file(tsl_file, tsl_lun)
      endif

   end subroutine close_logs
!>
!! \brief Sanitize a file name
!!
!! \warning This routine is used only after reading the restart file.
!! This may result in inconsistent file naming, when problem_name or run_id contains characters such as '+'.
!<

   function fix_string(str) result (outstr)

      implicit none

      character(len=*), intent(in) :: str

      character(len=len(str)) :: outstr

      integer :: i

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
