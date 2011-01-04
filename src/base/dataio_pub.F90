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
   use iso_fortran_env, only: error_unit, output_unit
   use types,           only: hdf, domlen, idlen

   implicit none

   public  ! QA_WARN most variables are not secrets here
   private :: colormessage, T_PLAIN, T_ERR, T_WARN, T_INFO, T_IO, T_SILENT, ansi_red, ansi_green, ansi_yellow, ansi_blue, ansi_magenta, ansi_cyan ! QA_WARN no need to use these symbols outside dataio_pub
   !mpisetup uses: ansi_white and ansi_black

   real, parameter    :: piernik_hdf5_version = 1.12   !< output version

   integer, parameter :: stdout = output_unit
   integer, parameter :: stderr = error_unit

   ! Buffer lengths for global use
   integer, parameter :: planelen = 2           !< length of planes names e.g. "xy","yz","rp" etc.
   integer, parameter :: varlen = 4             !< length of state variable names in hdf files
   integer, parameter :: fplen = 24             !< length of buffer for printed FP or integer number
   integer, parameter :: hnlen = 32             !< hostname length limit
   integer, parameter :: cwdlen = 512           !< allow for quite long CWD
   integer, parameter :: msglen = 1024          !< 1kB for a message ought to be enough for anybody ;-)
   integer, parameter :: cbuff_len=32           !< length for problem parameters
   integer, parameter :: ansilen = 7, ansirst=4 !< length of our ANSI color-strings
   integer, parameter :: maxparfilelen   = 128  !< max length of line in problem.par file
   integer, parameter :: maxparfilelines = 256  !< max number of lines in problem.par

   ! Simulation control
   character(len=cbuff_len) :: problem_name     !< The problem name
   character(len=idlen)     :: run_id           !< Auxiliary run identifier
   real               :: tend                   !< simulation time to end
   real               :: wend                   !< wall clock time to end (in hours)

   integer            :: nend                   !< number of the step to end simulation
   integer            :: nstep_start            !< number of start timestep
   integer            :: nhdf                   !< current number of hdf file
   integer            :: nres                   !< current number of restart file
   real               :: next_t_log             !< when to print statistics to the log file
   real               :: next_t_tsl             !< when to produce the timeslice file
   integer            :: nrestart               !< number of restart file to be read while restart is not set to ''
   integer            :: step_hdf               !< number of simulation timestep corresponding to values dumped in hdf file
   character(len=domlen) :: domain              !< string to choose if boundaries have to be dumped in hdf files

   ! Buffers for global use
   character(len=cwdlen) :: cwd                 !< path to the current working directory
   character(len=cwdlen) :: log_file            !< path to the current log file
   character(len=cwdlen) :: tmp_log_file        !< path to the temporary log file
   character(len=cwdlen) :: par_file            !< path to the parameter file

   ! Handy variables
   integer            :: ierrh                  !< variable for iostat
   type(hdf)          :: chdf                   !< container for some vital simulation parameters

   ! Simulation state
   integer, parameter :: PIERNIK_START       = 1                       ! before initialization
   integer, parameter :: PIERNIK_INITIALIZED = PIERNIK_START       + 1 ! initialized, running
   integer, parameter :: PIERNIK_FINISHED    = PIERNIK_INITIALIZED + 1
   integer, parameter :: PIERNIK_CLEANUP     = PIERNIK_FINISHED    + 1

   real               :: last_hdf_time                  !< time in simulation of the last resent hdf file dump
   integer            :: code_progress                  !< rough estimate of code execution progress

   logical, save      :: halfstep = .false.             !< true when X-Y-Z sweeps are done and Z-Y-X are not
   logical, save      :: skip_advection = .false.       !< .true. will instruct fluidupdate:make_3sweeps to skip sweeps (used by maclaurin problem, replaces precompiler symbol __NO_FLUID_STEP)
   logical, save      :: log_file_initialized = .false. !< \todo Comment me
   integer, save      :: require_init_prob = 0          !< 1 will call initproblem::init_prob on restart

   !! \todo
   !!  Currently to use PGPLOT you need to:
   !!   1. set one of i{x,y,z} to positive value and zero to the others
   !!   2. use only one-element vars array in problem.par
   !!
   logical            :: vizit = .false.        !< perform "live" visualization using pgplot (BEWARE: highly experimental)
   real               :: fmin                   !< minimum on pgplot scale
   real               :: fmax                   !< maximum on pgplot scale

   ! stdout and logfile messages
   integer, parameter :: T_PLAIN  = 0, &        !< enum for message types
        &                T_ERR    = T_PLAIN + 1, &
        &                T_WARN   = T_ERR   + 1, &
        &                T_INFO   = T_WARN  + 1, &
        &                T_IO     = T_INFO  + 1, &
        &                T_SILENT = T_IO    + 1
   character(len=msglen)  :: msg                !< buffer for messages
   character(len=ansirst) :: ansi_black
   character(len=ansilen) :: ansi_red, ansi_green, ansi_yellow, ansi_blue, ansi_magenta, ansi_cyan, ansi_white

   interface
      subroutine plt_hdf5(var,ij,xn,tab,ierrh)
         implicit none
         character(LEN=*), intent(in)        :: var   !< quantity to be plotted
         character(LEN=*), intent(in)        :: ij    !< plane of plot
         integer, intent(in)                 :: xn    !< no. of cell at which we are slicing the local block
         integer, intent(inout)              :: ierrh !< error handling
         real, dimension(:,:), intent(inout) :: tab   !< array  containing given quantity
      end subroutine plt_hdf5
   end interface

   interface
      subroutine vars_hdf5(var, tab, ierrh)
         implicit none
         character(len=*), intent(in)                    :: var
         real(kind=4), dimension(:,:,:), intent(inout)   :: tab
         integer, intent(inout)                          :: ierrh
      end subroutine vars_hdf5
   end interface

   interface
      subroutine tsl_out(user_vars, tsl_names)
         implicit none
         real, dimension(:), intent(inout), allocatable                       :: user_vars
         character(len=*), dimension(:), intent(inout), allocatable, optional :: tsl_names

      end subroutine tsl_out
   end interface

   procedure(plt_hdf5),  pointer :: user_plt_hdf5 => Null()
   procedure(vars_hdf5), pointer :: user_vars_hdf5 => Null()
   procedure(tsl_out),   pointer :: user_tsl => Null()

contains
!-----------------------------------------------------------------------------
   subroutine colormessage(nm, mode)
      use mpi,       only: MPI_COMM_WORLD
      implicit none

      character(len=*),  intent(in) :: nm
      integer,           intent(in) :: mode

      integer, parameter            :: log_lun = 3            !< luncher for log file
      character(len=ansilen)        :: ansicolor
      integer, parameter            :: msg_type_len = 7       !< length of the "Warning" word.
      character(len=msg_type_len)   :: msg_type_str
      integer                       :: proc
      integer                       :: outunit
      logical, save                 :: frun = .true.

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
         case (T_IO)
            ansicolor = ansi_blue
            outunit   = stdout
            msg_type_str = "I/O    "
         case (T_SILENT)
            ansicolor = ansi_black
            outunit   = stdout
            msg_type_str = ""
         case default ! T_PLAIN
            ansicolor = ansi_black
            outunit   = stdout
            msg_type_str = ""
      end select

      call MPI_Comm_rank(MPI_COMM_WORLD, proc, ierrh)

      if (mode /= T_SILENT) then
         if (mode == T_PLAIN) then
            write(outunit,'(a)') trim(nm)
         else
            write(outunit,'(a,a," @",a,i5,2a)') trim(ansicolor),msg_type_str,ansi_black, proc, ': ', trim(nm)
         endif
      endif

      if (log_file_initialized) then
         open(log_lun, file=log_file, position='append')
      else
         open(log_lun, file=tmp_log_file, status='unknown', position='append')
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
   subroutine printio(nm)

      implicit none

      character(len=*), intent(in) :: nm

      call colormessage(nm, T_IO)

   end subroutine printio
!-----------------------------------------------------------------------------
   subroutine warn(nm)

      implicit none

      character(len=*), intent(in) :: nm

      call colormessage(nm, T_WARN)

   end subroutine warn
!-----------------------------------------------------------------------------
   !! BEWARE: routine is not finished, it should kill PIERNIK gracefully
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
   subroutine get_container(nstep)

      implicit none

      integer, intent(out) :: nstep

      nstep         = chdf%nstep
      nhdf          = chdf%nhdf
      next_t_tsl    = chdf%next_t_tsl
      nres          = chdf%nres
      next_t_log    = chdf%next_t_log
      step_hdf      = chdf%step_hdf
      last_hdf_time = chdf%last_hdf_time
      nrestart      = chdf%nrestart
      domain        = chdf%domain

   end subroutine get_container
!-----------------------------------------------------------------------------
   subroutine set_container_chdf(nstep)

      implicit none

      integer, intent(in) ::  nstep

      chdf%nstep          = nstep
      chdf%nhdf           = nhdf
      chdf%next_t_tsl     = next_t_tsl
      chdf%nres           = nres
      chdf%next_t_log     = next_t_log
      chdf%step_hdf       = step_hdf
      chdf%last_hdf_time  = last_hdf_time
      chdf%nrestart       = nrestart
      chdf%domain         = domain

   end subroutine set_container_chdf
!-----------------------------------------------------------------------------
   subroutine namelist_errh(ierrh,nm)

      implicit none

      integer, intent(in) :: ierrh
      character(len=*), intent(in) :: nm

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
            write(msg,'(3a)') "Namelist: ",trim(nm)," not found in problem.par"
            call die(msg)
         case (239, 5010)
            write(msg,'(3a)') "One of the variables found in problem.par doesn't belong to the ",trim(nm)," namelist"
            call die(msg)
         case (0)
         case default
            write(msg, *)'Unknown error (', ierrh,') in namelist ',trim(nm)
            call die(msg)
      endselect

   end subroutine namelist_errh
!-----------------------------------------------------------------------------
   subroutine compare_namelist(nml_bef, nml_aft)
      use mpi,    only: MPI_COMM_WORLD
      implicit none

      character(len=*), intent(in)     :: nml_bef, nml_aft
      integer                          :: io
      character(len=maxparfilelen)     :: sa, sb
      integer, parameter               :: lun_bef=501, lun_aft=502
      integer                          :: proc

      call MPI_Comm_rank(MPI_COMM_WORLD, proc, ierrh)
      if (proc > 0) call die("[dataio_pub:compare_namelist] This routine must not be called by many threads at once. Make sure that diff_nml macro is called only from rank 0.")

      open(lun_bef, file=nml_bef, status='old')
      open(lun_aft, file=nml_aft, status='old')
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
end module dataio_pub
