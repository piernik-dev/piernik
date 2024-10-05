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
!! \brief MPI routines setup
!!
!! In this module following namelists of parameters are specified:
!! \copydetails mpisetup::init_mpi
!<
module mpisetup

   use constants, only: cbuff_len, INT4
   use MPIF,      only: MPI_ADDRESS_KIND
#ifdef MPIF08
   use MPIF, only: MPI_Comm
#endif /* MPIF08 */

   implicit none

   private
   public :: cleanup_mpi, init_mpi, bigbang, bigbang_shift, &
        &    buffer_dim, cbuff, ibuff, lbuff, rbuff, err_mpi, tag_ub, &
        &    master, slave, nproc, proc, FIRST, LAST, have_mpi, is_spawned, &
        &    report_to_master, report_string_to_master

   integer(kind=4), protected :: nproc          !< number of processes
   integer(kind=4), protected :: proc           !< rank of my process
   integer(kind=4), protected :: LAST           !< rank of the last process
   integer(kind=4) :: err_mpi                   !< error status
   integer(kind=INT4), parameter :: FIRST = 0   !< the rank of the master process
   real(kind=8), protected    :: bigbang        !< First result of MPI_Wtime()
   real(kind=8), protected    :: bigbang_shift  !< A correction applied to readouts of MPI_Wtime() if necessary
   real(kind=8), parameter    :: min_bigbang = 1e-6  !< Start all processes from time = 1 µs
   integer(kind=MPI_ADDRESS_KIND), protected :: tag_ub

   logical, protected :: master      !< .True. if proc == FIRST
   logical, protected :: slave       !< .True. if proc != FIRST
   logical, protected :: have_mpi    !< .True. when run on more than one processor
   logical, protected :: is_spawned  !< .True. if Piernik was run via MPI_Spawn

#ifdef MPIF08
   type(MPI_Comm), protected  :: intercomm  !< intercommunicator
#else /* !MPIF08 */
   integer(kind=4), protected :: intercomm  !< intercommunicator
#endif /* !MPIF08 */

   integer, parameter :: buffer_dim = 200                   !< size of [cilr]buff arrays used to exchange namelist parameters
   character(len=cbuff_len), dimension(buffer_dim) :: cbuff !< buffer for character parameters
   integer(kind=4),          dimension(buffer_dim) :: ibuff !< buffer for integer parameters
   real,                     dimension(buffer_dim) :: rbuff !< buffer for real parameters
   logical,                  dimension(buffer_dim) :: lbuff !< buffer for logical parameters

contains

!-----------------------------------------------------------------------------
!>
!! \brief Routine to start MPI
!<

   subroutine init_mpi

      use constants,     only: cwdlen, I_ONE, V_LOG, V_DEBUG, V_INFO
      use dataio_pub,    only: die, print_char_line, printinfo, msg, tmp_log_file, par_file, lun
      use MPIF,          only: MPI_COMM_WORLD, MPI_CHARACTER, MPI_INTEGER, MPI_COMM_NULL, MPI_TAG_UB, &
           &                   MPI_Wtime, MPI_Init, MPI_Comm_get_parent, MPI_Comm_rank, MPI_Comm_size
      use MPIFUN,        only: MPI_Gather, MPI_Comm_get_attr
      use signalhandler, only: SIGINT, register_sighandler
#if defined(__INTEL_COMPILER)
      use ifport,        only: getpid, getcwd, hostnm
#endif /* __INTEL_COMPILER */
#ifndef MPIF08
      use dataio_pub,    only: warn
#endif /* !MPIF08 */

      implicit none

      integer(kind=4), parameter :: hnlen = 32             !< hostname length limit
      character(len=cwdlen) :: cwd_proc
      character(len=hnlen)  :: host_proc
      integer(kind=4)       :: pid_proc
      character(len=cwdlen), allocatable, dimension(:) :: cwd_all
      character(len=hnlen) , allocatable, dimension(:) :: host_all
      integer(kind=4)      , allocatable, dimension(:) :: pid_all
#if !defined(__INTEL_COMPILER)
      integer(kind=1)       :: getcwd, hostnm
      integer(kind=4)       :: getpid
#endif /* !__INTEL_COMPILER */
      integer :: cwd_status, host_status
      logical :: par_file_exist
      logical :: tmp_log_exist
      integer :: iproc
      logical(kind=4) :: flag

      call MPI_Init( err_mpi )
      bigbang = MPI_Wtime()

#if defined(__INTEL_COMPILER) || defined(__GFORTRAN__)
      call register_sighandler(SIGINT, abort_sigint)
#endif /* ! __INTEL_COMPILER || __GFORTRAN__ */

      call MPI_Comm_get_parent(intercomm, err_mpi)
#ifdef MPIF08
      is_spawned = (intercomm%mpi_val /= MPI_COMM_NULL%mpi_val)
#else /* !MPIF08 */
      is_spawned = (intercomm /= MPI_COMM_NULL)
      call warn("[mpisetup:init_mpi] using old, deprecated Fortran interface to MPI library")
#endif /* !MPIF08 */

      call MPI_Comm_rank(MPI_COMM_WORLD, proc, err_mpi)
      call MPI_Comm_size(MPI_COMM_WORLD, nproc, err_mpi)
      call MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, tag_ub, flag, err_mpi)

      LAST = nproc - I_ONE
      master = (proc == FIRST)
      slave  = .not. master
      have_mpi = (LAST /= FIRST)

      !Assume that any tmp_log_file existed before Piernik was started and contains invalid/outdated/... data.
      !Delete it now and keep in mind that any warn, die, printinfo or printio messages issued before this point will be lost as well.
      if (master) then
         inquire(file = tmp_log_file, exist = tmp_log_exist)
         if (tmp_log_exist) then
            open(newunit=lun, file=tmp_log_file)
            close(lun, status="delete")
         endif
         call printinfo("[mpisetup:init_mpi]: commencing...", V_DEBUG)
         if (is_spawned) call printinfo("[mpisetup:init_mpi] Piernik was called via MPI_Spawn. Additional magic will happen!", V_INFO)
      endif

      if (allocated(cwd_all) .or. allocated(host_all) .or. allocated(pid_all)) &
           call die("[mpisetup:init_mpi] cwd_all, host_all or pid_all already allocated")
      !> \deprecated BEWARE on slave it is probably enough to allocate only one element or none at all (may depend on MPI implementation)
      allocate(cwd_all(FIRST:LAST), host_all(FIRST:LAST), pid_all(FIRST:LAST))

      pid_proc    = getpid()
      host_status = hostnm(host_proc)
      cwd_status  = getcwd(cwd_proc)

      if (cwd_status /= 0) call die("[mpisetup:init_mpi] problems accessing current working directory.")

!!$      write(msg,'(3a,i8,3a)') 'mpisetup: host="',trim(host_proc),'", PID=',pid_proc,' CWD="',trim(cwd_proc),'"'
!!$      call printinfo(msg, V_DEBUG)

      call MPI_Gather(cwd_proc,  cwdlen, MPI_CHARACTER, cwd_all,  cwdlen, MPI_CHARACTER, FIRST, MPI_COMM_WORLD, err_mpi)
      call MPI_Gather(host_proc, hnlen,  MPI_CHARACTER, host_all, hnlen,  MPI_CHARACTER, FIRST, MPI_COMM_WORLD, err_mpi)
      call MPI_Gather(pid_proc,  I_ONE,  MPI_INTEGER,   pid_all,  I_ONE,  MPI_INTEGER,   FIRST, MPI_COMM_WORLD, err_mpi)

      bigbang_shift = min_bigbang - bigbang

      if (master) then
         inquire(file=par_file, exist=par_file_exist)
         if (.not. par_file_exist) call die('[mpisetup:init_mpi] Cannot find "problem.par" in the working directory',0)

         call print_char_line("-")
         call printinfo("###############     Environment     ###############", V_LOG)
         call printinfo("", V_LOG)
         call printinfo("PROCESSES:", V_LOG)
         do iproc = FIRST, LAST
            write(msg,"(a,i4,a,i6,4a)") " proc=", iproc, ", pid= ",pid_all(iproc), " @",trim(host_all(iproc)), " cwd=",trim(cwd_all(iproc))
            call printinfo(msg, V_LOG)
         enddo
         call printinfo("", V_INFO, .true.)
         write(msg, "(a,i5,a)")'   Start of the PIERNIK code on ', nproc, " processes"
         call printinfo(msg, V_INFO, .true.)
         call printinfo("", V_INFO, .true.)
         call printinfo("###############     Namelist parameters     ###############", V_LOG)
      endif

      deallocate(host_all, pid_all, cwd_all)

   end subroutine init_mpi

!> \brief Prepare to clean exit from the code

   subroutine cleanup_mpi

      use dataio_pub,      only: print_char_line, close_logs
      use MPIF,            only: MPI_COMM_WORLD, MPI_Barrier, MPI_Comm_disconnect, MPI_Finalize
      use piernik_mpi_sig, only: sig
#if defined(__INTEL_COMPILER)
      use ifport,          only: sleep
#endif /* __INTEL_COMPILER */

      implicit none

      if (master) call print_char_line("+")
      call MPI_Barrier(MPI_COMM_WORLD,err_mpi)
      if (have_mpi) call sleep(1) ! Prevent random SIGSEGVs in openmpi's MPI_Finalize
      if (is_spawned) then
         call report_to_master(sig%clean_exit)
         call MPI_Comm_disconnect(intercomm, err_mpi)
      endif
      call close_logs
      call MPI_Finalize(err_mpi)

   end subroutine cleanup_mpi

!-----------------------------------------------------------------------------
!>
!! \brief Routine used to communicate events to master Python script
!! \todo just a proof of concept
   subroutine report_to_master(ivar4, only_master)

      use constants, only: I_ONE
      use MPIF,      only: MPI_INTEGER
      use MPIFUN,    only: MPI_Send

      implicit none

      integer(kind=4),   intent(in) :: ivar4 !< integer scalar that will be send to python
      !>
      !! if more than one process calls mpisetup::report_to_master, this
      !! variable lets the slaves skip sending message
      !<
      logical, optional, intent(in) :: only_master
      integer(kind=4)               :: tag  !< master scripts accepts ANY_TAG, so it can carry meaningful value too

      if (.not.is_spawned) return

      if (present(only_master)) then
         if (only_master .and. slave) return
      endif
      tag = proc ! use proc number as tag

      call MPI_Send(ivar4, I_ONE, MPI_INTEGER, FIRST, tag, intercomm, err_mpi)

   end subroutine report_to_master

   subroutine report_string_to_master(str, only_master)

      use constants, only: I_ONE
      use MPIF,      only: MPI_INTEGER, MPI_CHARACTER
      use MPIFUN,    only: MPI_Send

      implicit none
      character(len=*),  intent(in) :: str
      logical, optional, intent(in) :: only_master

      integer(kind=4) :: tag  !< master scripts accepts ANY_TAG, so it can carry meaningful value too
      integer(kind=4) :: buf

      if (.not.is_spawned) return

      if (present(only_master)) then
         if (only_master .and. slave) return
      endif
      tag = proc ! use proc number as tag
      buf = len(str, kind=4)
      call MPI_Send(buf, I_ONE, MPI_INTEGER, FIRST, tag, intercomm, err_mpi)
      call MPI_Send(str, buf, MPI_CHARACTER, FIRST, tag, intercomm, err_mpi)

   end subroutine report_string_to_master

   integer(kind=4) function abort_sigint(signum)

      use constants, only: I_ZERO
      use MPIF,      only: MPI_Abort, MPI_COMM_WORLD

      implicit none

      integer(kind=4), intent(in) :: signum !< signal identifier

      if (master) print *, "[mpisetup:abort_sigint] CTRL-C caught, calling abort"
      ! As per MPI documentation for MPI_Abort():
      !   "This routine should not be used from within a signal handler."
      call MPI_Abort(MPI_COMM_WORLD, I_ZERO, err_mpi) ! "I too like to live dangerously." -- Austin Powers
      abort_sigint = signum

   end function abort_sigint

end module mpisetup
