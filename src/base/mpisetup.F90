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
#include "macros.h"
!>
!! \brief (KK)
!!
!! In this module following namelists of parameters are specified:
!! \copydetails mpisetup::init_mpi
!<
module mpisetup

   use constants, only: cbuff_len, INT4

   implicit none

   private
   public :: cleanup_mpi, init_mpi, inflate_req, &
        &    buffer_dim, cbuff, ibuff, lbuff, rbuff, req, status, ierr, procmask, &
        &    master, slave, nproc, proc, FIRST, LAST, comm, have_mpi

   integer(kind=4), protected :: nproc, proc, LAST          !< number of processes, rank of my process, rank of last process
   integer(kind=4), protected :: comm, ierr                 !< global communicator, error status
   integer(kind=INT4), parameter :: FIRST = 0               !< the rank of the master process

   logical, protected :: master, slave      !< shortcuts for testing proc == FIRST
   logical, protected :: have_mpi           !< .true. when run on more than one processor

   integer, allocatable, dimension(:)   :: req    !< request array for MPI_Waitall
   integer, allocatable, dimension(:,:) :: status !< status array for MPI_Waitall
   !> \warning Because we use one centralized req(:) and status(:,:) arrays, the routines that are using them should not call each other to avoid any interference.
   !! If you want nested non-blocking communication, only one set of MPI transactions may use these arrays.
   !< All other sets of communication should define their own req(:) and status(:,:) arrays

   integer, dimension(:), allocatable :: procmask !< (FIRST:LAST)-sized auxiliary array for searching overlaps, neighbours etc. BEWARE: antiparallel

   integer, parameter :: buffer_dim = 200                   !< size of [cilr]buff arrays used to exchange namelist parameters
   character(len=cbuff_len), dimension(buffer_dim) :: cbuff !< buffer for character parameters
   integer,                  dimension(buffer_dim) :: ibuff !< buffer for integer parameters
   real,                     dimension(buffer_dim) :: rbuff !< buffer for real parameters
   logical,                  dimension(buffer_dim) :: lbuff !< buffer for logical parameters

contains

!-----------------------------------------------------------------------------
!>
!! \brief Routine to start MPI
!<

   subroutine init_mpi

      use constants,  only: cwdlen, I_ONE
      use mpi,        only: MPI_COMM_WORLD, MPI_CHARACTER, MPI_INTEGER
      use dataio_pub, only: die, printinfo, msg, cwd, ansi_white, ansi_black, tmp_log_file
      use dataio_pub, only: par_file, ierrh, namelist_errh, compare_namelist, cmdl_nml, lun  ! QA_WARN required for diff_nml

      implicit none

      integer, parameter :: hnlen = 32             !< hostname length limit
      character(len=cwdlen) :: cwd_proc
      character(len=hnlen)  :: host_proc
      integer(kind=4)       :: pid_proc
      character(len=cwdlen), allocatable, dimension(:) :: cwd_all
      character(len=hnlen) , allocatable, dimension(:) :: host_all
      integer(kind=4)      , allocatable, dimension(:) :: pid_all
      integer(kind=1)       :: getcwd, hostnm
      integer(kind=4)       :: getpid
      integer :: cwd_status, host_status
      logical :: par_file_exist
      logical :: tmp_log_exist
      integer :: iproc

      call MPI_Init( ierr )
      comm = MPI_COMM_WORLD

      call MPI_Comm_rank(comm, proc, ierr)
      call MPI_Comm_size(comm, nproc, ierr)

      LAST = nproc-I_ONE
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
#ifdef VERBOSE
         call printinfo("[mpisetup:init_mpi]: commencing...")
#endif /* VERBOSE */
      endif

      if (allocated(cwd_all) .or. allocated(host_all) .or. allocated(pid_all)) call die("[mpisetup:init_mpi] cwd_all, host_all or pid_all already allocated")
      !> \deprecated BEWARE on slave it is probably enough to allocate only one element or none at all (may depend on MPI implementation)
      allocate(cwd_all(FIRST:LAST), host_all(FIRST:LAST), pid_all(FIRST:LAST))

      pid_proc    = getpid()
      host_status = hostnm(host_proc)
      cwd_status  = getcwd(cwd_proc)

      if (cwd_status /= 0) call die("[mpisetup:init_mpi] problems accessing current working directory.")
#ifdef DEBUG
      write(msg,'(3a,i6,3a)') 'mpisetup: host="',trim(host_proc),'", PID=',pid_proc,' CWD="',trim(cwd_proc),'"'
      call printinfo(msg)
#endif /* DEBUG */

      call MPI_Gather(cwd_proc,  cwdlen, MPI_CHARACTER, cwd_all,  cwdlen, MPI_CHARACTER, FIRST, comm, ierr)
      call MPI_Gather(host_proc, hnlen,  MPI_CHARACTER, host_all, hnlen,  MPI_CHARACTER, FIRST, comm, ierr)
      call MPI_Gather(pid_proc,  I_ONE, MPI_INTEGER,   pid_all,  I_ONE, MPI_INTEGER,   FIRST, comm, ierr)

      if (master) then
         inquire(file=par_file, exist=par_file_exist)
         if (.not. par_file_exist) call die('[mpisetup:init_mpi] Cannot find "problem.par" in the working directory',0)

         call printinfo("------------------------------------------------------------------------------------------------------", .false.)
         call printinfo("###############     Environment     ###############", .false.)
         call printinfo("", .false.)
         call printinfo("PROCESSES:", .false.)
         do iproc = FIRST, LAST
            write(msg,"(a,i4,a,i6,4a)") " proc=", iproc, ", pid= ",pid_all(iproc), " @",trim(host_all(iproc)), " cwd=",trim(cwd_all(iproc))
            call printinfo(msg, .false.)
         enddo
         call printinfo("", .true.)
         write(msg,"(5a,i5)") 'Start of the',ansi_white,' PIERNIK ',ansi_black,'code. No. of procs = ', nproc
         call printinfo(msg, .true.)
         call printinfo("", .true.)
         call printinfo("###############     Namelist parameters     ###############", .false.)
      endif

      deallocate(host_all, pid_all, cwd_all)

      if (allocated(procmask)) call die("[mpisetup:init_mpi] procmask already allocated")
      allocate(procmask(FIRST:LAST))

   end subroutine init_mpi

!-----------------------------------------------------------------------------
!>
!! Increase size of req(:) and status(:,:) arrays for non-blocking communication on request.
!<
   subroutine inflate_req(nreq)

      use dataio_pub, only: warn, msg
      use mpi,        only: MPI_STATUS_SIZE

      implicit none

      integer, intent(in) :: nreq

      integer :: sreq

      if (allocated(req)) then
         sreq = size(req)
         if (sreq < nreq) then
            write(msg, '(2(a,i6))')"[mpisetup:inflate_req] reallocating req and status from ",sreq," to ",nreq
            if (master) call warn(msg)
            deallocate(req)
            if (allocated(status)) deallocate(status)
         endif
      else
         sreq = 0
      endif

      if (sreq < nreq) then
         allocate(req(nreq))
         allocate(status(MPI_STATUS_SIZE, nreq))
      endif

   end subroutine inflate_req

!-----------------------------------------------------------------------------
!>
!! \brief Prepare to clean exit from the code
!<
   subroutine cleanup_mpi

      use dataio_pub, only: printinfo

      implicit none

      if (allocated(procmask)) deallocate(procmask)
      if (allocated(req)) deallocate(req)
      if (allocated(status)) deallocate(status)

      if (master) call printinfo("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++", .false.)
      call MPI_Barrier(comm,ierr)
      if (have_mpi) call sleep(1) ! Prevent random SIGSEGVs in openmpi's MPI_Finalize
      call MPI_Finalize(ierr)

   end subroutine cleanup_mpi

end module mpisetup
