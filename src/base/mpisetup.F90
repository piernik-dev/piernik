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

   use constants, only: cbuff_len, INT4, pSUM, pLAND
   use MPIF,      only: MPI_ADDRESS_KIND
#ifdef MPIF08
   use MPIF, only: MPI_Comm, MPI_Request, MPI_Op
#endif /* MPIF08 */

   implicit none

   private
   public :: cleanup_mpi, init_mpi, inflate_req, bigbang, bigbang_shift, &
        &    buffer_dim, cbuff, ibuff, lbuff, rbuff, req, req2, err_mpi, tag_ub, &
        &    master, slave, nproc, proc, FIRST, LAST, have_mpi, is_spawned, &
        &    piernik_MPI_Allreduce, piernik_MPI_Barrier, piernik_MPI_Bcast, report_to_master, &
        &    report_string_to_master, extra_barriers

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

   logical, save :: extra_barriers = .false.  !< when changed to .true. additional MPI_Barriers will be called.

#ifdef MPIF08
   type(MPI_Request), allocatable, dimension(:), target :: req        !< request array for MPI_Waitall
   type(MPI_Request), allocatable, dimension(:), target :: req2       !< second request array for MPI_Waitall
   type(MPI_Comm), protected                            :: intercomm  !< intercommunicator
   type(MPI_Op), dimension(pSUM:pLAND)                  :: mpiop      !< translation between pSUM:pLAND and MPI_SUM:MPI_LAND
#else /* !MPIF08 */
   integer(kind=4), allocatable, dimension(:),   target :: req        !< request array for MPI_Waitall
   integer(kind=4), allocatable, dimension(:),   target :: req2       !< second request array for MPI_Waitall
   integer(kind=4), protected                           :: intercomm  !< intercommunicator
   integer(kind=4), dimension(pSUM:pLAND)               :: mpiop      !< translation between pSUM:pLAND and MPI_SUM:MPI_LAND
#endif /* !MPIF08 */

   !> \warning Because we use centralized req(:) and req2(:) arrays, the routines that are using them should not call each other to avoid any interference.
   !! If you want nested non-blocking communication, only one set of MPI transactions may use these arrays.
   !< All other sets of communication should define their own req(:) array

   integer, parameter :: buffer_dim = 200                   !< size of [cilr]buff arrays used to exchange namelist parameters
   character(len=cbuff_len), dimension(buffer_dim) :: cbuff !< buffer for character parameters
   integer(kind=4),          dimension(buffer_dim) :: ibuff !< buffer for integer parameters
   real,                     dimension(buffer_dim) :: rbuff !< buffer for real parameters
   logical,                  dimension(buffer_dim) :: lbuff !< buffer for logical parameters

   interface inflate_req
      module procedure doublesize_req
      module procedure setsize_req
   end interface inflate_req

   !! \todo expand this wrapper to make it more general, unlimited polymorphism will render this obsolete
   !! Switching to pure mpi_f08 interface should allow for great simplification of these routines.
   !! Meanwhile we keep that spaghetti to not break compatibility with systems where only mpi interface is available.
   interface piernik_MPI_Bcast
      module procedure MPI_Bcast_single_logical
      module procedure MPI_Bcast_single_string
      module procedure MPI_Bcast_single_real4
      module procedure MPI_Bcast_single_real8
      module procedure MPI_Bcast_single_int4
      module procedure MPI_Bcast_single_int8
      module procedure MPI_Bcast_vec_logical
      module procedure MPI_Bcast_vec_string
      module procedure MPI_Bcast_vec_int4
      module procedure MPI_Bcast_vec_int8
      module procedure MPI_Bcast_vec_real4
      module procedure MPI_Bcast_vec_real8
      module procedure MPI_Bcast_arr2d_int4
      module procedure MPI_Bcast_arr2d_int8
      module procedure MPI_Bcast_arr2d_real4
      module procedure MPI_Bcast_arr2d_real8
      module procedure MPI_Bcast_arr3d_int4
      module procedure MPI_Bcast_arr3d_int8
      module procedure MPI_Bcast_arr3d_real4
      module procedure MPI_Bcast_arr3d_real8
   end interface piernik_MPI_Bcast

   interface piernik_MPI_Allreduce
      module procedure MPI_Allreduce_single_logical
      module procedure MPI_Allreduce_single_real4
      module procedure MPI_Allreduce_single_real8
      module procedure MPI_Allreduce_single_int4
      module procedure MPI_Allreduce_single_int8
      module procedure MPI_Allreduce_vec_real4
      module procedure MPI_Allreduce_vec_real8
      module procedure MPI_Allreduce_vec_int4
      module procedure MPI_Allreduce_vec_int8
      module procedure MPI_Allreduce_arr2d_real4
      module procedure MPI_Allreduce_arr2d_real8
      module procedure MPI_Allreduce_arr3d_real8
   end interface piernik_MPI_Allreduce

contains

!-----------------------------------------------------------------------------
!>
!! \brief Routine to start MPI
!<

   subroutine init_mpi

      use constants,     only: cwdlen, I_ONE, V_LOG, V_DEBUG, V_INFO
      use dataio_pub,    only: die, print_char_line, printinfo, msg, tmp_log_file, par_file, lun
      use MPIF,          only: MPI_COMM_WORLD, MPI_CHARACTER, MPI_INTEGER, MPI_COMM_NULL, &
           &                   MPI_SUM, MPI_MIN, MPI_MAX, MPI_LOR, MPI_LAND, MPI_TAG_UB, &
           &                   MPI_Wtime, MPI_Init, MPI_Comm_get_parent, &
           &                   MPI_Comm_rank, MPI_Comm_size
      use MPIFUN,        only: MPI_Gather, MPI_Comm_get_attr
      use signalhandler, only: SIGINT, register_sighandler
#if defined(__INTEL_COMPILER)
      use ifport,        only: getpid, getcwd, hostnm
#endif /* __INTEL_COMPILER */

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
      mpiop(:) = [ MPI_SUM, MPI_MIN, MPI_MAX, MPI_LOR, MPI_LAND ]

#if defined(__INTEL_COMPILER) || defined(__GFORTRAN__)
      call register_sighandler(SIGINT, abort_sigint)
#endif /* ! __INTEL_COMPILER || __GFORTRAN__ */

      call MPI_Comm_get_parent(intercomm, err_mpi)
#ifdef MPIF08
      is_spawned = (intercomm%mpi_val /= MPI_COMM_NULL%mpi_val)
#else /* !MPIF08 */
      is_spawned = (intercomm /= MPI_COMM_NULL)
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

!> \brief Set size of req(:) or req2(:) array for non-blocking communication on request.

   subroutine setsize_req(nreq, use_req2)

      implicit none

      integer(kind=4),   intent(in) :: nreq      !< expected maximum number of concurrent MPI requests in non-blocking parts of the code
      logical, optional, intent(in) :: use_req2  !< use req2 if .true.

      integer :: sreq
      logical :: r2

      r2 = .false.
      if (present(use_req2)) r2 = use_req2

      ! warning: spaghetti
      if (r2) then
         if (allocated(req2)) then
            sreq = size(req2)
            if (sreq < nreq) deallocate(req2)
         else
            sreq = 0
         endif

         if (sreq < nreq) allocate(req2(nreq))
      else
         if (allocated(req)) then
            sreq = size(req)
            if (sreq < nreq) deallocate(req)
         else
            sreq = 0
         endif

         if (sreq < nreq) allocate(req(nreq))
      endif

   end subroutine setsize_req

!>
!! \brief Double size of req(:) or req2(:) array for non-blocking communication on request.
!!
!! \details Perform an emergency resize by a factor of 2. Save existing values stored in req(:).
!<

   subroutine doublesize_req(use_req2)

      use dataio_pub, only: warn, msg, die
#ifdef MPIF08
      use MPIF,       only: MPI_Request
#endif /* MPIF08 */

      implicit none

      logical, optional, intent(in) :: use_req2  !< use req2 if .true.

      integer :: sreq
#ifdef MPIF08
      type(MPI_Request), allocatable, dimension(:) :: new_req    !< new request array for MPI_Waitall
#else /* !MPIF08 */
      integer(kind=4), allocatable, dimension(:)   :: new_req    !< new request array for MPI_Waitall
#endif /* !MPIF08 */
      logical :: r2

      r2 = .false.
      if (present(use_req2)) r2 = use_req2

      if (r2) then
         if (.not. allocated(req2)) call die("[mpisetup:doublesize_req] req2 not allocated")
         sreq = size(req2)
      else
         if (.not. allocated(req)) call die("[mpisetup:doublesize_req] req not allocated")
         sreq = size(req)
      endif

      if (sreq <= 0) call die("[mpisetup:doublesize_req] req or req2 is a 0-sized array")

      write(msg, '(2(a,i6))')"[mpisetup:doublesize_req] Emergency doubling size of req or req2 from ", sreq, " to ", 2 * sreq
      if (master) call warn(msg)
      allocate(new_req(2*sreq))

      if (r2) then
         new_req(1:sreq) = req2(:)
         call move_alloc(from=new_req, to=req2)
      else
         new_req(1:sreq) = req(:)
         call move_alloc(from=new_req, to=req)
      endif

   end subroutine doublesize_req

!> \brief Prepare to clean exit from the code

   subroutine cleanup_mpi

      use dataio_pub,      only: print_char_line, close_logs
      use MPIF,            only: MPI_COMM_WORLD, MPI_Barrier, MPI_Comm_disconnect, MPI_Finalize
      use piernik_mpi_sig, only: sig
#if defined(__INTEL_COMPILER)
      use ifport,          only: sleep
#endif /* __INTEL_COMPILER */

      implicit none

      if (allocated(req)) deallocate(req)
      if (allocated(req2)) deallocate(req2)

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
!! \brief Wrapper for MPI_Barrier
!<
   subroutine piernik_MPI_Barrier

      use MPIF, only: MPI_COMM_WORLD, MPI_Barrier

      implicit none

      call MPI_Barrier(MPI_COMM_WORLD, err_mpi)

   end subroutine piernik_MPI_Barrier

!-----------------------------------------------------------------------------
!>
!! \brief Wrapper for MPI_Bcast
!! Broadcast single logical value from FIRST to all
!! \todo unlimited polymorphism will obsolete me
!<
   subroutine MPI_Bcast_single_logical(lvar)

      use constants, only: I_ONE
      use MPIF,      only: MPI_LOGICAL, MPI_COMM_WORLD
      use MPIFUN,    only: MPI_Bcast

      implicit none

      logical, intent(inout) :: lvar   !< logical scalar that will be broadcasted

      call MPI_Bcast(lvar, I_ONE, MPI_LOGICAL, FIRST, MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Bcast_single_logical
!-----------------------------------------------------------------------------
!>
!! \brief Wrapper for MPI_Bcast
!! Broadcast logical vector from FIRST to all
!! \todo unlimited polymorphism will obsolete me
!<
   subroutine MPI_Bcast_vec_logical(lvar)

      use MPIF,   only: MPI_LOGICAL, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Bcast

      implicit none

      logical, dimension(:), intent(inout) :: lvar   !< logical scalar that will be broadcasted

      call MPI_Bcast(lvar, size(lvar, kind=4), MPI_LOGICAL, FIRST, MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Bcast_vec_logical
!-----------------------------------------------------------------------------
!>
!! \brief Wrapper for MPI_Bcast
!! Broadcast single string value from FIRST to all
!! \todo unlimited polymorphism will obsolete me
!<
   subroutine MPI_Bcast_single_string(cvar, clen)

      use MPIF,   only: MPI_CHARACTER, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Bcast

      implicit none

      character(len=*), intent(inout) :: cvar   !< string that will be broadcasted
      integer(kind=4),  intent(in)    :: clen   !< length of the cvar

      call MPI_Bcast(cvar, clen, MPI_CHARACTER, FIRST, MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Bcast_single_string
!-----------------------------------------------------------------------------
!>
!! \brief Wrapper for MPI_Bcast
!! Broadcast vector of strings from FIRST to all
!! \todo unlimited polymorphism will obsolete me
!<
   subroutine MPI_Bcast_vec_string(cvar, clen)

      use MPIF,   only: MPI_CHARACTER, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Bcast

      implicit none

      character(len=*), dimension(:), intent(inout) :: cvar   !< vector of strings that will be broadcasted
      integer(kind=4),                intent(in)    :: clen   !< length of the cvar

      call MPI_Bcast(cvar, clen*size(cvar, kind=4), MPI_CHARACTER, FIRST, MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Bcast_vec_string
!-----------------------------------------------------------------------------
!>
!! \brief Wrapper for MPI_Bcast
!! Broadcast single integer(kind=4) value from FIRST to all
!! \todo unlimited polymorphism will obsolete me
!<
   subroutine MPI_Bcast_single_int4(ivar4)

      use constants, only: I_ONE
      use MPIF,      only: MPI_INTEGER, MPI_COMM_WORLD
      use MPIFUN,    only: MPI_Bcast

      implicit none

      integer(kind=4), intent(inout) :: ivar4   !< integer scalar that will be broadcasted

      call MPI_Bcast(ivar4, I_ONE, MPI_INTEGER, FIRST, MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Bcast_single_int4
!-----------------------------------------------------------------------------
!>
!! \brief Wrapper for MPI_Bcast
!! Broadcast single integer(kind=8) value from FIRST to all
!! \todo unlimited polymorphism will obsolete me
!<
   subroutine MPI_Bcast_single_int8(ivar8)

      use constants, only: I_ONE
      use MPIF,      only: MPI_INTEGER8, MPI_COMM_WORLD
      use MPIFUN,    only: MPI_Bcast

      implicit none

      integer(kind=8), intent(inout) :: ivar8   !< integer scalar that will be broadcasted

      call MPI_Bcast(ivar8, I_ONE, MPI_INTEGER8, FIRST, MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Bcast_single_int8
!-----------------------------------------------------------------------------
!>
!! \brief Wrapper for MPI_Bcast
!! Broadcast single real(kind=4) value from FIRST to all
!! \todo unlimited polymorphism will obsolete me
!<
   subroutine MPI_Bcast_single_real4(rvar4)

      use constants, only: I_ONE
      use MPIF,      only: MPI_REAL, MPI_COMM_WORLD
      use MPIFUN,    only: MPI_Bcast

      implicit none

      real(kind=4), intent(inout) :: rvar4   !< integer scalar that will be broadcasted

      call MPI_Bcast(rvar4, I_ONE, MPI_REAL, FIRST, MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Bcast_single_real4
!-----------------------------------------------------------------------------
!>
!! \brief Wrapper for MPI_Bcast
!! Broadcast single real(kind=8) value from FIRST to all
!! \todo unlimited polymorphism will obsolete me
!<
   subroutine MPI_Bcast_single_real8(rvar8)

      use constants, only: I_ONE
      use MPIF,      only: MPI_DOUBLE_PRECISION, MPI_COMM_WORLD
      use MPIFUN,    only: MPI_Bcast

      implicit none

      real(kind=8), intent(inout) :: rvar8   !< real scalar that will be broadcasted

      call MPI_Bcast(rvar8, I_ONE, MPI_DOUBLE_PRECISION, FIRST, MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Bcast_single_real8
!-----------------------------------------------------------------------------
!>
!! \brief Wrapper for MPI_Bcast
!! Broadcast real(kind=4) vector from FIRST to all
!! \todo unlimited polymorphism will obsolete me
!<
   subroutine MPI_Bcast_vec_real4(rvar4)

      use MPIF,   only: MPI_REAL, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Bcast

      implicit none

      real(kind=4), dimension(:), intent(inout) :: rvar4   !< real4 vector that will be broadcasted

      call MPI_Bcast(rvar4, size(rvar4, kind=4), MPI_REAL, FIRST, MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Bcast_vec_real4
!-----------------------------------------------------------------------------
!>
!! \brief Wrapper for MPI_Bcast
!! Broadcast real(kind=8) vector from FIRST to all
!! \todo unlimited polymorphism will obsolete me
!<
   subroutine MPI_Bcast_vec_real8(rvar8)

      use MPIF,   only: MPI_DOUBLE_PRECISION, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Bcast

      implicit none

      real(kind=8), dimension(:), intent(inout) :: rvar8   !< real8 vector that will be broadcasted

      call MPI_Bcast(rvar8, size(rvar8, kind=4), MPI_DOUBLE_PRECISION, FIRST, MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Bcast_vec_real8
!-----------------------------------------------------------------------------
!>
!! \brief Wrapper for MPI_Bcast
!! Broadcast integer(kind=4) vector from FIRST to all
!! \todo unlimited polymorphism will obsolete me
!<
   subroutine MPI_Bcast_vec_int4(ivar4)

      use MPIF,   only: MPI_INTEGER, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Bcast

      implicit none

      integer(kind=4), dimension(:), intent(inout) :: ivar4   !< int4 vector that will be broadcasted

      call MPI_Bcast(ivar4, size(ivar4, kind=4), MPI_INTEGER, FIRST, MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Bcast_vec_int4
!-----------------------------------------------------------------------------
!>
!! \brief Wrapper for MPI_Bcast
!! Broadcast integer(kind=8) vector from FIRST to all
!! \todo unlimited polymorphism will obsolete me
!<
   subroutine MPI_Bcast_vec_int8(ivar8)

      use MPIF,   only: MPI_INTEGER8, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Bcast

      implicit none

      integer(kind=8), dimension(:), intent(inout) :: ivar8   !< int8 vector that will be broadcasted

      call MPI_Bcast(ivar8, size(ivar8, kind=4), MPI_INTEGER8, FIRST, MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Bcast_vec_int8
!-----------------------------------------------------------------------------
!>
!! \brief Wrapper for MPI_Bcast
!! Broadcast real(kind=4) 2D array from FIRST to all
!! \todo unlimited polymorphism will obsolete me
!<
   subroutine MPI_Bcast_arr2d_real4(rvar4)

      use MPIF,   only: MPI_REAL, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Bcast

      implicit none

      real(kind=4), dimension(:,:), intent(inout) :: rvar4   !< real4 arr2d that will be broadcasted

      call MPI_Bcast(rvar4, size(rvar4, kind=4), MPI_REAL, FIRST, MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Bcast_arr2d_real4
!-----------------------------------------------------------------------------
!>
!! \brief Wrapper for MPI_Bcast
!! Broadcast real(kind=8) 2D array from FIRST to all
!! \todo unlimited polymorphism will obsolete me
!<
   subroutine MPI_Bcast_arr2d_real8(rvar8)

      use MPIF,   only: MPI_DOUBLE_PRECISION, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Bcast

      implicit none

      real(kind=8), dimension(:,:), intent(inout) :: rvar8   !< real8 arr2d that will be broadcasted

      call MPI_Bcast(rvar8, size(rvar8, kind=4), MPI_DOUBLE_PRECISION, FIRST, MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Bcast_arr2d_real8
!-----------------------------------------------------------------------------
!>
!! \brief Wrapper for MPI_Bcast
!! Broadcast integer(kind=4) 2D array from FIRST to all
!! \todo unlimited polymorphism will obsolete me
!<
   subroutine MPI_Bcast_arr2d_int4(ivar4)

      use MPIF,   only: MPI_INTEGER, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Bcast

      implicit none

      integer(kind=4), dimension(:, :), intent(inout) :: ivar4   !< int4 arr2d that will be broadcasted

      call MPI_Bcast(ivar4, size(ivar4, kind=4), MPI_INTEGER, FIRST, MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Bcast_arr2d_int4
!-----------------------------------------------------------------------------
!>
!! \brief Wrapper for MPI_Bcast
!! Broadcast integer(kind=8) 2D array from FIRST to all
!! \todo unlimited polymorphism will obsolete me
!<
   subroutine MPI_Bcast_arr2d_int8(ivar8)

      use MPIF,   only: MPI_INTEGER8, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Bcast

      implicit none

      integer(kind=8), dimension(:,:), intent(inout) :: ivar8   !< int8 arr2d that will be broadcasted

      call MPI_Bcast(ivar8, size(ivar8, kind=4), MPI_INTEGER8, FIRST, MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Bcast_arr2d_int8
!-----------------------------------------------------------------------------
!>
!! \brief Wrapper for MPI_Bcast
!! Broadcast real(kind=4) 3D array from FIRST to all
!! \todo unlimited polymorphism will obsolete me
!<
   subroutine MPI_Bcast_arr3d_real4(rvar4)

      use MPIF,   only: MPI_REAL, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Bcast

      implicit none

      real(kind=4), dimension(:,:,:), intent(inout) :: rvar4   !< real4 arr3d that will be broadcasted

      call MPI_Bcast(rvar4, size(rvar4, kind=4), MPI_REAL, FIRST, MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Bcast_arr3d_real4
!-----------------------------------------------------------------------------
!>
!! \brief Wrapper for MPI_Bcast
!! Broadcast real(kind=8) 3D array from FIRST to all
!! \todo unlimited polymorphism will obsolete me
!<
   subroutine MPI_Bcast_arr3d_real8(rvar8)

      use MPIF,   only: MPI_DOUBLE_PRECISION, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Bcast

      implicit none

      real(kind=8), dimension(:,:,:), intent(inout) :: rvar8   !< real8 arr3d that will be broadcasted

      call MPI_Bcast(rvar8, size(rvar8, kind=4), MPI_DOUBLE_PRECISION, FIRST, MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Bcast_arr3d_real8
!-----------------------------------------------------------------------------
!>
!! \brief Wrapper for MPI_Bcast
!! Broadcast integer(kind=4) 3D array from FIRST to all
!! \todo unlimited polymorphism will obsolete me
!<
   subroutine MPI_Bcast_arr3d_int4(ivar4)

      use MPIF,   only: MPI_INTEGER, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Bcast

      implicit none

      integer(kind=4), dimension(:,:, :), intent(inout) :: ivar4   !< int4 arr3d that will be broadcasted

      call MPI_Bcast(ivar4, size(ivar4, kind=4), MPI_INTEGER, FIRST, MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Bcast_arr3d_int4
!-----------------------------------------------------------------------------
!>
!! \brief Wrapper for MPI_Bcast
!! Broadcast integer(kind=8) 3D array from FIRST to all
!! \todo unlimited polymorphism will obsolete me
!<
   subroutine MPI_Bcast_arr3d_int8(ivar8)

      use MPIF,   only: MPI_INTEGER8, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Bcast

      implicit none

      integer(kind=8), dimension(:,:,:), intent(inout) :: ivar8   !< int8 arr3d that will be broadcasted

      call MPI_Bcast(ivar8, size(ivar8, kind=4), MPI_INTEGER8, FIRST, MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Bcast_arr3d_int8
!-----------------------------------------------------------------------------
!>
!! \brief Wrapper for MPI_Allreduce with MPI_IN_PLACE
!! \todo unlimited polymorphism will obsolete me
!<
   subroutine MPI_Allreduce_single_logical(lvar, reduction)

      use constants, only: I_ONE
      use MPIF,      only: MPI_LOGICAL, MPI_IN_PLACE, MPI_COMM_WORLD
      use MPIFUN,    only: MPI_Allreduce

      implicit none

      logical,         intent(inout) :: lvar      !< logical that will be reduced
      integer(kind=4), intent(in)    :: reduction !< integer to mark a reduction type

      call MPI_Allreduce(MPI_IN_PLACE, lvar, I_ONE, MPI_LOGICAL, mpiop(reduction), MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Allreduce_single_logical
!-----------------------------------------------------------------------------
!>
!! \brief Wrapper for MPI_Allreduce with MPI_IN_PLACE
!! \todo unlimited polymorphism will obsolete me
!<
   subroutine MPI_Allreduce_single_int4(ivar4, reduction)

      use constants, only: I_ONE
      use MPIF,      only: MPI_INTEGER, MPI_IN_PLACE, MPI_COMM_WORLD
      use MPIFUN,    only: MPI_Allreduce

      implicit none

      integer(kind=4), intent(inout) :: ivar4     !< int4 that will be reduced
      integer(kind=4), intent(in)    :: reduction !< integer to mark a reduction type

      call MPI_Allreduce(MPI_IN_PLACE, ivar4, I_ONE, MPI_INTEGER, mpiop(reduction), MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Allreduce_single_int4
!-----------------------------------------------------------------------------
!>
!! \brief Wrapper for MPI_Allreduce with MPI_IN_PLACE
!! \todo unlimited polymorphism will obsolete me
!<
   subroutine MPI_Allreduce_single_int8(ivar8, reduction)

      use constants, only: I_ONE
      use MPIF,      only: MPI_INTEGER8, MPI_IN_PLACE, MPI_COMM_WORLD
      use MPIFUN,    only: MPI_Allreduce

      implicit none

      integer(kind=8), intent(inout) :: ivar8     !< int8 that will be reduced
      integer(kind=4), intent(in)    :: reduction !< integer to mark a reduction type

      call MPI_Allreduce(MPI_IN_PLACE, ivar8, I_ONE, MPI_INTEGER8, mpiop(reduction), MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Allreduce_single_int8
!-----------------------------------------------------------------------------
!>
!! \brief Wrapper for MPI_Allreduce with MPI_IN_PLACE
!! \todo unlimited polymorphism will obsolete me
!<
   subroutine MPI_Allreduce_vec_int4(ivar4, reduction)

      use MPIF,   only: MPI_INTEGER, MPI_IN_PLACE, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Allreduce

      implicit none

      integer(kind=4), dimension(:), intent(inout) :: ivar4     !< int4 that will be reduced
      integer(kind=4),               intent(in)    :: reduction !< integer to mark a reduction type

      call MPI_Allreduce(MPI_IN_PLACE, ivar4, size(ivar4, kind=4), MPI_INTEGER, mpiop(reduction), MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Allreduce_vec_int4
!-----------------------------------------------------------------------------
!>
!! \brief Wrapper for MPI_Allreduce with MPI_IN_PLACE
!! \todo unlimited polymorphism will obsolete me
!<
   subroutine MPI_Allreduce_vec_int8(ivar8, reduction)

      use MPIF,   only: MPI_INTEGER8, MPI_IN_PLACE, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Allreduce

      implicit none

      integer(kind=8), dimension(:), intent(inout) :: ivar8     !< int8 that will be reduced
      integer(kind=4),               intent(in)    :: reduction !< integer to mark a reduction type

      call MPI_Allreduce(MPI_IN_PLACE, ivar8, size(ivar8, kind=4), MPI_INTEGER8, mpiop(reduction), MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Allreduce_vec_int8
!-----------------------------------------------------------------------------
!>
!! \brief Wrapper for MPI_Allreduce with MPI_IN_PLACE
!! \todo unlimited polymorphism will obsolete me
!<
   subroutine MPI_Allreduce_single_real4(rvar4, reduction)

      use constants, only: I_ONE
      use MPIF,      only: MPI_REAL, MPI_IN_PLACE, MPI_COMM_WORLD
      use MPIFUN,    only: MPI_Allreduce

      implicit none

      real(kind=4),    intent(inout) :: rvar4     !< real4 that will be reduced
      integer(kind=4), intent(in)    :: reduction !< integer to mark a reduction type

      call MPI_Allreduce(MPI_IN_PLACE, rvar4, I_ONE, MPI_REAL, mpiop(reduction), MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Allreduce_single_real4
!-----------------------------------------------------------------------------
!>
!! \brief Wrapper for MPI_Allreduce with MPI_IN_PLACE
!! \todo unlimited polymorphism will obsolete me
!<
   subroutine MPI_Allreduce_single_real8(rvar8, reduction)

      use constants, only: I_ONE
      use MPIF,      only: MPI_DOUBLE_PRECISION, MPI_IN_PLACE, MPI_COMM_WORLD
      use MPIFUN,    only: MPI_Allreduce

      implicit none

      real(kind=8),    intent(inout) :: rvar8     !< real8 that will be reduced
      integer(kind=4), intent(in)    :: reduction !< integer to mark a reduction type

      call MPI_Allreduce(MPI_IN_PLACE, rvar8, I_ONE, MPI_DOUBLE_PRECISION, mpiop(reduction), MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Allreduce_single_real8
!-----------------------------------------------------------------------------
!>
!! \brief Wrapper for MPI_Allreduce with MPI_IN_PLACE
!! \todo unlimited polymorphism will obsolete me
!<
   subroutine MPI_Allreduce_vec_real4(rvar4, reduction)

      use MPIF,   only: MPI_REAL, MPI_IN_PLACE, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Allreduce

      implicit none

      real(kind=4), dimension(:), intent(inout) :: rvar4     !< real4 that will be reduced
      integer(kind=4),            intent(in)    :: reduction !< integer to mark a reduction type

      call MPI_Allreduce(MPI_IN_PLACE, rvar4, size(rvar4, kind=4), MPI_REAL, mpiop(reduction), MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Allreduce_vec_real4
!-----------------------------------------------------------------------------
!>
!! \brief Wrapper for MPI_Allreduce with MPI_IN_PLACE
!! \todo unlimited polymorphism will obsolete me
!<
   subroutine MPI_Allreduce_vec_real8(rvar8, reduction)

      use MPIF,   only: MPI_DOUBLE_PRECISION, MPI_IN_PLACE, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Allreduce

      implicit none

      real(kind=8), dimension(:), intent(inout) :: rvar8     !< real8 that will be reduced
      integer(kind=4),            intent(in)    :: reduction !< integer to mark a reduction type

      call MPI_Allreduce(MPI_IN_PLACE, rvar8, size(rvar8, kind=4), MPI_DOUBLE_PRECISION, mpiop(reduction), MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Allreduce_vec_real8
!-----------------------------------------------------------------------------
!>
!! \brief Wrapper for MPI_Allreduce with MPI_IN_PLACE
!! \todo unlimited polymorphism will obsolete me
!<
   subroutine MPI_Allreduce_arr3d_real8(rvar8, reduction)

      use MPIF,   only: MPI_DOUBLE_PRECISION, MPI_IN_PLACE, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Allreduce

      implicit none

      real(kind=8), dimension(:,:,:), intent(inout) :: rvar8     !< real8 that will be reduced
      integer(kind=4),                intent(in)    :: reduction !< integer to mark a reduction type

      call MPI_Allreduce(MPI_IN_PLACE, rvar8, size(rvar8, kind=4), MPI_DOUBLE_PRECISION, mpiop(reduction), MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Allreduce_arr3d_real8
!-----------------------------------------------------------------------------
!>
!! \brief Wrapper for MPI_Allreduce with MPI_IN_PLACE
!! \todo unlimited polymorphism will obsolete me
!<
   subroutine MPI_Allreduce_arr2d_real8(rvar8, reduction)

      use MPIF,   only: MPI_DOUBLE_PRECISION, MPI_IN_PLACE, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Allreduce

      implicit none

      real(kind=8), dimension(:,:), intent(inout) :: rvar8     !< real8 that will be reduced
      integer(kind=4),              intent(in)    :: reduction !< integer to mark a reduction type

      call MPI_Allreduce(MPI_IN_PLACE, rvar8, size(rvar8, kind=4), MPI_DOUBLE_PRECISION, mpiop(reduction), MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Allreduce_arr2d_real8
!-----------------------------------------------------------------------------
!>
!! \brief Wrapper for MPI_Allreduce with MPI_IN_PLACE
!! \todo unlimited polymorphism will obsolete me
!<
   subroutine MPI_Allreduce_arr2d_real4(rvar4, reduction)

      use MPIF,   only: MPI_REAL, MPI_IN_PLACE, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Allreduce

      implicit none

      real(kind=4), dimension(:,:), intent(inout) :: rvar4     !< real4 that will be reduced
      integer(kind=4),              intent(in)    :: reduction !< integer to mark a reduction type

      call MPI_Allreduce(MPI_IN_PLACE, rvar4, size(rvar4, kind=4), MPI_REAL, mpiop(reduction), MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Allreduce_arr2d_real4
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
