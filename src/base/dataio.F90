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
!=====================================================================
!!
!!  dataio module responsible for data output
!!
!=====================================================================
!

!>
!! \brief [DW] Module containing all main routines  responsible for data output
!!
!!
!! In this module following namelists of parameters are specified:
!!
!! @b END_CONTROL
!!
!! \f[ \begin{tabular}{ | p{3cm} | p{3cm} | p{4cm} | p{8cm} | } \hline &&&\\
!! {\bf parameter} & {\bf default value} & {\bf possible values} & {\bf description} \\ &&&\\ \hline \hline &&&\\
!! nend & 1   & integer & step number to end simulation \\ &&&\\ \hline &&&\\
!! tend &-1.0 & real    & simulation time to end        \\ &&&\\ \hline
!! \end{tabular} \f]
!!
!! @b RESTART_CONTROL
!!
!! \f[ \begin{tabular}{ | p{3cm} | p{3cm} | p{4cm} | p{8cm} | } \hline &&&\\
!! {\bf parameter} & {\bf default value} & {\bf possible values} & {\bf description} \\ &&&\\ \hline \hline &&&\\
!! restart  & 'last' & 'last' or another string of characters & 'last': automatic choise of the last restart file regardless of "nrestart" value; if smth else is set: "nrestart" value is fixing \\ &&&\\ \hline &&&\\
!! new\_id  & ''     & string of characters & three character string to change run\_id when restarting simulation (e.g. to avoid overwriting of the output from the previous (pre-restart) simulation; if new\_id = '' then run\_id is still used) \\ &&&\\ \hline &&&\\
!! nrestart & 3      & integer & number of restart file to be read while restart is not set to '' \\ &&&\\ \hline &&&\\
!! resdel   & 0      & integer & number of recent restart dumps which should be saved; each n-resdel-1 restart file is supposed to be deleted while writing n restart file \\ &&&\\ \hline
!! \end{tabular} \f]
!!
!! @b OUTPUT_CONTROL
!!
!! \f[ \begin{tabular}{ | p{3cm} | p{3cm} | p{4cm} | p{8cm} | } \hline &&&\\
!! {\bf parameter} & {\bf default value} & {\bf possible values} & {\bf description} \\ &&&\\ \hline \hline &&&\\
!! dt\_hdf & 0.0 & real & time between successive hdf dumps          \\ &&&\\ \hline &&&\\
!! dt\_res & 0.0 & real & time between successive restart file dumps \\ &&&\\ \hline &&&\\
!! dt\_tsl & 0.0 & real & time between successive timeslice dumps    \\ &&&\\ \hline &&&\\
!! dt\_log & 0.0 & real & time between successive log dumps          \\ &&&\\ \hline &&&\\
!! dt\_plt & 0.0 & real & time between successive domain slices files dumps \\ &&&\\ \hline &&&\\
!! ix      &     & integer & index in x-direction of slice to dump in plt files \\ &&&\\ \hline &&&\\
!! iy      &     & integer & index in y-direction of slice to dump in plt files \\ &&&\\ \hline &&&\\
!! iz      &     & integer & index in z-direction of slice to dump in plt files \\ &&&\\ \hline &&&\\
!! domain  & 'phys\_domain' & 'phys\_domain' or 'full\_domain' & string to choose if boundaries have to be dumped in hdf files \\ &&&\\ \hline &&&\\
!! vars    & ''  & 'dens', 'velx', 'vely', 'velz', 'ener' and some more & array of 4-character strings standing for variables to dump in hdf files \\ &&&\\ \hline &&&\\
!! mag\_center & 'no' & 'yes'/'no' & choise to dump magnetic fields values from cell centers or not (if not then values from cell borders) \\ &&&\\ \hline
!! \end{tabular} \f] \f[
!! \begin{tabular}{ | p{3cm} | p{3cm} | p{4cm} | p{8cm} | } &&&\\
!! min\_disk\_space\_MB  & 100 & integer & minimum disk space in MB to continue simulation {\bf (currently not used)} \\ &&&\\ \hline &&&\\
!! sleep\_minutes        & 0   & integer & minutes of sleeping time before continue simulation \\ &&&\\ \hline &&&\\
!! sleep\_seconds        & 0   & integer & seconds of sleeping time before continue simulation \\ &&&\\ \hline &&&\\
!! user\_message\_file   & trim(cwd)//'/msg' & string similar to default value & path to possible user message file containing dt\_xxx changes or orders to dump/stop/end simulation \\ &&&\\ \hline &&&\\
!! system\_message\_file & '/tmp/piernik\_msg' & string of characters similar to default value  & path to possible system (UPS) messege file contaning orders to dump/stop/end simulation \\ &&&\\ \hline
!! \end{tabular} \f]
!!
!! \todo check and if necessary bring back usefullness of min_disk_space_MB parameter
!<
module dataio
!>
!! \brief (KK)
!!
!! \todo check the usefullness of wait logical variable
!<
! Written by G. Kowal
! Modified for this code and extended by M.Hanasz
   use types
   use mpisetup
   use initproblem
#ifdef SN_SRC
   use snsources
#endif /* SN_SRC */
#ifdef HDF5
   use dataio_hdf5, only : write_hdf5, write_restart_hdf5, read_restart_hdf5, &
       init_hdf5, cleanup_hdf5, write_plot
#endif /* HDF5 */
   implicit none

   type(hdf) :: chdf

   integer               :: nend                   !< number of the step to end simulation
   integer               :: nstep                  !< current number of simulation timestep
   integer               :: istep                  !< current number of substep (related to integration order)
   integer               :: nstep_start            !< number of start timestep
   real                  :: tend                   !< simulation time to end

   integer, parameter    :: nvarsmx = 16           !< maximum number of variables to dump in hdf files
   character(len=3)      :: new_id                 !< three character string to replace run_id when restarting simulation (if new_id = '' then run_id is still used)
   character(len=16)     :: restart                !< choise of restart file: if restart = 'last': automatic choise of the last restart file regardless of "nrestart" value; if smth else is set: "nrestart" value is fixing
   character(len=16)     :: domain                 !< string to choose if boundaries have to be dumped in hdf files
   character(len=16)     :: mag_center             !< choise to dump magnetic fields values from cell centers or not (if not then values from cell borders)
   integer               :: nrestart               !< number of restart file to be read while restart is not set to ''
   integer               :: resdel                 !< number of recent restart dumps which should be saved; each n-resdel-1 restart file is supposed to be deleted while writing n restart file
   real                  :: dt_hdf                 !< time between successive hdf dumps
   real                  :: dt_res                 !< time between successive restart file dumps
   real                  :: dt_tsl                 !< time between successive timeslice dumps
   real                  :: dt_log                 !< time between successive log dumps
   real                  :: dt_plt                 !< time between successive domain slices files dumps
   integer               :: min_disk_space_MB      !< minimum disk space in MB to continue simulation (currently not used)
   integer               :: sleep_minutes          !< minutes of sleeping time before continue simulation
   integer               :: sleep_seconds          !< seconds of sleeping time before continue simulation
   character(len=160)    :: user_message_file      !< path to possible user message file containing dt_xxx changes or orders to dump/stop/end simulation
   character(len=160)    :: system_message_file    !< path to possible system (UPS) messege file contaning orders to dump/stop/end simulation
   integer               :: ix                     !< index in x-direction of slice to dump in plt files
   integer               :: iy                     !< index in y-direction of slice to dump in plt files
   integer               :: iz                     !< index in z-direction of slice to dump in plt files
   integer               :: iv                     !< work index to count successive variables to dump in hdf files
   character(len=4), dimension(nvarsmx) :: vars    !< array of 4-character strings standing for variables to dump in hdf files

   integer               :: tsl_lun = 2            !< luncher for timeslice file
   integer               :: log_lun = 3            !< luncher for log file
   integer               :: nhdf                   !< current number of hdf file
   integer               :: nres                   !< current number of restart file
   integer               :: ntsl                   !< current number of timeslice file
   integer               :: nlog                   !< current number of log file
   integer               :: step_hdf               !< number of simulation timestep corresponding to values dumped in hdf file
   integer               :: step_res               !< number of simulation timestep corresponding to values dumped in restart file
   integer               :: nhdf_start             !< number of hdf file for the first hdf dump in simulation run
   integer               :: nres_start             !< number of restart file for the first restart dump in simulation run
   real                  :: t_start                !< time in simulation of start simulation run
   real                  :: last_hdf_time          !< time in simulation of the last resent hdf file dump
   character(len=128)    :: log_file               !< path to the current log file
   character(len=2)      :: pc1                    !< string of two characters to mark current block in x-direction in hdf4 files
   character(len=2)      :: pc2                    !< string of two characters to mark current block in y-direction in hdf4 files
   character(len=2)      :: pc3                    !< string of two characters to mark current block in z-direction in hdf4 files
   logical               :: tsl_firstcall          !< logical value to start a new timeslice file

   logical               :: wait                   !< logical value to have a break in simulation (currently not used)

   integer               :: nchar                  !< number of characters in a user/system message
   character             :: msg*16                 !< string of characters - content of a user/system message
   real                  :: msg_param              !< parameter changed by a user/system message

   character             :: hostfull*(80)
   character             :: host*8
   character             :: fhost*10
   character             :: fpid*10
   integer               :: pid
   integer               :: uid
   integer               :: ihost
   integer               :: scstatus
   real                  :: vx_max                 !< maximum value of velocity x-component
   real                  :: vy_max                 !< maximum value of velocity y-component
   real                  :: vz_max                 !< maximum value of velocity z-component
   real                  :: va2max
   real                  :: va_max                 !< maximum value of Alfven velocity
   real                  :: cs2max
   real                  :: cs_max                 !< maximum value of current sound speed distribution
   real                  :: dens_min               !< minimum value of current density distribution
   real                  :: dens_max               !< maximum value of current density distribution
   real                  :: pres_min               !< minimum value of current pressure distribution
   real                  :: pres_max               !< maximum value of current pressure distribution

#ifdef MAGNETIC
   real                  :: divb_max               !< maximum value of current divergence of magnetic induction distribution
   real                  :: b_min                  !< minimum value of current magnetic induction distribution
   real                  :: b_max                  !< maximum value of current magnetic induction distribution
#endif /* MAGNETIC */
#ifndef ISO
   real                  :: temp_min               !< minimum value of current temperature distribution
   real                  :: temp_max               !< maximum value of current temperature distribution
#endif /* ISO */
#ifdef COSM_RAYS
   real                  :: encr_min               !< minimum value of current cosmic ray energy density distribution
   real                  :: encr_max               !< maximum value of current cosmic ray energy density distribution
#endif /* COSM_RAYS */
#ifdef HDF5
    character(len=128) :: filename                 !< string of characters indicating currently used file
#endif /* HDF5 */

   namelist /END_CONTROL/ nend, tend
   namelist /RESTART_CONTROL/ restart, new_id, nrestart, resdel
   namelist /OUTPUT_CONTROL/ dt_hdf, dt_res, dt_tsl, dt_log, dt_plt, ix, iy, iz, &
                             domain, vars, mag_center, &
                             min_disk_space_MB, sleep_minutes, sleep_seconds, &
                             user_message_file, system_message_file

   contains

     subroutine set_container(chdf)
       use types
       implicit none
       type(hdf), intent(out) :: chdf

       chdf%nstep = nstep
       chdf%nhdf  = nhdf
       chdf%ntsl  = ntsl
       chdf%nres  = nres
       chdf%nlog  = nlog
       chdf%step_hdf = step_hdf
       chdf%log_lun = log_lun
       chdf%last_hdf_time = last_hdf_time
       chdf%log_file = log_file
       chdf%nrestart = nrestart
       chdf%domain  = domain

     end subroutine set_container

     subroutine get_container(chdf)
       use types
       implicit none
       type(hdf), intent(in) :: chdf

       nstep = chdf%nstep
       nhdf = chdf%nhdf
       ntsl = chdf%ntsl
       nres = chdf%nres
       nlog = chdf%nlog
       step_hdf =  chdf%step_hdf
       log_lun = chdf%log_lun
       last_hdf_time = chdf%last_hdf_time
       log_file = chdf%log_file
       nrestart = chdf%nrestart
       domain   = chdf%domain

     end subroutine get_container


!---------------------------------------------------------------------
!
! inititalizes dataio parameters
!
!---------------------------------------------------------------------
!
   subroutine init_dataio
      use errh,            only : namelist_errh
      use initproblem,     only : problem_name,run_id
      use version,         only : nenv,env
      use fluidboundaries, only : all_fluid_boundaries
      use fluidindex
#ifdef MAGNETIC
      use magboundaries,   only : all_mag_boundaries
#endif /* MAGNETIC */
      implicit none
      integer              :: ierrh
      integer(kind=1)      :: getpid
      integer(kind=1)      :: hostnm
      integer(kind=1)      :: system
      integer              :: system_status, i
      character(LEN=160)   :: system_command
      character(LEN=100)   :: par_file, tmp_log_file

      restart = 'last'   ! 'last': automatic choise of the last restart file
                         ! regardless of "nrestart" value;
                         ! if smth else is set: "nrestart" value is fixing
      new_id  = ''
      nrestart=  3
      resdel  = 0

      dt_hdf = 0.0
      dt_res = 0.0
      dt_tsl = 0.0
      dt_log = 0.0
      dt_plt = 0.0
      domain = 'phys_domain'
      vars(:)   = ''
      mag_center= 'no'
      min_disk_space_MB = 100
      sleep_minutes   = 0
      sleep_seconds   = 0
      user_message_file   = trim(cwd)//'/msg'
      system_message_file = '/tmp/piernik_msg'

      wait  = .false.
      tsl_firstcall = .true.

      nhdf  = 0
      ntsl  = 0
      nres  = 0
      nlog  = 0

      step_hdf  = -1
      step_res  = -1

      pc1 = '00'
      pc2 = '00'
      pc3 = '00'

      nend = 1
      tend = -1.0

      if(psize(1) .gt. 1) pc1 = '0x'
      if(psize(2) .gt. 1) pc2 = '0x'
      if(psize(3) .gt. 1) pc3 = '0x'

      pid = getpid()

      scstatus = hostnm(hostfull)
      ihost = index(hostfull,'.')
      if (ihost .eq. 0) ihost = index(hostfull,' ')
      host = hostfull(1:ihost-1)

      if(proc .eq. 0) then
         par_file = trim(cwd)//'/problem.par'
         tmp_log_file = trim(cwd)//'/tmp.log'
         open(1,file=par_file)
         read(unit=1,nml=OUTPUT_CONTROL,iostat=ierrh)
         call namelist_errh(ierrh,'OUTPUT_CONTROL')
         close(1)
         open(1,file=par_file)
         read(unit=1,nml=RESTART_CONTROL,iostat=ierrh)
         call namelist_errh(ierrh,'RESTART_CONTROL')
         close(1)
         open(1,file=par_file)
         read(unit=1,nml=END_CONTROL,iostat=ierrh)
         call namelist_errh(ierrh,'END_CONTROL')
         close(1)
         open(3, file=tmp_log_file, position='append')
         write(unit=3,nml=OUTPUT_CONTROL)
         write(unit=3,nml=RESTART_CONTROL)
         close(3)

!  namelist /END_CONTROL/ nend, tend
         ibuff(1)  = nend

         rbuff(1)  = tend

!  namelist /RESTART_CONTROL/ restart, new_id, nrestart, resdel

         cbuff(20) = restart
         cbuff(21) = new_id

         ibuff(20) = nrestart
         ibuff(21) = resdel

!  namelist /OUTPUT_CONTROL/ dt_hdf, dt_res, dt_tsl, domain, vars, mag_center, &
!                            min_disk_space_MB, sleep_minutes, ix, iy, iz&
!                            user_message_file, system_message_file

         ibuff(40) = min_disk_space_MB
         ibuff(41) = sleep_minutes
         ibuff(42) = sleep_seconds
         ibuff(43) = ix
         ibuff(44) = iy
         ibuff(45) = iz

         rbuff(40) = dt_hdf
         rbuff(41) = dt_res
         rbuff(42) = dt_tsl
         rbuff(43) = dt_log
         rbuff(44) = dt_plt

         cbuff(40) = domain

         do iv = 1, nvarsmx
            cbuff(40+iv) = vars(iv)
         enddo

         cbuff(90) = mag_center
         cbuff(91) = user_message_file(1:32)
         cbuff(92) = system_message_file(1:32)

         call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
         call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
         call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      else

         call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
         call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
         call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

!  namelist /END_CONTROL/ nend, tend
         nend                = ibuff(1)

         tend                = rbuff(1)
!  namelist /RESTART_CONTROL/ restart, new_id, nrestart, resdel

         restart             = trim(cbuff(20))
         new_id              = trim(cbuff(21))

         nrestart            = ibuff(20)
         resdel              = ibuff(21)

!  namelist /OUTPUT_CONTROL/ dt_hdf, dt_res, dt_tsl, domain, vars, mag_center, &
!                            min_disk_space_MB, sleep_minutes, ix, iy, iz &
!                            user_message_file, system_message_file

         min_disk_space_MB   = ibuff(40)
         sleep_minutes       = ibuff(41)
         sleep_seconds       = ibuff(42)
         ix                  = ibuff(43)
         iy                  = ibuff(44)
         iz                  = ibuff(45)

         dt_hdf              = rbuff(40)
         dt_res              = rbuff(41)
         dt_tsl              = rbuff(42)
         dt_log              = rbuff(43)
         dt_plt              = rbuff(44)

         domain              = trim(cbuff(40))
         do iv=1, nvarsmx
            vars(iv)         = trim(cbuff(40+iv))
         enddo

         mag_center          = trim(cbuff(90))

         user_message_file   = trim(cbuff(91))
         system_message_file = trim(cbuff(92))

      endif

      last_hdf_time = -dt_hdf

#ifdef HDF5
      call init_hdf5(vars,ix,iy,iz,dt_plt)
#endif /* HDF5 */

      if(proc == 0 .and. restart .eq. 'last') call find_last_restart(nrestart)
      call MPI_BARRIER(comm,ierr)
      call MPI_BCAST(nrestart, 1, MPI_INTEGER, 0, comm, ierr)
      if (nrestart == 0) then
         if (proc == 0) then
            write (log_file,'(a,a1,a3,a1,i3.3,a4)') &
                 trim(problem_name),'_', run_id,'_',nrestart,'.log'
            tmp_log_file = trim(cwd)//'/tmp.log'
            log_file = trim(cwd)//'/'//trim(log_file)
            system_command = 'mv '//trim(tmp_log_file)//' '//trim(log_file)
            system_status = SYSTEM(system_command)
            open(3, file=log_file, position='append')
            do i=1,nenv
               write(3,*) trim(env(i))
            enddo
            close(3)
         endif
         call set_container(chdf); chdf%nres = nrestart
         call write_data(output='all')

      else
         if (proc == 0) then
            write (log_file,'(a,a1,a3,a1,i3.3,a4)') &
               trim(problem_name),'_', run_id,'_',nrestart,'.log'
            tmp_log_file = trim(cwd)//'/tmp.log'
            log_file = trim(cwd)//'/'//log_file
            system_command = 'mv '//trim(tmp_log_file)//' '//trim(log_file)
            system_status = SYSTEM(system_command)
         endif
#ifdef HDF5
         call set_container(chdf); chdf%nres = nrestart
         call read_restart_hdf5(chdf)
         call get_container(chdf)
#else /* HDF5 */
         nres = nrestart
         call read_restart
#endif
         nstep_start = nstep
         t_start     = t
         nres_start  = nrestart
         nhdf_start  = nhdf-1
         if(proc==0) then
            open(3, file=log_file, position='append')
            do i=1,nenv
               write(3,*) trim(env(i))
            enddo
            close(3)
         endif
         if(new_id .ne. '') run_id=new_id
      endif
#ifdef HDF5
      call MPI_BCAST(log_file, 32, MPI_CHARACTER, 0, comm, ierr)
      call set_container(chdf)
#endif /* HDF5 */
      call all_fluid_boundaries
#ifdef MAGNETIC
      call all_mag_boundaries
#endif /* MAGNETIC */

   end subroutine init_dataio

   subroutine cleanup_dataio

      implicit none
#ifdef HDF5
      call cleanup_hdf5
#endif
   end subroutine cleanup_dataio

   subroutine user_msg_handler(end_sim)
      implicit none
      logical, intent(inout) :: end_sim
      integer :: tsleep
!      do while(disk_full)

!--- process 0 checks for messages

      if(proc == 0) then
         call read_file_msg
      endif
      call MPI_BCAST(msg,       16, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(msg_param, 16, MPI_DOUBLE_PRECISION, 0, comm, ierr)

!---  if a user message is received then:
      if (len_trim(msg) /= 0) then
         if(trim(msg) == 'res' .or. trim(msg) == 'dump' ) then
#ifndef HDF5
            nres = max(nres,1)
            call write_restart
            nres = nres + 1
            step_res = nstep
#else
            if(proc==0) then
              write (filename,'(a,a1,a3,a1,i4.4,a4)') &
                trim(problem_name),'_', run_id,'_',nres,'.res'
            endif
            call MPI_BCAST(filename, 128, MPI_CHARACTER, 0, comm, ierr)
            call set_container(chdf)
            call write_restart_hdf5(filename,chdf)
#endif
         endif
         if(trim(msg) .eq. 'hdf') then
#ifndef HDF5
            call write_hdf
#else /* HDF5 */
            call set_container(chdf)
            call write_hdf5(chdf)
#endif /* HDF5 */
            nhdf = nhdf + 1
            step_hdf = nstep
         endif
         if(trim(msg) .eq. 'log')    call write_log
         if(trim(msg) .eq. 'tsl')    call write_timeslice
         if(trim(msg) .eq. 'tend')   tend   = msg_param
         if(trim(msg) .eq. 'nend')   nend   = msg_param
         if(trim(msg) .eq. 'dtres')  dt_res = msg_param
         if(trim(msg) .eq. 'dthdf')  dt_hdf = msg_param
         if(trim(msg) .eq. 'dtlog')  dt_log = msg_param
         if(trim(msg) .eq. 'dttsl')  dt_tsl = msg_param

         if(trim(msg) .eq. 'sleep') then
            tsleep = 60*msg_param
            call sleep(tsleep)
         endif

         if(trim(msg) .eq. 'stop') end_sim = .true.
      endif
!  enddo ! while disk is full
   end subroutine user_msg_handler

!---------------------------------------------------------------------
!
! controls data dumping
!
!---------------------------------------------------------------------
!
   subroutine write_data(output)
#ifdef USER_IO
      use initproblem, only : user_io_routine
#endif /* USER_IO */
      implicit none
      character  :: output*3

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      if (output .eq. 'log' .or. output .eq. 'end') then
         call write_log
      endif

      if (output .eq. 'log' .or. output .eq. 'end') then
         call write_timeslice
      endif

!    CALL checkdf

#ifdef USER_IO
      call user_io_routine
#endif /* USER_IO */

#ifdef HDFSWEEP
      if (output .ne. 'gpt') then
         call write_log
         call write_timeslice
#else /* HDFSWEEP */
      if (dt_hdf .gt. 0.0 .and. nstep .gt. step_hdf .and. output .ne. 'gpt') then
#endif /* HDFSWEEP */

         if ((t-last_hdf_time) .ge. dt_hdf &
                .or. output .eq. 'hdf' .or. output .eq. 'end') then
#ifndef HDF5
            call write_hdf
#else /* HDF5 */
            call set_container(chdf)
            call write_hdf5(chdf)
#endif /* HDF5 */

            if((t-last_hdf_time) .ge. dt_hdf) last_hdf_time = last_hdf_time + dt_hdf
            if((t-last_hdf_time) .ge. dt_hdf) last_hdf_time = t ! additional control
                          ! in the case of changing dt_hdf into smaller value via msg
            nhdf = nhdf + 1
            step_hdf = nstep
         endif
      endif

      if (dt_res .gt. 0.0 .and. nstep .gt. step_res) then
         if ((nres-nres_start) .lt. (int((t-t_start) / dt_res) + 1) &
                .or. output .eq. 'res' .or. output .eq. 'end') then
            if (nres > 0) then
#ifdef HDF5
               if(proc==0) then
                  write (filename,'(a,a1,a3,a1,i4.4,a4)') &
                    trim(problem_name),'_', run_id,'_',nres,'.res'
               endif
               call MPI_BCAST(filename, 128, MPI_CHARACTER, 0, comm, ierr)
               call set_container(chdf)
               call write_restart_hdf5(filename,chdf)
#else /* HDF5 */
               call write_restart
#endif
            endif
            nres = nres + 1
            step_res = nstep
         endif
      endif
#ifdef HDF5
      call set_container(chdf)
      call write_plot(chdf)
#endif /* HDF5 */

   end subroutine write_data

!---------------------------------------------------------------------
!
! dumps data to hdf file
!
!---------------------------------------------------------------------
!
   subroutine write_hdf
      use mpisetup, only: cwd
      use arrays, only : wa,b,u
      use grid, only : dx,dy,dz,xmin,xmax,ymin,ymax,zmin,zmax,nxd,nyd,nzd,nb
      use grid, only : nx,ny,nz,nxb,nyb,nzb,x,y,z
      use initproblem, only : problem_name, run_id

      use fluidindex,  only : nfluid
      use fluidindex,  only : ibx,iby,ibz
      use fluidindex,  only : nvar, iarr_all_dn,iarr_all_mx,iarr_all_my,iarr_all_mz

#ifndef ISO
      use fluidindex, only : iarr_all_en
#endif /* ISO */

#ifdef COSM_RAYS
      use initcosmicrays, only : iecr
#endif /* COSM_RAYS */

#ifdef GRAV
      use arrays, only : gp
#endif /* GRAV */

      implicit none

      character(len=128) :: file_name_hdf,file_name_disp
      character(LEN=4)   :: varname

      integer :: sd_id, sds_id, dim_id, iostatus
      integer :: rank, comp_type
      integer, dimension(3) :: dims, istart, stride
      integer, dimension(1) :: comp_prm
      integer :: sfstart, sfend, sfsnatt, sfcreate, sfwdata, sfscompress, sfendacc &
               , sfdimid, sfsdmname, sfsdscale, sfsdmstr

      integer :: iv, ifl, iw
      integer :: nxo = 1, nyo = 1, nzo = 1, &
                 iso = 1, jso = 1, kso = 1, &
                 ieo = 1, jeo = 1, keo = 1

      real(kind=4), dimension(:,:,:), allocatable :: tmparr
      real(kind=4), dimension(:), allocatable     :: temp_scl
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      if(domain .eq. 'full_domain') then
         nxo = nx
         nyo = ny
         nzo = nz
         if(nxd /= 1) then
            iso = 1
            ieo = nx
         else
            iso = 1
            ieo = 1
         endif
         if(nyd /= 1) then
            jso = 1
            jeo = ny
         else
            jso = 1
            jeo = 1
         endif
         if(nzd /= 1) then
            kso = 1
            keo = nz
         else
            kso = 1
            keo = 1
         endif

      else if(domain .eq. 'phys_domain') then
         nxo = nxb
         nyo = nyb
         nzo = nzb
         if(nxd /= 1) then
            iso = nb+1
            ieo = nb+nxb
         else
            iso = 1
            ieo = 1
         endif
         if(nyd /= 1) then
            jso = nb+1
            jeo = nb+nyb
         else
            jso = 1
            jeo = 1
         endif
         if(nzd /= 1) then
            kso = nb+1
            keo = nb+nzb
         else
            kso = 1
            keo = 1
         endif
      endif

      allocate(tmparr(iso:ieo,jso:jeo,kso:keo))

      rank = 3
      comp_type = 4
      comp_prm(1) = 6
      dims(1) = nxo
      dims(2) = nyo
      dims(3) = nzo
      istart(:) = 0
      stride(:) = 1

!  generate filename
!

      write (file_name_hdf,'(a,a1,a,a1,a3,a1,i2.2,a1,i2.2,a1,i2.2,a1,i4.4,a4)') &
               trim(cwd),'/',trim(problem_name),'_', run_id, &
               '_',pcoords(1),'_',pcoords(2),'_',pcoords(3),'_',nhdf, '.hdf'

      write (file_name_disp,'(a,a1,a3,a1,a2,a1,a2,a1,a2,a1,i4.4,a4)') &
               trim(problem_name),'_', run_id,'_',pc1,'_',pc2,'_',pc3,'_',nhdf, '.hdf'

      sd_id = sfstart(trim(file_name_hdf), 4)

! write attributes
!
      iostatus = sfsnatt( sd_id, 'problem' , 4, 32, problem_name)
      iostatus = sfsnatt( sd_id, 'run_id'  , 4, 32, run_id      )
      iostatus = sfsnatt( sd_id, 'domain'  , 4, 32, domain      )

      iostatus = sfsnatt( sd_id, 'psize1'  , 23, 1, psize(1) )
      iostatus = sfsnatt( sd_id, 'psize2'  , 23, 1, psize(2) )
      iostatus = sfsnatt( sd_id, 'psize3'  , 23, 1, psize(3) )

      iostatus = sfsnatt( sd_id, 'pcoords1', 23, 1, pcoords(1) )
      iostatus = sfsnatt( sd_id, 'pcoords2', 23, 1, pcoords(2) )
      iostatus = sfsnatt( sd_id, 'pcoords3', 23, 1, pcoords(3) )

      iostatus = sfsnatt( sd_id, 'dims1'   , 23, 1, dims(1)    )
      iostatus = sfsnatt( sd_id, 'dims2'   , 23, 1, dims(2)    )
      iostatus = sfsnatt( sd_id, 'dims3'   , 23, 1, dims(3)    )

      iostatus = sfsnatt( sd_id, 'nxd'     , 23, 1, nxd     )
      iostatus = sfsnatt( sd_id, 'nyd'     , 23, 1, nyd     )
      iostatus = sfsnatt( sd_id, 'nzd'     , 23, 1, nzd     )

      iostatus = sfsnatt( sd_id, 'nxb'     , 23, 1, nxb     )
      iostatus = sfsnatt( sd_id, 'nyb'     , 23, 1, nyb     )
      iostatus = sfsnatt( sd_id, 'nzb'     , 23, 1, nzb     )
      iostatus = sfsnatt( sd_id, 'nb'      , 23, 1, nb      )

      iostatus = sfsnatt( sd_id, 'xmin'    , 6,  1, xmin    )
      iostatus = sfsnatt( sd_id, 'xmax'    , 6,  1, xmax    )
      iostatus = sfsnatt( sd_id, 'ymin'    , 6,  1, ymin    )
      iostatus = sfsnatt( sd_id, 'ymax'    , 6,  1, ymax    )
      iostatus = sfsnatt( sd_id, 'zmin'    , 6,  1, zmin    )
      iostatus = sfsnatt( sd_id, 'zmax'    , 6,  1, zmax    )

      iostatus = sfsnatt( sd_id, 'nstep'   , 23, 1, nstep          )
      iostatus = sfsnatt( sd_id, 'time'    ,  6, 1, t              )
      iostatus = sfsnatt( sd_id, 'timestep',  6, 1, dt             )

!    do ifl=1,nfluid
!    write(gammaifl,'(a5,i1)') 'gamma',ifl
!    iostatus = sfsnatt( sd_id, gammaifl   ,  6, 1, gamma(ifl)     )
!    enddo


! write selected problem dependent parameters
#ifdef SN_SRC
      iostatus = sfsnatt( sd_id, 'nsn'      , 24, 1, nsn     )
#endif /* SN_SRC */

      iv = 1
      iw = 0
      ifl = 1

      do while (len_trim(vars(iv)) .ne. 0)

         select case(vars(iv))
         case ('dens')
            write(varname,'(a3,i1)') 'den',ifl
            wa(iso:ieo,jso:jeo,kso:keo) = u(iarr_all_dn(ifl),iso:ieo,jso:jeo,kso:keo)
            call next_fluid_or_var(ifl,iw,nfluid)

         case ('velx')
            write(varname,'(a3,i1)') 'vlx',ifl
            wa(iso:ieo,jso:jeo,kso:keo) = u(iarr_all_mx(ifl),iso:ieo,jso:jeo,kso:keo) &
                                        / u(iarr_all_dn(ifl),iso:ieo,jso:jeo,kso:keo)
            call next_fluid_or_var(ifl,iw,nfluid)

         case ('vely')
            write(varname,'(a3,i1)') 'vly',ifl
            wa(iso:ieo,jso:jeo,kso:keo) = u(iarr_all_my(ifl),iso:ieo,jso:jeo,kso:keo) &
                                        / u(iarr_all_dn(ifl),iso:ieo,jso:jeo,kso:keo)
            call next_fluid_or_var(ifl,iw,nfluid)

         case ('velz')
            write(varname,'(a3,i1)') 'vlz',ifl
            wa(iso:ieo,jso:jeo,kso:keo) = u(iarr_all_mz(ifl),iso:ieo,jso:jeo,kso:keo) &
                                        / u(iarr_all_dn(ifl),iso:ieo,jso:jeo,kso:keo)
            call next_fluid_or_var(ifl,iw,nfluid)

#ifdef ISO
         case ('ener')
            write(varname,'(a3,i1)') 'ene',ifl
            wa(iso:ieo,jso:jeo,kso:keo) = 0.5*(u(iarr_all_mx(ifl),iso:ieo,jso:jeo,kso:keo)**2 &
                                              +u(iarr_all_my(ifl),iso:ieo,jso:jeo,kso:keo)**2 &
                                              +u(iarr_all_mz(ifl),iso:ieo,jso:jeo,kso:keo)**2)&
                                              /u(iarr_all_dn(ifl),iso:ieo,jso:jeo,kso:keo)
            call next_fluid_or_var(ifl,iw,nfluid)

#else /* ISO */
         case ('ener')
            write(varname,'(a3,i1)') 'ene',ifl
            wa(iso:ieo,jso:jeo,kso:keo) = u(iarr_all_en(ifl),iso:ieo,jso:jeo,kso:keo)
            call next_fluid_or_var(ifl,iw,nfluid)
#endif /* ISO */

#ifdef COSM_RAYS
         case ('encr')
            varname = 'encr'
            wa(iso:ieo,jso:jeo,kso:keo) = u(iecr, iso:ieo, jso:jeo, kso:keo)
#endif /* COSM_RAYS */

         case ('divb')
            wa(iso:ieo-1,jso:jeo-1,kso:keo-1) = &
             (b(ibx,iso+1:ieo,jso:jeo-1,kso:keo-1) - b(ibx,iso:ieo-1,jso:jeo-1,kso:keo-1))*dy*dz &
            +(b(iby,iso:ieo-1,jso+1:jeo,kso:keo-1) - b(iby,iso:ieo-1,jso:jeo-1,kso:keo-1))*dx*dz &
            +(b(ibz,iso:ieo-1,jso:jeo-1,kso+1:keo) - b(ibz,iso:ieo-1,jso:jeo-1,kso:keo-1))*dx*dy
            wa = abs(wa)
            wa(ieo,:,:) = wa(ieo-1,:,:)
            wa(:,jeo,:) = wa(:,jeo-1,:)
            wa(:,:,keo) = wa(:,:,keo-1)
#ifdef GRAV
         case ('gpot')
            varname = 'gpot'
            wa(iso:ieo,jso:jeo,kso:keo) = gp(iso:ieo,jso:jeo,kso:keo)
#endif /* GRAV */

         case ('magx')
            varname = 'magx'
            if(domain .eq. 'full_domain') then
               if(mag_center .eq. 'yes') then
                  wa(:,:,:) = 0.5*b(ibx,:,:,:)
                  wa(:,:,:) = wa(:,:,:)  + cshift(wa(:,:,:),shift=1,dim=1)
               else
                  wa(:,:,:) = b(ibx,:,:,:)
               endif
            else if(domain .eq. 'phys_domain') then
               if(mag_center .eq. 'yes') then
                  wa(iso:ieo,jso:jeo,kso:keo) = 0.5*(b(ibx,iso:ieo,jso:jeo,kso:keo))
                  wa(iso:ieo,jso:jeo,kso:keo) = wa(iso:ieo,jso:jeo,kso:keo) &
                                              + 0.5*(b(ibx,iso+1:ieo+1,jso:jeo,kso:keo))
               else
                  wa(iso:ieo,jso:jeo,kso:keo) = b(ibx,iso:ieo,jso:jeo,kso:keo)
               endif
            endif

         case ('magy')
            varname = 'magy'
            if(domain .eq. 'full_domain') then
               if(mag_center .eq. 'yes') then
                  wa(:,:,:) = 0.5*b(iby,:,:,:)
                  wa(:,:,:) = wa(:,:,:)  + cshift(wa(:,:,:),shift=1,dim=2)
               else
                  wa(:,:,:) = b(iby,:,:,:)
               endif
            else if(domain .eq. 'phys_domain') then
               if(mag_center .eq. 'yes') then
                  wa(iso:ieo,jso:jeo,kso:keo) = 0.5*(b(iby,iso:ieo,jso:jeo,kso:keo))
                  wa(iso:ieo,jso:jeo,kso:keo) = wa(iso:ieo,jso:jeo,kso:keo) &
                                              + 0.5*(b(iby,iso:ieo,jso+1:jeo+1,kso:keo))
               else
                  wa(iso:ieo,jso:jeo,kso:keo) = b(iby,iso:ieo,jso:jeo,kso:keo)
               endif
            endif

         case ('magz')
            varname = 'magz'
            if(domain .eq. 'full_domain') then
               if(mag_center .eq. 'yes') then
                  wa(:,:,:) = 0.5*b(ibz,:,:,:)
                  wa(:,:,:) = wa(:,:,:)  + cshift(wa(:,:,:),shift=1,dim=3)
               else
                  wa(:,:,:) = b(ibz,:,:,:)
               endif
            else if(domain .eq. 'phys_domain') then
               if(mag_center .eq. 'yes') then
                  wa(iso:ieo,jso:jeo,kso:keo) = 0.5*(b(ibz,iso:ieo,jso:jeo,kso:keo))
                  wa(iso:ieo,jso:jeo,kso:keo) = wa(iso:ieo,jso:jeo,kso:keo) &
                                              + 0.5*(b(ibz,iso:ieo,jso:jeo,kso+1:keo+1))
               else
                  wa(iso:ieo,jso:jeo,kso:keo) = b(ibz,iso:ieo,jso:jeo,kso:keo)
               endif
            endif

         case default
            print *, 'Variable ', vars(iv), ' is not defined! Skipping.'
         end select

! write data
!
         where( abs(wa) < 1.e-12 )  wa = sign(1e-12,wa)

         sds_id = sfcreate(sd_id, varname, 5, rank, dims)
         iostatus = sfscompress(sds_id, comp_type, comp_prm)
         tmparr = real(wa(iso:ieo,jso:jeo,kso:keo),4)
         iostatus = sfwdata(sds_id, istart, stride, dims, tmparr)

! write coords
!
         if(.not.allocated(temp_scl)) allocate(temp_scl(iso:ieo))
         temp_scl(iso:ieo) = real(x(iso:ieo),4)
         dim_id = sfdimid( sds_id, 0 )
         iostatus = sfsdmname( dim_id, 'xc' )
         iostatus = sfsdscale( dim_id, dims(1), 5, temp_scl(:))
         iostatus = sfsdmstr ( dim_id, 'X', 'pc', '' )
         if(allocated(temp_scl)) deallocate(temp_scl)

         if(.not.allocated(temp_scl)) allocate(temp_scl(jso:jeo))
         temp_scl(jso:jeo) = real(y(jso:jeo),4)
         dim_id = sfdimid( sds_id, 1 )
         iostatus = sfsdmname( dim_id, 'yc' )
         iostatus = sfsdscale( dim_id, dims(2), 5, temp_scl(:))
         iostatus = sfsdmstr ( dim_id, 'Y', 'pc', '' )
         if(allocated(temp_scl)) deallocate(temp_scl)

         if(.not.allocated(temp_scl)) allocate(temp_scl(kso:keo))
         temp_scl(kso:keo) = real(z(kso:keo),4)
         dim_id = sfdimid( sds_id, 2 )
         iostatus = sfsdmname( dim_id, 'zc' )
         iostatus = sfsdscale( dim_id, dims(3), 5, temp_scl(:))
         iostatus = sfsdmstr ( dim_id, 'Z', 'pc', '' )
         if(allocated(temp_scl)) deallocate(temp_scl)

         iostatus = sfendacc(sds_id)

         if(iw .eq. 0) iv = iv + 1
      end do

      iostatus = sfend(sd_id)

      if(proc.eq. 0) then
         open(log_lun, file=log_file, position='append')
         write(log_lun,*) 'Writing output   file: ', trim(file_name_disp)
         write(*,*)       'Writing output   file: ', trim(file_name_disp)
         close(log_lun)
      endif

      if(allocated(tmparr)) deallocate(tmparr)

   end subroutine write_hdf

   subroutine next_fluid_or_var(ifluid,ivar,nfluids)
      implicit none
      integer ifluid,ivar,nfluids
      if(ifluid .lt. nfluids) then
         ifluid=ifluid+1
         ivar=1
      else
         ivar=0
         ifluid=1
      endif
      return
   end subroutine next_fluid_or_var

!---------------------------------------------------------------------
!
! writes restart file
!
!---------------------------------------------------------------------
!
   subroutine write_restart

      use fluidindex,   only : nvar,nmag
      use arrays,       only : u,b
      use grid,         only : nxb,nyb,nzb,x,y,z,nx,ny,nz
      use grid,         only : xmin,xmax,ymin,ymax,zmin,zmax,nxd,nyd,nzd,nb
      use initproblem,  only : problem_name, run_id
#ifdef GRAV
      use arrays,       only : gp
#endif /* GRAV */

      implicit none

      character(len=128) :: file_name_res,file_name_disp,file_name_last
      character*160 syscom
      logical lastres_exist

      integer :: sd_id, sds_id, dim_id, scstatus
      integer :: iostatus, ranku, rankb, rank3d, comp_type
      integer, dimension(4) :: dimsu, dimsb, istart, stride
      integer, dimension(3) :: dims, dims3d
      integer, dimension(1) :: comp_prm
      integer :: sfstart, sfend, sfsnatt, sfcreate, sfwdata, sfscompress, sfendacc &
               , sfdimid, sfsdmname, sfsdscale

      integer(kind=1) :: system

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!  prepare data dimensions
!
      comp_type   = 4
      comp_prm(1) = 6
      istart(:)   = 0
      stride(:)   = 1

      dims(1)     = nx
      dims(2)     = ny
      dims(3)     = nz

      ranku       = 4
      dimsu(1)    = nvar
      dimsu(2)    = dims(1)
      dimsu(3)    = dims(2)
      dimsu(4)    = dims(3)

      rankb       = 4
      dimsb(1)    = nmag
      dimsb(2)    = dims(1)
      dimsb(3)    = dims(2)
      dimsb(4)    = dims(3)

      rank3d      = 3
      dims3d(1)    = dims(1)
      dims3d(2)    = dims(2)
      dims3d(3)    = dims(3)

!  generate filename
!
      write (file_name_res,'(a,a1,a,a1,a3,a1,3(i2.2,a1),i3.3,a4)') &
                    trim(cwd),'/',trim(problem_name),'_', run_id,  &
               '_',pcoords(1),'_',pcoords(2),'_',pcoords(3),'_',nres, '.res'

      write (file_name_disp,'(a,a1,a3,a1,3(a2,a1),i3.3,a4)') &
                    trim(problem_name),'_', run_id, &
               '_',pc1,'_',pc2,'_',pc3,'_',nres, '.res'

      if((resdel .ne. 0) .and. (nres .gt. resdel)) then
         write (file_name_last,'(a,a1,a,a1,a3,a1,3(i2.2,a1),i3.3,a4)') &
                        trim(cwd),'/',trim(problem_name),'_', run_id,  &
              '_',pcoords(1),'_',pcoords(2),'_',pcoords(3),'_',nres-resdel, '.res'
      endif

      sd_id = sfstart(file_name_res, 4)

!!  ATTRIBUTES
!!
! write config attributes
!
      iostatus = sfsnatt( sd_id, 'problem' , 4, 32, problem_name)
      iostatus = sfsnatt( sd_id, 'run_id'  , 4, 32, run_id      )
      iostatus = sfsnatt( sd_id, 'domain'  , 4, 32, domain      )

      iostatus = sfsnatt( sd_id, 'psize1'  , 23, 1, psize(1) )
      iostatus = sfsnatt( sd_id, 'psize2'  , 23, 1, psize(2) )
      iostatus = sfsnatt( sd_id, 'psize3'  , 23, 1, psize(3) )

      iostatus = sfsnatt( sd_id, 'pcoords1', 23, 1, pcoords(1) )
      iostatus = sfsnatt( sd_id, 'pcoords2', 23, 1, pcoords(2) )
      iostatus = sfsnatt( sd_id, 'pcoords3', 23, 1, pcoords(3) )

      iostatus = sfsnatt( sd_id, 'dimsu'   , 23, 1, nvar )
      iostatus = sfsnatt( sd_id, 'dims1'   , 23, 1, dims(1) )
      iostatus = sfsnatt( sd_id, 'dims2'   , 23, 1, dims(2) )
      iostatus = sfsnatt( sd_id, 'dims3'   , 23, 1, dims(3) )

      iostatus = sfsnatt( sd_id, 'nxd'     , 23,  1, nxd     )
      iostatus = sfsnatt( sd_id, 'nyd'     , 23,  1, nyd     )
      iostatus = sfsnatt( sd_id, 'nzd'     , 23,  1, nzd     )

      iostatus = sfsnatt( sd_id, 'nxb'     , 23,  1, nxb     )
      iostatus = sfsnatt( sd_id, 'nyb'     , 23,  1, nyb     )
      iostatus = sfsnatt( sd_id, 'nzb'     , 23,  1, nzb     )
      iostatus = sfsnatt( sd_id, 'nb'      , 23,  1, nb      )

      iostatus = sfsnatt( sd_id, 'xmin'    ,  6,  1, xmin   )
      iostatus = sfsnatt( sd_id, 'xmax'    ,  6,  1, xmax   )
      iostatus = sfsnatt( sd_id, 'ymin'    ,  6,  1, ymin   )
      iostatus = sfsnatt( sd_id, 'ymax'    ,  6,  1, ymax   )
      iostatus = sfsnatt( sd_id, 'zmin'    ,  6,  1, zmin   )
      iostatus = sfsnatt( sd_id, 'zmax'    ,  6,  1, zmax   )

! write evolution attributes
!
      iostatus = sfsnatt( sd_id, 'nstep'   , 24,  1, nstep   )
      iostatus = sfsnatt( sd_id, 'time'    ,  6,  1, t       )
      iostatus = sfsnatt( sd_id, 'timestep',  6,  1, dt      )

! write dataio attributes
!
      iostatus = sfsnatt( sd_id, 'nres'     , 24, 1, nres + 1 )
      iostatus = sfsnatt( sd_id, 'nhdf'     , 24, 1, nhdf     )
      iostatus = sfsnatt( sd_id, 'ntsl'     , 24, 1, ntsl     )
      iostatus = sfsnatt( sd_id, 'nlog'     , 24, 1, nlog     )
      iostatus = sfsnatt( sd_id, 'step_res' , 24, 1, nstep     )
      iostatus = sfsnatt( sd_id, 'step_hdf' , 24, 1, step_hdf )
      iostatus = sfsnatt( sd_id, 'last_hdf_time', 6, 1, last_hdf_time)

! write selected problem dependent parameters
#ifdef SN_SRC
      iostatus = sfsnatt( sd_id, 'nsn'      , 24, 1, nsn     )
#endif /* SN_SRC */

! write array of integer scalars
!
!    sds_id   = sfcreate(sd_id, 'intscal', 23, 1, nintscal)
!    iostatus = sfscompress(sds_id, comp_type, comp_prm)
!    iostatus = sfwdata(sds_id, istart, stride, nintscal, intscal)

! write array of real scalars
!
!    sds_id   = sfcreate(sd_id, 'rlscal', 6, 1, nrlscal)
!    iostatus = sfscompress(sds_id, comp_type, comp_prm)
!    iostatus = sfwdata(sds_id, istart, stride, nrlscal, rlscal)

! write initial vertical density profile
!
!    sds_id   = sfcreate(sd_id, 'dprof', 6, 1, nz)
!    iostatus = sfscompress(sds_id, comp_type, comp_prm)
!    iostatus = sfwdata(sds_id, istart, stride, nz, dprof)

! write fluid variables array
!
      sds_id = sfcreate(sd_id, 'fluid_vars', 6, ranku, dimsu)
      iostatus = sfscompress(sds_id, comp_type, comp_prm)
      iostatus = sfwdata(sds_id, istart, stride, dimsu, u)

! write magnetic field array
!
      sds_id = sfcreate(sd_id, 'mag_field', 6, rankb, dimsb)
      iostatus = sfscompress(sds_id, comp_type, comp_prm)
      iostatus = sfwdata(sds_id, istart, stride, dimsb, b)

#ifdef GRAV
! write gravitational potential
!
      sds_id = sfcreate(sd_id, 'grav_pot', 6, rank3d, dims3d)
      iostatus = sfscompress(sds_id, comp_type, comp_prm)
      iostatus = sfwdata(sds_id, istart, stride, dims3d, gp)
#endif /* GRAV */

! write coords
!
      dim_id = sfdimid( sds_id, 1 )
      iostatus = sfsdmname( dim_id, 'x' )
      iostatus = sfsdscale( dim_id, dimsu(2), 6, x)

      dim_id = sfdimid( sds_id, 2 )
      iostatus = sfsdmname( dim_id, 'y' )
      iostatus = sfsdscale( dim_id, dimsu(3), 6, y)

      dim_id = sfdimid( sds_id, 3 )
      iostatus = sfsdmname( dim_id, 'z' )
      iostatus = sfsdscale( dim_id, dimsu(4), 6, z)

      iostatus = sfendacc(sds_id)

      iostatus = sfend(sd_id)

      if(proc.eq. 0) then
         open(log_lun, file=log_file, position='append')
         write(log_lun,*) 'Writing restart  file: ', trim(file_name_disp)
         write(*,*)       'Writing restart  file: ', trim(file_name_disp)
         close(log_lun)
      endif

      if((resdel .ne. 0) .and. (nres .gt. resdel)) then
         inquire(file=file_name_last, exist=lastres_exist)
         if(lastres_exist) then
            syscom='rm -f'//file_name_last
            scstatus = system(syscom)
            write(*,*) trim(file_name_last),' removed'
         else
            write(*,*) trim(file_name_last),' does not exist'
         endif
      endif

   end subroutine write_restart

!---------------------------------------------------------------------
!
! read restart file
!
!---------------------------------------------------------------------
!
   subroutine read_restart !(all)
      use fluidindex,   only : nvar, nmag
      use arrays,       only : u,b
      use grid,         only : nx,ny,nz
      use initproblem,  only : problem_name, run_id
#ifdef GRAV
      use arrays,       only : gp
#endif /* GRAV */

      implicit none

      character(len=128) :: file_name_res,file_name_disp
      logical file_exist, log_exist

      integer :: sd_id, sds_id, attr_index, sds_index
      integer :: iostatus, ranku, rankb, rank3d
      integer, dimension(4) :: dimsu, dimsb, istart, stride
      integer, dimension(3) :: dims3d
      integer :: sfstart, sfend, sffattr, sfrnatt, sfn2index, sfselect, sfrdata, sfendacc

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!  prepare data dimensions
!
      istart(:) = 0
      stride(:) = 1

      ranku    = 4
      dimsu(1) = nvar
      dimsu(2) = nx
      dimsu(3) = ny
      dimsu(4) = nz

      rankb    = 4
      dimsb(1) = nmag
      dimsb(2) = nx
      dimsb(3) = ny
      dimsb(4) = nz

      rank3d   = 3
      dims3d(1) = nx
      dims3d(2) = ny
      dims3d(3) = nz

!  generate filename
!
      write (file_name_res,'(a,a1,a,a1,a3,a1,3(i2.2,a1),i3.3,a4)') &
              trim(cwd),'/',trim(problem_name),'_', run_id,  &
              '_',pcoords(1),'_',pcoords(2),'_',pcoords(3),'_',nres, '.res'

      write (file_name_disp,'(a,a1,a3,a1,3(a2,a1),i3.3,a4)') &
              trim(problem_name),'_', run_id, &
              '_',pc1,'_',pc2,'_',pc3,'_',nres, '.res'

      if(proc.eq. 0) then
         write(*,*)       'Reading restart  file: ', trim(file_name_disp)
      endif
      if(proc==0) then
         inquire(file =log_file , exist = log_exist)
         if(log_exist .eqv. .true.) then
            open(log_lun, file=log_file, position='append')
            write(log_lun,*) 'Reading restart  file: ', trim(file_name_disp)
            close(log_lun)
         endif
      endif

      inquire(file = file_name_res, exist = file_exist)
      if(file_exist .eqv. .false.) then
         if(log_exist) then
            open(log_lun, file=log_file, position='append')
            write(log_lun,*) 'Restart  file: ', trim(file_name_res), &
                         '               does not exist.  ABORTING !!! '
            close(log_lun)
         endif

         write(*,*)       'Restart  file: ', trim(file_name_res), &
                         '               does not exist.  ABORTING !!! '
         call MPI_BARRIER(comm,ierr)
         call mpistop
         stop
      endif


      sd_id = sfstart(file_name_res, 1)

! read evolution attributes
!
      attr_index = sffattr( sd_id, 'nstep'     )
      iostatus = sfrnatt( sd_id, attr_index, nstep   )
      attr_index = sffattr( sd_id, 'time'     )
      iostatus = sfrnatt( sd_id, attr_index, t   )
      attr_index = sffattr( sd_id, 'timestep' )
      iostatus = sfrnatt( sd_id, attr_index, dt     )

! read dataio attributes
!
      attr_index = sffattr( sd_id, 'nres'      )
      iostatus = sfrnatt( sd_id, attr_index, nres     )
      attr_index = sffattr( sd_id, 'nhdf'      )
      iostatus = sfrnatt( sd_id, attr_index, nhdf     )
      attr_index = sffattr( sd_id, 'ntsl'      )
      iostatus = sfrnatt( sd_id, attr_index, ntsl     )
      attr_index = sffattr( sd_id, 'nlog'      )
      iostatus = sfrnatt( sd_id, attr_index, nlog     )
      attr_index = sffattr( sd_id, 'step_res'  )
      iostatus = sfrnatt( sd_id, attr_index, step_res )
      attr_index = sffattr( sd_id, 'step_hdf'  )
      iostatus = sfrnatt( sd_id, attr_index, step_hdf )
      attr_index = sffattr( sd_id, 'last_hdf_time' )
      iostatus = sfrnatt( sd_id, attr_index, last_hdf_time )

! read selected problem dependent parameters
#ifdef SN_SRC
      attr_index = sffattr( sd_id, 'nsn'       )
      iostatus = sfrnatt( sd_id, attr_index, nsn_last    )
#endif /* SN_SRC */

! read variables array
!
      sds_index = sfn2index(sd_id, 'fluid_vars')
      sds_id = sfselect(sd_id, sds_index)
      iostatus = sfrdata(sds_id, istart, stride, dimsu, u)

! read magnetic field
!
      sds_index = sfn2index(sd_id, 'mag_field')
      sds_id = sfselect(sd_id, sds_index)
      iostatus = sfrdata(sds_id, istart, stride, dimsb, b)

#ifdef GRAV
! read gravitational potential
!
      sds_index = sfn2index(sd_id, 'grav_pot')
      sds_id = sfselect(sd_id, sds_index)
      iostatus = sfrdata(sds_id, istart, stride, dims3d, gp)
#endif /* GRAV */

      iostatus = sfendacc(sds_id)

      iostatus = sfend(sd_id)

   end subroutine read_restart

!------------------------------------------------------------------------

   subroutine find_last_restart(restart_number)

      use initproblem, only : problem_name, run_id

      implicit none

      character*120 file_name
      integer restart_number,nres
      logical exist
      character(len=128) :: file_name_base

      restart_number = 0

      write (file_name_base,'(a,a1,a3,a1)') trim(problem_name),'_',run_id,'_'

      call rm_file('restart_list.tmp')

      do nres =999,0,-1
#ifdef HDF5
         write (file_name,'(a,a1,a,a1,a3,a1,i4.4,a4)') &
               trim(cwd),'/',trim(problem_name),'_', run_id,'_',nres,'.res'
#else /* HDF5 */
         write (file_name,'(a,a1,a,a1,a3,a1,3(i2.2,a1),i3.3,a4)') &
               trim(cwd),'/',trim(problem_name),'_', run_id,'_',0,'_',0,'_',0,'_',nres,'.res'
#endif /* HDF5 */
         inquire(file = file_name, exist = exist)
         if(exist) then
            restart_number = nres
         return
         endif
      enddo

   end subroutine find_last_restart

!---------------------------------------------------------------------
!
! writes integrals to text file
!
!---------------------------------------------------------------------
!
   subroutine write_timeslice
      use types
      use fluidindex,      only : nfluid
      use fluidindex,      only : ibx,iby,ibz
      use fluidindex,      only : nvar, iarr_all_dn,iarr_all_mx,iarr_all_my,iarr_all_mz
      use grid,            only : dvol,dx,dy,dz,is,ie,js,je,ks,ke,x,y,z,nxd,nyd,nzd
      use arrays,          only : u,b,wa
      use initproblem,     only : problem_name, run_id
#ifdef IONIZED
      use initionized,     only : gamma_ion, cs_iso_ion,cs_iso_ion2
#endif /* IONIZED */
#ifdef NEUTRAL
      use initneutral,     only : gamma_neu, cs_iso_neu,cs_iso_neu2
#endif /* NEUTRAL */
#ifndef ISO
      use fluidindex,      only : iarr_all_en
#endif /* ISO */
#ifdef COSM_RAYS
      use initcosmicrays,  only : iecr
#endif /* COSM_RAYS */
#ifdef GRAV
      use arrays,          only : gp
#endif /* GRAV */
#ifdef ISO
#ifdef IONIZED
      use initionized,     only : cs_iso_ion2
#endif /* IONIZED */
#ifdef NEUTRAL
      use initneutral,     only : cs_iso_neu2
#endif /* NEUTRAL */
#endif /* ISO */
#ifdef RESISTIVE
      use resistivity
#endif /* RESISTIVE */
#ifdef COSM_RAYS
      use initcosmicrays,  only : iecr
#endif /* COSM_RAYS */
#ifdef SNE_DISTR
      use sndistr,         only : emagadd, tot_emagadd
#endif /* SNE_DISTR */

      implicit none

      character(len=128) :: tsl_file

      real :: mass = 0.0, momx = 0.0, momy = 0.0,  momz = 0.0, &
#ifndef ISO
              ener = 0.0, &
#endif /* ISO */
              ekin = 0.0,  emag = 0.0, &
              tot_mass = 0.0, tot_momx = 0.0, tot_momy = 0.0, tot_momz = 0.0, &
              tot_ener = 0.0, tot_eint = 0.0, tot_ekin = 0.0, tot_emag = 0.0, &
              tot_epot = 0.0, mflx = 0.0, mfly = 0.0, mflz = 0.0, &
              tot_mflx = 0.0, tot_mfly = 0.0, tot_mflz = 0.0

      type(tsl_container) :: tsl

#ifdef GRAV
      integer  :: i,j
      real     :: epot =0.0
#endif /* GRAV */
#ifdef COSM_RAYS
      real     :: encr = 0.0, tot_encr = 0.0
#endif /* COSM_RAYS */
      real     :: cs_iso2
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
#ifdef IONIZED
      cs_iso2 = cs_iso_ion2
#endif /* IONIZED */
#ifdef NEUTRAL
      cs_iso2 = cs_iso_neu2
#else
#ifndef IONIZED
      cs_iso2 = 0.0
#endif
#endif /* NEUTRAL */


      if (proc .eq. 0) then
         write (tsl_file,'(a,a1,a,a1,a3,a1,i3.3,a4)') &
              trim(cwd),'/',trim(problem_name),'_', run_id,'_',nrestart,'.tsl'

         if (tsl_firstcall) then
            open(tsl_lun, file=tsl_file)

            write (tsl_lun, '(a1,a8,50a16)') '#','nstep', 'time', 'timestep',  &
                                               'mass', 'momx', 'momy', 'momz', &
                                               'ener', 'epot', 'eint', 'ekin', &
#ifdef COSM_RAYS
                                               'encr_tot', 'encr_min', 'encr_max',&
#endif /* COSM_RAYS */

#ifdef MAGNETIC
                                               'emag', 'mflx', 'mfly', 'mflz', 'vai_max', &
                                               'b_min', 'b_max', 'divb_max', &
#ifdef RESISTIVE
                                               'eta_max', &
#endif /* RESISTIVE */
#endif /* MAGNETIC */
! some quantities computed in "write_log".One can add more, or change.
#ifdef IONIZED
                                               'vxi_max', 'vyi_max', 'vzi_max', 'csi_max', &
                                               'deni_min', 'deni_max', 'prei_min', 'prei_max', &
#ifndef ISO
                                               'temi_min', 'temi_max', &
#endif /* ISO */
#endif /* IONIZED */

#ifdef NEUTRAL
                                               'denn_min', 'denn_max', 'vxn_max', 'vyn_max', 'vzn_max', &
                                               'pren_min', 'pren_max', 'temn_min', 'temn_max', 'csn_max', &
#endif /* NEUTRAL */

#ifdef DUST
                                              'dend_min', 'dend_max', 'vxd_max', 'vyd_max', 'vzd_max', &
#endif /* DUST */
                                             'dummy'

            write (tsl_lun, '(a1)') '#'
            tsl_firstcall = .false.
         else
            open(tsl_lun, file=tsl_file, position='append')
         endif
      endif

      mass = sum(u(iarr_all_dn,is:ie,js:je,ks:ke)) * dvol
      call mpi_allreduce(mass, tot_mass, 1, mpi_real8, mpi_sum, comm3d, ierr)

      momx = sum(u(iarr_all_mx,is:ie,js:je,ks:ke)) * dvol
      call mpi_allreduce(momx, tot_momx, 1, mpi_real8, mpi_sum, comm3d, ierr)

      momy = sum(u(iarr_all_my,is:ie,js:je,ks:ke)) * dvol
      call mpi_allreduce(momy, tot_momy, 1, mpi_real8, mpi_sum, comm3d, ierr)

      momz = sum(u(iarr_all_mz,is:ie,js:je,ks:ke)) * dvol
      call mpi_allreduce(momz, tot_momz, 1, mpi_real8, mpi_sum, comm3d, ierr)

#ifdef GRAV
      epot = sum(u(iarr_all_dn(1),is:ie,js:je,ks:ke) *gp(is:ie,js:je,ks:ke)) * dvol
      call mpi_allreduce(epot, tot_epot, 1, mpi_real8, mpi_sum, comm3d, ierr)
#endif /* GRAV */

      wa(is:ie,js:je,ks:ke) &
          = 0.5 * (u(iarr_all_mx(1),is:ie,js:je,ks:ke)**2   &
                 + u(iarr_all_my(1),is:ie,js:je,ks:ke)**2   &
                 + u(iarr_all_mz(1),is:ie,js:je,ks:ke)**2)/ &
                   max(u(iarr_all_dn(1),is:ie,js:je,ks:ke),smalld)
      ekin = sum(wa(is:ie,js:je,ks:ke)) * dvol
      call mpi_allreduce(ekin, tot_ekin, 1, mpi_real8, mpi_sum, comm3d, ierr)

      wa(is:ie,js:je,ks:ke) &
         = 0.5 * (b(ibx,is:ie,js:je,ks:ke)**2 + &
                  b(iby,is:ie,js:je,ks:ke)**2 + &
                  b(ibz,is:ie,js:je,ks:ke)**2)
      emag = sum(wa(is:ie,js:je,ks:ke)) * dvol
      call mpi_allreduce(emag, tot_emag, 1, mpi_real8, mpi_sum, comm3d, ierr)

      wa(is:ie,js:je,ks:ke) = b(ibx,is:ie,js:je,ks:ke)
      mflx = sum(wa(is:ie,js:je,ks:ke)) * dy*dz/nxd
      call mpi_allreduce(mflx, tot_mflx, 1, mpi_real8, mpi_sum, comm3d, ierr)

      wa(is:ie,js:je,ks:ke) = b(iby,is:ie,js:je,ks:ke)
      mfly = sum(wa(is:ie,js:je,ks:ke)) * dx*dz/nyd
      call mpi_allreduce(mfly, tot_mfly, 1, mpi_real8, mpi_sum, comm3d, ierr)

      wa(is:ie,js:je,ks:ke) = b(ibz,is:ie,js:je,ks:ke)
      mflz = sum(wa(is:ie,js:je,ks:ke)) * dx*dy/nzd
      call mpi_allreduce(mflz, tot_mflz, 1, mpi_real8, mpi_sum, comm3d, ierr)
#ifdef ISO
      tot_eint = cs_iso2*tot_mass
      tot_ener = tot_eint+tot_ekin+tot_emag
#else /* ISO */
      ener = sum(u(iarr_all_en,is:ie,js:je,ks:ke)) * dvol
      call mpi_allreduce(ener, tot_ener, 1, mpi_real8, mpi_sum, comm3d, ierr)
      tot_eint = tot_ener - tot_ekin - tot_emag
#endif /* ISO */
#ifdef GRAV
      tot_ener = tot_ener + tot_epot
#endif /* GRAV */

#ifdef COSM_RAYS
      encr = sum(u(iecr,is:ie,js:je,ks:ke)) * dvol
      call mpi_allreduce(encr, tot_encr, 1, mpi_real8, mpi_sum, comm3d, ierr)
#endif /* COSM_RAYS */
#ifdef SNE_DISTR
      call mpi_allreduce(emagadd, sum_emagadd, 1, mpi_real8, mpi_sum, comm3d, ierr)
      tot_emagadd = tot_emagadd + sum_emagadd
#endif /* SNE_DISTR */

      call write_log(tsl)

      if (proc == 0) then
         write (tsl_lun, '(1x,i8,50(1x,1pe15.8))') &
                      nstep, &
                      t, dt, tot_mass, &
                      tot_momx, tot_momy, tot_momz, &
#ifdef GRAV
#endif /* GRAV */
                      tot_ener, tot_epot, tot_eint, tot_ekin,&
#ifdef MAGNETIC
                      tot_emag, tot_mflx, tot_mfly, tot_mflz, tsl%vai_max, &
                      tsl%b_min, tsl%b_max, tsl%divb_max, &
#ifdef RESISTIVE
                      tsl%etamax, &
#endif /* RESISTIVE */
#endif /* MAGNETIC */
#ifdef COSM_RAYS
                      tot_encr, tsl%encr_min, tsl%encr_max, &
#endif /* COSM_RAYS */

! some quantities computed in "write_log".One can add more, or change.
#ifdef IONIZED
                      tsl%vxi_max, tsl%vyi_max, tsl%vzi_max, tsl%csi_max, &
                      tsl%deni_min, tsl%deni_max, tsl%prei_min, tsl%prei_max, &
#ifndef ISO
                      tsl%temi_min, tsl%temi_max, &
#endif /* ISO */
#endif /* IONIZED */

#ifdef NEUTRAL
                      tsl%denn_min, tsl%denn_max, tsl%vxn_max, tsl%vyn_max, tsl%vzn_max, &
                      tsl%pren_min, tsl%pren_max, tsl%temn_min, tsl%temn_max, tsl%csn_max, &
#endif /* NEUTRAL */

#ifdef DUST
                      tsl%dend_min, tsl%dend_max, tsl%vxd_max, tsl%vyd_max, tsl%vzd_max, &
#endif /* DUST */
                      0.0
         close(tsl_lun)
      endif

   end subroutine write_timeslice

!---------------------------------------------------------------------
!
! writes timestep diagnostics to the logfile
!
!---------------------------------------------------------------------
!
   subroutine  write_log(tsl)
      use types
      use fluidindex,         only : ibx, iby, ibz, nfluid
      use arrays,             only : wa,u,b
      use grid,               only   : dx,dy,dz,dxmn,nb,is,ie,js,je,ks,ke,nx,ny,nz
      use constants,          only : small, hydro_mass, k_B
      use mpisetup,           only : smallei,cfl
#ifdef IONIZED
      use initionized,        only : gamma_ion, cs_iso_ion,cs_iso_ion2
      use initionized,        only : idni,imxi,imyi,imzi
#ifndef ISO
      use initionized,        only : ieni
#endif /* ISO */
#endif /* IONIZED */
#ifdef NEUTRAL
      use initneutral,        only : gamma_neu, cs_iso_neu,cs_iso_neu2
      use initneutral,        only : idnn,imxn,imyn,imzn
#ifndef ISO
      use initneutral,        only : ienn
#endif /* ISO */
#endif /* NEUTRAL */
#ifdef DUST
      use initdust,           only : idnd,imxd,imyd,imzd
#endif /* DUST */
#ifdef COSM_RAYS
      use timestepcosmicrays, only : dt_crs
      use initcosmicrays,     only : iecr
#endif /* COSM_RAYS */
#ifdef RESISTIVE
      use resistivity
#endif /* RESISTIVE */

      implicit none
      type(tsl_container), optional :: tsl

#ifdef MAGNETIC
      type(value) :: b_min, b_max, divb_max
#endif /* MAGNETIC */

#ifdef IONIZED
      type(value) :: deni_min, deni_max, vxi_max, vyi_max, vzi_max, &
                     prei_min, prei_max, temi_min, temi_max, vai_max, csi_max
#endif /* IONIZED */

#ifdef NEUTRAL
      type(value) :: denn_min, denn_max, vxn_max, vyn_max, vzn_max, &
                     pren_min, pren_max, temn_min, temn_max, csn_max
#endif /* NEUTRAL */

#ifdef DUST
      type(value) :: dend_min, dend_max, vxd_max, vyd_max, vzd_max
#endif /* DUST */

#ifdef COSM_RAYS
      type(value) :: encr_min, encr_max
#endif /* COSM_RAYS */

#ifdef RESISTIVE
      type(value) :: etamax
#endif /* RESISTIVE */

! Timestep diagnostics
#ifdef NEUTRAL
      wa            = u(idnn,:,:,:)
      denn_min%val  = minval(wa(is:ie,js:je,ks:ke))
      denn_min%loc  = minloc(wa(is:ie,js:je,ks:ke)) &
                    + (/nb,nb,nb/)
      call mpifind(denn_min%val, 'min', denn_min%loc, denn_min%proc)

      denn_max%val  = maxval(wa(is:ie,js:je,ks:ke))
      denn_max%loc  = maxloc(wa(is:ie,js:je,ks:ke)) &
                    + (/nb,nb,nb/)
      call mpifind(denn_max%val, 'max', denn_max%loc, denn_max%proc)

      wa          = abs(u(imxn,:,:,:)/u(idnn,:,:,:))
      vxn_max%val = maxval(wa(is:ie,js:je,ks:ke))
      vxn_max%loc = maxloc(wa(is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)
      call mpifind(vxn_max%val, 'max', vxn_max%loc, vxn_max%proc)

      wa          = abs(u(imyn,:,:,:)/u(idnn,:,:,:))
      vyn_max%val = maxval(wa(is:ie,js:je,ks:ke))
      vyn_max%loc = maxloc(wa(is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)
      call mpifind(vyn_max%val, 'max', vyn_max%loc, vyn_max%proc)

      wa           = abs(u(imzn,:,:,:)/u(idnn,:,:,:))
      vzn_max%val  = maxval(wa(is:ie,js:je,ks:ke))
      vzn_max%loc  = maxloc(wa(is:ie,js:je,ks:ke)) &
                   + (/nb,nb,nb/)
      call mpifind(vzn_max%val, 'max', vzn_max%loc, vzn_max%proc)
#ifdef ISO
      pren_min%val  = cs_iso_neu2*denn_min%val
      pren_min%loc  = denn_min%loc
      pren_min%proc = denn_min%proc
      pren_max%val  = cs_iso_neu2*denn_max%val
      pren_max%loc  = denn_max%loc
      pren_max%proc = denn_max%proc
      csn_max%val   = cs_iso_neu
      csn_max%loc   = 0
      csn_max%proc  = 0
      temn_min%val  = hydro_mass / k_B * cs_iso_neu2
      temn_min%loc  = 0
      temn_min%proc = 0
      temn_max%val  = hydro_mass / k_B * cs_iso_neu2
      temn_max%loc  = 0
      temn_max%proc = 0
#else /* ISO */
      wa(:,:,:) = (u(ienn,:,:,:) &                ! eint
                - 0.5*((u(imxn,:,:,:)**2 +u(imyn,:,:,:)**2 &
                      + u(imzn,:,:,:)**2)/u(idnn,:,:,:)))
      wa(:,:,:) = max(wa(:,:,:),smallei)
      wa(:,:,:) = (gamma_neu-1.0)*wa(:,:,:)           ! pres

      pren_min%val  = minval(wa(is:ie,js:je,ks:ke))
      pren_min%loc  = minloc(wa(is:ie,js:je,ks:ke)) + (/nb,nb,nb/)
      call mpifind(pren_min%val, 'min', pren_min%loc, pren_min%proc)

      pren_max%val  = maxval(wa(is:ie,js:je,ks:ke))
      pren_max%loc  = maxloc(wa(is:ie,js:je,ks:ke)) + (/nb,nb,nb/)
      call mpifind(pren_max%val, 'max', pren_max%loc, pren_max%proc)

      temn_max%val  = maxval( hydro_mass / k_B * wa(is:ie,js:je,ks:ke) &
                                             /u(idnn,is:ie,js:je,ks:ke))
      temn_max%loc  = maxloc(wa(is:ie,js:je,ks:ke)    &
                        /u(idnn,is:ie,js:je,ks:ke)  ) + (/nb,nb,nb/)
      call mpifind(temn_max%val, 'max', temn_max%loc, temn_max%proc)

      temn_min%val  = minval( hydro_mass / k_B * wa(is:ie,js:je,ks:ke) &
                                            /u(idnn,is:ie,js:je,ks:ke))
      temn_min%loc  = minloc(wa(is:ie,js:je,ks:ke) &
                        /u(idnn,is:ie,js:je,ks:ke)  ) &
                     + (/nb,nb,nb/)
      call mpifind(temn_min%val, 'min', temn_min%loc, temn_min%proc)

      wa(:,:,:) = gamma_neu*wa(:,:,:)
      csn_max%val    = sqrt(maxval(wa(is:ie,js:je,ks:ke) &
                              /u(idnn,is:ie,js:je,ks:ke)))
      csn_max%loc    = maxloc(wa(is:ie,js:je,ks:ke) &
                         /u(idnn,is:ie,js:je,ks:ke)) &
                     + (/nb,nb,nb/)
      call mpifind(csn_max%val, 'max', csn_max%loc, csn_max%proc)
#endif /* ISO */

#endif /* NEUTRAL */

#ifdef IONIZED
      wa            = u(idni,:,:,:)
      deni_min%val  = minval(wa(is:ie,js:je,ks:ke))
      deni_min%loc  = minloc(wa(is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)
      call mpifind(deni_min%val, 'min', deni_min%loc, deni_min%proc)

      deni_max%val  = maxval(wa(is:ie,js:je,ks:ke))
      deni_max%loc  = maxloc(wa(is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)
      call mpifind(deni_max%val, 'max', deni_max%loc, deni_max%proc)

      wa          = abs(u(imxi,:,:,:)/u(idni,:,:,:))
      vxi_max%val = maxval(wa(is:ie,js:je,ks:ke))
      vxi_max%loc = maxloc(wa(is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)
      call mpifind(vxi_max%val, 'max', vxi_max%loc, vxi_max%proc)

      wa          = abs(u(imyi,:,:,:)/u(idni,:,:,:))
      vyi_max%val = maxval(wa(is:ie,js:je,ks:ke))
      vyi_max%loc = maxloc(wa(is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)
      call mpifind(vyi_max%val, 'max', vyi_max%loc, vyi_max%proc)

      wa           = abs(u(imzi,:,:,:)/u(idni,:,:,:))
      vzi_max%val  = maxval(wa(is:ie,js:je,ks:ke))
      vzi_max%loc  = maxloc(wa(is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)
      call mpifind(vzi_max%val, 'max', vzi_max%loc, vzi_max%proc)

#ifdef MAGNETIC
      wa(:,:,:)  = b(1,:,:,:)*b(1,:,:,:) + b(2,:,:,:)*b(2,:,:,:) + &
                   b(3,:,:,:)*b(3,:,:,:)
      b_min%val  = sqrt(minval(wa(is:ie,js:je,ks:ke)))
      b_min%loc  = minloc(wa(is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)
      call mpifind(b_min%val, 'min', b_min%loc, b_min%proc)

      b_max%val  = sqrt(maxval(wa(is:ie,js:je,ks:ke)))
      b_max%loc  = maxloc(wa(is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)
      call mpifind(b_max%val, 'max', b_max%loc, b_max%proc)

      vai_max%val = sqrt(maxval(wa(is:ie,js:je,ks:ke) &
                           /u(idni,is:ie,js:je,ks:ke)))
      vai_max%loc = maxloc(wa(is:ie,js:je,ks:ke)     &
                      /u(idni,is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)
      call mpifind(vai_max%val, 'max', vai_max%loc, vai_max%proc)
#endif /* MAGNETIC */

#ifdef ISO
      prei_min%val  = cs_iso_ion2*deni_min%val
      prei_min%loc  = deni_min%loc
      prei_min%proc = deni_min%proc
      prei_max%val  = cs_iso_ion2*deni_max%val
      prei_max%loc  = deni_max%loc
      prei_max%proc = deni_max%proc
      csi_max%val   = cs_iso_ion
      csi_max%loc   = 0
      csi_max%proc  = 0
      temi_min%val  = hydro_mass / k_B * cs_iso_ion2
      temi_min%loc  = 0
      temi_min%proc = 0
      temi_max%val  = hydro_mass / k_B * cs_iso_ion2
      temi_max%loc  = 0
      temi_max%proc = 0
#else /* ISO */
      wa(:,:,:) = (u(ieni,:,:,:) &                ! eint
                - 0.5*((u(imxi,:,:,:)**2 +u(imyi,:,:,:)**2 &
                  + u(imzi,:,:,:)**2)/u(idni,:,:,:)))
#ifdef MAGNETIC
      wa(:,:,:) = wa(:,:,:) - 0.5*(b(ibx,:,:,:)**2 + b(iby,:,:,:)**2 + &
                 b(ibz,:,:,:)**2)
#endif /* MAGNETIC */
      wa(:,:,:) = max(wa(:,:,:),smallei)
      wa(:,:,:) = (gamma_ion-1.0)*wa(:,:,:)           ! pres

      prei_min%val  = minval(wa(is:ie,js:je,ks:ke))
      prei_min%loc  = minloc(wa(is:ie,js:je,ks:ke)) + (/nb,nb,nb/)
      call mpifind(prei_min%val, 'min', prei_min%loc, prei_min%proc)

      prei_max%val  = maxval(wa(is:ie,js:je,ks:ke))
      prei_max%loc  = maxloc(wa(is:ie,js:je,ks:ke)) + (/nb,nb,nb/)
      call mpifind(prei_max%val, 'max', prei_max%loc, prei_max%proc)

      temi_max%val  = maxval( hydro_mass / k_B * wa(is:ie,js:je,ks:ke) &
                                            /u(idni,is:ie,js:je,ks:ke))
      temi_max%loc  = maxloc(wa(is:ie,js:je,ks:ke)    &
                        /u(idni,is:ie,js:je,ks:ke)  ) + (/nb,nb,nb/)
      call mpifind(temi_max%val, 'max', temi_max%loc, temi_max%proc)

      temi_min%val  = minval( hydro_mass / k_B * wa(is:ie,js:je,ks:ke) &
                                            /u(idni,is:ie,js:je,ks:ke))
      temi_min%loc  = minloc(wa(is:ie,js:je,ks:ke) &
                        /u(idni,is:ie,js:je,ks:ke)  ) &
                     + (/nb,nb,nb/)
      call mpifind(temi_min%val, 'min', temi_min%loc, temi_min%proc)

      wa(:,:,:) = gamma_ion*wa(:,:,:)
      csi_max%val    = sqrt(maxval(wa(is:ie,js:je,ks:ke) &
                              /u(idni,is:ie,js:je,ks:ke)))
      csi_max%loc    = maxloc(wa(is:ie,js:je,ks:ke) &
                         /u(idni,is:ie,js:je,ks:ke)) &
                     + (/nb,nb,nb/)
      call mpifind(csi_max%val, 'max', csi_max%loc, csi_max%proc)
#endif /* ISO */

#endif /* IONIZED */

#ifdef DUST
      wa            = u(idnd,:,:,:)
      dend_min%val  = minval(wa(is:ie,js:je,ks:ke))
      dend_min%loc  = minloc(wa(is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)
      call mpifind(dend_min%val, 'min', dend_min%loc, dend_min%proc)

      dend_max%val  = maxval(wa(is:ie,js:je,ks:ke))
      dend_max%loc  = maxloc(wa(is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)
      call mpifind(dend_max%val, 'max', dend_max%loc, dend_max%proc)

      where(u(idnd,:,:,:) > 0.0)
         wa          = abs(u(imxd,:,:,:)/u(idnd,:,:,:))
      elsewhere
         wa          = 0.0
      endwhere
      vxd_max%val = maxval(wa(is:ie,js:je,ks:ke))
      vxd_max%loc = maxloc(wa(is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)
      call mpifind(vxd_max%val, 'max', vxd_max%loc, vxd_max%proc)

      where(u(idnd,:,:,:) > 0.0)
         wa          = abs(u(imyd,:,:,:)/u(idnd,:,:,:))
      elsewhere
         wa          = 0.0
      endwhere
      vyd_max%val = maxval(wa(is:ie,js:je,ks:ke))
      vyd_max%loc = maxloc(wa(is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)
      call mpifind(vyd_max%val, 'max', vyd_max%loc, vyd_max%proc)

      where(u(idnd,:,:,:) > 0.0)
         wa           = abs(u(imzd,:,:,:)/u(idnd,:,:,:))
      elsewhere
         wa          = 0.0
      endwhere
      vzd_max%val  = maxval(wa(is:ie,js:je,ks:ke))
      vzd_max%loc  = maxloc(wa(is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)
      call mpifind(vzd_max%val, 'max', vzd_max%loc, vzd_max%proc)
#endif /* DUST */

#ifdef RESISTIVE
      etamax%val = eta_max
      call mpifind(etamax%val, 'max', etamax%loc, etamax%proc)
#endif /* RESISTIVE */


#ifdef MAGNETIC
      wa(1:nx-1,1:ny-1,1:max(nz-1,1)) = &
                 (b(ibx,2:nx,1:ny-1,1:max(nz-1,1)) - b(ibx,1:nx-1,1:ny-1,1:max(nz-1,1)))*dy*dz &
                +(b(iby,1:nx-1,2:ny,1:max(nz-1,1)) - b(iby,1:nx-1,1:ny-1,1:max(nz-1,1)))*dx*dz &
                +(b(ibz,1:nx-1,1:ny-1,min(2,nz):nz) - b(ibz,1:nx-1,1:ny-1,1:max(nz-1,1)))*dx*dy
      wa = abs(wa)

      wa(ie,:,:) = wa(max(ie-1,1),:,:)
      wa(:,je,:) = wa(:,max(je-1,1),:)
      wa(:,:,ke) = wa(:,:,max(ke-1,1))

      divb_max%val  = maxval(wa(is:ie,js:je,ks:ke))
      divb_max%loc  = maxloc(wa(is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)
      call mpifind(divb_max%val, 'max', divb_max%loc, divb_max%proc)
#endif /* MAGNETIC */

#ifdef COSM_RAYS
      wa            = u(iecr,:,:,:)
      encr_min%val  = minval(wa(is:ie,js:je,ks:ke))
      encr_min%loc  = minloc(wa(is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)
      call mpifind(encr_min%val, 'min', encr_min%loc, encr_min%proc)

      encr_max%val  = maxval(wa(is:ie,js:je,ks:ke))
      encr_max%loc  = maxloc(wa(is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)
      call mpifind(encr_max%val, 'max', encr_max%loc, encr_max%proc)
#endif /* COSM_RAYS */

      if(proc == 0)  then
         if(.not.present(tsl)) then
            open(log_lun, file=log_file, position='append')
#ifdef IONIZED
            write(log_lun,771) 'min(dens)   ION  =', deni_min%val,  deni_min%proc,  deni_min%loc
            write(log_lun,771) 'max(dens)   ION  =', deni_max%val,  deni_max%proc,  deni_max%loc
#ifndef ISO
            write(log_lun,771) 'min(temp)   ION  =', temi_min%val,  temi_min%proc,  temi_min%loc
            write(log_lun,771) 'max(temp)   ION  =', temi_max%val,  temi_max%proc,  temi_max%loc
#endif /* ISO */
            write(log_lun,771) 'min(pres)   ION  =', prei_min%val,  prei_min%proc,  prei_min%loc
            write(log_lun,771) 'max(pres)   ION  =', prei_max%val,  prei_max%proc,  prei_max%loc
            write(log_lun,777) 'max(|vx|)   ION  =', vxi_max%val, 'dt=',cfl*dx/(vxi_max%val+small),   vxi_max%proc, vxi_max%loc
            write(log_lun,777) 'max(|vy|)   ION  =', vyi_max%val, 'dt=',cfl*dy/(vyi_max%val+small),   vyi_max%proc, vyi_max%loc
            write(log_lun,777) 'max(|vz|)   ION  =', vzi_max%val, 'dt=',cfl*dz/(vzi_max%val+small),   vzi_max%proc, vzi_max%loc
            write(log_lun,777) 'max(c_s )   ION  =', csi_max%val, 'dt=',cfl*dxmn/(csi_max%val+small), csi_max%proc, csi_max%loc
#ifdef MAGNETIC
            write(log_lun,777) 'max(c_f)    ION  =', sqrt(csi_max%val**2+vai_max%val**2),&
                                      'dt=',cfl*dxmn/sqrt(csi_max%val**2+vai_max%val**2)
            write(log_lun,777) 'max(v_a)    ION  =', vai_max%val, 'dt=',cfl*dxmn/(vai_max%val+small), vai_max%proc, vai_max%loc
            write(log_lun,771) 'min(|b|)    MAG  =', b_min%val,     b_min%proc,     b_min%loc
            write(log_lun,771) 'max(|b|)    MAG  =', b_max%val,     b_max%proc,     b_max%loc
            write(log_lun,771) 'max(|divb|) MAG  =', divb_max%val,  divb_max%proc,  divb_max%loc
#else /* MAGNETIC */
            write(log_lun,777) 'max(c_s)    ION  =', sqrt(csi_max%val**2), 'dt=',cfl*dxmn/sqrt(csi_max%val**2)
#endif /* MAGNETIC */
#endif /* IONIZED */

#ifdef NEUTRAL
            write(log_lun,771) 'min(dens)   NEU  =', denn_min%val,  denn_min%proc,  denn_min%loc
            write(log_lun,771) 'max(dens)   NEU  =', denn_max%val,  denn_max%proc,  denn_max%loc
#ifndef ISO
            write(log_lun,771) 'min(temp)   NEU  =', temn_min%val,  temn_min%proc,  temn_min%loc
            write(log_lun,771) 'max(temp)   NEU  =', temn_max%val,  temn_max%proc,  temn_max%loc
#endif /* ISO */
            write(log_lun,771) 'min(pres)   NEU  =', pren_min%val,  pren_min%proc,  pren_min%loc
            write(log_lun,771) 'max(pres)   NEU  =', pren_max%val,  pren_max%proc,  pren_max%loc
            write(log_lun,777) 'max(|vx|)   NEU  =', vxn_max%val, 'dt=',cfl*dx/(vxn_max%val+small),   vxn_max%proc, vxn_max%loc
            write(log_lun,777) 'max(|vy|)   NEU  =', vyn_max%val, 'dt=',cfl*dy/(vyn_max%val+small),   vyn_max%proc, vyn_max%loc
            write(log_lun,777) 'max(|vz|)   NEU  =', vzn_max%val, 'dt=',cfl*dz/(vzn_max%val+small),   vzn_max%proc, vzn_max%loc
            write(log_lun,777) 'max(c_s )   NEU  =', csn_max%val, 'dt=',cfl*dxmn/(csn_max%val+small), csn_max%proc, csn_max%loc
#endif /* NEUTRAL */
#ifdef DUST
            write(log_lun,771) 'min(dens)   DST  =', dend_min%val,  dend_min%proc,  dend_min%loc
            write(log_lun,771) 'max(dens)   DST  =', dend_max%val,  dend_max%proc,  dend_max%loc
            write(log_lun,777) 'max(|vx|)   DST  =', vxd_max%val, 'dt=',cfl*dx/(vxd_max%val+small),   vxd_max%proc, vxd_max%loc
            write(log_lun,777) 'max(|vy|)   DST  =', vyd_max%val, 'dt=',cfl*dy/(vyd_max%val+small),   vyd_max%proc, vyd_max%loc
            write(log_lun,777) 'max(|vz|)   DST  =', vzd_max%val, 'dt=',cfl*dz/(vzd_max%val+small),   vzd_max%proc, vzd_max%loc
#endif /* DUST */
#ifdef COSM_RAYS
            write(log_lun,771) 'min(encr)   CRS  =', encr_min%val,        encr_min%proc, encr_min%loc
            write(log_lun,777) 'max(encr)   CRS  =', encr_max%val,      'dt=',dt_crs,     encr_max%proc, encr_max%loc
#endif /* COSM_RAYS */
#ifdef RESISTIVE
            write(log_lun,777) 'max(eta)    RES  =', etamax%val ,      'dt=',dt_resist, etamax%proc,  etamax%loc
#endif /* RESISTIVE */

            write(log_lun,'(a80)') '================================================================================'
            write(log_lun,900) nstep,dt,t
            write(log_lun,'(a80)') '================================================================================'
            close(log_lun)
         else
#ifdef IONIZED
            tsl%deni_min = deni_min%val
            tsl%deni_max = deni_max%val
            tsl%vxi_max  = vxi_max%val
            tsl%vyi_max  = vyi_max%val
            tsl%vzi_max  = vzi_max%val
            tsl%prei_min = prei_min%val
            tsl%prei_max = prei_max%val
            tsl%temi_min = temi_min%val
            tsl%temi_max = temi_max%val
            tsl%vai_max  = vai_max%val
            tsl%csi_max  = csi_max%val
#endif /* IONIZED */
#ifdef MAGNETIC
            tsl%b_min = b_min%val
            tsl%b_max = b_max%val
            tsl%divb_max = divb_max%val
#endif /* MAGNETIC */
#ifdef NEUTRAL
            tsl%denn_min = denn_min%val
            tsl%denn_max = denn_max%val
            tsl%vxn_max  = vxn_max%val
            tsl%vyn_max  = vyn_max%val
            tsl%vzn_max  = vzn_max%val
            tsl%pren_min = pren_min%val
            tsl%pren_max = pren_max%val
            tsl%temn_min = temn_min%val
            tsl%temn_max = temn_max%val
            tsl%csn_max  = csn_max%val
#endif /* NEUTRAL */

#ifdef DUST
            tsl%dend_min = dend_min%val
            tsl%dend_max = dend_max%val
            tsl%vxd_max = vxd_max%val
            tsl%vyd_max = vyd_max%val
            tsl%vzd_max = vzd_max%val
#endif /* DUST */

#ifdef COSM_RAYS
            tsl%encr_min = encr_min%val
            tsl%encr_max = encr_max%val
#endif /* COSM_RAYS */

#ifdef RESISTIVE
            tsl%etamax = etamax%val
#endif /* RESISTIVE */
         endif
      endif

771 format(5x,a18,(1x,e15.9),16x,5(1x,i4))
777 format(5x,a18,(1x,e15.9),2x,a3,(1x,e10.4),5(1x,i4))

900 format('   nstep = ',i7,'   dt = ',f22.16,'   t = ',f22.16,2(1x,i4))
#ifdef RESISTIVE
776 format(5x,a18,(1x,e10.4),2x,a3,(1x,e10.4),4(1x,i4))
#endif /* RESISTIVE */

!#endif /* NOT_WORKING */

   end subroutine write_log


!
!=======================================================================
!
!                    B E G I N   S U B R O U T I N E
!                             R E A D F M S G
!
!=======================================================================
!
   subroutine read_file_msg

!     written by: michal hanasz
!     date:       26. june 2003
!
!-------------------------------------------------------------------------
!     configurable parameters: problem.par
!-------------------------------------------------------------------------
!      user_message_file           ! 1st (user) message file (eg.'./msg')
!      system_message_file         ! 2nd (ups)  message file (eg.'/etc/ups/user/msg')
!-------------------------------------------------------------------------
      implicit none
      character user_last_msg_file*80
      character system_last_msg_file*80

      character user_msg_time(10)*80,system_msg_time(10)*80

      character, save :: user_msg_time_old*80,system_msg_time_old*80
      character(len=265) :: syscom
      integer i
      integer(kind=1) :: system

      msg=''
      user_msg_time(9)=''
      nchar=0

      open(91,file=user_message_file,status='old',err=224)
      read(91,fmt=*,err=224,end=224) msg, msg_param

      close(91)
      goto 225
224   continue
      close(91)

      open(92,file=user_message_file,status='old',err=333)
      read(92,fmt=*,err=333,end=333) msg
      close(92)
      goto 225
333   continue
      close(92)
      goto 888

225   continue
      nchar=len_trim(msg)

      user_last_msg_file='./user_last_msg.tmp'
      call rm_file(user_last_msg_file)

      syscom='ls -l --full-time msg >'//user_last_msg_file
      scstatus = system(syscom)
      open(93,file=user_last_msg_file,status='old',err=888)
        read(93,fmt=*,err=888,end=888) (user_msg_time(i), i=1,10)
      close(93)

!---  do the requested action only once for a given user message file

      if(user_msg_time(7) .eq. user_msg_time_old) then
         msg=''
         nchar=0
      else
         system_msg_time_old= user_msg_time(7)
         return
      endif

888   continue

      call rm_file(user_message_file)

!------------------------------------------------------------------------


      open(96,file=system_message_file,status='old',err=424)
      read(96,fmt=*,err=424,end=424) msg, msg_param
      close(96)
      goto 425
424   continue
      close(96)

      open(97,file=system_message_file,status='old',err=533)
      read(97,fmt=*,err=533,end=533) msg
      close(97)
      goto 425
533   continue
      close(97)
      goto 999

425   continue
      nchar=len_trim(msg)

      system_last_msg_file='./system_last_msg.tmp'
      call rm_file(system_last_msg_file)

      syscom='ls -l --full-time '//system_message_file//' > '//system_last_msg_file
      scstatus = system(syscom)

      open(98,file=system_last_msg_file,status='old',err=999)
        read(98,fmt=*,err=999,end=999) (system_msg_time(i), i=1,10)
      close(98)

!---  do the requested action only once for a given system message file

      if(system_msg_time(7) .eq. system_msg_time_old) then
         msg=''
         nchar=0
      else
         system_msg_time_old= system_msg_time(7)
         return
      endif

999   continue


      return
   end subroutine read_file_msg

!------------------------------------------------------------------------

   subroutine rm_file(file_name)
      implicit none
      character*(*) file_name
      character(len=265) syscom
      logical exist
      integer(kind=1) :: system

!     delete file if exists
      inquire(file = file_name, exist = exist)
      if(exist) then
         syscom='rm -f '//trim(file_name)
         scstatus = system(syscom)
      endif

   end subroutine rm_file

end module dataio
