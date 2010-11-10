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
!! \copydetails dataio::init_dataio
!!
!! \todo check and if necessary bring back usefulness of min_disk_space_MB parameter
!<

module dataio
!>
!! \brief (KK)
!!
!! \todo check the usefulness of wait logical variable
!<
   use dataio_pub,    only: cwdlen, hnlen, varlen, domain, fmin, fmax, vizit, nend, tend, wend, nrestart
   use mpisetup,      only: cbuff_len
   use types,         only: idlen

   implicit none

   private
   public :: check_log, check_tsl, write_data, write_crashed, cleanup_dataio, init_dataio, user_msg_handler

   integer               :: istep                  !< current number of substep (related to integration order)

   integer, parameter    :: nvarsmx = 16           !< maximum number of variables to dump in hdf files
   character(len=idlen)  :: new_id                 !< three character string to change run_id when restarting simulation (e.g. to avoid overwriting of the output from the previous (pre-restart) simulation; if new_id = '' then run_id is still used)
   character(len=cbuff_len) :: restart             !< choice of restart file: if restart = 'last': automatic choice of the last restart file regardless of "nrestart" value; if something else is set: "nrestart" value is fixing
   character(len=cbuff_len) :: mag_center          !< choice to dump magnetic fields values from cell centers or not (if not then values from cell borders)
   integer               :: resdel                 !< number of recent restart dumps which should be saved; each n-resdel-1 restart file is supposed to be deleted while writing n restart file
   real                  :: dt_hdf                 !< time between successive hdf dumps
   real                  :: dt_res                 !< time between successive restart file dumps
   real                  :: dt_tsl                 !< time between successive timeslice dumps
   real                  :: dt_log                 !< time between successive log dumps
   real                  :: dt_plt                 !< time between successive domain slices files dumps
   integer               :: min_disk_space_MB      !< minimum disk space in MB to continue simulation <b>(currently not used)</b>
   integer               :: sleep_minutes          !< minutes of sleeping time before continue simulation
   integer               :: sleep_seconds          !< seconds of sleeping time before continue simulation
   character(len=cwdlen) :: user_message_file      !< path to possible user message file containing dt_xxx changes or orders to dump/stop/end simulation
   character(len=cwdlen) :: system_message_file    !< path to possible system (UPS) message file containing orders to dump/stop/end simulation
   integer               :: ix                     !< index in x-direction of slice to dump in plt files
   integer               :: iy                     !< index in y-direction of slice to dump in plt files
   integer               :: iz                     !< index in z-direction of slice to dump in plt files
   integer               :: iv                     !< work index to count successive variables to dump in hdf files
   character(len=varlen), dimension(nvarsmx) :: vars !< array of 4-character strings standing for variables to dump in hdf files

   integer, parameter    :: tsl_lun = 2            !< luncher for timeslice file
   integer               :: step_res               !< number of simulation timestep corresponding to values dumped in restart file
   integer               :: nhdf_start             !< number of hdf file for the first hdf dump in simulation run
   integer               :: nres_start             !< number of restart file for the first restart dump in simulation run
   real                  :: t_start                !< time in simulation of start simulation run
   logical               :: tsl_firstcall          !< logical value to start a new timeslice file

   logical               :: wait                   !< logical value to have a break in simulation (currently not used)

   integer               :: nchar                  !< number of characters in a user/system message
   integer, parameter    :: umsg_len = 16
   character(len=umsg_len) :: umsg                    !< string of characters - content of a user/system message
   real                  :: umsg_param              !< parameter changed by a user/system message

   character(len=hnlen)  :: hostfull, host
   integer               :: pid
   integer               :: uid
   integer               :: ihost
   integer               :: scstatus

   character(len=cwdlen) :: filename               !< string of characters indicating currently used file

   namelist /END_CONTROL/ nend, tend, wend
   namelist /RESTART_CONTROL/ restart, new_id, nrestart, resdel
   namelist /OUTPUT_CONTROL/ dt_hdf, dt_res, dt_tsl, dt_log, dt_plt, ix, iy, iz, &
                             domain, vars, mag_center, vizit, fmin, fmax, &
                             min_disk_space_MB, sleep_minutes, sleep_seconds, &
                             user_message_file, system_message_file

   interface mpi_addmul
      module procedure mpi_sum4d_and_multiply
      module procedure mpi_sum3d_and_multiply
   end interface mpi_addmul

   contains

      subroutine check_log

         use mpisetup,      only: t
         use dataio_pub,    only: nlog

         implicit none

         if (dt_log > 0.0) then
            if (nlog < (int(t / dt_log) + 1)) then
               call write_log
               nlog = nlog + 1
            endif
         endif

      end subroutine check_log

      subroutine check_tsl

         use mpisetup,      only: t
         use dataio_pub,    only: ntsl

         implicit none

         if (dt_tsl .gt. 0.0) then
            if (ntsl .lt. (int(t / dt_tsl) + 1)) then
               call write_timeslice
               ntsl = ntsl + 1
             endif
         endif
      end subroutine check_tsl

!---------------------------------------------------------------------
!
! initializes dataio parameters
!
!---------------------------------------------------------------------
!
!>
!! \brief Routine that sets the initial values of data input/output parameters from namelists @c END_CONTROL, @c RESTART_CONTROL and @c OUTPUT_CONTROL
!!
!! \n \n
!! @b END_CONTROL
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>nend</td><td>1    </td><td>integer</td><td>\copydoc dataio::nend</td></tr>
!! <tr><td>tend</td><td>-1.0 </td><td>real   </td><td>\copydoc dataio::tend</td></tr>
!! </table>
!! \n \n
!! @b RESTART_CONTROL
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>restart </td><td>'last'</td><td>'last' or another string of characters</td><td>\copydoc dataio::restart </td></tr>
!! <tr><td>new_id  </td><td>''    </td><td>string of characters                  </td><td>\copydoc dataio::new_id  </td></tr>
!! <tr><td>nrestart</td><td>3     </td><td>integer                               </td><td>\copydoc dataio::nrestart</td></tr>
!! <tr><td>resdel  </td><td>0     </td><td>integer                               </td><td>\copydoc dataio::resdel  </td></tr>
!! </table>
!! \n \n
!! @b OUTPUT_CONTROL
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>dt_hdf             </td><td>0.0                </td><td>real      </td><td>\copydoc dataio::dt_hdf           </td></tr>
!! <tr><td>dt_res             </td><td>0.0                </td><td>real      </td><td>\copydoc dataio::dt_res           </td></tr>
!! <tr><td>dt_tsl             </td><td>0.0                </td><td>real      </td><td>\copydoc dataio::dt_tsl           </td></tr>
!! <tr><td>dt_log             </td><td>0.0                </td><td>real      </td><td>\copydoc dataio::dt_log           </td></tr>
!! <tr><td>dt_plt             </td><td>0.0                </td><td>real      </td><td>\copydoc dataio::dt_plt           </td></tr>
!! <tr><td>ix                 </td><td>                   </td><td>integer   </td><td>\copydoc dataio::ix               </td></tr>
!! <tr><td>iy                 </td><td>                   </td><td>integer   </td><td>\copydoc dataio::iy               </td></tr>
!! <tr><td>iz                 </td><td>                   </td><td>integer   </td><td>\copydoc dataio::iz               </td></tr>
!! <tr><td>domain             </td><td>'phys_domain'     </td><td>'phys_domain' or 'full_domain'                     </td><td>\copydoc dataio::domain</td></tr>
!! <tr><td>vars               </td><td>''                 </td><td>'dens', 'velx', 'vely', 'velz', 'ener' and some more </td><td>\copydoc dataio::vars  </td></tr>
!! <tr><td>mag_center         </td><td>'no'               </td><td>'yes'/'no'</td><td>\copydoc dataio::mag_center       </td></tr>
!! <tr><td>min_disk_space_MB  </td><td>100                </td><td>integer   </td><td>\copydoc dataio::min_disk_space_MB</td></tr>
!! <tr><td>sleep_minutes      </td><td>0                  </td><td>integer   </td><td>\copydoc dataio::sleep_minutes    </td></tr>
!! <tr><td>sleep_seconds      </td><td>0                  </td><td>integer   </td><td>\copydoc dataio::sleep_seconds    </td></tr>
!! <tr><td>user_message_file  </td><td>trim(cwd)//'/msg'  </td><td>string similar to default value              </td><td>\copydoc dataio::user_message_file  </td></tr>
!! <tr><td>system_message_file</td><td>'/tmp/piernik_msg'</td><td>string of characters similar to default value</td><td>\copydoc dataio::system_message_file</td></tr>
!! </table>
!! \n \n
!<
   subroutine init_dataio

      use dataio_hdf5,     only: init_hdf5, read_restart_hdf5, parfile, parfilelines
      use dataio_pub,      only: chdf, nres, last_hdf_time, step_hdf, nlog, ntsl, dataio_initialized, log_file, cwdlen, maxparfilelines, cwd, &
           &                     tmp_log_file, msglen, printinfo, warn, nhdf, nstep_start, set_container_chdf, get_container
      use dataio_pub,      only: par_file, ierrh, namelist_errh, compare_namelist  ! QA_WARN required for diff_nml
      use fluidboundaries, only: all_fluid_boundaries
      use mpisetup,        only: lbuff, ibuff, rbuff, cbuff, proc, cbuff_len, comm, ierr, buffer_dim, &
           &                      MPI_CHARACTER, MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_LOGICAL, &
           &                      t, nstep, bnd_xl, bnd_xr, bnd_yl, bnd_yr, bnd_zl, bnd_zr
      use problem_pub,     only: problem_name, run_id
      use timer,           only: time_left
      use version,         only: nenv,env, init_version
#ifdef MAGNETIC
      use magboundaries,   only: all_mag_boundaries
#endif /* MAGNETIC */

      implicit none

      logical              :: tn
      integer(kind=1)      :: getpid
      integer(kind=1)      :: hostnm
      integer(kind=1)      :: system
      integer              :: system_status, i
      character(LEN=msglen):: system_command

      restart = 'last'   ! 'last': automatic choice of the last restart file
                         ! regardless of "nrestart" value;
                         ! if something else is set: "nrestart" value is fixing
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
      write(user_message_file,'(a,"/msg")') trim(cwd)
      system_message_file = "/tmp/piernik_msg"

      wait  = .false.
      tsl_firstcall = .true.

      nhdf  = 0
      ntsl  = 0
      nres  = 0
      nlog  = 0

      step_hdf  = -1
      step_res  = -1

      nend = 1
      tend = -1.0
      wend = huge(1.0)

      pid = getpid()

      scstatus = hostnm(hostfull)
      ihost = index(hostfull,'.')
      if (ihost .eq. 0) ihost = index(hostfull,' ')
      host = hostfull(1:ihost-1)

      if (proc == 0) then

         open(1,file=par_file)
         ierrh = 0
         do while (ierrh == 0 .and. parfilelines<maxparfilelines)
            read(unit=1, fmt='(a)', iostat=ierrh) parfile(parfilelines+1)
            if (ierrh == 0) then
               parfilelines = parfilelines + 1
               i = len_trim(parfile(parfilelines))
               if (i < len(parfile(parfilelines))) then
                  parfile(parfilelines)(i+1:i+1) = achar(0)
               else
                  call warn("[dataio:init_dataio] problem.par contains very long lines. The copy in the logfile and HDF dumps can be truncated.")
               endif
            endif
         enddo
         close(1)
         if (parfilelines == maxparfilelines) call warn("[dataio:init_dataio] problem.par has too many lines. The copy in the logfile and HDF dumps can be truncated.")

         diff_nml(OUTPUT_CONTROL)
         diff_nml(RESTART_CONTROL)
         diff_nml(END_CONTROL)

!  namelist /END_CONTROL/ nend, tend, wend
         ibuff(1)  = nend

         rbuff(1)  = tend
         rbuff(2)  = wend

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
         rbuff(45) = fmin
         rbuff(46) = fmax

         lbuff(1)  = vizit

         cbuff(40) = domain

         do iv = 1, nvarsmx
            cbuff(40+iv) = vars(iv)
         enddo

         cbuff(90) = mag_center
         cbuff(91) = user_message_file(1:cbuff_len)
         cbuff(92) = system_message_file(1:cbuff_len)

         call MPI_Bcast(cbuff, cbuff_len*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
         call MPI_Bcast(lbuff,    buffer_dim, MPI_LOGICAL,          0, comm, ierr)
         call MPI_Bcast(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
         call MPI_Bcast(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      else

         call MPI_Bcast(cbuff, cbuff_len*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
         call MPI_Bcast(lbuff,    buffer_dim, MPI_LOGICAL,          0, comm, ierr)
         call MPI_Bcast(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
         call MPI_Bcast(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

!  namelist /END_CONTROL/ nend, tend, wend
         nend                = ibuff(1)

         tend                = rbuff(1)
         wend                = rbuff(2)
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
         fmin                = rbuff(45)
         fmax                = rbuff(46)

         vizit               = lbuff(1)

         domain              = trim(cbuff(40))
         do iv=1, nvarsmx
            vars(iv)         = trim(cbuff(40+iv))
         enddo

         mag_center          = trim(cbuff(90))

         user_message_file   = trim(cbuff(91))
         system_message_file = trim(cbuff(92))

      endif

      tn = time_left(wend)

      last_hdf_time = -dt_hdf

      call init_hdf5(vars,ix,iy,iz,dt_plt)

      if (proc == 0 .and. restart .eq. 'last') call find_last_restart(nrestart)
      call MPI_Barrier(comm,ierr)
      call MPI_Bcast(nrestart, 1, MPI_INTEGER, 0, comm, ierr)

      call init_version
      if (proc == 0) then
         call printinfo("###############     Source configuration     ###############", .false.)
         do i=1,nenv
            call printinfo(env(i), .false.)
         enddo
         write(log_file,'(6a,i3.3,a)') trim(cwd),'/',trim(problem_name),'_',trim(run_id),'_',nrestart,'.log'
         ! ToDo: if the simulation is restarted then save previous log_file (if exists) under a different, unique name
         write(system_command, '("mv ",a," ",a)') trim(tmp_log_file), trim(log_file)
         system_status = SYSTEM(system_command)
      endif
      call set_container_chdf(nstep); chdf%nres = nrestart

      if (nrestart /= 0) then
         call read_restart_hdf5(chdf)
         call get_container(nstep)
         nstep_start = nstep
         t_start     = t
         nres_start  = nrestart
         nhdf_start  = nhdf-1
         if (new_id .ne. '') run_id=new_id
      endif
      call MPI_Bcast(log_file, cwdlen, MPI_CHARACTER, 0, comm, ierr)
      call set_container_chdf(nstep)
      if (all([bnd_xl,bnd_xr,bnd_yl,bnd_yr,bnd_zl,bnd_zr] /= "user")) then
         call all_fluid_boundaries
#ifdef MAGNETIC
         call all_mag_boundaries
#endif /* MAGNETIC */
      endif

      dataio_initialized = .true.

   end subroutine init_dataio

   subroutine cleanup_dataio
      use dataio_hdf5,     only: cleanup_hdf5
      implicit none

      call cleanup_hdf5
   end subroutine cleanup_dataio

   subroutine user_msg_handler(end_sim)

      use dataio_hdf5,   only: write_hdf5, write_restart_hdf5
      use dataio_pub,    only: chdf, step_hdf, msg, printinfo, warn, set_container_chdf
      use mpisetup,      only: MPI_CHARACTER, MPI_DOUBLE_PRECISION, comm, ierr, proc, nstep

      implicit none

      logical, intent(inout) :: end_sim
      integer :: tsleep

!--- process 0 checks for messages

      if (proc == 0) call read_file_msg

      call MPI_Bcast(umsg,       umsg_len, MPI_CHARACTER,        0, comm, ierr)
      call MPI_Bcast(umsg_param, 1,        MPI_DOUBLE_PRECISION, 0, comm, ierr)

!---  if a user message is received then:
      if (len_trim(umsg) /= 0) then
         select case (trim(umsg))
            case ('res','dump')
               call write_restart_hdf5
            case ('hdf')
               call set_container_chdf(nstep)
               call write_hdf5(chdf)
               step_hdf = nstep
            case ('log')
               call write_log
            case ('tsl')
               call write_timeslice
            case ('tend')
               tend   = umsg_param
            case ('nend')
               nend   = umsg_param
            case ('dtres')
               dt_res = umsg_param
            case ('dthdf')
               dt_hdf = umsg_param
            case ('dtlog')
               dt_log = umsg_param
            case ('dttsl')
               dt_tsl = umsg_param
            case ('dtplt')
               dt_plt = umsg_param
            case ('sleep')
               tsleep = 60*umsg_param
               call sleep(tsleep)
            case ('stop')
               end_sim = .true.
            case ('help')
               if (proc == 0) then
                  write(msg,*) "[dataio:user_msg_handler] Recognized messages:",char(10),&
                  &"  help     - prints this information",char(10),&
                  &"  stop     - finish the simulation",char(10),&
                  &"  res|dump - immediately dumps a restart file",char(10),&
                  &"  hdf      - dumps a plotfile",char(10),&
                  &"  log      - update logfile",char(10),&
                  &"  tsl      - write a timeslice",char(10),&
                  &"  sleep <number> - wait <number> seconds",char(10),&
                  &"  tend|nend|dtres|dthdf|dtlog|dttsl|dtplt <value> - update specified parameter with <value>",char(10),&
                  &"Note that only one line at a time is read."
                  call printinfo(msg)
               endif
            case default
               if (proc == 0) then
                  write(msg,*) "[dataio:user_msg_handler]: non-recognized message '",trim(umsg),"'. Use message 'help' for list of valid keys."
                  call warn(msg)
               endif
         end select
      endif

   end subroutine user_msg_handler

!---------------------------------------------------------------------
!
! Makes data dump on abnormal Piernik termination
!
!---------------------------------------------------------------------
!
   subroutine write_crashed(msg)

      use dataio_pub,    only: nres, die
      use problem_pub,   only: problem_name

      implicit none

      character(len=*), intent(in) :: msg

      ! force output for diagnostics
      problem_name = "crash"
      dt_hdf = tiny(1.0)
      nres = 1
      call write_data(output='end')

      call die(msg)

   end subroutine write_crashed

!---------------------------------------------------------------------
!
! controls data dumping
!
!---------------------------------------------------------------------
!
   subroutine write_data(output)

      use dataio_hdf5,   only: write_hdf5, write_restart_hdf5, write_plot
      use dataio_pub,    only: chdf, nres, last_hdf_time, step_hdf, set_container_chdf
      use mpisetup,      only: t, nstep

      implicit none

      character(len=*), intent(in) :: output

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      if (output .eq. 'log' .or. output .eq. 'end') call write_log
      if (output .eq. 'log' .or. output .eq. 'end') call write_timeslice

!    CALL checkdf

#ifdef HDFSWEEP
      if (output .ne. 'gpt') then
         call write_log
         call write_timeslice
#else /* !HDFSWEEP */
      if (dt_hdf .gt. 0.0 .and. nstep .gt. step_hdf .and. output .ne. 'gpt') then
#endif /* !HDFSWEEP */

         if ((t-last_hdf_time) .ge. dt_hdf &
                .or. output .eq. 'hdf' .or. output .eq. 'end') then
            call set_container_chdf(nstep)
            call write_hdf5(chdf)

            if ((t-last_hdf_time) .ge. dt_hdf) last_hdf_time = last_hdf_time + dt_hdf
            if ((t-last_hdf_time) .ge. dt_hdf) last_hdf_time = t ! additional control
                          ! in the case of changing dt_hdf into smaller value via msg
            step_hdf = nstep
         endif
      endif

      if (dt_res .gt. 0.0 .and. nstep .gt. step_res) then
         if ((nres-nres_start) .lt. (int((t-t_start) / dt_res) + 1) &
                .or. output .eq. 'res' .or. output .eq. 'end') then
            if (nres > 0) then
               call write_restart_hdf5
            else
               nres = 1
            endif
            step_res = nstep
         endif
      endif
      call set_container_chdf(nstep)
      call write_plot

   end subroutine write_data

!------------------------------------------------------------------------

   subroutine find_last_restart(restart_number)

      use dataio_pub,    only: cwdlen, cwd
      use problem_pub,   only: problem_name, run_id
#if defined(__INTEL_COMPILER)
      use ifport,        only: unlink
#endif /* __INTEL_COMPILER */

      implicit none

#if defined(__PGI)
      include "lib3f.h"
#endif /* __PGI */

      integer, intent(out) :: restart_number

      character(len=cwdlen) :: file_name
      integer               :: nres, unlink_stat
      logical               :: exist
      character(len=cwdlen) :: file_name_base

      restart_number = 0

      write(file_name_base,'(a,a1,a3,a1)') trim(problem_name),'_',run_id,'_'

      unlink_stat = unlink('restart_list.tmp')

      do nres =999,0,-1
         write(file_name,'(a,a1,a,a1,a3,a1,i4.4,a4)') &
               trim(cwd),'/',trim(problem_name),'_', run_id,'_',nres,'.res'
         inquire(file = file_name, exist = exist)
         if (exist) then
            restart_number = nres
         return
         endif
      enddo

   end subroutine find_last_restart


   function mpi_sum4d_and_multiply(tab,factor) result(output)
      use mpisetup, only: MPI_DOUBLE_PRECISION, MPI_SUM, comm3d, ierr
      implicit none
      real, dimension(:,:,:,:), intent(in) :: tab
      real, intent(in)                     :: factor
      real :: local, output
      local = sum(tab(:,:,:,:)) * factor
      call MPI_ALLREDUCE(local, output, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm3d, ierr)
   end function mpi_sum4d_and_multiply

   function mpi_sum3d_and_multiply(tab,factor) result(output)
      use mpisetup, only: MPI_DOUBLE_PRECISION, MPI_SUM, comm3d, ierr
      implicit none
      real, dimension(:,:,:), intent(in) :: tab
      real, intent(in)                   :: factor
      real :: local, output
      local = sum(tab(:,:,:)) * factor
      call MPI_ALLREDUCE(local, output, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm3d, ierr)
   end function mpi_sum3d_and_multiply
!---------------------------------------------------------------------
!
! writes integrals to text file
!
!---------------------------------------------------------------------
!
   subroutine write_timeslice

      use arrays,          only: u, b, wa
      use dataio_pub,      only: cwdlen, cwd, user_tsl
      use diagnostics,     only: pop_vector
      use fluidindex,      only: nvar, iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz, ibx, iby, ibz
      use grid,            only: dvol, dx, dy, dz, is, ie, js, je, ks, ke, x, y, z, nxd, nyd, nzd
      use mpisetup,        only: proc, t, dt, smalld, nstep
      use problem_pub,     only: problem_name, run_id
      use types,           only: tsl_container, phys_prop
#ifndef ISO
      use fluidindex,      only: iarr_all_en
#endif /* !ISO */
#ifdef COSM_RAYS
      use fluidindex,      only: iarr_all_crs
#endif /* COSM_RAYS */
#ifdef GRAV
      use arrays,          only: gpot
#endif /* GRAV */

      implicit none

      character(len=cwdlen)                           :: tsl_file, head_fmt
      character(len=hnlen), dimension(:), allocatable :: tsl_names
      real, allocatable, dimension(:)                 :: tsl_vars
      type(phys_prop), pointer                        :: sn

      real, save :: tot_mass = 0.0, tot_momx = 0.0, tot_momy = 0.0, tot_momz = 0.0, &
                    tot_ener = 0.0, tot_eint = 0.0, tot_ekin = 0.0, tot_emag = 0.0, &
                    tot_epot = 0.0, tot_mflx = 0.0, tot_mfly = 0.0, tot_mflz = 0.0

      type(tsl_container) :: tsl

#ifdef COSM_RAYS
      real, save :: tot_encr = 0.0
#endif /* COSM_RAYS */
      real     :: cs_iso2
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
#ifdef IONIZED
      cs_iso2 = nvar%ion%cs2
#endif /* IONIZED */
#ifdef NEUTRAL
      cs_iso2 = nvar%neu%cs2
#else /* !NEUTRAL */
#ifndef IONIZED
      cs_iso2 = 0.0
#endif /* !IONIZED */
#endif /* !NEUTRAL */


      if (proc == 0) then
         write(tsl_file,'(a,a1,a,a1,a3,a1,i3.3,a4)') &
              trim(cwd),'/',trim(problem_name),'_', run_id,'_',nrestart,'.tsl'

         if (tsl_firstcall) then
            call pop_vector(tsl_names, hnlen, ["nstep   ", "time    ", "timestep"])
            call pop_vector(tsl_names, hnlen, ["mass", "momx", "momy", "momz", "ener", "epot", "eint", "ekin"])

#ifdef COSM_RAYS
            call pop_vector(tsl_names, hnlen, ["encr_tot", "encr_min", "encr_max"])
#endif /* COSM_RAYS */
#ifdef MAGNETIC
            call pop_vector(tsl_names, hnlen, ["emag   ", "mflx   ", "mfly   ", "mflz   ", "vai_max", "b_min  ", "b_max  "])
            call pop_vector(tsl_names, hnlen, ["divb_max"])
#ifdef RESISTIVE
            call pop_vector(tsl_names, hnlen, ["eta_max"])
#endif /* RESISTIVE */
#endif /* MAGNETIC */
#ifdef IONIZED
            call pop_vector(tsl_names, hnlen, ["vxi_max ", "vyi_max ", "vzi_max " , "csi_max ", "deni_min", "deni_max", "prei_min", "prei_max"])
#ifndef ISO
            call pop_vector(tsl_names, hnlen, ["temi_min", "temi_max"])
#endif /* !ISO */
#endif /* IONIZED */
#ifdef NEUTRAL
            call pop_vector(tsl_names, hnlen, ["denn_min", "denn_max", "vxn_max ", "vyn_max ", "vzn_max ", "pren_min", &
               "pren_max", "temn_min", "temn_max", "csn_max "])
#endif /* NEUTRAL */
#ifdef DUST
            call pop_vector(tsl_names, hnlen, ["dend_min", "dend_max", "vxd_max ", "vyd_max ", "vzd_max "])
#endif /* DUST */
            if (associated(user_tsl)) call user_tsl(tsl_vars, tsl_names)
            write(head_fmt,'(A,I2,A)') "(a1,a8,",size(tsl_names)-1,"a16)"

            open(tsl_lun, file=tsl_file)
            write(tsl_lun,fmt=head_fmt) "#",tsl_names
            write(tsl_lun, '(a1)') '#'
            deallocate(tsl_names)
            tsl_firstcall = .false.
         else
            open(tsl_lun, file=tsl_file, position='append')
         endif
      endif

      tot_mass = mpi_addmul(u(iarr_all_dn,is:ie,js:je,ks:ke), dvol)
      tot_momx = mpi_addmul(u(iarr_all_mx,is:ie,js:je,ks:ke), dvol)
      tot_momy = mpi_addmul(u(iarr_all_my,is:ie,js:je,ks:ke), dvol)
      tot_momz = mpi_addmul(u(iarr_all_mz,is:ie,js:je,ks:ke), dvol)
#ifdef GRAV
      tot_epot = mpi_addmul(u(iarr_all_dn(1),is:ie,js:je,ks:ke) *gpot(is:ie,js:je,ks:ke),dvol)
#endif /* GRAV */

      wa(is:ie,js:je,ks:ke) &
          = 0.5 * (u(iarr_all_mx(1),is:ie,js:je,ks:ke)**2   &
                 + u(iarr_all_my(1),is:ie,js:je,ks:ke)**2   &
                 + u(iarr_all_mz(1),is:ie,js:je,ks:ke)**2)/ &
                   max(u(iarr_all_dn(1),is:ie,js:je,ks:ke),smalld)
      tot_ekin = mpi_addmul(wa(is:ie,js:je,ks:ke), dvol)

      wa(is:ie,js:je,ks:ke) &
         = 0.5 * (b(ibx,is:ie,js:je,ks:ke)**2 + &
                  b(iby,is:ie,js:je,ks:ke)**2 + &
                  b(ibz,is:ie,js:je,ks:ke)**2)
      tot_emag = mpi_addmul(wa(is:ie,js:je,ks:ke), dvol)

      tot_mflx = mpi_addmul(b(ibx,is:ie,js:je,ks:ke), dy*dz/nxd)
      tot_mfly = mpi_addmul(b(iby,is:ie,js:je,ks:ke), dx*dz/nyd)
      tot_mflz = mpi_addmul(b(ibz,is:ie,js:je,ks:ke), dx*dy/nzd)
#ifdef ISO
      tot_eint = cs_iso2*tot_mass
      tot_ener = tot_eint+tot_ekin+tot_emag
#else /* !ISO */
      tot_ener = mpi_addmul(u(iarr_all_en,is:ie,js:je,ks:ke), dvol)
      tot_eint = tot_ener - tot_ekin - tot_emag
#endif /* !ISO */
#ifdef GRAV
      tot_ener = tot_ener + tot_epot
#endif /* GRAV */

#ifdef COSM_RAYS
      tot_encr = mpi_addmul(u(iarr_all_crs,is:ie,js:je,ks:ke), dvol)
#endif /* COSM_RAYS */

      call write_log(tsl)

      if (proc == 0) then
         call pop_vector(tsl_vars, [t, dt, tot_mass, tot_momx, tot_momy, tot_momz, tot_ener, tot_epot, tot_eint, tot_ekin])
#ifdef MAGNETIC
         call pop_vector(tsl_vars, [tot_emag, tot_mflx, tot_mfly, tot_mflz, tsl%vai_max, tsl%b_min, tsl%b_max, tsl%divb_max])
#ifdef RESISTIVE
         call pop_vector(tsl_vars, [tsl%etamax])
#endif /* RESISTIVE */
#endif /* MAGNETIC */
#ifdef COSM_RAYS
         call pop_vector(tsl_vars, [tot_encr, tsl%encr_min, tsl%encr_max])
#endif /* COSM_RAYS */
#ifdef IONIZED
         sn=>nvar%ion%snap
         call pop_vector(tsl_vars, [sn%velx_max%val, sn%vely_max%val, sn%velz_max%val, sn%cs_max%val, &
                                    sn%dens_min%val, sn%dens_max%val, sn%pres_min%val, sn%pres_max%val])
#ifndef ISO
         call pop_vector(tsl_vars, [sn%temp_min%val, sn%temp_max%val])
#endif /* !ISO */
#endif /* IONIZED */

#ifdef NEUTRAL
         sn=>nvar%neu%snap
         call pop_vector(tsl_vars, [sn%dens_min%val, sn%dens_max%val, sn%velx_max%val, sn%vely_max%val, sn%velz_max%val, &
                                    sn%pres_min%val, sn%pres_max%val, sn%temp_min%val, sn%temp_max%val, sn%cs_max%val])
#endif /* NEUTRAL */
#ifdef DUST
         sn=>nvar%dst%snap
         call pop_vector(tsl_vars, [sn%dens_min%val, sn%dens_max%val, sn%velx_max%val, sn%vely_max%val, sn%velz_max%val])
#endif /* DUST */

         if (associated(user_tsl)) call user_tsl(tsl_vars)
         write(tsl_lun, '(1x,i8,50(1x,es15.8))') nstep, tsl_vars

! some quantities computed in "write_log".One can add more, or change.
         close(tsl_lun)
         deallocate(tsl_vars)
      endif

   end subroutine write_timeslice

   subroutine get_extremum(tab,minmax,prop)
      use dataio_pub,    only: msg, warn
      use types,         only: value
      use grid,          only: nb
      use mpisetup,      only: mpifind
      implicit none
      real, dimension(:,:,:), intent(in) :: tab
      character(len=*), intent(in)       :: minmax
      type(value), intent(out)           :: prop

      if (minmax=='min') then
         prop%val = minval(tab)
         prop%loc = minloc(tab) + [nb,nb,nb]
      elseif (minmax=='max') then
         prop%val = maxval(tab)
         prop%loc = maxloc(tab) + [nb,nb,nb]
      else
         write(msg,*) "[dataio:get_extremum]: I don't know what to do with minmax = ", trim(minmax)
         call warn(msg)
      endif

      call mpifind(prop%val, trim(minmax), prop%loc, prop%proc)
      return
   end subroutine get_extremum

   subroutine common_shout(pr, fluid, pres_tn, temp_tn, cs_tn)
      use constants,       only: small
      use dataio_pub,      only: msg, printinfo
      use grid,            only: dxmn, dx, dy, dz
      use mpisetup,        only: cfl
      use types,           only: phys_prop
      implicit none
      type(phys_prop), intent(in)  :: pr
      character(len=*), intent(in) :: fluid
      logical, intent(in)          :: pres_tn, temp_tn, cs_tn

      character(len=cwdlen), parameter :: fmt772 = "(5x,a12,2a3,(1x,es15.9),16x,5(1x,i4))"
      character(len=cwdlen), parameter :: fmt778 = "(5x,a12,2a3,(1x,es15.9),2x,a3,(1x,es10.4),5(1x,i4))"
      real :: dxmn_safe

      if (dxmn >= sqrt(huge(1.0))) then
         dxmn_safe = sqrt(huge(1.0))
      else
         dxmn_safe = dxmn
      endif

      write(msg, fmt772) 'min(dens)   ',fluid,'  =', pr%dens_min%val,  pr%dens_min%proc,  pr%dens_min%loc
      call printinfo(msg, .false.)
      write(msg, fmt772) 'max(dens)   ',fluid,'  =', pr%dens_max%val,  pr%dens_max%proc,  pr%dens_max%loc
      call printinfo(msg, .false.)
      if (temp_tn) then
         write(msg, fmt772) 'min(temp)   ',fluid,'  =', pr%temp_min%val,  pr%temp_min%proc,  pr%temp_min%loc
         call printinfo(msg, .false.)
         write(msg, fmt772) 'max(temp)   ',fluid,'  =', pr%temp_max%val,  pr%temp_max%proc,  pr%temp_max%loc
         call printinfo(msg, .false.)
      endif
      if (pres_tn) then
         write(msg, fmt772) 'min(pres)   ',fluid,'  =', pr%pres_min%val,  pr%pres_min%proc,  pr%pres_min%loc
         call printinfo(msg, .false.)
         write(msg, fmt772) 'max(pres)   ',fluid,'  =', pr%pres_max%val,  pr%pres_max%proc,  pr%pres_max%loc
         call printinfo(msg, .false.)
      endif

      write(msg, fmt778) 'max(|vx|)   ',fluid,'  =', pr%velx_max%val, 'dt=',cfl*dx/(pr%velx_max%val+small),   pr%velx_max%proc, pr%velx_max%loc
      call printinfo(msg, .false.)
      write(msg, fmt778) 'max(|vy|)   ',fluid,'  =', pr%vely_max%val, 'dt=',cfl*dy/(pr%vely_max%val+small),   pr%vely_max%proc, pr%vely_max%loc
      call printinfo(msg, .false.)
      write(msg, fmt778) 'max(|vz|)   ',fluid,'  =', pr%velz_max%val, 'dt=',cfl*dz/(pr%velz_max%val+small),   pr%velz_max%proc, pr%velz_max%loc
      call printinfo(msg, .false.)
      if (cs_tn) then
         write(msg, fmt778) 'max(c_s )   ',fluid,'  =', pr%cs_max%val, 'dt=',cfl*dxmn_safe/(pr%cs_max%val+small), pr%cs_max%proc, pr%cs_max%loc
         call printinfo(msg, .false.)
      endif
   end subroutine common_shout

   subroutine get_common_vars(fl)
      use arrays,     only: u, b, wa
      use mpisetup,   only: smallp
      use grid,       only: is, ie, js, je, ks, ke
      use types,      only: phys_prop, component_fluid
      use constants,  only: mH, kboltz, gasRconst
      implicit none
      type(component_fluid), intent(inout), target :: fl
      type(phys_prop), pointer                     :: pr

      pr => fl%snap
      wa = u(fl%idn,:,:,:)
      call get_extremum(wa(is:ie,js:je,ks:ke), 'max', pr%dens_max)
      call get_extremum(wa(is:ie,js:je,ks:ke), 'min', pr%dens_min)

      wa = abs(u(fl%imx,:,:,:)/u(fl%idn,:,:,:))
      call get_extremum(wa(is:ie,js:je,ks:ke), 'max', pr%velx_max)

      wa = abs(u(fl%imy,:,:,:)/u(fl%idn,:,:,:))
      call get_extremum(wa(is:ie,js:je,ks:ke), 'max', pr%vely_max)

      wa = abs(u(fl%imz,:,:,:)/u(fl%idn,:,:,:))
      call get_extremum(wa(is:ie,js:je,ks:ke), 'max', pr%velz_max)

#ifdef ISO
      pr%pres_min%val  = fl%cs2*pr%dens_min%val
      pr%pres_min%loc  = pr%dens_min%loc
      pr%pres_min%proc = pr%dens_min%proc
      pr%pres_max%val  = fl%cs2*pr%dens_max%val
      pr%pres_max%loc  = pr%dens_max%loc
      pr%pres_max%proc = pr%dens_max%proc
      pr%cs_max%val    = fl%cs
      pr%cs_max%loc    = 0
      pr%cs_max%proc   = 0
      pr%temp_min%val  = mH / kboltz * gasRconst/fl%gam * fl%cs2
      pr%temp_min%loc  = 0
      pr%temp_min%proc = 0
      pr%temp_max%val  = mH / kboltz * gasRconst/fl%gam * fl%cs2
      pr%temp_max%loc  = 0
      pr%temp_max%proc = 0
#else /* !ISO */
      if (fl%tag /= "DST") then
         wa(:,:,:) = (u(fl%ien,:,:,:) &                ! eint
                   - 0.5*((u(fl%imx,:,:,:)**2 +u(fl%imy,:,:,:)**2 + u(fl%imz,:,:,:)**2)/u(fl%idn,:,:,:)))
         if (fl%tag == "ION") wa(:,:,:) = wa(:,:,:) - 0.5*(sum(b**2,dim=1))

         wa(:,:,:) = max((fl%gam_1)*wa(:,:,:),smallp)  ! pres

         call get_extremum(wa(is:ie,js:je,ks:ke), 'max', pr%pres_max)
         call get_extremum(wa(is:ie,js:je,ks:ke), 'min', pr%pres_min)

         wa(:,:,:) = fl%gam*wa(:,:,:)
         call get_extremum(wa(is:ie,js:je,ks:ke), 'max', pr%cs_max)

         wa(:,:,:) = mH/kboltz / fl%gam * gasRconst * wa(:,:,:) / u(fl%idn,:,:,:)
         call get_extremum(wa(is:ie,js:je,ks:ke), 'max', pr%temp_max)
         call get_extremum(wa(is:ie,js:je,ks:ke), 'min', pr%temp_min)
      endif
#endif /* !ISO */
      end subroutine get_common_vars
   !---------------------------------------------------------------------
   !
   ! writes timestep diagnostics to the logfile
   !
   ! Quite costly routine due to extensive array searches
   !
   !---------------------------------------------------------------------
   !
      subroutine  write_log(tsl)

         use arrays,             only: wa, u, b
         use constants,          only: small
         use dataio_pub,         only: msg, printinfo
         use fluidindex,         only: ibx, iby, ibz, nvar
         use grid,               only: dx, dy, dz, dxmn, is, ie, js, je, ks, ke, nx, ny, nz
         use mpisetup,           only: cfl, t, dt, proc, mpifind
         use types,              only: tsl_container, value
#ifdef COSM_RAYS
         use fluidindex,         only: iarr_all_crs
         use timestepcosmicrays, only: dt_crs
#endif /* COSM_RAYS */
#ifdef RESISTIVE
         use resistivity,        only: dt_resist, eta_max
#endif /* RESISTIVE */
#ifdef VARIABLE_GP
         use arrays,             only: gp
#endif /* VARIABLE_GP */

         implicit none

#if (defined(MAGNETIC) || defined(COSM_RAYS)|| defined (VARIABLE_GP))
         character(len=cwdlen), parameter :: fmt771 = "(5x,a18,(1x,es15.9),16x,5(1x,i4))"
#endif /* MAGNETIC || COSM_RAYS || defined VARIABLE_GP */
#if (defined(MAGNETIC) || defined(COSM_RAYS) || defined(RESISTIVE))
         character(len=cwdlen), parameter :: fmt777 = "(5x,a18,(1x,es15.9),2x,a3,(1x,es10.4),5(1x,i4))"
#endif /* MAGNETIC || COSM_RAYS || RESISTIVE */

         type(tsl_container), optional  :: tsl
         real :: dxmn_safe

#ifdef MAGNETIC
         type(value) :: b_min, b_max, divb_max, vai_max
#endif /* MAGNETIC */
#ifdef COSM_RAYS
         type(value) :: encr_min, encr_max
#endif /* COSM_RAYS */
#ifdef RESISTIVE
         type(value) :: etamax
#endif /* RESISTIVE */
#ifdef VARIABLE_GP
         type(value) :: gpxmax, gpymax, gpzmax
#endif /* VARIABLE_GP */

         if (dxmn >= sqrt(huge(1.0))) then
            dxmn_safe = sqrt(huge(1.0))
         else
            dxmn_safe = dxmn
         endif

   ! Timestep diagnostics
#ifdef NEUTRAL
         call get_common_vars(nvar%neu)
#endif /* NEUTRAL */

#ifdef IONIZED
         call get_common_vars(nvar%ion)

#ifdef MAGNETIC
         wa(:,:,:)  = sqrt(b(1,:,:,:)*b(1,:,:,:) + b(2,:,:,:)*b(2,:,:,:) + b(3,:,:,:)*b(3,:,:,:))
         call get_extremum(wa(is:ie,js:je,ks:ke), 'max', b_max)
         call get_extremum(wa(is:ie,js:je,ks:ke), 'min', b_min)

         wa(:,:,:)  = wa(:,:,:) / sqrt(u(nvar%ion%idn,:,:,:))
         call get_extremum(wa(is:ie,js:je,ks:ke), 'max', vai_max)
#endif /* MAGNETIC */

#ifdef ISO
#ifdef ISO_LOCAL
!        wa            = cs_iso2_arr(:,:,:)*u(idni,:,:,:)
!        prei_min%val  = minval(wa(is:ie,js:je,ks:ke))
!        prei_min%loc  = minloc(wa(is:ie,js:je,ks:ke)) + [nb, nb, nb]
!        call mpifind(prei_min%val, 'min', prei_min%loc, prei_min%proc)
!
!        prei_max%val  = maxval(wa(is:ie,js:je,ks:ke))
!        prei_max%loc  = maxloc(wa(is:ie,js:je,ks:ke)) + [nb, nb, nb]
!        call mpifind(prei_max%val, 'max', prei_max%loc, prei_max%proc)
!
!        csi_max%val   = maxval(cs_iso2_arr(is:ie,js:je,ks:ke))
!        csi_max%loc   = maxloc(cs_iso2_arr(is:ie,js:je,ks:ke)) + [nb, nb, nb]
!        call mpifind(csi_max%val, 'max', csi_max%loc, csi_max%proc)
!
!        wa            = mH / kboltz * cs_iso2_arr(:,:,:)
!        temi_min%val  = minval(wa(is:ie,js:je,ks:ke))
!        temi_min%loc  = minloc(wa(is:ie,js:je,ks:ke)) + [nb, nb, nb]
!        call mpifind(temi_min%val, 'min', temi_min%loc, temi_min%proc)
!
!        temi_max%val  = maxval(wa(is:ie,js:je,ks:ke))
!        temi_max%loc  = maxloc(wa(is:ie,js:je,ks:ke)) + [nb, nb, nb]
!        call mpifind(temi_max%val, 'max', temi_max%loc, temi_max%proc)
#endif /* ISO_LOCAL */
#else /* !ISO */
!        wa(:,:,:) = (u(ieni,:,:,:) &                ! eint
!                  - 0.5*((u(imxi,:,:,:)**2 +u(imyi,:,:,:)**2 &
!                    + u(imzi,:,:,:)**2)/u(idni,:,:,:)))
#ifdef MAGNETIC
!        wa(:,:,:) = wa(:,:,:) - 0.5*(b(ibx,:,:,:)**2 + b(iby,:,:,:)**2 + &
!                   b(ibz,:,:,:)**2)
#endif /* MAGNETIC */
#endif /* !ISO */
#endif /* IONIZED */

#ifdef DUST
      call get_common_vars(nvar%dst)
#endif /* DUST */

#ifdef VARIABLE_GP
      wa(1:nx-1,:,:) = abs((gp(2:nx,:,:)-gp(1:nx-1,:,:))/dx)
      call get_extremum(wa(is:ie,js:je,ks:ke), 'max', gpxmax)
      wa(:,1:ny-1,:) = abs((gp(:,2:ny,:)-gp(:,1:ny-1,:))/dy)
      call get_extremum(wa(is:ie,js:je,ks:ke), 'max', gpymax)
      wa(:,:,1:nz-1) = abs((gp(:,:,2:nz)-gp(:,:,1:nz-1))/dz)
      call get_extremum(wa(is:ie,js:je,ks:ke), 'max', gpzmax)
#endif /* VARIABLE_GP */

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

      call get_extremum(wa(is:ie,js:je,ks:ke), 'max', divb_max)
#endif /* MAGNETIC */

#ifdef COSM_RAYS
      wa            = sum(u(iarr_all_crs,:,:,:),1)
      call get_extremum(wa(is:ie,js:je,ks:ke), 'max', encr_max)
      call get_extremum(wa(is:ie,js:je,ks:ke), 'min', encr_min)
#endif /* COSM_RAYS */

      if (proc == 0)  then
         if (.not.present(tsl)) then
            call printinfo('================================================================================', .false.)
#ifdef IONIZED
            call common_shout(nvar%ion%snap,'ION',.true.,.true.,.true.)
#ifdef MAGNETIC
            write(msg, fmt777) 'max(c_f)    ION  =', sqrt(nvar%ion%snap%cs_max%val**2+vai_max%val**2), &
               'dt=',cfl*dxmn_safe/sqrt(nvar%ion%snap%cs_max%val**2+vai_max%val**2+small)
            call printinfo(msg, .false.)
            write(msg, fmt777) 'max(v_a)    ION  =', vai_max%val, 'dt=',cfl*dxmn_safe/(vai_max%val+small), vai_max%proc, vai_max%loc
            call printinfo(msg, .false.)
            write(msg, fmt771) 'min(|b|)    MAG  =', b_min%val,     b_min%proc,     b_min%loc
            call printinfo(msg, .false.)
            write(msg, fmt771) 'max(|b|)    MAG  =', b_max%val,     b_max%proc,     b_max%loc
            call printinfo(msg, .false.)
            write(msg, fmt771) 'max(|divb|) MAG  =', divb_max%val,  divb_max%proc,  divb_max%loc
            call printinfo(msg, .false.)
#else /* !MAGNETIC */
!            if (csi_max%val > 0.) write(msg, fmt777) 'max(c_s )   ION  =', sqrt(csi_max%val**2), 'dt=',cfl*dxmn_safe/sqrt(csi_max%val**2)
!            call printinfo(msg, .false.)
#endif /* !MAGNETIC */
#endif /* IONIZED */
#ifdef NEUTRAL
            call common_shout(nvar%neu%snap,'NEU',.true.,.true.,.true.)
#endif /* NEUTRAL */
#ifdef DUST
            call common_shout(nvar%dst%snap,'DST',.false.,.false.,.false.)
#endif /* DUST */
#ifdef COSM_RAYS
            write(msg, fmt771) 'min(encr)   CRS  =', encr_min%val,        encr_min%proc, encr_min%loc
            call printinfo(msg, .false.)
            write(msg, fmt777) 'max(encr)   CRS  =', encr_max%val,      'dt=',dt_crs,     encr_max%proc, encr_max%loc
            call printinfo(msg, .false.)
#endif /* COSM_RAYS */
#ifdef RESISTIVE
            write(msg, fmt777) 'max(eta)    RES  =', etamax%val ,      'dt=',dt_resist, etamax%proc,  etamax%loc
            call printinfo(msg, .false.)
#endif /* RESISTIVE */
#ifdef VARIABLE_GP
            write(msg,fmt771) 'max(|gpx|)  GPT  =', gpxmax%val,   gpxmax%proc,     gpxmax%loc
            call printinfo(msg, .false.)
            write(msg,fmt771) 'max(|gpy|)  GPT  =', gpymax%val,   gpymax%proc,     gpymax%loc
            call printinfo(msg, .false.)
            write(msg,fmt771) 'max(|gpz|)  GPT  =', gpzmax%val,   gpzmax%proc,     gpzmax%loc
            call printinfo(msg, .false.)
#endif /* VARIABLE_GP */
            call printinfo('================================================================================', .false.)
         else
#ifdef MAGNETIC
            tsl%b_min = b_min%val
            tsl%b_max = b_max%val
            tsl%divb_max = divb_max%val
            tsl%vai_max  = vai_max%val
#endif /* MAGNETIC */

#ifdef COSM_RAYS
            tsl%encr_min = encr_min%val
            tsl%encr_max = encr_max%val
#endif /* COSM_RAYS */

#ifdef RESISTIVE
            tsl%etamax = etamax%val
#endif /* RESISTIVE */

#ifdef VARIABLE_GP
            tsl%gpxmax = gpxmax%val
            tsl%gpymax = gpymax%val
            tsl%gpzmax = gpzmax%val
#endif /* VARIABLE_GP */
         endif
      endif

   end subroutine write_log
!------------------------------------------------------------------------
   subroutine read_file_msg
!-------------------------------------------------------------------------
!     configurable parameters: problem.par
!-------------------------------------------------------------------------
!      user_message_file           ! 1st (user) message file (eg.'./msg')
!      system_message_file         ! 2nd (ups)  message file (eg.'/etc/ups/user/msg')
!-------------------------------------------------------------------------

!\todo: process multiple commands at once
      use dataio_pub,    only: cwdlen, msg, printinfo, warn
      use mpisetup,      only: proc
#if defined(__INTEL_COMPILER)
      use ifport,        only: unlink, stat
#endif /* __INTEL_COMPILER */

      implicit none

#if defined(__PGI)
      include "lib3f.h"
#endif /* __PGI */

      integer, parameter :: n_msg_origin = 2, msg_lun=91
      character(len=*), parameter, dimension(n_msg_origin) :: msg_origin = [ "user  ", "system" ]

      character(len=cwdlen), dimension(n_msg_origin), save :: fname
      integer ::  unlink_stat, io, sz, sts, i
      integer, dimension(13) :: stat_buff
      logical :: msg_param_read = .false., ex
      integer, dimension(n_msg_origin), save :: last_msg_stamp

      umsg=''
      umsg_param = 0.0
      sz = -1
      fname = [ user_message_file, system_message_file ]

      do i = 1, n_msg_origin
         inquire(FILE=fname(i), EXIST=ex)
         if (ex .and. sz<=0) then ! process only one message file at a time, user_message_file first
            ! I think system file is more important, but current logic prevents reading user_message_file when system_message_file is present

            sts = stat(fname(i), stat_buff)
            if (last_msg_stamp(i) == stat_buff(10)) exit
            last_msg_stamp(i) = stat_buff(10)

            open(msg_lun, file=fname(i), status='old')
            read(msg_lun, *, iostat=io) umsg, umsg_param
            if (io/=0) then
               rewind(msg_lun)
               read(msg_lun, *, iostat=io) umsg
               if (io/=0) then
                  write(msg, '(5a)'   )"[dataio:read_file_msg] ",trim(msg_origin(i))," message: '",trim(umsg),"'"
               else
                  write(msg, '(3a)'   )"[dataio:read_file_msg] Cannot read ",trim(msg_origin(i))," message."
                  call warn(msg)
                  msg=""
               endif
            else
               msg_param_read = .true.
               write(msg, '(5a,g15.7)')"[dataio:read_file_msg] ",trim(msg_origin(i))," message: '",trim(umsg),"', with parameter = ", umsg_param
            endif
            close(msg_lun)

            if (len_trim(msg) > 0 .and. proc==0) call printinfo(msg)

            sz = len_trim(msg)
            if (fname(i) == user_message_file) unlink_stat = unlink(user_message_file)

         endif
      enddo

   end subroutine read_file_msg
!------------------------------------------------------------------------
end module dataio
