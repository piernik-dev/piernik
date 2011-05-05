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
!! \brief (DW) Module containing all main routines  responsible for data output
!!
!!
!! In this module following namelists of parameters are specified:
!! \copydetails dataio::init_dataio
!!
!! \todo check and if necessary bring back usefulness of min_disk_space_MB parameter
!<

module dataio

   use dataio_pub,    only: domain, fmin, fmax, vizit, nend, tend, wend, nrestart, problem_name, run_id
   use constants,     only: cwdlen, cbuff_len, varlen, idlen

   implicit none

   private
   public :: check_log, check_tsl, write_data, write_crashed, cleanup_dataio, init_dataio, user_msg_handler

   integer               :: istep                  !< current number of substep (related to integration order)

   integer, parameter    :: nvarsmx = 16           !< maximum number of variables to dump in hdf files
   character(len=idlen)  :: new_id                 !< three character string to change run_id when restarting simulation (e.g. to avoid overwriting of the output from the previous (pre-restart) simulation; if new_id = '' then run_id is still used)
   character(len=cbuff_len) :: restart             !< choice of restart file: if restart = 'last': automatic choice of the last restart file regardless of "nrestart" value; if something else is set: "nrestart" value is fixing
   character(len=cbuff_len) :: mag_center          !< choice to dump magnetic fields values from cell centers or not (if not then values from cell borders, unused)
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

   integer               :: nchar                  !< number of characters in a user/system message
   integer, parameter    :: umsg_len = 16
   character(len=umsg_len) :: umsg                    !< string of characters - content of a user/system message
   real                  :: umsg_param              !< parameter changed by a user/system message

   character(len=cwdlen) :: filename               !< string of characters indicating currently used file
   character(len=cwdlen), parameter :: fmt_loc   = "(2x,a12,a3,'  = ',es16.9,16x,            5(1x,i4))"
   character(len=cwdlen), parameter :: fmt_dtloc = "(2x,a12,a3,'  = ',es16.9,'  dt=',es11.4, 5(1x,i4))"

   namelist /END_CONTROL/ nend, tend, wend
   namelist /RESTART_CONTROL/ restart, new_id, nrestart, resdel
   namelist /OUTPUT_CONTROL/ problem_name, run_id, dt_hdf, dt_res, dt_tsl, dt_log, dt_plt, ix, iy, iz, &
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
      use dataio_pub,    only: next_t_log

      implicit none

      if (dt_log > 0.0) then
         if (next_t_log <= t) then
            call write_log
            next_t_log = next_t_log + dt_log
         endif
      endif

   end subroutine check_log

   subroutine check_tsl

      use mpisetup,      only: t
      use dataio_pub,    only: next_t_tsl

      implicit none

      if (dt_tsl > 0.0) then
         if (next_t_tsl <= t) then
            call write_timeslice
            next_t_tsl = next_t_tsl + dt_tsl
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
!! <tr><td>nend</td><td>1        </td><td>integer</td><td>\copydoc dataio_pub::nend</td></tr>
!! <tr><td>tend</td><td>-1.0     </td><td>real   </td><td>\copydoc dataio_pub::tend</td></tr>
!! <tr><td>wend</td><td>huge(1.0)</td><td>real   </td><td>\copydoc dataio_pub::wend</td></tr>
!! </table>
!! \n \n
!! @b RESTART_CONTROL
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>restart </td><td>'last'</td><td>'last' or another string of characters</td><td>\copydoc dataio::restart </td></tr>
!! <tr><td>new_id  </td><td>''    </td><td>string of characters                  </td><td>\copydoc dataio::new_id  </td></tr>
!! <tr><td>nrestart</td><td>3     </td><td>integer                               </td><td>\copydoc dataio_pub::nrestart</td></tr>
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
!! <tr><td>domain             </td><td>'phys_domain'     </td><td>'phys_domain' or 'full_domain'                     </td><td>\copydoc dataio_pub::domain</td></tr>
!! <tr><td>vars               </td><td>''                 </td><td>'dens', 'velx', 'vely', 'velz', 'ener' and some more </td><td>\copydoc dataio::vars  </td></tr>
!! <tr><td>mag_center         </td><td>'no'               </td><td>'yes'/'no'</td><td>\copydoc dataio::mag_center       </td></tr>
!! <tr><td>vizit              </td><td>.false.            </td><td>logical   </td><td>\copydoc dataio_pub::vizit        </td></tr>
!! <tr><td>fmin               </td><td>                   </td><td>real      </td><td>\copydoc dataio_pub::fmin         </td></tr>
!! <tr><td>fmax               </td><td>                   </td><td>real      </td><td>\copydoc dataio_pub::fmax         </td></tr>
!! <tr><td>min_disk_space_MB  </td><td>100                </td><td>integer   </td><td>\copydoc dataio::min_disk_space_MB</td></tr>
!! <tr><td>sleep_minutes      </td><td>0                  </td><td>integer   </td><td>\copydoc dataio::sleep_minutes    </td></tr>
!! <tr><td>sleep_seconds      </td><td>0                  </td><td>integer   </td><td>\copydoc dataio::sleep_seconds    </td></tr>
!! <tr><td>user_message_file  </td><td>trim(cwd)//'/msg'  </td><td>string similar to default value              </td><td>\copydoc dataio::user_message_file  </td></tr>
!! <tr><td>system_message_file</td><td>'/tmp/piernik_msg'</td><td>string of characters similar to default value</td><td>\copydoc dataio::system_message_file</td></tr>
!! </table>
!! \n \n
!<
   subroutine init_dataio

      use constants,       only: small, cwdlen, cbuff_len, PIERNIK_INIT_IO_IC, BND_USER
      use dataio_hdf5,     only: init_hdf5, read_restart_hdf5, parfile, parfilelines
      use dataio_pub,      only: chdf, nres, last_hdf_time, step_hdf, next_t_log, next_t_tsl, log_file_initialized, log_file, maxparfilelines, cwd, &
           &                     tmp_log_file, msglen, printinfo, warn, msg, nhdf, nstep_start, set_container_chdf, get_container, die, code_progress
      use dataio_pub,      only: par_file, ierrh, namelist_errh, compare_namelist, cmdl_nml  ! QA_WARN required for diff_nml
      use fluidboundaries, only: all_fluid_boundaries
      use mpisetup,        only: lbuff, ibuff, rbuff, cbuff, master, slave, comm, ierr, buffer_dim, t, nstep
      use mpi,             only: MPI_CHARACTER, MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_LOGICAL
      use timer,           only: time_left
      use version,         only: nenv,env, init_version
      use grid,            only: cg
#ifdef MAGNETIC
      use magboundaries,   only: all_mag_boundaries
#endif /* MAGNETIC */

      implicit none

      logical              :: tn
      integer(kind=1)      :: system
      integer              :: system_status, i
      character(len=msglen):: system_command

      if (code_progress < PIERNIK_INIT_IO_IC) call die("[dataio:init_dataio] Some physics modules are not initialized.")

      problem_name = "nameless"
      run_id = "___"
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

      tsl_firstcall = .true.

      nhdf  = 0
      nres  = 0
      next_t_tsl  = -1.*small
      next_t_log  = -1.*small

      step_hdf  = -1
      step_res  = -1

      nend = 1
      tend = -1.0
      wend = huge(1.0)

      if (master) then

         open(1,file=par_file)
         ierrh = 0
         do while (ierrh == 0 .and. parfilelines<maxparfilelines)
            read(unit=1, fmt='(a)', iostat=ierrh) parfile(parfilelines+1)
            if (ierrh == 0) then
               parfilelines = parfilelines + 1
               i = len_trim(parfile(parfilelines))
               if (i >= len(parfile(parfilelines))) call warn("[dataio:init_dataio] problem.par contains very long lines. The copy in the logfile and HDF dumps can be truncated.")
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

         cbuff(31) = problem_name
         cbuff(32) = run_id
         cbuff(40) = domain

         do iv = 1, nvarsmx
            cbuff(40+iv) = vars(iv)
         enddo

         cbuff(90) = mag_center
         cbuff(91) = user_message_file(1:cbuff_len)
         cbuff(92) = system_message_file(1:cbuff_len)

      endif

      call MPI_Bcast(cbuff, cbuff_len*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_Bcast(lbuff,           buffer_dim, MPI_LOGICAL,          0, comm, ierr)
      call MPI_Bcast(ibuff,           buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_Bcast(rbuff,           buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      if (slave) then

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

         problem_name        = cbuff(31)
         run_id              = cbuff(32)(1:idlen)
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

      if (master .and. restart == 'last') call find_last_restart(nrestart)
      call MPI_Barrier(comm,ierr)
      call MPI_Bcast(nrestart, 1, MPI_INTEGER, 0, comm, ierr)

      call init_version
      if (master) then
         call printinfo("###############     Source configuration     ###############", .false.)
         do i=1,nenv
            call printinfo(env(i), .false.)
         enddo
         write(log_file,'(6a,i3.3,a)') trim(cwd),'/',trim(problem_name),'_',trim(run_id),'_',nrestart,'.log'
!> \todo if the simulation is restarted then save previous log_file (if exists) under a different, unique name
         write(system_command, '("mv ",a," ",a)') trim(tmp_log_file), trim(log_file)
         system_status = system(system_command)
         if (system_status /= 0) then
            write(msg,'(2a)')"[dataio:init_dataio] The log must be stored in ",tmp_log_file
            call warn(msg)
            log_file_initialized = .false.
         else
            log_file_initialized = .true.
         endif
      endif
      call MPI_Bcast(log_file, cwdlen, MPI_CHARACTER, 0, comm, ierr)          ! BEWARE: every msg issued by slaves before this sync may lead to race condition on tmp_log_file
      call MPI_Bcast(log_file_initialized, 1, MPI_LOGICAL, 0, comm, ierr)

      call set_container_chdf(nstep); chdf%nres = nrestart

      if (nrestart /= 0) then
         if (master) call printinfo("###############     Reading restart     ###############", .false.)
         call read_restart_hdf5(chdf)
         call get_container(nstep)
         nstep_start = nstep
         t_start     = t
         nres_start  = nrestart
         nhdf_start  = nhdf-1
         if (new_id /= '') run_id=new_id
         if (all(cg%bnd(:,:) /= BND_USER)) then
            call all_fluid_boundaries
#ifdef MAGNETIC
            call all_mag_boundaries
#endif /* MAGNETIC */
         endif
      endif
      call set_container_chdf(nstep)

   end subroutine init_dataio

   subroutine cleanup_dataio
      use dataio_hdf5,     only: cleanup_hdf5
      implicit none

      call cleanup_hdf5
   end subroutine cleanup_dataio

   subroutine user_msg_handler(end_sim)

      use dataio_hdf5,   only: write_hdf5, write_restart_hdf5
      use dataio_pub,    only: chdf, step_hdf, msg, printinfo, warn, set_container_chdf
      use mpisetup,      only: comm, ierr, master, nstep
      use mpi,           only: MPI_CHARACTER, MPI_DOUBLE_PRECISION

      implicit none

      logical, intent(inout) :: end_sim
      integer :: tsleep

!--- process 0 checks for messages

      if (master) call read_file_msg

      call MPI_Bcast(umsg,       umsg_len, MPI_CHARACTER,        0, comm, ierr)
      call MPI_Bcast(umsg_param, 1,        MPI_DOUBLE_PRECISION, 0, comm, ierr)

!---  if a user message is received then:
      if (len_trim(umsg) /= 0) then
         select case (trim(umsg))
            case ('res')
               call write_restart_hdf5
            case ('dump')
               call write_restart_hdf5(.true.)
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
               nend   = int(umsg_param)
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
               tsleep = int(60*umsg_param)
               call sleep(tsleep)
            case ('stop')
               end_sim = .true.
            case ('help')
               if (master) then
                  write(msg,*) "[dataio:user_msg_handler] Recognized messages:",char(10),&
                  &"  help     - prints this information",char(10),&
                  &"  stop     - finish the simulation",char(10),&
                  &"  res      - immediately dumps a restart file",char(10),&
                  &"  dump     - immediately dumps a restart file of full domain for all blocks",char(10),&
                  &"  hdf      - dumps a plotfile",char(10),&
                  &"  log      - update logfile",char(10),&
                  &"  tsl      - write a timeslice",char(10),&
                  &"  sleep <number> - wait <number> seconds",char(10),&
                  &"  tend|nend|dtres|dthdf|dtlog|dttsl|dtplt <value> - update specified parameter with <value>",char(10),&
                  &"Note that only one line at a time is read."
                  call printinfo(msg)
               endif
            case default
               if (master) then
                  write(msg,*) "[dataio:user_msg_handler]: non-recognized message '",trim(umsg),"'. Use message 'help' for list of valid keys."
                  call warn(msg)
               endif
         end select
      endif

   end subroutine user_msg_handler

!---------------------------------------------------------------------
!>
!! |brief Makes data dump on abnormal Piernik termination
!<
!---------------------------------------------------------------------
!
   subroutine write_crashed(msg)

      use dataio_pub,    only: nres, die

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
!>
!! \brief controls data dumping
!<
!---------------------------------------------------------------------
!
   subroutine write_data(output)

      use dataio_hdf5,   only: write_hdf5, write_restart_hdf5, write_plot
      use dataio_pub,    only: chdf, nres, last_hdf_time, step_hdf, set_container_chdf
      use mpisetup,      only: t, nstep

      implicit none

      character(len=*), intent(in) :: output

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      if (output == 'log' .or. output == 'end') call write_log
      if (output == 'tsl' .or. output == 'end') call write_timeslice

!    call checkdf

      if (dt_hdf > 0.0 .and. nstep > step_hdf .and. output /= 'gpt') then
         if ((t-last_hdf_time) >= dt_hdf .or. output == 'hdf' .or. output == 'end') then
            call set_container_chdf(nstep)
            call write_hdf5(chdf)

            if ((t-last_hdf_time) >= dt_hdf) last_hdf_time = last_hdf_time + dt_hdf
            if ((t-last_hdf_time) >= dt_hdf) last_hdf_time = t ! additional control
                          ! in the case of changing dt_hdf into smaller value via msg
            step_hdf = nstep
         endif
      endif

      if (dt_res > 0.0 .and. nstep > step_res) then
         if ((nres-nres_start) < (int((t-t_start) / dt_res) + 1) &
                .or. output == 'res' .or. output == 'end') then
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

      use dataio_pub,    only: cwd
      use constants,     only: cwdlen
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
      use mpisetup, only: comm, ierr
      use mpi,      only: MPI_DOUBLE_PRECISION, MPI_SUM
      implicit none
      real, dimension(:,:,:,:), intent(in) :: tab
      real, intent(in)                     :: factor
      real :: local, output
      local = sum(tab(:,:,:,:)) * factor
      call MPI_Allreduce(local, output, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
   end function mpi_sum4d_and_multiply

   function mpi_sum3d_and_multiply(tab,factor) result(output)
      use mpisetup, only: comm, ierr
      use mpi,      only: MPI_DOUBLE_PRECISION, MPI_SUM
      implicit none
      real, dimension(:,:,:), intent(in) :: tab
      real, intent(in)                   :: factor
      real :: local, output
      local = sum(tab(:,:,:)) * factor
      call MPI_Allreduce(local, output, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
   end function mpi_sum3d_and_multiply
!---------------------------------------------------------------------
!>
!! \brief writes integrals to text file
!<
!---------------------------------------------------------------------
!
   subroutine write_timeslice

      use arrays,          only: u, b, wa
      use dataio_pub,      only: cwd, user_tsl
      use fluids_pub,      only: has_ion, has_dst, has_neu
      use constants,       only: cwdlen, xdim, ydim, zdim
      use diagnostics,     only: pop_vector
      use fluidindex,      only: flind, iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz, ibx, iby, ibz
      use grid,            only: cg
      use mpisetup,        only: master, t, dt, smalld, nstep, dom
      use types,           only: tsl_container
      use fluidtypes,      only: phys_prop
#ifndef ISO
      use fluidindex,      only: iarr_all_en
#endif /* !ISO */
#ifdef COSM_RAYS
      use fluidindex,      only: iarr_all_crs
#endif /* COSM_RAYS */
#ifdef GRAV
      use arrays,          only: gpot
#endif /* GRAV */
#ifdef RESISTIVE
      use resistivity,     only: eta1_active
#endif /* RESISTIVE */

      implicit none

      character(len=cwdlen)                           :: tsl_file, head_fmt
      character(len=cbuff_len), dimension(:), allocatable :: tsl_names
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
      if (has_ion) then
         cs_iso2 = flind%ion%cs2
      elseif (has_neu) then
         cs_iso2 = flind%neu%cs2
      else
         cs_iso2 = 0.0
      endif

      if (master) then
         write(tsl_file,'(a,a1,a,a1,a3,a1,i3.3,a4)') &
              trim(cwd),'/',trim(problem_name),'_', run_id,'_',nrestart,'.tsl'

         if (tsl_firstcall) then
            call pop_vector(tsl_names, cbuff_len, ["nstep   ", "time    ", "timestep"])
            call pop_vector(tsl_names, cbuff_len, ["mass", "momx", "momy", "momz", "ener", "epot", "eint", "ekin"])

#ifdef MAGNETIC
            call pop_vector(tsl_names, cbuff_len, ["emag   ", "mflx   ", "mfly   ", "mflz   ", "vai_max", "b_min  ", "b_max  "])
            call pop_vector(tsl_names, cbuff_len, ["divb_max"])
#ifdef RESISTIVE
            if (eta1_active) call pop_vector(tsl_names, cbuff_len, ["eta_max"])
#endif /* RESISTIVE */
#endif /* MAGNETIC */
#ifdef COSM_RAYS
            call pop_vector(tsl_names, cbuff_len, ["encr_tot", "encr_min", "encr_max"])
#endif /* COSM_RAYS */
            ! \todo: replicated code, simplify me
            if (has_ion) call pop_vector(tsl_names, cbuff_len, ["deni_min", "deni_max", "vxi_max ", "vyi_max ", "vzi_max ", &
                                                   "prei_min", "prei_max", "temi_min", "temi_max", "csi_max "])
            if (has_neu) call pop_vector(tsl_names, cbuff_len, ["denn_min", "denn_max", "vxn_max ", "vyn_max ", "vzn_max ", &
                                                   "pren_min", "pren_max", "temn_min", "temn_max", "csn_max "])
            if (has_dst) call pop_vector(tsl_names, cbuff_len, ["dend_min", "dend_max", "vxd_max ", "vyd_max ", "vzd_max "])

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

      tot_mass = mpi_addmul(u(iarr_all_dn, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), cg%dvol)
      tot_momx = mpi_addmul(u(iarr_all_mx, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), cg%dvol)
      tot_momy = mpi_addmul(u(iarr_all_my, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), cg%dvol)
      tot_momz = mpi_addmul(u(iarr_all_mz, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), cg%dvol)
#ifdef GRAV
      tot_epot = mpi_addmul(u(iarr_all_dn(1), cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) *gpot(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), cg%dvol)
#endif /* GRAV */

      wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) = &
           & 0.5 * (u(iarr_all_mx(1), cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke)**2   &
           &      + u(iarr_all_my(1), cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke)**2   &
           &      + u(iarr_all_mz(1), cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke)**2)/ &
           & max(u(iarr_all_dn(1), cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke),smalld)
      tot_ekin = mpi_addmul(wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), cg%dvol)

      wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) = &
           & 0.5 * (b(ibx, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke)**2 + &
           &        b(iby, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke)**2 + &
           &        b(ibz, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke)**2)
      tot_emag = mpi_addmul(wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), cg%dvol)

      tot_mflx = mpi_addmul(b(ibx, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), cg%dy*cg%dz/dom%n_d(xdim))
      tot_mfly = mpi_addmul(b(iby, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), cg%dx*cg%dz/dom%n_d(ydim))
      tot_mflz = mpi_addmul(b(ibz, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), cg%dx*cg%dy/dom%n_d(zdim))
#ifdef ISO
      tot_eint = cs_iso2*tot_mass
      tot_ener = tot_eint+tot_ekin+tot_emag
#else /* !ISO */
      tot_ener = mpi_addmul(u(iarr_all_en, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), cg%dvol)
      tot_eint = tot_ener - tot_ekin - tot_emag
#endif /* !ISO */
#ifdef GRAV
      tot_ener = tot_ener + tot_epot
#endif /* GRAV */

#ifdef COSM_RAYS
      tot_encr = mpi_addmul(u(iarr_all_crs, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), cg%dvol)
#endif /* COSM_RAYS */

      call write_log(tsl)

      if (master) then
         call pop_vector(tsl_vars, [t, dt, tot_mass, tot_momx, tot_momy, tot_momz, tot_ener, tot_epot, tot_eint, tot_ekin])
#ifdef MAGNETIC
         call pop_vector(tsl_vars, [tot_emag, tot_mflx, tot_mfly, tot_mflz, tsl%vai_max, tsl%b_min, tsl%b_max, tsl%divb_max])
#ifdef RESISTIVE
         if (eta1_active) call pop_vector(tsl_vars, [tsl%etamax])
#endif /* RESISTIVE */
#endif /* MAGNETIC */
#ifdef COSM_RAYS
         call pop_vector(tsl_vars, [tot_encr, tsl%encr_min, tsl%encr_max])
#endif /* COSM_RAYS */

         ! \todo: replicated code, simplify me
         if (has_ion) then
            sn=>flind%ion%snap
            call pop_vector(tsl_vars, [sn%dens_min%val, sn%dens_max%val, sn%velx_max%val, sn%vely_max%val, sn%velz_max%val, &
                                    sn%pres_min%val, sn%pres_max%val, sn%temp_min%val, sn%temp_max%val, sn%cs_max%val])
         endif
         if (has_neu) then
            sn=>flind%neu%snap
            call pop_vector(tsl_vars, [sn%dens_min%val, sn%dens_max%val, sn%velx_max%val, sn%vely_max%val, sn%velz_max%val, &
                                    sn%pres_min%val, sn%pres_max%val, sn%temp_min%val, sn%temp_max%val, sn%cs_max%val])
         endif
         if (has_dst) then
            sn=>flind%dst%snap
            call pop_vector(tsl_vars, [sn%dens_min%val, sn%dens_max%val, sn%velx_max%val, sn%vely_max%val, sn%velz_max%val])
         endif

      endif

      if (associated(user_tsl)) call user_tsl(tsl_vars)

      if (master) then
         write(tsl_lun, '(1x,i8,50(1x,es15.8))') nstep, tsl_vars

! some quantities computed in "write_log".One can add more, or change.
         close(tsl_lun)
         deallocate(tsl_vars)
      endif

   end subroutine write_timeslice

   subroutine get_extremum(tab,minmax,prop)

      use constants,     only: MINL, MAXL
      use dataio_pub,    only: msg, warn
      use types,         only: value
      use grid,          only: cg
      use mpisetup,      only: mpifind

      implicit none

      real, dimension(:,:,:), intent(in) :: tab
      integer, intent(in)                :: minmax
      type(value), intent(out)           :: prop

      select case (minmax)
         case (MINL)
            prop%val = minval(tab)
            prop%loc = minloc(tab) + [cg%nb, cg%nb, cg%nb]
         case (MAXL)
            prop%val = maxval(tab)
            prop%loc = maxloc(tab) + [cg%nb, cg%nb, cg%nb]
         case default
            write(msg,*) "[dataio:get_extremum]: I don't know what to do with minmax = ", minmax
            call warn(msg)
      end select

      call mpifind(prop, minmax)

   end subroutine get_extremum

   subroutine common_shout(pr, fluid, pres_tn, temp_tn, cs_tn)

      use constants,       only: small
      use dataio_pub,      only: msg, printinfo
      use grid,            only: cg
      use mpisetup,        only: cfl
      use fluidtypes,      only: phys_prop

      implicit none

      type(phys_prop), intent(in)  :: pr
      character(len=*), intent(in) :: fluid
      logical, intent(in)          :: pres_tn, temp_tn, cs_tn

      real :: dxmn_safe

      if (cg%dxmn >= sqrt(huge(1.0))) then
         dxmn_safe = sqrt(huge(1.0))
      else
         dxmn_safe = cg%dxmn
      endif

      write(msg, fmt_loc)     'min(dens)   ',fluid, pr%dens_min%val, pr%dens_min%proc, pr%dens_min%loc
      call printinfo(msg, .false.)
      write(msg, fmt_loc)     'max(dens)   ',fluid, pr%dens_max%val, pr%dens_max%proc, pr%dens_max%loc
      call printinfo(msg, .false.)
      if (temp_tn) then
         write(msg, fmt_loc)  'min(temp)   ',fluid, pr%temp_min%val, pr%temp_min%proc, pr%temp_min%loc
         call printinfo(msg, .false.)
         write(msg, fmt_loc)  'max(temp)   ',fluid, pr%temp_max%val, pr%temp_max%proc, pr%temp_max%loc
         call printinfo(msg, .false.)
      endif
      if (pres_tn) then
         write(msg, fmt_loc)  'min(pres)   ',fluid, pr%pres_min%val, pr%pres_min%proc, pr%pres_min%loc
         call printinfo(msg, .false.)
         write(msg, fmt_loc)  'max(pres)   ',fluid, pr%pres_max%val, pr%pres_max%proc, pr%pres_max%loc
         call printinfo(msg, .false.)
      endif

      write(msg, fmt_dtloc)   'max(|vx|)   ',fluid, pr%velx_max%val, cfl*cg%dx/(pr%velx_max%val+small), pr%velx_max%proc, pr%velx_max%loc
      call printinfo(msg, .false.)
      write(msg, fmt_dtloc)   'max(|vy|)   ',fluid, pr%vely_max%val, cfl*cg%dy/(pr%vely_max%val+small), pr%vely_max%proc, pr%vely_max%loc
      call printinfo(msg, .false.)
      write(msg, fmt_dtloc)   'max(|vz|)   ',fluid, pr%velz_max%val, cfl*cg%dz/(pr%velz_max%val+small), pr%velz_max%proc, pr%velz_max%loc
      call printinfo(msg, .false.)
      if (cs_tn) then
         write(msg, fmt_dtloc)'max(c_s )   ',fluid, pr%cs_max%val, cfl*dxmn_safe/(pr%cs_max%val+small), pr%cs_max%proc, pr%cs_max%loc
         call printinfo(msg, .false.)
      endif
   end subroutine common_shout

   subroutine get_common_vars(fl)

      use arrays,     only: u, b, wa
      use constants,  only: ION, DST, MINL, MAXL
      use mpisetup,   only: smallp
      use grid,       only: cg
      use fluidtypes, only: phys_prop, component_fluid
      use units,      only: mH, kboltz

      implicit none

      type(component_fluid), intent(inout), target :: fl
      type(phys_prop), pointer                     :: pr

      pr => fl%snap
      wa = u(fl%idn,:,:,:)
      call get_extremum(wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), MAXL, pr%dens_max)
      call get_extremum(wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), MINL, pr%dens_min)

      wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) = abs(u(fl%imx,cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke)/u(fl%idn,cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke))
      call get_extremum(wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), MAXL, pr%velx_max)

      wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) = abs(u(fl%imy,cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke)/u(fl%idn,cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke))
      call get_extremum(wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), MAXL, pr%vely_max)

      wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) = abs(u(fl%imz,cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke)/u(fl%idn,cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke))
      call get_extremum(wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), MAXL, pr%velz_max)

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
      pr%temp_min%val  = (mH * fl%cs2)/ (kboltz * fl%gam)
      pr%temp_min%loc  = 0
      pr%temp_min%proc = 0
      pr%temp_max%val  = (mH * fl%cs2)/ (kboltz * fl%gam)
      pr%temp_max%loc  = 0
      pr%temp_max%proc = 0
#else /* !ISO */
      if (fl%tag /= DST) then
         wa(:,:,:) = (u(fl%ien,:,:,:) &                ! eint
                   - 0.5*((u(fl%imx,:,:,:)**2 +u(fl%imy,:,:,:)**2 + u(fl%imz,:,:,:)**2)/u(fl%idn,:,:,:)))
         if (fl%tag == ION) wa(:,:,:) = wa(:,:,:) - 0.5*(sum(b**2,dim=1))

         wa(:,:,:) = max(fl%gam_1*wa(:,:,:),smallp)  ! pres

         call get_extremum(wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), MAXL, pr%pres_max)
         call get_extremum(wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), MINL, pr%pres_min)

         wa(:,:,:) = fl%gam*wa(:,:,:)/u(fl%idn,:,:,:) ! sound speed squared
         call get_extremum(wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), MAXL, pr%cs_max)
         pr%cs_max%val = sqrt(pr%cs_max%val)

         wa(:,:,:) = (mH * wa(:,:,:))/ (kboltz * fl%gam) ! temperature
         call get_extremum(wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), MAXL, pr%temp_max)
         call get_extremum(wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), MINL, pr%temp_min)
      endif
#endif /* !ISO */
   end subroutine get_common_vars
!---------------------------------------------------------------------
!>
!! \brief writes timestep diagnostics to the logfile
!!
!! \deprecated Quite costly routine due to extensive array searches
!<
!---------------------------------------------------------------------
!
   subroutine  write_log(tsl)

      use arrays,             only: wa, u, b
      use constants,          only: small, MINL, MAXL
      use dataio_pub,         only: msg, printinfo
      use fluids_pub,         only: has_dst, has_ion, has_neu
      use fluidindex,         only: ibx, iby, ibz, flind
      use grid,               only: cg
      use mpisetup,           only: cfl, t, dt, master
      use types,              only: tsl_container, value
      use interactions,       only: has_interactions, collfaq
      use func,               only: L2norm
#ifdef COSM_RAYS
      use fluidindex,         only: iarr_all_crs
      use timestepcosmicrays, only: dt_crs
#endif /* COSM_RAYS */
#ifdef RESISTIVE
      use resistivity,        only: dt_resist, etamax, cu2max, eta1_active
#ifndef ISO
      use resistivity,        only: deimin
#endif /* !ISO */
#endif /* RESISTIVE */
#ifdef VARIABLE_GP
      use arrays,             only: gpot
#endif /* VARIABLE_GP */
#if defined VARIABLE_GP || defined MAGNETIC
      use grid,               only: D_x, D_y, D_z
#endif /* VARIABLE_GP || MAGNETIC */

      implicit none

      type(tsl_container), optional  :: tsl
      real :: dxmn_safe

      type(value) :: drag
#ifdef MAGNETIC
      type(value) :: b_min, b_max, divb_max, vai_max
#endif /* MAGNETIC */
#ifdef COSM_RAYS
      type(value) :: encr_min, encr_max
#endif /* COSM_RAYS */
#ifdef VARIABLE_GP
      type(value) :: gpxmax, gpymax, gpzmax
#endif /* VARIABLE_GP */
      character(len=idlen) :: id
#if defined VARIABLE_GP || defined MAGNETIC
      integer :: nxl, nyl, nzl, nxu, nyu, nzu !< shortcuts for indices to compute all gradients regardless of dimensionality of the simulation

      nxl = 1 + D_x
      nyl = 1 + D_y
      nzl = 1 + D_z
      nxu = cg%nx - D_x
      nyu = cg%ny - D_y
      nzu = cg%nz - D_z
#endif /* VARIABLE_GP || MAGNETIC */

      id = '' ! suppress compiler warnings if noe of the modules requiring the id variable are swithed on.
      if (cg%dxmn >= sqrt(huge(1.0))) then
         dxmn_safe = sqrt(huge(1.0))
      else
         dxmn_safe = cg%dxmn
      endif

   ! Timestep diagnostics
      if (has_neu) call get_common_vars(flind%neu)

      if (has_ion) then
         call get_common_vars(flind%ion)

#ifdef MAGNETIC
         wa(:,:,:)  = sqrt(b(1,:,:,:)*b(1,:,:,:) + b(2,:,:,:)*b(2,:,:,:) + b(3,:,:,:)*b(3,:,:,:))
         call get_extremum(wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), MAXL, b_max)
         call get_extremum(wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), MINL, b_min)

         wa(:,:,:)  = wa(:,:,:) / sqrt(u(flind%ion%idn,:,:,:))
         call get_extremum(wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), MAXL, vai_max)
#endif /* MAGNETIC */

#ifdef ISO
#ifdef ISO_LOCAL
!        wa            = cs_iso2_arr(:,:,:)*u(idni,:,:,:)
!        prei_min%val  = minval(wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke))
!        prei_min%loc  = minloc(wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke)) + [cg%nb, cg%nb, cg%nb]
!        call mpifind(prei_min, MINL)
!
!        prei_max%val  = maxval(wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke))
!        prei_max%loc  = maxloc(wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke)) + [cg%nb, cg%nb, cg%nb]
!        call mpifind(prei_max, MAXL)
!
!        csi_max%val   = maxval(cs_iso2_arr(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke))
!        csi_max%loc   = maxloc(cs_iso2_arr(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke)) + [cg%nb, cg%nb, cg%nb]
!        call mpifind(csi_max, MAXL)
!
!        wa            = mH / kboltz * cs_iso2_arr(:,:,:)
!        temi_min%val  = minval(wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke))
!        temi_min%loc  = minloc(wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke)) + [cg%nb, cg%nb, cg%nb]
!        call mpifind(temi_min, MINL)
!
!        temi_max%val  = maxval(wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke))
!        temi_max%loc  = maxloc(wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke)) + [cg%nb, cg%nb, cg%nb]
!        call mpifind(temi_max, MAXL)
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
      endif

      if (has_dst) call get_common_vars(flind%dst)

#ifdef VARIABLE_GP
      wa(1:nxu,:,:) = abs((gpot(nxl:cg%nx,:,:)-gpot(1:nxu,:,:))*cg%idx) ; wa(cg%nx,:,:) = wa(nxu,:,:)
      call get_extremum(wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), MAXL, gpxmax)
      wa(:,1:nyu,:) = abs((gpot(:,nyl:cg%ny,:)-gpot(:,1:nyu,:))*cg%idy) ; wa(:,cg%ny,:) = wa(:,nyu,:)
      call get_extremum(wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), MAXL, gpymax)
      wa(:,:,1:nzu) = abs((gpot(:,:,nzl:cg%nz)-gpot(:,:,1:nzu))*cg%idz) ; wa(:,:,cg%nz) = wa(:,:,nzu)
      call get_extremum(wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), MAXL, gpzmax)
#endif /* VARIABLE_GP */

#ifdef MAGNETIC
      wa(1:nxu,1:nyu,1:nzu) = &
                 (b(ibx,nxl:cg%nx,  1:nyu  ,  1:nzu  ) - b(ibx,1:nxu,1:nyu,1:nzu))*cg%dy*cg%dz &
               + (b(iby,  1:nxu  ,nyl:cg%ny,  1:nzu  ) - b(iby,1:nxu,1:nyu,1:nzu))*cg%dx*cg%dz &
               + (b(ibz,  1:nxu  ,  1:nyu  ,nzl:cg%nz) - b(ibz,1:nxu,1:nyu,1:nzu))*cg%dx*cg%dy
      wa = abs(wa)

      wa(cg%ie,:,:) = wa(cg%ie-D_x,:,:)
      wa(:,cg%je,:) = wa(:,cg%je-D_y,:)
      wa(:,:,cg%ke) = wa(:,:,cg%ke-D_z)

      call get_extremum(wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), MAXL, divb_max)
#endif /* MAGNETIC */

#ifdef COSM_RAYS
      wa            = sum(u(iarr_all_crs,:,:,:),1)
      call get_extremum(wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), MAXL, encr_max)
      call get_extremum(wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), MINL, encr_min)
#endif /* COSM_RAYS */

      if (has_interactions) then
         wa = L2norm(u(flind%dst%imx,:,:,:),u(flind%dst%imy,:,:,:),u(flind%dst%imz,:,:,:),u(flind%neu%imx,:,:,:),u(flind%neu%imy,:,:,:),u(flind%neu%imz,:,:,:) ) * u(flind%dst%idn,:,:,:)
         call get_extremum(wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), MAXL, drag)
      endif

      if (master)  then
         if (.not.present(tsl)) then
            call printinfo('===========================================================================', .false.)
            if (has_ion) then
               call common_shout(flind%ion%snap,'ION',.true.,.true.,.true.)
#ifdef MAGNETIC
               id = "ION"
               write(msg, fmt_dtloc) 'max(c_f)    ', id, sqrt(flind%ion%snap%cs_max%val**2+vai_max%val**2), &
                    &             cfl*dxmn_safe/sqrt(flind%ion%snap%cs_max%val**2+vai_max%val**2+small)
               call printinfo(msg, .false.)
               write(msg, fmt_dtloc) 'max(v_a)    ', id, vai_max%val, cfl*dxmn_safe/(vai_max%val+small), vai_max%proc, vai_max%loc
               call printinfo(msg, .false.)
               id = "MAG"
               write(msg, fmt_loc)   'min(|b|)    ', id, b_min%val,     b_min%proc,     b_min%loc
               call printinfo(msg, .false.)
               write(msg, fmt_loc)   'max(|b|)    ', id, b_max%val,     b_max%proc,     b_max%loc
               call printinfo(msg, .false.)
               write(msg, fmt_loc)   'max(|divb|) ', id, divb_max%val,  divb_max%proc,  divb_max%loc
               call printinfo(msg, .false.)
#else /* !MAGNETIC */
!               if (csi_max%val > 0.) write(msg, fmtff8) 'max(c_s )   ION  =', sqrt(csi_max%val**2), 'dt=',cfl*dxmn_safe/sqrt(csi_max%val**2)
!               call printinfo(msg, .false.)
#endif /* !MAGNETIC */
            endif
            if (has_neu) call common_shout(flind%neu%snap,'NEU',.true.,.true.,.true.)
            if (has_dst) call common_shout(flind%dst%snap,'DST',.false.,.false.,.false.)
            if (has_interactions) then
               write(msg, fmt_dtloc) 'max(drag)   ', "INT", drag%val, flind%neu%cs/(maxval(collfaq) * drag%val + small), drag%proc, drag%loc
               call printinfo(msg, .false.)
            endif
#ifdef COSM_RAYS
            id = "CRS"
            write(msg, fmt_loc)   'min(encr)   ', id, encr_min%val, encr_min%proc, encr_min%loc
            call printinfo(msg, .false.)
            write(msg, fmt_dtloc) 'max(encr)   ', id, encr_max%val, dt_crs, encr_max%proc, encr_max%loc
            call printinfo(msg, .false.)
#endif /* COSM_RAYS */
#ifdef RESISTIVE
            if (eta1_active) then
               id = "RES"
               write(msg, fmt_dtloc) 'max(eta)    ', id, etamax%val, dt_resist, etamax%proc, etamax%loc
               call printinfo(msg, .false.)
               write(msg, fmt_dtloc) 'max(cu2)    ', id, cu2max%val, dt_resist, cu2max%proc, cu2max%loc
               call printinfo(msg, .false.)
#ifndef ISO
               write(msg, fmt_dtloc) 'min(dei)    ', id, deimin%val, dt_resist, deimin%proc, deimin%loc
               call printinfo(msg, .false.)
#endif /* !ISO */
            endif
#endif /* RESISTIVE */
#ifdef VARIABLE_GP
            id = "GPT"
            write(msg, fmt_loc)   'max(|gpx|)  ', id, gpxmax%val, gpxmax%proc, gpxmax%loc
            call printinfo(msg, .false.)
            write(msg, fmt_loc)   'max(|gpy|)  ', id, gpymax%val, gpymax%proc, gpymax%loc
            call printinfo(msg, .false.)
            write(msg, fmt_loc)   'max(|gpz|)  ', id, gpzmax%val, gpzmax%proc, gpzmax%loc
            call printinfo(msg, .false.)
#endif /* VARIABLE_GP */
            call printinfo('===========================================================================', .false.)
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
            if (eta1_active) tsl%etamax = etamax%val
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

!> \todo process multiple commands at once
      use constants,     only: cwdlen
      use dataio_pub,    only: msg, printinfo, warn
      use mpisetup,      only: master
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

            if (len_trim(msg) > 0 .and. master) call printinfo(msg)

            sz = len_trim(msg)
            if (fname(i) == user_message_file) unlink_stat = unlink(user_message_file)

         endif
      enddo

   end subroutine read_file_msg
!------------------------------------------------------------------------
end module dataio
