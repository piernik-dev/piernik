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
!! \brief Module containing all main routines  responsible for data output
!!
!!
!! In this module following namelists of parameters are specified:
!! \copydetails dataio::init_dataio
!<

module dataio

   use dataio_pub, only: domain_dump, fmin, fmax, vizit, nend, tend, wend, new_id, nrestart, problem_name, run_id, multiple_h5files, use_v2_io, nproc_io, enable_compression, gzip_level
   use constants,  only: cwdlen, fmt_len, cbuff_len, dsetnamelen, RES, TSL, ndims

   implicit none

   private
   public :: check_log, check_tsl, dump, write_data, write_crashed, cleanup_dataio, init_dataio, init_dataio_parameters, user_msg_handler

   integer                  :: istep                 !< current number of substep (related to integration order)

   integer, parameter       :: nvarsmx = 20          !< maximum number of variables to dump in hdf files
   character(len=cbuff_len) :: restart               !< choice of restart %file: if restart = 'last': automatic choice of the last restart file regardless of "nrestart" value; if something else is set: "nrestart" value is fixing
   logical                  :: mag_center            !< choice to dump magnetic fields values from cell centers or not (if not then values from cell borders, unused)
   integer(kind=4)          :: resdel                !< number of recent restart dumps which should be saved; each n-resdel-1 restart file is supposed to be deleted while writing n restart file
   real                     :: dt_hdf                !< time between successive hdf dumps
   real                     :: dt_res                !< time between successive restart file dumps
   real                     :: dt_tsl                !< time between successive timeslice dumps
   real                     :: dt_log                !< time between successive log dumps
   real                     :: dt_plt                !< time between successive domain slices files dumps
   character(len=cwdlen)    :: user_message_file     !< path to possible user message file containing dt_xxx changes or orders to dump/stop/end simulation
   character(len=cwdlen)    :: system_message_file   !< path to possible system (UPS) message file containing orders to dump/stop/end simulation
   integer(kind=4), dimension(ndims) :: plt_plane    !< indices of cells that are sliced in plt files
   integer                  :: iv                    !< work index to count successive variables to dump in hdf files
   character(len=dsetnamelen), dimension(nvarsmx) :: vars !< array of 4-character strings standing for variables to dump in hdf files

   integer                  :: tsl_lun               !< logical unit number for timeslice file
   integer                  :: nhdf_start            !< number of hdf file for the first hdf dump in simulation run
   integer                  :: nres_start            !< number of restart file for the first restart dump in simulation run
   real                     :: t_start               !< time in simulation of start simulation run
   logical                  :: tsl_firstcall         !< logical value to start a new timeslice file
   logical                  :: initial_hdf_dump      !< force initial hdf dump
   logical, dimension(RES:TSL) :: dump = .false.     !< logical values for all dump types to restrict to only one dump of each type a step

   integer                  :: nchar                 !< number of characters in a user/system message
   integer, parameter       :: umsg_len = 16
   character(len=umsg_len)  :: umsg                  !< string of characters - content of a user/system message
   real                     :: umsg_param            !< parameter changed by a user/system message

   character(len=cwdlen)    :: filename              !< string of characters indicating currently used file
   character(len=fmt_len), protected, target :: fmt_loc, fmt_dtloc, fmt_vloc

   type :: tsl_container
      logical :: dummy
#ifdef COSM_RAYS
      real :: encr_min, encr_max
#endif /* COSM_RAYS */
#ifdef RESISTIVE
      real :: etamax
#endif /* RESISTIVE */
#ifdef MAGNETIC
      real :: b_min, b_max, divb_max, vai_max
#endif /* MAGNETIC */
#ifdef VARIABLE_GP
      real :: gpxmax, gpymax, gpzmax
#endif /* VARIABLE_GP */
   end type tsl_container

   namelist /END_CONTROL/     nend, tend, wend
   namelist /RESTART_CONTROL/ restart, new_id, nrestart, resdel
   namelist /OUTPUT_CONTROL/  problem_name, run_id, dt_hdf, dt_res, dt_tsl, dt_log, dt_plt, plt_plane, &
                              domain_dump, vars, mag_center, vizit, fmin, fmax, user_message_file, system_message_file, &
                              multiple_h5files, use_v2_io, nproc_io, enable_compression, gzip_level, initial_hdf_dump

contains

!>
!! \brief Routine that sets the initial values of data input/output parameters from namelists @c END_CONTROL, @c RESTART_CONTROL and @c OUTPUT_CONTROL.
!! Called as early as possible.
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
!! <tr><td>restart </td><td>'last'</td><td>'last' or another string of characters</td><td>\copydoc dataio::restart     </td></tr>
!! <tr><td>new_id  </td><td>''    </td><td>string of characters                  </td><td>\copydoc dataio_pub::new_id  </td></tr>
!! <tr><td>nrestart</td><td>3     </td><td>integer                               </td><td>\copydoc dataio_pub::nrestart</td></tr>
!! <tr><td>resdel  </td><td>0     </td><td>integer                               </td><td>\copydoc dataio::resdel      </td></tr>
!! </table>
!! \n \n
!! @b OUTPUT_CONTROL
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>problem_name       </td><td>"nameless"         </td><td>real      </td><td>\copydoc dataio_pub::problem_name </td></tr>
!! <tr><td>run_id             </td><td>"___"              </td><td>real      </td><td>\copydoc dataio_pub::run_id       </td></tr>
!! <tr><td>dt_hdf             </td><td>0.0                </td><td>real      </td><td>\copydoc dataio::dt_hdf           </td></tr>
!! <tr><td>dt_res             </td><td>0.0                </td><td>real      </td><td>\copydoc dataio::dt_res           </td></tr>
!! <tr><td>dt_tsl             </td><td>0.0                </td><td>real      </td><td>\copydoc dataio::dt_tsl           </td></tr>
!! <tr><td>dt_log             </td><td>0.0                </td><td>real      </td><td>\copydoc dataio::dt_log           </td></tr>
!! <tr><td>dt_plt             </td><td>0.0                </td><td>real      </td><td>\copydoc dataio::dt_plt           </td></tr>
!! <tr><td>plt_plane          </td><td>(nxd, nyd, nzd)/2  </td><td>integer(3)</td><td>\copydoc dataio::plt_plane        </td></tr>
!! <tr><td>domain_dump        </td><td>'phys_domain'      </td><td>'phys_domain' or 'full_domain'                       </td><td>\copydoc dataio_pub::domain_dump</td></tr>
!! <tr><td>vars               </td><td>''                 </td><td>'dens', 'velx', 'vely', 'velz', 'ener' and some more </td><td>\copydoc dataio::vars  </td></tr>
!! <tr><td>mag_center         </td><td>.false.            </td><td>logical   </td><td>\copydoc dataio::mag_center       </td></tr>
!! <tr><td>vizit              </td><td>.false.            </td><td>logical   </td><td>\copydoc dataio_pub::vizit        </td></tr>
!! <tr><td>fmin               </td><td>                   </td><td>real      </td><td>\copydoc dataio_pub::fmin         </td></tr>
!! <tr><td>fmax               </td><td>                   </td><td>real      </td><td>\copydoc dataio_pub::fmax         </td></tr>
!! <tr><td>user_message_file  </td><td>trim(wd_rd)//'/msg'</td><td>string similar to default value              </td><td>\copydoc dataio::user_message_file  </td></tr>
!! <tr><td>system_message_file</td><td>'/tmp/piernik_msg' </td><td>string of characters similar to default value</td><td>\copydoc dataio::system_message_file</td></tr>
!! <tr><td>multiple_h5files   </td><td>.false.            </td><td>logical   </td><td>\copydoc dataio_pub::multiple_h5files</td></tr>
!! <tr><td>use_v2_io          </td><td>.false.            </td><td>logical   </td><td>\copydoc dataio_pub::use_v2_io    </td></tr>
!! <tr><td>nproc_io           </td><td>1                  </td><td>integer   </td><td>\copydoc dataio_pub::nproc_io     </td></tr>
!! <tr><td>enable_compression </td><td>.false.            </td><td>logical   </td><td>\copydoc dataio_pub::enable_compression</td></tr>
!! <tr><td>gzip_level         </td><td>9                  </td><td>integer   </td><td>\copydoc dataio_pub::gzip_level   </td></tr>
!! </table>
!! \n \n
!<
   subroutine init_dataio_parameters

      use constants,  only: idlen, cwdlen, cbuff_len, PIERNIK_INIT_MPI, I_ONE, I_TWO
      use dataio_pub, only: nres, nrestart, last_hdf_time, last_plt_time, last_res_time, last_tsl_time, last_log_time, log_file_initialized, &
           &                tmp_log_file, printinfo, printio, warn, msg, nhdf, nimg, die, code_progress, wd_wr, wd_rd, &
           &                move_file, multiple_h5files, parfile, parfilelines, log_file, maxparfilelines, can_i_write
      use dataio_pub, only: par_file, ierrh, namelist_errh, compare_namelist, cmdl_nml, lun  ! QA_WARN required for diff_nml
      use domain,     only: dom
      use mpisetup,   only: lbuff, ibuff, rbuff, cbuff, master, slave, nproc, proc, piernik_MPI_Bcast, piernik_MPI_Barrier

      implicit none

      integer              :: system_status, i, par_lun

#ifdef VERBOSE
      if (master) call printinfo("[dataio:init_dataio_parameters] Commencing dataio module initialization")
#endif /* VERBOSE */

      if (code_progress < PIERNIK_INIT_MPI) call die("[dataio:init_dataio_parameters] Some physics modules are not initialized.")

      problem_name = "nameless"
      run_id       = "___"
      restart      = 'last'   ! 'last': automatic choice of the last restart file regardless of "nrestart" value;
                              ! if something else is set: "nrestart" value is fixing
      new_id       = ''
      nrestart     = 3
      resdel       = 0

      dt_hdf       = 0.0
      dt_res       = 0.0
      dt_tsl       = 0.0
      dt_log       = 0.0
      dt_plt       = 0.0

      plt_plane = max(I_ONE, dom%n_d(:)/I_TWO)

      initial_hdf_dump = .false.

      domain_dump       = 'phys_domain'
      vars(:)           = ''
      mag_center        = .false.
      write(user_message_file,'(a,"/msg")') trim(wd_rd)
      system_message_file = "/tmp/piernik_msg"

      tsl_firstcall      = .true.
      use_v2_io          = .false.
      nproc_io           = 1
      enable_compression = .false.
      gzip_level         = 9

      nhdf  = -1
      nimg  = -1
      nres  = 0

      nend = 1
      tend = -1.0
      wend = huge(1.0)

      if (master) then

         open(newunit=par_lun,file=par_file)
         ierrh = 0
         do while (ierrh == 0 .and. parfilelines<maxparfilelines)
            read(unit=par_lun, fmt='(a)', iostat=ierrh) parfile(parfilelines+1)
            if (ierrh == 0) then
               parfilelines = parfilelines + 1
               i = len_trim(parfile(parfilelines))
               if (i >= len(parfile(parfilelines))) call warn("[dataio:init_dataio_parameters] problem.par contains very long lines. The copy in the logfile and HDF dumps can be truncated.")
            endif
         enddo
         close(par_lun)
         if (parfilelines == maxparfilelines) call warn("[dataio:init_dataio_parameters] problem.par has too many lines. The copy in the logfile and HDF dumps can be truncated.")

         diff_nml(OUTPUT_CONTROL)
         diff_nml(RESTART_CONTROL)
         diff_nml(END_CONTROL)

         if (use_v2_io) then
            if (nproc_io <= 0 .or. nproc_io > nproc) nproc_io = nproc ! fully parallel v2 I/O

            if (nproc_io /= 1) &
               call warn("[dataio:init_dataio_parameters] Parallel v2 I/O is experimental feature")
         endif

         if (gzip_level < 1 .or. gzip_level > 9) then
            call warn("[dataio:init_dataio_parameters] invalid compression level")
            gzip_level = 9
         endif

!  namelist /END_CONTROL/ nend, tend, wend
         ibuff(1)  = nend

         rbuff(1)  = tend
         rbuff(2)  = wend

!  namelist /RESTART_CONTROL/ restart, new_id, nrestart, resdel

         cbuff(20) = restart
         cbuff(21) = new_id

         ibuff(20) = nrestart
         ibuff(21) = resdel

!   namelist /OUTPUT_CONTROL/ problem_name, run_id, dt_hdf, dt_res, dt_tsl, dt_log, dt_plt, plt_plane, &
!                             domain_dump, vars, mag_center, vizit, fmin, fmax, &
!                             user_message_file, system_message_file, multiple_h5files, use_v2_io, nproc_io

         ibuff(40:42) = plt_plane
         ibuff(43) = nproc_io
         ibuff(44) = gzip_level

         rbuff(40) = dt_hdf
         rbuff(41) = dt_res
         rbuff(42) = dt_tsl
         rbuff(43) = dt_log
         rbuff(44) = dt_plt
         rbuff(45) = fmin
         rbuff(46) = fmax

         lbuff(1)  = vizit
         lbuff(2)  = multiple_h5files
         lbuff(3)  = use_v2_io
         lbuff(4)  = mag_center
         lbuff(5)  = initial_hdf_dump

         cbuff(31) = problem_name
         cbuff(32) = run_id
         cbuff(40) = domain_dump

         do iv = 1, nvarsmx
            cbuff(40+iv) = vars(iv)
         enddo

         cbuff(90) = user_message_file(1:cbuff_len)
         cbuff(91) = system_message_file(1:cbuff_len)

      endif

      call piernik_MPI_Bcast(cbuff, cbuff_len)
      call piernik_MPI_Bcast(lbuff)
      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(rbuff)

      if (slave) then

!  namelist /END_CONTROL/ nend, tend, wend
         nend                = ibuff(1)

         tend                = rbuff(1)
         wend                = rbuff(2)
!  namelist /RESTART_CONTROL/ restart, new_id, nrestart, resdel

         restart             = trim(cbuff(20))
         new_id              = trim(cbuff(21))

         nrestart            = int(ibuff(20), kind=4)
         resdel              = ibuff(21)

!  namelist /OUTPUT_CONTROL/ dt_hdf, dt_res, dt_tsl, domain, vars, mag_center, &
!                            ix, iy, iz, user_message_file, system_message_file

         plt_plane           = ibuff(40:42)
         nproc_io            = int(ibuff(43), kind=4)
         gzip_level          = int(ibuff(44), kind=4)

         dt_hdf              = rbuff(40)
         dt_res              = rbuff(41)
         dt_tsl              = rbuff(42)
         dt_log              = rbuff(43)
         dt_plt              = rbuff(44)
         fmin                = rbuff(45)
         fmax                = rbuff(46)

         vizit               = lbuff(1)
         multiple_h5files    = lbuff(2)
         use_v2_io           = lbuff(3)
         mag_center          = lbuff(4)
         initial_hdf_dump    = lbuff(5)

         problem_name        = cbuff(31)
         run_id              = cbuff(32)(1:idlen)
         domain_dump         = trim(cbuff(40))
         do iv=1, nvarsmx
            vars(iv)         = trim(cbuff(40+iv))
         enddo

         user_message_file   = trim(cbuff(90))
         system_message_file = trim(cbuff(91))

      endif

      can_i_write = mod( proc*nproc_io, nproc) < nproc_io
      if (can_i_write) then
         write(msg,'(a,i6,a)')"Process ",proc," can write"
         call printio(msg)
      endif

      last_log_time = -dt_log
      last_tsl_time = -dt_tsl
      last_hdf_time = -dt_hdf
      last_plt_time = -dt_plt
      last_res_time = 0.0

      if (master .and. restart == 'last') call find_last_restart(nrestart)
      call piernik_MPI_Barrier
      call piernik_MPI_Bcast(nrestart)

      if (master) then
         write(log_file,'(6a,i3.3,a)') trim(wd_wr),'/',trim(problem_name),'_',trim(run_id),'_',nrestart,'.log'
!> \todo if the simulation is restarted then save previous log_file (if exists) under a different, unique name
         system_status = move_file(trim(tmp_log_file), trim(log_file))
         if (system_status /= 0) then
            write(msg,'(2a)')"[dataio:init_dataio_parameters] The log must be stored in ",tmp_log_file
            call warn(msg)
            log_file_initialized = .false.
         else
            log_file_initialized = .true.
         endif
      endif
      call piernik_MPI_Bcast(log_file, cwdlen)          ! BEWARE: every msg issued by slaves before this sync may lead to race condition on tmp_log_file
      call piernik_MPI_Bcast(log_file_initialized)

   end subroutine init_dataio_parameters

!> \brief Initialize these I/O variables that may depend on any other modules (called at the end of init_piernik)

   subroutine init_dataio

      use common_hdf5,  only: init_hdf5
      use constants,    only: PIERNIK_INIT_IO_IC
      use data_hdf5,    only: init_data
      use dataio_pub,   only: nres, nrestart, printinfo, nhdf, nstep_start, die, code_progress
      use domain,       only: dom
      use global,       only: t, nstep
      use mpisetup,     only: master
      use restart_hdf5, only: read_restart_hdf5
      use slice_hdf5,   only: init_plot
      use timer,        only: time_left
      use user_hooks,   only: user_vars_arr_in_restart
      use version,      only: nenv,env, init_version

      implicit none

      logical :: tn
      integer :: i

      if (code_progress < PIERNIK_INIT_IO_IC) call die("[dataio:init_dataio] Some physics modules are not initialized.")

      write(fmt_loc,  '(2(a,i1),a)') "(2x,a12,a3,'  = ',es16.9,16x,            ",dom%eff_dim+1,"(1x,i4),",dom%eff_dim,"(1x,f12.4))"
      write(fmt_dtloc,'(2(a,i1),a)') "(2x,a12,a3,'  = ',es16.9,'  dt=',es11.4, ",dom%eff_dim+1,"(1x,i4),",dom%eff_dim,"(1x,f12.4))"
      write(fmt_vloc, '(2(a,i1),a)') "(2x,a12,a3,'  = ',es16.9,'   v=',es11.4, ",dom%eff_dim+1,"(1x,i4),",dom%eff_dim,"(1x,f12.4))"

      if (master) tn = time_left(wend)

      call init_hdf5(vars)
      call init_data
      call init_plot( plt_plane, dt_plt)

      call init_version
      if (master) then
         call printinfo("###############     Source configuration     ###############", .false.)
         do i=1,nenv
            call printinfo(env(i), .false.)
         enddo
      endif

      if (associated(user_vars_arr_in_restart)) call user_vars_arr_in_restart

      nres = nrestart

      if (nrestart /= 0) then
         if (master) call printinfo("###############     Reading restart     ###############", .false.)
         call read_restart_hdf5
         nstep_start = nstep
         t_start     = t
         nres_start  = nrestart
         nhdf_start  = nhdf-1
         if (new_id /= '') run_id=new_id
      endif

#ifdef VERBOSE
      call printinfo("[dataio:init_dataio] finished. \o/")
#endif /* VERBOSE */

   end subroutine init_dataio

   subroutine cleanup_dataio
      use common_hdf5,     only: cleanup_hdf5
      implicit none

      call cleanup_hdf5
   end subroutine cleanup_dataio

   subroutine user_msg_handler(end_sim)

      use data_hdf5,    only: write_hdf5
      use dataio_pub,   only: msg, printinfo, warn
      use mpisetup,     only: master, piernik_MPI_Bcast
      use restart_hdf5, only: write_restart_hdf5
      use timer,        only: time_left

      implicit none

      logical, intent(inout) :: end_sim
      logical                :: tn
      integer                :: tsleep

!--- process 0 checks for messages

      if (master) call read_file_msg

      call piernik_MPI_Bcast(umsg, umsg_len)
      call piernik_MPI_Bcast(umsg_param)

!---  if a user message is received then:
      if (len_trim(umsg) /= 0) then
         select case (trim(umsg))
            case ('res', 'dump')
               call write_restart_hdf5
            case ('hdf')
               call write_hdf5
            case ('log')
               call write_log
            case ('tsl')
               call write_timeslice
            case ('wend')
               wend = umsg_param
               if (master) tn = time_left(wend)
            case ('wleft')
               if (master) tn = time_left(-1.0)
            case ('tend')
               tend   = umsg_param
            case ('nend')
               nend   = int(umsg_param, kind=4)
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
                  &"  wleft    - show how much walltime is left",char(10),&
                  &"  sleep <number> - wait <number> seconds",char(10),&
                  &"  wend|tend|nend|dtres|dthdf|dtlog|dttsl|dtplt <value> - update specified parameter with <value>",char(10),&
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

      use constants,  only: FINAL
      use dataio_pub, only: nres, die

      implicit none

      character(len=*), intent(in) :: msg

      ! force output for diagnostics
      problem_name = "crash"
      dt_hdf = tiny(1.0)
      nres = 1
      call write_data(output=FINAL)

      call die(msg)

   end subroutine write_crashed

!---------------------------------------------------------------------
!>
!! \brief controls data dumping
!<
!---------------------------------------------------------------------
!
   subroutine write_data(output)

      use constants,    only: FINAL, HDF, PLT, LOGF
      use data_hdf5,    only: write_hdf5
      use dataio_pub,   only: last_res_time, last_hdf_time, last_plt_time
      use dataio_user,  only: user_post_write_data
      use restart_hdf5, only: write_restart_hdf5
      use slice_hdf5,   only: write_plot

      implicit none

      integer(kind=4), intent(in) :: output

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      dump(LOGF) = (output == LOGF .or. output == FINAL) ; if (dump(LOGF)) call write_log
      dump(TSL)  = (output == TSL  .or. output == FINAL) ; if (dump(TSL))  call write_timeslice

      call determine_dump(dump(RES), last_res_time, dt_res, output, RES)
      if (dump(RES)) call write_restart_hdf5

      call determine_dump(dump(PLT), last_plt_time, dt_plt, output, PLT)
      if (dump(PLT)) call write_plot

      call determine_dump(dump(HDF), last_hdf_time, dt_hdf, output, HDF)
      call manage_hdf_dump(dump(HDF), output)
      if (dump(HDF)) call write_hdf5

      if (associated(user_post_write_data)) call user_post_write_data(output, dump)

   end subroutine write_data

   subroutine determine_dump(dmp, last_dump_time, dt_dump, output, dumptype)

      use constants, only: FINAL
      use global,    only: t

      implicit none

      integer(kind=4), intent(in)    :: output, dumptype
      real,            intent(in)    :: dt_dump
      real,            intent(inout) :: last_dump_time
      logical,         intent(inout) :: dmp

      dmp = ((.not.(dmp)) .and. (output == FINAL)) !! dmp==.true. means that this dump is already done due to another conditions and is stopped by dmp=.not.(dmp); important only for FINAL output
      dmp = (dmp .or. (t-last_dump_time) >= dt_dump)
      dmp = (dmp .and. dt_dump > 0.0)
      if (dmp) last_dump_time = last_dump_time + real(floor((t-last_dump_time)/dt_dump))*dt_dump
      dmp = (dmp .or. output == dumptype)

   end subroutine determine_dump

   subroutine manage_hdf_dump(dmp, output)

      use constants,    only: FINAL, INCEPTIVE

      implicit none

      integer(kind=4), intent(in)    :: output  !< type of output
      logical,         intent(inout) :: dmp     !< perform I/O if True

      if (output == FINAL .and. trim(problem_name) /= 'crash') write(problem_name, '(a,a6)') trim(problem_name), '_final'
      if ((output == INCEPTIVE) .and. initial_hdf_dump) dmp = .true.  !< \todo problem_name may be enhanced by '_initial', but this and nhdf should be reverted just after write_hdf5 is called

   end subroutine manage_hdf_dump

   subroutine check_log

      use constants,  only: CHK, LOGF
      use dataio_pub, only: last_log_time

      implicit none

      call determine_dump(dump(LOGF), last_log_time, dt_log, CHK, LOGF)
      if (dump(LOGF)) call write_log

   end subroutine check_log

   subroutine check_tsl

      use mpisetup,   only: report_to_master
      use mpisignals, only: sig
      use constants,  only: CHK
      use dataio_pub, only: last_tsl_time

      implicit none

      call determine_dump(dump(TSL), last_tsl_time, dt_tsl, CHK, TSL)
      if (dump(TSL)) then
         call write_timeslice
         call report_to_master(sig%tsl_updated, only_master=.True.)
      endif

   end subroutine check_tsl

!>
!! \brief Find the restart point with highest number
!!
!! \todo use restart_fname() function
!! \todo scan the 9999 .. 0 range somewhat smarter (get directory listing?)
!<

   subroutine find_last_restart(restart_number)

      use common_hdf5,   only: output_fname
      use constants,     only: RD
#if defined(__INTEL_COMPILER)
      use ifport,        only: unlink
#endif /* __INTEL_COMPILER */

      implicit none

#if defined(__PGI)
      include "lib3f.h"
#endif /* __PGI */

      integer(kind=4), intent(out) :: restart_number

      integer(kind=4)              :: nres
      integer                      :: unlink_stat
      logical                      :: exist

      restart_number = 0

      unlink_stat = unlink('restart_list.tmp')

      do nres = 999, 0, -1
         inquire(file = trim(output_fname(RD,'.res', nres)), exist = exist)
         if (exist) then
            restart_number = nres
            return
         endif
      enddo

   end subroutine find_last_restart

!>
!! \brief writes integrals to text file
!<

   subroutine write_timeslice

      use constants,      only: cwdlen, xdim, ydim, zdim, DST
      use dataio_pub,     only: wd_wr
      use dataio_user,    only: user_tsl
      use diagnostics,    only: pop_vector
      use domain,         only: dom
      use fluidindex,     only: flind, iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz
      use fluids_pub,     only: has_ion, has_dst, has_neu
      use fluidtypes,     only: phys_prop
      use func,           only: ekin, emag
      use cg_list,        only: cg_list_element
      use global,         only: t, dt, smalld, nstep
      use cg_leaves,      only: leaves
      use grid_cont,      only: grid_container
      use mass_defect,    only: update_magic_mass
      use mpi,            only: MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_SUM
      use mpisetup,       only: master, comm, mpi_err
      use named_array_list, only: wna
#ifdef GRAV
      use constants,      only: gpot_n
      use named_array_list, only: qna
#endif /* GRAV */
#ifndef ISO
      use fluidindex,     only: iarr_all_en
#endif /* !ISO */
#ifdef COSM_RAYS
      use fluidindex,     only: iarr_all_crs
#endif /* COSM_RAYS */
#ifdef RESISTIVE
      use resistivity,    only: eta1_active
#endif /* RESISTIVE */

      implicit none

      character(len=cwdlen)                               :: tsl_file, head_fmt
      character(len=cbuff_len), dimension(:), allocatable :: tsl_names
      real,                     dimension(:), allocatable :: tsl_vars
      real, dimension(:,:,:,:), pointer                   :: pu, pb
      type(phys_prop),          pointer                   :: sn
      type(tsl_container)                                 :: tsl
      type(grid_container),     pointer                   :: cg
      type(cg_list_element),    pointer                   :: cgl
      real                                                :: cs_iso2
      enum, bind(C)
         enumerator :: T_MASS                                  !< total mass
         enumerator :: T_MOMX, T_MOMY, T_MOMZ                  !< total momenta
         enumerator :: T_ENER, T_EINT, T_EKIN, T_EMAG, T_EPOT  !< total energies
         enumerator :: T_MFLX, T_MFLY, T_MFLZ                  !< total magnetic fluxes
#ifdef COSM_RAYS
         enumerator :: T_ENCR                                  !< total CR energy
#endif /* COSM_RAYS */
         enumerator :: T_LAST                                  !< DO NOT place any index behind this one
      end enum
      real, dimension(T_MASS:T_LAST-1), save :: tot_q          !< array of total quantities
      integer(kind=4)                        :: ifl

      if (has_ion) then
         cs_iso2 = flind%ion%cs2
      elseif (has_neu) then
         cs_iso2 = flind%neu%cs2
      else
         cs_iso2 = 0.0
      endif

      if (master) then
         write(tsl_file,'(a,a1,a,a1,a3,a1,i3.3,a4)') trim(wd_wr),'/',trim(problem_name),'_', run_id,'_',nrestart,'.tsl'

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
            if (has_ion) then
               call pop_vector(tsl_names, cbuff_len, ["deni_min", "deni_max", "vxi_max ", "vyi_max ", "vzi_max ", "prei_min", "prei_max", "temi_min", "temi_max", "csi_max "])
               call pop_vector(tsl_names, cbuff_len, ["ion_mmass_cur", "ion_mmass_cum"])
            endif
            if (has_neu) then
               call pop_vector(tsl_names, cbuff_len, ["denn_min", "denn_max", "vxn_max ", "vyn_max ", "vzn_max ", "pren_min", "pren_max", "temn_min", "temn_max", "csn_max "])
               call pop_vector(tsl_names, cbuff_len, ["neu_mmass_cur", "neu_mmass_cum"])
            endif
            if (has_dst) then
               call pop_vector(tsl_names, cbuff_len, ["dend_min", "dend_max", "vxd_max ", "vyd_max ", "vzd_max "])
               call pop_vector(tsl_names, cbuff_len, ["dst_mmass_cur", "dst_mmass_cum"])
            endif

            if (associated(user_tsl)) call user_tsl(tsl_vars, tsl_names)
            write(head_fmt,'(A,I2,A)') "(a1,a8,",size(tsl_names)-1,"a16)"

            open(newunit=tsl_lun, file=tsl_file)
            write(tsl_lun,fmt=head_fmt) "#",tsl_names
            write(tsl_lun, '(a1)') '#'
            deallocate(tsl_names)
            tsl_firstcall = .false.
         else
            open(newunit=tsl_lun, file=tsl_file, position='append')
         endif
      endif

      tot_q(:) = 0.
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         pu => cg%w(wna%fi)%span(cg%ijkse)
         pb => cg%w(wna%bi)%span(cg%ijkse)

         tot_q(T_MASS) = tot_q(T_MASS) + cg%dvol * sum(pu(iarr_all_dn,:,:,:))
         tot_q(T_MOMX) = tot_q(T_MOMX) + cg%dvol * sum(pu(iarr_all_mx,:,:,:))
         tot_q(T_MOMY) = tot_q(T_MOMY) + cg%dvol * sum(pu(iarr_all_my,:,:,:))
         tot_q(T_MOMZ) = tot_q(T_MOMZ) + cg%dvol * sum(pu(iarr_all_mz,:,:,:))
#ifdef GRAV
         tot_q(T_EPOT) = tot_q(T_EPOT) + cg%dvol * sum(sum(pu(iarr_all_dn(:),:,:,:),dim=1) * cg%q(qna%ind(gpot_n))%span(cg%ijkse))
#endif /* GRAV */

         tot_q(T_EKIN) = tot_q(T_EKIN) + cg%dvol * sum(ekin(pu(iarr_all_mx(:),:,:,:), pu(iarr_all_my(:),:,:,:), pu(iarr_all_mz(:),:,:,:), max(pu(iarr_all_dn(:),:,:,:),smalld)))
         tot_q(T_EMAG) = tot_q(T_EMAG) + cg%dvol * sum(emag(pb(xdim,:,:,:), pb(ydim,:,:,:), pb(zdim,:,:,:)))

         tot_q(T_MFLX) = tot_q(T_MFLX) + cg%dvol/dom%L_(xdim) * sum(pb(xdim,:,:,:)) !cg%dy*cg%dz/dom%n_d(xdim)
         tot_q(T_MFLY) = tot_q(T_MFLY) + cg%dvol/dom%L_(ydim) * sum(pb(ydim,:,:,:)) !cg%dx*cg%dz/dom%n_d(ydim)
         tot_q(T_MFLZ) = tot_q(T_MFLZ) + cg%dvol/dom%L_(zdim) * sum(pb(zdim,:,:,:)) !cg%dx*cg%dy/dom%n_d(zdim)
#ifndef ISO
         tot_q(T_ENER) = tot_q(T_ENER) + cg%dvol * sum(pu(iarr_all_en,:,:,:))
#endif /* !ISO */

#ifdef COSM_RAYS
         tot_q(T_ENCR) = tot_q(T_ENCR) + cg%dvol * sum(pu(iarr_all_crs,:,:,:))
         tot_q(T_ENER) = tot_q(T_ENER) + tot_q(T_ENCR)
#endif /* COSM_RAYS */

         cgl => cgl%nxt
      enddo

#ifdef ISO
      tot_q(T_EINT) = tot_q(T_EINT) + cs_iso2*tot_q(T_MASS)
      tot_q(T_ENER) = tot_q(T_ENER) + tot_q(T_EINT)+tot_q(T_EKIN)+tot_q(T_EMAG)
#else /* !ISO */
      tot_q(T_EINT) = tot_q(T_EINT) + tot_q(T_ENER) - tot_q(T_EKIN) - tot_q(T_EMAG)
#endif /* !ISO */
#ifdef GRAV
      tot_q(T_ENER) = tot_q(T_ENER) + tot_q(T_EPOT)
#endif /* GRAV */

      call MPI_Allreduce(MPI_IN_PLACE, tot_q(:), size(tot_q), MPI_DOUBLE_PRECISION, MPI_SUM, comm, mpi_err)

      call write_log(tsl)
      call update_magic_mass(tsl=.true.)

      if (master) then
         call pop_vector(tsl_vars, [t, dt, tot_q(T_MASS), tot_q(T_MOMX), tot_q(T_MOMY), tot_q(T_MOMZ), tot_q(T_ENER), tot_q(T_EPOT), tot_q(T_EINT), tot_q(T_EKIN)])
#ifdef MAGNETIC
         call pop_vector(tsl_vars, [tot_q(T_EMAG), tot_q(T_MFLX), tot_q(T_MFLY), tot_q(T_MFLZ), tsl%vai_max, tsl%b_min, tsl%b_max, tsl%divb_max])
#ifdef RESISTIVE
         if (eta1_active) call pop_vector(tsl_vars, [tsl%etamax])
#endif /* RESISTIVE */
#endif /* MAGNETIC */
#ifdef COSM_RAYS
         call pop_vector(tsl_vars, [tot_q(T_ENCR), tsl%encr_min, tsl%encr_max])
#endif /* COSM_RAYS */

         do ifl = lbound(flind%all_fluids, 1, kind=4), ubound(flind%all_fluids, 1, kind=4)
            sn => flind%all_fluids(ifl)%fl%snap
            call pop_vector(tsl_vars, [sn%dens_min%val, sn%dens_max%val, sn%velx_max%val, sn%vely_max%val, sn%velz_max%val])
            if (flind%all_fluids(ifl)%fl%tag /= DST) then
               call pop_vector(tsl_vars, [sn%pres_min%val, sn%pres_max%val, sn%temp_min%val, sn%temp_max%val, sn%cs_max%val])
            endif
            call pop_vector(tsl_vars, [sn%mmass_cur, sn%mmass_cum])
         enddo

      endif

      if (associated(user_tsl)) call user_tsl(tsl_vars)

      if (master) then
         write(tsl_lun, '(1x,i8,100(1x,es15.8))') nstep, tsl_vars

! some quantities computed in "write_log".One can add more, or change.
         close(tsl_lun)
         deallocate(tsl_vars)
      endif

   end subroutine write_timeslice

   subroutine common_shout(pr, fluid, pres_tn, temp_tn, cs_tn)

      use domain,      only: is_multicg
      use fluidtypes,  only: phys_prop

      implicit none

      type(phys_prop),  intent(in)  :: pr
      character(len=*), intent(in)  :: fluid
      logical,          intent(in)  :: pres_tn, temp_tn, cs_tn

      call cmnlog_s(fmt_loc, 'min(dens)   ', fluid, pr%dens_min)
      call cmnlog_s(fmt_loc, 'max(dens)   ', fluid, pr%dens_max)
      if (temp_tn) then
         call cmnlog_s(fmt_loc, 'min(temp)   ', fluid, pr%temp_min)
         call cmnlog_s(fmt_loc, 'max(temp)   ', fluid, pr%temp_max)
      endif
      if (pres_tn) then
