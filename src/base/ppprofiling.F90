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
!! \brief Module for obtaining precise parallel profiling data
!!
!! \details MPI_Wtime-based event log that will be useful for parallel profiling
!! of the Piernik code.
!!
!! General purpose profilers don't know what is important in Piernik and tend to
!! provide a lot of detailed information about irrelevant routines.
!! This sometimes adds too much overhead from instrumenting the code which can be
!! misleading. Here, the developer can choose which parts of the code need closer
!! look.
!!
!! Collecting the events should be as cheap as possible, thus we don't construct
!! any trees inside Piernik. All the collected data is meant to be used in
!! postprocessing. The developer can set up several logs, each focused on
!! different aspect of the code and not interfering with other logs.
!!
!! The profiling is turned off by default. One can turn it on by:
!! * enabling in problem.par or in CLI
!!     ./piernik -n "&PROFILING use_profiling = T/"
!!   to get the whole run profiled.
!! * enabling via user_message_file
!!      echo ppp > msg
!!   to enable profiling for just one step, or
!!      echo ppp 10 > msg
!!   to enable profiling for 10 steps.
!! Interacting with PPP via msg file will override the choice for use_profiling
!! namelist parameter and can be used to stop profiling and flush the data to the
!! output file and continue the simulation without interrupting. The profiling may
!! be later turned on when needed. For example:
!!    ./piernik -n "&PROFILING use_profiling = T/" &\
!!        sleep 15 && echo ppp > msg &&\
!!        sleep 10 && echo ppp 2 > msg &&\
!!        sleep 1 && echo stop > msg
!! may give profiling data containing "init_piernik", first few steps, last step and
!! "finalize" (exact list depends on how long the steps are executed in the actual run).
!!
!! To exclude some categories of events one can use watch_* parameters.
!! By default MPI, AMR, debug, single-cg and auxiliary categories are excluded.
!! If an event is tagged with multiple categories, it will be excluded if any
!! of them is excluded.
!!
!! Perhaps it is possible to get much smaller output file (currently it is plain ASCII)
!! by implementing HDF5 output. I don't think it is critical as PPP profiling is
!! a debugging-type activity, not a production thing. It would require upgrading
!! python/ppp_plot.py script.
!<

module ppp

   use ppp_eventlist, only: eventlist

   implicit none

   private
   public :: init_profiling, cleanup_profiling, update_profiling, ppp_main, umsg_request

   ! namelist parameters
   logical :: watch_io         !< watch timers related to I/O
   logical :: watch_multigrid  !< watch timers related to multigrid
   logical :: watch_gravity    !< watch timers related to gravity
   logical :: watch_cr         !< watch timers related to cosmic rays
   logical :: watch_particles  !< watch timers related to particles
   logical :: watch_MPI        !< watch timers related to communication
   logical :: watch_AMR        !< watch timers related to refinements
   logical :: watch_cg         !< watch timers related to single-cg operations
   logical :: watch_magnetic   !< watch timers related to magnetic field
   logical :: watch_problem    !< watch timers related to problem
   logical :: watch_debug      !< watch timers related to debugging
   logical :: watch_aux        !< watch auxiliary timers
   logical :: xxl              !< allow for significantly bigger number of collected events (at your own risk)

   integer, save :: umsg_request = 0  !< turn on profiling for next umsg_request steps (read from msg file)

   type(eventlist) :: ppp_main  ! main eventlist

contains

!>
!! \brief initialize profiling output according to parameters from namelist PROFILING
!!
!! \n \n
!! @b PROFILING
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>use_profiling  </td><td>.false.  </td><td>logical value  </td><td>\copydoc ppp::use_profiling   </td></tr>
!! <tr><td>watch_io       </td><td>.true.   </td><td>logical value  </td><td>\copydoc ppp::watch_io        </td></tr>
!! <tr><td>watch_multigrid</td><td>.true.   </td><td>logical value  </td><td>\copydoc ppp::watch_multigrid </td></tr>
!! <tr><td>watch_gravity  </td><td>.true.   </td><td>logical value  </td><td>\copydoc ppp::watch_gravity   </td></tr>
!! <tr><td>watch_cr       </td><td>.true.   </td><td>logical value  </td><td>\copydoc ppp::watch_cr        </td></tr>
!! <tr><td>watch_particles</td><td>.true.   </td><td>logical value  </td><td>\copydoc ppp::watch_particles </td></tr>
!! <tr><td>watch_MPI      </td><td>.false.  </td><td>logical value  </td><td>\copydoc ppp::watch_MPI       </td></tr>
!! <tr><td>watch_AMR      </td><td>.false.  </td><td>logical value  </td><td>\copydoc ppp::watch_AMR       </td></tr>
!! <tr><td>watch_cg       </td><td>.false.  </td><td>logical value  </td><td>\copydoc ppp::watch_cg        </td></tr>
!! <tr><td>watch_magnetic </td><td>.true.   </td><td>logical value  </td><td>\copydoc ppp::watch_magnetic  </td></tr>
!! <tr><td>watch_problem  </td><td>.true.   </td><td>logical value  </td><td>\copydoc ppp::watch_problem   </td></tr>
!! <tr><td>watch_debug    </td><td>.false.  </td><td>logical value  </td><td>\copydoc ppp::watch_debug     </td></tr>
!! <tr><td>watch_aux      </td><td>.false.  </td><td>logical value  </td><td>\copydoc ppp::watch_aux       </td></tr>
!! <tr><td>xxl            </td><td>.false.  </td><td>logical value  </td><td>\copydoc ppp::xxl             </td></tr>
!! </table>
!! \n \n
!<

   subroutine init_profiling

      use constants,     only: PPP_IO, PPP_MG, PPP_GRAV, PPP_CR, PPP_PART, PPP_MPI, &
           &                   PPP_AMR, PPP_CG, PPP_MAG, PPP_PROB, PPP_DEBUG, PPP_AUX
      use dataio_pub,    only: nh, log_wr, problem_name, run_id, nrestart
      use mpisetup,      only: lbuff, master, slave, piernik_MPI_Bcast
      use ppp_eventlist, only: use_profiling, disable_mask, profile_file

      implicit none

      namelist /PROFILING/ use_profiling, xxl, &
           &               watch_io, watch_multigrid, watch_gravity, watch_cr, &
           &               watch_particles, watch_MPI, watch_AMR, watch_cg, &
           &               watch_magnetic, watch_problem, watch_debug, watch_aux

      use_profiling   = .false.
      watch_io        = .true.
      watch_multigrid = .true.
      watch_gravity   = .true.
      watch_cr        = .true.
      watch_particles = .true.
      watch_MPI       = .false.
      watch_AMR       = .false.
      watch_cg        = .false.
      watch_magnetic  = .true.
      watch_problem   = .true.
      watch_debug     = .false.
      watch_aux       = .false.
      xxl             = .false.

      if (master) then

         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=PROFILING)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=""
         read(unit=nh%lun, nml=PROFILING, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "PROFILING")
         read(nh%cmdl_nml,nml=PROFILING, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "PROFILING", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=PROFILING)
         close(nh%lun)
         call nh%compare_namelist()

         lbuff(1)  = use_profiling
         lbuff(3)  = watch_io
         lbuff(4)  = watch_multigrid
         lbuff(5)  = watch_gravity
         lbuff(6)  = watch_cr
         lbuff(7)  = watch_particles
         lbuff(8)  = watch_MPI
         lbuff(9)  = watch_AMR
         lbuff(10) = watch_cg
         lbuff(11) = watch_magnetic
         lbuff(12) = watch_problem
         lbuff(13) = watch_debug
         lbuff(14) = watch_aux
         lbuff(15) = xxl

      endif

      call piernik_MPI_Bcast(lbuff)

      if (slave) then

         use_profiling   = lbuff(1)
         watch_io        = lbuff(3)
         watch_multigrid = lbuff(4)
         watch_gravity   = lbuff(5)
         watch_cr        = lbuff(6)
         watch_particles = lbuff(7)
         watch_MPI       = lbuff(8)
         watch_AMR       = lbuff(9)
         watch_cg        = lbuff(10)
         watch_magnetic  = lbuff(11)
         watch_problem   = lbuff(12)
         watch_debug     = lbuff(13)
         watch_aux       = lbuff(14)
         xxl             = lbuff(15)

      endif

      disable_mask = 0
      if (.not. watch_io       ) disable_mask = disable_mask + PPP_IO
      if (.not. watch_multigrid) disable_mask = disable_mask + PPP_MG
      if (.not. watch_gravity  ) disable_mask = disable_mask + PPP_GRAV
      if (.not. watch_cr       ) disable_mask = disable_mask + PPP_CR
      if (.not. watch_particles) disable_mask = disable_mask + PPP_PART
      if (.not. watch_MPI      ) disable_mask = disable_mask + PPP_MPI
      if (.not. watch_AMR      ) disable_mask = disable_mask + PPP_AMR
      if (.not. watch_cg       ) disable_mask = disable_mask + PPP_CG
      if (.not. watch_magnetic ) disable_mask = disable_mask + PPP_MAG
      if (.not. watch_problem  ) disable_mask = disable_mask + PPP_PROB
      if (.not. watch_debug    ) disable_mask = disable_mask + PPP_DEBUG
      if (.not. watch_aux      ) disable_mask = disable_mask + PPP_AUX

      write(profile_file, '(6a,i3.3,a)') trim(log_wr), '/', trim(problem_name), '_', trim(run_id), '_', nrestart, '.ppprofile.ascii'

      call ppp_main%init("main", xxl)

   end subroutine init_profiling

!> \brief Turn on and off profiling upon request

   subroutine update_profiling

      use dataio_pub,    only: printinfo
      use global,        only: nstep
      use mpisetup,      only: master
      use ppp_eventlist, only: use_profiling

      implicit none

      integer, save :: turn_off = huge(1)

      if (turn_off <= nstep) then  ! turn off ppp_main profiling
         call ppp_main%publish
         if (master) call printinfo("[ppprofiling:update_profiling] Stop PPP")
         use_profiling = .false.
         turn_off = huge(1)
      endif

      if (umsg_request < 1) return

      ! this allows for overriding previous umsg_request
      ! or to disable profiling enabled in problem.par and flush the data
      turn_off = nstep + umsg_request
      umsg_request = 0

      if (.not. use_profiling) then  ! turn on ppp_main profiling
         use_profiling = .true.
         call ppp_main%init("main", xxl)
         if (master) call printinfo("[ppprofiling:update_profiling] Start PPP")
      endif

   end subroutine update_profiling

!> \brief close the profile file

   subroutine cleanup_profiling

      use dataio_pub,    only: close_txt_file
      use mpisetup,      only: master
      use ppp_eventlist, only: profile_file, profile_lun

      implicit none

      if (master) call close_txt_file(profile_file, profile_lun)

   end subroutine cleanup_profiling

end module ppp
