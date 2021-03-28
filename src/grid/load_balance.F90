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
!! \brief This module contains variables and initialization routines related to
!!     * load balancing,
!!     * excluding poorly performing threads,
!!     * logging measured cg-related costs.
!<

module load_balance

   use cg_cost,   only: cost_labels
   use constants, only: cbuff_len

   implicit none

   private
   public :: init_load_balance, &
        &    auto_balance, cost_mask, balance_thr, verbosity_log, verbosity_stdout, verbosity_nstep, &
        &    enable_exclusion, watch_ind, nstep_exclusion, reset_exclusion, exclusion_thr

   ! namelist parameters
   ! balance:
   logical,                  protected :: auto_balance      !< When .true. the use cg-associated costs in rebalance routine
   character(len=cbuff_len), protected :: cost_to_balance   !< One of [ cg_cost:cost_labels, "all", "none" ], default: "MHD", ToDo: enable selected subset.
   real,                     protected :: balance_thr       !< Minimum tolerated imbalance

   ! verbosity
   integer(kind=4),          protected :: verbosity_log     !< Enumerated from 0: none, summary, detailed
   integer(kind=4),          protected :: verbosity_stdout  !< Enumerated from 0: none, summary, detailed
   integer(kind=4),          protected :: verbosity_nstep   !< How often to inform about cg costs

   ! thread exclusion
   logical,                  protected :: enable_exclusion  !< When .true. then threads detected as underperforming will be excluded for load balancing
   character(len=cbuff_len), protected :: watch_cost        !< Which cg cost to watch? One of [ cg_cost:cost_labels, "all", "none" ], default: "MHD"
   integer(kind=4),          protected :: nstep_exclusion   !< How often to check for underperforming CPUs
   integer(kind=4),          protected :: reset_exclusion   !< How often to routinely reset excluded status
   real,                     protected :: exclusion_thr     !< Exclusion threshold

   namelist /BALANCE/ auto_balance, cost_to_balance, balance_thr, verbosity_log, verbosity_stdout, verbosity_nstep, &
        &             enable_exclusion, watch_cost, nstep_exclusion, reset_exclusion, exclusion_thr

   logical, dimension(lbound(cost_labels,1):ubound(cost_labels,1)) :: cost_mask  !< translated cost_to_balance
   integer(kind=4) :: watch_ind  !< which index to watch for exclusion

   enum, bind(C)
      enumerator :: V_NONE = 0, V_SUM, V_DETAILED
   end enum

contains

!>
!! \brief Initialization of parameters of load balance mechanics
!!
!! @b BALANCE
!! \n \n
!! <table border="+1">
!!   <tr><td> auto_balance     </td><td> .false.  </td><td> logical     </td><td> \copydoc load_balance::auto_balance     </td></tr>
!!   <tr><td> cost_to_balance  </td><td> "MHD"    </td><td> character() </td><td> \copydoc load_balance::cost_to_balance  </td></tr>
!!   <tr><td> balance_thr      </td><td> 1.       </td><td> real        </td><td> \copydoc load_balance::balance_thr      </td></tr>
!!   <tr><td> verbosity_log    </td><td> 1        </td><td> integer     </td><td> \copydoc load_balance:verbosity_log:    </td></tr>
!!   <tr><td> verbosity_stdout </td><td> 0        </td><td> integer     </td><td> \copydoc load_balance::verbosity_stdout </td></tr>
!!   <tr><td> verbosity_nstep  </td><td> 10       </td><td> integer     </td><td> \copydoc load_balance::verbosity_nstep  </td></tr>
!!   <tr><td> enable_exclusion </td><td> .false.  </td><td> logical     </td><td> \copydoc load_balance::enable_exclusion </td></tr>
!!   <tr><td> watch_cost       </td><td> "MHD"    </td><td> character() </td><td> \copydoc load_balance::watch_cost       </td></tr>
!!   <tr><td> nstep_exclusion  </td><td> 1        </td><td> integer     </td><td> \copydoc load_balance::nstep_exclusion  </td></tr>
!!   <tr><td> reset_exclusion  </td><td> huge(1)  </td><td> integer     </td><td> \copydoc load_balance::reset_exclusion  </td></tr>
!!   <tr><td> exclusion_thr    </td><td> 3.       </td><td> real        </td><td> \copydoc load_balance::exclusion_thr    </td></tr>
!! </table>
!! \n \n
!<
   subroutine init_load_balance

      use constants,  only: I_ONE, INVALID, cbuff_len
      use dataio_pub, only: nh      ! QA_WARN required for diff_nml
      use dataio_pub, only: msg, printinfo, warn
      use mpisetup,   only: cbuff, ibuff, lbuff, rbuff, master, slave, piernik_MPI_Bcast

      implicit none

      integer, parameter :: verbosity_nstep_default = 10
      real,    parameter :: intolerable_perf = 3.
      real,    parameter :: tolerable_imbalance = 1.  ! A bit above 1. may prevent throwing cg back and forth
      integer            :: ind

      ! No code_progress dependencies apart from PIERNIK_INIT_MPI (obvious for our use of namelist)

      ! Namelist defaults
      auto_balance     = .false.
      cost_to_balance  = "MHD"
      balance_thr      = tolerable_imbalance
      verbosity_log    = V_SUM
      verbosity_stdout = V_NONE
      verbosity_nstep  = verbosity_nstep_default
      enable_exclusion = .false.
      watch_cost       = "MHD"
      nstep_exclusion  = I_ONE
      reset_exclusion  = huge(I_ONE)
      exclusion_thr    = intolerable_perf

      if (master) then

         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=BALANCE)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=""
         read(unit=nh%lun, nml=BALANCE, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "BALANCE")
         read(nh%cmdl_nml,nml=BALANCE, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "BALANCE", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=BALANCE)
         close(nh%lun)
         call nh%compare_namelist()

         cbuff(1) = cost_to_balance
         cbuff(2) = watch_cost

         ibuff(1) = verbosity_log
         ibuff(2) = verbosity_stdout
         ibuff(3) = verbosity_nstep
         ibuff(4) = nstep_exclusion
         ibuff(5) = reset_exclusion

         lbuff(1) = auto_balance
         lbuff(2) = enable_exclusion

         rbuff(1) = balance_thr
         rbuff(2) = exclusion_thr

      endif

      call piernik_MPI_Bcast(cbuff, cbuff_len)
      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(lbuff)
      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         cost_to_balance  = cbuff(1)
         watch_cost       = cbuff(2)

         verbosity_log    = ibuff(1)
         verbosity_stdout = ibuff(2)
         verbosity_nstep  = ibuff(3)
         nstep_exclusion  = ibuff(4)
         reset_exclusion  = ibuff(5)

         auto_balance     = lbuff(1)
         enable_exclusion = lbuff(2)

         balance_thr      = rbuff(1)
         exclusion_thr    = rbuff(2)

      endif

      ! decode strings
      watch_ind = decode_cost(watch_cost)

      cost_mask(:) = .false.
      ind = decode_cost(cost_to_balance)
      if (ind == INVALID) then
         if (trim(cost_to_balance) == "all") cost_mask(:) = .true.
      else
         cost_mask(ind) = .true.
      endif

      !sanitizing
      if (exclusion_thr <= 1.) then
         if (master) call warn("[load_balance] exclusion_thr <= 1. (disabling)")
         exclusion_thr = huge(1.)
      endif

      if (master) then
         if (any(cost_mask) .and. auto_balance) call printinfo("[load_balance] Auto-balance enabled")
         if (watch_ind /= INVALID .and. enable_exclusion) then
            write(msg, '(a,f4.1,a)')"[load_balance] Thread exclusion enabled (threshold = ", exclusion_thr, ")"
            call printinfo(msg)
         endif
      endif

   contains

      function decode_cost(str) result(cost_ind)

         use cg_cost,   only: cost_labels
         use constants, only: cbuff_len, INVALID

         implicit none

         character(len=cbuff_len), intent(in) :: str  !< string to decode

         integer(kind=4) :: cost_ind  !< index in costs related to str
         integer(kind=4) :: i

         cost_ind = INVALID
         do i = lbound(cost_labels, 1), ubound(cost_labels, 1)
            if (trim(cost_labels(i)) == trim(str)) then
               cost_ind = i
               exit
            endif
         enddo

      end function decode_cost

   end subroutine init_load_balance

end module load_balance
