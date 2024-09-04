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

   use cg_cost_data, only: cost_labels
   use constants,    only: cbuff_len

   implicit none

   public  ! QA_WARN no secrets are kept here

   ! namelist parameters
   !   balance
   real,                     protected :: balance_cg        !< Use the cg-associated costs in rebalance routine as weights. Value <= 0. means no cost weighting (assume all cg of equal weight), use 1. for full cost weighting.
   real,                     protected :: balance_host      !< Use averaged cg MHD costs to account for differing host speed. Value <=0 disables thread speed estimate (assume all hosts equally fast), use 1. for fully speed-weighted rebalance.
   logical,                  protected :: balance_thread    !< If .true. then use balance_host for each thread separately (CPU affinity has to be ensured outside of Piernik).
   character(len=cbuff_len), protected :: cost_to_balance   !< One of [ cg_cost_data:cost_labels, "all", "none" ], default: "all", ToDo: enable selected subset.
   real,                     protected :: balance_levels    !< If 1. then prioritize intra-level balancing for multigrid

   !   verbosity
   integer(kind=4),          protected :: verbosity         !< Enumerated from 0: none, summary, detailed, elaborate
   integer(kind=4),          protected :: verbosity_nstep   !< How often to inform about cg costs
   integer(kind=4)                     :: umsg_verbosity    !< Verbosity level requested by user msg

   !   thread exclusion
   logical,                  protected :: enable_exclusion  !< When .true. then threads detected as underperforming will be excluded for load balancing
   character(len=cbuff_len), protected :: watch_cost        !< Which cg cost to watch for anomalously slow hosts? One of [ cg_cost_data:cost_labels, "none" ], default: "MHD"
   real,                     protected :: exclusion_thr     !< Exclusion threshold
   logical,                  protected :: flexible_balance  !< can alter some parameters in this module if uneven load balancing is detected
   real,                     protected :: imbalance_tol     !< don't warn if load balance fits the range [imbalance_tol, 1./imbalance_tol]
   integer(kind=4),          protected :: n_rebalance       !< rebalance every n steps
   real,                     protected :: oop_thr           !< Maximum allowed ratio of Out-of-Place grid pieces (according to current ordering scheme)

   namelist /BALANCE/ balance_cg, balance_host, balance_thread, cost_to_balance, balance_levels, &
        &             verbosity, verbosity_nstep, flexible_balance, imbalance_tol, &
        &             enable_exclusion, watch_cost, exclusion_thr, n_rebalance, oop_thr

   logical, dimension(lbound(cost_labels,1):ubound(cost_labels,1)) :: cost_mask  !< translated from cost_to_balance
   integer(kind=4) :: watch_ind  !< which index to watch for exclusion in case of anomalously slow hosts

   enum, bind(C)
      enumerator :: VB_NONE = 0, VB_SUMMARY, VB_HOST, VB_DETAILED, VB_ELABORATE  !< verbosity levels
   end enum

   integer(kind=4), parameter :: r_rebalance = 2  !< routine rebalance 2 steps after start or restart
   logical :: rebalance_asap  !< flag to suggest rebalancing

contains

!>
!! \brief Initialization of parameters of load balance mechanics
!!
!! @b BALANCE
!! \n \n
!! <table border="+1">
!!   <tr><td> balance_cg       </td><td> 0.       </td><td> real        </td><td> \copydoc load_balance::balance_cg       </td></tr>
!!   <tr><td> balance_host     </td><td> 0.       </td><td> real        </td><td> \copydoc load_balance::balance_host     </td></tr>
!!   <tr><td> balance_thread   </td><td> .false.  </td><td> logical     </td><td> \copydoc load_balance::balance_thread   </td></tr>
!!   <tr><td> cost_to_balance  </td><td> "all"    </td><td> character() </td><td> \copydoc load_balance::cost_to_balance  </td></tr>
!!   <tr><td> balance_levels   </td><td> 0.       </td><td> real        </td><td> \copydoc load_balance::balance_levels   </td></tr>
!!   <tr><td> verbosity        </td><td> 2        </td><td> integer     </td><td> \copydoc load_balance::verbosity        </td></tr>
!!   <tr><td> verbosity_nstep  </td><td> 0        </td><td> integer     </td><td> \copydoc load_balance::verbosity_nstep  </td></tr>
!!   <tr><td> enable_exclusion </td><td> .false.  </td><td> logical     </td><td> \copydoc load_balance::enable_exclusion </td></tr>
!!   <tr><td> watch_cost       </td><td> "MHD"    </td><td> character() </td><td> \copydoc load_balance::watch_cost       </td></tr>
!!   <tr><td> exclusion_thr    </td><td> 3.       </td><td> real        </td><td> \copydoc load_balance::exclusion_thr    </td></tr>
!!   <tr><td> flexible_balance </td><td> .false.  </td><td> logical     </td><td> \copydoc load_balance::flexible_balance </td></tr>
!!   <tr><td> imbalance_tol    </td><td> 0.8      </td><td> real        </td><td> \copydoc load_balance::imbalance_tol    </td></tr>
!!   <tr><td> n_rebalance      </td><td> huge     </td><td> integer     </td><td> \copydoc load_balance::n_rebalance      </td></tr>
!!   <tr><td> oop_thr          </td><td> 0.05     </td><td> real        </td><td> \copydoc load_balance::oop_thr          </td></tr>
!! </table>
!! \n \n
!<
   subroutine init_load_balance

      use constants,  only: INVALID, cbuff_len, I_ONE
      use dataio_pub, only: nh      ! QA_WARN required for diff_nml
      use dataio_pub, only: msg, printinfo, warn
      use mpisetup,   only: cbuff, ibuff, lbuff, rbuff, master, slave, piernik_MPI_Bcast

      implicit none

      real,    parameter :: intolerable_perf = 3., insane_factor = 2.
      integer            :: ind

      ! No code_progress dependencies apart from PIERNIK_INIT_MPI (obvious for our use of namelist)

      ! Namelist defaults
      balance_cg       = 0.
      balance_host     = 0.
      balance_levels   = 0.
      balance_thread   = .false.
      cost_to_balance  = "all"
      verbosity        = VB_HOST
      verbosity_nstep  = 0
      enable_exclusion = .false.
      watch_cost       = "MHD"
      exclusion_thr    = intolerable_perf
      flexible_balance = .false.
      imbalance_tol    = 0.8
      n_rebalance      = huge(I_ONE)
      oop_thr          = 0.05

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

         ibuff(1) = verbosity
         ibuff(2) = verbosity_nstep
         ibuff(3) = n_rebalance

         lbuff(1) = balance_thread
         lbuff(2) = enable_exclusion
         lbuff(3) = flexible_balance

         rbuff(1) = exclusion_thr
         rbuff(2) = balance_cg
         rbuff(3) = balance_host
         rbuff(4) = balance_levels
         rbuff(5) = imbalance_tol
         rbuff(6) = oop_thr

      endif

      call piernik_MPI_Bcast(cbuff, cbuff_len)
      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(lbuff)
      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         cost_to_balance  = cbuff(1)
         watch_cost       = cbuff(2)

         verbosity        = ibuff(1)
         verbosity_nstep  = ibuff(2)
         n_rebalance      = ibuff(3)

         balance_thread   = lbuff(1)
         enable_exclusion = lbuff(2)
         flexible_balance = lbuff(3)

         exclusion_thr    = rbuff(1)
         balance_cg       = rbuff(2)
         balance_host     = rbuff(3)
         balance_levels   = rbuff(4)
         imbalance_tol    = rbuff(5)
         oop_thr          = rbuff(6)

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

      if (verbosity_nstep <= 0) then
         if (master) call warn("[load_balance] verbosity_nstep <= 0 (disabling)")
         verbosity_nstep = huge(1_4)
      endif

      if (balance_host > 1.) then
         if (balance_host > insane_factor) then
            if (master) call warn("[load_balance] balance_host >> 1. (reducing to 1.)")
            balance_host = 1.
         else
            if (master) call warn("[load_balance] balance_host > 1. may lead to unexpected behavior of load balancer")
         endif
      endif

      if (balance_cg > 1.) then
         if (balance_cg > insane_factor) then
            if (master) call warn("[load_balance] balance_cg >> 1. (reducing to 1.)")
            balance_cg = 1.
         else
            if (master) call warn("[load_balance] balance_cg > 1. may lead to unexpected behavior of load balancer")
         endif
      endif

      if (balance_levels > 1.) then
         if (balance_levels > insane_factor) then
            if (master) call warn("[load_balance] balance_levels >> 1. (reducing to 1.)")
            balance_levels = 1.
         else
            if (master) call warn("[load_balance] balance_levels > 1. may lead to unexpected behavior of load balancer")
         endif
      endif

      if (master) then
         if (any(cost_mask) .and. (balance_cg > 0. .or. balance_host > 0.)) then
            write(msg, '(a)')"[load_balance] Auto-balance enabled: "
            if (balance_host > 0.) write(msg(len_trim(msg)+1:), '(a,f5.2,a)') &
                 " balance_host = ", balance_host, " (" // trim(merge("thread", "host  ", balance_thread)) // "-based)"
            if (balance_cg > 0.) write(msg(len_trim(msg)+1:), '(a,f5.2)') " balance_cg = ", balance_cg
            call printinfo(msg)
         endif
         if (watch_ind /= INVALID .and. enable_exclusion) then
            write(msg, '(a,f4.1,a)')"[load_balance] Thread exclusion enabled (threshold = ", exclusion_thr, ")"
            call printinfo(msg)
         endif
      endif

      umsg_verbosity = VB_NONE
      if (imbalance_tol > 1.) imbalance_tol = 1. / imbalance_tol  ! normalize to [0. : 1.] range

      if (flexible_balance) then
         balance_thread = .true.
         if (balance_host > 0.) call warn("[load_balance] flexible_balance enforced balance_host = 1.")
         balance_host = 1.
#ifdef NBODY
         if (balance_cg <= epsilon(1.)) then
            ! Setting balance_cg to a non-0. value may help to minimize walltime for unbalanced particle setups.
            ! * The value 1.0 will try to equalize cg costs across all processes but since the particle operations don't
            !   overlap well with other computations, the wall time is not guaranteed to improve. Perhaps it may even
            !   degrade.
            ! * The value 0.5, that was usually close to the observable optimum in tests and may be treates as
            ! a conservative, safe choice.
            ! By running the code with balance_cg set to 1.0, 0.5 and 0.1 you should be able to estimate where to look
            ! for the minimum. It seems that 0.75 may beat the 0.5 and 1.0 by few percent.
            balance_cg = 0.5
            write(msg, '(a,f6.3,a)')"[load_balance] flexible_balance enforced balance_cg =", balance_cg, " but it is advised to check also other values up to 1.0 and estimate is the optimal value."
            if (master) call warn(msg)
            ! OPT: Generally, the problem with load balance of particles in Piernik is hard, at least for current
            ! implmentation. Possible ways to consider:
            ! * Rearrange particle computations in a way that will allow for overlap with other costly operations.
            !   Personally, I don't have idea to what extent it is possible if at all.
            ! * Remove the concept of ghost particles. This should significantly reduce the amount of work associated
            !   with particles and thus will reduce the imbalance. For proper communication of the TSC map, one would
            !   need to write modified boundary exchange routines that would compute sums on overlapping boundary
            !   regions. Special care is needed on fine-coarse boundaries to correctly account for projecting TSC onto
            !   grid of different cell size.
            ! * Write a dedicated load balancing and offload whole "ghost cg" with particles to particle-deficient
            !   threads and communicate the density and potential fields. This may work nicely when the number of
            !   particles on a cg exceeds the number of its cells. For most dense particle clumps consider offloading
            !   only subsets of particles. The ratios or thresholds at which the offloading should occur has to be
            !   determined, but it may require relatively small changes to the code.
            ! * Exploit shared memory benefits and maintain the list of particles as a property of a host (set of
            !   threads) rather than property of a cg. Requires big chamges in the codebase. If MPI-3 would be used
            !   as a shared memory enabler, then it should be possible to perform the migration to shared  memory use
            !   in parts and start from the particles as most benefit is expected here. It is also possible to perform
            !   a conditional migration and enable the shared memory only when the code runs on a single host.
            !   This will be still quite beneficial for medium-sized runs, that may fit single machine.
            !   After validating single-machine shared memory code it may be relatively easy to go for general case as
            !   the inter-host communication would be analogous to what we have in Piernik now with the communication
            !   between threads. Using shared memory should be also beneficial for other parts of the code:
            !   * a lot of boundary exchanges could be performed as a direct memory copy,
            !   * some memory copies related to global reductions may be even avoided (e.g. in the multipole algorithm),
            !   * it is possible to achieve better cache utilisation by parallel processing of a single cg, which may
            !     also result in less transfers between CPU and RAM memory.
         endif
#endif /* NBODY */
      endif

      rebalance_asap = .false.

   contains

      function decode_cost(str) result(cost_ind)

         use cg_cost_data, only: cost_labels
         use constants,    only: cbuff_len, INVALID

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
