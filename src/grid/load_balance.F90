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
   public :: init_load_balance, print_costs

   ! namelist parameters
   !   balance
   logical,                  protected :: auto_balance      !< When .true. the use cg-associated costs in rebalance routine
   character(len=cbuff_len), protected :: cost_to_balance   !< One of [ cg_cost:cost_labels, "all", "none" ], default: "MHD", ToDo: enable selected subset.
   real,                     protected :: balance_thr       !< Minimum tolerated imbalance
   real,                     protected :: avg_factor        !< Exponential moving average factor

   !   verbosity
   integer(kind=4),          protected :: verbosity         !< Enumerated from 0: none, summary, detailed, elaborate
   integer(kind=4),          protected :: verbosity_nstep   !< How often to inform about cg costs

   !   thread exclusion
   logical,                  protected :: enable_exclusion  !< When .true. then threads detected as underperforming will be excluded for load balancing
   character(len=cbuff_len), protected :: watch_cost        !< Which cg cost to watch? One of [ cg_cost:cost_labels, "all", "none" ], default: "MHD"
   real,                     protected :: exclusion_thr     !< Exclusion threshold

   namelist /BALANCE/ auto_balance, cost_to_balance, balance_thr, avg_factor, verbosity, verbosity_nstep, &
        &             enable_exclusion, watch_cost, exclusion_thr

   logical, dimension(lbound(cost_labels,1):ubound(cost_labels,1)) :: cost_mask  !< translated from cost_to_balance
   integer(kind=4) :: watch_ind  !< which index to watch for exclusion

   enum, bind(C)
      enumerator :: V_NONE = 0, V_SUMMARY, V_HOST, V_DETAILED, V_ELABORATE  !< verbosity levels
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
!!   <tr><td> avg_factor       </td><td> 0.2      </td><td> real        </td><td> \copydoc load_balance::avg_factor       </td></tr>
!!   <tr><td> verbosity        </td><td> 2        </td><td> integer     </td><td> \copydoc load_balance::verbosity        </td></tr>
!!   <tr><td> verbosity_nstep  </td><td> 20       </td><td> integer     </td><td> \copydoc load_balance::verbosity_nstep  </td></tr>
!!   <tr><td> enable_exclusion </td><td> .false.  </td><td> logical     </td><td> \copydoc load_balance::enable_exclusion </td></tr>
!!   <tr><td> watch_cost       </td><td> "MHD"    </td><td> character() </td><td> \copydoc load_balance::watch_cost       </td></tr>
!!   <tr><td> exclusion_thr    </td><td> 3.       </td><td> real        </td><td> \copydoc load_balance::exclusion_thr    </td></tr>
!! </table>
!! \n \n
!<
   subroutine init_load_balance

      use constants,  only: INVALID, cbuff_len
      use dataio_pub, only: nh      ! QA_WARN required for diff_nml
      use dataio_pub, only: msg, printinfo, warn
      use mpisetup,   only: cbuff, ibuff, lbuff, rbuff, master, slave, piernik_MPI_Bcast

      implicit none

      integer, parameter :: verbosity_nstep_default = 20
      real,    parameter :: intolerable_perf = 3.
      real,    parameter :: tolerable_imbalance = 1.  ! A bit above 1. may prevent throwing cg back and forth
      real,    parameter :: default_ema = 0.2  ! Exponential moving average factor (1.0 => no tail, 0.0 => no update)
      integer            :: ind

      ! No code_progress dependencies apart from PIERNIK_INIT_MPI (obvious for our use of namelist)

      ! Namelist defaults
      auto_balance     = .false.
      cost_to_balance  = "MHD"
      balance_thr      = tolerable_imbalance
      avg_factor       = default_ema
      verbosity        = V_HOST
      verbosity_nstep  = verbosity_nstep_default
      enable_exclusion = .false.
      watch_cost       = "MHD"
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

         ibuff(1) = verbosity
         ibuff(2) = verbosity_nstep

         lbuff(1) = auto_balance
         lbuff(2) = enable_exclusion

         rbuff(1) = balance_thr
         rbuff(2) = avg_factor
         rbuff(3) = exclusion_thr

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

         auto_balance     = lbuff(1)
         enable_exclusion = lbuff(2)

         balance_thr      = rbuff(1)
         avg_factor       = rbuff(2)
         exclusion_thr    = rbuff(3)

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

      if (abs(avg_factor) > 1.) then
         if (master) call warn("[load_balance] |avg_factor| > 1. (disabling)")
         avg_factor = 0.
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

!> \brief Report the measured costs

   subroutine print_costs(cglist)

      use cg_cost,       only: cost_labels
      use cg_cost_stats, only: cg_stats_t, stat_labels, I_MAX
      use cg_list,       only: cg_list_t, cg_list_element
      use constants,     only: I_ONE, base_level_id, PPP_AMR
      use dataio_pub,    only: msg, warn
      use global,        only: nstep
      use mpisetup,      only: master, FIRST, LAST, err_mpi
      use MPIF,          only: MPI_COMM_WORLD, MPI_DOUBLE_PRECISION, MPI_Wtime, MPI_Gather
      use ppp,           only: ppp_main

      implicit none

      class(cg_list_t), intent(in) :: cglist

      type(cg_list_element), pointer :: cgl
      type(cg_stats_t) :: leaves_stats, all_stats
      integer(kind=4) :: i
      real :: maxv, mul
      character(len=I_ONE) :: prefix
      real, save :: prev_time = -huge(1.)
      real, allocatable, dimension(:, :, :) :: all_proc_stats
      real, allocatable, dimension(:, :) :: send_stats
      integer, parameter :: N_STATS = size(stat_labels) + I_ONE
      character(len=*), parameter :: lbpc_label = "load_balance:print_costs"

      call ppp_main%start(lbpc_label, PPP_AMR)

      ! gather mean, standard deviation and extrema for the costs
      call leaves_stats%reset
      call all_stats%reset
      cgl => cglist%first
      do while (associated(cgl))
         if (cgl%cg%l%id >= base_level_id) call leaves_stats%add(cgl%cg%costs)
         call all_stats%add(cgl%cg%costs)
         cgl => cgl%nxt
      enddo

      if (master) then
         allocate(all_proc_stats(N_STATS, size(cost_labels), FIRST:LAST))
      else
         allocate(all_proc_stats(0, 0, 0))
      endif

      allocate(send_stats(N_STATS, size(cost_labels)))
      do i = lbound(cost_labels, 1), ubound(cost_labels, 1)
         send_stats(:, i - lbound(cost_labels, 1) + lbound(send_stats, 2)) = [ leaves_stats%get(i), all_stats%get_sum() ]
      enddo
      call MPI_Gather(send_stats,     size(send_stats, kind=4), MPI_DOUBLE_PRECISION, &
           &          all_proc_stats, size(send_stats, kind=4), MPI_DOUBLE_PRECISION, &
           &          FIRST, MPI_COMM_WORLD, err_mpi)

      if (master) then

         if (nstep >= 1) call update_costs

         if ((verbosity > V_NONE) .and. ((nstep <= 0 .and. verbosity_nstep <= 1) .or. (nstep > 0 .and. mod(nstep, verbosity_nstep) == 0))) then

            ! Choose between s and ms.
            maxv = maxval(all_proc_stats(:, I_MAX - lbound(stat_labels) + I_ONE, :))
            if (maxv > 0.999) then  ! let's prevent rounding up by write to 1000 ms
               mul = 1.
               prefix = " "
            else
               mul = 1e3
               prefix = "m"
            endif

            if (maxv <= 0.) then
               if (prev_time >= 0.) then
                  write(msg, '(a,f11.6,a)')"No cg cost data found anywhere. Walltime:", MPI_Wtime() - prev_time, " s"
                  call warn(msg)
               endif
            else
               select case (verbosity)
                  case (V_ELABORATE:)
                     call log_elaborate
                  case (V_DETAILED)
                     call log_detailed
                  case (V_HOST)
                  case (V_SUMMARY)
               end select
               if (verbosity >= V_HOST) call log_host
               if (verbosity >= V_SUMMARY) call log_summary
               if ((auto_balance .or. enable_exclusion) .and. (nstep >= 1)) call log_speed
            endif

         endif
      endif

      deallocate(all_proc_stats, send_stats)

      prev_time = MPI_Wtime()

      call ppp_main%stop(lbpc_label, PPP_AMR)

   contains

      subroutine update_costs

         use cg_cost,       only: cost_labels, I_MHD  ! Hardcoded MHD cost as the most reliable measure. ToDo: give some choice.
         use cg_cost_stats, only: I_AVG
         use procnames,     only: pnames

         implicit none

         logical, save :: firstcall = .true.
         integer :: p

         do p = FIRST, LAST
            associate (cur_speed => all_proc_stats(I_AVG, I_MHD - lbound(cost_labels, 1) + I_ONE, p), &
               &       pspeed => pnames%speed(p))
               if (firstcall .or. cur_speed <= 0.) then  ! first call or depleted thread
                  call pspeed%add(max(0., cur_speed), avg_factor)
               else
                  if (pnames%speed(p)%avg <= 0.) then  ! repopulated thread (e.g. after "unexclude" msg
                     call pspeed%add(cur_speed, avg_factor)
                  else
                     call pspeed%add(cur_speed)
                  endif
               endif
            end associate
         enddo

         if (firstcall) then
            call pnames%calc_hostspeed(avg_factor)
         else
            call pnames%calc_hostspeed
         endif

         call pnames%mark_for_exclusion(exclusion_thr)

         firstcall = .false.

      end subroutine update_costs

      subroutine log_elaborate

         use cg_cost_stats, only: I_SUM2
         use dataio_pub,    only: printinfo
         use procnames,     only: pnames

         implicit none

         integer :: j, p
         real, dimension(lbound(stat_labels,1):ubound(stat_labels,1)) :: stat

         do i = lbound(cost_labels, 1), ubound(cost_labels, 1)
            if (any(all_proc_stats(I_MAX - lbound(stat_labels) + I_ONE, i - lbound(cost_labels, 1) + I_ONE, :) > 0.)) then
               do p = FIRST, LAST  ! lbound(all_proc_stats, 3), ubound(all_proc_stats, 3)
                  stat = all_proc_stats(:N_STATS-I_ONE, i - lbound(cost_labels, 1) + I_ONE, p)
                  if (stat(I_MAX) > 0.) then
                     write(msg, '(2a,i5,3a)')"@", pnames%procnames(p)(:pnames%maxnamelen), p, " Cost('", cost_labels(i), "')"
                     do j = lbound(stat_labels, 1), ubound(stat_labels, 1)
                        if (j /= I_SUM2) &
                             write(msg(len_trim(msg)+1:), '(3a,f10.3,3a)') merge(": ", ", ", j == lbound(stat_labels, 1)), &
                             &                                             trim(stat_labels(j)), "= ", mul*stat(j), " ", trim(prefix), "s "
                     enddo
                     call printinfo(msg)
                  endif
               enddo
            endif
         enddo

      end subroutine log_elaborate

      subroutine log_detailed

         use cg_cost_stats, only: I_AVG, I_SIGMA
         use dataio_pub,    only: printinfo
         use procnames,     only: pnames

         implicit none

         integer :: p

         write(msg, '(2a,a5,a)')" host", repeat(" ", pnames%maxnamelen - 4), "rank", " :"
         do i = lbound(cost_labels, 1), ubound(cost_labels, 1)
            if (any(all_proc_stats(I_AVG, i - lbound(cost_labels, 1) + I_ONE, :) > 0.)) then
               write(msg(len_trim(msg)+1:), '(a24)') "avg(" // trim(cost_labels(i)) // ") ± σ"
            endif
         enddo
         write(msg(len_trim(msg)+1:), '(a)') " (per cg)"
         call printinfo(msg)

         do p = FIRST, LAST
            write(msg, '(2a,i5,a)')"@", pnames%procnames(p)(:pnames%maxnamelen), p, " :"
            do i = lbound(cost_labels, 1), ubound(cost_labels, 1)
               if (any(all_proc_stats(I_AVG, i - lbound(cost_labels, 1) + I_ONE, :) > 0.)) then
                  if (all_proc_stats(I_AVG,   i - lbound(cost_labels, 1) + I_ONE, p) > 0.) then
                     write(msg(len_trim(msg)+1:), '(f10.3,3a,f6.1,a)') &
                          &  mul*all_proc_stats(I_AVG,   i - lbound(cost_labels, 1) + I_ONE, p),  " ", prefix, "s ±", &
                          & 100.*all_proc_stats(I_SIGMA, i - lbound(cost_labels, 1) + I_ONE, p) / &
                          &      all_proc_stats(I_AVG,   i - lbound(cost_labels, 1) + I_ONE, p), "%"
                  else
                     write(msg(len_trim(msg)+1:), '(f10.3,3a,a6,a)') &
                          &  mul*all_proc_stats(I_AVG,   i - lbound(cost_labels, 1) + I_ONE, p),  " ", prefix, "s  ", "N/A", "%"
                  endif
               endif
            enddo
            call printinfo(msg)
         enddo

      end subroutine log_detailed

      subroutine log_host

         use cg_cost,       only: cg_cost_data_t
         use cg_cost_stats, only: I_AVG, I_SUM, I_SUM2
         use dataio_pub,    only: printinfo
         use procnames,     only: pnames

         implicit none

         integer :: host, p1, i
         type(cg_cost_data_t) :: n, sum, sum2

         write(msg, '(2a,a5,a)')" host", repeat(" ", pnames%maxnamelen - 4), " ", " :"
         do i = lbound(cost_labels, 1), ubound(cost_labels, 1)
            if (any(all_proc_stats(I_AVG, i - lbound(cost_labels, 1) + I_ONE, :) > 0.)) then
               write(msg(len_trim(msg)+1:), '(a24)') "avg(" // trim(cost_labels(i)) // ") ± σ"
            endif
         enddo
         write(msg(len_trim(msg)+1:), '(a)') " (per cg, on whole host)"
         call printinfo(msg)

         do host = lbound(pnames%proc_on_node, 1), ubound(pnames%proc_on_node, 1)
            ! this is a bit awkward but we need to find moments over a host from already computed moments on threads
            n%wtime = 0.
            sum%wtime = 0.
            sum2%wtime = 0.
            associate (ph => pnames%proc_on_node(host))
               do p1 = lbound(ph%proc, 1), ubound(ph%proc, 1)
                  associate (p => ph%proc(p1))
                     where (all_proc_stats(I_AVG, :, p) > 0.)  ! reconstruct number of cg, total time and sum of squared durations
                        n%wtime(:) = n%wtime(:) + all_proc_stats(I_SUM, :, p) / all_proc_stats(I_AVG, :, p)
                        sum%wtime(:) = sum%wtime(:) + all_proc_stats(I_SUM, :, p)
                        sum2%wtime(:) = sum2%wtime(:) + all_proc_stats(I_SUM2, :, p)
                     endwhere
                  end associate
               enddo
               write(msg, '(3a)')"@", ph%nodename(:pnames%maxnamelen), "      :"
               do i = lbound(cost_labels, 1), ubound(cost_labels, 1)
                  if (any(all_proc_stats(I_AVG, i - lbound(cost_labels, 1) + I_ONE, :) > 0.)) then
                     if (sum%wtime(i) > 0.) then
                        write(msg(len_trim(msg)+1:), '(f10.3,3a,f6.1,a)') &
                             &  mul * sum%wtime(i) / n%wtime(i),  " ", prefix, "s ±", &
                             & 100. * sqrt(sum2%wtime(i) * n%wtime(i) / sum%wtime(i)**2 - 1.), "%"
                     else
                        write(msg(len_trim(msg)+1:), '(f10.3,3a,a6,a)') 0.,  " ", prefix, "s  ", "N/A", "%"
                     endif
                  endif
               enddo
               call printinfo(msg)
            end associate
         enddo

      end subroutine log_host

      subroutine log_summary

         use dataio_pub, only: printinfo
         use procnames,  only: pnames

         implicit none

         integer :: p, host, lines, per_line, l
         real :: dt_wall
         integer, parameter :: max_one_line = 8, max_per_line = 12

         if (prev_time >= 0.) then
            dt_wall = MPI_Wtime() - prev_time

            write(msg, '(a,f11.3,3a)') "All accumulated cg costs out of ", mul*dt_wall, " ", trim(prefix), "s"
            call printinfo(msg)

            do host = lbound(pnames%proc_on_node, 1), ubound(pnames%proc_on_node, 1)
               associate (ph => pnames%proc_on_node(host))

                  write(msg, '(3a)')"@", ph%nodename(:pnames%maxnamelen), " : "

                  if (size(ph%proc) <= max_one_line) then
                     do p = lbound(ph%proc, 1), ubound(ph%proc, 1)
                        write(msg(len_trim(msg)+1:), '(f11.3)') mul*all_proc_stats(N_STATS, I_ONE, ph%proc(p))
                     enddo
                     write(msg(len_trim(msg)+1:), '(a)') " | "
                     do p = lbound(ph%proc, 1), ubound(ph%proc, 1)
                        write(msg(len_trim(msg)+1:), '(f6.1,a)') all_proc_stats(N_STATS, I_ONE, ph%proc(p))/dt_wall * 100., "%"
                     enddo
                     call printinfo(msg)
                  else
                     lines = ceiling(real(size(ph%proc)) / max_per_line)
                     per_line = ceiling(real(size(ph%proc)) / lines)  ! or max_per_line for extremely heterogenous set of nodes
                     do l = 1, lines
                        if (l > I_ONE) write(msg, '(2a)') repeat(" ", pnames%maxnamelen + 1), " : "
                        do p =   lbound(ph%proc, 1) + per_line * (l - 1), &
                             min(lbound(ph%proc, 1) + per_line *  l - 1 , &
                             &   ubound(ph%proc, 1))
                           write(msg(len_trim(msg)+1:), '(f11.3)') mul*all_proc_stats(N_STATS, I_ONE, ph%proc(p))
                        enddo
                        call printinfo(msg)
                     enddo
                     do l = 1, lines
                        write(msg, '(2a)') repeat(" ", pnames%maxnamelen + 1), " | "
                        do p =   lbound(ph%proc, 1) + per_line * (l - 1), &
                             min(lbound(ph%proc, 1) + per_line *  l - 1 , &
                             &   ubound(ph%proc, 1))
                           write(msg(len_trim(msg)+1:), '(f10.1,a)') all_proc_stats(N_STATS, I_ONE, ph%proc(p))/dt_wall * 100., "%"
                        enddo
                        call printinfo(msg)
                     enddo
                  endif

               end associate
            enddo

         endif

      end subroutine log_summary

      subroutine log_speed

         use constants,  only: fmt_len
         use dataio_pub, only: printinfo, warn
         use procnames,  only: pnames

         implicit none

         integer :: host, dec, ln
         character(len=fmt_len) :: fmt, header
         real :: mx
         integer, parameter :: mpl = 16, maxex = 128

         mx = maxval(1./pnames%speed(:)%avg, mask=(pnames%speed(:)%avg > 0.))
         dec = 3
         if (mx > 0.) dec = max(0, min(3, 3 - int(floor(log10(mx)))))

         do host = lbound(pnames%proc_on_node, 1), ubound(pnames%proc_on_node, 1)
            associate (ph => pnames%proc_on_node(host))

               do ln = 0, int((size(ph%proc) - 1)/ mpl)
                  associate (pb =>     lbound(ph%proc, 1) +  ln    * mpl, &
                       &     pe => min(lbound(ph%proc, 1) + (ln+1) * mpl - 1, ubound(ph%proc, 1)))

                     write(fmt, *)"(3a,f6.", dec, ",a)"
                     write(header, fmt) "@", ph%nodename(:pnames%maxnamelen), " <MHD speed> = ", &
                          merge(1. / ph%speed%avg, 0., ph%speed%avg > 0.), " blk/s ["
                     write(fmt,  *) "(a,", pe - pb + 1, "f6.", dec, ",a)"
                     write(msg, fmt) merge(trim(header), repeat(" ", len_trim(header)), ln == 0), &
                          merge(1. / pnames%speed(ph%proc(pb:pe))%avg, 0., pnames%speed(ph%proc(pb:pe))%avg > 0.), &
                          merge(" ]", "  ", ln == int((size(ph%proc) - 1)/ mpl))
                     call printinfo(msg)

                  end associate
               enddo

            end associate
         enddo

         if (any(pnames%exclude)) then
            call warn("[load_balance] There are threads marked for exclusion:")
            do host = lbound(pnames%proc_on_node, 1), ubound(pnames%proc_on_node, 1)
               associate (ph => pnames%proc_on_node(host))
                  if (any(pnames%exclude(ph%proc))) then
                     header = "@" // ph%nodename(:pnames%maxnamelen)
                     if (size(ph%proc) <= maxex) then
                        write(fmt, *)"(a,", size(ph%proc), "a2,a)"
                        write(msg, fmt) trim(header) // ": [", merge("X", ".", pnames%exclude(ph%proc)), " ]"
                     else
                        write(msg, '(a,2(i5,a))') trim(header), count(pnames%exclude(ph%proc)), " out of ", size(ph%proc), " threads"
                     endif
                     call warn(msg)
                  endif
               end associate
            enddo
         endif

      end subroutine log_speed

   end subroutine print_costs

end module load_balance
