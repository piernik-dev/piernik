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

!> \brief This module contains routines related to load balancing,

module lb_helpers

   implicit none

   private
   public :: costs_maintenance

contains

   subroutine costs_maintenance

      use cg_list_global, only: all_cg

      implicit none

      call print_costs
      call all_cg%reset_costs

   end subroutine costs_maintenance

!> \brief Report the measured costs

   subroutine print_costs

      use cg_cost_data,   only: cost_labels
      use cg_cost_stats,  only: cg_stats_t, stat_labels, I_MAX
      use cg_list,        only: cg_list_element
      use cg_list_global, only: all_cg
      use constants,      only: I_ONE, base_level_id, PPP_AMR
      use dataio_pub,     only: msg, warn
      use global,         only: nstep
      use load_balance,   only: balance_cg, balance_host, enable_exclusion, exclusion_thr, &
           &                    verbosity, verbosity_nstep, umsg_verbosity, &
           &                    V_NONE, V_SUMMARY, V_HOST, V_DETAILED, V_ELABORATE
      use mpisetup,       only: master, FIRST, LAST, err_mpi
      use MPIF,           only: MPI_COMM_WORLD, MPI_DOUBLE_PRECISION, MPI_Wtime
      use MPIFUN,         only: MPI_Gather
      use ppp,            only: ppp_main

      implicit none

      type(cg_list_element), pointer :: cgl
      type(cg_stats_t) :: leaves_stats, all_stats
      integer(kind=4) :: i, v
      real :: maxv, mul
      character(len=I_ONE) :: prefix
      real, save :: prev_time = -huge(1.)
      real, allocatable, dimension(:, :, :) :: all_proc_stats
      real, allocatable, dimension(:, :) :: send_stats
      integer, parameter :: N_STATS = size(stat_labels) + I_ONE
      character(len=*), parameter :: lbpc_label = "load_balance:print_costs"
      integer(kind=4), parameter :: h_l = len("host")

      call ppp_main%start(lbpc_label, PPP_AMR)

      ! gather mean, standard deviation and extrema for the costs
      call leaves_stats%reset
      call all_stats%reset
      cgl => all_cg%first
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

         if (((verbosity > V_NONE) .and. ((nstep <= 0 .and. verbosity_nstep <= 1) .or. (nstep > 0 .and. mod(nstep, verbosity_nstep) == 0))) &
              .or. (umsg_verbosity > V_NONE)) then

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
               v = verbosity
               if (umsg_verbosity > V_NONE) then
                  v = min(V_ELABORATE, umsg_verbosity)
                  umsg_verbosity = V_NONE
                  ! Intentionally the reset occurs only when umsg_verbosity triggers printing (maxv has to be > 0.)
               endif
               select case (v)
                  case (V_ELABORATE:)
                     call log_elaborate
                  case (V_DETAILED)
                     call log_detailed
                  case (V_HOST)
                  case (V_SUMMARY)
               end select
               if (v >= V_HOST) call log_host
               if (v >= V_SUMMARY) call log_summary
               if ((balance_cg > 0. .or. balance_host > 0. .or. enable_exclusion) .and. (nstep >= 1)) call log_speed
            endif

         endif

         call speed_check

      endif

      deallocate(all_proc_stats, send_stats)

      prev_time = MPI_Wtime()

      call ppp_main%stop(lbpc_label, PPP_AMR)

   contains

      ! It seems that most reliable mark of host speed is the MHD timer because of most even load there.

      subroutine update_costs

         use cg_cost_data,  only: cost_labels
         use cg_cost_stats, only: I_AVG
         use load_balance,  only: watch_ind
         use procnames,     only: pnames

         implicit none

         logical, save :: firstcall = .true.
         integer :: p

         do p = FIRST, LAST
            associate (cur_speed => all_proc_stats(I_AVG, watch_ind - lbound(cost_labels, 1) + I_ONE, p), &
               &       pspeed => pnames%wtime(p))
               if (firstcall .or. cur_speed <= 0.) then  ! first call or depleted thread
                  pspeed = max(0., cur_speed)
               else
                  pspeed = cur_speed
               endif
            end associate
         enddo

         call pnames%calc_hostspeed
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
         use constants,     only: I_ZERO
         use dataio_pub,    only: printinfo
         use procnames,     only: pnames

         implicit none

         integer :: p

         write(msg, '(2a,a5,a)')" host", repeat(" ", max(I_ZERO, pnames%maxnamelen - h_l)), "rank", " :"
         do i = lbound(cost_labels, 1), ubound(cost_labels, 1)
            if (any(all_proc_stats(I_AVG, i - lbound(cost_labels, 1) + I_ONE, :) > 0.)) then
               write(msg(len_trim(msg)+1:), '(a24)') "avg(" // trim(cost_labels(i)) // ") ± σ"
            endif
         enddo
         write(msg(len_trim(msg)+1:), '(a)') " (per cg)"
         call printinfo(msg)

         do p = FIRST, LAST
            write(msg, '(2a,i5,a)')"@", pnames%procnames(p)(:max(h_l, pnames%maxnamelen)), p, " :"
            do i = lbound(cost_labels, 1), ubound(cost_labels, 1)
               if (any(all_proc_stats(I_AVG, i - lbound(cost_labels, 1) + I_ONE, :) > 0.)) then
                  if (all_proc_stats(I_AVG,   i - lbound(cost_labels, 1) + I_ONE, p) > 0.) then
                     write(msg(len_trim(msg)+1:), '(f10.3,3a,f6.1,a)') &
                          &  mul*all_proc_stats(I_AVG,   i - lbound(cost_labels, 1) + I_ONE, p),  " ", prefix, "s ±", &
                          & 100.*all_proc_stats(I_SIGMA, i - lbound(cost_labels, 1) + I_ONE, p) / &
                          &      all_proc_stats(I_AVG,   i - lbound(cost_labels, 1) + I_ONE, p), "%"
                  else
                     write(msg(len_trim(msg)+1:), '(a22)') "-"
                  endif
               endif
            enddo
            call printinfo(msg)
         enddo

      end subroutine log_detailed

      subroutine log_host

         use cg_cost_data,  only: cg_cost_data_t
         use cg_cost_stats, only: I_AVG, I_SUM, I_SUM2
         use constants,     only: I_ZERO
         use dataio_pub,    only: printinfo
         use procnames,     only: pnames

         implicit none

         integer :: host, p1, i
         type(cg_cost_data_t) :: n, sum, sum2

         write(msg, '(2a,a5,a)')" host", repeat(" ", max(I_ZERO, pnames%maxnamelen - h_l)), " ", " :"
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
               write(msg, '(3a)')"@", ph%nodename(:max(h_l, pnames%maxnamelen)), "      :"
               do i = lbound(cost_labels, 1), ubound(cost_labels, 1)
                  if (any(all_proc_stats(I_AVG, i - lbound(cost_labels, 1) + I_ONE, :) > 0.)) then
                     if (sum%wtime(i) > 0.) then
                        write(msg(len_trim(msg)+1:), '(f10.3,3a,f6.1,a)') &
                             &  mul * sum%wtime(i) / n%wtime(i),  " ", prefix, "s ±", &
                             & 100. * sqrt(sum2%wtime(i) * n%wtime(i) / sum%wtime(i)**2 - 1.), "%"
                     else
                        write(msg(len_trim(msg)+1:), '(a22)') "-"
                     endif
                  endif
               enddo
               call printinfo(msg)
            end associate
         enddo

      end subroutine log_host

      subroutine log_summary

         use dataio_pub, only: printinfo
         use mpisetup,   only: nproc
         use procnames,  only: pnames

         implicit none

         integer :: p, host, lines, per_line, l
         real :: dt_wall
         integer, parameter :: max_one_line = 8, max_per_line = 12
         integer, parameter :: s_ind = I_ONE  ! can be anything within size(cost_labels) because all values are the same

         if (prev_time >= 0.) then
            dt_wall = MPI_Wtime() - prev_time

            write(msg, '(a,f11.3,3a,f5.1,a)') "All accumulated cg costs out of ", mul*dt_wall, " ", trim(prefix), "s, ", 100*sum(all_proc_stats(N_STATS, I_ONE, :))/nproc/dt_wall, "% spent on cg (averaged globally)"
            call printinfo(msg)

            do host = lbound(pnames%proc_on_node, 1), ubound(pnames%proc_on_node, 1)
               associate (ph => pnames%proc_on_node(host))

                  write(msg, '(3a)')"@", ph%nodename(:pnames%maxnamelen), " : "

                  if (size(ph%proc) <= max_one_line) then
                     do p = lbound(ph%proc, 1), ubound(ph%proc, 1)
                        if (all_proc_stats(N_STATS, s_ind, ph%proc(p)) > 0.) then
                           write(msg(len_trim(msg)+1:), '(f11.3)') mul*all_proc_stats(N_STATS, s_ind, ph%proc(p))
                        else
                           write(msg(len_trim(msg)+1:), '(a11)') "-"
                        endif
                     enddo
                     write(msg(len_trim(msg)+1:), '(a)') " | "
                     do p = lbound(ph%proc, 1), ubound(ph%proc, 1)
                        if (all_proc_stats(N_STATS, s_ind, ph%proc(p)) > 0.) then
                           write(msg(len_trim(msg)+1:), '(f6.1,a)') all_proc_stats(N_STATS, s_ind, ph%proc(p))/dt_wall * 100., "%"
                        else
                           write(msg(len_trim(msg)+1:), '(a7)') "-"
                        endif
                     enddo
                     call printinfo(msg)
                  else
                     lines = ceiling(real(size(ph%proc)) / max_per_line)
                     per_line = ceiling(real(size(ph%proc)) / lines)  ! or max_per_line for extremely heterogeneous set of nodes
                     do l = 1, lines
                        if (l > I_ONE) write(msg, '(2a)') repeat(" ", pnames%maxnamelen + 1), " : "
                        do p =   lbound(ph%proc, 1) + per_line * (l - 1), &
                             min(lbound(ph%proc, 1) + per_line *  l - 1 , &
                             &   ubound(ph%proc, 1))
                           if (all_proc_stats(N_STATS, s_ind, ph%proc(p)) > 0.) then
                              write(msg(len_trim(msg)+1:), '(f11.3)') mul*all_proc_stats(N_STATS, s_ind, ph%proc(p))
                           else
                              write(msg(len_trim(msg)+1:), '(a11)') "-"
                           endif
                        enddo
                        call printinfo(msg)
                     enddo
                     do l = 1, lines
                        write(msg, '(2a)') repeat(" ", pnames%maxnamelen + 1), " | "
                        do p =   lbound(ph%proc, 1) + per_line * (l - 1), &
                             min(lbound(ph%proc, 1) + per_line *  l - 1 , &
                             &   ubound(ph%proc, 1))
                           if (all_proc_stats(N_STATS, s_ind, ph%proc(p)) > 0.) then
                              write(msg(len_trim(msg)+1:), '(f10.1,a)') all_proc_stats(N_STATS, s_ind, ph%proc(p))/dt_wall * 100., "%"
                           else
                              write(msg(len_trim(msg)+1:), '(a11)') "-"
                           endif
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

         mx = maxval(1./pnames%wtime(:), mask=(pnames%wtime(:) > 0.))
         dec = 3
         if (mx > 0.) dec = max(0, min(3, 3 - int(floor(log10(mx)))))

         do host = lbound(pnames%proc_on_node, 1), ubound(pnames%proc_on_node, 1)
            associate (ph => pnames%proc_on_node(host))

               do ln = 0, int((size(ph%proc) - 1)/ mpl)
                  associate (pb =>     lbound(ph%proc, 1) +  ln    * mpl, &
                       &     pe => min(lbound(ph%proc, 1) + (ln+1) * mpl - 1, ubound(ph%proc, 1)))

                     if (ph%wtime > 0.) then
                        write(fmt, *)"(3a,f6.", dec, ",a)"
                        write(header, fmt) "@", ph%nodename(:pnames%maxnamelen), " <MHD speed> = ", 1./ph%wtime, " blk/s ["
                        write(fmt,  *) "(a,", pe - pb + 1, "f6.", dec, ",a)"
                        write(msg, fmt) merge(trim(header), repeat(" ", len_trim(header)), ln == 0), &
                             merge(1. / pnames%wtime(ph%proc(pb:pe)), 0., pnames%wtime(ph%proc(pb:pe)) > 0.), &
                             merge(" ]", "  ", ln == int((size(ph%proc) - 1)/ mpl))
                     else
                        write(msg, '(3a)') "@", ph%nodename(:pnames%maxnamelen), " <MHD speed> = N/A"
                     endif
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

      !>
      !! \brief Detects uneven load balance by checking the time spent on cgs and suggests a rebalance.
      !!
      !! Depends on prior call to print_costs.
      !<

      subroutine speed_check

         use load_balance, only: imbalance_tol, rebalance_asap, flexible_balance
         use mpisetup,     only: slave
         use procnames,    only: pnames

         implicit none

         integer :: p, n
         real :: min_s, max_s
         integer, parameter :: s_ind = I_ONE  ! can be anything within size(cost_labels) because all values are the same

         if (slave .or. .not. flexible_balance) return

         n = 0
         min_s = huge(1.)
         max_s = -min_s

         do p = lbound(all_proc_stats, 3), ubound(all_proc_stats, 3)
            associate (s => all_proc_stats(N_STATS, s_ind, p))
               if (pnames%wtime(p) > epsilon(1.)) then
                  if (s < min_s) min_s = s
                  if (s > max_s) max_s = s
                  n = n + 1
               endif
            end associate
         enddo

         if (n > 0) then
            if (max_s * imbalance_tol > min_s) rebalance_asap = .true.
         endif

      end subroutine speed_check

   end subroutine print_costs

end module lb_helpers
