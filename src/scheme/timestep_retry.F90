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
!! \brief Retry the timestep when Courant number or other problems are detected
!<

module timestep_retry

   use named_array_list, only: na_var_list

   implicit none

   private
   public :: repeat_fluidstep, reset_freezing_speed

   ! for simplicity create an array of pointers to qna and wna
   type :: na_p
      class(na_var_list), pointer :: p => null()
      integer, allocatable, dimension(:) :: indices
   end type na_p
   type(na_p), dimension(2) :: na_lists  ! currently we have only 2 such lists, fna may join in the future

contains


!<
!! \brief Find the name for backup field when repeating timestep is allowed.
!>

   function get_rname(name) result(rname)

      use constants,  only: dsetnamelen
      use dataio_pub, only: die, msg

      implicit none

      character(len=dsetnamelen), intent(in) :: name

      character(len=dsetnamelen), parameter :: retry_n = '.retry'  !< append to mark an auto-generated copy for timestep retrying
      character(len=dsetnamelen) :: rname

      if (len_trim(name) + len_trim(retry_n) > dsetnamelen) then
         write(msg,'(3a,i3,a)')"[timestep_retry:get_rname] field name '", name, "' is too long. Try to shorten it to ", dsetnamelen - len_trim(retry_n), " characters"
         call die(msg)
      endif
      write(rname, '(2a)') trim(name), trim(retry_n)

   end function get_rname

!<
!! \brief Save the current state after correct time-step or restore previously saved state and try a shorter timestep.
!!
!! \warning There might be other evolving variables (such as mass_defect::magic_mass) that should be added here
!!
!! \todo Move this routine somewhere else, because it should be available for all hydro schemes
!!
!! \todo Also take care of particles.
!>

   subroutine repeat_fluidstep

      use allreduce,        only: piernik_MPI_Allreduce
      use cg_cost_data,     only: I_OTHER
      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: pSUM, I_ZERO, I_ONE, dsetnamelen, AT_IGNORE, INVALID
      use dataio_pub,       only: warn, msg, die
      use global,           only: dt, dt_full, t, t_saved, max_redostep_attempts, nstep, nstep_saved, dt_cur_shrink, repetitive_steps, tstep_attempt
      use mass_defect,      only: downgrade_magic_mass
      use mpisetup,         only: master
      use named_array_list, only: qna, wna, na_var_list_q, na_var_list_w
      use ppp,              only: ppp_main
      use repeatstep,       only: repeat_step
      use user_hooks,       only: user_reaction_to_redo_step
#if defined(GRAV) && defined(NBODY)
      use domain,           only: dom
      use particle_func,    only: particle_in_area
      use particle_types,   only: particle, P_ID, P_MASS, P_POS_X, P_POS_Z, P_VEL_X, P_VEL_Z, P_ACC_X, P_ACC_Z, P_ENER, P_TFORM, P_TDYN, npf
      use particle_utils,   only: is_part_in_cg
#endif /* GRAV && NBODY */
#ifdef RANDOMIZE
      use randomization,    only: randoms_redostep
#endif /* RANDOMIZE */

      implicit none

      type(cg_list_element), pointer :: cgl
      integer(kind=4)                :: no_hist_count
      integer                        :: i, j
      character(len=dsetnamelen)     :: rname
      character(len=*), parameter :: rs_label = "repeat_step_"
#if defined(GRAV) && defined(NBODY)
      logical :: in, phy, out, fin, indomain
      integer :: p
      type(particle), pointer :: part
#endif /* GRAV && NBODY */

      if (.not. repetitive_steps) return

      call restart_arrays

      if (repeat_step()) then
         tstep_attempt = tstep_attempt + I_ONE
         if (tstep_attempt > max_redostep_attempts) then
            write(msg, '(a,i2,a)')"[timestep_retry:repeat_fluidstep] tstep_attempt > ", max_redostep_attempts, " (max_redostep_attempts)"
            call die(msg)
         endif
         write(msg, '(a,i2,a)') "[timestep_retry:repeat_fluidstep] Redoing previous step (", tstep_attempt, ")"
         if (master) call warn(msg)
         t = t_saved
         nstep = nstep_saved
         dt = dt_full * dt_cur_shrink
         call reset_freezing_speed
         call downgrade_magic_mass
         if (associated(user_reaction_to_redo_step)) call user_reaction_to_redo_step
      else
         tstep_attempt = I_ZERO
         nstep_saved = nstep
         t_saved = t
      endif
#ifdef RANDOMIZE
      call randoms_redostep(repeat_step())
#endif /* RANDOMIZE */

      ! refresh na_lists(:)%indices
      do j = lbound(na_lists, dim=1), ubound(na_lists, dim=1)
         associate (na => na_lists(j)%p)
            allocate(na_lists(j)%indices(size(na%lst)))
            do i = lbound(na%lst(:), dim=1), ubound(na%lst(:), dim=1)
               if (na%lst(i)%restart_mode > AT_IGNORE) then
                  rname = get_rname(na%lst(i)%name)
                  select type(na)
                     type is (na_var_list_q)
                        na_lists(j)%indices(i) = qna%ind(rname)
                     type is (na_var_list_w)
                        na_lists(j)%indices(i) = wna%ind(rname)
                     class default
                        call die("[timestep_retry:repeat_fluidstep] unknown named array list type (rname)")
                  end select
               else
                  na_lists(j)%indices(i) = INVALID
               endif
            enddo
         end associate
      enddo

      if (repeat_step()) then
         call ppp_main%start(rs_label // "reverting")
      else
         call ppp_main%start(rs_label // "saving")
      endif
      no_hist_count = 0
      cgl => leaves%first
      do while (associated(cgl))
         call cgl%cg%costs%start

         ! No need to take care of any cgl%cg%q arrays as long as gravity is extrapolated from the previous timestep.
         ! error checking should've been done in restart_arrays, called few lines earlier
         do j = lbound(na_lists, dim=1), ubound(na_lists, dim=1)
            associate (na => na_lists(j)%p)
               do i = lbound(na%lst(:), dim=1), ubound(na%lst(:), dim=1)
                  if (na%lst(i)%restart_mode > AT_IGNORE) then
                     rname = get_rname(na%lst(i)%name)
                     if (repeat_step()) then
                        if (cgl%cg%has_previous_timestep) then
                           select type(na)  !! ToDo: unify qna and wna somehow at least in the grid_container
                              type is (na_var_list_q)
                                 cgl%cg%q(i)%arr = cgl%cg%q(na_lists(j)%indices(i))%arr
                              type is (na_var_list_w)
                                 cgl%cg%w(i)%arr = cgl%cg%w(na_lists(j)%indices(i))%arr
                              class default
                                 call die("[timestep_retry:repeat_fluidstep] unknown named array list type ->")
                           end select
                        else
                           no_hist_count = no_hist_count + I_ONE
                        endif
                     else
                        select type(na)
                           type is (na_var_list_q)
                              cgl%cg%q(na_lists(j)%indices(i))%arr = cgl%cg%q(i)%arr
                           type is (na_var_list_w)
                              cgl%cg%w(na_lists(j)%indices(i))%arr = cgl%cg%w(i)%arr
                           class default
                              call die("[timestep_retry:repeat_fluidstep] unknown named array list type <-")
                        end select
                        cgl%cg%has_previous_timestep = .true.
                     endif
                  endif
               enddo
            end associate
         enddo
#if defined(GRAV) && defined(NBODY)
         if (repeat_step()) then
            if (.not. allocated(cgl%cg%psave)) call die("timestep_retry:repeat_fluidstep] .not. allocated(cgl%cg%psave)")
            call cgl%cg%pset%cleanup
            do p = lbound(cgl%cg%psave, 2, kind=4), ubound(cgl%cg%psave, 2, kind=4)
               indomain = particle_in_area(cgl%cg%psave(P_POS_X:P_POS_Z, p), dom%edge)
               call is_part_in_cg(cgl%cg, cgl%cg%psave(P_POS_X:P_POS_Z, p), indomain, in, phy, out, fin)
               call cgl%cg%pset%add(nint(cgl%cg%psave(P_ID, p), kind=4), cgl%cg%psave(P_MASS, p), &
                    cgl%cg%psave(P_POS_X:P_POS_Z, p), cgl%cg%psave(P_VEL_X:P_VEL_Z, p), &
                    cgl%cg%psave(P_ACC_X:P_ACC_Z, p), cgl%cg%psave(P_ENER, p), &
                    in, phy, out, fin, &
                    cgl%cg%psave(P_TFORM, p), cgl%cg%psave(P_TDYN, p))
            enddo
         else
            if (allocated(cgl%cg%psave)) deallocate(cgl%cg%psave)
            allocate(cgl%cg%psave(npf, cgl%cg%count_all_particles()))
            part => cgl%cg%pset%first
            p = 1
            do while (associated(part))
               cgl%cg%psave(P_ID, p)            = part%pdata%pid
               cgl%cg%psave(P_MASS, p)          = part%pdata%mass
               cgl%cg%psave(P_POS_X:P_POS_Z, p) = part%pdata%pos
               cgl%cg%psave(P_VEL_X:P_VEL_Z, p) = part%pdata%vel
               cgl%cg%psave(P_ACC_X:P_ACC_Z, p) = part%pdata%acc
               cgl%cg%psave(P_ENER, p)          = part%pdata%energy
               cgl%cg%psave(P_TFORM, p)         = part%pdata%tform
               cgl%cg%psave(P_TDYN, p)          = part%pdata%tdyn

               p = p + I_ONE
               part => part%nxt
            enddo
         endif
#endif /* GRAV && NBODY */

         call cgl%cg%costs%stop(I_OTHER)
         cgl => cgl%nxt
      enddo
      if (repeat_step()) then
         call ppp_main%stop(rs_label // "reverting")
      else
         call ppp_main%stop(rs_label // "saving")
      endif

      call piernik_MPI_Allreduce(no_hist_count, pSUM)
      if (master .and. no_hist_count/=0) then
         write(msg, '(a,i6,a)')"[timestep_retry:repeat_fluidstep] Error: not reverted: ", no_hist_count, " grid pieces."
         call die(msg)
         ! AMR domains require careful treatment of timestep retries.
         ! Going back past rebalancing or refinement change would require updating whole AMR structure, not just field values.
      endif

      do j = lbound(na_lists, dim=1), ubound(na_lists, dim=1)
         deallocate(na_lists(j)%indices)
      enddo

   end subroutine repeat_fluidstep

!>
!! \brief Check if the arrays exist and create them if needed
!<

   subroutine restart_arrays

      use cg_list_global,   only: all_cg
      use constants,        only: AT_BACKUP, AT_IGNORE, dsetnamelen, INVALID, V_VERBOSE
      use dataio_pub,       only: printinfo, msg
      use mpisetup,         only: master
      use named_array_list, only: qna, wna

      implicit none

      integer :: i, j
      character(len=dsetnamelen) :: rname
      integer(kind=4), dimension(:), allocatable, target :: pos_copy

      if (.not. associated(na_lists(1)%p)) na_lists(1)%p => qna
      if (.not. associated(na_lists(2)%p)) na_lists(2)%p => wna

      do j = lbound(na_lists, dim=1), ubound(na_lists, dim=1)
         associate (na => na_lists(j)%p)

            do i = lbound(na%lst(:), dim=1), ubound(na%lst(:), dim=1)
               if (na%lst(i)%restart_mode > AT_IGNORE) then
                  rname = get_rname(na%lst(i)%name)
                  if (.not. na%exists(rname)) then
                     if (master) then
                        write(msg,'(3a)')"[timestep_retry:restart_arrays] creating backup field '", rname, "'"
                        call printinfo(msg, V_VERBOSE)
                     endif
                     allocate(pos_copy(size(na%lst(i)%position)))
                     pos_copy = na%lst(i)%position
                     if (na%lst(i)%dim4 /= INVALID) then
                        call all_cg%reg_var(rname, dim4=na%lst(i)%dim4, position=pos_copy, multigrid=na%lst(i)%multigrid, restart_mode = AT_BACKUP)
                     else
                        call all_cg%reg_var(rname,                      position=pos_copy, multigrid=na%lst(i)%multigrid, restart_mode = AT_BACKUP)
                     endif
                     deallocate(pos_copy)
                  endif
               endif
            enddo
         end associate
      enddo

   end subroutine restart_arrays

   subroutine reset_freezing_speed

      use fluidindex, only: flind

      implicit none

      integer :: ifl

      do ifl = lbound(flind%all_fluids, dim=1), ubound(flind%all_fluids, dim=1)
         call flind%all_fluids(ifl)%fl%res_c
      enddo

   end subroutine reset_freezing_speed

end module timestep_retry
