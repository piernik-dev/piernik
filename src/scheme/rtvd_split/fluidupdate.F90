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
!! \brief COMMENT ME
!<

module fluidupdate   ! SPLIT
! pulled by RTVD

   use named_array_list, only: na_var_list

   implicit none

   private
   public :: fluid_update

   ! for simplicity create an array of pointers to qna and wna
   type :: na_p
      class(na_var_list), pointer :: p => null()
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
         write(msg,'(3a,i3,a)')"[fluidupdate:get_rname] field name '", name, "' is too long. Try to shorten it to ", dsetnamelen - len_trim(retry_n), " characters"
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
!>

   subroutine repeat_fluidstep

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: pSUM, I_ONE, dsetnamelen, AT_IGNORE
      use dataio_pub,       only: warn, msg, die
      use global,           only: dt, dtm, t, t_saved, cfl_violated, nstep, nstep_saved, dt_max_grow, repeat_step
      use mass_defect,      only: downgrade_magic_mass
      use mpisetup,         only: master, piernik_MPI_Allreduce
      use named_array_list, only: qna, wna, na_var_list_q, na_var_list_w
      use user_hooks,       only: user_reaction_to_redo_step

      implicit none

      type(cg_list_element), pointer :: cgl
      integer(kind=4) :: no_hist_count
      integer :: i, j
      character(len=dsetnamelen) :: rname

      if (.not.repeat_step) return

      call restart_arrays

      if (cfl_violated) then
         if (master) call warn("[fluidupdate:fluid_update] Redoing previous step...")
         t = t_saved
         nstep = nstep_saved
         dt = dtm/dt_max_grow**2
      else
         nstep_saved = nstep
         t_saved = t
      endif

      no_hist_count = 0
      cgl => leaves%first
      do while (associated(cgl))
         ! No need to take care of any cgl%cg%q arrays as long as graity is extrapolated from the prefious timestep.

         ! error checking should've been done in restart_arrays, called few lines earlier
         do j = lbound(na_lists, dim=1), ubound(na_lists, dim=1)
            associate (na => na_lists(j)%p)
               do i = lbound(na%lst(:), dim=1), ubound(na%lst(:), dim=1)
                  if (na%lst(i)%restart_mode /= AT_IGNORE) then
                     rname = get_rname(na%lst(i)%name)
                     if (cfl_violated) then
                        if (cgl%cg%has_previous_timestep) then
                           select type(na)  !! ToDo: unify qna and wna somehow at least in the grid_container
                              type is (na_var_list_q)
                                 cgl%cg%q(i)%arr = cgl%cg%q(qna%ind(rname))%arr
                              type is (na_var_list_w)
                                 cgl%cg%w(i)%arr = cgl%cg%w(wna%ind(rname))%arr
                              class default
                                 call die("[fluidupdate:repeat_fluidstep] unknown named array list type ->")
                           end select
                        else
                           no_hist_count = no_hist_count + I_ONE
                        endif
                        call downgrade_magic_mass
                        if (associated(user_reaction_to_redo_step)) call user_reaction_to_redo_step
                     else
                        select type(na)
                           type is (na_var_list_q)
                              cgl%cg%q(qna%ind(rname))%arr = cgl%cg%q(i)%arr
                           type is (na_var_list_w)
                              cgl%cg%w(wna%ind(rname))%arr = cgl%cg%w(i)%arr
                           class default
                              call die("[fluidupdate:repeat_fluidstep] unknown named array list type <-")
                        end select
                        cgl%cg%has_previous_timestep = .true.
                     endif
                  endif
               enddo
            end associate
         enddo
         cgl => cgl%nxt
      enddo
      call piernik_MPI_Allreduce(no_hist_count, pSUM)
      if (master .and. no_hist_count/=0) then
         write(msg, '(a,i6,a)')"[fluidupdate:repeat_fluidstep] Warning: not reverted: ", no_hist_count, " grid pieces."
         call die(msg)
         ! AMR domains require careful treatment of timestep retries.
         ! Going back past rebalancing or refinement change would require updating whole AMR structure, not just field values.
      endif

   end subroutine repeat_fluidstep

!>
!! \brief Check if the arrays exist and create them if needed
!<

   subroutine restart_arrays

      use cg_list_global,   only: all_cg
      use constants,        only: AT_IGNORE, dsetnamelen, INVALID
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
               if (na%lst(i)%restart_mode /= AT_IGNORE) then
                  rname = get_rname(na%lst(i)%name)
                  if (.not. na%exists(rname)) then
                     if (master) then
                        write(msg,'(3a)')"[fluidupdate:restart_arrays] creating backup field '", rname, "'"
                        call printinfo(msg)
                     endif
                     allocate(pos_copy(size(na%lst(i)%position)))
                     pos_copy = na%lst(i)%position
                     if (na%lst(i)%dim4 /= INVALID) then
                        call all_cg%reg_var(rname, dim4=na%lst(i)%dim4, position=pos_copy, multigrid=na%lst(i)%multigrid)
                     else
                        call all_cg%reg_var(rname,                      position=pos_copy, multigrid=na%lst(i)%multigrid)
                     endif
                     deallocate(pos_copy)
                  endif
               endif
            enddo
         end associate
      enddo

   end subroutine restart_arrays

!>
!! \brief Advance the solution by two timesteps using directional splitting
!<

   subroutine fluid_update

      use dataio_pub,  only: halfstep
      use global,      only: dt, dtm, t
      use mass_defect, only: update_magic_mass

      implicit none

      call repeat_fluidstep

      halfstep = .false.
      t=t+dt
      call make_3sweeps(.true.) ! X -> Y -> Z

! Sources should be hooked to problem_customize_solution with forward argument

      halfstep = .true.
      t=t+dt
      dtm = dt
      call make_3sweeps(.false.) ! Z -> Y -> X
      call update_magic_mass

   end subroutine fluid_update

!>
!! \brief Perform sweeps in all three directions plus sources that are calculated every timestep
!<
   subroutine make_3sweeps(forward)

      use cg_list_dataop,      only: expanded_domain
      use constants,           only: xdim, ydim, zdim, I_ONE
      use global,              only: skip_sweep, use_fargo
      use fargo,               only: make_fargosweep
      use sweeps,              only: sweep
      use user_hooks,          only: problem_customize_solution
#ifdef GRAV
      use global,              only: t, dt
      use gravity,             only: source_terms_grav
      use particle_pub,        only: pset, psolver
#endif /* GRAV */
#if defined(COSM_RAYS) && defined(MULTIGRID)
      use all_boundaries,      only: all_fluid_boundaries
      use initcosmicrays,      only: use_split
      use multigrid_diffusion, only: multigrid_solve_diff
#endif /* COSM_RAYS && MULTIGRID */
#ifdef SHEAR
      use shear,               only: shear_3sweeps
#endif /* SHEAR */

      implicit none

      logical, intent(in) :: forward  !< If .true. then do X->Y->Z sweeps, if .false. then reverse that order

      integer(kind=4) :: s

#ifdef SHEAR
      call shear_3sweeps
#endif /* SHEAR */

#ifdef GRAV
      call source_terms_grav
#endif /* GRAV */

#if defined(COSM_RAYS) && defined(MULTIGRID)
      if (.not. use_split) then
         call multigrid_solve_diff
         call all_fluid_boundaries
      endif
#endif /* COSM_RAYS && MULTIGRID */

      call expanded_domain%delete ! at this point everything should be initialized after domain expansion and we no longer need this list
      if (use_fargo) then
         if (.not.skip_sweep(zdim)) call make_sweep(zdim, forward)
         if (.not.skip_sweep(xdim)) call make_sweep(xdim, forward)
         if (.not.skip_sweep(ydim)) call make_fargosweep
      else
         if (forward) then
            do s = xdim, zdim
               if (.not.skip_sweep(s)) call make_sweep(s, forward)
            enddo
         else
            do s = zdim, xdim, -I_ONE
               if (.not.skip_sweep(s)) call make_sweep(s, forward)
            enddo
         endif
      endif
#ifdef GRAV
      if (associated(psolver)) call pset%evolve(psolver, t-dt, dt)
#endif /* GRAV */
      if (associated(problem_customize_solution)) call problem_customize_solution(forward)

   end subroutine make_3sweeps

!>
!! \brief Perform single sweep in forward or backward direction
!<
   subroutine make_sweep(dir, forward)

      use constants,      only: DIVB_CT
      use dataio_pub,     only: die
      use domain,         only: dom
      use global,         only: geometry25D, divB_0_method
      use sweeps,         only: sweep
#ifdef COSM_RAYS
      use crdiffusion,    only: cr_diff
      use initcosmicrays, only: use_split
#endif /* COSM_RAYS */
#ifdef DEBUG
      use piernikiodebug,   only: force_dumps
#endif /* DEBUG */

      implicit none

      integer(kind=4), intent(in) :: dir      !< direction, one of xdim, ydim, zdim
      logical,         intent(in) :: forward  !< if .false. then reverse operation order in the sweep

      if (divB_0_method /= DIVB_CT) call die("[fluidupdate:make_sweep] only CT is implemented in RTVD")

      if (dom%has_dir(dir)) then
         if (.not. forward) then
#ifdef COSM_RAYS
            if (use_split) call cr_diff(dir)
#endif /* COSM_RAYS */
#ifdef MAGNETIC
            call magfield(dir)
#endif /* MAGNETIC */
         endif

         call sweep(dir)

         if (forward) then
#ifdef MAGNETIC
            call magfield(dir)
#endif /* MAGNETIC */
#ifdef COSM_RAYS
            if (use_split) call cr_diff(dir)
#endif /* COSM_RAYS */
         endif
      else
         if (geometry25D) call sweep(dir)
      endif

#ifdef DEBUG
      call force_dumps
#endif /* DEBUG */

   end subroutine make_sweep

#ifdef MAGNETIC
   subroutine magfield(dir)

      use advects,     only: advectb
      use constants,   only: ndims, I_ONE
#ifdef RESISTIVE
      use resistivity, only: diffuseb
#endif /* RESISTIVE */

      implicit none

      integer(kind=4), intent(in) :: dir

      integer(kind=4)             :: bdir, dstep

      do dstep = 0, 1
         bdir  = I_ONE + mod(dir+dstep,ndims)
         call advectb(bdir, dir)
#ifdef RESISTIVE
         call diffuseb(bdir, dir)
#endif /* RESISTIVE */
         call mag_add(dir, bdir)
      enddo

   end subroutine magfield

!------------------------------------------------------------------------------------------

   subroutine mag_add(dim1, dim2)

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use grid_cont,        only: grid_container
      use all_boundaries,   only: all_mag_boundaries
      use user_hooks,       only: custom_emf_bnd
#ifdef RESISTIVE
      use constants,        only: wcu_n
      use dataio_pub,       only: die
      use domain,           only: is_multicg
      use named_array_list, only: qna
#endif /* RESISTIVE */

      implicit none

      integer(kind=4), intent(in)    :: dim1, dim2

      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
#ifdef RESISTIVE
      real, dimension(:,:,:), pointer :: wcu
#endif /* RESISTIVE */

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
#ifdef RESISTIVE
! DIFFUSION FULL STEP
         wcu => cg%q(qna%ind(wcu_n))%arr
         if (is_multicg) call die("[fluidupdate:mag_add] multiple grid pieces per processor not implemented yet") ! not tested custom_emf_bnd
         if (associated(custom_emf_bnd)) call custom_emf_bnd(wcu)
         cg%b(dim2,:,:,:) = cg%b(dim2,:,:,:) -              wcu*cg%idl(dim1)
         cg%b(dim2,:,:,:) = cg%b(dim2,:,:,:) + pshift(wcu,dim1)*cg%idl(dim1)
         cg%b(dim1,:,:,:) = cg%b(dim1,:,:,:) +              wcu*cg%idl(dim2)
         cg%b(dim1,:,:,:) = cg%b(dim1,:,:,:) - pshift(wcu,dim2)*cg%idl(dim2)
#endif /* RESISTIVE */
! ADVECTION FULL STEP
         if (associated(custom_emf_bnd)) call custom_emf_bnd(cg%wa)
         cg%b(dim2,:,:,:) = cg%b(dim2,:,:,:) - cg%wa*cg%idl(dim1)
         cg%wa = mshift(cg%wa,dim1)
         cg%b(dim2,:,:,:) = cg%b(dim2,:,:,:) + cg%wa*cg%idl(dim1)
         cg%b(dim1,:,:,:) = cg%b(dim1,:,:,:) - cg%wa*cg%idl(dim2)
         cg%wa = pshift(cg%wa,dim2)
         cg%b(dim1,:,:,:) = cg%b(dim1,:,:,:) + cg%wa*cg%idl(dim2)
         cgl => cgl%nxt
      enddo

      call all_mag_boundaries

   end subroutine mag_add
!------------------------------------------------------------------------------------------

!>
!! \brief Function pshift makes one-cell, forward circular shift of 3D array in any direction
!! \param tab input array
!! \param d direction of the shift, where 1,2,3 corresponds to \a x,\a y,\a z respectively
!! \return real, dimension(size(tab,1),size(tab,2),size(tab,3))
!!
!! The function was written in order to significantly improve
!! the performance at the cost of the flexibility of original \p CSHIFT.
!<
   function pshift(tab, d)

      use dataio_pub,    only: warn

      implicit none

      real, dimension(:,:,:), intent(inout) :: tab
      integer(kind=4),        intent(in)    :: d

      integer :: ll
      real, dimension(size(tab,1),size(tab,2),size(tab,3)) :: pshift

      ll = size(tab,d)

      if (ll==1) then
         pshift = tab
         return
      endif

      if (d==1) then
         pshift(1:ll-1,:,:) = tab(2:ll,:,:); pshift(ll,:,:) = tab(1,:,:)
      else if (d==2) then
         pshift(:,1:ll-1,:) = tab(:,2:ll,:); pshift(:,ll,:) = tab(:,1,:)
      else if (d==3) then
         pshift(:,:,1:ll-1) = tab(:,:,2:ll); pshift(:,:,ll) = tab(:,:,1)
      else
         call warn('[fluidupdate:pshift]: Dim ill defined in pshift!')
      endif

      return
   end function pshift

!>
!! \brief Function mshift makes one-cell, backward circular shift of 3D array in any direction
!! \param tab input array
!! \param d direction of the shift, where 1,2,3 corresponds to \a x,\a y,\a z respectively
!! \return real, dimension(size(tab,1),size(tab,2),size(tab,3))
!!
!! The function was written in order to significantly improve
!! the performance at the cost of the flexibility of original \p CSHIFT.
!<
   function mshift(tab,d)

      use dataio_pub,    only: warn

      implicit none

      real, dimension(:,:,:), intent(inout) :: tab
      integer(kind=4),        intent(in)    :: d

      integer :: ll
      real, dimension(size(tab,1) , size(tab,2) , size(tab,3)) :: mshift

      ll = size(tab,d)

      if (ll==1) then
         mshift = tab
         return
      endif

      if (d==1) then
         mshift(2:ll,:,:) = tab(1:ll-1,:,:); mshift(1,:,:) = tab(ll,:,:)
      else if (d==2) then
         mshift(:,2:ll,:) = tab(:,1:ll-1,:); mshift(:,1,:) = tab(:,ll,:)
      else if (d==3) then
         mshift(:,:,2:ll) = tab(:,:,1:ll-1); mshift(:,:,1) = tab(:,:,ll)
      else
         call warn('[fluidupdate:mshift]: Dim ill defined in mshift!')
      endif

      return
   end function mshift

#endif /* MAGNETIC */

end module fluidupdate
