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
module fluidupdate   ! SPLIT
! pulled by RTVD
   implicit none

   private
   public :: fluid_update

contains

   subroutine repeat_fluidstep

      use constants,  only: I_ONE, u0_n, b0_n
      use dataio_pub, only: warn
      use gc_list,    only: cg_list_element
      use global,     only: dt, dtm, t, cfl_violated, nstep, dt_max_grow, repeat_step
      use grid,       only: all_cg
      use grid_cont,  only: grid_container
      use mpisetup,   only: master

      implicit none

      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      if (.not.repeat_step) return

      cgl => all_cg%first
      do while (associated(cgl))
         cg => cgl%cg

         if (cfl_violated) then
            t = t-2.0*dtm
            cg%u = cg%w(cg%get_na_ind_4d(u0_n))%arr
            cg%b = cg%w(cg%get_na_ind_4d(b0_n))%arr
            dt = dtm/dt_max_grow**2
            nstep = nstep - I_ONE
            if (master) call warn("[fluidupdate:fluid_update] Redoing previous step...")
         else
            cg%w(cg%get_na_ind_4d(u0_n))%arr = cg%u
            cg%w(cg%get_na_ind_4d(b0_n))%arr = cg%b
         endif

         cgl => cgl%nxt
      enddo

   end subroutine repeat_fluidstep

   subroutine fluid_update

      use dataio_pub, only: halfstep
      use global,     only: dt, dtm, t
#ifdef SN_SRC
      use snsources,  only: random_sn
#endif /* SN_SRC */

      implicit none

      call repeat_fluidstep

      halfstep = .false.
      t=t+dt
      call make_3sweeps(.true.) ! X -> Y -> Z

! Sources ----------------------------------------

#ifdef SN_SRC
      call random_sn !> \todo hook this to problem_customize_solution
#endif /* SN_SRC */

      halfstep = .true.
      t=t+dt
      dtm = dt
      call make_3sweeps(.false.) ! Z -> Y -> X

   end subroutine fluid_update

!>
!! \brief Perform sweeps in all three directions plus sources that are calculated every timestep
!<
   subroutine make_3sweeps(forward)

      use constants,           only: xdim, ydim, zdim, I_ONE
      use user_hooks,          only: problem_customize_solution
      use global,              only: skip_sweep
#ifdef SHEAR
      use dataio_pub,          only: die
      use domain,              only: dom, is_multicg
      use fluidboundaries,     only: bnd_u
      use global,              only: t, dt
      use grid,                only: all_cg
      use grid_cont,           only: grid_container
      use shear,               only: yshift
#endif /* SHEAR */
#ifdef GRAV
      use gravity,             only: source_terms_grav
#endif /* GRAV */
#ifdef COSM_RAYS
      use initcosmicrays,      only: use_split
#ifdef MULTIGRID
      use multigrid_diffusion, only: multigrid_solve_diff
      use fluidboundaries,     only: all_fluid_boundaries
#endif /* MULTIGRID */
#endif /* COSM_RAYS */

      implicit none

      logical, intent(in) :: forward  !< If .true. then do X->Y->Z sweeps, if .false. then reverse that order

      integer(kind=4) :: s
#ifdef SHEAR
      type(grid_container), pointer :: cg

      cg => all_cg%first%cg
      if (is_multicg) call die("[fluidupdate:make_3sweeps] multiple grid pieces per procesor not implemented yet") !nontrivial SHEAR

      if (dom%has_dir(ydim)) call yshift(t, dt)
      if (dom%has_dir(xdim)) call bnd_u(xdim, cg)
      if (dom%has_dir(ydim)) call bnd_u(ydim, cg)
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

      if (forward) then
         do s = xdim, zdim
            if (.not.skip_sweep(s)) call make_sweep(s, forward)
         enddo
      else
         do s = zdim, xdim, -I_ONE
            if (.not.skip_sweep(s)) call make_sweep(s, forward)
         enddo
      endif
      if (associated(problem_customize_solution)) call problem_customize_solution

   end subroutine make_3sweeps

!>
!! \brief Perform single sweep in forward or backward direction
!<
   subroutine make_sweep(dir, forward)

      use constants,      only: ydim
      use domain,         only: dom
      use sweeps,         only: sweep
#if defined SHEAR && defined FLUID_INTERACTIONS
      use sweeps,         only: source_terms_y
#endif /* SHEAR && FLUID_INTERACTIONS */
#ifdef COSM_RAYS
      use crdiffusion,    only: cr_diff
      use initcosmicrays, only: use_split
#endif /* COSM_RAYS */
#ifdef DEBUG
      use piernikdebug,   only: force_dumps
#endif /* DEBUG */

      implicit none

      integer(kind=4), intent(in) :: dir      !< direction, one of xdim, ydim, zdim
      logical, intent(in) :: forward  !< if .false. then reverse operation order in the sweep

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
#if defined SHEAR && defined FLUID_INTERACTIONS
         if (dir == ydim) call source_terms_y
#endif /* SHEAR && FLUID_INTERACTIONS */
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

      use func,          only: pshift, mshift
      use grid,          only: all_cg
      use gc_list,       only: cg_list_element
      use grid_cont,     only: grid_container
      use magboundaries, only: all_mag_boundaries
      use user_hooks,    only: custom_emf_bnd
#ifdef RESISTIVE
      use dataio_pub,    only: die
      use domain,        only: is_multicg
      use resistivity,   only: wcu
#endif /* RESISTIVE */

      implicit none

      integer(kind=4), intent(in)    :: dim1, dim2
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg

      cgl => all_cg%first
      do while (associated(cgl))
         cg => cgl%cg
#ifdef RESISTIVE
! DIFFUSION FULL STEP
         if (is_multicg) call die("[fluidupdate:mag_add] multiple grid pieces per procesor not implemented yet") ! move wcu into cg
         if (associated(custom_emf_bnd)) call custom_emf_bnd(wcu%arr)
         cg%b(dim2,:,:,:) = cg%b(dim2,:,:,:) - wcu%arr*cg%idl(dim1)
         wcu%arr = pshift(wcu%arr,dim1)
         cg%b(dim2,:,:,:) = cg%b(dim2,:,:,:) + wcu%arr*cg%idl(dim1)
         wcu%arr = mshift(wcu%arr,dim1)
         cg%b(dim1,:,:,:) = cg%b(dim1,:,:,:) + wcu%arr*cg%idl(dim2)
         wcu%arr = pshift(wcu%arr,dim2)
         cg%b(dim1,:,:,:) = cg%b(dim1,:,:,:) - wcu%arr*cg%idl(dim2)
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
#endif /* MAGNETIC */
!------------------------------------------------------------------------------------------

end module fluidupdate
