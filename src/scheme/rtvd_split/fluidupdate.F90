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
! pulled by ANY
   implicit none

   private
   public :: fluid_update

contains

   subroutine repeat_fluidstep

      use dataio_pub, only: warn
      use global,     only: dt, dtm, t, cfl_violated, nstep, dt_max_grow, repeat_step
      use grid,       only: cga
      use grid_cont,  only: cg_list_element, grid_container
      use mpisetup,   only: master


      implicit none

      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      if (.not.repeat_step) return

      call cga%get_root(cgl)
      do while (associated(cgl))
         cg => cgl%cg

         if (cfl_violated) then
            t = t-2.0*dtm
            cg%u%arr = cg%u0%arr
            cg%b%arr = cg%b0%arr
            dt = dtm/dt_max_grow**2
            nstep = nstep-1
            if (master) call warn("[fluidupdate:fluid_update] Redoing previous step...")
         else
            cg%u0%arr = cg%u%arr
            cg%b0%arr = cg%b%arr
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

      use constants,           only: xdim, ydim, zdim
      use types,               only: problem_customize_solution
#ifdef SHEAR
      use dataio_pub,          only: die
      use domain,              only: has_dir
      use fluidboundaries,     only: bnd_u
      use global,              only: t, dt
      use grid,                only: cga
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

      integer :: s
#ifdef SHEAR
      type(grid_container), pointer :: cg

      cg => cga%cg_all(1)
      if (ubound(cga%cg_all(:), dim=1) > 1) call die("[fluidupdate:make_3sweeps] multiple grid pieces per procesor not implemented yet") !nontrivial SHEAR

      if (has_dir(ydim)) call yshift(t, dt)
      if (has_dir(xdim)) call bnd_u(xdim, cg)
      if (has_dir(ydim)) call bnd_u(ydim, cg)
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
            call make_sweep(s, forward)
         enddo
      else
         do s = zdim, xdim, -1
            call make_sweep(s, forward)
         enddo
      endif
      if (associated(problem_customize_solution)) call problem_customize_solution

   end subroutine make_3sweeps

!>
!! \brief Perform single sweep in forward or backward direction
!<
   subroutine make_sweep(dir, forward)

      use constants,      only: xdim, ydim, zdim
      use dataio_pub,     only: msg, die
      use domain,         only: has_dir
      use grid,           only: cga
      use grid_cont,      only: cg_list_element, grid_container
      use sweeps,         only: sweepx, sweepy, sweepz
#if defined SHEAR && defined FLUID_INTERACTIONS
      use sweeps,         only: source_terms_y
#endif /* SHEAR && FLUID_INTERACTIONS */
#ifdef COSM_RAYS
      use crdiffusion,    only: cr_diff_x, cr_diff_y, cr_diff_z
      use initcosmicrays, only: use_split
#endif /* COSM_RAYS */
#ifdef DEBUG
      use dataio,         only: write_data
      use dataio_hdf5,    only: write_hdf5, write_restart_hdf5
      use dataio_pub,     only: chdf, set_container_chdf
      use global,         only: nstep
      use piernikdebug,   only: force_hdf5_dump, force_log_dump, force_res_dump, force_allbnd_dump
#endif /* DEBUG */

      implicit none

      integer, intent(in) :: dir      !< direction, one of xdim, ydim, zdim
      logical, intent(in) :: forward  !< if .false. then reverse operation order in the sweep

      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      call cga%get_root(cgl)
      do while (associated(cgl))
         cg => cgl%cg

         select case (dir)

            case (xdim)
               if (has_dir(xdim)) then
                  if (.not. forward) then
#ifdef COSM_RAYS
                     if (use_split) call cr_diff_x
#endif /* COSM_RAYS */
#ifdef MAGNETIC
                     call magfieldbyzx
#endif /* MAGNETIC */
                  endif

                  call sweepx(cg)

                  if (forward) then
#ifdef MAGNETIC
                     call magfieldbyzx
#endif /* MAGNETIC */
#ifdef COSM_RAYS
                     if (use_split) call cr_diff_x
#endif /* COSM_RAYS */
                  endif
               endif

            case (ydim)
               if (has_dir(ydim)) then
                  if (.not. forward) then
#ifdef COSM_RAYS
                     if (use_split) call cr_diff_y
#endif /* COSM_RAYS */
#ifdef MAGNETIC
                     call magfieldbzxy
#endif /* MAGNETIC */
                  endif
                  call sweepy(cg)

                  if (forward) then
#ifdef MAGNETIC
                     call magfieldbzxy
#endif /* MAGNETIC */
#ifdef COSM_RAYS
                     if (use_split) call cr_diff_y
#endif /* COSM_RAYS */
                  endif
               else
#if defined SHEAR && defined FLUID_INTERACTIONS
                  call source_terms_y
#endif /* SHEAR && FLUID_INTERACTIONS */
               endif

            case (zdim)
               if (has_dir(zdim)) then
                  if (.not. forward) then
#ifdef COSM_RAYS
                     if (use_split) call cr_diff_z
#endif /* COSM_RAYS */
#ifdef MAGNETIC
                     call magfieldbxyz
#endif /* MAGNETIC */
                  endif

                  call sweepz(cg)

                  if (forward) then
#ifdef MAGNETIC
                     call magfieldbxyz
#endif /* MAGNETIC */
#ifdef COSM_RAYS
                     if (use_split) call cr_diff_z
#endif /* COSM_RAYS */
                  endif
               endif

            case default
               write(msg,'(a,i10)')"[fluidupdate:make_sweep] Illegal direction ",dir
               call die(msg)

         end select

         cgl => cgl%nxt
      enddo

#ifdef DEBUG
      call set_container_chdf(nstep)
      if (force_hdf5_dump)   call write_hdf5(chdf)
      if (force_res_dump)    call write_restart_hdf5
      if (force_allbnd_dump) call write_restart_hdf5(debug_res=.true.)
      if (force_log_dump)    call write_data(output='log')
#endif /* DEBUG */

   end subroutine make_sweep

#ifdef MAGNETIC
   subroutine magfieldbyzx

      use advects,     only: advectby_x, advectbz_x
      use fluidindex,  only: ibx, iby, ibz
      use constants,   only: xdim, ydim, zdim
#ifdef RESISTIVE
      use resistivity, only: diffuseb
#endif /* RESISTIVE */

      implicit none

      call advectby_x

#ifdef RESISTIVE
      call diffuseb(iby,xdim,zdim,'emfz',ydim,zdim)
#endif /* RESISTIVE */

      call mag_add(iby,xdim,ibx,ydim)

      call advectbz_x

#ifdef RESISTIVE
      call diffuseb(ibz,xdim,ydim,'emfy',ydim,zdim)
#endif /* RESISTIVE */

      call mag_add(ibz,xdim,ibx,zdim)

   end subroutine magfieldbyzx

!------------------------------------------------------------------------------------------

   subroutine magfieldbzxy

      use advects,     only: advectbx_y, advectbz_y
      use fluidindex,  only: ibx, iby, ibz
      use constants,   only: xdim, ydim, zdim
#ifdef RESISTIVE
      use resistivity, only: diffuseb
#endif /* RESISTIVE */

      implicit none

      call advectbz_y

#ifdef RESISTIVE
      call diffuseb(ibz,ydim,xdim,'emfx',zdim,xdim)
#endif /* RESISTIVE */

      call mag_add(ibz,ydim,iby,zdim)

      call advectbx_y

#ifdef RESISTIVE
      call diffuseb(ibx,ydim,zdim,'emfz',zdim,xdim)
#endif /* RESISTIVE */

      call mag_add(ibx,ydim,iby,xdim)

   end subroutine magfieldbzxy

!------------------------------------------------------------------------------------------

   subroutine magfieldbxyz

      use advects,     only: advectbx_z, advectby_z
      use fluidindex,  only: ibx, iby, ibz
      use constants,   only: xdim, ydim, zdim
#ifdef RESISTIVE
      use resistivity, only: diffuseb
#endif /* RESISTIVE */

      implicit none

      call advectbx_z
#ifdef RESISTIVE
      call diffuseb(ibx,zdim,ydim,'emfy',xdim,ydim)
#endif /* RESISTIVE */
      call mag_add(ibx,zdim,ibz,xdim)

      call advectby_z
#ifdef RESISTIVE
      call diffuseb(iby,zdim,xdim,'emfx',xdim,ydim)
#endif /* RESISTIVE */

      call mag_add(iby,zdim,ibz,ydim)
   end subroutine magfieldbxyz

!------------------------------------------------------------------------------------------

   subroutine mag_add(ib1,dim1,ib2,dim2)

      use dataio_pub,    only: die
      use func,          only: pshift, mshift
      use grid,          only: cga
      use grid_cont,     only: cg_list_element, grid_container
      use magboundaries, only: all_mag_boundaries
      use types,         only: custom_emf_bnd
#ifdef RESISTIVE
      use resistivity,   only: wcu
#endif /* RESISTIVE */

      implicit none

      integer, intent(in)            :: ib1,ib2,dim1,dim2
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg

      if (ubound(cga%cg_all(:), dim=1) > 1) call die("[fluidupdate:mag_add] multiple grid pieces per procesor not implemented yet") !nontrivial not really checked

      call cga%get_root(cgl)
      do while (associated(cgl))
         cg => cgl%cg
#ifdef RESISTIVE
! DIFFUSION FULL STEP
         if (associated(custom_emf_bnd)) call custom_emf_bnd(wcu%arr)
         cg%b%arr(ib1,:,:,:) = cg%b%arr(ib1,:,:,:) - wcu%arr*cg%idl(dim1)
         wcu%arr = pshift(wcu%arr,dim1)
         cg%b%arr(ib1,:,:,:) = cg%b%arr(ib1,:,:,:) + wcu%arr*cg%idl(dim1)
         wcu%arr = mshift(wcu%arr,dim1)
         cg%b%arr(ib2,:,:,:) = cg%b%arr(ib2,:,:,:) + wcu%arr*cg%idl(dim2)
         wcu%arr = pshift(wcu%arr,dim2)
         cg%b%arr(ib2,:,:,:) = cg%b%arr(ib2,:,:,:) - wcu%arr*cg%idl(dim2)
#endif /* RESISTIVE */
! ADVECTION FULL STEP
         if (associated(custom_emf_bnd)) call custom_emf_bnd(cg%wa%arr)
         cg%b%arr(ib1,:,:,:) = cg%b%arr(ib1,:,:,:) - cg%wa%arr*cg%idl(dim1)
         cg%wa%arr = mshift(cg%wa%arr,dim1)
         cg%b%arr(ib1,:,:,:) = cg%b%arr(ib1,:,:,:) + cg%wa%arr*cg%idl(dim1)
         cg%b%arr(ib2,:,:,:) = cg%b%arr(ib2,:,:,:) - cg%wa%arr*cg%idl(dim2)
         cg%wa%arr = pshift(cg%wa%arr,dim2)
         cg%b%arr(ib2,:,:,:) = cg%b%arr(ib2,:,:,:) + cg%wa%arr*cg%idl(dim2)
         cgl => cgl%nxt
      enddo

      call all_mag_boundaries

   end subroutine mag_add
#endif /* MAGNETIC */
!------------------------------------------------------------------------------------------

end module fluidupdate
