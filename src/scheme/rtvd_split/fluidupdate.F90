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
#include "piernik.def"
module fluidupdate   ! SPLIT

   implicit none

   private
   public :: fluid_update

   integer, parameter :: DIR_X = 1, DIR_Y = DIR_X + 1, DIR_Z = DIR_Y + 1

contains

   subroutine fluid_update

      use dataio,        only: check_log, check_tsl
      use dataio_pub,    only: cwdlen, halfstep, msg, printinfo
      use mpisetup,      only: proc, dt, dtm, t, nstep
      use timer,         only: timer_
      use timestep,      only: time_step
#ifdef SN_SRC
      use snsources,     only: random_sn
#endif /* SN_SRC */

      implicit none
      character(len=cwdlen), parameter :: fmt900 = "('   nstep = ',i7,'   dt = ',es22.16,'   t = ',es22.16,'   dWallClock = ',f7.2,' s')"
      logical, save                    :: first_run = .true.
      real                             :: ts   ! Timestep wallclock

      halfstep = .false.

      if (first_run) then
         dtm = 0.0
         ts=timer_("fluid_update",.true.)
         ts = 0.0
      else
         dtm = dt
         ts=timer_("fluid_update")
      endif
      call time_step

#ifdef RESISTIVE
      if (first_run) then
         dtm = 0.0
         dt  = 0.0
      endif
#endif /* RESISTIVE */

      call check_log
      call check_tsl
      if (first_run .and. proc == 0) then
         write(msg, fmt900) 0,dt,t,ts
         call printinfo(msg, .true.)
      endif

      t=t+dt

      call make_3sweeps(.true.) ! X -> Y -> Z

! Sources ----------------------------------------

#ifdef SN_SRC
      call random_sn
#endif /* SN_SRC */

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      t=t+dt
      dtm = dt
      halfstep = .true.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call make_3sweeps(.false.) ! Z -> Y -> X

      if (first_run) first_run = .false.

      ts=timer_("fluid_update")
      if (proc == 0) then
         write(msg, fmt900) nstep,dt,t,ts
         call printinfo(msg, .true.)
      endif

   end subroutine fluid_update

!------------------------------------------------------------------------------------------
!
! Perform sweeps in all three directions plus sources that are calculated every timestep
!

   subroutine make_3sweeps(forward)

      use dataio_pub,      only: skip_advection
      use types,           only: problem_customize_solution
#ifdef SHEAR
      use fluidboundaries, only: bnd_u
      use grid,            only: nxd, nyd
      use mpisetup,        only: t, dt
      use shear,           only: yshift
#endif /* SHEAR */
#ifdef GRAV
      use gravity,         only: source_terms_grav
#endif /* GRAV */
#ifdef COSM_RAYS
      use initcosmicrays,     only: use_split
#ifdef MULTIGRID
      use multigrid_diffusion, only: multigrid_solve_diff
      use fluidboundaries,     only: all_fluid_boundaries
#endif /* MULTIGRID */
#endif /* COSM_RAYS */

      implicit none

      logical, intent(in) :: forward  !< If .true. then do X->Y->Z sweeps, if .false. then reverse that order

      integer :: s

#ifdef SHEAR
      if (nyd /= 1) call yshift(t, dt)
      if (nxd /= 1) call bnd_u('xdim')
      if (nyd /= 1) call bnd_u('ydim')
#endif /* SHEAR */

#ifdef GRAV
      call source_terms_grav
#endif /* GRAV */

#ifdef COSM_RAYS
#ifdef MULTIGRID
      if (.not. use_split) then
         call multigrid_solve_diff
         call all_fluid_boundaries
      endif
#endif /* MULTIGRID */
#endif /* COSM_RAYS */

      if (.not. skip_advection) then
         if (forward) then
            do s = DIR_X, DIR_Z
               call make_sweep(s, forward)
            enddo
         else
            do s = DIR_Z, DIR_X, -1
               call make_sweep(s, forward)
            enddo
         endif
         if (associated(problem_customize_solution)) call problem_customize_solution
      endif

   end subroutine make_3sweeps

!------------------------------------------------------------------------------------------
!
! Perform single sweep in forward or backward direction
!

   subroutine make_sweep(dir, forward)

      use dataio_pub,     only: msg, die
      use grid,           only: nxd, nyd, nzd
      use sweeps,         only: sweepx, sweepy, sweepz
#if defined SHEAR && defined FLUID_INTERACTIONS
      use sweeps,         only: source_terms_y
#endif /* SHEAR && FLUID_INTERACTIONS */
#ifdef COSM_RAYS
      use crdiffusion,    only: cr_diff_x, cr_diff_y, cr_diff_z
      use initcosmicrays, only: use_split
#endif /* COSM_RAYS */
#ifdef DEBUG
      use dataio_hdf5,    only: write_hdf5
      use dataio_pub,     only: chdf
#endif /* DEBUG */

      implicit none

      integer, intent(in) :: dir      !< direction, one of DIR_X, DIR_Y, DIR_Z
      logical, intent(in) :: forward  !< if .false. then reverse operation order in the sweep

      select case (dir)

         case (DIR_X)
            if (nxd /= 1) then
               if (.not. forward) then
#ifdef COSM_RAYS
                  if (use_split) call cr_diff_x
#endif /* COSM_RAYS */
#ifdef MAGNETIC
                  call magfieldbyzx
#endif /* MAGNETIC */
               endif

               call sweepx

               if (forward) then
#ifdef MAGNETIC
                  call magfieldbyzx
#endif /* MAGNETIC */
#ifdef COSM_RAYS
                  if (use_split) call cr_diff_x
#endif /* COSM_RAYS */
               endif
            endif

         case (DIR_Y)
            if (nyd /= 1) then
               if (.not. forward) then
#ifdef COSM_RAYS
                  if (use_split) call cr_diff_y
#endif /* COSM_RAYS */
#ifdef MAGNETIC
                  call magfieldbzxy
#endif /* MAGNETIC */
               endif
               call sweepy

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
#endif /* SHEAR */
            endif

         case (DIR_Z)
            if (nzd /= 1) then
               if (.not. forward) then
#ifdef COSM_RAYS
                  if (use_split) call cr_diff_z
#endif /* COSM_RAYS */
#ifdef MAGNETIC
                  call magfieldbxyz
#endif /* MAGNETIC */
               endif

               call sweepz

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

#ifdef DEBUG
      call write_hdf5(chdf)
#endif /* DEBUG */

   end subroutine make_sweep

#ifdef MAGNETIC
   subroutine magfieldbyzx

      use advects,     only: advectby_x, advectbz_x
      use arrays,      only: b
      use fluidindex,  only: ibx, iby, ibz
      use grid,        only: xdim, ydim, zdim
#ifdef RESISTIVE
      use resistivity, only: diffuseby_x, diffusebz_x
#endif /* RESISTIVE */

      implicit none

      call advectby_x

#ifdef RESISTIVE
      call diffuseby_x
#endif /* RESISTIVE */

      call mag_add(iby,xdim,ibx,ydim)

      call advectbz_x

#ifdef RESISTIVE
      call diffusebz_x
#endif /* RESISTIVE */

      call mag_add(ibz,xdim,ibx,zdim)

   end subroutine magfieldbyzx

!------------------------------------------------------------------------------------------

   subroutine magfieldbzxy

      use advects,     only: advectbx_y, advectbz_y
      use arrays,      only: b
      use fluidindex,  only: ibx, iby, ibz
      use grid,        only: xdim, ydim, zdim
#ifdef RESISTIVE
      use resistivity, only: diffusebx_y, diffusebz_y
#endif /* RESISTIVE */

      implicit none

      call advectbz_y

#ifdef RESISTIVE
      call diffusebz_y
#endif /* RESISTIVE */

      call mag_add(ibz,ydim,iby,zdim)

      call advectbx_y

#ifdef RESISTIVE
      call diffusebx_y
#endif /* RESISTIVE */

      call mag_add(ibx,ydim,iby,xdim)

   end subroutine magfieldbzxy

!------------------------------------------------------------------------------------------

   subroutine magfieldbxyz

      use advects,     only: advectbx_z, advectby_z
      use arrays,      only: b
      use fluidindex,  only: ibx, iby, ibz
      use grid,        only: xdim, ydim, zdim
#ifdef RESISTIVE
      use resistivity, only: diffusebx_z, diffuseby_z
#endif /* RESISTIVE */

      implicit none

      call advectbx_z
#ifdef RESISTIVE
      call diffusebx_z
#endif /* RESISTIVE */
      call mag_add(ibx,zdim,ibz,xdim)

      call advectby_z
#ifdef RESISTIVE
      call diffuseby_z
#endif /* RESISTIVE */

      call mag_add(iby,zdim,ibz,ydim)
   end subroutine magfieldbxyz

!------------------------------------------------------------------------------------------

   subroutine mag_add(ib1,dim1,ib2,dim2)

      use arrays,        only: b, wa
      use func,          only: pshift, mshift
      use grid,          only: idl
      use magboundaries, only: all_mag_boundaries
      use types,         only: custom_emf_bnd
#ifdef RESISTIVE
      use arrays,        only: wcu
#endif /* RESISTIVE */

      implicit none

      integer             :: ib1,ib2,dim1,dim2

#ifdef RESISTIVE
! DIFFUSION FULL STEP
      if (associated(custom_emf_bnd)) call custom_emf_bnd(wcu)
      b(ib1,:,:,:) = b(ib1,:,:,:) - wcu*idl(dim1)
      wcu = pshift(wcu,dim1)
      b(ib1,:,:,:) = b(ib1,:,:,:) + wcu*idl(dim1)
      wcu = mshift(wcu,dim1)
      b(ib2,:,:,:) = b(ib2,:,:,:) + wcu*idl(dim2)
      wcu = pshift(wcu,dim2)
      b(ib2,:,:,:) = b(ib2,:,:,:) - wcu*idl(dim2)
#endif /* RESISTIVE */
! ADVECTION FULL STEP
      if (associated(custom_emf_bnd)) call custom_emf_bnd(wa)
      b(ib1,:,:,:) = b(ib1,:,:,:) - wa*idl(dim1)
      wa = mshift(wa,dim1)
      b(ib1,:,:,:) = b(ib1,:,:,:) + wa*idl(dim1)
      b(ib2,:,:,:) = b(ib2,:,:,:) - wa*idl(dim2)
      wa = pshift(wa,dim2)
      b(ib2,:,:,:) = b(ib2,:,:,:) + wa*idl(dim2)

      call all_mag_boundaries

   end subroutine mag_add
#endif /* MAGNETIC */
!------------------------------------------------------------------------------------------

end module fluidupdate
