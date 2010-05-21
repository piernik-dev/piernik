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
!    Initial implemetation of PIERNIK code was based on TVD split MHD code by
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

   integer, parameter :: DIR_X = 1, DIR_Y = DIR_X + 1, DIR_Z = DIR_Y + 1

contains

   subroutine fluid_update

      use timer,         only : timer_
      use dataio,        only : check_log, check_tsl
      use timestep,      only : time_step
      use mpisetup,      only : proc, dt, dtm, t, nstep
      use dataio_public, only : halfstep
#ifdef SN_SRC
      use snsources,     only : random_sn
#endif /* SN_SRC */
#ifdef SNE_DISTR
      use sndistr,       only : supernovae_distribution
#endif /* SNE_DISTR */
#ifdef DEBUG
      use dataio_public, only : nhdf
      use dataio,        only : write_hdf
#endif /* DEBUG */

      implicit none

      logical, save :: first_run = .true.
      real          :: ts   ! Timestep wallclock
#ifdef DEBUG
      integer       :: system, syslog
#endif /* DEBUG */

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

      if(proc.eq.0) write(*,900) nstep,dt,t,ts
900   format('   nstep = ',i7,'   dt = ',es22.16,'   t = ',es22.16,'   dWallClock = ',f7.2,' s')

      t=t+dt

      call make_3sweeps(.true.) ! X -> Y -> Z

! Sources ----------------------------------------

#ifdef SN_SRC
#ifndef SNE_DISTR
      call random_sn
#endif /* SNE_DISTR */
#endif /* SN_SRC */

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      t=t+dt
      dtm = dt
      halfstep = .true.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call make_3sweeps(.false.) ! Z -> Y -> X

#ifdef SNE_DISTR
      call supernovae_distribution
#ifdef DEBUG
      syslog = system('echo -n sne_dis')
      call write_hdf
      nhdf = nhdf + 1
#endif /* DEBUG */
#endif /* SNE_DISTR */

      if (first_run) first_run = .false.

   end subroutine fluid_update

!------------------------------------------------------------------------------------------
!
! Perform sweeps in all three directions plus sources that are calculated every timestep
!

   subroutine make_3sweeps(forward)

      use types,           only : problem_customize_solution
      use dataio_public,   only : skip_advection
#ifdef SHEAR
      use shear,           only : yshift
      use fluidboundaries, only : bnd_u
      use mpisetup,        only : t, dt
      use grid,            only : nxd, nyd, nzd
#endif /* SHEAR */
#ifdef GRAV
      use gravity,         only : source_terms_grav
#endif /* GRAV */

      implicit none

      logical, intent(in) :: forward  !< If .true. then do X->Y->Z sweeps, if .false. then reverse that order

      integer :: s

#ifdef SHEAR
      call yshift(t, dt)
      if(nxd /= 1) call bnd_u('xdim')
      if(nyd /= 1) call bnd_u('ydim')
#endif /* SHEAR */

#ifdef GRAV
      call source_terms_grav
#endif /* GRAV */

      if (.not. skip_advection) then
         if (forward) then
            do s = DIR_X, DIR_Z
               call make_sweep(s, forward)
            end do
         else
            do s = DIR_Z, DIR_X, -1
               call make_sweep(s, forward)
            end do
         end if
         if (associated(problem_customize_solution)) call problem_customize_solution
      end if

   end subroutine make_3sweeps

!------------------------------------------------------------------------------------------
!
! Perform single sweep in forward or backward direction
!

   subroutine make_sweep(dir, forward)

      use sweeps,        only : sweepx, sweepy, sweepz
      use grid,          only : nxd, nyd, nzd
      use errh,          only : die
#if defined SHEAR && defined FLUID_INTERACTIONS
      use sweeps,        only : source_terms_y
#endif /* SHEAR */
#ifdef COSM_RAYS
      use crdiffusion,   only : cr_diff_x, cr_diff_y, cr_diff_z
#endif /* COSM_RAYS */
#ifdef DEBUG
      use dataio_public, only : nhdf
      use dataio,        only : write_hdf
#endif /* DEBUG */

      implicit none

      integer, intent(in) :: dir      !< direction, one of DIR_X, DIR_Y, DIR_Z
      logical, intent(in) :: forward  !< if .false. then reverse operation order in the sweep

      select case (dir)

         case(DIR_X)
            if (nxd /= 1) then
               if (.not. forward) then
#ifdef COSM_RAYS
                  call cr_diff_x
#endif /* COSM_RAYS */
#ifdef MAGNETIC
                  call magfieldbyzx
#endif /* MAGNETIC */
               end if

               call sweepx

               if (forward) then
#ifdef MAGNETIC
                  call magfieldbyzx
#endif /* MAGNETIC */
#ifdef COSM_RAYS
                  call cr_diff_x
#endif /* COSM_RAYS */
               end if
            endif

         case(DIR_Y)
            if (nyd /= 1) then
               if (.not. forward) then
#ifdef COSM_RAYS
                  call cr_diff_y
#endif /* COSM_RAYS */
#ifdef MAGNETIC
                  call magfieldbzxy
#endif /* MAGNETIC */
               end if
               call sweepy

               if (forward) then
#ifdef MAGNETIC
                  call magfieldbzxy
#endif /* MAGNETIC */
#ifdef COSM_RAYS
                  call cr_diff_y
#endif /* COSM_RAYS */
               end if
            else
#if defined SHEAR && defined FLUID_INTERACTIONS
               call source_terms_y
#endif /* SHEAR */
            endif

         case(DIR_Z)
            if (nzd /= 1) then
               if (.not. forward) then
#ifdef COSM_RAYS
                  call cr_diff_z
#endif /* COSM_RAYS */
#ifdef MAGNETIC
                  call magfieldbxyz
#endif /* MAGNETIC */
               end if

               call sweepz

               if (forward) then
#ifdef MAGNETIC
                  call magfieldbxyz
#endif /* MAGNETIC */
#ifdef COSM_RAYS
                  call cr_diff_z
#endif  /* COSM_RAYS */
               end if
            endif

         case default
            write(*,'(a,i10)')"[fluidupdate:make_sweep] Illegal direction ",dir
            call die("[fluidupdate:make_sweep] Illegal direction.")

      end select

#ifdef DEBUG
      ! syslog = system('echo -n sweep z')
      call write_hdf
      nhdf = nhdf + 1 !\todo should go inside write_hdf
#endif /* DEBUG */

   end subroutine make_sweep

#ifdef MAGNETIC
   subroutine magfieldbyzx

      use fluidindex,  only : ibx,iby,ibz
      use arrays,      only : b
      use grid,        only : xdim,ydim,zdim,nyd,nzd
      use advects,     only : advectby_x,advectbz_x
#ifdef RESISTIVE
      use resistivity, only : diffuseby_x,diffusebz_x
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

      use fluidindex,  only : ibx,iby,ibz
      use arrays,      only : b
      use grid,        only : xdim,ydim,zdim,nzd,nxd
      use advects,     only : advectbx_y,advectbz_y
#ifdef RESISTIVE
      use resistivity, only : diffusebx_y,diffusebz_y
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

      use fluidindex,  only : ibx,iby,ibz
      use arrays,      only : b
      use grid,        only : xdim,ydim,zdim,nxd,nyd
      use advects,     only : advectbx_z,advectby_z
#ifdef RESISTIVE
      use resistivity, only : diffusebx_z,diffuseby_z
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

      use func,          only : pshift, mshift
      use arrays,        only : b, wa, wcu
      use grid,          only : dl
      use magboundaries, only : all_mag_boundaries

      implicit none

      integer             :: ib1,ib2,dim1,dim2

#ifdef RESISTIVE
! DIFFUSION FULL STEP

      b(ib1,:,:,:) = b(ib1,:,:,:) - wcu/dl(dim1)
!   wcu = cshift(wcu,shift= 1,dim=dim1)
      wcu = pshift(wcu,dim1)
      b(ib1,:,:,:) = b(ib1,:,:,:) + wcu/dl(dim1)
!   wcu = cshift(wcu,shift=-1,dim=dim1)
      wcu = mshift(wcu,dim1)
      b(ib2,:,:,:) = b(ib2,:,:,:) + wcu/dl(dim2)
!   wcu = cshift(wcu,shift= 1,dim=dim2)
      wcu = pshift(wcu,dim2)
      b(ib2,:,:,:) = b(ib2,:,:,:) - wcu/dl(dim2)
#endif /* RESISTIVE */
! ADVECTION FULL STEP

      b(ib1,:,:,:) = b(ib1,:,:,:) - wa/dl(dim1)
!   wa = cshift(wa,shift=-1,dim=dim1)
      wa = mshift(wa,dim1)
      b(ib1,:,:,:) = b(ib1,:,:,:) + wa/dl(dim1)
      b(ib2,:,:,:) = b(ib2,:,:,:) - wa/dl(dim2)
!   wa = cshift(wa,shift=1,dim=dim2)
      wa = pshift(wa,dim2)
      b(ib2,:,:,:) = b(ib2,:,:,:) + wa/dl(dim2)

      call all_mag_boundaries

   end subroutine mag_add
#endif /* MAGNETIC */
!------------------------------------------------------------------------------------------

end module fluidupdate
