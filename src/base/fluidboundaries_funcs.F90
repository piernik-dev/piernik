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

module fluidboundaries_funcs
! pulled by ANY

   implicit none

   private
   public :: bnd_null, bnd_xl_per, bnd_xr_per, bnd_xl_ref, bnd_xr_ref, bnd_xl_out, bnd_xr_out, bnd_xl_outd, bnd_xr_outd, &
        &    user_bnd_xl, user_bnd_xr, user_bnd_yl, user_bnd_yr, user_bnd_zl, user_bnd_zr, func_bnd_xl, func_bnd_xr, &
        &    init_default_fluidboundaries

   interface

      subroutine user_bnd(cg)

         use grid_cont, only: grid_container

         implicit none

         type(grid_container), pointer, intent(inout) :: cg

      end subroutine user_bnd

   end interface

   procedure(user_bnd), pointer :: user_bnd_xl, user_bnd_xr, user_bnd_yl, user_bnd_yr, user_bnd_zl, user_bnd_zr
   procedure(user_bnd), pointer :: func_bnd_xl, func_bnd_xr

contains

!--------------------------------------------------------------------------------------------------
   subroutine default_bnd(cg)

      use grid_cont,  only: grid_container
      use dataio_pub, only: die

      implicit none

      type(grid_container), pointer, intent(inout) :: cg

      call die("User boundaries are not defined")

      if (.true. .or. cg%grid_id >=0) return ! suppress compiler warnings

   end subroutine default_bnd
!--------------------------------------------------------------------------------------------------
   subroutine init_default_fluidboundaries

      use constants,  only: PIERNIK_INIT_MPI
      use dataio_pub, only: code_progress, die
#ifdef VERBOSE
      use dataio_pub, only: printinfo
#endif /* VERBOSE */

      implicit none

      if (code_progress < PIERNIK_INIT_MPI) call die("[fluidboundaries_funcs:init_default_fluidboundaries] MPI not initialized.")

#ifdef VERBOSE
      call printinfo("[fluidboundaries_funcs:init_default_fluidboundaries]: commencing...")
#endif /* VERBOSE */

      user_bnd_xl => default_bnd
      user_bnd_xr => default_bnd
      user_bnd_yl => default_bnd
      user_bnd_yr => default_bnd
      user_bnd_zl => default_bnd
      user_bnd_zr => default_bnd

#ifdef VERBOSE
      call printinfo("[fluidboundaries_funcs:init_default_fluidboundaries]: finished. \o/")
#endif /* VERBOSE */

   end subroutine init_default_fluidboundaries

   subroutine bnd_null(cg)

      use grid_cont, only: grid_container

      implicit none

      type(grid_container), pointer, intent(inout) :: cg

      if (.true. .or. cg%grid_id >= 0) return ! suppress compiler warnings

   end subroutine bnd_null

   subroutine bnd_xl_per(cg)

      use domain,    only: dom
      use grid_cont, only: grid_container

      implicit none

      type(grid_container), pointer, intent(inout) :: cg

      cg%u(:,1:dom%nb,:,:) = cg%u(:, cg%ieb:cg%ie,:,:)

   end subroutine bnd_xl_per

   subroutine bnd_xl_ref(cg)

      use domain,     only: dom
      use grid_cont,  only: grid_container
      use fluidindex, only: iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz
#ifndef ISO
      use fluidindex, only: iarr_all_en
#endif /* !ISO */
#ifdef COSM_RAYS
      use fluidindex, only: iarr_all_crs
#endif /* COSM_RAYS */

      implicit none

      type(grid_container), pointer, intent(inout) :: cg

      integer :: ib

      do ib=1, dom%nb
         cg%u([iarr_all_dn,iarr_all_my,iarr_all_mz], cg%is-ib,:,:)  = cg%u([iarr_all_dn,iarr_all_my,iarr_all_mz], dom%nb+ib,:,:)
         cg%u(iarr_all_mx, cg%is-ib,:,:)  = -cg%u(iarr_all_mx, dom%nb+ib,:,:)
#ifndef ISO
         cg%u(iarr_all_en, cg%is-ib,:,:)  =  cg%u(iarr_all_en, dom%nb+ib,:,:)
#endif /* !ISO */
#ifdef COSM_RAYS
         cg%u(iarr_all_crs, cg%is-ib,:,:) =  cg%u(iarr_all_crs, dom%nb+ib,:,:)
#endif /* COSM_RAYS */
      enddo

   end subroutine bnd_xl_ref

   subroutine bnd_xl_out(cg)

      use domain,         only: dom
      use grid_cont,      only: grid_container
#ifdef COSM_RAYS
      use initcosmicrays, only: smallecr
      use fluidindex,     only: iarr_all_crs
#endif /* COSM_RAYS */

      implicit none

      type(grid_container), pointer, intent(inout) :: cg

      integer :: ib

      do ib = 1, dom%nb
         cg%u(:,ib,:,:)            = cg%u(:, cg%is,:,:)
#ifdef COSM_RAYS
         cg%u(iarr_all_crs,ib,:,:) = smallecr
#endif /* COSM_RAYS */
      enddo

   end subroutine bnd_xl_out

   subroutine bnd_xl_outd(cg)

      use domain,         only: dom
      use grid_cont,      only: grid_container
      use fluidindex,     only: iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz
#ifndef ISO
      use fluidindex,     only: iarr_all_en
#endif /* !ISO */
#ifdef COSM_RAYS
      use initcosmicrays, only: smallecr
      use fluidindex,     only: iarr_all_crs
#endif /* COSM_RAYS */

      implicit none

      type(grid_container), pointer, intent(inout) :: cg

      integer :: ib

      do ib = 1, dom%nb
         cg%u([iarr_all_dn,iarr_all_my,iarr_all_mz],ib,:,:)  = cg%u([iarr_all_dn,iarr_all_my,iarr_all_mz], cg%is,:,:)
         cg%u(iarr_all_mx,ib,:,:)  = min(cg%u(iarr_all_mx, cg%is,:,:),0.0)
#ifndef ISO
         cg%u(iarr_all_en,ib,:,:)  = cg%u(iarr_all_en, cg%is,:,:)
#endif /* !ISO */
#ifdef COSM_RAYS
         cg%u(iarr_all_crs,ib,:,:) = smallecr
#endif /* COSM_RAYS */
      enddo

   end subroutine bnd_xl_outd

   subroutine bnd_xr_per(cg)

      use constants, only: xdim
      use grid_cont, only: grid_container

      implicit none

      type(grid_container), pointer, intent(inout) :: cg

      cg%u(:, cg%ie+1:cg%n_(xdim),:,:) = cg%u(:, cg%is:cg%isb,:,:)

   end subroutine bnd_xr_per

   subroutine bnd_xr_ref(cg)

      use domain,     only: dom
      use grid_cont,  only: grid_container
      use fluidindex, only: iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz
#ifndef ISO
      use fluidindex, only: iarr_all_en
#endif /* !ISO */
#ifdef COSM_RAYS
      use fluidindex, only: iarr_all_crs
#endif /* COSM_RAYS */

      implicit none

      type(grid_container), pointer, intent(inout) :: cg

      integer :: ib

      do ib=1, dom%nb
         cg%u([iarr_all_dn,iarr_all_my,iarr_all_mz], cg%ie+ib,:,:) = cg%u([iarr_all_dn,iarr_all_my,iarr_all_mz], cg%ie+1-ib,:,:)
         cg%u(iarr_all_mx, cg%ie+ib,:,:)  = -cg%u(iarr_all_mx, cg%ie+1-ib,:,:)
#ifndef ISO
         cg%u(iarr_all_en, cg%ie+ib,:,:)  =  cg%u(iarr_all_en, cg%ie+1-ib,:,:)
#endif /* !ISO */
#ifdef COSM_RAYS
         cg%u(iarr_all_crs, cg%ie+ib,:,:) =  cg%u(iarr_all_crs, cg%ie+1-ib,:,:)
#endif /* COSM_RAYS */
      enddo

   end subroutine bnd_xr_ref

   subroutine bnd_xr_out(cg)

      use domain,         only: dom
      use grid_cont,      only: grid_container
#ifdef COSM_RAYS
      use initcosmicrays, only: smallecr
      use fluidindex,     only: iarr_all_crs
#endif /* COSM_RAYS */

      implicit none

      type(grid_container), pointer, intent(inout) :: cg

      integer :: ib

      do ib = 1, dom%nb
         cg%u(:, cg%ie+ib,:,:)            = cg%u(:, cg%ie,:,:)
#ifdef COSM_RAYS
         cg%u(iarr_all_crs, cg%ie+ib,:,:) = smallecr
#endif /* COSM_RAYS */
      enddo

   end subroutine bnd_xr_out

   subroutine bnd_xr_outd(cg)

      use domain,         only: dom
      use grid_cont,      only: grid_container
      use fluidindex,     only: iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz
#ifndef ISO
      use fluidindex,     only: iarr_all_en
#endif /* !ISO */
#ifdef COSM_RAYS
      use initcosmicrays, only: smallecr
      use fluidindex,     only: iarr_all_crs
#endif /* COSM_RAYS */

      implicit none

      type(grid_container), pointer, intent(inout) :: cg

      integer :: ib

      do ib = 1, dom%nb
         cg%u([iarr_all_dn,iarr_all_my,iarr_all_mz], cg%ie+ib,:,:) = cg%u([iarr_all_dn,iarr_all_my,iarr_all_mz], cg%ie,:,:)
         cg%u(iarr_all_mx, cg%ie+ib,:,:)  = max(cg%u(iarr_all_mx, cg%ie,:,:),0.0)
#ifndef ISO
         cg%u(iarr_all_en, cg%ie+ib,:,:)  = cg%u(iarr_all_en, cg%ie,:,:)
#endif /* !ISO */
#ifdef COSM_RAYS
         cg%u(iarr_all_crs, cg%ie+ib,:,:) = smallecr
#endif /* COSM_RAYS */
      enddo

   end subroutine bnd_xr_outd

end module fluidboundaries_funcs
