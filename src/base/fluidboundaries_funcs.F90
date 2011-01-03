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
   public :: &
      & bnd_null, bnd_xl_per, bnd_xr_per, bnd_xl_ref, bnd_xr_ref, bnd_xl_out, bnd_xr_out, bnd_xl_outd, bnd_xr_outd

   contains

      subroutine bnd_null
         implicit none
         return
      end subroutine bnd_null

      subroutine bnd_xl_per

         use arrays, only: u
         use grid,   only: cg

         implicit none

         u(:,1:cg%nb,:,:) = u(:, cg%nxb+1:cg%ie,:,:)

      end subroutine bnd_xl_per

      subroutine bnd_xl_ref

         use arrays,         only: u
         use grid,           only: cg
         use fluidindex,     only: iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz
#ifndef ISO
         use fluidindex,     only: iarr_all_en
#endif /* !ISO */
#ifdef COSM_RAYS
         use fluidindex,     only: iarr_all_crs
#endif /* COSM_RAYS */

         implicit none

         integer :: ib

         do ib=1, cg%nb
            u([iarr_all_dn,iarr_all_my,iarr_all_mz], cg%is-ib,:,:)  = u([iarr_all_dn,iarr_all_my,iarr_all_mz], cg%nb+ib,:,:)
            u(iarr_all_mx, cg%is-ib,:,:)  = -u(iarr_all_mx, cg%nb+ib,:,:)
#ifndef ISO
            u(iarr_all_en, cg%is-ib,:,:)  =  u(iarr_all_en, cg%nb+ib,:,:)
#endif /* !ISO */
#ifdef COSM_RAYS
            u(iarr_all_crs, cg%is-ib,:,:) =  u(iarr_all_crs, cg%nb+ib,:,:)
#endif /* COSM_RAYS */
         enddo

      end subroutine bnd_xl_ref

      subroutine bnd_xl_out

         use arrays,         only: u
         use grid,           only: cg
#ifdef COSM_RAYS
         use initcosmicrays, only: smallecr
         use fluidindex,     only: iarr_all_crs
#endif /* COSM_RAYS */

         implicit none

         integer :: ib

         do ib = 1, cg%nb
            u(:,ib,:,:)            = u(:, cg%is,:,:)
#ifdef COSM_RAYS
            u(iarr_all_crs,ib,:,:) = smallecr
#endif /* COSM_RAYS */
         enddo

      end subroutine bnd_xl_out

      subroutine bnd_xl_outd

         use arrays,         only: u
         use grid,           only: cg
         use fluidindex,     only: iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz
#ifndef ISO
         use fluidindex,     only: iarr_all_en
#endif /* !ISO */
#ifdef COSM_RAYS
         use initcosmicrays, only: smallecr
         use fluidindex,     only: iarr_all_crs
#endif /* COSM_RAYS */

         implicit none

         integer :: ib

         do ib = 1, cg%nb
            u([iarr_all_dn,iarr_all_my,iarr_all_mz],ib,:,:)  = u([iarr_all_dn,iarr_all_my,iarr_all_mz], cg%is,:,:)
            u(iarr_all_mx,ib,:,:)  = min(u(iarr_all_mx, cg%is,:,:),0.0)
#ifndef ISO
            u(iarr_all_en,ib,:,:)  = u(iarr_all_en, cg%is,:,:)
#endif /* !ISO */
#ifdef COSM_RAYS
            u(iarr_all_crs,ib,:,:) = smallecr
#endif /* COSM_RAYS */
         enddo

      end subroutine bnd_xl_outd

      subroutine bnd_xr_per

         use arrays, only: u
         use grid,   only: cg

         implicit none

         u(:, cg%ie+1:cg%nx,:,:) = u(:, cg%is:cg%isb,:,:)

      end subroutine bnd_xr_per

      subroutine bnd_xr_ref

         use arrays,         only: u
         use grid,           only: cg
         use fluidindex,     only: iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz
#ifndef ISO
         use fluidindex,     only: iarr_all_en
#endif /* !ISO */
#ifdef COSM_RAYS
         use fluidindex,     only: iarr_all_crs
#endif /* COSM_RAYS */

         implicit none

         integer :: ib

         do ib=1, cg%nb
            u([iarr_all_dn,iarr_all_my,iarr_all_mz], cg%ie+ib,:,:) = u([iarr_all_dn,iarr_all_my,iarr_all_mz], cg%ie+1-ib,:,:)
            u(iarr_all_mx, cg%ie+ib,:,:)  = -u(iarr_all_mx, cg%ie+1-ib,:,:)
#ifndef ISO
            u(iarr_all_en, cg%ie+ib,:,:)  =  u(iarr_all_en, cg%ie+1-ib,:,:)
#endif /* !ISO */
#ifdef COSM_RAYS
            u(iarr_all_crs, cg%ie+ib,:,:) =  u(iarr_all_crs, cg%ie+1-ib,:,:)
#endif /* COSM_RAYS */
         enddo

      end subroutine bnd_xr_ref

      subroutine bnd_xr_out

         use arrays,         only: u
         use grid,           only: cg
#ifdef COSM_RAYS
         use initcosmicrays, only: smallecr
         use fluidindex,     only: iarr_all_crs
#endif /* COSM_RAYS */

         implicit none

         integer :: ib

         do ib = 1, cg%nb
            u(:, cg%ie+ib,:,:)            = u(:, cg%ie,:,:)
#ifdef COSM_RAYS
            u(iarr_all_crs, cg%ie+ib,:,:) = smallecr
#endif /* COSM_RAYS */
         enddo

      end subroutine bnd_xr_out

      subroutine bnd_xr_outd

         use arrays,         only: u
         use grid,           only: cg
         use fluidindex,     only: iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz
#ifndef ISO
         use fluidindex,     only: iarr_all_en
#endif /* !ISO */
#ifdef COSM_RAYS
         use initcosmicrays, only: smallecr
         use fluidindex,     only: iarr_all_crs
#endif /* COSM_RAYS */

         implicit none

         integer :: ib

         do ib = 1, cg%nb
            u([iarr_all_dn,iarr_all_my,iarr_all_mz], cg%ie+ib,:,:) = u([iarr_all_dn,iarr_all_my,iarr_all_mz], cg%ie,:,:)
            u(iarr_all_mx, cg%ie+ib,:,:)  = max(u(iarr_all_mx, cg%ie,:,:),0.0)
#ifndef ISO
            u(iarr_all_en, cg%ie+ib,:,:)  = u(iarr_all_en, cg%ie,:,:)
#endif /* !ISO */
#ifdef COSM_RAYS
            u(iarr_all_crs, cg%ie+ib,:,:) = smallecr
#endif /* COSM_RAYS */
         enddo

      end subroutine bnd_xr_outd

end module fluidboundaries_funcs
