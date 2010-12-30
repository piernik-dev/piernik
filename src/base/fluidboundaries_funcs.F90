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
         use grid,   only: nb, nxb

         implicit none

         u(:,1:nb,:,:) = u(:,nxb+1:nxb+nb,:,:)

      end subroutine bnd_xl_per

      subroutine bnd_xl_ref

         use arrays,         only: u
         use grid,           only: nb
         use fluidindex,     only: iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz
#ifndef ISO
         use fluidindex,     only: iarr_all_en
#endif /* !ISO */
#ifdef COSM_RAYS
         use fluidindex,     only: iarr_all_crs
#endif /* COSM_RAYS */

         implicit none

         integer :: ib

         do ib=1,nb
            u([iarr_all_dn,iarr_all_my,iarr_all_mz],nb+1-ib,:,:)  = u([iarr_all_dn,iarr_all_my,iarr_all_mz],nb+ib,:,:)
            u(iarr_all_mx,nb+1-ib,:,:)  = -u(iarr_all_mx,nb+ib,:,:)
#ifndef ISO
            u(iarr_all_en,nb+1-ib,:,:)  =  u(iarr_all_en,nb+ib,:,:)
#endif /* !ISO */
#ifdef COSM_RAYS
            u(iarr_all_crs,nb+1-ib,:,:) =  u(iarr_all_crs,nb+ib,:,:)
#endif /* COSM_RAYS */
         enddo

      end subroutine bnd_xl_ref

      subroutine bnd_xl_out

         use arrays,         only: u
         use grid,           only: nb
#ifdef COSM_RAYS
         use initcosmicrays, only: smallecr
         use fluidindex,     only: iarr_all_crs
#endif /* COSM_RAYS */

         implicit none

         integer :: ib

         do ib = 1, nb
            u(:,ib,:,:)            = u(:,nb+1,:,:)
#ifdef COSM_RAYS
            u(iarr_all_crs,ib,:,:) = smallecr
#endif /* COSM_RAYS */
         enddo

      end subroutine bnd_xl_out

      subroutine bnd_xl_outd

         use arrays,         only: u
         use grid,           only: nb
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

         do ib = 1, nb
            u([iarr_all_dn,iarr_all_my,iarr_all_mz],ib,:,:)  = u([iarr_all_dn,iarr_all_my,iarr_all_mz],nb+1,:,:)
            u(iarr_all_mx,ib,:,:)  = min(u(iarr_all_mx,nb+1,:,:),0.0)
#ifndef ISO
            u(iarr_all_en,ib,:,:)  = u(iarr_all_en,nb+1,:,:)
#endif /* !ISO */
#ifdef COSM_RAYS
            u(iarr_all_crs,ib,:,:) = smallecr
#endif /* COSM_RAYS */
         enddo

      end subroutine bnd_xl_outd

      subroutine bnd_xr_per

         use arrays, only: u
         use grid,   only: nb, nxb

         implicit none

         u(:,nxb+nb+1:nxb+2*nb,:,:) = u(:,nb+1:2*nb,:,:)

      end subroutine bnd_xr_per

      subroutine bnd_xr_ref

         use arrays,         only: u
         use grid,           only: nb, nxb
         use fluidindex,     only: iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz
#ifndef ISO
         use fluidindex,     only: iarr_all_en
#endif /* !ISO */
#ifdef COSM_RAYS
         use fluidindex,     only: iarr_all_crs
#endif /* COSM_RAYS */

         implicit none

         integer :: ib

         do ib=1,nb
            u([iarr_all_dn,iarr_all_my,iarr_all_mz],nb+nxb+ib,:,:) = u([iarr_all_dn,iarr_all_my,iarr_all_mz],nb+nxb+1-ib,:,:)
            u(iarr_all_mx,nb+nxb+ib,:,:)  = -u(iarr_all_mx,nb+nxb+1-ib,:,:)
#ifndef ISO
            u(iarr_all_en,nb+nxb+ib,:,:)  =  u(iarr_all_en,nb+nxb+1-ib,:,:)
#endif /* !ISO */
#ifdef COSM_RAYS
            u(iarr_all_crs,nb+nxb+ib,:,:) =  u(iarr_all_crs,nb+nxb+1-ib,:,:)
#endif /* COSM_RAYS */
         enddo

      end subroutine bnd_xr_ref

      subroutine bnd_xr_out

         use arrays,         only: u
         use grid,           only: nb, nxb
#ifdef COSM_RAYS
         use initcosmicrays, only: smallecr
         use fluidindex,     only: iarr_all_crs
#endif /* COSM_RAYS */

         implicit none

         integer :: ib

         do ib = 1, nb
            u(:,nb+nxb+ib,:,:)            = u(:,nb+nxb,:,:)
#ifdef COSM_RAYS
            u(iarr_all_crs,nb+nxb+ib,:,:) = smallecr
#endif /* COSM_RAYS */
         enddo

      end subroutine bnd_xr_out

      subroutine bnd_xr_outd

         use arrays,         only: u
         use grid,           only: nb, nxb
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

         do ib = 1, nb
            u([iarr_all_dn,iarr_all_my,iarr_all_mz],nb+nxb+ib,:,:) = u([iarr_all_dn,iarr_all_my,iarr_all_mz],nb+nxb,:,:)
            u(iarr_all_mx,nb+nxb+ib,:,:)  = max(u(iarr_all_mx,nb+nxb,:,:),0.0)
#ifndef ISO
            u(iarr_all_en,nb+nxb+ib,:,:)  = u(iarr_all_en,nb+nxb,:,:)
#endif /* !ISO */
#ifdef COSM_RAYS
            u(iarr_all_crs,nb+nxb+ib,:,:) = smallecr
#endif /* COSM_RAYS */
         enddo

      end subroutine bnd_xr_outd

end module fluidboundaries_funcs
