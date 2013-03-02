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

!>
!! \brief Module of boundary conditions for fluids
!<
module fluidboundaries
! pulled by ANY

   implicit none

   private
   public :: bnd_u, all_fluid_boundaries

contains

   subroutine init_fluidboundaries(cg)

      use constants,             only: PIERNIK_INIT_DOMAIN, xdim, zdim, LO, HI, &
           &                           BND_MPI, BND_FC, BND_MPI_FC, BND_PER, BND_REF, BND_OUT, BND_OUTD, BND_OUTH, BND_OUTHD, BND_COR, BND_SHE, BND_USER
      use dataio_pub,            only: msg, warn, die, code_progress
      use domain,                only: is_multicg
      use grid_cont,             only: grid_container

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer(kind=4)                           :: dir, side

      if (code_progress < PIERNIK_INIT_DOMAIN) call die("[fluidboundaries:init_fluidboundaries] MPI not initialized.") ! bnd_xl, bnd_xr

      do dir = xdim, zdim
         do side = LO, HI

            select case (cg%bnd(dir, side))
               case (BND_MPI, BND_REF, BND_OUT, BND_OUTD, BND_USER, BND_PER)
                  ! Do nothing
               case (BND_FC, BND_MPI_FC)
                  call die("[fluidboundaries:init_fluidboundaries] fine-coarse interfaces not implemented yet")
               case (BND_COR)
                  if (dir == zdim) then
                     write(msg,'("[fluidboundaries:init_fluidboundaries] corner ",i1," boundary condition ",i3," not implemented in ",i1,"-direction")') side, cg%bnd(dir, side), dir
                     call warn(msg)
                  endif
               case (BND_SHE)
                  if (dir /= xdim) then
                     write(msg,'("[fluidboundaries:init_fluidboundaries] shear ",i1," boundary condition ",i3," not implemented in ",i1,"-direction")') side, cg%bnd(dir, side), dir
                     call warn(msg)
                  endif
               case (BND_OUTH)
                  if (dir == zdim) then
                     if (is_multicg) call die("[fluidboundaries:init_fluidboundaries] hydrostatic:outh_bnd with multiple grid pieces per processor not implemented yet") !nontrivial not really checked
                  else
                     write(msg,'("[fluidboundaries:init_fluidboundaries] outflow hydrostatic ",i1," boundary condition ",i3," not implemented in ",i1,"-direction")') side, cg%bnd(dir, side), dir
                     call warn(msg)
                  endif
               case (BND_OUTHD)
                  if (dir == zdim) then
                     if (is_multicg) call die("[fluidboundaries:init_fluidboundaries] hydrostatic:outh_bnd with multiple grid pieces per processor not implemented yet") !nontrivial not really checked
                  else
                     write(msg,'("[fluidboundaries:init_fluidboundaries] outflow hydrostatic ",i1," boundary condition ",i3," not implemented in ",i1,"-direction")') side, cg%bnd(dir, side), dir
                     call warn(msg)
                  endif
               case default
                  write(msg,'("[fluidboundaries:init_fluidboundaries] unknown ",i1," boundary condition ",i3," not implemented in ",i1,"-direction")') side, cg%bnd(dir, side), dir
                  call warn(msg)
            end select
         enddo
      enddo


   end subroutine init_fluidboundaries

   subroutine bnd_u(dir, cg)

      use constants,             only: ndims, xdim, ydim, zdim, LO, HI, INT4, &
           &                           BND_MPI, BND_FC, BND_MPI_FC, BND_PER, BND_REF, BND_OUT, BND_OUTD, BND_COR, BND_SHE, BND_USER
      use dataio_pub,            only: msg, warn, die
      use domain,                only: dom
      use fluidboundaries_funcs, only: user_fluidbnd
      use fluidindex,            only: iarr_all_dn
      use grid_cont,             only: grid_container
      use named_array_list,      only: wna
#ifdef COSM_RAYS
      use initcosmicrays,        only: smallecr
      use fluidindex,            only: iarr_all_crs
#endif /* COSM_RAYS */
#ifdef GRAV
      use constants,             only: BND_OUTH, BND_OUTHD
      use hydrostatic,           only: outh_bnd
#endif /* GRAV */

      implicit none

      integer(kind=4),               intent(in)    :: dir
      type(grid_container), pointer, intent(inout) :: cg

      integer(kind=4), dimension(ndims,LO:HI)      :: l, r
      logical, save                                :: frun = .true.
      integer(kind=4)                              :: side, ssign, ib

      if (.not. any([xdim, ydim, zdim] == dir)) call die("[fluidboundaries:bnd_u] Invalid direction.")

      if (frun) then
         call init_fluidboundaries(cg)
         frun = .false.
      endif

!===============================================================

! Non-MPI boundary conditions
      l = cg%lhn ; r = l
      do side = LO, HI

         select case (cg%bnd(dir, side))
         case (BND_MPI, BND_COR, BND_SHE, BND_FC, BND_MPI_FC, BND_PER)
            ! Do nothing
         case (BND_USER)
            call user_fluidbnd(dir, side, cg, wn=wna%fi)
         case (BND_REF)
            ssign = 2_INT4*side-3_INT4
            do ib=1_INT4, dom%nb
               l(dir,:) = cg%ijkse(dir,side)+ssign*ib ; r(dir,:) = cg%ijkse(dir,side)+ssign*(1_INT4-ib)
               cg%u(:,l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) = cg%u(:,r(xdim,LO):r(xdim,HI),r(ydim,LO):r(ydim,HI),r(zdim,LO):r(zdim,HI))
               cg%u(iarr_all_dn+dir, l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) = -cg%u(iarr_all_dn+dir,l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI))
            enddo
         case (BND_OUT)
            r(dir,:) = cg%ijkse(dir,side)
            ssign = 2_INT4*side-3_INT4
            do ib=1_INT4, dom%nb
               l(dir,:) = cg%ijkse(dir,side)+ssign*ib
               cg%u(:,l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) = cg%u(:,r(xdim,LO):r(xdim,HI),r(ydim,LO):r(ydim,HI),r(zdim,LO):r(zdim,HI))
#ifdef COSM_RAYS
               cg%u(iarr_all_crs,l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) = smallecr
#endif /* COSM_RAYS */
            enddo
         case (BND_OUTD)
            r(dir,:) = cg%ijkse(dir,side)
            ssign = 2_INT4*side-3_INT4
            do ib=1_INT4, dom%nb
               l(dir,:) = cg%ijkse(dir,side)+ssign*ib
               cg%u(:,l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) = cg%u(:,r(xdim,LO):r(xdim,HI),r(ydim,LO):r(ydim,HI),r(zdim,LO):r(zdim,HI))
!> \deprecated BEWARE: use of uninitialized value on first call (a side effect of r1726)
#ifdef COSM_RAYS
               cg%u(iarr_all_crs,l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) = smallecr
#endif /* COSM_RAYS */
            enddo
            l(dir,:) = [1_INT4, dom%nb] + cg%ijkse(dir,side)*(side-1_INT4)
            if (side == LO) then
               cg%u(iarr_all_dn+dir,l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) = min(cg%u(iarr_all_dn+dir,l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)),0.0)
            else
               cg%u(iarr_all_dn+dir,l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) = max(cg%u(iarr_all_dn+dir,l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)),0.0)
            endif
#ifdef GRAV
         case (BND_OUTH)
            if (dir == zdim) call outh_bnd(side, cg, .false.)
         case (BND_OUTHD)
            if (dir == zdim) call outh_bnd(side, cg, .true.)
#endif /* GRAV */
         case default
            write(msg,'("[fluidboundaries:bnd_u]: Unrecognized ",i1," boundary condition ",i3," not implemented in ",i1,"-direction")') side, cg%bnd(dir, side), dir
            call warn(msg)
         end select

      enddo

   end subroutine bnd_u

   subroutine all_fluid_boundaries

      use cg_level_base,      only: base
      use cg_level_connected, only: cg_level_connected_T
!      use cg_level_finest,    only: finest
      use cg_list,            only: cg_list_element
      use constants,          only: xdim, zdim
      use domain,             only: dom
      use named_array_list,   only: wna

      implicit none

      type(cg_level_connected_T), pointer :: curl
      type(cg_list_element), pointer :: cgl
      integer(kind=4)                :: dir

      curl => base%level

!      call finest%level%restrict_to_base

      ! should be more selective (modified leaves?)
      do while (associated(curl))
         call curl%arr4d_boundaries(wna%fi)
         cgl => curl%first
         do while (associated(cgl))
            do dir = xdim, zdim
               if (dom%has_dir(dir)) call bnd_u(dir, cgl%cg)
            enddo
            cgl => cgl%nxt
         enddo
         curl => curl%finer
      enddo

   end subroutine all_fluid_boundaries

end module fluidboundaries
