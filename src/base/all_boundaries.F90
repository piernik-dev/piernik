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
!! \brief Global update of boundary routines
!!
!! \details moved to a separate file due to concerns about circular dependencies
!<
module all_boundaries

   implicit none

   private
   public :: all_bnd, all_fluid_boundaries
#ifdef MAGNETIC
   public :: all_mag_boundaries
#endif /* MAGNETIC */

contains

!>
!! Subroutine calling all type boundaries after initialization of new run or restart reading
!! \todo make sure that all_fluid_boundaries and all_mag_boundaries can handle BND_USER boundaries right now, or do the boundaries later
!<

   subroutine all_bnd

      implicit none

!      if (all(cg%bnd(:,:) /= BND_USER)) then
      call all_fluid_boundaries
#ifdef MAGNETIC
      call all_mag_boundaries
#endif /* MAGNETIC */
!      endif

   end subroutine all_bnd

   subroutine all_fluid_boundaries(dir)

      use cg_leaves,          only: leaves
!      use cg_level_finest,    only: finest
      use constants,          only: xdim, zdim
      use domain,             only: dom
      use named_array_list,   only: wna

      implicit none

      integer(kind=4), optional, intent(in) :: dir !< select only this direction

      integer(kind=4)                     :: d

      if (present(dir)) then
         if (.not. dom%has_dir(dir)) return
      endif

!      call finest%level%restrict_to_base

      ! should be more selective (modified leaves?)
      call leaves%leaf_arr4d_boundaries(wna%fi, dir=dir)
      if (present(dir)) then
         call leaves%bnd_u(dir)
      else
         do d = xdim, zdim
            if (dom%has_dir(d)) call leaves%bnd_u(d)
         enddo
      endif

   end subroutine all_fluid_boundaries

#ifdef MAGNETIC
   subroutine all_mag_boundaries

      use cg_leaves,        only: leaves
      use cg_list_global,   only: all_cg
      use constants,        only: xdim, zdim
      use domain,           only: dom
      use named_array_list, only: wna

      implicit none

      integer(kind=4) :: dir

      do dir = xdim, zdim
         if (dom%has_dir(dir)) call all_cg%internal_boundaries_4d(wna%bi, dir=dir) ! should be more selective (modified leaves?)
      enddo

      ! Do not fuse these loops
      do dir = xdim, zdim
         if (dom%has_dir(dir)) call leaves%bnd_b(dir)
      enddo

   end subroutine all_mag_boundaries
#endif /* MAGNETIC */

end module all_boundaries
