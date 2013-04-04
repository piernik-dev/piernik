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
!! \brief This module handles refinement criteria that are provided by simple geometric shapes.
!!
!! \todo Implement empty spheres (shells), add radius to points.
!<

module refinement_primitives

   implicit none

   private
   public :: mark_all_primitives

contains

!> \brief Call routines related to all primitive shapes

   subroutine mark_all_primitives

      implicit none

      call mark_points
      call mark_boxes

   end subroutine mark_all_primitives

!> \brief Mark all refinement points

   subroutine mark_points

      use cg_level_base,      only: base
      use cg_level_connected, only: cg_level_connected_T
      use cg_level_finest,    only: finest
      use cg_list,            only: cg_list_element
      use constants,          only: LO, HI, ndims, xdim, ydim, zdim
      use domain,             only: dom
      use refinement,         only: refine_points

      implicit none

      type(cg_level_connected_T), pointer :: curl
      type(cg_list_element), pointer :: cgl
      integer :: ip
      integer, dimension(ndims) :: ijk

      do ip = lbound(refine_points, dim=1), ubound(refine_points, dim=1)

         curl => finest%level
         do while (associated(curl) .and. curl%level_id>=base%level%level_id)
            if (curl%level_id <= refine_points(ip)%level) then
               cgl => curl%first
               do while (associated(cgl))
                  if (all((cgl%cg%fbnd(:, LO)<=refine_points(ip)%coords(:) .and. cgl%cg%fbnd(:, HI)>=refine_points(ip)%coords(:)) .or. .not. dom%has_dir(:))) then
                     if (curl%level_id < refine_points(ip)%level) then
                        ijk = cgl%cg%ijkse(:, LO)
                        where (dom%has_dir) ijk = cgl%cg%ijkse(:, LO) + int((refine_points(ip)%coords(:) - cgl%cg%fbnd(:, LO)) / cgl%cg%dl(:))
                        where (ijk > cgl%cg%ijkse(:, HI)) ijk = cgl%cg%ijkse(:, HI)
                        if (cgl%cg%leafmap(ijk(xdim), ijk(ydim), ijk(zdim))) cgl%cg%refine_flags%refine = .true.
                     endif
                     cgl%cg%refine_flags%derefine = .false.
                  endif
                  cgl => cgl%nxt
               enddo
            endif
            curl => curl%coarser
         enddo

      enddo

   end subroutine mark_points

!> \brief Mark all refinement boxes

   subroutine mark_boxes

      use cg_level_base,      only: base
      use cg_level_connected, only: cg_level_connected_T
      use cg_level_finest,    only: finest
      use cg_list,            only: cg_list_element
      use constants,          only: LO, HI, ndims, xdim, ydim, zdim
      use domain,             only: dom
      use refinement,         only: refine_boxes

      implicit none

      type(cg_level_connected_T), pointer :: curl
      type(cg_list_element), pointer :: cgl
      integer :: ip
      integer, dimension(ndims) :: ijk_l, ijk_h

      do ip = lbound(refine_boxes, dim=1), ubound(refine_boxes, dim=1)

         curl => finest%level
         do while (associated(curl) .and. curl%level_id>=base%level%level_id)
            if (curl%level_id <= refine_boxes(ip)%level) then
               cgl => curl%first
               do while (associated(cgl))
                  if (all((cgl%cg%fbnd(:, LO)<=refine_boxes(ip)%hcoords(:) .and. cgl%cg%fbnd(:, HI)>=refine_boxes(ip)%lcoords(:)) .or. .not. dom%has_dir(:))) then
                     if (curl%level_id < refine_boxes(ip)%level) then
                        ijk_l = cgl%cg%ijkse(:, LO)
                        ijk_h = cgl%cg%ijkse(:, HI)
                        where (dom%has_dir)
                           ijk_l = cgl%cg%ijkse(:, LO) + int((refine_boxes(ip)%lcoords(:) - cgl%cg%fbnd(:, LO)) / cgl%cg%dl(:))
                           ijk_h = cgl%cg%ijkse(:, LO) + int((refine_boxes(ip)%hcoords(:) - cgl%cg%fbnd(:, LO)) / cgl%cg%dl(:))
                        endwhere
                        where (ijk_l > cgl%cg%ijkse(:, HI)) ijk_l = cgl%cg%ijkse(:, HI)
                        where (ijk_h > cgl%cg%ijkse(:, HI)) ijk_h = cgl%cg%ijkse(:, HI)
                        where (ijk_l < cgl%cg%ijkse(:, LO)) ijk_l = cgl%cg%ijkse(:, LO)
                        where (ijk_h < cgl%cg%ijkse(:, LO)) ijk_h = cgl%cg%ijkse(:, LO)
                        if (any(cgl%cg%leafmap(ijk_l(xdim):ijk_h(xdim), ijk_l(ydim):ijk_h(ydim), ijk_l(zdim):ijk_h(zdim)))) cgl%cg%refine_flags%refine = .true.
                     endif
                     cgl%cg%refine_flags%derefine = .false.
                  endif
                  cgl => cgl%nxt
               enddo
            endif
            curl => curl%coarser
         enddo

      enddo

   end subroutine mark_boxes

end module refinement_primitives
