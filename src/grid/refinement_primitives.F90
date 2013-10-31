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
      integer(kind=8), dimension(ndims) :: ip_ijk

      do ip = lbound(refine_points, dim=1), ubound(refine_points, dim=1)

         curl => finest%level
         do while (associated(curl) .and. curl%level_id>=base%level%level_id)
            if (curl%level_id <= refine_points(ip)%level) then

               ip_ijk(:) = curl%off(:)
               where (dom%has_dir) ip_ijk(:) = curl%off(:) + floor((refine_points(ip)%coords(:) - dom%edge(:, LO))/dom%L_(:)*curl%n_d)
               !BEWARE: ip_ijk can contain indices outside the domain

               cgl => curl%first
               do while (associated(cgl))
                  cgl%cg%refine_flags%derefine = .true.
                  if (all(ip_ijk >= cgl%cg%ijkse(:, LO)) .and. all(ip_ijk <= cgl%cg%ijkse(:, HI))) then
                     if (curl%level_id < refine_points(ip)%level) then
                        if (cgl%cg%leafmap(ip_ijk(xdim), ip_ijk(ydim), ip_ijk(zdim))) cgl%cg%refine_flags%refine = .true.
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
      integer(kind=8), dimension(ndims, LO:HI) :: ip_ijk, cg_ijk  ! cell coordinates of given box on current level

      do ip = lbound(refine_boxes, dim=1), ubound(refine_boxes, dim=1)

         curl => finest%level
         do while (associated(curl) .and. curl%level_id>=base%level%level_id)

            if (curl%level_id <= refine_boxes(ip)%level) then

               ip_ijk(:, LO) = curl%off(:)
               ip_ijk(:, HI) = curl%off(:)
               where (dom%has_dir)
                  ip_ijk(:, LO) = curl%off(:) + floor((refine_boxes(ip)%coords(:, LO) - dom%edge(:, LO))/dom%L_(:)*curl%n_d)
                  ip_ijk(:, HI) = curl%off(:) + floor((refine_boxes(ip)%coords(:, HI) - dom%edge(:, LO))/dom%L_(:)*curl%n_d)
               end where
               !BEWARE: ip_ijk can contain indices outside the domain

               cgl => curl%first
               do while (associated(cgl))
                  cgl%cg%refine_flags%derefine = .true.
                  if (all(ip_ijk(:, HI) >= cgl%cg%ijkse(:, LO)) .and. all(ip_ijk(:, LO) <= cgl%cg%ijkse(:, HI))) then
                     if (curl%level_id < refine_boxes(ip)%level) then
                        cg_ijk(:, LO) = min(max(int(ip_ijk(:, LO), kind=4), cgl%cg%ijkse(:, LO)), cgl%cg%ijkse(:, HI))
                        cg_ijk(:, HI) = min(max(int(ip_ijk(:, HI), kind=4), cgl%cg%ijkse(:, LO)), cgl%cg%ijkse(:, HI))
                        if (any(cgl%cg%leafmap(cg_ijk(xdim, LO):cg_ijk(xdim, HI), &
                             &                 cg_ijk(ydim, LO):cg_ijk(ydim, HI), &
                             &                 cg_ijk(zdim, LO):cg_ijk(zdim, HI)))) cgl%cg%refine_flags%refine = .true.
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
