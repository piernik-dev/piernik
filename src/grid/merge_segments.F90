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

!> \brief Module that merges individual segments for internal boundaries into large clumps to optimize MPI communication on AMR grids

module merge_segments

   use sort_segment_list, only: sort_segment_list_T

   implicit none

   private
   public :: merge_segments_T

   type :: merge_segments_T
      type(sort_segment_list_T), dimension(:, :), allocatable :: sl ! segment list (FIRST:LAST, IN:OUT)
      logical :: valid
   contains
      procedure          :: merge     ! Merge segments
      procedure, private :: populate  ! Initialize segment list with fresh data
   end type merge_segments_T

   enum, bind(C)
      enumerator :: IN, OUT
   end enum

contains

!> \brief Merge segments

   subroutine merge(this, list)

      use cg_list,  only: cg_list_T
      use mpisetup, only: FIRST, LAST !, proc

      implicit none

      class(merge_segments_T), intent(inout) :: this
      class(cg_list_T),        intent(in)    :: list

      integer :: p

      if (allocated(this%sl)) deallocate(this%sl) !> \todo check if it properly frees this%sl(:)%list
      allocate(this%sl(FIRST:LAST, IN:OUT))

      call this%populate(list)
      do p = FIRST, LAST
         ! technically we don't need to aggregate on p==proc, but it is safer to do it anyway
         call this%sl(p, IN )%sort
         call this%sl(p, OUT)%sort
         !> \todo OPT: we can avoid one call to sort if we put both incoming and outgoing segments in
         !! this%sl(p)%list(:). Note that in this%populate we cannot assume that o_bnd will be sorted in the same
         !! way as i_bnd was. Then we'll be able to drop the IN|OUT index as well.
      enddo

      !this%valid = .true.

   end subroutine merge

!> \brief Initialize segment list with fresh data

   subroutine populate(this, list)

      use cg_list,           only: cg_list_T, cg_list_element
      use constants,         only: xdim, cor_dim

      implicit none

      class(merge_segments_T), intent(inout) :: this
      class(cg_list_T),        intent(in)    :: list

      type(cg_list_element), pointer :: cgl
      integer :: d, i

      cgl => list%first
      do while (associated(cgl))
         do d = xdim, cor_dim
            if (allocated(cgl%cg%i_bnd(d)%seg)) then
               do i = lbound(cgl%cg%i_bnd(d)%seg, dim=1), ubound(cgl%cg%i_bnd(d)%seg, dim=1)
                  call this%sl(cgl%cg%i_bnd(d)%seg(i)%proc, IN)%add(cgl%cg%i_bnd(d)%seg(i)%tag, cgl%cg%i_bnd(d)%seg(i)%se)
               enddo
            endif

            if (allocated(cgl%cg%o_bnd(d)%seg)) then
               do i = lbound(cgl%cg%o_bnd(d)%seg, dim=1), ubound(cgl%cg%o_bnd(d)%seg, dim=1)
                  call this%sl(cgl%cg%o_bnd(d)%seg(i)%proc, OUT)%add(cgl%cg%o_bnd(d)%seg(i)%tag, cgl%cg%o_bnd(d)%seg(i)%se)
               enddo
            endif
         enddo
         cgl => cgl%nxt
      enddo

   end subroutine populate

end module merge_segments
