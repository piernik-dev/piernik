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
      type(sort_segment_list_T), dimension(:, :), allocatable :: sl ! array of sortable segment lists (FIRST:LAST, IN:OUT)
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

      integer :: p, i

      if (allocated(this%sl)) deallocate(this%sl) !> \todo check if it properly frees this%sl(:)%list
      allocate(this%sl(FIRST:LAST, IN:OUT))

      this%valid = .true.
      call this%populate(list)
      if (this%valid) then
         do p = FIRST, LAST
            do i = IN, OUT
               ! technically we don't need to aggregate on p==proc, but it is safer to do it anyway
               call this%sl(p, i)%sort
               !> \todo OPT: we can avoid one call to sort if we put both incoming and outgoing segments in
               !! this%sl(p)%list(:). Note that in this%populate we cannot assume that o_bnd will be sorted in the same
               !! way as i_bnd was. Then we'll be able to drop the IN|OUT index as well.
               call this%sl(p, i)%find_offsets
            enddo
         enddo
      else
         deallocate(this%sl)
      endif

   end subroutine merge

!> \brief Initialize segment list with fresh data

   subroutine populate(this, list)

      use cg_list,    only: cg_list_T, cg_list_element
      use constants,  only: xdim, cor_dim
      use dataio_pub, only: warn
      use mpisetup,   only: proc

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
                  call this%sl(cgl%cg%i_bnd(d)%seg(i)%proc, IN)%add( &
                       &       cgl%cg%i_bnd(d)%seg(i)%tag, &
                       &       cgl%cg%i_bnd(d)%seg(i)%se, cgl%cg, d)
                  if (cgl%cg%i_bnd(d)%seg(i)%proc == proc .and. .not. associated(cgl%cg%i_bnd(d)%seg(i)%local)) then
                     this%valid = .false.
                     call warn("[merge_segments:populate] local i_bnd without pointer set. Cannot safely use merged messages.") ! Or perhaps it is better to die here?
                  endif
               enddo
            endif

            if (allocated(cgl%cg%o_bnd(d)%seg)) then
               do i = lbound(cgl%cg%o_bnd(d)%seg, dim=1), ubound(cgl%cg%o_bnd(d)%seg, dim=1)
                  call this%sl(cgl%cg%o_bnd(d)%seg(i)%proc, OUT)%add( &
                       &       cgl%cg%o_bnd(d)%seg(i)%tag, &
                       &       cgl%cg%o_bnd(d)%seg(i)%se, cgl%cg, d)
               enddo
            endif
         enddo
         cgl => cgl%nxt
      enddo

   end subroutine populate

end module merge_segments
