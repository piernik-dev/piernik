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
!! \brief A library that converts a 3D logical map into a set of boxes
!!
!! \details We want to produce as small set of non overlapping boxes that covers a given map.
!! The map can be e.g a list of grid regions that is covered by intra-level communication. Then the boxes are the part that need to be prolonged from the coarse grid.
!! Finding best (optimal) solution can be computationally expensive, so we  are satisfied by "good enough" solutions.
!! In our approach we first look for all corners that are convex as some those for sure should belong to the biggest boxes that fit into the map.
!! Then we arbitrarily decide that the most important are the boxes that have corners near Center of Mass, which often, but not always tends to be true.
!! Next, look for a biggest box with that corner by "growing" it like a grain of salt.
!! Iterate until all flags get covered by boxes.
!!
!! Alternative approaches:
!! * Find all the boxes possible at any time (note that even if a box connect a set of corners, another box can be found for the other corners),
!! then pick up the biggest one and iterate. This is much more costly, especially for complicated maps with a lot of corners (like plotted by concentric spheres).
!! * Coarsen the map as long as there exist at least one coarse cell that is fully covered by flags. Pic it as a fine box and try to grow. Iterate as in other approaches.
!! Coarsening may be done in directionally-split fashion for better flexibility
!!
!! \todo Use it to merge smaller grids into larger ones.
!<

module mergebox

   use box_list,    only: box_list_T
   use corner_list, only: corner_list_T
   use constants,   only: ndims

   implicit none

   private
   public :: wmap

   type :: wmap
      logical, allocatable, dimension(:,:,:) :: map !< Logical map to be processed
      type(corner_list_T) :: clist                  !< List of convec corners
      type(box_list_T) :: blist                     !< List of boxes
      real, dimension(ndims) :: CoM                 !< Center of Mass
   contains
      procedure :: init          !< Initialize
      procedure :: print         !< Print the map
      procedure :: make_clist    !< Create convect corner list
      procedure :: grow_cuboids  !< Pick up which box to grow
      procedure :: cleanup       !< Free the memory
      procedure :: calcCoM       !< Compute Center of Mass
      procedure :: grow_box      !< Grow a single box
      procedure :: update_map    !< Clear part of the map that is covered by the boxes already found
      procedure :: find_boxes    !< Find a list of boxes
      procedure :: set           !< Set flags inside a cuboid
      procedure :: clear         !< Clear flags inside a cuboid
   end type wmap

contains

!> \brief Initialize: Allocate map with 1 guardcell in each side and mark the map as .false.

   subroutine init(this, bnd)

      use constants,  only: xdim, ydim, zdim, LO, HI
      use dataio_pub, only: die

      implicit none

      class(wmap),                                  intent(inout) :: this
      integer(kind=8), dimension(xdim:zdim, LO:HI), intent(in)    :: bnd

      if (allocated(this%map)) call die("[mergeboxes:init] Already allocated")
      allocate(this%map(bnd(xdim, LO)-1:bnd(xdim, HI)+1, bnd(ydim, LO)-1:bnd(ydim, HI)+1, bnd(zdim, LO)-1:bnd(zdim, HI)+1))
      allocate(this%clist%clist(0))
      allocate(this%blist%blist(0))
      call this%clear

   end subroutine init

!> \brief Set flags inside a cuboid or set the whole map (except for guardcells)

   subroutine set(this, bnd)

      use constants, only: xdim, ydim, zdim, LO, HI

      implicit none

      class(wmap),                                            intent(inout) :: this
      integer(kind=8), dimension(xdim:zdim, LO:HI), optional, intent(in)    :: bnd

      if (present(bnd)) then
         this%map(bnd(xdim, LO):bnd(xdim, HI), bnd(ydim, LO):bnd(ydim, HI), bnd(zdim, LO):bnd(zdim, HI)) = .true.
      else ! set everything except for guardcells
         this%map(lbound(this%map, dim=1)+1:ubound(this%map, dim=1)-1, &
              &   lbound(this%map, dim=2)+1:ubound(this%map, dim=2)-1, &
              &   lbound(this%map, dim=3)+1:ubound(this%map, dim=3)-1) = .true.
      endif

   end subroutine set

!> \brief Clear flags inside a cuboid or clear the whole map

   subroutine clear(this, bnd)

      use constants, only: xdim, ydim, zdim, LO, HI

      implicit none

      class(wmap),                                            intent(inout) :: this
      integer(kind=8), dimension(xdim:zdim, LO:HI), optional, intent(in)    :: bnd

      if (present(bnd)) then
         this%map(bnd(xdim, LO):bnd(xdim, HI), bnd(ydim, LO):bnd(ydim, HI), bnd(zdim, LO):bnd(zdim, HI)) = .false.
      else
         this%map = .false.
      endif

   end subroutine clear

!> \brief Print the map (useful for debugging)

   subroutine print(this)

      use dataio_pub, only: msg, printinfo

      implicit none

      class(wmap), intent(inout) :: this

      integer(kind=8) :: i, j, k
      real :: a

      do k = lbound(this%map, dim=3),ubound(this%map, dim=3)
         do j = lbound(this%map, dim=2),ubound(this%map, dim=2)
            do i = lbound(this%map, dim=1),ubound(this%map, dim=1)
               a = 0.
               if (this%map(i,j,k)) a = 1.
               write(msg,'(a,3i3,f5.2)')"[mergeboxes:print] ",i,j,k,a
               call printinfo(msg)
            enddo
         enddo
      enddo

   end subroutine print

!> \brief Create convect corner list and set their vectors

   subroutine make_clist(this)

      use corner_list, only: corner
      use constants,   only: xdim, ydim, zdim
      use dataio_pub,  only: die

      implicit none

      class(wmap), intent(inout) :: this

      integer(kind=8) :: i, j, k
      type(corner) :: c

      call this%clist%cleanup
      allocate(this%clist%clist(0))

      do k = lbound(this%map, dim=3)+1,ubound(this%map, dim=3)-1
         do j = lbound(this%map, dim=2)+1,ubound(this%map, dim=2)-1
            do i = lbound(this%map, dim=1)+1,ubound(this%map, dim=1)-1
               if (this%map(i, j ,k)) then
                  if (.not. (this%map(i-1, j ,  k  ) .and. this%map(i+1, j,   k  ))) then
                     if (.not. (this%map(i,   j-1, k  ) .and. this%map(i,   j+1, k  ))) then
                        if (.not. (this%map(i,   j,   k-1) .and. this%map(i,   j,   k+1))) then
                           c%pos = [ i, j, k ]
                           c%ivec = [ 0, 0, 0 ]
                           c%dist = huge(1.)
                           if (this%map(i+1, j, k)) c%ivec(xdim) =  1
                           if (this%map(i-1, j, k)) c%ivec(xdim) = -1
                           if (this%map(i, j+1, k)) c%ivec(ydim) =  1
                           if (this%map(i, j-1, k)) c%ivec(ydim) = -1
                           if (this%map(i, j, k+1)) c%ivec(zdim) =  1
                           if (this%map(i, j, k-1)) c%ivec(zdim) = -1
                           call this%clist%add_c(c)
                        endif
                     endif
                  endif
               endif
            enddo
         enddo
      enddo

      if (count(this%map) > 0 .and. size(this%clist%clist) == 0) call die("[mergeboxes:make_clist] no corners detected on an nonempty map")

   end subroutine make_clist

!> \brief Pick up which box to grow

   subroutine grow_cuboids(this)

      use box_list,  only: box
      use constants, only: LO, HI

      implicit none

      class(wmap), intent(inout) :: this

      integer :: i, ci
      real :: closest
      integer(kind=8), dimension(ndims) :: b
      type(box) :: bx

      ci = 0
      closest = huge(1.)
      do i = lbound(this%clist%clist, dim=1), ubound(this%clist%clist, dim=1)
         if (this%clist%clist(i)%dist < closest) then
            closest = this%clist%clist(i)%dist
            ci = i
         endif
      enddo

      call this%grow_box(ci, b)
      bx%b(:, LO) = this%clist%clist(ci)%pos
      bx%b(:, HI) = b
      call this%blist%add_b(bx)

   end subroutine grow_cuboids

!> \brief Free the memory

   subroutine cleanup(this)

      implicit none

      class(wmap), intent(inout) :: this

      if (allocated(this%map)) deallocate(this%map)
      call this%clist%cleanup
      call this%blist%cleanup

   end subroutine cleanup

!> \brief Compute Center of Mass

   subroutine calcCoM(this)

      use constants, only: xdim, zdim, LONG

      implicit none

      class(wmap), intent(inout) :: this

      integer, parameter :: MASS=xdim-1
      real, dimension(MASS:zdim) :: mx
      integer(kind=8) :: i, j, k

      mx = 0
      do k = lbound(this%map, dim=3)+1,ubound(this%map, dim=3)-1
         do j = lbound(this%map, dim=2)+1,ubound(this%map, dim=2)-1
            do i = lbound(this%map, dim=1)+1,ubound(this%map, dim=1)-1
               if (this%map(i,j,k)) then
                  mx = mx + [ 1_LONG, i, j, k ]
               endif
            enddo
         enddo
      enddo

      if (mx(MASS) /= 0) then
         this%CoM(xdim:zdim) = mx(xdim:zdim) / mx(MASS)
      else
         mx = -huge(1.)
      endif

   end subroutine calcCoM

!> \brief Grow a single box

   subroutine grow_box(this, ci, bx)

      use box_list,  only: box
      use constants, only: LO, HI, xdim, ydim, zdim

      implicit none

      class(wmap),                       intent(inout) :: this
      integer,                           intent(in)    :: ci
      integer(kind=8), dimension(ndims), intent(out)   :: bx

      integer(kind=8), dimension(ndims) :: vec
      type(box) :: b
      integer :: d
      integer(kind=8) :: i, j, k
      logical :: flag

      associate ( pos => this%clist%clist(ci)%pos)
      bx = pos
      vec = this%clist%clist(ci)%ivec

      do while (any(vec /= 0))
         do d = xdim, zdim
            if (vec(d) /= 0) then
               b%b(:, LO) = pos
               b%b(:, HI) = bx
               call b%sanitize
               if (vec(d) > 0) then
                  b%b(d, LO) = b%b(d, HI) + 1
                  b%b(d, HI) = b%b(d, LO)
               else
                  b%b(d, HI) = b%b(d, LO) - 1
                  b%b(d, LO) = b%b(d, HI)
               endif
               flag = .true.
               do k = b%b(zdim, LO), b%b(zdim, HI)
                  do j = b%b(ydim, LO), b%b(ydim, HI)
                     do i = b%b(xdim, LO), b%b(xdim, HI)
                        flag = flag .and. this%map(i, j, k)
                        if (.not. flag) exit
                     enddo
                     if (.not. flag) exit
                  enddo
                  if (.not. flag) exit
               enddo
               if (flag) then
                  bx(d) = bx(d)+vec(d)
               else
                  vec(d) = 0
               endif
            endif
         enddo
      enddo

      end associate

   end subroutine grow_box

!> \brief Clear part of the map that is covered by the boxes already found

   subroutine update_map(this)

      implicit none

      class(wmap), intent(inout), target :: this

      call this%clear(this%blist%blist(ubound(this%blist%blist, dim=1))%b)

   end subroutine update_map

!> \brief Find a list of boxes: main iteration loop

   subroutine find_boxes(this)

      implicit none

      class(wmap), intent(inout), target :: this

!      call this%print
      do while (count(this%map) > 0)
         call this%make_clist
         call this%calcCoM
         call this%clist%update_dists(this%CoM)
!         call this%clist%print
         call this%grow_cuboids
         call this%update_map
      enddo

!      call this%blist%print

   end subroutine find_boxes

end module mergebox

!!$program mergeboxes
!!$
!!$   use mergebox, only: wmap
!!$
!!$   implicit none
!!$
!!$   integer, parameter :: cubesize = 8, smallcube = cubesize/2
!!$   type(wmap) :: cube
!!$
!!$   integer :: i, j, k
!!$
!!$   call cube%init([cubesize, cubesize, cubesize])
!!$
!!$   do k = 1, cubesize
!!$      do j = 1, cubesize
!!$         do i = 1, cubesize
!!$            cube%map(i, j, k) = (sum(( [i, j, k] - (cubesize+1.)/2.)**2) < (cubesize/2.)**2)
!!$         enddo
!!$      enddo
!!$   enddo
!!$
!!$!   cube%map = .true.
!!$!   cube%map(2:smallcube, 2:smallcube, 2:smallcube) = .false.
!!$
!!$!   cube%map = .false.
!!$!   cube%map(smallcube, smallcube, smallcube) = .true.
!!$
!!$   call cube%find_boxes
!!$
!!$   call cube%cleanup
!!$
!!$end program mergeboxes
