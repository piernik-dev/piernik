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

!> \brief This module contains grid container level list and related methods

module cg_list_lev

   use constants, only: ndims
   use gc_list,   only: cg_list
   use domain,    only: cuboids

   implicit none

   private
   public :: cg_list_level, cg_list_patch

   !>
   !! \brief A list of all cg of the same resolution.
   !!
   !! \details For positive refinement levels the list may be composed of several disconnected subsets of cg ("islands: made of one or more cg: cg_list_patch).
   !<
   type, extends(cg_list) :: cg_list_level
      integer :: lev                           !< level number (relative to base level). For printing, debug, and I/O use only. No arithmetic should depend on it.
      integer(kind=8), dimension(ndims) :: n_d !< maximum number of grid cells in each direction (size of fully occupied level)
      type(cuboids), dimension(:), allocatable :: pse  !< lists of grid chunks on each process (FIRST:LAST); Use with care, because this is an antiparallel thing
      type(cg_list_level), pointer :: coarser  !< coarser level cg set or null()
      type(cg_list_level), pointer :: finer    !< finer level cg set or null()
    contains
      procedure :: prolong_level0              !< interpolate the grid data to this%finer level
      procedure :: restrict_level              !< interpolate the grid data from this%coarser level
      procedure :: print_segments              !< print detailed information about current level decomposition
      ! fine-coarse boundary exchanges may also belong to this type
   end type cg_list_level

   !>
   !! \brief A list of grid containers that cover single box (or rectangle) on a certain resolution level
   !!
   !! \details This set would be a result of base domain or patch decomposition
   !<
   type, extends(cg_list) :: cg_list_patch
      integer(kind=4), dimension(ndims) :: n_d                !< number of grid cells
      integer(kind=8), dimension(ndims) :: off                !< offset (with respect to the base level, counted on own level)
      type(cg_list_patch), pointer :: parent                  !< Parent patch (or null()). \todo Consider relaxing this restriction and allow multi-parent patches
      type(cg_list_patch), dimension(:), pointer :: children  !< refined patches
      type(cg_list_level), pointer :: list_level              !< all cg on the same level
      !> \todo consider creating neigbour list (or (ndims, LO:HI) lists)
   end type cg_list_patch

contains

!>
!! \brief Simplest restriction (averaging).
!!
!! \todo implement high order restriction and test its influence on V-cycle convergence rate
!!
!! \details Some data can be locally copied without MPI, but this seems to have really little impact on the performance.
!! Some tests show that purely MPI code without local copies is marginally faster.
!<

   subroutine restrict_level(this, iv)

      use constants,  only: xdim, ydim, zdim, LO, HI, LONG, I_ONE
      use dataio_pub, only: msg, warn!, die
      use domain,     only: dom
      use gc_list,    only: cg_list_element
      use grid_cont,  only: grid_container
      use mpisetup,   only: proc, comm, ierr, req, status
      use mpi,        only: MPI_DOUBLE_PRECISION

      implicit none

      class(cg_list_level), target, intent(inout) :: this !< object invoking type-bound procedure
      integer,                      intent(in)    :: iv   !< variable to be restricted

      integer(kind=4), parameter :: tag1 = 1
      type(cg_list_level), pointer :: coarse
      integer :: g
      integer(kind=8), dimension(:,:), pointer :: fse, cse ! shortcuts for fine segment and coarse segment
      integer(kind=8) :: i, j, k, ic, jc, kc
      integer(kind=8), dimension(xdim:zdim) :: off1
      real :: norm
      integer(kind=4) :: nr
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg            !< current grid container

!      if (iv < lbound(this%first%cg%q, dim=1) .or. iv > ubound(this%first%cg%q, dim=1)) call die("[gc_list:restrict_level] Invalid variable index.")

      coarse => this%coarser
      if (.not. associated(coarse)) then
         write(msg,'(a,i3)')"[gc_list:restrict_level] no coarse level here: ", this%lev
         call warn(msg) ! can't restrict base level
         return
      endif

      nr = 0
      cgl => coarse%first
      do while (associated(cgl))
         cg => cgl%cg

         !OPT find a way to reduce this to areas with nonlocal incoming restriction
         cg%q(iv)%arr(:,:,:) = 0. ! disables check_dirty

         ! be ready to receive everything into right buffers
         if (allocated(cg%mg%f_tgt%seg)) then
            do g = 1, ubound(cg%mg%f_tgt%seg(:), dim=1)
               nr = nr + I_ONE
               call MPI_Irecv(cg%mg%f_tgt%seg(g)%buf(1, 1, 1), size(cg%mg%f_tgt%seg(g)%buf(:, :, :)), MPI_DOUBLE_PRECISION, cg%mg%f_tgt%seg(g)%proc, tag1, comm, req(nr), ierr)
            enddo
         endif

         cgl => cgl%nxt
      enddo

      ! interpolate to coarse buffer and send it
      norm = 1./2**dom%eff_dim
      cgl => this%first
      do while (associated(cgl))
         cg => cgl%cg
         do g = 1, ubound(cg%mg%c_tgt%seg(:), dim=1)

            fse => cg%mg%c_tgt%seg(g)%se

            cg%mg%c_tgt%seg(g)%buf(:, :, :) = 0.
            off1(:) = mod(cg%off(:)+fse(:, LO) - cg%ijkse(:, LO), 2_LONG)
            do k = fse(zdim, LO), fse(zdim, HI)
               kc = (k-fse(zdim, LO)+off1(zdim))/2 + 1
               do j = fse(ydim, LO), fse(ydim, HI)
                  jc = (j-fse(ydim, LO)+off1(ydim))/2 + 1
                  do i = fse(xdim, LO), fse(xdim, HI)
                     ic = (i-fse(xdim, LO)+off1(xdim))/2 + 1
                     cg%mg%c_tgt%seg(g)%buf(ic, jc, kc) = cg%mg%c_tgt%seg(g)%buf(ic, jc, kc) + cg%q(iv)%arr(i, j, k) * norm
                  enddo
               enddo
            enddo
            nr = nr + I_ONE
            call MPI_Isend(cg%mg%c_tgt%seg(g)%buf(1, 1, 1), size(cg%mg%c_tgt%seg(g)%buf(:, :, :)), MPI_DOUBLE_PRECISION, cg%mg%c_tgt%seg(g)%proc, tag1, comm, req(nr), ierr)
         enddo
         cgl => cgl%nxt
      enddo

      if (nr>0) call MPI_Waitall(nr, req(:nr), status(:,:nr), ierr)

      ! copy the received buffers to the right places
      cgl => coarse%first
      do while (associated(cgl))
         cg => cgl%cg
         if (allocated(cg%mg%f_tgt%seg)) then
            do g = 1, ubound(cg%mg%f_tgt%seg(:), dim=1)
               cse => cg%mg%f_tgt%seg(g)%se
               cg%q(iv)%arr     (cse(xdim, LO):cse(xdim, HI), cse(ydim, LO):cse(ydim, HI), cse(zdim, LO):cse(zdim, HI)) = &
                    cg%q(iv)%arr(cse(xdim, LO):cse(xdim, HI), cse(ydim, LO):cse(ydim, HI), cse(zdim, LO):cse(zdim, HI)) + cg%mg%f_tgt%seg(g)%buf(:, :, :)
            enddo
         endif
         cgl => cgl%nxt
      enddo

   end subroutine restrict_level

!> \brief 0th order prolongation : injection

   subroutine prolong_level0(this, iv)

      use constants,  only: xdim, ydim, zdim, LO, HI, LONG, I_ONE
      use dataio_pub, only: msg, warn!, die
      use gc_list,    only: cg_list_element
      use grid_cont,  only: grid_container
      use mpisetup,   only: proc, comm, ierr, req, status
      use mpi,        only: MPI_DOUBLE_PRECISION

      implicit none

      class(cg_list_level), target, intent(inout) :: this !< object invoking type-bound procedure
      integer,                      intent(in)    :: iv   !< variable to be prolonged

      integer(kind=4), parameter :: tag1 = 1
      type(cg_list_level), pointer :: fine
      integer :: g
      integer(kind=8), dimension(:,:), pointer :: fse, cse ! shortcuts for fine segment and coarse segment
      integer(kind=8) :: i, j, k, ic, jc, kc
      integer(kind=8), dimension(xdim:zdim) :: off1
      integer(kind=4) :: nr
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg            !< current grid container

!      if (iv < lbound(this%first%cg%q, dim=1) .or. iv > ubound(this%first%cg%q, dim=1)) call die("[gc_list:prolong_level0] Invalid variable index.")

      fine => this%finer
      if (.not. associated(fine)) then
         write(msg,'(a,i3)')"[gc_list:restrict_level] no fine level here: ", this%lev
         call warn(msg) ! can't prolong finest level
         return
      endif

      cgl => fine%first
      do while (associated(cgl))
         cg => cgl%cg

         ! OPT: try to remove or limit this  ~20% Ir, ~50% Dw_m
         fine%first%cg%q(iv)%arr(:,:,:) = 0. ! disables check_dirty

         ! be ready to receive everything into right buffers
         nr = 0
         if (allocated(cg%mg%c_tgt%seg)) then
            do g = 1, ubound(cg%mg%c_tgt%seg(:), dim=1)
               nr = nr + I_ONE
               call MPI_Irecv(cg%mg%c_tgt%seg(g)%buf(1, 1, 1), size(cg%mg%c_tgt%seg(g)%buf(:, :, :)), MPI_DOUBLE_PRECISION, cg%mg%c_tgt%seg(g)%proc, tag1, comm, req(nr), ierr)
            enddo
         endif
         cgl => cgl%nxt
      enddo

      ! send coarse data
      cgl => this%first
      do while (associated(cgl))
         cg => cgl%cg
         off1(:) = mod(fine%first%cg%off(:), 2_LONG) !!!
         do g = 1, ubound(cg%mg%f_tgt%seg(:), dim=1)

            cse => cg%mg%f_tgt%seg(g)%se

            nr = nr + I_ONE
            cg%mg%f_tgt%seg(g)%buf(:, :, :) = cg%q(iv)%arr(cse(xdim, LO):cse(xdim, HI), cse(ydim, LO):cse(ydim, HI), cse(zdim, LO):cse(zdim, HI))
            call MPI_Isend(cg%mg%f_tgt%seg(g)%buf(1, 1, 1), size(cg%mg%f_tgt%seg(g)%buf(:, :, :)), MPI_DOUBLE_PRECISION, cg%mg%f_tgt%seg(g)%proc, tag1, comm, req(nr), ierr)

         enddo
         cgl => cgl%nxt
      enddo

      if (nr>0) call MPI_Waitall(nr, req(:nr), status(:,:nr), ierr)

      ! interpolate received coarse data to the right place
      cgl => fine%first
      do while (associated(cgl))
         cg => cgl%cg
         if (allocated(cg%mg%c_tgt%seg)) then
            do g = 1, ubound(cg%mg%c_tgt%seg(:), dim=1)
               fse => cg%mg%c_tgt%seg(g)%se
               if (associated(this%first)) then
                  off1(:) = mod(cg%off(:)+fse(:, LO) - this%first%cg%ijkse(:, LO), 2_LONG) !!!
               else
                  off1(:) = mod(cg%off(:)+fse(:, LO) - 1, 2_LONG) !!!
               endif
               do k = fse(zdim, LO), fse(zdim, HI)
                  kc = (k-fse(zdim, LO)+off1(zdim))/2 + 1
                  do j = fse(ydim, LO), fse(ydim, HI)
                     jc = (j-fse(ydim, LO)+off1(ydim))/2 + 1
                     do i = fse(xdim, LO), fse(xdim, HI)
                        ic = (i-fse(xdim, LO)+off1(xdim))/2 + 1
                        cg%q(iv)%arr(i, j, k) = cg%mg%c_tgt%seg(g)%buf(ic, jc, kc)
                     enddo
                  enddo
               enddo
            enddo
         endif
         cgl => cgl%nxt
      enddo

   end subroutine prolong_level0

!> \brief Print detailed information about current level decomposition

   subroutine print_segments(this)

      use constants,  only: LO, HI
      use dataio_pub, only: printinfo, msg
      use mpisetup,   only: proc, FIRST, LAST, master, nproc

      implicit none

      class(cg_list_level), intent(in) :: this

      integer :: p, i, hl
      integer(kind=8) :: ccnt
      real, allocatable, dimension(:) :: maxcnt
      character(len=len(msg)) :: header

      if (.not. master) return

      !call dom%print_me
      allocate(maxcnt(FIRST:LAST))
      maxcnt(:) = 0
      do p = FIRST, LAST
         hl = 0
         do i= 1, ubound(this%pse(p)%sel(:, :, :), dim=1)
            ccnt = product(this%pse(p)%sel(i, :, HI) - this%pse(p)%sel(i, :, LO) + 1)
            maxcnt(p) = maxcnt(p) + ccnt
            if (i == 1) then
               write(header, '(a,i4)')"[gc_list:print_segments] segment @", p
               hl = len_trim(header)
            else
               header = repeat(" ", hl)
            endif
            write(msg,'(2a,2(3i6,a),i8,a)') header(:hl), " : [", this%pse(p)%sel(i, :, LO), "] : [", this%pse(p)%sel(i, :, HI), "] #", ccnt, " cells"
            call printinfo(msg)
         enddo
      enddo
      write(msg,'(a,i3,a,f8.5)')"[gc_list:print_segments] Load balance at level ", this%lev," : ",product(real(this%n_d(:)))/(nproc*maxval(maxcnt(:)))
      !> \todo add calculation of total internal boundary surface in cells
      call printinfo(msg)
      deallocate(maxcnt)

   end subroutine print_segments

end module cg_list_lev
