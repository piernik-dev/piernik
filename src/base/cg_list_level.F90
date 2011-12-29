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

!!$ ============================================================================
!>
!! \brief Simplest restriction (averaging).
!! \todo implement high order restriction and test its influence on V-cycle convergence rate
!<

   subroutine restrict_level(this, iv)

      use constants,  only: xdim, ydim, zdim, LO, HI, LONG, I_ONE
      use dataio_pub, only: msg, warn, die
      use domain,     only: dom
      use mpisetup,   only: proc, comm, ierr, req, status
      use mpi,        only: MPI_DOUBLE_PRECISION

      implicit none

      class(cg_list_level), intent(inout), target  :: this !< object invoking type-bound procedure
      integer(kind=4), intent(in)      :: iv               !< variable to be restricted

      integer(kind=4), parameter :: tag1 = 1
      type(cg_list_level), pointer :: coarse
      integer :: g, g1, d
      integer(kind=8), dimension(:,:), pointer :: fse, cse ! shortcuts for fine segment and coarse segment
      integer(kind=8) :: i, j, k, ic, jc, kc
      integer(kind=8), dimension(xdim:zdim) :: off1
      real :: norm
      integer(kind=4) :: nr

      if (iv < lbound(this%first%cg%q, dim=1) .or. iv > ubound(this%first%cg%q, dim=1)) call die("[gc_list:restrict_level] Invalid variable index.")

      coarse => this%coarser
      if (.not. associated(coarse)) then
         write(msg,'(a,i3)')"[gc_list:restrict_level] no coarse level here: ", this%lev
         call warn(msg) ! can't restrict base level
      else
         !OPT find a way to reduce this to areas with nonlocal incoming restriction
         coarse%first%cg%q(iv)%arr(:,:,:) = 0. ! disables check_dirty
      endif

      nr = 0
      if (allocated(coarse%first%cg%mg%f_tgt%seg)) then
         do g = 1, ubound(coarse%first%cg%mg%f_tgt%seg(:), dim=1)
            if (coarse%first%cg%mg%f_tgt%seg(g)%proc /= proc) then
               nr = nr + I_ONE
               call MPI_Irecv(coarse%first%cg%mg%f_tgt%seg(g)%buf(1, 1, 1), size(coarse%first%cg%mg%f_tgt%seg(g)%buf(:, :, :)), MPI_DOUBLE_PRECISION, coarse%first%cg%mg%f_tgt%seg(g)%proc, tag1, comm, req(nr), ierr)
            endif
         enddo
      endif

      do g = 1, ubound(this%first%cg%mg%c_tgt%seg(:), dim=1)

         fse => this%first%cg%mg%c_tgt%seg(g)%se

         off1(:) = mod(this%first%cg%off(:), 2_LONG)
         norm = 1./2**dom%eff_dim
         if (this%first%cg%mg%c_tgt%seg(g)%proc == proc) then

            nullify(cse)
            do g1 = 1, ubound(coarse%first%cg%mg%f_tgt%seg(:), dim=1) !> \todo should be set up in init_multigrid
               if (coarse%first%cg%mg%f_tgt%seg(g1)%proc == proc) then
                  if (.not. associated(cse)) then
                     cse => coarse%first%cg%mg%f_tgt%seg(g1)%se
                  else
                     call die("[gc_list:restrict_level] multiple local coarse grid targets")
                  endif
               endif
            enddo
            if (.not. associated(cse)) call die("[gc_list:restrict_level] cannot find local coarse grid")

            do d = xdim, zdim ! debug
               if (dom%has_dir(d)) then
                  if ((fse(d, LO)-off1(d)-dom%nb-1)/2+cse(d, LO) < cse(d, LO)) call die("mv:rl <cse")
                  if ((fse(d, HI)-off1(d)-dom%nb-1)/2+cse(d, LO) > cse(d, HI)) call die("mv:rl >cse")
               endif
            enddo
            ! OPT: completely unoptimized,
            ! note that e.g. 10 cells on fine grid may contribute to 5 or 6 cells on coarse grid, depending on offset and both cases require a bit different code
            ! \todo use array sections to restrict the fully covered interior of coarse segment and finish the boundaries, where necessary
            ! OPT: old code operations: Ir:Dr:Dw:Dr_m:Dw_m = 13:8:1:1:0.1, this code 130:20:8:1.4:0.2
            ! OPT: computation of ic consumes ~30% Ir
            ! OPT: convert the loop over i (at least) to array section operation
            ! previous implementation of restriction (sum of array sections without loops) was probably faster and did produce slightly different results
            ! (with typical relative difference at level 1e-15 due to different numerical roundoffs)
            ! It is possible that future optimisations will affect the results and bit-to-bit comparisions (with h5diff) will fail
            do k = fse(zdim, LO), fse(zdim, HI)
               kc = (k-off1(zdim)-this%first%cg%ks)/2+cse(zdim, LO)
               do j = fse(ydim, LO), fse(ydim, HI)
                  jc = (j-off1(ydim)-this%first%cg%js)/2+cse(ydim, LO)
                  do i = fse(xdim, LO), fse(xdim, HI)
                     ic = (i-off1(xdim)-this%first%cg%is)/2+cse(xdim, LO)
                     coarse%first%cg%q(iv)%arr(ic, jc, kc) = coarse%first%cg%q(iv)%arr(ic, jc, kc) + this%first%cg%q(iv)%arr(i, j, k) * norm
                  enddo
               enddo
            enddo
            !\todo add geometrical terms to improve convergence on cylindrical grids
         else
            ! OPT: see the notes above
            this%first%cg%mg%c_tgt%seg(g)%buf(:, :, :) = 0.
            off1(:) = mod(this%first%cg%off(:)+fse(:, LO) - this%first%cg%ijkse(:, LO), 2_LONG)
            do k = fse(zdim, LO), fse(zdim, HI)
               kc = (k-fse(zdim, LO)+off1(zdim))/2 + 1
               do j = fse(ydim, LO), fse(ydim, HI)
                  jc = (j-fse(ydim, LO)+off1(ydim))/2 + 1
                  do i = fse(xdim, LO), fse(xdim, HI)
                     ic = (i-fse(xdim, LO)+off1(xdim))/2 + 1
                     this%first%cg%mg%c_tgt%seg(g)%buf(ic, jc, kc) = this%first%cg%mg%c_tgt%seg(g)%buf(ic, jc, kc) + this%first%cg%q(iv)%arr(i, j, k) * norm
                  enddo
               enddo
            enddo
            nr = nr + I_ONE
            call MPI_Isend(this%first%cg%mg%c_tgt%seg(g)%buf(1, 1, 1), size(this%first%cg%mg%c_tgt%seg(g)%buf(:, :, :)), MPI_DOUBLE_PRECISION, this%first%cg%mg%c_tgt%seg(g)%proc, tag1, comm, req(nr), ierr)
         endif
      enddo

      if (nr>0) call MPI_Waitall(nr, req(:nr), status(:,:nr), ierr)

      if (allocated(coarse%first%cg%mg%f_tgt%seg)) then
         do g = 1, ubound(coarse%first%cg%mg%f_tgt%seg(:), dim=1)
            if (coarse%first%cg%mg%f_tgt%seg(g)%proc /= proc) then
               cse => coarse%first%cg%mg%f_tgt%seg(g)%se
               coarse%first%cg%q(iv)%arr     (cse(xdim, LO):cse(xdim, HI), cse(ydim, LO):cse(ydim, HI), cse(zdim, LO):cse(zdim, HI)) = &
                    coarse%first%cg%q(iv)%arr(cse(xdim, LO):cse(xdim, HI), cse(ydim, LO):cse(ydim, HI), cse(zdim, LO):cse(zdim, HI)) + coarse%first%cg%mg%f_tgt%seg(g)%buf(:, :, :)
            endif
         enddo
      endif

   end subroutine restrict_level

!!$ ============================================================================
!>
!! \brief 0th order prolongation : injection
!<

   subroutine prolong_level0(this, iv)

      use constants,  only: xdim, ydim, zdim, LO, HI, LONG, I_ONE
      use dataio_pub, only: msg, warn, die
      use domain,     only: dom
      use mpisetup,   only: proc, comm, ierr, req, status
      use mpi,        only: MPI_DOUBLE_PRECISION

      implicit none

      class(cg_list_level), intent(inout), target  :: this !< object invoking type-bound procedure
      integer(kind=4), intent(in)         :: iv    !< variable to be prolonged

      integer(kind=4), parameter :: tag1 = 1
      type(cg_list_level), pointer :: fine
      integer :: g, g1, d
      integer(kind=8), dimension(:,:), pointer :: fse, cse ! shortcuts for fine segment and coarse segment
      integer(kind=8) :: i, j, k, ic, jc, kc
      integer(kind=8), dimension(xdim:zdim) :: off1
      integer(kind=4) :: nr

      if (iv < lbound(this%first%cg%q, dim=1) .or. iv > ubound(this%first%cg%q, dim=1)) call die("[gc_list:prolong_level0] Invalid variable index.")

      fine => this%finer
      if (.not. associated(fine)) then
         write(msg,'(a,i3)')"[gc_list:restrict_level] no fine level here: ", this%lev
         call warn(msg) ! can't prolong finest level
      else
         ! OPT: try to remove or limit this  ~20% Ir, ~50% Dw_m
         fine%first%cg%q(iv)%arr(:,:,:) = 0. ! disables check_dirty
      endif

      nr = 0
      if (allocated(fine%first%cg%mg%c_tgt%seg)) then
         do g = 1, ubound(fine%first%cg%mg%c_tgt%seg(:), dim=1)
            if (fine%first%cg%mg%c_tgt%seg(g)%proc /= proc) then
               nr = nr + I_ONE
               call MPI_Irecv(fine%first%cg%mg%c_tgt%seg(g)%buf(1, 1, 1), size(fine%first%cg%mg%c_tgt%seg(g)%buf(:, :, :)), MPI_DOUBLE_PRECISION, fine%first%cg%mg%c_tgt%seg(g)%proc, tag1, comm, req(nr), ierr)
            endif
         enddo
      endif

      off1(:) = mod(fine%first%cg%off(:), 2_LONG)
      do g = 1, ubound(this%first%cg%mg%f_tgt%seg(:), dim=1)

         cse => this%first%cg%mg%f_tgt%seg(g)%se

         if (this%first%cg%mg%f_tgt%seg(g)%proc == proc) then

            nullify(fse)
            do g1 = 1, ubound(fine%first%cg%mg%c_tgt%seg(:), dim=1) !> \todo should be set up in init_multigrid
               if (fine%first%cg%mg%c_tgt%seg(g1)%proc == proc) then
                  if (.not. associated(fse)) then
                     fse => fine%first%cg%mg%c_tgt%seg(g1)%se
                  else
                     call die("[gc_list:prolong_level0] multiple local fine grid targets")
                  endif
               endif
            enddo
            if (.not. associated(fse)) call die("[gc_list:prolong_level0] cannot find local fine grid")

            ! No guardcells required here

            ! Possible optimization candidate: reduce L1 and L2 cache misses on both read and write (RBGS only, secondary importance)
            do d = xdim, zdim ! debug
               if (dom%has_dir(d)) then
                  if ((fse(d, LO)-off1(d)-dom%nb-1)/2+cse(d, LO) < cse(d, LO)) call die("mv:rl <cse")
                  if ((fse(d, HI)-off1(d)-dom%nb-1)/2+cse(d, LO) > cse(d, HI)) call die("mv:rl >cse")
               endif
            enddo
            ! OPT: completely unoptimized
            ! OPT: old code operations: Ir:Dr:Dw:Dr_m:Dw_m = 37:9:11:1.5:1.5, new code = 72:10:9:0.1:1.7
            ! OPT: computation of ic consumes ~25% Ir
            do k = fse(zdim, LO), fse(zdim, HI)
               kc = (k-off1(zdim)-this%first%cg%ks)/2+cse(zdim, LO)
               do j = fse(ydim, LO), fse(ydim, HI)
                  jc = (j-off1(ydim)-this%first%cg%js)/2+cse(ydim, LO)
                  do i = fse(xdim, LO), fse(xdim, HI)
                     ic = (i-off1(xdim)-this%first%cg%is)/2+cse(xdim, LO)
                     fine%first%cg%q(iv)%arr(i, j, k) = this%first%cg%q(iv)%arr(ic, jc, kc)
                  enddo
               enddo
            enddo
         else
            nr = nr + I_ONE
            this%first%cg%mg%f_tgt%seg(g)%buf(:, :, :) = this%first%cg%q(iv)%arr(cse(xdim, LO):cse(xdim, HI), cse(ydim, LO):cse(ydim, HI), cse(zdim, LO):cse(zdim, HI))
            call MPI_Isend(this%first%cg%mg%f_tgt%seg(g)%buf(1, 1, 1), size(this%first%cg%mg%f_tgt%seg(g)%buf(:, :, :)), MPI_DOUBLE_PRECISION, this%first%cg%mg%f_tgt%seg(g)%proc, tag1, comm, req(nr), ierr)
         endif

      enddo

      if (nr>0) call MPI_Waitall(nr, req(:nr), status(:,:nr), ierr)

      if (allocated(fine%first%cg%mg%c_tgt%seg)) then
         do g = 1, ubound(fine%first%cg%mg%c_tgt%seg(:), dim=1)
            if (fine%first%cg%mg%c_tgt%seg(g)%proc /= proc) then
               fse => fine%first%cg%mg%c_tgt%seg(g)%se
               off1(:) = mod(fine%first%cg%off(:)+fse(:, LO) - this%first%cg%ijkse(:, LO), 2_LONG)
               do k = fse(zdim, LO), fse(zdim, HI)
                  kc = (k-fse(zdim, LO)+off1(zdim))/2 + 1
                  do j = fse(ydim, LO), fse(ydim, HI)
                     jc = (j-fse(ydim, LO)+off1(ydim))/2 + 1
                     do i = fse(xdim, LO), fse(xdim, HI)
                        ic = (i-fse(xdim, LO)+off1(xdim))/2 + 1
                        fine%first%cg%q(iv)%arr(i, j, k) = fine%first%cg%mg%c_tgt%seg(g)%buf(ic, jc, kc)
                     enddo
                  enddo
               enddo
            endif
         enddo
      endif

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
