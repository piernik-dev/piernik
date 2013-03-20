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
!! \brief Update internal boundaries
!!
!! \details This module contains subroutines that are responsible for preparing internal and external boundary cells for cg%q(:) and cg%w(:) arrays.
!! For simplicity, let's assume that all boundaries that rely on MPI communication are internal.
!! This implies that periodic, corner shear and fine-coarse boundaries are also "internal"
!!
!! Note that this routine may not properly update some layers of guardcells when number of guardcell layers exceeds number of active cells.
!! Appropriate checks should be made in decompose_patch routine.
!!
!! \todo integrate here as much stuff from fluidboundaries, magboundaries, etc. as possible.
!!
!! \todo implement a way to ensure proper block pairing for leaves list.
!<

module cg_list_bnd

! pulled by ANY

   use cg_list_dataop, only: cg_list_dataop_T

   implicit none

   private
   public :: cg_list_bnd_T

   !>
   !! \brief Lists of grid containers with boundary update
   !!
   !! \details Procedure internal_boundaries requires existence of matching mpi data structures, which is done in cg_level type.
   !! Thus this type is usably only when the lost consist of one or moge full cg levels.
   !<

   type, extends(cg_list_dataop_T) :: cg_list_bnd_T
    contains
      procedure          :: internal_boundaries_3d     !< A wrapper that calls internal_boundaries for 3D arrays stored in cg%q(:)
      procedure          :: internal_boundaries_4d     !< A wrapper that calls internal_boundaries for 4D arrays stored in cg%w(:)
      procedure, private :: internal_boundaries        !< Exchanges guardcells for BND_MPI and BND_PER boundaries (internal and periodic external boundaries)
      procedure          :: external_boundaries        !< Set up external boundary values
      procedure          :: clear_boundaries           !< Clear (set to 0) all boundaries
      procedure          :: dirty_boundaries           !< Put dirty values to all boundaries
      procedure          :: level_3d_boundaries        !< Perform internal boundary exchanges and external boundary extrapolations on 3D named arrays
      procedure          :: level_4d_boundaries        !< Perform internal boundary exchanges and external boundary extrapolations on 4D named arrays
      !> \todo move routines for external guardcells for rank-4 arrays here as well (fluidboundaries and magboundaries)
   end type cg_list_bnd_T

contains

!> \brief A wrapper that calls internal_boundaries for 3D arrays stored in cg%q(:)

   subroutine internal_boundaries_3d(this, ind, dim)

      implicit none

      class(cg_list_bnd_T),      intent(in) :: this   !< the list on which to perform the boundary exchange
      integer(kind=4),           intent(in) :: ind    !< index of cg%q(:) 3d array
      integer(kind=4), optional, intent(in) :: dim    !< do the internal boundaries only in the specified dimension

      call internal_boundaries(this, ind, .true., dim)

   end subroutine internal_boundaries_3d

!> \brief A wrapper that calls internal_boundaries for 4D arrays stored in cg%w(:)

   subroutine internal_boundaries_4d(this, ind, dim)

      implicit none

      class(cg_list_bnd_T),      intent(in) :: this !< the list on which to perform the boundary exchange
      integer(kind=4),           intent(in) :: ind  !< index of cg%w(:) 4d array
      integer(kind=4), optional, intent(in) :: dim  !< do the internal boundaries only in the specified dimension

      call internal_boundaries(this, ind, .false., dim)

   end subroutine internal_boundaries_4d

!>
!! \brief This routine exchanges guardcells for BND_MPI and BND_PER boundaries on rank-3 and rank-4 arrays
!!
!! \details This routine should not be called directly. Appropriate wrappers for rank-3 and rank-4 arrays are provided above.
!! The corners should be properly updated if this%[io]_bnd(:, ind) was set up appropriately and this routine is called separately for each dimension.
!!
!! OPT: Since r6414 (see also r6406 .. 6409) we don't use MPI types for boundary exchanges.
!! The implementation with MPI_Types was a bit faster, but it had a severe limitation of total number of MPI types declared (in my case it was 261888)
!! That amount was exceeded much faster than expected because for each grid container we declared number_of_neighbours*dom%nb*(1+size(wna%lst)) MPI types.
!! Note that some of them were never used.
!! \todo Try to define MPI_types for communication right before MPI_Isend/MPI_Irecv calls and release just after use. Then compare performance.
!!
!! \todo Check how much performance is lost due to using MPI calls even for local copies. Decide whether it is worth to convert local MPI calls to direct memory copies.
!! For othes suggestions on performance optimisation see decription of cg_level::mpi_bnd_types.
!!
!! \warning this == leaves could be unsafe: need to figure out how to handle unneeded edges; this == all_cg or base%level or other concatenation of whole levels should work well
!<

   subroutine internal_boundaries(this, ind, tgt3d, dim)

      use cg_list,          only: cg_list_element
      use constants,        only: xdim, ydim, zdim, LO, HI, I_ONE, I_TWO
      use dataio_pub,       only: die, warn
      use domain,           only: dom
      use grid_cont,        only: grid_container, segment
      use mpi,              only: MPI_DOUBLE_PRECISION, MPI_STATUS_SIZE
      use mpisetup,         only: comm, mpi_err, req, inflate_req
      use named_array_list, only: wna

      implicit none

      class(cg_list_bnd_T),      intent(in)        :: this   !< the list on which to perform the boundary exchange
      integer(kind=4),           intent(in)        :: ind    !< index of cg%q(:) 3d array or cg%w(:) 4d array
      logical,                   intent(in)        :: tgt3d  !< .true. for cg%q, .false. for cg%w
      integer(kind=4), optional, intent(in)        :: dim    !< do the internal boundaries only in the specified dimension

      integer                                      :: g, d
      integer(kind=4)                              :: nr     !< index of first free slot in req and status arrays
      logical, dimension(xdim:zdim)                :: dmask
      type(grid_container),     pointer            :: cg
      type(cg_list_element),    pointer            :: cgl
      real, dimension(:,:,:),   pointer            :: pa3d
      real, dimension(:,:,:,:), pointer            :: pa4d
      logical                                      :: active
      type(segment), pointer                       :: i_seg, o_seg !< shortcuts
      integer(kind=4), allocatable, dimension(:,:) :: mpistatus !< status array for MPI_Waitall

      dmask(:) = dom%has_dir(:)
      if (present(dim)) then
         dmask(:) = .false.
         dmask(dim) = dom%has_dir(dim)
      endif

      nr = 0
      cgl => this%first
      do while (associated(cgl))
         cg => cgl%cg

         ! exclude non-multigrid variables below base level
         if (tgt3d) then
            active = associated(cg%q(ind)%arr)
         else
            active = associated(cg%w(ind)%arr)
         endif

         do d = xdim, zdim
            if (dmask(d) .and. active) then
               if (allocated(cg%i_bnd(d)%seg)) then
                  if (.not. allocated(cg%o_bnd(d)%seg)) call die("[cg_list_bnd:internal_boundaries] cg%i_bnd without cg%o_bnd")
                  if (ubound(cg%i_bnd(d)%seg(:), dim=1) /= ubound(cg%o_bnd(d)%seg(:), dim=1)) &
                       call die("[cg_list_bnd:internal_boundaries] cg%i_bnd differs in number of entries from cg%o_bnd")
                  do g = lbound(cg%i_bnd(d)%seg(:), dim=1), ubound(cg%i_bnd(d)%seg(:), dim=1)

                     if (nr+I_TWO >  ubound(req(:), dim=1)) call inflate_req
                     i_seg => cg%i_bnd(d)%seg(g)
                     o_seg => cg%o_bnd(d)%seg(g)

                     !> \deprecated: A lot of semi-duplicated code below
                     if (tgt3d) then
                        if (ind > ubound(cg%q(:), dim=1) .or. ind < lbound(cg%q(:), dim=1)) call die("[cg_list_bnd:internal_boundaries] wrong 3d index")

                        if (allocated(i_seg%buf)) then
                           call warn("clb:ib allocated i-buf")
                           deallocate(i_seg%buf) !> \todo check shape and recycle if possible
                        endif
                        allocate(i_seg%buf(i_seg%se(xdim, HI) - i_seg%se(xdim, LO) + 1, &
                             &             i_seg%se(ydim, HI) - i_seg%se(ydim, LO) + 1, &
                             &             i_seg%se(zdim, HI) - i_seg%se(zdim, LO) + 1))
                        pa3d => cg%q(ind)%span(i_seg%se(:,:))
                        i_seg%buf(:,:,:) = pa3d(:,:,:)
                        call MPI_Irecv(i_seg%buf, size(i_seg%buf), MPI_DOUBLE_PRECISION, i_seg%proc, i_seg%tag, comm, req(nr+I_ONE), mpi_err)

                        if (allocated(o_seg%buf)) then
                           call warn("clb:ib allocated o-buf")
                           deallocate(o_seg%buf) !> \todo check shape and recycle if possible
                        endif
                        allocate(o_seg%buf(o_seg%se(xdim, HI) - o_seg%se(xdim, LO) + 1, &
                             &             o_seg%se(ydim, HI) - o_seg%se(ydim, LO) + 1, &
                             &             o_seg%se(zdim, HI) - o_seg%se(zdim, LO) + 1))
                        pa3d => cg%q(ind)%span(o_seg%se(:,:))
                        o_seg%buf(:,:,:) = pa3d(:,:,:)
                        call MPI_Isend(o_seg%buf, size(o_seg%buf), MPI_DOUBLE_PRECISION, o_seg%proc, o_seg%tag, comm, req(nr+I_TWO), mpi_err)

                     else
                        if (ind > ubound(cg%w(:), dim=1) .or. ind < lbound(cg%w(:), dim=1)) call die("[cg_list_bnd:internal_boundaries] wrong 4d index")

                        if (allocated(i_seg%buf4)) then
                           call warn("clb:ib allocated i-buf")
                           deallocate(i_seg%buf4) !> \todo check shape and recycle if possible
                        endif
                        allocate(i_seg%buf4(wna%lst(ind)%dim4, &
                             &              i_seg%se(xdim, HI) - i_seg%se(xdim, LO) + 1, &
                             &              i_seg%se(ydim, HI) - i_seg%se(ydim, LO) + 1, &
                             &              i_seg%se(zdim, HI) - i_seg%se(zdim, LO) + 1))
                        !>
                        !! \todo optimize me
                        !! Following 2 lines (along with other occurences of
                        !! similar constructs) cause ~10% perfomance drop wrt
                        !! in-place communication using MPI_Types.
                        !! They can be optimized by using explicit loop over
                        !! last index:
                        !!    do ni = lbound(i_seg%buf4, 4), ubound(i_seg%buf4, 4)
                        !!       hhi = i_seg%se(zdim,LO) - 1 + ni
                        !!       i_seg%buf4(:,:,:,ni) = &
                        !!          cg%w(ind)%arr(:,i_seg%se(xdim,LO):i_seg%se(xdim,HI),i_seg%se(ydim,LO):i_seg%se(ydim,HI),hhi)
                        !!    enddo
                        !<
                        pa4d => cg%w(ind)%span(i_seg%se(:,:))
                        i_seg%buf4(:,:,:,:) = pa4d(:,:,:,:)
                        call MPI_Irecv(i_seg%buf4, size(i_seg%buf4), MPI_DOUBLE_PRECISION, i_seg%proc, i_seg%tag, comm, req(nr+I_ONE), mpi_err)

                        if (allocated(o_seg%buf4)) then
                           call warn("clb:ib allocated o-buf")
                           deallocate(o_seg%buf4) !> \todo check shape and recycle if possible
                        endif
                        allocate(o_seg%buf4(wna%lst(ind)%dim4, &
                             &              o_seg%se(xdim, HI) - o_seg%se(xdim, LO) + 1, &
                             &              o_seg%se(ydim, HI) - o_seg%se(ydim, LO) + 1, &
                             &              o_seg%se(zdim, HI) - o_seg%se(zdim, LO) + 1))
                        !>
                        !! \todo optimize me
                        !! do ni = lbound(o_seg%buf4, 4), ubound(o_seg%buf4, 4)
                        !!    hhi = o_seg%se(zdim,LO) - 1 + ni
                        !!    o_seg%buf4(:,:,:,ni) = &
                        !!       cg%w(ind)%arr(:,o_seg%se(xdim,LO):o_seg%se(xdim,HI),o_seg%se(ydim,LO):o_seg%se(ydim,HI),hhi)
                        !! enddo
                        !<
                        pa4d => cg%w(ind)%span(o_seg%se(:,:))
                        o_seg%buf4(:,:,:,:) = pa4d(:,:,:,:)
                        call MPI_Isend(o_seg%buf4, size(o_seg%buf4), MPI_DOUBLE_PRECISION, o_seg%proc, o_seg%tag, comm, req(nr+I_TWO), mpi_err)

                     endif
                     nr = nr + I_TWO
                  enddo
               else
                  if (allocated(cg%o_bnd(d)%seg)) call die("[cg_list_bnd:internal_boundaries] cg%o_bnd without cg%i_bnd")
               endif
            endif
         enddo

         cgl => cgl%nxt
      enddo

      allocate(mpistatus(MPI_STATUS_SIZE, nr))
      call MPI_Waitall(nr, req(:nr), mpistatus, mpi_err)
      deallocate(mpistatus)

      ! Move the received data from buffers to the right place. Deallocate buffers
      cgl => this%first
      do while (associated(cgl))
         cg => cgl%cg
         ! exclude non-multigrid variables below base level
         if (tgt3d) then
            active = associated(cg%q(ind)%arr)
         else
            active = associated(cg%w(ind)%arr)
         endif

         do d = xdim, zdim
            if (dmask(d) .and. active) then
               if (allocated(cg%i_bnd(d)%seg)) then
                  ! sanity checks are already done
                  do g = lbound(cg%i_bnd(d)%seg(:), dim=1), ubound(cg%i_bnd(d)%seg(:), dim=1)

                     i_seg => cg%i_bnd(d)%seg(g)
                     o_seg => cg%o_bnd(d)%seg(g)

                     if (tgt3d) then
                        pa3d => cg%q(ind)%span(i_seg%se(:,:))
                        pa3d(:,:,:) = i_seg%buf(:,:,:)
                        deallocate(i_seg%buf)
                        deallocate(o_seg%buf)
                     else
                        !>
                        !! \todo optimize me
                        !! do ni = lbound(i_seg%buf4, 4), ubound(i_seg%buf4, 4)
                        !!    hhi = i_seg%se(zdim,LO) - 1 + ni
                        !!    cg%w(ind)%arr(:,i_seg%se(xdim,LO):i_seg%se(xdim,HI),i_seg%se(ydim,LO):i_seg%se(ydim,HI),hhi) = &
                        !!       i_seg%buf4(:,:,:,ni)
                        !! enddo
                        !<
                        pa4d => cg%w(ind)%span(i_seg%se(:,:))
                        pa4d(:,:,:,:) = i_seg%buf4(:,:,:,:)
                        deallocate(i_seg%buf4)
                        deallocate(o_seg%buf4)
                     endif
                  enddo
               endif
            endif
         enddo

         cgl => cgl%nxt
      enddo

   end subroutine internal_boundaries

!>
!! \brief Set up external boundary values
!!
!! \details This is purely local operation. We don't bother here to limit the number of layers to be filled, because it is cheap.
!! Can be moved to cg_list or cg_list_dataop.
!<

   subroutine external_boundaries(this, ind, area_type, bnd_type)

      use cg_list,    only: cg_list_element
      use constants,  only: ndims, xdim, ydim, zdim, LO, HI, AT_NO_B, I_ONE, I_TWO, I_THREE, pLOR, &
           &                BND_PER, BND_MPI, BND_FC, BND_MPI_FC, BND_SHE, BND_COR, BND_REF, BND_NEGREF, BND_ZERO, BND_XTRAP, BND_NONE
      use dataio_pub, only: die, msg
      use domain,     only: dom
      use grid_cont,  only: grid_container
      use mpisetup,   only: piernik_MPI_Allreduce

      implicit none

      class(cg_list_bnd_T),      intent(in)   :: this       !< the list on which to perform the boundary exchange
      integer(kind=4),           intent(in)   :: ind        !< Negative value: index of cg%q(:) 3d array
      integer(kind=4), optional, intent(in)   :: area_type  !< defines how do we treat boundaries
      integer(kind=4), optional, intent(in)   :: bnd_type   !< Override default boundary type on external boundaries (useful in multigrid solver).
                                                            !< Note that BND_PER, BND_MPI, BND_SHE and BND_COR aren't external and cannot be overridden

      integer(kind=4)                         :: lh, clh, d, b_type, i
      integer(kind=4), dimension(ndims,LO:HI) :: l, r, rh
      logical                                 :: dodie
      type(cg_list_element),  pointer         :: cgl
      type(grid_container),   pointer         :: cg
      real, dimension(:,:,:), pointer         :: pa3d

      dodie = .false.

      cgl => this%first
      do while (associated(cgl))
         cg => cgl%cg

         if (ind > ubound(cg%q(:), dim=1) .or. ind < lbound(cg%q(:), dim=1)) call die("[cg_list_bnd:external_boundaries] wrong 3d index")
         pa3d => cg%q(ind)%arr

         do d = xdim, zdim
            if (dom%has_dir(d)) then
               l = reshape([lbound(pa3d, kind=4),ubound(pa3d, kind=4)],shape=[ndims,HI]) ; r = l
               do lh = LO, HI

                  select case (cg%bnd(d, lh))
                     case (BND_PER, BND_MPI, BND_FC, BND_MPI_FC) ! Already done in internal_bnd or arr3d_boundaries
                     case (BND_SHE) !> \todo move appropriate code from poissonsolver::poisson_solve or do nothing. or die until someone really needs SHEAR.
                        write(msg,*) "[cg_list_bnd:external_boundaries] 'she' not implemented"
                        dodie = .true.
                     case (BND_COR)
                        if (present(area_type)) then
                           if (area_type /= AT_NO_B) cycle
                        endif
                        write(msg,*) "[cg_list_bnd:external_boundaries] 'cor' not implemented"
                        dodie = .true.
                     case default ! Set gradient == 0 on the external boundaries
                        if (present(area_type)) then
                           if (area_type /= AT_NO_B) cycle
                        endif
                        b_type = cg%bnd(d, lh)
                        if (present(bnd_type)) b_type = bnd_type
                        select case (b_type)
                           case (BND_REF)  ! reflecting BC (e.g. homogeneous Neumamnn)
                              ! there will be special rules for vector fields (velocity, magnetic) perpendicular to the given boundary (like BND_NEGREF)
                              do i = 1, dom%nb
                                 l(d,:) = cg%ijkse(d,lh)   -i     *(I_THREE-I_TWO*lh)
                                 r(d,:) = cg%ijkse(d,lh)+(i-I_ONE)*(I_THREE-I_TWO*lh)
                                 pa3d(l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) = pa3d(r(xdim,LO):r(xdim,HI),r(ydim,LO):r(ydim,HI),r(zdim,LO):r(zdim,HI))
                              enddo
                           case (BND_NEGREF)  ! reflecting BC (e.g. homogeneous Dirichlet BC with 0 at domain face)
                              do i = 1, dom%nb
                                 l(d,:) = cg%ijkse(d,lh)   -i     *(I_THREE-I_TWO*lh)
                                 r(d,:) = cg%ijkse(d,lh)+(i-I_ONE)*(I_THREE-I_TWO*lh)
                                 pa3d(l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) = - pa3d(r(xdim,LO):r(xdim,HI),r(ydim,LO):r(ydim,HI),r(zdim,LO):r(zdim,HI))
                              enddo
                           case (BND_ZERO)  ! zero BC (e.g. homogenous Dirichlet BC with 0 at first layer of cells)
                              clh = LO + HI - lh ; l(d,HI) = ubound(pa3d, dim=d, kind=4) ! restore after lh==LO case
                              l(d,clh) = cg%ijkse(d,lh)
                              pa3d(l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) = 0.
                           case (BND_NONE) ! remember to initialize everything first!
                           case (BND_XTRAP) !> \deprecated mixed-type BC: free flux; BEWARE: it is not protected from inflow
                              rh(:,:) = r(:,:) ; rh(d,:) = cg%ijkse(d,lh) ; r(d,:) = cg%ijkse(d,lh)+I_THREE-I_TWO*lh
                              do i = 1, dom%nb
                                 l(d,:) = cg%ijkse(d,lh)-i*(I_THREE-I_TWO*lh)
                                 pa3d(l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) = &
                                       (1+i) * pa3d(rh(xdim,LO):rh(xdim,HI),rh(ydim,LO):rh(ydim,HI),rh(zdim,LO):rh(zdim,HI)) &
                                         - i * pa3d( r(xdim,LO): r(xdim,HI), r(ydim,LO): r(ydim,HI), r(zdim,LO): r(zdim,HI))
                              enddo
                           case default ! BND_OUT, BND_OUTD, BND_OUTH, BND_OUTHD
                              r(d,:) = cg%ijkse(d,lh)
                              do i = 1, dom%nb
                                 l(d,:) = cg%ijkse(d,lh)   -i     *(I_THREE-I_TWO*lh)
                                 pa3d(l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) = pa3d(r(xdim,LO):r(xdim,HI),r(ydim,LO):r(ydim,HI),r(zdim,LO):r(zdim,HI))
                              enddo
                        end select
                  end select

               enddo
            endif
         enddo

         cgl => cgl%nxt
      enddo

      call piernik_MPI_Allreduce(dodie, pLOR)
      if (dodie) call die(msg)

   end subroutine external_boundaries

!> \brief Set zero to all boundaries (will defeat any attemts of use of dirty checks on boundaries)

   subroutine clear_boundaries(this, ind, value)

      use cg_list,   only: cg_list_element
      use constants, only: ndims, xdim, zdim, LO, HI
      use domain,    only: dom
      use grid_cont, only: grid_container

      implicit none

      class(cg_list_bnd_T), intent(in)        :: this  !< the list on which to perform the action
      integer(kind=4),      intent(in)        :: ind   !< Negative value: index of cg%q(:) 3d array
      real, optional,       intent(in)        :: value !< Value to be put in the boundaries (could be dirty)

      integer(kind=4)                         :: lh, clh, d
      integer(kind=4), dimension(ndims,LO:HI) :: l
      type(cg_list_element),  pointer         :: cgl
      type(grid_container),   pointer         :: cg
      real, dimension(:,:,:), pointer         :: pa3d
      real :: v

      !> \todo fill corners with big_float ?

      v = 0.
      if (present(value)) v=value
      cgl => this%first
      do while (associated(cgl))
         cg => cgl%cg
         do d = xdim, zdim
            if (dom%has_dir(d)) then
               do lh = LO, HI
                  l = reshape([lbound(cg%q(ind)%arr, kind=4),ubound(cg%q(ind)%arr, kind=4)],shape=[ndims,HI])
                  clh = LO + HI - lh
                  l(d,clh) = cg%ijkse(d,lh) + HI*lh-(LO+HI) ! -1 for LO, +1 for HI
                  pa3d => cg%q(ind)%span(l)
                  pa3d = v
               enddo
            endif
         enddo
         cgl => cgl%nxt
      enddo

   end subroutine clear_boundaries

!> \brief Put dirty values to all boundaries

   subroutine dirty_boundaries(this, ind)

      use constants, only: dirtyH

      implicit none

      class(cg_list_bnd_T), intent(in) :: this  !< the list on which to perform the action
      integer(kind=4),      intent(in) :: ind   !< Negative value: index of cg%q(:) 3d array

      call this%clear_boundaries(ind, value=dirtyH)

   end subroutine dirty_boundaries

!>
!! \brief This routine sets up all guardcells (internal and external) for given rank-3 arrays.
!!
!! \details No fine-coarse exchanges can be done here, see cg_level_connected::arr3d_boundaries for that feature
!<

   subroutine level_3d_boundaries(this, ind, area_type, bnd_type)

      use constants, only: AT_NO_B

      implicit none

      class(cg_list_bnd_T),      intent(in) :: this       !< the list on which to perform the boundary exchange
      integer(kind=4),           intent(in) :: ind        !< index of cg%q(:) 3d array
      integer(kind=4), optional, intent(in) :: area_type  !< defines how do we treat boundaries
      integer(kind=4), optional, intent(in) :: bnd_type   !< Override default boundary type on external boundaries (useful in multigrid solver).
                                                          !< Note that BND_PER, BND_MPI, BND_SHE and BND_COR aren't external and cannot be overridden

      logical                               :: do_permpi

      !> \todo fill corners with big_float ?

      do_permpi = .true.
      if (present(area_type)) then
         if (area_type /= AT_NO_B) do_permpi = .false.
      endif

      if (do_permpi) call this%internal_boundaries_3d(ind)

      call this%external_boundaries(ind, area_type, bnd_type)

   end subroutine level_3d_boundaries

!>
!! \brief This routine sets up all guardcells (internal and external) for given rank-4 arrays.
!!
!! \details No fine-coarse exchanges can be done here, see cg_level_connected::arr4d_boundaries for that feature
!<

   subroutine level_4d_boundaries(this, ind, area_type)

      use constants, only: AT_NO_B

      implicit none

      class(cg_list_bnd_T),      intent(in) :: this       !< the list on which to perform the boundary exchange
      integer(kind=4),           intent(in) :: ind        !< index of cg%w(:) 4d array
      integer(kind=4), optional, intent(in) :: area_type  !< defines how do we treat boundaries

      logical                               :: do_permpi

      !> \todo fill corners with big_float ?

      do_permpi = .true.
      if (present(area_type)) then
         if (area_type /= AT_NO_B) do_permpi = .false.
      endif

      if (do_permpi) call this%internal_boundaries_4d(ind)

!      call this%external_boundaries(ind, area_type, bnd_type) ! should call fluidboundaries:bnd_u, but that depends on hydrostatic too much

   end subroutine level_4d_boundaries

end module cg_list_bnd
