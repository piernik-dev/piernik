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
!! \brief Update internal and external boundaries
!!
!! \details This module contains subroutines that are responsible for preparing internal and external boundary cells for cg%q(:) and cg%w(:) arrays.
!! For simplicity, let's assume that all boundaries that rely on MPI communication are internal.
!! This implies that periodic, corner shear and fine-coarse boundaries are also "internal"
!!
!! Note that this routine may not properly update some layers of guardcells when number of guardcell layers exceeds number of active cells.
!! Appropriate checks should be made in decompose_patch routine.
!!
!! \todo integrate here as much stuff from magboundaries, etc. as possible.
!!
!! \todo implement a way to ensure proper block pairing for leaves list.
!<

module cg_list_bnd

! pulled by ANY

   use cg_list_dataop, only: cg_list_dataop_t
   use merge_segments, only: merge_segments_t

   implicit none

   private
   public :: cg_list_bnd_t

   !>
   !! \brief Lists of grid containers with boundary update
   !!
   !! \details Procedure internal_boundaries requires existence of matching mpi data structures, which is done in cg_level type.
   !! Thus this type is usable only when the list consist of one or more full cg levels.
   !<

   type, extends(cg_list_dataop_t), abstract :: cg_list_bnd_t
      type(merge_segments_t) :: ms                         !< merged segments
   contains
      procedure          :: level_3d_boundaries            !< Perform internal boundary exchanges and external boundary extrapolations on 3D named arrays
      procedure          :: level_4d_boundaries            !< Perform internal boundary exchanges and external boundary extrapolations on 4D named arrays
      procedure, private :: internal_boundaries            !< Exchanges guardcells for BND_MPI and BND_PER boundaries (internal and periodic external boundaries)
      procedure, private :: internal_boundaries_local      !< Exchanges guardcells between local grid containers
      procedure, private :: internal_boundaries_MPI_1by1   !< Exchanges guardcells with remote grid containers, one-by-one
      procedure, private :: internal_boundaries_MPI_merged !< Exchanges guardcells with remote grid containers with merged MPI messages
      procedure          :: clear_boundaries               !< Clear (set to 0) all boundaries
      procedure          :: dirty_boundaries               !< Put dirty values to all boundaries
      procedure          :: external_boundaries            !< Set up external boundary values
      procedure          :: bnd_u                          !< External (Non-MPI) boundary conditions for the fluid array: cg%u
      procedure          :: bnd_b                          !< External (Non-MPI) boundary conditions for the magnetic field array: cg%b
      !> \todo move routines for external guardcells for rank-4 arrays here as well (fluidboundaries and magboundaries)
   end type cg_list_bnd_t

contains

!>
!! \brief This routine sets up all guardcells (internal and external) for given rank-3 arrays.
!!
!! \details No fine-coarse exchanges can be done here, see cg_level_connected::arr3d_boundaries for that feature
!<

   subroutine level_3d_boundaries(this, ind, area_type, bnd_type, dir, nocorners)

      use constants, only: AT_NO_B

      implicit none

      class(cg_list_bnd_t),      intent(inout) :: this       !< the list on which to perform the boundary exchange
      integer(kind=4),           intent(in)    :: ind        !< index of cg%q(:) 3d array
      integer(kind=4), optional, intent(in)    :: area_type  !< defines how do we treat boundaries
      integer(kind=4), optional, intent(in)    :: bnd_type   !< Override default boundary type on external boundaries (useful in multigrid solver).
                                                          !< Note that BND_PER, BND_MPI, BND_SHE and BND_COR aren't external and cannot be overridden
      integer(kind=4), optional, intent(in)    :: dir        !< select only this direction
      logical,         optional, intent(in)    :: nocorners  !< .when .true. then don't care about proper edge and corner update

      logical :: do_permpi

      !> \todo fill corners with sth*dirtyH1 ?

      do_permpi = .true.
      if (present(area_type)) then
         if (area_type /= AT_NO_B) do_permpi = .false.
      endif

      if (do_permpi) call this%internal_boundaries(ind, .true., dir=dir, nocorners=nocorners)

      call this%external_boundaries(ind, area_type=area_type, bnd_type=bnd_type)

   end subroutine level_3d_boundaries

!>
!! \brief This routine sets up all guardcells (internal and external) for given rank-4 arrays.
!!
!! \details No fine-coarse exchanges can be done here, see cg_level_connected::arr4d_boundaries for that feature
!<

   subroutine level_4d_boundaries(this, ind, area_type, dir, nocorners)

      use constants, only: AT_NO_B

      implicit none

      class(cg_list_bnd_t),      intent(inout) :: this       !< the list on which to perform the boundary exchange
      integer(kind=4),           intent(in)    :: ind        !< index of cg%w(:) 4d array
      integer(kind=4), optional, intent(in)    :: area_type  !< defines how do we treat boundaries
      integer(kind=4), optional, intent(in)    :: dir        !< select only this direction
      logical,         optional, intent(in)    :: nocorners  !< .when .true. then don't care about proper edge and corner update

      logical :: do_permpi

      !> \todo fill corners with sth*dirtyH1 ?

      do_permpi = .true.
      if (present(area_type)) then
         if (area_type /= AT_NO_B) do_permpi = .false.
      endif

      if (do_permpi) call this%internal_boundaries(ind, .false., dir=dir, nocorners=nocorners)

!      call this%external_boundaries(ind, area_type=area_type, bnd_type=bnd_type) ! should call cg_list_bnd:bnd_u, but that depends on hydrostatic too much

   end subroutine level_4d_boundaries

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
!! For other suggestions on performance optimisation see description of cg_list_neighbors::find_neighbors_bruteforce
!!
!! \warning this == leaves could be unsafe: need to figure out how to handle unneeded edges; this == all_cg or base%level or other concatenation of whole levels should work well
!<

   subroutine internal_boundaries(this, ind, tgt3d, dir, nocorners)

      use constants, only: xdim, zdim, cor_dim, PPP_AMR
      use domain,    only: dom
      use global,    only: prefer_merged_MPI
      use ppp,       only: ppp_main

      implicit none

      class(cg_list_bnd_t),      intent(inout) :: this      !< the list on which to perform the boundary exchange
      integer(kind=4),           intent(in)    :: ind       !< index of cg%q(:) 3d array or cg%w(:) 4d array
      logical,                   intent(in)    :: tgt3d     !< .true. for cg%q, .false. for cg%w
      integer(kind=4), optional, intent(in)    :: dir       !< do the internal boundaries only in the specified dimension
      logical,         optional, intent(in)    :: nocorners !< .when .true. then don't care about proper edge and corner update

      logical, dimension(xdim:cor_dim) :: dmask
      character(len=*), parameter :: ib_label = "internal_boundaries", ibl_label = "internal_boundaries_local", ibm_label = "internal_boundaries_MPI_merged", ib1_label = "internal_boundaries_MPI_1by1"

      call ppp_main%start(ib_label)
      dmask(xdim:zdim) = dom%has_dir
      if (present(dir)) then
         dmask(xdim:zdim) = .false.
         dmask(dir) = dom%has_dir(dir)
      endif

      dmask(cor_dim) = .true.
      if (present(nocorners)) dmask(cor_dim) = .not. nocorners

      ! for CRESP (many components, big rank-4 arrays) it is definitely better to use internal_boundaries_MPI_1by1
      ! for multigrid (rank-3 arrays) it is better to use internal_boundaries_MPI_merged
      ! for non-CRESP (few components, rank-4 arrays, small bsize) it is a bit better to use internal_boundaries_MPI_merged
      !
      ! OPT: at what size of cg%w array (number of components, size of block) one approach wins over another?
      ! it seems that at [5, 16, 16, 16] internal_boundaries_MPI_merged and internal_boundaries_MPI_1by1 have similar performance

      if (this%ms%valid .and. (prefer_merged_MPI .or. tgt3d)) then
         call ppp_main%start(ibl_label)
         call internal_boundaries_local(this, ind, tgt3d, dmask)
         call ppp_main%stop(ibl_label)
         call ppp_main%start(ibm_label, PPP_AMR)
         call internal_boundaries_MPI_merged(this, ind, tgt3d, dmask)
         call ppp_main%stop(ibm_label, PPP_AMR)
      else
         call ppp_main%start(ib1_label, PPP_AMR)
         call internal_boundaries_MPI_1by1(this, ind, tgt3d, dmask)
         call ppp_main%stop(ib1_label, PPP_AMR)
      endif
      call ppp_main%stop(ib_label)

   end subroutine internal_boundaries

!> \brief This routine exchanges guardcells between local blocks for BND_MPI and BND_PER boundaries on rank-3 and rank-4 arrays.

   subroutine internal_boundaries_local(this, ind, tgt3d, dmask)

      use cg_cost_data,  only: I_OTHER
      use cg_list,       only: cg_list_element
      use constants,     only: xdim, ydim, zdim, LO, HI, cor_dim, INVALID
      use dataio_pub,    only: die
      use grid_cont,     only: grid_container
      use grid_cont_bnd, only: segment

      implicit none

      class(cg_list_bnd_t),             intent(in) :: this  !< the list on which to perform the boundary exchange
      integer(kind=4),                  intent(in) :: ind   !< index of cg%q(:) 3d array or cg%w(:) 4d array
      logical,                          intent(in) :: tgt3d !< .true. for cg%q, .false. for cg%w
      logical, dimension(xdim:cor_dim), intent(in) :: dmask !< .true. for the directions we want to exchange

      integer                           :: g, d, g_o, i
      integer(kind=8)                   :: j
      type(grid_container),     pointer :: cg
      type(cg_list_element),    pointer :: cgl
      real, dimension(:,:,:),   pointer :: pa3d, pa3d_o
      type(segment), pointer            :: i_seg, o_seg !< shortcuts

      cgl => this%first
      do while (associated(cgl))
         cg => cgl%cg
         call cg%costs%start

         do d = lbound(cg%i_bnd, dim=1), ubound(cg%i_bnd, dim=1)
            if (dmask(d) .and. is_active(cg, ind, tgt3d)) then
               if (allocated(cg%i_bnd(d)%seg)) then
                  if (.not. allocated(cg%o_bnd(d)%seg)) call die("[cg_list_bnd:internal_boundaries_local] cg%i_bnd without cg%o_bnd")
                  if (ubound(cg%i_bnd(d)%seg(:), dim=1) /= ubound(cg%o_bnd(d)%seg(:), dim=1)) &
                       call die("[cg_list_bnd:internal_boundaries_local] cg%i_bnd differs in number of entries from cg%o_bnd")
                  do g = lbound(cg%i_bnd(d)%seg(:), dim=1), ubound(cg%i_bnd(d)%seg(:), dim=1)

                     if (associated(cg%i_bnd(d)%seg(g)%local)) then
                        i_seg => cg%i_bnd(d)%seg(g)
                        ! find the right segment on the other grid container. OPT: can be done while searching for the pointer i_seg%local
                        g_o = INVALID
                        do i = lbound(i_seg%local%o_bnd(d)%seg, dim=1), ubound(i_seg%local%o_bnd(d)%seg, dim=1)
                           if (i_seg%tag == i_seg%local%o_bnd(d)%seg(i)%tag) then
                              g_o = i
                              exit
                           endif
                        enddo
                        if (g_o == INVALID) call die("[cg_list_bnd:internal_boundaries_local] cannot find the other grid chunk")
                        o_seg => i_seg%local%o_bnd(d)%seg(g_o)
                        if (tgt3d) then
                           pa3d   =>          cg%q(ind)%span(i_seg%se(:,:))
                           pa3d_o => i_seg%local%q(ind)%span(o_seg%se(:,:))
                           pa3d(:,:,:) = pa3d_o(:,:,:)
                        else
                           ! BEWARE: manual optimisation ... but it works (at least in gfortran)
                           do j = o_seg%se(zdim, LO), o_seg%se(zdim, HI)
                              cg%w(ind)%arr(:, i_seg%se(xdim, LO):i_seg%se(xdim, HI), &
                                   &           i_seg%se(ydim, LO):i_seg%se(ydim, HI), &
                                   &           j - o_seg%se(zdim, LO) + i_seg%se(zdim, LO)) = &
                                   i_seg%local%w(ind)%arr(:, o_seg%se(xdim, LO):o_seg%se(xdim, HI), o_seg%se(ydim, LO):o_seg%se(ydim, HI), j)
                           enddo
                        endif
                     endif
                  enddo
               else
                  if (allocated(cg%o_bnd(d)%seg)) call die("[cg_list_bnd:internal_boundaries_local] cg%o_bnd without cg%i_bnd")
               endif
            endif
         enddo

         call cg%costs%stop(I_OTHER)
         cgl => cgl%nxt
      enddo

   end subroutine internal_boundaries_local
!>
!! \brief This routine exchanges guardcells with remote blocks for BND_MPI and BND_PER boundaries on rank-3 and
!! rank-4 arrays. Pieces of boundaries that have to be sent between same pair of processes are merged
!<

   subroutine internal_boundaries_MPI_merged(this, ind, tgt3d, dmask)

      use constants,        only: xdim, ydim, zdim, cor_dim, I_ONE, LO, HI
      use dataio_pub,       only: die
      use merge_segments,   only: IN, OUT
      use MPIF,             only: MPI_DOUBLE_PRECISION, MPI_COMM_WORLD
      use MPIFUN,           only: MPI_Irecv, MPI_Isend, MPI_Comm_dup, MPI_Comm_free
      use mpisetup,         only: FIRST, LAST, proc, err_mpi, req, req2, inflate_req, nproc
      use named_array_list, only: wna
      use ppp_mpi,          only: piernik_Waitall
#ifdef MPIF08
      use MPIF,             only: MPI_Comm
#endif /* MPIF08 */

      implicit none

      class(cg_list_bnd_t),             intent(inout) :: this  !< the list on which to perform the boundary exchange
      integer(kind=4),                  intent(in)    :: ind   !< index of cg%q(:) 3d array or cg%w(:) 4d array
      logical,                          intent(in)    :: tgt3d !< .true. for cg%q, .false. for cg%w
      logical, dimension(xdim:cor_dim), intent(in)    :: dmask !< .true. for the directions we want to exchange

      integer :: i
      integer(kind=4) :: p
      integer(kind=4) :: nr !< index of first free slot in req array
#ifdef MPIF08
      type(MPI_Comm)  :: ibmpi_comm
#else /* !MPIF08 */
      integer(kind=4) :: ibmpi_comm
#endif /* !MPIF08 */

      if (.not. this%ms%valid) call die("[cg_list_bnd:internal_boundaries_MPI_merged] this%ms%valid .eqv. .false.")

      call inflate_req(nproc, .true.)

      call MPI_Comm_dup(MPI_COMM_WORLD, ibmpi_comm, err_mpi)
      nr = 0
      do p = FIRST, LAST
         if (p /= proc) then
            call this%ms%sl(p, IN )%find_offsets(dmask)
            call this%ms%sl(p, OUT)%find_offsets(dmask)
            if (this%ms%sl(p, IN)%total_size /= this%ms%sl(p, OUT)%total_size) &
                 call die("[cg_list_bnd:internal_boundaries_MPI_merged] this%ms%sl(p, :)%total_size /=")

            if (this%ms%sl(p, IN)%total_size /= 0) then ! we have something to communicate with process p
               if (tgt3d) then
                  allocate(this%ms%sl(p, IN )%buf(this%ms%sl(p, IN )%total_size))
                  allocate(this%ms%sl(p, OUT)%buf(this%ms%sl(p, OUT)%total_size))
                  do i = lbound(this%ms%sl(p, OUT)%list, dim=1), this%ms%sl(p, OUT)%cur_last
                     if (dmask( this%ms%sl(p, OUT)%list(i)%dir)) then
                        this     %ms%sl(p, OUT)%buf( &
                             this%ms%sl(p, OUT)%list(i)%offset: &
                             this%ms%sl(p, OUT)%list(i)%off_ceil) = reshape( &
                             this%ms%sl(p, OUT)%list(i)%cg%q(ind)%arr( &
                             this%ms%sl(p, OUT)%list(i)%se(xdim, LO): &
                             this%ms%sl(p, OUT)%list(i)%se(xdim, HI), &
                             this%ms%sl(p, OUT)%list(i)%se(ydim, LO): &
                             this%ms%sl(p, OUT)%list(i)%se(ydim, HI), &
                             this%ms%sl(p, OUT)%list(i)%se(zdim, LO): &
                             this%ms%sl(p, OUT)%list(i)%se(zdim, HI)), [ &
                             this%ms%sl(p, OUT)%list(i)%off_ceil - &
                             this%ms%sl(p, OUT)%list(i)%offset + I_ONE ] )
                     endif
                  enddo
               else
                  allocate(this%ms%sl(p, IN )%buf(this%ms%sl(p, IN )%total_size*wna%lst(ind)%dim4))
                  allocate(this%ms%sl(p, OUT)%buf(this%ms%sl(p, OUT)%total_size*wna%lst(ind)%dim4))
                  do i = lbound(this%ms%sl(p, OUT)%list, dim=1), this%ms%sl(p, OUT)%cur_last
                     if (dmask( this%ms%sl(p, OUT)%list(i)%dir)) then
                        this     %ms%sl(p, OUT)%buf( &
                            (this%ms%sl(p, OUT)%list(i)%offset - I_ONE) * wna%lst(ind)%dim4 + I_ONE : &
                             this%ms%sl(p, OUT)%list(i)%off_ceil        * wna%lst(ind)%dim4 ) = reshape( &
                             this%ms%sl(p, OUT)%list(i)%cg%w(ind)%arr( &
                             1:wna%lst(ind)%dim4, &
                             this%ms%sl(p, OUT)%list(i)%se(xdim, LO): &
                             this%ms%sl(p, OUT)%list(i)%se(xdim, HI), &
                             this%ms%sl(p, OUT)%list(i)%se(ydim, LO): &
                             this%ms%sl(p, OUT)%list(i)%se(ydim, HI), &
                             this%ms%sl(p, OUT)%list(i)%se(zdim, LO): &
                             this%ms%sl(p, OUT)%list(i)%se(zdim, HI)), [ &
                            (this%ms%sl(p, OUT)%list(i)%off_ceil - &
                             this%ms%sl(p, OUT)%list(i)%offset + I_ONE) * wna%lst(ind)%dim4 ] )
                     endif
                  enddo
               endif
               if (nr+I_ONE >  ubound(req(:), dim=1)) call inflate_req
               if (nr+I_ONE >  ubound(req2(:), dim=1)) call inflate_req(.true.)
               call MPI_Irecv(this%ms%sl(p, IN )%buf, size(this%ms%sl(p, IN )%buf, kind=4), MPI_DOUBLE_PRECISION, p, p,    ibmpi_comm, req( nr+I_ONE), err_mpi)
               call MPI_Isend(this%ms%sl(p, OUT)%buf, size(this%ms%sl(p, OUT)%buf, kind=4), MPI_DOUBLE_PRECISION, p, proc, ibmpi_comm, req2(nr+I_ONE), err_mpi)
               nr = nr + I_ONE
            endif

         endif
      enddo

      call piernik_Waitall(nr, "int_bnd_merged_R")

      do p = FIRST, LAST
         if (p /= proc) then
            if (this%ms%sl(p, IN)%total_size /= 0) then ! we have something received from process p
               if (tgt3d) then
                  do i = lbound(this%ms%sl(p, IN)%list, dim=1), this%ms%sl(p, IN)%cur_last
                     if (dmask( this%ms%sl(p, IN)%list(i)%dir)) then
                        this     %ms%sl(p, IN)%list(i)%cg%q(ind)%arr( &
                             this%ms%sl(p, IN)%list(i)%se(xdim, LO): &
                             this%ms%sl(p, IN)%list(i)%se(xdim, HI), &
                             this%ms%sl(p, IN)%list(i)%se(ydim, LO): &
                             this%ms%sl(p, IN)%list(i)%se(ydim, HI), &
                             this%ms%sl(p, IN)%list(i)%se(zdim, LO): &
                             this%ms%sl(p, IN)%list(i)%se(zdim, HI)) = reshape ( &
                             this%ms%sl(p, IN)%buf( &
                             this%ms%sl(p, IN)%list(i)%offset: &
                             this%ms%sl(p, IN)%list(i)%off_ceil), [ &
                             this%ms%sl(p, IN)%list(i)%se(xdim, HI) - &
                             this%ms%sl(p, IN)%list(i)%se(xdim, LO) + I_ONE, &
                             this%ms%sl(p, IN)%list(i)%se(ydim, HI) - &
                             this%ms%sl(p, IN)%list(i)%se(ydim, LO) + I_ONE, &
                             this%ms%sl(p, IN)%list(i)%se(zdim, HI) - &
                             this%ms%sl(p, IN)%list(i)%se(zdim, LO) + I_ONE ] )
                     endif
                  enddo
               else
                  do i = lbound(this%ms%sl(p, IN)%list, dim=1), this%ms%sl(p, IN)%cur_last
                     if (dmask( this%ms%sl(p, IN)%list(i)%dir)) then
                        this     %ms%sl(p, IN)%list(i)%cg%w(ind)%arr( &
                             1:wna%lst(ind)%dim4, &
                             this%ms%sl(p, IN)%list(i)%se(xdim, LO): &
                             this%ms%sl(p, IN)%list(i)%se(xdim, HI), &
                             this%ms%sl(p, IN)%list(i)%se(ydim, LO): &
                             this%ms%sl(p, IN)%list(i)%se(ydim, HI), &
                             this%ms%sl(p, IN)%list(i)%se(zdim, LO): &
                             this%ms%sl(p, IN)%list(i)%se(zdim, HI)) = reshape ( &
                             this%ms%sl(p, IN)%buf( &
                            (this%ms%sl(p, IN)%list(i)%offset - I_ONE) * wna%lst(ind)%dim4 + I_ONE : &
                             this%ms%sl(p, IN)%list(i)%off_ceil * wna%lst(ind)%dim4 ), [ &
                             int(wna%lst(ind)%dim4, kind=8), &
                             this%ms%sl(p, IN)%list(i)%se(xdim, HI) - &
                             this%ms%sl(p, IN)%list(i)%se(xdim, LO) + I_ONE, &
                             this%ms%sl(p, IN)%list(i)%se(ydim, HI) - &
                             this%ms%sl(p, IN)%list(i)%se(ydim, LO) + I_ONE, &
                             this%ms%sl(p, IN)%list(i)%se(zdim, HI) - &
                             this%ms%sl(p, IN)%list(i)%se(zdim, LO) + I_ONE ] )
                     endif
                  enddo
               endif
            endif
         endif
      enddo

      call piernik_Waitall(nr, "int_bnd_merged_S", use_req2 = .true.)
      call MPI_Comm_free(ibmpi_comm, err_mpi)

      do p = FIRST, LAST
         if (p /= proc) then
            if (allocated(this%ms%sl(p, IN )%buf)) deallocate(this%ms%sl(p, IN )%buf)
            if (allocated(this%ms%sl(p, OUT)%buf)) deallocate(this%ms%sl(p, OUT)%buf)
         endif
      enddo

   end subroutine internal_boundaries_MPI_merged

!>
!! \brief This routine exchanges guardcells with remote blocks for BND_MPI and BND_PER boundaries on rank-3 and
!! rank-4 arrays. There is one message per each piece of boundary.
!!
!! \details This routine will exchange local blocks as well (which would degrade the performance a bit) where the
!! pointer in cg%i_bnd(:)%seg(:)%local is not set.
!<

   subroutine internal_boundaries_MPI_1by1(this, ind, tgt3d, dmask)

      use cg_cost_data,     only: I_OTHER
      use cg_list,          only: cg_list_element
      use constants,        only: xdim, cor_dim, LO, HI, I_ONE, I_TWO, I_THREE, I_FOUR
      use dataio_pub,       only: die
      use grid_cont,        only: grid_container
      use grid_cont_bnd,    only: segment
      use MPIF,             only: MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, MPI_ORDER_FORTRAN, &
           &                      MPI_Type_create_subarray, MPI_Type_commit, MPI_Type_free
      use MPIFUN,           only: MPI_Irecv, MPI_Isend, MPI_Comm_dup, MPI_Comm_free
      use mpisetup,         only: err_mpi, req, inflate_req
      use named_array_list, only: wna
      use ppp_mpi,          only: piernik_Waitall
#ifdef MPIF08
      use MPIF,             only: MPI_Comm
#endif /* MPIF08 */

      implicit none

      class(cg_list_bnd_t),             intent(in) :: this  !< the list on which to perform the boundary exchange
      integer(kind=4),                  intent(in) :: ind   !< index of cg%q(:) 3d array or cg%w(:) 4d array
      logical,                          intent(in) :: tgt3d !< .true. for cg%q, .false. for cg%w
      logical, dimension(xdim:cor_dim), intent(in) :: dmask !< .true. for the directions we want to exchange

      integer                                      :: g, d
      integer(kind=4)                              :: nr     !< index of first free slot in req array
      type(grid_container),     pointer            :: cg
      type(cg_list_element),    pointer            :: cgl
      type(segment), pointer                       :: i_seg, o_seg !< shortcuts

      integer(kind=4), parameter :: rank3 = I_THREE, rank4 = I_FOUR
      integer(kind=4), dimension(rank3) :: b3sz, b3su, b3st
      integer(kind=4), dimension(rank4) :: b4sz, b4su, b4st
#ifdef MPIF08
      type(MPI_Comm)  :: ib1by1_comm
#else /* !MPIF08 */
      integer(kind=4) :: ib1by1_comm
#endif /* !MPIF08 */

      call MPI_Comm_dup(MPI_COMM_WORLD, ib1by1_comm, err_mpi)
      nr = 0
      cgl => this%first
      do while (associated(cgl))
         cg => cgl%cg
         call cg%costs%start

         ! exclude non-multigrid variables below base level
         if (tgt3d) then
            if (ind > ubound(cg%q(:), dim=1) .or. ind < lbound(cg%q(:), dim=1)) call die("[cg_list_bnd:internal_boundaries_MPI_1by1] wrong 3d index")
            b3sz = shape(cg%q(ind)%arr, kind=4)
         else
            if (ind > ubound(cg%w(:), dim=1) .or. ind < lbound(cg%w(:), dim=1)) call die("[cg_list_bnd:internal_boundaries_MPI_1by1] wrong 4d index")
            b4sz = shape(cg%w(ind)%arr, kind=4)
         endif

         do d = lbound(cg%i_bnd, dim=1), ubound(cg%i_bnd, dim=1)
            if (dmask(d) .and. is_active(cg, ind, tgt3d)) then
               if (allocated(cg%i_bnd(d)%seg)) then
                  if (.not. allocated(cg%o_bnd(d)%seg)) call die("[cg_list_bnd:internal_boundaries_MPI_1by1] cg%i_bnd without cg%o_bnd")
                  if (ubound(cg%i_bnd(d)%seg(:), dim=1) /= ubound(cg%o_bnd(d)%seg(:), dim=1)) &
                       call die("[cg_list_bnd:internal_boundaries_MPI_1by1] cg%i_bnd differs in number of entries from cg%o_bnd")
                  do g = lbound(cg%i_bnd(d)%seg(:), dim=1), ubound(cg%i_bnd(d)%seg(:), dim=1)

                     if (nr+I_TWO > ubound(req(:), dim=1)) call inflate_req

                     i_seg => cg%i_bnd(d)%seg(g)
                     o_seg => cg%o_bnd(d)%seg(g)

                     !> \deprecated: A lot of semi-duplicated code below
                     ! array_of_starts has to be C-like, so b3st(:) = 0  points to lbound(cg%q(ind)%arr)
                     if (tgt3d) then

                        b3su = int(i_seg%se(:, HI) - i_seg%se(:, LO) + I_ONE, kind=4)
                        b3st = int(i_seg%se(:, LO), kind=4) - lbound(cg%q(ind)%arr, kind=4)
                        call MPI_Type_create_subarray(rank3, b3sz, b3su, b3st, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, i_seg%sub_type, err_mpi)
                        call MPI_Type_commit(i_seg%sub_type, err_mpi)
                        call MPI_Irecv(cg%q(ind)%arr(:,:,:), I_ONE, i_seg%sub_type, i_seg%proc, i_seg%tag, ib1by1_comm, req(nr+I_ONE), err_mpi)

                        b3su = int(o_seg%se(:, HI) - o_seg%se(:, LO) + I_ONE, kind=4)
                        b3st = int(o_seg%se(:, LO), kind=4) - lbound(cg%q(ind)%arr, kind=4)
                        call MPI_Type_create_subarray(rank3, b3sz, b3su, b3st, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, o_seg%sub_type, err_mpi)
                        call MPI_Type_commit(o_seg%sub_type, err_mpi)
                        call MPI_Isend(cg%q(ind)%arr(:,:,:), I_ONE, o_seg%sub_type, o_seg%proc, o_seg%tag, ib1by1_comm, req(nr+I_TWO), err_mpi)

                     else

                        b4su = [ int(wna%lst(ind)%dim4, kind=4), int(i_seg%se(:, HI) - i_seg%se(:, LO) + I_ONE, kind=4) ]
                        b4st = [ I_ONE, int(i_seg%se(:, LO), kind=4) ] - lbound(cg%w(ind)%arr, kind=4)
                        call MPI_Type_create_subarray(rank4, b4sz, b4su, b4st, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, i_seg%sub_type, err_mpi)
                        call MPI_Type_commit(i_seg%sub_type, err_mpi)
                        call MPI_Irecv(cg%w(ind)%arr(:,:,:,:), I_ONE, i_seg%sub_type, i_seg%proc, i_seg%tag, ib1by1_comm, req(nr+I_ONE), err_mpi)

                        b4su = [ int(wna%lst(ind)%dim4, kind=4), int(o_seg%se(:, HI) - o_seg%se(:, LO) + I_ONE, kind=4) ]
                        b4st = [ I_ONE, int(o_seg%se(:, LO), kind=4) ] - lbound(cg%w(ind)%arr, kind=4)
                        call MPI_Type_create_subarray(rank4, b4sz, b4su, b4st, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, o_seg%sub_type, err_mpi)
                        call MPI_Type_commit(o_seg%sub_type, err_mpi)
                        call MPI_Isend(cg%w(ind)%arr(:,:,:,:), I_ONE, o_seg%sub_type, o_seg%proc, o_seg%tag, ib1by1_comm, req(nr+I_TWO), err_mpi)

                     endif
                     nr = nr + I_TWO
                  enddo
               else
                  if (allocated(cg%o_bnd(d)%seg)) call die("[cg_list_bnd:internal_boundaries_MPI_1by1] cg%o_bnd without cg%i_bnd")
               endif
            endif
         enddo

         call cg%costs%stop(I_OTHER)
         cgl => cgl%nxt
      enddo

      call piernik_Waitall(nr, "int_bnd_1by1")
      call MPI_Comm_free(ib1by1_comm, err_mpi)

      cgl => this%first
      do while (associated(cgl))
         cg => cgl%cg
         call cg%costs%start

         do d = lbound(cg%i_bnd, dim=1), ubound(cg%i_bnd, dim=1)
            if (dmask(d) .and. is_active(cg, ind, tgt3d)) then
               if (allocated(cg%i_bnd(d)%seg)) then
                  ! sanity checks are already done
                  do g = lbound(cg%i_bnd(d)%seg(:), dim=1), ubound(cg%i_bnd(d)%seg(:), dim=1)
                     call MPI_Type_free(cg%i_bnd(d)%seg(g)%sub_type, err_mpi)
                     call MPI_Type_free(cg%o_bnd(d)%seg(g)%sub_type, err_mpi)
                  enddo
               endif
            endif
         enddo

         call cg%costs%stop(I_OTHER)
         cgl => cgl%nxt
      enddo

   end subroutine internal_boundaries_MPI_1by1

!> \brief exclude non-multigrid variables below base level

   pure logical function is_active(cg, ind, tgt3d)

      use grid_cont, only: grid_container

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer(kind=4),               intent(in) :: ind   !< index of cg%q(:) 3d array or cg%w(:) 4d array
      logical,                       intent(in) :: tgt3d !< .true. for cg%q, .false. for cg%w

      if (tgt3d) then ! cannot use merge() here
         is_active = associated(cg%q(ind)%arr)
      else
         is_active = associated(cg%w(ind)%arr)
      endif

   end function is_active

!> \brief Set zero to all boundaries (will defeat any attempts of use of dirty checks on boundaries)

   subroutine clear_boundaries(this, ind, value)

      use cg_list,   only: cg_list_element
      use constants, only: ndims, xdim, zdim, LO, HI, BND_MPI, BND_FC, BND_MPI_FC
      use domain,    only: dom
      use grid_cont, only: grid_container

      implicit none

      class(cg_list_bnd_t), intent(in)        :: this  !< the list on which to perform the action
      integer(kind=4),      intent(in)        :: ind   !< Negative value: index of cg%q(:) 3d array
      real, optional,       intent(in)        :: value !< Value to be put in the boundaries (could be dirty)

      integer(kind=4)                         :: lh, clh, d
      integer(kind=4), dimension(ndims,LO:HI) :: l
      type(cg_list_element),  pointer         :: cgl
      type(grid_container),   pointer         :: cg
      real, dimension(:,:,:), pointer         :: pa3d
      real :: v

      !> \todo fill corners with sth*dirtyH1 ?

      v = 0.
      if (present(value)) v=value
      cgl => this%first
      do while (associated(cgl))
         cg => cgl%cg
         do d = xdim, zdim
            if (dom%has_dir(d)) then
               do lh = LO, HI
                  if (any(cg%bnd(d, lh) == [ BND_MPI, BND_FC, BND_MPI_FC ])) then !> \warning there should be no exceptions but something in crtest still needs this
                     !> \todo Find out what is the problem with crtest
                     l = reshape([lbound(cg%q(ind)%arr, kind=4),ubound(cg%q(ind)%arr, kind=4)],shape=[ndims,HI])
                     clh = LO + HI - lh
                     l(d,clh) = cg%ijkse(d,lh) + HI*lh-(LO+HI) ! -1 for LO, +1 for HI
                     pa3d => cg%q(ind)%span(l)
                     pa3d = v
                  endif
               enddo
            endif
         enddo
         cgl => cgl%nxt
      enddo

   end subroutine clear_boundaries

!> \brief Put dirty values to all boundaries

   subroutine dirty_boundaries(this, ind)

      use constants, only: dirtyH1

      implicit none

      class(cg_list_bnd_t), intent(in) :: this  !< the list on which to perform the action
      integer(kind=4),      intent(in) :: ind   !< Negative value: index of cg%q(:) 3d array

      call this%clear_boundaries(ind, value=0.87*dirtyH1)

   end subroutine dirty_boundaries

!>
!! \brief Set up external boundary values
!!
!! \details This is purely local operation. We don't bother here to limit the number of layers to be filled, because it is cheap.
!! Can be moved to cg_list or cg_list_dataop.
!<

   subroutine external_boundaries(this, ind, area_type, bnd_type)

      use cg_cost_data, only: I_OTHER
      use cg_list,      only: cg_list_element
      use constants,    only: ndims, xdim, ydim, zdim, LO, HI, AT_NO_B, I_ONE, I_TWO, I_THREE, &
           &                  BND_PER, BND_MPI, BND_FC, BND_MPI_FC, BND_SHE, BND_COR, BND_REF, BND_NEGREF, BND_ZERO, BND_XTRAP, BND_NONE
      use dataio_pub,   only: die
      use domain,       only: dom
      use grid_cont,    only: grid_container

      implicit none

      class(cg_list_bnd_t),      intent(in)   :: this       !< the list on which to perform the boundary exchange
      integer(kind=4),           intent(in)   :: ind        !< Negative value: index of cg%q(:) 3d array
      integer(kind=4), optional, intent(in)   :: area_type  !< defines how do we treat boundaries
      integer(kind=4), optional, intent(in)   :: bnd_type   !< Override default boundary type on external boundaries (useful in multigrid solver).
                                                            !< Note that BND_PER, BND_MPI, BND_SHE and BND_COR aren't external and cannot be overridden

      integer(kind=4)                         :: lh, clh, d, b_type, i
      integer(kind=4), dimension(ndims,LO:HI) :: l, r, rh
      type(cg_list_element),  pointer         :: cgl
      type(grid_container),   pointer         :: cg
      real, dimension(:,:,:), pointer         :: pa3d

      cgl => this%first
      do while (associated(cgl))
         cg => cgl%cg
         call cg%costs%start

         if (ind > ubound(cg%q(:), dim=1) .or. ind < lbound(cg%q(:), dim=1)) call die("[cg_list_bnd:external_boundaries] wrong 3d index")
         pa3d => cg%q(ind)%arr

         do d = xdim, zdim
            if (dom%has_dir(d)) then
               l = reshape([lbound(pa3d, kind=4),ubound(pa3d, kind=4)],shape=[ndims,HI]) ; r = l
               do lh = LO, HI

                  select case (cg%bnd(d, lh))
                     case (BND_PER, BND_MPI, BND_FC, BND_MPI_FC) ! Already done in internal_bnd or arr3d_boundaries
                     case (BND_SHE) !> die until someone really needs SHEAR.
                        call die("[cg_list_bnd:external_boundaries] 'she' not implemented")
                     case (BND_COR)
                        if (present(area_type)) then
                           if (area_type /= AT_NO_B) cycle
                        endif
                        call die("[cg_list_bnd:external_boundaries] 'cor' not implemented")
                     case default ! Set gradient == 0 on the external boundaries
                        if (present(area_type)) then
                           if (area_type /= AT_NO_B) cycle
                        endif
                        b_type = cg%bnd(d, lh)
                        if (present(bnd_type)) b_type = bnd_type
                        select case (b_type)
                           case (BND_REF)  ! reflecting BC (e.g. homogeneous Neumann)
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
                           case (BND_ZERO)  ! zero BC (e.g. homogeneous Dirichlet BC with 0 at first layer of cells)
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

         call cg%costs%stop(I_OTHER)
         cgl => cgl%nxt
      enddo

   end subroutine external_boundaries

!> \brief External (Non-MPI) boundary conditions for the fluid array: cg%u

   subroutine bnd_u(this, dir)

      use cg_cost_data,          only: I_OTHER
      use cg_list,               only: cg_list_element
      use constants,             only: ndims, xdim, ydim, zdim, LO, HI, INT4, I_ONE, &
           &                           BND_MPI, BND_FC, BND_MPI_FC, BND_PER, BND_REF, BND_OUT, BND_OUTD, BND_COR, BND_SHE, BND_USER
      use dataio_pub,            only: msg, warn, die
      use domain,                only: dom, vel_outd
      use fluidboundaries_funcs, only: user_fluidbnd
      use fluidindex,            only: iarr_all_dn
      use grid_cont,             only: grid_container
      use named_array_list,      only: wna
      use ppp,                   only: ppp_main
#ifdef COSM_RAYS
      use fluidindex,            only: iarr_all_crn
      use initcosmicrays,        only: smallecr
#endif /* COSM_RAYS */
#ifdef CRESP
      use initcrspectrum,        only: smallcree, smallcren
      use initcosmicrays,        only: iarr_cre_e, iarr_cre_n
#endif /* CRESP */
#ifdef GRAV
      use constants,             only: BND_OUTH, BND_OUTHD, I_ZERO
      use fluidboundaries_funcs, only: outh_fluidbnd
#endif /* GRAV */

      implicit none

      class(cg_list_bnd_t), intent(in) :: this !< the list on which to perform the boundary exchange
      integer(kind=4),      intent(in) :: dir  !< the direction in which we perform fluid boundary update (xdim, ydim or zdim)

      type(grid_container), pointer           :: cg
      integer(kind=4), dimension(ndims,LO:HI) :: l, r
      logical, save                           :: frun = .true.
      integer(kind=4)                         :: side, ssign, ib
      type(cg_list_element), pointer          :: cgl
      character(len=*), parameter             :: bu_label = "bnd_u"

      if (.not. any([xdim, ydim, zdim] == dir)) call die("[cg_list_bnd:bnd_u] Invalid direction.")

      call ppp_main%start(bu_label)

      if (frun) then
         call init_fluidboundaries
         frun = .false.
         if (HI-LO /= I_ONE) call die("[cg_list_bnd:bnd_u] HI-LO /= I_ONE")
      endif

      cgl => this%first
      do while (associated(cgl))
         cg => cgl%cg
         call cg%costs%start

         l = cg%lhn ; r = l
         do side = LO, HI

            select case (cg%bnd(dir, side))
               case (BND_MPI, BND_COR, BND_SHE, BND_FC, BND_MPI_FC, BND_PER)
                  ! Do nothing
               case (BND_USER)
                  call user_fluidbnd(dir, side, cg, wn=wna%fi)
               case (BND_REF)
                  ssign = lh_2_pm1(side)
                  do ib=1_INT4, dom%nb
                     l(dir,:) = cg%ijkse(dir,side)+ssign*ib ; r(dir,:) = cg%ijkse(dir,side)+ssign*(1_INT4-ib)
                     cg%u(:,l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) = cg%u(:,r(xdim,LO):r(xdim,HI),r(ydim,LO):r(ydim,HI),r(zdim,LO):r(zdim,HI))
                     cg%u(iarr_all_dn+dir, l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) = -cg%u(iarr_all_dn+dir,l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI))
                  enddo
               case (BND_OUT)
                  r(dir,:) = cg%ijkse(dir,side)
                  ssign = lh_2_pm1(side)
                  do ib=1_INT4, dom%nb
                     l(dir,:) = cg%ijkse(dir,side)+ssign*ib
                     cg%u(:,l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) = cg%u(:,r(xdim,LO):r(xdim,HI),r(ydim,LO):r(ydim,HI),r(zdim,LO):r(zdim,HI))
#ifdef COSM_RAYS
                     cg%u(iarr_all_crn, l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) = smallecr
#endif /* COSM_RAYS */
#ifdef CRESP
                     cg%u(iarr_cre_n,   l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) = smallcren   !< CRESP number density component
                     cg%u(iarr_cre_e,   l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) = smallcree   !< CRESP energy density component
#endif /* CRESP */
                  enddo
               case (BND_OUTD)
                  r(dir,:) = cg%ijkse(dir,side)
                  ssign = lh_2_pm1(side)
                  do ib=1_INT4, dom%nb
                     l(dir,:) = cg%ijkse(dir,side)+ssign*ib
                     cg%u(:,l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) = cg%u(:,r(xdim,LO):r(xdim,HI),r(ydim,LO):r(ydim,HI),r(zdim,LO):r(zdim,HI))
                     !> \deprecated BEWARE: use of uninitialized value on first call (a side effect of r1726)
#ifdef COSM_RAYS
                     cg%u(iarr_all_crn, l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) = smallecr
#endif /* COSM_RAYS */
#ifdef CRESP
                     cg%u(iarr_cre_n,   l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) = smallcren   !< CRESP number density component
                     cg%u(iarr_cre_e,   l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) = smallcree   !< CRESP energy density component
#endif /* CRESP */
                  enddo
                  l(dir,:) = cg%ijkse(dir,side) - [dom%nb, 1_INT4] +(dom%nb+1_INT4)*(side-LO)
                  if (side == LO) then
                     cg%u(iarr_all_dn+dir,l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) = min(cg%u(iarr_all_dn+dir,l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)), -vel_outd * cg%u(iarr_all_dn,l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)))
                  else
                     cg%u(iarr_all_dn+dir,l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) = max(cg%u(iarr_all_dn+dir,l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)),  vel_outd * cg%u(iarr_all_dn,l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)))
                  endif
#ifdef GRAV
               case (BND_OUTH)
                  call outh_fluidbnd(dir, side, cg, wn=I_ZERO)
               case (BND_OUTHD)
                  call outh_fluidbnd(dir, side, cg, wn=I_ONE)
#endif /* GRAV */
               case default
                  write(msg,'("[cg_list_bnd:bnd_u]: Unrecognized ",i1," boundary condition ",i3," not implemented in ",i1,"-direction")') side, cg%bnd(dir, side), dir
                  call warn(msg)
            end select

         enddo

         call cg%costs%stop(I_OTHER)
         cgl => cgl%nxt
      enddo

      call ppp_main%stop(bu_label)

   contains

!> \brief convert [ LO, HI ] to [ -1, 1 ]. Assumes HI-LO == 1

      elemental function lh_2_pm1(lh)

         use constants, only: LO, HI, I_TWO

         implicit none

         integer(kind=4), intent(in) :: lh

         integer(kind=4) :: lh_2_pm1

         lh_2_pm1 = I_TWO*lh - (LO+HI)

      end function lh_2_pm1

!> \brief Perform some checks

      subroutine init_fluidboundaries

         use constants,  only: PIERNIK_INIT_DOMAIN, xdim, zdim, LO, HI, &
              &                BND_MPI, BND_FC, BND_MPI_FC, BND_PER, BND_REF, BND_OUT, BND_OUTD, BND_OUTH, BND_OUTHD, BND_COR, BND_SHE, BND_USER
         use dataio_pub, only: msg, warn, die, code_progress
         use domain,     only: dom, is_multicg

         implicit none

         integer(kind=4) :: dir, side

         if (code_progress < PIERNIK_INIT_DOMAIN) call die("[cg_list_bnd:init_fluidboundaries] MPI not initialized.") ! bnd_xl, bnd_xr

         do dir = xdim, zdim
            do side = LO, HI

               select case (dom%bnd(dir, side))
                  case (BND_MPI, BND_REF, BND_OUT, BND_OUTD, BND_USER, BND_PER)
                     ! Do nothing
                  case (BND_FC, BND_MPI_FC)
                     call die("[cg_list_bnd:init_fluidboundaries] fine-coarse interfaces not implemented yet")
                  case (BND_COR)
                     if (dir == zdim) then
                        write(msg,'("[cg_list_bnd:init_fluidboundaries] corner ",i1," boundary condition ",i3," not implemented in ",i1,"-direction")') &
                             side, dom%bnd(dir, side), dir
                        call warn(msg)
                     endif
                  case (BND_SHE)
                     if (dir /= xdim) then
                        write(msg,'("[cg_list_bnd:init_fluidboundaries] shear ",i1," boundary condition ",i3," not implemented in ",i1,"-direction")') &
                             side, dom%bnd(dir, side), dir
                        call warn(msg)
                     endif
                  case (BND_OUTH)
                     if (dir == zdim) then
                        if (is_multicg) call die("[cg_list_bnd:init_fluidboundaries] hydrostatic:outh_bnd with multiple grid pieces per processor not implemented yet")
                        !nontrivial not really checked
                     else
                        write(msg,'("[cg_list_bnd:init_fluidboundaries] outflow hydrostatic ",i1," boundary condition ",i3," not implemented in ",i1,"-direction")') &
                             side, dom%bnd(dir, side), dir
                        call warn(msg)
                     endif
                  case (BND_OUTHD)
                     if (dir == zdim) then
                        if (is_multicg) call die("[cg_list_bnd:init_fluidboundaries] hydrostatic:outh_bnd with multiple grid pieces per processor not implemented yet")
                        !nontrivial not really checked
                     else
                        write(msg,'("[cg_list_bnd:init_fluidboundaries] outflow hydrostatic ",i1," boundary condition ",i3," not implemented in ",i1,"-direction")') &
                             side, dom%bnd(dir, side), dir
                        call warn(msg)
                     endif
                  case default
                     write(msg,'("[cg_list_bnd:init_fluidboundaries] unknown ",i1," boundary condition ",i3," not implemented in ",i1,"-direction")') &
                          side, dom%bnd(dir, side), dir
                     call warn(msg)
               end select
            enddo
         enddo

      end subroutine init_fluidboundaries

   end subroutine bnd_u

!> \brief External (Non-MPI) boundary conditions for the magnetic field array: cg%b

   subroutine bnd_b(this, dir)

      use cg_cost_data,          only: I_OTHER
      use cg_list,               only: cg_list_element
      use constants,             only: ndims, xdim, ydim, zdim, LO, HI, I_TWO, I_THREE, &
           &                           BND_MPI, BND_FC, BND_MPI_FC, BND_PER, BND_REF, BND_OUT, BND_OUTD, BND_OUTH, BND_OUTHD, BND_COR, BND_SHE, BND_USER
      use dataio_pub,            only: msg, warn, die
      use domain,                only: dom
      use fluidboundaries_funcs, only: user_fluidbnd
      use global,                only: cc_mag
      use grid_cont,             only: grid_container
      use mpisetup,              only: master
      use named_array_list,      only: wna
      use ppp,                   only: ppp_main

      implicit none

      class(cg_list_bnd_t), intent(in) :: this !< the list on which to perform the boundary exchange
      integer(kind=4),      intent(in) :: dir  !< the direction in which we perform magnetic boundary update (xdim, ydim or zdim)

      type(grid_container), pointer           :: cg
      integer(kind=4)                         :: side
      logical, save                           :: frun = .true.
      logical, save,   dimension(ndims,LO:HI) :: bnd_not_provided = .false.
      integer(kind=4), dimension(ndims,LO:HI) :: l, r
      type(cg_list_element), pointer          :: cgl
      character(len=*), parameter             :: bb_label = "bnd_b"

      call ppp_main%start(bb_label)

! Non-MPI boundary conditions
      if (frun) then
         !> \deprecated This may not work as intended when there are many grid containers per process. PLEASE CHECK IT
         bnd_not_provided(:,         :) = (dom%bnd(:,:) == BND_REF)      .or. (dom%bnd(:,         :) == BND_MPI)
         bnd_not_provided(xdim:ydim, :) = bnd_not_provided(xdim:ydim, :) .or. (dom%bnd(xdim:ydim, :) == BND_COR)
         bnd_not_provided(xdim,      :) = bnd_not_provided(xdim,      :) .or. (dom%bnd(xdim,      :) == BND_SHE)
      endif

      if (bnd_not_provided(dir,LO) .and. bnd_not_provided(dir,HI)) return  ! avoid triple case

      cgl => this%first
      do while (associated(cgl))
         cg => cgl%cg
         call cg%costs%start

         l = cg%lhn ; r = l

         do side = LO, HI
            select case (cg%bnd(dir, side))
               case (BND_MPI, BND_REF)
                  ! Do nothing
               case (BND_USER)
                  call user_fluidbnd(dir, side, cg, wn=wna%bi)
               case (BND_FC, BND_MPI_FC)
                  if (.not. cc_mag) &
                       call die("[cg_list_bnd:bnd_b] fine-coarse interfaces not implemented yet for face-centered B field.")
               case (BND_COR)
                  if (dir == zdim) then
                     write(msg,'(2(a,i3))') "[cg_list_bnd:bnd_b]: Boundary condition ",cg%bnd(dir, side)," not implemented in ",dir
                     if (master) call warn(msg)
                  endif
               case (BND_SHE)
                  if (dir /= xdim) then
                     write(msg,'(2(a,i3))') "[cg_list_bnd:bnd_b]: Boundary condition ",cg%bnd(dir, side)," not implemented in ",dir
                     if (master) call warn(msg)
                  endif
               case (BND_PER)
               case (BND_OUT, BND_OUTD, BND_OUTH, BND_OUTHD)
                  call outflow_b(cg, dir, side)
               case default
                  write(msg,'(2(a,i3))') "[cg_list_bnd:bnd_b]: Boundary condition ",cg%bnd(dir, side)," not implemented in ",dir
                  if (master) call warn(msg)
            end select
         enddo

         call cg%costs%stop(I_OTHER)
         cgl => cgl%nxt
      enddo

      call ppp_main%stop(bb_label)

   contains

      subroutine outflow_b(cg, dir, side)

         ! use global,                only: cc_mag
         use grid_cont,             only: grid_container

         implicit none

         type(grid_container), pointer    :: cg
         integer(kind=4),      intent(in) :: dir
         integer(kind=4),      intent(in) :: side

         integer :: i, it
         integer :: pm_one   !< +1 for LO and -1 for HI
         integer :: pm_two   !< +2 for LO and -2 for HI

         pm_one = I_THREE - I_TWO * side
         pm_two = 2 * pm_one

         ! Apparently this is already written for cell-centered magnetic field.

         ! Simulations with Constrained Transport may exhibit slight asymmetries because
         ! rightmost face is reset here while leftmost is not. Use expressions like
         !
         !   it = cg%ijkse(dir, side) - pm_one * i + (side - LO)
         !
         ! when cc_mag is .false. in evaluation of dir-component of magnetic field
         ! for more strict external boundary treatment.

         ! BEWARE: this kind of boundaries does not guarantee div(B) == 0 .
         ! Expect div(B) growing proportionally to the distance from the domain boundary.

         select case (dir)
            case (xdim)
               do i = 1, dom%nb
                  it = cg%ijkse(dir, side) - pm_one * i
                  cg%b(xdim, it, :, :) = 2.0 * cg%b(xdim, it + pm_one, :, :) - cg%b(xdim, it + pm_two, :, :)
                  cg%b(ydim, it, :, :) = cg%b(ydim, it + pm_one, :, :)
                  cg%b(zdim, it, :, :) = cg%b(zdim, it + pm_one, :, :)
               enddo
            case (ydim)
               do i = 1, dom%nb
                  it = cg%ijkse(dir, side) - pm_one * i
                  cg%b(ydim, :, it, :) = 2.0 * cg%b(ydim, :, it + pm_one, :) - cg%b(ydim, :, it + pm_two, :)
                  cg%b(xdim, :, it, :) = cg%b(xdim, :, it + pm_one, :)
                  cg%b(zdim, :, it, :) = cg%b(zdim, :, it + pm_one, :)
               enddo
            case (zdim)
               do i = 1, dom%nb
                  it = cg%ijkse(dir, side) - pm_one * i
                  cg%b(zdim, :, :, it) = 2.0 * cg%b(zdim, :, :, it + pm_one) - cg%b(zdim, :, :, it + pm_two)
                  cg%b(xdim, :, :, it) = cg%b(xdim, :, :, it + pm_one)
                  cg%b(ydim, :, :, it) = cg%b(ydim, :, :, it + pm_one)
               enddo
         end select

      end subroutine outflow_b

   end subroutine bnd_b

end module cg_list_bnd
