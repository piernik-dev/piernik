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

!> \brief A module with an abstract type created to take out neighbor finding code from cg_level

module cg_list_neighbors

   use cg_list_rebalance, only: cg_list_rebalance_T

   implicit none

   private
   public :: cg_list_neighbors_T

   !>
   !! \brief An abstract type created to take out neighbor finding code from cg_level
   !!
   !! \details
   !! OPT: Searching through this%dot%gse for neighbours, prolongation/restriction overlaps etc is quite costly.
   !! The cost is O(this%cnt^2). Provide a list, sorted according to Morton/Hilbert id's and do a bisection search
   !! instead of checking against all grids on AMR-style 'blocky' grids or cartesian decomposition with equal-sized
   !! grids. It will result in massive speedups on cg_list_neighbors_T%find_neighbors and
   !! cg_level_connected_T%{vertical_prep,vertical_b_prep). It may also simplify the process of fixing refinement
   !! structure in refinement_update::fix_refinement. Grids which are larger than AMR_bsize (merged grids, non-block
   !! decompositions, both not implemented yet) may be referred by several id's that correspond with AMR_bsize-d virtual
   !! grid pieces. Non-cartesian decompositions should be handled with the bruteforce way. It is possible to optimize
   !! them slightly if we save the numbers found during decomposition, but I don't think it is really important.
   !! Unequal cartesian decompositions should be handled with the bruteforce way. It can be optimized too, but the
   !! impact of optimization would be similar to optimization of non-cartesian decompositions.
   !! Alternatively, construct a searchable binary tree or oct-tree and provide fast routines for searching grid
   !! pieces covering specified position.
   !<
   type, extends(cg_list_rebalance_T), abstract :: cg_list_neighbors_T
   contains
      procedure          :: find_neighbors            !< Choose between more general nad fast routine for neighbor searching
      procedure, private :: find_neighbors_SFC        !< Make full description of intra-level communication with neighbors. Approach exploiting strict SFC distribution.
      procedure, private :: find_neighbors_bruteforce !< Make full description of intra-level communication with neighbors. Brute-force approach.
   end type cg_list_neighbors_T

contains

!>
!! \brief Choose between more general and fast routine for neighbor searching
!!
!! \details
!!
!! Possible improvements of performance
!! * merge smaller blocks into larger ones,
!!
!! \todo Write variant of find_neighbors_* routine to achieve previous (pre-a27c945a) performance and maintain
!! correctness on corners on complicated topologies:
!! * Divide the descriptions of communicated regions into 4 categories: X-faces, Y-faces + XY-corners, Z-faces +
!!   [XY]Z-corners, other corners. The other corners would be non-empty only for some refinement local topologies,
!!   it would certainly be empty on an uniform grid.
!! * When no corners are required, perform simultaneous exchange described by the three directional categories.
!!   Some corners might be set up correctly by a chance, some might not.
!! * When corners are required, perform sequential exchange described by the three directional categories and
!!   supplement it with communication of "other corners". The sequence of Isend/Irecv should be as follows: Isend
!!   X-faces, Irecv X-faces, Waitall, Isend Y-faces, Irecv Y-faces, Waitall, Isend Z-faces, Irecv Z-faces, Waitall
!!   "Other corners" can be Isend at any time and must be Irecv after Z-faces are copied to the right place.
!!
!! \todo Provide a way to merge single-boundary messages into large clumps. On AMR grids it may outperform current
!! cartesian decompositions of uniform grids and may get close to pre-a27c945a performance.
!<

   subroutine find_neighbors(this)

      use refinement, only: prefer_n_bruteforce

      implicit none

      class(cg_list_neighbors_T), intent(inout) :: this !< object invoking type bound procedure

      this%ms%valid = .false.
      if (this%dot%is_blocky .and. .not. prefer_n_bruteforce) then
         call this%find_neighbors_SFC
         call this%ms%merge(this)
      else
         call this%find_neighbors_bruteforce
         ! calling this%ms%merge(this) here makes sense only for such setups, so in periodic boundaries 2 blocks
         ! covers the domain in at least one direction.
      endif

   end subroutine find_neighbors

!>
!! \brief Make full description of intra-level communication with neighbors. Approach exploiting strict SFC
!! distribution.
!!
!! \details Assume that cuboids don't collide (no overlapping grid pieces on same refinement level are allowed)
!! Should produce the same set of blocks to be communicated as find_neighbors_bruteforce, but should be way faster,
!! especially in massively parallel runs.
!!
!! This approach works best if the grid containers are distributing strictly according to the SFC curve and are sorted.
!! If each of p processes has g grid containers on current level (giving n = p * g grids on the level),
!! the cost should be proportional to (log_2(p)+log_2(g))*g, assuming that we have already sorted array
!! containing most critical information from this%dot%gse
!!
!! OPT: there is an easy way to determine here if the corner updates can be performed by face updates done in a
!! sequence of "sweeps" or not. It is possible even to isolate only those corners that can't be updated with face
!! communication. The question is: will it really matter for AMR communication, when we manage to aggregate separate
!! messages for pairs of processes to a single message? IMO not that much.
!!
!! OPT: To achieve pre-a27c945a performance in uniform-grid simulations we may also implement another variant of
!! find_neighbors_* that looks only for face neighbors (typically 6, in worst case 12, still far less than 26) for
!! levels that are fully covered (or even for levels that are covered by strictly convect, well separated groups of
!! grids).
!<

   subroutine find_neighbors_SFC(this)

      use cg_list,    only: cg_list_element
      use constants,  only: xdim, ydim, zdim, cor_dim, ndims, LO, HI, INVALID, BND_FC
      use dataio_pub, only: die
      use domain,     only: dom
      use gcpa,       only: gcpa_T
      use grid_cont,  only: grid_container
      use mpisetup,   only: proc
      use ordering,   only: SFC_order

      implicit none

      class(cg_list_neighbors_T), intent(inout) :: this !< object invoking type bound procedure

      type(grid_container),  pointer    :: cg      !< grid container that we are currently working on
      type(cg_list_element), pointer    :: cgl
      integer(kind=4)                   :: ix, iy, iz
      integer                           :: lh
      integer(kind=8), dimension(ndims) :: n_off     !< neighbor's offset
      integer(kind=8)                   :: n_id      !< neighbor's id
      integer                           :: n_p       !< neighbor's process
      integer                           :: n_grid_id !< neighbor's grid_id on n_p
      integer                           :: n_dd      !< neighbor's direction
      integer(kind=4)                   :: tag
      integer(kind=8), dimension(xdim:zdim, LO:HI) :: overlap
      type(gcpa_T) :: l_pse

      if (.not. this%dot%is_blocky) call die("[cg_list_neighbors:find_neighbors_SFC] Can work only on regular cartesian decompositions")

      call l_pse%init(this)

      cgl => this%first
      do while (associated(cgl))
         cg => cgl%cg

         if (allocated(cg%i_bnd)) deallocate(cg%i_bnd)
         if (allocated(cg%o_bnd)) deallocate(cg%o_bnd)
         allocate(cg%i_bnd(xdim:cor_dim), cg%o_bnd(xdim:cor_dim))

         ! for all potential neighbors:
         do iz = -dom%D_z, dom%D_z
            do iy = -dom%D_y, dom%D_y
               do ix = -dom%D_x, dom%D_x
                  if (any( [ ix, iy, iz ] /= 0)) then
                     ! find their SFC_id (take care about periodicity)
                     n_off = cg%my_se(:, LO) + [ ix, iy, iz ] * cg%n_b
                     where (dom%periodic) n_off = mod(n_off + this%n_d - this%off, this%n_d) + this%off
                     n_id = INVALID
                     if ( all(n_off >= this%off          .or. .not. dom%has_dir) .and. &
                          all(n_off <  this%off+this%n_d .or. .not. dom%has_dir)) then ! it is internal boundary

                        n_dd = INVALID
                        if (count([ix, iy, iz] /= 0) > 1) then
                           n_dd = cor_dim
                        else if (ix /= 0) then
                           n_dd = xdim
                        else if (iy /= 0) then
                           n_dd = ydim
                        else if (iz /= 0) then
                           n_dd = zdim
                        endif
                        if (n_dd == INVALID) call die("[cg_list_neighbors:find_neighbors_SFC] undefined direction")

                        n_id = SFC_order(n_off-this%off)
                        call this%dot%find_grid(n_id, n_p, n_grid_id) ! find on what process they may reside
                        if (n_grid_id == INVALID) then ! find if they really occur on that process
                           ! if it not occurs set cg%bnd(d, lh) to BND_FC or BND_MPI_FC
                           lh = LO + (ix + iy + iz + 1)/2 ! -1 => LO, +1 => HI, only valid for count([ix, iy, iz] /= 0) == 1) which implies n_dd /= cor_dim
                           if (n_dd >=xdim .and. n_dd <=zdim) cg%bnd(n_dd, lh) = BND_FC
                        else
                           ! incoming part:
                           tag = uniq_tag([-ix, -iy, -iz], n_grid_id)
                           overlap(:, LO) = max(cg%lhn(:, LO), cg%ijkse(:, LO) + [ ix, iy, iz ] * cg%n_b)
                           overlap(:, HI) = min(cg%lhn(:, HI), cg%ijkse(:, HI) + [ ix, iy, iz ] * cg%n_b)
                           call cg%i_bnd(n_dd)%add_seg(n_p, overlap, tag)
                           if (n_p == proc) cg%i_bnd(n_dd)%seg(ubound(cg%i_bnd(n_dd)%seg, dim=1))%local => l_pse%l_pse(n_grid_id)%p

                           ! outgoing part:
                           tag = uniq_tag([ix, iy, iz], cg%grid_id)
                           overlap(:, LO) = max(cg%ijkse(:, LO), cg%lhn(:, LO) + [ ix, iy, iz ] * cg%n_b)
                           overlap(:, HI) = min(cg%ijkse(:, HI), cg%lhn(:, HI) + [ ix, iy, iz ] * cg%n_b)
                           call cg%o_bnd(n_dd)%add_seg(n_p, overlap, tag)

                        endif
                     endif
                  endif
               enddo
            enddo
         enddo

         cgl => cgl%nxt
      enddo

      call l_pse%cleanup

   contains

      !>
      !! \brief Create unique tag for cg - cg exchange
      !!
      !! \details If we put a constraint that a grid piece can not be smaller than dom%nb, then total number of
      !! neighbours that affect local guardcells is exactly 3 for AMR, cartesian decomposition with equal-size
      !! blocks
      !! Thus, in each direction we can describe realtive position as one of three cases:
      !! * LEFT, RIGHT - corner neighbours
      !! * FACE - face neighbour
      !<
      pure function uniq_tag(ixyz, grid_id)

         use constants, only: xdim, ydim, zdim, I_ONE

         implicit none

         integer(kind=4), dimension(xdim:zdim), intent(in) :: ixyz    ! offset in whole grid blocks
         integer,                               intent(in) :: grid_id ! grid piece id

         integer(kind=4) :: uniq_tag
         integer(kind=4), dimension(xdim:zdim) :: r
         integer(kind=4), parameter :: N_POS=3 ! -1, 0, +1

         r = ixyz + I_ONE ! -1 => LEFT, 0 => FACE, +1 => RIGHT
         uniq_tag = int(((grid_id*N_POS+r(zdim))*N_POS+r(ydim))*N_POS+r(xdim), kind=4)

      end function uniq_tag

   end subroutine find_neighbors_SFC

!>
!! \brief Make full description of intra-level communication with neighbors. Brute-force approach.
!!
!! \details Assume that cuboids don't collide (no overlapping grid pieces on same refinement level are allowed)
!!
!! This is very general but also quite slow approach. If each of p processes has g grid containers on current level
!! (giving n = p * g grids on the level), the cost is proportional to p*p*g or n*n/p.
!!
!! Current implementation (commit a27c945a) implies correct update of all corners, even on complicated refinement
!! topologies (concave fine region - convect coarse region or fine regions touching each other only by corners).
!! Previous implementation could correctly fill the corners only on an uniform grid and when boundary exchange was
!! called for x, y and z directions separately. Warning: that change introduces measurable performance degradation!
!! This is caused by the fact that in 3D it is required to make 26 separate exchanges to fill all guardcells (in
!! cg_list_bnd::internal_boundaries), while in previous approach only 6 exchanges were required. Unfortunately
!! the previous approach did not work properly for complicated refinements.
!!
!! \todo consider going back to sweeped boundary exchanges (only 6 neighbours to communicate with) as soon as
!! find_neighbors_SFC is tested enough to be chosen as the only option for AMR 'blocky' grids.
!<

   subroutine find_neighbors_bruteforce(this)

      use cg_list,    only: cg_list_element
      use constants,  only: xdim, ydim, zdim, cor_dim, ndims, LO, HI, BND_MPI_FC, BND_FC
      use domain,     only: dom
      use gcpa,       only: gcpa_T
      use grid_cont,  only: grid_container, is_overlap
      use mpisetup,   only: FIRST, LAST, proc

      implicit none

      class(cg_list_neighbors_T), intent(inout) :: this !< object invoking type bound procedure

      type(grid_container),  pointer                  :: cg      !< grid container that we are currently working on
      type(cg_list_element), pointer                  :: cgl
      integer                                         :: j, b, id, ix, iy, iz
      integer(kind=8)                                 :: n_lbnd_face_cells
      integer(kind=4)                                 :: d, dd, hl, lh, tag
      integer(kind=8), dimension(xdim:zdim)           :: per
      integer(kind=8), dimension(xdim:zdim, LO:HI)    :: b_layer, poff, aux
      type :: fmap
         logical, dimension(:,:,:), allocatable       :: map
         integer(kind=8), dimension(xdim:zdim, LO:HI) :: b_layer
      end type fmap
      type(fmap), dimension(xdim:zdim, LO:HI)         :: f
      integer(kind=8), dimension(ndims, LO:HI)        :: box_8   !< temporary storage
      type(gcpa_T) :: l_pse

      call l_pse%init(this)

      cgl => this%first
      do while (associated(cgl))
         cg => cgl%cg

         if (allocated(cg%i_bnd)) deallocate(cg%i_bnd)
         if (allocated(cg%o_bnd)) deallocate(cg%o_bnd)
         allocate(cg%i_bnd(xdim:cor_dim), cg%o_bnd(xdim:cor_dim))

         per(:) = 0
         where (dom%periodic(:)) per(:) = this%n_d(:)

         ! Create maps to mark neighbouring face cells
         do d = xdim, zdim
            if (.not. allocated(cg%i_bnd(d)%seg)) allocate(cg%i_bnd(d)%seg(0))
            if (.not. allocated(cg%o_bnd(d)%seg)) allocate(cg%o_bnd(d)%seg(0))
            if (dom%has_dir(d)) then
               do lh = LO, HI
                  hl = LO+HI-lh ! HI for LO, LO for HI
                  f(d, lh)%b_layer(:,:) = cg%my_se(:, :)
                  f(d, lh)%b_layer(d, hl) = f(d, lh)%b_layer(d, lh)          ! interior cell layer, 1 cell thick, without corners
                  allocate(f(d, lh)%map(f(d, lh)%b_layer(xdim,LO):f(d, lh)%b_layer(xdim,HI), &
                       &                f(d, lh)%b_layer(ydim,LO):f(d, lh)%b_layer(ydim,HI), &
                       &                f(d, lh)%b_layer(zdim,LO):f(d, lh)%b_layer(zdim,HI)))
                  f(d, lh)%map = .false.
               enddo
            endif
         enddo
         allocate(cg%i_bnd(cor_dim)%seg(0), cg%o_bnd(cor_dim)%seg(0))

         do j = FIRST, LAST
            do b = lbound(this%dot%gse(j)%c(:), dim=1), ubound(this%dot%gse(j)%c(:), dim=1)
               box_8 = int(cg%lhn, kind=8)
               if (is_overlap(box_8, this%dot%gse(j)%c(b)%se(:,:), per(:))) then                ! identify processes with interesting neighbour data

                  do d = xdim, zdim
                     if (dom%has_dir(d)) then
                        do lh = LO, HI
                           hl = LO+HI-lh ! HI for LO, LO for HI

                           ! First, update the map of faces with neighbours

                           ! create 1-layer thick map of neighbours
                           b_layer = this%dot%gse(j)%c(b)%se
                           b_layer(d, hl) = b_layer(d, hl) - lh+hl  ! move the opposite boundary
                           b_layer(d, lh) = b_layer(d, hl)

                           do id = -1, 1 ! scan through periodic images of the domain
                              if (id == 0 .or. per(d)>0) then
                                 poff = b_layer
                                 poff(d, :) = poff(d, :) + id*per(d)
                                 poff(:, LO) = max(poff(:, LO), f(d, lh)%b_layer(:, LO))
                                 poff(:, HI) = min(poff(:, HI), f(d, lh)%b_layer(:, HI))
                                 ! construct the layer to be send to the _interior_ of neighbouring grid and set the flag map
                                 if (is_overlap(f(d, lh)%b_layer, poff)) &
                                      f(d, lh)%map(poff(xdim,LO):poff(xdim,HI), poff(ydim,LO):poff(ydim,HI), poff(zdim,LO):poff(zdim,HI)) = .true.
                              endif
                           enddo

                           ! Second, describe incoming data
                           b_layer = cg%my_se
                           b_layer(d, hl) = b_layer(d, lh) + (lh-hl)
                           b_layer(d, lh) = b_layer(d, lh) + (lh-hl)*dom%nb ! dom%nb thick layer without corners
                           b_layer(:d-1, LO) = b_layer(:d-1, LO) - dom%nb*dom%D_(:d-1)
                           b_layer(:d-1, HI) = b_layer(:d-1, HI) + dom%nb*dom%D_(:d-1) ! corners added in only one way
                           ! faces and corners are included in y and z direction to minimize number of pieces in non-cartesian grid decompositions

                           ! set up segments to be received
                           do iz = -1, 1 ! scan through all periodic possibilities
                              if (iz == 0 .or. per(zdim)>0) then
                                 do iy = -1, 1
                                    if (iy == 0 .or. per(ydim)>0) then
                                       do ix = -1, 1
                                          if (ix == 0 .or. per(xdim)>0) then
                                             poff = this%dot%gse(j)%c(b)%se
                                             poff(:, LO) = poff(:, LO) + [ ix, iy, iz ] * per(:)
                                             poff(:, HI) = poff(:, HI) + [ ix, iy, iz ] * per(:)
                                             if (is_overlap(b_layer, poff)) then
                                                poff(:, LO) = max(b_layer(:, LO), poff(:, LO))
                                                poff(:, HI) = min(b_layer(:, HI), poff(:, HI))
                                                aux = this%dot%gse(j)%c(b)%se
                                                aux(:, LO) = aux(:, LO) + [ ix, iy, iz ] * per(:)
                                                aux(:, HI) = aux(:, HI) + [ ix, iy, iz ] * per(:)
                                                tag = uniq_tag(cg%my_se, aux, b)
                                                aux = poff
                                                aux(d, :) = aux(d, :) + [ -1, 1 ]
                                                if (is_overlap(cg%my_se, aux)) then
                                                   dd = d
                                                else
                                                   dd = cor_dim
                                                endif
                                                call cg%i_bnd(dd)%add_seg(j, poff, tag)
                                                if (j == proc) cg%i_bnd(dd)%seg(ubound(cg%i_bnd(dd)%seg, dim=1))%local => l_pse%l_pse(b)%p
                                             endif
                                          endif
                                       enddo
                                    endif
                                 enddo
                              endif
                           enddo

                           ! Third, describe outgoing data
                           !> \warning replicated code, see above
                           b_layer = this%dot%gse(j)%c(b)%se
                           b_layer(d, hl) = b_layer(d, lh) + (lh-hl)
                           b_layer(d, lh) = b_layer(d, lh) + (lh-hl)*dom%nb ! dom%nb thick layer without corners
                           b_layer(:d-1, LO) = b_layer(:d-1, LO) - dom%nb*dom%D_(:d-1)
                           b_layer(:d-1, HI) = b_layer(:d-1, HI) + dom%nb*dom%D_(:d-1) ! corners added in only one way
                           ! faces and corners are included in y and z direction to minimize number of pieces in non-cartesian grid decompositions

                           ! set up segments to be send
                           do iz = -1, 1 ! scan through all periodic possibilities
                              if (iz == 0 .or. per(zdim)>0) then
                                 do iy = -1, 1
                                    if (iy == 0 .or. per(ydim)>0) then
                                       do ix = -1, 1
                                          if (ix == 0 .or. per(xdim)>0) then
                                             poff = b_layer
                                             poff(:, LO) = poff(:, LO) + [ ix, iy, iz ] * per(:)
                                             poff(:, HI) = poff(:, HI) + [ ix, iy, iz ] * per(:)
                                             if (is_overlap(poff(:,:), cg%my_se)) then
                                                poff(:, LO) = max(cg%my_se(:, LO), poff(:, LO))
                                                poff(:, HI) = min(cg%my_se(:, HI), poff(:, HI))
                                                aux = cg%my_se
                                                aux(:, LO) = aux(:, LO) - [ ix, iy, iz ] * per(:)
                                                aux(:, HI) = aux(:, HI) - [ ix, iy, iz ] * per(:)
                                                tag = uniq_tag(this%dot%gse(j)%c(b)%se, aux, cg%grid_id)
                                                aux = poff
                                                aux(:, LO) = aux(:, LO) - [ ix, iy, iz ] * per(:)
                                                aux(:, HI) = aux(:, HI) - [ ix, iy, iz ] * per(:)
                                                aux(d, :) = aux(d, :) + [ -1, 1 ]
                                                if (is_overlap(this%dot%gse(j)%c(b)%se, aux)) then
                                                   dd = d
                                                else
                                                   dd = cor_dim
                                                endif
                                                call cg%o_bnd(dd)%add_seg(j, poff, tag)
                                             endif
                                          endif
                                       enddo
                                    endif
                                 enddo
                              endif
                           enddo

                        enddo
                     endif
                  enddo

               endif

            enddo
         enddo

         ! Detect fine-coarse boundaries and update boundary types. When not all mapped cells are facing neighbours, then we may deal with fine/coarse boundary (full or partial)
         do d = xdim, zdim
            if (dom%has_dir(d)) then
               do lh = LO, HI
                  n_lbnd_face_cells = count(f(d,lh)%map(:,:,:))
                  if (.not. cg%ext_bnd(d, lh)) then
                     if (n_lbnd_face_cells < size(f(d,lh)%map(:,:,:))) cg%bnd(d, lh) = BND_MPI_FC
                     if (n_lbnd_face_cells == 0)                       cg%bnd(d, lh) = BND_FC
                  endif
                  deallocate(f(d,lh)%map)
               enddo
            endif
         enddo

         cgl => cgl%nxt
      enddo

      call l_pse%cleanup

   contains

      !>
      !! \brief Create unique tag for cg - cg exchange
      !!
      !! \details If we put a constraint that a grid piece can not be smaller than dom%nb, then total number of neighbours that affect local guardcells is
      !! * 3 for cartesian decomposition (including AMR with equal-size blocks)
      !! * 2 to 4 for noncartesian decomposition
      !! * more for AMR with consolidated blocks (unimplemented yet, not compatible with current approach)
      !! Thus, in each direction we can describe realtive position as one of four cases, or a bit easier one of five cases:
      !! * FAR_LEFT, FAR_RIGHT - corner neighbours, either touching corner or a bit further away
      !! * LEFT, RIGHT - partially face, partially corner neighbours
      !! * FACE - face neighbour (may cover also some corners)
      !<
      pure function uniq_tag(se, nb_se, grid_id)

         use constants, only: LO, HI, xdim, ydim, zdim, INVALID

         implicit none

         integer(kind=8), dimension(xdim:zdim, LO:HI), intent(in) :: se       ! a grid piece
         integer(kind=8), dimension(xdim:zdim, LO:HI), intent(in) :: nb_se    ! neighboring grid piece
         integer,                                      intent(in) :: grid_id  ! grid piece id

         integer(kind=4) :: uniq_tag
         integer, dimension(xdim:zdim) :: r
         integer :: d
         enum, bind(C)
            enumerator :: FAR_LEFT=0, LEFT, FACE, RIGHT, FAR_RIGHT, N_POS
         end enum

         r = INVALID
         do d = xdim, zdim
            if (nb_se(d, LO) > se(d, HI)) then
               r(d) = FAR_RIGHT
            else if (nb_se(d, HI) < se(d, LO)) then
               r(d) = FAR_LEFT
            else if ((nb_se(d, LO) < se(d, LO)) .and. (nb_se(d, HI) < se(d, HI)) .and. (nb_se(d, HI) >= se(d, LO))) then
               r(d) = LEFT
            else if ((nb_se(d, HI) > se(d, HI)) .and. (nb_se(d, LO) > se(d, LO)) .and. (nb_se(d, LO) <= se(d, HI))) then
               r(d) = RIGHT
            else
               r(d) = FACE
            endif
         enddo
         uniq_tag = int(((grid_id*N_POS+r(zdim))*N_POS+r(ydim))*N_POS+r(xdim), kind=4)

      end function uniq_tag

   end subroutine find_neighbors_bruteforce

end module cg_list_neighbors
