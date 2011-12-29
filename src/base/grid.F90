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
#include "macros.h"
!>
!! \brief (DW) Module containing lists of grid containers for computational mesh and initialization and cleanup routines
!<
module grid

   use gc_list,     only: cg_list_global, cg_list
   use cg_list_lev, only: cg_list_level, cg_list_patch

   implicit none

   private
   public :: init_grid, cleanup_grid, all_cg, base_lev, leaves, mpi_bnd_types, set_q_mbc, is_overlap

   type(cg_list_global) :: all_cg                                     !< all grid containers; \todo restore protected
   type(cg_list_level), target  :: base_lev                !< base level grid containers \todo restore "protected"
   type(cg_list), protected  :: leaves                                !< grid containers not fully covered by finer grid containers
   integer, parameter :: NBD = 1                                      !< at the moment the base domain may be composed of only one patch
   type(cg_list_patch), dimension(NBD), target, protected :: base_dom !< base level patches; \todo relax the NBD=1 restriction if we want something like L-shaped or more complex domains

   interface is_overlap
      module procedure is_overlap_simple, is_overlap_per
   end interface

contains

!>
!! \brief Routine that allocates all grid containers and most important field arrays inside gcs
!<
   subroutine init_grid

      use constants,  only: PIERNIK_INIT_DOMAIN, AT_NO_B, AT_OUT_B, AT_IGNORE, INVALID, &
           &                ndims, xdim, zdim, fluid_n, uh_n, mag_n, wa_n, u0_n, b0_n, cs_i2_n, base_level_id, base_level_offset
      use dataio_pub, only: printinfo, die, code_progress
      use domain,     only: dom, pdom, is_multicg
      use fluidindex, only: flind
      use gc_list,    only: cg_list_element
      use global,     only: repeat_step
      use grid_cont,  only: grid_container
      use mpisetup,   only: proc, inflate_req, FIRST

      implicit none

      integer :: nrq, d
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
      type(cg_list_patch), pointer :: pbd

      if (code_progress < PIERNIK_INIT_DOMAIN) call die("[grid:init_grid] domain not initialized.")

#ifdef VERBOSE
      call printinfo("[grid:init_grid]: commencing...")
#endif /* VERBOSE */

      ! Create the empty main lists for base level only.
      ! Refinement lists will be added by iterating the initproblem::init_prob routine, in restart_hdf5::read_restart_hdf5 or in not_yet_implemented::refinement_update
      ! Underground levels will be added in multigrid::init_multigrid
      call all_cg%init
      call base_lev%init
      call leaves%init
      do d = lbound(base_dom, dim=1), ubound(base_dom, dim=1) ! currently we have only one base patch
         call base_dom(d)%init
      enddo

      pbd => base_dom(NBD)
      call dom2cg(pdom%n_d(:), base_level_offset, base_level_id, pbd)

!#ifdef VERBOSE
      call base_lev%print_segments
!#endif /* VERBOSE */

#ifdef VERBOSE
      call printinfo("[grid:init_grid]: all_cg finished. \o/")
#endif /* VERBOSE */

      call all_cg%reg_var(wa_n, AT_IGNORE)                 !! BEWARE: magic string across multiple files
      call all_cg%reg_var(fluid_n, AT_NO_B, flind%all)     !! Main array of all fluids' components, "u"
      call all_cg%reg_var(uh_n, AT_IGNORE, flind%all)      !! Main array of all fluids' components (for t += dt/2)
      call all_cg%reg_var(mag_n, AT_OUT_B, ndims)          !! Main array of magnetic field's components, "b"
      if (repeat_step) then
         call all_cg%reg_var(u0_n, AT_IGNORE, flind%all)   !! Copy of main array of all fluids' components
         call all_cg%reg_var(b0_n, AT_IGNORE, ndims)       !! Copy of main array of magnetic field's components
      endif

      nrq = 0
      cgl => all_cg%first
      do while (associated(cgl))
         cg => cgl%cg

         cg%u  => cg%w(all_cg%ind_4d(fluid_n))%arr
         cg%b  => cg%w(all_cg%ind_4d(mag_n))%arr
         cg%wa => cg%q(all_cg%ind(wa_n))%arr

         if (allocated(cg%w)) then
            do d = xdim, zdim
               if (allocated(cg%w(1)%w_i_mbc(d, pdom%nb)%mbc)) nrq = nrq + 2 * count(cg%w(1)%w_i_mbc(d, pdom%nb)%mbc(:) /= INVALID) ! w(1) is probably fluid, but it can be any registered field
            enddo
         endif

         cgl => cgl%nxt
      enddo
      call inflate_req(nrq)

#ifdef ISO
      if (is_multicg) call die("[grid:init_cs_iso2] multiple grid pieces per procesor not fully implemented yet") !nontrivial maxval

      call all_cg%reg_var(cs_i2_n, AT_NO_B) ! BEWARE: magic string across multiple files

      cgl => all_cg%first
      do while (associated(cgl))
         cgl%cg%cs_iso2 => cgl%cg%q(all_cg%ind(cs_i2_n))%arr
         cgl%cg%cs_iso2(:,:,:) = maxval(flind%all_fluids(:)%cs2)   ! set cs2 with sane values
         cgl => cgl%nxt
      enddo
#endif /* ISO */

#ifdef VERBOSE
      call printinfo("[grid:init_grid]: cg finished. \o/")
#endif /* VERBOSE */

   end subroutine init_grid

!> \brief Create new cg according to domain decomposition data and add them to appropriate lists
   subroutine dom2cg(n_d, offset, level, patch)

      use constants,     only: ndims, I_ONE, base_level_id
      use decomposition, only: divide_domain
      use dataio_pub,    only: warn, die
      use domain,        only: is_mpi_noncart, is_multicg, is_refined, is_uneven
      use mpi,           only: MPI_IN_PLACE, MPI_COMM_NULL, MPI_LOGICAL, MPI_LAND, MPI_LOR
      use mpisetup,      only: proc, comm, ierr, LAST, master
      use types,         only: cdd

      implicit none

      integer(kind=4), dimension(ndims), intent(in) :: n_d
      integer(kind=8), dimension(ndims), intent(in) :: offset
      integer, intent(in) :: level
      type(cg_list_patch), pointer, intent(inout) :: patch

      integer :: g
      logical :: dom_divided

      patch%n_d = n_d
      patch%off = offset
      patch%parent => null()
      patch%children => null()

      if (level == base_level_id) then
         base_lev%lev = level
         base_lev%n_d = n_d
         patch%list_level => base_lev
      else
         patch%list_level => null()
         call die("[grid:dom2cg] Only base level is currently supported")
      endif


      dom_divided = divide_domain(patch)
      call MPI_Allreduce(MPI_IN_PLACE, dom_divided, I_ONE, MPI_LOGICAL, MPI_LAND, comm, ierr)
      if (.not. dom_divided .and. master) call die("[grid:dom2cg] Domain decomposition failed")

      !\todo Analyze the decomposition and set up [ is_uneven, is_mpi_noncart, is_refined, ... ]
      is_multicg = (ubound(base_lev%pse(proc)%sel(:, :, :), dim=1) > 1)
      call MPI_Allreduce(MPI_IN_PLACE, is_multicg, I_ONE, MPI_LOGICAL, MPI_LOR, comm, ierr)
      if (is_multicg .and. cdd%comm3d /= MPI_COMM_NULL) call die("[grid:dom2cg] is_multicg cannot be used with comm3d")
      if (is_refined) then
         is_mpi_noncart = .true.
         is_multicg = .true.
      endif
      if (is_mpi_noncart) is_uneven = .true.

      ! bnd_[xyz][lr] now become altered according to local topology of processes
      if (is_multicg .and. master) call warn("[grid:dom2cg] Multiple blocks per process not fully implemented yet")
      if (is_refined) call die("[grid:dom2cg] Refinements are not implemented")

      do g = lbound(base_lev%pse(proc)%sel(:,:,:), dim=1), ubound(base_lev%pse(proc)%sel(:,:,:), dim=1)

         ! create the new element in the patch and initialize it
         call patch%add
         call patch%last%cg%init(base_lev%n_d, base_lev%pse(proc)%sel(g, :, :), g, level)

         call mpi_bnd_types(patch%last%cg)
         call set_q_mbc(patch%last%cg)

         ! add to the other lists
         call all_cg%add(patch%last%cg)
         if (level == base_level_id) call base_lev%add(patch%last%cg)
         call leaves%add(patch%last%cg)

      enddo

   end subroutine dom2cg

! \brief Update the list of leaves
! subroutine update leaves
! end subroutine update leaves

!> \brief deallocate everything
   subroutine cleanup_grid

      use gc_list, only: cg_list_element

      implicit none

      integer :: d
      type(cg_list_element), pointer :: cgle

      call leaves%delete
      call base_lev%delete
      if (allocated(base_lev%pse)) deallocate(base_lev%pse)
      do d = lbound(base_dom, dim=1), ubound(base_dom, dim=1) ! currently we have only one base patch
         call base_dom(d)%delete
      enddo

      ! manually deallocate all grid containers first
      cgle => all_cg%first
      do while (associated(cgle))
         deallocate(cgle%cg)
         cgle => cgle%nxt
      enddo
      if (allocated(all_cg%q_lst)) deallocate(all_cg%q_lst)
      if (allocated(all_cg%w_lst)) deallocate(all_cg%w_lst)
      call all_cg%delete

   end subroutine cleanup_grid

!> \brief Initialize the communicators for q even if there are no q arrays at the moment

   subroutine set_q_mbc(this)

      use constants, only: ndims, xdim, zdim
      use domain,    only: dom
      use grid_cont, only: grid_container, set_mpi_types

      implicit none

      class(grid_container), intent(inout) :: this
      integer :: d, ib, g

      allocate(this%q_i_mbc(ndims, dom%nb), this%q_o_mbc(ndims, dom%nb))
      do d = xdim, zdim
         do ib = 1, dom%nb
            if (allocated(this%i_bnd(d, ib)%seg)) then
               allocate(this%q_i_mbc(d, ib)%mbc(lbound(this%i_bnd(d, ib)%seg, dim=1):ubound(this%i_bnd(d, ib)%seg, dim=1)), &
                    &   this%q_o_mbc(d, ib)%mbc(lbound(this%i_bnd(d, ib)%seg, dim=1):ubound(this%i_bnd(d, ib)%seg, dim=1)))
               do g = lbound(this%i_bnd(d, ib)%seg, dim=1), ubound(this%i_bnd(d, ib)%seg, dim=1)
                  call set_mpi_types(this%n_(:), this%i_bnd(d, ib)%seg(g)%se(:,:), this%q_i_mbc(d, ib)%mbc(g))
                  call set_mpi_types(this%n_(:), this%o_bnd(d, ib)%seg(g)%se(:,:), this%q_o_mbc(d, ib)%mbc(g))
               enddo
            endif
         enddo
      enddo

   end subroutine set_q_mbc

!> \brief Create MPI types for boundary exchanges

   subroutine mpi_bnd_types(this)

      use constants,  only: FLUID, ARR, xdim, zdim, ndims, LO, HI, BND, BLK, I_ONE
      use dataio_pub, only: die
      use domain,     only: dom
      use fluidindex, only: flind
      use grid_cont,  only: grid_container
      use mpi,        only: MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, MPI_COMM_NULL
      use mpisetup,   only: ierr, proc, FIRST, LAST, procmask
      use types,      only: cdd

      implicit none

      class(grid_container), intent(inout) :: this

      integer(kind=4), dimension(:), allocatable :: sizes, subsizes, starts
      integer :: t, g, j, b
      integer(kind=4) :: d, dd, hl, lh, ib
      integer(kind=4), parameter, dimension(FLUID:ARR) :: dims = [ I_ONE+ndims, I_ONE+ndims, I_ONE+ndims, ndims ] !< dimensionality of arrays
      integer(kind=4), dimension(FLUID:ARR) :: nc
      integer(kind=8), dimension(xdim:zdim) :: ijks, per
      integer(kind=8), dimension(xdim:zdim, LO:HI) :: b_layer, bp_layer, poff
      type(cg_list_level), pointer :: my_lev

      nc = [ flind%all, ndims, max(flind%crs%all,I_ONE), I_ONE ]      !< number of fluids, magnetic field components, CRs, and 1 for a rank-3 array

      my_lev => base_lev
      do while (this%level_id /= my_lev%lev)
         if (this%level_id > my_lev%lev) then
            my_lev => my_lev%finer
         else
            my_lev => my_lev%coarser
         endif
         if (.not. associated(my_lev)) call die("[grid:mpi_bnd_types] Cannot find my refinement level")
      enddo
      if (this%level_id /= my_lev%lev) call die("[grid:mpi_bnd_types] Found wrong level")

      if (allocated(this%i_bnd) .or. allocated(this%o_bnd)) call die("[grid:mpi_bnd_types] this%i_bnd or this%o_bnd already allocated")
      allocate(this%i_bnd(xdim:zdim, dom%nb), this%o_bnd(xdim:zdim, dom%nb))

      ! assume that cuboids fill the domain and don't collide

      ijks(:) = this%ijkse(:, LO) - this%off(:)
      per(:) = 0
      where (dom%periodic(:)) per(:) = my_lev%n_d(:)

      if (cdd%comm3d == MPI_COMM_NULL) then

         do d = xdim, zdim
            if (dom%has_dir(d) .and. .not. this%empty) then

               ! identify processes with interesting neighbour data
               procmask(:) = 0
               do lh = LO, HI
                  hl = LO+HI-lh ! HI for LO, LO for HI
                  b_layer(:,:) = this%my_se(:, :)
                  b_layer(d, lh) = b_layer(d, lh) + lh-hl ! -1 for LO, +1 for HI
                  b_layer(d, hl) = b_layer(d, lh) ! adjacent boundary layer, 1 cell wide, without corners
                  do j = FIRST, LAST
                     do b = lbound(my_lev%pse(j)%sel(:, :, :), dim=1), ubound(my_lev%pse(j)%sel(:, :, :), dim=1)
                        if (is_overlap(b_layer(:,:), my_lev%pse(j)%sel(b, :, :), per(:))) procmask(j) = procmask(j) + 1 ! count how many boundaries we share with that process
                     enddo
                  enddo
               enddo
               do ib = 1, dom%nb
                  allocate(this%i_bnd(d, ib)%seg(sum(procmask(:))), this%o_bnd(d, ib)%seg(sum(procmask(:))))
               enddo

               ! set up segments to be sent or received
               g = 0
               do j = FIRST, LAST
                  if (procmask(j) /= 0) then
                     do lh = LO, HI
                        hl = LO+HI-lh
                        do b = lbound(my_lev%pse(j)%sel(:, :, :), dim=1), ubound(my_lev%pse(j)%sel(:, :, :), dim=1)
                           b_layer(:,:) = this%my_se(:, :)
                           b_layer(d, lh) = b_layer(d, lh) + lh-hl
                           b_layer(d, hl) = b_layer(d, lh) ! adjacent boundary layer, 1 cell wide, without corners
                           bp_layer(:, :) = b_layer(:, :)
                           where (per(:) > 0)
                              bp_layer(:, LO) = mod(b_layer(:, LO) + per(:), per(:))
                              bp_layer(:, HI) = mod(b_layer(:, HI) + per(:), per(:)) ! adjacent boundary layer, 1 cell wide, without corners, corrected for periodicity
                           endwhere
                           !> \todo save b_layer(:,:) and bp_layer(:,:) and move above calculations outside the b loop

                           if (is_overlap(bp_layer(:,:), my_lev%pse(j)%sel(b, :, :))) then
                              poff(:,:) = bp_layer(:,:) - b_layer(:,:) ! displacement due to periodicity
                              bp_layer(:, LO) = max(bp_layer(:, LO), my_lev%pse(j)%sel(b, :, LO))
                              bp_layer(:, HI) = min(bp_layer(:, HI), my_lev%pse(j)%sel(b, :, HI))
                              b_layer(:,:) = bp_layer(:,:) - poff(:,:)
                              g = g + 1
                              do ib = 1, dom%nb
                                 this%i_bnd(d, ib)%seg(g)%proc = j
                                 this%i_bnd(d, ib)%seg(g)%se(:,LO) = b_layer(:, LO) + ijks(:)
                                 this%i_bnd(d, ib)%seg(g)%se(:,HI) = b_layer(:, HI) + ijks(:)
                                 if (any(this%i_bnd(d, ib)%seg(g)%se(d, :) < 0)) &
                                      this%i_bnd(d, ib)%seg(g)%se(d, :) = this%i_bnd(d, ib)%seg(g)%se(d, :) + my_lev%n_d(d)
                                 if (any(this%i_bnd(d, ib)%seg(g)%se(d, :) > this%n_b(d) + 2*dom%nb)) &
                                      this%i_bnd(d, ib)%seg(g)%se(d, :) = this%i_bnd(d, ib)%seg(g)%se(d, :) - my_lev%n_d(d)

                                 ! expand to cover corners (requires separate MPI_Waitall for each direction)
                                 !! \todo create separate %mbc for corner-less exchange with one MPI_Waitall (can scale better)
                                 !! \warning edges and corners will be filled multiple times
                                 do dd = xdim, zdim
                                    if (dd /= d .and. dom%has_dir(dd)) then
                                       this%i_bnd(d, ib)%seg(g)%se(dd, LO) = this%i_bnd(d, ib)%seg(g)%se(dd, LO) - ib
                                       this%i_bnd(d, ib)%seg(g)%se(dd, HI) = this%i_bnd(d, ib)%seg(g)%se(dd, HI) + ib
                                    endif
                                 enddo
                                 this%o_bnd(d, ib)%seg(g) = this%i_bnd(d, ib)%seg(g)
                                 this%i_bnd(d, ib)%seg(g)%tag = int(HI*ndims*b           + (HI*d+lh-LO), kind=4) ! Assume that we won't mix communication with different ib
                                 this%o_bnd(d, ib)%seg(g)%tag = int(HI*ndims*this%grid_id + (HI*d+hl-LO), kind=4)
                                 select case (lh)
                                    case (LO)
                                       this%i_bnd(d, ib)%seg(g)%se(d, LO) = this%i_bnd(d, ib)%seg(g)%se(d, HI) - (ib - 1)
                                       this%o_bnd(d, ib)%seg(g)%se(d, LO) = this%i_bnd(d, ib)%seg(g)%se(d, HI) + 1
                                       this%o_bnd(d, ib)%seg(g)%se(d, HI) = this%o_bnd(d, ib)%seg(g)%se(d, LO) + (ib - 1)
                                    case (HI)
                                       this%i_bnd(d, ib)%seg(g)%se(d, HI) = this%i_bnd(d, ib)%seg(g)%se(d, LO) + (ib - 1)
                                       this%o_bnd(d, ib)%seg(g)%se(d, HI) = this%i_bnd(d, ib)%seg(g)%se(d, LO) - 1
                                       this%o_bnd(d, ib)%seg(g)%se(d, LO) = this%o_bnd(d, ib)%seg(g)%se(d, HI) - (ib - 1)
                                 end select

                              enddo
                           endif
                        enddo
                     enddo
                  endif
               enddo
            endif
         enddo

      else

         do d = xdim, zdim
            if (dom%has_dir(d)) then
               do t = FLUID, ARR  ! fluid, Bfield, wcr, grav

                  allocate(sizes(dims(t)), subsizes(dims(t)), starts(dims(t)))
                  if (dims(t) == 1+ndims) sizes(1) = nc(t)
                  sizes(dims(t)-zdim+xdim:dims(t)) = this%n_(:)

                  do ib = 1, dom%nb

                     subsizes(:) = sizes(:)
                     subsizes(dims(t)-zdim+d) = ib
                     starts(:) = 0

                     starts(dims(t)-zdim+d) = dom%nb-ib
                     call MPI_Type_create_subarray(dims(t), sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, this%mbc(t, d, LO, BND, ib), ierr)
                     call MPI_Type_commit(this%mbc(t, d, LO, BND, ib), ierr)

                     starts(dims(t)-zdim+d) = dom%nb
                     call MPI_Type_create_subarray(dims(t), sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, this%mbc(t, d, LO, BLK, ib), ierr)
                     call MPI_Type_commit(this%mbc(t, d, LO, BLK, ib), ierr)

                     starts(dims(t)-zdim+d) = this%n_b(d) + dom%nb - ib
                     call MPI_Type_create_subarray(dims(t), sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, this%mbc(t, d, HI, BLK, ib), ierr)
                     call MPI_Type_commit(this%mbc(t, d, HI, BLK, ib), ierr)

                     starts(dims(t)-zdim+d) = this%ijkse(d, HI)
                     call MPI_Type_create_subarray(dims(t), sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, this%mbc(t, d, HI, BND, ib), ierr)
                     call MPI_Type_commit(this%mbc(t, d, HI, BND, ib), ierr)

                  enddo

                  deallocate(sizes, subsizes, starts)

               enddo
            endif
         enddo

      endif

   end subroutine mpi_bnd_types

!>
!! \brief is_overlap_per checks if two given blocks placed within a periodic domain are overlapping.
!!
!! \details to handle shearing box which is divided in y-direction at the edges, one has to provide another subroutine (is_overlap_per_shear) and add it to interface is_overlap
!<
   logical function is_overlap_per(this, other, periods) result(share)

      use constants,  only: xdim, ydim, zdim, ndims, LO, HI
      use domain,     only: dom

      implicit none

      integer(kind=8), dimension(xdim:zdim, LO:HI), intent(in) :: this !< object invoking type-bound procedure
      integer(kind=8), dimension(xdim:zdim, LO:HI), intent(in) :: other !< the other box
      integer(kind=8), dimension(xdim:zdim), intent(in)        :: periods !< where >0 then the direction is periodic with the given number of cells

      integer :: i, j, k
      integer(kind=8), dimension(xdim:zdim, LO:HI) :: oth

      share = .false.
      do i = -1, 1
         if ((dom%has_dir(xdim) .or. periods(xdim)>0) .or. i==0) then
            do j = -1, 1
               if ((dom%has_dir(ydim) .or. periods(ydim)>0) .or. j==0) then
                  do k = -1, 1
                     if ((dom%has_dir(zdim) .or. periods(zdim)>0) .or. k==0) then
                        oth(:,:) = other(:,:) + reshape([i*periods(xdim), j*periods(ydim), k*periods(zdim), i*periods(xdim), j*periods(ydim), k*periods(zdim)], [ndims, HI])
                        share = share .or. is_overlap_simple(this, oth)
                     endif
                  enddo
               endif
            enddo
         endif
      enddo

   end function is_overlap_per

!!
!> \brief is_overlap_simple checks if two given blocks placed within a nonperiodic domain are overlapping.
!! This routine is not supposed to take care of periodic domain - use is_overlap_per when you check overlap for boxes that cross the periodic domain boundary
!<
   logical function is_overlap_simple(this, other) result(share)

      use constants,  only: xdim, zdim, LO, HI
      use domain,     only: dom

      implicit none

      integer(kind=8), dimension(xdim:zdim, LO:HI), intent(in) :: this !< object invoking type-bound procedure
      integer(kind=8), dimension(xdim:zdim, LO:HI), intent(in) :: other !< the other box

      integer :: d

      share = .true.
      do d = xdim, zdim
         if (dom%has_dir(d)) share = share .and. (other(d, LO) <= this(d, HI)) .and. (other(d, HI) >= this(d, LO))
      enddo

   end function is_overlap_simple

end module grid
