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
!! \brief (DW) Module containing routines to specify required computational mesh.
!! \date January/February 2006
!!
!!
!! In this module a following namelist of parameters is specified:
!! \copydetails grid::init_grid
!<
module grid

   use constants, only: ndims
   use grid_cont, only: cg_set

   implicit none

   private
   public :: init_grid, init_arrays, grid_mpi_boundaries_prep, arr3d_boundaries, cleanup_grid
   public :: total_ncells, cga, D_x, D_y, D_z, D_

   integer(kind=8), protected :: total_ncells !< total number of %grid cells
   integer, dimension(ndims), protected :: D_!< set to 1 for existing directions, 0 otherwise. Useful for dimensionally-safe indices for difference operators on arrays,
   integer, protected :: D_x          !< set to 1 when x-direction exists, 0 otherwise
   integer, protected :: D_y          !< set to 1 when y-direction exists, 0 otherwise.
   integer, protected :: D_z          !< set to 1 when z-direction exists, 0 otherwise.
   type(cg_set), target :: cga        !< A container for all grids.

contains

!>
!! \brief Routine which sets numbers of cells for the domain, MPI blocks and initializes direction meshes (x,y,z).
!! Also compute domain maximum and minimum of coordinates, lengths of cells and coordinates of zone centers and left/right zone boundaries.
!!
!<
   subroutine init_grid

      use constants,  only: PIERNIK_INIT_DOMAIN, xdim, ydim, zdim
      use dataio_pub, only: printinfo, die, code_progress
      use domain,     only: dom, has_dir
      use grid_cont,  only: cg_list_element
      use mpisetup,   only: proc

      implicit none

      type(cg_list_element), pointer :: cgl
      integer :: g

      if (code_progress < PIERNIK_INIT_DOMAIN) call die("[grid:init_grid] MPI not initialized.")

#ifdef VERBOSE
      call printinfo("[grid:init_grid]: commencing...")
#endif /* VERBOSE */

      if (ubound(dom%pse(proc)%sel(:,:,:), dim=1) > 1) call die("[grid:init_grid] Multiple blocks per process not implemented yet")

      allocate(cga%cg_all(1)) !At the moment we use only one grid container per thread, but this will change in AMR
      allocate(cga%cg_base%cg_l(size(cga%cg_all)))  ! We have no refinement yet, so everything is on the base level
      allocate(cga%cg_leafs%cg_l(size(cga%cg_all))) ! We have no refinement yet, so everything is the finest grid as well
      allocate(cga%cg_levels(1))                   ! We have no refinement yet, so there is only one level
      do g = lbound(cga%cg_levels(:), dim=1), ubound(cga%cg_levels(:), dim=1)
         allocate(cga%cg_levels(g)%cg_l(size(cga%cg_all)))
      enddo

      !for an uniform grid set up trivial lists of grid containers
      do g = lbound(cga%cg_all(:), dim=1), ubound(cga%cg_all(:), dim=1)
         cga%cg_base%cg_l(g)%cg      => cga%cg_all(g)
         cga%cg_leafs%cg_l(g)%cg     => cga%cg_all(g)
         cga%cg_levels(1)%cg_l(g)%cg => cga%cg_all(g)

         if (g /= lbound(cga%cg_all(:), dim=1)) then
            cga%cg_base%cg_l(g)%prv      => cga%cg_base%cg_l(g-1)
            cga%cg_leafs%cg_l(g)%prv     => cga%cg_leafs%cg_l(g-1)
            cga%cg_levels(1)%cg_l(g)%prv => cga%cg_levels(1)%cg_l(g-1)
         else
            cga%cg_base%cg_l(g)%prv      => null()
            cga%cg_leafs%cg_l(g)%prv     => null()
            cga%cg_levels(1)%cg_l(g)%prv => null()
         endif

         if (g /= ubound(cga%cg_all(:), dim=1)) then
            cga%cg_base%cg_l(g)%nxt      => cga%cg_base%cg_l(g+1)
            cga%cg_leafs%cg_l(g)%nxt     => cga%cg_leafs%cg_l(g+1)
            cga%cg_levels(1)%cg_l(g)%nxt => cga%cg_levels(1)%cg_l(g+1)
         else
            cga%cg_base%cg_l(g)%nxt      => null()
            cga%cg_leafs%cg_l(g)%nxt     => null()
            cga%cg_levels(1)%cg_l(g)%nxt => null()
         endif
      enddo

      cgl => cga%cg_base%cg_l(1)
      do while (associated(cgl))
         call cgl%cg%init(dom)
         cgl => cgl%nxt
      enddo

      where (has_dir(:))
         D_(:) = 1
      elsewhere
         D_(:) = 0
      endwhere

      ! shortcuts
      D_x = D_(xdim)
      D_y = D_(ydim)
      D_z = D_(zdim)

      total_ncells = product(int(dom%n_d(:), kind=8))
      if (any(total_ncells < dom%n_d(:))) call die("[grid:init_grid] Integer overflow: too many cells")

#ifdef VERBOSE
      call printinfo("[grid:init_grid]: finished. \o/")
#endif /* VERBOSE */

   end subroutine init_grid

!>
!! Routine that allocates all arrays
!<

   subroutine init_arrays(flind)

      use constants,   only: PIERNIK_INIT_BASE, ndims, zdim
      use diagnostics, only: my_allocate
      use dataio_pub,  only: die, code_progress
      use global,      only: repeat_step
      use fluidtypes,  only: var_numbers
      use grid_cont,   only: cg_list_element, grid_container

      implicit none

      type(var_numbers), intent(in) :: flind !< fluid database; cannot use fluidindex::flind here due to circular dependencies in some setups
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      if (code_progress < PIERNIK_INIT_BASE) call die("[arrays:init_arrays] grid or fluids not initialized.")

      call cga%get_root(cgl)
      do while (associated(cgl))
         cg => cgl%cg
         call cg%u%init( [ flind%all, cg%n_(:) ] )
         call cg%uh%init( [ flind%all, cg%n_(:) ] )
         if (repeat_step) call cg%u0%init( [ flind%all, cg%n_(:) ] )

         call cg%b%init( [ ndims, cg%n_(:) ] )
         if (repeat_step) call cg%b0%init( [ ndims, cg%n_(:) ] )

         call cg%wa%init(cg%n_(:))

#ifdef GRAV
         call my_allocate(cg%dprof, [cg%n_(zdim)], "dprof")
#endif /* GRAV */
         cgl => cgl%nxt
      enddo

   end subroutine init_arrays

!>
!! \brief deallocate everything
!<

   subroutine cleanup_grid

      use grid_cont,  only: cg_list_element

      implicit none

      type(cg_list_element), pointer :: cgl
      integer :: g

      cgl => cga%cg_base%cg_l(1)
      do while (associated(cgl))

         call cgl%cg%u%clean()
         call cgl%cg%u0%clean()
         call cgl%cg%uh%clean()

         call cgl%cg%b%clean()
         call cgl%cg%b0%clean()

         call cgl%cg%wa%clean()

#ifdef GRAV
         if (allocated(cgl%cg%dprof)) deallocate(cgl%cg%dprof)
#endif /* GRAV */

         call cgl%cg%cleanup
         cgl => cgl%nxt
      enddo

      do g = lbound(cga%cg_levels(:), dim=1), ubound(cga%cg_levels(:), dim=1)
         deallocate(cga%cg_levels(g)%cg_l)
      enddo
      deallocate(cga%cg_levels)
      deallocate(cga%cg_leafs%cg_l)
      deallocate(cga%cg_base%cg_l)
      deallocate(cga%cg_all)

   end subroutine cleanup_grid

!>
!! \brief Set up subsets of u,b and sgp arrays for MPI communication
!!
!! \todo this should be called from cg%init only
!<

   subroutine grid_mpi_boundaries_prep(numfluids, numcrs)

      use constants,  only: PIERNIK_INIT_BASE, FLUID, ARR, xdim, zdim, ndims, LO, HI, BND, BLK, INVALID, I_ONE
      use dataio_pub, only: die, code_progress
      use domain,     only: has_dir, dom, is_overlap, cdd
      use grid_cont,  only: cg_list_element, grid_container
      use mpi,        only: MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, MPI_COMM_NULL
      use mpisetup,   only: ierr, proc, FIRST, LAST, procmask

      implicit none

      integer(kind=4), intent(in) :: numfluids !< expect flind%all,     here, cannot grab it directly because of cyclic deps in CR-based setups
      integer(kind=4), intent(in) :: numcrs    !< expect flind%crs%all, here, cannot grab it directly because of cyclic deps in CR-based setups

      integer(kind=4), dimension(:), allocatable :: sizes, subsizes, starts
      integer :: d, t, g, j
      integer(kind=4) :: hl, lh
      integer(kind=4), dimension(FLUID:ARR) :: nc
      integer(kind=4), parameter, dimension(FLUID:ARR) :: dims = [ I_ONE+ndims, I_ONE+ndims, I_ONE+ndims, ndims ] !< dimensionality of arrays
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg
      integer(kind=8), dimension(xdim:zdim) :: ijks, per
      integer(kind=8), dimension(xdim:zdim, LO:HI) :: b_layer, bp_layer, poff
      logical :: sharing

      if (code_progress < PIERNIK_INIT_BASE) call die("[grid:grid_mpi_boundaries_prep] grid or fluids not initialized.")
      if (ubound(dom%pse(proc)%sel(:,:,:), dim=1) > 1) call die("[grid:grid_mpi_boundaries_prep] Multiple blocks per process not implemented yet")

      nc = [ numfluids, ndims, max(numcrs,I_ONE), I_ONE ]      !< number of fluids, magnetic field components, CRs, and 1 for rank-3 array

      call cga%get_root(cgl)
      do while (associated(cgl))
         cg => cgl%cg
         ! find neighbours and set up the MPI containers
         if (cdd%comm3d == MPI_COMM_NULL) then

            if (allocated(cg%i_bnd) .or. allocated(cg%o_bnd)) call die("[grid:grid_mpi_boundaries_prep] cg%i_bnd or cg%o_bnd already allocated")
            allocate(cg%i_bnd(xdim:zdim, FLUID:ARR), cg%o_bnd(xdim:zdim, FLUID:ARR))

            ! assume that cuboids fill the domain and don't collide

            ijks(:) = cg%ijkse(:, LO) - cg%off(:)
            per(:) = 0
            where (dom%periodic(:)) per(:) = dom%n_d(:)

            do d = xdim, zdim
               if (has_dir(d) .and. .not. cg%empty) then

                  ! identify processes with interesting neighbour data
                  procmask(:) = 0
                  do lh = LO, HI
                     hl = LO+HI-lh ! HI for LO, LO for HI
                     b_layer(:,:) = dom%pse(proc)%sel(1, :, :)
                     b_layer(d, lh) = b_layer(d, lh) + lh-hl ! -1 for LO, +1 for HI
                     b_layer(d, hl) = b_layer(d, lh) ! boundary layer without corners
                     do j = FIRST, LAST
                        call is_overlap(b_layer(:,:), dom%pse(j)%sel(1, :, :), sharing, per(:))
                        if (sharing) procmask(j) = procmask(j) + 1
                     enddo
                  enddo
                  do j = FLUID, ARR
                     allocate(cg%i_bnd(d, j)%seg(sum(procmask(:))))
                     allocate(cg%o_bnd(d, j)%seg(sum(procmask(:))))
                  enddo

                  ! set up segments to be sent or received
                  g = 0
                  do j = FIRST, LAST
                     if (procmask(j) /= 0) then
                        do lh = LO, HI
                           hl = LO+HI-lh
                           b_layer(:,:) = dom%pse(proc)%sel(1, :, :)
                           b_layer(d, lh) = b_layer(d, lh) + lh-hl
                           b_layer(d, hl) = b_layer(d, lh)

                           bp_layer(:, :) = b_layer(:, :)
                           where (per(:) > 0)
                              bp_layer(:, LO) = mod(b_layer(:, LO) + per(:), per(:))
                              bp_layer(:, HI) = mod(b_layer(:, HI) + per(:), per(:))
                           endwhere
                           call is_overlap(bp_layer(:,:), dom%pse(j)%sel(1, :, :), sharing)

                           if (sharing) then
                              poff(:,:) = bp_layer(:,:) - b_layer(:,:) ! displacement due to periodicity
                              bp_layer(:, LO) = max(bp_layer(:, LO), dom%pse(j)%sel(1, :, LO))
                              bp_layer(:, HI) = min(bp_layer(:, HI), dom%pse(j)%sel(1, :, HI))

                              b_layer(:,:) = bp_layer(:,:) - poff(:,:)
                              g = g + 1
                              do t = FLUID, ARR
                                 cg%i_bnd(d, t)%seg(g)%mbc = INVALID
                                 cg%i_bnd(d, t)%seg(g)%proc = j
                                 cg%i_bnd(d, t)%seg(g)%se(:,LO) = b_layer(:, LO) + ijks(:)
                                 cg%i_bnd(d, t)%seg(g)%se(:,HI) = b_layer(:, HI) + ijks(:)
                                 if (any(cg%i_bnd(d, t)%seg(g)%se(d, :) < 0)) &
                                      cg%i_bnd(d, t)%seg(g)%se(d, :) = cg%i_bnd(d, t)%seg(g)%se(d, :) + dom%n_d(d)
                                 if (any(cg%i_bnd(d, t)%seg(g)%se(d, :) > cg%n_b(d) + 2*cg%nb)) &
                                      cg%i_bnd(d, t)%seg(g)%se(d, :) = cg%i_bnd(d, t)%seg(g)%se(d, :) - dom%n_d(d)

                                 ! expand to cover corners (requires separate MPI_Waitall for each direction)
                                 ! \todo create separate %mbc for corner-less exchange with one MPI_Waitall (can scale better)
                                 where (has_dir(:d-1))
                                    cg%i_bnd(d, t)%seg(g)%se(:d-1, LO) = cg%i_bnd(d, t)%seg(g)%se(:d-1, LO) - cg%nb
                                    cg%i_bnd(d, t)%seg(g)%se(:d-1, HI) = cg%i_bnd(d, t)%seg(g)%se(:d-1, HI) + cg%nb
                                 endwhere
                                 cg%o_bnd(d, t)%seg(g) = cg%i_bnd(d, t)%seg(g)
                                 cg%i_bnd(d, t)%seg(g)%lh = lh
                                 cg%o_bnd(d, t)%seg(g)%lh = hl
                                 select case (lh)
                                    case (LO)
                                       cg%i_bnd(d, t)%seg(g)%se(d, LO) = cg%i_bnd(d, t)%seg(g)%se(d, HI) - (cg%nb - 1)
                                       cg%o_bnd(d, t)%seg(g)%se(d, LO) = cg%i_bnd(d, t)%seg(g)%se(d, HI) + 1
                                       cg%o_bnd(d, t)%seg(g)%se(d, HI) = cg%o_bnd(d, t)%seg(g)%se(d, LO) + (cg%nb - 1)
                                    case (HI)
                                       cg%i_bnd(d, t)%seg(g)%se(d, HI) = cg%i_bnd(d, t)%seg(g)%se(d, LO) + (cg%nb - 1)
                                       cg%o_bnd(d, t)%seg(g)%se(d, HI) = cg%i_bnd(d, t)%seg(g)%se(d, LO) - 1
                                       cg%o_bnd(d, t)%seg(g)%se(d, LO) = cg%o_bnd(d, t)%seg(g)%se(d, HI) - (cg%nb - 1)
                                 end select
                                 ! set MPI type only for non-local transfers

                                 allocate(sizes(dims(t)), subsizes(dims(t)), starts(dims(t)))

                                 starts(:) = 0
                                 if (dims(t) == 1+ndims) then
                                    sizes(1) = nc(t)
                                    subsizes(1) = sizes(1)
                                 endif
                                 sizes   (dims(t)-zdim+xdim:dims(t)) = cg%n_(:)
                                 subsizes(dims(t)-zdim+xdim:dims(t)) = int(cg%i_bnd(d, t)%seg(g)%se(:, HI) - cg%i_bnd(d, t)%seg(g)%se(:, LO) + 1, kind=4)
                                 starts  (dims(t)-zdim+xdim:dims(t)) = int(cg%i_bnd(d, t)%seg(g)%se(:, LO) - 1, kind=4)
                                 call MPI_Type_create_subarray(dims(t), sizes, subsizes, starts,  MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, cg%i_bnd(d, t)%seg(g)%mbc, ierr)
                                 call MPI_Type_commit(cg%i_bnd(d, t)%seg(g)%mbc, ierr)

                                 subsizes(dims(t)-zdim+xdim:dims(t)) = int(cg%o_bnd(d, t)%seg(g)%se(:, HI) - cg%o_bnd(d, t)%seg(g)%se(:, LO) + 1, kind=4)
                                 starts  (dims(t)-zdim+xdim:dims(t)) = int(cg%o_bnd(d, t)%seg(g)%se(:, LO) - 1, kind=4)
                                 call MPI_Type_create_subarray(dims(t), sizes, subsizes, starts,  MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, cg%o_bnd(d, t)%seg(g)%mbc, ierr)
                                 call MPI_Type_commit(cg%o_bnd(d, t)%seg(g)%mbc, ierr)

                                 deallocate(sizes, subsizes, starts)
                              enddo
                           endif
                        enddo
                     endif
                  enddo
               endif
            enddo

         else

            do d = xdim, zdim
               if (has_dir(d)) then
                  do t = FLUID, ARR  ! fluid, Bfield, wcr, grav

                     allocate(sizes(dims(t)), subsizes(dims(t)), starts(dims(t)))

                     if (dims(t) == 1+ndims) sizes(1) = nc(t)
                     sizes(dims(t)-zdim+xdim:dims(t)) = cg%n_(:)
                     subsizes(:) = sizes(:)
                     subsizes(dims(t)-zdim+d) = cg%nb

                     starts(:) = 0
                     call MPI_Type_create_subarray(dims(t), sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, cg%mbc(t, d, LO, BND), ierr)
                     call MPI_Type_commit(cg%mbc(t, d, LO, BND), ierr)

                     starts(dims(t)-zdim+d) = cg%nb
                     call MPI_Type_create_subarray(dims(t), sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, cg%mbc(t, d, LO, BLK), ierr)
                     call MPI_Type_commit(cg%mbc(t, d, LO, BLK), ierr)

                     starts(dims(t)-zdim+d) = cg%n_b(d)
                     call MPI_Type_create_subarray(dims(t), sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, cg%mbc(t, d, HI, BLK), ierr)
                     call MPI_Type_commit(cg%mbc(t, d, HI, BLK), ierr)

                     starts(dims(t)-zdim+d) = cg%ijkse(d, HI)
                     call MPI_Type_create_subarray(dims(t), sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, cg%mbc(t, d, HI, BND), ierr)
                     call MPI_Type_commit(cg%mbc(t, d, HI, BND), ierr)

                     deallocate(sizes, subsizes, starts)

                  enddo
               endif
            enddo
         endif

         cgl => cgl%nxt
      enddo

   end subroutine grid_mpi_boundaries_prep

!-----------------------------------------------------------------------------
!
! This routine sets up all guardcells (internal and external) for given rank-3 arrays.
!

   subroutine arr3d_boundaries(pa3d, area_type, dname)

      use constants,  only: ARR, xdim, ydim, zdim, LO, HI, BND, BLK, BND_PER, BND_MPI, BND_SHE, BND_COR, AT_NO_B, I_ONE
      use dataio_pub, only: die, msg
      use domain,     only: has_dir, cdd
      use grid_cont,  only: cg_list_element, grid_container
      use mpi,        only: MPI_REQUEST_NULL, MPI_IN_PLACE, MPI_LOGICAL, MPI_LOR, MPI_COMM_NULL
      use mpisetup,   only: ierr, comm, proc, req, status

      implicit none

      real, dimension(:,:,:), pointer, intent(inout) :: pa3d
      integer(kind=4), intent(in), optional          :: area_type
      character(len=*), intent(in), optional         :: dname

      integer :: i, d
      integer(kind=4) :: lh
      logical :: dodie, do_permpi
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      dodie = .false.

      !> \todo fill corners with big_float ?

      if (ubound(cga%cg_all(:), dim=1) > 1) call die("[grid:arr3d_boundaries] multiple grid pieces per procesor not implemented yet") !nontrivial MPI_Waitall should be outside do while (associated(cgl)) loop

      call cga%get_root(cgl)
      do while (associated(cgl))
         cg => cgl%cg

         if (cdd%comm3d == MPI_COMM_NULL) then

            do_permpi = .true.
            if (present(area_type)) then
               if (area_type /= AT_NO_B) do_permpi = .false.
            endif

            if (do_permpi) call cg%internal_boundaries(ARR, pa3d=pa3d)

         endif

         req(:) = MPI_REQUEST_NULL

         do d = xdim, zdim
            if (has_dir(d)) then
               do lh = LO, HI

                  select case (cg%bnd(d, lh))
                     case (BND_PER)
                        if (cdd%comm3d /= MPI_COMM_NULL) then
                           if (present(area_type)) then
                              if (area_type /= AT_NO_B) cycle
                           endif
                           do i = 1, ceiling(cg%nb/real(cg%n_b(d))) ! Repeating is important for domains that are narrower than their guardcells (e.g. cg%n_b(d) = 2)
                              select case (2*d+lh)
                                 case (2*xdim+LO)
                                    pa3d(1:cg%nb, :, :) = pa3d(cg%ieb:cg%ie, :, :)
                                 case (2*ydim+LO)
                                    pa3d(:, 1:cg%nb, :) = pa3d(:, cg%jeb:cg%je, :)
                                 case (2*zdim+LO)
                                    pa3d(:, :, 1:cg%nb) = pa3d(:, :, cg%keb:cg%ke)
                                 case (2*xdim+HI)
                                    pa3d(cg%ie+1:cg%n_(xdim), :, :) = pa3d(cg%is:cg%isb, :, :)
                                 case (2*ydim+HI)
                                    pa3d(:, cg%je+1:cg%n_(ydim), :) = pa3d(:, cg%js:cg%jsb, :)
                                 case (2*zdim+HI)
                                    pa3d(:, :, cg%ke+1:cg%n_(zdim)) = pa3d(:, :, cg%ks:cg%ksb)
                              end select
                           enddo
                        endif
                     case (BND_MPI)
                        if (cdd%comm3d /= MPI_COMM_NULL) then
                           if (cdd%psize(d) > 1) then
                              call MPI_Isend(pa3d(1, 1, 1), I_ONE, cg%mbc(ARR, d, lh, BLK), cdd%procn(d, lh), int(2*d+(LO+HI-lh), kind=4), cdd%comm3d, req(4*(d-xdim)+1+2*(lh-LO)), ierr)
                              call MPI_Irecv(pa3d(1, 1, 1), I_ONE, cg%mbc(ARR, d, lh, BND), cdd%procn(d, lh), int(2*d+       lh,  kind=4), cdd%comm3d, req(4*(d-xdim)+2+2*(lh-LO)), ierr)
                           else
                              call die("[grid:arr3d_boundaries] bnd_[xyz][lr] == 'mpi' && cdd%psize([xyz]dim) <= 1")
                           endif
                        endif
                     case (BND_SHE) !> \todo move appropriate code from poissonsolver::poisson_solve or do nothing. or die until someone really needs SHEAR.
                        write(msg,*) "[grid:arr3d_boundaries] 'she' not implemented for ",dname
                        dodie = .true.
                     case (BND_COR)
                        if (present(area_type)) then
                           if (area_type /= AT_NO_B) cycle
                        endif
                        write(msg,*) "[grid:arr3d_boundaries] 'cor' not implemented for ", dname
                        dodie = .true.
                     case default ! Set gradient == 0 on the external boundaries
                        if (present(area_type)) then
                           if (area_type /= AT_NO_B) cycle
                        endif
                        do i = 1, cg%nb
                           select case (2*d+lh)
                              case (2*xdim+LO)
                                 pa3d(i, :, :) = pa3d(cg%is, :, :)
                              case (2*ydim+LO)
                                 pa3d(:, i, :) = pa3d(:, cg%js, :)
                              case (2*zdim+LO)
                                 pa3d(:, :, i) = pa3d(:, :, cg%ks)
                              case (2*xdim+HI)
                                 pa3d(cg%ie+i, :, :) = pa3d(cg%ie, :, :)
                              case (2*ydim+HI)
                                 pa3d(:, cg%je+i, :) = pa3d(:, cg%je, :)
                              case (2*zdim+HI)
                                 pa3d(:, :, cg%ke+i) = pa3d(:, :, cg%ke)
                           end select
                        enddo
                  end select

               enddo
            endif
            !> \warning outside xdim-zdim loop MPI_Waitall may change the operations order and as a result may leave mpi-corners uninitiallized
            if (cdd%comm3d /= MPI_COMM_NULL) call MPI_Waitall(size(req(:)), req(:), status(:,:), ierr)
         enddo

         cgl => cgl%nxt
      enddo

      call MPI_Allreduce(MPI_IN_PLACE, dodie, I_ONE, MPI_LOGICAL, MPI_LOR, comm, ierr)
      if (dodie) call die(msg)

   end subroutine arr3d_boundaries

end module grid
