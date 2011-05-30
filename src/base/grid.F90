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

   use grid_cont, only: grid_container

   implicit none

   private
   public :: init_grid, grid_mpi_boundaries_prep, arr3d_boundaries
   public :: total_ncells, cg, D_x, D_y, D_z

   integer, protected :: total_ncells                   !< total number of %grid cells
   integer, protected :: D_x                            !< set to 1 when x-direction exists, 0 otherwise. Use to construct dimensionally-safe indices for arrays
   integer, protected :: D_y                            !< set to 1 when y-direction exists, 0 otherwise.
   integer, protected :: D_z                            !< set to 1 when z-direction exists, 0 otherwise.
   type(grid_container), target :: cg        !< A container for the grid. For AMR this will be a dynamically resized array

contains

!>
!! \brief Routine which sets numbers of cells for the domain, MPI blocks and initializes direction meshes (x,y,z).
!! Also compute domain maximum and minimum of coordinates, lengths of cells and coordinates of zone centers and left/right zone boundaries.
!!
!<
   subroutine init_grid

      use constants,  only: PIERNIK_INIT_MPI, xdim, ydim, zdim
      use dataio_pub, only: printinfo, die, code_progress
      use mpisetup,   only: dom, has_dir

      implicit none

      if (code_progress < PIERNIK_INIT_MPI) call die("[grid:init_grid] MPI not initialized.")

#ifdef VERBOSE
      call printinfo("[grid:init_grid]: commencing...")
#endif /* VERBOSE */

      call cg%init(dom)

      D_x = 0; D_y = 0; D_z = 0
      if (has_dir(xdim)) D_x = 1
      if (has_dir(ydim)) D_y = 1
      if (has_dir(zdim)) D_z = 1
      total_ncells = product(dom%n_d(:))

#ifdef VERBOSE
      call printinfo("[grid:init_grid]: finished. \o/")
#endif /* VERBOSE */

   end subroutine init_grid

!>
!! \brief Set up subsets of u,b and sgp arrays for MPI communication
!!
!! \todo this should be called from cg%init only
!<

   subroutine grid_mpi_boundaries_prep(numfluids, numcrs)

      use dataio_pub, only: die, code_progress
      use constants,  only: PIERNIK_INIT_BASE, FLUID, ARR, xdim, zdim, ndims, LO, HI, BND, BLK, INVALID
      use mpi,        only: MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, MPI_COMM_NULL
      use mpisetup,   only: ierr, has_dir, comm3d, dom, proc, nproc, is_overlap, procmask

      implicit none

      integer, intent(in) :: numfluids !< expect flind%all,     here, cannot grab it directly because of cyclic deps in CR-based setups
      integer, intent(in) :: numcrs    !< expect flind%crs%all, here, cannot grab it directly because of cyclic deps in CR-based setups

      integer, dimension(:), allocatable :: sizes, subsizes, starts
      integer :: d, t, g, hl, lh, j
      integer, dimension(FLUID:ARR) :: nc
      integer, parameter, dimension(FLUID:ARR) :: dims = [ 1+ndims, 1+ndims, 1+ndims, ndims ] !< dimensionality of arrays

      integer(kind=8), dimension(xdim:zdim) :: ijks, per
      integer(kind=8), dimension(xdim:zdim, LO:HI) :: b_layer, bp_layer, poff
      logical :: sharing

      if (code_progress < PIERNIK_INIT_BASE) call die("[grid:grid_mpi_boundaries_prep] grid or fluids not initialized.")

      nc = [ numfluids, ndims, max(numcrs,1), 1 ]      !< number of fluids, magnetic field components, CRs, and 1 for rank-3 array

      ! find neighbours and set up the MPI containers
      if (comm3d == MPI_COMM_NULL) then

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
                  b_layer(:,:) = dom%se(proc, :, :)
                  b_layer(d, lh) = b_layer(d, lh) + lh-hl ! -1 for LO, +1 for HI
                  b_layer(d, hl) = b_layer(d, lh) ! boundary layer without corners
                  do j = 0, nproc-1
                     call is_overlap(b_layer(:,:), dom%se(j, :, :), sharing, per(:))
                     if (sharing) procmask(j) = procmask(j) + 1
                  enddo
               enddo
               do j = FLUID, ARR
                  allocate(cg%i_bnd(d, j)%seg(sum(procmask(:))))
                  allocate(cg%o_bnd(d, j)%seg(sum(procmask(:))))
               enddo

               ! set up segments to be sent or received
               g = 0
               do j = 0, nproc-1
                  if (procmask(j) /= 0) then
                     do lh = LO, HI
                        hl = LO+HI-lh
                        b_layer(:,:) = dom%se(proc, :, :)
                        b_layer(d, lh) = b_layer(d, lh) + lh-hl
                        b_layer(d, hl) = b_layer(d, lh)

                        bp_layer(:, :) = b_layer(:, :)
                        where (per(:) > 0)
                           bp_layer(:, LO) = mod(b_layer(:, LO) + per(:), per(:))
                           bp_layer(:, HI) = mod(b_layer(:, HI) + per(:), per(:))
                        endwhere
                        call is_overlap(bp_layer(:,:), dom%se(j, :, :), sharing)

                        if (sharing) then
                           poff(:,:) = bp_layer(:,:) - b_layer(:,:) ! displacement due to periodicity
                           bp_layer(:, LO) = max(bp_layer(:, LO), dom%se(j, :, LO))
                           bp_layer(:, HI) = min(bp_layer(:, HI), dom%se(j, :, HI))

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
                              sizes   (dims(t)-zdim+xdim:dims(t)) = [ cg%nx, cg%ny, cg%nz ]
                              subsizes(dims(t)-zdim+xdim:dims(t)) = int(cg%i_bnd(d, t)%seg(g)%se(:, HI) - cg%i_bnd(d, t)%seg(g)%se(:, LO) + 1, kind=4)
                              starts  (dims(t)-zdim+xdim:dims(t)) = int(cg%i_bnd(d, t)%seg(g)%se(:, LO), kind=4)-1
                              call MPI_Type_create_subarray(dims(t), sizes, subsizes, starts,  MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, cg%i_bnd(d, t)%seg(g)%mbc, ierr)
                              call MPI_Type_commit(cg%i_bnd(d, t)%seg(g)%mbc, ierr)

                              subsizes(dims(t)-zdim+xdim:dims(t)) = int(cg%o_bnd(d, t)%seg(g)%se(:, HI) - cg%o_bnd(d, t)%seg(g)%se(:, LO) + 1, kind=4)
                              starts  (dims(t)-zdim+xdim:dims(t)) = int(cg%o_bnd(d, t)%seg(g)%se(:, LO), kind=4)-1
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
                  sizes(dims(t)-zdim+xdim:dims(t)) = [ cg%nx, cg%ny, cg%nz ]
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

   end subroutine grid_mpi_boundaries_prep

!-----------------------------------------------------------------------------
!
! This routine sets up all guardcells (internal and external) for given rank-3 arrays.
!

   subroutine arr3d_boundaries(pa3d, area_type, dname)

      use constants,  only: ARR, xdim, ydim, zdim, LO, HI, BND, BLK, BND_PER, BND_MPI, BND_SHE, BND_COR, AT_NO_B
      use dataio_pub, only: die, msg
      use mpi,        only: MPI_REQUEST_NULL, MPI_IN_PLACE, MPI_LOGICAL, MPI_LOR, MPI_COMM_NULL
      use mpisetup,   only: ierr, has_dir, psize, procn, comm, comm3d, proc, req, status

      implicit none

      real, dimension(:,:,:), pointer, intent(inout) :: pa3d
      integer, intent(in), optional                  :: area_type
      character(len=*), intent(in), optional         :: dname

      integer :: i, d, lh
      logical :: dodie, do_permpi

      dodie = .false.

      !> \todo fill corners with big_float ?

      if (comm3d == MPI_COMM_NULL) then

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
                     if (comm3d /= MPI_COMM_NULL) then
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
                                 pa3d(cg%ie+1:cg%nx, :, :) = pa3d(cg%is:cg%isb, :, :)
                              case (2*ydim+HI)
                                 pa3d(:, cg%je+1:cg%ny, :) = pa3d(:, cg%js:cg%jsb, :)
                              case (2*zdim+HI)
                                 pa3d(:, :, cg%ke+1:cg%nz) = pa3d(:, :, cg%ks:cg%ksb)
                           end select
                        enddo
                     endif
                  case (BND_MPI)
                     if (comm3d /= MPI_COMM_NULL) then
                        if (psize(d) > 1) then
                           call MPI_Isend(pa3d(1, 1, 1), 1, cg%mbc(ARR, d, lh, BLK), procn(d, lh), 2*d+(LO+HI-lh), comm3d, req(4*(d-xdim)+1+2*(lh-LO)), ierr)
                           call MPI_Irecv(pa3d(1, 1, 1), 1, cg%mbc(ARR, d, lh, BND), procn(d, lh), 2*d+       lh,  comm3d, req(4*(d-xdim)+2+2*(lh-LO)), ierr)
                        else
                           call die("[grid:arr3d_boundaries] bnd_[xyz][lr] == 'mpi' && psize([xyz]dim) <= 1")
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
         if (comm3d /= MPI_COMM_NULL) call MPI_Waitall(size(req(:)), req(:), status(:,:), ierr)
      enddo

      call MPI_Allreduce(MPI_IN_PLACE, dodie, 1, MPI_LOGICAL, MPI_LOR, comm, ierr)
      if (dodie) call die(msg)

   end subroutine arr3d_boundaries

end module grid
