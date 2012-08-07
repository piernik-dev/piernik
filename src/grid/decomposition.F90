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
!! \brief Domain decomposition routines and variables
!!
!! \details This module contains everything closely related to grid redistribution over processes
!<

module decomposition

   use constants, only: ndims
   use primes,    only: primes_T
   use cart_comm, only: cart_decomposition

   implicit none

   private
   public :: cleanup_decomposition, init_decomposition, set_pse_sel, box_T, cuboids, cdd

   type(cart_decomposition) :: cdd !< Cartesian Domain Decomposition stuff

   type :: cuboids
      integer(kind=8), dimension(:,:,:), allocatable :: sel !< list of grid chunks (:, xdim:zdim, LO:HI)
   end type cuboids

   !> \brief A box (or rectangle) within a certain refinement level to be decomposed into smaller pieces
   type :: box_T
      integer(kind=4), dimension(ndims) :: n_d          !< number of grid cells
      integer(kind=8), dimension(ndims) :: off          !< offset (with respect to the base level, counted on own level)
      type(cuboids), dimension(:), allocatable :: pse   !< lists of grid chunks on each process (FIRST:LAST); Use with care, because this is an antiparallel thing

    contains

      procedure          :: decompose_patch             !< Main wrapper for a block decomposer
      procedure, private :: decompose_patch_int         !< Compute allowed domain decomposition

      procedure, private :: cartesian_tiling            !< Decomposes the box into a topologically cartesian grid
      procedure, private :: choppy_tiling               !< Less structured box decomposition
      procedure, private :: stamp_cg                    !< Divide the box into lots of identical blocks
      procedure :: allocate_pse                !< allocate the segment list
      procedure :: deallocate_pse              !< deallocate the segment list
   end type box_T

   ! Private variables
   type(primes_T) :: primes
   real :: ideal_bsize

contains

!> \brief Initialize the cdd structure and find some prime numbers

   subroutine init_decomposition

      use mpisetup,  only: nproc

      implicit none

      call cdd%init            ! cdd% will contain valid values if and only if comm3d becomes valid communicator
      call primes%sieve(nproc) ! it is possible to use primes only to sqrt(nproc), but it is easier to have the full table. Cheap for any reasonable nproc.

   end subroutine init_decomposition

!> \brief Free the resources

   subroutine cleanup_decomposition

      implicit none

      call cdd%cleanup
      call primes%erase

   end subroutine cleanup_decomposition

!> \brief Main wrapper for a block decomposer

   logical function decompose_patch(this, n_d, off, n_pieces) result(patch_divided)

      use constants, only: I_ONE
      use mpi,       only: MPI_IN_PLACE, MPI_LOGICAL, MPI_LAND
      use mpisetup,  only: comm, mpi_err

      implicit none

      class(box_T),                      intent(inout) :: this     !< the patch, which we want to be chopped into pieces
      integer(kind=4), dimension(ndims), intent(in)    :: n_d      !< number of grid cells
      integer(kind=8), dimension(ndims), intent(in)    :: off      !< offset (with respect to the base level, counted on own level), \todo make use of it
      integer, optional,                 intent(in)    :: n_pieces !< how many pieces the patch should be divided to?

      this%n_d(:) = n_d(:)
      this%off(:) = off(:)

      call this%decompose_patch_int(patch_divided, n_pieces)
      if (patch_divided) patch_divided = is_not_too_small(this, "not catched anywhere")
      call MPI_Allreduce(MPI_IN_PLACE, patch_divided, I_ONE, MPI_LOGICAL, MPI_LAND, comm, mpi_err)

   end function decompose_patch

!> \brief This routine computes (hopefully close-to-optimal) allowed domain decomposition

   subroutine decompose_patch_int(patch, patch_divided, n_pieces)

      use constants,  only: ndims
      use dataio_pub, only: warn, printinfo, msg
      use domain,     only: dom, psize, bsize, allow_noncart, allow_uneven, dd_rect_quality, dd_unif_quality, minsize, use_comm3d
      use mpisetup,   only: nproc, master, have_mpi

      implicit none

      class(box_T),      intent(inout) :: patch         !< the patch, which we want to be chopped into pieces
      logical,           intent(out)   :: patch_divided !< Set to .true. after a successful decomposition
      integer, optional, intent(in)    :: n_pieces      !< how many pieces the patch should be divided to?

      real :: quality
      integer(kind=4), dimension(ndims) :: p_size
      integer :: pieces, ml

      pieces = nproc
      if (present(n_pieces)) pieces = n_pieces
      call primes%sieve(pieces)

      patch_divided = .false.

      ! Try the decomposition into same-size blocks
      if (all(bsize(:) > 0 .or. .not. dom%has_dir(:)) .and. .not. present(n_pieces)) then
         call patch%stamp_cg
         patch_divided = allocated(patch%pse)
         if (patch_divided) patch_divided = is_not_too_small(patch, "stamp_cg")
         if (patch_divided) return
      endif

      ! Try the cartesian decomposition, specified in problem.par
      if (product(psize(:)) == pieces) then
         if (all(mod(patch%n_d(:), int(psize(:), kind=4)) == 0)) then
            if (master .and. have_mpi) then
               write(msg,'(a,3i4,a,3i6,a)')"[decomposition:decompose_patch_int] Domain divided to [",psize(:)," ] pieces, each of [",patch%n_d(:)/psize(:)," ] cells."
               call printinfo(msg)
            endif
            call patch%cartesian_tiling(psize(:), pieces)
            patch_divided = is_not_too_small(patch, "cartesian_tiling")
            if (patch_divided) return
         else
            write(msg,'(a,3i6,a,3i4,a)')"[decomposition:decompose_patch_int] Cannot divide domain with [",patch%n_d(:)," ] cells to [",psize(:)," ] piecess. "
            if (master) call warn(msg)
         endif
      endif

      ! this is the minimal total area of internal boundaries (periodic case), achievable for some perfect domain divisions
      ideal_bsize = dom%eff_dim * (pieces * product(real(patch%n_d(:)))**(dom%eff_dim-1))**(1./dom%eff_dim)

      ! Try to find a close-to-optimal cartesian decomposition into same-sized blocks
      call decompose_patch_uniform(p_size(:), patch%n_d, pieces)
      if (product(p_size(:)) == pieces) then
         quality = ideal_bsize / sum(p_size(:)/real(patch%n_d(:)) * product(real(patch%n_d(:))), MASK = patch%n_d(:) > 1)
         if (quality >= dd_unif_quality .or. .not. (allow_uneven .or. allow_noncart)) then
            call patch%cartesian_tiling(p_size(:), pieces)
            patch_divided = is_not_too_small(patch, "decompose_patch_uniform")
            if (patch_divided) return
         else
            write(msg,'(2(a,f6.3),a)')"[decomposition:decompose_patch_int] Quality of uniform division = ",quality," is below threshold ",dd_unif_quality, ", trying harder ..."
            if (master) call warn(msg)
         endif
      endif

      ! Try to find a close-to-optimal cartesian decomposition into similar-sized blocks
      if (allow_uneven) then
         call decompose_patch_rectlinear(p_size(:), patch%n_d, pieces)
         quality = 1 !< \todo make an estimate
         if (product(p_size(:)) == pieces) then
            if (quality > dd_rect_quality .or. .not. allow_noncart) then
               call patch%cartesian_tiling(p_size(:), pieces)
               patch_divided = is_not_too_small(patch, "decompose_patch_rectlinear")
               if (patch_divided) return
            endif
         endif
         if (master) call warn("[decomposition:decompose_patch_int] decompose_patch_rectlinear failed")
      else
         if (master) call warn("[decomposition:decompose_patch_int] Did not try uneven domain division")
      endif

      ! Try to find a close-to-optimal decomposition into similar-volume blocks, minimize the boundary area
      if (allow_noncart) then
         p_size(:) = psize(:)
         call decompose_patch_slices(p_size(:), patch%n_d, pieces)
         ! if good_enough then return
         call patch%choppy_tiling(p_size(:), pieces)
         patch_divided = is_not_too_small(patch, "decompose_patch_slices")
         if (patch_divided) return
         if (master) call warn("[decomposition:decompose_patch_int] decompose_patch_slices failed")
      else
         if (master) call warn("[decomposition:decompose_patch_int] Did not try non-cartesian domain division")
      endif

      ! The domain is probably too small for given number of processes, decompose the domain into smallest possible pieces and leave some processes unemployed
      if (use_comm3d) then
         ! We have only single cartesian communicator, so we divide either to nproc or put everything on master
         ! Putting the only piece of the grid on master is valid as long as local boundaries (periodic case) are done by copying memory, not mpi communication
         p_size(:) = 1
      else
         p_size(:) = patch%n_d(:) / minsize(:)
         do while (product(p_size(:)) > nproc)
            ml = maxloc(p_size(:), dim=1)
            if (p_size(ml) > 1) p_size(ml) = p_size(ml) - 1
         enddo
      endif
      call patch%cartesian_tiling(p_size(:), product(p_size(:)))
      patch_divided = is_not_too_small(patch, "decompose_patch_cartesian_less_than_nproc")
      if (patch_divided) return

      ! Everything failed
      write(msg,'(a,3i6,a,i4,a)') "[decomposition:decompose_patch_int] Cannot divide domain with [",patch%n_d(:)," ] cells to ",pieces," piecess. "
      if (master) call warn(msg) ! should die

   end subroutine decompose_patch_int

!>
!! \brief Decomposes the domain into topologically Cartesian grid
!!
!! \details Each process gets a single piece of the grid. There are at most 6 neighbours; each boundary is either external boundary or is shared with one neighbour.
!<

   subroutine cartesian_tiling(patch, p_size, pieces)

      use constants,  only: xdim, ydim, ndims, LO, HI, I_ONE
      use dataio_pub, only: printinfo, die
      use domain,     only: dom, use_comm3d
      use mpi,        only: MPI_COMM_NULL
      use mpisetup,   only: master, proc, nproc, mpi_err

      implicit none

      class(box_T),                      intent(inout) :: patch   !< the patch, which we want to be chopped into pieces
      integer(kind=4), dimension(ndims), intent(in)    :: p_size
      integer,                           intent(in)    :: pieces      !< number of pieces

      integer(kind=4) :: p
      integer(kind=4), dimension(ndims) :: pc
      integer, dimension(nproc) :: n_cg ! BEWARE: antiparallel

      if (product(p_size(:)) /= pieces) call die("[decomposition:cartesian_tiling] product(p_size(:)) /= pieces")

      if (pieces > nproc) call die("[decomposition:cartesian_tiling] cartesian decomposition into more pieces than processes not implemented yet")
      n_cg(:) = 0
      n_cg(1:pieces) = 1
      call patch%allocate_pse(n_cg(:))

      if (use_comm3d) then
         if (pieces > 1 .or. nproc == 1) then
            if (pieces /= nproc) call die("[decomposition:cartesian_tiling] use_comm3d is incompatible with decomposition to a number of pieces different than number of processes")
            call cdd%init_cart(p_size)
            !> \warning Perhaps cdd%init_cart should be called only once.
            if (master) call printinfo("[decomposition:cartesian_tiling] Cartesian decomposition with comm3d")
         else
            if (master) call printinfo("[decomposition:cartesian_tiling] Whole grid on the master PE")
         endif
      else
         if (master) call printinfo("[decomposition:cartesian_tiling] Cartesian decomposition without comm3d")
      endif

      do p = 0, pieces-1
         if (cdd%comm3d == MPI_COMM_NULL) then
            if (use_comm3d .and. nproc == pieces) call die("[decomposition:cartesian_tiling] MPI_Cart_create failed")
            pc(:) = [ mod(p, p_size(xdim)), mod(p/p_size(xdim), p_size(ydim)), p/product(p_size(xdim:ydim)) ]
         else
            if (pieces == nproc) then
               call MPI_Cart_coords(cdd%comm3d, p, ndims, pc, mpi_err)
            else
               if (pieces /= 1) call die("decomposition:cartesian_tiling] pieces /= nproc and pieces /= 1 and cdd%comm3d set")
               pc(:) = 0
            endif
         endif
         where (dom%has_dir(:))
            patch%pse(p)%sel(I_ONE, :, LO) = (patch%n_d(:) *  pc(:) ) / p_size(:)     ! offset of low boundaries of the local domain (0 at low external boundaries)
            patch%pse(p)%sel(I_ONE, :, HI) = (patch%n_d(:) * (pc(:)+1))/p_size(:) - 1 ! offset of high boundaries of the local domain (n_d(:) - 1 at right external boundaries)
         endwhere
      enddo

      if (ubound(patch%pse(proc)%sel(:,:,:), dim=1) > 1) call die("[decomposition:cartesian_tiling] No multiblock support.")

   end subroutine cartesian_tiling

!>
!! \brief Less structured domain decomposition.
!!
!! \details Each process gets a single piece of the grid. Each non-external boundary can be shared wit one or more processes.
!! When pieces == product(p_size(:)), the decomposition is identical to the result of cartesian_tiling.
!>

   subroutine choppy_tiling(patch, p_size, pieces)

      use constants,  only: ndims, xdim, ydim, zdim, LO, HI, I_ZERO, I_ONE
      use dataio_pub, only: printinfo, die
      use mpisetup,   only: master, proc

      implicit none

      class(box_T),                      intent(inout) :: patch   !< the patch, which we want to be chopped into pieces
      integer(kind=4), dimension(ndims), intent(in)    :: p_size
      integer,                           intent(in)    :: pieces      !< number of pieces

      integer(kind=4) :: p, px, py
      integer(kind=4), dimension(:), allocatable :: pz_slab, py_slab

      call patch%allocate_pse

      if (master) call printinfo("[decomposition:choppy_tiling] Non-cartesian decomposition (no comm3d)")
      allocate(pz_slab(p_size(zdim) + 1))
      pz_slab(1) = I_ZERO
      do p = I_ONE, p_size(zdim)
         pz_slab(p+1) = pz_slab(p) + pieces / p_size(zdim)
         if (p <= mod(pieces, p_size(zdim))) pz_slab(p+1) = pz_slab(p+1) + I_ONE ! longer slabs go first
      enddo
      do p = I_ONE, p_size(zdim)
         do px = pz_slab(p), pz_slab(p+1)-I_ONE
            patch%pse(px)%sel(I_ONE, zdim, LO) = nint((patch%n_d(zdim) *  pz_slab(p)   ) / real(pieces))
            patch%pse(px)%sel(I_ONE, zdim, HI) = nint((patch%n_d(zdim) *  pz_slab(p+1) ) / real(pieces)) - 1
         enddo
         allocate(py_slab(p_size(ydim) + 1))
         py_slab(1) = I_ZERO
         do py = I_ONE, p_size(ydim)
            py_slab(py+1) = py_slab(py) + (pz_slab(p+1)-pz_slab(p)) / p_size(ydim)
            if (py <= mod((pz_slab(p+1)-pz_slab(p)), p_size(ydim))) py_slab(py+1) = py_slab(py+1) + I_ONE ! longer slabs go first
         enddo
         do py = I_ONE, p_size(ydim)
            do px = pz_slab(p)+py_slab(py), pz_slab(p)+py_slab(py+1) - I_ONE
               patch%pse(px)%sel(I_ONE, ydim, LO) = nint((patch%n_d(ydim) *  py_slab(py)   ) / real(pz_slab(p+1)-pz_slab(p)))
               patch%pse(px)%sel(I_ONE, ydim, HI) = nint((patch%n_d(ydim) *  py_slab(py+1) ) / real(pz_slab(p+1)-pz_slab(p))) - I_ONE
            enddo
            do px = I_ZERO, py_slab(py+1)-py_slab(py) - I_ONE
               patch%pse(pz_slab(p)+py_slab(py)+px)%sel(I_ONE, xdim, LO) = (patch%n_d(xdim) *  px    ) / (py_slab(py+1)-py_slab(py))
               patch%pse(pz_slab(p)+py_slab(py)+px)%sel(I_ONE, xdim, HI) = (patch%n_d(xdim) * (px+1) ) / (py_slab(py+1)-py_slab(py)) - 1 ! no need to sort lengths here
            enddo
         enddo
         if (allocated(py_slab)) deallocate(py_slab)
      enddo
      if (allocated(pz_slab)) deallocate(pz_slab)

      if (ubound(patch%pse(proc)%sel(:,:,:), dim=1) > 1) call die("[decomposition:choppy_tiling] No multiblock support.")

   end subroutine choppy_tiling

!>
!! \brief This routine tries to divide the computational domain into local domains.
!!
!! \details The goal is to minimize the ratio of longest to shortest edge to minimize the amount of inter-process communication.
!! If the benchmarks show that some direction should be partitioned in more pieces than other directions, implement appropriate weighting in j1, j2 and j3 calculations.
!!
!! For some weird domains and PE counts this routine may find tiling that does not satisfy multigrid restrictions even if there exists an acceptable tiling.
!! In such case the user must divide domain manually by providing psize(:) parameters through problem.par.
!!
!! \attention Must be called by all procs to avoid communication and ensure that every proc has proper psize
!<
   subroutine decompose_patch_uniform(p_size, n_d, pieces)

      use constants,  only: xdim, zdim, ndims
      use dataio_pub, only: warn, printinfo, msg
      use domain,     only: dom
      use mpisetup,   only: master

      implicit none

      integer(kind=4), dimension(ndims), intent(out) :: p_size
      integer(kind=4), dimension(ndims), intent(in)  :: n_d         !< size of the box to be divided
      integer,                           intent(in)  :: pieces      !< number of pieces

      integer(kind=4) :: n
      integer :: j1, j2, j3, jj, p
      integer(kind=4), dimension(ndims) :: ldom, tmp

      ldom(xdim:zdim) = int(n_d(zdim:xdim:-1), kind=4) ! Maxloc returns first occurrence of max, reversing direction order (to ZYX) gives better cache utilization.
      n = pieces
      p_size(:) = 1
      if (n == 1) return

      do p = size(primes%tab), 1, -1 ! start from largest defined primes, continue down to 2
         do while (mod(n, primes%tab(p))==0)

            jj = 0
            j1 = sum(maxloc(ldom), 1) ! First try the longest edge; note the trick to make a scalar from 1-element vector without assignment to another variable
            if (mod(ldom(j1), primes%tab(p))==0) then
               jj = j1
            else
               j2 = 1 + mod(j1 + 0, int(ndims))
               j3 = 1 + mod(j1 + ndims -2, int(ndims))
               if (ldom(j2) > ldom(j3)) then
                  j2 = 1 + mod(j1 + ndims -2, int(ndims))
                  j3 = 1 + mod(j1 + 0, int(ndims))
               endif
               if (mod(ldom(j2), primes%tab(p))==0) jj = j2 ! middle edge ...
               if (jj == 0 .and. mod(ldom(j3), primes%tab(p))==0) jj = j3 ! try the shortest edge on last resort
            endif

            if (jj == 0) then
               if (master) call warn("[decomposition:decompose_patch_uniform]: Can't find divisible edge")
               p_size(:) = 1
               return
            else
               p_size(jj) = p_size(jj) * primes%tab(p)
               n          = n          / primes%tab(p)
               ldom(jj)   = ldom(jj)   / primes%tab(p)
            endif

         enddo
      enddo

      if (any(ldom(:) < dom%nb .and. dom%has_dir(:)) .or. n /= 1) then
         if (master) call warn("[decomposition:decompose_patch_uniform]: I am not that intelligent") ! pieces has too big prime factors
         p_size(:) = 1
         return
      endif

      tmp(xdim:zdim) = p_size(zdim:xdim:-1) ! directions were reverted at ldom assignment
      p_size(:) = tmp(:)

      if (master) then
         write(msg,'(a,3i4,a,3i6,a)')"[decomposition:decompose_patch_uniform] Domain divided to [",p_size(:)," ] pieces, each of [",ldom(zdim:xdim:-1)," ] cells."
         call printinfo(msg)
      endif

   end subroutine decompose_patch_uniform

!>
!! \brief Divide the computational domain into local domains. Allow their size to change by +/- 1 depending on CPU rank (this will introduce some load imbalance)
!! if it is not possible to divide an edge evenly. Try to minimize the imbalance and total internal boundaries size.
!<
   subroutine decompose_patch_rectlinear(p_size, n_d, pieces)

      use constants,  only: xdim, ydim, ndims, I_ZERO, I_ONE
      use dataio_pub, only: printinfo, msg
      use domain,     only: dom, is_uneven
      use mpisetup,   only: master

      implicit none

      integer(kind=4), dimension(ndims), intent(out) :: p_size
      integer(kind=4), dimension(ndims), intent(in)  :: n_d         !< size of the box to be divided
      integer,                           intent(in)  :: pieces      !< number of pieces

      real, parameter :: b_load_fac = 0.25 ! estimated increase of execution time after doubling the total size of internal boundaries.
      ! \todo estimate this factor for massively parallel runs and for Intel processors

      integer(kind=4), allocatable, dimension(:) :: ppow
      integer(kind=4), allocatable, dimension(:,:) :: fac
      integer(kind=4), dimension(ndims) :: ldom
      integer(kind=4) :: p, i, j, k, nf
      integer :: n, ii, bsize
      real :: load_balance, best, quality

      p_size(:) = 1
      if (pieces == 1) return

      allocate(ppow(size(primes%tab)))

      p = pieces
      do i = I_ONE, int(size(primes%tab), kind=4)
         ppow(i) = 0
         do while (mod(p, primes%tab(i)) == 0)
            ppow(i) = ppow(i) + I_ONE
            p = p / primes%tab(i)
         enddo
      enddo

      nf = int(count(ppow(:) > 0), kind=4)
      allocate(fac(nf,3))
      j = I_ONE
      do i = I_ONE, int(size(primes%tab), kind=4)
         if (ppow(i)>0) then
            ! prime, its power and number of different decompositions in three dimensions for this prime
            fac(j,:) = [ primes%tab(i), ppow(i), int((ppow(i)+1)*(ppow(i)+2)/2, kind=4) ]
            j = j + I_ONE
         endif
      enddo
      deallocate(ppow)

      best = 0.
      ii = 0
      do while (all(fac(:,3) > 0))
         ldom(:) = 1
         do n = 1, nf ! find an unique decomposition of fac(n,2) into [i,j,k], all([i,j,k] >= 0) && i+j+k = fac(n,2). The decompositions are enumerated with fac(n,3).
            i = int(sqrt(1./4.+2.*(fac(n,3)-1)) - 1./2., kind=4) ! i and k enumerate a point in a triangle: (i>=0 && k>=0 && i+k<=fac(n,2))
            k = fac(n,3) - int(1 + i*(i+1)/2, kind=4)
            i = fac(n,2) - i
            j = fac(n,2) - (i + k)
            ldom(:) = ldom(:) * fac(n,1)**[i, j, k]
         enddo

         bsize = int(sum(ldom(:)/real(n_d(:), kind=8) * product(int(n_d(:), kind=8)), MASK = n_d(:) > 1)) !ldom(1)*n_d(2)*n_d(3) + ldom(2)*n_d(1)*n_d(3) + ldom(3)*n_d(1)*n_d(2)
         load_balance = product(real(n_d(:))) / ( real(pieces) * product( int((n_d(:)-1)/ldom(:)) + 1 ) )

         quality = load_balance/ (1 + b_load_fac*(bsize/ideal_bsize - 1.))
         ! \todo add a factor that estimates lower cost when x-direction is not chopped too much
         quality = quality * (1. - (0.001 * ldom(xdim) + 0.0001 * ldom(ydim))/pieces) ! \deprecated estimate these magic numbers

         if (any(ldom(:) > n_d(:))) quality = 0
         if (any(n_d(:)/ldom(:) < dom%nb .and. dom%has_dir(:))) quality = 0

#ifdef DEBUG
         if (quality > 0 .and. master) then
            ii = ii + 1
            write(msg,'(a,i3,a,3i4,a,i10,2(a,f10.7))')"m:ddr ",ii," p_size= [",ldom(:)," ], bcells= ", bsize, ", balance = ", load_balance, ", est_quality = ", quality
            call printinfo(msg)
         endif
#endif /* DEBUG */
         if (quality > best) then
            best = quality
            p_size(:) = ldom(:)
         endif
         do j = I_ONE, nf ! search for next unique combination
            if (fac(j,3) > 1) then
               fac(j,3) = fac(j,3) - I_ONE
               exit
            else
               if (j<nf) then
                  fac(j,3) = int((fac(j,2)+1)*(fac(j,2)+2)/2, kind=4)
               else
                  fac(:,3) = I_ZERO ! no more combinations to try
               endif
            endif
         enddo
      enddo

      deallocate(fac)

      is_uneven = any(mod(n_d(:), int(p_size(:), kind=4)) /= 0)

      if (master) then
#ifdef DEBUG
         write(msg,'(a,3f10.2,a,i10)')"m:ddr id p_size = [",(pieces/product(real(n_d(:), kind=8)))**(1./dom%eff_dim)*n_d(:),"], bcells= ", int(ideal_bsize)
         call printinfo(msg)
#endif /* DEBUG */
         write(msg,'(a,3i4,a)')      "[decomposition:decompose_patch_rectlinear] Domain divided to [",p_size(:)," ] pieces"
         call printinfo(msg)
         if (is_uneven) then
            write(msg,'(2(a,3i5),a)')"                                         Sizes are from [", int(n_d(:)/p_size(:))," ] to [",int((n_d(:)-1)/p_size(:))+1," ] cells."
            call printinfo(msg)
            write(msg,'(a,f8.5)')    "                                         Load balance is ",product(int(n_d(:), kind=8)) / ( real(pieces, kind=8) * product( int((n_d(:)-1)/p_size(:)) + 1 ) )
         else
            write(msg,'(a,3i5,a)')   "                                         Size is [", int(n_d(:)/p_size(:))," ] cells."
         endif
         call printinfo(msg)
      endif

   end subroutine decompose_patch_rectlinear

!>
!! \brief Divide the computational domain into local domains. Allow their size to depend significantly on CPU rank and allow for more than one neighbour on a single boundary.
!! Try to minimize the imbalance and total internal boundaries size.
!<
   subroutine decompose_patch_slices(p_size, n_d, pieces)

      use constants,  only: xdim, ydim, zdim, ndims, I_ONE
      use dataio_pub, only: msg, printinfo, warn
      use domain,     only: dom, is_mpi_noncart, is_uneven
      use mpisetup,   only: master

      implicit none

      integer(kind=4), dimension(ndims), intent(inout) :: p_size
      integer(kind=4), dimension(ndims), intent(in)    :: n_d         !< size of the box to be divided
      integer,                           intent(in)    :: pieces      !< number of pieces

      real, parameter :: minfac = 1.3 ! prevent domain division to halves if cell count in a given direction is too low. (not verified for optimality)
      real :: optc

      !> \todo Try to make an intelligent guess for slicing, then go down to the local minimum and explore neighbourhood. Exploring all possibilities is an O(pieces)**2 task
      ! The best solution is probably near (pieces/product(real(n_d(:), kind=8)))**(1./dom%eff_dim)*n_d(:)

      is_mpi_noncart = .true.
      is_uneven = .true.

      if (all(p_size(ydim:zdim) == 1)) then
         if (dom%has_dir(zdim)) then
            optc = max(real(dom%nb), (product(int(n_d(:), kind=8))/real(pieces)) ** (1./dom%eff_dim)) ! number of cells for ideal cubes
            if (n_d(zdim) > minfac*optc) p_size(zdim) = int(ceiling(n_d(zdim)/optc), kind=4)
         endif
         if (dom%has_dir(ydim)) then
            optc = max(real(dom%nb), (product(int(n_d(xdim:ydim), kind=8))*p_size(zdim)/real(pieces)) ** (1./count(dom%has_dir(xdim:ydim))))
            if (n_d(ydim) > minfac*optc) p_size(ydim) = int(ceiling(n_d(ydim)/optc), kind=4)
         endif
      endif
      if (dom%has_dir(xdim)) p_size(xdim) = (pieces - I_ONE)/(p_size(ydim)*p_size(zdim)) + I_ONE !sometimes it might be less by 1

      where (.not. dom%has_dir(:)) p_size(:) = 1 ! just in case
      do while (product(p_size(:)) < pieces)
         write(msg,'(a,3i4,a)') "[decomposition:decompose_patch_slices] imperfect noncartesian division to [",p_size(:)," ] pieces"
         if (master) call warn(msg)
         if (dom%has_dir(xdim)) then
            p_size(xdim) = p_size(xdim) + I_ONE
         else if (dom%has_dir(ydim)) then
            p_size(ydim) = p_size(ydim) + I_ONE
         else
            p_size(zdim) = p_size(zdim) + I_ONE
         endif
      enddo
      write(msg,'(a,3i4,a)') "[decomposition:decompose_patch_slices] performed noncartesian division to [",p_size(:)," ] pieces"
      if (master) call printinfo(msg)

   end subroutine decompose_patch_slices

!> \brief Divide the domain into lots of identical blocks

   subroutine stamp_cg(patch)

      use constants,  only: xdim, zdim, LO, HI
      use dataio_pub, only: warn, msg, die
      use domain,     only: dom, bsize
      use mpisetup,   only: master, proc, nproc

      implicit none

      class(box_T), intent(inout) :: patch  !< the patch, which we want to be chopped into pieces

      integer, dimension(xdim:zdim) :: n_bl, b_loc
      integer, dimension(xdim:zdim, LO:HI) :: se
      integer :: tot_bl, loc_bl, bl_s, bl_e
      integer, dimension(nproc) :: n_cg ! BEWARE: antiparallel
      integer, allocatable, dimension(:) :: pb ! BEWARE: antiparallel
      integer :: p, b

      if (any(bsize(xdim:zdim) <=0)) then
         if (master) call warn("[initproblem:stamp_cg] some(bsize(1:3)) <=0")
         return
      endif

      if (any(mod(patch%n_d(:), int(bsize(xdim:zdim), kind=4)) /= 0 .and. dom%has_dir(:))) then
         write(msg,'(a,3f10.3,a)')"[initproblem:stamp_cg] Fractional number of blocks: n_d(:)/bsize(1:3) = [",patch%n_d(:)/real(bsize(xdim:zdim)),"]"
         if (master) call warn(msg)
         return
      endif

      where (dom%has_dir(:))
         n_bl(:) = int(patch%n_d(:) / bsize(xdim:zdim), kind=4)
      elsewhere
         n_bl(:) = 1
      endwhere
      tot_bl = product(n_bl(:), mask=dom%has_dir(:))
      if (allocated(pb)) call die("[initproblem:stamp_cg] pb already allocated")
      allocate(pb(tot_bl))

      if (tot_bl < nproc) then
         if (master) call warn("[initproblem:stamp_cg] Found more processes than available blocks. Reverting to some other method.")
         return
      endif

      do p = 0, nproc-1
         n_cg(p+1) = (((p + 1) * tot_bl)/nproc) - ((p * tot_bl)/nproc) ! number of blocks on process p
      enddo
      call patch%allocate_pse(n_cg(:))

      bl_s = (proc * tot_bl)/nproc + 1
      bl_e = ((proc + 1) * tot_bl)/nproc
      loc_bl = bl_e - bl_s + 1
      if (loc_bl /= n_cg(proc + 1)) call die("[initproblem:stamp_cg] loc_bl /= n_cg(proc + 1)")

      do p = 0, nproc-1
         n_cg(p+1) = (p * tot_bl)/nproc + 1 ! first block on process p
         if (p>0) pb(n_cg(p):n_cg(p+1)-1) = p-1
      enddo
      pb(n_cg(nproc):) = nproc-1

      do p = 0, nproc-1
         if (size(patch%pse(p)%sel(:,:,:), dim=1) /= count(pb(:) == p)) call die("[initproblem:stamp_cg] size(pse(p)%sel(:,:,:), dim=1) /= count(pb(:) == p)")
      enddo

      se(:,:) = 0
      do b = 1, tot_bl
         ! pb(b) = min((b*nproc)/tot_bl, nproc-1)
         call simple_ordering(b, n_bl(:), b_loc(:)) !> \todo implement Morton and Hilbert ordering

         where (dom%has_dir(:))
            se(:, LO) = b_loc(:) * bsize(xdim:zdim)
            se(:, HI) = se(:, LO) + bsize(xdim:zdim) - 1
         endwhere
         patch%pse(pb(b))%sel(b-n_cg(pb(b)+1)+1,:,:) = se(:,:) ! equivalent to call set_pse_sel(pb(b), b-n_cg(pb(b)+1)+1, se)
      enddo
      !> \todo implement merging to larger cuboids

      if (master) call warn("[initproblem:stamp_cg] Experimental implementation")

      deallocate(pb)

   end subroutine stamp_cg

!> \brief Simple block ordering

   subroutine simple_ordering(b, n_bl, b_loc)

      use constants, only: xdim, ydim, zdim

      implicit none

      integer, intent(in) :: b                            !< block number
      integer, dimension(xdim:zdim), intent(in)  :: n_bl  !< extents of block array
      integer, dimension(xdim:zdim), intent(out) :: b_loc !< location of block b

      b_loc(:) = [ mod(b-1, n_bl(xdim)), mod((b-1)/n_bl(xdim), n_bl(ydim)), (b-1)/product(n_bl(xdim:ydim)) ]

   end subroutine simple_ordering

!>
!! \brief Prevent domain decompositions into pieces that are narrower than number of guardcells
!!
!! \details When a piece of grid is narrower than number of guardcells we may expect the following problems:
!!  * Complicated boundary exchange routines because a single neighbour cannot provide valid boundary data in one step.
!!  * Huge memory overhead because number of guardcells is much larger than number of active cells.
!!  * Huge performance penalty because everything becomes dominated by guardcell operations.
!! If this routine prevents you running a simulation then probably you either try to use too many processors or you use wrong domain decomposition scheme.
!<

   logical function is_not_too_small(patch, label) result(patch_divided)

      use constants,  only: I_ONE, LO, HI
      use dataio_pub, only: warn, msg
      use domain,     only: dom, minsize
      use mpi,        only: MPI_IN_PLACE, MPI_LOGICAL, MPI_LAND
      use mpisetup,   only: proc, comm, mpi_err

      implicit none

      type(box_T), intent(inout) :: patch
      character(len=*), intent(in) :: label

      integer :: p

      patch_divided = .true.

      if (allocated(patch%pse)) then
         if (allocated(patch%pse(proc)%sel)) then
            do p = lbound(patch%pse(proc)%sel(:, :, :), dim=1), ubound(patch%pse(proc)%sel(:, :, :), dim=1)
               patch_divided = patch_divided .and. all(patch%pse(proc)%sel(p, :, HI) - patch%pse(proc)%sel(p, :, LO) >= minsize(:) - 1 .or. .not. dom%has_dir(:))
               if (.not. all(patch%pse(proc)%sel(p, :, HI) - patch%pse(proc)%sel(p, :, LO) >= minsize(:) - 1 .or. .not. dom%has_dir(:))) then
                  write(msg,'(3a,2(3i6,a))')"[decomposition:is_not_too_small] ",label," [",patch%pse(proc)%sel(p, :, LO),"]:[",patch%pse(proc)%sel(p, :, HI), "]"
                  call warn(msg)
               endif
            enddo
         else
            patch_divided = .false.
            write(msg,'(3a)')"[decomposition:is_not_too_small] ",label," no pse(proc)%sel"
            call warn(msg)
         endif
      else
         patch_divided = .false.
         write(msg,'(3a)')"[decomposition:is_not_too_small] ",label," no pse"
         call warn(msg)
      endif
      call MPI_Allreduce(MPI_IN_PLACE, patch_divided, I_ONE, MPI_LOGICAL, MPI_LAND, comm, mpi_err)

      if (allocated(patch%pse) .and. .not. patch_divided) call patch%deallocate_pse

   end function is_not_too_small

!> \brief public routine for setting user decompositions

   subroutine set_pse_sel(this, p, n, se)

      use constants,   only: xdim, zdim, LO, HI

      implicit none

      class(box_T),                         intent(inout) :: this
      integer(kind=4),                      intent(in)    :: p            !< process
      integer,                              intent(in)    :: n            !< block number
      integer, dimension(xdim:zdim, LO:HI), intent(in)    :: se           !< segment

      this%pse(p)%sel(n,:,:) = se(:,:)

   end subroutine set_pse_sel

!>
!! \brief allocate the segment list
!!
!! \details Allocate one cuboid spec per process by default or the amount passed in n_cg argument
!<

   subroutine allocate_pse(patch, n_cg)

      use constants,   only: xdim, zdim, LO, HI, I_ONE
      use dataio_pub,  only: die
      use mpisetup,    only: FIRST, LAST, nproc

      implicit none

      class(box_T),                        intent(inout) :: patch
      integer, dimension(nproc), optional, intent(in)    :: n_cg        !< how many segments per process?

      integer :: p

      if (allocated(patch%pse)) call die("[decomposition:allocate_pse] pse already allocated")
      allocate(patch%pse(FIRST:LAST))
      do p = FIRST, LAST
         if (present(n_cg)) then
            !if (n_cg(p+1) > 0)
            allocate(patch%pse(p)%sel(n_cg(p+1), xdim:zdim, LO:HI))
         else
            allocate(patch%pse(p)%sel(I_ONE, xdim:zdim, LO:HI))
         endif
         patch%pse(p)%sel(:, :, :) = 0
      enddo

   end subroutine allocate_pse

!> \brief deallocate the segment list

   subroutine deallocate_pse(patch)

      use mpisetup, only: FIRST, LAST

      implicit none

      class(box_T), intent(inout) :: patch

      integer :: p

      if (allocated(patch%pse)) then
         do p = FIRST, LAST
            deallocate(patch%pse(p)%sel)
         enddo
         deallocate(patch%pse)
      endif

   end subroutine deallocate_pse

end module decomposition
