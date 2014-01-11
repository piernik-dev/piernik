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
!! \brief Depiction of Topology
!!
!! \details The structure that contains most important information of all blocks on all processes on a given level
!! and information on decomposition as well.
!<

module dot

   use decomposition,   only: cuboid

   implicit none

   private
   public :: dot_T

   !> \brief cuboid with SFC_id
   type, extends(cuboid) :: c_id
      integer(kind=8) :: SFCid
   end type C_id

   !> \brief A list of grid pieces (typically used as a list of all grids residing on a given process)
   type :: cuboids
      type(c_id),  allocatable, dimension(:) :: c      !< an array of grid piece
      logical                                :: sorted !< .true. when this%c%SFC_id was sorted
      !> \todo consider adding some integer translation array that allows use of unsorted this%c arrays
   end type cuboids

   !> \brief Depiction of global Topology of a level. Use with care, because this is an antiparallel thing
   type :: dot_T
      type(cuboids),   dimension(:),   allocatable :: gse           !< lists of grid chunks on each process (FIRST:LAST)
      integer(kind=8), dimension(:,:), allocatable :: SFC_id_range  !< min and max SFC id on processes
      integer                                      :: tot_se        !< global number of grids on the level
      logical                                      :: is_blocky     !< .true. when all grid pieces on this level on all processes have same shape and size
      logical                                      :: is_strict_SFC !< .true. when all grid pieces are distributd along SFC curve
   contains
      procedure :: cleanup              !< Deallocate everything
      procedure :: update_global        !< Gather updated information about the level and overwrite it to this%gse
      procedure :: update_local         !< Copy info on local blocks from list of blocks to this%gse
      procedure :: update_tot_se        !< Count all cg on current level for computing tags in vertical_prep
      procedure :: is_consitent         !< Check local consistency
      procedure :: check_blocky         !< Check if all blocks in the domain have same size and shape
      procedure :: update_SFC_id_range  !< Update SFC_id_range array and set this%is_strict_SFC
      procedure :: find_grid            !< Find process and grid_id using SFC_id
   end type dot_T

contains

!> \brief Deallocate everything

   subroutine cleanup(this)

      implicit none

      class(dot_T), intent(inout) :: this

      if (allocated(this%gse)) deallocate(this%gse) ! this%gse(:)%c should be deallocated automagically
      if (allocated(this%SFC_id_range)) deallocate(this%SFC_id_range)

   end subroutine cleanup

!>
!! \brief Gather information on cg's currently present on local level, and write new this%dot%gse array
!!
!! OPT: For strict SFC distribution it is possible to determine complete list of neighbors (on the same level and
!! also one level up and down) and exchange only that data. It might be a bit faster for massively parallel runs.
!<

   subroutine update_global(this, first_cgl, cnt, off)

      use cg_list,    only: cg_list_element
      use constants,  only: I_ZERO, I_ONE, ndims, LO, HI
      use dataio_pub, only: die
      use mpi,        only: MPI_IN_PLACE, MPI_DATATYPE_NULL, MPI_INTEGER
      use mpisetup,   only: FIRST, LAST, proc, comm, mpi_err
      use ordering,   only: SFC_order

      implicit none

      class(dot_T),                      intent(inout) :: this       !< object invoking type bound procedure
      type(cg_list_element), pointer,    intent(in)    :: first_cgl  !< first grid on the list belonging to given level
      integer(kind=4),                   intent(in)    :: cnt        !< number of grids on given level
      integer(kind=8), dimension(ndims), intent(in)    :: off        !< offset of the level

      integer(kind=4), dimension(FIRST:LAST) :: allcnt, alloff, ncub_allcnt, ncub_alloff
      integer(kind=4), allocatable, dimension(:) :: allse
      type(cg_list_element), pointer :: cgl
      integer :: i, p
      integer, parameter :: ncub = ndims*HI ! the number of integers in each cuboid
      integer(kind=8) :: prev_id, cur_id

      ! get the count of grid pieces on each process
      ! Beware: int(this%cnt, kind=4) is not properly updated after calling this%distribute.
      ! Use size(this%dot%gse(proc)%c) if you want to propagate gse before the grid containers are actually added to the level
      ! OPT: this call can be quite long to complete
      call MPI_Allgather(cnt, I_ONE, MPI_INTEGER, allcnt, I_ONE, MPI_INTEGER, comm, mpi_err)

      ! compute offsets for  a composite table of all grid pieces
      alloff(FIRST) = I_ZERO
      do i = FIRST+I_ONE, LAST
         alloff(i) = alloff(i-1) + allcnt(i-1)
      enddo

      allocate(allse(ncub*sum(allcnt(:))))

      ! Collect definitions of own grid pieces
      allse = 0
      cgl => first_cgl
      do i = alloff(proc), alloff(proc) + allcnt(proc) - 1
         if (.not. associated(cgl)) call die("[dot:update_global]] Run out of cg.")
         allse(ncub*i      +1:ncub*i+   ndims) = int(cgl%cg%my_se(:, LO), kind=4)
         allse(ncub*i+ndims+1:ncub*i+HI*ndims) = int(cgl%cg%my_se(:, HI), kind=4) ! we do it in low-level way here. Is it worth using reshape() or something?
         if (any(cgl%cg%my_se > huge(allse(1)))) call die("[dot:update_global]] Implement 8-byte integers in MPI transactions for such huge refinements")
         cgl => cgl%nxt
      enddo
      if (associated(cgl)) call die("[dot:update_global]] Not all cg were read.")

      ! First use of MPI_Allgatherv in the Piernik Code!
      ncub_allcnt(:) = int(ncub * allcnt(:), kind=4)
      ncub_alloff(:) = int(ncub * alloff(:), kind=4)
      call MPI_Allgatherv(MPI_IN_PLACE, I_ZERO, MPI_DATATYPE_NULL, allse, ncub_allcnt, ncub_alloff, MPI_INTEGER, comm, mpi_err)

      ! Rewrite the gse array, forget about past.
      if (.not. allocated(this%gse)) allocate(this%gse(FIRST:LAST))
      do p = FIRST, LAST
         if (allocated(this%gse(p)%c)) deallocate(this%gse(p)%c)
         allocate(this%gse(p)%c(allcnt(p)))
         this%gse(p)%sorted = .true.
         prev_id = -huge(1_8)
         do i = alloff(p), alloff(p) + allcnt(p) - 1
            this%gse(p)%c(i-alloff(p)+1)%se(:, LO) = allse(ncub*i      +1:ncub*i+   ndims) ! we do it in low-level way here again.
            this%gse(p)%c(i-alloff(p)+1)%se(:, HI) = allse(ncub*i+ndims+1:ncub*i+HI*ndims)
            cur_id = SFC_order(this%gse(p)%c(i-alloff(p)+1)%se(:, LO) - off)
            this%gse(p)%c(i-alloff(p)+1)%SFCid = cur_id
            if (prev_id > cur_id) this%gse(p)%sorted = .false.
            prev_id = cur_id
         enddo
      enddo

   end subroutine update_global

!>
!! \brief Copy info on local blocks from list of blocks to this%gse
!!
!! \details Recreate local gse in case anything was derefined, refresh grid_id and make room for new pieces in the
!! gse array
!<

   subroutine update_local(this, first_cgl, cnt)

      use cg_list,  only: cg_list_element
      use mpisetup, only: FIRST, LAST, proc

      implicit none

      class(dot_T),                   intent(inout) :: this       !< object invoking type bound procedure
      type(cg_list_element), pointer, intent(in)    :: first_cgl  !< first grid on the list belonging to given level
      integer(kind=4),                intent(in)    :: cnt        !< number of grids on given level to be initialized
      integer                        :: i
      type(cg_list_element), pointer :: cgl

      if (.not. allocated(this%gse)) allocate(this%gse(FIRST:LAST))
      if (allocated(this%gse(proc)%c)) deallocate(this%gse(proc)%c)
      allocate(this%gse(proc)%c(cnt))
      this%gse(proc)%sorted = .false.
      i = 0
      cgl => first_cgl
      do while (associated(cgl))
         i = i + 1
         this%gse(proc)%c(i)%se(:,:) = cgl%cg%my_se(:,:)
         cgl%cg%grid_id = i
         cgl => cgl%nxt
      enddo
      do while (i<ubound(this%gse(proc)%c, dim=1))
         i = i + 1
         this%gse(proc)%c(i)%se = -huge(1)
      enddo

   end subroutine update_local

!> \brief Count all cg on current level. Useful for computing tags in vertical_prep

   subroutine update_tot_se(this)

      use mpisetup, only: FIRST, LAST

      implicit none

      class(dot_T), intent(inout) :: this  !< object invoking type bound procedure

      integer :: p

      this%tot_se = 0
      do p = FIRST, LAST
         if (allocated(this%gse)) this%tot_se = this%tot_se + ubound(this%gse(p)%c(:), dim=1)
      enddo

   end subroutine update_tot_se

!> \brief Check local consistency

   subroutine is_consitent(this, first_cgl)

      use cg_list,    only: cg_list_element
      use dataio_pub, only: die
      use mpisetup,   only: proc

      implicit none

      class(dot_T),                   intent(inout) :: this       !< object invoking type bound procedure
      type(cg_list_element), pointer, intent(in)    :: first_cgl  !< first grid on the list belonging to given level

      type(cg_list_element), pointer :: cgl
      logical                        :: found_id
      integer                        :: i

      cgl => first_cgl
      do while (associated(cgl))
         found_id = .false.
         do i = lbound(this%gse(proc)%c, dim=1), ubound(this%gse(proc)%c, dim=1)
            if (all(this%gse(proc)%c(i)%se(:,:) == cgl%cg%my_se(:,:))) then
               if (found_id) call die("[dot:is_consitent] multiple occurrences")
               found_id = .true.
            endif
         enddo
         if (.not. found_id) call die("[dot:is_consitent] no occurrences")
         cgl => cgl%nxt
      enddo

   end subroutine is_consitent

   subroutine check_blocky(this)

      use constants,  only: ndims, LO, HI, pLAND, I_ONE
      use mpi,        only: MPI_INTEGER, MPI_REQUEST_NULL
      use mpisetup,   only: proc, req, status, comm, mpi_err, LAST, inflate_req, slave, piernik_MPI_Allreduce

      implicit none

      class(dot_T), intent(inout) :: this       !< object invoking type bound procedure

      integer(kind=4), dimension(ndims) :: shape, shape1
      integer(kind=4), parameter :: sh_tag = 7
      integer, parameter :: nr = 2
      integer :: i

      call inflate_req(nr)
      this%is_blocky = .true.
      shape = 0
      shape1 = 0

      if (allocated(this%gse(proc)%c)) then
         if (size(this%gse(proc)%c, dim=1) > 0) then
            i = lbound(this%gse(proc)%c, dim=1)
            shape = int(this%gse(proc)%c(i)%se(:, HI) - this%gse(proc)%c(i)%se(:, LO), kind=4)
            do i = lbound(this%gse(proc)%c, dim=1) + 1, ubound(this%gse(proc)%c, dim=1)-1
               if (any((this%gse(proc)%c(i)%se(:, HI) - this%gse(proc)%c(i)%se(:, LO)) /= shape)) &
                    this%is_blocky = .false.
            enddo
         endif
      endif
      req = MPI_REQUEST_NULL
      if (slave)     call MPI_Irecv(shape1, size(shape1), MPI_INTEGER, proc-I_ONE, sh_tag, comm, req(1 ), mpi_err)
      if (proc<LAST) call MPI_Isend(shape,  size(shape),  MPI_INTEGER, proc+I_ONE, sh_tag, comm, req(nr), mpi_err)
      call MPI_Waitall(nr, req(:nr), status(:, :nr), mpi_err)
      if (any(shape /= 0) .and. any(shape1 /= 0)) then
         if (any(shape /= shape1)) this%is_blocky = .false.
      endif
      call piernik_MPI_Allreduce(this%is_blocky, pLAND)

   end subroutine check_blocky

!>
!! \brief Update SFC_id_range array and check if the level is decomposed into processes strictly along currently
!! used space-filling curve. Update this%gse(:)%c(:)%SFCid as well.
!<

   subroutine update_SFC_id_range(this, off)

      use constants,  only: LO, HI, ndims
      use dataio_pub, only: die
      use mpi,        only: MPI_INTEGER8
      use mpisetup,   only: FIRST, LAST, proc, comm, mpi_err
      use ordering,   only: SFC_order

      implicit none

      class(dot_T),                      intent(inout) :: this !< object invoking type bound procedure
      integer(kind=8), dimension(ndims), intent(in)    :: off  !< offset of the level

      integer(kind=8) :: SFC_id
      integer(kind=8), dimension(:), allocatable :: id_buf
      integer :: i

      if (.not. allocated(this%SFC_id_range)) allocate(this%SFC_id_range(FIRST:LAST, LO:HI))
      if (any(lbound(this%SFC_id_range) /= [ FIRST, LO ]) .or. any(ubound(this%SFC_id_range) /= [ LAST, HI ])) &
           call die("[dot:update_SFC_id_range] bogus this%SFC_id_range dimensions")

      this%SFC_id_range(proc, :) = [ huge(1), -huge(1) ]
      if (allocated(this%gse)) then
         if (allocated(this%gse(proc)%c)) then
            if (size(this%gse(proc)%c, dim=1) > 0) then
               do i = lbound(this%gse(proc)%c, dim=1), ubound(this%gse(proc)%c, dim=1)
                  this%gse(proc)%c(i)%SFCid = SFC_order(this%gse(proc)%c(i)%se(:, LO)-off)
                  if (this%SFC_id_range(proc, LO) > this%gse(proc)%c(i)%SFCid) this%SFC_id_range(proc, LO) = this%gse(proc)%c(i)%SFCid
                  if (this%SFC_id_range(proc, HI) < this%gse(proc)%c(i)%SFCid) this%SFC_id_range(proc, HI) = this%gse(proc)%c(i)%SFCid
               enddo
            endif
         endif
      endif

      allocate(id_buf(size(this%SFC_id_range)))
      call MPI_Allgather(this%SFC_id_range(proc, :), HI-LO+1, MPI_INTEGER8, id_buf, HI-LO+1, MPI_INTEGER8, comm, mpi_err)
      this%SFC_id_range(:, LO) = id_buf(1::2)
      this%SFC_id_range(:, HI) = id_buf(2::2)

      deallocate(id_buf)

      SFC_id = -huge(1_8)
      this%is_strict_SFC = .true.
      do i = FIRST, LAST
         if (this%SFC_id_range(i, LO) < huge(1)) then ! skip processes that have no grids
            this%is_strict_SFC = this%is_strict_SFC .and. (SFC_id < this%SFC_id_range(i, LO))
            SFC_id = this%SFC_id_range(i, HI)
         endif
      enddo

   end subroutine update_SFC_id_range

!> \brief Find process and grid_id using SFC_id

   subroutine find_grid(this, SFC_id, p, grid_id)

      use constants, only: INVALID, LO, HI
      use mpisetup,  only: FIRST, LAST

      implicit none

      class(dot_T),    intent(inout) :: this
      integer(kind=8), intent(in)    :: SFC_id
      integer,         intent(out)   :: p
      integer,         intent(out)   :: grid_id

      integer, dimension(:), allocatable :: pp
      integer :: ip, il, iu, j

      ! search for p
      p = INVALID
      if (this%is_strict_SFC) then
         ! binary search
         allocate(pp(1))
         pp = INVALID
         il = lbound(this%SFC_id_range, dim=1)
         iu = ubound(this%SFC_id_range, dim=1)
         do while (iu-il > 1)
            j = (il+iu)/2
            if (this%SFC_id_range(j, LO)  == huge(1)) exit ! stop bisection when some processes don't own anything
            ! OPT: this may require somewhat better solution for massively parallel runs.
            if (this%SFC_id_range(j, LO) <= SFC_id) then
               il = j
            else
               iu = j
            endif
         enddo
         do j = il, iu
            if (this%SFC_id_range(j, LO) <= SFC_id .and. this%SFC_id_range(j, HI) >= SFC_id) pp(1) = j
         enddo
      else
         ! sequential search, return many possible p's
         allocate(pp(FIRST:LAST))
         pp = INVALID
         il = lbound(pp, dim=1)
         do ip = lbound(pp, dim=1), ubound(pp, dim=1)
            if (this%SFC_id_range(ip, LO) <= SFC_id .and. this%SFC_id_range(ip, HI) >= SFC_id) then
               pp(il) = ip
               il = il + 1
            endif
         enddo
      endif

      ! search for grid_id
      grid_id = INVALID
      do ip = lbound(pp, dim=1), ubound(pp, dim=1)
         if (pp(ip) /= INVALID) then
            if (this%gse(pp(ip))%sorted) then
               ! binary search
               il = lbound(this%gse(pp(ip))%c, dim=1)
               iu = ubound(this%gse(pp(ip))%c, dim=1)
               do while (iu-il > 1)
                  j = (il+iu)/2
                  if (this%gse(pp(ip))%c(j)%SFCid <= SFC_id) then
                     il = j
                  else
                     iu = j
                  endif
               enddo
               call seq_search_g(il, iu)
            else
               ! sequential search
               call seq_search_g(lbound(this%gse(pp(ip))%c, dim=1), ubound(this%gse(pp(ip))%c, dim=1))
            endif
         else
            exit
         endif
         if (grid_id /= INVALID) exit
      enddo

      if (grid_id == INVALID) p = INVALID
      deallocate(pp)

   contains

      subroutine seq_search_g(i1, i2)

         implicit none

         integer, intent(in) :: i1, i2

         do j = i1, i2
            if (this%gse(pp(ip))%c(j)%SFCid == SFC_id) then
               grid_id = j
               p = pp(ip)
               exit
            endif
         enddo

      end subroutine seq_search_g

   end subroutine find_grid

end module dot
