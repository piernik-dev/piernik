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

   !> \brief A list of grid pieces (typically used as a list of all grids residing on a given process)
   type :: cuboids
      type(cuboid), allocatable, dimension(:) :: c !< an array of grid piece
   end type cuboids

   !> \brief Depiction of global Topology of a level. Use with care, because this is an antiparallel thing
   type :: dot_T
      type(cuboids), dimension(:), allocatable :: gse !< lists of grid chunks on each process (FIRST:LAST)
   contains
      procedure :: cleanup        !< Deallocate everything
      procedure :: update_global  !< Gather updated information about the level and overwrite it to this%dot%gse
   end type dot_T

contains

!> \brief Deallocate everything

   subroutine cleanup(this)

      implicit none

      class(dot_T), intent(inout) :: this

      if (allocated(this%gse)) deallocate(this%gse) ! this%gse(:)%c should be deallocated automagically

   end subroutine cleanup

!>
!! \brief Gather information on cg's currently present on local level, and write new this%dot%gse array
!!
!! OPT: For strict SFC distribution it is possible to determine complete list of neighbors (on the same level and
!! also one level up and down) and exchange only that data. It might be a bit faster for massively parallel runs.
!<

   subroutine update_global(this, first_cgl, cnt)

      use cg_list,    only: cg_list_element
      use constants,  only: I_ZERO, I_ONE, ndims, LO, HI
      use dataio_pub, only: die
      use mpi,        only: MPI_IN_PLACE, MPI_DATATYPE_NULL, MPI_INTEGER
      use mpisetup,   only: FIRST, LAST, proc, comm, mpi_err

      implicit none

      class(dot_T),                   intent(inout) :: this       !< object invoking type bound procedure
      type(cg_list_element), pointer, intent(in)    :: first_cgl  !< first grid on the list belonging to given level
      integer(kind=4),                intent(in)    :: cnt        !< number of grids on given level

      integer(kind=4), dimension(FIRST:LAST) :: allcnt, alloff, ncub_allcnt, ncub_alloff
      integer(kind=4), allocatable, dimension(:) :: allse
      type(cg_list_element), pointer :: cgl
      integer :: i, p
      integer, parameter :: ncub = ndims*HI ! the number of integers in each cuboid

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
         do i = alloff(p), alloff(p) + allcnt(p) - 1
            this%gse(p)%c(i-alloff(p)+1)%se(:, LO) = allse(ncub*i      +1:ncub*i+   ndims) ! we do it in low-level way here again.
            this%gse(p)%c(i-alloff(p)+1)%se(:, HI) = allse(ncub*i+ndims+1:ncub*i+HI*ndims)
         enddo
      enddo

   end subroutine update_global

end module dot
