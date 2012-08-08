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
!! \brief Definition of the type that contains cartesian MPI communicator and related stuff
!!
!! \warning This approach (Cartesian communicator) can not be used on AMR grids and its use in deep multigrid levels is very limited and troublesome.
!! The reason we still keep this quite old feature are pieces of code like SHEAR or corner boundary conditions, which aren't implemented on more general communication yet.
!! It may happen in the fututr that this whole module will disappear.
!<
module cart_comm

   use constants, only: ndims, LO, HI

   implicit none

   private
   public :: cdd

   type cart_decomposition
      integer(kind=4)                          :: comm3d  !< cartesian communicator
      integer(kind=4), dimension(ndims)        :: psize   !< number of divisions in each direction
      integer(kind=4), dimension(ndims)        :: pcoords !< own process coordinates within psize(:)-shaped array of processes
      integer(kind=4), dimension(ndims, LO:HI) :: procn   !< array of neighbours proc numbers
      integer(kind=4)                          :: procxyl !< neighbour in corner boundaries
      integer(kind=4)                          :: procyxl !< neighbour in corner boundaries
    contains
       procedure :: init      !< initialize to invalid values
       procedure :: cleanup   !< free the resources
       procedure :: init_cart !< initialize Cartesian communicator
   end type cart_decomposition

   type(cart_decomposition) :: cdd !< Cartesian Domain Decomposition stuff for global decompositions

contains

   !> \brief Initialize to invalid values

   subroutine init(this)

      use constants, only: I_ONE
      use mpi,       only: MPI_COMM_NULL, MPI_PROC_NULL

      implicit none

      class(cart_decomposition), intent(inout) :: this

      this%procn(:,:) = MPI_PROC_NULL
      this%psize(:)   = -I_ONE
      this%pcoords(:) = -I_ONE
      this%procxyl    = MPI_PROC_NULL
      this%procyxl    = MPI_PROC_NULL
      this%comm3d     = MPI_COMM_NULL

   end subroutine init

   !> \brief Free the resources

   subroutine cleanup(this)

      use mpi,      only: MPI_COMM_NULL
      use mpisetup, only: mpi_err

      implicit none

      class(cart_decomposition), intent(inout) :: this

      if (this%comm3d /= MPI_COMM_NULL) call MPI_Comm_free(this%comm3d, mpi_err)

   end subroutine cleanup

!> \brief Initialize Cartesian communicator

   subroutine init_cart(this, p_size)

      use constants,  only: xdim, ydim, zdim, ndims, LO, HI, BND_COR, I_ONE
      use dataio_pub, only: die
      use domain,     only: dom, is_mpi_noncart, is_refined, reorder
      use mpisetup,   only: comm, proc, mpi_err

      implicit none

      class(cart_decomposition),         intent(inout) :: this
      integer(kind=4), dimension(ndims), intent(in)    :: p_size !< how many cuts are there in each direction?

      integer(kind=4) :: p
      integer(kind=4), dimension(ndims) :: pc

      this%psize(:) = p_size(:)
      if (is_mpi_noncart .or. is_refined) call die("[decomposition:cartesian_tiling] MPI_Cart_create cannot be used for non-rectilinear or AMR domains")

      call MPI_Cart_create(comm, ndims, this%psize, dom%periodic, reorder, this%comm3d, mpi_err)
      call MPI_Cart_coords(this%comm3d, proc, ndims, this%pcoords, mpi_err)

      ! Compute neighbors
      do p = xdim, zdim
         call MPI_Cart_shift(this%comm3d, p-xdim, I_ONE, this%procn(p, LO), this%procn(p, HI), mpi_err)
      enddo

      if (any(dom%bnd(xdim:ydim, LO) == BND_COR)) then
         if (this%pcoords(xdim) == 0 .and. this%pcoords(ydim) > 0) then
            pc(:) = [ this%pcoords(ydim), this%pcoords(xdim), this%pcoords(zdim) ]
            call MPI_Cart_rank(this%comm3d, pc, this%procxyl, mpi_err)
         endif
         if (this%pcoords(ydim) == 0 .and. this%pcoords(xdim) > 0 ) then
            pc(:) = [ this%pcoords(ydim), this%pcoords(xdim), this%pcoords(zdim) ]
            call MPI_Cart_rank(this%comm3d, pc, this%procyxl, mpi_err)
         endif
      endif

      if (any(dom%bnd(xdim:ydim, HI) == BND_COR)) call die("[decomposition:cartesian_tiling] Corner boundary on the right side not implemented anywhere")

   end subroutine init_cart

end module cart_comm
