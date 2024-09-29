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

!> \brief The segment type: contaier for various internal boundary exchanges

module grid_cont_bseg

   use constants,    only: xdim, zdim, LO, HI
   use grid_cont_na, only: grid_container_na_t
#ifdef MPIF08
   use MPIF,         only: MPI_Datatype
#endif /* MPIF08 */

   implicit none

   private
   public :: tgt_list, segment, grid_container_bseg_t

   !> \brief Specification of segment of data for boundary exchange
   type :: segment
      integer(kind=4) :: proc                              !< target process
      integer(kind=8), dimension(xdim:zdim, LO:HI) :: se   !< range
      integer(kind=4) :: tag                               !< unique tag for data exchange
      real, allocatable, dimension(:,:,:)   :: buf         !< buffer for the 3D (scalar) data to be sent or received
      real, allocatable, dimension(:,:,:,:) :: buf4        !< buffer for the 4D (vector) data to be sent or received
      integer :: ireq  !< request index (for MPI_Test in fc_fluxes)
#ifdef MPIF08
      type(MPI_Datatype) :: sub_type                       !< MPI type related to this segment
#else /* !MPIF08 */
      integer(kind=4) :: sub_type                          !< MPI type related to this segment
#endif /* !MPIF08 */

      integer(kind=8), dimension(xdim:zdim, LO:HI) :: se2  !< auxiliary range, used in cg_level_connected:vertical_bf_prep
      class(grid_container_bseg_t), pointer :: local       !< set this pointer to non-null when the exchange is local
   contains
      procedure :: send_buf   !< Perform piernik_Isend, when buf, proc and tag are ready
      procedure :: recv_buf   !< Perform piernik_Irecv, when buf, proc and tag are ready
      procedure :: send_buf4  !< Perform piernik_Isend, when buf4, proc and tag are ready
      procedure :: recv_buf4  !< Perform piernik_Irecv, when buf4, proc and tag are ready
   end type segment

   !< \brief target list container for prolongations, restrictions and boundary exchanges
   type :: tgt_list
      type(segment), dimension(:), allocatable :: seg  !< segments of data to be received or sent
   contains
      procedure :: add_seg  !< Add an new segment, reallocate if necessary
      procedure :: cleanup  !< Deallocate internals
   end type tgt_list

   type, extends(grid_container_na_t), abstract :: grid_container_bseg_t
      ! External boundary conditions
      type(tgt_list), dimension(:), allocatable :: i_bnd  !< description of incoming boundary data
      type(tgt_list), dimension(:), allocatable :: o_bnd  !< description of outgoing boundary data
      ! initialization of i_bnd and o_bnd is done in cg_list_neighbors because we don't have access to cg_level%dot here
   contains
      procedure :: cleanup_bseg  !< Deallocate all internals
   end type grid_container_bseg_t

contains

!> \brief Perform piernik_Isend, when buf, proc and tag are ready

   subroutine send_buf(this, req)

      use isend_irecv,  only: piernik_Isend
      use MPIF,         only: MPI_DOUBLE_PRECISION
      use req_array,    only: req_arr

      implicit none

      class(segment), intent(inout) :: this  !< object invoking type-bound procedure
      class(req_arr), intent(inout) :: req   !< array for requests

      ! explicit buf(lbound(buf, ...), ...) needed to prevent valgrind complains on "Invalid read of size 8", at least with gfortran 12.3
      call piernik_Isend(this%buf(lbound(this%buf, 1):, lbound(this%buf, 2):, lbound(this%buf, 3):), size(this%buf, kind=4), MPI_DOUBLE_PRECISION, this%proc, this%tag, req)
      this%ireq = req%n

   end subroutine send_buf

!< \brief Perform piernik_Irecv, when buf, proc and tag are ready

   subroutine recv_buf(this, req)

      use isend_irecv,  only: piernik_Irecv
      use MPIF,         only: MPI_DOUBLE_PRECISION
      use req_array,    only: req_arr

      implicit none

      class(segment), intent(inout) :: this  !< object invoking type-bound procedure
      class(req_arr), intent(inout) :: req   !< array for requests

      call piernik_Irecv(this%buf(lbound(this%buf, 1):, lbound(this%buf, 2):, lbound(this%buf, 3):), size(this%buf, kind=4), MPI_DOUBLE_PRECISION, this%proc, this%tag, req)
      this%ireq = req%n

   end subroutine recv_buf

!> \brief Perform piernik_Isend, when buf4, proc and tag are ready

   subroutine send_buf4(this, req)

      use isend_irecv,  only: piernik_Isend
      use MPIF,         only: MPI_DOUBLE_PRECISION
      use req_array,    only: req_arr

      implicit none

      class(segment), intent(inout) :: this  !< object invoking type-bound procedure
      class(req_arr), intent(inout) :: req   !< array for requests

      call piernik_Isend(this%buf4(lbound(this%buf4, 1):, lbound(this%buf4, 2):, lbound(this%buf4, 3):, lbound(this%buf4, 4):), size(this%buf4, kind=4), MPI_DOUBLE_PRECISION, this%proc, this%tag, req)
      this%ireq = req%n

   end subroutine send_buf4

!< \brief Perform piernik_Irecv, when buf4, proc and tag are ready

   subroutine recv_buf4(this, req)

      use isend_irecv,  only: piernik_Irecv
      use MPIF,         only: MPI_DOUBLE_PRECISION
      use req_array,    only: req_arr

      implicit none

      class(segment), intent(inout) :: this  !< object invoking type-bound procedure
      class(req_arr), intent(inout) :: req   !< array for requests

      call piernik_Irecv(this%buf4(lbound(this%buf4, 1):, lbound(this%buf4, 2):, lbound(this%buf4, 3):, lbound(this%buf4, 4):), size(this%buf4, kind=4), MPI_DOUBLE_PRECISION, this%proc, this%tag, req)
      this%ireq = req%n

   end subroutine recv_buf4

!> \brief Add a new segment, reallocate if necessary

   subroutine add_seg(this, proc, se, tag)

      use dataio_pub, only: die

      implicit none

      class(tgt_list),                              intent(inout) :: this !< object invoking type-bound procedure
      integer(kind=4),                              intent(in)    :: proc !< process to be communicated
      integer(kind=8), dimension(xdim:zdim, LO:HI), intent(in)    :: se   !< segment definition
      integer(kind=4),                              intent(in)    :: tag  !< tag for MPI calls

      type(segment), dimension(:), allocatable :: tmp
      integer :: g

      if (tag < 0) call die("[grid_container_bseg:add_seg] tag<0")

      if (allocated(this%seg)) then
         allocate(tmp(lbound(this%seg, dim=1):ubound(this%seg, dim=1)+1))
         tmp(:ubound(this%seg, dim=1)) = this%seg
         call move_alloc(from=tmp, to=this%seg)
      else
         allocate(this%seg(1))
      endif

      g = ubound(this%seg, dim=1)

      this%seg(g)%proc = proc
      this%seg(g)%se = se
      this%seg(g)%tag = tag
      nullify(this%seg(g)%local)

   end subroutine add_seg

!> \brief Deallocate all internals

   subroutine cleanup(this)

      implicit none

      class(tgt_list), intent(inout) :: this  !< object invoking type-bound procedure

      integer :: g

      if (allocated(this%seg)) then
         do g = lbound(this%seg, dim=1), ubound(this%seg, dim=1)
            if (allocated(this%seg(g)%buf )) deallocate(this%seg(g)%buf )
            if (allocated(this%seg(g)%buf4)) deallocate(this%seg(g)%buf4)
         enddo
         deallocate(this%seg)
      endif

   end subroutine cleanup

!> \brief Deallocate all internals

   subroutine cleanup_bseg(this)

      implicit none

      class(grid_container_bseg_t), intent(inout) :: this  !< object invoking type-bound procedure

      integer :: d

      if (allocated(this%i_bnd)) then
         do d = lbound(this%i_bnd, dim=1), ubound(this%i_bnd, dim=1)
            call this%i_bnd(d)%cleanup
         enddo
         deallocate(this%i_bnd)
      endif

      if (allocated(this%o_bnd)) then
         do d = lbound(this%o_bnd, dim=1), ubound(this%o_bnd, dim=1)
            call this%o_bnd(d)%cleanup
         enddo
         deallocate(this%o_bnd)
      endif

   end subroutine cleanup_bseg

end module grid_cont_bseg
