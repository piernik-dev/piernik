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
!! \brief MPI wrappers for MPI_Isend and MPI_Irecv.
!!
!! These are supposed to be used only where we expect 1-to-1 correspondence
!! between send and receive, i.e. no tag duplicates are allowed between any pair
!! of processes.
!<

module isend_irecv

   use cnt_array,  only: arrsum

   implicit none

   private
   public :: piernik_Isend, piernik_Irecv, init_sr, cleanup_sr

   integer(kind=4) :: err_mpi  !< error status
   type(arrsum) :: size_sr
   enum, bind(C)
      enumerator :: I_S = 1, I_R
   end enum
   integer, parameter :: max_op = I_R  ! separate counters for send and recaive operations

contains

!> \brief Initialize MPI_Isend + MPI_Irecv stat counters

   subroutine init_sr

      use mpi_wrappers, only: max_rank, T_LAST

      implicit none

      call size_sr%init([max_op, max_rank + 1, T_LAST])

   end subroutine init_sr

!> \brief Print log and clean up MPI_Isend + MPI_Irecv stat counter

   subroutine cleanup_sr

      use constants,    only: V_DEBUG
      use mpi_wrappers, only: MPI_wrapper_stats, row_descr, col_descr
      use mpisetup,     only: master

      implicit none

      if (master .and. MPI_wrapper_stats) &
           call size_sr%print("Isend/Irecv total elements(calls). " // trim(col_descr) // " " // trim(row_descr), V_DEBUG)

      call size_sr%cleanup

   end subroutine cleanup_sr

!>
!! \brief Polymorphic wrapper for MPI_Isend
!!
!! Requires Fortran 2018 for DIMENSION(..)
!<

   subroutine piernik_Isend(buf, count, datatype, dest, tag, req)

      use MPIFUN,       only: MPI_Isend
      use req_array,    only: req_arr
#ifdef MPIF08
      use MPIF,         only: MPI_Datatype, MPI_Request
#endif /* MPIF08 */
#ifndef NO_F2018
      use mpi_wrappers, only: MPI_wrapper_stats, what_type
#endif /* !NO_F2018 */

      implicit none

      class(*), dimension(..), target, intent(inout) :: buf       !< buffer that will be sent
      integer(kind=4),                 intent(in)    :: count     !< number of elements
      integer(kind=4),                 intent(in)    :: dest      !< data target
      integer(kind=4),                 intent(in)    :: tag       !< the tag
#ifdef MPIF08
      type(MPI_Datatype),              intent(in)    :: datatype  !< MPI data type
#else /* !MPIF08 */
      integer(kind=4),                 intent(in)    :: datatype
#endif /* !MPIF08 */
      class(req_arr),                  intent(inout) :: req       !< array for requests

#ifdef MPIF08
      type(MPI_Request), pointer :: r
#else /* !MPIF08 */
      integer(kind=4),   pointer :: r
#endif /* !MPIF08 */

#ifndef NO_F2018
      if (MPI_wrapper_stats) &
           call size_sr%add([int(I_S), rank(buf)+1, what_type(buf)], size(buf, kind=8))
#endif /* !NO_F2018 */

      ! No idea why MPI_Irecv here doesn't work correctly when req%nxt() is passed directly.
      r => req%nxt()
      call MPI_Isend(buf, count, datatype, dest, tag, req%comm, r, err_mpi)
      call req%store_tag(tag, dest, recv = .false.)

   end subroutine piernik_Isend

!>
!! \brief Polymorphic wrapper for MPI_Irecv
!!
!! Requires Fortran 2018 for DIMENSION(..)
!<

   subroutine piernik_Irecv(buf, count, datatype, source, tag, req)

      use MPIFUN,       only: MPI_Irecv
      use req_array,    only: req_arr
#ifdef MPIF08
      use MPIF,         only: MPI_Datatype, MPI_Request
#endif /* MPIF08 */
#ifndef NO_F2018
      use mpi_wrappers, only: MPI_wrapper_stats, what_type
#endif /* !NO_F2018 */

      implicit none

      class(*), dimension(..), target, intent(inout) :: buf       !< buffer that will be received
      integer(kind=4),                 intent(in)    :: count     !< number of elements
      integer(kind=4),                 intent(in)    :: source    !< data origin
      integer(kind=4),                 intent(in)    :: tag       !< the tag
#ifdef MPIF08
      type(MPI_Datatype),              intent(in)    :: datatype  !< MPI data type
#else /* !MPIF08 */
      integer(kind=4),                 intent(in)    :: datatype
#endif /* !MPIF08 */
      class(req_arr),                  intent(inout) :: req       !< array for requests

#ifdef MPIF08
      type(MPI_Request), pointer :: r
#else /* !MPIF08 */
      integer(kind=4),   pointer :: r
#endif /* !MPIF08 */

#ifndef NO_F2018
      if (MPI_wrapper_stats) &
           call size_sr%add([int(I_R), rank(buf)+1, what_type(buf)], size(buf, kind=8))
#endif /* !NO_F2018 */

      ! No idea why MPI_Irecv here doesn't work correctly when req%nxt() is passed directly.
      r => req%nxt()
      call MPI_Irecv(buf, count, datatype, source, tag, req%comm, r, err_mpi)
      call req%store_tag(tag, source, recv = .true.)

   end subroutine piernik_Irecv

end module isend_irecv
