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

!> \brief MPI wrappers for MPI_Allreduce - in place and most simplified.

module isend_irecv

   use cnt_array,  only: arrsum

   implicit none

   private
   public :: piernik_Isend, piernik_Irecv, init_sr, cleanup_sr

   integer(kind=4) :: err_mpi  !< error status

   type(arrsum) :: size_sr
   enum, bind(C)
      enumerator :: T_BOO = 1, T_I4, T_I8, T_R4, T_R8
   end enum
   enum, bind(C)
      enumerator :: I_S = 1, I_R
   end enum
   integer, parameter :: max_dtype = T_R8  ! required for stats gathering
   integer, parameter :: max_rank = 4
   integer, parameter :: max_op = I_R  ! separate counters for send and recaive operations

contains

!> \brief Initialize MPI_Isend + MPI_Irecv stat counters

   subroutine init_sr

      implicit none

      call size_sr%init([max_op, max_rank + 1, max_dtype])

   end subroutine init_sr

!> \brief Print log and clean up MPI_Isend + MPI_Irecv stat counter

   subroutine cleanup_sr

      use constants,    only: V_DEBUG
      use mpi_wrappers, only: MPI_wrapper_stats
      use mpisetup,     only: master

      implicit none

      if (master .and. MPI_wrapper_stats) &
           call size_sr%print("Isend/Irecv total elements(calls). Columns: logical, int32, int64, real32, real64. Rows from scalars to rank-4.", V_DEBUG)

      call size_sr%cleanup

   end subroutine cleanup_sr

!>
!! \brief Polymorphic wrapper for MPI_Isend
!!
!! Requires Fortran 2018 for SELECT RANK and DIMENSION(..)
!<

   subroutine piernik_Isend(buf, count, datatype, dest, tag, req)

      use dataio_pub,   only: die
      use mpi_wrappers, only: MPI_wrapper_stats
      use MPIF,         only: MPI_LOGICAL, MPI_INTEGER, MPI_INTEGER8, MPI_REAL, MPI_DOUBLE_PRECISION
      use MPIFUN,       only: MPI_Isend
      use req_array,    only: req_arr
#ifdef MPIF08
      use MPIF,         only: MPI_Datatype, MPI_Request
#endif /* MPIF08 */

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
      type(MPI_Datatype) :: dtype
      type(MPI_Request), pointer :: r
#else /* !MPIF08 */
      integer(kind=4) :: dtype
      integer(kind=4),   pointer :: r
#endif /* !MPIF08 */
      class(*), pointer :: pbuf
      logical, target :: dummy
      integer :: it

      ! duplicated code: see bcast
      select rank (buf)
         rank (0)
            pbuf => buf
         rank (1)
            pbuf => buf(lbound(buf, 1))
         rank (2)  ! most likely never used
            pbuf => buf(lbound(buf, 1), lbound(buf, 2))
         rank (3)  ! most likely never used
            pbuf => buf(lbound(buf, 1), lbound(buf, 2), lbound(buf, 3))
         rank (max_rank)  ! most likely never used
            pbuf => buf(lbound(buf, 1), lbound(buf, 2), lbound(buf, 3), lbound(buf, 4))
         rank default
            call die("[isend_irecv:piernik_Isend] non-implemented rank")
            pbuf => dummy  ! suppress compiler warning
      end select
      ! In Fortran 2023 it should be possible to reduce the above construct to just
      !     pbuf => buf(@lbound(buf))

      select type (pbuf)
         type is (logical)
            dtype = MPI_LOGICAL
            it = T_BOO
         type is (integer(kind=4))
            dtype = MPI_INTEGER
            it = T_I4
         type is (integer(kind=8))  ! most likely never used
            dtype = MPI_INTEGER8
            it = T_I8
         type is (real(kind=4))  ! most likely never used
            dtype = MPI_REAL
            it = T_R4
         type is (real(kind=8))
            dtype = MPI_DOUBLE_PRECISION
            it = T_R8
         class default
            call die("[isend_irecv:piernik_Isend] non-implemented type")
            it = huge(1)
      end select

      if (MPI_wrapper_stats) call size_sr%add([int(I_S), rank(buf)+1, it], size(buf, kind=8))

      ! No idea why MPI_Irecv here doesn't work correctly when req%nxt() is passed directly.
      r => req%nxt()
      call MPI_Isend(buf, count, datatype, dest, tag, req%comm, r, err_mpi)

   end subroutine piernik_Isend

   subroutine piernik_Irecv(buf, count, datatype, source, tag, req)

      use dataio_pub,   only: die
      use mpi_wrappers, only: MPI_wrapper_stats
      use MPIF,         only: MPI_LOGICAL, MPI_INTEGER, MPI_INTEGER8, MPI_REAL, MPI_DOUBLE_PRECISION
      use MPIFUN,       only: MPI_Irecv
      use req_array,    only: req_arr
#ifdef MPIF08
      use MPIF,         only: MPI_Datatype, MPI_Request
#endif /* MPIF08 */

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
      type(MPI_Datatype) :: dtype
      type(MPI_Request), pointer :: r
#else /* !MPIF08 */
      integer(kind=4) :: dtype
      integer(kind=4),   pointer :: r
#endif /* !MPIF08 */
      class(*), pointer :: pbuf
      logical, target :: dummy
      integer :: it

      ! duplicated code: see bcast
      select rank (buf)
         rank (0)
            pbuf => buf
         rank (1)
            pbuf => buf(lbound(buf, 1))
         rank (2)  ! most likely never used
            pbuf => buf(lbound(buf, 1), lbound(buf, 2))
         rank (3)  ! most likely never used
            pbuf => buf(lbound(buf, 1), lbound(buf, 2), lbound(buf, 3))
         rank (max_rank)  ! most likely never used
            pbuf => buf(lbound(buf, 1), lbound(buf, 2), lbound(buf, 3), lbound(buf, 4))
         rank default
            call die("[isend_irecv:piernik_Irecv] non-implemented rank")
            pbuf => dummy  ! suppress compiler warning
      end select
      ! In Fortran 2023 it should be possible to reduce the above construct to just
      !     pbuf => buf(@lbound(buf))

      select type (pbuf)
         type is (logical)
            dtype = MPI_LOGICAL
            it = T_BOO
         type is (integer(kind=4))
            dtype = MPI_INTEGER
            it = T_I4
         type is (integer(kind=8))  ! most likely never used
            dtype = MPI_INTEGER8
            it = T_I8
         type is (real(kind=4))  ! most likely never used
            dtype = MPI_REAL
            it = T_R4
         type is (real(kind=8))
            dtype = MPI_DOUBLE_PRECISION
            it = T_R8
         class default
            call die("[isend_irecv:piernik_Irecv] non-implemented type")
            it = huge(1)
      end select

      if (MPI_wrapper_stats) call size_sr%add([int(I_R), rank(buf)+1, it], size(buf, kind=8))

      ! No idea why MPI_Irecv here doesn't work correctly when req%nxt() is passed directly.
      r => req%nxt()
      call MPI_Irecv(buf, count, datatype, source, tag, req%comm, r, err_mpi)

   end subroutine piernik_Irecv

end module isend_irecv
