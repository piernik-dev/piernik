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

!> \brief MPI wrappers - helpers, parameters and variables

module mpi_wrappers

   implicit none

   private
   public :: mpi_type, MPI_wrapper_stats, &
        &    max_rank, row_descr, T_BOO, T_STR, T_I4, T_I8, T_R4, T_R8, T_LAST, col_descr
#ifndef NO_F2018
   public :: first_element, what_type
#endif /* !NO_F2018 */

   logical, save :: MPI_wrapper_stats = .false.  !< collect usage statistics in piernik_MPI_* wrappers

   integer, parameter :: descr_len = 128
   integer, parameter :: max_rank = 4
   character(len=descr_len), parameter :: row_descr = "Rows from scalars to rank-4."

   enum, bind(C)
      enumerator :: T_BOO = 1, T_STR, T_I4, T_I8, T_R4, T_R8
   end enum
   integer, parameter :: T_LAST = T_R8
   character(len=descr_len), parameter :: col_descr = "Columns: logical, character, int32, int64, real32, real64."

contains

#ifndef NO_F2018
#if defined(__INTEL_COMPILER)
#  warning If Intel oneAPI cannot compile this module, try to add -dNO_F2018 to the setup call
#endif /* __INTEL_COMPILER */
!>
!! \brief Return a pointer to the first element of the array
!!
!! Requires Fortran 2018 for SELECT RANK and DIMENSION(..)
!<

   function first_element(buf)

      use dataio_pub, only: die

      implicit none

      class(*), dimension(..), target, intent(in) :: buf !< buffer that will be examined

      class(*), pointer :: first_element
      logical, target :: dummy

      select rank (buf)
         rank (0)
            first_element => buf
         rank (1)
            first_element => buf(lbound(buf, 1))
         rank (2)
            first_element => buf(lbound(buf, 1), lbound(buf, 2))
         rank (3)
            first_element => buf(lbound(buf, 1), lbound(buf, 2), lbound(buf, 3))
         rank (max_rank)
            first_element => buf(lbound(buf, 1), lbound(buf, 2), lbound(buf, 3), lbound(buf, 4))
         rank default
            call die("[mpi_wrappers:first_element] non-implemented rank")
            first_element => dummy  ! suppress compiler warning
      end select
      ! In Fortran 2023 it should be possible to reduce the above construct to just
      !     first_element => buf(@lbound(buf))

   end function first_element

!>
!! \brief Return datatype
!!
!! Requires Fortran 2018 for DIMENSION(..)
!<

   function what_type(buf)

      use dataio_pub, only: die

      implicit none

      class(*), dimension(..), target, intent(in) :: buf !< buffer that will be examined

      integer :: what_type
      class(*), pointer :: pbuf

      pbuf => first_element(buf)
      select type (pbuf)
         type is (logical)
            what_type = T_BOO
         type is (character(len=*))
            what_type = T_STR
         type is (integer(kind=4))
            what_type = T_I4
         type is (integer(kind=8))
            what_type = T_I8
         type is (real(kind=4))
            what_type = T_R4
         type is (real(kind=8))
            what_type = T_R8
         class default
            call die("[mpi_wrappers:what_type] non-implemented type")
            what_type = huge(1)
      end select

   end function what_type
#endif /* !NO_F2018 */

!> \brief convert what_type output to MPI types

   function mpi_type(it)

      use dataio_pub, only: die
      use MPIF,       only: MPI_LOGICAL, MPI_CHARACTER, MPI_INTEGER, MPI_INTEGER8, MPI_REAL, MPI_DOUBLE_PRECISION
#ifdef MPIF08
      use MPIF,       only: MPI_Datatype
#endif /* MPIF08 */

      implicit none

      integer, intent(in) :: it  !< input type

#ifdef MPIF08
      type(MPI_Datatype) :: mpi_type
#else /* !MPIF08 */
      integer(kind=4) :: mpi_type
#endif /* !MPIF08 */

      select case (it)
         case (T_BOO)
            mpi_type = MPI_LOGICAL
         case (T_STR)
            mpi_type = MPI_CHARACTER
         case (T_I4)
            mpi_type = MPI_INTEGER
         case (T_I8)
            mpi_type = MPI_INTEGER8
         case (T_R4)
            mpi_type = MPI_REAL
         case (T_R8)
            mpi_type = MPI_DOUBLE_PRECISION
         case default
            call die("[mpi_wrappers:mpi_type] non-implemented type")
#ifndef MPIF08
            mpi_type = -huge(1)
#endif /* !MPIF08 */
      end select

   end function mpi_type

end module mpi_wrappers
