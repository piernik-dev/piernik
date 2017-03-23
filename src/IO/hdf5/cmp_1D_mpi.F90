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
!! \brief Module that contains routines for comparing 1D arrays across processes.
!!
!! \todo Move it somewhere else if it becomes useful outside reading restart routine
!<

module cmp_1D_mpi
! pulled by HDF5

   implicit none

   private
   public :: compare_array1D

!> \brief Compare a 1D array or string across all processed. Die if any difference is found

   interface compare_array1D
      module procedure compare_int_array1D
      module procedure compare_real_array1D
      module procedure compare_string
   end interface

contains

!>
!! \brief Compare a 1D real array across all processes. Die if any difference is found
!!
!! \todo Add an optional parameter (string), which would identify which comparison failed
!<

   subroutine compare_real_array1D(arr)

      use constants,  only: I_ONE
      use dataio_pub, only: die
      use func,       only: operator(.notequals.)
      use mpi,        only: MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE
      use mpisetup,   only: slave, LAST, comm, proc, mpi_err

      implicit none

      real, dimension(:), intent(in)  :: arr

      real, dimension(size(arr)) :: aux
      integer(kind=4), parameter :: tag = 10

      if (proc /= LAST) call MPI_Send(arr(:), size(arr(:), kind=4), MPI_DOUBLE_PRECISION, proc+I_ONE, tag, comm, mpi_err)
      if (slave) then
         call MPI_Recv(aux(:), size(aux(:), kind=4), MPI_DOUBLE_PRECISION, proc-I_ONE, tag, comm, MPI_STATUS_IGNORE, mpi_err)
         if (any(aux(:).notequals.arr(:))) call die("[cmp_1D_mpi:compare_real_array1D] Inconsistency found.")
      endif

   end subroutine compare_real_array1D

!> \brief Compare a 1D integer array across all processes. Die if any difference is found

   subroutine compare_int_array1D(arr)

      use constants,  only: I_ONE
      use dataio_pub, only: die
      use mpi,        only: MPI_INTEGER, MPI_STATUS_IGNORE
      use mpisetup,   only: slave, LAST, comm, proc, mpi_err

      implicit none

      integer(kind=4), dimension(:), intent(in)  :: arr

      integer(kind=4), dimension(size(arr)) :: aux
      integer(kind=4), parameter            :: tag = 10

      if (proc /= LAST) call MPI_Send(arr(:), size(arr(:), kind=4), MPI_INTEGER, proc+I_ONE, tag, comm, mpi_err)
      if (slave) then
         call MPI_Recv(aux(:), size(aux(:), kind=4), MPI_INTEGER, proc-I_ONE, tag, comm, MPI_STATUS_IGNORE, mpi_err)
         if (any(aux(:) /= arr(:))) call die("[cmp_1D_mpi:compare_int_array1D] Inconsistency found.")
      endif

   end subroutine compare_int_array1D

!> \brief Compare a 1D integer array across all processes. Die if any difference is found

   subroutine compare_string(str)

      use constants,  only: I_ONE
      use dataio_pub, only: die
      use mpi,        only: MPI_CHARACTER, MPI_STATUS_IGNORE
      use mpisetup,   only: slave, LAST, comm, proc, mpi_err

      implicit none

      character(len=*), intent(in) :: str

      character(len=len(str))      :: aux
      integer(kind=4), parameter   :: tag = 10

      if (proc /= LAST) call MPI_Send(str, len(str, kind=4), MPI_CHARACTER, proc+I_ONE, tag, comm, mpi_err)
      if (slave) then
         call MPI_Recv(aux, len(aux, kind=4), MPI_CHARACTER, proc-I_ONE, tag, comm, MPI_STATUS_IGNORE, mpi_err)
         if (aux /= str) call die("[cmp_1D_mpi:compare_string] Inconsistency found.")
      endif

   end subroutine compare_string

end module cmp_1D_mpi
