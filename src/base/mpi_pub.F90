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
!! \brief (KK)
!!
!! This modules is a middleware between Piernik and host mpi implementation
!! IMPORTANT: other modules beside this one should *not* use host mpi module
!! directly
!! REMARK: it required >=ifort-12.1 or >=gcc-FUTURE
!<
module mpi_pub
! pulled by NONE
   implicit none

   private
   public :: mpi_container

   type :: mpi_container
      integer(kind=4) :: float_type
      integer(kind=4) :: logical_type
      integer(kind=4) :: int_type
      integer(kind=4) :: char_type
      integer(kind=4) :: err
      integer(kind=4) :: comm_world

!      contains
!        procedure :: init => init_mpi_container
   end type mpi_container

contains

   subroutine mpi_container(this)
      use mpi, only: MPI_COMM_WORLD, MPI_LOGICAL, MPI_CHARACTER

      implicit none

      class(mpi_container), intent(out) :: this !< object invoking type-bound procedure
      integer, target :: dummy_int
      real, target :: dummy_real
      class(*), pointer :: p

      this%comm_world = MPI_COMM_WORLD
      this%err = 0
      this%logical_type = MPI_LOGICAL
      this%char_type = MPI_CHARACTER

      p => dummy_int; this%int_type = get_mpi_type(p)
      p => dummy_real; this%float_type = get_mpi_type(p)

      contains

         function get_mpi_type(dummy) result (mpi_type)
            use mpi, only: MPI_REAL, MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_INTEGER8

            implicit none
            class(*), pointer :: dummy

            select type (dummy)
               type is (real(kind=4))
                  mpi_type = MPI_REAL
               type is (real(kind=8))
                  mpi_type = MPI_DOUBLE_PRECISION
               type is (integer(kind=4))
                  mpi_type = MPI_INTEGER
               type is (integer(kind=8))
                  mpi_type = MPI_INTEGER8
            end select
         end function get_mpi_type

   end subroutine mpi_container

end module mpi_pub
