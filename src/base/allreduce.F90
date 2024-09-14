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

!> \brief MPI wrappers for MPI_Allreduce in place and most simplified.

module allreduce

   use constants, only: pSUM, pLAND
#ifdef MPIF08
   use MPIF,      only: MPI_Op
#endif /* MPIF08 */

   implicit none

   private
   public :: piernik_MPI_Allreduce, init_allreduce

   integer(kind=4) :: err_mpi  !< error status

   !< translation between pSUM:pLAND and MPI_SUM:MPI_LAND
#ifdef MPIF08
   type(MPI_Op),    &
#else /* !MPIF08 */
   integer(kind=4), &
#endif /* !MPIF08 */
                    & dimension(pSUM:pLAND) :: mpiop

   interface piernik_MPI_Allreduce
      module procedure MPI_Allreduce_single_logical
      module procedure MPI_Allreduce_single_real4
      module procedure MPI_Allreduce_single_real8
      module procedure MPI_Allreduce_single_int4
      module procedure MPI_Allreduce_single_int8
      module procedure MPI_Allreduce_vec_real4
      module procedure MPI_Allreduce_vec_real8
      module procedure MPI_Allreduce_vec_int4
      module procedure MPI_Allreduce_vec_int8
      module procedure MPI_Allreduce_arr2d_real4
      module procedure MPI_Allreduce_arr2d_real8
      module procedure MPI_Allreduce_arr3d_real8
   end interface piernik_MPI_Allreduce

contains

!> \brief Initialize mpiop translation array as some MPI implementations create MPI_SUM, MPI_MIN, etc. as non-constants

   subroutine init_allreduce

      use MPIF, only: MPI_SUM, MPI_MIN, MPI_MAX, MPI_LOR, MPI_LAND

      implicit none

      mpiop = [ MPI_SUM, MPI_MIN, MPI_MAX, MPI_LOR, MPI_LAND ]  ! this must match the ordering in constants enum

   end subroutine init_allreduce

!>
!! \brief Wrapper for MPI_Allreduce with MPI_IN_PLACE
!! \todo unlimited polymorphism will obsolete me
!<
   subroutine MPI_Allreduce_single_logical(lvar, reduction)

      use constants, only: I_ONE
      use MPIF,      only: MPI_LOGICAL, MPI_IN_PLACE, MPI_COMM_WORLD
      use MPIFUN,    only: MPI_Allreduce

      implicit none

      logical,         intent(inout) :: lvar      !< logical that will be reduced
      integer(kind=4), intent(in)    :: reduction !< integer to mark a reduction type

      call MPI_Allreduce(MPI_IN_PLACE, lvar, I_ONE, MPI_LOGICAL, mpiop(reduction), MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Allreduce_single_logical

!>
!! \brief Wrapper for MPI_Allreduce with MPI_IN_PLACE
!! \todo unlimited polymorphism will obsolete me
!<
   subroutine MPI_Allreduce_single_int4(ivar4, reduction)

      use constants, only: I_ONE
      use MPIF,      only: MPI_INTEGER, MPI_IN_PLACE, MPI_COMM_WORLD
      use MPIFUN,    only: MPI_Allreduce

      implicit none

      integer(kind=4), intent(inout) :: ivar4     !< int4 that will be reduced
      integer(kind=4), intent(in)    :: reduction !< integer to mark a reduction type

      call MPI_Allreduce(MPI_IN_PLACE, ivar4, I_ONE, MPI_INTEGER, mpiop(reduction), MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Allreduce_single_int4

!>
!! \brief Wrapper for MPI_Allreduce with MPI_IN_PLACE
!! \todo unlimited polymorphism will obsolete me
!<
   subroutine MPI_Allreduce_single_int8(ivar8, reduction)

      use constants, only: I_ONE
      use MPIF,      only: MPI_INTEGER8, MPI_IN_PLACE, MPI_COMM_WORLD
      use MPIFUN,    only: MPI_Allreduce

      implicit none

      integer(kind=8), intent(inout) :: ivar8     !< int8 that will be reduced
      integer(kind=4), intent(in)    :: reduction !< integer to mark a reduction type

      call MPI_Allreduce(MPI_IN_PLACE, ivar8, I_ONE, MPI_INTEGER8, mpiop(reduction), MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Allreduce_single_int8

!>
!! \brief Wrapper for MPI_Allreduce with MPI_IN_PLACE
!! \todo unlimited polymorphism will obsolete me
!<
   subroutine MPI_Allreduce_vec_int4(ivar4, reduction)

      use MPIF,   only: MPI_INTEGER, MPI_IN_PLACE, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Allreduce

      implicit none

      integer(kind=4), dimension(:), intent(inout) :: ivar4     !< int4 that will be reduced
      integer(kind=4),               intent(in)    :: reduction !< integer to mark a reduction type

      call MPI_Allreduce(MPI_IN_PLACE, ivar4, size(ivar4, kind=4), MPI_INTEGER, mpiop(reduction), MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Allreduce_vec_int4

!>
!! \brief Wrapper for MPI_Allreduce with MPI_IN_PLACE
!! \todo unlimited polymorphism will obsolete me
!<
   subroutine MPI_Allreduce_vec_int8(ivar8, reduction)

      use MPIF,   only: MPI_INTEGER8, MPI_IN_PLACE, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Allreduce

      implicit none

      integer(kind=8), dimension(:), intent(inout) :: ivar8     !< int8 that will be reduced
      integer(kind=4),               intent(in)    :: reduction !< integer to mark a reduction type

      call MPI_Allreduce(MPI_IN_PLACE, ivar8, size(ivar8, kind=4), MPI_INTEGER8, mpiop(reduction), MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Allreduce_vec_int8

!>
!! \brief Wrapper for MPI_Allreduce with MPI_IN_PLACE
!! \todo unlimited polymorphism will obsolete me
!<
   subroutine MPI_Allreduce_single_real4(rvar4, reduction)

      use constants, only: I_ONE
      use MPIF,      only: MPI_REAL, MPI_IN_PLACE, MPI_COMM_WORLD
      use MPIFUN,    only: MPI_Allreduce

      implicit none

      real(kind=4),    intent(inout) :: rvar4     !< real4 that will be reduced
      integer(kind=4), intent(in)    :: reduction !< integer to mark a reduction type

      call MPI_Allreduce(MPI_IN_PLACE, rvar4, I_ONE, MPI_REAL, mpiop(reduction), MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Allreduce_single_real4

!>
!! \brief Wrapper for MPI_Allreduce with MPI_IN_PLACE
!! \todo unlimited polymorphism will obsolete me
!<
   subroutine MPI_Allreduce_single_real8(rvar8, reduction)

      use constants, only: I_ONE
      use MPIF,      only: MPI_DOUBLE_PRECISION, MPI_IN_PLACE, MPI_COMM_WORLD
      use MPIFUN,    only: MPI_Allreduce

      implicit none

      real(kind=8),    intent(inout) :: rvar8     !< real8 that will be reduced
      integer(kind=4), intent(in)    :: reduction !< integer to mark a reduction type

      call MPI_Allreduce(MPI_IN_PLACE, rvar8, I_ONE, MPI_DOUBLE_PRECISION, mpiop(reduction), MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Allreduce_single_real8

!>
!! \brief Wrapper for MPI_Allreduce with MPI_IN_PLACE
!! \todo unlimited polymorphism will obsolete me
!<
   subroutine MPI_Allreduce_vec_real4(rvar4, reduction)

      use MPIF,   only: MPI_REAL, MPI_IN_PLACE, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Allreduce

      implicit none

      real(kind=4), dimension(:), intent(inout) :: rvar4     !< real4 that will be reduced
      integer(kind=4),            intent(in)    :: reduction !< integer to mark a reduction type

      call MPI_Allreduce(MPI_IN_PLACE, rvar4, size(rvar4, kind=4), MPI_REAL, mpiop(reduction), MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Allreduce_vec_real4

!>
!! \brief Wrapper for MPI_Allreduce with MPI_IN_PLACE
!! \todo unlimited polymorphism will obsolete me
!<
   subroutine MPI_Allreduce_vec_real8(rvar8, reduction)

      use MPIF,   only: MPI_DOUBLE_PRECISION, MPI_IN_PLACE, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Allreduce

      implicit none

      real(kind=8), dimension(:), intent(inout) :: rvar8     !< real8 that will be reduced
      integer(kind=4),            intent(in)    :: reduction !< integer to mark a reduction type

      call MPI_Allreduce(MPI_IN_PLACE, rvar8, size(rvar8, kind=4), MPI_DOUBLE_PRECISION, mpiop(reduction), MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Allreduce_vec_real8

!>
!! \brief Wrapper for MPI_Allreduce with MPI_IN_PLACE
!! \todo unlimited polymorphism will obsolete me
!<
   subroutine MPI_Allreduce_arr3d_real8(rvar8, reduction)

      use MPIF,   only: MPI_DOUBLE_PRECISION, MPI_IN_PLACE, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Allreduce

      implicit none

      real(kind=8), dimension(:,:,:), intent(inout) :: rvar8     !< real8 that will be reduced
      integer(kind=4),                intent(in)    :: reduction !< integer to mark a reduction type

      call MPI_Allreduce(MPI_IN_PLACE, rvar8, size(rvar8, kind=4), MPI_DOUBLE_PRECISION, mpiop(reduction), MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Allreduce_arr3d_real8

!>
!! \brief Wrapper for MPI_Allreduce with MPI_IN_PLACE
!! \todo unlimited polymorphism will obsolete me
!<
   subroutine MPI_Allreduce_arr2d_real8(rvar8, reduction)

      use MPIF,   only: MPI_DOUBLE_PRECISION, MPI_IN_PLACE, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Allreduce

      implicit none

      real(kind=8), dimension(:,:), intent(inout) :: rvar8     !< real8 that will be reduced
      integer(kind=4),              intent(in)    :: reduction !< integer to mark a reduction type

      call MPI_Allreduce(MPI_IN_PLACE, rvar8, size(rvar8, kind=4), MPI_DOUBLE_PRECISION, mpiop(reduction), MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Allreduce_arr2d_real8

!>
!! \brief Wrapper for MPI_Allreduce with MPI_IN_PLACE
!! \todo unlimited polymorphism will obsolete me
!<
   subroutine MPI_Allreduce_arr2d_real4(rvar4, reduction)

      use MPIF,   only: MPI_REAL, MPI_IN_PLACE, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Allreduce

      implicit none

      real(kind=4), dimension(:,:), intent(inout) :: rvar4     !< real4 that will be reduced
      integer(kind=4),              intent(in)    :: reduction !< integer to mark a reduction type

      call MPI_Allreduce(MPI_IN_PLACE, rvar4, size(rvar4, kind=4), MPI_REAL, mpiop(reduction), MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Allreduce_arr2d_real4

end module allreduce
