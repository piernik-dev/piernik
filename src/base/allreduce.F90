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

   use constants,  only: pSUM, pLAND
   use dataio_pub, only: die  ! QA_WARN this is needed only when NO_F2018 is not set, which is determined at compile time, so Makefile should not depend on it.
#ifdef MPIF08
   use MPIF,       only: MPI_Op
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

#ifdef NO_F2018
   ! We keep that spaghetti to not break compatibility with systems where only
   ! the old mpi interface or compilers are available.
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
#endif /* NO_F2018 */

contains

!> \brief Initialize mpiop translation array as some MPI implementations create MPI_SUM, MPI_MIN, etc. as non-constants

   subroutine init_allreduce

      use MPIF, only: MPI_SUM, MPI_MIN, MPI_MAX, MPI_LOR, MPI_LAND

      implicit none

      mpiop = [ MPI_SUM, MPI_MIN, MPI_MAX, MPI_LOR, MPI_LAND ]  ! this must match the ordering in constants enum

   end subroutine init_allreduce

#ifndef NO_F2018
!>
!! \brief Polymorphic wrapper for MPI_Allreduce with MPI_IN_PLACE
!!
!! Unlimited polymorphism here still requires some spaghetti code but much less than in non-f08 way
!!
!! Requires Fortran 2018 for SELECT RANK and DIMENSION(..)
!<

   subroutine piernik_MPI_Allreduce(var, reduction)

      use dataio_pub, only: die
      use MPIF,       only: MPI_LOGICAL, MPI_INTEGER, MPI_INTEGER8, MPI_REAL, MPI_DOUBLE_PRECISION, &
           &                MPI_IN_PLACE, MPI_COMM_WORLD
      use MPIFUN,     only: MPI_Allreduce
#ifdef MPIF08
      use MPIF,       only: MPI_Datatype
#endif /* MPIF08 */

      implicit none

      class(*), dimension(..), target, intent(inout) :: var       !< variable that will be reduced
      integer(kind=4),                 intent(in)    :: reduction !< integer to mark a reduction type

#ifdef MPIF08
      type(MPI_Datatype) :: dtype
#else /* !MPIF08 */
      integer(kind=4) :: dtype
#endif /* !MPIF08 */

      class(*), pointer :: pvar
      logical, target :: dummy

      ! duplicated code: see bcast
      select rank (var)
         rank (0)
            pvar => var
         rank (1)
            pvar => var(lbound(var, 1))
         rank (2)  ! most likely never used
            pvar => var(lbound(var, 1), lbound(var, 2))
         rank (3)  ! most likely never used
            pvar => var(lbound(var, 1), lbound(var, 2), lbound(var, 3))
         rank (4)  ! most likely never used
            pvar => var(lbound(var, 1), lbound(var, 2), lbound(var, 3), lbound(var, 4))
         rank default
            call die("[allreduce:piernik_MPI_Allreduce] non-implemented rank")
            pvar => dummy  ! suppress compiler warning
      end select
      ! In Fortran 2023 it should be possible to reduce the above construct to just
      !     pvar => var(@lbound(var))

      select type (pvar)
         type is (logical)
            dtype = MPI_LOGICAL
         type is (integer(kind=4))
            dtype = MPI_INTEGER
         type is (integer(kind=8))  ! most likely never used
            dtype = MPI_INTEGER8
         type is (real(kind=4))  ! most likely never used
            dtype = MPI_REAL
         type is (real(kind=8))
            dtype = MPI_DOUBLE_PRECISION
         class default
            call die("[allreduce:piernik_MPI_Allreduce] non-implemented type")
      end select

      call MPI_Allreduce(MPI_IN_PLACE, var, size(var, kind=4), dtype, mpiop(reduction), MPI_COMM_WORLD, err_mpi)

   end subroutine piernik_MPI_Allreduce

#else /* !NO_F2018 */

!> \brief Wrapper for MPI_Allreduce with MPI_IN_PLACE for logical

   subroutine MPI_Allreduce_single_logical(lvar, reduction)

      use constants, only: I_ONE
      use MPIF,      only: MPI_LOGICAL, MPI_IN_PLACE, MPI_COMM_WORLD
      use MPIFUN,    only: MPI_Allreduce

      implicit none

      logical,         intent(inout) :: lvar      !< logical that will be reduced
      integer(kind=4), intent(in)    :: reduction !< integer to mark a reduction type

      call MPI_Allreduce(MPI_IN_PLACE, lvar, I_ONE, MPI_LOGICAL, mpiop(reduction), MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Allreduce_single_logical

!> \brief Wrapper for MPI_Allreduce with MPI_IN_PLACE for integer(kind=4)

   subroutine MPI_Allreduce_single_int4(ivar4, reduction)

      use constants, only: I_ONE
      use MPIF,      only: MPI_INTEGER, MPI_IN_PLACE, MPI_COMM_WORLD
      use MPIFUN,    only: MPI_Allreduce

      implicit none

      integer(kind=4), intent(inout) :: ivar4     !< int4 that will be reduced
      integer(kind=4), intent(in)    :: reduction !< integer to mark a reduction type

      call MPI_Allreduce(MPI_IN_PLACE, ivar4, I_ONE, MPI_INTEGER, mpiop(reduction), MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Allreduce_single_int4

!> \brief Wrapper for MPI_Allreduce with MPI_IN_PLACE for integer(kind=8) most likely never used

   subroutine MPI_Allreduce_single_int8(ivar8, reduction)

      use constants, only: I_ONE
      use MPIF,      only: MPI_INTEGER8, MPI_IN_PLACE, MPI_COMM_WORLD
      use MPIFUN,    only: MPI_Allreduce

      implicit none

      integer(kind=8), intent(inout) :: ivar8     !< int8 that will be reduced
      integer(kind=4), intent(in)    :: reduction !< integer to mark a reduction type

      call MPI_Allreduce(MPI_IN_PLACE, ivar8, I_ONE, MPI_INTEGER8, mpiop(reduction), MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Allreduce_single_int8

!> \brief Wrapper for MPI_Allreduce with MPI_IN_PLACE for integer(kind=4), dimension(:) most likely never used

   subroutine MPI_Allreduce_vec_int4(ivar4, reduction)

      use MPIF,   only: MPI_INTEGER, MPI_IN_PLACE, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Allreduce

      implicit none

      integer(kind=4), dimension(:), intent(inout) :: ivar4     !< int4 that will be reduced
      integer(kind=4),               intent(in)    :: reduction !< integer to mark a reduction type

      call MPI_Allreduce(MPI_IN_PLACE, ivar4, size(ivar4, kind=4), MPI_INTEGER, mpiop(reduction), MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Allreduce_vec_int4

!> \brief Wrapper for MPI_Allreduce with MPI_IN_PLACE for integer(kind=8), dimension(:) most likely never used

   subroutine MPI_Allreduce_vec_int8(ivar8, reduction)

      use MPIF,   only: MPI_INTEGER8, MPI_IN_PLACE, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Allreduce

      implicit none

      integer(kind=8), dimension(:), intent(inout) :: ivar8     !< int8 that will be reduced
      integer(kind=4),               intent(in)    :: reduction !< integer to mark a reduction type

      call MPI_Allreduce(MPI_IN_PLACE, ivar8, size(ivar8, kind=4), MPI_INTEGER8, mpiop(reduction), MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Allreduce_vec_int8

!> \brief Wrapper for MPI_Allreduce with MPI_IN_PLACE for real(kind=4) most likely never used

   subroutine MPI_Allreduce_single_real4(rvar4, reduction)

      use constants, only: I_ONE
      use MPIF,      only: MPI_REAL, MPI_IN_PLACE, MPI_COMM_WORLD
      use MPIFUN,    only: MPI_Allreduce

      implicit none

      real(kind=4),    intent(inout) :: rvar4     !< real4 that will be reduced
      integer(kind=4), intent(in)    :: reduction !< integer to mark a reduction type

      call MPI_Allreduce(MPI_IN_PLACE, rvar4, I_ONE, MPI_REAL, mpiop(reduction), MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Allreduce_single_real4

!> \brief Wrapper for MPI_Allreduce with MPI_IN_PLACE for real(kind=8)

   subroutine MPI_Allreduce_single_real8(rvar8, reduction)

      use constants, only: I_ONE
      use MPIF,      only: MPI_DOUBLE_PRECISION, MPI_IN_PLACE, MPI_COMM_WORLD
      use MPIFUN,    only: MPI_Allreduce

      implicit none

      real(kind=8),    intent(inout) :: rvar8     !< real8 that will be reduced
      integer(kind=4), intent(in)    :: reduction !< integer to mark a reduction type

      call MPI_Allreduce(MPI_IN_PLACE, rvar8, I_ONE, MPI_DOUBLE_PRECISION, mpiop(reduction), MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Allreduce_single_real8

!> \brief Wrapper for MPI_Allreduce with MPI_IN_PLACE for real(kind=4), dimension(:) most likely never used

   subroutine MPI_Allreduce_vec_real4(rvar4, reduction)

      use MPIF,   only: MPI_REAL, MPI_IN_PLACE, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Allreduce

      implicit none

      real(kind=4), dimension(:), intent(inout) :: rvar4     !< real4 that will be reduced
      integer(kind=4),            intent(in)    :: reduction !< integer to mark a reduction type

      call MPI_Allreduce(MPI_IN_PLACE, rvar4, size(rvar4, kind=4), MPI_REAL, mpiop(reduction), MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Allreduce_vec_real4

!> \brief Wrapper for MPI_Allreduce with MPI_IN_PLACE for real(kind=8), dimension(:)

   subroutine MPI_Allreduce_vec_real8(rvar8, reduction)

      use MPIF,   only: MPI_DOUBLE_PRECISION, MPI_IN_PLACE, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Allreduce

      implicit none

      real(kind=8), dimension(:), intent(inout) :: rvar8     !< real8 that will be reduced
      integer(kind=4),            intent(in)    :: reduction !< integer to mark a reduction type

      call MPI_Allreduce(MPI_IN_PLACE, rvar8, size(rvar8, kind=4), MPI_DOUBLE_PRECISION, mpiop(reduction), MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Allreduce_vec_real8

!> \brief Wrapper for MPI_Allreduce with MPI_IN_PLACE for real(kind=8), dimension(:,:,:) most likely never used

   subroutine MPI_Allreduce_arr3d_real8(rvar8, reduction)

      use MPIF,   only: MPI_DOUBLE_PRECISION, MPI_IN_PLACE, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Allreduce

      implicit none

      real(kind=8), dimension(:,:,:), intent(inout) :: rvar8     !< real8 that will be reduced
      integer(kind=4),                intent(in)    :: reduction !< integer to mark a reduction type

      call MPI_Allreduce(MPI_IN_PLACE, rvar8, size(rvar8, kind=4), MPI_DOUBLE_PRECISION, mpiop(reduction), MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Allreduce_arr3d_real8

!> \brief Wrapper for MPI_Allreduce with MPI_IN_PLACE for real(kind=8), dimension(:,:) most likely never used

   subroutine MPI_Allreduce_arr2d_real8(rvar8, reduction)

      use MPIF,   only: MPI_DOUBLE_PRECISION, MPI_IN_PLACE, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Allreduce

      implicit none

      real(kind=8), dimension(:,:), intent(inout) :: rvar8     !< real8 that will be reduced
      integer(kind=4),              intent(in)    :: reduction !< integer to mark a reduction type

      call MPI_Allreduce(MPI_IN_PLACE, rvar8, size(rvar8, kind=4), MPI_DOUBLE_PRECISION, mpiop(reduction), MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Allreduce_arr2d_real8

!> \brief Wrapper for MPI_Allreduce with MPI_IN_PLACE for real(kind=4), dimension(:,:) most likely never used

   subroutine MPI_Allreduce_arr2d_real4(rvar4, reduction)

      use MPIF,   only: MPI_REAL, MPI_IN_PLACE, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Allreduce

      implicit none

      real(kind=4), dimension(:,:), intent(inout) :: rvar4     !< real4 that will be reduced
      integer(kind=4),              intent(in)    :: reduction !< integer to mark a reduction type

      call MPI_Allreduce(MPI_IN_PLACE, rvar4, size(rvar4, kind=4), MPI_REAL, mpiop(reduction), MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Allreduce_arr2d_real4

#endif /* !NO_F2018 */

end module allreduce
