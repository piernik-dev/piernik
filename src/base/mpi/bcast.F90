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

!> \brief MPI wrappers for MPI_Bcast from FIRST to all, most simplified.

module bcast

   use cnt_array,  only: arrsum
   use dataio_pub, only: die  ! QA_WARN this is needed for correct calculation of dependencies when NO_F2018 is not set. NO_F2018 is determined at compile time.

   implicit none

   private
   public :: piernik_MPI_Bcast, init_bcast, cleanup_bcast

   integer(kind=4) :: err_mpi  !< error status
   type(arrsum) :: size_bcast

#ifdef NO_F2018
   ! We keep that spaghetti to not break compatibility with systems where only
   ! the old mpi interface or compilers are available.
   interface piernik_MPI_Bcast
      module procedure MPI_Bcast_single_logical
      module procedure MPI_Bcast_single_string
      module procedure MPI_Bcast_single_real4
      module procedure MPI_Bcast_single_real8
      module procedure MPI_Bcast_single_int4
      module procedure MPI_Bcast_single_int8
      module procedure MPI_Bcast_vec_logical
      module procedure MPI_Bcast_vec_string
      module procedure MPI_Bcast_vec_int4
      module procedure MPI_Bcast_vec_int8
      module procedure MPI_Bcast_vec_real4
      module procedure MPI_Bcast_vec_real8
      module procedure MPI_Bcast_arr2d_int4
      module procedure MPI_Bcast_arr2d_int8
      module procedure MPI_Bcast_arr2d_real4
      module procedure MPI_Bcast_arr2d_real8
      module procedure MPI_Bcast_arr3d_int4
      module procedure MPI_Bcast_arr3d_int8
      module procedure MPI_Bcast_arr3d_real4
      module procedure MPI_Bcast_arr3d_real8
   end interface piernik_MPI_Bcast
#endif /* NO_F2018 */

contains

!> \brief Initialize MPI_Bcast stat counter

   subroutine init_bcast

      use mpi_wrappers, only: max_rank, T_LAST

      implicit none

      call size_bcast%init([max_rank + 1, T_LAST])

   end subroutine init_bcast

!> \brief Print log and clean up MPI_Bcast stat counter

   subroutine cleanup_bcast

      use constants,    only: V_DEBUG
      use mpi_wrappers, only: MPI_wrapper_stats, row_descr, col_descr
      use mpisetup,     only: master

      implicit none

      if (master .and. MPI_wrapper_stats) &
           call size_bcast%print("Bcast total elements(calls). " // trim(col_descr) // " " // trim(row_descr), V_DEBUG)

      call size_bcast%cleanup

   end subroutine cleanup_bcast

#ifndef NO_F2018

!>
!! \brief Polymorphic wrapper for MPI_Bcast
!!
!! Unlimited polymorphism here still requires some spaghetti code but much less than in non-f08 way
!!
!! Requires Fortran 2018 for DIMENSION(..)
!<

   subroutine piernik_MPI_Bcast(var, clen)

      use dataio_pub,   only: die
      use mpi_wrappers, only: MPI_wrapper_stats, first_element, what_type, mpi_type, T_STR
      use mpisetup,     only: FIRST
      use MPIF,         only: MPI_COMM_WORLD
      use MPIFUN,       only: MPI_Bcast
#ifdef MPIF08
      use MPIF,         only: MPI_Datatype
#endif /* MPIF08 */

      implicit none

      class(*), dimension(..), target, intent(inout) :: var   !< variable that will be broadcasted
      integer(kind=4), optional,       intent(in)    :: clen  !< length of the string (only when type(var) is character)

#ifdef MPIF08
      type(MPI_Datatype) :: dtype
#else /* !MPIF08 */
      integer(kind=4) :: dtype
#endif /* !MPIF08 */
      class(*), pointer :: pvar
      integer(kind=4) :: lenmul
      integer :: it

      pvar => first_element(var)

      it = what_type(pvar)
      dtype = mpi_type(it)
      lenmul = 1
      select type (pvar)
         type is (character(len=*))
            if (.not. present(clen)) call die("[bcast:piernik_MPI_Bcast] clen is required for strings")
            lenmul = clen
      end select

      if (present(clen) .and. (it /= T_STR)) call die("[bcast:piernik_MPI_Bcast] clen makes no sense for non-strings")

      if (MPI_wrapper_stats) call size_bcast%add([rank(var)+1, it], lenmul*size(var, kind=8))

      call MPI_Bcast(var, lenmul*size(var, kind=4), dtype, FIRST, MPI_COMM_WORLD, err_mpi)

   end subroutine piernik_MPI_Bcast

#else /* !NO_F2018 */

!> \brief Wrapper for MPI_Bcast for logical

   subroutine MPI_Bcast_single_logical(lvar)

      use constants, only: I_ONE
      use mpisetup,  only: FIRST
      use MPIF,      only: MPI_LOGICAL, MPI_COMM_WORLD
      use MPIFUN,    only: MPI_Bcast

      implicit none

      logical, intent(inout) :: lvar  !< logical scalar that will be broadcasted

      call MPI_Bcast(lvar, I_ONE, MPI_LOGICAL, FIRST, MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Bcast_single_logical

!> \brief Wrapper for MPI_Bcast for logical, dimension(:)

   subroutine MPI_Bcast_vec_logical(lvar)

      use mpisetup,  only: FIRST
      use MPIF,   only: MPI_LOGICAL, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Bcast

      implicit none

      logical, dimension(:), intent(inout) :: lvar  !< logical scalar that will be broadcasted

      call MPI_Bcast(lvar, size(lvar, kind=4), MPI_LOGICAL, FIRST, MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Bcast_vec_logical

!> \brief Wrapper for MPI_Bcast for character(len=clen)

   subroutine MPI_Bcast_single_string(cvar, clen)

      use mpisetup,  only: FIRST
      use MPIF,   only: MPI_CHARACTER, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Bcast

      implicit none

      character(len=*), intent(inout) :: cvar  !< string that will be broadcasted
      integer(kind=4),  intent(in)    :: clen  !< length of the cvar

      call MPI_Bcast(cvar, clen, MPI_CHARACTER, FIRST, MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Bcast_single_string

!> \brief Wrapper for MPI_Bcast for character(len=clen), dimension(:)

   subroutine MPI_Bcast_vec_string(cvar, clen)

      use mpisetup,  only: FIRST
      use MPIF,   only: MPI_CHARACTER, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Bcast

      implicit none

      character(len=*), dimension(:), intent(inout) :: cvar  !< vector of strings that will be broadcasted
      integer(kind=4),                intent(in)    :: clen  !< length of the cvar

      call MPI_Bcast(cvar, clen*size(cvar, kind=4), MPI_CHARACTER, FIRST, MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Bcast_vec_string

!> \brief Wrapper for MPI_Bcast for integer(kind=4)

   subroutine MPI_Bcast_single_int4(ivar4)

      use constants, only: I_ONE
      use mpisetup,  only: FIRST
      use MPIF,      only: MPI_INTEGER, MPI_COMM_WORLD
      use MPIFUN,    only: MPI_Bcast

      implicit none

      integer(kind=4), intent(inout) :: ivar4  !< integer scalar that will be broadcasted

      call MPI_Bcast(ivar4, I_ONE, MPI_INTEGER, FIRST, MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Bcast_single_int4

!> \brief Wrapper for MPI_Bcast for integer(kind=8)

   subroutine MPI_Bcast_single_int8(ivar8)

      use constants, only: I_ONE
      use mpisetup,  only: FIRST
      use MPIF,      only: MPI_INTEGER8, MPI_COMM_WORLD
      use MPIFUN,    only: MPI_Bcast

      implicit none

      integer(kind=8), intent(inout) :: ivar8  !< integer scalar that will be broadcasted

      call MPI_Bcast(ivar8, I_ONE, MPI_INTEGER8, FIRST, MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Bcast_single_int8

!> \brief Wrapper for MPI_Bcast for real(kind=4)

   subroutine MPI_Bcast_single_real4(rvar4)

      use constants, only: I_ONE
      use mpisetup,  only: FIRST
      use MPIF,      only: MPI_REAL, MPI_COMM_WORLD
      use MPIFUN,    only: MPI_Bcast

      implicit none

      real(kind=4), intent(inout) :: rvar4  !< integer scalar that will be broadcasted

      call MPI_Bcast(rvar4, I_ONE, MPI_REAL, FIRST, MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Bcast_single_real4

!> \brief Wrapper for MPI_Bcast for real(kind=8)

   subroutine MPI_Bcast_single_real8(rvar8)

      use constants, only: I_ONE
      use mpisetup,  only: FIRST
      use MPIF,      only: MPI_DOUBLE_PRECISION, MPI_COMM_WORLD
      use MPIFUN,    only: MPI_Bcast

      implicit none

      real(kind=8), intent(inout) :: rvar8  !< real scalar that will be broadcasted

      call MPI_Bcast(rvar8, I_ONE, MPI_DOUBLE_PRECISION, FIRST, MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Bcast_single_real8

!> \brief Wrapper for MPI_Bcast for real(kind=4), dimension(:)

   subroutine MPI_Bcast_vec_real4(rvar4)

      use mpisetup,  only: FIRST
      use MPIF,   only: MPI_REAL, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Bcast

      implicit none

      real(kind=4), dimension(:), intent(inout) :: rvar4  !< real4 vector that will be broadcasted

      call MPI_Bcast(rvar4, size(rvar4, kind=4), MPI_REAL, FIRST, MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Bcast_vec_real4

!> \brief Wrapper for MPI_Bcast for real(kind=8), dimension(:)

   subroutine MPI_Bcast_vec_real8(rvar8)

      use mpisetup,  only: FIRST
      use MPIF,   only: MPI_DOUBLE_PRECISION, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Bcast

      implicit none

      real(kind=8), dimension(:), intent(inout) :: rvar8  !< real8 vector that will be broadcasted

      call MPI_Bcast(rvar8, size(rvar8, kind=4), MPI_DOUBLE_PRECISION, FIRST, MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Bcast_vec_real8

!> \brief Wrapper for MPI_Bcast fot integer(kind=4), dimension(:)

   subroutine MPI_Bcast_vec_int4(ivar4)

      use mpisetup,  only: FIRST
      use MPIF,   only: MPI_INTEGER, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Bcast

      implicit none

      integer(kind=4), dimension(:), intent(inout) :: ivar4  !< int4 vector that will be broadcasted

      call MPI_Bcast(ivar4, size(ivar4, kind=4), MPI_INTEGER, FIRST, MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Bcast_vec_int4

!> \brief Wrapper for MPI_Bcast for integer(kind=8), dimension(:)

   subroutine MPI_Bcast_vec_int8(ivar8)

      use mpisetup,  only: FIRST
      use MPIF,   only: MPI_INTEGER8, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Bcast

      implicit none

      integer(kind=8), dimension(:), intent(inout) :: ivar8  !< int8 vector that will be broadcasted

      call MPI_Bcast(ivar8, size(ivar8, kind=4), MPI_INTEGER8, FIRST, MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Bcast_vec_int8

!> \brief Wrapper for MPI_Bcast for real(kind=4), dimension(:,:)

   subroutine MPI_Bcast_arr2d_real4(rvar4)

      use mpisetup,  only: FIRST
      use MPIF,   only: MPI_REAL, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Bcast

      implicit none

      real(kind=4), dimension(:,:), intent(inout) :: rvar4  !< real4 arr2d that will be broadcasted

      call MPI_Bcast(rvar4, size(rvar4, kind=4), MPI_REAL, FIRST, MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Bcast_arr2d_real4

!> \brief Wrapper for MPI_Bcast for real(kind=8), dimension(:,:)

   subroutine MPI_Bcast_arr2d_real8(rvar8)

      use mpisetup,  only: FIRST
      use MPIF,   only: MPI_DOUBLE_PRECISION, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Bcast

      implicit none

      real(kind=8), dimension(:,:), intent(inout) :: rvar8  !< real8 arr2d that will be broadcasted

      call MPI_Bcast(rvar8, size(rvar8, kind=4), MPI_DOUBLE_PRECISION, FIRST, MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Bcast_arr2d_real8

!> \brief Wrapper for MPI_Bcast for integer(kind=4), dimension(:, :)

   subroutine MPI_Bcast_arr2d_int4(ivar4)

      use mpisetup,  only: FIRST
      use MPIF,   only: MPI_INTEGER, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Bcast

      implicit none

      integer(kind=4), dimension(:, :), intent(inout) :: ivar4  !< int4 arr2d that will be broadcasted

      call MPI_Bcast(ivar4, size(ivar4, kind=4), MPI_INTEGER, FIRST, MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Bcast_arr2d_int4

!> \brief Wrapper for MPI_Bcast for integer(kind=8), dimension(:,:)

   subroutine MPI_Bcast_arr2d_int8(ivar8)

      use mpisetup,  only: FIRST
      use MPIF,   only: MPI_INTEGER8, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Bcast

      implicit none

      integer(kind=8), dimension(:,:), intent(inout) :: ivar8  !< int8 arr2d that will be broadcasted

      call MPI_Bcast(ivar8, size(ivar8, kind=4), MPI_INTEGER8, FIRST, MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Bcast_arr2d_int8

!> \brief Wrapper for MPI_Bcast for real(kind=4), dimension(:,:,:)

   subroutine MPI_Bcast_arr3d_real4(rvar4)

      use mpisetup,  only: FIRST
      use MPIF,   only: MPI_REAL, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Bcast

      implicit none

      real(kind=4), dimension(:,:,:), intent(inout) :: rvar4  !< real4 arr3d that will be broadcasted

      call MPI_Bcast(rvar4, size(rvar4, kind=4), MPI_REAL, FIRST, MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Bcast_arr3d_real4

!> \brief Wrapper for MPI_Bcast for real(kind=8), dimension(:,:,:)

   subroutine MPI_Bcast_arr3d_real8(rvar8)

      use mpisetup,  only: FIRST
      use MPIF,   only: MPI_DOUBLE_PRECISION, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Bcast

      implicit none

      real(kind=8), dimension(:,:,:), intent(inout) :: rvar8  !< real8 arr3d that will be broadcasted

      call MPI_Bcast(rvar8, size(rvar8, kind=4), MPI_DOUBLE_PRECISION, FIRST, MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Bcast_arr3d_real8

!> \brief Wrapper for MPI_Bcast for integer(kind=4), dimension(:,:,:)

   subroutine MPI_Bcast_arr3d_int4(ivar4)

      use mpisetup,  only: FIRST
      use MPIF,   only: MPI_INTEGER, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Bcast

      implicit none

      integer(kind=4), dimension(:,:, :), intent(inout) :: ivar4  !< int4 arr3d that will be broadcasted

      call MPI_Bcast(ivar4, size(ivar4, kind=4), MPI_INTEGER, FIRST, MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Bcast_arr3d_int4

!> \brief Wrapper for MPI_Bcast for integer(kind=8), dimension(:,:,:)

   subroutine MPI_Bcast_arr3d_int8(ivar8)

      use mpisetup,  only: FIRST
      use MPIF,   only: MPI_INTEGER8, MPI_COMM_WORLD
      use MPIFUN, only: MPI_Bcast

      implicit none

      integer(kind=8), dimension(:,:,:), intent(inout) :: ivar8  !< int8 arr3d that will be broadcasted

      call MPI_Bcast(ivar8, size(ivar8, kind=4), MPI_INTEGER8, FIRST, MPI_COMM_WORLD, err_mpi)

   end subroutine MPI_Bcast_arr3d_int8

#endif /* !NO_F2018 */

end module bcast
