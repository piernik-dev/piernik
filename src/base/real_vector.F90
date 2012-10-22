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

!> \brief Definition of tyhe vector of real values

module real_vector

   implicit none

   private
   public :: real_vec_T

   !> \brief a type for vector of real data
   type :: real_vec_T
      real, dimension(:), pointer :: r => null()
   contains
      procedure       :: real_vec_T_allocate
      procedure       :: real_vec_T_allocate2
      generic, public :: allocate   => real_vec_T_allocate, real_vec_T_allocate2
      procedure       :: associated => real_vec_T_associated
      procedure       :: deallocate => real_vec_T_deallocate
   end type real_vec_T

contains

!> \brief Check if it is associated

   logical function real_vec_T_associated(this)

      implicit none

      class(real_vec_T), intent(in) :: this

      real_vec_T_associated = associated(this%r)

   end function real_vec_T_associated

!> \brief Simple allocation of given number of elements

   subroutine real_vec_T_allocate(this, n)

      implicit none

      class(real_vec_T), intent(inout) :: this
      integer(kind=4),   intent(in)    :: n

      allocate(this%r(n))

   end subroutine real_vec_T_allocate

!> \brief Allocation of given range of indices

   subroutine real_vec_T_allocate2(this, n1, n2)

      implicit none

      class(real_vec_T), intent(inout) :: this
      integer(kind=4),   intent(in)    :: n1
      integer(kind=4),   intent(in)    :: n2

      allocate(this%r(n1:n2))

   end subroutine real_vec_T_allocate2

!> free the memory

   subroutine real_vec_T_deallocate(this)

      implicit none

      class(real_vec_T), intent(inout) :: this

      deallocate(this%r)
      nullify(this%r)

   end subroutine real_vec_T_deallocate

end module real_vector
