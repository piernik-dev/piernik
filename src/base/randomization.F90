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

!> \brief This module serves an array of prime numbers. Used primarily in decomposition module.

module randomization
! pulled by RANDOMIZE

   implicit none

   private
   public :: initseed, redoseed, seed_size, init_randomization, cleanup_randomization, randoms_redostep
#ifdef HDF5
   public :: read_current_seed_from_restart, write_current_seed_to_restart
#endif /* HDF5 */

   integer, dimension(:), allocatable :: initseed, redoseed
   integer                            :: seed_size

contains

!>
!! \brief Routine to prepare seed for generating random numbers
!<
   subroutine init_randomization

      use constants,  only: INVALID
      use dataio_pub, only: nrestart

      implicit none

      integer :: i

      call random_seed(size=seed_size)
      if (.not.allocated(initseed)) allocate(initseed(seed_size))
      if (.not.allocated(redoseed)) allocate(redoseed(seed_size))

      if (nrestart == INVALID) then
         do i= 1,seed_size; initseed(i)= 2**i; enddo   ! Simple seed 2,4,8,16,..
      endif
      call random_seed(PUT=initseed)

   end subroutine init_randomization

   subroutine cleanup_randomization

      implicit none

      if (allocated(initseed)) deallocate(initseed)
      if (allocated(redoseed)) deallocate(redoseed)


   end subroutine cleanup_randomization

   subroutine randoms_redostep(dorepeat)

      implicit none

      logical, intent(in) :: dorepeat

      if (dorepeat) then
         call random_seed(put=redoseed)
      else
         call random_seed(get=redoseed)
      endif

   end subroutine randoms_redostep

#ifdef HDF5
   subroutine write_current_seed_to_restart(file_id)

      use hdf5, only: HID_T, SIZE_T
      use h5lt, only: h5ltset_attribute_int_f

      implicit none

      integer(HID_T), intent(in)    :: file_id
      integer(SIZE_T)               :: bufsize
      integer(kind=4)               :: error
      integer, dimension(seed_size) :: resseed

      bufsize = 1
      call random_seed(get=resseed)
      bufsize = seed_size
      ! warning: this doesn't work with -fdefault-integer-8 because we don't have proper Fortran H5LT interface for 64-bit integer attributes and random_seed works with the default integer kind
      call h5ltset_attribute_int_f(file_id, "/", "current_seed", resseed, bufsize, error)
      ! There is no h5ltset_attribute_long_f in the library, so it can't compile with -fdefault-integer-8.
      ! This is not something critical, so we can:
      ! * wait for full support,
      ! * rewrite LT calls with more primitive ones,
      ! * exclude 64-bit tests with RANDOMIZE from Jenkins (current workaround).

   end subroutine write_current_seed_to_restart

!-----------------------------------------------------------------------------

   subroutine read_current_seed_from_restart(file_id)

      use cmp_1D_mpi, only: compare_array1D
      use hdf5,       only: HID_T
      use h5lt,       only: h5ltget_attribute_int_f

      implicit none

      integer(HID_T), intent(in)         :: file_id
      integer(kind=4)                    :: error
      integer, dimension(:), allocatable :: resseed

      if (.not.allocated(resseed)) allocate(resseed(seed_size))
      ! warning: this doesn't work with -fdefault-integer-8 because we don't have proper Fortran H5LT interface for 64-bit integer attributes and random_seed works with the default integer kind
      call h5ltget_attribute_int_f(file_id, "/", "current_seed", resseed, error)
      ! There is no h5ltget_attribute_long_f in the library, so it can't compile with -fdefault-integer-8, see comments in write_current_seed_to_restart.

      call compare_array1D(resseed)
      call random_seed(put=resseed)
      deallocate(resseed)

   end subroutine read_current_seed_from_restart
#endif /* HDF5 */

end module randomization
