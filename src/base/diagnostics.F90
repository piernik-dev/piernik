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
module diagnostics
   implicit none

   interface my_allocate
      module procedure allocate_array_1D_int
      module procedure allocate_array_2D_int
      module procedure allocate_array_3D_int
      module procedure allocate_array_1D_real
      module procedure allocate_array_2D_real
      module procedure allocate_array_3D_real
      module procedure allocate_array_4D_real
      module procedure allocate_array_5D_real
   end interface my_allocate

   interface my_deallocate
      module procedure deallocate_array_1D_int
      module procedure deallocate_array_2D_int
      module procedure deallocate_array_3D_int
      module procedure deallocate_array_1D_real
      module procedure deallocate_array_2D_real
      module procedure deallocate_array_3D_real
      module procedure deallocate_array_4D_real
      module procedure deallocate_array_5D_real
   end interface my_deallocate

   interface incr_vec
      module procedure increase_char_vector
      module procedure increase_real_vector
   end interface incr_vec

   interface pop_vector
      module procedure pop_char_vector
      module procedure pop_real_vector
   end interface pop_vector

   private
   public :: diagnose_arrays, ma1d, ma2d, ma3d, ma4d, ma5d, my_allocate, my_deallocate, pop_vector

   real,    parameter :: MiB = 8./1048576.  ! sizeof(double) / 2**20
   integer, parameter :: an_len = 64
   real,    parameter :: big_float =  huge(real(1.0,4))
   integer, parameter :: big_int   =  huge(int(1,4))
   real, save         :: used_memory = 0.0

   real, dimension(:), allocatable                  :: array_sizes
   character(len=an_len), dimension(:), allocatable :: array_names
   integer, dimension(1) :: ma1d
   integer, dimension(2) :: ma2d
   integer, dimension(3) :: ma3d
   integer, dimension(4) :: ma4d
   integer, dimension(5) :: ma5d

contains
   subroutine diagnose_arrays
      use dataio_pub,    only: printinfo, warn, msg
      implicit none
      integer :: i

      if (allocated(array_names)) then
         write(msg,'(a,I3,a)') "[diagnostics:diagnose_arrays]: I am aware of ",size(array_names)," arrays..."
         call printinfo(msg)

         do i = LBOUND(array_names,1), UBOUND(array_names,1)
            write(msg,'(3a,F7.3,a)') "Array ",trim(array_names(i))," has ",array_sizes(i)," MiB"
            call printinfo(msg)
         enddo
         write(msg,'(a,F8.3,a)') "[diagnostics:diagnose_arrays]: Total memory used = ", used_memory," MiB"
         call printinfo(msg)
      else
         call warn("[diagnostics:diagnose_arrays]: I am not aware of any arrays :( ")
      endif

   end subroutine diagnose_arrays

   subroutine keep_track_of_arrays(new_size,new_name)
      implicit none
      real, intent(in)             :: new_size
      character(len=*), intent(in) :: new_name

      call incr_vec(array_sizes,1)
      call incr_vec(array_names,1,an_len)

      array_sizes(ubound(array_sizes)) = new_size
      array_names(ubound(array_names)) = new_name

   end subroutine keep_track_of_arrays

   subroutine pop_char_vector(vec,lensize,words)
      use dataio_pub,    only: die
      implicit none
      integer, intent(in) :: lensize
      character(len=*), intent(in), dimension(:) :: words
      character(len=lensize), dimension(:), allocatable, intent(inout) :: vec
      integer :: old , i

      old = 0
      do i = 1, size(words)
         if (len_trim(words(i)) > lensize) call die("[diagnostics:pop_char_vector] word > lensize")
      enddo

      if (allocated(vec)) old = size(vec)
      call incr_vec(vec,size(words),lensize)
      vec(old+1:old+size(words)) = words(:)
      return

   end subroutine pop_char_vector

   subroutine pop_real_vector(vec,words)
      implicit none
      real, intent(in), dimension(:) :: words
      real, dimension(:), allocatable, intent(inout) :: vec
      integer :: old

      old = 0
      if (allocated(vec)) old = size(vec)
      call incr_vec(vec,size(words))
      vec(old+1:old+size(words)) = words(:)
      return

   end subroutine pop_real_vector

   subroutine increase_char_vector(vec,addlen,lensize)
!> \todo get rid of lensize
      implicit none
      integer, intent(in) :: lensize, addlen
      character(len=lensize), dimension(:), allocatable, intent(inout) :: vec
      character(len=lensize), dimension(:), allocatable :: temp
      integer :: old_size

      if (.not.allocated(vec)) then
         allocate(vec(addlen))
         vec = ''
      else
         old_size = size(vec)
         allocate(temp(old_size))
         temp = vec
         deallocate(vec)
         allocate(vec(old_size+addlen)) !> \deprecated BEWARE: vec not deallocated
         vec = ''
         vec(:old_size) = temp
         deallocate(temp)
      endif
   end subroutine increase_char_vector

   subroutine increase_real_vector(vec,addlen)
      implicit none
      integer, intent(in) :: addlen
      real, dimension(:), allocatable, intent(inout) :: vec
      real, dimension(:), allocatable :: temp
      integer :: old_size

      if (.not.allocated(vec)) then
         allocate(vec(addlen))
         vec = 0.0
      else
         old_size = size(vec)
         allocate(temp(old_size))
         temp = vec
         deallocate(vec)
         allocate(vec(old_size+addlen)) !> \deprecated BEWARE: vec not deallocated
         vec = 0.0
         vec(:old_size) = temp
         deallocate(temp)
      endif
   end subroutine increase_real_vector

   ! GOD I NEED TEMPLATES IN FORTRAN!!!!

   subroutine deallocate_array_1D_int(array)
      implicit none
      integer, dimension(:), allocatable, intent(inout)  :: array

      used_memory = used_memory - size(array)*MiB*0.5
      if (allocated(array)) deallocate(array)

   end subroutine deallocate_array_1D_int

   subroutine deallocate_array_2D_int(array)
      implicit none
      integer, dimension(:,:), allocatable, intent(inout)  :: array

      used_memory = used_memory - size(array)*MiB*0.5
      if (allocated(array)) deallocate(array)

   end subroutine deallocate_array_2D_int

   subroutine deallocate_array_3D_int(array)
      implicit none
      integer, dimension(:,:,:), allocatable, intent(inout)  :: array

      used_memory = used_memory - size(array)*MiB*0.5
      if (allocated(array)) deallocate(array)

   end subroutine deallocate_array_3D_int

   subroutine deallocate_array_1D_real(array)
      implicit none
      real, dimension(:), allocatable, intent(inout)  :: array

      used_memory = used_memory - size(array)*MiB
      if (allocated(array)) deallocate(array)

   end subroutine deallocate_array_1D_real

   subroutine deallocate_array_2D_real(array)
      implicit none
      real, dimension(:,:), allocatable, intent(inout)  :: array

      used_memory = used_memory - size(array)*MiB
      if (allocated(array)) deallocate(array)

   end subroutine deallocate_array_2D_real

   subroutine deallocate_array_3D_real(array)
      implicit none
      real, dimension(:,:,:), allocatable, intent(inout)  :: array

      used_memory = used_memory - size(array)*MiB
      if (allocated(array)) deallocate(array)

   end subroutine deallocate_array_3D_real

   subroutine deallocate_array_4D_real(array)
      implicit none
      real, dimension(:,:,:,:), allocatable, intent(inout)  :: array

      used_memory = used_memory - size(array)*MiB
      if (allocated(array)) deallocate(array)

   end subroutine deallocate_array_4D_real

   subroutine deallocate_array_5D_real(array)
      implicit none
      real, dimension(:,:,:,:,:), allocatable, intent(inout)  :: array

      used_memory = used_memory - size(array)*MiB
      if (allocated(array)) deallocate(array)

   end subroutine deallocate_array_5D_real

   subroutine allocate_array_1D_int(array,as,aname)
      implicit none
      integer, dimension(:), allocatable, intent(inout)  :: array
      integer, dimension(1), intent(in)                  :: as
      character(len=*), intent(in), optional             :: aname

      if (.not.allocated(array)) then
         allocate( array(as(1)) )
         array = big_int
      endif
      used_memory = used_memory + size(array)*MiB*0.5
      if (present(aname)) call keep_track_of_arrays(size(array)*MiB*0.5,aname)

   end subroutine allocate_array_1D_int

   subroutine allocate_array_2D_int(array,as,aname)
      implicit none
      integer, dimension(:,:), allocatable, intent(inout)  :: array
      integer, dimension(2), intent(in)                    :: as
      character(len=*), intent(in), optional               :: aname

      if (.not.allocated(array)) then
         allocate( array(as(1),as(2)) )
         array = big_int
      endif
      used_memory = used_memory + size(array)*MiB*0.5
      if (present(aname)) call keep_track_of_arrays(size(array)*MiB*0.5,aname)

   end subroutine allocate_array_2D_int

   subroutine allocate_array_3D_int(array,as,aname)
      implicit none
      integer, dimension(:,:,:), allocatable, intent(inout)  :: array
      integer, dimension(3), intent(in)                      :: as
      character(len=*), intent(in), optional                 :: aname

      if (.not.allocated(array)) then
         allocate( array(as(1),as(2),as(3)) )
         array = big_int
      endif
      used_memory = used_memory + size(array)*MiB*0.5
      if (present(aname)) call keep_track_of_arrays(size(array)*MiB*0.5,aname)

   end subroutine allocate_array_3D_int

   subroutine allocate_array_1D_real(array,as,aname)
      implicit none
      real, dimension(:), allocatable, intent(inout)  :: array
      integer, dimension(1), intent(in)               :: as
      character(len=*), intent(in), optional          :: aname

      if (.not.allocated(array)) then
         allocate( array(as(1)) )
         array = big_float
      endif
      used_memory = used_memory + size(array)*MiB
      if (present(aname)) call keep_track_of_arrays(size(array)*MiB,aname)

   end subroutine allocate_array_1D_real

   subroutine allocate_array_2D_real(array,as,aname)
      implicit none
      real, dimension(:,:), allocatable, intent(inout)  :: array
      integer, dimension(2), intent(in)                 :: as
      character(len=*), intent(in), optional            :: aname

      if (.not.allocated(array)) then
         allocate( array(as(1),as(2)) )
         array = big_float
      endif
      used_memory = used_memory + size(array)*MiB
      if (present(aname)) call keep_track_of_arrays(size(array)*MiB,aname)

   end subroutine allocate_array_2D_real

   subroutine allocate_array_3D_real(array,as,aname)
      implicit none
      real, dimension(:,:,:), allocatable, intent(inout)  :: array
      integer, dimension(3), intent(in)                   :: as
      character(len=*), intent(in), optional              :: aname

      if (.not.allocated(array)) then
         allocate( array(as(1),as(2),as(3)) )
         array = big_float
      endif
      used_memory = used_memory + size(array)*MiB
      if (present(aname)) call keep_track_of_arrays(size(array)*MiB,aname)

   end subroutine allocate_array_3D_real

   subroutine allocate_array_4D_real(array,as,aname)
      implicit none
      real, dimension(:,:,:,:), allocatable, intent(inout)  :: array
      integer, dimension(4), intent(in)                     :: as
      character(len=*), intent(in), optional                :: aname

      if (.not.allocated(array)) then
         allocate( array(as(1),as(2),as(3),as(4)) )
         array = big_float
      endif
      used_memory = used_memory + size(array)*MiB
      if (present(aname)) call keep_track_of_arrays(size(array)*MiB,aname)

   end subroutine allocate_array_4D_real

   subroutine allocate_array_5D_real(array,as,aname)
      implicit none
      real, dimension(:,:,:,:,:), allocatable, intent(inout)  :: array
      integer, dimension(5), intent(in)                       :: as
      character(len=*), intent(in), optional                  :: aname

      if (.not.allocated(array)) then
         allocate( array(as(1),as(2),as(3),as(4),as(5)) )
         array = big_float
      endif
      used_memory = used_memory + size(array)*MiB
      if (present(aname)) call keep_track_of_arrays(size(array)*MiB,aname)

   end subroutine allocate_array_5D_real

end module diagnostics
