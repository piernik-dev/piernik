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
      module procedure allocate_array_1D_int4
      module procedure allocate_array_2D_int4
      module procedure allocate_array_3D_int4
      module procedure allocate_array_1D_real
      module procedure allocate_array_2D_real
      module procedure allocate_array_3D_real
      module procedure allocate_array_4D_real
      module procedure allocate_array_5D_real
   end interface my_allocate

   interface my_deallocate
      module procedure deallocate_array_1D_int4
      module procedure deallocate_array_2D_int4
      module procedure deallocate_array_3D_int4
      module procedure deallocate_array_1D_real
      module procedure deallocate_array_2D_real
      module procedure deallocate_array_3D_real
      module procedure deallocate_array_4D_real
      module procedure deallocate_array_5D_real
   end interface my_deallocate

   interface incr_vec
      module procedure increase_char_vector
      module procedure increase_real_vector
      module procedure increase_int8_vector
   end interface incr_vec

   interface pop_vector
      module procedure pop_char_vector
      module procedure pop_real_vector
   end interface pop_vector

   private
   public :: diagnose_arrays, ma1d, ma2d, ma3d, ma4d, ma5d, my_allocate, my_deallocate, pop_vector, check_environment, cleanup_diagnostics

   integer(kind=8), parameter :: i4_s=4, r8_s=8 ! sizeof(int(kind=4)), sizeof(double)
!   real,    parameter :: MiB = 8./1048576.  ! sizeof(double) / 2**20
   integer, parameter :: an_len = 64
   integer(kind=8), save :: used_memory = 0

   integer(kind=8), dimension(:), allocatable :: array_sizes
   character(len=an_len), dimension(:), allocatable :: array_names
   integer(kind=4), dimension(1) :: ma1d
   integer(kind=4), dimension(2) :: ma2d
   integer(kind=4), dimension(3) :: ma3d
   integer(kind=4), dimension(4) :: ma4d
   integer(kind=4), dimension(5) :: ma5d

contains

   subroutine cleanup_diagnostics

      implicit none

      ! Unfotrunately deallocating through my_deallocate does not erase the entry in array_names
      !call diagnose_arrays

      if (allocated(array_names)) deallocate(array_names)
      if (allocated(array_sizes)) deallocate(array_sizes)

   end subroutine cleanup_diagnostics

   subroutine check_environment
      implicit none

      call check_lhs_realloc

   end subroutine check_environment

   subroutine check_lhs_realloc
      use dataio_pub,   only: die
      implicit none
      integer, dimension(:), allocatable :: ivec
      integer, parameter :: n = 1
      integer :: i, nsize

      allocate(ivec(n))
      ivec  = [(i,i=1,n+1)]       ! valgrind: Invalid write of size 4
      nsize = abs(n-size(ivec))
      deallocate(ivec)

      ! \todo Decide are we using lhs realloc or not
      !   Against:
      !     * the whole code needs to be reviewed for possible unintentional usage
      !     * requires heightened awareness, i.e. BE C4R3FUL WH4T Y0U TYP3 !!
      !   In favour:
      !     * it can be handy (especially in mc module)
      !     * supported by both gcc and intel

!      if (nsize == 0) call die("[diagnostics:check_lhs_realloc]: No lhs realloc!")
      if (nsize /= 0) call die("[diagnostics:check_lhs_realloc]: LHS realloc detected.")

   end subroutine check_lhs_realloc

   subroutine diagnose_arrays

      use constants,  only: fplen
      use dataio_pub, only: printinfo, msg

      implicit none

      integer :: i
      character(len=fplen) :: sstr

      if (allocated(array_names)) then

         write(msg,'(a,i4,a)') "[diagnostics:diagnose_arrays]: I am aware of ",size(array_names)," arrays..."
         call printinfo(msg)

         do i = lbound(array_names,1), ubound(array_names,1)
            call size2str(array_sizes(i), sstr)
            write(msg,'(4a)') "Array ",trim(array_names(i))," has ", trim(sstr)
            call printinfo(msg)
         enddo

         call size2str(used_memory, sstr)
         write(msg,'(2a)') "[diagnostics:diagnose_arrays]: Total memory used = ", trim(sstr)
         call printinfo(msg)

      endif

   end subroutine diagnose_arrays

!>
!! \brief Print the size in bytes together with a binary prefix, when applicable
!!
!! \details When size is less than 2**1 no prefix is printed, only "B" that stands for bytes.
!! KiB prefix is used for sizes from 2**11 to 2**21-1, then MiB, GiB, etc. fo larger sizes.
!! ZiB and YiB are out of 8-byte integer range.
!<

   subroutine size2str(size, str)

      implicit none

      integer(kind=8), intent(in) :: size
      character(len=*), intent(out) :: str

      real, parameter :: over = 2.
      real, parameter :: Ki = 2.**10
      character(len=*), parameter :: prefixes = "KMGTPEZY"
      integer ::  bin

      if (abs(size) < over*Ki) then
         write(str,'(i5,a)') size," B"
      else
         bin = int(log(abs(size)/over)/log(Ki))
         if (bin < len(prefixes)) then
            write(str,'(F8.2,3a)') size/Ki**bin, " ",prefixes(bin:bin),"iB"
         else
            write(str, *) size, " B"
         endif
      endif

   end subroutine size2str

   subroutine keep_track_of_arrays(new_size,new_name)

      implicit none

      integer(kind=8),  intent(in) :: new_size
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

   end subroutine pop_real_vector

   subroutine increase_char_vector(vec,addlen,lensize)

      implicit none

      integer, intent(in) :: lensize, addlen !> \todo get rid of lensize
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

   subroutine increase_int8_vector(vec,addlen)

      implicit none

      integer, intent(in) :: addlen
      integer(kind=8), dimension(:), allocatable, intent(inout) :: vec
      integer(kind=8), dimension(:), allocatable :: temp
      integer :: old_size

      if (.not.allocated(vec)) then
         allocate(vec(addlen))
         vec = 0
      else
         old_size = size(vec)
         allocate(temp(old_size))
         temp = vec
         deallocate(vec)
         allocate(vec(old_size+addlen)) !> \deprecated BEWARE: vec not deallocated
         vec = 0
         vec(:old_size) = temp
         deallocate(temp)
      endif
   end subroutine increase_int8_vector

   ! GOD I NEED TEMPLATES IN FORTRAN!!!!

   subroutine deallocate_array_1D_int4(array)
      implicit none
      integer(kind=4), dimension(:), allocatable, intent(inout)  :: array

      used_memory = used_memory - size(array)*i4_s
      if (allocated(array)) deallocate(array)

   end subroutine deallocate_array_1D_int4

   subroutine deallocate_array_2D_int4(array)
      implicit none
      integer(kind=4), dimension(:,:), allocatable, intent(inout)  :: array

      used_memory = used_memory - size(array)*i4_s
      if (allocated(array)) deallocate(array)

   end subroutine deallocate_array_2D_int4

   subroutine deallocate_array_3D_int4(array)
      implicit none
      integer(kind=4), dimension(:,:,:), allocatable, intent(inout)  :: array

      used_memory = used_memory - size(array)*i4_s
      if (allocated(array)) deallocate(array)

   end subroutine deallocate_array_3D_int4

   subroutine deallocate_array_1D_real(array)
      implicit none
      real, dimension(:), allocatable, intent(inout)  :: array

      used_memory = used_memory - size(array)*r8_s
      if (allocated(array)) deallocate(array)

   end subroutine deallocate_array_1D_real

   subroutine deallocate_array_2D_real(array)
      implicit none
      real, dimension(:,:), allocatable, intent(inout)  :: array

      used_memory = used_memory - size(array)*r8_s
      if (allocated(array)) deallocate(array)

   end subroutine deallocate_array_2D_real

   subroutine deallocate_array_3D_real(array)
      implicit none
      real, dimension(:,:,:), allocatable, intent(inout)  :: array

      used_memory = used_memory - size(array)*r8_s
      if (allocated(array)) deallocate(array)

   end subroutine deallocate_array_3D_real

   subroutine deallocate_array_4D_real(array)
      implicit none
      real, dimension(:,:,:,:), allocatable, intent(inout)  :: array

      used_memory = used_memory - size(array)*r8_s
      if (allocated(array)) deallocate(array)

   end subroutine deallocate_array_4D_real

   subroutine deallocate_array_5D_real(array)
      implicit none
      real, dimension(:,:,:,:,:), allocatable, intent(inout)  :: array

      used_memory = used_memory - size(array)*r8_s
      if (allocated(array)) deallocate(array)

   end subroutine deallocate_array_5D_real

   subroutine allocate_array_1D_int4(array,as,aname)

      use constants, only: big_int

      implicit none

      integer(kind=4), dimension(:), allocatable, intent(inout)  :: array
      integer(kind=4), dimension(1), intent(in)          :: as
      character(len=*), intent(in), optional             :: aname

      if (.not.allocated(array)) then
         allocate( array(as(1)) )
         array = big_int
      endif
      used_memory = used_memory + size(array)*i4_s
      if (present(aname)) call keep_track_of_arrays(size(array)*i4_s,aname)

   end subroutine allocate_array_1D_int4

   subroutine allocate_array_2D_int4(array,as,aname)

      use constants, only: big_int

      implicit none

      integer(kind=4), dimension(:,:), allocatable, intent(inout)  :: array
      integer(kind=4), dimension(2), intent(in)            :: as
      character(len=*), intent(in), optional               :: aname

      if (.not.allocated(array)) then
         allocate( array(as(1),as(2)) )
         array = big_int
      endif
      used_memory = used_memory + size(array)*i4_s
      if (present(aname)) call keep_track_of_arrays(size(array)*i4_s,aname)

   end subroutine allocate_array_2D_int4

   subroutine allocate_array_3D_int4(array,as,aname)

      use constants, only: big_int

      implicit none

      integer(kind=4), dimension(:,:,:), allocatable, intent(inout)  :: array
      integer(kind=4), dimension(3), intent(in)              :: as
      character(len=*), intent(in), optional                 :: aname

      if (.not.allocated(array)) then
         allocate( array(as(1),as(2),as(3)) )
         array = big_int
      endif
      used_memory = used_memory + size(array)*i4_s
      if (present(aname)) call keep_track_of_arrays(size(array)*i4_s,aname)

   end subroutine allocate_array_3D_int4

   subroutine allocate_array_1D_real(array,as,aname)

      use constants, only: big_float

      implicit none

      real, dimension(:), allocatable, intent(inout)  :: array
      integer(kind=4), dimension(1), intent(in)       :: as
      character(len=*), intent(in), optional          :: aname

      if (.not.allocated(array)) then
         allocate( array(as(1)) )
         array = big_float
      endif
      used_memory = used_memory + size(array)*r8_s
      if (present(aname)) call keep_track_of_arrays(size(array)*r8_s,aname)

   end subroutine allocate_array_1D_real

   subroutine allocate_array_2D_real(array,as,aname)

      use constants, only: big_float

      implicit none

      real, dimension(:,:), allocatable, intent(inout)  :: array
      integer(kind=4), dimension(2), intent(in)         :: as
      character(len=*), intent(in), optional            :: aname

      if (.not.allocated(array)) then
         allocate( array(as(1),as(2)) )
         array = big_float
      endif
      used_memory = used_memory + size(array)*r8_s
      if (present(aname)) call keep_track_of_arrays(size(array)*r8_s,aname)

   end subroutine allocate_array_2D_real

   subroutine allocate_array_3D_real(array,as,aname)

      use constants, only: big_float

      implicit none

      real, dimension(:,:,:), allocatable, intent(inout)  :: array
      integer(kind=4), dimension(3), intent(in)           :: as
      character(len=*), intent(in), optional              :: aname

      if (.not.allocated(array)) then
         allocate( array(as(1),as(2),as(3)) )
         array = big_float
      endif
      used_memory = used_memory + size(array)*r8_s
      if (present(aname)) call keep_track_of_arrays(size(array)*r8_s,aname)

   end subroutine allocate_array_3D_real

   subroutine allocate_array_4D_real(array,as,aname)

      use constants, only: big_float

      implicit none

      real, dimension(:,:,:,:), allocatable, intent(inout)  :: array
      integer(kind=4), dimension(4), intent(in)             :: as
      character(len=*), intent(in), optional                :: aname

      if (.not.allocated(array)) then
         allocate( array(as(1),as(2),as(3),as(4)) )
         array = big_float
      endif
      used_memory = used_memory + size(array)*r8_s
      if (present(aname)) call keep_track_of_arrays(size(array)*r8_s,aname)

   end subroutine allocate_array_4D_real

   subroutine allocate_array_5D_real(array,as,aname)

      use constants, only: big_float

      implicit none

      real, dimension(:,:,:,:,:), allocatable, intent(inout)  :: array
      integer(kind=4), dimension(5), intent(in)               :: as
      character(len=*), intent(in), optional                  :: aname

      if (.not.allocated(array)) then
         allocate( array(as(1),as(2),as(3),as(4),as(5)) )
         array = big_float
      endif
      used_memory = used_memory + size(array)*r8_s
      if (present(aname)) call keep_track_of_arrays(size(array)*r8_s,aname)

   end subroutine allocate_array_5D_real

end module diagnostics
