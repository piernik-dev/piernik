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
!    Initial implemetation of PIERNIK code was based on TVD split MHD code by
!    Ue-Li Pen
!        see: Pen, Arras & Wong (2003) for algorithm and
!             http://www.cita.utoronto.ca/~pen/MHD
!             for original source code "mhd.f90"
!
!    For full list of developers see $PIERNIK_HOME/license/pdt.txt
!
#include "piernik.def"
!>
!! \brief (KK) Module that contains unclassified functions
!!
!! This module should be empty. Every function or subroutine placed here belong
!! elsewhere. We are yet unsure where to put them.
!! \todo Move all structures elsewhere
!! \warning Procedures \a dipol and \a rn_angles were moved to sn_sources.F90
!<
module func

   implicit none
   contains

!>
!! \brief Function pshift makes one-cell,foreward circular shift of 3D array in any direction
!! \param tab input array
!! \param d direction of the shift, where 1,2,3 corresponds to \a x,\a y,\a z respectively
!! \return real, dimension(SIZE(tab,1),SIZE(tab,2),SIZE(tab,3))
!!
!! The function was written in order to significantly improve
!! the performance at the cost of the flexibility of original \p CSHIFT.
!<
   function pshift(tab,d)
      use dataio_public, only: warn
      implicit none
      real, dimension(:,:,:) :: tab
      integer :: d
      integer :: lx,ly,lz
      real, dimension(SIZE(tab,1),SIZE(tab,2),SIZE(tab,3)) :: pshift

      lx = SIZE(tab,1)
      ly = SIZE(tab,2)
      lz = SIZE(tab,3)

      if (d==1 .and. lx==1) then
        pshift = tab
        return
      endif
      if (d==2 .and. ly==1) then
        pshift = tab
        return
      endif
      if (d==3 .and. lz==1) then
        pshift = tab
        return
      endif

      if (d==1) then
         pshift(1:lx-1,:,:) = tab(2:lx,:,:); pshift(lx,:,:) = tab(1,:,:)
      else if (d==2) then
         pshift(:,1:ly-1,:) = tab(:,2:ly,:); pshift(:,ly,:) = tab(:,1,:)
      else if (d==3) then
         pshift(:,:,1:lz-1) = tab(:,:,2:lz); pshift(:,:,lz) = tab(:,:,1)
      else
         call warn('[func:pshift]: Dim ill defined in pshift!')
      endif

      return
   end function pshift

!>
!! \brief Function mshift makes one-cell, backward circular shift of 3D array in any direction
!! \param tab input array
!! \param d direction of the shift, where 1,2,3 corresponds to \a x,\a y,\a z respectively
!! \return real, dimension(SIZE(tab,1),SIZE(tab,2),SIZE(tab,3))
!!
!! The function was written in order to significantly improve
!! the performance at the cost of the flexibility of original \p CSHIFT.
!<
   function mshift(tab,d)
      use dataio_public, only: warn
      implicit none
      real, dimension(:,:,:) :: tab
      integer :: d
      integer :: lx,ly,lz
      real, dimension(SIZE(tab,1) , SIZE(tab,2) , SIZE(tab,3)) :: mshift

      lx = SIZE(tab,1)
      ly = SIZE(tab,2)
      lz = SIZE(tab,3)

      if (d==1 .and. lx==1) then
        mshift = tab
        return
      endif
      if (d==2 .and. ly==1) then
        mshift = tab
        return
      endif
      if (d==3 .and. lz==1) then
        mshift = tab
        return
      endif

      if (d==1) then
         mshift(2:lx,:,:) = tab(1:lx-1,:,:); mshift(1,:,:) = tab(lx,:,:)
      else if (d==2) then
         mshift(:,2:ly,:) = tab(:,1:ly-1,:); mshift(:,1,:) = tab(:,ly,:)
      else if (d==3) then
         mshift(:,:,2:lz) = tab(:,:,1:lz-1); mshift(:,:,1) = tab(:,:,lz)
      else
         call warn('[func:mshift]: Dim ill defined in mshift!')
      endif

      return
   end function mshift

!-----------------------------------------------------------------------------
   subroutine compare_namelist(nml_bef, nml_aft)

      use dataio_public, only: cwd, msg, maxparfilelen, printinfo

      implicit none

      character(len=*), intent(in)     :: nml_bef, nml_aft
      integer                          :: io
      character(len=maxparfilelen)     :: sa, sb
      integer, parameter               :: lun_bef=501, lun_aft=502

      open(lun_bef, file=nml_bef, status='old')
      open(lun_aft, file=nml_aft, status='old')
      io = 0
      do
         read(lun_bef,'(a)', iostat=io) sa
         read(lun_aft,'(a)', iostat=io) sb
         if (io/=0) exit
         if ((sa/=sb)) then
            write(msg,'(a1,a)') '*',trim(sb)
         else
            write(msg,'(a1,a)') ' ',trim(sb)
         endif
         call printinfo(msg, .false.)
      enddo
      close(lun_aft, status="delete")
      close(lun_bef, status="delete")

   end subroutine compare_namelist

   function fix_string(str) result (outstr)
      implicit none
      character(len=*), intent(in)  :: str
      character(len=len(str)) :: outstr

      integer          :: i
      character(len=1) :: c

      do i=1, len(str)
         outstr(i:i) = " "
      enddo

      do i=1, len(str)
         c = str(i:i)
         outstr(i:i) = ''
         if ( is_lowercase(c) .or. is_uppercase(c) .or. is_digit(c) .or. c=='_' .or. c=='-' ) then
            outstr(i:i) = c
         else
            return
         endif
      enddo
      return
   end function fix_string

   logical function is_lowercase(c)
      implicit none
      character(len=1), intent(in) :: c

      is_lowercase = .false.

      if (ichar(c) >= 97 .and. ichar(c) <= 122) is_lowercase=.true.
   end function is_lowercase

   logical function is_uppercase(c)
      implicit none
      character(len=1), intent(in) :: c

      is_uppercase = .false.

      if (ichar(c) >= 65 .and. ichar(c) <= 90) is_uppercase=.true.
   end function is_uppercase

   logical function is_digit(c)
      implicit none
      character(len=1), intent(in) :: c

      is_digit = .false.

      if (ichar(c) >= 48 .and. ichar(c) <= 57) is_digit=.true.
   end function is_digit

end module func
