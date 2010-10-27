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
#include "piernik.def"
!>
!! \brief (KK) Module that contains unclassified functions
!!
!! This module should be empty. Every function or subroutine placed here belong
!! elsewhere. We are yet unsure where to put them.
!! \todo Move all structures elsewhere
!! \warning Procedures \a dipole and \a rn_angles were moved to sn_sources.F90
!<
module func

   implicit none
   contains

!>
!! \brief Function pshift makes one-cell, forward circular shift of 3D array in any direction
!! \param tab input array
!! \param d direction of the shift, where 1,2,3 corresponds to \a x,\a y,\a z respectively
!! \return real, dimension(SIZE(tab,1),SIZE(tab,2),SIZE(tab,3))
!!
!! The function was written in order to significantly improve
!! the performance at the cost of the flexibility of original \p CSHIFT.
!<
   function pshift(tab,d)
      use dataio_pub,    only: warn
      implicit none
      real, dimension(:,:,:) :: tab
      integer :: d
      integer :: ll
      real, dimension(SIZE(tab,1),SIZE(tab,2),SIZE(tab,3)) :: pshift

      ll = SIZE(tab,d)

      if (ll==1) then
        pshift = tab
        return
      endif

      if (d==1) then
         pshift(1:ll-1,:,:) = tab(2:ll,:,:); pshift(ll,:,:) = tab(1,:,:)
      else if (d==2) then
         pshift(:,1:ll-1,:) = tab(:,2:ll,:); pshift(:,ll,:) = tab(:,1,:)
      else if (d==3) then
         pshift(:,:,1:ll-1) = tab(:,:,2:ll); pshift(:,:,ll) = tab(:,:,1)
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
      use dataio_pub,    only: warn
      implicit none
      real, dimension(:,:,:) :: tab
      integer :: d
      integer :: ll
      real, dimension(SIZE(tab,1) , SIZE(tab,2) , SIZE(tab,3)) :: mshift

      ll = SIZE(tab,d)

      if (ll==1) then
        mshift = tab
        return
      endif

      if (d==1) then
         mshift(2:ll,:,:) = tab(1:ll-1,:,:); mshift(1,:,:) = tab(ll,:,:)
      else if (d==2) then
         mshift(:,2:ll,:) = tab(:,1:ll-1,:); mshift(:,1,:) = tab(:,ll,:)
      else if (d==3) then
         mshift(:,:,2:ll) = tab(:,:,1:ll-1); mshift(:,:,1) = tab(:,:,ll)
      else
         call warn('[func:mshift]: Dim ill defined in mshift!')
      endif

      return
   end function mshift

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
