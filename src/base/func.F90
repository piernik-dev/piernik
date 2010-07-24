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
!! \param d shift's direction, where 1,2,3 corresponds to \a x,\a y,\a z respectively
!! \return real, dimension(SIZE(tab,1),SIZE(tab,2),SIZE(tab,3))
!!
!! The function was written in order to significantly improve
!! the performance at the cost of the flexibility of original \p CSHIFT.
!<
   function pshift(tab,d)
      implicit none
      real, dimension(:,:,:) :: tab
      integer :: d
      integer :: lx,ly,lz
      real, dimension(SIZE(tab,1),SIZE(tab,2),SIZE(tab,3)) :: pshift

      lx = SIZE(tab,1)
      ly = SIZE(tab,2)
      lz = SIZE(tab,3)

      if(d==1 .and. lx==1) then
        pshift = tab
        return
      endif
      if(d==2 .and. ly==1) then
        pshift = tab
        return
      endif
      if(d==3 .and. lz==1) then
        pshift = tab
        return
      endif

      if(d==1) then
         pshift(1:lx-1,:,:) = tab(2:lx,:,:); pshift(lx,:,:) = tab(1,:,:)
      else if(d==2) then
         pshift(:,1:ly-1,:) = tab(:,2:ly,:); pshift(:,ly,:) = tab(:,1,:)
      else if(d==3) then
         pshift(:,:,1:lz-1) = tab(:,:,2:lz); pshift(:,:,lz) = tab(:,:,1)
      else
         write(*,*) 'Dim ill defined in pshift!'
      endif

      return
   end function pshift

!>
!! \brief Function mshift makes one-cell, backward circular shift of 3D array in any direction
!! \param tab input array
!! \param d shift's direction, where 1,2,3 corresponds to \a x,\a y,\a z respectively
!! \return real, dimension(SIZE(tab,1),SIZE(tab,2),SIZE(tab,3))
!!
!! The function was written in order to significantly improve
!! the performance at the cost of the flexibility of original \p CSHIFT.
!<
   function mshift(tab,d)
      implicit none
      real, dimension(:,:,:) :: tab
      integer :: d
      integer :: lx,ly,lz
      real, dimension(SIZE(tab,1) , SIZE(tab,2) , SIZE(tab,3)) :: mshift

      lx = SIZE(tab,1)
      ly = SIZE(tab,2)
      lz = SIZE(tab,3)

      if(d==1 .and. lx==1) then
        mshift = tab
        return
      endif
      if(d==2 .and. ly==1) then
        mshift = tab
        return
      endif
      if(d==3 .and. lz==1) then
        mshift = tab
        return
      endif

      if(d==1) then
         mshift(2:lx,:,:) = tab(1:lx-1,:,:); mshift(1,:,:) = tab(lx,:,:)
      else if(d==2) then
         mshift(:,2:ly,:) = tab(:,1:ly-1,:); mshift(:,1,:) = tab(:,ly,:)
      else if(d==3) then
         mshift(:,:,2:lz) = tab(:,:,1:lz-1); mshift(:,:,1) = tab(:,:,lz)
      else
         write(*,*) 'Dim ill defined in mshift!'
      endif

      return
   end function mshift

!-----------------------------------------------------------------------------
#ifdef COSM_RAYS
   subroutine div_v(ifluid)
      use initionized, only : idni,imxi,imyi,imzi
      use fluidindex,  only : nvar
      use fluidindex,  only : iarr_all_dn,iarr_all_mx,iarr_all_my,iarr_all_mz
      use grid,        only : nx,ny,nz
      use grid,        only : dx,dy,dz,nxd,nyd,nzd
      use arrays,      only : u,divvel
      implicit none
      real, dimension(nx) :: vx
      real, dimension(ny) :: vy
      real, dimension(nz) :: vz
      integer             :: i,j,k,ifluid
      integer             :: idnf,imxf,imyf,imzf

      idnf = iarr_all_dn(ifluid)
      imxf = iarr_all_mx(ifluid)
      imyf = iarr_all_my(ifluid)
      imzf = iarr_all_mz(ifluid)

      divvel(:,:,:) = 0.0

      if(nxd /= 1) then
         do k = 1, nz
            do j = 1, ny
               vx = u(imxf,:,j,k) / u(idnf,:,j,k)
               divvel(2:nx-1,j,k) = ( vx(3:nx) - vx(1:nx-2) )  / (2.*dx)
            enddo
         enddo
         divvel(1,:,:) = divvel(2,:,:); divvel(nx,:,:) = divvel(nx-1,:,:) ! for sanity
      endif

      if(nyd /= 1) then
         do k = 1, nz
            do i = 1, nx
               vy = u(imyf,i,:,k) / u(idnf,i,:,k)
               divvel(i,2:ny-1,k) = divvel(i,2:ny-1,k)+( vy(3:ny) - vy(1:ny-2) )  / (2.*dy)
            enddo
         enddo
         divvel(:,1,:) = divvel(:,2,:); divvel(:,ny,:) = divvel(:,ny-1,:) ! for sanity
      endif

      if(nzd /= 1) then
         do j = 1, ny
            do i = 1, nx
               vz = u(imzf,i,j,:) / u(idnf,i,j,:)
               divvel(i,j,2:nz-1) = divvel(i,j,2:nz-1)+( vz(3:nz) - vz(1:nz-2) )  / (2.*dz)
            enddo
         enddo
         divvel(:,:,1) = divvel(:,:,2); divvel(:,:,nz) = divvel(:,:,nz-1) ! for sanity
      endif

   end subroutine div_v


   subroutine div_vx(k,j)
      use grid,        only : nx
      use arrays,      only : divvel
      implicit none
      real,dimension(nx) :: divv
      integer :: j,k

      divv = divvel(:,j,k)

   end subroutine div_vx

   subroutine div_vy(k,i)
      use grid,        only : ny
      use arrays, only      : divvel
      implicit none
      real,dimension(ny) :: divv
      integer :: i,k

      divv = divvel(i,:,k)

   end subroutine div_vy

   subroutine  div_vz(j,i)
      use grid,        only : nz
      use arrays, only      : divvel
      implicit none
      real,dimension(nz) :: divv
      integer :: i,j

      divv = divvel(i,j,:)

   end subroutine div_vz

   subroutine whichfaq(faq,i,j,n)
      implicit none
      real :: faq
      integer :: i,j,n

      faq = 0.5
      if(i == 0) then
         i   = 1
         faq = 1.0
      endif
      if(j-1 == n) then
         j   = n
         faq = 1.0
      endif

   end subroutine whichfaq

#endif /* COSM_RAYS  */

   subroutine compare_namelist(nml_bef, nml_aft, output)
      implicit none
      character(len=*), intent(in)     :: nml_bef, nml_aft, output
      integer                          :: io
      character(len=256)               :: sa, sb

      open(501,file=nml_bef, status='old')
      open(502,file=nml_aft, status='old')
      open(503,file="temp3.dat", status='unknown')
      io = 0
      do
         read(501,'(a)', iostat=io) sa
         read(502,'(a)', iostat=io) sb
         if(io/=0) exit
         if((sa/=sb)) then
            write(503,'(a1,a)') '*',trim(sb)
         else
            write(503,'(a1,a)') ' ',trim(sb)
         endif
      enddo
      close(502,status="delete")
      close(501,status="delete")
      io = 0
      rewind(503)
      open(3, file=output, position='append')
         do while(io == 0)
            read(503,'(a)',iostat=io) sb
            write(3,'(a)') trim(sb)
         enddo
      close(3)
      close(503,status="delete")
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
         if( is_lowercase(c) .or. is_uppercase(c) .or. is_digit(c) .or. c=='_' .or. c=='-' ) then
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

      if(ichar(c) >= 97 .and. ichar(c) <= 122) is_lowercase=.true.
   end function is_lowercase

   logical function is_uppercase(c)
      implicit none
      character(len=1), intent(in) :: c

      is_uppercase = .false.

      if(ichar(c) >= 65 .and. ichar(c) <= 90) is_uppercase=.true.
   end function is_uppercase

   logical function is_digit(c)
      implicit none
      character(len=1), intent(in) :: c

      is_digit = .false.

      if(ichar(c) >= 48 .and. ichar(c) <= 57) is_digit=.true.
   end function is_digit

end module func
