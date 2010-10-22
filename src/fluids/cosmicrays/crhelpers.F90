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
module crhelpers
   implicit none
   private
   public :: div_v, div_vx, div_vy, div_vz, whichfaq

   contains

   subroutine div_v(ifluid)
      use fluidindex,  only: nvar
      use fluidindex,  only: iarr_all_dn,iarr_all_mx,iarr_all_my,iarr_all_mz
      use grid,        only: nx,ny,nz
      use grid,        only: dx,dy,dz,nxd,nyd,nzd
      use arrays,      only: u,divvel
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

      if (nxd /= 1) then
         do k = 1, nz
            do j = 1, ny
               vx = u(imxf,:,j,k) / u(idnf,:,j,k)
               divvel(2:nx-1,j,k) = ( vx(3:nx) - vx(1:nx-2) )  / (2.*dx)
            enddo
         enddo
         divvel(1,:,:) = divvel(2,:,:); divvel(nx,:,:) = divvel(nx-1,:,:) ! for sanity
      endif

      if (nyd /= 1) then
         do k = 1, nz
            do i = 1, nx
               vy = u(imyf,i,:,k) / u(idnf,i,:,k)
               divvel(i,2:ny-1,k) = divvel(i,2:ny-1,k)+( vy(3:ny) - vy(1:ny-2) )  / (2.*dy)
            enddo
         enddo
         divvel(:,1,:) = divvel(:,2,:); divvel(:,ny,:) = divvel(:,ny-1,:) ! for sanity
      endif

      if (nzd /= 1) then
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
      use grid,        only: nx
      use arrays,      only: divvel
      implicit none
      real,dimension(nx) :: divv
      integer :: j,k

      divv = divvel(:,j,k)

   end subroutine div_vx

   subroutine div_vy(k,i)
      use grid,        only: ny
      use arrays, only: divvel
      implicit none
      real,dimension(ny) :: divv
      integer :: i,k

      divv = divvel(i,:,k)

   end subroutine div_vy

   subroutine  div_vz(j,i)
      use grid,        only: nz
      use arrays, only: divvel
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
      if (i == 0) then
         i   = 1
         faq = 1.0
      endif
      if (j-1 == n) then
         j   = n
         faq = 1.0
      endif

   end subroutine whichfaq

end module crhelpers
