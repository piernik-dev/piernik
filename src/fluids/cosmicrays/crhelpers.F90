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

module crhelpers
! pulled by COSM_RAYS
   implicit none

   private
   public :: div_v, div_vx, div_vy, div_vz, whichfaq

contains

   subroutine div_v(ifluid)

      use arrays,      only: u, divvel
      use fluidindex,  only: iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz
      use grid,        only: cg
      use mpisetup,    only: xdim, ydim, zdim, has_dir

      implicit none

      real, dimension(cg%nx) :: vx
      real, dimension(cg%ny) :: vy
      real, dimension(cg%nz) :: vz
      integer                :: i, j, k, ifluid
      integer                :: idnf, imxf, imyf, imzf

      idnf = iarr_all_dn(ifluid)
      imxf = iarr_all_mx(ifluid)
      imyf = iarr_all_my(ifluid)
      imzf = iarr_all_mz(ifluid)

      divvel(:,:,:) = 0.0

      if (has_dir(xdim)) then
         do k = 1, cg%nz
            do j = 1, cg%ny
               vx = u(imxf,:,j,k) / u(idnf,:,j,k)
               divvel(2:cg%nx-1,j,k) = ( vx(3:cg%nx) - vx(1:cg%nx-2) )  / (2.*cg%dx)
            enddo
         enddo
         divvel(1,:,:) = divvel(2,:,:); divvel(cg%nx,:,:) = divvel(cg%nx-1,:,:) ! for sanity
      endif

      if (has_dir(ydim)) then
         do k = 1, cg%nz
            do i = 1, cg%nx
               vy = u(imyf,i,:,k) / u(idnf,i,:,k)
               divvel(i,2:cg%ny-1,k) = divvel(i,2:cg%ny-1,k)+( vy(3:cg%ny) - vy(1:cg%ny-2) )  / (2.*cg%dy)
            enddo
         enddo
         divvel(:,1,:) = divvel(:,2,:); divvel(:, cg%ny,:) = divvel(:, cg%ny-1,:) ! for sanity
      endif

      if (has_dir(zdim)) then
         do j = 1, cg%ny
            do i = 1, cg%nx
               vz = u(imzf,i,j,:) / u(idnf,i,j,:)
               divvel(i,j,2:cg%nz-1) = divvel(i,j,2:cg%nz-1)+( vz(3:cg%nz) - vz(1:cg%nz-2) )  / (2.*cg%dz)
            enddo
         enddo
         divvel(:,:,1) = divvel(:,:,2); divvel(:,:, cg%nz) = divvel(:,:, cg%nz-1) ! for sanity
      endif

   end subroutine div_v

   subroutine div_vx(k,j)

      use arrays, only: divvel
      use grid,   only: cg

      implicit none

      real, dimension(cg%nx) :: divv
      integer                :: j, k

      divv = divvel(:,j,k)

   end subroutine div_vx

   subroutine div_vy(k,i)

      use arrays, only: divvel
      use grid,   only: cg

      implicit none

      real, dimension(cg%ny) :: divv
      integer                :: i, k

      divv = divvel(i,:,k)

   end subroutine div_vy

   subroutine  div_vz(j,i)

      use arrays, only: divvel
      use grid,   only: cg

      implicit none

      real, dimension(cg%nz) :: divv
      integer                :: i, j

      divv = divvel(i,j,:)

   end subroutine div_vz

   subroutine whichfaq(faq,i,j,n)
      implicit none
      real    :: faq
      integer :: i, j, n

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
