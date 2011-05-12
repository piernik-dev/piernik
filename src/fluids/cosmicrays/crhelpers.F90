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
   public :: div_v, set_div_v1d, cleanup_crhelpers

   real, dimension(:,:,:), allocatable, target :: divvel

contains

   subroutine cleanup_crhelpers
      use diagnostics, only: my_deallocate
      implicit none
      call my_deallocate(divvel)
   end subroutine cleanup_crhelpers

   subroutine set_div_v1d(p, dir, i1, i2)
      use constants,   only: xdim, ydim, zdim
      implicit none
      real, dimension(:), pointer, intent(inout) :: p
      integer, intent(in) :: dir, i1, i2

      select case (dir)
         case (xdim)
            p => divvel(:, i1, i2)
         case (ydim)
            p => divvel(i2, :, i1)
         case (zdim)
            p => divvel(i1, i2, :)
      end select
      return
   end subroutine set_div_v1d

   subroutine div_v(ifluid)

      use diagnostics, only: ma3d, my_allocate
      use arrays,      only: u
      use fluidindex,  only: iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz
      use grid,        only: cg
      use constants,   only: xdim, ydim, zdim
      use mpisetup,    only: has_dir

      implicit none

      real, dimension(cg%nx) :: vx
      real, dimension(cg%ny) :: vy
      real, dimension(cg%nz) :: vz
      integer                :: i, j, k, ifluid
      integer                :: idnf, imxf, imyf, imzf

      if (.not.allocated(divvel)) then
         ma3d = [cg%nx, cg%ny, cg%nz]
         call my_allocate(divvel, ma3d, "divvel")
      endif

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

end module crhelpers
