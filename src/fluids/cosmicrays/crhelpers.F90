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
   public :: div_v, set_div_v1d

contains

   subroutine set_div_v1d(p, dir, i1, i2, cg)

      use constants,   only: xdim, ydim, zdim
      use cr_data,     only: divv_n
      use dataio_pub,  only: die
      use grid_cont,   only: grid_container

      implicit none

      integer(kind=4), intent(in) :: dir
      integer, intent(in) :: i1, i2
      real, dimension(:), pointer, intent(inout) :: p
      type(grid_container), pointer, intent(in) :: cg

      real, dimension(:,:,:), pointer :: divvel

      divvel => cg%ptr(divv_n)
      if (.not. associated(divvel)) call die("[crhelpers:set_div_v1d] cannot get divvel")

      select case (dir)
         case (xdim)
            p => divvel(:, i1, i2)
         case (ydim)
            p => divvel(i2, :, i1)
         case (zdim)
            p => divvel(i1, i2, :)
      end select

   end subroutine set_div_v1d

!>
!! \brief Compute divergence of velocity
!!
!! \details This routine requires a single layer of valid guardcells uin cg%u arrays
!<

   subroutine div_v(ifluid, cg)

      use constants,   only: xdim, ydim, zdim, half
      use cr_data,     only: divv_n
      use dataio_pub,  only: die
      use domain,      only: dom, is_multicg
      use fluidindex,  only: iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz
      use grid_cont,   only: grid_container

      implicit none

      integer(kind=4), intent(in) :: ifluid
      type(grid_container), pointer, intent(inout) :: cg

      real, dimension(:), allocatable :: vx, vy, vz
      integer                         :: i, j, k
      integer                         :: idnf, imxf, imyf, imzf
      real, dimension(:,:,:), pointer :: divvel

      divvel => cg%ptr(divv_n)
      if (.not. associated(divvel)) call die("[crhelpers:div_v] cannot get divvel")

      if (any([allocated(vx), allocated(vy), allocated(vz)])) call die("[crhelpers:div_v] v[xyz] already allocated")
      allocate(vx(cg%n_(xdim)), vy(cg%n_(ydim)), vz(cg%n_(zdim)))

      idnf = iarr_all_dn(ifluid)
      imxf = iarr_all_mx(ifluid)
      imyf = iarr_all_my(ifluid)
      imzf = iarr_all_mz(ifluid)

      divvel(:,:,:) = 0.0

      if (dom%has_dir(xdim)) then
         do k = 1, cg%n_(zdim)
            do j = 1, cg%n_(ydim)
               vx = cg%u(imxf,:,j,k) / cg%u(idnf,:,j,k)
               divvel(2:cg%n_(xdim)-1,j,k) = ( vx(3:cg%n_(xdim)) - vx(1:cg%n_(xdim)-2) )  * (half*cg%idx)
            enddo
         enddo
         divvel(1,:,:) = divvel(2,:,:); divvel(cg%n_(xdim),:,:) = divvel(cg%n_(xdim)-1,:,:) ! for sanity
      endif

      if (dom%has_dir(ydim)) then
         do k = 1, cg%n_(zdim)
            do i = 1, cg%n_(xdim)
               vy = cg%u(imyf,i,:,k) / cg%u(idnf,i,:,k)
               divvel(i,2:cg%n_(ydim)-1,k) = divvel(i,2:cg%n_(ydim)-1,k)+( vy(3:cg%n_(ydim)) - vy(1:cg%n_(ydim)-2) )  * (half*cg%idy)
            enddo
         enddo
         divvel(:,1,:) = divvel(:,2,:); divvel(:, cg%n_(ydim),:) = divvel(:, cg%n_(ydim)-1,:) ! for sanity
      endif

      if (dom%has_dir(zdim)) then
         do j = 1, cg%n_(ydim)
            do i = 1, cg%n_(xdim)
               vz = cg%u(imzf,i,j,:) / cg%u(idnf,i,j,:)
               divvel(i,j,2:cg%n_(zdim)-1) = divvel(i,j,2:cg%n_(zdim)-1)+( vz(3:cg%n_(zdim)) - vz(1:cg%n_(zdim)-2) )  * (half*cg%idz)
            enddo
         enddo
         divvel(:,:,1) = divvel(:,:,2); divvel(:,:, cg%n_(zdim)) = divvel(:,:, cg%n_(zdim)-1) ! for sanity
      endif

      deallocate(vx, vy, vz)

   end subroutine div_v

end module crhelpers
