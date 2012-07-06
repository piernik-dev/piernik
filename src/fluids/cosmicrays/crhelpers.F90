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
   use constants, only: dsetnamelen
   implicit none

   private
   public :: div_v, set_div_v1d, divv_n

   character(len=dsetnamelen), parameter :: divv_n = "divvel" !< divergence of velocity

contains

   subroutine set_div_v1d(p, dir, i1, i2, cg)

      use dataio_pub,  only: die
      use grid_cont,   only: grid_container
      use named_array, only: qna

      implicit none

      integer(kind=4), intent(in) :: dir
      integer, intent(in) :: i1, i2
      real, dimension(:), pointer, intent(inout) :: p
      type(grid_container), pointer, intent(in) :: cg

      if (.not. qna%exists(divv_n)) call die("[crhelpers:set_div_v1d] cannot get divvel")
      p => cg%q(qna%ind(divv_n))%get_sweep(dir, i1, i2)

   end subroutine set_div_v1d

!>
!! \brief Compute divergence of velocity
!!
!! The divergence of velocity computed with the aid of 6-th order finite
!! differencing based on the Lagendre Polynomial interpolation.
!!
!! \details This routine requires a three layers of valid guardcells uin cg%u arrays
!<

! Should be moved in future to a dedicated module containing general purpose interpolation
! and derivation routines, and placed together with the other useful scheems described
! in particular in http://turbulence.pha.jhu.edu/Database-functions.pdf
!


   subroutine div_v(ifluid, cg)

      use constants,   only: xdim, ydim, zdim, big
      use dataio_pub,  only: die
      use domain,      only: dom
      use fluidindex,  only: iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz
      use grid_cont,   only: grid_container
      use named_array, only: qna

      implicit none

      integer(kind=4), intent(in) :: ifluid
      type(grid_container), pointer, intent(inout) :: cg

      real, dimension(:), allocatable :: vx, vy, vz
      integer                         :: i, j, k
      integer                         :: idnf, imxf, imyf, imzf
      real, dimension(:,:,:), pointer :: divvel
      integer                         :: nx, ny, nz
      real                            :: idx, idy, idz
      real, parameter                 :: p3_4 = 3./4., m3_20 = -3./20., p1_60 = 1./60.

      divvel => cg%q(qna%ind(divv_n))%arr
      if (.not. associated(divvel)) call die("[crhelpers:div_v] cannot get divvel")

      if (any([allocated(vx), allocated(vy), allocated(vz)])) call die("[crhelpers:div_v] v[xyz] already allocated")
      allocate(vx(cg%n_(xdim)), vy(cg%n_(ydim)), vz(cg%n_(zdim)))

      idnf = iarr_all_dn(ifluid)
      imxf = iarr_all_mx(ifluid)
      imyf = iarr_all_my(ifluid)
      imzf = iarr_all_mz(ifluid)
      nx   = cg%n_(xdim)  ; idx = cg%idx
      ny   = cg%n_(ydim)  ; idy = cg%idy
      nz   = cg%n_(zdim)  ; idy = cg%idy

      divvel(:,:,:) = 0.0

      if (dom%has_dir(xdim)) then
         do k = 1, cg%n_(zdim)
            do j = 1, cg%n_(ydim)
               vx = cg%u(imxf,:,j,k) / cg%u(idnf,:,j,k)
               divvel(4:nx-3,j,k) = divvel(4:nx-3,j,k) + ( vx(5:nx) - vx(3:nx-4) )  * (p3_4 *idx)
               divvel(4:nx-3,j,k) = divvel(4:nx-3,j,k) + ( vx(6:nx) - vx(2:nx-5) )  * (m3_20*idx)
               divvel(4:nx-3,j,k) = divvel(4:nx-3,j,k) + ( vx(7:nx) - vx(1:nx-6) )  * (p1_60*idx)
            enddo
         enddo
         divvel(1:3,:,:) = big; divvel(nx-2:nx,:,:) = big ! for sanity
      endif

      if (dom%has_dir(ydim)) then
         do k = 1, cg%n_(zdim)
            do i = 1, cg%n_(xdim)
               vy = cg%u(imyf,i,:,k) / cg%u(idnf,i,:,k)
               divvel(i,4:ny-3,k) = divvel(i,4:ny-3,k) + ( vy(5:ny) - vy(3:ny-4) )  * (p3_4 *idy)
               divvel(i,4:ny-3,k) = divvel(i,4:ny-3,k) + ( vy(6:ny) - vy(2:ny-5) )  * (m3_20*idy)
               divvel(i,4:ny-3,k) = divvel(i,4:ny-3,k) + ( vy(7:ny) - vy(1:ny-6) )  * (p1_60*idy)
            enddo
         enddo
         divvel(:,1:3,:) = big; divvel(:,ny-2:ny,:) = big ! for sanity
      endif

      if (dom%has_dir(zdim)) then
         do j = 1, cg%n_(ydim)
            do i = 1, cg%n_(xdim)
               vz = cg%u(imzf,i,j,:) / cg%u(idnf,i,j,:)
               divvel(i,j,4:nz-3) = divvel(i,j,4:nz-3) + ( vz(5:nz) - vz(3:nz-4) )  * (p3_4 *idz)
               divvel(i,j,4:nz-3) = divvel(i,j,4:nz-3) + ( vz(6:nz) - vz(2:nz-5) )  * (m3_20*idz)
               divvel(i,j,4:nz-3) = divvel(i,j,4:nz-3) + ( vz(7:nz) - vz(1:nz-6) )  * (p1_60*idz)
            enddo
         enddo
         divvel(:,:,1:3) = big; divvel(:,:,nz-2:nz) = big ! for sanity
      endif

      deallocate(vx, vy, vz)

   end subroutine div_v

end module crhelpers
