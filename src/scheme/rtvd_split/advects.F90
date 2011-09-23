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

!>
!! \todo remove workaround for http://gcc.gnu.org/bugzilla/show_bug.cgi?id=48955
!<
module advects
! pulled by MAGNETIC && RTVD
   implicit none

   private
   public  :: advectb
   real, dimension(:), pointer :: vibj => null()

contains

!>
!!   advectby_x --> advectb(bdir=ydim, vdir=xdim, emf='vxby')
!!   advectbz_x --> advectb(bdir=zdim, vdir=xdim, emf='vxbz')
!!   advectbx_y --> advectb(bdir=xdim, vdir=ydim, emf='vybx')
!!   advectbz_y --> advectb(bdir=zdim, vdir=ydim, emf='vybz')
!!   advectbx_z --> advectb(bdir=xdim, vdir=zdim, emf='vzbx')
!!   advectby_z --> advectb(bdir=ydim, vdir=zdim, emf='vzby')

!<
   subroutine advectb(bdir, vdir)

      use constants, only: xdim, ydim, zdim

      implicit none

      integer, intent(in) :: bdir, vdir

      select case ( bdir)
         case (xdim)
            if (vdir == ydim) call advectbx_y
            if (vdir == zdim) call advectbx_z
         case (ydim)
            if (vdir == xdim) call advectby_x
            if (vdir == zdim) call advectby_z
         case (zdim)
            if (vdir == xdim) call advectbz_x
            if (vdir == ydim) call advectbz_y
      end select

   end subroutine advectb

   subroutine advectby_x

      use constants,     only: xdim, ydim, zdim
      use dataio_pub,    only: die
      use domain,        only: has_dir, D_y
      use fluidindex,    only: iby, flind
      use global,        only: dt
      use grid,          only: all_cg
      use gc_list,       only: cg_list_element
      use grid_cont,     only: grid_container
      use magboundaries, only: bnd_emf
      use rtvd,          only: tvdb

      implicit none

      integer                         :: j, k, jm, j_s, j_e
      real, dimension(:), allocatable :: vx
      real, dimension(:), allocatable :: vx0 !< \todo workaround for bug in gcc-4.6, REMOVE ME
      type(cg_list_element), pointer  :: cgl
      type(grid_container),  pointer  :: cg
      integer :: i_wa

      cgl => all_cg%first
      do while (associated(cgl))
         cg => cgl%cg
         i_wa = cg%get_na_ind("wa") ! BEWARE: magic strings across multiple files

         if (any([allocated(vx), allocated(vx0)])) call die("[advects:advectby_x] vx or vx0 already allocated")
         allocate(vx(cg%n_(xdim)), vx0(cg%n_(xdim)))

         j_s = 1 + D_y
         j_e = cg%n_(ydim)

         do k = cg%ks, cg%ke
            do j = j_s, j_e         ! cg%n_(ydim)b cg%is /= 1 in by_x
               jm=j-1
               vx=0.0
               if (jm == 0) then
                  !vx = cg%u%arr(flind%ion%imx,:,1,k) / cg%u%arr(flind%ion%idn,:,1,k)
                  vx0 = cg%u%arr(flind%ion%imx,:,1,k) / cg%u%arr(flind%ion%idn,:,1,k) !< \todo workaround for bug in gcc-4.6, REMOVE ME
               else
                  !vx =(cg%u%arr(flind%ion%imx,:,jm,k)+cg%u%arr(flind%ion%imx,:,j,k))/(cg%u%arr(flind%ion%idn,:,jm,k)+cg%u%arr(flind%ion%idn,:,j,k))
                  vx0 =(cg%u%arr(flind%ion%imx,:,jm,k)+cg%u%arr(flind%ion%imx,:,j,k))/(cg%u%arr(flind%ion%idn,:,jm,k)+cg%u%arr(flind%ion%idn,:,j,k)) !< \todo workaround for bug in gcc-4.6, REMOVE ME
               endif
               !vx(2:cg%n_(xdim)-1)=(vx(1:cg%n_(xdim)-2) + vx(3:cg%n_(xdim)) + 2.0*vx(2:cg%n_(xdim)-1))*0.25
               vx(2:cg%n_(xdim)-1)=(vx0(1:cg%n_(xdim)-2) + vx0(3:cg%n_(xdim)) + 2.0*vx0(2:cg%n_(xdim)-1))*0.25 !< \todo workaround for bug in gcc-4.6, REMOVE ME
               vx(1)  = vx(2)
               vx(cg%n_(xdim)) = vx(cg%n_(xdim)-1)

               vibj => cg%q(i_wa)%get_sweep(xdim,j,k)
               call tvdb(vibj, cg%b%get_sweep(xdim,iby,j,k), vx, cg%n_(xdim),dt, cg%idx)

            enddo
         enddo

         do j = xdim, zdim
            if (has_dir(j)) call bnd_emf(cg%wa,'vxby',j)
         enddo

         deallocate(vx)
         deallocate(vx0)

         cgl => cgl%nxt
      enddo

   end subroutine advectby_x

   subroutine advectbz_x

      use constants,     only: xdim, zdim
      use dataio_pub,    only: die
      use domain,        only: has_dir, D_z
      use fluidindex,    only: ibz, flind
      use global,        only: dt
      use grid,          only: all_cg
      use gc_list,       only: cg_list_element
      use grid_cont,     only: grid_container
      use magboundaries, only: bnd_emf
      use rtvd,          only: tvdb

      implicit none

      integer                         :: j, k, km, k_s, k_e
      real, dimension(:), allocatable :: vx
      real, dimension(:), allocatable :: vx0 !< \todo workaround for bug in gcc-4.6, REMOVE ME
      type(cg_list_element), pointer  :: cgl
      type(grid_container),  pointer  :: cg
      integer :: i_wa

      cgl => all_cg%first
      do while (associated(cgl))
         cg => cgl%cg
         i_wa = cg%get_na_ind("wa")

         if (any([allocated(vx), allocated(vx0)])) call die("[advects:advectby_x] vx or vx0 already allocated")
         allocate(vx(cg%n_(xdim)), vx0(cg%n_(xdim)))

         k_s = 1 + D_z
         k_e = cg%n_(zdim)

         do k=k_s,k_e
            km=k-1
            do j=cg%js, cg%je
               vx=0.0
               if (km == 0) then
                  !vx = cg%u%arr(flind%ion%imx,:,j,1) / cg%u%arr(flind%ion%idn,:,j,1)
                  vx0 = cg%u%arr(flind%ion%imx,:,j,1) / cg%u%arr(flind%ion%idn,:,j,1) !< \todo workaround for bug in gcc-4.6, REMOVE ME
               else
                  !vx = (cg%u%arr(flind%ion%imx,:,j,km)+cg%u%arr(flind%ion%imx,:,j,k))/(cg%u%arr(flind%ion%idn,:,j,km)+cg%u%arr(flind%ion%idn,:,j,k))
                  vx0 = (cg%u%arr(flind%ion%imx,:,j,km)+cg%u%arr(flind%ion%imx,:,j,k))/(cg%u%arr(flind%ion%idn,:,j,km)+cg%u%arr(flind%ion%idn,:,j,k)) !< \todo workaround for bug in gcc-4.6, REMOVE ME
               endif
               !vx(2:cg%n_(xdim)-1)=(vx(1:cg%n_(xdim)-2) + vx(3:cg%n_(xdim)) + 2.0*vx(2:cg%n_(xdim)-1))*0.25
               vx(2:cg%n_(xdim)-1)=(vx0(1:cg%n_(xdim)-2) + vx0(3:cg%n_(xdim)) + 2.0*vx0(2:cg%n_(xdim)-1))*0.25 !< \todo workaround for bug in gcc-4.6, REMOVE ME
               vx(1)  = vx(2)
               vx(cg%n_(xdim)) = vx(cg%n_(xdim)-1)

               vibj => cg%q(i_wa)%get_sweep(xdim,j,k)
               call tvdb(vibj, cg%b%get_sweep(xdim,ibz,j,k), vx, cg%n_(xdim),dt, cg%idx)

            enddo
         enddo

         do j = xdim, zdim
            if (has_dir(j)) call bnd_emf(cg%wa,'vxbz',j)
         enddo

         deallocate(vx)
         deallocate(vx0)

        cgl => cgl%nxt
      enddo

   end subroutine advectbz_x

   subroutine advectbz_y

      use constants,     only: xdim, ydim, zdim
      use dataio_pub,    only: die
      use domain,        only: has_dir, D_z
      use fluidindex,    only: ibz, flind
      use global,        only: dt
      use grid,          only: all_cg
      use gc_list,       only: cg_list_element
      use grid_cont,     only: grid_container
      use magboundaries, only: bnd_emf
      use rtvd,          only: tvdb

      implicit none

      integer                         :: i, k, km, k_s, k_e
      real, dimension(:), allocatable :: vy
      real, dimension(:), allocatable :: vy0 !< \todo workaround for bug in gcc-4.6, REMOVE ME
      type(cg_list_element), pointer  :: cgl
      type(grid_container),  pointer  :: cg
      integer :: i_wa

      cgl => all_cg%first
      do while (associated(cgl))
         cg => cgl%cg
         i_wa = cg%get_na_ind("wa")

         if (any([allocated(vy), allocated(vy0)])) call die("[advects:advectby_y] vy or vy0 already allocated")
         allocate(vy(cg%n_(ydim)), vy0(cg%n_(ydim)))

         k_s = 1 + D_z
         k_e = cg%n_(zdim)

         do k=k_s,k_e
            km=k-1
            do i=cg%is, cg%ie
               vy=0.0
               if (km == 0) then
                  !vy = cg%u%arr(flind%ion%imy,i,:,1) / cg%u%arr(flind%ion%idn,i,:,1)
                  vy0 = cg%u%arr(flind%ion%imy,i,:,1) / cg%u%arr(flind%ion%idn,i,:,1) !< \todo workaround for bug in gcc-4.6, REMOVE ME
               else
                  !vy=(cg%u%arr(flind%ion%imy,i,:,km)+cg%u%arr(flind%ion%imy,i,:,k))/(cg%u%arr(flind%ion%idn,i,:,km)+cg%u%arr(flind%ion%idn,i,:,k))
                  vy0=(cg%u%arr(flind%ion%imy,i,:,km)+cg%u%arr(flind%ion%imy,i,:,k))/(cg%u%arr(flind%ion%idn,i,:,km)+cg%u%arr(flind%ion%idn,i,:,k)) !< \todo workaround for bug in gcc-4.6, REMOVE ME
               endif
               !vy(2:cg%n_(ydim)-1)=(vy(1:cg%n_(ydim)-2) + vy(3:cg%n_(ydim)) + 2.0*vy(2:cg%n_(ydim)-1))*0.25
               vy(2:cg%n_(ydim)-1)=(vy0(1:cg%n_(ydim)-2) + vy0(3:cg%n_(ydim)) + 2.0*vy0(2:cg%n_(ydim)-1))*0.25 !< \todo workaround for bug in gcc-4.6, REMOVE ME
               vy(1)  = vy(2)
               vy(cg%n_(ydim)) = vy(cg%n_(ydim)-1)

               vibj => cg%q(i_wa)%get_sweep(ydim,k,i)
               call tvdb(vibj, cg%b%get_sweep(ydim,ibz,k,i), vy, cg%n_(ydim),dt, cg%idy)

            enddo
         enddo

         do i = xdim, zdim
            if (has_dir(i)) call bnd_emf(cg%wa,'vybz',i)
         enddo

         deallocate(vy)
         deallocate(vy0)

         cgl => cgl%nxt
      enddo

   end subroutine advectbz_y

   subroutine advectbx_y

      use constants,     only: xdim, ydim, zdim
      use dataio_pub,    only: die
      use domain,        only: has_dir, D_x
      use fluidindex,    only: ibx, flind
      use global,        only: dt
      use grid,          only: all_cg
      use gc_list,       only: cg_list_element
      use grid_cont,     only: grid_container
      use magboundaries, only: bnd_emf
      use rtvd,          only: tvdb

      implicit none

      integer                         :: k, i, im, i_s, i_e
      real, dimension(:), allocatable :: vy
      real, dimension(:), allocatable :: vy0 !< \todo workaround for bug in gcc-4.6, REMOVE ME
      type(cg_list_element), pointer  :: cgl
      type(grid_container),  pointer  :: cg
      integer :: i_wa

      cgl => all_cg%first
      do while (associated(cgl))
         cg => cgl%cg
         i_wa = cg%get_na_ind("wa")

         if (any([allocated(vy), allocated(vy0)])) call die("[advects:advectby_y] vy or vy0 already allocated")
         allocate(vy(cg%n_(ydim)), vy0(cg%n_(ydim)))

         i_s = 1 + D_x
         i_e = cg%n_(xdim)

         do k=cg%ks, cg%ke
            do i=i_s,i_e
               im=i-1
               vy=0.0
               if (im == 0) then
                  !vy = cg%u%arr(flind%ion%imy,1,:,k) / cg%u%arr(flind%ion%idn,1,:,k)
                  vy0 = cg%u%arr(flind%ion%imy,1,:,k) / cg%u%arr(flind%ion%idn,1,:,k) !< \todo workaround for bug in gcc-4.6, REMOVE ME
               else
                  !vy=(cg%u%arr(flind%ion%imy,im,:,k)+cg%u%arr(flind%ion%imy,i,:,k))/(cg%u%arr(flind%ion%idn,im,:,k)+cg%u%arr(flind%ion%idn,i,:,k))
                  vy0=(cg%u%arr(flind%ion%imy,im,:,k)+cg%u%arr(flind%ion%imy,i,:,k))/(cg%u%arr(flind%ion%idn,im,:,k)+cg%u%arr(flind%ion%idn,i,:,k)) !< \todo workaround for bug in gcc-4.6, REMOVE ME
               endif
               !vy(2:cg%n_(ydim)-1)=(vy(1:cg%n_(ydim)-2) + vy(3:cg%n_(ydim)) + 2.0*vy(2:cg%n_(ydim)-1))*0.25
               vy(2:cg%n_(ydim)-1)=(vy0(1:cg%n_(ydim)-2) + vy0(3:cg%n_(ydim)) + 2.0*vy0(2:cg%n_(ydim)-1))*0.25 !< \todo workaround for bug in gcc-4.6, REMOVE ME
               vy(1)  = vy(2)
               vy(cg%n_(ydim)) = vy(cg%n_(ydim)-1)

               vibj => cg%q(i_wa)%get_sweep(ydim,k,i)
               call tvdb(vibj, cg%b%get_sweep(ydim,ibx,k,i), vy, cg%n_(ydim),dt, cg%idy)

            enddo
         enddo

         do i = xdim, zdim
            if (has_dir(i)) call bnd_emf(cg%wa,'vybx',i)
         enddo

         deallocate(vy)
         deallocate(vy0)

         cgl => cgl%nxt
      enddo

   end subroutine advectbx_y

   subroutine advectbx_z

      use constants,     only: xdim, zdim
      use dataio_pub,    only: die
      use domain,        only: has_dir, D_x
      use fluidindex,    only: ibx, flind
      use global,        only: dt
      use grid,          only: all_cg
      use gc_list,       only: cg_list_element
      use grid_cont,     only: grid_container
      use magboundaries, only: bnd_emf
      use rtvd,          only: tvdb

      implicit none

      integer                         :: j, i, im, i_s, i_e
      real, dimension(:), allocatable :: vz
      real, dimension(:), allocatable :: vz0 !< \todo workaround for bug in gcc-4.6, REMOVE ME
      type(cg_list_element), pointer  :: cgl
      type(grid_container),  pointer  :: cg
      integer :: i_wa

      cgl => all_cg%first
      do while (associated(cgl))
         cg => cgl%cg
         i_wa = cg%get_na_ind("wa")

         if (any([allocated(vz), allocated(vz0)])) call die("[advects:advectby_z] vz or vz0 already allocated")
         allocate(vz(cg%n_(zdim)), vz0(cg%n_(zdim)))

         i_s = 1 + D_x
         i_e = cg%n_(xdim)

         do j=cg%js, cg%je
            do i=i_s,i_e
               im=i-1
               vz=0.0
               if (im == 0) then
                  !vz = cg%u%arr(flind%ion%imz,1,j,:) / cg%u%arr(flind%ion%idn,1,j,:)
                  vz0 = cg%u%arr(flind%ion%imz,1,j,:) / cg%u%arr(flind%ion%idn,1,j,:) !< \todo workaround for bug in gcc-4.6, REMOVE ME
               else
                  !vz = (cg%u%arr(flind%ion%imz,im,j,:)+cg%u%arr(flind%ion%imz,i,j,:))/(cg%u%arr(flind%ion%idn,im,j,:)+cg%u%arr(flind%ion%idn,i,j,:))
                  vz0 = (cg%u%arr(flind%ion%imz,im,j,:)+cg%u%arr(flind%ion%imz,i,j,:))/(cg%u%arr(flind%ion%idn,im,j,:)+cg%u%arr(flind%ion%idn,i,j,:)) !< \todo workaround for bug in gcc-4.6, REMOVE ME
               endif
               !vz(2:cg%n_(zdim)-1)=(vz(1:cg%n_(zdim)-2) + vz(3:cg%n_(zdim)) + 2.0*vz(2:cg%n_(zdim)-1))*0.25
               vz(2:cg%n_(zdim)-1)=(vz0(1:cg%n_(zdim)-2) + vz0(3:cg%n_(zdim)) + 2.0*vz0(2:cg%n_(zdim)-1))*0.25 !< \todo workaround for bug in gcc-4.6, REMOVE ME
               vz(1)  = vz(2)
               vz(cg%n_(zdim)) = vz(cg%n_(zdim)-1)

               vibj => cg%q(i_wa)%get_sweep(zdim,i,j)
               call tvdb(vibj, cg%b%get_sweep(zdim,ibx,i,j), vz, cg%n_(zdim),dt, cg%idz)

            enddo
         enddo

         do i = xdim, zdim
            if (has_dir(i)) call bnd_emf(cg%wa,'vzbx',i)
         enddo

         deallocate(vz)
         deallocate(vz0)

         cgl => cgl%nxt
      enddo

   end subroutine advectbx_z

   subroutine advectby_z

      use constants,     only: xdim, ydim, zdim
      use dataio_pub,    only: die
      use domain,        only: has_dir, D_y
      use fluidindex,    only: iby, flind
      use global,        only: dt
      use grid,          only: all_cg
      use gc_list,       only: cg_list_element
      use grid_cont,     only: grid_container
      use magboundaries, only: bnd_emf
      use rtvd,          only: tvdb

      implicit none

      integer                         :: i, j, jm, j_s, j_e
      real, dimension(:), allocatable :: vz
      real, dimension(:), allocatable :: vz0 !< \todo workaround for bug in gcc-4.6, REMOVE ME
      type(cg_list_element), pointer  :: cgl
      type(grid_container),  pointer  :: cg
      integer :: i_wa

      cgl => all_cg%first
      do while (associated(cgl))
         cg => cgl%cg
         i_wa = cg%get_na_ind("wa")

         if (any([allocated(vz), allocated(vz0)])) call die("[advects:advectby_z] vz or vz0 already allocated")
         allocate(vz(cg%n_(zdim)), vz0(cg%n_(zdim)))

         j_s = 1 + D_y
         j_e = cg%n_(ydim)

         do j=j_s,j_e
            jm=j-1
            do i=cg%is, cg%ie
               vz=0.0
               if (jm == 0) then
                  !vz = cg%u%arr(flind%ion%imz,i,1,:) / cg%u%arr(flind%ion%idn,i,1,:)
                  vz0 = cg%u%arr(flind%ion%imz,i,1,:) / cg%u%arr(flind%ion%idn,i,1,:) !< \todo workaround for bug in gcc-4.6, REMOVE ME
               else
                  !vz=(cg%u%arr(flind%ion%imz,i,jm,:)+cg%u%arr(flind%ion%imz,i,j,:))/(cg%u%arr(flind%ion%idn,i,jm,:)+cg%u%arr(flind%ion%idn,i,j,:))
                  vz0=(cg%u%arr(flind%ion%imz,i,jm,:)+cg%u%arr(flind%ion%imz,i,j,:))/(cg%u%arr(flind%ion%idn,i,jm,:)+cg%u%arr(flind%ion%idn,i,j,:)) !< \todo workaround for bug in gcc-4.6, REMOVE ME
               endif
               !vz(2:cg%n_(zdim)-1)=(vz(1:cg%n_(zdim)-2) + vz(3:cg%n_(zdim)) + 2.0*vz(2:cg%n_(zdim)-1))*0.25
               vz(2:cg%n_(zdim)-1)=(vz0(1:cg%n_(zdim)-2) + vz0(3:cg%n_(zdim)) + 2.0*vz0(2:cg%n_(zdim)-1))*0.25 !< \todo workaround for bug in gcc-4.6, REMOVE ME
               vz(1)  = vz(2)
               vz(cg%n_(zdim)) = vz(cg%n_(zdim)-1)

               vibj => cg%q(i_wa)%get_sweep(zdim,i,j)
               call tvdb(vibj, cg%b%get_sweep(zdim,iby,i,j), vz, cg%n_(zdim),dt, cg%idz)

            enddo
         enddo

         do i = xdim, zdim
            if (has_dir(i)) call bnd_emf(cg%wa,'vzby',i)
         enddo

         deallocate(vz)
         deallocate(vz0)

         cgl => cgl%nxt
      enddo

   end subroutine advectby_z

end module advects
