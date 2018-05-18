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
module ct
! pulled by MAGNETIC && RIEMANN
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

!---------------------------------------------------------------------------------------------------------------------

  subroutine ctb(vibj, b, vg, n, dt, idi)

    use constants, only: big, half

    implicit none

    integer(kind=4),               intent(in)    :: n       !< array size
      real,                        intent(in)    :: dt      !< time step
      real,                        intent(in)    :: idi     !< cell length, depends on direction x, y or z
      real, dimension(:), pointer, intent(inout) :: vibj    !< face-centered electromotive force components (b*vg)
      real, dimension(:), pointer, intent(in)    :: b       !< magnetic field
      real, dimension(n),          intent(in)    :: vg      !< velocity in the center of cell boundary
! locals
      real, dimension(n)                         :: b1      !< magnetic field
      real, dimension(n)                         :: vibj1   !< face-centered electromotive force (EMF) components (b*vg)
      real, dimension(n)                         :: vh      !< velocity interpolated to the cell edges
      real                                       :: dti     !< dt/di
      real                                       :: v       !< auxiliary variable to compute EMF
      real                                       :: w       !< EMF component
      real                                       :: dw      !< The second-order correction to EMF component
      real                                       :: dwm     !< face centered EMF interpolated to left cell-edge
      real                                       :: dwp     !< face centered EMF interpolated to right cell-edge
      integer                                    :: i       !< auxiliary array indicator
      integer                                    :: ip      !< i+1
      integer                                    :: ipp     !< i+2
      integer                                    :: im      !< i-1

! unlike the B field, the vibj lives on the right cell boundary
      vh = 0.0

! velocity interpolation to the cell boundaries

      vh(1:n-1) =(vg(1:n-1)+ vg(2:n))*half;     vh(n) = vh(n-1)

      dti = dt*idi

! face-centered EMF components computation, depending on the sign of vh, the components are upwinded  to cell edges, leading to 1st order EMF

      where (vh > 0.)
         vibj1=b*vg
      elsewhere
         vibj1=eoshift(b*vg,1,boundary=big)
      endwhere

! values of magnetic field computation in Runge-Kutta half step

      ! TODO: GEOFACTOR missing
      b1(2:n) = b(2:n) - (vibj1(2:n) - vibj1(1:n-1)) * dti * half;    b1(1) = b(2)

      do i = 3, n-3
         ip  = i  + 1
         ipp = ip + 1
         im  = i  - 1
         v   = vh(i)

! recomputation of EMF components (w) with b1 and face centered EMF interpolation to cell-edges (dwp, dwm), depending on the sign of v.

         if (v > 0.0) then
            w   = vg(i) * b1(i)
            dwp = (vg(ip) * b1(ip) - w) * half
            dwm = (w - vg(im) * b1(im)) * half
         else
            w   = vg(ip) * b1(ip)
            dwp = (w - vg(ipp) * b1(ipp)) * half
            dwm = (vg(i) * b1(i) - w) * half
         endif

! the second-order corrections to the EMF components computation with the aid of the van Leer monotonic interpolation and 2nd order EMF computation

         dw=0.0
         if (dwm * dwp > 0.0) dw = 2.0 * dwm * dwp / (dwm + dwp)
         vibj(i) = (w + dw) * dt
      enddo
      return

    end subroutine ctb

!---------------------------------------------------------------------------------------------------------------------
   subroutine advectb(bdir, vdir)

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: xdim, ydim, zdim, LO, HI, ndims, INT4
      use dataio_pub,       only: die
      use domain,           only: dom
      use fluidindex,       only: flind
      use global,           only: dt
      use grid_cont,        only: grid_container
      use magboundaries,    only: bnd_emf
      use named_array_list, only: qna, wna
      !use rtvd,             only: tvdb

      implicit none

      integer(kind=4),       intent(in) :: bdir, vdir
      integer, pointer                  :: i1, i2, i1m, i2m
      integer, dimension(ndims), target :: ii, im
      integer                           :: rdir, i, j
      integer(kind=4)                   :: imom                   !< index of vdir momentum
      integer(kind=4)                   :: dir
      integer(kind=4), dimension(ndims) :: emf
      real, dimension(:), allocatable   :: vv, vv0 !< \todo workaround for bug in gcc-4.6, REMOVE ME
      real, dimension(:),    pointer    :: pm1, pm2, pd1, pd2
      type(cg_list_element), pointer    :: cgl
      type(grid_container),  pointer    :: cg

      imom = flind%ion%idn + int(vdir, kind=4)
      rdir = sum([xdim,ydim,zdim]) - bdir - vdir
      emf(vdir) = 1_INT4 ; emf(bdir) = 2_INT4 ; emf(rdir) = 3_INT4

      if (mod(3+vdir-bdir,3) == 1) then     ! even permutation
         i1 => ii(rdir) ; i1m => ii(rdir) ; i2 => ii(bdir) ; i2m => im(bdir)
      elseif (mod(3+vdir-bdir,3) == 2) then !  odd permutation
         i1 => ii(bdir) ; i1m => im(bdir) ; i2 => ii(rdir) ; i2m => ii(rdir)
      else
         call die('[advects:advectb] neither even nor odd permutation.')
         i1 => ii(rdir) ; i1m => ii(rdir) ; i2 => ii(rdir) ; i2m => ii(rdir) ! suppress compiler warnings
      endif

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         if (any([allocated(vv), allocated(vv0)])) call die("[advects:advectb] vv or vv0 already allocated")
         allocate(vv(cg%n_(vdir)), vv0(cg%n_(vdir)))

         im(bdir) = cg%lhn(bdir, LO)
         do i = cg%lhn(bdir, LO) + dom%D_(bdir), cg%lhn(bdir, HI)
            ii(bdir) = i
            do j = cg%ijkse(rdir,LO), cg%ijkse(rdir,HI)
               ii(rdir) = j
               im(rdir) = ii(rdir)
               vv=0.0
               pm1 => cg%w(wna%fi)%get_sweep(vdir,imom,i1m,i2m)
               pm2 => cg%w(wna%fi)%get_sweep(vdir,imom,i1 ,i2 )
               pd1 => cg%w(wna%fi)%get_sweep(vdir,flind%ion%idn,i1m,i2m)
               pd2 => cg%w(wna%fi)%get_sweep(vdir,flind%ion%idn,i1 ,i2 )
               vv0(:) = (pm1+pm2)/(pd1+pd2) !< \todo workaround for bug in gcc-4.6, REMOVE ME
               !vv(2:cg%n_(vdir)-1)=(vv(1:cg%n_(vdir)-2) + vv(3:cg%n_(vdir)) + 2.0*vv(2:cg%n_(vdir)-1))*0.25
               vv(2:cg%n_(vdir)-1)=(vv0(1:cg%n_(vdir)-2) + vv0(3:cg%n_(vdir)) + 2.0*vv0(2:cg%n_(vdir)-1))*0.25 !< \todo workaround for bug in gcc-4.6, REMOVE ME
               vv(1)  = vv(2)
               vv(cg%n_(vdir)) = vv(cg%n_(vdir)-1)

               vibj => cg%q(qna%wai)%get_sweep(vdir,i1,i2)
               call ctb(vibj, cg%w(wna%bi)%get_sweep(vdir,bdir,i1,i2), vv, cg%n_(vdir),dt, cg%idl(vdir))
               NULLIFY(pm1, pm2, pd1, pd2)

            enddo
            im(bdir) = ii(bdir)
         enddo

         do dir = xdim, zdim
            if (dom%has_dir(dir)) call bnd_emf(qna%wai,emf(dir),dir, cg)
         enddo

         deallocate(vv, vv0)

         cgl => cgl%nxt
      enddo
      NULLIFY(i1, i1m, i2, i2m)

   end subroutine advectb

 end module ct

