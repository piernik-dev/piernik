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

      use constants,     only: xdim, ydim, zdim, LO, HI, ndims, varlen, fluid_n, mag_n, wa_n
      use dataio_pub,    only: die
      use domain,        only: dom
      use fluidindex,    only: flind
      use global,        only: dt
      use grid,          only: leaves
      use gc_list,       only: cg_list_element
      use grid_cont,     only: grid_container
      use magboundaries, only: bnd_emf
      use rtvd,          only: tvdb

      implicit none

      integer(kind=4),       intent(in) :: bdir, vdir
      integer, pointer                  :: i1, i2, i1m, i2m
      integer, dimension(ndims), target :: ii, im
      integer                           :: rdir, i, j, i_wa, ui, bi
      integer(kind=4)                   :: imom                   !< index of vdir momentum
      real, dimension(:), allocatable   :: vv, vv0 !< \todo workaround for bug in gcc-4.6, REMOVE ME
      real, dimension(:),    pointer    :: pm1, pm2, pd1, pd2
      type(cg_list_element), pointer    :: cgl
      type(grid_container),  pointer    :: cg
      character(len=varlen)             :: emf
      character(len=varlen), dimension(ndims,ndims) :: v_b_          !< \deprecated v_b_ should be moved to a init routine

      v_b_(xdim,:) = [ 'vxbx', 'vxby', 'vxbz' ]
      v_b_(ydim,:) = [ 'vybx', 'vyby', 'vybz' ]
      v_b_(zdim,:) = [ 'vzbx', 'vzby', 'vzbz' ]

      emf = v_b_(vdir, bdir)

      imom = flind%ion%idn + int(vdir, kind=4)
      rdir = sum([xdim,ydim,zdim]) - bdir - vdir
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
         i_wa = cg%ind(wa_n)
         ui   = cg%ind_4d(fluid_n)
         bi   = cg%ind_4d(mag_n)

         if (any([allocated(vv), allocated(vv0)])) call die("[advects:advectb] vv or vv0 already allocated")
         allocate(vv(cg%n_(vdir)), vv0(cg%n_(vdir)))

         im(bdir) = 1
         do i = 1 + dom%D_(bdir), cg%n_(bdir)
            ii(bdir) = i
            do j = cg%ijkse(rdir,LO), cg%ijkse(rdir,HI)
               ii(rdir) = j
               im(rdir) = ii(rdir)
               vv=0.0
               pm1 => cg%w(ui)%get_sweep(vdir,imom,i1m,i2m)
               pm2 => cg%w(ui)%get_sweep(vdir,imom,i1 ,i2 )
               pd1 => cg%w(ui)%get_sweep(vdir,flind%ion%idn,i1m,i2m)
               pd2 => cg%w(ui)%get_sweep(vdir,flind%ion%idn,i1 ,i2 )
               vv0 = (pm1+pm2)/(pd1+pd2) !< \todo workaround for bug in gcc-4.6, REMOVE ME
               !vv(2:cg%n_(vdir)-1)=(vv(1:cg%n_(vdir)-2) + vv(3:cg%n_(vdir)) + 2.0*vv(2:cg%n_(vdir)-1))*0.25
               vv(2:cg%n_(vdir)-1)=(vv0(1:cg%n_(vdir)-2) + vv0(3:cg%n_(vdir)) + 2.0*vv0(2:cg%n_(vdir)-1))*0.25 !< \todo workaround for bug in gcc-4.6, REMOVE ME
               vv(1)  = vv(2)
               vv(cg%n_(vdir)) = vv(cg%n_(vdir)-1)

               vibj => cg%q(i_wa)%get_sweep(vdir,i1,i2)
               call tvdb(vibj, cg%w(bi)%get_sweep(vdir,bdir,i1,i2), vv, cg%n_(vdir),dt, cg%idl(vdir))
               NULLIFY(pm1, pm2, pd1, pd2)

            enddo
            im(bdir) = ii(bdir)
         enddo

         do i = xdim, zdim
            if (dom%has_dir(i)) call bnd_emf(cg%wa,emf,i, cg)
         enddo

         deallocate(vv, vv0)

         cgl => cgl%nxt
      enddo
      NULLIFY(i1, i1m, i2, i2m)

   end subroutine advectb

end module advects
