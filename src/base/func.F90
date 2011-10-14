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
!! \brief (KK) Module that contains unclassified functions
!!
!! This module should be empty. Every function or subroutine placed here belong
!! elsewhere. We are yet unsure where to put them.
!! \todo Move all structures elsewhere
!! \warning Procedures \a dipole and \a rn_angles were moved to sn_sources.F90
!<
module func

   implicit none

   private
   public :: ekin, emag, L2norm, get_extremum

contains

!> \todo move to a better place (possibly dataio_pub since stop types module stops using die routine)
   subroutine get_extremum(tab, minmax, prop, cg)

      use constants,  only: MINL, MAXL, I_ONE, ndims, xdim, ydim, zdim
      use dataio_pub, only: msg, warn, die
      use domain,     only: dom, is_multicg
      use grid_cont,  only: grid_container
      use mpi,        only: MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_STATUS_IGNORE, MPI_2DOUBLE_PRECISION, MPI_MINLOC, MPI_MAXLOC, MPI_IN_PLACE
      use mpisetup,   only: comm, ierr, master, proc, FIRST
      use types,      only: value

      implicit none

      real, dimension(:,:,:), intent(in), pointer  :: tab
      integer(kind=4),        intent(in)  :: minmax
      type(value),            intent(out) :: prop
      type(grid_container), pointer, intent(in) :: cg

      integer, parameter :: tag1 = 11
      integer, parameter :: tag2 = 12
      real, dimension(2)  :: v_red
      integer, dimension(MINL:MAXL), parameter :: op = [ MPI_MINLOC, MPI_MAXLOC ]

      if (is_multicg) call die("[func:get_extremum] multiple grid pieces per procesor not implemented yet") !nontrivial

      select case (minmax)
         case (MINL)
            prop%val = minval(tab)
            prop%loc = minloc(tab) + [cg%nb, cg%nb, cg%nb]
         case (MAXL)
            prop%val = maxval(tab)
            prop%loc = maxloc(tab) + [cg%nb, cg%nb, cg%nb]
         case default
            write(msg,*) "[func:get_extremum]: I don't know what to do with minmax = ", minmax
            call warn(msg)
      end select

      v_red(:) = [ prop%val, real(proc) ]

      call MPI_Allreduce(MPI_IN_PLACE, v_red, I_ONE, MPI_2DOUBLE_PRECISION, op(minmax), comm, ierr)

      prop%val = v_red(1)
      prop%proc = int(v_red(2))

      if (proc == prop%proc) then
         where (.not. dom%has_dir(:)) prop%coords(:) = 0.
         if (dom%has_dir(xdim)) prop%coords(xdim) = cg%x(prop%loc(xdim))
         if (dom%has_dir(ydim)) prop%coords(ydim) = cg%y(prop%loc(ydim))
         if (dom%has_dir(zdim)) prop%coords(zdim) = cg%z(prop%loc(zdim))
      endif

      if (prop%proc /= 0) then
         if (proc == prop%proc) then ! slave
            call MPI_Send (prop%loc,    ndims, MPI_INTEGER,          FIRST, tag1, comm, ierr)
            call MPI_Send (prop%coords, ndims, MPI_DOUBLE_PRECISION, FIRST, tag2, comm, ierr)
         endif
         if (master) then
            call MPI_Recv (prop%loc,    ndims, MPI_INTEGER,          prop%proc, tag1, comm, MPI_STATUS_IGNORE, ierr)
            call MPI_Recv (prop%coords, ndims, MPI_DOUBLE_PRECISION, prop%proc, tag2, comm, MPI_STATUS_IGNORE, ierr)
         endif
      endif

   end subroutine get_extremum

   elemental real function L2norm(x1,x2,x3,y1,y2,y3)
      implicit none
      real, intent(in) :: x1, x2, x3
      real, intent(in) :: y1, y2, y3

      L2norm = sqrt( (x1 - y1)**2 + (x2 - y2)**2 + (x3 - y3)**2 )
   end function L2norm

   elemental real function emag(bx,by,bz)
      use constants,  only: half
      implicit none
      real, intent(in) :: bx, by, bz

      emag = half*(bx**2 + by**2 + bz**2)

   end function emag

   elemental real function ekin(mx,my,mz,dn)
      implicit none
      real, intent(in) :: mx, my, mz, dn

      ekin = emag(mx,my,mz)/dn

   end function ekin

end module func
