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
!! \brief Routines for converting face-centered magnetic field to cell-centered values.
!<

module bfc_bcc
! pulled by ANY

   implicit none

   private
   public :: interpolate_mag_field

contains

!>
!! \brief return 1D vector of magnetic field components from cg in direction cdim at indices i1 and i2.
!<

   function interpolate_mag_field(cdim, cg, i1, i2, ind) result (b)

      use constants,        only: pdims, xdim, ydim, zdim, half, ORTHO1, ORTHO2
      use domain,           only: dom
      use fluidindex,       only: iarr_mag_swp, nmag
      use grid_cont,        only: grid_container

      implicit none

      integer(kind=4),               intent(in) :: cdim
      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: i1, i2
      integer(kind=4),               intent(in) :: ind

      real, dimension(cg%n_(cdim), nmag)        :: b
      real, dimension(:), pointer               :: pb, pb1
      integer(kind=4)                           :: ibx, iby, ibz
      integer                                   :: i1p, i2p

      !> OPTIMIZE ME

      ibx = iarr_mag_swp(cdim, xdim)
      iby = iarr_mag_swp(cdim, ydim)
      ibz = iarr_mag_swp(cdim, zdim)

      i1p = i1 + dom%D_(pdims(cdim, ORTHO1))
      i2p = i2 + dom%D_(pdims(cdim, ORTHO2))

      pb => cg%w(ind)%get_sweep(cdim, ibx, i1, i2)
      b(1:cg%n_(cdim) - 1, ibx) = half * (pb(1:cg%n_(cdim)-1) + pb(2:cg%n_(cdim)))
      b(cg%n_(cdim),       ibx) = b(cg%n_(cdim)-1, ibx)

      pb  => cg%w(ind)%get_sweep(cdim, iby, i1, i2)
      if (cdim == xdim) then
         pb1 => cg%w(ind)%get_sweep(cdim, iby, i1p, i2)
      else
         pb1 => cg%w(ind)%get_sweep(cdim, iby, i1, i2p)
      endif
      b(:, iby) = half * (pb + pb1)

      pb  => cg%w(ind)%get_sweep(cdim, ibz, i1, i2)
      if (cdim == xdim) then
         pb1 => cg%w(ind)%get_sweep(cdim, ibz, i1, i2p)
      else
         pb1 => cg%w(ind)%get_sweep(cdim, ibz, i1p, i2)
      endif
      b(:, ibz) = half * (pb + pb1)

      b(:, iarr_mag_swp(cdim,:)) = b(:,:)
      nullify(pb,pb1)

   end function interpolate_mag_field

end module bfc_bcc
