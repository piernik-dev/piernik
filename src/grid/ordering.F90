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
!! \brief Functions that convert coordinates into single integer (position on a space-filling curve)
!!
!! \todo implement inverse functions: convert id to offsets
!!
!! \todo implement also Hilbert or Peano-Hilbert ordering
!<

module ordering

   implicit none

   private
   public :: SFC_order

contains

!>
!! \brief Wrapper for space-filling curve ID routine
!!
!! \details Currently only Morton_id is implemented, but in the future we may implement also Hilbert or Peano-Hilbert curve.
!! When we implement other SFC, change call to Morton_id with a call to pointer that is set to the right type of SFC.
!<

   function SFC_order(off) result(id)

      use constants,  only: ndims, xdim, ydim, zdim, INVALID
      use dataio_pub, only: die
      use domain,     only: dom

      implicit none

      integer(kind=8), dimension(ndims), intent(in) :: off

      integer(kind=8) :: id
      integer :: j1, j2

      id = INVALID
      select case (dom%eff_dim)
         case (0) ! Special case
            id = 1
         case (1) ! No need to process coordinate
            do j1 = xdim, zdim
               if (dom%has_dir(j1)) id = off(j1)
            enddo
         case (2) ! select coordinates only in existing dimensions
            if (dom%has_dir(xdim)) then
               j1 = xdim
            else
               j1 = ydim
            endif
            if (dom%has_dir(zdim)) then
               j2 = zdim
            else
               j2 = ydim
            endif
            id = Morton_id([ off(j1), off(j2) ])
         case (3) ! just do the conversion
            id = Morton_id(off)
         case default
            call die("[ordering:SFC_order] invalid dimensionality")
      end select

      if (id == INVALID .or. id < 0) call die("[ordering:SFC_order] invalid id")

   end function SFC_order

!>
!! \brief Convert contiguous vector of coordinates into its Morton identifier
!!
!! \details 64-bit integers are in use so offsets up to 2**21 in 3D and up to 2**31 in 2D should be supported.
!! For a 64**3 domain that corresponds to a 15 levels of refinement.
!! If more levels are required one can:
!! * use AMR::bsize divisible by some 2*k and ignore lowest k bits of offset. This way additional k levels of refinement can be supported.
!! * divide off by AMR::bsize to increase maximum effective resolution to AMR::bsize*2**21 in 3D
!! * use 128-bit integers (if they're available) to increase allowed range by 21 levels of refinement
!! * emulate 128-bit integer in two 64-bit integers (note that we don't need arithmetic here, just bitwise operations and comparison)
!<

   function Morton_id(off) result(id)

      use constants,  only: ndims
      use dataio_pub, only: die

      implicit none

      integer(kind=8), dimension(:), intent(in) :: off

      integer(kind=8) :: id
      integer(kind=8), dimension(ndims) :: o
      integer :: i
      integer(kind=8) :: mask

      if (any(off < 0)) call die("[ordering:Morton_id] only non-negative offsets are allowed")

      o(:size(off)) = off
      id = 0
      mask = 1
      do while (any(o(:size(off)) /= 0))
         do i = lbound(off, dim=1), ubound(off, dim=1)
            if (btest(o(i), 0)) then
               id = ior(id, mask)
               if (mask <= 0) call die("[ordering:Morton_id] mask overflow")
            endif
            mask = ishft(mask, 1)
         enddo
         o(:size(off)) = ishft(o(:size(off)), -1)
      enddo

      !id = ieor(id, ishft(id, -1)) Gray code

   end function Morton_id

end module ordering
