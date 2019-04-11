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
!! \brief Diagnostics of particle module
!!
!! This module contains particle diagnostics related things that are needed by I/O and debug routines.
!<

module particle_utils
! pulled by NBODY

   implicit none

   private
   public :: max_pvel_1d, max_pacc_3d

contains

   subroutine max_pvel_1d(max_v)

      use constants,      only: ndims, xdim, zdim
      use particle_types, only: pset

      implicit none

      real, dimension(ndims), intent(out) :: max_v
      integer(kind=4)                     :: cdim

      do cdim = xdim, zdim
         max_v(cdim)  = maxval(abs(pset%p(:)%vel(cdim)))
      enddo

   end subroutine max_pvel_1d

   subroutine max_pacc_3d(cg, max_pacc)

      use constants,      only: big, CENTER, half, LO, xdim, zdim, zero
      use grid_cont,      only: grid_container
      use mpisetup,       only: proc
      use particle_types, only: pset
      use types,          only: value

      implicit none

      type(grid_container), pointer, intent(in)  :: cg
      type(value),                   intent(out) :: max_pacc
      real                                       :: acc2, max_acc
      integer                                    :: i, n_part
      integer(kind=4) :: cdim

      max_pacc%assoc = big

      max_acc  = zero
      n_part = size(pset%p, dim=1)

      do i = 1, n_part
         acc2 = sum(pset%p(i)%acc(:)**2)
         if (acc2 > max_acc) then
            max_acc = acc2
            max_pacc%coords(:) = pset%p(i)%pos(:)
            !max_pacc%proc = i !> \todo it might be an information about extremum particle, but the scheme of log file is to print the process number
         endif
      enddo
      max_pacc%val = sqrt(max_acc)
      do cdim = xdim, zdim
         max_pacc%loc(cdim) = int( half + (max_pacc%coords(cdim) - cg%coord(CENTER,cdim)%r(cg%ijkse(cdim, LO))) * cg%idl(cdim) )
      enddo
      max_pacc%proc = proc

   end subroutine max_pacc_3d

end module particle_utils
