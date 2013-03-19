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
!! \brief This module contains problem variables intended to use across the code.
!!
!! \details It is strongly discouraged to use this feature for purposes other than debugging.
!<

module problem_pub

#ifdef MACLAURIN_PROBLEM
   use constants, only: ndims
#endif /* MACLAURIN_PROBLEM */

   implicit none

   public  ! QA_WARN nothing to hide here

   ! hack for tests
#ifdef JEANS_PROBLEM
   real :: jeans_d0
   integer :: jeans_mode
#endif /* JEANS_PROBLEM */
#ifdef MACLAURIN_PROBLEM
   real, dimension(ndims) :: xs ! position of the sphere
   real :: as ! factor used to calculate the potential

contains

!> \brief Calculate point-like potential for the maclaurin problem

   real function ap_potential(x, y, z) result(phi)

      use constants, only: xdim, ydim, zdim

      implicit none

      real, intent(in) :: x, y, z

      phi = as / sqrt((x-xs(xdim))**2 + (y-xs(ydim))**2 + (z-xs(zdim))**2)

   end function ap_potential

   subroutine maclaurin2bnd_potential

      use cg_leaves,  only: leaves
      use cg_list,    only: cg_list_element
      use constants,  only: xdim, ydim, zdim, LO, HI, GEO_XYZ
      use dataio_pub, only: die
      use domain,     only: dom
      use grid_cont,  only: grid_container
      use units,      only: fpiG

      implicit none

      integer :: i, j, k
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      if (dom%geometry_type /= GEO_XYZ) call die("[problem_pub:maclaurin2bnd_potential] only cartesian geometry implemented")

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         if (any(cg%ext_bnd(xdim, :))) then
            do j = cg%js, cg%je
               do k = cg%ks, cg%ke
                  if (cg%ext_bnd(xdim, LO)) cg%mg%bnd_x(j, k, LO) = ap_potential(cg%fbnd(xdim, LO), cg%y(j), cg%z(k)) * fpiG
                  if (cg%ext_bnd(xdim, HI)) cg%mg%bnd_x(j, k, HI) = ap_potential(cg%fbnd(xdim, HI), cg%y(j), cg%z(k)) * fpiG
               enddo
            enddo
         endif

         if (any(cg%ext_bnd(ydim, :))) then
            do i = cg%is, cg%ie
               do k = cg%ks, cg%ke
                  if (cg%ext_bnd(ydim, LO)) cg%mg%bnd_y(i, k, LO) = ap_potential(cg%x(i), cg%fbnd(ydim, LO), cg%z(k)) * fpiG
                  if (cg%ext_bnd(ydim, HI)) cg%mg%bnd_y(i, k, HI) = ap_potential(cg%x(i), cg%fbnd(ydim, HI), cg%z(k)) * fpiG
               enddo
            enddo
         endif

         if (any(cg%ext_bnd(zdim, :))) then
            do i = cg%is, cg%ie
               do j = cg%js, cg%je
                  if (cg%ext_bnd(zdim, LO)) cg%mg%bnd_z(i, j, LO) = ap_potential(cg%x(i), cg%y(j), cg%fbnd(zdim, LO)) * fpiG
                  if (cg%ext_bnd(zdim, HI)) cg%mg%bnd_z(i, j, HI) = ap_potential(cg%x(i), cg%y(j), cg%fbnd(zdim, HI)) * fpiG
               enddo
            enddo
         endif
         cgl => cgl%nxt
      enddo


   end subroutine maclaurin2bnd_potential

#endif /* MACLAURIN_PROBLEM */

end module problem_pub
