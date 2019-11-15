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

!>  \brief Auxiliary functions used by particle integrators

module particle_func
! pulled by GRAV

   implicit none

   private
   public :: particle_in_area, check_ord, df_d_p, d2f_d2_p, d2f_dd_p, df_d_o2, d2f_d2_o2, d2f_dd_o2

   interface
      function dxi(cell, cg, ig, dir)
         use constants, only: ndims
         use grid_cont, only: grid_container
         implicit none
         integer, dimension(ndims),     intent(in) :: cell
         type(grid_container), pointer, intent(in) :: cg
         integer(kind=4),               intent(in) :: ig, dir
         real                                      :: dxi
      end function dxi

      function d2dxi2(cell, cg, ig, dir)
         use constants, only: ndims
         use grid_cont, only: grid_container
         implicit none
         integer, dimension(ndims),     intent(in) :: cell
         type(grid_container), pointer, intent(in) :: cg
         integer(kind=4),               intent(in) :: ig, dir
         real                                      :: d2dxi2
      end function d2dxi2

      function d2dxixj(cell, cg, ig, dir1, dir2)
         use constants, only: ndims
         use grid_cont, only: grid_container
         implicit none
         integer, dimension(ndims),     intent(in) :: cell
         type(grid_container), pointer, intent(in) :: cg
         integer(kind=4),               intent(in) :: ig, dir1, dir2
         real                                      :: d2dxixj
      end function d2dxixj

   end interface

   procedure(dxi),     pointer :: df_d_p   => NULL()
   procedure(d2dxi2),  pointer :: d2f_d2_p => NULL()
   procedure(d2dxixj), pointer :: d2f_dd_p => NULL()

contains

!> \brief check if the particle locates inside given area

   function particle_in_area(pos, area) result(itis)

      use constants, only: ndims, LO, HI

      implicit none

      real, dimension(ndims),       intent(in) :: pos
      real, dimension(ndims,LO:HI), intent(in) :: area
      logical                                  :: itis

      itis = (all(pos >= area(:,LO)) .and. all(pos <= area(:,HI)))

   end function particle_in_area

   subroutine check_ord(order)

      implicit none

      integer, intent(in)    :: order

      if (order == 2) then
         df_d_p   => df_d_o2
         d2f_d2_p => d2f_d2_o2
         d2f_dd_p => d2f_dd_o2
      else
         df_d_p   => df_d_o4
         d2f_d2_p => d2f_d2_o4
         d2f_dd_p => d2f_dd_o4
      endif

   end subroutine check_ord

   function df_d_o2(cell, cg, ig, dir)

      use constants, only: idm, ndims, half
      use grid_cont, only: grid_container

      implicit none

      integer, dimension(ndims),     intent(in) :: cell
      type(grid_container), pointer, intent(in) :: cg
      integer(kind=4),               intent(in) :: ig, dir
      real, target                              :: df_d_o2

      !o(R^2)
      df_d_o2 = (cg%q(ig)%point(cell+idm(dir,:)) - cg%q(ig)%point(cell-idm(dir,:)) ) * half *cg%idl(dir)

   end function df_d_o2

   function d2f_d2_o2(cell, cg, ig, dir)

      use constants, only: idm, ndims, two
      use grid_cont, only: grid_container

      implicit none

      integer, dimension(ndims),     intent(in) :: cell
      type(grid_container), pointer, intent(in) :: cg
      integer(kind=4),               intent(in) :: ig, dir
      real, target                              :: d2f_d2_o2

      !o(R^2)
      d2f_d2_o2 = (cg%q(ig)%point(cell+idm(dir,:)) - two*cg%q(ig)%point(cell) + cg%q(ig)%point(cell-idm(dir,:)) ) * cg%idl(dir)**2

   end function d2f_d2_o2

   function d2f_dd_o2(cell, cg, ig, dir1, dir2)

      use constants, only: idm, ndims, oneq
      use grid_cont, only: grid_container

      implicit none

      integer, dimension(ndims),     intent(in) :: cell
      type(grid_container), pointer, intent(in) :: cg
      integer(kind=4),               intent(in) :: ig, dir1, dir2
      real, target                              :: d2f_dd_o2

      !o(R^2)
      d2f_dd_o2 = (cg%q(ig)%point(cell+idm(dir1,:)+idm(dir2,:)) - cg%q(ig)%point(cell+idm(dir1,:)-idm(dir2,:)) + &
                   cg%q(ig)%point(cell-idm(dir1,:)-idm(dir2,:)) - cg%q(ig)%point(cell-idm(dir1,:)+idm(dir2,:)) ) * oneq*cg%idl(dir1)*cg%idl(dir2)

   end function d2f_dd_o2

   function df_d_o4(cell, cg, ig, dir)

      use constants, only: idm, ndims, onet
      use grid_cont, only: grid_container

      implicit none

      integer, dimension(ndims),     intent(in) :: cell
      type(grid_container), pointer, intent(in) :: cg
      integer(kind=4),               intent(in) :: ig, dir
      real, target                              :: df_d_o4

      !o(R^4)
      df_d_o4 = 2.0 * (cg%q(ig)%point(cell +   idm(dir,:)) - cg%q(ig)%point(cell -   idm(dir,:)) ) * cg%idl(dir) * onet - &
                      (cg%q(ig)%point(cell + 2*idm(dir,:)) - cg%q(ig)%point(cell - 2*idm(dir,:)) ) * cg%idl(dir) / 12.0

   end function df_d_o4

   function d2f_d2_o4(cell, cg, ig, dir)

      use constants, only: idm, ndims, two, onet
      use grid_cont, only: grid_container

      implicit none

      integer, dimension(ndims),     intent(in) :: cell
      type(grid_container), pointer, intent(in) :: cg
      integer(kind=4),               intent(in) :: ig, dir
      real, target                              :: d2f_d2_o4

      !o(R^4)
      d2f_d2_o4 = 4.0 * (cg%q(ig)%point(cell +   idm(dir,:)) + cg%q(ig)%point(cell -   idm(dir,:)) - two *  cg%q(ig)%point(cell)) * cg%idl(dir)**2 * onet - &
                        (cg%q(ig)%point(cell + 2*idm(dir,:)) + cg%q(ig)%point(cell - 2*idm(dir,:)) - two *  cg%q(ig)%point(cell)) * cg%idl(dir)**2 / 12.0

   end function d2f_d2_o4

   function d2f_dd_o4(cell, cg, ig, dir1, dir2)

      use constants, only: idm, ndims, onet
      use grid_cont, only: grid_container

      implicit none

      integer, dimension(ndims),     intent(in) :: cell
      type(grid_container), pointer, intent(in) :: cg
      integer(kind=4),               intent(in) :: ig, dir1, dir2
      real, target                              :: d2f_dd_o4

      !o(R^4)
      d2f_dd_o4 = (cg%q(ig)%point(cell +   idm(dir1,:) +   idm(dir2,:)) + cg%q(ig)%point(cell -   idm(dir1,:) -   idm(dir2,:)) - &
                   cg%q(ig)%point(cell +   idm(dir1,:) -   idm(dir2,:)) - cg%q(ig)%point(cell -   idm(dir1,:) +   idm(dir2,:)) ) * cg%idl(dir1)*cg%idl(dir2) * onet - &
                  (cg%q(ig)%point(cell + 2*idm(dir1,:) + 2*idm(dir2,:)) + cg%q(ig)%point(cell - 2*idm(dir1,:) - 2*idm(dir2,:)) - &
                   cg%q(ig)%point(cell + 2*idm(dir1,:) - 2*idm(dir2,:)) - cg%q(ig)%point(cell - 2*idm(dir1,:) + 2*idm(dir2,:)) ) * cg%idl(dir1)*cg%idl(dir2) / 48.0

   end function d2f_dd_o4


end module particle_func
