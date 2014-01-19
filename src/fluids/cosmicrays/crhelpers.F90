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
!! \brief Module containing routines for divergence computing which are used (so far only) in cosmic rays integration schemes
!<
module crhelpers
! pulled by COSM_RAYS
   use constants, only: dsetnamelen
   implicit none

   private
   public :: div_v, set_div_v1d, divv_n
#if defined(__INTEL_COMPILER) || defined(_CRAYFTN)
      !! \deprecated remove this clause as soon as Intel Compiler gets required
      !! features and/or bug fixes
   public :: init_div_v
#endif /* __INTEL_COMPILER || _CRAYFTN */

   interface
      subroutine div_v_func(ifluid, cg)
         use grid_cont, only: grid_container
         implicit none
         integer(kind=4),               intent(in)    :: ifluid
         type(grid_container), pointer, intent(inout) :: cg
      end subroutine div_v_func
   end interface

   character(len=dsetnamelen), parameter :: divv_n = "divvel" !< divergence of velocity

#if defined(__INTEL_COMPILER) || defined(_CRAYFTN)
      !! \deprecated remove this clause as soon as Intel Compiler gets required
      !! features and/or bug fixes
   procedure(div_v_func), pointer :: div_v
#else /* !__INTEL_COMPILER && !_CRAYFTN */
   procedure(div_v_func), pointer :: div_v => init_divv
#endif /* !__INTEL_COMPILER && !_CRAYFTN */

contains

#if defined(__INTEL_COMPILER) || defined(_CRAYFTN)
   !! \deprecated remove this clause as soon as Intel Compiler gets required
   !! features and/or bug fixes
   subroutine init_div_v

      use initcosmicrays, only: divv_scheme

      implicit none

      select case (trim(divv_scheme))
         case ("6lp", "6th_order_legandre")
            div_v => div_v_6th_lp
         case default
            div_v => div_v_1st
      end select

      return
   end subroutine init_div_v
#endif /* __INTEL_COMPILER || _CRAYFTN */

   subroutine init_divv(ifluid, cg)

      use grid_cont,      only: grid_container
      use initcosmicrays, only: divv_scheme

      implicit none

      integer(kind=4),               intent(in)    :: ifluid
      type(grid_container), pointer, intent(inout) :: cg

      select case (trim(divv_scheme))
         case ("6lp", "6th_order_legandre")
            div_v => div_v_6th_lp
         case default
            div_v => div_v_1st
      end select

      call div_v(ifluid, cg)
      return
   end subroutine init_divv

   subroutine set_div_v1d(p, dir, i1, i2, cg)

      use dataio_pub,       only: die
      use grid_cont,        only: grid_container
      use named_array_list, only: qna

      implicit none

      integer(kind=4),               intent(in)    :: dir
      integer,                       intent(in)    :: i1, i2
      real, dimension(:),   pointer, intent(inout) :: p
      type(grid_container), pointer, intent(inout) :: cg

      if (.not. qna%exists(divv_n)) call die("[crhelpers:set_div_v1d] cannot get divvel")
      p => cg%q(qna%ind(divv_n))%get_sweep(dir, i1, i2)

   end subroutine set_div_v1d

!>
!! \brief Compute divergence of velocity
!!
!! \details This routine requires a single layer of valid guardcells in cg%u arrays
!!
!! The divergence of velocity computed with the aid of 6-th order finite
!! differencing based on the Lagendre Polynomial interpolation.
!!
!! \todo Should be moved to a dedicated module containing general purpose interpolation
!! and derivation routines, and placed together with the other useful scheems described
!! in particular in http://turbulence.pha.jhu.edu/Database-functions.pdf
!<
   subroutine div_v_6th_lp(ifluid, cg)

      use constants,        only: xdim, zdim, pdims, big, LO, HI, ORTHO1, ORTHO2
      use domain,           only: dom
      use fluidindex,       only: iarr_all_dn
      use grid_cont,        only: grid_container
      use named_array_list, only: qna, wna

      implicit none

      integer(kind=4),               intent(in)    :: ifluid
      type(grid_container), pointer, intent(inout) :: cg
      real, dimension(:),   pointer                :: divvel, mom, dens
      integer(kind=4)                              :: dir
      integer                                      :: i2, i3
      real, parameter                              :: p3_4 = 3./4., m3_20 = -3./20., p1_60 = 1./60.

      cg%q(qna%ind(divv_n))%arr(:,:,:) = 0.0

      do dir = xdim, zdim
         if (.not. dom%has_dir(dir)) cycle
         do i2 = cg%lhn(pdims(dir, ORTHO1), LO), cg%lhn(pdims(dir, ORTHO1), HI)
            do i3 = cg%lhn(pdims(dir, ORTHO2), LO), cg%lhn(pdims(dir, ORTHO2), HI)
               divvel => cg%q(qna%ind(divv_n))%get_sweep(dir, i2, i3)
               mom  => cg%w(wna%fi)%get_sweep(dir, iarr_all_dn(ifluid) + dir, i2, i3)
               dens => cg%w(wna%fi)%get_sweep(dir, iarr_all_dn(ifluid)      , i2, i3)
               associate( &
#ifndef _CRAYFTN
                  vv => mom(:) / dens(:), &
#endif /* _CRAYFTN */
                  nn => cg%n_(dir), &
                  idl => cg%idl(dir) &
               )
#ifndef _CRAYFTN
                  divvel(4:nn-3)  = divvel(4:nn-3) + (vv(5:nn-2) - vv(3:nn-4)) * (p3_4  * idl)
                  divvel(4:nn-3)  = divvel(4:nn-3) + (vv(6:nn-1) - vv(2:nn-5)) * (m3_20 * idl)
                  divvel(4:nn-3)  = divvel(4:nn-3) + (vv(7:nn  ) - vv(1:nn-6)) * (p1_60 * idl)
#else /* _CRAFTN */
                  divvel(4:nn-3)  = divvel(4:nn-3) + (mom(5:nn-2)/dens(5:nn-2) - mom(3:nn-4)/dens(3:nn-4)) * (p3_4  * idl)
                  divvel(4:nn-3)  = divvel(4:nn-3) + (mom(6:nn-1)/dens(6:nn-1) - mom(2:nn-5)/dens(2:nn-5)) * (m3_20 * idl)
                  divvel(4:nn-3)  = divvel(4:nn-3) + (mom(7:nn  )/dens(7:nn  ) - mom(1:nn-6)/dens(1:nn-6)) * (p1_60 * idl)
#endif /* _CRAYFTN */
                  divvel(1:3)     = big
                  divvel(nn-2:nn) = big
               end associate
            enddo
         enddo
      enddo

   end subroutine div_v_6th_lp
!>
!! \brief Compute divergence of velocity
!!
!! \details This routine requires a single layer of valid guardcells uin cg%u arrays
!<

   subroutine div_v_1st(ifluid, cg)

      use constants,        only: xdim, zdim, pdims, half, LO, HI, ORTHO1, ORTHO2
      use domain,           only: dom
      use fluidindex,       only: iarr_all_dn
      use grid_cont,        only: grid_container
      use named_array_list, only: qna, wna

      implicit none

      integer(kind=4),               intent(in)    :: ifluid
      type(grid_container), pointer, intent(inout) :: cg
      real, dimension(:),   pointer                :: divvel, mom, dn
      integer(kind=4)                              :: dir
      integer                                      :: i2, i3

      cg%q(qna%ind(divv_n))%arr(:,:,:) = 0.0

      do dir = xdim, zdim
         if (.not.dom%has_dir(dir)) cycle
         do i2 = cg%lhn(pdims(dir, ORTHO1), LO), cg%lhn(pdims(dir, ORTHO1), HI)
            do i3 = cg%lhn(pdims(dir, ORTHO2), LO), cg%lhn(pdims(dir, ORTHO2), HI)
               divvel => cg%q(qna%ind(divv_n))%get_sweep(dir, i2, i3)
               mom    => cg%w(wna%fi)%get_sweep(dir, iarr_all_dn(ifluid)+dir, i2, i3)
               dn     => cg%w(wna%fi)%get_sweep(dir, iarr_all_dn(ifluid)    , i2, i3)
               associate( &
#ifndef _CRAYFTN
                  vv => mom(:) / dn(:), &
#endif /* _CRAYFTN */
                  nn => cg%n_(dir), &
                  idl => cg%idl(dir) &
               )
#ifndef _CRAYFTN
                  divvel(2:nn-1) = divvel(2:nn-1) + (vv(3:nn) - vv(1:nn-2)) * (half * idl)
#else /* _CRAYFTN */
                  divvel(2:nn-1) = divvel(2:nn-1) + (mom(3:nn) / dn(3:nn) - mom(1:nn-2) / dn(1:nn-2)) * (half * idl)
#endif /* _CRAYFTN */
                  divvel(1) = divvel(2)
                  divvel(nn) = divvel(nn-1) ! for sanity
               end associate
            enddo
         enddo
      enddo

   end subroutine div_v_1st

end module crhelpers
