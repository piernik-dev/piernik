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
!    HLLD Riemann solver for ideal magnetohydrodynamics
!    Varadarajan Parthasarathy, CAMK, Warszawa. 2015.
!    Dr. Artur Gawryszczak, CAMK, Warszawa.
!
!    Energy fix up routines for CT and its related comments are not used in the current version.
!    The algorithm is simply present for experimental purposes.
!--------------------------------------------------------------------------------------------------------------

#include "piernik.def"

module fluxlimiters
! pulled by RIEMANN

   implicit none
   private
   public :: set_limiters, limiter, flimiter, blimiter

   interface
      function limiter(q) result(dq)

         implicit none

         real, dimension(:,:), intent(in) :: q

         real, dimension(size(q,1), size(q,2)) :: dq

      end function limiter
   end interface

   procedure(limiter), pointer :: flimiter => null(), blimiter => null()

contains

!!!-------------------------------------------------------------------------------

   subroutine set_limiters(flim_str, blim_str)

      use dataio_pub, only: die

      implicit none

      character(len=*), intent(in) :: flim_str, blim_str

      if (associated(flimiter)) call die("[fluxlimiters:set_limiters] flimiter already associated")
      call set_limiter(flim_str, flimiter)

      if (associated(blimiter)) call die("[fluxlimiters:set_limiters] blimiter already associated")
      call set_limiter(blim_str, blimiter)

   end subroutine set_limiters

!!!-------------------------------------------------------------------------------

   subroutine set_limiter(lim_str, lim)

      use dataio_pub, only: msg, die

      implicit none

      character(len=*),            intent(in)  :: lim_str
      procedure(limiter), pointer, intent(out) :: lim

      select case (lim_str)
         case ('vanleer', 'VANLEER')
            lim => calculate_slope_vanleer
         case ('minmod', 'MINMOD')
            lim => slope_limiter_minmod
         case ('moncen', 'MONCEN')
            lim => slope_limiter_moncen
         case ('superbee', 'SUPERBEE')
            lim => slope_limiter_superbee
         case default
            write(msg,'(2a)') "[fluxlimiters:set_limiter] unknown limiter ", lim_str
            call die(msg)
            lim => null()
      end select

   end subroutine set_limiter

!!!-------------------------------------------------------------------------------

   function calculate_slope_vanleer(u) result(dq)

      implicit none

      real, dimension(:,:), intent(in)     :: u

      real, dimension(size(u,1),size(u,2)) :: dlft, drgt, dcen, dq
      integer :: n

      n = size(u,2)

      dlft(:,2:n)   = (u(:,2:n) - u(:,1:n-1)) ; dlft(:,1) = dlft(:,2)    ! (14.38)
      drgt(:,1:n-1) = dlft(:,2:n) ;             drgt(:,n) = drgt(:,n-1)

      dcen = dlft*drgt

      where (dcen>0.0)
         dq = 2.0*dcen / (dlft+drgt)       ! (14.54)
      elsewhere
         dq = 0.0
      endwhere

   end function calculate_slope_vanleer

!-----------------------------------------------------------------------------------------------------------------------

! Minmod slope limiter method. Rezolla book, Pg 428 Eq. 9.44, 9.45

   function slope_limiter_minmod(u) result(dq)

      use constants, only: half, one

      implicit none

      real, dimension(:,:), intent(in)     :: u

      real, dimension(size(u,1),size(u,2)) :: dlft, drgt, dq
      integer :: n

      n = size(u,2)

      dlft(:,2:n)   = (u(:,2:n) - u(:,1:n-1)) ; dlft(:,1) = dlft(:,2)
      drgt(:,1:n-1) = dlft(:,2:n) ;             drgt(:,n) = drgt(:,n-1)

      dq = (sign(one, dlft) + sign(one, drgt))*min(abs(dlft),abs(drgt))*half

   end function slope_limiter_minmod

!-----------------------------------------------------------------------------------------------------------------------

! Monotonized central limiter

   function slope_limiter_moncen(u) result(dq)

      use constants, only: half, one, two

      implicit none

      real, dimension(:,:), intent(in)     :: u

      real, dimension(size(u,1),size(u,2)) :: dlft, drgt, dq
      integer :: n

      n = size(u,2)

      dlft(:,2:n)   = (u(:,2:n) - u(:,1:n-1)) ; dlft(:,1) = dlft(:,2)
      drgt(:,1:n-1) = dlft(:,2:n) ;             drgt(:,n) = drgt(:,n-1)

      dq = (sign(one,dlft)+sign(one,drgt))*min(two*abs(dlft),two*abs(drgt),half*abs(dlft+drgt))*half

   end function slope_limiter_moncen
!-----------------------------------------------------------------------------------------------------------------------

! Super-bee limiter

   function slope_limiter_superbee(u) result(dq)

      use constants, only: half, one, two

      implicit none

      real, dimension(:,:), intent(in)     :: u

      real, dimension(size(u,1),size(u,2)) :: dlft, drgt, dq
      integer :: n

      n = size(u,2)

      dlft(:,2:n)   = (u(:,2:n) - u(:,1:n-1)) ; dlft(:,1) = dlft(:,2)
      drgt(:,1:n-1) = dlft(:,2:n) ;             drgt(:,n) = drgt(:,n-1)

      where (abs(dlft) > abs(drgt))
         dq = (sign(one,dlft)+sign(one,drgt))*min(abs(dlft), abs(two*drgt))*half
      elsewhere
         dq = (sign(one,dlft)+sign(one,drgt))*min(abs(two*dlft), abs(drgt))*half
      endwhere

    end function slope_limiter_superbee

end module fluxlimiters
