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
!! \brief Routines used to interpolate grid fields onto particle positions
!<

module particle_interpolation
! pulled by ANY
   implicit none

   private
   public :: quadratic_spline

contains

!>
!! \brief Quadratic spline interpolation
!<
    function quadratic_spline(f, dp, ijkp) result(gp)

      use constants,    only: xdim, ydim, zdim, IM, I0, IP, ndims
      use domain,       only: dom

      implicit none

      real, dimension (:,:,:,:), pointer, intent(in) :: f   !< field
      real, dimension(ndims), intent(in)    :: dp           !< point position in [0,1] range
      integer, dimension(ndims), intent(in) :: ijkp         !< nearest grid point where dp resides
      real, dimension(size(f,1)) :: gp                      !< interpolated value

      ! locals
      real, dimension(ndims, IM:IP) :: fac

!     This should be moved to particle specific modules, to avoid depending on  cg
!     do cdim = xdim, zdim
!        dp(cdim) = (xp(cdim) - cg%coord(cdim)%r(ijkp(cdim))) * cg%idl(ijkp(cdim))
!     enddo

      !  Interpolation coefficients
      fac(:, :) = 0.0
      where (dom%has_dir)
         fac(:, IM) = 0.5 * (0.5 - dp(:))**2
         fac(:, I0) = 0.75 - dp(:)**2
         fac(:, IP) = 0.5 * (0.5 + dp(:))**2
      elsewhere
         fac(:, IM) = 0.0
         fac(:, I0) = 1.0
         fac(:, IP) = 0.0
      endwhere

      ! Some computations can be avoided by spliting following statement into pieces
      ! and using ifs
      gp(:) = fac(xdim, I0)*fac(ydim, I0)*fac(zdim, I0)*f(ijkp(xdim),ijkp(ydim),ijkp(zdim),:) + &
              fac(xdim, I0)*fac(ydim, I0)*( f(ijkp(xdim)             ,ijkp(ydim)             ,ijkp(zdim)+dom%D_(zdim),:)*fac(zdim, IP) + &
                                            f(ijkp(xdim)             ,ijkp(ydim)             ,ijkp(zdim)-dom%D_(zdim),:)*fac(zdim, IM) ) + &
              fac(xdim, I0)*fac(zdim, I0)*( f(ijkp(xdim)             ,ijkp(ydim)+dom%D_(ydim),ijkp(zdim)             ,:)*fac(ydim, IP) + &
                                            f(ijkp(xdim)             ,ijkp(ydim)-dom%D_(ydim),ijkp(zdim)             ,:)*fac(ydim, IM) ) + &
              fac(ydim, I0)*fac(zdim, I0)*( f(ijkp(xdim)+dom%D_(xdim),ijkp(ydim)             ,ijkp(zdim)             ,:)*fac(xdim, IP) + &
                                            f(ijkp(xdim)-dom%D_(xdim),ijkp(ydim)             ,ijkp(zdim)             ,:)*fac(xdim, IM) ) + &
              fac(xdim, IP)*fac(ydim, IP)*( f(ijkp(xdim)+dom%D_(xdim),ijkp(ydim)+dom%D_(ydim),ijkp(zdim)+dom%D_(zdim),:)*fac(zdim, IP) + &
                                            f(ijkp(xdim)+dom%D_(xdim),ijkp(ydim)+dom%D_(ydim),ijkp(zdim)-dom%D_(zdim),:)*fac(zdim, IM) ) + &
              fac(xdim, IP)*fac(ydim, IM)*( f(ijkp(xdim)+dom%D_(xdim),ijkp(ydim)-dom%D_(ydim),ijkp(zdim)+dom%D_(zdim),:)*fac(zdim, IP) + &
                                            f(ijkp(xdim)+dom%D_(xdim),ijkp(ydim)-dom%D_(ydim),ijkp(zdim)-dom%D_(zdim),:)*fac(zdim, IM) ) + &
              fac(xdim, IM)*fac(ydim, IP)*( f(ijkp(xdim)-dom%D_(xdim),ijkp(ydim)+dom%D_(ydim),ijkp(zdim)+dom%D_(zdim),:)*fac(zdim, IP) + &
                                            f(ijkp(xdim)-dom%D_(xdim),ijkp(ydim)+dom%D_(ydim),ijkp(zdim)-dom%D_(zdim),:)*fac(zdim, IM) ) + &
              fac(xdim, IM)*fac(ydim, IM)*( f(ijkp(xdim)-dom%D_(xdim),ijkp(ydim)-dom%D_(ydim),ijkp(zdim)+dom%D_(zdim),:)*fac(zdim, IP) + &
                                            f(ijkp(xdim)-dom%D_(xdim),ijkp(ydim)-dom%D_(ydim),ijkp(zdim)-dom%D_(zdim),:)*fac(zdim, IM) ) + &
              fac(xdim, I0)*fac(ydim, IP)*( f(ijkp(xdim)             ,ijkp(ydim)+dom%D_(ydim),ijkp(zdim)+dom%D_(zdim),:)*fac(zdim, IP) + &
                                            f(ijkp(xdim)             ,ijkp(ydim)+dom%D_(ydim),ijkp(zdim)-dom%D_(zdim),:)*fac(zdim, IM) ) + &
              fac(xdim, I0)*fac(ydim, IM)*( f(ijkp(xdim)             ,ijkp(ydim)-dom%D_(ydim),ijkp(zdim)+dom%D_(zdim),:)*fac(zdim, IP) + &
                                            f(ijkp(xdim)             ,ijkp(ydim)-dom%D_(ydim),ijkp(zdim)-dom%D_(zdim),:)*fac(zdim, IM) ) + &
              fac(ydim, I0)*fac(zdim, IP)*( f(ijkp(xdim)+dom%D_(xdim),ijkp(ydim)             ,ijkp(zdim)+dom%D_(zdim),:)*fac(xdim, IP) + &
                                            f(ijkp(xdim)-dom%D_(xdim),ijkp(ydim)             ,ijkp(zdim)+dom%D_(zdim),:)*fac(xdim, IM) ) + &
              fac(ydim, I0)*fac(zdim, IM)*( f(ijkp(xdim)+dom%D_(xdim),ijkp(ydim)             ,ijkp(zdim)-dom%D_(zdim),:)*fac(xdim, IP) + &
                                            f(ijkp(xdim)-dom%D_(xdim),ijkp(ydim)             ,ijkp(zdim)-dom%D_(zdim),:)*fac(xdim, IM) ) + &
              fac(zdim, I0)*fac(xdim, IP)*( f(ijkp(xdim)+dom%D_(xdim),ijkp(ydim)+dom%D_(ydim),ijkp(zdim)             ,:)*fac(ydim, IP) + &
                                            f(ijkp(xdim)+dom%D_(xdim),ijkp(ydim)-dom%D_(ydim),ijkp(zdim)             ,:)*fac(ydim, IM) ) + &
              fac(zdim, I0)*fac(xdim, IM)*( f(ijkp(xdim)-dom%D_(xdim),ijkp(ydim)+dom%D_(ydim),ijkp(zdim)             ,:)*fac(ydim, IP) + &
                                            f(ijkp(xdim)-dom%D_(xdim),ijkp(ydim)-dom%D_(ydim),ijkp(zdim)             ,:)*fac(ydim, IM) )
    end function quadratic_spline
end module
