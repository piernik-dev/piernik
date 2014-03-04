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

!> \brief This module contains standard routines for finding where to refine and where to derefine

module refinement_filters

   implicit none

   private
   public :: refine_on_gradient, refine_on_relative_gradient, ref_crit

   type :: ref_crit
      integer              :: iv                      !< field index in cg%q or cg%w array
      integer              :: ic                      !< component index (cg%w(iv)%arr(ic,:,:,:)) or INVALID for 3D arrays
      real                 :: ref_thr                 !< refinement threshold
      real                 :: deref_thr               !< derefinement threshold
      real                 :: aux                     !< auxiliary parameter (can be smoother or filter strength)
      procedure(refine_crit), pass, pointer :: refine !< refinement routine
   end type ref_crit

   ! all routines that are public in this module should confotm to this interface
   interface
      subroutine refine_crit(this, cg, p3d)

         use grid_cont, only: grid_container
         import ref_crit

         implicit none

         class(ref_crit),                 intent(in)    :: this !< this contains refinement parameters
         type(grid_container), pointer,   intent(inout) :: cg   !< current grid piece
         real, dimension(:,:,:), pointer, intent(in)    :: p3d  !< pointer to array to be examined for (de)refinement
      end subroutine refine_crit
   end interface

contains

!>
!! \brief Refine/derefine based on ||grad u||
!! This is sensitive to gradients, but the thresholds must be rescaled, when you change units of the problem.
!<

   subroutine refine_on_gradient(this, cg, p3d)

      use domain,    only: dom
      use grid_cont, only: grid_container

      implicit none

      class(ref_crit),                 intent(in)    :: this !< this contains refinement parameters
      type(grid_container), pointer,   intent(inout) :: cg   !< current grid piece
      real, dimension(:,:,:), pointer, intent(in)    :: p3d  !< pointer to array to be examined for (de)refinement needs (should contain at least one layer of updated guardcells)

      integer :: i, j, k
      real :: r, max_r

      max_r = -huge(1.)

      !> \todo implement how far we should look for (de)refinements

      do k = cg%ks, cg%ke
         do j = cg%js, cg%je
            do i = cg%is, cg%ie
               r = grad2(i, j, k)
               max_r = max(max_r, r)
               cg%refinemap(i, j, k) = cg%refinemap(i, j, k) .or. (r >= this%ref_thr**2)
               ! we can avoid calculating square root here
            enddo
         enddo
      enddo

      ! check additional 1 perimeter of cells for derefinement
      max_r = max(max_r, &
           maxgradoverarea(cg%is-dom%D_x, cg%is-dom%D_x, cg%js-dom%D_y, cg%je+dom%D_y, cg%ks-dom%D_z, cg%ke+dom%D_z), &
           maxgradoverarea(cg%ie+dom%D_x, cg%ie+dom%D_x, cg%js-dom%D_y, cg%je+dom%D_y, cg%ks-dom%D_z, cg%ke+dom%D_z), &
           maxgradoverarea(cg%is, cg%ie, cg%js-dom%D_y, cg%js-dom%D_y, cg%ks-dom%D_z, cg%ke+dom%D_z), &
           maxgradoverarea(cg%is, cg%ie, cg%je+dom%D_y, cg%je+dom%D_y, cg%ks-dom%D_z, cg%ke+dom%D_z), &
           maxgradoverarea(cg%is, cg%ie, cg%js, cg%je, cg%ks-dom%D_z, cg%ks-dom%D_z), &
           maxgradoverarea(cg%is, cg%ie, cg%js, cg%je, cg%ke+dom%D_z, cg%ke+dom%D_z) )

      cg%refine_flags%derefine = cg%refine_flags%derefine .or. (max_r < this%deref_thr**2)

   contains

      !> \brief scan given area of indices for maximum value of specified function

      elemental real function maxgradoverarea(i1, i2, j1, j2, k1, k2)

         implicit none

         integer, intent(in) :: i1, i2, j1, j2, k1, k2

         integer :: i, j, k

         maxgradoverarea = -huge(1.)

         do k = k1, k2
            do j = j1, j2
               do i = i1, i2
                  maxgradoverarea = max(maxgradoverarea, grad2(i, j, k))
               enddo
            enddo
         enddo

      end function maxgradoverarea

      !>
      !! \brief Square of a gradient without taking into account cell sizes
      !!
      !! \details Perhaps it is not gradient but rather a norm of 3D difference
      !<

      elemental real function grad2(i, j, k)

         use domain, only: dom

         implicit none

         integer,                         intent(in) :: i, j, k !< indices

         grad2 = (p3d(i+dom%D_x, j, k) - p3d(i-dom%D_x, j, k))**2 + &
              &  (p3d(i, j+dom%D_y, k) - p3d(i, j-dom%D_y, k))**2 + &
              &  (p3d(i, j, k+dom%D_z) - p3d(i, j, k-dom%D_z))**2

      end function grad2

   end subroutine refine_on_gradient

!>
!! \brief Refine/derefine based on ||grad u||/||u||
!! This is sensitive to changes of sign or strong gradients (such as fast approaching 0.)
!<

   subroutine refine_on_relative_gradient(this, cg, p3d)

      use domain,    only: dom
      use grid_cont, only: grid_container

      implicit none

      class(ref_crit),                 intent(in)    :: this !< this contains refinement parameters
      type(grid_container), pointer,   intent(inout) :: cg   !< current grid piece
      real, dimension(:,:,:), pointer, intent(in)    :: p3d  !< pointer to array to be examined for (de)refinement needs (should contain at least one layer of updated guardcells)

      integer :: i, j, k
      real :: r, max_r

      max_r = -huge(1.)

      !> \todo implement how far we should look for (de)refinements

      do k = cg%ks, cg%ke
         do j = cg%js, cg%je
            do i = cg%is, cg%ie
               r = rel_grad2(i, j, k)
               max_r = max(max_r, r)
               cg%refinemap(i, j, k) = cg%refinemap(i, j, k) .or. (r >= this%ref_thr**2)
               ! we can avoid calculating square root here
            enddo
         enddo
      enddo

      ! check additional 1 perimeter of cells for derefinement
      max_r = max(max_r, &
           maxrelgradoverarea(cg%is-dom%D_x, cg%is-dom%D_x, cg%js-dom%D_y, cg%je+dom%D_y, cg%ks-dom%D_z, cg%ke+dom%D_z), &
           maxrelgradoverarea(cg%ie+dom%D_x, cg%ie+dom%D_x, cg%js-dom%D_y, cg%je+dom%D_y, cg%ks-dom%D_z, cg%ke+dom%D_z), &
           maxrelgradoverarea(cg%is, cg%ie, cg%js-dom%D_y, cg%js-dom%D_y, cg%ks-dom%D_z, cg%ke+dom%D_z), &
           maxrelgradoverarea(cg%is, cg%ie, cg%je+dom%D_y, cg%je+dom%D_y, cg%ks-dom%D_z, cg%ke+dom%D_z), &
           maxrelgradoverarea(cg%is, cg%ie, cg%js, cg%je, cg%ks-dom%D_z, cg%ks-dom%D_z), &
           maxrelgradoverarea(cg%is, cg%ie, cg%js, cg%je, cg%ke+dom%D_z, cg%ke+dom%D_z) )

      cg%refine_flags%derefine = cg%refine_flags%derefine .or. (max_r < this%deref_thr**2)

   contains

      !> \brief scan given area of indices for maximum value of specified function

      elemental real function maxrelgradoverarea(i1, i2, j1, j2, k1, k2)

         implicit none

         integer, intent(in) :: i1, i2, j1, j2, k1, k2

         integer :: i, j, k

         maxrelgradoverarea = -huge(1.)

         do k = k1, k2
            do j = j1, j2
               do i = i1, i2
                  maxrelgradoverarea = max(maxrelgradoverarea, rel_grad2(i, j, k))
               enddo
            enddo
         enddo

      end function maxrelgradoverarea

      !>
      !! \brief Square of a 'relative gradient' without taking into account cell sizes
      !!
      !! \details This routine should return value between 0 and 1.
      !<

      elemental real function rel_grad2(i, j, k)

         use constants, only: xdim, ydim, zdim
         use domain,    only: dom

         implicit none

         integer,                         intent(in) :: i, j, k !< indices

         rel_grad2 = 0.

         if (dom%has_dir(xdim)) rel_grad2 = sumsq01( &
              rel_grad_1pair(p3d(i, j, k), p3d(i+dom%D_x, j, k)), &
              rel_grad_1pair(p3d(i, j, k), p3d(i-dom%D_x, j, k)) )

         if (dom%has_dir(ydim)) rel_grad2 = sumsq01( rel_grad2, sumsq01( &
              rel_grad_1pair(p3d(i, j, k), p3d(i, j+dom%D_y, k)), &
              rel_grad_1pair(p3d(i, j, k), p3d(i, j-dom%D_y, k)) ) )

         if (dom%has_dir(zdim)) rel_grad2 = sumsq01( rel_grad2, sumsq01( &
              rel_grad_1pair(p3d(i, j, k), p3d(i, j, k+dom%D_z)), &
              rel_grad_1pair(p3d(i, j, k), p3d(i, j, k-dom%D_z)) ) )

      end function rel_grad2

      !>
      !! \brief Return 'relative gradient' defined as |a-b|/(|a| + |b|)
      !!
      !! \details This routine should return value between 0 and 1.
      !! rel_grad_1pair == 0 when a == b, also when a == b == 0
      !! rel_grad_1pair == 1 when one argument ==0 and the other /= 0
      !! rel_grad_1pair == 1 when a*b < 0
      !<

      elemental real function rel_grad_1pair(a, b)

         use func, only: operator(.notequals.)

         implicit none

         real, intent(in) :: a, b

         if ((abs(a) .notequals. 0.) .and. (abs(b) .notequals. 0.)) then
            rel_grad_1pair = abs(a-b) / (abs(a) + abs(b))
         else
            rel_grad_1pair = 0.
         endif

      end function rel_grad_1pair

      !>
      !! \brief square of an 'addition' of two numbers from [0, 1] range.
      !!
      !! \details Let us define a+b so it obeys the following rules:
      !! a (+) b = b (+) a
      !! a (+) 0 = a
      !! a (+) 1 = 1
      !! if a << 1 and b << 1, a (+) b \simeq sqrt(a**2+ b**2)
      !!
      !! One of possible formulas is:
      !! a (+) b = sqrt(a**2 - a**2 * b**2 + b**2)
      !! Note that: (a (+) b) (+) c = a (+) (b (+) c) = sqrt( (a (+) b)**2 * (1 - c**2) + c**2)
      !! That's why this routine returns (a (+) b)**2 for best performance when 'adding' multiple numbers.
      !<

      elemental real function sumsq01(a, b)

         implicit none

         real, intent(in) :: a, b

         sumsq01 = a**2 * (1-b**2) + b**2

      end function sumsq01

   end subroutine refine_on_relative_gradient

end module refinement_filters
