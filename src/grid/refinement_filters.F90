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
   public :: refine_on_gradient, refine_on_relative_gradient, refine_on_second_derivative, ref_crit

   type :: ref_crit
      integer              :: iv                      !< field index in cg%q or cg%w array
      integer              :: ic                      !< component index (cg%w(iv)%arr(ic,:,:,:)) or INVALID for 3D arrays
      integer              :: iplot                   !< field index for storing the refinement criterion or INVALID
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
!! \brief R. Loechner criterion
!! Original paper: https://www.researchgate.net/publication/222452974_An_adaptive_finite_element_scheme_for_transient_problems_in_CFD
!! Cartesian grid implementation: http://flash.uchicago.edu/~jbgallag/2012/flash4_ug/node14.html#SECTION05163100000000000000 (note that some indices in the left part of denominator seem to be slightly messed up)
!!
!! this%aux is the noise filter (epsilon in Loechner's paper)
!<

   subroutine refine_on_second_derivative(this, cg, p3d)

      use constants,  only: xdim, ydim, zdim, GEO_XYZ, INVALID
      use dataio_pub, only: die
      use domain,     only: dom
      use grid_cont,  only: grid_container

      implicit none

      class(ref_crit),                 intent(in)    :: this !< this contains refinement parameters
      type(grid_container), pointer,   intent(inout) :: cg   !< current grid piece
      real, dimension(:,:,:), pointer, intent(in)    :: p3d  !< pointer to array to be examined for (de)refinement needs (should contain at least two layers of updated guardcells)

      integer :: i, j, k
      real :: sn, sd, r, max_r
      integer, parameter :: how_far = 2

      if (dom%geometry_type /= GEO_XYZ) call die("[refinement_filters:refine_on_second_derivative] noncartesian geometry not supported yet")
      if (dom%nb < how_far+1) call die("[refinement_filters:refine_on_second_derivative] at east 2 guardcells are required")

      ! Perhaps it will be a bit faster with arrays for storing first-order differences
      ! but let's see if it works first and then how expensive it is.

      max_r = 0.
      do k = cg%ks - how_far*dom%D_z, cg%ke + how_far*dom%D_z
         do j = cg%js - how_far*dom%D_y, cg%je + how_far*dom%D_y
            do i = cg%is - how_far*dom%D_x, cg%ie + how_far*dom%D_x

               sn = 0.
               sd = 0.

               if (dom%has_dir(xdim)) then
                  ! d/dx d/dx
                  sn = sn + ((dux(i+1, j, k) - dux(i-1, j, k)) * cg%idx)**2
                  sd = sd + ((abs(dux(i+1, j, k)) + abs(dux(i-1, j, k)) + this%aux * (daux(i+1, j, k) + daux(i-1, j, k))) * cg%idx)**2
                  if (dom%has_dir(ydim)) then
                     ! d/dy d/dx
                     sn = sn + 2 * ((dux(i, j+1, k) - dux(i, j-1, k)) * cg%idy)**2  ! == d/dx d/dy
                     sd = sd + ((abs(dux(i, j+1, k)) + abs(dux(i, j-1, k)) + this%aux * (daux(i, j+1, k) + daux(i, j-1, k))) * cg%idy)**2 + &
                          &    ((abs(duy(i+1, j, k)) + abs(duy(i-1, j, k)) + this%aux * (dauy(i+1, j, k) + dauy(i-1, j, k))) * cg%idx)**2     ! d/dx d/dy
                     ! to be exploited: (dauy(i+1, j, k) + dauy(i-1, j, k)) * cg%idx == (daux(i, j+1, k) + daux(i, j-1, k)) * cg%idy
                  endif
                  if (dom%has_dir(zdim)) then
                     ! d/dz d/dx
                     sn = sn + 2 * ((dux(i, j, k+1) - dux(i, j, k-1)) * cg%idz)**2  ! == d/dx d/dz
                     sd = sd + ((abs(dux(i, j, k+1)) + abs(dux(i, j, k-1)) + this%aux * (daux(i, j, k+1) + daux(i, j, k-1))) * cg%idz)**2 + &
                          &    ((abs(duz(i+1, j, k)) + abs(duz(i-1, j, k)) + this%aux * (dauz(i+1, j, k) + dauz(i-1, j, k))) * cg%idx)**2     ! d/dx d/dz
                     ! to be exploited: dauz(i+1, j, k) + dauz(i-1, j, k)) * cg%idx == (daux(i, j, k+1) + daux(i, j, k-1)) * cg%idz
                  endif
               endif

               if (dom%has_dir(ydim)) then
                  ! d/dy d/dy
                  sn = sn + ((duy(i, j+1, k) - duy(i, j-1, k)) * cg%idy)**2
                  sd = sd + ((abs(duy(i, j+1, k)) + abs(duy(i, j-1, k)) + this%aux * (dauy(i, j+1, k) + dauy(i, j-1, k))) * cg%idy)**2
                  if (dom%has_dir(zdim)) then
                     ! d/dz d/dy
                     sn = sn + 2 * ((duy(i, j, k+1) - duy(i, j, k-1)) * cg%idz)**2  ! == d/dy d/dz
                     sd = sd + ((abs(duy(i, j, k+1)) + abs(duy(i, j, k-1)) + this%aux * (dauy(i, j, k+1) + dauy(i, j, k-1))) * cg%idz)**2 + &
                          &    ((abs(duz(i, j+1, k)) + abs(duz(i, j-1, k)) + this%aux * (dauz(i, j+1, k) + dauz(i, j-1, k))) * cg%idy)**2     ! d/dy d/dz
                     ! to be exploited: (dauz(i, j+1, k) + dauz(i, j-1, k)) * cg%idy == (dauy(i, j, k+1) + dauy(i, j, k-1)) * cg%idz
                  endif
               endif

               if (dom%has_dir(zdim)) then
                  ! d/dz d/dz
                  sn = sn + ((duz(i, j, k+1) - duz(i, j, k-1)) * cg%idz)**2
                  sd = sd + ((abs(duz(i, j, k+1)) + abs(duz(i, j, k-1)) + this%aux * (dauz(i, j, k+1) + dauz(i, j, k-1))) * cg%idz)**2
               endif

               if (sd > 0.) then
                  r = sn / sd
               else  ! sd == 0 because it should never be < 0.
                  if (sn > 0.) then
                     r = 1.  ! strange case
                  else
                     r = 0.  ! most likely constant field == 0.
                  endif
               endif
               if (this%iplot /= INVALID) cg%q(this%iplot)%arr(i, j, k) = r
               max_r = max(max_r, r)

               if (i >= cg%is .and. i <= cg%ie .and. j >= cg%js .and. j <= cg%je .and. k >= cg%ks .and. k <= cg%ke) &
                    cg%refinemap(i, j, k) = cg%refinemap(i, j, k) .or. (r >= this%ref_thr)

            enddo
         enddo
      enddo

      cg%refine_flags%derefine = cg%refine_flags%derefine .or. (max_r < this%deref_thr)

   contains

      elemental real function dux(i, j, k)
         implicit none
         integer, intent(in) :: i, j, k !< (x, y, z)-indices
         dux = (p3d(i+1, j, k) - p3d(i-1, j, k)) * cg%idx
      end function dux

      elemental real function daux(i, j, k)
         implicit none
         integer, intent(in) :: i, j, k !< (x, y, z)-indices
         daux = (abs(p3d(i+1, j, k)) + abs(p3d(i-1, j, k))) * cg%idx
      end function daux

      elemental real function duy(i, j, k)
         implicit none
         integer, intent(in) :: i, j, k !< (x, y, z)-indices
         duy = (p3d(i, j+1, k) - p3d(i, j-1, k)) * cg%idy
      end function duy

      elemental real function dauy(i, j, k)
         implicit none
         integer, intent(in) :: i, j, k !< (x, y, z)-indices
         dauy = (abs(p3d(i, j+1, k)) + abs(p3d(i, j-1, k))) * cg%idy
      end function dauy

      elemental real function duz(i, j, k)
         implicit none
         integer, intent(in) :: i, j, k !< (x, y, z)-indices
         duz = (p3d(i, j, k+1) - p3d(i, j, k-1)) * cg%idz
      end function duz

      elemental real function dauz(i, j, k)
         implicit none
         integer, intent(in) :: i, j, k !< (x, y, z)-indices
         dauz = (abs(p3d(i, j, k+1)) + abs(p3d(i, j, k-1))) * cg%idz
      end function dauz

   end subroutine refine_on_second_derivative

!>
!! \brief Refine/derefine based on ||grad u||
!! This is sensitive to gradients, but the thresholds must be rescaled, when you change units of the problem.
!<

   subroutine refine_on_gradient(this, cg, p3d)

      use constants, only: INVALID
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
               if (this%iplot /= INVALID) cg%q(this%iplot)%arr(i, j, k) = r
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

         integer(kind=4), intent(in) :: i1, i2, j1, j2, k1, k2

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

         integer, intent(in) :: i !< x-index
         integer, intent(in) :: j !< y-index
         integer, intent(in) :: k !< z-index

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

      use constants, only: INVALID
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
               if (this%iplot /= INVALID) cg%q(this%iplot)%arr(i, j, k) = r
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

         integer(kind=4), intent(in) :: i1, i2, j1, j2, k1, k2

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

         integer, intent(in) :: i !< x-index
         integer, intent(in) :: j !< y-index
         integer, intent(in) :: k !< z-index

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
