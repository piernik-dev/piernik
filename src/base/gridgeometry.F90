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
#include "piernik.def"
!>
!! \brief Module that sets coordinate system
!! \date December 2010
!!
!! This module contains routines that calculates grid coefficients and source terms used during flux update in rtvd::relaxing_tvd
!!
!<
module gridgeometry

   implicit none

   private
   public   :: gc, GC1, GC2, GC3, init_geometry, set_geo_coeffs, geometry_source_terms

   real, dimension(:,:,:), pointer :: gc            !< array of geometrical coefficients, \todo move it to the grid container
   enum, bind(C)
      enumerator :: GC1, GC2, GC3                   !< \todo change to a some meaningful names
   end enum

   interface
      !>
      !! \brief interface for routine setting grid coefficients
      !<
      subroutine set_gc(sweep, flind, i1, i2, cg)

         use fluidtypes, only: var_numbers
         use grid_cont,  only: grid_container

         implicit none

         integer(kind=4),  intent(in)  :: sweep         !< direction (xdim, ydim or zdim) we are doing calculations for
         type(var_numbers), intent(in) :: flind         !< \copydoc fluidindex::flind
         integer, intent(in)           :: i1            !< coordinate of sweep in the 1st remaining direction
         integer, intent(in)           :: i2            !< coordinate of sweep in the 2st remaining direction
         type(grid_container), pointer :: cg            !< current grid container

      end subroutine set_gc

      !>
      !! \brief interface for routine returning grid dependent source terms
      !!
      !! Currently, gsrc function returns accelerations
      !<
      function gsrc(u, p, sweep, cg) result(res)

         use grid_cont, only: grid_container

         implicit none

         integer(kind=4), intent(in)      :: sweep !< direction (x, y or z) we are doing calculations for
         real, dimension(:,:), intent(in) :: u     !< sweep of fluid conservative variables
         real, dimension(:,:), intent(in) :: p     !< sweep of pressure
         type(grid_container), intent(in) :: cg    !< current grid container

         real, dimension(size(p,1),size(p,2)) :: res   !< output sweep of accelerations

      end function gsrc

   end interface

   procedure(set_gc), pointer :: set_geo_coeffs          !< generic pointer for routine setting geometrical coefficients
   procedure(gsrc),   pointer :: geometry_source_terms   !< generic pointer for routine calculating source terms

contains
!>
!! \brief Routine for module initialization
!!
!! \details Routine associates gridgeometry::set_geo_coeffs() and gridgeometry::geometry_source_terms() pointers.
!<
   subroutine init_geometry

      use dataio_pub, only: msg, die, code_progress
      use constants,  only: PIERNIK_INIT_DOMAIN, GEO_XYZ, GEO_RPZ
      use domain,     only: geometry_type

      implicit none

      if (code_progress < PIERNIK_INIT_DOMAIN) call die("[gridgeometry:init_geometry] domain not initialized")

      select case (geometry_type)
         case (GEO_XYZ)
            set_geo_coeffs          => set_cart_coeffs
            geometry_source_terms   => cart_geometry_source_terms
         case (GEO_RPZ)
            set_geo_coeffs          => set_cyl_coeffs
            geometry_source_terms   => cyl_geometry_source_terms
         case default
            write(msg,'(a,i3)') "[gridgeometry:init_geometry] Unknown system of coordinates ", geometry_type
            call die(msg)
      end select

   end subroutine init_geometry

!>
!! \brief Routine allocating auxiliary arrays
!!
!! We need to allocate those arrays later to have flind%all
!<
   subroutine geo_coeffs_arrays(flind, cg)

      use constants,  only: xdim, ydim, zdim
      use dataio_pub, only: die
      use fluidtypes, only: var_numbers
      use grid_cont,  only: grid_container

      implicit none

      type(var_numbers), intent(in) :: flind
      type(grid_container), pointer :: cg

      if ( any( [allocated(cg%gc_xdim), allocated(cg%gc_ydim), allocated(cg%gc_zdim)] ) ) then
         call die("[gridgeometry:geo_coeffs_arrays] double allocation")
      else
         allocate(cg%gc_xdim(GC1:GC3, flind%all, cg%n_(xdim)))
         allocate(cg%gc_ydim(GC1:GC3, flind%all, cg%n_(ydim)))
         allocate(cg%gc_zdim(GC1:GC3, flind%all, cg%n_(zdim)))
      endif

   end subroutine geo_coeffs_arrays
!>
!! \brief A dummy routine. There is no need to set any cartesian coefficients, because all of them are equal to 1. so we don't use them in rtvd
!<
   subroutine set_cart_coeffs(sweep, flind, i1, i2, cg)

      use fluidtypes, only: var_numbers
      use grid_cont,  only: grid_container

      implicit none

      integer(kind=4),   intent(in) :: sweep
      type(var_numbers), intent(in) :: flind
      integer, intent(in)           :: i1, i2
      type(grid_container), pointer :: cg

      return
      if (.false.) write(0,*) cg%u(flind%all, sweep, i1, i2) ! suppress compiler warnings

   end subroutine set_cart_coeffs
!>
!! \brief routine setting geometrical coefficients for cylindrical grid
!<
   subroutine set_cyl_coeffs(sweep, flind, i1, i2, cg)

      use constants,  only: xdim, ydim, zdim
      use dataio_pub, only: die, msg
      use domain,     only: is_multicg
      use fluidtypes, only: var_numbers
      use grid_cont,  only: grid_container

      implicit none

      integer(kind=4), intent(in)   :: sweep
      type(var_numbers), intent(in) :: flind
      integer, intent(in)           :: i1, i2
      type(grid_container), pointer :: cg

      integer                        :: i
      logical, save                  :: frun = .true.

      if (is_multicg) call die("[gridgeometry:set_cyl_coeffs] multiple grid pieces per procesor not implemented yet") ! move gc to grid_container%, fix initialization

      !> \todo This should be probably called from cg%init (beware of cyclic dependencies) or init_grid
      if (frun) then
         call geo_coeffs_arrays(flind, cg)

         cg%gc_xdim(GC1,:,:) = spread( cg%inv_x(:), 1, flind%all)
         cg%gc_xdim(GC2,:,:) = spread( cg%xr(:),    1, flind%all)
         cg%gc_xdim(GC3,:,:) = spread( cg%xl(:),    1, flind%all)

         do i = lbound(flind%all_fluids,1), ubound(flind%all_fluids,1)
            cg%gc_xdim(GC1, flind%all_fluids(i)%imy, :) = cg%gc_xdim(GC1, flind%all_fluids(i)%imy, :) * cg%inv_x(:)
            cg%gc_xdim(GC2, flind%all_fluids(i)%imy, :) = cg%gc_xdim(GC2, flind%all_fluids(i)%imy, :) * cg%xr(:)
            cg%gc_xdim(GC3, flind%all_fluids(i)%imy, :) = cg%gc_xdim(GC3, flind%all_fluids(i)%imy, :) * cg%xl(:)
         enddo
         cg%gc_ydim(GC2:GC3,:,:) = 1.0     ! [ 1/r , 1 , 1]

         cg%gc_zdim(:,:,:) = 1.0           ! [ 1, 1, 1]

         frun = .false.
      endif

      select case (sweep)
         case (xdim)
            gc => cg%gc_xdim
         case (ydim)
            cg%gc_ydim(GC1,:,:)   = cg%inv_x(i2)
            gc => cg%gc_ydim
         case (zdim)
            gc => cg%gc_zdim
         case default
            write(msg,'(a,i2)') "[gridgeometry:set_cyl_coeffs] Unknown sweep : ",sweep
            call die(msg)
            if (.false.) frun = (i1==i2) ! suppress compiler warnings
      end select
      frun = .false.

   end subroutine set_cyl_coeffs
!>
!! \brief routine calculating geometrical source term for cartesian grid
!<
   function cart_geometry_source_terms(u, p, sweep, cg) result(res)

      use grid_cont,  only: grid_container

      implicit none

      integer(kind=4), intent(in)      :: sweep
      real, dimension(:,:), intent(in) :: u, p
      type(grid_container), intent(in) :: cg

      real, dimension(size(p,1),size(p,2))   :: res

      res = 0.0
      return
      if (.false.) write(0,*) sweep, u, p, cg%inv_x(1)

   end function cart_geometry_source_terms
!>
!! \brief routine calculating geometrical source term for cylindrical grid
!<
   function cyl_geometry_source_terms(u, p, sweep, cg) result(res)

      use constants,  only: xdim
      use fluidindex, only: iarr_all_dn, iarr_all_my
      use grid_cont,  only: grid_container

      implicit none

      integer(kind=4), intent(in)      :: sweep
      real, dimension(:,:), intent(in) :: u, p
      type(grid_container), intent(in) :: cg

      real, dimension(size(p,1),size(p,2)) :: res
      integer :: i

      if (sweep == xdim) then
         do i = 1, size(iarr_all_dn)
            res(i,:) = cg%inv_x(:) * (u(iarr_all_my(i),:)**2 / u(iarr_all_dn(i),:) + p(i,:))
         enddo
      else
         ! Note that there's no additional source term for angular momentum since we're using
         ! modified equation of motion following  Mignone et al. (2007), ApJS 170:228- and
         ! Skinner & Ostriker (2010), ApJSS 188:290-311
         ! The conservative implementation uses the gc(:,:,:) array to modify the azimuthal momentum flux to mimic angular momentum flux.
         res = 0.0
      endif

   end function cyl_geometry_source_terms

end module gridgeometry
