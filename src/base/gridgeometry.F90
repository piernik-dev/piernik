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
   public   :: gc, init_geometry, cleanup_geometry, set_geo_coeffs, geometry_source_terms

   real, dimension(:,:,:), pointer :: gc            !< array of geometrical coefficients
   integer, parameter              :: ngc = 3       !< number of geometrical coefficients

   real, allocatable, dimension(:,:,:), target  :: gc_xdim !< array of geometrical coefficients in x-direction
   real, allocatable, dimension(:,:,:), target  :: gc_ydim !< array of geometrical coefficients in y-direction
   real, allocatable, dimension(:,:,:), target  :: gc_zdim !< array of geometrical coefficients in z-direction

   interface
      !>
      !! \brief interface for routine setting grid coefficients
      !<
      subroutine set_gc(sweep,nvar,i1,i2)
         use types,         only: var_numbers
         implicit none
         character(len=*), intent(in)  :: sweep         !< direction (x, y or z) we are doing calculations for
         type(var_numbers), intent(in) :: nvar          !< \copydoc fluidindex::nvar
         integer, intent(in)           :: i1            !< coordinate of sweep in the 1st remaining direction
         integer, intent(in)           :: i2            !< coordinate of sweep in the 2st remaining direction
      end subroutine set_gc

      !>
      !! \brief interface for routine returning grid dependent source terms
      !!
      !! Currently, gsrc function returns accelerations
      !<
      function gsrc(u,p,sweep) result(res)
         implicit none
         character(len=*), intent(in)           :: sweep !< direction (x, y or z) we are doing calculations for
         real, dimension(:,:), intent(in)       :: u     !< sweep of fluid conservative variables
         real, dimension(:,:), intent(in)       :: p     !< sweep of pressure
         real, dimension(size(p,1),size(p,2))   :: res   !< output sweep of accelerations
      end function gsrc
   end interface

   procedure(set_gc), pointer :: set_geo_coeffs          !< generic pointer for routine setting geometrical coefficients
   procedure(gsrc),   pointer :: geometry_source_terms   !< generic pointer for routine calculating source terms
contains
!>
!! \brief Generic routine for module initialization
!<
   subroutine init_geometry
      use grid, only: geometry
      implicit none

      call set_geometry(geometry)
   end subroutine init_geometry

!>
!! \brief Generic routine for module cleanup
!!
!! Currently, deallocates auxiliary arrays for geometrical coefficients
!<
   subroutine cleanup_geometry
      implicit none

      if (allocated(gc_xdim)) deallocate(gc_xdim)
      if (allocated(gc_ydim)) deallocate(gc_ydim)
      if (allocated(gc_zdim)) deallocate(gc_zdim)

   end subroutine cleanup_geometry

!>
!! \brief Routine associating generic pointers
!!
!! Routine associates gridgeometry::set_geo_coeffs() and gridgeometry::geometry_source_terms()
!<
   subroutine set_geometry(geometry)
      use dataio_pub,    only: die, msg
      implicit none
      character(len=*), intent(in) :: geometry  !< string denoting geometry, "cartesian" or "cylindrical"

      select case (geometry)
         case ("cartesian")
            set_geo_coeffs          => set_cart_coeffs
            geometry_source_terms   => cart_geometry_source_terms
         case ("cylindrical")
            set_geo_coeffs          => set_cyl_coeffs
            geometry_source_terms   => cyl_geometry_source_terms
         case default
            write(msg,'(2a)') "[gridgeometry:set_geometry] Unknown system of coordinates ", geometry
            call die(msg)
      end select
   end subroutine set_geometry
!>
!! \brief Routine allocating auxiliary arrays
!!
!! We need to allocate those arrays later to have nvar%all
!<
   subroutine geo_coeffs_arrays(nvar)

      use types,         only: var_numbers
      use dataio_pub,    only: die
      use grid,          only: cg

      implicit none

      type(var_numbers), intent(in) :: nvar

      if ( any( [allocated(gc_xdim), allocated(gc_ydim), allocated(gc_zdim)] ) ) then
         call die("[gridgeometry:geo_coeffs_arrays] double allocation")
      else
         allocate(gc_xdim(ngc,nvar%all, cg%nx))
         allocate(gc_ydim(ngc,nvar%all, cg%ny))
         allocate(gc_zdim(ngc,nvar%all, cg%nz))
      endif

   end subroutine geo_coeffs_arrays
!>
!! \brief routine setting geometrical coefficients for cartesian grid
!<
   subroutine set_cart_coeffs(sweep,nvar,i1,i2)
      use types,         only: var_numbers
      use dataio_pub,    only: die, msg
      implicit none
      character(len=*), intent(in)  :: sweep
      type(var_numbers), intent(in) :: nvar
      integer, intent(in)           :: i1, i2
      logical, save                 :: frun = .true.

      if (frun) then
         call geo_coeffs_arrays(nvar)
         gc_xdim = 1.0
         gc_ydim = 1.0
         gc_zdim = 1.0
         frun = .false.
      endif

      select case (sweep)
         case ("xsweep")
            gc => gc_xdim
         case ("ysweep")
            gc => gc_ydim
         case ("zsweep")
            gc => gc_zdim
         case default
            write(msg,'(2a)') "[gridgeometry:set_cart_coeffs] Unknown sweep : ",sweep
            call die(msg)
            write(6,*) i1,i2   ! fool the compiler    QA_WARN
      end select
      frun = .false.
   end subroutine set_cart_coeffs
!>
!! \brief routine setting geometrical coefficients for cylindrical grid
!<
   subroutine set_cyl_coeffs(sweep,nvar,i1,i2)

      use types,         only: var_numbers
      use dataio_pub,    only: die, msg
      use grid,          only: cg

      implicit none

      character(len=*), intent(in)  :: sweep
      type(var_numbers), intent(in) :: nvar
      integer, intent(in)           :: i1, i2
      integer                       :: i
      logical, save                 :: frun = .true.

      if (frun) then
         call geo_coeffs_arrays(nvar)

         gc_xdim(1,:,:) = spread( cg%inv_x(:), 1, nvar%all)
         gc_xdim(2,:,:) = spread( cg%xr(:), 1, nvar%all)
         gc_xdim(3,:,:) = spread( cg%xl(:), 1, nvar%all)

         do i = LBOUND(nvar%all_fluids,1), UBOUND(nvar%all_fluids,1)
            gc_xdim(1,nvar%all_fluids(i)%imy,:) = gc_xdim(1,nvar%all_fluids(i)%imy,:) * cg%inv_x(:)
            gc_xdim(2,nvar%all_fluids(i)%imy,:) = gc_xdim(2,nvar%all_fluids(i)%imy,:) * cg%xr(:)
            gc_xdim(3,nvar%all_fluids(i)%imy,:) = gc_xdim(3,nvar%all_fluids(i)%imy,:) * cg%xl(:)
         enddo
         gc_ydim(2:3,:,:) = 1.0     ! [ 1/r , 1 , 1]

         gc_zdim = 1.0              ! [ 1, 1, 1]

         frun = .false.
      endif

      select case (sweep)
         case ("xsweep")
            gc => gc_xdim
         case ("ysweep")
            gc_ydim(1,:,:)   = cg%inv_x(i2)
            gc => gc_ydim
         case ("zsweep")
            gc => gc_zdim
         case default
            write(msg,'(2a)') "[gridgeometry:set_cyl_coeffs] Unknown sweep : ",sweep
            call die(msg)
            write(6,*) i1,i2   ! fool the compiler    QA_WARN
      end select
      frun = .false.

   end subroutine set_cyl_coeffs
!>
!! \brief routine calculating geometrical source term for cartesian grid
!<
   function cart_geometry_source_terms(u,p,sweep) result(res)

      implicit none

      character(len=*), intent(in)           :: sweep
      real, dimension(:,:), intent(in)       :: u, p
      real, dimension(size(p,1),size(p,2))   :: res

      res = 0.0
      return
      if (.false.) write(0,*) sweep, u, p

   end function cart_geometry_source_terms
!>
!! \brief routine calculating geometrical source term for cylindrical grid
!<
   function cyl_geometry_source_terms(u,p,sweep) result(res)

      use fluidindex, only: iarr_all_dn, iarr_all_my
      use grid,       only: cg

      implicit none

      character(len=*), intent(in)           :: sweep
      real, dimension(:,:), intent(in)       :: u, p
      real, dimension(size(p,1),size(p,2))   :: res
      integer :: i

      if (sweep == 'xsweep') then
         do i = 1, SIZE(iarr_all_dn)
            res(i,:) = cg%inv_x(:) * (u(iarr_all_my(i),:)**2 / u(iarr_all_dn(i),:) + p(i,:))
         enddo
      else
         res = 0.0
      endif

   end function cyl_geometry_source_terms

end module gridgeometry
