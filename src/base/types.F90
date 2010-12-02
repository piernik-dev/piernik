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
!! \brief (KK)
!<
module types
   private
   public :: hdf, value, grid_container, tsl_container, phys_prop, &
   & problem_customize_solution, finalize_problem, domlen, idlen, &
   & component_fluid, var_numbers, custom_emf_bnd

   integer, parameter :: domlen = 16 ! should be <= mpisetup::cbuff_len
   integer, parameter :: idlen  = 3
   integer, parameter :: dims   = 3

   type :: hdf
      integer :: nhdf, ntsl, nres, nlog, step_hdf, step_res, nstep, nrestart
      real    :: last_hdf_time
      character(len=domlen)  :: domain
      character(len=idlen)   :: new_id
   end type hdf

   type :: value
      real                     :: val
      integer, dimension(dims) :: loc
      integer                  :: proc
   end type value

   type :: phys_prop
      type(value) :: dens_min
      type(value) :: dens_max
      type(value) :: velx_max
      type(value) :: vely_max
      type(value) :: velz_max
      type(value) :: pres_max
      type(value) :: pres_min
      type(value) :: temp_max
      type(value) :: temp_min
      type(value) :: cs_max
      real        :: c
      real        :: dt
   end type phys_prop

   type :: grid_container
      real    :: dx, dy, dz, dxmn, dvol
      integer :: nb, nx, ny, nz
      integer :: nxb, nyb, nzb
      integer :: nxt, nyt, nzt
      integer :: is, ie, js, je, ks, ke
      integer :: maxxyz

      real    :: xmin, xmax, ymin, ymax, zmin, zmax
      real    :: xminb, xmaxb, yminb, ymaxb, zminb, zmaxb
      real    :: Lx, Ly, Lz
      integer, pointer  :: xdim, ydim, zdim

      real, dimension(:), pointer :: dl
      real, dimension(:), pointer  :: x, xl, xr
      real, dimension(:), pointer  :: y, yl, yr
      real, dimension(:), pointer  :: z, zl, zr
   end type grid_container

   ! ToDo: Fix qa2.py script rather than adding QA_WARN for type components

   type :: component
      !>
      !! number of all variables in fluid/component
      !<
      integer :: all = 0              ! QA_WARN implicit save does not work in types
      !>
      !! beginning number of variables in fluid/component
      !<
      integer :: beg = 0              ! QA_WARN implicit save does not work in types
      !>
      !! end number of variables in fluid/component
      !<
      integer :: end = 0              ! QA_WARN implicit save does not work in types
      !>
      !! index denoting position of the fluid in the row of fluids
      !<
      integer :: pos = 0              ! QA_WARN implicit save does not work in types
   end type component

   type, extends(component) :: component_fluid
      !>
      !! index denoting position of the fluid density in array arrays::u
      !<
      integer :: idn = -1              ! QA_WARN implicit save does not work in types
      !>
      !! index denoting position of the fluid x-momentum in array arrays::u
      !<
      integer :: imx = -1              ! QA_WARN implicit save does not work in types
      !>
      !! index denoting position of the fluid y-momentum in array arrays::u
      !<
      integer :: imy = -1              ! QA_WARN implicit save does not work in types
      !>
      !! index denoting position of the fluid z-momentum in array arrays::u
      !<
      integer :: imz = -1              ! QA_WARN implicit save does not work in types
      !>
      !! index denoting position of the fluid energy in array arrays::u
      !<
      integer :: ien = -1              ! QA_WARN implicit save does not work in types

      !>
      !! fluid's isothermal sound speed
      !<
      real    :: cs    = 0.0           ! QA_WARN implicit save does not work in types
      !>
      !! fluid's isothermal sound speed squared
      !<
      real    :: cs2   = 0.0           ! QA_WARN implicit save does not work in types
      !>
      !! fluid's adiabatic index
      !<
      real    :: gam   = -1.0          ! QA_WARN implicit save does not work in types
      !>
      !! fluid's adiabatic index minus one
      !<
      real    :: gam_1 = -1.0          ! QA_WARN implicit save does not work in types
      !>
      !! True if fluid is selfgravitating
      !<
      logical :: sg  = .false.         ! QA_WARN implicit save does not work in types
      !>
      !! Human readable tag describing fluid
      !<
      character(len=idlen) :: tag = '' ! QA_WARN implicit save does not work in types

      integer, allocatable, dimension(:)  :: iarr
      integer, allocatable, dimension(:)  :: iarr_swpx
      integer, allocatable, dimension(:)  :: iarr_swpy
      integer, allocatable, dimension(:)  :: iarr_swpz

      type(phys_prop) :: snap
   end type component_fluid

   type :: var_numbers
      !>
      !! total number of fluid variables = the size of array \a u(:,:,;,:) in the first index
      !<
      integer :: all         = 0     ! QA_WARN implicit save does not work in types
      !>
      !! number of fluids (ionized gas, neutral gas, dust)
      !<
      integer :: fluids      = 0     ! QA_WARN implicit save does not work in types
      !>
      !! number of adiabatic fluids (indicating the presence of energy density in the vector of conservative variables)
      !<
      integer :: adiab       = 0     ! QA_WARN implicit save does not work in types
      !>
      !! number of components, such as CRs, tracers, magnetic helicity (in future), whose formal description does not involve [???]
      !<
      integer :: components  = 0     ! QA_WARN implicit save does not work in types
      !>
      !! number of selfgravitating fluids (ionized gas, neutral gas, dust)
      !<
      integer :: fluids_sg   = 0     ! QA_WARN implicit save does not work in types

      type(component_fluid), dimension(:), pointer :: all_fluids

      type(component_fluid), pointer :: ion         !< numbers of variables for the ionized fluid
      type(component_fluid), pointer :: neu         !< numbers of variables for the neutral fluid
      type(component_fluid), pointer :: dst         !< numbers of variables for the dust fluid

      !! ToDo: those vars should be converted to pointers
      type(component) :: crs         !< numbers of variables in all cosmic ray components
      type(component) :: crn         !< numbers of variables in cosmic ray nuclear components
      type(component) :: cre         !< numbers of variables in cosmic ray electron components

   end type var_numbers

   type :: tsl_container
      logical :: dummy
#ifdef COSM_RAYS
      real :: encr_min, encr_max
#endif /* COSM_RAYS */

#ifdef RESISTIVE
      real :: etamax
#endif /* RESISTIVE */

#ifdef MAGNETIC
      real :: b_min, b_max, divb_max
#endif /* MAGNETIC */

#ifdef IONIZED
      real :: vai_max
#endif /* IONIZED */

#ifdef VARIABLE_GP
      real :: gpxmax, gpymax, gpzmax
#endif /* VARIABLE_GP */

   end type tsl_container

   interface
      subroutine no_args
         implicit none
      end subroutine no_args
      subroutine tab_args(tab)
         implicit none
         real, dimension(:,:,:), intent(inout) :: tab
      end subroutine tab_args
   end interface

   procedure(no_args),  pointer :: problem_customize_solution => NULL()
   procedure(no_args),  pointer :: finalize_problem           => NULL()
   procedure(tab_args), pointer :: custom_emf_bnd             => NULL()

end module types
