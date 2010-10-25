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
!! \brief (KK)
!<
module types
   private
   public :: indx, hdf, value, grid_container, tsl_container, phys_prop, &
   & problem_customize_solution, finalize_problem, domlen, idlen, &
   & component_fluid, var_numbers

   integer, parameter :: domlen = 16 ! should be <= mpisetup::cbuff_len
   integer, parameter :: idlen  = 3
   integer, parameter :: dims   = 3

   type :: indx
      integer :: dnd = -1, dnn = -1, dni = -1
      integer :: mxd = -1, mxn = -1, mxi = -1
      integer :: myd = -1, myn = -1, myi = -1
      integer :: mzd = -1, mzn = -1, mzi = -1
      integer :: enn = -1, eni = -1
      integer :: bx = -1, by = -1, bz = -1
      integer, dimension(:), pointer :: arr_crs
   end type indx

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
      integer :: nxd, nyd, nzd, nb
      integer :: nx, ny, nz
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

   type :: component
      integer :: all = 0              !< number of all variables in fluid/component
      integer :: beg = 0              !< beginning number of variables in fluid/component
      integer :: end = 0              !< end number of variables in fluid/component
      integer :: pos = 0              !< index denoting position of the fluid in the row of fluids
   end type component

   type :: component_fluid
      integer :: all = 0              !< number of all variables in fluid/component
      integer :: beg = 0              !< beginning number of variables in fluid/component
      integer :: end = 0              !< end number of variables in fluid/component
      integer :: pos = 0              !< index denoting position of the fluid in the row of fluids
      integer :: idn = -1
      integer :: imx = -1
      integer :: imy = -1
      integer :: imz = -1
      integer :: ien = -1

      real    :: cs    = 0.0
      real    :: cs2   = 0.0
      real    :: gam   = -1.0
      real    :: gam_1 = -1.0

      logical :: sg  = .false.

      character(len=idlen) :: tag = ''

      integer, allocatable, dimension(:)  :: iarr
      integer, allocatable, dimension(:)  :: iarr_swpx
      integer, allocatable, dimension(:)  :: iarr_swpy
      integer, allocatable, dimension(:)  :: iarr_swpz

      type(phys_prop) :: snap
   end type component_fluid

   type :: var_numbers
      integer :: all         = 0     !< total number of fluid variables = the size of array \a u(:,:,;,:) in the first index
      integer :: fluids      = 0     !< number of fluids (ionized gas, neutral gas, dust)
      integer :: adiab       = 0     !< number of adiabatic fluids (indicating the presence of energy density in the vector of conservative variables)
      integer :: components  = 0     !< number of components, such as CRs, tracers, magnetic helicity (in future), whose formal description does not involve [???]
      integer :: fluids_sg   = 0     !< number of selfgravitating fluids (ionized gas, neutral gas, dust)

      type(component_fluid), dimension(:), pointer :: all_fluids

      type(component_fluid), pointer :: ion         !< numbers of variables for the ionized fluid
      type(component_fluid), pointer :: neu         !< numbers of variables for the neutral fluid
      type(component_fluid), pointer :: dst         !< numbers of variables for the dust fluid

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
   end type tsl_container

   interface
      subroutine no_args
         implicit none
      end subroutine no_args
   end interface

   procedure(no_args), pointer :: problem_customize_solution => NULL()
   procedure(no_args), pointer :: finalize_problem           => NULL()

end module types
