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
   implicit none
   private
   public :: hdf, value, grid_container, tsl_container, phys_prop, &
   & problem_customize_solution, finalize_problem, domlen, idlen, &
   & component_fluid, var_numbers, custom_emf_bnd, cleanup_problem, &
   & problem_grace_passed

   integer, parameter :: domlen = 16 ! should be <= mpisetup::cbuff_len
   integer, parameter :: idlen  = 3
   integer, parameter :: dims   = 3

   type :: hdf
      integer :: nhdf, nres, step_hdf, step_res, nstep, nrestart
      real    :: last_hdf_time, next_t_tsl,  next_t_log
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

      real    :: dx                             !< length of the %grid cell in x-direction
      real    :: dy                             !< length of the %grid cell in y-direction
      real    :: dz                             !< length of the %grid cell in z-direction
      real    :: idx                            !< inverted length of the %grid cell in x-direction
      real    :: idy                            !< inverted length of the %grid cell in y-direction
      real    :: idz                            !< inverted length of the %grid cell in z-direction
      real    :: dxmn                           !< the smallest length of the %grid cell (among dx, dy, and dz)
      real    :: dvol                           !< volume of one %grid cell
      real    :: xmin                           !< physical domain left x-boundary position
      real    :: xmax                           !< physical domain right x-boundary position
      real    :: ymin                           !< physical domain left y-boundary position
      real    :: ymax                           !< physical domain right y-boundary position
      real    :: zmin                           !< physical domain left z-boundary position
      real    :: zmax                           !< physical domain right z-boundary position
      real    :: xminb                          !< current block left x-boundary position
      real    :: xmaxb                          !< current block right x-boundary position
      real    :: yminb                          !< current block left y-boundary position
      real    :: ymaxb                          !< current block right y-boundary position
      real    :: zminb                          !< current block left z-boundary position
      real    :: zmaxb                          !< current block right z-boundary position

      real    :: Lx                             !< span of the physical domain in x-direction (xmax-xmin)
      real    :: Ly                             !< span of the physical domain in y-direction (ymax-ymin)
      real    :: Lz                             !< span of the physical domain in z-direction (zmax-zmin)
      real    :: Vol                            !< total volume of the physical domain

      real, allocatable, dimension(:) :: dl     !< array of %grid cell sizes in all directions
      real, allocatable, dimension(:) :: idl    !< array of inverted %grid cell sizes in all directions
      real, allocatable, dimension(:) :: x      !< array of x-positions of %grid cells centers
      real, allocatable, dimension(:) :: inv_x  !< array of invert x-positions of %grid cells centers
      real, allocatable, dimension(:) :: y      !< array of y-positions of %grid cells centers
      real, allocatable, dimension(:) :: inv_y  !< array of invert y-positions of %grid cells centers
      real, allocatable, dimension(:) :: z      !< array of z-positions of %grid cells centers
      real, allocatable, dimension(:) :: inv_z  !< array of invert z-positions of %grid cells centers
      real, allocatable, dimension(:) :: xl     !< array of x-positions of %grid cells left borders
      real, allocatable, dimension(:) :: yl     !< array of y-positions of %grid cells left borders
      real, allocatable, dimension(:) :: zl     !< array of z-positions of %grid cells left borders
      real, allocatable, dimension(:) :: xr     !< array of x-positions of %grid cells right borders
      real, allocatable, dimension(:) :: yr     !< array of y-positions of %grid cells right borders
      real, allocatable, dimension(:) :: zr     !< array of z-positions of %grid cells right borders

      integer :: nx                             !< number of %grid cells in one block in x-direction
      integer :: ny                             !< number of %grid cells in one block in y-direction
      integer :: nz                             !< number of %grid cells in one block in z-direction
      integer :: nxb                            !< number of physical domain %grid cells in one block (without boundary cells) in x-direction
      integer :: nyb                            !< number of physical domain %grid cells in one block (without boundary cells) in y-direction
      integer :: nzb                            !< number of physical domain %grid cells in one block (without boundary cells) in z-direction
      integer :: nxt                            !< total number of %grid cells in the whole domain in x-direction
      integer :: nyt                            !< total number of %grid cells in the whole domain in y-direction
      integer :: nzt                            !< total number of %grid cells in the whole domain in z-direction
      integer :: is                             !< index of the first %grid cell of physical domain in x-direction
      integer :: ie                             !< index of the last %grid cell of physical domain in x-direction
      integer :: js                             !< index of the first %grid cell of physical domain in y-direction
      integer :: je                             !< index of the last %grid cell of physical domain in y-direction
      integer :: ks                             !< index of the first %grid cell of physical domain in z-direction
      integer :: ke                             !< index of the last %grid cell of physical domain in z-direction
      integer :: nb                             !< number of boundary cells surrounding the physical domain, same for all directions
      integer :: isb, ieb, jsb, jeb, ksb, keb   !< auxiliary indices for exchanging boundary data, (e.g. is:isb -> ie+1:nx, ieb:ie -> 1:nb)
      integer :: maxxyz                         !< maximum number of %grid cells in any direction

      ! a pointer to the next cg
      ! grid level
      ! a list of MPI types for communication with targets

   end type grid_container

   type :: component
      integer :: all = 0   !< number of all variables in fluid/component
      integer :: beg = 0   !< beginning number of variables in fluid/component
      integer :: end = 0   !< end number of variables in fluid/component
      integer :: pos = 0   !< index denoting position of the fluid in the row of fluids
   end type component

   type, extends(component) :: component_fluid

      integer :: idn = -1      !< index denoting position of the fluid density in array arrays::u
      integer :: imx = -1      !< index denoting position of the fluid x-momentum in array arrays::u
      integer :: imy = -1      !< index denoting position of the fluid y-momentum in array arrays::u
      integer :: imz = -1      !< index denoting position of the fluid z-momentum in array arrays::u
      integer :: ien = -1      !< index denoting position of the fluid energy in array arrays::u

      real    :: cs    = 0.0   !< fluid's isothermal sound speed
      real    :: cs2   = 0.0   !< fluid's isothermal sound speed squared
      real    :: gam   = -1.0  !< fluid's adiabatic index
      real    :: gam_1 = -1.0  !< fluid's adiabatic index minus one

      logical :: is_selfgrav   = .false. !< True if fluid is selfgravitating
      logical :: is_magnetized = .false. !< True if fluid is magnetized
      logical :: has_energy    = .false. !< True if fluid has additional energy array

      character(len=idlen) :: tag = '' !< Human readable tag describing fluid

      integer, allocatable, dimension(:)  :: iarr
      integer, allocatable, dimension(:)  :: iarr_swpx
      integer, allocatable, dimension(:)  :: iarr_swpy
      integer, allocatable, dimension(:)  :: iarr_swpz

      type(phys_prop) :: snap
   end type component_fluid

   type :: var_numbers

      integer :: all         = 0      !< total number of fluid variables = the size of array \a u(:,:,;,:) in the first index
      integer :: fluids      = 0      !< number of fluids (ionized gas, neutral gas, dust)
      integer :: energ       = 0      !< number of non-isothermal fluids (indicating the presence of energy density in the vector of conservative variables)
      integer :: components  = 0      !< number of components, such as CRs, tracers, magnetic helicity (in future), whose formal description does not involve [???]
      integer :: fluids_sg   = 0      !< number of selfgravitating fluids (ionized gas, neutral gas, dust)

      type(component_fluid), dimension(:), pointer :: all_fluids

      type(component_fluid), pointer :: ion         !< numbers of variables for the ionized fluid
      type(component_fluid), pointer :: neu         !< numbers of variables for the neutral fluid
      type(component_fluid), pointer :: dst         !< numbers of variables for the dust fluid

      !> \todo those vars should be converted to pointers
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
   procedure(no_args),  pointer :: problem_grace_passed       => NULL()
   procedure(no_args),  pointer :: finalize_problem           => NULL()
   procedure(no_args),  pointer :: cleanup_problem            => NULL()
   procedure(tab_args), pointer :: custom_emf_bnd             => NULL()

end module types
