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
!! \brief Definitions of compound types and subroutine templates for user hooks
!<
module types

   use constants, only: domlen, idlen, bndlen, ndims

   implicit none

   private
   public :: hdf, axes, domain_container, tsl_container, value, problem_customize_solution, problem_grace_passed, finalize_problem, cleanup_problem, custom_emf_bnd

   type :: hdf
      integer :: nhdf, nres, step_hdf, step_res, nstep, nrestart
      real    :: last_hdf_time, next_t_tsl,  next_t_log
      character(len=domlen)  :: domain
      character(len=idlen)   :: new_id
   end type hdf

   type :: value
      real                      :: val
      integer, dimension(ndims) :: loc
      integer                   :: proc
   end type value

! AMR: There will be at least one domain container for the base grid. It will be possible to host one or more refined domains on the base container and on the refined containers.
!      The refined domains may cover whole parent domain or only a tiny part of it.
!      The multigrid solver will operate also on a stack of coarser domains - parents of the base domain. The coarser domains must be no smaller than the base domain.
   type :: domain_container
      ! primary parameters, read from /DOMAIN_SIZES/, /BOUNDARIES/ and /DOMAIN_LIMITS/ namelists
      real    :: xmin                           !< physical domain left x-boundary position
      real    :: xmax                           !< physical domain right x-boundary position
      real    :: ymin                           !< physical domain left y-boundary position
      real    :: ymax                           !< physical domain right y-boundary position
      real    :: zmin                           !< physical domain left z-boundary position
      real    :: zmax                           !< physical domain right z-boundary position
      integer, dimension(ndims) :: n_d          !< number of grid cells in physical domain in x-, y- and z-direction (where equal to 1, the dimension is reduced to a point with no boundary cells)
      integer                   :: nb           !< number of boundary cells surrounding the physical domain, same for all directions
      character(len=bndlen) :: bnd_xl_dom       !< low-x computational domain boundary
      character(len=bndlen) :: bnd_xr_dom       !< high-x computational domain boundary
      character(len=bndlen) :: bnd_yl_dom       !< low-y computational domain boundary
      character(len=bndlen) :: bnd_yr_dom       !< high-y computational domain boundary
      character(len=bndlen) :: bnd_zl_dom       !< low-z computational domain boundary
      character(len=bndlen) :: bnd_zr_dom       !< high-z computational domain boundary
      !\todo use integer, dimension(ndims, LO:HI) :: bnd , like cg%

      ! derived parameters
      real    :: Lx                             !< span of the physical domain in x-direction (xmax-xmin)
      real    :: Ly                             !< span of the physical domain in y-direction (ymax-ymin)
      real    :: Lz                             !< span of the physical domain in z-direction (zmax-zmin)
      real    :: x0                             !< center of the physical domain in x-direction (xmax+xmin)/2.
      real    :: y0                             !< center of the physical domain in y-direction (ymax+ymin)/2.
      real    :: z0                             !< center of the physical domain in z-direction (zmax+zmin)/2.
      real    :: Vol                            !< total volume of the physical domain

      ! Do not use n[xyz]t components in the Piernik source tree. Avoid using them in your problem-specific files because it seriuosly limits paralelization
      ! \todo move them to another type, that extends domain_container
      integer :: nxt                            !< total number of %grid cells in the whole domain in x-direction
      integer :: nyt                            !< total number of %grid cells in the whole domain in y-direction
      integer :: nzt                            !< total number of %grid cells in the whole domain in z-direction

      !! \todo add hooks to parent and a list of children domains
   end type domain_container

   type :: axes
      real, allocatable, dimension(:) :: x      !< array of x-positions of %grid cells centers
      real, allocatable, dimension(:) :: y      !< array of y-positions of %grid cells centers
      real, allocatable, dimension(:) :: z      !< array of z-positions of %grid cells centers
   end type axes

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

   ! User hooks

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
