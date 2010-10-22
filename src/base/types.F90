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
   & problem_customize_solution, finalize_problem, domlen, idlen

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

   type :: tsl_container
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
