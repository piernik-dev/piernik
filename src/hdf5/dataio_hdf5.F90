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

#define RNG cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke
#include "piernik.h"

!>
!! \brief (KK) Module that contains I/O routines using HDF5 library
!!
!! Modules contains routines for creating HDF5 output such as
!! plots, snapshots, restart files.
!!
!<
module dataio_hdf5
! pulled by ANY
   use dataio_pub,    only: maxparfilelen, maxparfilelines
   use list_hdf5,     only: S_LEN

   implicit none

   private
   public :: init_hdf5, read_restart_hdf5, cleanup_hdf5, write_hdf5, write_restart_hdf5, write_plot, write_3darr_to_restart, read_3darr_from_restart
   public :: parfile, parfilelines, maxparfilelines

   integer, parameter :: dnamelen=10
   character(len=dnamelen), dimension(2) :: dname = (/"fluid     ","mag       "/)  !< dataset names for restart files
   character(len=S_LEN), allocatable, dimension(:) :: hdf_vars  !< dataset names for hdf files
   integer :: nhdf_vars !< number of quantities plotted to hdf files
   integer :: ix !< no. of cell (1 <= ix < nxd) for YZ slice in plt files
   integer :: iy !< no. of cell (1 <= iy < nyd) for XZ slice in plt files
   integer :: iz !< no. of cell (1 <= iz < nzd) for XY slice in plt files
   real    :: dt_plt !< frequency of plt output

   ! storage for the problem.par
   character(len=maxparfilelen), dimension(maxparfilelines) :: parfile !< contents of the parameter file
   integer, save                             :: parfilelines = 0       !< number of lines in the parameter file

contains

!>
!! \brief Empty routine assigned to additional_attrs unless the problem provides something else (see maclaurin or wt4 for examples)
!<

   subroutine null_attrs(file_id)

      use hdf5, only: HID_T

      implicit none

      integer(HID_T), intent(in) :: file_id

      integer :: dummy ! just for suppressing compiler warning
      if (.false.) dummy = kind(file_id)

   end subroutine null_attrs

!>
!! \brief Procedure initializing HDF5 module
!<

   subroutine init_hdf5(vars,tix,tiy,tiz,tdt_plt)

      use fluidindex,    only: iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz
      use list_hdf5,     only: additional_attrs, problem_write_restart, problem_read_restart
      use dataio_pub,    only: varlen
#ifdef COSM_RAYS
      use fluidindex,    only: iarr_all_crs
      use dataio_pub,    only: warn, msg
#endif /* COSM_RAYS */

      implicit none

      character(len=varlen), dimension(:), intent(in) :: vars  !< quantities to be plotted, see dataio::vars
      integer, intent(in) :: tix     !< local copy of dataio::ix
      integer, intent(in) :: tiy     !< local copy of dataio::iy
      integer, intent(in) :: tiz     !< local copy of dataio::iz
      real, intent(in)    :: tdt_plt !< local copy of dataio::dt_plt

      integer :: nvars, i, j
#if defined COSM_RAYS && !defined NEW_HDF5
      integer :: k
      character(len=varlen) :: aux
#endif /* COSM_RAYS && !NEW_HDF5 */

      ix = tix; iy = tiy; iz = tiz; dt_plt = tdt_plt

      nvars = 1
      do while ( len(trim(vars(nvars))) > 1)
         nvars = nvars + 1
      enddo
      nvars = nvars - 1

      nhdf_vars = 0
      do i = 1, nvars
         select case (vars(i))
            case ('dens')
               nhdf_vars = nhdf_vars + size(iarr_all_dn,1)
            case ('velx')
               nhdf_vars = nhdf_vars + size(iarr_all_mx,1)
            case ('vely')
               nhdf_vars = nhdf_vars + size(iarr_all_my,1)
            case ('velz')
               nhdf_vars = nhdf_vars + size(iarr_all_mz,1)
            case ('ener')
               nhdf_vars = nhdf_vars + size(iarr_all_mz,1)
#ifdef DUST
               nhdf_vars = nhdf_vars - 1
#endif /* DUST */
#ifdef GRAV
            case ('gpot')
               nhdf_vars = nhdf_vars + 1
#ifdef MULTIGRID
            case ('mgso') ! multigrid solution
               nhdf_vars = nhdf_vars + 1
#endif /* MULTIGRID */
#endif /* GRAV */
            case ('magx')
               nhdf_vars = nhdf_vars + 1
            case ('magy')
               nhdf_vars = nhdf_vars + 1
            case ('magz')
               nhdf_vars = nhdf_vars + 1
#ifdef COSM_RAYS
            case ('encr')
#ifndef NEW_HDF5
               nhdf_vars = nhdf_vars + size(iarr_all_crs,1)
#endif /* !NEW_HDF5 */
#endif /* COSM_RAYS */
            case ('pres')
               nhdf_vars = nhdf_vars + 1
            case default
               nhdf_vars = nhdf_vars + 1
         end select
      enddo
      allocate(hdf_vars(nhdf_vars)); j = 1
      do i = 1, nvars
         select case (vars(i))
            case ('dens')
#ifdef DUST
               hdf_vars(j) = 'dend' ; j = j + 1
#endif /* DUST */
#ifdef NEUTRAL
               hdf_vars(j) = 'denn' ; j = j + 1
#endif /* NEUTRAL */
#ifdef IONIZED
               hdf_vars(j) = 'deni' ; j = j + 1
#endif /* IONIZED */
            case ('velx')
#ifdef DUST
               hdf_vars(j) = 'vlxd' ; j = j + 1
#endif /* DUST */
#ifdef NEUTRAL
               hdf_vars(j) = 'vlxn' ; j = j + 1
#endif /* NEUTRAL */
#ifdef IONIZED
               hdf_vars(j) = 'vlxi' ; j = j + 1
#endif /* IONIZED */
            case ('vely')
#ifdef DUST
               hdf_vars(j) = 'vlyd' ; j = j + 1
#endif /* DUST */
#ifdef NEUTRAL
               hdf_vars(j) = 'vlyn' ; j = j + 1
#endif /* NEUTRAL */
#ifdef IONIZED
               hdf_vars(j) = 'vlyi' ; j = j + 1
#endif /* IONIZED */
            case ('velz')
#ifdef DUST
               hdf_vars(j) = 'vlzd' ; j = j + 1
#endif /* DUST */
#ifdef NEUTRAL
               hdf_vars(j) = 'vlzn' ; j = j + 1
#endif /* NEUTRAL */
#ifdef IONIZED
               hdf_vars(j) = 'vlzi' ; j = j + 1
#endif /* IONIZED */
            case ('ener')
#ifdef NEUTRAL
               hdf_vars(j) = 'enen' ; j = j + 1
#endif /* NEUTRAL */
#ifdef IONIZED
               hdf_vars(j) = 'enei' ; j = j + 1
#endif /* IONIZED */
            case ('magx')
               hdf_vars(j) = 'magx' ; j = j + 1
            case ('magy')
               hdf_vars(j) = 'magy' ; j = j + 1
            case ('magz')
               hdf_vars(j) = 'magz' ; j = j + 1
#ifdef COSM_RAYS
            case ('encr')
#ifndef NEW_HDF5
               do k = 1, size(iarr_all_crs,1)
                  if (k<=9) then
                     write(aux,'(A3,I1)') 'ecr',k
                     hdf_vars(j) = aux ; j = j + 1
                  else
                     write(msg, '(a,i3)')"[dataio_hdf5:init_hdf5] Cannot create name for CR energy component #",k
                     call warn(msg)
                  endif
               enddo
#endif /* !NEW_HDF5 */
#endif /* COSM_RAYS */
#ifdef GRAV
            case ('gpot')
               hdf_vars(j) = 'gpot' ; j = j + 1
#ifdef MULTIGRID
            case ('mgso')
               hdf_vars(j) = 'mgso' ; j = j + 1
#endif /* MULTIGRID */
#endif /* GRAV */
            case ('pres')
#ifdef NEUTRAL
               hdf_vars(j) = 'pren' ; j = j + 1
#endif /* NEUTRAL */
#ifdef IONIZED
               hdf_vars(j) = 'prei' ; j = j + 1
#endif /* IONIZED */
            case default
               hdf_vars(j) = trim(vars(i)) ; j = j + 1
         end select
      enddo

      if ( .not. associated(additional_attrs))      additional_attrs      => null_attrs
      if ( .not. associated(problem_write_restart)) problem_write_restart => null_attrs
      if ( .not. associated(problem_read_restart))  problem_read_restart  => null_attrs

   end subroutine init_hdf5

!>
!! \brief Procedure finalizing HDF5 module
!<
   subroutine cleanup_hdf5
      implicit none

      if (allocated(hdf_vars)) deallocate(hdf_vars)

   end subroutine cleanup_hdf5

!>
!! \brief Routine calculating quantities for plot files
!<
   subroutine common_plt_hdf5(var,ij,xn,tab,ierrh)
      use arrays,        only: u, b
      use dataio_pub,    only: varlen, planelen
      use grid,          only: cg
#ifdef GRAV
      use arrays,        only: gpot
#endif /* GRAV */
      use fluidindex,    only: flind, ibx, iby, ibz
#ifdef COSM_RAYS
      use fluidindex,    only: iarr_all_crs
#endif /* COSM_RAYS */

      implicit none
      character(len=varlen)   :: var !< quantity to be plotted
      character(len=planelen) :: ij  !< plane of plot
      integer                 :: xn  !< no. of cell at which we are slicing the local block
      integer                 :: ierrh !< error handling
      real, dimension(:,:)    :: tab !< array containing given quantity
#ifdef COSM_RAYS
      integer                 :: i
#endif /* COSM_RAYS */

      ierrh = 0
      select case (var)
         case ("dend")
            if (ij=="yz") tab(:,:) = u(flind%dst%idn, xn, cg%js:cg%je, cg%ks:cg%ke)
            if (ij=="xz") tab(:,:) = u(flind%dst%idn, cg%is:cg%ie, xn, cg%ks:cg%ke)
            if (ij=="xy") tab(:,:) = u(flind%dst%idn, cg%is:cg%ie, cg%js:cg%je, xn)
         case ("denn")
            if (ij=="yz") tab(:,:) = u(flind%neu%idn, xn, cg%js:cg%je, cg%ks:cg%ke)
            if (ij=="xz") tab(:,:) = u(flind%neu%idn, cg%is:cg%ie, xn, cg%ks:cg%ke)
            if (ij=="xy") tab(:,:) = u(flind%neu%idn, cg%is:cg%ie, cg%js:cg%je, xn)
         case ("deni")
            if (ij=="yz") tab(:,:) = u(flind%ion%idn, xn, cg%js:cg%je, cg%ks:cg%ke)
            if (ij=="xz") tab(:,:) = u(flind%ion%idn, cg%is:cg%ie, xn, cg%ks:cg%ke)
            if (ij=="xy") tab(:,:) = u(flind%ion%idn, cg%is:cg%ie, cg%js:cg%je, xn)
         case ("vlxd")
            if (ij=="yz") tab(:,:) = u(flind%dst%imx, xn, cg%js:cg%je, cg%ks:cg%ke) / &
                           u(flind%dst%idn, xn, cg%js:cg%je, cg%ks:cg%ke)
            if (ij=="xz") tab(:,:) = u(flind%dst%imx, cg%is:cg%ie, xn, cg%ks:cg%ke) / &
                           u(flind%dst%idn, cg%is:cg%ie, xn, cg%ks:cg%ke)
            if (ij=="xy") tab(:,:) = u(flind%dst%imx, cg%is:cg%ie, cg%js:cg%je, xn) / &
                           u(flind%dst%idn, cg%is:cg%ie, cg%js:cg%je, xn)
         case ("vlxn")
            if (ij=="yz") tab(:,:) = u(flind%neu%imx, xn, cg%js:cg%je, cg%ks:cg%ke) / &
                           u(flind%neu%idn, xn, cg%js:cg%je, cg%ks:cg%ke)
            if (ij=="xz") tab(:,:) = u(flind%neu%imx, cg%is:cg%ie, xn, cg%ks:cg%ke) / &
                           u(flind%neu%idn, cg%is:cg%ie, xn, cg%ks:cg%ke)
            if (ij=="xy") tab(:,:) = u(flind%neu%imx, cg%is:cg%ie, cg%js:cg%je, xn) / &
                           u(flind%neu%idn, cg%is:cg%ie, cg%js:cg%je, xn)
         case ("vlxi")
            if (ij=="yz") tab(:,:) = u(flind%ion%imx, xn, cg%js:cg%je, cg%ks:cg%ke) / &
                           u(flind%ion%idn, xn, cg%js:cg%je, cg%ks:cg%ke)
            if (ij=="xz") tab(:,:) = u(flind%ion%imx, cg%is:cg%ie, xn, cg%ks:cg%ke) / &
                           u(flind%ion%idn, cg%is:cg%ie, xn, cg%ks:cg%ke)
            if (ij=="xy") tab(:,:) = u(flind%ion%imx, cg%is:cg%ie, cg%js:cg%je, xn) / &
                           u(flind%ion%idn, cg%is:cg%ie, cg%js:cg%je, xn)
         case ("vlyd")
            if (ij=="yz") tab(:,:) = u(flind%dst%imy, xn, cg%js:cg%je, cg%ks:cg%ke) / &
                           u(flind%dst%idn, xn, cg%js:cg%je, cg%ks:cg%ke)
            if (ij=="xz") tab(:,:) = u(flind%dst%imy, cg%is:cg%ie, xn, cg%ks:cg%ke) / &
                           u(flind%dst%idn, cg%is:cg%ie, xn, cg%ks:cg%ke)
            if (ij=="xy") tab(:,:) = u(flind%dst%imy, cg%is:cg%ie, cg%js:cg%je, xn) / &
                           u(flind%dst%idn, cg%is:cg%ie, cg%js:cg%je, xn)
         case ("vlyn")
            if (ij=="yz") tab(:,:) = u(flind%neu%imy, xn, cg%js:cg%je, cg%ks:cg%ke) / &
                           u(flind%neu%idn, xn, cg%js:cg%je, cg%ks:cg%ke)
            if (ij=="xz") tab(:,:) = u(flind%neu%imy, cg%is:cg%ie, xn, cg%ks:cg%ke) / &
                           u(flind%neu%idn, cg%is:cg%ie, xn, cg%ks:cg%ke)
            if (ij=="xy") tab(:,:) = u(flind%neu%imy, cg%is:cg%ie, cg%js:cg%je, xn) / &
                           u(flind%neu%idn, cg%is:cg%ie, cg%js:cg%je, xn)
         case ("vlyi")
            if (ij=="yz") tab(:,:) = u(flind%ion%imy, xn, cg%js:cg%je, cg%ks:cg%ke) / &
                           u(flind%ion%idn, xn, cg%js:cg%je, cg%ks:cg%ke)
            if (ij=="xz") tab(:,:) = u(flind%ion%imy, cg%is:cg%ie, xn, cg%ks:cg%ke) / &
                           u(flind%ion%idn, cg%is:cg%ie, xn, cg%ks:cg%ke)
            if (ij=="xy") tab(:,:) = u(flind%ion%imy, cg%is:cg%ie, cg%js:cg%je, xn) / &
                           u(flind%ion%idn, cg%is:cg%ie, cg%js:cg%je, xn)
         case ("vlzd")
            if (ij=="yz") tab(:,:) = u(flind%dst%imz, xn, cg%js:cg%je, cg%ks:cg%ke) / &
                           u(flind%dst%idn, xn, cg%js:cg%je, cg%ks:cg%ke)
            if (ij=="xz") tab(:,:) = u(flind%dst%imz, cg%is:cg%ie, xn, cg%ks:cg%ke) / &
                           u(flind%dst%idn, cg%is:cg%ie, xn, cg%ks:cg%ke)
            if (ij=="xy") tab(:,:) = u(flind%dst%imz, cg%is:cg%ie, cg%js:cg%je, xn) / &
                           u(flind%dst%idn, cg%is:cg%ie, cg%js:cg%je, xn)
         case ("vlzn")
            if (ij=="yz") tab(:,:) = u(flind%neu%imz, xn, cg%js:cg%je, cg%ks:cg%ke) / &
                           u(flind%neu%idn, xn, cg%js:cg%je, cg%ks:cg%ke)
            if (ij=="xz") tab(:,:) = u(flind%neu%imz, cg%is:cg%ie, xn, cg%ks:cg%ke) / &
                           u(flind%neu%idn, cg%is:cg%ie, xn, cg%ks:cg%ke)
            if (ij=="xy") tab(:,:) = u(flind%neu%imz, cg%is:cg%ie, cg%js:cg%je, xn) / &
                           u(flind%neu%idn, cg%is:cg%ie, cg%js:cg%je, xn)
         case ("vlzi")
            if (ij=="yz") tab(:,:) = u(flind%ion%imz, xn, cg%js:cg%je, cg%ks:cg%ke) / &
                           u(flind%ion%idn, xn, cg%js:cg%je, cg%ks:cg%ke)
            if (ij=="xz") tab(:,:) = u(flind%ion%imz, cg%is:cg%ie, xn, cg%ks:cg%ke) / &
                           u(flind%ion%idn, cg%is:cg%ie, xn, cg%ks:cg%ke)
            if (ij=="xy") tab(:,:) = u(flind%ion%imz, cg%is:cg%ie, cg%js:cg%je, xn) / &
                           u(flind%ion%idn, cg%is:cg%ie, cg%js:cg%je, xn)
         case ("enen")
#ifdef ISO
            if (ij=="yz") tab(:,:) = 0.5 * (                     &
                          u(flind%neu%imx, xn, cg%js:cg%je, cg%ks:cg%ke)**2 &
                        + u(flind%neu%imy, xn, cg%js:cg%je, cg%ks:cg%ke)**2 &
                        + u(flind%neu%imz, xn, cg%js:cg%je, cg%ks:cg%ke)**2) / &
                             u(flind%neu%idn, xn, cg%js:cg%je, cg%ks:cg%ke)
            if (ij=="xz") tab(:,:) = 0.5 * (                     &
                          u(flind%neu%imx, cg%is:cg%ie, xn, cg%ks:cg%ke)**2 &
                         +u(flind%neu%imy, cg%is:cg%ie, xn, cg%ks:cg%ke)**2 &
                         +u(flind%neu%imz, cg%is:cg%ie, xn, cg%ks:cg%ke)**2) / &
                             u(flind%neu%idn, cg%is:cg%ie, xn, cg%ks:cg%ke)
            if (ij=="xy") tab(:,:) = 0.5 * (                     &
                          u(flind%neu%imx, cg%is:cg%ie, cg%js:cg%je, xn)**2 &
                         +u(flind%neu%imy, cg%is:cg%ie, cg%js:cg%je, xn)**2 &
                         +u(flind%neu%imz, cg%is:cg%ie, cg%js:cg%je, xn)**2) / &
                             u(flind%neu%idn, cg%is:cg%ie, cg%js:cg%je, xn)
#else /* !ISO */
            if (ij=="yz") tab(:,:) = u(flind%neu%ien, xn, cg%js:cg%je, cg%ks:cg%ke)
            if (ij=="xz") tab(:,:) = u(flind%neu%ien, cg%is:cg%ie, xn, cg%ks:cg%ke)
            if (ij=="xy") tab(:,:) = u(flind%neu%ien, cg%is:cg%ie, cg%js:cg%je, xn)
#endif /* !ISO */
         case ('prei')
            tab(:,:) = 0.0
#ifdef NEUTRAL
         case ("pren")
#ifdef ISO
            tab = 0.0
#else /* !ISO */
            if (ij=="xy") then
               tab(:,:) = real( u(flind%neu%ien, cg%is:cg%ie, cg%js:cg%je, xn) - &
                 0.5 *( u(flind%neu%imx, cg%is:cg%ie, cg%js:cg%je, xn)**2 + u(flind%neu%imy, cg%is:cg%ie, cg%js:cg%je, xn)**2 + &
                        u(flind%neu%imz, cg%is:cg%ie, cg%js:cg%je, xn)**2 ) / u(flind%neu%idn, cg%is:cg%ie, cg%js:cg%je, xn),4)*(flind%neu%gam_1)
            endif
#endif /* !ISO */
#endif /* NEUTRAL */
         case ("enei")
#ifdef ISO
            if (ij=="yz") tab(:,:) = 0.5 * (                     &
                          u(flind%ion%imx, xn, cg%js:cg%je, cg%ks:cg%ke)**2 &
                        + u(flind%ion%imy, xn, cg%js:cg%je, cg%ks:cg%ke)**2 &
                        + u(flind%ion%imz, xn, cg%js:cg%je, cg%ks:cg%ke)**2) / &
                             u(flind%ion%idn, xn, cg%js:cg%je, cg%ks:cg%ke)
            if (ij=="xz") tab(:,:) = 0.5 * (                     &
                          u(flind%ion%imx, cg%is:cg%ie, xn, cg%ks:cg%ke)**2 &
                         +u(flind%ion%imy, cg%is:cg%ie, xn, cg%ks:cg%ke)**2 &
                         +u(flind%ion%imz, cg%is:cg%ie, xn, cg%ks:cg%ke)**2) / &
                             u(flind%ion%idn, cg%is:cg%ie, xn, cg%ks:cg%ke)
            if (ij=="xy") tab(:,:) = 0.5 * (                     &
                          u(flind%ion%imx, cg%is:cg%ie, cg%js:cg%je, xn)**2 &
                         +u(flind%ion%imy, cg%is:cg%ie, cg%js:cg%je, xn)**2 &
                         +u(flind%ion%imz, cg%is:cg%ie, cg%js:cg%je, xn)**2) / &
                             u(flind%ion%idn, cg%is:cg%ie, cg%js:cg%je, xn)
#else /* !ISO */
            if (ij=="yz") tab(:,:) = u(flind%ion%ien, xn, cg%js:cg%je, cg%ks:cg%ke)
            if (ij=="xz") tab(:,:) = u(flind%ion%ien, cg%is:cg%ie, xn, cg%ks:cg%ke)
            if (ij=="xy") tab(:,:) = u(flind%ion%ien, cg%is:cg%ie, cg%js:cg%je, xn)
#endif /* !ISO */

         case ("magx")
            if (ij=="yz") tab(:,:) = b(ibx, xn, cg%js:cg%je, cg%ks:cg%ke)
            if (ij=="xz") tab(:,:) = b(ibx, cg%is:cg%ie, xn, cg%ks:cg%ke)
            if (ij=="xy") tab(:,:) = b(ibx, cg%is:cg%ie, cg%js:cg%je, xn)
         case ("magy")
            if (ij=="yz") tab(:,:) = b(iby, xn, cg%js:cg%je, cg%ks:cg%ke)
            if (ij=="xz") tab(:,:) = b(iby, cg%is:cg%ie, xn, cg%ks:cg%ke)
            if (ij=="xy") tab(:,:) = b(iby, cg%is:cg%ie, cg%js:cg%je, xn)
         case ("magz")
            if (ij=="yz") tab(:,:) = b(ibz, xn, cg%js:cg%je, cg%ks:cg%ke)
            if (ij=="xz") tab(:,:) = b(ibz, cg%is:cg%ie, xn, cg%ks:cg%ke)
            if (ij=="xy") tab(:,:) = b(ibz, cg%is:cg%ie, cg%js:cg%je, xn)
#ifdef GRAV
         case ("gpot")
            if (ij=="yz") tab(:,:) = gpot(xn, cg%js:cg%je, cg%ks:cg%ke)
            if (ij=="xz") tab(:,:) = gpot(cg%is:cg%ie, xn, cg%ks:cg%ke)
            if (ij=="xy") tab(:,:) = gpot(cg%is:cg%ie, cg%js:cg%je, xn)
#endif /* GRAV */
         case default
#ifdef COSM_RAYS
            if (var(1:3) == 'ecr') then
               i = iarr_all_crs(ichar(var(4:4))-48)
               if (ij=="yz") tab(:,:) = u(i, xn, cg%js:cg%je, cg%ks:cg%ke)
               if (ij=="xz") tab(:,:) = u(i, cg%is:cg%ie, xn, cg%ks:cg%ke)
               if (ij=="xy") tab(:,:) = u(i, cg%is:cg%ie, cg%js:cg%je, xn)
            else
#endif /* COSM_RAYS */
            ierrh = -1
#ifdef COSM_RAYS
            endif
#endif /* COSM_RAYS */
      end select

   end subroutine common_plt_hdf5

!>
!! \brief Routine calculating quantities for .hdf files
!<
   subroutine common_vars_hdf5(var,tab, ierrh)

      use arrays,        only: u, b
      use dataio_pub,    only: varlen
      use fluidindex,    only: flind, ibx, iby, ibz
      use grid,          only: cg  ! QA_WARN required for RNG
#ifdef GRAV
      use arrays,        only: gpot
#ifdef MULTIGRID
      use arrays,        only: sgp
#endif /* MULTIGRID */
#endif /* GRAV */

      implicit none
      character(len=varlen), intent(in) :: var
      real(kind=4), dimension(:,:,:)    :: tab
      integer, intent(out)              :: ierrh
#ifdef COSM_RAYS
      integer :: i
      integer, parameter    :: auxlen = varlen - 1
      character(len=auxlen) :: aux
#endif /* COSM_RAYS */

      ierrh = 0
      select case (var)
#ifdef COSM_RAYS
         case ("ecr1" : "ecr9")
            read(var,'(A3,I1)') aux,i !> \deprecated BEWARE 0 <= i <= 9, no other indices can be dumped to hdf file
            tab(:,:,:) = real(u(flind%crs%beg+i-1,RNG),4)
#endif /* COSM_RAYS */
         case ("dend")
            tab(:,:,:) = real(u(flind%dst%idn,RNG),4)
         case ("denn")
            tab(:,:,:) = real(u(flind%neu%idn,RNG),4)
         case ("deni")
            tab(:,:,:) = real(u(flind%ion%idn,RNG),4)
         case ("vlxd")
            tab(:,:,:) = real(u(flind%dst%imx,RNG) / u(flind%dst%idn,RNG),4)
         case ("vlxn")
            tab(:,:,:) = real(u(flind%neu%imx,RNG) / u(flind%neu%idn,RNG),4)
         case ("vlxi")
            tab(:,:,:) = real(u(flind%ion%imx,RNG) / u(flind%ion%idn,RNG),4)
         case ("vlyd")
            tab(:,:,:) = real(u(flind%dst%imy,RNG) / u(flind%dst%idn,RNG),4)
         case ("vlyn")
            tab(:,:,:) = real(u(flind%neu%imy,RNG) / u(flind%neu%idn,RNG),4)
         case ("vlyi")
            tab(:,:,:) = real(u(flind%ion%imy,RNG) / u(flind%ion%idn,RNG),4)
         case ("vlzd")
            tab(:,:,:) = real(u(flind%dst%imz,RNG) / u(flind%dst%idn,RNG),4)
         case ("vlzn")
            tab(:,:,:) = real(u(flind%neu%imz,RNG) / u(flind%neu%idn,RNG),4)
         case ("vlzi")
            tab(:,:,:) = real(u(flind%ion%imz,RNG) / u(flind%ion%idn,RNG),4)
#ifdef NEUTRAL
         case ("pren")
#ifdef ISO
            tab = 0.0
#else /* !ISO */
            tab(:,:,:) = real( u(flind%neu%ien,RNG) - &
              0.5 *( u(flind%neu%imx,RNG)**2 + u(flind%neu%imy,RNG)**2 + u(flind%neu%imz,RNG)**2 ) / u(flind%neu%idn,RNG),4)*(flind%neu%gam_1)
#endif /* !ISO */
#endif /* NEUTRAL */
         case ("enen")
#ifdef ISO
            tab(:,:,:) = real(0.5 *( u(flind%neu%imx,RNG)**2 + &
                                     u(flind%neu%imy,RNG)**2 + &
                                     u(flind%neu%imz,RNG)**2 ) &
                              / u(flind%neu%idn,RNG),4)
#else /* !ISO */
            tab(:,:,:) = real(u(flind%neu%ien,RNG),4)
#endif /* !ISO */
#ifdef IONIZED
         case ("enei")
#ifdef ISO
            tab(:,:,:) = real(0.5 *( u(flind%ion%imx,RNG)**2 + &
                                     u(flind%ion%imy,RNG)**2 + &
                                     u(flind%ion%imz,RNG)**2 ) &
                              / u(flind%ion%idn,RNG),4)
#else /* !ISO */
            tab(:,:,:) = real(u(flind%ion%ien,RNG),4)
#endif /* !ISO */
         case ("prei")
#ifdef ISO
            tab = 0.0
#else /* !ISO */
            tab(:,:,:) = real( u(flind%ion%ien,RNG) - &
              0.5 *( u(flind%ion%imx,RNG)**2 + u(flind%ion%imy,RNG)**2 + u(flind%ion%imz,RNG)**2 ) / u(flind%ion%idn,RNG),4)*(flind%ion%gam_1)
            tab(:,:,:) = tab(:,:,:) - real( 0.5*(flind%ion%gam_1)*(b(ibx,RNG)**2 + &
               b(iby,RNG)**2 + b(ibz,RNG)**2),4)
#endif /* !ISO */
#endif /* IONIZED */
         case ("magx")
            tab(:,:,:) = real(b(ibx,RNG),4)
         case ("magy")
            tab(:,:,:) = real(b(iby,RNG),4)
         case ("magz")
            tab(:,:,:) = real(b(ibz,RNG),4)
#ifdef GRAV
         case ("gpot")
            tab(:,:,:) = real(gpot(RNG),4)
#ifdef MULTIGRID
         case ("mgso")
            tab(:,:,:) = real(sgp(RNG),4)
#endif /* MULTIGRID */
#endif /* GRAV */

         case default
            ierrh = -1
      end select

   end subroutine common_vars_hdf5

   subroutine write_plot

      use dataio_pub,    only: cwdlen, log_file
      use hdf5,          only: HID_T, H5open_f, H5Fcreate_f, H5Gcreate_f, H5F_ACC_TRUNC_F, H5Gclose_f, H5close_f, h5fclose_f
      use mpisetup,      only: t, comm3d, ierr, master

      implicit none

      integer, save     :: nimg = 0, error = 0
      real, save        :: last_plt_time = 0.0
      character(len=cwdlen) :: fname
      integer           :: i, fe
      logical, save     :: first_entry = .true.
      integer(HID_T)    :: file_id                 !> File identifier
      integer(HID_T)    :: gr_id, gr2_id           !> Groups identifier

      if ( ((t-last_plt_time > dt_plt) .or. first_entry ) .and. dt_plt > 0.0 ) then
         fe = len_trim(log_file)
         write(fname,'(2a)') trim(log_file(1:fe-3)),"plt"
         call H5open_f(error)

         if (master .and. first_entry) then
            call H5Fcreate_f(fname, H5F_ACC_TRUNC_F, file_id, error)
            call H5Gcreate_f(file_id,"xy",gr_id,error)
            do i=1,nhdf_vars
               call H5Gcreate_f(gr_id,hdf_vars(i),gr2_id,error)
               call H5Gclose_f(gr2_id,error)
            enddo
            call H5Gclose_f(gr_id,error)
            call H5Gcreate_f(file_id,"yz",gr_id,error)
            do i = 1, nhdf_vars
               call H5Gcreate_f(gr_id,hdf_vars(i),gr2_id,error)
               call H5Gclose_f(gr2_id,error)
            enddo
            call H5Gclose_f(gr_id,error)
            call H5Gcreate_f(file_id,"xz",gr_id,error)
            do i = 1, nhdf_vars
               call H5Gcreate_f(gr_id,hdf_vars(i),gr2_id,error)
               call H5Gclose_f(gr2_id,error)
            enddo
            call H5Gclose_f(gr_id,error)
            call H5Fclose_f(file_id,error)
            first_entry = .false.
         endif
         call MPI_Barrier(comm3d,ierr)
         do i = 1,nhdf_vars
            if (ix > 0) call write_plot_hdf5(hdf_vars(i),"yz",nimg)
            if (iy > 0) call write_plot_hdf5(hdf_vars(i),"xz",nimg)
            if (iz > 0) call write_plot_hdf5(hdf_vars(i),"xy",nimg)
         enddo

         nimg = nimg+1
         first_entry = .false.
         call H5close_f(error)

         last_plt_time = t

      endif

   end subroutine write_plot

   subroutine write_plot_hdf5(var,plane,nimg)

      use arrays,        only: u
      use dataio_pub,    only: vizit, fmin, fmax, cwdlen, log_file, msg, varlen, die, warn, user_plt_hdf5, planelen
      use grid,          only: cg
      use hdf5,          only: HID_T, HSIZE_T, SIZE_T, H5F_ACC_RDWR_F, h5fopen_f, h5gopen_f, h5gclose_f, h5fclose_f
      use h5lt,          only: h5ltmake_dataset_double_f, h5ltset_attribute_double_f
      use mpisetup,      only: comm3d, ierr, psize, t, xdim, ydim, zdim, has_dir, dom
      use mpi,           only: MPI_CHARACTER, MPI_DOUBLE_PRECISION
#ifdef PGPLOT
      use viz,           only: draw_me
#endif /* PGPLOT */

      implicit none

      character(len=planelen), intent(in) :: plane
      character(len=varlen), intent(in)   :: var                           !> not yet implemented
      integer, intent(in)                 :: nimg
      logical, dimension(3)               :: remain
      logical                             :: ok_plt_var
      real, dimension(:), allocatable     :: buff
      real, dimension(:,:), allocatable   :: send, img
      real, dimension(:,:,:), allocatable :: temp
      real                                :: imax, imin, di
      integer                             :: ierrh, i, j
      integer                             :: comm2d, lp, ls, xn, error
      integer, parameter                  :: suffixed_planelen = planelen + 1  !< plane name + suffix e.g "xy_"
      character(len=suffixed_planelen)    :: pij
      character(len=cwdlen)               :: fname
      integer, parameter                  :: vdn_len = 12
      character(len=vdn_len)              :: vdname

      integer(HID_T)                      :: file_id                       !> File identifier
      integer(HID_T)                      :: gr_id,gr2_id                  !> Group identifier
      integer(HSIZE_T), dimension(2)      :: dims
      integer(SIZE_T)                     :: bufsize
      integer                             :: rank

      integer                             :: nib, nid, njb, njd, nkb
      integer                             :: pisize, pjsize, fe

      real, dimension(1)                  :: timebuffer

      rank = 2
      fe = len_trim(log_file)
      write(fname,'(2a)') trim(log_file(1:fe-3)),"plt"
      call MPI_Bcast(fname, cwdlen, MPI_CHARACTER, 0, comm3d, ierr)

      nib = 0; nid = 0; njb = 0; njd = 0; nkb = 0; pisize = 0; pjsize = 0
      select case (plane)
         case ("yz")
            if (has_dir(xdim)) then
               xn     = ix + cg%nb - cg%off(xdim)
            else
               xn     = 1
            endif
            remain = (/.false.,.true.,.true./)
            pij    = "yz_"
            nib    = cg%nyb
            nid    = dom%n_d(ydim)
            njb    = cg%nzb
            njd    = dom%n_d(zdim)
            nkb    = cg%nxb
            pisize = psize(ydim)
            pjsize = psize(zdim)
         case ("xz")
            if (has_dir(ydim)) then
               xn     = iy + cg%nb - cg%off(ydim)
            else
               xn     = 1
            endif
            remain = (/.true.,.false.,.true./)
            pij    = "xz_"
            nib    = cg%nxb
            nid    = dom%n_d(xdim)
            njb    = cg%nzb
            njd    = dom%n_d(zdim)
            nkb    = cg%nyb
            pisize = psize(xdim)
            pjsize = psize(zdim)
         case ("xy")
            if (has_dir(zdim)) then
               xn = iz + cg%nb - cg%off(zdim)
            else
               xn = 1
            endif
            remain = (/.true.,.true.,.false./)
            pij    = "xy_"
            nib    = cg%nxb
            nid    = dom%n_d(xdim)
            njb    = cg%nyb
            njd    = dom%n_d(ydim)
            nkb    = cg%nzb
            pisize = psize(xdim)
            pjsize = psize(ydim)
         case default
            call die("[dataio_hdf5:write_plot_hdf5] nonrecognized plane")
      end select

      dims(1) = nid
      dims(2) = njd
      call MPI_Barrier(comm3d,ierr)
      call MPI_Cart_sub(comm3d,remain,comm2d,ierr)
      call MPI_Comm_size(comm2d, ls, ierr)
      call MPI_Comm_rank(comm2d, lp, ierr)
      if ((xn > cg%nb .and. xn <= nkb+cg%nb).or.xn == 1) then
         allocate(temp(nib,njb,pisize*pjsize),img(nid,njd))
         allocate(buff(nid*njd))
         allocate(send(nib,njb))

         ok_plt_var = .false.
         call common_plt_hdf5(var,plane,xn,send,ierrh)
         if (associated(user_plt_hdf5) .and. ierrh /= 0) call user_plt_hdf5(var,plane,xn,send,ierrh)
         if (ierrh==0) ok_plt_var = .true.
         if (.not.ok_plt_var) then
            write(msg,'(2a)') var, " is not defined in common_plt_hdf5, neither in user_plt_hdf5 !!!"
            call die(msg)
         endif

         temp = -1.0
         call MPI_Gather(send, nib*njb, MPI_DOUBLE_PRECISION, temp, nib*njb, MPI_DOUBLE_PRECISION, 0, comm2d,ierr)

         if (lp == 0) then
            imax = maxval(temp); imin = minval(temp)
            di = imax-imin
            do i = 0, pisize-1
               do j = 0, pjsize-1
                  img(i*nib+1:(i+1)*nib,j*njb+1:(j+1)*njb) = temp(:,:,(j+1)+i*pjsize)
               enddo
            enddo
            if (vizit) then
#ifdef PGPLOT
               call draw_me(real(img,4), real(fmin,4), real(fmax,4))
#else /* !PGPLOT */
               call warn("[dataio_hdf5:write_plot_hdf5] vizit used without PGPLOT")
#endif /* !PGPLOT */
            else
               call H5Fopen_f(fname, H5F_ACC_RDWR_F, file_id, error)
               call H5Gopen_f(file_id,plane,gr_id,error)
               call H5Gopen_f(gr_id,var,gr2_id,error)
               write(vdname,'(a3,a4,a1,i4.4)') pij,var,"_",nimg
               call h5ltmake_dataset_double_f(gr2_id, vdname, rank, dims, img, error)
               bufsize = 1
               timebuffer = (/t/)
               call h5ltset_attribute_double_f(gr2_id,vdname,"time",timebuffer,bufsize,error)
               call H5Gclose_f(gr2_id,error)
               call H5Gclose_f(gr_id,error)
               call H5Fclose_f(file_id,error)
            endif
         endif
         if (allocated(send)) deallocate(send)
         if (allocated(temp)) deallocate(temp)
         if (allocated(img))  deallocate(img)
      endif
      call MPI_Barrier(comm3d,ierr)

   end subroutine write_plot_hdf5

!>
!! \brief Routine to set parameters and dimensions of arrays in restart file
!! \param area_type case name; possibilities: all - whole domain with mpi boundaries, outbnd - physical domain with outer boundaries, nobnd - only physical domain without any boundaries
!! \param area grid dimensions in the file
!! \param chnk dimensions of data array dumped by this process
!! \param lleft left limits of data from array to be dumped
!! \param lright right limits of data from array to be dumped
!! \loffs offset in area for this process
!<
   subroutine set_dims_to_write(area_type, area, chnk, lleft, lright, loffs)

      use grid,     only: cg
      use mpisetup, only: dom, has_dir, psize, pcoords

      implicit none

      character(len=*),      intent(in)  :: area_type
      integer, dimension(3), intent(out) :: area, lleft, lright, loffs, chnk

      select case (area_type)
         case ('allbnd')                           ! whole domain with mpi boundaries
            chnk(:)   = [cg%nx, cg%ny, cg%nz]
            area(:)   = chnk*psize
            lleft(:)  = 1
            lright(:) = chnk
            loffs(:)  = cg%off(:) + 2 * cg%nb * pcoords(:) !\todo invent something better
         case ('outbnd')                           ! physical domain with outer boundaries
            area(:)   = [dom%nxt, dom%nyt, dom%nzt]
            lleft(:)  = [cg%is,   cg%js,   cg%ks  ]
            lright(:) = [cg%ie,   cg%je,   cg%ke  ]
            chnk(:)   = cg%n_b(:)
            where (cg%off(:) == 0 .and. has_dir(:))
               lleft(:)  = lleft(:)  - cg%nb
               chnk(:)   = chnk(:)   + cg%nb
            endwhere
            where (cg%off(:) + cg%n_b(:) == dom%n_d(:) .and. has_dir(:))
               lright(:) = lright(:) + cg%nb
               chnk(:)   = chnk(:)   + cg%nb
            endwhere
            loffs(:)  = cg%off(:)
            where (loffs(:)>0) loffs(:) = loffs(:) + cg%nb ! the block adjacent to the left boundary are cg%nb cells wider than cg%n[xyz]b
         case ('no_bnd')                           ! only physical domain without any boundaries
            area(:)   = dom%n_d(:)
            lleft(:)  = [cg%is,   cg%js,   cg%ks  ]
            lright(:) = [cg%ie,   cg%je,   cg%ke  ]
            chnk(:)   = cg%n_b(:)
            loffs(:)  = cg%off(:)
      endselect

   end subroutine set_dims_to_write

!>
!! \brief This routine writes restart dump and updates restart counter
!<

   subroutine write_restart_hdf5(debug_res)

      use arrays,        only: u, b
      use dataio_pub,    only: chdf, nres, set_container_chdf, cwdlen, problem_name, run_id
      use fluidindex,    only: flind
      use grid,          only: cg
      use hdf5,          only: HID_T, HSIZE_T, HSSIZE_T, H5P_FILE_ACCESS_F, H5F_ACC_TRUNC_F, H5P_DATASET_CREATE_F, H5S_SELECT_SET_F, &
           &                   H5P_DATASET_XFER_F, H5FD_MPIO_INDEPENDENT_F, H5T_NATIVE_DOUBLE, &
           &                   h5open_f, h5close_f, h5fcreate_f, h5fclose_f, &
           &                   h5dcreate_f, h5dwrite_f, h5dclose_f, h5dget_space_f, &
           &                   h5pcreate_f, h5pclose_f, h5pset_fapl_mpio_f, h5pset_chunk_f, h5pset_dxpl_mpio_f, &
           &                   h5screate_simple_f, h5sclose_f, h5sselect_hyperslab_f
      use list_hdf5,     only: problem_write_restart
      use mpisetup,      only: comm3d, comm, info, ierr, master, nstep, dom, xdim, ydim, zdim
      use mpi,           only: MPI_CHARACTER
      use types,         only: hdf
#ifdef ISO_LOCAL
      use arrays,        only: cs_iso2_arr
#endif /* ISO_LOCAL */
#ifdef GRAV
      use arrays,        only: gp
#endif /* GRAV */

      implicit none
      logical, optional, intent(in) :: debug_res
      integer            :: nu
      character(len=cwdlen) :: filename  !> HDF File name

      integer(HID_T) :: file_id       !> File identifier
      integer(HID_T) :: dset_id       !> Dataset identifier
      integer(HID_T) :: plist_id      !> Property list identifier
      integer(HID_T) :: filespace     !> Dataspace identifier in file
      integer(HID_T) :: memspace      !> Dataspace identifier in memory

      integer(HSIZE_T),  DIMENSION(:), allocatable :: count
      integer(HSSIZE_T), DIMENSION(:), allocatable :: offset
      integer(HSIZE_T),  DIMENSION(:), allocatable :: stride
      integer(HSIZE_T),  DIMENSION(:), allocatable :: block
      integer(HSIZE_T),  DIMENSION(:), allocatable :: dimsf, dimsfi, chunk_dims

      integer :: error, rank
      real, pointer, dimension(:,:,:,:) :: p

      integer, parameter    :: atlen = 6
      integer, parameter    :: extlen = 4
      character(len=atlen)  :: area_type
      character(len=extlen) :: file_extension
      integer, dimension(3) :: area, lleft, lright, loffs, chnk

      if (present(debug_res)) then
         file_extension = '.dmp'
      else
         file_extension = '.res'
      endif

      if (master) write(filename, '(a,a1,a3,a1,i4.4,a4)') trim(problem_name), '_', run_id, '_', nres, file_extension
      call MPI_Bcast(filename, cwdlen, MPI_CHARACTER, 0, comm, ierr)
      call set_container_chdf(nstep)

      nu = flind%all

      call h5open_f(error)
      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
      call h5pset_fapl_mpio_f(plist_id, comm3d, info, error)

      call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)
      call h5pclose_f(plist_id, error)

      if (associated(problem_write_restart)) call problem_write_restart(file_id)

#ifdef ISO_LOCAL
      call write_3darr_to_restart(cs_iso2_arr(:,:,:),file_id,"cs_iso2", cg%nx, cg%ny, cg%nz)
#endif /* ISO_LOCAL */
#ifdef GRAV
      call write_3darr_to_restart(gp(:,:,:),file_id, "gp", cg%nx, cg%ny, cg%nz)
#endif /* GRAV */

      rank = 4
      allocate(dimsf(rank),dimsfi(rank),chunk_dims(rank))
      allocate(count(rank),offset(rank),stride(rank),block(rank))
      !> \deprecated The code for writing fluid data, mag field and write_3darr_to_restart is almost replicated.
      !> \todo Try to put it in a separate subroutine.
      !----------------------------------------------------------------------------------
      !  WRITE FLUID VARIABLES
      !

      area_type = 'no_bnd'
      if (present(debug_res)) area_type = 'allbnd'
      call set_dims_to_write(area_type, area, chnk, lleft, lright, loffs)
      dimsf = [nu, area(:)]      ! Dataset dimensions
      dimsfi = dimsf
      chunk_dims = [nu, chnk(:)] ! Chunks dimensions

      ! Create the data space for the  dataset.
      call h5screate_simple_f(rank, dimsf, filespace, error)
      call h5screate_simple_f(rank, chunk_dims, memspace, error)

      ! Create chunked dataset.
      call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
      if (area_type /= 'outbnd') call h5pset_chunk_f(plist_id, rank, chunk_dims, error)
      call h5dcreate_f(file_id, dname(1), H5T_NATIVE_DOUBLE, filespace, dset_id, error, plist_id)
      call h5sclose_f(filespace, error)
      call h5pclose_f(plist_id, error)

      ! Each process defines dataset in memory and writes it to the hyperslab
      ! in the file.
      stride(:) = 1
      count(:)  = 1
      block(:)  = chunk_dims(:)

      offset(:)   = [0, loffs(:)]

      ! Select hyperslab in the file.
      call h5dget_space_f(dset_id, filespace, error)
      call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error, stride, block)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
      p => u(:,lleft(xdim):lright(xdim),lleft(ydim):lright(ydim),lleft(zdim):lright(zdim))
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, p, dimsfi, error, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

      call h5sclose_f(filespace, error)
      call h5sclose_f(memspace, error)
      call h5dclose_f(dset_id, error)
      call h5pclose_f(plist_id, error)
      !----------------------------------------------------------------------------------

      !----------------------------------------------------------------------------------
      !  WRITE MAG VARIABLES
      !
      area_type = 'outbnd' ! unlike fluids, we need magnetic field boundaries values. Then chunks might be non-uniform
      if (present(debug_res)) area_type = 'allbnd'
      call set_dims_to_write(area_type, area, chnk, lleft, lright, loffs)
      dimsf = [3, area(:)] ! Dataset dimensions
      dimsfi = dimsf

      chunk_dims = [3, chnk(:)]                 ! Chunks dimensions

      ! Create the data space for the  dataset.
      call h5screate_simple_f(rank, dimsf, filespace, error)
      call h5screate_simple_f(rank, chunk_dims, memspace, error)

      ! Create chunked dataset.
      call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
      if (area_type /= 'outbnd') call h5pset_chunk_f(plist_id, rank, chunk_dims, error)
      call h5dcreate_f(file_id, dname(2), H5T_NATIVE_DOUBLE, filespace, dset_id, error, plist_id)
      call h5sclose_f(filespace, error)
      call h5pclose_f(plist_id, error)

      ! Each process defines dataset in memory and writes it to the hyperslab
      ! in the file.
      stride(:) = 1
      count(:)  = 1
      block(:)  = chunk_dims(:)

      offset(:) = [0, loffs(:)]

      ! Select hyperslab in the file.
      call h5dget_space_f(dset_id, filespace, error)
      call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error, stride, block)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
      p => b(:,lleft(xdim):lright(xdim),lleft(ydim):lright(ydim),lleft(zdim):lright(zdim))
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, p, dimsfi, error, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

      call h5sclose_f(filespace, error)
      call h5sclose_f(memspace, error)
      call h5dclose_f(dset_id, error)
      call h5pclose_f(plist_id, error)
      !----------------------------------------------------------------------------------
      if (allocated(dimsf))      deallocate(dimsf)
      if (allocated(dimsfi))     deallocate(dimsfi)
      if (allocated(chunk_dims)) deallocate(chunk_dims)
      if (allocated(count))      deallocate(count)
      if (allocated(offset))     deallocate(offset)
      if (allocated(stride))     deallocate(stride)
      if (allocated(block))      deallocate(block)

      rank = 1
      allocate(dimsf(rank),dimsfi(rank),chunk_dims(rank))
      allocate(count(rank),offset(rank),stride(rank),block(rank))

      !> \deprecated The code for writing axes is almost replicated.
      !> \todo Try to put them in a separate subroutine.
      !> \deprecated Some values are written multiple times, eg. each cg%x(:) is written by psize(ydim)*psize(zdim) CPUs.
      !----------------------------------------------------------------------------------
      !  WRITE X Axis
      !
      dimsf  = [dom%n_d(xdim)] ! Dataset dimensions
      dimsfi = dimsf
      chunk_dims = [cg%nxb]    ! Chunks dimensions

      ! Create the data space for the  dataset.
      call h5screate_simple_f(rank, dimsf, filespace, error)
      call h5screate_simple_f(rank, chunk_dims, memspace, error)

      ! Create chunked dataset.
      call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
      call h5pset_chunk_f(plist_id, rank, chunk_dims, error)
      call h5dcreate_f(file_id, "X axis", H5T_NATIVE_DOUBLE, filespace, dset_id, error, plist_id)
      call h5sclose_f(filespace, error)
      call h5pclose_f(plist_id, error)

      ! Each process defines dataset in memory and writes it to the hyperslab
      ! in the file.
      stride(:) = 1
      count(:)  = 1
      block(:)  = chunk_dims(:)

      offset(1) = cg%off(xdim)

      ! Select hyperslab in the file.
      call h5dget_space_f(dset_id, filespace, error)
      call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error, stride, block)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, cg%x(cg%is:cg%ie), dimsfi, error, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

      call h5sclose_f(filespace, error)
      call h5sclose_f(memspace, error)
      call h5dclose_f(dset_id, error)
      call h5pclose_f(plist_id, error)
      !----------------------------------------------------------------------------------

      !----------------------------------------------------------------------------------
      !  WRITE Y Axis
      !
      dimsf  = [dom%n_d(ydim)] ! Dataset dimensions
      dimsfi = dimsf
      chunk_dims = [cg%nyb]    ! Chunks dimensions

      ! Create the data space for the  dataset.
      call h5screate_simple_f(rank, dimsf, filespace, error)
      call h5screate_simple_f(rank, chunk_dims, memspace, error)

      ! Create chunked dataset.
      call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
      call h5pset_chunk_f(plist_id, rank, chunk_dims, error)
      call h5dcreate_f(file_id, "Y axis", H5T_NATIVE_DOUBLE, filespace, dset_id, error, plist_id)
      call h5sclose_f(filespace, error)
      call h5pclose_f(plist_id, error)

      ! Each process defines dataset in memory and writes it to the hyperslab
      ! in the file.
      stride(:) = 1
      count(:)  = 1
      block(:)  = chunk_dims(:)

      offset(1) = cg%off(ydim)

      ! Select hyperslab in the file.
      call h5dget_space_f(dset_id, filespace, error)
      call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error, stride, block)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, cg%y(cg%js:cg%je), dimsfi, error, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

      call h5sclose_f(filespace, error)
      call h5sclose_f(memspace, error)
      call h5dclose_f(dset_id, error)
      call h5pclose_f(plist_id, error)
      !----------------------------------------------------------------------------------

      !----------------------------------------------------------------------------------
      !  WRITE Z Axis
      !
      dimsf  = [dom%n_d(zdim)] ! Dataset dimensions
      dimsfi = dimsf
      chunk_dims = [cg%nzb]    ! Chunks dimensions

      ! Create the data space for the  dataset.
      call h5screate_simple_f(rank, dimsf, filespace, error)
      call h5screate_simple_f(rank, chunk_dims, memspace, error)

      ! Create chunked dataset.
      call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
      call h5pset_chunk_f(plist_id, rank, chunk_dims, error)
      call h5dcreate_f(file_id, "Z axis", H5T_NATIVE_DOUBLE, filespace, dset_id, error, plist_id)
      call h5sclose_f(filespace, error)
      call h5pclose_f(plist_id, error)

      ! Each process defines dataset in memory and writes it to the hyperslab
      ! in the file.
      stride(:) = 1
      count(:)  = 1
      block(:)  = chunk_dims(:)

      offset(1) = cg%off(zdim)

      ! Select hyperslab in the file.
      call h5dget_space_f(dset_id, filespace, error)
      call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error, stride, block)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, cg%z(cg%ks:cg%ke), dimsfi, error, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

      call h5sclose_f(filespace, error)
      call h5sclose_f(memspace, error)
      call h5dclose_f(dset_id, error)
      !----------------------------------------------------------------------------------
      call h5pclose_f(plist_id, error)
      call h5fclose_f(file_id, error)
      if (allocated(dimsf))      deallocate(dimsf)
      if (allocated(dimsfi))     deallocate(dimsfi)
      if (allocated(chunk_dims)) deallocate(chunk_dims)
      if (allocated(count))      deallocate(count)
      if (allocated(offset))     deallocate(offset)
      if (allocated(stride))     deallocate(stride)
      if (allocated(block))      deallocate(block)

      call set_common_attributes(filename, chdf, "restart")

      call h5close_f(error)

      nres = nres + 1
      if (associated(p)) nullify(p)

   end subroutine write_restart_hdf5

   subroutine write_3darr_to_restart(tab,file_id,dname,nx,ny,nz)

      use grid,        only: cg
      use hdf5,        only: HID_T, HSIZE_T, HSSIZE_T, H5P_DATASET_CREATE_F, H5S_SELECT_SET_F, &
           &                 H5P_DATASET_XFER_F, H5FD_MPIO_INDEPENDENT_F, H5T_NATIVE_DOUBLE, &
           &                 h5dcreate_f, h5dwrite_f, h5dclose_f, h5dget_space_f, &
           &                 h5pcreate_f, h5pset_chunk_f, h5pset_dxpl_mpio_f, &
           &                 h5screate_simple_f, h5sclose_f, h5sselect_hyperslab_f
      use mpisetup,    only: psize

      implicit none

      integer, intent(in)                      :: nx,ny,nz
      real, dimension(nx,ny,nz), intent(in)    :: tab
      integer(HID_T), intent(in)               :: file_id       !> File identifier
      character(len=*), intent(in)             :: dname

      integer(HID_T) :: dset_id       !> Dataset identifier
      integer(HID_T) :: plist_id      !> Property list identifier
      integer(HID_T) :: filespace     !> Dataspace identifier in file
      integer(HID_T) :: memspace      !> Dataspace identifier in memory

      integer(HSIZE_T),  DIMENSION(:), allocatable :: count
      integer(HSSIZE_T), DIMENSION(:), allocatable :: offset
      integer(HSIZE_T),  DIMENSION(:), allocatable :: stride
      integer(HSIZE_T),  DIMENSION(:), allocatable :: block
      integer(HSIZE_T),  DIMENSION(:), allocatable :: dimsf, dimsfi, chunk_dims

      integer :: error, rank

      rank = 3
      allocate(dimsf(rank),dimsfi(rank),chunk_dims(rank))
      allocate(count(rank),offset(rank),stride(rank),block(rank))
      !----------------------------------------------------------------------------------
      !  WRITE TAB
      !
      dimsf = [nx, ny, nz] * psize(:)           ! Dataset dimensions
      dimsfi = dimsf
      chunk_dims = [nx, ny, nz]                 ! Chunks dimensions

      ! Create the data space for the  dataset.
      call h5screate_simple_f(rank, dimsf, filespace, error)
      call h5screate_simple_f(rank, chunk_dims, memspace, error)

      ! Create chunked dataset.
      call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
      call h5pset_chunk_f(plist_id, rank, chunk_dims, error)
      call h5dcreate_f(file_id, trim(dname), H5T_NATIVE_DOUBLE, filespace, dset_id, error, plist_id)
      call h5sclose_f(filespace, error)

      ! Each process defines dataset in memory and writes it to the hyperslab
      ! in the file.
      stride(:) = 1
      count(:)  = 1
      block(:)  = chunk_dims(:)

      offset(:) = cg%off(:)

      ! Select hyperslab in the file.
      call h5dget_space_f(dset_id, filespace, error)
      call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error, stride, block)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, tab, dimsfi, error, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

      call h5sclose_f(filespace, error)
      call h5sclose_f(memspace, error)
      call h5dclose_f(dset_id, error)

      if (allocated(dimsf))      deallocate(dimsf)
      if (allocated(dimsfi))     deallocate(dimsfi)
      if (allocated(chunk_dims)) deallocate(chunk_dims)
      if (allocated(count))      deallocate(count)
      if (allocated(offset))     deallocate(offset)
      if (allocated(stride))     deallocate(stride)
      if (allocated(block))      deallocate(block)
      !----------------------------------------------------------------------------------

   end subroutine write_3darr_to_restart

   subroutine read_3darr_from_restart(file_id,dname,p3d,nx,ny,nz)

      use grid,         only: cg
      use hdf5,         only: HID_T, HSIZE_T, HSSIZE_T, SIZE_T, H5T_NATIVE_DOUBLE, &
           &                  H5S_SELECT_SET_F, H5FD_MPIO_INDEPENDENT_F, H5P_DATASET_XFER_F, &
           &                  h5pcreate_f, h5pclose_f, h5screate_simple_f, h5dopen_f, &
           &                  h5dget_space_f, h5sget_simple_extent_ndims_f, h5dget_create_plist_f, &
           &                  h5sselect_hyperslab_f, h5dread_f, h5sclose_f, h5pset_dxpl_mpio_f, h5dclose_f
      use mpisetup,     only: psize

      implicit none

      integer, intent(in)                  :: nx,ny,nz
      real, dimension(:,:,:), pointer      :: p3d
      integer(HID_T), intent(in)           :: file_id       ! File identifier
      character(len=*), intent(in)         :: dname

      integer(HID_T)        :: dset_id       ! Dataset identifier
      integer(HID_T)        :: plist_id      ! Property list identifier
      integer(HID_T)        :: filespace     ! Dataspace identifier in file
      integer(HID_T)        :: memspace      ! Dataspace identifier in memory

      integer(HSIZE_T),  DIMENSION(:), allocatable :: count
      integer(HSSIZE_T), DIMENSION(:), allocatable :: offset
      integer(HSIZE_T),  DIMENSION(:), allocatable :: stride
      integer(HSIZE_T),  DIMENSION(:), allocatable :: block
      integer(HSIZE_T),  DIMENSION(:), allocatable :: dimsf, dimsfi, chunk_dims

      integer               :: error, rank

      !----------------------------------------------------------------------------------
      !  READ TAB
      !
      rank = 3
      allocate(dimsf(rank), dimsfi(rank), chunk_dims(rank))
      allocate(block(rank), offset(rank), count(rank), stride(rank))
      dimsf = [nx, ny, nz] * psize(:)           ! Dataset dimensions
      dimsfi = dimsf
      chunk_dims = [nx, ny, nz]

      ! Create chunked dataset.
      call h5dopen_f(file_id, dname, dset_id, error)

      call h5dget_space_f(dset_id, filespace, error)
      call H5sget_simple_extent_ndims_f (filespace,rank,error)
      call H5dget_create_plist_f (dset_id,plist_id,error)

      ! Each process defines dataset in memory and writes it to the hyperslab
      ! in the file.
      stride(:) = 1
      count(:)  = 1
      block(:)  = chunk_dims(:)

      offset(:) = cg%off(:)

      ! Select hyperslab in the file.
      call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error, stride, block)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
      call h5screate_simple_f(rank, chunk_dims, memspace, error)
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, p3d, dimsfi, error, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

      call h5sclose_f(filespace, error)
      call h5sclose_f(memspace, error)
      call h5pclose_f(plist_id, error)
      call h5dclose_f(dset_id, error)
      !----------------------------------------------------------------------------------
      if (allocated(dimsf))      deallocate(dimsf)
      if (allocated(dimsfi))     deallocate(dimsfi)
      if (allocated(chunk_dims)) deallocate(chunk_dims)
      if (allocated(count))      deallocate(count)
      if (allocated(offset))     deallocate(offset)
      if (allocated(stride))     deallocate(stride)
      if (allocated(block))      deallocate(block)

   end subroutine read_3darr_from_restart

   subroutine read_restart_hdf5(chdf)

      use arrays,        only: u, b
      use dataio_pub,    only: cwdlen, msg, printio, warn, die, require_init_prob, problem_name, run_id, piernik_hdf5_version
      use fluidindex,    only: flind
      use func,          only: fix_string
      use grid,          only: cg
      use hdf5,          only: HID_T, HSIZE_T, HSSIZE_T, SIZE_T, H5P_FILE_ACCESS_F, H5T_NATIVE_DOUBLE, &
           &                   H5S_SELECT_SET_F, H5F_ACC_RDONLY_F, H5FD_MPIO_INDEPENDENT_F, H5P_DATASET_XFER_F, &
           &                   h5open_f, h5pcreate_f, h5pset_fapl_mpio_f, h5fopen_f, h5pclose_f, h5dopen_f, &
           &                   h5dget_space_f, h5sget_simple_extent_ndims_f, h5dget_create_plist_f, &
           &                   h5sselect_hyperslab_f, h5dread_f, h5sclose_f, h5pset_dxpl_mpio_f, h5dclose_f, &
           &                   h5screate_simple_f, h5fclose_f, h5close_f
      use h5lt,          only: h5ltget_attribute_double_f, h5ltget_attribute_int_f, h5ltget_attribute_string_f
      use list_hdf5,     only: problem_read_restart
      use mpisetup,      only: comm, ierr, magic_mass, master, t, info, comm3d, dt, cbuff_len, dom, has_dir, xdim, ydim, zdim
      use mpi,           only: MPI_CHARACTER, MPI_INTEGER, MPI_DOUBLE_PRECISION
      use types,         only: hdf, domlen, idlen
#ifdef ISO_LOCAL
      use arrays,        only: cs_iso2_arr
#endif /* ISO_LOCAL */
#ifdef GRAV
      use arrays,        only: gp
#endif /* GRAV */

      implicit none

      type(hdf)             :: chdf
      integer               :: nu
      character(len=cwdlen) :: filename  ! File name

      integer(HID_T)        :: file_id       ! File identifier
      integer(HID_T)        :: dset_id       ! Dataset identifier
      integer(HID_T)        :: plist_id      ! Property list identifier
      integer(HID_T)        :: filespace     ! Dataspace identifier in file
      integer(HID_T)        :: memspace      ! Dataspace identifier in memory
      integer(SIZE_T)       :: bufsize

      integer(HSIZE_T),  DIMENSION(:), allocatable :: count
      integer(HSSIZE_T), DIMENSION(:), allocatable :: offset
      integer(HSIZE_T),  DIMENSION(:), allocatable :: stride
      integer(HSIZE_T),  DIMENSION(:), allocatable :: block
      integer(HSIZE_T),  DIMENSION(:), allocatable :: dimsf, dimsfi, chunk_dims

      integer               :: error, rank
      logical               :: file_exist

      real, dimension(1)    :: rbuf
      integer, dimension(1) :: ibuf

      real, dimension(:,:,:,:), pointer :: p4d
#if defined(ISO_LOCAL) || defined(GRAV)
      real, dimension(:,:,:), pointer :: p3d
#endif /* ISO_LOCAL || GRAV */

      real                  :: restart_hdf5_version

      nu = flind%all

      if (master) then
         write(filename,'(a,a1,a3,a1,i4.4,a4)') trim(problem_name),'_', run_id,'_',chdf%nres,'.res'
         write(msg, '(2a)') 'Reading restart file: ',trim(filename)
         call printio(msg)
      endif
      call MPI_Bcast(filename, cwdlen, MPI_CHARACTER, 0, comm, ierr)

      inquire(file = filename, exist = file_exist)
      if (.not. file_exist) then
         write(msg,'(3a)') '[dataio_hdf5:read_restart_hdf5]: Restart file: ',trim(filename),' does not exist'
         call die(msg)
      endif

      call h5open_f(error)
      if (master) then
         call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, error)

         call h5ltget_attribute_double_f(file_id,"/","piernik", rbuf, error)
         if (error /= 0) call die("[dataio_hdf5:read_restart_hdf5] Cannot read 'piernik' attribute from the restart file. The file may be either damaged or incompatible")
         if (rbuf(1) > piernik_hdf5_version) call die("[dataio_hdf5:read_restart_hdf5] Cannot read future versions of the restart file. Upgrade your Piernik code and try again.")
         if (int(rbuf(1)) < int(piernik_hdf5_version)) call die("[dataio_hdf5:read_restart_hdf5] The restart file is too ancient. Try to temporarily disable this statement if you are sure that it can be interpreted correctly.")
         if (rbuf(1) < piernik_hdf5_version) call warn("[dataio_hdf5:read_restart_hdf5] Old versions of the restart file may not always work fully correctly.")

         call h5ltget_attribute_int_f(file_id,"/","nxd", ibuf,error)
         if (ibuf(1) /= dom%n_d(xdim) .or. error /= 0) call die("[dataio_hdf5:read_restart_hdf5] nxd does not match")
         if (has_dir(xdim)) then
            call h5ltget_attribute_double_f(file_id,"/","xmin", rbuf, error)
            if (rbuf(1) /= dom%xmin .or. error /= 0) call die("[dataio_hdf5:read_restart_hdf5] xmin does not match")
            call h5ltget_attribute_double_f(file_id,"/","xmax", rbuf, error)
            if (rbuf(1) /= dom%xmax .or. error /= 0) call die("[dataio_hdf5:read_restart_hdf5] xmax does not match")
         endif

         call h5ltget_attribute_int_f(file_id,"/","nyd", ibuf,error)
         if (ibuf(1) /= dom%n_d(ydim) .or. error /= 0) call die("[dataio_hdf5:read_restart_hdf5] nyd does not match")
         if (has_dir(ydim)) then
            call h5ltget_attribute_double_f(file_id,"/","ymin", rbuf, error)
            if (rbuf(1) /= dom%ymin .or. error /= 0) call die("[dataio_hdf5:read_restart_hdf5] ymin does not match")
            call h5ltget_attribute_double_f(file_id,"/","ymax", rbuf, error)
            if (rbuf(1) /= dom%ymax .or. error /= 0) call die("[dataio_hdf5:read_restart_hdf5] ymax does not match")
         endif

         call h5ltget_attribute_int_f(file_id,"/","nzd", ibuf,error)
         if (ibuf(1) /= dom%n_d(zdim) .or. error /= 0) call die("[dataio_hdf5:read_restart_hdf5] nzd does not match")
         if (has_dir(zdim)) then
            call h5ltget_attribute_double_f(file_id,"/","zmin", rbuf, error)
            if (rbuf(1) /= dom%zmin .or. error /= 0) call die("[dataio_hdf5:read_restart_hdf5] zmin does not match")
            call h5ltget_attribute_double_f(file_id,"/","zmax", rbuf, error)
            if (rbuf(1) /= dom%zmax .or. error /= 0) call die("[dataio_hdf5:read_restart_hdf5] zmax does not match")
         endif

         call h5fclose_f(file_id, error)
      endif

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
      call h5pset_fapl_mpio_f(plist_id, comm3d, info, error)

      call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, error, access_prp = plist_id)
      call h5pclose_f(plist_id, error)

      if (associated(problem_read_restart)) call problem_read_restart(file_id)

#ifdef ISO_LOCAL
      if (.not.associated(p3d)) p3d => cs_iso2_arr(:,:,:)
      call read_3darr_from_restart(file_id, "cs_iso2", p3d, cg%nx, cg%ny, cg%nz)
      if (associated(p3d)) nullify(p3d)
#endif /* ISO_LOCAL */
#ifdef GRAV
      if (.not.associated(p3d)) p3d => gp(:,:,:)
      call read_3darr_from_restart(file_id, "gp", p3d, cg%nx, cg%ny, cg%nz)
      if (associated(p3d)) nullify(p3d)
#endif /* GRAV */
      !> \deprecated The code for reading fluid data, mag field and read_3darr_from_restart is almost replicated.
      !> \todo Try to put it in a separate subroutine.
      !----------------------------------------------------------------------------------
      !  READ FLUID VARIABLES
      !

      rank = 4
      allocate(dimsf(rank), dimsfi(rank), chunk_dims(rank))
      allocate(block(rank), offset(rank), count(rank), stride(rank))
      dimsf = [nu, dom%n_d(:)] ! Dataset dimensions
      dimsfi = dimsf
      chunk_dims = [nu, cg%n_b(:)]
      if (.not.associated(p4d)) p4d => u(:,cg%is:cg%ie,cg%js:cg%je,cg%ks:cg%ke)

      ! Create chunked dataset.
      call h5dopen_f(file_id, dname(1), dset_id, error)

      call h5dget_space_f(dset_id, filespace, error)
      call H5sget_simple_extent_ndims_f (filespace,rank,error)
      call H5dget_create_plist_f (dset_id,plist_id,error)

      ! Each process defines dataset in memory and writes it to the hyperslab
      ! in the file.
      stride(:) = 1
      count(:)  = 1
      block(:)  = chunk_dims(:)

      offset(:)   = [0, cg%off(:)]

      ! Select hyperslab in the file.
      call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error, stride, block)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
      call h5screate_simple_f(rank, chunk_dims, memspace, error)
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, p4d, dimsfi, error, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

      call h5sclose_f(filespace, error)
      call h5sclose_f(memspace, error)
      call h5pclose_f(plist_id, error)
      call h5dclose_f(dset_id, error)
      if (associated(p4d)) nullify(p4d)
      !----------------------------------------------------------------------------------

      !----------------------------------------------------------------------------------
      !  READ MAG VARIABLES
      !
      dimsf = [3, dom%nxt, dom%nyt, dom%nzt] ! Dataset dimensions - we do need the whole domain of mag data
      dimsfi = dimsf
      chunk_dims = [3, cg%nx, cg%ny, cg%nz]
      if (.not.associated(p4d)) p4d => b(:,:,:,:)

      ! Create chunked dataset.
      call h5dopen_f(file_id, dname(2), dset_id, error)

      call h5dget_space_f(dset_id, filespace, error)
      call H5Sget_simple_extent_ndims_f (filespace,rank,error)
      call H5Dget_create_plist_f (dset_id,plist_id,error)

      ! Each process defines dataset in memory and writes it to the hyperslab
      ! in the file.
      stride(:) = 1
      count(:)  = 1
      block(:)  = chunk_dims(:)

      offset(:) = [0, cg%off(:)]

      ! Select hyperslab in the file.
      call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error, stride, block)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
      call h5screate_simple_f(rank, chunk_dims, memspace, error)
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, p4d, dimsfi, error, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

      call h5sclose_f(filespace, error)
      call h5sclose_f(memspace, error)
      call h5pclose_f(plist_id, error)
      call h5dclose_f(dset_id, error)
      if (associated(p4d)) nullify(p4d)
      !----------------------------------------------------------------------------------
      call h5fclose_f(file_id, error)
      if (allocated(dimsf))      deallocate(dimsf)
      if (allocated(dimsfi))     deallocate(dimsfi)
      if (allocated(chunk_dims)) deallocate(chunk_dims)
      if (allocated(count))      deallocate(count)
      if (allocated(offset))     deallocate(offset)
      if (allocated(stride))     deallocate(stride)
      if (allocated(block))      deallocate(block)
      if (master) then
         call h5fopen_f (filename, H5F_ACC_RDONLY_F, file_id, error)
         bufsize = 1
         call h5ltget_attribute_double_f(file_id,"/","piernik", rbuf,error)
         restart_hdf5_version = rbuf(1)
         call h5ltget_attribute_double_f(file_id,"/","time", rbuf,error)
         t = rbuf(1)
         call h5ltget_attribute_double_f(file_id,"/","timestep", rbuf,error)
         dt = rbuf(1)
         call h5ltget_attribute_double_f(file_id,"/","magic_mass", rbuf,error)
         magic_mass = rbuf(1)
         call h5ltget_attribute_int_f(file_id,"/","nstep", ibuf,error)
         chdf%nstep = ibuf(1)
         call h5ltget_attribute_int_f(file_id,"/","nres", ibuf,error)
         chdf%nres = ibuf(1)
         call h5ltget_attribute_int_f(file_id,"/","nhdf", ibuf,error)
         chdf%nhdf = ibuf(1)
         call h5ltget_attribute_int_f(file_id,"/","step_res", ibuf,error)
         chdf%step_res = ibuf(1)
         call h5ltget_attribute_int_f(file_id,"/","step_hdf", ibuf,error)
         chdf%step_hdf = ibuf(1)
         call h5ltget_attribute_double_f(file_id,"/","next_t_tsl", rbuf,error)
         chdf%next_t_tsl = rbuf(1)
         call h5ltget_attribute_double_f(file_id,"/","next_t_log", rbuf,error)
         chdf%next_t_log = rbuf(1)
         call h5ltget_attribute_double_f(file_id,"/","last_hdf_time", rbuf,error)
         chdf%last_hdf_time = rbuf(1)

         call h5ltget_attribute_string_f(file_id,"/","problem_name", problem_name,error)
         call h5ltget_attribute_string_f(file_id,"/","domain", chdf%domain,error)
         call h5ltget_attribute_string_f(file_id,"/","run_id", chdf%new_id,error)

         if (restart_hdf5_version > 1.11) then
            call h5ltget_attribute_int_f(file_id,"/","require_init_prob", ibuf,error)
            require_init_prob = ibuf(1)
         endif

         problem_name = fix_string(problem_name)   !> \deprecated BEWARE: >=HDF5-1.8.4 has weird issues with strings
         chdf%new_id  = fix_string(chdf%new_id)    !> \deprecated   this bit hacks it around
         chdf%domain  = fix_string(chdf%domain)

         call h5fclose_f(file_id, error)

         write(msg,'(2a)') 'Done reading restart file: ',trim(filename)
         call printio(msg)
      endif
      call h5close_f(error)

      call MPI_Bcast(restart_hdf5_version,    1, MPI_DOUBLE_PRECISION, 0, comm3d, ierr)

      call MPI_Bcast(chdf%nstep,    1, MPI_INTEGER, 0, comm3d, ierr)
      call MPI_Bcast(chdf%nres,     1, MPI_INTEGER, 0, comm3d, ierr)
      call MPI_Bcast(chdf%nhdf,     1, MPI_INTEGER, 0, comm3d, ierr)
      call MPI_Bcast(chdf%step_res, 1, MPI_INTEGER, 0, comm3d, ierr)
      call MPI_Bcast(chdf%step_hdf, 1, MPI_INTEGER, 0, comm3d, ierr)
      if (restart_hdf5_version > 1.11) call MPI_Bcast(require_init_prob, 1, MPI_INTEGER, 0, comm3d, ierr)

      call MPI_Bcast(chdf%next_t_tsl,    1, MPI_DOUBLE_PRECISION, 0, comm3d, ierr)
      call MPI_Bcast(chdf%next_t_log,    1, MPI_DOUBLE_PRECISION, 0, comm3d, ierr)
      call MPI_Bcast(chdf%last_hdf_time, 1, MPI_DOUBLE_PRECISION, 0, comm3d, ierr)
      call MPI_Bcast(t,                  1, MPI_DOUBLE_PRECISION, 0, comm3d, ierr)
      call MPI_Bcast(dt,                 1, MPI_DOUBLE_PRECISION, 0, comm3d, ierr)

      call MPI_Bcast(problem_name, cbuff_len, MPI_CHARACTER, 0, comm3d, ierr)
      call MPI_Bcast(chdf%domain,  domlen,    MPI_CHARACTER, 0, comm3d, ierr)
      call MPI_Bcast(chdf%new_id,  idlen,     MPI_CHARACTER, 0, comm3d, ierr)

   end subroutine read_restart_hdf5
!
! ------------------------------------------------------------------------------------
!
   subroutine write_hdf5(chdf)
      use dataio_pub,    only: cwdlen, msg, die, user_vars_hdf5, nhdf, varlen, problem_name, run_id
      use grid,          only: cg
      use hdf5,          only: HID_T, H5F_ACC_TRUNC_F, H5P_FILE_ACCESS_F, H5P_DEFAULT_F, &
           &                   h5open_f, h5close_f, h5fcreate_f, h5fclose_f, h5pcreate_f, h5pclose_f, h5pset_fapl_mpio_f
      use mpisetup,      only: comm3d, ierr, info
      use types,         only: hdf
#ifdef NEW_HDF5
      use list_hdf5,     only: iterate_lhdf5
#endif /* NEW_HDF5 */

      implicit none

      type(hdf), intent(in)   :: chdf
      integer(HID_T)          :: file_id       ! File identifier
      integer(HID_T)          :: plist_id      ! Property list identifier
      integer                 :: ierrh, error, i
      logical                 :: ok_var
      character(len=varlen)   :: dd
      character(len=cwdlen)   :: fname

      real(kind=4), allocatable :: data (:,:,:)  ! Data to write

      ! Initialize HDF5 library and Fortran interfaces.
      !
      write(dd,'(i4.4)') chdf%nhdf
      write(fname, '(6a)') trim(problem_name),"_",trim(run_id),"_",dd,".h5"

      call h5open_f(error)
      !
      ! Setup file access property list with parallel I/O access.
      !
      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
      call h5pset_fapl_mpio_f(plist_id, comm3d, info, error)
      !
      ! Create the file collectively.
      !
      call h5fcreate_f(trim(fname), H5F_ACC_TRUNC_F, file_id, error, creation_prp = H5P_DEFAULT_F, access_prp = plist_id)
      call h5pclose_f(plist_id, error)
      if (.not.allocated(data)) allocate(data(cg%nxb, cg%nyb, cg%nzb))
      do i = 1, nhdf_vars
         ierrh = 0; ok_var = .false.
         call common_vars_hdf5(hdf_vars(i),data,ierrh)
         if (associated(user_vars_hdf5) .and. ierrh /= 0) call user_vars_hdf5(hdf_vars(i),data,ierrh)
         if (ierrh>=0) ok_var = .true.
         if (.not.ok_var) then
            write(msg,'(3a)') "[dataio_hdf5:write_hdf5]: ",hdf_vars(i)," is not defined in common_vars_hdf5, neither in user_vars_hdf5."
            call die(msg)
         endif
         call write_arr(data,hdf_vars(i),file_id)
      enddo
      if (allocated(data)) deallocate(data)
      !
      ! Close the property list.
      !
#ifdef NEW_HDF5
      call iterate_lhdf5(file_id)
#endif /* NEW_HDF5 */

      !
      ! Close the property list.
      !
      call h5fclose_f(file_id, error)

      call set_common_attributes(fname, chdf, "output")

      call MPI_Barrier(comm3d,ierr)
      call h5close_f(error)

      nhdf = nhdf + 1

   end subroutine write_hdf5

   subroutine write_arr(data,dsetname,file_id)

      use dataio_pub,    only: varlen
      use grid,          only: cg
      use hdf5,          only: HID_T, HSIZE_T, HSSIZE_T, h5screate_simple_f, h5pcreate_f, h5pset_chunk_f, &
           &                   h5sclose_f, h5pset_dxpl_mpio_f, h5dwrite_f, h5dclose_f, H5P_DATASET_XFER_F, H5P_DATASET_CREATE_F, &
           &                   H5T_NATIVE_REAL, H5S_SELECT_SET_F, H5FD_MPIO_INDEPENDENT_F, h5dcreate_f, h5dget_space_f, &
           &                   h5sselect_hyperslab_f
      use mpisetup,      only: dom

      implicit none
      real(kind=4), dimension(:,:,:) :: data

      integer, parameter :: rank = 3       ! Dataset rank

      character(len=varlen) :: dsetname       ! Dataset name

      integer(HID_T)        :: file_id        ! Dataset identifier
      integer(HID_T)        :: dset_id        ! Dataset identifier
      integer(HID_T)        :: plist_id       ! Dataset identifier
      integer(HID_T)        :: filespace      ! Dataspace identifier in file
      integer(HID_T)        :: memspace       ! Dataspace identifier in memory

      integer, parameter :: ndims = 3
      integer(HSIZE_T),  DIMENSION(ndims) :: count
      integer(HSSIZE_T), DIMENSION(ndims) :: offset
      integer(HSIZE_T),  DIMENSION(ndims) :: stride
      integer(HSIZE_T),  DIMENSION(ndims) :: block
      integer(HSIZE_T),  DIMENSION(ndims) :: dimsf, dimsfi, chunk_dims
      integer :: error

      chunk_dims = cg%n_b(:) ! Chunks dimensions
      dimsf  = dom%n_d(:)    ! Dataset dimensions
      dimsfi = dimsf
      !
      ! Create the data space for the  dataset.
      !
      call h5screate_simple_f(rank, dimsf, filespace, error)
      call h5screate_simple_f(rank, chunk_dims, memspace, error)

      !
      ! Create chunked dataset.
      !
      call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
      call h5pset_chunk_f(plist_id, rank, chunk_dims, error)
      call h5dcreate_f(file_id, dsetname, H5T_NATIVE_REAL, filespace, dset_id, error, plist_id)
      call h5sclose_f(filespace, error)

      !
      ! Each process defines dataset in memory and writes it to the hyperslab
      ! in the file.
      !
      stride(:) = 1
      count(:) =  1
      block(:) = chunk_dims(:)

      offset(:) = cg%off(:)
      !
      ! Select hyperslab in the file.
      !
      call h5dget_space_f(dset_id, filespace, error)
      call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error, stride, block)

      !
      ! Create property list for collective dataset write
      !
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)

      !
      ! Write the dataset collectively.
      !
      call h5dwrite_f(dset_id, H5T_NATIVE_REAL, data, dimsfi, error, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

      !
      ! Close dataspaces.
      !
      call h5sclose_f(filespace, error)
      call h5sclose_f(memspace, error)
      !
      ! Close the dataset.
      !
      call h5dclose_f(dset_id, error)

   end subroutine write_arr

!>
!! \brief This routine writes all attributes that are common to restart and output files.
!! \details Other common elements may also be moved here.
!<
   subroutine set_common_attributes(filename, chdf, stype)

      use dataio_pub,    only: msg, printio, require_init_prob, piernik_hdf5_version, problem_name, run_id
      use grid,          only: cg
      use hdf5,          only: HID_T, SIZE_T, HSIZE_T, H5F_ACC_RDWR_F, H5T_NATIVE_CHARACTER, H5Z_FILTER_DEFLATE_F, H5P_DATASET_CREATE_F, &
           &                   h5fopen_f, h5fclose_f, H5Zfilter_avail_f, H5Pcreate_f, H5Pset_deflate_f, H5Pset_chunk_f, &
           &                   h5tcopy_f, h5tset_size_f, h5screate_simple_f, H5Dcreate_f, H5Dwrite_f, H5Dclose_f, H5Sclose_f, H5Tclose_f, H5Pclose_f
      use h5lt,          only: h5ltset_attribute_double_f, h5ltset_attribute_int_f, h5ltset_attribute_string_f
      use list_hdf5,     only: additional_attrs
      use mpisetup,      only: slave, t, dt, cbuff_len, local_magic_mass, comm, ierr, magic_mass, dom
      use mpi,           only: MPI_DOUBLE_PRECISION, MPI_SUM
      use types,         only: hdf
      use version,       only: env, nenv

      implicit none

      character(len=*), intent(in) :: filename  !> HDF File name
      character(len=*), intent(in) :: stype     !> "output" or "restart"
      type(hdf), intent(in)        :: chdf

      integer(HID_T)                 :: file_id       !> File identifier
      integer(HID_T)                 :: type_id, dspace_id, dset_id, prp_id
      integer(HSIZE_T), dimension(1) :: dimstr
      logical                        :: Z_avail
      integer                        :: fe, i
      integer(SIZE_T)                :: bufsize, maxlen
      integer                        :: error
      real                           :: magic_mass0
      integer, parameter             :: buf_len = 50
      integer, dimension(buf_len)    :: ibuffer
      real,    dimension(buf_len)    :: rbuffer
      character(len=cbuff_len), dimension(buf_len) :: ibuffer_name = ''
      character(len=cbuff_len), dimension(buf_len) :: rbuffer_name = ''

      call MPI_Reduce(local_magic_mass, magic_mass0, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
      magic_mass       = magic_mass + magic_mass0
      local_magic_mass = 0.0

      if (slave) return ! This data need not be written in parallel.

      !! The rr1 marks critical attributes that are read by read_restart_hdf5 and compared against value read from the problem.par file.
      !! The rr2 marks runtime values that are read by read_restart_hdf5 and assigned to something in the code.
      !> \todo Set up an universal table(s) of attribute names for use by both set_common_attributes and read_restart_hdf5.
      !! Provide indices for critical attributes (rr1) and for runtime attributes (rr2).
      !<

      call h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, error)

      rbuffer(1)  = t                        ; rbuffer_name(1)  = "time" !rr2
      rbuffer(2)  = dt                       ; rbuffer_name(2)  = "timestep" !rr2
      rbuffer(3)  = chdf%last_hdf_time       ; rbuffer_name(3)  = "last_hdf_time" !rr2
      rbuffer(4)  = dom%xmin                 ; rbuffer_name(4)  = "xmin" !rr1
      rbuffer(5)  = dom%xmax                 ; rbuffer_name(5)  = "xmax" !rr1
      rbuffer(6)  = dom%ymin                 ; rbuffer_name(6)  = "ymin" !rr1
      rbuffer(7)  = dom%ymax                 ; rbuffer_name(7)  = "ymax" !rr1
      rbuffer(8)  = dom%zmin                 ; rbuffer_name(8)  = "zmin" !rr1
      rbuffer(9)  = dom%zmax                 ; rbuffer_name(9)  = "zmax" !rr1
      rbuffer(10) = piernik_hdf5_version     ; rbuffer_name(10) = "piernik" !rr1,rr2
      rbuffer(11) = magic_mass               ; rbuffer_name(11) = "magic_mass" !rr2
      rbuffer(12) = chdf%next_t_tsl          ; rbuffer_name(12) = "next_t_tsl" !rr2
      rbuffer(13) = chdf%next_t_log          ; rbuffer_name(13) = "next_t_log" !rr2

      ibuffer(1)   = chdf%nstep              ; ibuffer_name(1)   = "nstep" !rr2
      ibuffer(2)   = chdf%nres+1             ; ibuffer_name(2)   = "nres" !rr2
      ibuffer(3)   = chdf%nhdf               ; ibuffer_name(3)   = "nhdf" !rr2
      ibuffer(4)   = chdf%nstep              ; ibuffer_name(4)   = "step_res" !rr2
      ibuffer(5)   = chdf%step_hdf           ; ibuffer_name(5)   = "step_hdf" !rr2
      ibuffer(6:8) = dom%n_d(:)              ; ibuffer_name(6:8) = [ "nxd", "nyd", "nzd" ] !rr1
      ibuffer(9)   = cg%nb                   ; ibuffer_name(9)   = "nb"
      ibuffer(10)  = require_init_prob       ; ibuffer_name(10)  = "require_init_prob" !rr2

      bufsize = 1

      i = 1
      do while (rbuffer_name(i) /= "")
         call h5ltset_attribute_double_f(file_id, "/", rbuffer_name(i), rbuffer(i), bufsize, error)
         i = i+1
      enddo

      i = 1
      do while (ibuffer_name(i) /= "")
         call h5ltset_attribute_int_f(file_id, "/", ibuffer_name(i), ibuffer(i), bufsize, error)
         i = i+1
      enddo

      ! Store a compressed copy of the problem.par file.
      maxlen = maxval(len_trim(parfile(:parfilelines)))
      dimstr = [parfilelines]
      call H5Zfilter_avail_f(H5Z_FILTER_DEFLATE_F, Z_avail, error)
      ! call H5Zget_filter_info_f ! everything should be always fine for gzip
      call H5Pcreate_f(H5P_DATASET_CREATE_F, prp_id, error)
      if (Z_avail) then
         call H5Pset_deflate_f(prp_id, 9, error)
         call H5Pset_chunk_f(prp_id, 1, dimstr, error)
      endif
      call H5Tcopy_f(H5T_NATIVE_CHARACTER, type_id, error)
      call H5Tset_size_f(type_id, maxlen, error)
      call H5Screate_simple_f(1, dimstr, dspace_id, error)
      call H5Dcreate_f(file_id, "problem.par", type_id,  dspace_id, dset_id, error, dcpl_id = prp_id)
      call H5Dwrite_f(dset_id, type_id, parfile(:)(:maxlen), dimstr, error)
      call H5Dclose_f(dset_id, error)
      call H5Sclose_f(dspace_id, error)

      ! Store a compressed copy of the piernik.def file and Id lines from source files.
      ! We recycle type_id and prp_id, so we don't close them yet.
      maxlen = maxval(len_trim(env(:nenv)))
      dimstr = [nenv]
      if (Z_avail) call H5Pset_chunk_f(prp_id, 1, dimstr, error)
      call H5Tset_size_f(type_id, maxlen, error)
      call H5Screate_simple_f(1, dimstr, dspace_id, error)
      call H5Dcreate_f(file_id, "env", type_id,  dspace_id, dset_id, error, dcpl_id = prp_id)
      call H5Dwrite_f(dset_id, type_id, env(:)(:maxlen), dimstr, error)
      call H5Dclose_f(dset_id, error)
      call H5Sclose_f(dspace_id, error)
      call H5Tclose_f(type_id, error)
      call H5Pclose_f(prp_id, error)

      !bufsize = 3
      !call h5ltset_attribute_int_f(file_id, "/", "psize", psize, bufsize, error) ! unused and will be obsolete soon

      fe = len_trim(problem_name)
      call h5ltset_attribute_string_f(file_id, "/", "problem_name", problem_name(1:fe), error) !rr2
      fe = len_trim(chdf%domain)
      call h5ltset_attribute_string_f(file_id, "/", "domain", chdf%domain(1:fe), error) !rr2
      fe = len_trim(run_id)
      call h5ltset_attribute_string_f(file_id, "/", "run_id", run_id(1:fe), error) !rr2

      call additional_attrs(file_id)

      call h5fclose_f(file_id, error)

      write(msg,'(4a)') 'Writing ',stype,' file: ',trim(filename)
      call printio(msg)

      ! only master process exits here

   end subroutine set_common_attributes

end module dataio_hdf5
