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

#define RNG is:ie, js:je, ks:ke
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

   real, parameter    :: piernik_hdf5_version = 1.1   !< output version
   integer, parameter :: dnamelen=10
   character(LEN=dnamelen), dimension(2) :: dname = (/"fluid     ","mag       "/)  !< dataset names for restart files
   character(len=S_LEN), allocatable, dimension(:) :: hdf_vars  !< dataset names for hdf files
   integer :: nhdf_vars !< number of quantities plotted to hdf files
   integer :: ix !< no. of cell (1 <= ix < nxd) for YZ slice in plt files
   integer :: iy !< no. of cell (1 <= iy < nyd) for XZ slice in plt files
   integer :: iz !< no. of cell (1 <= iz < nzd) for XY slice in plt files
   integer :: is !< COMMENT ME
   integer :: ie !< COMMENT ME
   integer :: js !< COMMENT ME
   integer :: je !< COMMENT ME
   integer :: ks !< COMMENT ME
   integer :: ke !< COMMENT ME
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
      use grid,          only: nx, ny, nz, has_dir, xdim, ydim, zdim, nb
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

      if (has_dir(xdim)) then
         is = nb+1; ie = nx-nb
      else
         is = 1; ie = 1
      endif

      if (has_dir(ydim)) then
         js = nb+1; je = ny-nb
      else
         js = 1; je = 1
      endif

      if (has_dir(zdim)) then
         ks = nb+1; ke = nz-nb
      else
         ks = 1; ke = 1
      endif

      nvars = 1
      do while ( LEN(TRIM(vars(nvars))) > 1)
        nvars = nvars + 1
      enddo
      nvars = nvars - 1

      nhdf_vars = 0
      do i = 1, nvars
         select case (vars(i))
            case ('dens')
               nhdf_vars = nhdf_vars + SIZE(iarr_all_dn,1)
            case ('velx')
               nhdf_vars = nhdf_vars + SIZE(iarr_all_mx,1)
            case ('vely')
               nhdf_vars = nhdf_vars + SIZE(iarr_all_my,1)
            case ('velz')
               nhdf_vars = nhdf_vars + SIZE(iarr_all_mz,1)
            case ('ener')
               nhdf_vars = nhdf_vars + SIZE(iarr_all_mz,1)
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
               nhdf_vars = nhdf_vars + SIZE(iarr_all_crs,1)
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

!<
!! \brief Routine calculating quantities for plot files
!<
   subroutine common_plt_hdf5(var,ij,xn,tab,ierrh)
      use arrays,        only: u, b
      use dataio_pub,    only: varlen, planelen
      use grid,          only: nb, nyb, nzb, nxb
#ifdef GRAV
      use arrays,        only: gpot
#endif /* GRAV */
      use fluidindex,    only: nvar, ibx, iby, ibz
#ifdef COSM_RAYS
      use fluidindex,    only: iarr_all_crs
#endif /* COSM_RAYS */

      implicit none
      character(LEN=varlen)   :: var !< quantity to be plotted
      character(LEN=planelen) :: ij  !< plane of plot
      integer                 :: xn  !< no. of cell at which we are slicing the local block
      integer                 :: ierrh !< error handling
      real, dimension(:,:)    :: tab !< array containing given quantity
#ifdef COSM_RAYS
      integer                 :: i
#endif /* COSM_RAYS */

      ierrh = 0
      select case (var)
         case ("dend")
            if (ij=="yz") tab(:,:) = u(nvar%dst%idn,xn,nb+1:nyb+nb,nb+1:nzb+nb)
            if (ij=="xz") tab(:,:) = u(nvar%dst%idn,nb+1:nxb+nb,xn,nb+1:nzb+nb)
            if (ij=="xy") tab(:,:) = u(nvar%dst%idn,nb+1:nxb+nb,nb+1:nyb+nb,xn)
         case ("denn")
            if (ij=="yz") tab(:,:) = u(nvar%neu%idn,xn,nb+1:nyb+nb,nb+1:nzb+nb)
            if (ij=="xz") tab(:,:) = u(nvar%neu%idn,nb+1:nxb+nb,xn,nb+1:nzb+nb)
            if (ij=="xy") tab(:,:) = u(nvar%neu%idn,nb+1:nxb+nb,nb+1:nyb+nb,xn)
         case ("deni")
            if (ij=="yz") tab(:,:) = u(nvar%ion%idn,xn,nb+1:nyb+nb,nb+1:nzb+nb)
            if (ij=="xz") tab(:,:) = u(nvar%ion%idn,nb+1:nxb+nb,xn,nb+1:nzb+nb)
            if (ij=="xy") tab(:,:) = u(nvar%ion%idn,nb+1:nxb+nb,nb+1:nyb+nb,xn)
         case ("vlxd")
            if (ij=="yz") tab(:,:) = u(nvar%dst%imx,xn,nb+1:nyb+nb,nb+1:nzb+nb) / &
                           u(nvar%dst%idn,xn,nb+1:nyb+nb,nb+1:nzb+nb)
            if (ij=="xz") tab(:,:) = u(nvar%dst%imx,nb+1:nxb+nb,xn,nb+1:nzb+nb) / &
                           u(nvar%dst%idn,nb+1:nxb+nb,xn,nb+1:nzb+nb)
            if (ij=="xy") tab(:,:) = u(nvar%dst%imx,nb+1:nxb+nb,nb+1:nyb+nb,xn) / &
                           u(nvar%dst%idn,nb+1:nxb+nb,nb+1:nyb+nb,xn)
         case ("vlxn")
            if (ij=="yz") tab(:,:) = u(nvar%neu%imx,xn,nb+1:nyb+nb,nb+1:nzb+nb) / &
                           u(nvar%neu%idn,xn,nb+1:nyb+nb,nb+1:nzb+nb)
            if (ij=="xz") tab(:,:) = u(nvar%neu%imx,nb+1:nxb+nb,xn,nb+1:nzb+nb) / &
                           u(nvar%neu%idn,nb+1:nxb+nb,xn,nb+1:nzb+nb)
            if (ij=="xy") tab(:,:) = u(nvar%neu%imx,nb+1:nxb+nb,nb+1:nyb+nb,xn) / &
                           u(nvar%neu%idn,nb+1:nxb+nb,nb+1:nyb+nb,xn)
         case ("vlxi")
            if (ij=="yz") tab(:,:) = u(nvar%ion%imx,xn,nb+1:nyb+nb,nb+1:nzb+nb) / &
                           u(nvar%ion%idn,xn,nb+1:nyb+nb,nb+1:nzb+nb)
            if (ij=="xz") tab(:,:) = u(nvar%ion%imx,nb+1:nxb+nb,xn,nb+1:nzb+nb) / &
                           u(nvar%ion%idn,nb+1:nxb+nb,xn,nb+1:nzb+nb)
            if (ij=="xy") tab(:,:) = u(nvar%ion%imx,nb+1:nxb+nb,nb+1:nyb+nb,xn) / &
                           u(nvar%ion%idn,nb+1:nxb+nb,nb+1:nyb+nb,xn)
         case ("vlyd")
            if (ij=="yz") tab(:,:) = u(nvar%dst%imy,xn,nb+1:nyb+nb,nb+1:nzb+nb) / &
                           u(nvar%dst%idn,xn,nb+1:nyb+nb,nb+1:nzb+nb)
            if (ij=="xz") tab(:,:) = u(nvar%dst%imy,nb+1:nxb+nb,xn,nb+1:nzb+nb) / &
                           u(nvar%dst%idn,nb+1:nxb+nb,xn,nb+1:nzb+nb)
            if (ij=="xy") tab(:,:) = u(nvar%dst%imy,nb+1:nxb+nb,nb+1:nyb+nb,xn) / &
                           u(nvar%dst%idn,nb+1:nxb+nb,nb+1:nyb+nb,xn)
         case ("vlyn")
            if (ij=="yz") tab(:,:) = u(nvar%neu%imy,xn,nb+1:nyb+nb,nb+1:nzb+nb) / &
                           u(nvar%neu%idn,xn,nb+1:nyb+nb,nb+1:nzb+nb)
            if (ij=="xz") tab(:,:) = u(nvar%neu%imy,nb+1:nxb+nb,xn,nb+1:nzb+nb) / &
                           u(nvar%neu%idn,nb+1:nxb+nb,xn,nb+1:nzb+nb)
            if (ij=="xy") tab(:,:) = u(nvar%neu%imy,nb+1:nxb+nb,nb+1:nyb+nb,xn) / &
                           u(nvar%neu%idn,nb+1:nxb+nb,nb+1:nyb+nb,xn)
         case ("vlyi")
            if (ij=="yz") tab(:,:) = u(nvar%ion%imy,xn,nb+1:nyb+nb,nb+1:nzb+nb) / &
                           u(nvar%ion%idn,xn,nb+1:nyb+nb,nb+1:nzb+nb)
            if (ij=="xz") tab(:,:) = u(nvar%ion%imy,nb+1:nxb+nb,xn,nb+1:nzb+nb) / &
                           u(nvar%ion%idn,nb+1:nxb+nb,xn,nb+1:nzb+nb)
            if (ij=="xy") tab(:,:) = u(nvar%ion%imy,nb+1:nxb+nb,nb+1:nyb+nb,xn) / &
                           u(nvar%ion%idn,nb+1:nxb+nb,nb+1:nyb+nb,xn)
         case ("vlzd")
            if (ij=="yz") tab(:,:) = u(nvar%dst%imz,xn,nb+1:nyb+nb,nb+1:nzb+nb) / &
                           u(nvar%dst%idn,xn,nb+1:nyb+nb,nb+1:nzb+nb)
            if (ij=="xz") tab(:,:) = u(nvar%dst%imz,nb+1:nxb+nb,xn,nb+1:nzb+nb) / &
                           u(nvar%dst%idn,nb+1:nxb+nb,xn,nb+1:nzb+nb)
            if (ij=="xy") tab(:,:) = u(nvar%dst%imz,nb+1:nxb+nb,nb+1:nyb+nb,xn) / &
                           u(nvar%dst%idn,nb+1:nxb+nb,nb+1:nyb+nb,xn)
         case ("vlzn")
            if (ij=="yz") tab(:,:) = u(nvar%neu%imz,xn,nb+1:nyb+nb,nb+1:nzb+nb) / &
                           u(nvar%neu%idn,xn,nb+1:nyb+nb,nb+1:nzb+nb)
            if (ij=="xz") tab(:,:) = u(nvar%neu%imz,nb+1:nxb+nb,xn,nb+1:nzb+nb) / &
                           u(nvar%neu%idn,nb+1:nxb+nb,xn,nb+1:nzb+nb)
            if (ij=="xy") tab(:,:) = u(nvar%neu%imz,nb+1:nxb+nb,nb+1:nyb+nb,xn) / &
                           u(nvar%neu%idn,nb+1:nxb+nb,nb+1:nyb+nb,xn)
         case ("vlzi")
            if (ij=="yz") tab(:,:) = u(nvar%ion%imz,xn,nb+1:nyb+nb,nb+1:nzb+nb) / &
                           u(nvar%ion%idn,xn,nb+1:nyb+nb,nb+1:nzb+nb)
            if (ij=="xz") tab(:,:) = u(nvar%ion%imz,nb+1:nxb+nb,xn,nb+1:nzb+nb) / &
                           u(nvar%ion%idn,nb+1:nxb+nb,xn,nb+1:nzb+nb)
            if (ij=="xy") tab(:,:) = u(nvar%ion%imz,nb+1:nxb+nb,nb+1:nyb+nb,xn) / &
                           u(nvar%ion%idn,nb+1:nxb+nb,nb+1:nyb+nb,xn)
         case ("enen")
#ifdef ISO
            if (ij=="yz") tab(:,:) = 0.5 * (                     &
                          u(nvar%neu%imx,xn,nb+1:nyb+nb,nb+1:nzb+nb)**2 &
                        + u(nvar%neu%imy,xn,nb+1:nyb+nb,nb+1:nzb+nb)**2 &
                        + u(nvar%neu%imz,xn,nb+1:nyb+nb,nb+1:nzb+nb)**2) / &
                             u(nvar%neu%idn,xn,nb+1:nyb+nb,nb+1:nzb+nb)
            if (ij=="xz") tab(:,:) = 0.5 * (                     &
                          u(nvar%neu%imx,nb+1:nxb+nb,xn,nb+1:nzb+nb)**2 &
                         +u(nvar%neu%imy,nb+1:nxb+nb,xn,nb+1:nzb+nb)**2 &
                         +u(nvar%neu%imz,nb+1:nxb+nb,xn,nb+1:nzb+nb)**2) / &
                             u(nvar%neu%idn,nb+1:nxb+nb,xn,nb+1:nzb+nb)
            if (ij=="xy") tab(:,:) = 0.5 * (                     &
                          u(nvar%neu%imx,nb+1:nxb+nb,nb+1:nyb+nb,xn)**2 &
                         +u(nvar%neu%imy,nb+1:nxb+nb,nb+1:nyb+nb,xn)**2 &
                         +u(nvar%neu%imz,nb+1:nxb+nb,nb+1:nyb+nb,xn)**2) / &
                             u(nvar%neu%idn,nb+1:nxb+nb,nb+1:nyb+nb,xn)
#else /* !ISO */
            if (ij=="yz") tab(:,:) = u(nvar%neu%ien,xn,nb+1:nyb+nb,nb+1:nzb+nb)
            if (ij=="xz") tab(:,:) = u(nvar%neu%ien,nb+1:nxb+nb,xn,nb+1:nzb+nb)
            if (ij=="xy") tab(:,:) = u(nvar%neu%ien,nb+1:nxb+nb,nb+1:nyb+nb,xn)
#endif /* !ISO */
         case ('prei')
            tab(:,:) = 0.0
#ifdef NEUTRAL
         case ("pren")
#ifdef ISO
            tab = 0.0
#else /* !ISO */
            if (ij=="xy") then
               tab(:,:) = real( u(nvar%neu%ien,nb+1:nxb+nb,nb+1:nyb+nb,xn) - &
                 0.5 *( u(nvar%neu%imx,nb+1:nxb+nb,nb+1:nyb+nb,xn)**2 + u(nvar%neu%imy,nb+1:nxb+nb,nb+1:nyb+nb,xn)**2 + &
                        u(nvar%neu%imz,nb+1:nxb+nb,nb+1:nyb+nb,xn)**2 ) / u(nvar%neu%idn,nb+1:nxb+nb,nb+1:nyb+nb,xn),4)*(nvar%neu%gam_1)
            endif
#endif /* !ISO */
#endif /* NEUTRAL */
         case ("enei")
#ifdef ISO
            if (ij=="yz") tab(:,:) = 0.5 * (                     &
                          u(nvar%ion%imx,xn,nb+1:nyb+nb,nb+1:nzb+nb)**2 &
                        + u(nvar%ion%imy,xn,nb+1:nyb+nb,nb+1:nzb+nb)**2 &
                        + u(nvar%ion%imz,xn,nb+1:nyb+nb,nb+1:nzb+nb)**2) / &
                             u(nvar%ion%idn,xn,nb+1:nyb+nb,nb+1:nzb+nb)
            if (ij=="xz") tab(:,:) = 0.5 * (                     &
                          u(nvar%ion%imx,nb+1:nxb+nb,xn,nb+1:nzb+nb)**2 &
                         +u(nvar%ion%imy,nb+1:nxb+nb,xn,nb+1:nzb+nb)**2 &
                         +u(nvar%ion%imz,nb+1:nxb+nb,xn,nb+1:nzb+nb)**2) / &
                             u(nvar%ion%idn,nb+1:nxb+nb,xn,nb+1:nzb+nb)
            if (ij=="xy") tab(:,:) = 0.5 * (                     &
                          u(nvar%ion%imx,nb+1:nxb+nb,nb+1:nyb+nb,xn)**2 &
                         +u(nvar%ion%imy,nb+1:nxb+nb,nb+1:nyb+nb,xn)**2 &
                         +u(nvar%ion%imz,nb+1:nxb+nb,nb+1:nyb+nb,xn)**2) / &
                             u(nvar%ion%idn,nb+1:nxb+nb,nb+1:nyb+nb,xn)
#else /* !ISO */
            if (ij=="yz") tab(:,:) = u(nvar%ion%ien,xn,nb+1:nyb+nb,nb+1:nzb+nb)
            if (ij=="xz") tab(:,:) = u(nvar%ion%ien,nb+1:nxb+nb,xn,nb+1:nzb+nb)
            if (ij=="xy") tab(:,:) = u(nvar%ion%ien,nb+1:nxb+nb,nb+1:nyb+nb,xn)
#endif /* !ISO */

         case ("magx")
            if (ij=="yz") tab(:,:) = b(ibx,xn,nb+1:nyb+nb,nb+1:nzb+nb)
            if (ij=="xz") tab(:,:) = b(ibx,nb+1:nxb+nb,xn,nb+1:nzb+nb)
            if (ij=="xy") tab(:,:) = b(ibx,nb+1:nxb+nb,nb+1:nyb+nb,xn)
         case ("magy")
            if (ij=="yz") tab(:,:) = b(iby,xn,nb+1:nyb+nb,nb+1:nzb+nb)
            if (ij=="xz") tab(:,:) = b(iby,nb+1:nxb+nb,xn,nb+1:nzb+nb)
            if (ij=="xy") tab(:,:) = b(iby,nb+1:nxb+nb,nb+1:nyb+nb,xn)
         case ("magz")
            if (ij=="yz") tab(:,:) = b(ibz,xn,nb+1:nyb+nb,nb+1:nzb+nb)
            if (ij=="xz") tab(:,:) = b(ibz,nb+1:nxb+nb,xn,nb+1:nzb+nb)
            if (ij=="xy") tab(:,:) = b(ibz,nb+1:nxb+nb,nb+1:nyb+nb,xn)
#ifdef GRAV
         case ("gpot")
            if (ij=="yz") tab(:,:) = gpot(xn,nb+1:nyb+nb,nb+1:nzb+nb)
            if (ij=="xz") tab(:,:) = gpot(nb+1:nxb+nb,xn,nb+1:nzb+nb)
            if (ij=="xy") tab(:,:) = gpot(nb+1:nxb+nb,nb+1:nyb+nb,xn)
#endif /* GRAV */
         case default
#ifdef COSM_RAYS
            if (var(1:3) == 'ecr') then
               i = iarr_all_crs(ichar(var(4:4))-48)
               if (ij=="yz") tab(:,:) = u(i,xn,nb+1:nyb+nb,nb+1:nzb+nb)
               if (ij=="xz") tab(:,:) = u(i,nb+1:nxb+nb,xn,nb+1:nzb+nb)
               if (ij=="xy") tab(:,:) = u(i,nb+1:nxb+nb,nb+1:nyb+nb,xn)
            else
#endif /* COSM_RAYS */
            ierrh = -1
#ifdef COSM_RAYS
            endif
#endif /* COSM_RAYS */
      end select

   end subroutine common_plt_hdf5

!<
!! \brief Routine calculating quantities for .hdf files
!<
   subroutine common_vars_hdf5(var,tab, ierrh)
      use arrays,        only: u, b
      use dataio_pub,    only: varlen
      use fluidindex,    only: nvar, ibx, iby, ibz
#ifdef GRAV
      use arrays,        only: gpot
#ifdef MULTIGRID
      use arrays,        only: sgp
#endif /* MULTIGRID */
#endif /* GRAV */

      implicit none
      character(LEN=varlen), intent(in) :: var
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
            read(var,'(A3,I1)') aux,i !BEWARE 0 <= i <= 9, no other indices can be dumped to hdf file
            tab(:,:,:) = real(u(nvar%crs%beg+i-1,RNG),4)
#endif /* COSM_RAYS */
         case ("dend")
            tab(:,:,:) = real(u(nvar%dst%idn,RNG),4)
         case ("denn")
            tab(:,:,:) = real(u(nvar%neu%idn,RNG),4)
         case ("deni")
            tab(:,:,:) = real(u(nvar%ion%idn,RNG),4)
         case ("vlxd")
            tab(:,:,:) = real(u(nvar%dst%imx,RNG) / u(nvar%dst%idn,RNG),4)
         case ("vlxn")
            tab(:,:,:) = real(u(nvar%neu%imx,RNG) / u(nvar%neu%idn,RNG),4)
         case ("vlxi")
            tab(:,:,:) = real(u(nvar%ion%imx,RNG) / u(nvar%ion%idn,RNG),4)
         case ("vlyd")
            tab(:,:,:) = real(u(nvar%dst%imy,RNG) / u(nvar%dst%idn,RNG),4)
         case ("vlyn")
            tab(:,:,:) = real(u(nvar%neu%imy,RNG) / u(nvar%neu%idn,RNG),4)
         case ("vlyi")
            tab(:,:,:) = real(u(nvar%ion%imy,RNG) / u(nvar%ion%idn,RNG),4)
         case ("vlzd")
            tab(:,:,:) = real(u(nvar%dst%imz,RNG) / u(nvar%dst%idn,RNG),4)
         case ("vlzn")
            tab(:,:,:) = real(u(nvar%neu%imz,RNG) / u(nvar%neu%idn,RNG),4)
         case ("vlzi")
            tab(:,:,:) = real(u(nvar%ion%imz,RNG) / u(nvar%ion%idn,RNG),4)
#ifdef NEUTRAL
         case ("pren")
#ifdef ISO
            tab = 0.0
#else /* !ISO */
            tab(:,:,:) = real( u(nvar%neu%ien,RNG) - &
              0.5 *( u(nvar%neu%imx,RNG)**2 + u(nvar%neu%imy,RNG)**2 + u(nvar%neu%imz,RNG)**2 ) / u(nvar%neu%idn,RNG),4)*(nvar%neu%gam_1)
#endif /* !ISO */
#endif /* NEUTRAL */
         case ("enen")
#ifdef ISO
            tab(:,:,:) = real(0.5 *( u(nvar%neu%imx,RNG)**2 + &
                                     u(nvar%neu%imy,RNG)**2 + &
                                     u(nvar%neu%imz,RNG)**2 ) &
                              / u(nvar%neu%idn,RNG),4)
#else /* !ISO */
            tab(:,:,:) = real(u(nvar%neu%ien,RNG),4)
#endif /* !ISO */
#ifdef IONIZED
         case ("enei")
#ifdef ISO
            tab(:,:,:) = real(0.5 *( u(nvar%ion%imx,RNG)**2 + &
                                     u(nvar%ion%imy,RNG)**2 + &
                                     u(nvar%ion%imz,RNG)**2 ) &
                              / u(nvar%ion%idn,RNG),4)
#else /* !ISO */
            tab(:,:,:) = real(u(nvar%ion%ien,RNG),4)
#endif /* !ISO */
         case ("prei")
#ifdef ISO
            tab = 0.0
#else /* !ISO */
            tab(:,:,:) = real( u(nvar%ion%ien,RNG) - &
              0.5 *( u(nvar%ion%imx,RNG)**2 + u(nvar%ion%imy,RNG)**2 + u(nvar%ion%imz,RNG)**2 ) / u(nvar%ion%idn,RNG),4)*(nvar%ion%gam_1)
            tab(:,:,:) = tab(:,:,:) - real( 0.5*(nvar%ion%gam_1)*(b(ibx,RNG)**2 + &
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
      use mpisetup,      only: t, comm3d, ierr, proc

      implicit none

      integer, save     :: nimg = 0, error = 0
      real, save        :: last_plt_time = 0.0
      character(LEN=cwdlen) :: fname
      integer           :: i, fe
      logical, save     :: first_entry = .true.
      integer(HID_T)    :: file_id                 !> File identifier
      integer(HID_T)    :: gr_id, gr2_id           !> Groups identifier

      if ( ((t-last_plt_time > dt_plt) .or. first_entry ) .and. dt_plt > 0.0 ) then
         fe = len_trim(log_file)
         write(fname,'(2a)') trim(log_file(1:fe-3)),"plt"
         call H5open_f(error)

         if (proc==0 .and. first_entry) then
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
      use grid,          only: nxb, nyb, nzb, nb, has_dir, xdim, ydim, zdim
      use hdf5,          only: HID_T, HSIZE_T, SIZE_T, H5F_ACC_RDWR_F, h5fopen_f, h5gopen_f, h5gclose_f, h5fclose_f
      use h5lt,          only: h5ltmake_dataset_double_f, h5ltset_attribute_double_f
      use mpisetup,      only: comm3d, ierr, pxsize, pysize, pzsize, t, pcoords
      use mpi,           only: MPI_CHARACTER, MPI_DOUBLE_PRECISION
#ifdef PGPLOT
      use viz,           only: draw_me
#endif /* PGPLOT */

      implicit none

      character(LEN=planelen), intent(in) :: plane
      character(LEN=varlen), intent(in)   :: var                           !> not yet implemented
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
      character(LEN=suffixed_planelen)    :: pij
      character(LEN=cwdlen)               :: fname
      integer, parameter                  :: vdn_len = 12
      character(LEN=vdn_len)              :: vdname

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
               xn     = ix + nb - pcoords(1)*nxb
            else
               xn     = 1
            endif
            remain = (/.false.,.true.,.true./)
            pij    = "yz_"
            nib    = nyb
            nid    = nyb*pysize
            njb    = nzb
            njd    = nzb*pzsize
            nkb    = nxb
            pisize = pysize
            pjsize = pzsize
         case ("xz")
            if (has_dir(ydim)) then
               xn     = iy + nb - pcoords(2)*nyb
            else
               xn     = 1
            endif
            remain = (/.true.,.false.,.true./)
            pij    = "xz_"
            nib    = nxb
            nid    = nxb*pxsize
            njb    = nzb
            njd    = nzb*pzsize
            nkb    = nyb
            pisize = pxsize
            pjsize = pzsize
         case ("xy")
            if (has_dir(zdim)) then
               xn = iz + nb - pcoords(3)*nzb
            else
               xn = 1
            endif
            remain = (/.true.,.true.,.false./)
            pij    = "xy_"
            nib    = nxb
            nid    = nxb*pxsize
            njb    = nyb
            njd    = nyb*pysize
            nkb    = nzb
            pisize = pxsize
            pjsize = pysize
         case default
            call die("[dataio_hdf5:write_plot_hdf5] nonrecognized plane")
      end select

      dims(1) = nid
      dims(2) = njd
      call MPI_Barrier(comm3d,ierr)
      call MPI_Cart_sub(comm3d,remain,comm2d,ierr)
      call MPI_Comm_size(comm2d, ls, ierr)
      call MPI_Comm_rank(comm2d, lp, ierr)
      if ((xn > nb .and. xn <= nkb+nb).or.xn == 1) then
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
         call MPI_Gather(send, nib*njb, MPI_DOUBLE_PRECISION, &
                         temp, nib*njb, MPI_DOUBLE_PRECISION, &
                         0, comm2d,ierr)

         if (lp == 0) then
            imax = maxval(temp); imin = minval(temp)
            di = imax-imin
            do i = 0, pisize-1
               do j = 0, pjsize-1
                  img(i*nib+1:(i+1)*nib,j*njb+1:(j+1)*njb) = &
                    temp(:,:,(j+1)+i*pjsize)
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

!<
!! \brief This routine writes restart dump and updates restart counter
!<

   subroutine write_restart_hdf5

      use arrays,        only: u, b
      use dataio_pub,    only: chdf, nres, set_container_chdf, cwdlen
      use fluidindex,    only: nvar
      use grid,          only: nxb, nyb, nzb, x, y, z, nx, ny, nz
      use hdf5,          only: HID_T, HSIZE_T, HSSIZE_T, H5P_FILE_ACCESS_F, H5F_ACC_TRUNC_F, H5P_DATASET_CREATE_F, H5S_SELECT_SET_F, &
           &                   H5P_DATASET_XFER_F, H5FD_MPIO_INDEPENDENT_F, H5T_NATIVE_DOUBLE, &
           &                   h5open_f, h5close_f, h5fcreate_f, h5fclose_f, &
           &                   h5dcreate_f, h5dwrite_f, h5dclose_f, h5dget_space_f, &
           &                   h5pcreate_f, h5pclose_f, h5pset_fapl_mpio_f, h5pset_chunk_f, h5pset_dxpl_mpio_f, &
           &                   h5screate_simple_f, h5sclose_f, h5sselect_hyperslab_f
      use list_hdf5,     only: problem_write_restart
      use mpisetup,      only: pcoords, pxsize, pysize, pzsize, comm3d, comm, info, ierr, proc, nstep
      use mpi,           only: MPI_CHARACTER
      use problem_pub,   only: problem_name, run_id
      use types,         only: hdf
#ifdef ISO_LOCAL
      use arrays,        only: cs_iso2_arr
#endif /* ISO_LOCAL */

      implicit none

      integer            :: nu
      CHARACTER(len=cwdlen) :: filename  !> HDF File name

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

      if (proc==0) write(filename, '(a,a1,a3,a1,i4.4,a4)') trim(problem_name), '_', run_id, '_', nres, '.res'
      call MPI_Bcast(filename, cwdlen, MPI_CHARACTER, 0, comm, ierr)
      call set_container_chdf(nstep)

      nu = nvar%all

      CALL h5open_f(error)
      CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
      CALL h5pset_fapl_mpio_f(plist_id, comm3d, info, error)

      CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)
      CALL h5pclose_f(plist_id, error)

      if (associated(problem_write_restart)) call problem_write_restart(file_id)

#ifdef ISO_LOCAL
      call write_3darr_to_restart(cs_iso2_arr(:,:,:),file_id,"cs_iso2",nx,ny,nz)
#endif /* ISO_LOCAL */

      rank = 4
      allocate(dimsf(rank),dimsfi(rank),chunk_dims(rank))
      allocate(count(rank),offset(rank),stride(rank),block(rank))
      !----------------------------------------------------------------------------------
      !  WRITE FLUID VARIABLES
      !
      dimsf = [nu, nxb*pxsize, nyb*pysize, nzb*pzsize] ! Dataset dimensions
      dimsfi = dimsf
      chunk_dims = [nu, nxb, nyb, nzb]                    ! Chunks dimensions

      ! Create the data space for the  dataset.
      CALL h5screate_simple_f(rank, dimsf, filespace, error)
      CALL h5screate_simple_f(rank, chunk_dims, memspace, error)

      ! Create chunked dataset.
      CALL h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
      CALL h5pset_chunk_f(plist_id, rank, chunk_dims, error)
      CALL h5dcreate_f(file_id, dname(1), H5T_NATIVE_DOUBLE, filespace, &
                       dset_id, error, plist_id)
      CALL h5sclose_f(filespace, error)
      CALL h5pclose_f(plist_id, error)

      ! Each process defines dataset in memory and writes it to the hyperslab
      ! in the file.
      stride(:) = 1
      count(:)  = 1
      block(:)  = chunk_dims(:)

      offset(1)   = 0
      offset(2:4) = pcoords(1:3)*chunk_dims(2:4)

      ! Select hyperslab in the file.
      CALL h5dget_space_f(dset_id, filespace, error)
      CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error, &
                                  stride, block)
      CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, u(:,is:ie,js:je,ks:ke), dimsfi, error, &
                      file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

      CALL h5sclose_f(filespace, error)
      CALL h5sclose_f(memspace, error)
      CALL h5dclose_f(dset_id, error)
      CALL h5pclose_f(plist_id, error)
      !----------------------------------------------------------------------------------

      !----------------------------------------------------------------------------------
      !  WRITE MAG VARIABLES
      !
      dimsf = [3, nxb*pxsize, nyb*pysize, nzb*pzsize] ! Dataset dimensions
      dimsfi = dimsf
      chunk_dims = [3, nxb, nyb, nzb]                 ! Chunks dimensions

      ! Create the data space for the  dataset.
      CALL h5screate_simple_f(rank, dimsf, filespace, error)
      CALL h5screate_simple_f(rank, chunk_dims, memspace, error)

      ! Create chunked dataset.
      CALL h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
      CALL h5pset_chunk_f(plist_id, rank, chunk_dims, error)
      CALL h5dcreate_f(file_id, dname(2), H5T_NATIVE_DOUBLE, filespace, &
                       dset_id, error, plist_id)
      CALL h5sclose_f(filespace, error)
      CALL h5pclose_f(plist_id, error)

      ! Each process defines dataset in memory and writes it to the hyperslab
      ! in the file.
      stride(:) = 1
      count(:)  = 1
      block(:)  = chunk_dims(:)

      offset(1)   = 0
      offset(2:4) = pcoords(1:3)*chunk_dims(2:4)

      ! Select hyperslab in the file.
      CALL h5dget_space_f(dset_id, filespace, error)
      CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error, &
                                  stride, block)
      CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, b(:,is:ie,js:je,ks:ke), dimsfi, error, &
                      file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

      CALL h5sclose_f(filespace, error)
      CALL h5sclose_f(memspace, error)
      CALL h5dclose_f(dset_id, error)
      CALL h5pclose_f(plist_id, error)
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

      !----------------------------------------------------------------------------------
      !  WRITE X Axis
      !
      dimsf  = [nxb*pxsize] ! Dataset dimensions
      dimsfi = dimsf
      chunk_dims = [nxb]    ! Chunks dimensions

      ! Create the data space for the  dataset.
      CALL h5screate_simple_f(rank, dimsf, filespace, error)
      CALL h5screate_simple_f(rank, chunk_dims, memspace, error)

      ! Create chunked dataset.
      CALL h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
      CALL h5pset_chunk_f(plist_id, rank, chunk_dims, error)
      CALL h5dcreate_f(file_id, "X axis", H5T_NATIVE_DOUBLE, filespace, &
                       dset_id, error, plist_id)
      CALL h5sclose_f(filespace, error)
      CALL h5pclose_f(plist_id, error)

      ! Each process defines dataset in memory and writes it to the hyperslab
      ! in the file.
      stride(:) = 1
      count(:)  = 1
      block(:)  = chunk_dims(:)

      offset(1) = pcoords(1)*chunk_dims(1)

      ! Select hyperslab in the file.
      CALL h5dget_space_f(dset_id, filespace, error)
      CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error, &
                                  stride, block)
      CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, x(is:ie), dimsfi, error, &
                      file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

      CALL h5sclose_f(filespace, error)
      CALL h5sclose_f(memspace, error)
      CALL h5dclose_f(dset_id, error)
      CALL h5pclose_f(plist_id, error)
      !----------------------------------------------------------------------------------

      !----------------------------------------------------------------------------------
      !  WRITE Y Axis
      !
      dimsf  = [nyb*pysize] ! Dataset dimensions
      dimsfi = dimsf
      chunk_dims = [nyb]    ! Chunks dimensions

      ! Create the data space for the  dataset.
      CALL h5screate_simple_f(rank, dimsf, filespace, error)
      CALL h5screate_simple_f(rank, chunk_dims, memspace, error)

      ! Create chunked dataset.
      CALL h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
      CALL h5pset_chunk_f(plist_id, rank, chunk_dims, error)
      CALL h5dcreate_f(file_id, "Y axis", H5T_NATIVE_DOUBLE, filespace, &
                       dset_id, error, plist_id)
      CALL h5sclose_f(filespace, error)
      CALL h5pclose_f(plist_id, error)

      ! Each process defines dataset in memory and writes it to the hyperslab
      ! in the file.
      stride(:) = 1
      count(:)  = 1
      block(:)  = chunk_dims(:)

      offset(1) = pcoords(2)*chunk_dims(1)

      ! Select hyperslab in the file.
      CALL h5dget_space_f(dset_id, filespace, error)
      CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error, &
                                  stride, block)
      CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, y(js:je), dimsfi, error, &
                      file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

      CALL h5sclose_f(filespace, error)
      CALL h5sclose_f(memspace, error)
      CALL h5dclose_f(dset_id, error)
      CALL h5pclose_f(plist_id, error)
      !----------------------------------------------------------------------------------

      !----------------------------------------------------------------------------------
      !  WRITE Z Axis
      !
      dimsf  = [nzb*pzsize] ! Dataset dimensions
      dimsfi = dimsf
      chunk_dims = [nzb]    ! Chunks dimensions

      ! Create the data space for the  dataset.
      CALL h5screate_simple_f(rank, dimsf, filespace, error)
      CALL h5screate_simple_f(rank, chunk_dims, memspace, error)

      ! Create chunked dataset.
      CALL h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
      CALL h5pset_chunk_f(plist_id, rank, chunk_dims, error)
      CALL h5dcreate_f(file_id, "Z axis", H5T_NATIVE_DOUBLE, filespace, &
                       dset_id, error, plist_id)
      CALL h5sclose_f(filespace, error)
      CALL h5pclose_f(plist_id, error)

      ! Each process defines dataset in memory and writes it to the hyperslab
      ! in the file.
      stride(:) = 1
      count(:)  = 1
      block(:)  = chunk_dims(:)

      offset(1) = pcoords(3)*chunk_dims(1)

      ! Select hyperslab in the file.
      CALL h5dget_space_f(dset_id, filespace, error)
      CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error, &
                                  stride, block)
      CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, z(ks:ke), dimsfi, error, &
                      file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

      CALL h5sclose_f(filespace, error)
      CALL h5sclose_f(memspace, error)
      CALL h5dclose_f(dset_id, error)
      !----------------------------------------------------------------------------------
      CALL h5pclose_f(plist_id, error)
      CALL h5fclose_f(file_id, error)
      if (allocated(dimsf))      deallocate(dimsf)
      if (allocated(dimsfi))     deallocate(dimsfi)
      if (allocated(chunk_dims)) deallocate(chunk_dims)
      if (allocated(count))      deallocate(count)
      if (allocated(offset))     deallocate(offset)
      if (allocated(stride))     deallocate(stride)
      if (allocated(block))      deallocate(block)

      call set_common_attributes(filename, chdf, "restart")

      CALL h5close_f(error)

      nres = nres + 1

   end subroutine write_restart_hdf5

   subroutine write_3darr_to_restart(tab,file_id,dname,nx,ny,nz)
      use hdf5,        only: HID_T, HSIZE_T, HSSIZE_T, H5P_DATASET_CREATE_F, H5S_SELECT_SET_F, &
           &                 H5P_DATASET_XFER_F, H5FD_MPIO_INDEPENDENT_F, H5T_NATIVE_DOUBLE, &
           &                 h5dcreate_f, h5dwrite_f, h5dclose_f, h5dget_space_f, &
           &                 h5pcreate_f, h5pset_chunk_f, h5pset_dxpl_mpio_f, &
           &                 h5screate_simple_f, h5sclose_f, h5sselect_hyperslab_f
      use mpisetup,    only: pcoords, pxsize, pysize, pzsize
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
      dimsf = [nx*pxsize, ny*pysize, nz*pzsize] ! Dataset dimensions
      dimsfi = dimsf
      chunk_dims = [nx, ny, nz]                 ! Chunks dimensions

      ! Create the data space for the  dataset.
      CALL h5screate_simple_f(rank, dimsf, filespace, error)
      CALL h5screate_simple_f(rank, chunk_dims, memspace, error)

      ! Create chunked dataset.
      CALL h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
      CALL h5pset_chunk_f(plist_id, rank, chunk_dims, error)
      CALL h5dcreate_f(file_id, trim(dname), H5T_NATIVE_DOUBLE, filespace, &
                       dset_id, error, plist_id)
      CALL h5sclose_f(filespace, error)

      ! Each process defines dataset in memory and writes it to the hyperslab
      ! in the file.
      stride(:) = 1
      count(:)  = 1
      block(:)  = chunk_dims(:)

      offset(1)   = 0
      offset(1:3) = pcoords(1:3)*chunk_dims(1:3)

      ! Select hyperslab in the file.
      CALL h5dget_space_f(dset_id, filespace, error)
      CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error, &
                                  stride, block)
      CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, tab, dimsfi, error, &
                     file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

      CALL h5sclose_f(filespace, error)
      CALL h5sclose_f(memspace, error)
      CALL h5dclose_f(dset_id, error)

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
      use hdf5,         only: HID_T, HSIZE_T, HSSIZE_T, SIZE_T, H5T_NATIVE_DOUBLE, &
                              H5S_SELECT_SET_F, H5FD_MPIO_INDEPENDENT_F, H5P_DATASET_XFER_F, &
                              h5pcreate_f, h5pclose_f, h5screate_simple_f, h5dopen_f, &
                              h5dget_space_f, h5sget_simple_extent_ndims_f, h5dget_create_plist_f, h5pget_chunk_f, &
                              h5sselect_hyperslab_f, h5dread_f, h5sclose_f, h5pset_dxpl_mpio_f, h5dclose_f
      use mpisetup,     only: pcoords, pxsize, pysize, pzsize

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
      integer(HSIZE_T),  DIMENSION(:), allocatable :: dimsf, dimsfi, chunk_dims, old_chunk

      integer               :: error, rank

      !----------------------------------------------------------------------------------
      !  READ TAB
      !
      rank = 3
      allocate(dimsf(rank), dimsfi(rank), chunk_dims(rank), old_chunk(rank))
      allocate(block(rank), offset(rank), count(rank), stride(rank))
      dimsf = [nx*pxsize, ny*pysize, nz*pzsize] ! Dataset dimensions
      dimsfi = dimsf
      chunk_dims = [nx, ny, nz]

      ! Create chunked dataset.
      CALL h5dopen_f(file_id, dname, dset_id, error)

      call h5dget_space_f(dset_id, filespace, error)
      call H5sget_simple_extent_ndims_f (filespace,rank,error)
      call H5dget_create_plist_f (dset_id,plist_id,error)
      call h5pget_chunk_f(plist_id, rank, old_chunk, error)

      ! Each process defines dataset in memory and writes it to the hyperslab
      ! in the file.
      stride(:) = 1
      count(:)  = 1
      block(:)  = chunk_dims(:)

      offset(1:3) = pcoords(1:3)*chunk_dims(1:3)

      ! Select hyperslab in the file.
      CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error, &
                                  stride, block)
      CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
      CALL h5screate_simple_f(rank, chunk_dims, memspace, error)
      CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, p3d, dimsfi, error, &
                     file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

      CALL h5sclose_f(filespace, error)
      CALL h5sclose_f(memspace, error)
      CALL h5pclose_f(plist_id, error)
      CALL h5dclose_f(dset_id, error)
      !----------------------------------------------------------------------------------
      if (allocated(dimsf))      deallocate(dimsf)
      if (allocated(dimsfi))     deallocate(dimsfi)
      if (allocated(chunk_dims)) deallocate(chunk_dims)
      if (allocated(old_chunk))  deallocate(old_chunk)
      if (allocated(count))      deallocate(count)
      if (allocated(offset))     deallocate(offset)
      if (allocated(stride))     deallocate(stride)
      if (allocated(block))      deallocate(block)

   end subroutine read_3darr_from_restart

   subroutine read_restart_hdf5(chdf)

      use arrays,        only: u, b
      use dataio_pub,    only: cwdlen, msg, printio, die
      use fluidindex,    only: nvar
      use func,          only: fix_string
      use grid,          only: nx, ny, nz, x, y, z, nxb, nyb, nzb
      use hdf5,          only: HID_T, HSIZE_T, HSSIZE_T, SIZE_T, H5P_FILE_ACCESS_F, H5T_NATIVE_DOUBLE, &
                               H5S_SELECT_SET_F, H5F_ACC_RDONLY_F, H5FD_MPIO_INDEPENDENT_F, H5P_DATASET_XFER_F, &
                               h5open_f, h5pcreate_f, h5pset_fapl_mpio_f, h5fopen_f, h5pclose_f, h5dopen_f, &
                               h5dget_space_f, h5sget_simple_extent_ndims_f, h5dget_create_plist_f, h5pget_chunk_f, &
                               h5sselect_hyperslab_f, h5dread_f, h5sclose_f, h5pset_dxpl_mpio_f, h5dclose_f, &
                               h5screate_simple_f, h5fclose_f, h5close_f
      use h5lt,          only: h5ltget_attribute_double_f, h5ltget_attribute_int_f, h5ltget_attribute_string_f
      use list_hdf5,     only: problem_read_restart
      use mpisetup,      only: comm, ierr, pcoords, pxsize, pysize, pzsize, magic_mass, proc, t, info, comm3d, dt, cbuff_len
      use mpi,           only: MPI_CHARACTER, MPI_INTEGER, MPI_DOUBLE_PRECISION
      use problem_pub,   only: problem_name, run_id
      use types,         only: hdf, domlen, idlen
#ifdef ISO_LOCAL
      use arrays,        only: cs_iso2_arr
#endif /* ISO_LOCAL */

      IMPLICIT NONE

      type(hdf)             :: chdf
      integer               :: nu
      CHARACTER(LEN=cwdlen) :: filename  ! File name

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
      integer(HSIZE_T),  DIMENSION(:), allocatable :: dimsf, dimsfi, chunk_dims, old_chunk

      integer               :: error, rank
      logical               :: file_exist

      real, dimension(1)    :: rbuf
      integer, dimension(1) :: ibuf

      real, dimension(:,:,:,:), pointer :: p4d
#ifdef ISO_LOCAL
      real, dimension(:,:,:), pointer :: p3d
#endif /* ISO_LOCAL */

      nu = nvar%all

      if (proc==0) then
         write(filename,'(a,a1,a3,a1,i4.4,a4)') trim(problem_name),'_', run_id,'_',chdf%nres,'.res'
         write(msg, '(2a)') 'Reading restart file: ',trim(filename)
         call printio(msg)
      endif
      call MPI_Bcast(filename, cwdlen, MPI_CHARACTER, 0, comm, ierr)

      inquire(file = filename, exist = file_exist)
      if (file_exist .eqv. .false.) then
         write(msg,'(3a)') '[dataio_hdf5:read_restart_hdf5]: Restart file: ',trim(filename),' does not exist'
         call die(msg)
      endif

      CALL h5open_f(error)
      CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
      CALL h5pset_fapl_mpio_f(plist_id, comm3d, info, error)

      CALL h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, error, access_prp = plist_id)
      CALL h5pclose_f(plist_id, error)

      if (associated(problem_read_restart)) call problem_read_restart(file_id)

#ifdef ISO_LOCAL
      if (.not.associated(p3d)) p3d => cs_iso2_arr(:,:,:)
      call read_3darr_from_restart(file_id,"cs_iso2",p3d,nx,ny,nz)
      if (associated(p3d)) nullify(p3d)
#endif /* ISO_LOCAL */
      !----------------------------------------------------------------------------------
      !  READ FLUID VARIABLES
      !
      rank = 4
      allocate(dimsf(rank), dimsfi(rank), chunk_dims(rank), old_chunk(rank))
      allocate(block(rank), offset(rank), count(rank), stride(rank))
      dimsf = [nu, nxb*pxsize, nyb*pysize, nzb*pzsize] ! Dataset dimensions
      dimsfi = dimsf
      chunk_dims = [nu, nxb, nyb, nzb]
      if (.not.associated(p4d)) p4d => u(:,is:ie,js:je,ks:ke)

      ! Create chunked dataset.
      CALL h5dopen_f(file_id, dname(1), dset_id, error)

      call h5dget_space_f(dset_id, filespace, error)
      call H5sget_simple_extent_ndims_f (filespace,rank,error)
      call H5dget_create_plist_f (dset_id,plist_id,error)
      call h5pget_chunk_f(plist_id, rank, old_chunk, error)
!      if ( any(old_chunk /= chunk_dims) ) call die("[dataio_hdf5:read_restart_hdf5]: Chunk_dims differs from chunk in restart file")

      ! Each process defines dataset in memory and writes it to the hyperslab
      ! in the file.
      stride(:) = 1
      count(:)  = 1
      block(:)  = chunk_dims(:)

      offset(1)   = 0
      offset(2:4) = pcoords(1:3)*chunk_dims(2:4)

      ! Select hyperslab in the file.
      CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error, &
                                  stride, block)
      CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
      CALL h5screate_simple_f(rank, chunk_dims, memspace, error)
      CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, p4d, dimsfi, error, &
                     file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

      CALL h5sclose_f(filespace, error)
      CALL h5sclose_f(memspace, error)
      CALL h5pclose_f(plist_id, error)
      CALL h5dclose_f(dset_id, error)
      if (associated(p4d)) nullify(p4d)
      !----------------------------------------------------------------------------------

      !----------------------------------------------------------------------------------
      !  READ MAG VARIABLES
      !
      dimsf = [3, nxb*pxsize, nyb*pysize, nzb*pzsize] ! Dataset dimensions
      dimsfi = dimsf
      chunk_dims = [3, nxb, nyb, nzb]
      if (.not.associated(p4d)) p4d => b(:,is:ie,js:je,ks:ke)

      ! Create chunked dataset.
      CALL h5dopen_f(file_id, dname(2), dset_id, error)

      call h5dget_space_f(dset_id, filespace, error)
      call H5Sget_simple_extent_ndims_f (filespace,rank,error)
      call H5Dget_create_plist_f (dset_id,plist_id,error)
      call h5pget_chunk_f(plist_id, rank, old_chunk, error)
!      if ( any(old_chunk /= chunk_dims) ) call die("[dataio_hdf5:read_restart_hdf5]: Chunk_dims differs from chunk in restart file")

      ! Each process defines dataset in memory and writes it to the hyperslab
      ! in the file.
      stride(:) = 1
      count(:)  = 1
      block(:)  = chunk_dims(:)

      offset(1)   = 0
      offset(2:4) = pcoords(1:3)*chunk_dims(2:4)

      ! Select hyperslab in the file.
      CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error, &
                                  stride, block)
      CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
      CALL h5screate_simple_f(rank, chunk_dims, memspace, error)
      CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, p4d, dimsfi, error, &
                     file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

      CALL h5sclose_f(filespace, error)
      CALL h5sclose_f(memspace, error)
      CALL h5pclose_f(plist_id, error)
      CALL h5dclose_f(dset_id, error)
      if (associated(p4d)) nullify(p4d)
      !----------------------------------------------------------------------------------
      CALL h5fclose_f(file_id, error)
      if (allocated(dimsf))      deallocate(dimsf)
      if (allocated(dimsfi))     deallocate(dimsfi)
      if (allocated(chunk_dims)) deallocate(chunk_dims)
      if (allocated(old_chunk))  deallocate(old_chunk)
      if (allocated(count))      deallocate(count)
      if (allocated(offset))     deallocate(offset)
      if (allocated(stride))     deallocate(stride)
      if (allocated(block))      deallocate(block)
      if (proc == 0) then
         CALL h5fopen_f (filename, H5F_ACC_RDONLY_F, file_id, error)
         bufsize = 1
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
         call h5ltget_attribute_int_f(file_id,"/","ntsl", ibuf,error)
         chdf%ntsl = ibuf(1)
         call h5ltget_attribute_int_f(file_id,"/","nlog", ibuf,error)
         chdf%nlog = ibuf(1)
         call h5ltget_attribute_int_f(file_id,"/","step_res", ibuf,error)
         chdf%step_res = ibuf(1)
         call h5ltget_attribute_int_f(file_id,"/","step_hdf", ibuf,error)
         chdf%step_hdf = ibuf(1)
         call h5ltget_attribute_double_f(file_id,"/","last_hdf_time", rbuf,error)
         chdf%last_hdf_time = rbuf(1)

         call h5ltget_attribute_string_f(file_id,"/","problem_name", problem_name,error)
         call h5ltget_attribute_string_f(file_id,"/","domain", chdf%domain,error)
         call h5ltget_attribute_string_f(file_id,"/","run_id", chdf%new_id,error)

         problem_name = fix_string(problem_name)   ! BEWARE: >=HDF5-1.8.4 has weird issues with strings
         chdf%new_id  = fix_string(chdf%new_id)    !   this bit hacks it around
         chdf%domain  = fix_string(chdf%domain)

         CALL h5fclose_f(file_id, error)

         write(msg,'(2a)') 'Done reading restart file: ',trim(filename)
         call printio(msg)
      endif

      call MPI_Bcast(chdf%nstep,    1, MPI_INTEGER, 0, comm3d, ierr)
      call MPI_Bcast(chdf%nres,     1, MPI_INTEGER, 0, comm3d, ierr)
      call MPI_Bcast(chdf%nhdf,     1, MPI_INTEGER, 0, comm3d, ierr)
      call MPI_Bcast(chdf%ntsl,     1, MPI_INTEGER, 0, comm3d, ierr)
      call MPI_Bcast(chdf%nlog,     1, MPI_INTEGER, 0, comm3d, ierr)
      call MPI_Bcast(chdf%step_res, 1, MPI_INTEGER, 0, comm3d, ierr)
      call MPI_Bcast(chdf%step_hdf, 1, MPI_INTEGER, 0, comm3d, ierr)

      call MPI_Bcast(chdf%last_hdf_time, 1, MPI_DOUBLE_PRECISION, 0, comm3d, ierr)
      call MPI_Bcast(t,                  1, MPI_DOUBLE_PRECISION, 0, comm3d, ierr)
      call MPI_Bcast(dt,                 1, MPI_DOUBLE_PRECISION, 0, comm3d, ierr)

      CALL MPI_Bcast(problem_name, cbuff_len, MPI_CHARACTER, 0, comm3d, ierr)
      CALL MPI_Bcast(chdf%domain,  domlen,    MPI_CHARACTER, 0, comm3d, ierr)
      CALL MPI_Bcast(chdf%new_id,  idlen,     MPI_CHARACTER, 0, comm3d, ierr)
      CALL h5close_f(error)

   end subroutine read_restart_hdf5
!
! ------------------------------------------------------------------------------------
!
   subroutine write_hdf5(chdf)
      use dataio_pub,    only: cwdlen, msg, die, user_vars_hdf5, nhdf, varlen
      use grid,          only: nxb, nyb, nzb
      use hdf5,          only: HID_T, H5F_ACC_TRUNC_F, H5P_FILE_ACCESS_F, H5P_DEFAULT_F, &
           &                   h5open_f, h5close_f, h5fcreate_f, h5fclose_f, h5pcreate_f, h5pclose_f, h5pset_fapl_mpio_f
      use mpisetup,      only: comm3d, ierr, info
      use problem_pub,   only: problem_name, run_id
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
      CHARACTER(LEN=cwdlen)   :: fname

      real(kind=4), allocatable :: data (:,:,:)  ! Data to write

      ! Initialize HDF5 library and Fortran interfaces.
      !
      write(dd,'(i4.4)') chdf%nhdf
      write(fname, '(6a)') trim(problem_name),"_",trim(run_id),"_",dd,".h5"

      CALL h5open_f(error)
      !
      ! Setup file access property list with parallel I/O access.
      !
      CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
      CALL h5pset_fapl_mpio_f(plist_id, comm3d, info, error)
      !
      ! Create the file collectively.
      !
      CALL h5fcreate_f(trim(fname), H5F_ACC_TRUNC_F, file_id, error, creation_prp = H5P_DEFAULT_F, access_prp = plist_id)
      CALL h5pclose_f(plist_id, error)
      if (.not.allocated(data)) allocate(data(nxb,nyb,nzb))
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
      CALL h5fclose_f(file_id, error)

      call set_common_attributes(fname, chdf, "output")

      call MPI_Barrier(comm3d,ierr)
      CALL h5close_f(error)

      nhdf = nhdf + 1

   end subroutine write_hdf5

   subroutine write_arr(data,dsetname,file_id)
      use dataio_pub,    only: varlen
      use grid,          only: nxb, nyb, nzb
      use hdf5,          only: HID_T, HSIZE_T, HSSIZE_T, h5screate_simple_f, h5pcreate_f, h5pset_chunk_f, &
                           h5sclose_f, h5pset_dxpl_mpio_f, h5dwrite_f, h5dclose_f, H5P_DATASET_XFER_F, H5P_DATASET_CREATE_F, &
                           H5T_NATIVE_REAL, H5S_SELECT_SET_F, H5FD_MPIO_INDEPENDENT_F, h5dcreate_f, h5dget_space_f, &
                           h5sselect_hyperslab_f
      use mpisetup,      only: pcoords, pxsize, pysize, pzsize

      implicit none
      real(kind=4), dimension(:,:,:) :: data

      integer, parameter :: rank = 3       ! Dataset rank

      CHARACTER(LEN=varlen) :: dsetname       ! Dataset name

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

      chunk_dims = [nxb,nyb,nzb]                     ! Chunks dimensions
      dimsf  = chunk_dims * [pxsize, pysize, pzsize]  ! Dataset dimensions
      dimsfi = dimsf
      !
      ! Create the data space for the  dataset.
      !
      CALL h5screate_simple_f(rank, dimsf, filespace, error)
      CALL h5screate_simple_f(rank, chunk_dims, memspace, error)

      !
      ! Create chunked dataset.
      !
      CALL h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
      CALL h5pset_chunk_f(plist_id, rank, chunk_dims, error)
      CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_REAL, filespace, &
                       dset_id, error, plist_id)
      CALL h5sclose_f(filespace, error)

      !
      ! Each process defines dataset in memory and writes it to the hyperslab
      ! in the file.
      !
      stride(:) = 1
      count(:) =  1
      block(:) = chunk_dims(:)

      offset(:) = pcoords(:)*chunk_dims(:)
      !
      ! Select hyperslab in the file.
      !
      CALL h5dget_space_f(dset_id, filespace, error)
      CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error, &
                                  stride, block)
      !
      ! Create property list for collective dataset write
      !
      CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)

      !
      ! Write the dataset collectively.
      !
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, data, dimsfi, error, &
                      file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

      !
      ! Close dataspaces.
      !
      CALL h5sclose_f(filespace, error)
      CALL h5sclose_f(memspace, error)
      !
      ! Close the dataset.
      !
      CALL h5dclose_f(dset_id, error)

   end subroutine write_arr

! This routine writes all attributes that are common to restart and output files.
! Other common elements may also be moved here.

   subroutine set_common_attributes(filename, chdf, stype)

      use dataio_pub,    only: msg, printio
      use grid,          only: nxb, nyb, nzb, nb, xmin, xmax, ymin, ymax, zmin, zmax
      use hdf5,          only: HID_T, SIZE_T, H5F_ACC_RDWR_F, h5fopen_f, h5fclose_f, h5gcreate_f, h5gclose_f
      use h5lt,          only: h5ltset_attribute_double_f, h5ltset_attribute_int_f, h5ltmake_dataset_string_f, h5ltset_attribute_string_f
      use list_hdf5,     only: additional_attrs
      use mpisetup,      only: proc, t, dt, psize, cbuff_len, pxsize, pysize, pzsize, local_magic_mass, comm, ierr, magic_mass
      use mpi,           only: MPI_DOUBLE_PRECISION, MPI_SUM
      use problem_pub,   only: problem_name, run_id
      use types,         only: hdf
      use version,       only: env, nenv

      implicit none

      character(len=*), intent(in) :: filename  !> HDF File name
      character(len=*), intent(in) :: stype     !> "output" or "restart"
      type(hdf), intent(in)        :: chdf

      integer(HID_T) :: file_id       !> File identifier
      integer(HID_T) :: gr_id         !> Group identifier
      character(len=cbuff_len) :: dset_name
      character(len=maxparfilelen) :: trim_env
      integer            :: fe, i
      integer(SIZE_T) :: bufsize
      integer :: error
      real    :: magic_mass0
      integer, parameter          :: buf_len = 50
      integer, dimension(buf_len) :: ibuffer
      real,    dimension(buf_len) :: rbuffer
      character(len=cbuff_len), dimension(buf_len) :: ibuffer_name = ''
      character(len=cbuff_len), dimension(buf_len) :: rbuffer_name = ''

      call MPI_REDUCE(local_magic_mass, magic_mass0, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
      magic_mass       = magic_mass + magic_mass0
      local_magic_mass = 0.0

      if (proc == 0) then

         call h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, error)

         rbuffer(1)  = t                        ; rbuffer_name(1)  = "time"
         rbuffer(2)  = dt                       ; rbuffer_name(2)  = "timestep"
         rbuffer(3)  = chdf%last_hdf_time       ; rbuffer_name(3)  = "last_hdf_time"
         rbuffer(4)  = xmin                     ; rbuffer_name(4)  = "xmin"
         rbuffer(5)  = xmax                     ; rbuffer_name(5)  = "xmax"
         rbuffer(6)  = ymin                     ; rbuffer_name(6)  = "ymin"
         rbuffer(7)  = ymax                     ; rbuffer_name(7)  = "ymax"
         rbuffer(8)  = zmin                     ; rbuffer_name(8)  = "zmin"
         rbuffer(9)  = zmax                     ; rbuffer_name(9)  = "zmax"
         rbuffer(10) = piernik_hdf5_version     ; rbuffer_name(10) = "piernik"
         rbuffer(11) = magic_mass               ; rbuffer_name(11) = "magic_mass"

         ibuffer(1)  = chdf%nstep               ; ibuffer_name(1)  = "nstep"
         ibuffer(2)  = chdf%nres+1              ; ibuffer_name(2)  = "nres"
         ibuffer(3)  = chdf%nhdf                ; ibuffer_name(3)  = "nhdf"
         ibuffer(4)  = chdf%ntsl                ; ibuffer_name(4)  = "ntsl"
         ibuffer(5)  = chdf%nlog                ; ibuffer_name(5)  = "nlog"
         ibuffer(6)  = chdf%nstep               ; ibuffer_name(6)  = "step_res"
         ibuffer(7)  = chdf%step_hdf            ; ibuffer_name(7)  = "step_hdf"
         ibuffer(8)  = nxb*pxsize               ; ibuffer_name(8)  = "nxd"
         ibuffer(9)  = nyb*pysize               ; ibuffer_name(9)  = "nyd"
         ibuffer(10) = nzb*pzsize               ; ibuffer_name(10) = "nzd"
         ibuffer(11) = nxb                      ; ibuffer_name(11) = "nxb"
         ibuffer(12) = nyb                      ; ibuffer_name(12) = "nyb"
         ibuffer(13) = nzb                      ; ibuffer_name(13) = "nzb"
         ibuffer(14) = nb                       ; ibuffer_name(14) = "nb"

         !BEWARE: A memory leak was detected here. h5lt calls use HD5f2cstring and probably sometimes don't free the allocated buffer

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

         call H5Gcreate_f(file_id,"file_versions",gr_id,error)
         do i = 1, nenv
            write(dset_name,'(A4,I3.3)') "env_",i
            trim_env=trim(env(i))
            if (len_trim(trim_env) < len(trim_env)) trim_env(len_trim(trim_env)+1:len_trim(trim_env)+1) = achar(0)
            ! I don't understand why passing trim(env(i)) or trim(trim_env) to h5ltmake_dataset_string_f sometimes results in improper detection of the string length.
            call h5ltmake_dataset_string_f (gr_id, dset_name, trim_env, error)
         enddo
         call H5Gclose_f(gr_id, error)

         call H5Gcreate_f(file_id,"problem.par",gr_id,error)
         do i = 1, parfilelines
            write(dset_name,'(A8,I4.4)') "parfile_",i
            call h5ltmake_dataset_string_f (gr_id, trim(dset_name), trim(parfile(i)), error)
         enddo
         call H5Gclose_f(gr_id, error)

         bufsize = 3
         call h5ltset_attribute_int_f(file_id, "/", "psize", psize, bufsize, error)

         fe = len_trim(problem_name)
         call h5ltset_attribute_string_f(file_id, "/", "problem_name", problem_name(1:fe), error)
         fe = len_trim(chdf%domain)
         call h5ltset_attribute_string_f(file_id, "/", "domain", chdf%domain(1:fe), error)
         fe = len_trim(run_id)
         call h5ltset_attribute_string_f(file_id, "/", "run_id", run_id(1:fe), error)

         call additional_attrs(file_id)

         call h5fclose_f(file_id, error)

         write(msg,'(4a)') 'Writing ',stype,' file: ',trim(filename)
         call printio(msg)
      endif

   end subroutine set_common_attributes

end module dataio_hdf5
