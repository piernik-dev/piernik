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
!    Initial implemetation of PIERNIK code was based on TVD split MHD code by
!    Ue-Li Pen
!        see: Pen, Arras & Wong (2003) for algorithm and
!             http://www.cita.utoronto.ca/~pen/MHD
!             for original source code "mhd.f90"
!
!    For full list of developers see $PIERNIK_HOME/license/pdt.txt
!

#define RNG is:ie, js:je, ks:ke
#include "piernik.def"

!>
!! \brief (KK) Module that contains I/O routines using HDF5 library
!!
!! Modules contains routines for creating HDF5 output such as
!! plots, snapshots, restart files.
!!
!<
module dataio_hdf5
   use list_hdf5, only : S_LEN

   implicit none

   private
   public :: init_hdf5, read_restart_hdf5, cleanup_hdf5, write_hdf5, write_restart_hdf5, write_plot

   character(LEN=10), dimension(3) :: dname = (/"fluid     ","mag       ","dinit     "/)  !< dataset names for restart files
   character(len=S_LEN), allocatable, dimension(:) :: hdf_vars  !< dataset names for hdf files
   integer :: nhdf_vars !< number of quantities ploted to hdf files
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

   contains

!>
!! \brief Procedure initializing HDF5 module
!!
!<

   subroutine init_hdf5(vars,tix,tiy,tiz,tdt_plt)
      use fluidindex, only : iarr_all_crs, iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz
      use grid, only : nx,ny,nz,nxd,nyd,nzd,nb
      implicit none
      character(len=4), dimension(:), intent(in) :: vars  !< quantities to be plotted, see dataio::vars
      integer,intent(in) :: tix     !< local copy of dataio::ix
      integer,intent(in) :: tiy     !< local copy of dataio::iy
      integer,intent(in) :: tiz     !< local copy of dataio::iz
      real,intent(in)    :: tdt_plt !< local copy of dataio::dt_plt
      integer :: nvars,i,j,k
      character(len=4) :: aux

      dname(1) = "fluid"
      dname(2) = "mag"
      dname(3) = "dinit"

      ix = tix; iy = tiy; iz = tiz; dt_plt = tdt_plt

      if(nxd == 1) then
         is = 1; ie = 1
      else
         is = nb+1; ie = nx-nb
      endif

      if(nyd == 1) then
         js = 1; je = 1
      else
         js = nb+1; je = ny-nb
      endif

      if(nzd == 1) then
         ks = 1; ke = 1
      else
         ks = nb+1; ke = nz-nb
      endif

      nvars = 1
      do while( LEN(TRIM(vars(nvars))) > 1)
        nvars = nvars + 1
      enddo
      nvars = nvars - 1

      nhdf_vars = 0
      do i = 1, nvars
         select case(vars(i))
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
            case ('magx')
               nhdf_vars = nhdf_vars + 1
            case ('magy')
               nhdf_vars = nhdf_vars + 1
            case ('magz')
               nhdf_vars = nhdf_vars + 1
#ifdef COSM_RAYS
            case ('encr')
               nhdf_vars = nhdf_vars + SIZE(iarr_all_crs,1)
#endif /* COSM_RAYS */
         end select
      enddo
      allocate(hdf_vars(nhdf_vars)); j = 1
      do i = 1, nvars
         select case(vars(i))
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
               do k = 1, size(iarr_all_crs,1)
                  write(aux,'(A3,I1)') 'ecr',k
                  hdf_vars(j) = aux ; j = j + 1
               enddo
#endif /* COSM_RAYS */
         end select
      enddo

   end subroutine init_hdf5

!>
!! \brief Procedure finalizing HDF5 module
!<
   subroutine cleanup_hdf5
      implicit none

      if(allocated(hdf_vars)) deallocate(hdf_vars)

   end subroutine cleanup_hdf5

!<
!! \brief Routine calculating quantities for plot files
!<
   subroutine common_plt_hdf5(var,ij,xn,tab,ierrh)
      use grid,   only : nb, nyb, nzb, nxb
      use arrays, only : u, b
      use fluidindex, only : ind

      implicit none
      character(LEN=4)     :: var !< quantity to be plotted
      character(LEN=2)     :: ij  !< plane of plot
      integer              :: xn  !< no. of cell at which we are slicing the local block
      integer              :: ierrh !< error handling
      real, dimension(:,:) :: tab !< array containing given quantity

      ierrh = 0
      select case(var)
         case ("dend")
            if(ij=="yz") tab(:,:) = u(ind%dnd,xn,nb+1:nyb+nb,nb+1:nzb+nb)
            if(ij=="xz") tab(:,:) = u(ind%dnd,nb+1:nxb+nb,xn,nb+1:nzb+nb)
            if(ij=="xy") tab(:,:) = u(ind%dnd,nb+1:nxb+nb,nb+1:nyb+nb,xn)
         case ("denn")
            if(ij=="yz") tab(:,:) = u(ind%dnn,xn,nb+1:nyb+nb,nb+1:nzb+nb)
            if(ij=="xz") tab(:,:) = u(ind%dnn,nb+1:nxb+nb,xn,nb+1:nzb+nb)
            if(ij=="xy") tab(:,:) = u(ind%dnn,nb+1:nxb+nb,nb+1:nyb+nb,xn)
         case ("deni")
            if(ij=="yz") tab(:,:) = u(ind%dni,xn,nb+1:nyb+nb,nb+1:nzb+nb)
            if(ij=="xz") tab(:,:) = u(ind%dni,nb+1:nxb+nb,xn,nb+1:nzb+nb)
            if(ij=="xy") tab(:,:) = u(ind%dni,nb+1:nxb+nb,nb+1:nyb+nb,xn)
         case ("vlxd")
            if(ij=="yz") tab(:,:) = u(ind%mxd,xn,nb+1:nyb+nb,nb+1:nzb+nb) / &
                           u(ind%dnd,xn,nb+1:nyb+nb,nb+1:nzb+nb)
            if(ij=="xz") tab(:,:) = u(ind%mxd,nb+1:nxb+nb,xn,nb+1:nzb+nb) / &
                           u(ind%dnd,nb+1:nxb+nb,xn,nb+1:nzb+nb)
            if(ij=="xy") tab(:,:) = u(ind%mxd,nb+1:nxb+nb,nb+1:nyb+nb,xn) / &
                           u(ind%dnd,nb+1:nxb+nb,nb+1:nyb+nb,xn)
         case ("vlxn")
            if(ij=="yz") tab(:,:) = u(ind%mxn,xn,nb+1:nyb+nb,nb+1:nzb+nb) / &
                           u(ind%dnn,xn,nb+1:nyb+nb,nb+1:nzb+nb)
            if(ij=="xz") tab(:,:) = u(ind%mxn,nb+1:nxb+nb,xn,nb+1:nzb+nb) / &
                           u(ind%dnn,nb+1:nxb+nb,xn,nb+1:nzb+nb)
            if(ij=="xy") tab(:,:) = u(ind%mxn,nb+1:nxb+nb,nb+1:nyb+nb,xn) / &
                           u(ind%dnn,nb+1:nxb+nb,nb+1:nyb+nb,xn)
         case ("vlxi")
            if(ij=="yz") tab(:,:) = u(ind%mxi,xn,nb+1:nyb+nb,nb+1:nzb+nb) / &
                           u(ind%dni,xn,nb+1:nyb+nb,nb+1:nzb+nb)
            if(ij=="xz") tab(:,:) = u(ind%mxi,nb+1:nxb+nb,xn,nb+1:nzb+nb) / &
                           u(ind%dni,nb+1:nxb+nb,xn,nb+1:nzb+nb)
            if(ij=="xy") tab(:,:) = u(ind%mxi,nb+1:nxb+nb,nb+1:nyb+nb,xn) / &
                           u(ind%dni,nb+1:nxb+nb,nb+1:nyb+nb,xn)
         case ("vlyd")
            if(ij=="yz") tab(:,:) = u(ind%myd,xn,nb+1:nyb+nb,nb+1:nzb+nb) / &
                           u(ind%dnd,xn,nb+1:nyb+nb,nb+1:nzb+nb)
            if(ij=="xz") tab(:,:) = u(ind%myd,nb+1:nxb+nb,xn,nb+1:nzb+nb) / &
                           u(ind%dnd,nb+1:nxb+nb,xn,nb+1:nzb+nb)
            if(ij=="xy") tab(:,:) = u(ind%myd,nb+1:nxb+nb,nb+1:nyb+nb,xn) / &
                           u(ind%dnd,nb+1:nxb+nb,nb+1:nyb+nb,xn)
         case ("vlyn")
            if(ij=="yz") tab(:,:) = u(ind%myn,xn,nb+1:nyb+nb,nb+1:nzb+nb) / &
                           u(ind%dnn,xn,nb+1:nyb+nb,nb+1:nzb+nb)
            if(ij=="xz") tab(:,:) = u(ind%myn,nb+1:nxb+nb,xn,nb+1:nzb+nb) / &
                           u(ind%dnn,nb+1:nxb+nb,xn,nb+1:nzb+nb)
            if(ij=="xy") tab(:,:) = u(ind%myn,nb+1:nxb+nb,nb+1:nyb+nb,xn) / &
                           u(ind%dnn,nb+1:nxb+nb,nb+1:nyb+nb,xn)
         case ("vlyi")
            if(ij=="yz") tab(:,:) = u(ind%myi,xn,nb+1:nyb+nb,nb+1:nzb+nb) / &
                           u(ind%dni,xn,nb+1:nyb+nb,nb+1:nzb+nb)
            if(ij=="xz") tab(:,:) = u(ind%myi,nb+1:nxb+nb,xn,nb+1:nzb+nb) / &
                           u(ind%dni,nb+1:nxb+nb,xn,nb+1:nzb+nb)
            if(ij=="xy") tab(:,:) = u(ind%myi,nb+1:nxb+nb,nb+1:nyb+nb,xn) / &
                           u(ind%dni,nb+1:nxb+nb,nb+1:nyb+nb,xn)
         case ("vlzd")
            if(ij=="yz") tab(:,:) = u(ind%mzd,xn,nb+1:nyb+nb,nb+1:nzb+nb) / &
                           u(ind%dnd,xn,nb+1:nyb+nb,nb+1:nzb+nb)
            if(ij=="xz") tab(:,:) = u(ind%mzd,nb+1:nxb+nb,xn,nb+1:nzb+nb) / &
                           u(ind%dnd,nb+1:nxb+nb,xn,nb+1:nzb+nb)
            if(ij=="xy") tab(:,:) = u(ind%mzd,nb+1:nxb+nb,nb+1:nyb+nb,xn) / &
                           u(ind%dnd,nb+1:nxb+nb,nb+1:nyb+nb,xn)
         case ("vlzn")
            if(ij=="yz") tab(:,:) = u(ind%mzn,xn,nb+1:nyb+nb,nb+1:nzb+nb) / &
                           u(ind%dnn,xn,nb+1:nyb+nb,nb+1:nzb+nb)
            if(ij=="xz") tab(:,:) = u(ind%mzn,nb+1:nxb+nb,xn,nb+1:nzb+nb) / &
                           u(ind%dnn,nb+1:nxb+nb,xn,nb+1:nzb+nb)
            if(ij=="xy") tab(:,:) = u(ind%mzn,nb+1:nxb+nb,nb+1:nyb+nb,xn) / &
                           u(ind%dnn,nb+1:nxb+nb,nb+1:nyb+nb,xn)
         case ("vlzi")
            if(ij=="yz") tab(:,:) = u(ind%mzi,xn,nb+1:nyb+nb,nb+1:nzb+nb) / &
                           u(ind%dni,xn,nb+1:nyb+nb,nb+1:nzb+nb)
            if(ij=="xz") tab(:,:) = u(ind%mzi,nb+1:nxb+nb,xn,nb+1:nzb+nb) / &
                           u(ind%dni,nb+1:nxb+nb,xn,nb+1:nzb+nb)
            if(ij=="xy") tab(:,:) = u(ind%mzi,nb+1:nxb+nb,nb+1:nyb+nb,xn) / &
                           u(ind%dni,nb+1:nxb+nb,nb+1:nyb+nb,xn)
         case ("enen")
#ifndef ISO
            if(ij=="yz") tab(:,:) = u(ind%enn,xn,nb+1:nyb+nb,nb+1:nzb+nb)
            if(ij=="xz") tab(:,:) = u(ind%enn,nb+1:nxb+nb,xn,nb+1:nzb+nb)
            if(ij=="xy") tab(:,:) = u(ind%enn,nb+1:nxb+nb,nb+1:nyb+nb,xn)
#else
            if(ij=="yz") tab(:,:) = 0.5 * (                     &
                          u(ind%mxn,xn,nb+1:nyb+nb,nb+1:nzb+nb)**2 &
                        + u(ind%myn,xn,nb+1:nyb+nb,nb+1:nzb+nb)**2 &
                        + u(ind%mzn,xn,nb+1:nyb+nb,nb+1:nzb+nb)**2) / &
                             u(ind%dnn,xn,nb+1:nyb+nb,nb+1:nzb+nb)
            if(ij=="xz") tab(:,:) = 0.5 * (                     &
                          u(ind%mxn,nb+1:nxb+nb,xn,nb+1:nzb+nb)**2 &
                         +u(ind%myn,nb+1:nxb+nb,xn,nb+1:nzb+nb)**2 &
                         +u(ind%mzn,nb+1:nxb+nb,xn,nb+1:nzb+nb)**2) / &
                             u(ind%dnn,nb+1:nxb+nb,xn,nb+1:nzb+nb)
            if(ij=="xy") tab(:,:) = 0.5 * (                     &
                          u(ind%mxn,nb+1:nxb+nb,nb+1:nyb+nb,xn)**2 &
                         +u(ind%myn,nb+1:nxb+nb,nb+1:nyb+nb,xn)**2 &
                         +u(ind%mzn,nb+1:nxb+nb,nb+1:nyb+nb,xn)**2) / &
                             u(ind%dnn,nb+1:nxb+nb,nb+1:nyb+nb,xn)
#endif /* ISO */
         case ("enei")
#ifndef ISO
            if(ij=="yz") tab(:,:) = u(ind%eni,xn,nb+1:nyb+nb,nb+1:nzb+nb)
            if(ij=="xz") tab(:,:) = u(ind%eni,nb+1:nxb+nb,xn,nb+1:nzb+nb)
            if(ij=="xy") tab(:,:) = u(ind%eni,nb+1:nxb+nb,nb+1:nyb+nb,xn)
#else
            if(ij=="yz") tab(:,:) = 0.5 * (                     &
                          u(ind%mxi,xn,nb+1:nyb+nb,nb+1:nzb+nb)**2 &
                        + u(ind%myi,xn,nb+1:nyb+nb,nb+1:nzb+nb)**2 &
                        + u(ind%mzi,xn,nb+1:nyb+nb,nb+1:nzb+nb)**2) / &
                             u(ind%dni,xn,nb+1:nyb+nb,nb+1:nzb+nb)
            if(ij=="xz") tab(:,:) = 0.5 * (                     &
                          u(ind%mxi,nb+1:nxb+nb,xn,nb+1:nzb+nb)**2 &
                         +u(ind%myi,nb+1:nxb+nb,xn,nb+1:nzb+nb)**2 &
                         +u(ind%mzi,nb+1:nxb+nb,xn,nb+1:nzb+nb)**2) / &
                             u(ind%dni,nb+1:nxb+nb,xn,nb+1:nzb+nb)
            if(ij=="xy") tab(:,:) = 0.5 * (                     &
                          u(ind%mxi,nb+1:nxb+nb,nb+1:nyb+nb,xn)**2 &
                         +u(ind%myi,nb+1:nxb+nb,nb+1:nyb+nb,xn)**2 &
                         +u(ind%mzi,nb+1:nxb+nb,nb+1:nyb+nb,xn)**2) / &
                             u(ind%dni,nb+1:nxb+nb,nb+1:nyb+nb,xn)
#endif /* ISO */

         case ("magx")
            if(ij=="yz") tab(:,:) = b(ind%bx,xn,nb+1:nyb+nb,nb+1:nzb+nb)
            if(ij=="xz") tab(:,:) = b(ind%bx,nb+1:nxb+nb,xn,nb+1:nzb+nb)
            if(ij=="xy") tab(:,:) = b(ind%bx,nb+1:nxb+nb,nb+1:nyb+nb,xn)
         case ("magy")
            if(ij=="yz") tab(:,:) = b(ind%by,xn,nb+1:nyb+nb,nb+1:nzb+nb)
            if(ij=="xz") tab(:,:) = b(ind%by,nb+1:nxb+nb,xn,nb+1:nzb+nb)
            if(ij=="xy") tab(:,:) = b(ind%by,nb+1:nxb+nb,nb+1:nyb+nb,xn)
         case ("magz")
            if(ij=="yz") tab(:,:) = b(ind%bz,xn,nb+1:nyb+nb,nb+1:nzb+nb)
            if(ij=="xz") tab(:,:) = b(ind%bz,nb+1:nxb+nb,xn,nb+1:nzb+nb)
            if(ij=="xy") tab(:,:) = b(ind%bz,nb+1:nxb+nb,nb+1:nyb+nb,xn)
         case default
            ierrh = -1
      end select


   end subroutine common_plt_hdf5

!<
!! \brief Routine calculating quantities for .hdf files
!<
   subroutine common_vars_hdf5(var,tab, ierrh)
      use fluidindex, only : ind
      use grid,   only : nb,nx,ny,nz
      use arrays, only : u,b

      implicit none
      character(LEN=4)     :: var
      real(kind=4), dimension(:,:,:) :: tab
      integer :: ierrh, i
      character(len=3) :: aux

      ierrh = 0

#ifdef COSM_RAYS
      select case(var(1:3))
         case("ecr")
           read(var,'(A3,I1)') aux,i
           tab(:,:,:) = real(u(ind%arr_crs(i),RNG),4)
      end select
#endif /* COSM_RAYS */
      select case(var)
         case("dend")
            tab(:,:,:) = real(u(ind%dnd,RNG),4)
         case("denn")
            tab(:,:,:) = real(u(ind%dnn,RNG),4)
         case("deni")
            tab(:,:,:) = real(u(ind%dni,RNG),4)
         case("vlxd")
            tab(:,:,:) = real(u(ind%mxd,RNG) / u(ind%dnd,RNG),4)
         case("vlxn")
            tab(:,:,:) = real(u(ind%mxn,RNG) / u(ind%dnn,RNG),4)
         case("vlxi")
            tab(:,:,:) = real(u(ind%mxi,RNG) / u(ind%dni,RNG),4)
         case("vlyd")
            tab(:,:,:) = real(u(ind%myd,RNG) / u(ind%dnd,RNG),4)
         case("vlyn")
            tab(:,:,:) = real(u(ind%myn,RNG) / u(ind%dnn,RNG),4)
         case("vlyi")
            tab(:,:,:) = real(u(ind%myi,RNG) / u(ind%dni,RNG),4)
         case("vlzd")
            tab(:,:,:) = real(u(ind%mzd,RNG) / u(ind%dnd,RNG),4)
         case("vlzn")
            tab(:,:,:) = real(u(ind%mzn,RNG) / u(ind%dnn,RNG),4)
         case("vlzi")
            tab(:,:,:) = real(u(ind%mzi,RNG) / u(ind%dni,RNG),4)
         case("enen")
#ifdef ISO
            tab(:,:,:) = real(0.5 *( u(ind%mxn,RNG)**2 + &
                                     u(ind%myn,RNG)**2 + &
                                     u(ind%mzn,RNG)**2 ) &
                              / u(ind%dnn,RNG),4)
#else
            tab(:,:,:) = real(u(ind%enn,RNG),4)
#endif
         case("enei")
#ifdef ISO
            tab(:,:,:) = real(0.5 *( u(ind%mxi,RNG)**2 + &
                                     u(ind%myi,RNG)**2 + &
                                     u(ind%mzi,RNG)**2 ) &
                              / u(ind%dni,RNG),4)
#else
            tab(:,:,:) = real(u(ind%eni,RNG),4)
#endif
         case("magx")
            tab(:,:,:) = real(b(ind%bx,RNG),4)
         case("magy")
            tab(:,:,:) = real(b(ind%by,RNG),4)
         case("magz")
            tab(:,:,:) = real(b(ind%bz,RNG),4)
         case default
            ierrh = -1
      end select

   end subroutine common_vars_hdf5

   subroutine write_plot(chdf)
      use hdf5, only  : HID_T, H5open_f, H5Fcreate_f, H5Gcreate_f, H5F_ACC_TRUNC_F, H5Gclose_f, &
         H5close_f, h5fclose_f
      use types, only : hdf
      use mpisetup, only : t, comm3d, ierr, proc
      implicit none
      type(hdf)     :: chdf !< Container for all necessary variables from dataio
      integer, save :: nimg = 0
      real, save    :: last_plt_time = 0.0
      character(LEN=32) ::fname
      integer :: i,error,fe
      logical, save :: first_entry = .true.
      integer(HID_T) :: file_id       !> File identifier
      integer(HID_T) :: gr_id,gr2_id  !> Group indentifier

      if( ((t-last_plt_time) > dt_plt) .and. dt_plt > 0.0 .or. first_entry) then
      fe = LEN(trim(chdf%log_file))
      fname = trim(chdf%log_file(1:fe-3)//"plt")
      call H5open_f(error)

      if(proc==0 .and. first_entry) then
         call H5Fcreate_f(fname, H5F_ACC_TRUNC_F, file_id, error)
         call H5Gcreate_f(file_id,"xy",gr_id,error)
         do i=1,nhdf_vars
            call H5Gcreate_f(gr_id,hdf_vars(i),gr2_id,error)
            call H5Gclose_f(gr2_id,error)
         enddo
         call H5Gclose_f(gr_id,error)
         call H5Gcreate_f(file_id,"yz",gr_id,error)
         do i=1,nhdf_vars
            call H5Gcreate_f(gr_id,hdf_vars(i),gr2_id,error)
            call H5Gclose_f(gr2_id,error)
         enddo
         call H5Gclose_f(gr_id,error)
         call H5Gcreate_f(file_id,"xz",gr_id,error)
         do i=1,nhdf_vars
            call H5Gcreate_f(gr_id,hdf_vars(i),gr2_id,error)
            call H5Gclose_f(gr2_id,error)
         enddo
         call H5Gclose_f(gr_id,error)
         call H5Fclose_f(file_id,error)
         first_entry = .false.
      endif
      call MPI_BARRIER(comm3d,ierr)
      do i = 1,nhdf_vars
        if(ix > 0) call write_plot_hdf5(chdf,hdf_vars(i),"yz",nimg)
        if(iy > 0) call write_plot_hdf5(chdf,hdf_vars(i),"xz",nimg)
        if(iz > 0) call write_plot_hdf5(chdf,hdf_vars(i),"xy",nimg)
      enddo

      nimg = nimg+1
      first_entry = .false.
      call H5close_f(error)

      last_plt_time = t

      endif

   end subroutine write_plot

   subroutine write_plot_hdf5(chdf,var,plane,nimg)
      use types,    only : hdf
      use mpisetup, only : MPI_CHARACTER, comm3d, ierr, pxsize, pysize, pzsize, MPI_DOUBLE_PRECISION, t, &
         pcoords, mpistop
      use hdf5,     only : HID_T, HSIZE_T, SIZE_T, H5F_ACC_RDWR_F, h5fopen_f, h5gopen_f, h5gclose_f, &
         h5fclose_f
      use h5lt,     only : h5ltmake_dataset_double_f, h5ltset_attribute_double_f
      use arrays,   only : u
      use grid,     only : nxb,nyb,nzb,nxd,nyd,nzd,nb
      use errh,     only : die

      implicit none
      character(LEN=2) :: plane
      logical, dimension(3) :: remain
      integer :: comm2d,lp,ls,xn,error
      logical :: ok_plt_var
      real, dimension(:,:), allocatable :: send
      real, dimension(:,:,:), allocatable :: temp
      real, dimension(:,:), allocatable :: img
      real, dimension(:), allocatable :: buff
      real :: imax,imin,di
      integer :: nimg, ierrh,i,j
      type(hdf) :: chdf !< Container for all necessary variables from dataio
      character(LEN=3) :: pij
      character(LEN=4) :: var    !> not yet implemented
      character(LEN=32) ::fname
      character(LEN=12) :: dname

      integer(HID_T) :: file_id       !> File identifier
      integer(HID_T) :: gr_id,gr2_id  !> Group indentifier
      integer(HSIZE_T), dimension(2) :: dims
      integer :: rank

      integer :: nib,nid,njb,njd,nkb,pisize,pjsize, fe
      integer(SIZE_T) :: bufsize

      rank = 2
      fe = LEN(trim(chdf%log_file))
      fname = trim(chdf%log_file(1:fe-3)//"plt")
      call MPI_BCAST(fname, 32, MPI_CHARACTER, 0, comm3d, ierr)

      nib = 0; nid = 0; njb = 0; njd = 0; nkb = 0; pisize = 0; pjsize = 0
      select case(plane)
         case("yz")
            xn     = ix + nb - pcoords(1)*nxb
            remain = (/.false.,.true.,.true./)
            pij    = "yz_"
            nib    = nyb
            nid    = nyd
            njb    = nzb
            njd    = nzd
            nkb    = nxb
            pisize = pysize
            pjsize = pzsize
         case("xz")
            xn     = iy + nb - pcoords(2)*nyb
            remain = (/.true.,.false.,.true./)
            pij    = "xz_"
            nib    = nxb
            nid    = nxd
            njb    = nzb
            njd    = nzd
            nkb    = nyb
            pisize = pxsize
            pjsize = pzsize
         case("xy")
            xn     = iz + nb - pcoords(3)*nzb
            remain = (/.true.,.true.,.false./)
            pij    = "xy_"
            nib    = nxb
            nid    = nxd
            njb    = nyb
            njd    = nyd
            nkb    = nzb
            pisize = pxsize
            pjsize = pysize
         case default
            call die("[dataio_hdf5:write_plot_hdf5] nonrecognized plane")
      end select

      dims(1) = nid
      dims(2) = njd
      call MPI_BARRIER(comm3d,ierr)
      call MPI_CART_SUB(comm3d,remain,comm2d,ierr)
      call MPI_COMM_SIZE(comm2d, ls, ierr)
      call MPI_COMM_RANK(comm2d, lp, ierr)
      if(xn > nb .and. xn <= nkb+nb) then
         allocate(temp(nib,njb,pisize*pjsize),img(nid,njd))
         allocate(buff(nid*njd))
         allocate(send(nib,njb))

         ok_plt_var = .false.
         call common_plt_hdf5(var,plane,xn,send,ierrh)
         if(ierrh==0) ok_plt_var = .true.
         if(ierrh==0) ok_plt_var = .true.
         if(.not.ok_plt_var) call die(var//" is not defined in common_plt_hdf5, neither in user_plt !!!")

         temp = -1.0
         call MPI_GATHER(send, nib*njb, MPI_DOUBLE_PRECISION, &
                         temp, nib*njb, MPI_DOUBLE_PRECISION, &
                         0, comm2d,ierr)

         if(lp == 0) then
            imax = maxval(temp); imin = minval(temp)
            di = imax-imin
            do i = 0, pisize-1
               do j = 0, pjsize-1
                  img(i*nib+1:(i+1)*nib,j*njb+1:(j+1)*njb) = &
                    temp(:,:,(j+1)+i*pjsize)
               enddo
            enddo
            call H5Fopen_f(fname, H5F_ACC_RDWR_F, file_id, error)
            call H5Gopen_f(file_id,plane,gr_id,error)
            call H5Gopen_f(gr_id,var,gr2_id,error)
            write(dname,'(a3,a4,a1,i4.4)') pij,var,"_",nimg
            call h5ltmake_dataset_double_f(gr2_id, dname, rank, dims, img, error)
            bufsize = 1
            call h5ltset_attribute_double_f(gr2_id,dname,"time",(/t/),bufsize,error)
            call H5Gclose_f(gr2_id,error)
            call H5Gclose_f(gr_id,error)
            call H5Fclose_f(file_id,error)
         endif
         if(allocated(send)) deallocate(send)
         if(allocated(temp)) deallocate(temp)
         if(allocated(img))  deallocate(img)
      endif
      call MPI_BARRIER(comm3d,ierr)

   end subroutine write_plot_hdf5

   subroutine write_restart_hdf5(filename,chdf)
      use hdf5, only : HID_T, HSIZE_T, HSSIZE_T, SIZE_T, H5F_ACC_TRUNC_F, h5open_f, h5pcreate_f, H5P_FILE_ACCESS_F, &
         h5pset_fapl_mpio_f, h5fcreate_f, h5pclose_f, H5P_DATASET_XFER_F, H5F_ACC_RDWR_F, H5S_SELECT_SET_F, &
         H5FD_MPIO_INDEPENDENT_F, H5T_NATIVE_DOUBLE, H5P_DATASET_CREATE_F, h5screate_simple_f, h5pset_chunk_f, &
         h5dcreate_f, h5sclose_f, h5dget_space_f, h5sselect_hyperslab_f, h5pset_dxpl_mpio_f, h5dwrite_f, &
         h5pset_chunk_f, h5dclose_f, h5fclose_f, h5close_f, h5fopen_f
      use h5lt, only : h5ltset_attribute_double_f, h5ltset_attribute_int_f, h5ltset_attribute_string_f
      use types, only : hdf
      use mpisetup, only : pcoords, comm3d, info, dt, t, proc, pxsize, pysize, pzsize, psize
      use grid, only : nx,ny,nz, x,y,z, nxb,nyb,nzb
      use arrays, only : u,b
      use grid, only : nxd,nyd,nzd,nb, xmin,xmax, &
          ymin,ymax, zmin,zmax
      use initproblem, only : problem_name, run_id
      use fluidindex, only : nvar
      IMPLICIT NONE
      type(hdf) :: chdf
      integer   :: llun,fe,nu
      CHARACTER(len=128) :: filename !> HDF File name
      character(len=128) :: lfile

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

      integer, dimension(3) :: dims
      integer :: error, rank = 4

      integer(SIZE_T) :: bufsize = 1

      nu = nvar%all

      llun  = chdf%log_lun
      lfile = chdf%log_file

      dims(1) = nx
      dims(2) = ny
      dims(3) = nz

      CALL h5open_f(error)
      CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
      CALL h5pset_fapl_mpio_f(plist_id, comm3d, info, error)

      CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)
      CALL h5pclose_f(plist_id, error)

#ifdef MASS_COMPENS
      rank = 3
      allocate(dimsf(rank),dimsfi(rank),chunk_dims(rank))
      allocate(count(rank),offset(rank),stride(rank),block(rank))
      !----------------------------------------------------------------------------------
      !  WRITE DINIT
      !
      dimsf = (/nx*pxsize,ny*pysize,nz*pzsize/) ! Dataset dimensions
      dimsfi = dimsf
      chunk_dims = (/nx,ny,nz/)                   ! Chunks dimensions

      ! Create the data space for the  dataset.
      CALL h5screate_simple_f(rank, dimsf, filespace, error)
      CALL h5screate_simple_f(rank, chunk_dims, memspace, error)

      ! Create chunked dataset.
      CALL h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
      CALL h5pset_chunk_f(plist_id, rank, chunk_dims, error)
      CALL h5dcreate_f(file_id, dname(3), H5T_NATIVE_DOUBLE, filespace, &
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
      CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dinit(:,:,:), dimsfi, error, &
                     file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

      CALL h5sclose_f(filespace, error)
      CALL h5sclose_f(memspace, error)
      CALL h5dclose_f(dset_id, error)

      if(allocated(dimsf))      deallocate(dimsf)
      if(allocated(dimsfi))     deallocate(dimsfi)
      if(allocated(chunk_dims)) deallocate(chunk_dims)
      if(allocated(count))      deallocate(count)
      if(allocated(offset))     deallocate(offset)
      if(allocated(stride))     deallocate(stride)
      if(allocated(block))      deallocate(block)
      !----------------------------------------------------------------------------------
#endif /* MASS_COMPENS */
      rank = 4
      allocate(dimsf(rank),dimsfi(rank),chunk_dims(rank))
      allocate(count(rank),offset(rank),stride(rank),block(rank))
      !----------------------------------------------------------------------------------
      !  WRITE FLUID VARIABLES
      !
      dimsf = (/nu,nx*pxsize,ny*pysize,nz*pzsize/) ! Dataset dimensions
      dimsfi = dimsf
      chunk_dims = (/nu,nx,ny,nz/)                   ! Chunks dimensions

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
      CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, u(:,:,:,:), dimsfi, error, &
                      file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

      CALL h5sclose_f(filespace, error)
      CALL h5sclose_f(memspace, error)
      CALL h5dclose_f(dset_id, error)
      CALL h5pclose_f(plist_id, error)
      !----------------------------------------------------------------------------------

      !----------------------------------------------------------------------------------
      !  WRITE MAG VARIABLES
      !
      dimsf = (/3,nx*pxsize,ny*pysize,nz*pzsize/) ! Dataset dimensions
      dimsfi = dimsf
      chunk_dims = (/3,nx,ny,nz/)                 ! Chunks dimensions

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
      CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, b(:,:,:,:), dimsfi, error, &
                      file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

      CALL h5sclose_f(filespace, error)
      CALL h5sclose_f(memspace, error)
      CALL h5dclose_f(dset_id, error)
      CALL h5pclose_f(plist_id, error)
      !----------------------------------------------------------------------------------
      if(allocated(dimsf))      deallocate(dimsf)
      if(allocated(dimsfi))     deallocate(dimsfi)
      if(allocated(chunk_dims)) deallocate(chunk_dims)
      if(allocated(count))      deallocate(count)
      if(allocated(offset))     deallocate(offset)
      if(allocated(stride))     deallocate(stride)
      if(allocated(block))      deallocate(block)

      rank = 1
      allocate(dimsf(rank),dimsfi(rank),chunk_dims(rank))
      allocate(count(rank),offset(rank),stride(rank),block(rank))

      !----------------------------------------------------------------------------------
      !  WRITE X Axis
      !
      dimsf  = (/nx*pxsize/) ! Dataset dimensions
      dimsfi = dimsf
      chunk_dims = (/nx/)    ! Chunks dimensions

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
      CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, x(:), dimsfi, error, &
                      file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

      CALL h5sclose_f(filespace, error)
      CALL h5sclose_f(memspace, error)
      CALL h5dclose_f(dset_id, error)
      CALL h5pclose_f(plist_id, error)
      !----------------------------------------------------------------------------------

      !----------------------------------------------------------------------------------
      !  WRITE Y Axis
      !
      dimsf  = (/ny*pysize/) ! Dataset dimensions
      dimsfi = dimsf
      chunk_dims = (/ny/)    ! Chunks dimensions

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
      CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, y(:), dimsfi, error, &
                      file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

      CALL h5sclose_f(filespace, error)
      CALL h5sclose_f(memspace, error)
      CALL h5dclose_f(dset_id, error)
      CALL h5pclose_f(plist_id, error)
      !----------------------------------------------------------------------------------

      !----------------------------------------------------------------------------------
      !  WRITE Z Axis
      !
      dimsf  = (/nz*pzsize/) ! Dataset dimensions
      dimsfi = dimsf
      chunk_dims = (/nz/)    ! Chunks dimensions

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
      CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, z(:), dimsfi, error, &
                      file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

      CALL h5sclose_f(filespace, error)
      CALL h5sclose_f(memspace, error)
      CALL h5dclose_f(dset_id, error)
      !----------------------------------------------------------------------------------
      CALL h5pclose_f(plist_id, error)
      CALL h5fclose_f(file_id, error)
      if(allocated(dimsf))      deallocate(dimsf)
      if(allocated(dimsfi))     deallocate(dimsfi)
      if(allocated(chunk_dims)) deallocate(chunk_dims)
      if(allocated(count))      deallocate(count)
      if(allocated(offset))     deallocate(offset)
      if(allocated(stride))     deallocate(stride)
      if(allocated(block))      deallocate(block)
      if(proc == 0) then
         CALL h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, error)

         bufsize = 1
         call h5ltset_attribute_double_f(file_id,"/","time", &
            (/t/),bufsize,error)
         call h5ltset_attribute_double_f(file_id,"/","timestep", &
            (/dt/),bufsize,error)
         call h5ltset_attribute_int_f(file_id,"/","nstep", &
            (/chdf%nstep/),bufsize,error)

         call h5ltset_attribute_int_f(file_id,"/","nres", &
            (/chdf%nres+1/),bufsize,error)
         call h5ltset_attribute_int_f(file_id,"/","nhdf", &
            (/chdf%nhdf/),bufsize,error)
         call h5ltset_attribute_int_f(file_id,"/","ntsl", &
            (/chdf%ntsl/),bufsize,error)
         call h5ltset_attribute_int_f(file_id,"/","nlog", &
            (/chdf%nlog/),bufsize,error)
         call h5ltset_attribute_int_f(file_id,"/","step_res", &
            (/chdf%nstep/),bufsize,error)
         call h5ltset_attribute_int_f(file_id,"/","step_hdf", &
            (/chdf%step_hdf/),bufsize,error)
         call h5ltset_attribute_double_f(file_id,"/","last_hdf_time", &
            (/chdf%last_hdf_time/),bufsize,error)

         call h5ltset_attribute_int_f(file_id,"/","nxd", &
            (/nxd/),bufsize,error)
         call h5ltset_attribute_int_f(file_id,"/","nyd", &
            (/nyd/),bufsize,error)
         call h5ltset_attribute_int_f(file_id,"/","nzd", &
            (/nzd/),bufsize,error)
         call h5ltset_attribute_int_f(file_id,"/","nxb", &
            (/nxb/),bufsize,error)
         call h5ltset_attribute_int_f(file_id,"/","nyb", &
            (/nyb/),bufsize,error)
         call h5ltset_attribute_int_f(file_id,"/","nzb", &
            (/nzb/),bufsize,error)
         call h5ltset_attribute_int_f(file_id,"/","nb", &
            (/nb/),bufsize,error)

         call h5ltset_attribute_double_f(file_id,"/","xmin", &
            (/xmin/),bufsize,error)
         call h5ltset_attribute_double_f(file_id,"/","xmax", &
            (/xmax/),bufsize,error)
         call h5ltset_attribute_double_f(file_id,"/","ymin", &
            (/ymin/),bufsize,error)
         call h5ltset_attribute_double_f(file_id,"/","ymax", &
            (/ymax/),bufsize,error)
         call h5ltset_attribute_double_f(file_id,"/","zmin", &
            (/zmin/),bufsize,error)
         call h5ltset_attribute_double_f(file_id,"/","zmax", &
            (/zmax/),bufsize,error)

         bufsize = 3
         call h5ltset_attribute_int_f(file_id,"/","psize", &
            psize,bufsize,error)

         fe = len(problem_name)
         call h5ltset_attribute_string_f(file_id,"/","problem name", &
            problem_name(1:fe),error)
         fe = len(chdf%domain)
         call h5ltset_attribute_string_f(file_id,"/","domain", &
            chdf%domain(1:fe),error)
         fe = len(run_id)
         call h5ltset_attribute_string_f(file_id,"/","run id", &
            run_id(1:fe),error)

         CALL h5fclose_f(file_id, error)

         open(llun, file=lfile, position='append')
         write(llun,*) 'Writing restart file: ',trim(filename)
         write(*,*)    'Writing restart file: ',trim(filename)
         close(llun)
      endif
      CALL h5close_f(error)

   end subroutine write_restart_hdf5

   subroutine read_restart_hdf5(chdf)
      use types, only : hdf
      use hdf5, only : HID_T, HSIZE_T, HSSIZE_T, SIZE_T, H5P_FILE_ACCESS_F, H5T_NATIVE_DOUBLE, &
          H5S_SELECT_SET_F, H5F_ACC_RDONLY_F, H5FD_MPIO_INDEPENDENT_F, H5P_DATASET_XFER_F, &
          h5open_f, h5pcreate_f, h5pset_fapl_mpio_f, h5fopen_f, h5pclose_f, h5dopen_f, &
          h5dget_space_f, h5sget_simple_extent_ndims_f, h5dget_create_plist_f, h5pget_chunk_f, &
          h5sselect_hyperslab_f, h5dread_f, h5sclose_f, h5pset_dxpl_mpio_f, h5dclose_f, &
          h5screate_simple_f, h5fclose_f, h5close_f
      use h5lt, only : h5ltget_attribute_double_f, h5ltget_attribute_int_f, h5ltget_attribute_string_f
      use mpisetup, only : MPI_CHARACTER, comm, ierr, pcoords, pxsize, pysize, pzsize, &
          MPI_CHARACTER, MPI_INTEGER, MPI_DOUBLE_PRECISION, proc, t, info, comm3d, dt, mpistop
      use fluidindex, only : nvar
      use grid, only : nx,ny,nz,x,y,z, nxb,nyb,nzb
      use arrays, only : u,b
      use grid, only : nxd,nyd,nzd,nb, xmin,xmax, &
          ymin,ymax, zmin,zmax
      use initproblem, only : problem_name, run_id

      IMPLICIT NONE
      type(hdf) :: chdf
      integer :: log_lun, nu
      CHARACTER(LEN=128) :: log_file  ! File name
      CHARACTER(LEN=128) :: filename  ! File name

      integer(HID_T) :: file_id       ! File identifier
      integer(HID_T) :: dset_id       ! Dataset identifier
      integer(HID_T) :: plist_id      ! Property list identifier
      integer(HID_T) :: filespace     ! Dataspace identifier in file
      integer(HID_T) :: memspace      ! Dataspace identifier in memory

      integer(HSIZE_T),  DIMENSION(:), allocatable :: count
      integer(HSSIZE_T), DIMENSION(:), allocatable :: offset
      integer(HSIZE_T),  DIMENSION(:), allocatable :: stride
      integer(HSIZE_T),  DIMENSION(:), allocatable :: block
      integer(HSIZE_T),  DIMENSION(:), allocatable :: dimsf, dimsfi, chunk_dims

      integer :: error, rank
      logical :: file_exist, log_exist

      real, dimension(1) :: rbuf
      integer, dimension(1) :: ibuf
      integer(SIZE_T) :: bufsize = 1

      nu = nvar%all

      log_lun  = chdf%log_lun
      log_file = chdf%log_file

      if(proc==0) then
         write (filename,'(a,a1,a3,a1,i4.4,a4)') &
            trim(problem_name),'_', run_id,'_',chdf%nres,'.res'
         write(*,*) 'Reading restart  file: ', trim(filename)
         inquire(file=log_file , exist = log_exist)
         if(log_exist .eqv. .true.) then
            open(log_lun, file=log_file, position='append')
            write(log_lun,*) 'Reading restart  file: ',trim(filename)
            close(log_lun)
         endif
      endif
      call MPI_BCAST(filename, 128,MPI_CHARACTER, 0, comm, ierr)

      inquire(file = filename, exist = file_exist)
      if(file_exist .eqv. .false.) then
         if(log_exist) then
            open(log_lun, file=log_file, position='append')
            write(log_lun,*) 'Restart  file: ', trim(filename), &
               ' does not exist.  ABORTING !!! '
            close(log_lun)
         endif

         write(*,*)       'Restart  file: ', trim(filename), &
            ' does not exist.  ABORTING !!! '
         call MPI_BARRIER(comm3d,ierr)
         call mpistop
         stop
      endif

      CALL h5open_f(error)
      CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
      CALL h5pset_fapl_mpio_f(plist_id, comm3d, info, error)

      CALL h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, error, access_prp = plist_id)
      CALL h5pclose_f(plist_id, error)

#ifdef MASS_COMPENS
      !----------------------------------------------------------------------------------
      !  READ DINIT
      !
      rank = 3
      allocate(dimsf(rank), dimsfi(rank), chunk_dims(rank))
      allocate(block(rank), offset(rank), count(rank), stride(rank))
      dimsf = (/nx*pxsize,ny*pysize,nz*pzsize/) ! Dataset dimensions
      dimsfi = dimsf

      ! Create chunked dataset.
      CALL h5dopen_f(file_id, dname(3), dset_id, error)

      call h5dget_space_f(dset_id, filespace, error)
      call H5sget_simple_extent_ndims_f (filespace,rank,error)
      call H5dget_create_plist_f (dset_id,plist_id,error)
      call h5pget_chunk_f(plist_id, rank, chunk_dims, error)

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
      CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, dinit, dimsfi, error, &
                     file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

      CALL h5sclose_f(filespace, error)
      CALL h5sclose_f(memspace, error)
      CALL h5pclose_f(plist_id, error)
      CALL h5dclose_f(dset_id, error)
      !----------------------------------------------------------------------------------
      if(allocated(dimsf))      deallocate(dimsf)
      if(allocated(dimsfi))     deallocate(dimsfi)
      if(allocated(chunk_dims)) deallocate(chunk_dims)
      if(allocated(count))      deallocate(count)
      if(allocated(offset))     deallocate(offset)
      if(allocated(stride))     deallocate(stride)
      if(allocated(block))      deallocate(block)
#endif /* MASS_COMPENS */
      !----------------------------------------------------------------------------------
      !  READ FLUID VARIABLES
      !
      rank = 4
      allocate(dimsf(rank), dimsfi(rank), chunk_dims(rank))
      allocate(block(rank), offset(rank), count(rank), stride(rank))
      dimsf = (/nu,nx*pxsize,ny*pysize,nz*pzsize/) ! Dataset dimensions
      dimsfi = dimsf

      ! Create chunked dataset.
      CALL h5dopen_f(file_id, dname(1), dset_id, error)

      call h5dget_space_f(dset_id, filespace, error)
      call H5sget_simple_extent_ndims_f (filespace,rank,error)
      call H5dget_create_plist_f (dset_id,plist_id,error)
      call h5pget_chunk_f(plist_id, rank, chunk_dims, error)

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
      CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, u, dimsfi, error, &
                     file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

      CALL h5sclose_f(filespace, error)
      CALL h5sclose_f(memspace, error)
      CALL h5pclose_f(plist_id, error)
      CALL h5dclose_f(dset_id, error)
      !----------------------------------------------------------------------------------

      !----------------------------------------------------------------------------------
      !  READ MAG VARIABLES
      !
      dimsf = (/3,nx*pxsize,ny*pysize,nz*pzsize/) ! Dataset dimensions
      dimsfi = dimsf

      ! Create chunked dataset.
      CALL h5dopen_f(file_id, dname(2), dset_id, error)

      call h5dget_space_f(dset_id, filespace, error)
      call H5Sget_simple_extent_ndims_f (filespace,rank,error)
      call H5Dget_create_plist_f (dset_id,plist_id,error)
      call h5pget_chunk_f(plist_id, rank, chunk_dims, error)

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
      CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, b, dimsfi, error, &
                     file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

      CALL h5sclose_f(filespace, error)
      CALL h5sclose_f(memspace, error)
      CALL h5pclose_f(plist_id, error)
      CALL h5dclose_f(dset_id, error)
      !----------------------------------------------------------------------------------
      CALL h5fclose_f(file_id, error)
      if(allocated(dimsf))      deallocate(dimsf)
      if(allocated(dimsfi))     deallocate(dimsfi)
      if(allocated(chunk_dims)) deallocate(chunk_dims)
      if(allocated(count))      deallocate(count)
      if(allocated(offset))     deallocate(offset)
      if(allocated(stride))     deallocate(stride)
      if(allocated(block))      deallocate(block)
      if(proc == 0) then
         CALL h5fopen_f (filename, H5F_ACC_RDONLY_F, file_id, error)
         bufsize = 1
         call h5ltget_attribute_double_f(file_id,"/","time", rbuf,error)
         t = rbuf(1)
         call h5ltget_attribute_double_f(file_id,"/","timestep", rbuf,error)
         dt = rbuf(1)
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

         call h5ltget_attribute_string_f(file_id,"/","problem name", problem_name,error)
         call h5ltget_attribute_string_f(file_id,"/","domain", chdf%domain,error)
         call h5ltget_attribute_string_f(file_id,"/","run id", chdf%new_id,error)

         CALL h5fclose_f(file_id, error)

         open(log_lun, file=log_file, position='append')
         write(log_lun,*) 'Done reading restart file: ',trim(filename)
         write(*,*)    'Done reading restart file: ',trim(filename)
         close(log_lun)
      endif

      call MPI_BCAST(chdf%nstep, 1, MPI_INTEGER, 0, comm3d, ierr)
      call MPI_BCAST(chdf%nres, 1, MPI_INTEGER, 0, comm3d, ierr)
      call MPI_BCAST(chdf%nhdf, 1, MPI_INTEGER, 0, comm3d, ierr)
      call MPI_BCAST(chdf%ntsl, 1, MPI_INTEGER, 0, comm3d, ierr)
      call MPI_BCAST(chdf%nlog, 1, MPI_INTEGER, 0, comm3d, ierr)
      call MPI_BCAST(chdf%step_res, 1, MPI_INTEGER, 0, comm3d, ierr)
      call MPI_BCAST(chdf%step_hdf, 1, MPI_INTEGER, 0, comm3d, ierr)

      call MPI_BCAST(chdf%last_hdf_time, 1, MPI_DOUBLE_PRECISION, 0, comm3d, ierr)
      call MPI_BCAST(t, 1, MPI_DOUBLE_PRECISION, 0, comm3d, ierr)
      call MPI_BCAST(dt, 1, MPI_DOUBLE_PRECISION, 0, comm3d, ierr)

      CALL MPI_BCAST(problem_name, 32, MPI_CHARACTER, 0, comm3d, ierr)
      CALL MPI_BCAST(chdf%domain, 16, MPI_CHARACTER, 0, comm3d, ierr)
      CALL MPI_BCAST(chdf%new_id, 3, MPI_CHARACTER, 0, comm3d, ierr)
      CALL h5close_f(error)

   end subroutine read_restart_hdf5
!
! ------------------------------------------------------------------------------------
!
   subroutine write_hdf5(chdf)
      use hdf5, only : HID_T, SIZE_T, H5F_ACC_RDWR_F, H5P_FILE_ACCESS_F, h5open_f, h5pcreate_f, h5pset_fapl_mpio_f, &
         h5pclose_f, h5close_f, H5F_ACC_TRUNC_F, H5P_DEFAULT_F, h5fcreate_f, h5fclose_f, h5fopen_f
      use h5lt, only : h5ltset_attribute_double_f, h5ltset_attribute_int_f, h5ltset_attribute_string_f
      use types, only : hdf
      use mpisetup, only: pcoords, comm3d, proc, info, psize,ierr, mpistop, t, dt
      use grid, only : nxb,nyb,nzb,nx,ny,nz,nxd,nyd,nzd,nb, xmin,xmax, &
         ymin,ymax, zmin,zmax
     use initproblem, only : problem_name, run_id
     use list_hdf5, only : write_arr
#ifdef NEW_HDF5
     use list_hdf5, only : iterate_lhdf5
#endif /* NEW_HDF5 */

      IMPLICIT NONE
      type(hdf) :: chdf
      integer(HID_T) :: file_id       ! File identifier
      integer(HID_T) :: plist_id      ! Property list identifier
      integer :: llun,fe, ierrh
      logical :: ok_var
      CHARACTER(len=128) :: lfile     ! File name
      character(len=4)   :: dd
      CHARACTER(LEN=32)  :: fname

      real(kind=4), allocatable :: data (:,:,:)  ! Data to write
      integer :: error, i
      integer(SIZE_T) :: bufsize = 1

      ! Initialize HDF5 library and Fortran interfaces.
      !
      llun = chdf%log_lun
      lfile = chdf%log_file
      write(dd,'(i4.4)') chdf%nhdf
      fname = trim(problem_name)//"_"//trim(run_id)//"_"//dd//".h5"

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
      if(.not.allocated(data)) allocate(data(nxb,nyb,nzb))
      ierrh = 0; ok_var = .false.
      do i = 1, nhdf_vars
         call common_vars_hdf5(hdf_vars(i),data,ierrh);  if(ierrh == 0) ok_var = .true.
         if(.not.ok_var) then
            write(*,*) hdf_vars(i),' is not defined in common_vars_hdf5, neither in user_hdf5 !!!'
            call MPI_BARRIER(comm3d,ierr)
            call mpistop
            stop
         endif
         call write_arr(data,hdf_vars(i),file_id)
      enddo
     if(allocated(data)) deallocate(data)
#ifdef NEW_HDF5
     call iterate_lhdf5(file_id)
#endif /* NEW_HDF5 */
     !
     ! Close the property list.
     !
     CALL h5fclose_f(file_id, error)
     if(proc == 0) then
         CALL h5fopen_f (fname, H5F_ACC_RDWR_F, file_id, error)

         bufsize = 1
         call h5ltset_attribute_double_f(file_id,"/","time",     (/t/),bufsize,error)
         call h5ltset_attribute_double_f(file_id,"/","timestep", (/dt/),bufsize,error)
         call h5ltset_attribute_int_f(file_id,"/","nstep",       (/chdf%nstep/),bufsize,error)

         call h5ltset_attribute_int_f(file_id,"/","nxd", (/nxd/),bufsize,error)
         call h5ltset_attribute_int_f(file_id,"/","nyd", (/nyd/),bufsize,error)
         call h5ltset_attribute_int_f(file_id,"/","nzd", (/nzd/),bufsize,error)
         call h5ltset_attribute_int_f(file_id,"/","nxb", (/nxb/),bufsize,error)
         call h5ltset_attribute_int_f(file_id,"/","nyb", (/nyb/),bufsize,error)
         call h5ltset_attribute_int_f(file_id,"/","nzb", (/nzb/),bufsize,error)
         call h5ltset_attribute_int_f(file_id,"/","nb",  (/nb/),bufsize,error)

         call h5ltset_attribute_double_f(file_id,"/","xmin", (/xmin/),bufsize,error)
         call h5ltset_attribute_double_f(file_id,"/","xmax", (/xmax/),bufsize,error)
         call h5ltset_attribute_double_f(file_id,"/","ymin", (/ymin/),bufsize,error)
         call h5ltset_attribute_double_f(file_id,"/","ymax", (/ymax/),bufsize,error)
         call h5ltset_attribute_double_f(file_id,"/","zmin", (/zmin/),bufsize,error)
         call h5ltset_attribute_double_f(file_id,"/","zmax", (/zmax/),bufsize,error)

         bufsize = 3
         call h5ltset_attribute_int_f(file_id,"/","psize", psize,bufsize,error)

         fe = len(problem_name)
         call h5ltset_attribute_string_f(file_id,"/","problem name", problem_name(1:fe),error)
         fe = len(chdf%domain)
         call h5ltset_attribute_string_f(file_id,"/","domain", chdf%domain(1:fe),error)
         call h5ltset_attribute_string_f(file_id,"/","run id", run_id(1:3),error)

         CALL h5fclose_f(file_id, error)
         open(llun, file=lfile, position='append')
         write(llun,*) 'Writing output   file: ',trim(fname)
         write(*,*)       'Writing output   file: ',trim(fname)
         close(llun)
      endif
      call MPI_BARRIER(comm3d,ierr)
      CALL h5close_f(error)

   end subroutine write_hdf5

   subroutine write_arr(data,dsetname,file_id)
      use hdf5, only : HID_T, HSIZE_T, HSSIZE_T, h5screate_simple_f, h5pcreate_f, h5pset_chunk_f, &
         h5sclose_f, h5pset_dxpl_mpio_f, h5dwrite_f, h5dclose_f, H5P_DATASET_XFER_F, H5P_DATASET_CREATE_F, &
         H5T_NATIVE_REAL, H5S_SELECT_SET_F, H5FD_MPIO_INDEPENDENT_F, h5dcreate_f, h5dget_space_f, &
         h5sselect_hyperslab_f
      use grid, only : nxd,nyd,nzd,nb
      use grid, only : nxb,nyb,nzb
      use mpisetup, only: pcoords

      implicit none
      real(kind=4), dimension(:,:,:) :: data

      integer :: rank = 3 ! Dataset rank

      CHARACTER(LEN=4) :: dsetname    ! Dataset name

      integer(HID_T) :: file_id       ! Dataset identifier
      integer(HID_T) :: dset_id       ! Dataset identifier
      integer(HID_T) :: plist_id       ! Dataset identifier
      integer(HID_T) :: filespace     ! Dataspace identifier in file
      integer(HID_T) :: memspace      ! Dataspace identifier in memory

      integer, parameter :: ndims = 3
      integer(HSIZE_T),  DIMENSION(ndims) :: count
      integer(HSSIZE_T), DIMENSION(ndims) :: offset
      integer(HSIZE_T),  DIMENSION(ndims) :: stride
      integer(HSIZE_T),  DIMENSION(ndims) :: block
      integer(HSIZE_T), DIMENSION(ndims) :: dimsf, dimsfi, chunk_dims
      integer :: error

      dimsf = (/nxd,nyd,nzd/) ! Dataset dimensions
      dimsfi = dimsf
      chunk_dims = (/nxb,nyb,nzb/) ! Chunks dimensions
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

end module dataio_hdf5

