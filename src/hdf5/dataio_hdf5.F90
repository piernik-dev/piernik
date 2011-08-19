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

   use constants,     only: FLUID, MAG, xdim, zdim
   use dataio_pub,    only: maxparfilelen, maxparfilelines
   use list_hdf5,     only: S_LEN

   implicit none

   private
   public :: init_hdf5, read_restart_hdf5, cleanup_hdf5, write_hdf5, write_restart_hdf5, write_plot, write_arr_to_restart, read_arr_from_restart
   public :: parfile, parfilelines

   integer, parameter :: dnamelen=5
   integer, parameter :: planelen = 2           !< length of plane names e.g. "xy", "yz", "rp" etc.
   character(len=planelen), dimension(xdim:zdim), parameter :: pl_id = [ "yz", "xz", "xy" ]
   character(len=dnamelen), dimension(FLUID:MAG) :: dname = [ "fluid", "mag  " ]  !< dataset names for restart files

   character(len=S_LEN), allocatable, dimension(:) :: hdf_vars  !< dataset names for hdf files
   integer :: nhdf_vars !< number of quantities plotted to hdf files
   integer, dimension(xdim:zdim) :: pl_i !< no. of cell ( 1 <= pl_i(:) < dom%n_d(:) ) for YZ, XZ and XY slices in plt files
   real :: dt_plt !< frequency of plt output

   ! storage for the problem.par
   character(len=maxparfilelen), dimension(maxparfilelines) :: parfile !< contents of the parameter file
   integer, save                             :: parfilelines = 0       !< number of lines in the parameter file

   interface write_arr_to_restart
      module procedure write_4darr_to_restart, write_3darr_to_restart
   end interface

   interface read_arr_from_restart
      module procedure read_4darr_from_restart, read_3darr_from_restart
   end interface

contains

!>
!! \brief Procedure initializing HDF5 module
!<

   subroutine init_hdf5(vars, tix, tiy, tiz, tdt_plt)

      use constants,   only: varlen
      use fluids_pub,  only: has_ion, has_dst, has_neu
      use fluidindex,  only: iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz
#ifdef COSM_RAYS
      use dataio_pub,  only: warn, msg
      use fluidindex,  only: iarr_all_crs
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

      pl_i(:) = [ tix, tiy, tiz ]
      dt_plt = tdt_plt

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
               if (has_dst) nhdf_vars = nhdf_vars - 1
#ifdef GRAV
            case ('gpot')
               nhdf_vars = nhdf_vars + 1
#ifdef MULTIGRID
            case ('mgso') ! multigrid solution
               nhdf_vars = nhdf_vars + 1
#endif /* MULTIGRID */
#endif /* GRAV */
            case ('magx', 'magy', 'magz', 'pres')
               nhdf_vars = nhdf_vars + 1
#ifdef COSM_RAYS
            case ('encr')
#ifndef NEW_HDF5
               nhdf_vars = nhdf_vars + size(iarr_all_crs,1)
#endif /* !NEW_HDF5 */
#endif /* COSM_RAYS */
            case default
               nhdf_vars = nhdf_vars + 1
         end select
      enddo
      allocate(hdf_vars(nhdf_vars)); j = 1
      do i = 1, nvars
         select case (vars(i))
            case ('dens')
               if (has_dst) then ; hdf_vars(j) = 'dend' ; j = j + 1 ; endif
               if (has_neu) then ; hdf_vars(j) = 'denn' ; j = j + 1 ; endif
               if (has_ion) then ; hdf_vars(j) = 'deni' ; j = j + 1 ; endif
            case ('velx')
               if (has_dst) then ; hdf_vars(j) = 'vlxd' ; j = j + 1 ; endif
               if (has_neu) then ; hdf_vars(j) = 'vlxn' ; j = j + 1 ; endif
               if (has_ion) then ; hdf_vars(j) = 'vlxi' ; j = j + 1 ; endif
            case ('vely')
               if (has_dst) then ; hdf_vars(j) = 'vlyd' ; j = j + 1 ; endif
               if (has_neu) then ; hdf_vars(j) = 'vlyn' ; j = j + 1 ; endif
               if (has_ion) then ; hdf_vars(j) = 'vlyi' ; j = j + 1 ; endif
            case ('velz')
               if (has_dst) then ; hdf_vars(j) = 'vlzd' ; j = j + 1 ; endif
               if (has_neu) then ; hdf_vars(j) = 'vlzn' ; j = j + 1 ; endif
               if (has_ion) then ; hdf_vars(j) = 'vlzi' ; j = j + 1 ; endif
            case ('ener')
               if (has_neu) then ; hdf_vars(j) = 'enen' ; j = j + 1 ; endif
               if (has_ion) then ; hdf_vars(j) = 'enei' ; j = j + 1 ; endif
            case ("magx", "magy", "magz")
               hdf_vars(j) = vars(i) ; j = j + 1
#ifdef COSM_RAYS
            case ('encr')
#ifndef NEW_HDF5
               do k = 1, size(iarr_all_crs,1)
                  if (k<=9) then
                     write(aux,'(A3,I1)') 'ecr', k
                     hdf_vars(j) = aux ; j = j + 1
                  else
                     write(msg, '(a,i3)')"[dataio_hdf5:init_hdf5] Cannot create name for CR energy component #", k
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
               if (has_neu) then ; hdf_vars(j) = 'pren' ; j = j + 1 ; endif
               if (has_ion) then ; hdf_vars(j) = 'prei' ; j = j + 1 ; endif
            case default
               hdf_vars(j) = trim(vars(i)) ; j = j + 1
         end select
      enddo

   end subroutine init_hdf5

!>
!! \brief Procedure finalizing HDF5 module
!<
   subroutine cleanup_hdf5
      implicit none

      if (allocated(hdf_vars)) deallocate(hdf_vars)

   end subroutine cleanup_hdf5

!>
!! \brief decode some useful indices from variable name, if possible
!<
   subroutine common_shortcuts(var, fl_dni, i_xyz)

      use constants,  only: varlen, singlechar
      use fluidindex, only: flind
      use fluidtypes, only: component_fluid

      implicit none

      character(len=varlen), intent(in) :: var
      type(component_fluid), pointer, intent(inout) :: fl_dni
      integer, intent(out) :: i_xyz

      character(len=singlechar) :: dc

      nullify(fl_dni)
      if (any([ "den", "vlx", "vly", "vlz", "ene" ] == var(1:3))) then
         select case (var(4:4))
            case ("d")
               fl_dni => flind%dst
            case ("n")
               fl_dni => flind%neu
            case ("i")
               fl_dni => flind%ion
         end select
      endif

      i_xyz = huge(1)
      if (var(1:2) == "vl") then
         dc = var(3:3)
      else if (var(1:3) == "mag") then
         dc = var(4:4)
      else
         dc = '_'
      endif
      if (any([ "x", "y", "z" ] == dc)) i_xyz = ichar(dc) - ichar("x")

   end subroutine common_shortcuts

!>
!! \brief Routine calculating quantities for plot files
!<
   subroutine common_plt_hdf5(var, ij, xn, tab, ierrh, cg)

      use constants,  only: varlen, xdim, ydim, zdim
      use grid_cont,  only: grid_container
      use fluidindex, only: flind, ibx
      use fluidtypes, only: component_fluid
#ifdef COSM_RAYS
      use fluidindex, only: iarr_all_crs
#endif /* COSM_RAYS */

      implicit none

      character(len=varlen), intent(in) :: var !< quantity to be plotted
      integer, intent(in)               :: ij  !< direction perpendicular to the plane of plot, xdim means "yz" plane
      integer(kind=8), intent(in)       :: xn  !< no. of cell at which we are slicing the local block
      integer, intent(out)              :: ierrh !< error handling
      real, dimension(:,:), intent(out) :: tab !< array containing given quantity
      type(grid_container), pointer, intent(in) :: cg

      integer :: is, ie, js, je, ks, ke, i_xyz
      type(component_fluid), pointer :: fl_dni
#ifdef COSM_RAYS
      integer :: i
#endif /* COSM_RAYS */

      ierrh = 0

      is=cg%is; ie=cg%ie
      js=cg%js; je=cg%je
      ks=cg%ks; ke=cg%ke
      select case (ij)
         case (xdim)
            is = int(xn, kind=4)
            ie = int(xn, kind=4)
         case (ydim)
            js = int(xn, kind=4)
            je = int(xn, kind=4)
         case (zdim)
            ks = int(xn, kind=4)
            ke = int(xn, kind=4)
      end select

      call common_shortcuts(var, fl_dni, i_xyz)

      select case (var)
         case ("dend", "deni", "denn")
            tab(:,:) = reshape(cg%u%arr(fl_dni%idn, is:ie, js:je, ks:ke), shape(tab))
         case ("vlxd", "vlxn", "vlxi", "vlyd", "vlyn", "vlyi", "vlzd", "vlzn", "vlzi")
            tab(:,:) = reshape(cg%u%arr(fl_dni%imx + i_xyz, is:ie, js:je, ks:ke) / cg%u%arr(fl_dni%idn, is:ie, js:je, ks:ke), shape(tab))
         case ("enen", "enei")
#ifdef ISO
            tab(:,:) = reshape(0.5 * ( cg%u%arr(fl_dni%imx, is:ie, js:je, ks:ke)**2 + &
                 &                     cg%u%arr(fl_dni%imy, is:ie, js:je, ks:ke)**2 + &
                 &                     cg%u%arr(fl_dni%imz, is:ie, js:je, ks:ke)**2 ) / &
                 &                     cg%u%arr(fl_dni%idn, is:ie, js:je, ks:ke), shape(tab))
#else /* !ISO */
            tab(:,:) = reshape(cg%u%arr(fl_dni%ien, is:ie, js:je, ks:ke), shape(tab))
#endif /* !ISO */
         case ('prei')
            tab(:,:) = 0.0
         case ("pren")
#ifdef ISO
            tab = 0.0
#else /* !ISO */
            ! BEWARE: Why there is only one case here?
            if (ij==zdim) then
               tab(:,:) = real( cg%u%arr(flind%neu%ien, is:ie, js:je, xn) - &
                    0.5 *( cg%u%arr(flind%neu%imx, is:ie, js:je, xn)**2 + &
                    &      cg%u%arr(flind%neu%imy, is:ie, js:je, xn)**2 + &
                    &      cg%u%arr(flind%neu%imz, is:ie, js:je, xn)**2 ) / cg%u%arr(flind%neu%idn, is:ie, js:je, xn), kind=4)*(flind%neu%gam_1)
            endif
#endif /* !ISO */
         case ("magx", "magy", "magz")
            tab(:,:) = reshape(cg%b%arr(ibx + i_xyz, is:ie, js:je, ks:ke), shape(tab))
#ifdef GRAV
         case ("gpot")
            tab(:,:) = reshape(cg%gpot%arr(is:ie, js:je, ks:ke), shape(tab))
#endif /* GRAV */
#ifdef COSM_RAYS
         case ("ecr*")
            i = iarr_all_crs(ichar(var(4:4))-ichar('0'))
            tab(:,:) = reshape(cg%u%arr(i, is:ie, js:je, ks:ke), shape(tab))
#endif /* COSM_RAYS */
         case default
            ierrh = -1
      end select

   end subroutine common_plt_hdf5

!>
!! \brief Routine calculating quantities for .hdf files
!<
   subroutine common_vars_hdf5(var, tab, ierrh, cg)

      use constants,  only: varlen
      use fluidindex, only: flind, ibx, iby, ibz
      use fluidtypes, only: component_fluid
      use grid_cont,  only: grid_container

      implicit none

      character(len=varlen), intent(in) :: var
      real(kind=4), dimension(:,:,:)    :: tab
      integer, intent(out)              :: ierrh
      type(grid_container), pointer, intent(in) :: cg

      type(component_fluid), pointer :: fl_dni
      integer :: i_xyz
#ifdef COSM_RAYS
      integer :: i
      integer, parameter    :: auxlen = varlen - 1
      character(len=auxlen) :: aux
#endif /* COSM_RAYS */

      call common_shortcuts(var, fl_dni, i_xyz)

      ierrh = 0
      tab = 0.0
      select case (var)
#ifdef COSM_RAYS
         case ("ecr1" : "ecr9")
            read(var,'(A3,I1)') aux, i !> \deprecated BEWARE 0 <= i <= 9, no other indices can be dumped to hdf file
            tab(:,:,:) = real(cg%u%arr(flind%crs%beg+i-1, RNG), kind=4)
#endif /* COSM_RAYS */
         case ("dend", "deni", "denn")
            tab(:,:,:) = real(cg%u%arr(fl_dni%idn, RNG), kind=4)
         case ("vlxd", "vlxn", "vlxi", "vlyd", "vlyn", "vlyi", "vlzd", "vlzn", "vlzi")
            tab(:,:,:) = real(cg%u%arr(fl_dni%imx + i_xyz, RNG) / cg%u%arr(fl_dni%idn, RNG), kind=4)
         case ("enen", "enei")
#ifdef ISO
            tab(:,:,:) = real(0.5 *( cg%u%arr(fl_dni%imx, RNG)**2 + &
                 &                   cg%u%arr(fl_dni%imy, RNG)**2 + &
                 &                   cg%u%arr(fl_dni%imz, RNG)**2 ) / cg%u%arr(fl_dni%idn, RNG), kind=4)
#else /* !ISO */
            tab(:,:,:) = real(cg%u%arr(fl_dni%ien, RNG), kind=4)
#endif /* !ISO */
#ifdef NEUTRAL
         case ("pren")
#ifndef ISO
            tab(:,:,:) = real( cg%u%arr(flind%neu%ien, RNG) - 0.5 * ( &
                 &             cg%u%arr(flind%neu%imx, RNG)**2 + &
                 &             cg%u%arr(flind%neu%imy, RNG)**2 + &
                 &             cg%u%arr(flind%neu%imz, RNG)**2 ) / cg%u%arr(flind%neu%idn, RNG), kind=4) * real(flind%neu%gam_1, kind=4)
#endif /* !ISO */
#endif /* NEUTRAL */
         case ("prei")
#ifndef ISO
            tab(:,:,:) = real( cg%u%arr(flind%ion%ien, RNG) - 0.5 *( &
                 &             cg%u%arr(flind%ion%imx, RNG)**2 + &
                 &             cg%u%arr(flind%ion%imy, RNG)**2 + &
                 &             cg%u%arr(flind%ion%imz, RNG)**2 ) / cg%u%arr(flind%ion%idn, RNG), kind=4) * real(flind%ion%gam_1, kind=4) - &
                 &       real( 0.5*(flind%ion%gam_1)*(cg%b%arr(ibx, RNG)**2 + cg%b%arr(iby, RNG)**2 + cg%b%arr(ibz, RNG)**2), kind=4)
#endif /* !ISO */
         case ("magx", "magy", "magz")
            tab(:,:,:) = real(cg%b%arr(ibx + i_xyz, RNG), kind=4)
         case ("gpot")
            if (associated(cg%gpot%arr)) tab(:,:,:) = real(cg%gpot%arr(RNG), kind=4)
         case ("mgso")
            if (associated(cg%sgp%arr))  tab(:,:,:) = real(cg%sgp%arr(RNG),  kind=4)
         case default
            ierrh = -1
      end select

   end subroutine common_vars_hdf5

   subroutine write_plot

      use constants,  only: cwdlen, xdim, zdim
      use dataio_pub, only: log_file
      use global,     only: t
      use hdf5,       only: HID_T, H5open_f, H5Fcreate_f, H5Gcreate_f, H5F_ACC_TRUNC_F, H5Gclose_f, H5close_f, h5fclose_f
      use mpisetup,   only: comm, ierr, master

      implicit none

      integer, save     :: nimg = 0
      integer(kind=4)   :: error
      real, save        :: last_plt_time = 0.0
      character(len=cwdlen) :: fname
      integer           :: i, d
      logical, save     :: first_entry = .true.
      integer(HID_T)    :: file_id                 !> File identifier
      integer(HID_T)    :: gr_id, gr2_id           !> Groups identifier

      if ( ((t-last_plt_time > dt_plt) .or. first_entry ) .and. dt_plt > 0.0 ) then

         write(fname,'(2a)') trim(log_file(1:len_trim(log_file)-3)),"plt"
         call H5open_f(error)

         if (master .and. first_entry) then
            call H5Fcreate_f(fname, H5F_ACC_TRUNC_F, file_id, error)

            do d = xdim, zdim
               if (pl_i(d) > 0) then
                  call H5Gcreate_f(file_id, pl_id(d), gr_id, error)
                  do i=1, nhdf_vars
                     call H5Gcreate_f(gr_id, hdf_vars(i), gr2_id, error)
                     call H5Gclose_f(gr2_id, error)
                  enddo
                  call H5Gclose_f(gr_id, error)
               endif
            enddo

            call H5Fclose_f(file_id, error)
         endif

         call MPI_Barrier(comm, ierr)

         do i = 1, nhdf_vars
            do d = xdim, zdim
               if (pl_i(d) > 0) call write_plot_hdf5(hdf_vars(i), d, nimg)
            enddo
         enddo

         nimg = nimg+1
         first_entry = .false.
         call H5close_f(error)

         last_plt_time = t

      endif

   end subroutine write_plot

! /details This routine appends a nimg-th slice of data, taken perpendiclularly to the plane direction, of the var field
!
! In AMR it will probably work only on base level (or any other regular, uniform grid) unless someone decides to complicate the format
!
! I/O is not paralelized here because it is not straightforwart to paralelize, when vizit is .true.
!
! /todo non-blocking receive ?

   subroutine write_plot_hdf5(var, plane, nimg)

      use constants,   only: xdim, ydim, zdim, varlen, cwdlen, LO, HI
      use dataio_pub,  only: vizit, fmin, fmax, log_file, msg, die, warn
      use dataio_user, only: user_plt_hdf5
      use domain,      only: dom, has_dir
      use global,      only: t
      use grid,        only: cga
      use grid_cont,   only: grid_container!, cg_list_element
      use hdf5,        only: HID_T, HSIZE_T, SIZE_T, H5F_ACC_RDWR_F, h5fopen_f, h5gopen_f, h5gclose_f, h5fclose_f
      use h5lt,        only: h5ltmake_dataset_double_f, h5ltset_attribute_double_f
      use mpi,         only: MPI_DOUBLE_PRECISION
      use mpisetup,    only: comm, ierr, proc, FIRST, LAST, status, master
#ifdef PGPLOT
      use viz,         only: draw_me
#endif /* PGPLOT */

      implicit none

      integer, intent(in)                 :: plane  ! xdim means "yz" and so on
      character(len=varlen), intent(in)   :: var                           !> not yet implemented
      integer, intent(in)                 :: nimg

      real, dimension(:,:), allocatable   :: send, img, recv
      integer                             :: ierrh, p
      integer(kind=4)                     :: error
      integer(kind=8)                     :: xn, xn_r
      character(len=cwdlen)               :: fname
      integer, parameter                  :: vdn_len = 12
      character(len=vdn_len)              :: vdname
      integer(HID_T)                      :: file_id                       !> File identifier
      integer(HID_T)                      :: gr_id, gr2_id                  !> Group identifier
      integer(kind=4), parameter                  :: rank = 2
      integer(HSIZE_T), dimension(rank)   :: dims
      integer, dimension(xdim:zdim), parameter :: d1 = [ ydim, xdim, xdim ] , d2 = [ zdim, zdim, ydim ] ! d1(d) and d2(2) are perpendicular to direction d
      integer(SIZE_T), parameter          :: bufsize = 1
      real, dimension(bufsize)            :: timebuffer
      integer, parameter                  :: tag = 101
      type(grid_container), pointer :: cg

      cg => cga%cg_all(1)
      if (ubound(cga%cg_all(:), dim=1) > 1) call die("[dataio_hdf5:write_plot_hdf5] multiple grid pieces per procesor not implemented yet") !nontrivial message tagging

      xn = 1
      if (has_dir(plane)) xn = pl_i(plane) + cg%nb - cg%off(plane)

      if ((xn > cg%nb .and. xn <= cg%n_b(plane)+cg%nb) .or. (xn == 1 .and. .not. has_dir(plane))) then
         allocate(send(cg%n_b(d1(plane)), cg%n_b(d2(plane))))
         call common_plt_hdf5(var, plane, xn, send, ierrh, cg)
         if (associated(user_plt_hdf5) .and. ierrh /= 0) call user_plt_hdf5(var, plane, xn, send, ierrh, cg)
         if (ierrh /= 0) then
            write(msg,'(2a)') var, " is not defined in common_plt_hdf5, neither in user_plt_hdf5 !!!"
            call die(msg)
         endif
      endif

      if (master) then

         write(fname,'(2a)') trim(log_file(1:len_trim(log_file)-3)),"plt"

         dims(:) = [ dom%n_d(d1(plane)), dom%n_d(d2(plane)) ] ! Dataset dimensions
         allocate(img(dims(1), dims(2)))

         do p = FIRST, LAST
            xn_r = 1
            if (has_dir(plane)) xn_r = pl_i(plane) + cg%nb - dom%pse(p)%sel(1, plane, LO)
            if ((xn_r > cg%nb .and. xn_r <= int(dom%pse(p)%sel(1, plane, HI) - dom%pse(p)%sel(1, plane, LO) + 1, 4) + cg%nb) .or. (xn_r == 1 .and. .not. has_dir(plane))) then
               if (p == proc) then
                  img(1+dom%pse(p)%sel(1, d1(plane), LO):1+dom%pse(p)%sel(1, d1(plane), HI), 1+dom%pse(p)%sel(1, d2(plane), LO):1+dom%pse(p)%sel(1, d2(plane), HI)) = send(:,:)
               else
                  allocate(recv(dom%pse(p)%sel(1, d1(plane), HI)-dom%pse(p)%sel(1, d1(plane), LO)+1, dom%pse(p)%sel(1, d2(plane), HI)-dom%pse(p)%sel(1, d2(plane), LO)+1))
                  call MPI_Recv(recv, size(recv), MPI_DOUBLE_PRECISION, p, tag, comm, status(:,p), ierr)
                  img(1+dom%pse(p)%sel(1, d1(plane), LO):1+dom%pse(p)%sel(1, d1(plane), HI), 1+dom%pse(p)%sel(1, d2(plane), LO):1+dom%pse(p)%sel(1, d2(plane), HI)) = recv(:,:)
                  deallocate(recv)
               endif
            endif
         enddo

         if (vizit) then
#ifdef PGPLOT
            call draw_me(real(img, kind=4), real(fmin, kind=4), real(fmax, kind=4))
#else /* !PGPLOT */
            call warn("[dataio_hdf5:write_plot_hdf5] vizit used without PGPLOT")
#endif /* !PGPLOT */
         else
            call H5Fopen_f(fname, H5F_ACC_RDWR_F, file_id, error)
            call H5Gopen_f(file_id, pl_id(plane), gr_id, error)
            call H5Gopen_f(gr_id, var, gr2_id, error)
            write(vdname,'(a2,"_",a4,"_",i4.4)') pl_id(plane), var, nimg
            call h5ltmake_dataset_double_f(gr2_id, vdname, rank, dims, img, error)
            timebuffer(:) = [ t ]
            call h5ltset_attribute_double_f(gr2_id, vdname, "time", timebuffer, bufsize, error)
            call H5Gclose_f(gr2_id, error)
            call H5Gclose_f(gr_id, error)
            call H5Fclose_f(file_id, error)
         endif

         if (allocated(img))  deallocate(img)
      else
         if ((xn > cg%nb .and. xn <= cg%n_b(plane)+cg%nb) .or. (xn == 1 .and. .not. has_dir(plane))) &
              call MPI_Send(send, size(send), MPI_DOUBLE_PRECISION, FIRST, tag, comm, ierr)
      endif

      call MPI_Barrier(comm, ierr) ! We must synchronize everyone before we reuse buffers and variables

      if (allocated(send)) deallocate(send)

   end subroutine write_plot_hdf5

!>
!! \brief Routine to set parameters and dimensions of arrays in restart file
!! \param area_type case name; possibilities:
!!   AT_ALL_B - whole domain with mpi boundaries,
!!   AT_OUT_B - physical domain with outer boundaries,
!    AT_NO_B  - only physical domain without any boundaries
!! \param area grid dimensions in the file
!! \param chnk dimensions of data array dumped by this process
!! \param lleft left limits of data from array to be dumped
!! \param lright right limits of data from array to be dumped
!! \loffs offset in area for this process
!<
   subroutine set_dims_to_write(area_type, area, chnk, lleft, lright, loffs, cg)

      use constants,  only: ndims, AT_ALL_B, AT_OUT_B, AT_NO_B, AT_USER, LO, HI
      use dataio_pub, only: die, warn
      use domain,     only: dom, has_dir, cdd, is_uneven, is_mpi_noncart
      use grid_cont,  only: grid_container
      use mpi,        only: MPI_COMM_NULL
      use types,      only: at_user_settings

      implicit none

      integer(kind=4),                   intent(in)  :: area_type
      integer,         dimension(ndims), intent(out) :: area, lleft, lright, chnk
      integer(kind=8), dimension(ndims), intent(out) :: loffs
      type(grid_container), pointer,     intent(in)  :: cg

      select case (area_type)
         case (AT_ALL_B)                           ! whole domain with mpi boundaries
            if (is_mpi_noncart) call die("[dataio_hdf5:set_dims_to_write] allbnd dump is too hard to implement with noncartesian domain division") !psize, pcoords
            if (cdd%comm3d == MPI_COMM_NULL) call die("[dataio_hdf5:set_dims_to_write] allbnd dump requires cdd%comm3d and cdd%")
            if (is_uneven) call warn("[dataio_hdf5:set_dims_to_write] allbnd dump with uneven domain division")
            chnk(:)   = cg%n_
            area(:)   = dom%n_d(:) + 2 * cg%nb * cdd%psize(:) ! \todo invent something better
            lleft(:)  = 1
            lright(:) = chnk
            loffs(:)  = cg%off(:) + 2 * cg%nb * cdd%pcoords(:) !\todo invent something better
         case (AT_OUT_B)                                   ! physical domain with outer boundaries
            area(:)   = dom%n_t(:)
            lleft(:)  = cg%ijkse(:, LO)
            lright(:) = cg%ijkse(:, HI)
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
            where (loffs(:)>0) loffs(:) = loffs(:) + cg%nb ! the block adjacent to the left boundary are cg%nb cells wider than cg%n_b(:)
         case (AT_NO_B)                                    ! only physical domain without any boundaries
            area(:)   = dom%n_d(:)
            lleft(:)  = cg%ijkse(:, LO)
            lright(:) = cg%ijkse(:, HI)
            chnk(:)   = cg%n_b(:)
            loffs(:)  = cg%off(:)
         case (AT_USER)                                    ! user defined domain (with no reference to simulations domain)
            if (associated(at_user_settings)) then
               call at_user_settings(area, lleft, lright, chnk, loffs)
            else
               call die("[dataio_hdf5:set_dims_to_write] Routine at_user_settings not associated")
            endif
         case default
            call die("[dataio_hdf5:set_dims_to_write] Non-recognized area_type.")
            area(:) = 0 ! suppres compiler warnings
      end select

   end subroutine set_dims_to_write

!>
!! \brief This routine writes restart dump and updates restart counter
!<

   subroutine write_restart_hdf5(debug_res)

      use constants,   only: cwdlen, AT_ALL_B, AT_OUT_B, AT_NO_B, INT4
      use dataio_pub,  only: chdf, nres, set_container_chdf, problem_name, run_id, msg, printio, hdf
      use global,      only: nstep
      use grid,        only: cga
      use grid_cont,   only: cg_list_element, grid_container
      use hdf5,        only: HID_T, H5P_FILE_ACCESS_F, H5F_ACC_TRUNC_F, h5open_f, h5close_f, h5fcreate_f, h5fclose_f, h5pcreate_f, h5pclose_f, h5pset_fapl_mpio_f
      !, H5P_DATASET_XFER_F, h5pset_preserve_f
      use dataio_user, only: problem_write_restart
      use mpi,         only: MPI_CHARACTER
      use mpisetup,    only: comm, info, ierr, master, FIRST

      implicit none

      logical, optional, intent(in) :: debug_res

      integer, parameter    :: extlen = 4
      character(len=extlen) :: file_extension
      character(len=cwdlen) :: filename  !> HDF File name
      integer(HID_T)        :: file_id       !> File identifier
      integer(HID_T)        :: plist_id      !> Property list identifier
      integer(kind=4)       :: area_type
      integer(kind=4)       :: error
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      ! Construct the filename
      if (present(debug_res)) then
         file_extension = '.dmp'
      else
         file_extension = '.res'
      endif

      if (master) then
         write(filename, '(a,a1,a3,a1,i4.4,a4)') trim(problem_name), '_', run_id, '_', nres, file_extension
         write(msg,'(3a)') 'Writing restart ', trim(filename), " ... "
         call printio(msg, .true.)
      endif
      call MPI_Bcast(filename, cwdlen, MPI_CHARACTER, FIRST, comm, ierr)
      call set_container_chdf(nstep)

      ! Set up a new HDF5 file for parallel write
      call h5open_f(error)
      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, error)
      call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)

      ! Write all data in parallel
      call cga%get_root(cgl)
      do while (associated(cgl))
         cg => cgl%cg

         if (associated(problem_write_restart)) call problem_write_restart(file_id, cg)

         if (associated(cg%cs_iso2%arr)) call write_arr_to_restart(file_id, cg%cs_iso2%arr, AT_NO_B, "cs_iso2", cg)
         if (associated(cg%gp%arr))      call write_arr_to_restart(file_id, cg%gp%arr, AT_OUT_B, "gp", cg)

         ! Write fluids
         area_type = AT_NO_B
         if (present(debug_res)) area_type = AT_ALL_B
         if (associated(cg%u%arr)) call write_arr_to_restart(file_id, cg%u%arr, area_type, dname(FLUID), cg)

         ! Write magnetic field
         area_type = AT_OUT_B ! unlike fluids, we need magnetic field boundaries values. Then chunks might be non-uniform
         if (present(debug_res)) area_type = AT_ALL_B
         if (associated(cg%b%arr)) call write_arr_to_restart(file_id, cg%b%arr, area_type, dname(MAG), cg)
         cgl => cgl%nxt
      enddo

!     \todo writing axes using collective I/O takes order of magnitude more than
!        dumping U and B arrays alltogether, since XYZ-axis is not even read
!        back during restart I'm commenting this out. Rewrite or punt.
!      call write_axes_to_restart(file_id)

      ! End of parallel writing (close the HDF5 file stuff)
      call h5pclose_f(plist_id, error)

      ! dump cg
      !call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      !call h5pset_preserve_f(plist_id, .true., error)
      !call write_grid_containter(cg, file_id, plist_id)
      !call h5pclose_f(plist_id, error)
      !
      call h5fclose_f(file_id, error)
      call h5close_f(error)

      ! Write some global variables
      call set_common_attributes(filename, chdf)

      nres = nres + 1_INT4

   end subroutine write_restart_hdf5

   !----------------------------------------------------------------------------------
   ! Write fluid, mag or other variables (4-D and 3-D arrays)
   !
   ! Having both rank-3 array pointer and rank-4 array pointer doesn;t look elegant, but works.
   ! Is there a way to pass only one, "universal" array pointer in Fortran?

   subroutine prep_arr_write(rank, ir, area_type, loffs, chunk_dims, dimsf, file_id, dname, memspace, plist_id, filespace, dset_id, dplist_id, dfilespace)

      use hdf5,       only: HID_T, HSIZE_T, H5T_NATIVE_DOUBLE, &
           &                H5P_DATASET_CREATE_F, H5S_SELECT_SET_F, H5P_DATASET_XFER_F, H5FD_MPIO_COLLECTIVE_F, &
           &                h5screate_simple_f, h5pcreate_f, h5dcreate_f, h5dget_space_f, &
           &                h5pset_chunk_f, h5pset_dxpl_mpio_f, h5sselect_hyperslab_f
      use domain,     only: is_uneven
      use constants,  only: ndims, AT_OUT_B, LONG

      implicit none

      integer, parameter                             :: rank4 = 1 + ndims
      integer(kind=4), intent(in)                    :: rank
      integer, intent(in)                            :: ir
      integer(kind=4), intent(in)                    :: area_type !> no boundaries, only outer boundaries or all boundaries
      integer(HSIZE_T), dimension(rank4), intent(in) :: chunk_dims, dimsf
      integer(kind=8),  dimension(ndims), intent(in) :: loffs
      integer(HID_T), intent(in)                     :: file_id   !> File identifier
      character(len=*), intent(in)                   :: dname

      integer(HID_T), intent(out) :: dset_id    !> Dataset identifier
      integer(HID_T), intent(out) :: filespace  !>
      integer(HID_T), intent(out) :: dfilespace !>
      integer(HID_T), intent(out) :: memspace   !> Dataspace identifier in memory
      integer(HID_T), intent(out) :: plist_id   !>
      integer(HID_T), intent(out) :: dplist_id  !>

      integer(HSIZE_T), dimension(rank4) :: count, offset, stride, block
      integer(kind=4) :: error

      ! Create the file space for the dataset and make it chunked if possible
      call h5screate_simple_f(rank, dimsf(ir:), dfilespace, error)
      call h5pcreate_f(H5P_DATASET_CREATE_F, dplist_id, error)
      if (.not. (is_uneven .or. (area_type == AT_OUT_B))) call h5pset_chunk_f(dplist_id, rank, chunk_dims(ir:), error)
      call h5dcreate_f(file_id, dname, H5T_NATIVE_DOUBLE, dfilespace, dset_id, error, dplist_id)

      ! Each process defines dataset in memory and writes it to the hyperslab in the file.
      call h5dget_space_f(dset_id, filespace, error)

      stride(:) = 1
      count(:)  = 1
      block(:)  = chunk_dims(:)
      offset(:) = [0_LONG, loffs(:)]
      call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset(ir:), count(ir:), error, stride(ir:), block(ir:))

      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

      call h5screate_simple_f(rank, chunk_dims(ir:), memspace, error)
      return
   end subroutine prep_arr_write

   subroutine clean_arr_write(memspace, plist_id, filespace, dset_id, dplist_id, dfilespace)
      use hdf5,       only: HID_T, h5sclose_f, h5pclose_f, h5dclose_f
      implicit none
      integer(HID_T), intent(inout) :: memspace, plist_id, filespace, dset_id, dplist_id, dfilespace
      integer(kind=4) :: error

      call h5sclose_f(memspace, error)
      call h5pclose_f(plist_id, error)
      call h5sclose_f(filespace, error)
      call h5dclose_f(dset_id, error)
      call h5pclose_f(dplist_id, error)
      call h5sclose_f(dfilespace, error)
   end subroutine clean_arr_write

   subroutine write_4darr_to_restart(file_id, pa4d, area_type, dname, cg)

      use constants,  only: xdim, ydim, zdim, ndims
      use dataio_pub, only: die
      use grid_cont,  only: grid_container
      use hdf5,       only: HID_T, HSIZE_T, H5T_NATIVE_DOUBLE, h5dwrite_f

      implicit none

      integer(HID_T), intent(in)                    :: file_id   !> File identifier
      real, pointer, dimension(:,:,:,:), intent(in) :: pa4d      !> 4-D array pointer
      integer(kind=4), intent(in)                   :: area_type !> no boundaries, only outer boundaries or all boundaries
      character(len=*), intent(in)                  :: dname
      type(grid_container), pointer, intent(in)     :: cg

      integer, parameter :: rank4 = 1 + ndims
      integer(HSIZE_T), dimension(rank4) :: dimsf, chunk_dims
      integer(HID_T) :: dset_id               !> Dataset identifier
      integer(HID_T) :: dplist_id, plist_id   !> Property list identifiers
      integer(HID_T) :: dfilespace, filespace !> Dataspace identifiers in file
      integer(HID_T) :: memspace              !> Dataspace identifier in memory

      integer,         dimension(ndims) :: area, lleft, lright, chnk
      integer(kind=8), dimension(ndims) :: loffs

      integer(kind=4) :: rank
      integer(kind=4) :: error

      call set_dims_to_write(area_type, area, chnk, lleft, lright, loffs, cg)

      if (.not.associated(pa4d)) call die("[dataio_hdf5:write_4darr_to_restart] Null pointer given.")
      dimsf      = [size(pa4d,1), area(:)]   ! Dataset dimensions
      chunk_dims = [size(pa4d,1), chnk(:)]   ! Chunks dimensions
      rank = rank4

      call prep_arr_write(rank, 1, area_type, loffs, chunk_dims, dimsf, file_id, dname, memspace, plist_id, filespace, dset_id, dplist_id, dfilespace)

      ! write data
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, pa4d(:, lleft(xdim):lright(xdim), lleft(ydim):lright(ydim), lleft(zdim):lright(zdim)), &
           &          dimsf(:), error, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

      call clean_arr_write(memspace, plist_id, filespace, dset_id, dplist_id, dfilespace)

   end subroutine write_4darr_to_restart

   subroutine write_3darr_to_restart(file_id, pa3d, area_type, dname, cg)

      use constants,  only: xdim, ydim, zdim, ndims
      use dataio_pub, only: die
      use grid_cont,  only: grid_container
      use hdf5,       only: HID_T, HSIZE_T, H5T_NATIVE_DOUBLE, h5dwrite_f

      implicit none

      integer(HID_T), intent(in)                    :: file_id   !> File identifier
      real, pointer, dimension(:,:,:), intent(in)   :: pa3d      !> 3-D array pointer, mutually exclusive with pa4d
      integer(kind=4), intent(in)                   :: area_type !> no boundaries, only outer boundaries or all boundaries
      character(len=*), intent(in)                  :: dname
      type(grid_container), pointer, intent(in)     :: cg

      integer, parameter :: rank4 = 1 + ndims
      integer(HSIZE_T), dimension(rank4) :: dimsf, chunk_dims
      integer(HID_T) :: dset_id               !> Dataset identifier
      integer(HID_T) :: dplist_id, plist_id   !> Property list identifiers
      integer(HID_T) :: dfilespace, filespace !> Dataspace identifiers in file
      integer(HID_T) :: memspace              !> Dataspace identifier in memory

      integer,         dimension(ndims) :: area, lleft, lright, chnk
      integer(kind=8), dimension(ndims) :: loffs

      integer(kind=4) :: rank
      integer(kind=4) :: error

      if (.not. associated(pa3d)) call die("[dataio_hdf5:write_3darr_to_restart] Null pointer given.")

      call set_dims_to_write(area_type, area, chnk, lleft, lright, loffs, cg)

      dimsf = [1, area(:)]      ! Dataset dimensions
      chunk_dims = [1, chnk(:)] ! Chunks dimensions

      rank = ndims

      call prep_arr_write(rank, 2, area_type, loffs, chunk_dims, dimsf, file_id, dname, memspace, plist_id, filespace, dset_id, dplist_id, dfilespace)

      ! write data
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, pa3d(lleft(xdim):lright(xdim), lleft(ydim):lright(ydim), lleft(zdim):lright(zdim)), &
           &          dimsf(2:), error, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

      call clean_arr_write(memspace, plist_id, filespace, dset_id, dplist_id, dfilespace)

   end subroutine write_3darr_to_restart

!!$   !----------------------------------------------------------------------------------
!!$   !
!!$   !  WRITE Axes
!!$   !
!!$   ! AMR: the axes should be associated with fluid and datasets
!!$   subroutine write_axes_to_restart(file_id)
!!$
!!$      use constants,  only: xdim, ydim, zdim, ndims
!!$      use hdf5,       only: HSIZE_T, HID_T, H5T_NATIVE_DOUBLE, H5FD_MPIO_COLLECTIVE_F, H5P_DATASET_CREATE_F, H5P_DATASET_XFER_F, H5S_SELECT_SET_F, &
!!$           &                h5screate_simple_f, h5pcreate_f, h5pset_chunk_f, h5pset_dxpl_mpio_f, h5dcreate_f, h5dget_space_f, h5dwrite_f, &
!!$           &                h5sclose_f, h5pclose_f, h5dclose_f, h5sselect_hyperslab_f
!!$      use domain,     only: dom, is_uneven
!!$      use grid,       only: cg
!!$
!!$      implicit none
!!$
!!$      integer(HID_T), intent(in) :: file_id !> File identifier
!!$
!!$      integer :: dir
!!$      integer :: error
!!$      integer, parameter :: rank=1
!!$      integer(HSIZE_T), dimension(rank) :: offset, count, stride, block, dimsf, chunk_dims
!!$      integer(HID_T) :: dset_id                !> Dataset identifier
!!$      integer(HID_T) :: dplist_id, plist_id    !> Property list identifiers
!!$      integer(HID_T) :: dfilespace, filespace  !> Dataspace identifiers in file
!!$      integer(HID_T) :: memspace               !> Dataspace identifier in memory
!!$      integer, parameter :: asis_n_len = 6
!!$      character(len=asis_n_len) :: dset_axis_n !> Dataspace name
!!$      character(len=ndims), parameter :: axis_n = "XYZ"
!!$
!!$      do dir = xdim, zdim
!!$
!!$         write(dset_axis_n,'(a,"-axis")')axis_n(dir:dir)
!!$         dimsf      = [dom%n_d(dir)] ! Dataset dimensions
!!$         chunk_dims = [cg%n_b(dir)]  ! Chunk dimensions
!!$
!!$         ! Create the file space for the dataset and make it chunked if possible
!!$         call h5screate_simple_f(rank, dimsf, dfilespace, error)
!!$         call h5pcreate_f(H5P_DATASET_CREATE_F, dplist_id, error)
!!$         if (.not. is_uneven) call h5pset_chunk_f(dplist_id, rank, chunk_dims, error)
!!$         call h5dcreate_f(file_id, dset_axis_n, H5T_NATIVE_DOUBLE, dfilespace, dset_id, error, dplist_id)
!!$
!!$         ! It is sufficient that only one CPU writes each piece of axis data.
!!$         ! The other CPUs also need to call at least everything up to h5dcreate_f, and then just close the stuff.
!!$         if (cg%off(mod(dir, ndims)+1) == 0 .and. cg%off(mod(dir+ndims-2, ndims)+1) == 0) then
!!$
!!$            call h5dget_space_f(dset_id, filespace, error)
!!$
!!$            ! Each contributing process defines dataset in memory and writes it to the hyperslab in the file.
!!$            stride(:) = 1
!!$            count(:)  = 1
!!$            block(:)  = chunk_dims(:)
!!$            offset(:) = [ cg%off(dir) ]
!!$            call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error, stride, block)
!!$
!!$            call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
!!$            call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
!!$
!!$            ! Create the memory space for the dataset.
!!$            call h5screate_simple_f(rank, chunk_dims, memspace, error)
!!$            select case (dir)
!!$               case (xdim)
!!$                  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, cg%x(cg%is:cg%ie), dimsf, error, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
!!$               case (ydim)
!!$                  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, cg%y(cg%js:cg%je), dimsf, error, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
!!$               case (zdim)
!!$                  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, cg%z(cg%ks:cg%ke), dimsf, error, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
!!$            end select
!!$            call h5sclose_f(memspace, error)
!!$            call h5pclose_f(plist_id, error)
!!$            call h5sclose_f(filespace, error)
!!$
!!$         endif
!!$
!!$         call h5dclose_f(dset_id, error)
!!$         call h5pclose_f(dplist_id, error)
!!$         call h5sclose_f(dfilespace, error)
!!$
!!$      enddo
!!$
!!$   end subroutine write_axes_to_restart

! This routine reads only interior cells of a pointed array from the restart file.
! Boundary cells are exchanged with the neughbours. Corner boundary cells are not guaranteed to be correct (area_type = AT_ALL_B not implemented yet).
! External boundary cells are not stored in the restart file and thus all of them are lost (area_type = AT_OUT_B not implemented yet).

   subroutine prep_arr_read(rank, ir, loffs, chunk_dims, file_id, dname, memspace, plist_id, filespace, dset_id)

      use hdf5,       only: HID_T, HSIZE_T, H5S_SELECT_SET_F, H5P_DATASET_XFER_F, H5FD_MPIO_COLLECTIVE_F, &
           &                h5dopen_f, h5sget_simple_extent_ndims_f, h5dget_space_f, &
           &                h5pcreate_f, h5pset_dxpl_mpio_f, h5sselect_hyperslab_f, h5screate_simple_f
      use constants,  only: ndims, LONG
      use dataio_pub, only: msg, die

      implicit none

      integer, parameter                             :: rank4 = 1 + ndims
      integer(kind=4), intent(in)                    :: rank
      integer, intent(in)                            :: ir
      integer(HSIZE_T), dimension(rank4), intent(in) :: chunk_dims
      integer(kind=8),  dimension(ndims), intent(in) :: loffs
      integer(HID_T), intent(in)                     :: file_id   !> File identifier
      character(len=*), intent(in)                   :: dname

      integer(HID_T), intent(out) :: dset_id    !> Dataset identifier
      integer(HID_T), intent(out) :: filespace  !>
      integer(HID_T), intent(out) :: memspace   !> Dataspace identifier in memory
      integer(HID_T), intent(out) :: plist_id   !>

      integer(HSIZE_T), dimension(rank4) :: count, offset, stride, block
      integer(kind=4) :: error
      integer(kind=4) :: rankf

      ! Create dataset.and filespace
      call h5dopen_f(file_id, dname, dset_id, error)
      if (error /= 0) then
         write(msg, '(3a)') "[dataio_hdf5:prep_arr_read] Opening dataset '", dname,"' failed."
         call die(msg)
      endif

      call h5dget_space_f(dset_id, filespace, error)
      call h5sget_simple_extent_ndims_f (filespace, rankf, error)
      if (rank /= rankf) then
         write(msg,'(3a,2(i2,a))')"[dataio_hdf5:prep_arr_read] Rank mismatch in array '", dname, "' (", rank, " /= ", rankf, ")"
         call die(msg)
      endif

      ! Select hyperslab in the file.
      stride(:) = 1
      count(:)  = 1
      block(:)  = chunk_dims(:)
      offset(:) = [ 0_LONG, loffs(:) ]
      call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset(ir:), count(ir:), error, stride(ir:), block(ir:))

      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
      call h5screate_simple_f(rank, chunk_dims(ir:), memspace, error)

   end subroutine prep_arr_read

   subroutine clean_arr_read(memspace, plist_id, filespace, dset_id)
      use hdf5,       only: HID_T, h5sclose_f, h5pclose_f, h5dclose_f
      implicit none
      integer(HID_T), intent(inout) :: memspace, plist_id, filespace, dset_id
      integer(kind=4) :: error

      call h5sclose_f(memspace, error)
      call h5pclose_f(plist_id, error)
      call h5sclose_f(filespace, error)
      call h5dclose_f(dset_id, error)

   end subroutine clean_arr_read

   subroutine read_4darr_from_restart(file_id, pa4d, area_type, dname, cg)

      use constants,  only: xdim, ydim, zdim, ndims
      use dataio_pub, only: msg, die
      use grid_cont,  only: grid_container
      use hdf5,       only: HID_T, HSIZE_T, SIZE_T, H5T_NATIVE_DOUBLE, h5dread_f

      implicit none

      integer(HID_T), intent(in)                       :: file_id   ! File identifier
      real, dimension(:,:,:,:), pointer, intent(inout) :: pa4d      ! pointer to (:, 1:cg%nx, 1:cg%ny, 1:cg%nz)-sized array
      integer(kind=4), intent(in)                      :: area_type !> no boundaries, only outer boundaries or all boundaries
      character(len=*), intent(in)                     :: dname
      type(grid_container), pointer, intent(in)        :: cg

      integer(HID_T)        :: dset_id       ! Dataset identifier
      integer(HID_T)        :: plist_id      ! Property list identifier
      integer(HID_T)        :: filespace     ! Dataspace identifier in file
      integer(HID_T)        :: memspace      ! Dataspace identifier in memory

      integer, parameter :: rank4 = 1 + ndims
      integer(HSIZE_T), dimension(rank4) :: dimsf, chunk_dims

      integer,         dimension(ndims) :: area, lleft, lright, chnk
      integer(kind=8), dimension(ndims) :: loffs
      integer :: ir
      integer(kind=4) :: rank
      integer(kind=4) :: error

      call set_dims_to_write(area_type, area, chnk, lleft, lright, loffs, cg)

      dimsf = [size(pa4d,1), area(:)]      ! Dataset dimensions
      chunk_dims = [size(pa4d,1), chnk(:)] ! Chunks dimensions
      rank = rank4
      ir = rank4 - rank + 1 ! 1 for 4-D arrays, 2 for 3-D arrays (to simplify use of count(:), offset(:), stride(:), block(:), dimsf(:) and chunk_dims(:)

      call prep_arr_read(rank, ir, loffs, chunk_dims, file_id, dname, memspace, plist_id, filespace, dset_id)

      ! Read the array
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, pa4d(:, lleft(xdim):lright(xdim), lleft(ydim):lright(ydim), lleft(zdim):lright(zdim)), &
           &         dimsf(ir:), error, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

      if (error /= 0) then
         write(msg, '(3a)') "[dataio_hdf5:read_4darr_from_restart] Reading dataset '", dname,"' failed."
         call die(msg)
      endif

      call clean_arr_read(memspace, plist_id, filespace, dset_id)
      ! rank-4 arrays (cg%u%arr(:,:,:,:) and b(:,:,:,:)) have their own guardcell-exchange routines, which can also be called here

   end subroutine read_4darr_from_restart

   subroutine read_3darr_from_restart(file_id, pa3d, area_type, dname, cg)

      use constants,  only: xdim, ydim, zdim, ndims
      use dataio_pub, only: msg, die
      use grid,       only: arr3d_boundaries
      use grid_cont,  only: grid_container
      use hdf5,       only: HID_T, HSIZE_T, SIZE_T, h5dread_f, H5T_NATIVE_DOUBLE

      implicit none

      integer(HID_T), intent(in)                       :: file_id   ! File identifier
      real, dimension(:,:,:), pointer, intent(inout)   :: pa3d      ! pointer to (1:cg%nx, 1:cg%ny, 1:cg%nz)-sized array
      integer(kind=4), intent(in)                      :: area_type !> no boundaries, only outer boundaries or all boundaries
      character(len=*), intent(in)                     :: dname
      type(grid_container), pointer, intent(in)        :: cg

      integer(HID_T)        :: dset_id       ! Dataset identifier
      integer(HID_T)        :: plist_id      ! Property list identifier
      integer(HID_T)        :: filespace     ! Dataspace identifier in file
      integer(HID_T)        :: memspace      ! Dataspace identifier in memory

      integer, parameter :: rank4 = 1 + ndims
      integer(HSIZE_T), dimension(rank4) :: dimsf, chunk_dims

      integer,         dimension(ndims) :: area, lleft, lright, chnk
      integer(kind=8), dimension(ndims) :: loffs
      integer :: ir
      integer(kind=4) :: rank
      integer(kind=4) :: error

      call set_dims_to_write(area_type, area, chnk, lleft, lright, loffs, cg)

      dimsf = [1, area(:)]      ! Dataset dimensions
      chunk_dims = [1, chnk(:)] ! Chunks dimensions
      if (.not. associated(pa3d)) call die("[dataio_hdf5:read_3darr_from_restart] Null pointer given.")
      rank = ndims
      ir = rank4 - rank + 1 ! 1 for 4-D arrays, 2 for 3-D arrays (to simplify use of count(:), offset(:), stride(:), block(:), dimsf(:) and chunk_dims(:)

      call prep_arr_read(rank, ir, loffs, chunk_dims, file_id, dname, memspace, plist_id, filespace, dset_id)

      ! Read the array
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, pa3d(lleft(xdim):lright(xdim), lleft(ydim):lright(ydim), lleft(zdim):lright(zdim)), &
           &         dimsf(ir:), error, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

      if (error /= 0) then
         write(msg, '(3a)') "[dataio_hdf5:read_3darr_from_restart] Reading dataset '", dname,"' failed."
         call die(msg)
      endif

      call clean_arr_read(memspace, plist_id, filespace, dset_id)

      ! Originally the pa3d array was written with the guardcells. The internal guardcells will be exchanged but the external ones are lost.
      call arr3d_boundaries(pa3d, area_type, dname)

   end subroutine read_3darr_from_restart

   subroutine read_restart_hdf5(chdf)

      use constants,   only: cwdlen, cbuff_len, domlen, idlen, xdim, ydim, zdim, AT_NO_B, AT_OUT_B, LO, HI
      use dataio_pub,  only: msg, printio, warn, die, require_init_prob, problem_name, run_id, piernik_hdf5_version, hdf
      use dataio_user, only: problem_read_restart
      use domain,      only: dom, has_dir
      use fluidindex,  only: flind
      use func,        only: fix_string
      use global,      only: magic_mass, t, dt
      use grid,        only: cga
      use grid_cont,   only: cg_list_element, grid_container
      use hdf5,        only: HID_T, SIZE_T, H5P_FILE_ACCESS_F, H5F_ACC_RDONLY_F, &
           &                 h5open_f, h5pcreate_f, h5pset_fapl_mpio_f, h5fopen_f, h5pclose_f, h5fclose_f, h5close_f
      use h5lt,        only: h5ltget_attribute_double_f, h5ltget_attribute_int_f, h5ltget_attribute_string_f
      use mpi,         only: MPI_CHARACTER, MPI_INTEGER, MPI_DOUBLE_PRECISION
      use mpisetup,    only: comm, ierr, info, comm, master, FIRST

      implicit none

      type(hdf)             :: chdf
      integer               :: nu
      character(len=cwdlen) :: filename  ! File name

      integer(HID_T)        :: file_id       ! File identifier
      integer(HID_T)        :: plist_id      ! Property list identifier
      integer(SIZE_T)       :: bufsize

      integer(kind=4)       :: error
      logical               :: file_exist

      real, dimension(1)    :: rbuf
      integer(kind=4), dimension(1) :: ibuf

      real                  :: restart_hdf5_version
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      nu = flind%all

      if (master) then
         write(filename,'(a,a1,a3,a1,i4.4,a4)') trim(problem_name),'_', run_id,'_', chdf%nres,'.res'
         write(msg, '(2a)') 'Reading restart file: ', trim(filename)
         call printio(msg)
      endif
      call MPI_Bcast(filename, cwdlen, MPI_CHARACTER, FIRST, comm, ierr)

      inquire(file = filename, exist = file_exist)
      if (.not. file_exist) then
         write(msg,'(3a)') '[dataio_hdf5:read_restart_hdf5]: Restart file: ', trim(filename),' does not exist'
         call die(msg)
      endif

      call h5open_f(error)
      if (master) then
         call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, error)

         call h5ltget_attribute_double_f(file_id,"/","piernik", rbuf, error)
         if (error /= 0) call die("[dataio_hdf5:read_restart_hdf5] Cannot read 'piernik' attribute from the restart file. The file may be either damaged or incompatible")
         if (rbuf(1) > piernik_hdf5_version) then
            write(msg,'(2(a,f5.2))')"[dataio_hdf5:read_restart_hdf5] Cannot read future versions of the restart file: ", rbuf(1)," > ", piernik_hdf5_version
            call die(msg)
         else if (int(rbuf(1)) < int(piernik_hdf5_version)) then
            write(msg,'(2(a,f5.2))')"[dataio_hdf5:read_restart_hdf5] The restart file is too ancient. It is unlikely that it could work correctly: ", rbuf(1)," << ", piernik_hdf5_version
            call die(msg)
         else if (rbuf(1) < piernik_hdf5_version) then
            write(msg,'(2(a,f5.2))')"[dataio_hdf5:read_restart_hdf5] Old versions of the restart file may not always work fully correctly: ", rbuf(1)," < ", piernik_hdf5_version
            call warn(msg)
         endif

         call h5ltget_attribute_int_f(file_id,"/","nxd", ibuf, error)
         if (ibuf(1) /= dom%n_d(xdim) .or. error /= 0) call die("[dataio_hdf5:read_restart_hdf5] nxd does not match")
         if (has_dir(xdim)) then
            call h5ltget_attribute_double_f(file_id,"/","xmin", rbuf, error)
            if (rbuf(1) /= dom%edge(xdim, LO) .or. error /= 0) call die("[dataio_hdf5:read_restart_hdf5] xmin does not match")
            call h5ltget_attribute_double_f(file_id,"/","xmax", rbuf, error)
            if (rbuf(1) /= dom%edge(xdim, HI) .or. error /= 0) call die("[dataio_hdf5:read_restart_hdf5] xmax does not match")
         endif

         call h5ltget_attribute_int_f(file_id,"/","nyd", ibuf, error)
         if (ibuf(1) /= dom%n_d(ydim) .or. error /= 0) call die("[dataio_hdf5:read_restart_hdf5] nyd does not match")
         if (has_dir(ydim)) then
            call h5ltget_attribute_double_f(file_id,"/","ymin", rbuf, error)
            if (rbuf(1) /= dom%edge(ydim, LO) .or. error /= 0) call die("[dataio_hdf5:read_restart_hdf5] ymin does not match")
            call h5ltget_attribute_double_f(file_id,"/","ymax", rbuf, error)
            if (rbuf(1) /= dom%edge(ydim, HI) .or. error /= 0) call die("[dataio_hdf5:read_restart_hdf5] ymax does not match")
         endif

         call h5ltget_attribute_int_f(file_id,"/","nzd", ibuf, error)
         if (ibuf(1) /= dom%n_d(zdim) .or. error /= 0) call die("[dataio_hdf5:read_restart_hdf5] nzd does not match")
         if (has_dir(zdim)) then
            call h5ltget_attribute_double_f(file_id,"/","zmin", rbuf, error)
            if (rbuf(1) /= dom%edge(zdim, LO) .or. error /= 0) call die("[dataio_hdf5:read_restart_hdf5] zmin does not match")
            call h5ltget_attribute_double_f(file_id,"/","zmax", rbuf, error)
            if (rbuf(1) /= dom%edge(zdim, HI) .or. error /= 0) call die("[dataio_hdf5:read_restart_hdf5] zmax does not match")
         endif

         call h5fclose_f(file_id, error)
      endif

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, error)

      call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, error, access_prp = plist_id)
      call h5pclose_f(plist_id, error)

      call cga%get_root(cgl)
      do while (associated(cgl))
         cg => cgl%cg

         if (associated(problem_read_restart)) call problem_read_restart(file_id, cg)

         if (associated(cg%cs_iso2%arr)) call read_arr_from_restart(file_id, cg%cs_iso2%arr, AT_NO_B, "cs_iso2", cg)
         if (associated(cg%gp%arr))      call read_arr_from_restart(file_id, cg%gp%arr, AT_OUT_B, "gp", cg)

         !  READ FLUID VARIABLES
         if (associated(cg%u%arr)) call read_arr_from_restart(file_id, cg%u%arr, AT_NO_B, dname(FLUID), cg)

         !  READ MAG VARIABLES
         if (associated(cg%b%arr)) call read_arr_from_restart(file_id, cg%b%arr, AT_OUT_B, dname(MAG), cg)
         cgl => cgl%nxt
      enddo

      call h5fclose_f(file_id, error)

      if (master) then
         call h5fopen_f (filename, H5F_ACC_RDONLY_F, file_id, error)
         bufsize = 1
         call h5ltget_attribute_double_f(file_id,"/","piernik", rbuf, error)
         restart_hdf5_version = rbuf(1)
         call h5ltget_attribute_double_f(file_id,"/","time", rbuf, error)
         t = rbuf(1)
         call h5ltget_attribute_double_f(file_id,"/","timestep", rbuf, error)
         dt = rbuf(1)
         call h5ltget_attribute_double_f(file_id,"/","magic_mass", rbuf, error)
         magic_mass = rbuf(1)
         call h5ltget_attribute_int_f(file_id,"/","nstep", ibuf, error)
         chdf%nstep = ibuf(1)
         call h5ltget_attribute_int_f(file_id,"/","nres", ibuf, error)
         chdf%nres = ibuf(1)
         call h5ltget_attribute_int_f(file_id,"/","nhdf", ibuf, error)
         chdf%nhdf = ibuf(1)
         call h5ltget_attribute_int_f(file_id,"/","step_res", ibuf, error)
         chdf%step_res = ibuf(1)
         call h5ltget_attribute_int_f(file_id,"/","step_hdf", ibuf, error)
         chdf%step_hdf = ibuf(1)
         call h5ltget_attribute_double_f(file_id,"/","next_t_tsl", rbuf, error)
         chdf%next_t_tsl = rbuf(1)
         call h5ltget_attribute_double_f(file_id,"/","next_t_log", rbuf, error)
         chdf%next_t_log = rbuf(1)
         call h5ltget_attribute_double_f(file_id,"/","last_hdf_time", rbuf, error)
         chdf%last_hdf_time = rbuf(1)

         call h5ltget_attribute_string_f(file_id,"/","problem_name", problem_name, error)
         call h5ltget_attribute_string_f(file_id,"/","domain", chdf%domain_dump, error)
         call h5ltget_attribute_string_f(file_id,"/","run_id", chdf%new_id, error)

         if (restart_hdf5_version > 1.11) then
            call h5ltget_attribute_int_f(file_id,"/","require_init_prob", ibuf, error)
            require_init_prob = ibuf(1)
         endif

         problem_name = fix_string(problem_name)   !> \deprecated BEWARE: >=HDF5-1.8.4 has weird issues with strings
         chdf%new_id  = fix_string(chdf%new_id)    !> \deprecated   this bit hacks it around
         chdf%domain_dump  = fix_string(chdf%domain_dump)

         call h5fclose_f(file_id, error)

         write(msg,'(2a)') 'Done reading restart file: ', trim(filename)
         call printio(msg)
      endif
      call h5close_f(error)

      call MPI_Bcast(restart_hdf5_version,    1, MPI_DOUBLE_PRECISION, FIRST, comm, ierr)

      call MPI_Bcast(chdf%nstep,    1, MPI_INTEGER, FIRST, comm, ierr)
      call MPI_Bcast(chdf%nres,     1, MPI_INTEGER, FIRST, comm, ierr)
      call MPI_Bcast(chdf%nhdf,     1, MPI_INTEGER, FIRST, comm, ierr)
      call MPI_Bcast(chdf%step_res, 1, MPI_INTEGER, FIRST, comm, ierr)
      call MPI_Bcast(chdf%step_hdf, 1, MPI_INTEGER, FIRST, comm, ierr)
      if (restart_hdf5_version > 1.11) call MPI_Bcast(require_init_prob, 1, MPI_INTEGER, FIRST, comm, ierr)

      call MPI_Bcast(chdf%next_t_tsl,    1, MPI_DOUBLE_PRECISION, FIRST, comm, ierr)
      call MPI_Bcast(chdf%next_t_log,    1, MPI_DOUBLE_PRECISION, FIRST, comm, ierr)
      call MPI_Bcast(chdf%last_hdf_time, 1, MPI_DOUBLE_PRECISION, FIRST, comm, ierr)
      call MPI_Bcast(t,                  1, MPI_DOUBLE_PRECISION, FIRST, comm, ierr)
      call MPI_Bcast(dt,                 1, MPI_DOUBLE_PRECISION, FIRST, comm, ierr)

      call MPI_Bcast(problem_name, cbuff_len, MPI_CHARACTER, FIRST, comm, ierr)
      call MPI_Bcast(chdf%domain_dump,domlen, MPI_CHARACTER, FIRST, comm, ierr)
      call MPI_Bcast(chdf%new_id,  idlen,     MPI_CHARACTER, FIRST, comm, ierr)

   end subroutine read_restart_hdf5
!
! ------------------------------------------------------------------------------------
!
   subroutine write_hdf5(chdf)

      use constants,   only: cwdlen, INT4
      use dataio_pub,  only: printio, msg, die, nhdf, problem_name, run_id, hdf
      use dataio_user, only: user_vars_hdf5
      use grid,        only: cga
      use grid_cont,   only: cg_list_element, grid_container
      use hdf5,        only: HID_T, H5F_ACC_TRUNC_F, H5P_FILE_ACCESS_F, H5P_DEFAULT_F, &
           &                 h5open_f, h5close_f, h5fcreate_f, h5fclose_f, h5pcreate_f, h5pclose_f, h5pset_fapl_mpio_f
      use mpisetup,    only: comm, ierr, info, master, FIRST
      use mpi,         only: MPI_CHARACTER
#ifdef NEW_HDF5
      use list_hdf5,   only: iterate_lhdf5
#endif /* NEW_HDF5 */
      use list_hdf5,   only: write_arr

      implicit none

      type(hdf), intent(in)   :: chdf
      integer(HID_T)          :: file_id       ! File identifier
      integer(HID_T)          :: plist_id      ! Property list identifier
      integer                 :: ierrh, i
      integer(kind=4)         :: error
      logical                 :: ok_var
      character(len=cwdlen)   :: fname
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg
      real(kind=4), allocatable :: data (:,:,:)  ! Data to write

      ! Initialize HDF5 library and Fortran interfaces.
      !
      if (master) then
         write(fname, '(a,a1,a3,a1,i4.4,a3)') trim(problem_name),"_", trim(run_id),"_", chdf%nhdf,".h5"
         write(msg,'(3a)') 'Writing datafile ', trim(fname), " ... "
         call printio(msg, .true.)
      endif
      call MPI_Bcast(fname, cwdlen, MPI_CHARACTER, FIRST, comm, ierr)

      call h5open_f(error)
      !
      ! Setup file access property list with parallel I/O access.
      !
      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
      call h5pset_fapl_mpio_f(plist_id, comm, info, error)
      !
      ! Create the file collectively.
      !
      call h5fcreate_f(trim(fname), H5F_ACC_TRUNC_F, file_id, error, creation_prp = H5P_DEFAULT_F, access_prp = plist_id)
      call h5pclose_f(plist_id, error)

      call cga%get_root(cgl)
      do while (associated(cgl))
         cg => cgl%cg

         if (.not.allocated(data)) allocate(data(cg%nxb, cg%nyb, cg%nzb))
         do i = 1, nhdf_vars
            ierrh = 0; ok_var = .false.
            call common_vars_hdf5(hdf_vars(i), data, ierrh, cg)
            if (associated(user_vars_hdf5) .and. ierrh /= 0) call user_vars_hdf5(hdf_vars(i), data, ierrh, cg)
            if (ierrh>=0) ok_var = .true.
            if (.not.ok_var) then
               write(msg,'(3a)') "[dataio_hdf5:write_hdf5]: ", hdf_vars(i)," is not defined in common_vars_hdf5, neither in user_vars_hdf5."
               call die(msg)
            endif
            call write_arr(data, hdf_vars(i), file_id)
         enddo
         if (allocated(data)) deallocate(data)

         cgl => cgl%nxt
      enddo

#ifdef NEW_HDF5
      call iterate_lhdf5(file_id)
#endif /* NEW_HDF5 */

      !
      ! Close the property list.
      !
      call h5fclose_f(file_id, error)
      call h5close_f(error)

      call set_common_attributes(fname, chdf)

      nhdf = nhdf + 1_INT4

   end subroutine write_hdf5

!>
!! \brief This routine writes all attributes that are common to restart and output files.
!! \details Other common elements may also be moved here.
!<
   subroutine set_common_attributes(filename, chdf)

      use constants,   only: cbuff_len, xdim, ydim, zdim, INT4
      use dataio_pub,  only: msg, printio, require_init_prob, piernik_hdf5_version, problem_name, run_id, hdf
      use dataio_user, only: additional_attrs
      use domain,      only: dom
      use global,      only: magic_mass, t, dt, local_magic_mass
      use grid,        only: cga
      use hdf5,        only: HID_T, SIZE_T, HSIZE_T, H5F_ACC_RDWR_F, H5T_NATIVE_CHARACTER, H5Z_FILTER_DEFLATE_F, H5P_DATASET_CREATE_F, &
           &                 h5open_f, h5fopen_f, h5fclose_f, H5Zfilter_avail_f, H5Pcreate_f, H5Pset_deflate_f, H5Pset_chunk_f, &
           &                 h5tcopy_f, h5tset_size_f, h5screate_simple_f, H5Dcreate_f, H5Dwrite_f, H5Dclose_f, H5Sclose_f, H5Tclose_f, H5Pclose_f, h5close_f
      use h5lt,        only: h5ltset_attribute_double_f, h5ltset_attribute_int_f, h5ltset_attribute_string_f
      use mpi,         only: MPI_DOUBLE_PRECISION, MPI_SUM
      use mpisetup,    only: slave, comm, ierr, FIRST
      use version,     only: env, nenv

      implicit none

      character(len=*), intent(in) :: filename  !> HDF File name
      type(hdf), intent(in)        :: chdf

      integer(HID_T)                 :: file_id       !> File identifier
      integer(HID_T)                 :: type_id, dspace_id, dset_id, prp_id
      integer(HSIZE_T), dimension(1) :: dimstr
      logical(kind=4)                :: Z_avail
      integer                        :: fe, i
      integer(SIZE_T)                :: bufsize, maxlen
      integer(kind=4)                :: error
      real                           :: magic_mass0
      integer, parameter             :: buf_len = 50
      integer(kind=4), dimension(buf_len) :: ibuffer
      real,    dimension(buf_len)    :: rbuffer
      character(len=cbuff_len), dimension(buf_len) :: ibuffer_name = ''
      character(len=cbuff_len), dimension(buf_len) :: rbuffer_name = ''

      call MPI_Reduce(local_magic_mass, magic_mass0, 1, MPI_DOUBLE_PRECISION, MPI_SUM, FIRST, comm, ierr)
      magic_mass       = magic_mass + magic_mass0
      local_magic_mass = 0.0

      if (slave) return ! This data need not be written in parallel.

      !! The rr1 marks critical attributes that are read by read_restart_hdf5 and compared against value read from the problem.par file.
      !! The rr2 marks runtime values that are read by read_restart_hdf5 and assigned to something in the code.
      !> \todo Set up an universal table(s) of attribute names for use by both set_common_attributes and read_restart_hdf5.
      !! Provide indices for critical attributes (rr1) and for runtime attributes (rr2).
      !<

      call h5open_f(error)
      call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)

      rbuffer(1)   = t                       ; rbuffer_name(1)   = "time" !rr2
      rbuffer(2)   = dt                      ; rbuffer_name(2)   = "timestep" !rr2
      rbuffer(3)   = chdf%last_hdf_time      ; rbuffer_name(3)   = "last_hdf_time" !rr2
      rbuffer(4:5) = dom%edge(xdim, :)       ; rbuffer_name(4:5) = [ "xmin", "xmax" ] !rr1
      rbuffer(6:7) = dom%edge(ydim, :)       ; rbuffer_name(6:7) = [ "ymin", "ymax" ] !rr1
      rbuffer(8:9) = dom%edge(zdim, :)       ; rbuffer_name(8:9) = [ "zmin", "zmax" ] !rr1
      rbuffer(10)  = piernik_hdf5_version    ; rbuffer_name(10)  = "piernik" !rr1, rr2
      rbuffer(11)  = magic_mass              ; rbuffer_name(11)  = "magic_mass" !rr2
      rbuffer(12)  = chdf%next_t_tsl         ; rbuffer_name(12)  = "next_t_tsl" !rr2
      rbuffer(13)  = chdf%next_t_log         ; rbuffer_name(13)  = "next_t_log" !rr2

      ibuffer(1)   = chdf%nstep              ; ibuffer_name(1)   = "nstep" !rr2
      ibuffer(2)   = chdf%nres+1_INT4        ; ibuffer_name(2)   = "nres" !rr2
      ibuffer(3)   = chdf%nhdf               ; ibuffer_name(3)   = "nhdf" !rr2
      ibuffer(4)   = chdf%nstep              ; ibuffer_name(4)   = "step_res" !rr2
      ibuffer(5)   = chdf%step_hdf           ; ibuffer_name(5)   = "step_hdf" !rr2
      ibuffer(6:8) = dom%n_d(:)              ; ibuffer_name(6:8) = [ "nxd", "nyd", "nzd" ] !rr1
      ibuffer(9)   = cga%cg_all(1)%nb        ; ibuffer_name(9)   = "nb"                             ! BEWARE: assuming cga%cg_all(:)%nb equal everywhere
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
      maxlen = int(maxval(len_trim(parfile(:parfilelines))), kind=4)
      dimstr = [parfilelines]
      call H5Zfilter_avail_f(H5Z_FILTER_DEFLATE_F, Z_avail, error)
      ! call H5Zget_filter_info_f ! everything should be always fine for gzip
      call H5Pcreate_f(H5P_DATASET_CREATE_F, prp_id, error)
      if (Z_avail) then
         call H5Pset_deflate_f(prp_id, 9_INT4, error)
         call H5Pset_chunk_f(prp_id, 1_INT4, dimstr, error)
      endif
      call H5Tcopy_f(H5T_NATIVE_CHARACTER, type_id, error)
      call H5Tset_size_f(type_id, maxlen, error)
      call H5Screate_simple_f(1_INT4, dimstr, dspace_id, error)
      call H5Dcreate_f(file_id, "problem.par", type_id,  dspace_id, dset_id, error, dcpl_id = prp_id)
      call H5Dwrite_f(dset_id, type_id, parfile(:)(:maxlen), dimstr, error)
      call H5Dclose_f(dset_id, error)
      call H5Sclose_f(dspace_id, error)

      ! Store a compressed copy of the piernik.def file and Id lines from source files.
      ! We recycle type_id and prp_id, so we don't close them yet.
      maxlen = int(maxval(len_trim(env(:nenv))), kind=4)
      dimstr = [nenv]
      if (Z_avail) call H5Pset_chunk_f(prp_id, 1_INT4, dimstr, error)
      call H5Tset_size_f(type_id, maxlen, error)
      call H5Screate_simple_f(1_INT4, dimstr, dspace_id, error)
      call H5Dcreate_f(file_id, "env", type_id,  dspace_id, dset_id, error, dcpl_id = prp_id)
      call H5Dwrite_f(dset_id, type_id, env(:)(:maxlen), dimstr, error)
      call H5Dclose_f(dset_id, error)
      call H5Sclose_f(dspace_id, error)
      call H5Tclose_f(type_id, error)
      call H5Pclose_f(prp_id, error)

      !> \todo store full domain decomposition for all procs here [cg%n_b(:), cg%off(:)] @(FIRST:LAST) ! MPI_Gather them?

      fe = len_trim(problem_name)
      call h5ltset_attribute_string_f(file_id, "/", "problem_name", problem_name(1:fe), error) !rr2
      fe = len_trim(chdf%domain_dump)
      call h5ltset_attribute_string_f(file_id, "/", "domain", chdf%domain_dump(1:fe), error) !rr2
      fe = len_trim(run_id)
      call h5ltset_attribute_string_f(file_id, "/", "run_id", run_id(1:fe), error) !rr2

      if (associated(additional_attrs)) call additional_attrs(file_id)

      call h5fclose_f(file_id, error)
      call h5close_f(error)

      write(msg,'(a)') 'done'
      call printio(msg)

      ! only master process exits here

   end subroutine set_common_attributes

   subroutine write_grid_containter(cg, file_id, plist_id)

      use constants, only: xdim, ydim, zdim, ndims, LO, HI, INT4
      use grid_cont, only: grid_container
      use hdf5,      only: HID_T, SIZE_T, HSIZE_T, H5T_NATIVE_INTEGER, H5T_STD_I8LE, H5T_NATIVE_DOUBLE, H5T_COMPOUND_F, &
           &               h5screate_simple_f, h5tarray_create_f, h5tget_size_f, h5tcreate_f, h5tinsert_f, h5dwrite_f, h5sclose_f, h5tclose_f, h5dclose_f, h5dcreate_f

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer(HID_T), intent(in)                :: file_id, plist_id

      integer, parameter :: dsetnamelen = 10
      integer(SIZE_T), parameter :: n_int4 = 19, n_r8 = 14, n_nxarr_r8 = 4, n_nyarr_r8 = 4, n_nzarr_r8 = 4, &
         & n_ndims_r8 = 2, n_ndims_i4 =1, n_ndims_i8 = 1, n_ndims_lohi_i4 = 2

      integer(SIZE_T) :: n_arr3d_r8, n_ndims_arr4d_r8, n_u_arr4d_r8, n_stub
      integer :: total_no

      integer(SIZE_T) :: int4_ts, r8_ts, nxarr_r8_ts, nyarr_r8_ts, nzarr_r8_ts, arr3d_r8_ts, ndims_r8_ts, &
         & ndims_i4_ts, ndims_i8_ts, ndims_lohi_i4_ts, ndims_arr4d_r8_ts, u_arr4d_r8_ts, type_size, offset
      integer(HID_T)  :: ndims_r8_t, ndims_i4_t, ndims_i8_t, nxarr_r8_t, nyarr_r8_t, nzarr_r8_t, arr3d_r8_t, ndims_lohi_i4_t, &
         & ndims_arr4d_r8_t, u_arr4d_r8_t, dtype_id, dspace_id, dset_id
      integer(HSIZE_T),  dimension(1) :: dims
      integer(HID_T),    dimension(:), allocatable :: types, dmem_id
      integer(SIZE_T),   dimension(:), allocatable :: types_sizes
      character(len=dsetnamelen), dimension(:), allocatable :: types_names
      character(len=dsetnamelen) :: dset_name

      integer(kind=4) :: error
      integer :: i

      dims = 1
      call h5screate_simple_f(1_INT4, dims, dspace_id, error)

      call h5tarray_create_f(H5T_NATIVE_DOUBLE,  1_INT4, [integer(HSIZE_T):: ndims],              ndims_r8_t, error)
      call h5tarray_create_f(H5T_NATIVE_INTEGER, 1_INT4, [integer(HSIZE_T):: ndims],              ndims_i4_t, error)
      call h5tarray_create_f(H5T_NATIVE_INTEGER, 2_INT4, [integer(HSIZE_T):: ndims, HI-LO+1],     ndims_lohi_i4_t, error)
      call h5tarray_create_f(H5T_STD_I8LE,       1_INT4, [integer(HSIZE_T):: ndims],              ndims_i8_t, error)
      call h5tarray_create_f(H5T_NATIVE_DOUBLE,  1_INT4, [integer(HSIZE_T):: cg%n_(xdim)],        nxarr_r8_t, error)
      call h5tarray_create_f(H5T_NATIVE_DOUBLE,  1_INT4, [integer(HSIZE_T):: cg%n_(ydim)],        nyarr_r8_t, error)
      call h5tarray_create_f(H5T_NATIVE_DOUBLE,  1_INT4, [integer(HSIZE_T):: cg%n_(zdim)],        nzarr_r8_t, error)
      call h5tarray_create_f(H5T_NATIVE_DOUBLE,  3_INT4, [integer(HSIZE_T):: cg%n_(:)   ],        arr3d_r8_t, error)
      call h5tarray_create_f(H5T_NATIVE_DOUBLE,  4_INT4, [integer(HSIZE_T):: ndims, cg%n_(:)],    ndims_arr4d_r8_t, error)
      call h5tarray_create_f(H5T_NATIVE_DOUBLE,  4_INT4, [integer(HSIZE_T):: size(cg%u%arr,1), cg%n_(:)], u_arr4d_r8_t, error)

      n_arr3d_r8 = 10  ! gc_{x,y,z}dim
      n_stub     = 0
      if (associated(cg%cs_iso2%arr))  n_stub = n_stub + 1_INT4
      if (associated(cg%wa%arr))       n_stub = n_stub + 1_INT4
      if (associated(cg%gpot%arr))     n_stub = n_stub + 1_INT4
      if (associated(cg%hgpot%arr))    n_stub = n_stub + 1_INT4
      if (associated(cg%gp%arr))       n_stub = n_stub + 1_INT4
      if (associated(cg%sgp%arr))      n_stub = n_stub + 1_INT4
      if (associated(cg%sgpm%arr))     n_stub = n_stub + 1_INT4
      n_arr3d_r8 = n_arr3d_r8 - n_stub

      n_ndims_arr4d_r8 = 1  ! b
      if (associated(cg%b0%arr)) then
         n_ndims_arr4d_r8 = n_ndims_arr4d_r8 + 1_INT4
      else
         n_stub = n_stub + 1_INT4
      endif

      n_u_arr4d_r8 = 2 ! u,uh
      if (associated(cg%u0%arr)) then
         n_u_arr4d_r8 = n_u_arr4d_r8 + 1_INT4
      else
         n_stub = n_stub + 1_INT4
      endif

      call h5tget_size_f(H5T_NATIVE_INTEGER, int4_ts,error)
      call h5tget_size_f(H5T_NATIVE_DOUBLE,  r8_ts, error)
      call h5tget_size_f(ndims_r8_t,         ndims_r8_ts, error)
      call h5tget_size_f(ndims_i4_t,         ndims_i4_ts, error)
      call h5tget_size_f(ndims_i8_t,         ndims_i8_ts, error)
      call h5tget_size_f(ndims_lohi_i4_t,    ndims_lohi_i4_ts, error)
      call h5tget_size_f(nxarr_r8_t,         nxarr_r8_ts, error)
      call h5tget_size_f(nyarr_r8_t,         nyarr_r8_ts, error)
      call h5tget_size_f(nzarr_r8_t,         nzarr_r8_ts, error)
      call h5tget_size_f(arr3d_r8_t,         arr3d_r8_ts, error)
      call h5tget_size_f(ndims_arr4d_r8_t,   ndims_arr4d_r8_ts, error)
      call h5tget_size_f(u_arr4d_r8_t,       u_arr4d_r8_ts, error)

      type_size = (n_int4+n_stub)*int4_ts + n_r8*r8_ts + n_ndims_r8*ndims_r8_ts + n_ndims_i4*ndims_i4_ts + n_ndims_i8*ndims_i8_ts + n_ndims_lohi_i4*ndims_lohi_i4_ts &
            & + n_nxarr_r8*nxarr_r8_ts + n_nyarr_r8*nyarr_r8_ts + n_nzarr_r8*nzarr_r8_ts + n_arr3d_r8*arr3d_r8_ts + n_ndims_arr4d_r8*ndims_arr4d_r8_ts + n_u_arr4d_r8*u_arr4d_r8_ts

      call h5tcreate_f(H5T_COMPOUND_F, type_size, dtype_id, error)

      total_no =  int(n_int4 + n_r8 + n_ndims_r8 + n_ndims_i4  + n_ndims_i8 + n_ndims_lohi_i4 + n_nxarr_r8 + n_nyarr_r8 + &
         n_nzarr_r8 + n_arr3d_r8 + n_ndims_arr4d_r8 + n_u_arr4d_r8 + n_stub, kind(total_no))

      allocate(types(total_no), types_sizes(total_no), types_names(total_no), dmem_id(total_no))

      types(1)  = H5T_NATIVE_INTEGER;  types_sizes(1)  = int4_ts;        types_names(1)  = "nx"
      types(2)  = H5T_NATIVE_INTEGER;  types_sizes(2)  = int4_ts;        types_names(2)  = "ny"
      types(3)  = H5T_NATIVE_INTEGER;  types_sizes(3)  = int4_ts;        types_names(3)  = "nz"
      types(4)  = H5T_NATIVE_INTEGER;  types_sizes(4)  = int4_ts;        types_names(4)  = "nxb"
      types(5)  = H5T_NATIVE_INTEGER;  types_sizes(5)  = int4_ts;        types_names(5)  = "nyb"
      types(6)  = H5T_NATIVE_INTEGER;  types_sizes(6)  = int4_ts;        types_names(6)  = "nzb"
      types(7)  = H5T_NATIVE_INTEGER;  types_sizes(7)  = int4_ts;        types_names(7)  = "is"
      types(8)  = H5T_NATIVE_INTEGER;  types_sizes(8)  = int4_ts;        types_names(8)  = "ie"
      types(9)  = H5T_NATIVE_INTEGER;  types_sizes(9)  = int4_ts;        types_names(9)  = "js"
      types(10) = H5T_NATIVE_INTEGER;  types_sizes(10) = int4_ts;        types_names(10) = "je"
      types(11) = H5T_NATIVE_INTEGER;  types_sizes(11) = int4_ts;        types_names(11) = "ks"
      types(12) = H5T_NATIVE_INTEGER;  types_sizes(12) = int4_ts;        types_names(12) = "ke"
      types(13) = H5T_NATIVE_INTEGER;  types_sizes(13) = int4_ts;        types_names(13) = "maxxyz"
      types(14) = H5T_NATIVE_INTEGER;  types_sizes(14) = int4_ts;        types_names(14) = "isb"
      types(15) = H5T_NATIVE_INTEGER;  types_sizes(15) = int4_ts;        types_names(15) = "ieb"
      types(16) = H5T_NATIVE_INTEGER;  types_sizes(16) = int4_ts;        types_names(16) = "jsb"
      types(17) = H5T_NATIVE_INTEGER;  types_sizes(17) = int4_ts;        types_names(17) = "jeb"
      types(18) = H5T_NATIVE_INTEGER;  types_sizes(18) = int4_ts;        types_names(18) = "ksb"
      types(19) = H5T_NATIVE_INTEGER;  types_sizes(19) = int4_ts;        types_names(19) = "keb"

      types(20) = H5T_NATIVE_DOUBLE;   types_sizes(20) = r8_ts;          types_names(20) = "dx"
      types(21) = H5T_NATIVE_DOUBLE;   types_sizes(21) = r8_ts;          types_names(21) = "dy"
      types(22) = H5T_NATIVE_DOUBLE;   types_sizes(22) = r8_ts;          types_names(22) = "dz"
      types(23) = H5T_NATIVE_DOUBLE;   types_sizes(23) = r8_ts;          types_names(23) = "idx"
      types(24) = H5T_NATIVE_DOUBLE;   types_sizes(24) = r8_ts;          types_names(24) = "idy"
      types(25) = H5T_NATIVE_DOUBLE;   types_sizes(25) = r8_ts;          types_names(25) = "idz"
      types(26) = H5T_NATIVE_DOUBLE;   types_sizes(26) = r8_ts;          types_names(26) = "dxmn"
      types(27) = H5T_NATIVE_DOUBLE;   types_sizes(27) = r8_ts;          types_names(27) = "dvol"
      types(28) = H5T_NATIVE_DOUBLE;   types_sizes(28) = r8_ts;          types_names(28) = "xminb"
      types(29) = H5T_NATIVE_DOUBLE;   types_sizes(29) = r8_ts;          types_names(29) = "xmaxb"
      types(30) = H5T_NATIVE_DOUBLE;   types_sizes(30) = r8_ts;          types_names(30) = "yminb"
      types(31) = H5T_NATIVE_DOUBLE;   types_sizes(31) = r8_ts;          types_names(31) = "ymaxb"
      types(32) = H5T_NATIVE_DOUBLE;   types_sizes(32) = r8_ts;          types_names(32) = "zminb"
      types(33) = H5T_NATIVE_DOUBLE;   types_sizes(33) = r8_ts;          types_names(33) = "zmaxb"

      types(34) = ndims_r8_t;          types_sizes(34) = ndims_r8_ts;    types_names(34) = "dl"
      types(35) = ndims_r8_t;          types_sizes(35) = ndims_r8_ts;    types_names(35) = "idl"

      types(36) = nxarr_r8_t;          types_sizes(36) = nxarr_r8_ts;    types_names(36) = "x"
      types(37) = nxarr_r8_t;          types_sizes(37) = nxarr_r8_ts;    types_names(37) = "xl"
      types(38) = nxarr_r8_t;          types_sizes(38) = nxarr_r8_ts;    types_names(38) = "xr"
      types(39) = nxarr_r8_t;          types_sizes(39) = nxarr_r8_ts;    types_names(39) = "inv_x"

      types(40) = nyarr_r8_t;          types_sizes(40) = nyarr_r8_ts;    types_names(40) = "y"
      types(41) = nyarr_r8_t;          types_sizes(41) = nyarr_r8_ts;    types_names(41) = "yl"
      types(42) = nyarr_r8_t;          types_sizes(42) = nyarr_r8_ts;    types_names(42) = "yr"
      types(43) = nyarr_r8_t;          types_sizes(43) = nyarr_r8_ts;    types_names(43) = "inv_y"

      types(44) = nzarr_r8_t;          types_sizes(44) = nzarr_r8_ts;    types_names(44) = "z"
      types(45) = nzarr_r8_t;          types_sizes(45) = nzarr_r8_ts;    types_names(45) = "zl"
      types(46) = nzarr_r8_t;          types_sizes(46) = nzarr_r8_ts;    types_names(46) = "zr"
      types(47) = nzarr_r8_t;          types_sizes(47) = nzarr_r8_ts;    types_names(47) = "inv_z"

      types(48) = ndims_i4_t;          types_sizes(48) = ndims_i4_ts;    types_names(48) = "n_b"

      types(49) = ndims_i8_t;          types_sizes(49) = ndims_i8_ts;    types_names(49) = "off"

      types(50) = ndims_lohi_i4_t;     types_sizes(50) = ndims_lohi_i4_ts; types_names(50) = "ijkse"
      types(51) = ndims_lohi_i4_t;     types_sizes(51) = ndims_lohi_i4_ts; types_names(51) = "bnd"

      types(52) = arr3d_r8_t;          types_sizes(52) = arr3d_r8_ts;    types_names(52) = "gc_xdim"
      types(53) = arr3d_r8_t;          types_sizes(53) = arr3d_r8_ts;    types_names(53) = "gc_ydim"
      types(54) = arr3d_r8_t;          types_sizes(54) = arr3d_r8_ts;    types_names(54) = "gc_zdim"

      types_names(55) = "wa"
      if (associated(cg%wa%arr)) then
         types(55) = arr3d_r8_t;          types_sizes(55) = arr3d_r8_ts
      else
         types(55) = H5T_NATIVE_INTEGER;  types_sizes(55) = int4_ts
      endif
      types_names(56) = "gpot"
      if (associated(cg%gpot%arr)) then
         types(56) = arr3d_r8_t;          types_sizes(56) = arr3d_r8_ts
      else
         types(56) = H5T_NATIVE_INTEGER;  types_sizes(56) = int4_ts
      endif
      types_names(57) = "hgpot"
      if (associated(cg%hgpot%arr)) then
         types(57) = arr3d_r8_t;          types_sizes(57) = arr3d_r8_ts
      else
         types(57) = H5T_NATIVE_INTEGER;  types_sizes(57) = int4_ts
      endif
      types_names(58) = "gp"
      if (associated(cg%gp%arr)) then
         types(58) = arr3d_r8_t;          types_sizes(58) = arr3d_r8_ts
      else
         types(58) = H5T_NATIVE_INTEGER;  types_sizes(58) = int4_ts
      endif
      types_names(59) = "sgp"
      if (associated(cg%sgp%arr)) then
         types(59) = arr3d_r8_t;          types_sizes(59) = arr3d_r8_ts
      else
         types(59) = H5T_NATIVE_INTEGER;  types_sizes(59) = int4_ts
      endif
      types_names(60) = "sgpm"
      if (associated(cg%sgpm%arr)) then
         types(60) = arr3d_r8_t;          types_sizes(60) = arr3d_r8_ts
      else
         types(60) = H5T_NATIVE_INTEGER;  types_sizes(60) = int4_ts
      endif
      types_names(61) = "cs_iso2"
      if (associated(cg%cs_iso2%arr)) then
         types(61) = arr3d_r8_t;          types_sizes(61) = arr3d_r8_ts
      else
         types(61) = H5T_NATIVE_INTEGER;  types_sizes(61) = int4_ts
      endif

      types(62) = ndims_arr4d_r8_t;    types_sizes(62) = ndims_arr4d_r8_ts;    types_names(62) = "b"
      types_names(63) = "b0"
      if (associated(cg%b0%arr)) then
         types(63) = ndims_arr4d_r8_t;    types_sizes(63) = ndims_arr4d_r8_ts
      else
         types(63) = H5T_NATIVE_INTEGER;  types_sizes(63) = int4_ts
      endif

      types(64) = u_arr4d_r8_t;        types_sizes(64) = u_arr4d_r8_ts;        types_names(64) = "u"
      types(65) = u_arr4d_r8_t;        types_sizes(65) = u_arr4d_r8_ts;        types_names(65) = "uh"
      types_names(66) = "u0"
      if (associated(cg%u0%arr)) then
         types(66) = u_arr4d_r8_t;        types_sizes(66) = u_arr4d_r8_ts
      else
         types(66) = H5T_NATIVE_INTEGER;  types_sizes(66) = int4_ts
      endif

      offset = 0
      do i = 1, total_no
         call h5tinsert_f(dtype_id, types_names(i),  offset, types(i), error)
         offset = offset + types_sizes(i)
      enddo

      write(dset_name,'("cg",i4.4)') 1
      call h5dcreate_f(file_id, dset_name, dtype_id, dspace_id, dset_id, error)

      do i = 1, total_no
         call h5tcreate_f(H5T_COMPOUND_F, types_sizes(i), dmem_id(i), error)
         offset = 0
         call h5tinsert_f(dmem_id(i), types_names(i), offset, types(i), error)
      enddo

      dims = 1
      call h5dwrite_f(dset_id, dmem_id(1),  int(cg%n_(xdim), kind=4), dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(2),  int(cg%n_(ydim), kind=4), dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(3),  int(cg%n_(zdim), kind=4), dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(4),  int(cg%nxb, kind=4),    dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(5),  int(cg%nyb, kind=4),    dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(6),  int(cg%nzb, kind=4),    dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(7),  int(cg%is, kind=4),     dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(8),  int(cg%ie, kind=4),     dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(9),  int(cg%js, kind=4),     dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(10), int(cg%je, kind=4),     dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(11), int(cg%ks, kind=4),     dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(12), int(cg%ke, kind=4),     dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(13), int(cg%maxxyz, kind=4), dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(14), int(cg%isb, kind=4),    dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(15), int(cg%ieb, kind=4),    dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(16), int(cg%jsb, kind=4),    dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(17), int(cg%jeb, kind=4),    dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(18), int(cg%ksb, kind=4),    dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(19), int(cg%keb, kind=4),    dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(20), cg%dx,     dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(21), cg%dy,     dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(22), cg%dz,     dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(23), cg%idx,    dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(24), cg%idy,    dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(25), cg%idz,    dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(26), cg%dxmn,   dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(27), cg%dvol,   dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(28), cg%fbnd(xdim, LO), dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(29), cg%fbnd(xdim, HI), dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(30), cg%fbnd(ydim, LO), dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(31), cg%fbnd(ydim, HI), dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(32), cg%fbnd(zdim, LO), dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(33), cg%fbnd(zdim, HI), dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(34), cg%dl,     dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(35), cg%idl,    dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(36), cg%x,      dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(37), cg%xl,     dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(38), cg%xr,     dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(39), cg%inv_x,  dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(40), cg%y,      dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(41), cg%yl,     dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(42), cg%yr,     dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(43), cg%inv_y,  dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(44), cg%z,      dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(45), cg%zl,     dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(46), cg%zr,     dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(47), cg%inv_z,  dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(48), int(cg%n_b, kind=4),    dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(49), int(cg%n_b, kind=4),    dims, error, xfer_prp=plist_id) !!!
      call h5dwrite_f(dset_id, dmem_id(50), int(cg%ijkse, kind=4),  dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(51), int(cg%bnd, kind=4),    dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(52), int(cg%gc_xdim, kind=4),dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(53), int(cg%gc_ydim, kind=4),dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(54), int(cg%gc_zdim, kind=4),dims, error, xfer_prp=plist_id)
      if (associated(cg%wa%arr)) then
         call h5dwrite_f(dset_id, dmem_id(55), cg%wa%arr     ,dims, error, xfer_prp=plist_id)
      else
         call h5dwrite_f(dset_id, dmem_id(55), -999_INT4      ,dims, error, xfer_prp=plist_id)
      endif
      if (associated(cg%gpot%arr)) then
         call h5dwrite_f(dset_id, dmem_id(56), cg%gpot%arr   ,dims, error, xfer_prp=plist_id)
      else
         call h5dwrite_f(dset_id, dmem_id(56), -999_INT4      ,dims, error, xfer_prp=plist_id)
      endif
      if (associated(cg%hgpot%arr)) then
         call h5dwrite_f(dset_id, dmem_id(57), cg%hgpot%arr  ,dims, error, xfer_prp=plist_id)
      else
         call h5dwrite_f(dset_id, dmem_id(57), -999_INT4      ,dims, error, xfer_prp=plist_id)
      endif
      if (associated(cg%gp%arr)) then
         call h5dwrite_f(dset_id, dmem_id(58), cg%gp%arr     ,dims, error, xfer_prp=plist_id)
      else
         call h5dwrite_f(dset_id, dmem_id(58), -999_INT4      ,dims, error, xfer_prp=plist_id)
      endif
      if (associated(cg%sgp%arr)) then
         call h5dwrite_f(dset_id, dmem_id(59), cg%sgp%arr    ,dims, error, xfer_prp=plist_id)
      else
         call h5dwrite_f(dset_id, dmem_id(59), -999_INT4      ,dims, error, xfer_prp=plist_id)
      endif
      if (associated(cg%sgpm%arr)) then
         call h5dwrite_f(dset_id, dmem_id(60), cg%sgpm%arr   ,dims, error, xfer_prp=plist_id)
      else
         call h5dwrite_f(dset_id, dmem_id(60), -999_INT4      ,dims, error, xfer_prp=plist_id)
      endif
      if (associated(cg%cs_iso2%arr)) then
         call h5dwrite_f(dset_id, dmem_id(61), cg%cs_iso2%arr,dims, error, xfer_prp=plist_id)
      else
         call h5dwrite_f(dset_id, dmem_id(61), -999_INT4      ,dims, error, xfer_prp=plist_id)
      endif
      call h5dwrite_f(dset_id, dmem_id(62), cg%b%arr,dims, error, xfer_prp=plist_id)
      if (associated(cg%b0%arr)) then
         call h5dwrite_f(dset_id, dmem_id(63), cg%b0%arr     ,dims, error, xfer_prp=plist_id)
      else
         call h5dwrite_f(dset_id, dmem_id(63), -999_INT4      ,dims, error, xfer_prp=plist_id)
      endif
      call h5dwrite_f(dset_id, dmem_id(64), cg%u%arr,  dims, error, xfer_prp=plist_id)
      call h5dwrite_f(dset_id, dmem_id(65), cg%uh%arr, dims, error, xfer_prp=plist_id)
      if (associated(cg%u0%arr)) then
         call h5dwrite_f(dset_id, dmem_id(66), cg%u0%arr     ,dims, error, xfer_prp=plist_id)
      else
         call h5dwrite_f(dset_id, dmem_id(66), -999_INT4      ,dims, error, xfer_prp=plist_id)
      endif

      call h5dclose_f(dset_id, error)
      do i = 1, total_no
         call h5tclose_f(dmem_id(i),error)
      enddo
      CALL h5sclose_f(dspace_id, error)
      CALL h5tclose_f(dtype_id, error)
      call h5tclose_f(ndims_r8_t, error)
      call h5tclose_f(ndims_i4_t, error)
      call h5tclose_f(ndims_i8_t, error)
      call h5tclose_f(ndims_lohi_i4_t, error)
      call h5tclose_f(nxarr_r8_t, error)
      call h5tclose_f(nyarr_r8_t, error)
      call h5tclose_f(nzarr_r8_t, error)
      call h5tclose_f(arr3d_r8_t, error)
      call h5tclose_f(ndims_arr4d_r8_t, error)
      call h5tclose_f(u_arr4d_r8_t, error)

      deallocate(types,types_sizes,types_names,dmem_id)

   end subroutine write_grid_containter

end module dataio_hdf5
