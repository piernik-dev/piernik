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
!! \brief Module that contains HDF5 I/O routines for writing plotfiles (2D slices)
!<
module slice_hdf5

! pulled by ANY

   use constants, only: xdim, zdim

   implicit none

   private
   public :: init_plot, write_plot

   integer, parameter :: planelen = 2           !< length of plane names e.g. "xy", "yz", "rp" etc.
   character(len=planelen), dimension(xdim:zdim), parameter :: pl_id = [ "yz", "xz", "xy" ]
   integer, dimension(xdim:zdim) :: pl_i !< no. of cell ( 1 <= pl_i(:) < dom%n_d(:) ) for YZ, XZ and XY slices in plt files
   real :: dt_plt !< frequency of plt output

contains

   subroutine init_plot(ti, tdt_plt)

      use constants, only: ndims

      implicit none

      integer, dimension(ndims), intent(in) :: ti  !< local copy of [ dataio::ix, dataio::iy, dataio::iz ]
      real, intent(in)    :: tdt_plt !< local copy of dataio::dt_plt

      pl_i(:) = ti(:)
      dt_plt = tdt_plt

   end subroutine init_plot

!>
!! \brief Routine calculating quantities for plot files
!<
   subroutine common_plt_hdf5(var, ij, xn, tab, ierrh, cg)

      use constants,   only: varlen, xdim, ydim, zdim, half
      use common_hdf5, only: common_shortcuts
      use grid_cont,   only: grid_container
      use fluidindex,  only: flind
      use fluidtypes,  only: component_fluid
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
            tab(:,:) = reshape(half * ( cg%u%arr(fl_dni%imx, is:ie, js:je, ks:ke)**2 + &
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
                    half *( cg%u%arr(flind%neu%imx, is:ie, js:je, xn)**2 + &
                    &      cg%u%arr(flind%neu%imy, is:ie, js:je, xn)**2 + &
                    &      cg%u%arr(flind%neu%imz, is:ie, js:je, xn)**2 ) / cg%u%arr(flind%neu%idn, is:ie, js:je, xn), kind=4)*(flind%neu%gam_1)
            endif
#endif /* !ISO */
         case ("magx", "magy", "magz")
            tab(:,:) = reshape(cg%b%arr(xdim + i_xyz, is:ie, js:je, ks:ke), shape(tab))
#ifdef GRAV
         case ("gpot")
            tab(:,:) = reshape(cg%gpot(is:ie, js:je, ks:ke), shape(tab))
#endif /* GRAV */
#ifdef COSM_RAYS
         case ("cr*")
            i = iarr_all_crs(ichar(var(3:3))-ichar('0'))
            tab(:,:) = reshape(cg%u%arr(i, is:ie, js:je, ks:ke), shape(tab))
#endif /* COSM_RAYS */
         case default
            ierrh = -1
      end select

   end subroutine common_plt_hdf5

   subroutine write_plot

      use constants,   only: cwdlen, xdim, zdim
      use common_hdf5, only: nhdf_vars, hdf_vars
      use dataio_pub,  only: log_file
      use global,      only: t
      use hdf5,        only: HID_T, H5open_f, H5Fcreate_f, H5Gcreate_f, H5F_ACC_TRUNC_F, H5Gclose_f, H5close_f, h5fclose_f
      use mpisetup,    only: comm, ierr, master

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
      use domain,      only: dom, has_dir, is_multicg
      use global,      only: t
      use grid,        only: all_cg
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

      cg => all_cg%first%cg
      if (is_multicg) call die("[slice_hdf5:write_plot_hdf5] multiple grid pieces per procesor not implemented yet") !nontrivial message tagging

      xn = 1
      if (has_dir(plane)) xn = pl_i(plane) + cg%nb - cg%off(plane)

      if ((xn > cg%nb .and. xn <= cg%n_b(plane)+cg%nb) .or. (xn == 1 .and. .not. has_dir(plane))) then
         allocate(send(cg%n_b(d1(plane)), cg%n_b(d2(plane))))
         call common_plt_hdf5(var, plane, xn, send, ierrh, cg)
         if (associated(user_plt_hdf5) .and. ierrh /= 0) then
            ierrh = 0
            call user_plt_hdf5(var, plane, xn, send, ierrh, cg)
         endif
         if (ierrh /= 0) then
            write(msg,'(3a)')"[slice_hdf5:write_plot_hdf5]", var, " is not defined in common_plt_hdf5, neither in user_plt_hdf5 !!!"
            call warn(msg)
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
            call warn("[slice_hdf5:write_plot_hdf5] vizit used without PGPLOT")
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

end module slice_hdf5
