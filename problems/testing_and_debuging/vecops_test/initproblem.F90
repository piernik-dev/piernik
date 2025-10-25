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

module initproblem

!>
!! \brief Initial condition for testing the type bound procedure on cg for divergence/curl/cross product/dot product
!! The first vector is set by the fluid velocity :  vx = cosy sinx , vy = - cosx siny , vz = 0.0 and density = 1
!! The second vector is set by the magnetic field : Bx = - siny cosx , By = - sinx cosy , Bz = 0.0
!! The dump_vec_gnuplot subroutine have been written with the help of LLM (ChatGPT)
!<

   use constants, only: dsetnamelen

   implicit none



   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   ! namelist parameters

   integer :: order

   namelist /PROBLEM_CONTROL/  order

   ! other private data
   character(len=dsetnamelen), parameter :: divvnum = "divvnum"
   character(len=dsetnamelen), parameter :: divvana = "divvana"
   character(len=dsetnamelen), parameter :: divbnum = "divbnum"
   character(len=dsetnamelen), parameter :: divbana = "divbana"
   character(len=dsetnamelen), parameter :: curlvnum = "curlvnum"
   character(len=dsetnamelen), parameter :: curlbnum = "curlbnum"
   character(len=dsetnamelen), parameter :: vcrossb  = "vcrossb"
   character(len=dsetnamelen), parameter :: vdotb    = "vdotb"


contains

!> \brief Set up custom pointers to tweak the code execution according to our needs

   subroutine problem_pointers

      use dataio_user,     only: user_vars_hdf5
      use user_hooks,      only: finalize_problem

      implicit none

      user_vars_hdf5             => user_out
      finalize_problem           => dump_vec_gnuplot

   end subroutine problem_pointers

!> \brief Read the runtime parameters specified in the namelist

   subroutine read_problem_par

      use bcast,          only: piernik_MPI_Bcast
      use cg_list_global, only: all_cg
      use constants,      only: AT_IGNORE
      use dataio_pub,     only: nh
      use mpisetup,       only: ibuff, master, slave

      implicit none

      ! the default values
      order  = 2


      if (master) then

         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=PROBLEM_CONTROL)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=""
         read(unit=nh%lun, nml=PROBLEM_CONTROL, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "PROBLEM_CONTROL")
         read(nh%cmdl_nml,nml=PROBLEM_CONTROL, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "PROBLEM_CONTROL", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=PROBLEM_CONTROL)
         close(nh%lun)
         call nh%compare_namelist()

         ibuff(1) = order

      endif

      call piernik_MPI_Bcast(ibuff)

      if (slave) then

         order = ibuff(1)

      endif


      call all_cg%reg_var(divvnum,  restart_mode = AT_IGNORE )
      call all_cg%reg_var(divvana,  restart_mode = AT_IGNORE )
      call all_cg%reg_var(divbnum,  restart_mode = AT_IGNORE )
      call all_cg%reg_var(divbana,  restart_mode = AT_IGNORE )
      call all_cg%reg_var(curlvnum, dim4 = 3, restart_mode = AT_IGNORE)
      call all_cg%reg_var(curlbnum, dim4 = 3, restart_mode = AT_IGNORE)
      call all_cg%reg_var(vcrossb,  dim4 = 3, restart_mode = AT_IGNORE)
      call all_cg%reg_var(vdotb,              restart_mode = AT_IGNORE)


   end subroutine read_problem_par


   subroutine problem_initial_conditions

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: xdim, ydim, zdim, LO, HI
      use fluidindex,       only: flind
      use fluidtypes,       only: component_fluid
      use func,             only: ekin, emag
      use grid_cont,        only: grid_container
      use named_array_list, only: wna, qna
      implicit none

      class(component_fluid), pointer :: fl
      integer                         :: i, j, k
      real                            :: xi, yj, zk, sx, sy, cx, cy
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg
      integer                         :: p
      integer :: vmap(3)

      do p = 1, flind%fluids

         fl => flind%all_fluids(p)%fl

         cgl => leaves%first
         do while (associated(cgl))
            cg => cgl%cg

            do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
               yj = cg%y(j)
               do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
                  xi = cg%x(i)
                  do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
                     zk = cg%z(k)
                     sx = sin(xi);  cx = cos(xi)
                     sy = sin(yj);  cy = cos(yj)
                     cg%u(fl%idn, i, j, k) = 1.0
                     cg%u(fl%imx,i,j,k) =  cy*sx
                     cg%u(fl%imy,i,j,k) = -cx*sy
                     cg%u(fl%imz, i, j, k) = 0.0

                     cg%q(qna%ind(divvana))%arr(i,j,k) = 0.0


                     if (fl%has_energy) then
                        cg%u(fl%ien,i,j,k) = 1.0
                        cg%u(fl%ien,i,j,k) = cg%u(fl%ien,i,j,k) + ekin(cg%u(fl%imx,i,j,k), cg%u(fl%imy,i,j,k), cg%u(fl%imz,i,j,k), cg%u(fl%idn,i,j,k))

                        if (fl%is_magnetized) then
                           cg%q(qna%ind(divbana))%arr(i, j, k) = 0.0
                           cg%b(xdim,i,j,k) = -sy*cx
                           cg%b(ydim,i,j,k) =  sx*cy
                           cg%b(zdim,i,j,k) =  0.0
                           cg%u(fl%ien,i,j,k) = cg%u(fl%ien,i,j,k) + emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))
                        endif
                     endif

                  enddo
               enddo
            enddo

            cgl => cgl%nxt
         enddo
      enddo
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         cg%q(qna%ind(divvnum))%arr = cg%get_divergence(ord=order, iw=wna%fi, vec=[2,3,4])
         cg%q(qna%ind(divbnum))%arr = cg%get_divergence(ord=order, iw=wna%bi, vec=[1,2,3])
         cg%w(wna%ind(curlbnum))%arr = cg%get_curl(ord=order, iw=wna%bi, vec=[1,2,3])
         cg%w(wna%ind(curlvnum))%arr = cg%get_curl(ord=order, iw=wna%fi, vec=[2,3,4])
         cg%w(wna%ind(vcrossb))%arr = cg%cross(wna%fi,wna%bi,[2,3,4],[1,2,3])
         cg%q(qna%ind(vdotb))%arr =  cg%dot(wna%fi,wna%bi,[2,3,4],[1,2,3])

         cgl => cgl%nxt
      enddo
   end subroutine problem_initial_conditions

#ifdef HDF5
   subroutine user_out(var, tab, ierrh, cg)

      use dataio_pub,       only: die
      use grid_cont,        only: grid_container
      use named_array_list, only: qna, wna

      implicit none

      character(len=*),               intent(in)    :: var
      real, dimension(:,:,:),         intent(inout) :: tab
      integer,                        intent(inout) :: ierrh
      type(grid_container), pointer,  intent(in)    :: cg


      if (.not. qna%exists(divvana)) call die("[initproblem:user_out] cannot find divvana")
      if (.not. qna%exists(divvana)) call die("[initproblem:user_out] cannot find divvnum")

      ierrh = 0
      select case (trim(var))
         case ("divv_ana")
            tab(:,:,:) = real(cg%q(qna%ind(divvana))%span(cg%ijkse), kind(tab))
         case ("divv_num")
            tab(:,:,:) = real(cg%q(qna%ind(divvnum))%span(cg%ijkse), kind(tab))
         case ("divb_ana")
            tab(:,:,:) = real(cg%q(qna%ind(divbana))%span(cg%ijkse), kind(tab))
         case ("divb_num")
            tab(:,:,:) = real(cg%q(qna%ind(divbnum))%span(cg%ijkse), kind(tab))

         case ("curlv_num_x")
            tab(:,:,:) = real( cg%w(wna%ind(curlvnum))%arr(1, RNG), kind(tab) )
         case ("curlv_num_y")
            tab(:,:,:) = real( cg%w(wna%ind(curlvnum))%arr(2, RNG), kind(tab) )
         case ("curlv_num_z")
            tab(:,:,:) = real( cg%w(wna%ind(curlvnum))%arr(3, RNG), kind(tab) )

         case ("curlb_num_x")
            tab(:,:,:) = real( cg%w(wna%ind(curlbnum))%arr(1, RNG), kind(tab) )
         case ("curlb_num_y")
            tab(:,:,:) = real( cg%w(wna%ind(curlbnum))%arr(2, RNG), kind(tab) )
         case ("curlb_num_z")
            tab(:,:,:) = real( cg%w(wna%ind(curlbnum))%arr(3, RNG), kind(tab) )

         case ("vcrossb_x")
            tab(:,:,:) = real( cg%w(wna%ind(vcrossb))%arr(1,RNG), kind(tab) )
         case ("vcrossb_y")
            tab(:,:,:) = real( cg%w(wna%ind(vcrossb))%arr(2,RNG), kind(tab) )
         case ("vcrossb_z")
            tab(:,:,:) = real( cg%w(wna%ind(vcrossb))%arr(3,RNG), kind(tab) )
         case ("vdotb")
            tab(:,:,:) = real( cg%q(qna%ind(vdotb))%arr(RNG), kind(tab) )
         case default
            ierrh = -1
      end select

   end subroutine user_out

#endif /* HDF5 */
   subroutine verify_test
      use allreduce,        only: piernik_MPI_Allreduce
      use cg_list,          only: cg_list_element
      use cg_leaves,        only: leaves
      use constants,        only: V_ESSENTIAL, V_INFO, pSUM, pMAX
      use dataio_pub,       only: msg, printinfo
      use grid_cont,        only: grid_container
      use mpisetup,         only: master
      use named_array_list, only: qna, wna
      implicit none
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg

      real :: ssum, smax

      integer :: i, j, k
      real :: xi, yj, sx, cx, sy, cy, dv
      real :: val, ana, ex, ey, ez, err2, lsum, lmax

      if (master) call printinfo('--- Diagnostics (global volume-weighted L2 errors) ---', V_INFO)

      ssum = 0.0; smax = 0.0
      cgl => leaves%first
      do while (associated(cgl))
         cg  => cgl%cg
         dv  = cg%dvol
         lsum = 0.0; lmax = 0.0
         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               do i = cg%is, cg%ie
                  val  = cg%q(qna%ind(divvnum))%arr(i,j,k)
                  err2 = val*val
                  lsum = lsum + err2
                  lmax = max(lmax, abs(val))
               enddo
            enddo
         enddo
         ssum = ssum + lsum*dv
         smax = max(smax, lmax)
         cgl => cgl%nxt
      enddo
      call piernik_MPI_Allreduce(ssum, pSUM)
      call piernik_MPI_Allreduce(smax, pMAX)
      if (master) then
         write(msg,'("||div(v)||_err:  L2_vw=",1pe12.5,"  Linf=",1pe12.5)') sqrt(ssum), smax
         call printinfo(msg, V_ESSENTIAL)
      endif

      ssum = 0.0; smax = 0.0
      cgl => leaves%first
      do while (associated(cgl))
         cg  => cgl%cg
         dv  = cg%dvol
         lsum = 0.0; lmax = 0.0
         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               do i = cg%is, cg%ie
                  val  = cg%q(qna%ind(divbnum))%arr(i,j,k)
                  err2 = val*val
                  lsum = lsum + err2
                  lmax = max(lmax, abs(val))
               enddo
            enddo
         enddo
         ssum = ssum + lsum*dv
         smax = max(smax, lmax)
         cgl => cgl%nxt
      enddo
      call piernik_MPI_Allreduce(ssum, pSUM)
      call piernik_MPI_Allreduce(smax, pMAX)
      if (master) then
         write(msg,'("||div(B)||_err:  L2_vw=",1pe12.5,"  Linf=",1pe12.5)') sqrt(ssum), smax
         call printinfo(msg, V_ESSENTIAL)
      endif

      ssum = 0.0; smax = 0.0
      cgl => leaves%first
      do while (associated(cgl))
         cg  => cgl%cg
         dv  = cg%dvol
         lsum = 0.0; lmax = 0.0
         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               yj = cg%y(j);  sy = sin(yj);  cy = cos(yj)
               do i = cg%is, cg%ie
                  xi = cg%x(i); sx = sin(xi); cx = cos(xi)
                  ex = cg%w(wna%ind(curlvnum))%arr(1,i,j,k) - 0.0
                  ey = cg%w(wna%ind(curlvnum))%arr(2,i,j,k) - 0.0
                  ez = cg%w(wna%ind(curlvnum))%arr(3,i,j,k) - (2.0*sx*sy)
                  err2 = ex*ex + ey*ey + ez*ez
                  lsum = lsum + err2
                  lmax = max(lmax, sqrt(err2))
               enddo
            enddo
         enddo
         ssum = ssum + lsum*dv
         smax = max(smax, lmax)
         cgl => cgl%nxt
      enddo
      call piernik_MPI_Allreduce(ssum, pSUM)
      call piernik_MPI_Allreduce(smax, pMAX)
      if (master) then
         write(msg,'("||curl(v)||_err: L2_vw=",1pe12.5,"  Linf=",1pe12.5)') sqrt(ssum), smax
         call printinfo(msg, V_ESSENTIAL)
      endif

      ssum = 0.0; smax = 0.0
      cgl => leaves%first
      do while (associated(cgl))
         cg  => cgl%cg
         dv  = cg%dvol
         lsum = 0.0; lmax = 0.0
         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               yj = cg%y(j);  sy = sin(yj);  cy = cos(yj)
               do i = cg%is, cg%ie
                  xi = cg%x(i); sx = sin(xi); cx = cos(xi)
                  ex = cg%w(wna%ind(curlbnum))%arr(1,i,j,k) - 0.0
                  ey = cg%w(wna%ind(curlbnum))%arr(2,i,j,k) - 0.0
                  ez = cg%w(wna%ind(curlbnum))%arr(3,i,j,k) - (2.0*cx*cy)
                  err2 = ex*ex + ey*ey + ez*ez
                  lsum = lsum + err2
                  lmax = max(lmax, sqrt(err2))
               enddo
            enddo
         enddo
         ssum = ssum + lsum*dv
         smax = max(smax, lmax)
         cgl => cgl%nxt
      enddo
      call piernik_MPI_Allreduce(ssum, pSUM)
      call piernik_MPI_Allreduce(smax, pMAX)
      if (master) then
         write(msg,'("||curl(B)||_err: L2_vw=",1pe12.5,"  Linf=",1pe12.5)') sqrt(ssum), smax
         call printinfo(msg, V_ESSENTIAL)
      endif

      ssum = 0.0; smax = 0.0
      cgl => leaves%first
      do while (associated(cgl))
         cg  => cgl%cg
         dv  = cg%dvol
         lsum = 0.0; lmax = 0.0
         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               yj = cg%y(j);  sy = sin(yj);  cy = cos(yj)
               do i = cg%is, cg%ie
                  xi = cg%x(i); sx = sin(xi); cx = cos(xi)
                  ex = cg%w(wna%ind(vcrossb))%arr(1,i,j,k) - 0.0
                  ey = cg%w(wna%ind(vcrossb))%arr(2,i,j,k) - 0.0
                  ez = cg%w(wna%ind(vcrossb))%arr(3,i,j,k) - (sx*sx*cy*cy - cx*cx*sy*sy)
                  err2 = ex*ex + ey*ey + ez*ez
                  lsum = lsum + err2
                  lmax = max(lmax, sqrt(err2))
               enddo
            enddo
         enddo
         ssum = ssum + lsum*dv
         smax = max(smax, lmax)
         cgl => cgl%nxt
      enddo
      call piernik_MPI_Allreduce(ssum, pSUM)
      call piernik_MPI_Allreduce(smax, pMAX)
      if (master) then
         write(msg,'("||v×B||_err  : L2_vw=",1pe12.5,"  Linf=",1pe12.5)') sqrt(ssum), smax
         call printinfo(msg, V_ESSENTIAL)
      endif

      ssum = 0.0; smax = 0.0
      cgl => leaves%first
      do while (associated(cgl))
         cg  => cgl%cg
         dv  = cg%dvol
         lsum = 0.0; lmax = 0.0
         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               yj = cg%y(j);  sy = sin(yj);  cy = cos(yj)
               do i = cg%is, cg%ie
                  xi = cg%x(i); sx = sin(xi); cx = cos(xi)
                  ana = -2.0*sx*cx*sy*cy
                  val = cg%q(qna%ind(vdotb))%arr(i,j,k)
                  err2 = (val - ana)**2
                  lsum = lsum + err2
                  lmax = max(lmax, abs(val - ana))
               enddo
            enddo
         enddo
         ssum = ssum + lsum*dv
         smax = max(smax, lmax)
         cgl => cgl%nxt
      enddo
      call piernik_MPI_Allreduce(ssum, pSUM)
      call piernik_MPI_Allreduce(smax, pMAX)
      if (master) then
         write(msg,'("||v·B||_err  : L2_vw=",1pe12.5,"  Linf=",1pe12.5)') sqrt(ssum), smax
         call printinfo(msg, V_ESSENTIAL)
      endif
   end subroutine verify_test

   subroutine dump_vec_gnuplot
      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use grid_cont,        only: grid_container
      use named_array_list, only: qna, wna
      use constants,        only: V_INFO
      use dataio_pub,       only: printinfo
      use mpisetup,         only: master
      implicit none

      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
      integer :: i, j, j0, k0, lun
      real :: x, y, z0, sx, cx, sy, cy
      logical :: wrote_files
      call verify_test
      if (.not. master) return

      call printinfo('Writing .dat files for the gnuplot. 1D is along x for y = midpoint of the domain', V_INFO)
      wrote_files = .false.

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         ! Mid-line in y, mid-plane in z
         j0 = (cg%js + cg%je)/2
         k0 = (cg%ks + cg%ke)/2
         y  = cg%y(j0);  sy = sin(y);  cy = cos(y)
         z0 = cg%z(k0)

         ! ----------------------- 1D: div(v) -----------------------
         open(newunit=lun, file='divv_line.dat', status='unknown')
         do i = cg%is, cg%ie
            x  = cg%x(i)
            write(lun,'(4(1pe23.15))') x, 0.0, cg%q(qna%ind(divvnum))%arr(i,j0,k0), &
                                       (cg%q(qna%ind(divvnum))%arr(i,j0,k0) - 0.0)
         enddo
         close(lun)

         ! ----------------------- 1D: div(B) -----------------------
         open(newunit=lun, file='divb_line.dat', status='unknown')
         do i = cg%is, cg%ie
            x  = cg%x(i)
            write(lun,'(4(1pe23.15))') x, 0.0, cg%q(qna%ind(divbnum))%arr(i,j0,k0), &
                                       (cg%q(qna%ind(divbnum))%arr(i,j0,k0) - 0.0)
         enddo
         close(lun)

         ! ---------------------- 1D: curl(v)_z (2 sin x sin y) -----
         open(newunit=lun, file='curlv_z_line.dat', status='unknown')
         do i = cg%is, cg%ie
            x  = cg%x(i);  sx = sin(x)
            write(lun,'(4(1pe23.15))') x, 2.0*sx*sy, cg%w(wna%ind(curlvnum))%arr(3,i,j0,k0), &
                                       (cg%w(wna%ind(curlvnum))%arr(3,i,j0,k0) - 2.0*sx*sy)
         enddo
         close(lun)

         ! ---------------------- 1D: curl(B)_z (2 cos x cos y) -----
         open(newunit=lun, file='curlb_z_line.dat', status='unknown')
         do i = cg%is, cg%ie
            x  = cg%x(i);  cx = cos(x)
            write(lun,'(4(1pe23.15))') x, 2.0*cx*cy, cg%w(wna%ind(curlbnum))%arr(3,i,j0,k0), &
                                       (cg%w(wna%ind(curlbnum))%arr(3,i,j0,k0) - 2.0*cx*cy)
         enddo
         close(lun)

         ! ----------------------- 1D: (v x B)_z --------------------
         open(newunit=lun, file='vcrossb_z_line.dat', status='unknown')
         do i = cg%is, cg%ie
            x  = cg%x(i);  sx = sin(x); cx = cos(x)
            write(lun,'(4(1pe23.15))') x, (sx*sx*cy*cy - cx*cx*sy*sy), cg%w(wna%ind(vcrossb))%arr(3,i,j0,k0), &
                                       (cg%w(wna%ind(vcrossb))%arr(3,i,j0,k0) - (sx*sx*cy*cy - cx*cx*sy*sy))
         enddo
         close(lun)

         ! ------------------------- 1D: v·B ------------------------
         open(newunit=lun, file='vdotb_line.dat', status='unknown')
         do i = cg%is, cg%ie
            x  = cg%x(i);  sx = sin(x); cx = cos(x)
            write(lun,'(4(1pe23.15))') x, (-2.0*sx*cx*sy*cy), cg%q(qna%ind(vdotb))%arr(i,j0,k0), &
                                       (cg%q(qna%ind(vdotb))%arr(i,j0,k0) - (-2.0*sx*cx*sy*cy))
         enddo
         close(lun)

         ! ----------------------- 2D slice dumps (x,y at fixed z) --
         ! columns: x y ana num
         open(newunit=lun, file='divv_2d.dat', status='unknown')
         do j = cg%js, cg%je
            y = cg%y(j)
            do i = cg%is, cg%ie
               x  = cg%x(i)
               write(lun,'(4(1pe23.15))') x, y, 0.0, cg%q(qna%ind(divvnum))%arr(i,j,k0)
            enddo
            write(lun,'(a)') ''
         enddo
         close(lun)

         open(newunit=lun, file='divb_2d.dat', status='unknown')
         do j = cg%js, cg%je
            y = cg%y(j)
            do i = cg%is, cg%ie
               x  = cg%x(i)
               write(lun,'(4(1pe23.15))') x, y, 0.0, cg%q(qna%ind(divbnum))%arr(i,j,k0)
            enddo
            write(lun,'(a)') ''
         enddo
         close(lun)

         open(newunit=lun, file='curlv_z_2d.dat', status='unknown')
         do j = cg%js, cg%je
            y = cg%y(j); sy = sin(y)
            do i = cg%is, cg%ie
               x  = cg%x(i); sx = sin(x)
               write(lun,'(4(1pe23.15))') x, y, 2.0*sx*sy, cg%w(wna%ind(curlvnum))%arr(3,i,j,k0)
            enddo
            write(lun,'(a)') ''
         enddo
         close(lun)

         open(newunit=lun, file='curlb_z_2d.dat', status='unknown')
         do j = cg%js, cg%je
            y = cg%y(j); cy = cos(y)
            do i = cg%is, cg%ie
               x  = cg%x(i); cx = cos(x)
               write(lun,'(4(1pe23.15))') x, y, 2.0*cx*cy, cg%w(wna%ind(curlbnum))%arr(3,i,j,k0)
            enddo
            write(lun,'(a)') ''
         enddo
         close(lun)

         open(newunit=lun, file='vcrossb_z_2d.dat', status='unknown')
         do j = cg%js, cg%je
            y = cg%y(j); sy = sin(y); cy = cos(y)
            do i = cg%is, cg%ie
               x  = cg%x(i); sx = sin(x); cx = cos(x)
               write(lun,'(4(1pe23.15))') x, y, (sx*sx*cy*cy - cx*cx*sy*sy), cg%w(wna%ind(vcrossb))%arr(3,i,j,k0)
            enddo
            write(lun,'(a)') ''
         enddo
         close(lun)

         open(newunit=lun, file='vdotb_2d.dat', status='unknown')
         do j = cg%js, cg%je
            y = cg%y(j); sy = sin(y); cy = cos(y)
            do i = cg%is, cg%ie
               x  = cg%x(i); sx = sin(x); cx = cos(x)
               write(lun,'(4(1pe23.15))') x, y, (-2.0*sx*cx*sy*cy), cg%q(qna%ind(vdotb))%arr(i,j,k0)
            enddo
            write(lun,'(a)') ''
         enddo
         close(lun)

         wrote_files = .true.
         exit
      enddo

      if (master .and. wrote_files) then
         open(newunit=lun, file='test_plot_1d.gp', status='unknown')
         write(lun,'(a)') 'set encoding utf8'
         write(lun,'(a)') 'set term pngcairo size 1100,700 font "Helvetica,10" enhanced'
         write(lun,'(a)') 'set grid xtics ytics mxtics mytics back lc rgb "#cccccc"'
         write(lun,'(a)') 'set mxtics 2; set mytics 2'
         write(lun,'(a)') 'set key at graph 0.5,0.98 center top horizontal opaque samplen 2 spacing 1.2'
         write(lun,'(a)') 'set tics out; set format y "%1.2g"; set format y2 "%1.2g"; set y2tics'
         write(lun,'(a)') 'set style line 1 lc rgb "#7b3294" lw 2.2'                  ! Analytical
         write(lun,'(a)') 'set style line 2 lc rgb "#008837" lw 2.2'                  ! Numerical
         write(lun,'(a)') 'set style line 3 lc rgb "#E69F00" lw 2.6 dt 1 pt 5 ps 0.8' ! Residual
         write(lun,'(a)') 'set style line 4 lc rgb "#E69F00" lw 2.6 dt 3 pt 5 ps 0.8' ! 100×Residual
         write(lun,'(a)') 'set style line 5 lc rgb "#0571B0" lw 3.0 dt 2'             ! |err| (y2)
         write(lun,'(a)') 'set style line 9 lc rgb "#BBBBBB" lw 1.0 dt 3'             ! zero line
         write(lun,'(a)') 'sf_fixed = 100.0'                                          ! scale for 100×Residual
         write(lun,'(a)') ''

         ! -------- div(v) --------
         write(lun,'(a)') 'file="divv_line.dat"'
         write(lun,'(a)') 'stats file using 4 prefix "E"  nooutput'
         write(lun,'(a)') 'stats file using (column(4)**2) prefix "E2" nooutput'
         write(lun,'(a)') 'L2 = sprintf("%.3e", sqrt(E2_sum/E2_records))'
         write(lun,'(a)') 'Linf = sprintf("%.3e", E_max)'
         write(lun,'(a)') 'set output "divv.png"'
         write(lun,'(a)') 'set title "div(v) - mid-line slice"'
         write(lun,'(a)') 'set xlabel "x"; set ylabel "value"; set y2label "|error|"'
         write(lun,'(a)') 'unset label; set label 1 sprintf("L2 = %s   Linf = %s", L2, Linf) at graph 0.5,0.92 center front'
         write(lun,'(a)') 'plot file u 1:2 w l  ls 1 t "Analytical", ' // &
      &                 'file u 1:3 w l  ls 2 t "Numerical", '       // &
      &                 'file u 1:($3-$2)          w lp ls 3 t "Residual", ' // &
      &                 'file u 1:(sf_fixed*($3-$2)) w lp ls 4 t "100*Residual", ' // &
      &                 'file u 1:(abs($4)) axes x1y2 w l  ls 5 t "|err|", ' // &
      &                 '0 w l ls 9 notitle'
         write(lun,'(a)') ''

         ! -------- div(B) --------
         write(lun,'(a)') 'file="divb_line.dat"'
         write(lun,'(a)') 'stats file using 4 prefix "E"  nooutput'
         write(lun,'(a)') 'stats file using (column(4)**2) prefix "E2" nooutput'
         write(lun,'(a)') 'L2 = sprintf("%.3e", sqrt(E2_sum/E2_records))'
         write(lun,'(a)') 'Linf = sprintf("%.3e", E_max)'
         write(lun,'(a)') 'set output "divb.png"; set title "div(B) - mid-line slice"'
         write(lun,'(a)') 'unset label; set label 1 sprintf("L2 = %s   Linf = %s", L2, Linf) at graph 0.5,0.92 center front'
         write(lun,'(a)') 'plot file u 1:2 w l  ls 1 t "Analytical", ' // &
      &                 'file u 1:3 w l  ls 2 t "Numerical", '       // &
      &                 'file u 1:($3-$2)          w lp ls 3 t "Residual", ' // &
      &                 'file u 1:(sf_fixed*($3-$2)) w lp ls 4 t "100*Residual", ' // &
      &                 'file u 1:(abs($4)) axes x1y2 w l  ls 5 t "|err|", ' // &
      &                 '0 w l ls 9 notitle'
         write(lun,'(a)') ''

         ! -------- curl(v)_z --------
         write(lun,'(a)') 'file="curlv_z_line.dat"'
         write(lun,'(a)') 'stats file using 4 prefix "E"  nooutput'
         write(lun,'(a)') 'stats file using (column(4)**2) prefix "E2" nooutput'
         write(lun,'(a)') 'L2 = sprintf("%.3e", sqrt(E2_sum/E2_records))'
         write(lun,'(a)') 'Linf = sprintf("%.3e", E_max)'
         write(lun,'(a)') 'set output "curlv_z.png"; set title "curl(v)_z - mid-line slice"'
         write(lun,'(a)') 'unset label; set label 1 sprintf("L2 = %s   Linf = %s", L2, Linf) at graph 0.5,0.92 center front'
         write(lun,'(a)') 'plot file u 1:2 w l  ls 1 t "Analytical", ' // &
      &                 'file u 1:3 w l  ls 2 t "Numerical", '       // &
      &                 'file u 1:($3-$2)          w lp ls 3 t "Residual", ' // &
      &                 'file u 1:(sf_fixed*($3-$2)) w lp ls 4 t "100*Residual", ' // &
      &                 'file u 1:(abs($4)) axes x1y2 w l  ls 5 t "|err|", ' // &
      &                 '0 w l ls 9 notitle'
         write(lun,'(a)') ''

         ! -------- curl(B)_z --------
         write(lun,'(a)') 'file="curlb_z_line.dat"'
         write(lun,'(a)') 'stats file using 4 prefix "E"  nooutput'
         write(lun,'(a)') 'stats file using (column(4)**2) prefix "E2" nooutput'
         write(lun,'(a)') 'L2 = sprintf("%.3e", sqrt(E2_sum/E2_records))'
         write(lun,'(a)') 'Linf = sprintf("%.3e", E_max)'
         write(lun,'(a)') 'set output "curlb_z.png"; set title "curl(B)_z - mid-line slice"'
         write(lun,'(a)') 'unset label; set label 1 sprintf("L2 = %s   Linf = %s", L2, Linf) at graph 0.5,0.92 center front'
         write(lun,'(a)') 'plot file u 1:2 w l  ls 1 t "Analytical", ' // &
      &                 'file u 1:3 w l  ls 2 t "Numerical", '       // &
      &                 'file u 1:($3-$2)          w lp ls 3 t "Residual", ' // &
      &                 'file u 1:(sf_fixed*($3-$2)) w lp ls 4 t "100*Residual", ' // &
      &                 'file u 1:(abs($4)) axes x1y2 w l  ls 5 t "|err|", ' // &
      &                 '0 w l ls 9 notitle'
         write(lun,'(a)') ''

         ! -------- (v x B)_z --------
         write(lun,'(a)') 'file="vcrossb_z_line.dat"'
         write(lun,'(a)') 'stats file using 4 prefix "E"  nooutput'
         write(lun,'(a)') 'stats file using (column(4)**2) prefix "E2" nooutput'
         write(lun,'(a)') 'L2 = sprintf("%.3e", sqrt(E2_sum/E2_records))'
         write(lun,'(a)') 'Linf = sprintf("%.3e", E_max)'
         write(lun,'(a)') 'set output "vcrossb_z.png"; set title "(v x B)_z - mid-line slice"'
         write(lun,'(a)') 'unset label; set label 1 sprintf("L2 = %s   Linf = %s", L2, Linf) at graph 0.5,0.92 center front'
         write(lun,'(a)') 'plot file u 1:2 w l  ls 1 t "Analytical", ' // &
      &                 'file u 1:3 w l  ls 2 t "Numerical", '       // &
      &                 'file u 1:($3-$2)          w lp ls 3 t "Residual", ' // &
      &                 'file u 1:(sf_fixed*($3-$2)) w lp ls 4 t "100*Residual", ' // &
      &                 'file u 1:(abs($4)) axes x1y2 w l  ls 5 t "|err|", ' // &
      &                 '0 w l ls 9 notitle'
         write(lun,'(a)') ''

         ! -------- v·B --------
         write(lun,'(a)') 'file="vdotb_line.dat"'
         write(lun,'(a)') 'stats file using 4 prefix "E"  nooutput'
         write(lun,'(a)') 'stats file using (column(4)**2) prefix "E2" nooutput'
         write(lun,'(a)') 'L2 = sprintf("%.3e", sqrt(E2_sum/E2_records))'
         write(lun,'(a)') 'Linf = sprintf("%.3e", E_max)'
         write(lun,'(a)') 'set output "vdotb.png"; set title "v·B - mid-line slice"'
         write(lun,'(a)') 'unset label; set label 1 sprintf("L2 = %s   Linf = %s", L2, Linf) at graph 0.5,0.92 center front'
         write(lun,'(a)') 'plot file u 1:2 w l  ls 1 t "Analytical", ' // &
      &                 'file u 1:3 w l  ls 2 t "Numerical", '       // &
      &                 'file u 1:($3-$2)          w lp ls 3 t "Residual", ' // &
      &                 'file u 1:(sf_fixed*($3-$2)) w lp ls 4 t "100*Residual", ' // &
      &                 'file u 1:(abs($4)) axes x1y2 w l  ls 5 t "|err|", ' // &
      &                 '0 w l ls 9 notitle'
      ! ------------------ 2D per-figure gnuplot ------------------
      open(newunit=lun, file='test_plot_2d.gp', status='replace')
      write(lun,'(a)') 'set encoding utf8'
      write(lun,'(a)') 'set term pngcairo size 1200,520 font "Helvetica,10" enhanced'
      write(lun,'(a)') 'set pm3d map'
      write(lun,'(a)') 'set palette rgbformulae 33,13,10'
      write(lun,'(a)') 'set view map'
      write(lun,'(a)') 'set colorbox vertical; set cblabel "value"'
      write(lun,'(a,f0.6)') 'z0 = ', z0
      write(lun,'(a)') ''

      ! helper: two panes with shared cbrange + L2/Linf in the overall title
      ! Each pane gets its own title: Analytical (left), Numerical (right)

      ! ---- div(v) ----
      write(lun,'(a)') 'fname = "divv_2d.dat"'
      write(lun,'(a)') 'stats fname u 3 prefix "A"  nooutput'
      write(lun,'(a)') 'stats fname u 4 prefix "N"  nooutput'
      write(lun,'(a)') 'stats fname u (abs($3-$4))   prefix "E"  nooutput'
      write(lun,'(a)') 'stats fname u (($3-$4)**2)   prefix "E2" nooutput'
      write(lun,'(a)') 'cmin = (A_min < N_min) ? A_min : N_min'
      write(lun,'(a)') 'cmax = (A_max > N_max) ? A_max : N_max'
      write(lun,'(a)') 'L2   = sprintf("%.3e", sqrt(E2_sum/E2_records))'
      write(lun,'(a)') 'Linf = sprintf("%.3e", E_max)'
      write(lun,'(a)') 'set cbrange [cmin' // ':' // 'cmax]'
      write(lun,'(a)') 'set output "divv_2d.png"'
      write(lun,'(a)') 'set multiplot layout 1,2 title sprintf("div(v) - 2D slice z=%g  (L2=%s, Linf=%s)", z0, L2, Linf)'
      write(lun,'(a)') 'set xlabel "x"; set ylabel "y"'
      write(lun,'(a)') 'set title "Analytical"'
      write(lun,'(a)') 'splot fname u 1:2:3 w pm3d notitle'
      write(lun,'(a)') 'set title "Numerical"'
      write(lun,'(a)') 'splot fname u 1:2:4 w pm3d notitle'
      write(lun,'(a)') 'unset multiplot'
      write(lun,'(a)') ''

      ! ---- div(B) ----
      write(lun,'(a)') 'fname = "divb_2d.dat"'
      write(lun,'(a)') 'stats fname u 3 prefix "A"  nooutput'
      write(lun,'(a)') 'stats fname u 4 prefix "N"  nooutput'
      write(lun,'(a)') 'stats fname u (abs($3-$4))   prefix "E"  nooutput'
      write(lun,'(a)') 'stats fname u (($3-$4)**2)   prefix "E2" nooutput'
      write(lun,'(a)') 'cmin = (A_min < N_min) ? A_min : N_min'
      write(lun,'(a)') 'cmax = (A_max > N_max) ? A_max : N_max'
      write(lun,'(a)') 'L2   = sprintf("%.3e", sqrt(E2_sum/E2_records))'
      write(lun,'(a)') 'Linf = sprintf("%.3e", E_max)'
      write(lun,'(a)') 'set cbrange [cmin' // ':' // 'cmax]'
      write(lun,'(a)') 'set output "divb_2d.png"'
      write(lun,'(a)') 'set multiplot layout 1,2 title sprintf("div(B) - 2D slice z=%g  (L2=%s, Linf=%s)", z0, L2, Linf)'
      write(lun,'(a)') 'set xlabel "x"; set ylabel "y"'
      write(lun,'(a)') 'set title "Analytical"'
      write(lun,'(a)') 'splot fname u 1:2:3 w pm3d notitle'
      write(lun,'(a)') 'set title "Numerical"'
      write(lun,'(a)') 'splot fname u 1:2:4 w pm3d notitle'
      write(lun,'(a)') 'unset multiplot'
      write(lun,'(a)') ''

      ! ---- curl(v)_z ----
      write(lun,'(a)') 'fname = "curlv_z_2d.dat"'
      write(lun,'(a)') 'stats fname u 3 prefix "A"  nooutput'
      write(lun,'(a)') 'stats fname u 4 prefix "N"  nooutput'
      write(lun,'(a)') 'stats fname u (abs($3-$4))   prefix "E"  nooutput'
      write(lun,'(a)') 'stats fname u (($3-$4)**2)   prefix "E2" nooutput'
      write(lun,'(a)') 'cmin = (A_min < N_min) ? A_min : N_min'
      write(lun,'(a)') 'cmax = (A_max > N_max) ? A_max : N_max'
      write(lun,'(a)') 'L2   = sprintf("%.3e", sqrt(E2_sum/E2_records))'
      write(lun,'(a)') 'Linf = sprintf("%.3e", E_max)'
      write(lun,'(a)') 'set cbrange [cmin' // ':' // 'cmax]'
      write(lun,'(a)') 'set output "curlv_z_2d.png"'
      write(lun,'(a)') 'set multiplot layout 1,2 title sprintf("curl(v)_z - 2D slice z=%g  (L2=%s, Linf=%s)", z0, L2, Linf)'
      write(lun,'(a)') 'set xlabel "x"; set ylabel "y"'
      write(lun,'(a)') 'set title "Analytical"'
      write(lun,'(a)') 'splot fname u 1:2:3 w pm3d notitle'
      write(lun,'(a)') 'set title "Numerical"'
      write(lun,'(a)') 'splot fname u 1:2:4 w pm3d notitle'
      write(lun,'(a)') 'unset multiplot'
      write(lun,'(a)') ''

      ! ---- curl(B)_z ----
      write(lun,'(a)') 'fname = "curlb_z_2d.dat"'
      write(lun,'(a)') 'stats fname u 3 prefix "A"  nooutput'
      write(lun,'(a)') 'stats fname u 4 prefix "N"  nooutput'
      write(lun,'(a)') 'stats fname u (abs($3-$4))   prefix "E"  nooutput'
      write(lun,'(a)') 'stats fname u (($3-$4)**2)   prefix "E2" nooutput'
      write(lun,'(a)') 'cmin = (A_min < N_min) ? A_min : N_min'
      write(lun,'(a)') 'cmax = (A_max > N_max) ? A_max : N_max'
      write(lun,'(a)') 'L2   = sprintf("%.3e", sqrt(E2_sum/E2_records))'
      write(lun,'(a)') 'Linf = sprintf("%.3e", E_max)'
      write(lun,'(a)') 'set cbrange [cmin' // ':' // 'cmax]'
      write(lun,'(a)') 'set output "curlb_z_2d.png"'
      write(lun,'(a)') 'set multiplot layout 1,2 title sprintf("curl(B)_z - 2D slice z=%g  (L2=%s, Linf=%s)", z0, L2, Linf)'
      write(lun,'(a)') 'set xlabel "x"; set ylabel "y"'
      write(lun,'(a)') 'set title "Analytical"'
      write(lun,'(a)') 'splot fname u 1:2:3 w pm3d notitle'
      write(lun,'(a)') 'set title "Numerical"'
      write(lun,'(a)') 'splot fname u 1:2:4 w pm3d notitle'
      write(lun,'(a)') 'unset multiplot'
      write(lun,'(a)') ''

      ! ---- (v x B)_z ----
      write(lun,'(a)') 'fname = "vcrossb_z_2d.dat"'
      write(lun,'(a)') 'stats fname u 3 prefix "A"  nooutput'
      write(lun,'(a)') 'stats fname u 4 prefix "N"  nooutput'
      write(lun,'(a)') 'stats fname u (abs($3-$4))   prefix "E"  nooutput'
      write(lun,'(a)') 'stats fname u (($3-$4)**2)   prefix "E2" nooutput'
      write(lun,'(a)') 'cmin = (A_min < N_min) ? A_min : N_min'
      write(lun,'(a)') 'cmax = (A_max > N_max) ? A_max : N_max'
      write(lun,'(a)') 'L2   = sprintf("%.3e", sqrt(E2_sum/E2_records))'
      write(lun,'(a)') 'Linf = sprintf("%.3e", E_max)'
      write(lun,'(a)') 'set cbrange [cmin' // ':' // 'cmax]'
      write(lun,'(a)') 'set output "vcrossb_z_2d.png"'
      write(lun,'(a)') 'set multiplot layout 1,2 title sprintf("(v x B)_z - 2D slice z=%g  (L2=%s, Linf=%s)", z0, L2, Linf)'
      write(lun,'(a)') 'set xlabel "x"; set ylabel "y"'
      write(lun,'(a)') 'set title "Analytical"'
      write(lun,'(a)') 'splot fname u 1:2:3 w pm3d notitle'
      write(lun,'(a)') 'set title "Numerical"'
      write(lun,'(a)') 'splot fname u 1:2:4 w pm3d notitle'
      write(lun,'(a)') 'unset multiplot'
      write(lun,'(a)') ''

      ! ---- v·B ----
      write(lun,'(a)') 'fname = "vdotb_2d.dat"'
      write(lun,'(a)') 'stats fname u 3 prefix "A"  nooutput'
      write(lun,'(a)') 'stats fname u 4 prefix "N"  nooutput'
      write(lun,'(a)') 'stats fname u (abs($3-$4))   prefix "E"  nooutput'
      write(lun,'(a)') 'stats fname u (($3-$4)**2)   prefix "E2" nooutput'
      write(lun,'(a)') 'cmin = (A_min < N_min) ? A_min : N_min'
      write(lun,'(a)') 'cmax = (A_max > N_max) ? A_max : N_max'
      write(lun,'(a)') 'L2   = sprintf("%.3e", sqrt(E2_sum/E2_records))'
      write(lun,'(a)') 'Linf = sprintf("%.3e", E_max)'
      write(lun,'(a)') 'set cbrange [cmin' // ':' // 'cmax]'
      write(lun,'(a)') 'set output "vdotb_2d.png"'
      write(lun,'(a)') 'set multiplot layout 1,2 title sprintf("v·B - 2D slice z=%g  (L2=%s, Linf=%s)", z0, L2, Linf)'
      write(lun,'(a)') 'set xlabel "x"; set ylabel "y"'
      write(lun,'(a)') 'set title "Analytical"'
      write(lun,'(a)') 'splot fname u 1:2:3 w pm3d notitle'
      write(lun,'(a)') 'set title "Numerical"'
      write(lun,'(a)') 'splot fname u 1:2:4 w pm3d notitle'
      write(lun,'(a)') 'unset multiplot'

         close(lun)

         call printinfo('Wrote: *_line.dat, *_2d.dat, test_plot_1d.gp, test_plot_2d.gp', V_INFO)
         call printinfo('1D: gnuplot test_plot_1d.gp  -> divv.png, divb.png, curlv_z.png, curlb_z.png, vcrossb_z.png, vdotb.png', V_INFO)
         call printinfo('2D: gnuplot test_plot_2d.gp     -> *_2d.png (Analytical | Numerical)', V_INFO)
      endif

   end subroutine dump_vec_gnuplot

end module initproblem
