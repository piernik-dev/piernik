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

   use constants, only: cbuff_len

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

   enum, bind(C)
      enumerator :: DIVV, DIVB, CURLV_Z, CURLB_Z, V_DOT_B, VCROSSB_Z
   end enum
   character(len=cbuff_len), dimension(DIVV:VCROSSB_Z), parameter :: titles = [ 'div(v)   ', 'div(B)   ', 'curl(v)_z', 'curl(B)_z', 'v·B     ', '(v×B)_z ' ]

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

      call all_cg%reg_var(divvnum,            restart_mode = AT_IGNORE)
      call all_cg%reg_var(divvana,            restart_mode = AT_IGNORE)
      call all_cg%reg_var(divbnum,            restart_mode = AT_IGNORE)
      call all_cg%reg_var(divbana,            restart_mode = AT_IGNORE)
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

      do p = 1, flind%fluids

         fl => flind%all_fluids(p)%fl

         cgl => leaves%first
         do while (associated(cgl))
            cg => cgl%cg

            do j = cg%lhn(ydim, LO), cg%lhn(ydim, HI)
               yj = cg%y(j)
               do i = cg%lhn(xdim, LO), cg%lhn(xdim, HI)
                  xi = cg%x(i)
                  do k = cg%lhn(zdim, LO), cg%lhn(zdim, HI)
                     zk = cg%z(k)
                     sx = sin(xi);  cx = cos(xi)
                     sy = sin(yj);  cy = cos(yj)
                     cg%u(fl%idn, i, j, k) = 1.0
                     cg%u(fl%imx, i, j, k) =  cy*sx
                     cg%u(fl%imy, i, j, k) = -cx*sy
                     cg%u(fl%imz, i, j, k) = 0.0

                     cg%q(qna%ind(divvana))%arr(i, j, k) = 0.0

                     if (fl%has_energy) then
                        cg%u(fl%ien, i, j, k) = 1.0
                        cg%u(fl%ien, i, j, k) = cg%u(fl%ien, i, j, k) + ekin(cg%u(fl%imx, i, j, k), cg%u(fl%imy, i, j, k), cg%u(fl%imz, i, j, k), cg%u(fl%idn, i, j, k))

                        if (fl%is_magnetized) then
                           cg%q(qna%ind(divbana))%arr(i, j, k) = 0.0
                           cg%b(xdim, i, j, k) = -sy*cx
                           cg%b(ydim, i, j, k) =  sx*cy
                           cg%b(zdim, i, j, k) =  0.0
                           cg%u(fl%ien, i, j, k) = cg%u(fl%ien, i, j, k) + emag(cg%b(xdim, i, j, k), cg%b(ydim, i, j, k), cg%b(zdim, i, j, k))
                        endif
                     endif

                  enddo
               enddo
            enddo

            cgl => cgl%nxt
         enddo
      enddo

      fl => flind%all_fluids(1)%fl
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         cg%q(qna%ind(divvnum))%arr  = cg%get_divergence(ord=order, iw=wna%fi, vec=[fl%imx, fl%imy, fl%imz])
         cg%q(qna%ind(divbnum))%arr  = cg%get_divergence(ord=order, iw=wna%bi, vec=[xdim, ydim, zdim])
         cg%w(wna%ind(curlbnum))%arr = cg%get_curl      (ord=order, iw=wna%bi, vec=[xdim, ydim, zdim])
         cg%w(wna%ind(curlvnum))%arr = cg%get_curl      (ord=order, iw=wna%fi, vec=[fl%imx, fl%imy, fl%imz])
         cg%w(wna%ind(vcrossb))%arr  = cg%cross(wna%fi, wna%bi, [fl%imx, fl%imy, fl%imz], [xdim, ydim, zdim])
         cg%q(qna%ind(vdotb))%arr    = cg%dot  (wna%fi, wna%bi, [fl%imx, fl%imy, fl%imz], [xdim, ydim, zdim])

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

      if (.not. qna%exists(divvana))  call die("[initproblem:user_out] cannot find divvana")
      if (.not. qna%exists(divvnum))  call die("[initproblem:user_out] cannot find divvnum")
      if (.not. qna%exists(divbana))  call die("[initproblem:user_out] cannot find divbana")
      if (.not. qna%exists(divbnum))  call die("[initproblem:user_out] cannot find divbnum")
      if (.not. wna%exists(curlvnum)) call die("[initproblem:user_out] cannot find curlvnum")
      if (.not. wna%exists(curlbnum)) call die("[initproblem:user_out] cannot find curlbnum")
      if (.not. wna%exists(vcrossb))  call die("[initproblem:user_out] cannot find vcrossb")
      if (.not. qna%exists(vdotb))    call die("[initproblem:user_out] cannot find vdotb")

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

      use constants,  only: V_ESSENTIAL, V_INFO
      use dataio_pub, only: msg, printinfo
      use mpisetup,   only: master

      implicit none

      real :: ssum, smax
      integer :: comp

      if (master) then
         write(msg,'("--- Diagnostics (global volume-weighted L2 errors, ", i1, a2, "-order of operators) ---")') order, merge("nd", "th", order == 2)
         call printinfo(msg, V_INFO)
      endif

      ! Main execution
      do comp = DIVV, VCROSSB_Z
         call calc_error(ssum, smax)
         if (master) then
            write(msg,'("||",a10,"||_err : L2_vw=",1pe12.5,"  Linf=",1pe12.5)') trim(titles(comp)), sqrt(ssum), smax
            call printinfo(msg, V_ESSENTIAL)
         endif
      enddo

   contains

      subroutine calc_error(err2, err_max)

         use allreduce,        only: piernik_MPI_Allreduce
         use cg_list,          only: cg_list_element
         use cg_leaves,        only: leaves
         use constants,        only: pSUM, pMAX, xdim, zdim
         use dataio_pub,       only: die
         use grid_cont,        only: grid_container
         use named_array_list, only: qna, wna

         implicit none

         real, intent(out) :: err2, err_max

         type(cg_list_element), pointer :: cgl
         type(grid_container),  pointer :: cg
         real :: lsum, lmax, xi, yj, val
         integer :: i, j, k
         real, dimension(xdim:zdim) :: vec

         err2 = 0.0
         err_max = 0.0

         cgl => leaves%first
         do while (associated(cgl))
            cg => cgl%cg

            lsum = 0.0
            lmax = 0.0
            do k = cg%ks, cg%ke
               do j = cg%js, cg%je
                  yj = cg%y(j)
                  associate(sy => sin(yj), cy => cos(yj))
                  do i = cg%is, cg%ie
                     xi = cg%x(i)
                     associate (sx => sin(xi), cx => cos(xi))
                     select case (comp)
                        case (DIVV)
                           val = cg%q(qna%ind(divvnum))%arr(i, j, k)
                           lmax = max(lmax, abs(val))
                           lsum = lsum + val*val
                        case (DIVB)
                           val = cg%q(qna%ind(divbnum))%arr(i, j, k)
                           lmax = max(lmax, abs(val))
                           lsum = lsum + val*val
                        case (CURLV_Z)
                           vec(:) = cg%w(wna%ind(curlvnum))%arr(xdim:zdim, i, j, k) - [ 0., 0., 2.0*sx*sy ]
                           lmax = max(lmax, norm2(vec))
                           lsum = lsum + norm2(vec)**2
                        case (CURLB_Z)
                           vec(:) = cg%w(wna%ind(curlbnum))%arr(xdim:zdim, i, j, k) - [ 0., 0., 2.0*cx*cy ]
                           lmax = max(lmax, norm2(vec))
                           lsum = lsum + norm2(vec)**2
                        case (V_DOT_B)
                           val = cg%q(qna%ind(vdotb))%arr(i, j, k) + 2.0*sx*cx*sy*cy
                           lmax = max(lmax, abs(val))
                           lsum = lsum + val*val
                        case (VCROSSB_Z)
                           vec(:) = cg%w(wna%ind(vcrossb))%arr(xdim:zdim, i, j, k) - [ 0., 0., (sx*sx*cy*cy - cx*cx*sy*sy) ]
                           lmax = max(lmax, norm2(vec))
                           lsum = lsum + norm2(vec)**2
                        case default
                           call die("[initproblem:verify_test:calc_error] Unknown case")
                     end select
                     end associate
                  enddo
                  end associate
               enddo
            enddo
            err2 = err2 + lsum * cg%dvol
            err_max = max(err_max, lmax)

            cgl => cgl%nxt
         enddo

         call piernik_MPI_Allreduce(err2, pSUM)
         call piernik_MPI_Allreduce(err_max, pMAX)

      end subroutine calc_error

   end subroutine verify_test

   subroutine dump_vec_gnuplot

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: V_INFO, zdim, cbuff_len
      use dataio_pub,       only: printinfo, die, warn, msg
      use grid_cont,        only: grid_container
      use mpisetup,         only: master
      use named_array_list, only: qna, wna

      implicit none

      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
      integer :: i, j, j0, k0, lun, dset
      real :: x, y, z0, va, vn
      logical :: first_cg

      character(len=cbuff_len), dimension(DIVV:VCROSSB_Z), parameter :: fnames     = [ "divv     ", "divb     ", "curlv_z  ", "curlb_z  ", "vdotb    ", "vcrossb_z" ]

      call verify_test

      if (.not. master) then
         call warn("[initproblem:dump_vec_gnuplot] This is serial routine. The dumps will be incomplete – limited to master process")
         return
      endif

      call printinfo('Writing .dat files for the gnuplot. 1D is along x for y = midpoint of the domain', V_INFO)

      first_cg = .true.
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         ! Mid-line in y, mid-plane in z (works nice only on single-block domain)
         j0 = (cg%js + cg%je)/2
         k0 = (cg%ks + cg%ke)/2
         y  = cg%y(j0)
         associate (sy => sin(y), cy => cos(y))
         z0 = cg%z(k0)

         ! 1D slices
         do dset = DIVV, VCROSSB_Z
            ! delete old file before writing first cg
            open(newunit=lun, file=trim(fnames(dset)) // '_line.dat', status='unknown', access=trim(merge('sequential', 'append    ', first_cg)))
            do i = cg%is, cg%ie
               x  = cg%x(i)
               associate (sx => sin(x), cx => cos(x))
               select case (dset)
                  case (DIVV)
                     va = 0.0
                     vn = cg%q(qna%ind(divvnum))%arr(i, j0, k0)
                  case (DIVB)
                     va = 0.0
                     vn = cg%q(qna%ind(divbnum))%arr(i, j0, k0)
                  case (CURLV_Z)
                     va = 2.0 * sx * sy
                     vn = cg%w(wna%ind(curlvnum))%arr(zdim, i, j0, k0)
                  case (CURLB_Z)
                     va = 2.0 * cx * cy
                     vn = cg%w(wna%ind(curlbnum))%arr(zdim, i, j0, k0)
                  case (V_DOT_B)
                     va = -2.0 * sx * cx * sy * cy
                     vn = cg%q(qna%ind(vdotb))%arr(i, j0, k0)
                  case (VCROSSB_Z)
                     va = sx * sx * cy * cy - cx * cx * sy * sy
                     vn = cg%w(wna%ind(vcrossb))%arr(zdim, i, j0, k0)
                  case default
                     call die("[initproblem:dump_vec_gnuplot] unknown 1D case")
               end select
               end associate
               write(lun,'(4(1pe23.15))') x, va, vn, vn - va
            enddo
            write(lun,'(a)') ''
            close(lun)
         enddo
         end associate

         ! 2D maps (x,y at fixed z) --
         ! columns: x y ana num
         do dset = DIVV, VCROSSB_Z
            open(newunit=lun, file=trim(fnames(dset)) // '_2d.dat', status='unknown', access=trim(merge('sequential', 'append    ', first_cg)))
            do j = cg%js, cg%je
               y = cg%y(j)
               associate (sy => sin(y), cy => cos(y))
               do i = cg%is, cg%ie
                  x  = cg%x(i)
                  associate (sx => sin(x), cx => cos(x))
                  select case (dset)
                     case (DIVV)
                        va = 0.0
                        vn = cg%q(qna%ind(divvnum))%arr(i, j, k0)
                     case (DIVB)
                        va = 0.0
                        vn = cg%q(qna%ind(divbnum))%arr(i, j, k0)
                     case (CURLV_Z)
                        va = 2.0 * sx * sy
                        vn = cg%w(wna%ind(curlvnum))%arr(zdim, i, j, k0)
                     case (CURLB_Z)
                        va = 2.0 * cx * cy
                        vn = cg%w(wna%ind(curlbnum))%arr(zdim, i, j, k0)
                     case (V_DOT_B)
                        va = -2.0 * sx * cx * sy * cy
                        vn = cg%q(qna%ind(vdotb))%arr(i, j, k0)
                     case (VCROSSB_Z)
                        va = sx * sx * cy * cy - cx * cx * sy * sy
                        vn = cg%w(wna%ind(vcrossb))%arr(zdim, i, j, k0)
                     case default
                        call die("[initproblem:dump_vec_gnuplot] unknown 2D case")
                  end select
                  end associate
                  write(lun,'(4(1pe23.15))') x, y, va, vn
               enddo
               end associate
               write(lun,'(a)') ''
            enddo
            write(lun,'(a)') ''
            close(lun)
         enddo

         cgl => cgl%nxt
         first_cg = .false.
      enddo

      ! 1D script
      open(newunit=lun, file='test_plot_1d.gp', status='unknown')
      write(lun,'(a)') 'set encoding utf8'
      write(lun,'(a)') 'print "Creating 1D plots for VecOps test"'
      write(lun,'(a)') 'set term pngcairo size 1100,700 font "Helvetica,10" enhanced'
      write(lun,'(a)') 'set grid xtics ytics mxtics mytics back lc rgb "#cccccc"'
      write(lun,'(a)') 'set mxtics 2; set mytics 2'
      write(lun,'(a)') 'set key at graph 0.5,0.98 center top horizontal opaque samplen 2 spacing 1.2'
      write(lun,'(a)') 'set tics out; set format y "%1.2g"; set format y2 "%1.2g"; set y2tics'
      write(lun,'(a)') 'set style line 1 lc rgb "#7b3294" lw 2.2'                  ! Analytical
      write(lun,'(a)') 'set style line 2 lc rgb "#008837" lw 2.2'                  ! Numerical
      write(lun,'(a)') 'set style line 3 lc rgb "#E69F00" lw 2.6 dt 1 pt 5 ps 0.8' ! Residual
      write(lun,'(a)') 'set style line 4 lc rgb "#E69F00" lw 2.6 dt 3 pt 5 ps 0.8' ! 100×Residual
      write(lun,'(a)') 'set style line 5 lc rgb "#0571B0" lw 1.0'                  ! |err| (y2)
      write(lun,'(a)') 'set style line 9 lc rgb "#BBBBBB" lw 1.0 dt 3'             ! zero line
      write(lun,'(a)') 'sf_fixed = 100.0'                                            ! scale for 100×Residual
      write(lun,'(a)') ''

      do dset = DIVV, VCROSSB_Z
         write(lun,'(a)') 'file="' // trim(fnames(dset)) // '_line.dat"'
         write(lun,'(a)') 'stats file using 4 prefix "E"  nooutput'
         write(lun,'(a)') 'stats file using (column(4)**2) prefix "E2" nooutput'
         write(lun,'(a)') 'L2 = sprintf("%.3e", sqrt(E2_sum/E2_records))'
         write(lun,'(a)') 'Linf = sprintf("%.3e", E_max)'
         write(lun,'(a)') 'set output "' // trim(fnames(dset)) // '.png"'
         write(lun,'(a)') 'set title "' // trim(titles(dset)) // ' – mid-line slice"'
         write(lun,'(a)') 'set xlabel "x"; set ylabel "value"; set y2label "|error|"'
         write(lun,'(a)') 'unset label; set label 1 sprintf("L_2 = %s   L_∞ = %s", L2, Linf) at graph 0.5,0.92 center front'
         write(lun,'(a)') 'plot file u 1:2 w l  ls 1 t "Analytical", ' // &
      &                 'file u 1:3 w l  ls 2 t "Numerical", '       // &
      &                 'file u 1:($3-$2)          w lp ls 3 t "Residual", ' // &
      &                 'file u 1:(sf_fixed*($3-$2)) w lp ls 4 t "100*Residual", ' // &
      &                 'file u 1:(abs($4)) axes x1y2 w l  ls 5 t "|err|", ' // &
      &                 '0 w l ls 9 notitle'
         write(lun,'(a)') ''
      enddo
      close(lun)

      ! 2D script
      open(newunit=lun, file='test_plot_2d.gp', status='replace')
      write(lun,'(a)') 'set encoding utf8'
      write(lun,'(a)') 'set term pngcairo size 1200,520 font "Helvetica,10" enhanced'
      write(lun,'(a)') 'set pm3d map'
      write(lun,'(a)') 'set palette rgbformulae 33,13,10'
      write(lun,'(a)') 'set view map'
      write(lun,'(a)') 'set colorbox vertical; #set cblabel "value"'  ! cblabel would need an offset and rotation, otherwise it collides with ylabel of the right panel
      write(lun,'(a,f0.6)') 'z0 = ', z0
      write(lun,'(a)') ''

      write(lun,'(a)') 'print "Creating 2D plots for VecOps test"'
      do dset = DIVV, VCROSSB_Z
         write(lun,'(a)') 'fname = "' // trim(fnames(dset)) // '_2d.dat"'
         write(lun,'(a)') 'stats fname u 3 prefix "A"  nooutput'
         write(lun,'(a)') 'stats fname u 4 prefix "N"  nooutput'
         write(lun,'(a)') 'stats fname u (abs($3-$4))   prefix "E"  nooutput'
         write(lun,'(a)') 'stats fname u (($3-$4)**2)   prefix "E2" nooutput'
         write(lun,'(a)') 'cmin = (A_min < N_min) ? A_min : N_min'
         write(lun,'(a)') 'cmax = (A_max > N_max) ? A_max : N_max'
         write(lun,'(a)') 'L2   = sprintf("%.3e", sqrt(E2_sum/E2_records))'
         write(lun,'(a)') 'Linf = sprintf("%.3e", E_max)'
         write(lun,'(a)') '#set cbrange [cmin' // ':' // 'cmax]'
         write(lun,'(a)') 'set output "' // trim(fnames(dset)) // '_2d.png"'
         write(lun,'(a)') 'set multiplot layout 1,2 title sprintf("' // trim(titles(dset)) // ' – 2D slice z = %g  (L_2 = %s, L_∞ = %s)", z0, L2, Linf)'
         write(lun,'(a)') 'set xlabel "x"; set ylabel "y"'
         write(lun,'(a)') 'set title "Numerical - Analytical"'
         write(lun,'(a)') 'splot fname u 1:2:($4-$3) w pm3d notitle'
         write(lun,'(a)') 'set title "Numerical"'
         write(lun,'(a)') 'splot fname u 1:2:4 w pm3d notitle'
         write(lun,'(a)') 'unset multiplot'
         write(lun,'(a)') ''
      enddo
      close(lun)

      call printinfo('Wrote: *_line.dat, *_2d.dat, test_plot_1d.gp, test_plot_2d.gp', V_INFO)
      msg = '1D: gnuplot test_plot_1d.gp  ->'
      do dset = DIVV, VCROSSB_Z
         write(msg, '(a)') trim(msg) // ' ' // trim(fnames(dset)) // '.png'
      enddo
      call printinfo(trim(msg), V_INFO)
      call printinfo('2D: gnuplot test_plot_2d.gp  ->  *_2d.png (Numerical - Analytical | Numerical)', V_INFO)

   end subroutine dump_vec_gnuplot

end module initproblem
