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
!! \brief Calculate the lovely shape of the Mandelbrot set. Refine the
!! set to follow the interesting details.
!!
!! The Mandelbrot problem is intended for stress testing of the AMR
!! subsystem. There are more efficient fractal generators, but
!! one may consider some fun ideas:
!! * Use multiprecision or implement a custom fixed-point representation
!!   optimized for these calculations.
!! * Use speedup tricks like these in fast deep zoom programs.
!! * Detect the interiors of minibrots for further speedups. This already
!!   works to some extent because they are not refined too much.
!<

module initproblem

   use constants, only: dsetnamelen, cbuff_len, FP_REAL, FP_DOUBLE, FP_EXT, FP_QUAD

   implicit none

   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   ! namelist parameters
   integer(kind=4) :: maxiter             !< Maximum number of iterations
   logical :: smooth_map                  !< Try continuous colouring
   real :: ref_thr                        !< threshold for refining a grid
   character(len=cbuff_len) :: precision  !< precision of Mandelbrot calculations
   logical :: log_polar                   !< Use polar mapping around x_center + i * y_center
   character(len=cbuff_len) :: x_center   !< x-coordinate of the center (crucial for polar mode)
   character(len=cbuff_len) :: y_center   !< y-coordinate of the center (crucial for polar mode)
   real :: c_polar                        !< correct colouring with x-coordinate multiplied by this factor

   namelist /PROBLEM_CONTROL/  maxiter, smooth_map, log_polar, x_center, y_center, c_polar, ref_thr, precision

   ! other private data
   character(len=dsetnamelen), parameter :: mand_n = "mand", re_n = "real", imag_n = "imag"

   integer :: prec  ! decoded precision level
   real(kind=FP_QUAD) :: xcq, ycq  ! decoded central coordinates
   enum, bind(C)
      enumerator :: FP_QCMPLX = maxval([FP_REAL, FP_DOUBLE, FP_EXT, FP_QUAD]) + 1 !, FP_2FLOAT, FP_2DOUBLE
   end enum

contains

!> Set up some user hooks

   subroutine problem_pointers

      use dataio_user, only: user_reg_var_restart
#ifdef HDF5
      use dataio_user, only: user_vars_hdf5
#endif /* HDF5 */
      use user_hooks,  only: problem_refine_derefine

      implicit none

      user_reg_var_restart    => register_user_var
#ifdef HDF5
      user_vars_hdf5          => mand_vars
#endif /* HDF5 */
      problem_refine_derefine => mark_the_set

   end subroutine problem_pointers

!> \brief Read runtime parameters

   subroutine read_problem_par

      use bcast,      only: piernik_MPI_Bcast
      use constants,  only: ydim, LO, HI, dpi, INVALID
      use dataio_pub, only: warn, die, nh
      use domain,     only: dom
      use mpisetup,   only: cbuff, lbuff, ibuff, rbuff, master, slave

      implicit none

      ! namelist default parameter values
      maxiter = 100
      smooth_map = .true.
      log_polar = .false.
      x_center = "0."
      y_center = "0."
      precision = "double"
      c_polar = 0.
      ref_thr = 1.

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

         ibuff(1) = maxiter

         lbuff(1) = smooth_map
         lbuff(2) = log_polar

         rbuff(1) = ref_thr
         rbuff(2) = c_polar

         cbuff(1) = x_center
         cbuff(2) = y_center
         cbuff(3) = precision

      endif

      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(lbuff)
      call piernik_MPI_Bcast(rbuff)
      call piernik_MPI_Bcast(cbuff, cbuff_len)

      if (slave) then

         maxiter    = ibuff(1)

         smooth_map = lbuff(1)
         log_polar  = lbuff(2)

         ref_thr    = rbuff(1)
         c_polar    = rbuff(2)

         x_center   = cbuff(1)
         y_center   = cbuff(2)
         precision  = cbuff(3)

      endif

      if (any(dom%has_dir(:) .neqv. [ .true., .true., .false. ])) &
           call die("[initproblem:read_problem_par] Mandelbrot is supposed to be run only in the XY plane without a Z direction")

      select case (trim(precision))
         case ("single", "float")
            prec = FP_REAL
         case ("double")
            prec = FP_DOUBLE
         case ("extended", "long")
            prec = FP_EXT
         case ("quad")
            prec = FP_QUAD
         case ("quad complex")
            prec = FP_QCMPLX
         case default
            call die("[initproblem:read_problem_par] precision must be single, double, extended, or quad")
            prec = INVALID
      end select

      if (log_polar .and. master) then
         if (dom%edge(ydim, HI) - dom%edge(ydim, LO) < 0.999*dpi) call warn("[initproblem:read_problem_par] not covering full angle")
         if (dom%edge(ydim, HI) - dom%edge(ydim, LO) > 1.001*dpi) call warn("[initproblem:read_problem_par] covering more than full angle")
      endif

      read(x_center, *) xcq
      read(y_center, *) ycq

      call register_user_var

   end subroutine read_problem_par

!> \brief Calculate the Mandelbrot iterations for all leaves

   subroutine problem_initial_conditions

      use cg_list,          only: cg_list_element
      use cg_leaves,        only: leaves
      use dataio_pub,       only: warn
      use grid_cont,        only: grid_container
      use fluidindex,       only: iarr_all_dn
      use mpisetup,         only: master
      use named_array_list, only: qna, wna

      implicit none

      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg
      integer :: i, j, k
      real, dimension(:,:,:), pointer :: mand, r__l, imag

      ! Create the initial density arrays
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         if (.not. cg%is_old) then

            call cg%set_constant_b_field([0., 0., 0.])
            cg%u(:, :, :, :) = 0.

            mand => cg%q(qna%ind(mand_n))%arr
            r__l => cg%q(qna%ind(re_n  ))%arr
            imag => cg%q(qna%ind(imag_n))%arr
            if (.not. associated(mand) .or. .not. associated(r__l) .or. .not. associated(imag)) then
               if (master) call warn("[initproblem:problem_initial_conditions] Cannot store the set")
               return
            endif

            do k = cg%ks, cg%ke
               do j = cg%js, cg%je
                  do i = cg%is, cg%ie
                     call calculate_mandelbrot(cg%x(i), cg%y(j), mand(i,j,k), r__l(i,j,k), imag(i,j,k))
                  enddo
               enddo
            enddo

         endif
         cgl => cgl%nxt
      enddo

      call leaves%qw_copy(qna%ind(mand_n), wna%fi, iarr_all_dn(1)) ! prevent spurious FP exceptions

   end subroutine problem_initial_conditions

! Fortran 2018 has no assumed-kind polymorphism for this calculation. The
! explicit blocks below keep the selected kind visible for teaching purposes.

! TODO: Implement a double-double approach for comparison with built-in quad.

   subroutine calculate_mandelbrot(xx, yy, mand, real_z, imag_z)

      use constants,  only: FP_REAL, FP_DOUBLE, FP_EXT, FP_QUAD
      use dataio_pub, only: die

      implicit none

      real, intent(in) :: xx, yy
      real, intent(out) :: mand, real_z, imag_z

      integer :: nit
      real :: rnit, r, x, y
      real, parameter :: bailout2 = 10., min_log_mand = 0.1

      nit = 1
      if (log_polar) then
         r = 10.**xx
         x = r*cos(yy)
         y = r*sin(yy)
      else
         x = xx
         y = yy
      endif

      select case (prec)
         case (FP_REAL)
            block
               real(kind=FP_REAL) :: zx, zy, zt, cx, cy

               cx = real(xcq + x, kind=kind(cx))
               cy = real(ycq + y, kind=kind(cy))

               zx = cx
               zy = cy
               do while (zx*zx + zy*zy < bailout2 .and. nit < maxiter)
                  zt = zx*zx - zy*zy + cx
                  zy = 2._FP_REAL*zx*zy + cy
                  zx = zt
                  nit = nit + 1
               enddo
               x = real(zx, kind=kind(x))
               y = real(zy, kind=kind(y))

            end block
         case (FP_DOUBLE)
            block
               real(kind=FP_DOUBLE) :: zx, zy, zt, cx, cy

               cx = real(xcq + x, kind=kind(cx))
               cy = real(ycq + y, kind=kind(cy))

               zx = cx
               zy = cy
               do while (zx*zx + zy*zy < bailout2 .and. nit < maxiter)
                  zt = zx*zx - zy*zy + cx
                  zy = 2.*zx*zy + cy
                  zx = zt
                  nit = nit + 1
               enddo
               x = real(zx, kind=kind(x))
               y = real(zy, kind=kind(y))

            end block
         case (FP_EXT)
            block
               real(kind=FP_EXT) :: zx, zy, zt, cx, cy

               cx = real(xcq + x, kind=kind(cx))
               cy = real(ycq + y, kind=kind(cy))

               zx = cx
               zy = cy
               do while (zx*zx + zy*zy < bailout2 .and. nit < maxiter)
                  zt = zx*zx - zy*zy + cx
                  zy = 2.*zx*zy + cy
                  zx = zt
                  nit = nit + 1
               enddo
               x = real(zx, kind=kind(x))
               y = real(zy, kind=kind(y))

            end block
         case (FP_QUAD)
            block
               real(kind=FP_QUAD) :: zx, zy, zt, cx, cy

               cx = real(xcq + x, kind=kind(cx))
               cy = real(ycq + y, kind=kind(cy))

               zx = cx
               zy = cy
               do while (zx*zx + zy*zy < bailout2 .and. nit < maxiter)
                  zt = zx*zx - zy*zy + cx
                  zy = 2.*zx*zy + cy
                  zx = zt
                  nit = nit + 1
               enddo
               x = real(zx, kind=kind(x))
               y = real(zy, kind=kind(y))

            end block
         case (FP_QCMPLX)
            ! Using native complex type is more compact but also about 10% slower
            block
               complex(kind=FP_QUAD) :: z, c

               c = complex(xcq + x, ycq + y)
               z = c
               do while (real(z)**2 + imag(z)**2 < bailout2 .and. nit < maxiter)
                  z = z*z + c
                  nit = nit + 1
               enddo
               x = real(z, kind=kind(x))
               y = real(imag(z), kind=kind(y))

            end block
         case default
            call die("[initproblem:calculate_mandelbrot] non-implemented precision")
            x = 0.
            y = 0.
            nit = 0
      end select

      rnit = nit
      if (smooth_map .and. x*x + y*y > bailout2) rnit = rnit + 1 - log(log(sqrt(x*x + y*y)))/log(2.)
      if (nit >= maxiter) then
         mand = min_log_mand
      else
         mand = max(min_log_mand, log(max(rnit, min_log_mand)) + c_polar * xx)
      endif
      real_z = x
      imag_z = y

   end subroutine calculate_mandelbrot

!> \brief Add fields for iteration count and real and imaginary coordinate of the point at the end of iterations

   subroutine register_user_var

      use cg_list_global, only: all_cg
      use constants,      only: AT_NO_B

      implicit none

      call all_cg%reg_var(mand_n, restart_mode = AT_NO_B)
      call all_cg%reg_var(re_n,   restart_mode = AT_NO_B)
      call all_cg%reg_var(imag_n, restart_mode = AT_NO_B)

   end subroutine register_user_var

!> \brief Make some derived variables available for output

   subroutine mand_vars(var, tab, ierrh, cg)

      use grid_cont,   only: grid_container
      use named_array_list, only: qna

      implicit none

      character(len=*),              intent(in)    :: var
      real, dimension(:,:,:),        intent(inout) :: tab
      integer,                       intent(inout) :: ierrh
      type(grid_container), pointer, intent(in)    :: cg

      ierrh = 0
      select case (trim(var))
         case ("distance", "dist") ! Supply the alternative name to comply with the old 4-letter limit
            tab(:,:,:) = log(sqrt(cg%q(qna%ind(re_n))%span(cg%ijkse)**2 + cg%q(qna%ind(imag_n))%span(cg%ijkse)**2))
         case ("angle", "ang")
            tab(:,:,:) = atan2(cg%q(qna%ind(imag_n))%span(cg%ijkse), cg%q(qna%ind(re_n))%span(cg%ijkse))
         case default
            ierrh = -1
      end select

   end subroutine mand_vars

!> \brief Mark interesting blocks for refinement

   subroutine mark_the_set

      use cg_list,    only: cg_list_element
      use cg_leaves,  only: leaves
      use domain,     only: dom
      use fluidindex, only: iarr_all_dn

      implicit none

      type(cg_list_element), pointer :: cgl
      real :: nitd, diffmax
      integer(kind=4) :: i, j, k

      cgl => leaves%first
      do while (associated(cgl))
         associate (cg => cgl%cg)
            ! Cannot use mand_n as long as it stays uninitialized during second call to problem_refine_derefine in update_refinement
            ! wna%fi is vital and is therefore automatically prolonged
            ! Possible fixes:
            ! * make mand_n vital
            ! * do not call problem_refine_derefine twice in update_refinement
!!$            nitd = exp(maxval(cg%q(qna%ind(mand_n))%span(cg%ijkse), mask=cg%leafmap)) - &
!!$                 & exp(minval(cg%q(qna%ind(mand_n))%span(cg%ijkse), mask=cg%leafmap))

            diffmax = -huge(1.)
            ! Look one cell beyond boundary to prevent unnecessary derefinements
            do i = cg%is-dom%D_x, cg%ie+dom%D_x
               do j = cg%js-dom%D_y, cg%je+dom%D_y
                  do k = cg%ks-dom%D_z, cg%ke+dom%D_z
                     nitd = maxval(abs(exp(cg%u(iarr_all_dn(1), i, j, k)) - [ &
                          &            exp(cg%u(iarr_all_dn(1), i+dom%D_x, j, k)), &
                          &            exp(cg%u(iarr_all_dn(1), i-dom%D_x, j, k)), &
                          &            exp(cg%u(iarr_all_dn(1), i, j+dom%D_y, k)), &
                          &            exp(cg%u(iarr_all_dn(1), i, j-dom%D_y, k)), &
                          &            exp(cg%u(iarr_all_dn(1), i, j, k+dom%D_z)), &
                          &            exp(cg%u(iarr_all_dn(1), i, j, k-dom%D_z)) ] ) )
                     if (nitd >= ref_thr) call cg%flag%set(i, j, k)
                     diffmax = max(diffmax, nitd)
                  enddo
               enddo
            enddo

         end associate
         cgl => cgl%nxt
      enddo

   end subroutine mark_the_set

end module initproblem
