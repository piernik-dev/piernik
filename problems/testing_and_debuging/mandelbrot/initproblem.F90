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
!! subsystem. It can be also used for presenting how precision of different
!! floating-point types affect the results.
!!
!! It is not intended for serious research. There are more efficient fractal
!! generators, but one may consider some fun ideas:
!! * Use multiprecision or implement a custom fixed-point representation
!!   optimized for these calculations.
!! * Use speedup tricks like these in fast deep zoom programs.
!! * Detect the interiors of minibrots for further speedups. This already
!!   works to some extent because they are not refined too much.
!<

module initproblem

   use constants, only: dsetnamelen, cwdlen, FP_REAL, FP_DOUBLE, FP_EXT, FP_QUAD

   implicit none

   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   ! namelist parameters
   integer(kind=4)       :: maxiter     !< Maximum number of iterations
   logical               :: smooth_map  !< Try continuous colouring
   real                  :: ref_thr     !< threshold for refining a grid
   character(len=cwdlen) :: precision   !< precision of Mandelbrot calculations, long string used only for simplicity here
   logical               :: log_polar   !< Use polar mapping around x_center + i * y_center
   character(len=cwdlen) :: x_center    !< x-coordinate of the center (crucial for polar mode)
   character(len=cwdlen) :: y_center    !< y-coordinate of the center (crucial for polar mode)
   real                  :: c_polar     !< correct colouring with x-coordinate multiplied by this factor (useful for polar mode)

   namelist /PROBLEM_CONTROL/  maxiter, smooth_map, log_polar, x_center, y_center, c_polar, ref_thr, precision

   ! other private data
   character(len=dsetnamelen), parameter :: mand_n = "mand", re_n = "real", imag_n = "imag"

   integer :: prec  ! decoded precision level
   real(kind=FP_QUAD) :: xcq, ycq  ! decoded central coordinates
   enum, bind(C)
      enumerator :: FP_QCMPLX = maxval([FP_REAL, FP_DOUBLE, FP_EXT, FP_QUAD]) + 1, FP_2FLOAT, FP_2DOUBLE !, FP_2EXT, FP_2QUAD
   end enum

   type :: double_float
      real(kind=FP_REAL) :: hi, lo
   end type double_float

   type :: double_double
      real(kind=FP_DOUBLE) :: hi, lo
   end type double_double

contains

   pure function df_from_quad(value) result(pair)

      implicit none

      real(kind=FP_QUAD), intent(in) :: value
      type(double_float) :: pair

      pair%hi = real(value, kind=FP_REAL)
      pair%lo = real(value - real(pair%hi, kind=FP_QUAD), kind=FP_REAL)

   end function df_from_quad

   pure subroutine df_two_sum(a, b, sum, error)

      implicit none

      real(kind=FP_REAL), intent(in) :: a, b
      real(kind=FP_REAL), intent(out) :: sum, error
      real(kind=FP_REAL) :: b_virtual

      sum = a + b
      b_virtual = sum - a
      error = (a - (sum - b_virtual)) + (b - b_virtual)

   end subroutine df_two_sum

   pure subroutine df_two_product(a, b, product, error)

      implicit none

      real(kind=FP_REAL), intent(in) :: a, b
      real(kind=FP_REAL), intent(out) :: product, error
      real(kind=FP_REAL), parameter :: splitter = 4097._FP_REAL
      real(kind=FP_REAL) :: a_split, a_hi, a_lo, b_split, b_hi, b_lo
      real(kind=FP_REAL) :: error1, error2, error3

      product = a * b
      a_split = splitter * a
      a_hi = a_split - (a_split - a)
      a_lo = a - a_hi
      b_split = splitter * b
      b_hi = b_split - (b_split - b)
      b_lo = b - b_hi
      error1 = product - a_hi * b_hi
      error2 = error1 - a_lo * b_hi
      error3 = error2 - a_hi * b_lo
      error = a_lo * b_lo - error3

   end subroutine df_two_product

   pure function df_add(a, b) result(sum)

      implicit none

      type(double_float), intent(in) :: a, b
      type(double_float) :: sum
      real(kind=FP_REAL) :: hi_sum, hi_error, lo_sum, lo_error
      real(kind=FP_REAL) :: correction, correction_error, result_hi, result_error

      call df_two_sum(a%hi, b%hi, hi_sum, hi_error)
      call df_two_sum(a%lo, b%lo, lo_sum, lo_error)
      call df_two_sum(hi_error, lo_sum, correction, correction_error)
      call df_two_sum(hi_sum, correction, result_hi, result_error)
      result_error = result_error + correction_error + lo_error
      call df_two_sum(result_hi, result_error, sum%hi, sum%lo)

   end function df_add

   pure function df_subtract(a, b) result(difference)

      implicit none

      type(double_float), intent(in) :: a, b
      type(double_float) :: difference, negative_b

      negative_b%hi = -b%hi
      negative_b%lo = -b%lo
      difference = df_add(a, negative_b)

   end function df_subtract

   pure function df_multiply(a, b) result(product)

      implicit none

      type(double_float), intent(in) :: a, b
      type(double_float) :: product, cross

      call df_two_product(a%hi, b%hi, product%hi, product%lo)
      call df_two_product(a%hi, b%lo, cross%hi, cross%lo)
      product = df_add(product, cross)
      call df_two_product(a%lo, b%hi, cross%hi, cross%lo)
      product = df_add(product, cross)
      call df_two_product(a%lo, b%lo, cross%hi, cross%lo)
      product = df_add(product, cross)

   end function df_multiply

   pure function dd_from_quad(value) result(pair)

      implicit none

      real(kind=FP_QUAD), intent(in) :: value

      type(double_double) :: pair

      pair%hi = real(value, kind=FP_DOUBLE)
      pair%lo = real(value - real(pair%hi, kind=FP_QUAD), kind=FP_DOUBLE)

   end function dd_from_quad

   pure subroutine two_sum(a, b, sum, error)

      implicit none

      real(kind=FP_DOUBLE), intent(in) :: a, b
      real(kind=FP_DOUBLE), intent(out) :: sum, error

      real(kind=FP_DOUBLE) :: b_virtual

      sum = a + b
      b_virtual = sum - a
      error = (a - (sum - b_virtual)) + (b - b_virtual)

   end subroutine two_sum

   pure subroutine two_product(a, b, product, error)

      implicit none

      real(kind=FP_DOUBLE), intent(in) :: a, b
      real(kind=FP_DOUBLE), intent(out) :: product, error

      real(kind=FP_DOUBLE), parameter :: splitter = 134217729._FP_DOUBLE
      real(kind=FP_DOUBLE) :: a_split, a_hi, a_lo, b_split, b_hi, b_lo
      real(kind=FP_DOUBLE) :: error1, error2, error3

      product = a * b
      a_split = splitter * a
      a_hi = a_split - (a_split - a)
      a_lo = a - a_hi
      b_split = splitter * b
      b_hi = b_split - (b_split - b)
      b_lo = b - b_hi
      error1 = product - a_hi * b_hi
      error2 = error1 - a_lo * b_hi
      error3 = error2 - a_hi * b_lo
      error = a_lo * b_lo - error3

   end subroutine two_product

   pure function dd_add(a, b) result(sum)

      implicit none

      type(double_double), intent(in) :: a, b

      type(double_double) :: sum
      real(kind=FP_DOUBLE) :: hi_sum, hi_error, lo_sum, lo_error
      real(kind=FP_DOUBLE) :: correction, correction_error, result_hi, result_error

      call two_sum(a%hi, b%hi, hi_sum, hi_error)
      call two_sum(a%lo, b%lo, lo_sum, lo_error)
      call two_sum(hi_error, lo_sum, correction, correction_error)
      call two_sum(hi_sum, correction, result_hi, result_error)
      result_error = result_error + correction_error + lo_error
      call two_sum(result_hi, result_error, sum%hi, sum%lo)

   end function dd_add

   pure function dd_subtract(a, b) result(difference)

      implicit none

      type(double_double), intent(in) :: a, b

      type(double_double) :: difference, negative_b

      negative_b%hi = -b%hi
      negative_b%lo = -b%lo
      difference = dd_add(a, negative_b)

   end function dd_subtract

   pure function dd_multiply(a, b) result(product)

      implicit none

      type(double_double), intent(in) :: a, b

      type(double_double) :: product, cross

      call two_product(a%hi, b%hi, product%hi, product%lo)
      call two_product(a%hi, b%lo, cross%hi, cross%lo)
      product = dd_add(product, cross)
      call two_product(a%lo, b%hi, cross%hi, cross%lo)
      product = dd_add(product, cross)
      call two_product(a%lo, b%lo, cross%hi, cross%lo)
      product = dd_add(product, cross)

   end function dd_multiply

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
      use constants,  only: ydim, LO, HI, dpi, INVALID, I_THREE
      use dataio_pub, only: warn, die, nh
      use domain,     only: dom
      use mpisetup,   only: lbuff, ibuff, rbuff, master, slave

      implicit none

      character(len=cwdlen), dimension(I_THREE) :: lcbuff !< buffer for long string parameters

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

         lcbuff(1) = x_center
         lcbuff(2) = y_center
         lcbuff(3) = precision

      endif

      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(lbuff)
      call piernik_MPI_Bcast(rbuff)
      call piernik_MPI_Bcast(lcbuff, cwdlen)

      if (slave) then

         maxiter    = ibuff(1)

         smooth_map = lbuff(1)
         log_polar  = lbuff(2)

         ref_thr    = rbuff(1)
         c_polar    = rbuff(2)

         x_center   = lcbuff(1)
         y_center   = lcbuff(2)
         precision  = lcbuff(3)

      endif

      if (any(dom%has_dir(:) .neqv. [ .true., .true., .false. ])) &
           call die("[initproblem:read_problem_par] Mandelbrot is supposed to be run only in the XY plane without a Z direction")

      select case (trim(precision))
         case ("single", "float")
            ! Human-friendly aliases for the internal scalar precision levels.
            prec = FP_REAL
         case ("double")
            prec = FP_DOUBLE
         case ("extended", "long")
            prec = FP_EXT
         case ("quad")
            prec = FP_QUAD
         case ("quad complex")
            ! ``quad complex'' is not a scalar kind here; it uses the native complex type.
            prec = FP_QCMPLX
         case ("double-float")
            ! ``double-float'' uses two single-precision components.
            prec = FP_2FLOAT
         case ("double-double")
            ! ``double-double'' is a demonstration of Dekker (1971) double-double arithmetic.
            prec = FP_2DOUBLE
         case default
            call die("[initproblem:read_problem_par] unsupported precision; use 'single', 'double', 'extended', 'quad', 'quad complex', 'double-float', or 'double-double'")
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

! The double-double arithmetic is kept separate from the built-in quad path.

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
         case (FP_2FLOAT)
            ! This is a demonstration of Dekker (1971) double-float arithmetic.
            block
               type(double_float) :: zx, zy, cx, cy, zt

               cx = df_from_quad(xcq + real(x, kind=FP_QUAD))
               cy = df_from_quad(ycq + real(y, kind=FP_QUAD))
               zx = cx
               zy = cy
               do while (zx%hi*zx%hi + zy%hi*zy%hi < bailout2 .and. nit < maxiter)
                  zt = df_add(df_subtract(df_multiply(zx, zx), df_multiply(zy, zy)), cx)
                  zy = df_add(df_add(df_multiply(zx, zy), df_multiply(zx, zy)), cy)
                  zx = zt
                  nit = nit + 1
               enddo
               x = real(zx%hi, kind=kind(x))
               y = real(zy%hi, kind=kind(y))

            end block
         case (FP_2DOUBLE)
            ! This is a demonstration of Dekker (1971) double-double arithmetic.
            block
               type(double_double) :: zx, zy, cx, cy, zt

               cx = dd_from_quad(xcq + real(x, kind=FP_QUAD))
               cy = dd_from_quad(ycq + real(y, kind=FP_QUAD))
               zx = cx
               zy = cy
               do while (zx%hi*zx%hi + zy%hi*zy%hi < bailout2 .and. nit < maxiter)
                  zt = dd_add(dd_subtract(dd_multiply(zx, zx), dd_multiply(zy, zy)), cx)
                  zy = dd_add(dd_add(dd_multiply(zx, zy), dd_multiply(zx, zy)), cy)
                  zx = zt
                  nit = nit + 1
               enddo
               x = real(zx%hi, kind=kind(x))
               y = real(zy%hi, kind=kind(y))

            end block
         case default
            call die("[initproblem:calculate_mandelbrot] non-implemented precision")
            x = 0.
            y = 0.
            nit = 0
      end select

      rnit = nit
      ! The stored field is a smooth, log-like escape-time measure rather than the
      ! raw iteration count. This produces a continuous color map without changing
      ! the actual escape test.
      if (smooth_map .and. x*x + y*y > bailout2) rnit = rnit + 1 - log(log(sqrt(x*x + y*y)))/log(2.)
      if (nit >= maxiter) then
         mand = min_log_mand
      else
         ! The polar bias is useful for making logarithmic radial features stand out,
         ! while the lower bound keeps the colouring finite for interior points.
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
            ! Check a one-cell halo before deciding to derefine. This avoids flagging
            ! blocks only because the escape metric changes sharply at a coarse/fine
            ! interface, which would otherwise cause spurious AMR oscillations.
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
