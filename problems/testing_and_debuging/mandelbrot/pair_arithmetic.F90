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

!> \brief Compensated arithmetic for values represented by two floating-point components.
!!
!! A pair represents the unevaluated sum hi + lo. The leading component uses
!! either single or double precision; the low component retains rounding
!! information that would otherwise be lost. Error-free transforms recover
!! residuals from basic additions and products and carry them in the low
!! component; the resulting pair still has finite precision.
!!
!! The transforms assume IEEE round-to-nearest arithmetic. Their guarantees
!! require finite inputs, no overflowing intermediate operations, and
!! representable residuals. They do not protect against overflow or underflow.
!!
!! New pair formats require a pair type, conversion routine, and arithmetic
!! overloads. A double-extended or double-quad pair may carry more precision
!! than FP_QUAD, so coordinates for those formats must be parsed directly into
!! the pair or through a wider intermediate rather than first rounded to quad.
!<

module pair_arithmetic

   use constants, only: FP_REAL, FP_DOUBLE

   implicit none

   private
   public :: double_float, double_double, float_from_quad, double_from_quad
   public :: pair_add, pair_subtract, pair_multiply

   !> A two-component value based on the FP_REAL kind.
   !!
   !! The represented value is approximately hi + lo.
   type :: double_float
      real(kind=FP_REAL) :: hi !< Leading component
      real(kind=FP_REAL) :: lo !< Low component carrying the rounding residual
   end type double_float

   !> A two-component value based on the FP_DOUBLE kind.
   !!
   !! The represented value is approximately hi + lo.
   type :: double_double
      real(kind=FP_DOUBLE) :: hi !< Leading component
      real(kind=FP_DOUBLE) :: lo !< Low component carrying the rounding residual
   end type double_double

   !> \brief Add two values represented by the same pair type.
   !! \details Overloaded for double_float and double_double. The result is
   !! renormalized so its leading component contains the rounded sum and its
   !! low component retains the remaining correction.
   interface pair_add
      module procedure df_add, dd_add
   end interface pair_add

   !> \brief Subtract the second pair from the first.
   !! \details Overloaded for double_float and double_double.
   interface pair_subtract
      module procedure df_subtract, dd_subtract
   end interface pair_subtract

   !> \brief Multiply two values represented by the same pair type.
   !! \details Overloaded for double_float and double_double. The four
   !! component products are accumulated with their product residuals.
   interface pair_multiply
      module procedure df_multiply, dd_multiply
   end interface pair_multiply

contains

   !> \brief Split a quad-precision value into two FP_REAL components.
   !! \param[in] value Value to convert.
   !! \return A double_float pair whose components retain the leading value
   !! and representable residual.
   pure function float_from_quad(value) result(pair)

      use constants, only: FP_QUAD

      implicit none

      real(kind=FP_QUAD), intent(in) :: value

      type(double_float) :: pair

      pair%hi = real(value, kind=FP_REAL)
      pair%lo = real(value - real(pair%hi, kind=FP_QUAD), kind=FP_REAL)

   end function float_from_quad

   !> \brief Split a quad-precision value into two FP_DOUBLE components.
   !! \param[in] value Value to convert.
   !! \return A double_double pair whose components retain the leading value
   !! and representable residual.
   pure function double_from_quad(value) result(pair)

      use constants, only: FP_QUAD

      implicit none

      real(kind=FP_QUAD), intent(in) :: value

      type(double_double) :: pair

      pair%hi = real(value, kind=FP_DOUBLE)
      pair%lo = real(value - real(pair%hi, kind=FP_QUAD), kind=FP_DOUBLE)

   end function double_from_quad

   !> \brief Compute a rounded sum and its addition residual for FP_REAL values.
   !! \details Knuth's TwoSum transform. With round-to-nearest arithmetic, no
   !! overflowing intermediates, and a representable residual, sum + error
   !! equals the exact sum.
   !! \param[in] a First addend.
   !! \param[in] b Second addend.
   !! \param[out] sum Rounded sum.
   !! \param[out] error Addition residual.
   pure subroutine df_two_sum(a, b, sum, error)

      implicit none

      real(kind=FP_REAL), intent(in) :: a, b
      real(kind=FP_REAL), intent(out) :: sum, error

      real(kind=FP_REAL) :: b_virtual

      sum = a + b
      b_virtual = sum - a
      error = (a - (sum - b_virtual)) + (b - b_virtual)

   end subroutine df_two_sum

   !> \brief Compute a rounded sum and its addition residual for FP_DOUBLE values.
   !! \details This is the double-precision form of Knuth's TwoSum transform.
   !! Under the same rounding and range conditions as df_two_sum, sum + error
   !! equals the exact sum.
   !! \param[in] a First addend.
   !! \param[in] b Second addend.
   !! \param[out] sum Rounded sum.
   !! \param[out] error Addition residual.
   pure subroutine dd_two_sum(a, b, sum, error)

      implicit none

      real(kind=FP_DOUBLE), intent(in) :: a, b
      real(kind=FP_DOUBLE), intent(out) :: sum, error

      real(kind=FP_DOUBLE) :: b_virtual

      sum = a + b
      b_virtual = sum - a
      error = (a - (sum - b_virtual)) + (b - b_virtual)

   end subroutine dd_two_sum

   !> \brief Recover the rounding residual of an FP_REAL product.
   !! \details Dekker splitting separates each operand into leading and low
   !! parts before reconstructing the product error.
   !! \param[in] a First factor.
   !! \param[in] b Second factor.
   !! \param[out] product Rounded product.
   !! \param[out] error Product residual.
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

   !> \brief Recover the rounding residual of an FP_DOUBLE product.
   !! \details This is the double-precision form of Dekker's splitting
   !! algorithm; the splitter is chosen for the FP_DOUBLE significand width.
   !! \param[in] a First factor.
   !! \param[in] b Second factor.
   !! \param[out] product Rounded product.
   !! \param[out] error Product residual.
   pure subroutine dd_two_product(a, b, product, error)

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

   end subroutine dd_two_product

   !> \brief Add two double_float values and renormalize the result.
   !! \param[in] a First addend.
   !! \param[in] b Second addend.
   !! \return Sum represented by leading and low components.
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

   !> \brief Add two double_double values and renormalize the result.
   !! \param[in] a First addend.
   !! \param[in] b Second addend.
   !! \return Sum represented by leading and low components.
   pure function dd_add(a, b) result(sum)

      implicit none

      type(double_double), intent(in) :: a, b
      type(double_double) :: sum
      real(kind=FP_DOUBLE) :: hi_sum, hi_error, lo_sum, lo_error
      real(kind=FP_DOUBLE) :: correction, correction_error, result_hi, result_error

      call dd_two_sum(a%hi, b%hi, hi_sum, hi_error)
      call dd_two_sum(a%lo, b%lo, lo_sum, lo_error)
      call dd_two_sum(hi_error, lo_sum, correction, correction_error)
      call dd_two_sum(hi_sum, correction, result_hi, result_error)
      result_error = result_error + correction_error + lo_error
      call dd_two_sum(result_hi, result_error, sum%hi, sum%lo)

   end function dd_add

   !> \brief Subtract one double_float value from another.
   !! \param[in] a Minuend.
   !! \param[in] b Subtrahend.
   !! \return Difference represented by leading and low components.
   pure function df_subtract(a, b) result(difference)

      implicit none

      type(double_float), intent(in) :: a, b

      type(double_float) :: difference, negative_b

      negative_b%hi = -b%hi
      negative_b%lo = -b%lo
      difference = df_add(a, negative_b)

   end function df_subtract

   !> \brief Subtract one double_double value from another.
   !! \param[in] a Minuend.
   !! \param[in] b Subtrahend.
   !! \return Difference represented by leading and low components.
   pure function dd_subtract(a, b) result(difference)

      implicit none

      type(double_double), intent(in) :: a, b

      type(double_double) :: difference, negative_b

      negative_b%hi = -b%hi
      negative_b%lo = -b%lo
      difference = dd_add(a, negative_b)

   end function dd_subtract

   !> \brief Multiply two double_float values.
   !! \param[in] a First factor.
   !! \param[in] b Second factor.
   !! \return Product represented by leading and low components.
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

   !> \brief Multiply two double_double values.
   !! \param[in] a First factor.
   !! \param[in] b Second factor.
   !! \return Product represented by leading and low components.
   pure function dd_multiply(a, b) result(product)

      implicit none

      type(double_double), intent(in) :: a, b

      type(double_double) :: product, cross

      call dd_two_product(a%hi, b%hi, product%hi, product%lo)
      call dd_two_product(a%hi, b%lo, cross%hi, cross%lo)
      product = dd_add(product, cross)
      call dd_two_product(a%lo, b%hi, cross%hi, cross%lo)
      product = dd_add(product, cross)
      call dd_two_product(a%lo, b%lo, cross%hi, cross%lo)
      product = dd_add(product, cross)

   end function dd_multiply

end module pair_arithmetic
