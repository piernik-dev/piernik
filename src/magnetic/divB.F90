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
!! \brief Module responsible for calculation of magnetic field divergence
!<

module div_B
! pulled by ANY
! originally it was 'by MAGNETIC' but I think I the future we may choose to drop most #define, including MAGNETIC

   use constants, only: dsetnamelen, INVALID

   implicit none

   character(len=dsetnamelen), parameter :: divB_n = "div_B"
   integer(kind=4), protected            :: idivB = INVALID

!! d/dx factors for cell-centered field:                                       ; positive half of the stencil from Maxima:
!! 2nd order: [                      -1/2, 0, 1/2                       ] / dx ; linsolve_by_lu(matrix([2]),                                                                                       matrix([1]));
!! 4th order: [                1/12, -2/3, 0, 2/3, -1/12                ] / dx ; linsolve_by_lu(matrix([2,2*2],        [2,2*2**3]),                                                                matrix([1],[0]));
!! 6th order: [        -1/60,  3/20, -3/4, 0, 3/4, -3/20, 1/60          ] / dx ; linsolve_by_lu(matrix([2,2*2,2*3],    [2,2*2**3,2*3**3],       [2,2*2**5,2*3**5]),                                matrix([1],[0],[0]));
!! 8th order: [ 1/280, -4/105,  1/5, -4/5, 0, 4/5, -1/5,  4/105, -1/280 ] / dx ; linsolve_by_lu(matrix([2,2*2,2*3,2*4],[2,2*2**3,2*3**3,2*4**3],[2,2*2**5,2*3**5,2*4**5],[2,2*2**7,2*3**7,2*4**7]),matrix([1],[0],[0],[0]));
!!
!! higher order formula generator:
!! echo 10 | awk '{o=$1/2; printf("linsolve_by_lu(matrix("); for (i=1;i<=o;i++) { printf("[2"); for (j=2;j<=o;j++) { printf(",2*%d**%d",j,(2*i-1))} printf("]"); if (i<o) printf(","); } printf("),matrix([1]"); for (i=1;i<o;i++) printf(",[0]"); printf("));\n")}' | maxima
!!
!! d/dx factors for face-centered field:                                                              ; positive half of the stencil from Maxima:
!! 2nd order: [                                -1,         1                                   ] / dx ; linsolve_by_lu(matrix([1]),                                                               matrix([1]));
!! 4th order: [                     1/24       -9/8        9/8       -1/24                     ] / dx ; linsolve_by_lu(matrix([1,3],    [1,3**3]),                                                matrix([1],[0]));
!! 6th order: [          -3/640,   25/384,    -75/64,     75/64,    -25/384,   3/640           ] / dx ; linsolve_by_lu(matrix([1,3,5],  [1,3**3,5**3],     [1,3**5,5**5]),                        matrix([1],[0],[0]));
!! 8th order: [ 5/7168, -49/5120, 245/3072, -1125/1024, 1225/1024, -245/3072, 49/5120, -5/7168 ] / dx ; linsolve_by_lu(matrix([1,3,5,7],[1,3**3,5**3,7**3],[1,3**5,5**5,7**5],[1,3**7,5**7,7**7]),matrix([1],[0],[0],[0]));
!!
!! higher order formula generator:
!! echo 10 | awk '{o=$1/2; printf("linsolve_by_lu(matrix("); for (i=1;i<=o;i++) { printf("[1"); for (j=2;j<=o;j++) { printf(",%d**%d",(2*j-1),(2*i-1))} printf("]"); if (i<o) printf(","); } printf("),matrix([1]"); for (i=1;i<o;i++) printf(",[0]"); printf("));\n")}' | maxima
!!
!! For Ubuntu one may need maxima-share package

   integer(kind=4), parameter :: max_c = 4
   real, dimension(max_c, max_c), parameter :: coeff_c = reshape( [ [ 1./2.,  0.,     0.,       0.      ], &
        &                                                           [ 2./3., -1./12., 0.,       0.      ], &
        &                                                           [ 3./4., -3./20., 1./60.,   0.      ], &
        &                                                           [ 4./5., -1./5.,  4./105., -1./280. ] ], [max_c, max_c] ) ! it is sort of stupid that we can't directly initialize this as a 2D array here
   real, dimension(max_c, max_c), parameter :: coeff_f = reshape( [ [    1.,          0.,        0.,        0.       ], &
        &                                                           [    9./8.,      -1./24.,    0.,        0.       ], &
        &                                                           [   75./64.,    -25./384.,   3./640.,   0.       ], &
        &                                                           [ 1225./1024., -245./3072., 49./5120., -5./7168. ] ], [max_c, max_c] )

   private
   public :: divB, idivB, divB_c_IO, print_divB_norm

contains

!>
!! \brief Print estimates of div(B) norm on demand
!<

   subroutine print_divB_norm

      use cg_leaves,        only: leaves
      use constants,        only: I_TWO, psi_n, V_INFO
      use dataio_pub,       only: printinfo, msg
      use domain,           only: dom
      use mpisetup,         only: master
      use named_array_list, only: qna

      implicit none

      integer(kind=4) :: i

#ifndef MAGNETIC
      return
#endif /* !MAGNETIC */

      write(msg,'(a)')"[divB:print_divB_norm]"
      do i = I_TWO, I_TWO*min(max_c, dom%nb), I_TWO  ! only even-order norms;  tricky: I_TWO*max_c
         call divB(i)
         write(msg(len_trim(msg)+1:), '(a,i1,a,g12.5)') " |divB|_", i, "= ", leaves%norm_sq(idivB) / sqrt(dom%Vol)
      enddo
      if (qna%exists(psi_n)) write(msg(len_trim(msg)+1:), '(a,g12.5)') " |psi|= ", leaves%norm_sq(qna%ind(psi_n)) / sqrt(dom%Vol)
      if (master) call printinfo(msg, V_INFO)

   end subroutine print_divB_norm

!>
!! \brief Allocate extra space for divB
!!
!! ToDo: allow for use of cg%wa in memory-critical applications (beware of delayed read of divB).
!<

   subroutine divB_init

      use cg_list_global,   only: all_cg
      use constants,        only: INVALID
      use dataio_pub,       only: die
      use named_array_list, only: qna

      implicit none

#ifndef MAGNETIC
      return
#endif /* !MAGNETIC */

      if (qna%exists(divB_n)) then
         if (idivB /= qna%ind(divB_n)) call die ("[divB:divB_init] qna%exists(divB_n) .and. idivB /= qna%ind(divB_n)")
      else
         if (idivB /= INVALID) call die ("[divB:divB_init] .not. qna%exists(divB_n) .and. idivB /= INVALID")
         call all_cg%reg_var(divB_n)
         idivB = qna%ind(divB_n)
      endif

   end subroutine divB_init

!>
!! \brief Calculate divB of requested order
!!
!! BEWARE: do not use this routine in datafields_hdf5 calls because it may result in terribly slow IO in heavy-AMR simulations. Use divB_c_IO instead.
!! BEWARE: magic integers
!<

   subroutine divB(order, cell_centered)

      use constants,  only: GEO_XYZ
      use dataio_pub, only: msg, die, warn
      use domain,     only: dom
      use global,     only: cc_mag

      implicit none

      integer(kind=4), optional, intent(in) :: order
      logical,         optional, intent(in) :: cell_centered

      integer(kind=4) :: ord
      logical         :: ccB
      integer, parameter :: max_ord = 8

#ifndef MAGNETIC
      call warn("[[divB:divB] No magnetic field!")
      return
#endif /* !MAGNETIC */

      if (dom%geometry_type /= GEO_XYZ) call die("[divB:divB] non-cartesian geometry not implemented yet.")

      ord = 2
      if (present(order)) then
         if (order > max_ord) call die("[divB:divB] Highest allowed order of div(B) is currently equal to 8.")
         select case (order)
            case (7, 8)
               ord = 8
            case (5, 6)
               ord = 6
            case (3, 4)
               ord = 4
            case (:2)
               ord = 2
            case default
               call die("[divB:divB] invalid order of div(B)")
         end select
         if (order /= ord) then
            write(msg, '(a, i3)')"[divB:divB] order of div(B) was upgraded to ", ord
            call warn(msg)
         endif
         if (dom%nb < ord/2) call die("[divB:divB] Not enough guardcells for requested approximation order")  ! assumed stencils on uniform grid
      endif

      ccB = cc_mag
      if (present(cell_centered)) ccB = cell_centered

      call divB_init  ! it should be BOTH safe and cheap to call it multiple times

      call divB_c(ord, ccB)
      if (ord > max_ord) call warn("[divB:divB] only 8th order of div(B) is currently implemented for magnetic field.")  !! BEWARE: order hardcoded in the string

   end subroutine divB

!>
!! \brief Calculate 2nd order of divB for cell- and face-centered field and put it in the proper array.
!! If required, calculate higher order corrections and add it to the proper array.
!!
!! OPT: this may not be the best approach performance-wise.
!<

   subroutine divB_c(ord, ccB)

      use cg_cost_data, only: I_OTHER
      use cg_leaves,    only: leaves
      use cg_list,      only: cg_list_element
      use constants,    only: I_ONE, I_TWO
      use dataio_pub,   only: die
      use func,         only: operator(.notequals.)

      implicit none

      integer(kind=4), intent(in) :: ord
      logical,         intent(in) :: ccB

      type(cg_list_element), pointer :: cgl
      integer(kind=4) :: o, i
      real, dimension(max_c) :: coeff

      o = ord/I_TWO  ! BEWARE: tricky, assumed stencils on uniform grid
      if (o < I_ONE .or. o > max_c .or. I_TWO*o /= ord) call die("[divB:divB_c] cannot find coefficient") ! no odd order allowed here just in case

      if (ccB) then
         coeff = coeff_c(:, o)
      else
         coeff = coeff_f(:, o)
      endif

      cgl => leaves%first
      do while (associated(cgl))
         call cgl%cg%costs%start

         cgl%cg%q(idivB)%arr(                                cgl%cg%is:cgl%cg%ie, cgl%cg%js:cgl%cg%je, cgl%cg%ks:cgl%cg%ke) = sixpoint(cgl%cg, coeff(I_ONE), I_ONE, ccB)
         do i = I_TWO, max_c
            if (coeff(i) .notequals. 0.) cgl%cg%q(idivB)%arr(cgl%cg%is:cgl%cg%ie, cgl%cg%js:cgl%cg%je, cgl%cg%ks:cgl%cg%ke) = &
                 &                       cgl%cg%q(idivB)%arr(cgl%cg%is:cgl%cg%ie, cgl%cg%js:cgl%cg%je, cgl%cg%ks:cgl%cg%ke) + sixpoint(cgl%cg, coeff(i),     i,     ccB)
         enddo

         call cgl%cg%costs%stop(I_OTHER)
         cgl => cgl%nxt
      enddo

   end subroutine divB_c

!>
!! \brief single-cg variant of the routine divB_c intended for I/O  usage.
!!
!! The output is already converted to single-precision.
!<

   function divB_c_IO(cg, ord, cell_centered)

      use constants,  only: I_ONE, I_TWO
      use dataio_pub, only: die
      use func,       only: operator(.notequals.)
      use grid_cont,  only: grid_container

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer(kind=4),               intent(in) :: ord
      logical,                       intent(in) :: cell_centered

      real, dimension(RNG)                      :: divB_c_IO
      integer(kind=4)                           :: o, i
      real, dimension(max_c)                    :: coeff

      o = ord/I_TWO  ! BEWARE: tricky, assumed stencils on uniform grid
      if (o < I_ONE .or. o > max_c .or. I_TWO*o /= ord) call die("[divB:divB_c] cannot find coefficient") ! no odd order allowed here just in case

      if (cell_centered) then
         coeff = coeff_c(:, o)
      else
         coeff = coeff_f(:, o)
      endif

      divB_c_IO = sixpoint(cg, coeff(I_ONE), I_ONE, cell_centered)
      do i = I_TWO, max_c
         if (coeff(i) .notequals. 0.) divB_c_IO = divB_c_IO + sixpoint(cg, coeff(i), i, cell_centered)
      enddo

   end function divB_c_IO

!>
!! \brief Return estimate of derivative taken at given span (in cells) and multiplied by given coefficient for a given grid container.
!! This is a basic piece useful to construct various order derivatives and takes the advantage of stencil symmetry on an uniform grid.
!!
!! OPT: this may not be the best approach performance-wise but it is compact and we don't expect to evaluate it foo often.
!<

   function sixpoint(cg, coeff, span, cell_centered)

      use constants,  only: xdim, ydim, zdim
      use dataio_pub, only: die
      use domain,     only: dom
      use func,       only: operator(.equals.)
      use grid_cont,  only: grid_container

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      real,                          intent(in) :: coeff
      integer(kind=4),               intent(in) :: span
      logical,                       intent(in) :: cell_centered

      real, dimension(RNG)                      :: sixpoint
      integer                                   :: spm1

      if ((coeff .equals. 0.) .or. span <=0 .or. span > dom%nb) call die("[divB:sixpoint] coeff == 0. or unacceptable span")

      if (.not. associated(cg%b)) call die("[divB:sixpoint] no cg%b")

      if (cell_centered) then
         sixpoint = coeff * ( &
                 & (cg%b(xdim, cg%is+span*dom%D_x:cg%ie+span*dom%D_x, cg%js             :cg%je,              cg%ks             :cg%ke             ) - &
                 &  cg%b(xdim, cg%is-span*dom%D_x:cg%ie-span*dom%D_x, cg%js             :cg%je,              cg%ks             :cg%ke             )   )/cg%dx + &
                 & (cg%b(ydim, cg%is             :cg%ie,              cg%js+span*dom%D_y:cg%je+span*dom%D_y, cg%ks             :cg%ke             ) - &
                 &  cg%b(ydim, cg%is             :cg%ie,              cg%js-span*dom%D_y:cg%je-span*dom%D_y, cg%ks             :cg%ke             )   )/cg%dy + &
                 & (cg%b(zdim, cg%is             :cg%ie,              cg%js             :cg%je,              cg%ks+span*dom%D_z:cg%ke+span*dom%D_z) - &
                 &  cg%b(zdim, cg%is             :cg%ie,              cg%js             :cg%je,              cg%ks-span*dom%D_z:cg%ke-span*dom%D_z)   )/cg%dz )
      else
         spm1 = span - 1
         sixpoint = coeff * ( &
                 & (cg%b(xdim, cg%is+span*dom%D_x:cg%ie+span*dom%D_x, cg%js             :cg%je,              cg%ks             :cg%ke             ) - &
                 &  cg%b(xdim, cg%is-spm1*dom%D_x:cg%ie-spm1*dom%D_x, cg%js             :cg%je,              cg%ks             :cg%ke             )   )/cg%dx + &
                 & (cg%b(ydim, cg%is             :cg%ie,              cg%js+span*dom%D_y:cg%je+span*dom%D_y, cg%ks             :cg%ke             ) - &
                 &  cg%b(ydim, cg%is             :cg%ie,              cg%js-spm1*dom%D_y:cg%je-spm1*dom%D_y, cg%ks             :cg%ke             )   )/cg%dy + &
                 & (cg%b(zdim, cg%is             :cg%ie,              cg%js             :cg%je,              cg%ks+span*dom%D_z:cg%ke+span*dom%D_z) - &
                 &  cg%b(zdim, cg%is             :cg%ie,              cg%js             :cg%je,              cg%ks-spm1*dom%D_z:cg%ke-spm1*dom%D_z)   )/cg%dz )
      endif

   end function sixpoint

end module div_B
