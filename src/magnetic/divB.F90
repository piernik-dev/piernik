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
! pulled by MAGNETIC

   use constants, only: dsetnamelen, INVALID

   implicit none

   character(len=dsetnamelen), parameter :: divB_n = "div_B"
   integer, protected                    :: idivB = INVALID

   private
   public :: divB, idivB

contains

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

      if (qna%exists(divB_n)) then
         if (idivB /= qna%ind(divB_n)) call die ("[divB::divB_init] qna%exists(divB_n) .and. idivB /= qna%ind(divB_n)")
      else
         if (idivB /= INVALID) call die ("[divB::divB_init] .not. qna%exists(divB_n) .and. idivB /= INVALID")
         call all_cg%reg_var(divB_n)
         idivB = qna%ind(divB_n)
      endif

   end subroutine divB_init

!>
!! \brief Calculate divB of requested order
!!
!! BEWARE: do not use this routie in datafields_hdf5 calls because it may result in terribly slow IO in heavy-AMR simulations.
!! In order to remove code duplication between this file and data_hdf5.F90, some additional assumptions are required, such as:
!! * use only one div(B) field and don't compare different orders/placements), or
!! * use separate qna fields for each divB field, or
!! * reorganize the routines so a single call per cg% are possible.
!! The third option seems to be most favorable.
!!
!! BEWARE: magic integers
!!
!! d/dx factors for cell-centered field:                                       ; positive half of the stencil from Maxima:
!! 2nd order: [                      -1/2, 0, 1/2                       ] / dx ; linsolve_by_lu(matrix([2]),                                                                                       matrix([1]));
!! 4th order: [                1/12, -2/3, 0, 2/3, -1/12                ] / dx ; linsolve_by_lu(matrix([2,2*2],        [2,2*2**3]),                                                                matrix([1],[0]));
!! 6th order: [        -1/60,  3/20, -3/4, 0, 3/4, -3/20, 1/60          ] / dx ; linsolve_by_lu(matrix([2,2*2,2*3],    [2,2*2**3,2*3**3],       [2,2*2**5,2*3**5]),                                matrix([1],[0],[0]));
!! 8th order: [ 1/280, -4/105,  1/5, -4/5, 0, 4/5, -1/5,  4/105, -1/280 ] / dx ; linsolve_by_lu(matrix([2,2*2,2*3,2*4],[2,2*2**3,2*3**3,2*4**3],[2,2*2**5,2*3**5,2*4**5],[2,2*2**7,2*3**7,2*4**7]),matrix([1],[0],[0],[0]));
!!
!! higher order formula generator:
!! echo 10 | awk '{o=$1/2; printf("linsolve_by_lu(matrix("); for (i=1;i<=o;i++) { printf("[2"); for (j=2;j<=o;j++) { printf(",2*%d**%d",j,(2*i-1))} printf("]"); if (i<o) printf(","); } printf("),matrix([1]"); for (i=1;i<o;i++) printf(",[0]"); printf("));\n")}'
!!
!! d/dx factors for face-centered field:                                                              ; positive half of the stencil from Maxima:
!! 2nd order: [                                -1,         1                                   ] / dx ; linsolve_by_lu(matrix([1]),                                                               matrix([1]));
!! 4th order: [                     1/24       -9/8        9/8       -1/24                     ] / dx ; linsolve_by_lu(matrix([1,3],    [1,3**3]),                                                matrix([1],[0]));
!! 6th order: [          -3/640,   25/384,    -75/64,     75/64,    -25/384,   3/640           ] / dx ; linsolve_by_lu(matrix([1,3,5],  [1,3**3,5**3],     [1,3**5,5**5]),                        matrix([1],[0],[0]));
!! 8th order: [ 5/7168, -49/5120, 245/3072, -1125/1024, 1225/1024, -245/3072, 49/5120, -5/7168 ] / dx ; linsolve_by_lu(matrix([1,3,5,7],[1,3**3,5**3,7**3],[1,3**5,5**5,7**5],[1,3**7,5**7,7**7]),matrix([1],[0],[0],[0]));
!!
!! higher order formula generator:
!! echo 10 | awk '{o=$1/2; printf("linsolve_by_lu(matrix("); for (i=1;i<=o;i++) { printf("[1"); for (j=2;j<=o;j++) { printf(",%d**%d",(2*j-1),(2*i-1))} printf("]"); if (i<o) printf(","); } printf("),matrix([1]"); for (i=1;i<o;i++) printf(",[0]"); printf("));\n")}'
!!
!! Note that for some strange reasons Maxima cannot properly handle autogenerated formulas for large orders, such as 50 :-)
!<

   subroutine divB(order, cell_centered)

      use constants,  only: GEO_XYZ
      use dataio_pub, only: msg, die, warn
      use domain,     only: dom

      implicit none

      integer, optional, intent(in) :: order
      logical, optional, intent(in) :: cell_centered

      integer :: ord
      logical :: ccB

      if (dom%geometry_type /= GEO_XYZ) call die("[divB::divB] non-cartesian geometry not implemented yet.")

      ord = 2
      if (present(order)) then
         if (order > 6) call die("[divB::divB] Highest allowed order of div(B) is currently equal to 6.")
         select case (order)
            case (5, 6)
               ord = 6
            case (3, 4)
               ord = 4
            case (:2)
               ord = 2
            case default
               call die("[divB::divB] invalid order of div(B)")
         end select
         if (order /= ord) then
            write(msg, '(a, i3)')"[divB::divB] order of div(B) was upgraded to ", ord
            call warn(msg)
         endif
         if (dom%nb < ord/2) call die("[divB::divB] Not enough guardcells for requested approximation order")  ! assumed stencils on uniform grid
      endif

      ccB = .false.  ! ToDo detect automatically
      if (present(cell_centered)) ccB = cell_centered

      call divB_init  ! it should be BOTH safe and cheap to call it multiple times
      
      if (ccB) then
         call divB_c2(ord)
         if (ord > 2) call divB_c4(ord)
         if (ord > 4) call divB_c6(ord)
         if (ord > 6) call warn("[divB::divB] only 6th order of div(B) is currently implemented for cell-centered field.")
      else  ! face-centered B
         call divB_f2(ord)
         if (ord > 2) call divB_f4(ord)
         if (ord > 4) call divB_f6(ord)
         if (ord > 6) call warn("[divB::divB] only 2nd order of div(B) is currently implemented for face-centered field.")
      endif

   end subroutine divB

!>
!! \brief Calculate 2nd order of divB for cell-centered field and put it in the proper array.
!<

   subroutine divB_c2(ord)

      use cg_leaves,  only: leaves
      use cg_list,    only: cg_list_element
      use constants,  only: xdim, ydim, zdim
      use dataio_pub, only: die
      use domain,     only: dom
      use func,       only: operator(.equals.)

      implicit none

      integer, intent(in) :: ord

      type(cg_list_element), pointer :: cgl
      real, dimension(4), parameter :: coeff = [ 1./2., 2./3., 3./4., 4./5. ]
      integer :: o

      o = ord/2  ! BEWARE: tricky, assumed stencils on uniform grid
      if (o < lbound(coeff, dim=1) .or. o > ubound(coeff, dim=1) .or. 2*o /= ord .or. (coeff(o) .equals. 0.)) call die("[divB::divB_c2] cannot find coefficient") ! no odd order allowed here just in case

      cgl => leaves%first
      do while (associated(cgl))
         associate (is => cgl%cg%is, ie => cgl%cg%ie, &
              &     js => cgl%cg%js, je => cgl%cg%je, &
              &     ks => cgl%cg%ks, ke => cgl%cg%ke)
            cgl%cg%q(idivB)%arr(   is        :ie,         js        :je,         ks        :ke        ) = coeff(o) * ( &
                 & (cgl%cg%b(xdim, is+dom%D_x:ie+dom%D_x, js        :je,         ks        :ke        ) - &
                 &  cgl%cg%b(xdim, is-dom%D_x:ie-dom%D_x, js        :je,         ks        :ke        )   )/cgl%cg%dx + &
                 & (cgl%cg%b(ydim, is        :ie,         js+dom%D_y:je+dom%D_y, ks        :ke        ) - &
                 &  cgl%cg%b(ydim, is        :ie,         js-dom%D_y:je-dom%D_y, ks        :ke        )   )/cgl%cg%dy + &
                 & (cgl%cg%b(zdim, is        :ie,         js        :je,         ks+dom%D_z:ke+dom%D_z) - &
                 &  cgl%cg%b(zdim, is        :ie,         js        :je,         ks-dom%D_z:ke-dom%D_z)   )/cgl%cg%dz )
         end associate
         cgl => cgl%nxt
      enddo

   end subroutine divB_c2

!>
!! \brief Calculate 4th order correction for 2nd order estimate of divB for cell-centered field and add it to the proper array.
!!
!! OPT: this may not be the best approach performance-wise.
!!
!! Spaghetti warning: this routime may be generalized and unified with divB_c2
!<

   subroutine divB_c4(ord)

      use cg_leaves,  only: leaves
      use cg_list,    only: cg_list_element
      use constants,  only: xdim, ydim, zdim
      use dataio_pub, only: die
      use domain,     only: dom
      use func,       only: operator(.equals.)

      implicit none

      integer, intent(in) :: ord

      type(cg_list_element), pointer :: cgl
      real, dimension(4), parameter :: coeff = [ 0., -1./12., -3./20., -1./5. ]
      integer :: o

      o = ord/2  ! BEWARE: tricky, assumed stencils on uniform grid
      if (o < lbound(coeff, dim=1) .or. o > ubound(coeff, dim=1) .or. 2*o /= ord .or. (coeff(o) .equals. 0.)) call die("[divB::divB_c4] cannot find coefficient") ! no odd order allowed here just in case

      cgl => leaves%first
      do while (associated(cgl))
         associate (is => cgl%cg%is, ie => cgl%cg%ie, &
              &     js => cgl%cg%js, je => cgl%cg%je, &
              &     ks => cgl%cg%ks, ke => cgl%cg%ke)
            cgl%cg%q(idivB)%arr(     is        :ie,             js        :je,             ks        :ke            ) = &
                 cgl%cg%q(idivB)%arr(is        :ie,             js        :je,             ks        :ke            ) + coeff(o) * ( &
                 & (cgl%cg%b(xdim,   is+2*dom%D_x:ie+2*dom%D_x, js        :je,             ks        :ke            ) - &
                 &  cgl%cg%b(xdim,   is-2*dom%D_x:ie-2*dom%D_x, js        :je,             ks        :ke            )   )/cgl%cg%dx + &
                 & (cgl%cg%b(ydim,   is        :ie,             js+2*dom%D_y:je+2*dom%D_y, ks        :ke            ) - &
                 &  cgl%cg%b(ydim,   is        :ie,             js-2*dom%D_y:je-2*dom%D_y, ks        :ke            )   )/cgl%cg%dy + &
                 & (cgl%cg%b(zdim,   is        :ie,             js        :je,             ks+2*dom%D_z:ke+2*dom%D_z) - &
                 &  cgl%cg%b(zdim,   is        :ie,             js        :je,             ks-2*dom%D_z:ke-2*dom%D_z)   )/cgl%cg%dz )
         end associate
         cgl => cgl%nxt
      enddo

   end subroutine divB_c4

!>
!! \brief Calculate 6th order correction for 4th order estimate of divB for cell-centered field and add it to the proper array.
!!
!! OPT: this may not be the best approach performance-wise.
!!
!! Spaghetti warning: this routime may be generalized and unified with divB_c4
!<

   subroutine divB_c6(ord)

      use cg_leaves,  only: leaves
      use cg_list,    only: cg_list_element
      use constants,  only: xdim, ydim, zdim
      use dataio_pub, only: die
      use domain,     only: dom
      use func,       only: operator(.equals.)

      implicit none

      integer, intent(in) :: ord

      type(cg_list_element), pointer :: cgl
      real, dimension(4), parameter :: coeff = [ 0., 0., 1./60., 4./105. ]
      integer :: o

      o = ord/2  ! BEWARE: tricky, assumed stencils on uniform grid
      if (o < lbound(coeff, dim=1) .or. o > ubound(coeff, dim=1) .or. 2*o /= ord .or. (coeff(o) .equals. 0.)) call die("[divB::divB_c6] cannot find coefficient") ! no odd order allowed here just in case

      cgl => leaves%first
      do while (associated(cgl))
         associate (is => cgl%cg%is, ie => cgl%cg%ie, &
              &     js => cgl%cg%js, je => cgl%cg%je, &
              &     ks => cgl%cg%ks, ke => cgl%cg%ke)
            cgl%cg%q(idivB)%arr(     is        :ie,             js        :je,             ks        :ke            ) = &
                 cgl%cg%q(idivB)%arr(is        :ie,             js        :je,             ks        :ke            ) + coeff(o) * ( &
                 & (cgl%cg%b(xdim,   is+3*dom%D_x:ie+3*dom%D_x, js        :je,             ks        :ke            ) - &
                 &  cgl%cg%b(xdim,   is-3*dom%D_x:ie-3*dom%D_x, js        :je,             ks        :ke            )   )/cgl%cg%dx + &
                 & (cgl%cg%b(ydim,   is        :ie,             js+3*dom%D_y:je+3*dom%D_y, ks        :ke            ) - &
                 &  cgl%cg%b(ydim,   is        :ie,             js-3*dom%D_y:je-3*dom%D_y, ks        :ke            )   )/cgl%cg%dy + &
                 & (cgl%cg%b(zdim,   is        :ie,             js        :je,             ks+3*dom%D_z:ke+3*dom%D_z) - &
                 &  cgl%cg%b(zdim,   is        :ie,             js        :je,             ks-3*dom%D_z:ke-3*dom%D_z)   )/cgl%cg%dz )
         end associate
         cgl => cgl%nxt
      enddo

   end subroutine divB_c6

!>
!! \brief Calculate 2nd order of divB for face-centered field and put it in the proper array.
!<

   subroutine divB_f2(ord)

      use cg_leaves,  only: leaves
      use cg_list,    only: cg_list_element
      use constants,  only: xdim, ydim, zdim
      use dataio_pub, only: die
      use domain,     only: dom
      use func,       only: operator(.equals.)

      implicit none

      integer, intent(in) :: ord

      type(cg_list_element), pointer :: cgl
      real, dimension(4), parameter :: coeff = [ 1., 9./8., 75./64., 1225./1024. ]
      integer :: o

      o = ord/2  ! BEWARE: tricky, assumed stencils on uniform grid
      if (o < lbound(coeff, dim=1) .or. o > ubound(coeff, dim=1) .or. 2*o /= ord .or. (coeff(o) .equals. 0.)) call die("[divB::divB_f2] cannot find coefficient") ! no odd order allowed here just in case

      cgl => leaves%first
      do while (associated(cgl))
         associate (is => cgl%cg%is, ie => cgl%cg%ie, &
              &     js => cgl%cg%js, je => cgl%cg%je, &
              &     ks => cgl%cg%ks, ke => cgl%cg%ke)
            cgl%cg%q(idivB)%arr(   is        :ie,         js        :je,         ks        :ke        ) = coeff(o) * ( &
                 & (cgl%cg%b(xdim, is+dom%D_x:ie+dom%D_x, js        :je,         ks        :ke        ) - &
                 &  cgl%cg%b(xdim, is        :ie,         js        :je,         ks        :ke        )   )/cgl%cg%dx + &
                 & (cgl%cg%b(ydim, is        :ie,         js+dom%D_y:je+dom%D_y, ks        :ke        ) - &
                 &  cgl%cg%b(ydim, is        :ie,         js        :je,         ks        :ke        )   )/cgl%cg%dy + &
                 & (cgl%cg%b(zdim, is        :ie,         js        :je,         ks+dom%D_z:ke+dom%D_z) - &
                 &  cgl%cg%b(zdim, is        :ie,         js        :je,         ks        :ke        )   )/cgl%cg%dz )
         end associate
         cgl => cgl%nxt
      enddo

   end subroutine divB_f2

!>
!! \brief Calculate 4th order correction for 2nd order estimate of divB for face-centered field and add it to the proper array.
!!
!! OPT: this may not be the best approach performance-wise.
!!
!! Spaghetti warning: this routime may be generalized and unified with divB_f2
!<

   subroutine divB_f4(ord)

      use cg_leaves,  only: leaves
      use cg_list,    only: cg_list_element
      use constants,  only: xdim, ydim, zdim
      use dataio_pub, only: die
      use domain,     only: dom
      use func,       only: operator(.equals.)

      implicit none

      integer, intent(in) :: ord

      type(cg_list_element), pointer :: cgl
      real, dimension(4), parameter :: coeff = [ 0., -1./24., -25./284., -245./3072. ]
      integer :: o

      o = ord/2  ! BEWARE: tricky, assumed stencils on uniform grid
      if (o < lbound(coeff, dim=1) .or. o > ubound(coeff, dim=1) .or. 2*o /= ord .or. (coeff(o) .equals. 0.)) call die("[divB::divB_f4] cannot find coefficient") ! no odd order allowed here just in case

      cgl => leaves%first
      do while (associated(cgl))
         associate (is => cgl%cg%is, ie => cgl%cg%ie, &
              &     js => cgl%cg%js, je => cgl%cg%je, &
              &     ks => cgl%cg%ks, ke => cgl%cg%ke)
            cgl%cg%q(idivB)%arr(     is          :ie,           js          :je,           ks          :ke          ) = &
                 cgl%cg%q(idivB)%arr(is          :ie,           js          :je,           ks          :ke          ) + coeff(o) * ( &
                 & (cgl%cg%b(xdim,   is+2*dom%D_x:ie+2*dom%D_x, js          :je,           ks          :ke          ) - &
                 &  cgl%cg%b(xdim,   is          :ie,           js          :je,           ks          :ke          )   )/cgl%cg%dx + &
                 & (cgl%cg%b(ydim,   is          :ie,           js+2*dom%D_y:je+2*dom%D_y, ks          :ke          ) - &
                 &  cgl%cg%b(ydim,   is          :ie,           js          :je,           ks          :ke          )   )/cgl%cg%dy + &
                 & (cgl%cg%b(zdim,   is          :ie,           js          :je,           ks+2*dom%D_z:ke+2*dom%D_z) - &
                 &  cgl%cg%b(zdim,   is          :ie,           js          :je,           ks          :ke          )   )/cgl%cg%dz )
         end associate
         cgl => cgl%nxt
      enddo

   end subroutine divB_f4

!>
!! \brief Calculate 4th order correction for 2nd order estimate of divB for face-centered field and add it to the proper array.
!!
!! OPT: this may not be the best approach performance-wise.
!!
!! Spaghetti warning: this routime may be generalized and unified with divB_f4
!<

   subroutine divB_f6(ord)

      use cg_leaves,  only: leaves
      use cg_list,    only: cg_list_element
      use constants,  only: xdim, ydim, zdim
      use dataio_pub, only: die
      use domain,     only: dom
      use func,       only: operator(.equals.)

      implicit none

      integer, intent(in) :: ord

      type(cg_list_element), pointer :: cgl
      real, dimension(4), parameter :: coeff = [ 0., 0., 3./640., 49./5120. ]
      integer :: o

      o = ord/2  ! BEWARE: tricky, assumed stencils on uniform grid
      if (o < lbound(coeff, dim=1) .or. o > ubound(coeff, dim=1) .or. 2*o /= ord .or. (coeff(o) .equals. 0.)) call die("[divB::divB_f6] cannot find coefficient") ! no odd order allowed here just in case

      cgl => leaves%first
      do while (associated(cgl))
         associate (is => cgl%cg%is, ie => cgl%cg%ie, &
              &     js => cgl%cg%js, je => cgl%cg%je, &
              &     ks => cgl%cg%ks, ke => cgl%cg%ke)
            cgl%cg%q(idivB)%arr(     is          :ie,           js          :je,           ks          :ke          ) = &
                 cgl%cg%q(idivB)%arr(is          :ie,           js          :je,           ks          :ke          ) + coeff(o) * ( &
                 & (cgl%cg%b(xdim,   is+3*dom%D_x:ie+3*dom%D_x, js          :je,           ks          :ke          ) - &
                 &  cgl%cg%b(xdim,   is          :ie,           js          :je,           ks          :ke          )   )/cgl%cg%dx + &
                 & (cgl%cg%b(ydim,   is          :ie,           js+3*dom%D_y:je+3*dom%D_y, ks          :ke          ) - &
                 &  cgl%cg%b(ydim,   is          :ie,           js          :je,           ks          :ke          )   )/cgl%cg%dy + &
                 & (cgl%cg%b(zdim,   is          :ie,           js          :je,           ks+3*dom%D_z:ke+3*dom%D_z) - &
                 &  cgl%cg%b(zdim,   is          :ie,           js          :je,           ks          :ke          )   )/cgl%cg%dz )
         end associate
         cgl => cgl%nxt
      enddo

   end subroutine divB_f6

end module div_B
