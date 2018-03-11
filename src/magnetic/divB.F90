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
      endif

      ccB = .false.  ! ToDo detect automatically
      if (present(cell_centered)) ccB = cell_centered

      call divB_init  ! it should be BOTH safe and cheap to call it multiple times
      
      if (ccB) then
         call divB_c2
         if (ord > 2) call warn("[divB::divB] only 2nd order of div(B) is currently implemented for cell-centered field.")
      else  ! face-centered B
         call divB_f2
         if (ord > 2) call warn("[divB::divB] only 2nd order of div(B) is currently implemented for face-centered field.")
      endif

   end subroutine divB

!>
!! \brief Calculate 2nd order of divB for cell-centered field and put it in the proper array.
!<

   subroutine divB_c2

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: xdim, ydim, zdim, half
      use domain,           only: dom

      implicit none

      type(cg_list_element), pointer :: cgl

      cgl => leaves%first
      do while (associated(cgl))
         associate (is => cgl%cg%is, ie => cgl%cg%ie, &
              &     js => cgl%cg%js, je => cgl%cg%je, &
              &     ks => cgl%cg%ks, ke => cgl%cg%ke)
            cgl%cg%q(idivB)%arr(   is        :ie,         js        :je,         ks        :ke        ) = half * ( &
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
!! \brief Calculate 2nd order of divB for cell-centered field and put it in the proper array.
!<

   subroutine divB_f2

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: xdim, ydim, zdim
      use domain,           only: dom

      implicit none

      type(cg_list_element), pointer :: cgl

      cgl => leaves%first
      do while (associated(cgl))
         associate (is => cgl%cg%is, ie => cgl%cg%ie, &
              &     js => cgl%cg%js, je => cgl%cg%je, &
              &     ks => cgl%cg%ks, ke => cgl%cg%ke)
            cgl%cg%q(idivB)%arr(   is        :ie,         js        :je,         ks        :ke        ) = ( &
                 & (cgl%cg%b(xdim, is+dom%D_x:ie+dom%D_x, js        :je,         ks        :ke        ) - &
                 &  cgl%cg%b(xdim, is        :ie        , js        :je,         ks        :ke        )   )/cgl%cg%dx + &
                 & (cgl%cg%b(ydim, is        :ie,         js+dom%D_y:je+dom%D_y, ks        :ke        ) - &
                 &  cgl%cg%b(ydim, is        :ie,         js        :je        , ks        :ke        )   )/cgl%cg%dy + &
                 & (cgl%cg%b(zdim, is        :ie,         js        :je,         ks+dom%D_z:ke+dom%D_z) - &
                 &  cgl%cg%b(zdim, is        :ie,         js        :je,         ks        :ke        )   )/cgl%cg%dz )
         end associate
         cgl => cgl%nxt
      enddo

   end subroutine divB_f2

end module div_B
